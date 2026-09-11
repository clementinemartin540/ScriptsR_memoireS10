# ============================================================================
#  calib_fit_cv.R  -  CALAGE avec VALIDATION CROISEE (jours x postes, par classe).
#
#  Cale (tau, fdry, beta) au KGE median journalier sur un TRAIN (sous-ensemble de
#  jours ET de postes, tire stratifie par classe), puis produit précipitations modélisées (P_mod) sur TOUS
#  les postes/jours avec ces parametres. L'evaluation train vs test est faite
#  dans diagnostics_cv.R a partir des objets sauvegardes ici.
#
#  Principe : le split agit UNIQUEMENT sur les observations montrees au calage
#  (P_obs masque). Le moteur optimize_tau_ortho_streaming n'est pas modifie : il
#  ne "voit" que le bloc train. c_j (daily) se calcule sur le domaine continental
#  -> identique train/test (projection physique, non ajustee sur obs) ; la calibration/validation
#  valide donc la transferabilite de (tau, fdry, beta).
#
#  Sorties (dans OUT_DIR, suffixe _cv<tag>) :
#    optimum_ortho_cv_<tag>.rds   : opt (parametres cales sur le train)
#    pmod_stations_cv_<tag>.rds   : P_mod/P_obs/P_era5 sur TOUS postes/jours
#    cv_split_<tag>.rds           : masques train/test (jours, postes) + reglages
# ============================================================================

source("C:/Users/cmartin/ownCloud/partage_LPM_20260824/calib_ortho/calib_cube_auxDEF.R")
library(ncdf4); library(sf); library(terra); library(lubridate); library(dplyr)

# ------------------------------------------------------------------- chemins
POSTES      <- "C:/Users/cmartin/ownCloud/scriptsR/calib/stations_MF.gpkg"
RR_CSV      <- "C:/Users/cmartin/ownCloud/scriptsR/QUOT/RR_data_light.csv"
ERA5_DAILY  <- "C:/Users/cmartin/ownCloud/daily/total_precip_daily_06utc_1940_2025.nc"
CLASSES_RDS <- "C:/Users/cmartin/ownCloud/scriptsR/classif_jours_tot.RDS"
CUBE_ROOT   <- "C:/Users/cmartin/ownCloud/partage_LPM_20260824/calib_ortho"

# ======= REGLAGES DU CUBE (doivent CORRESPONDRE a calib_build_cube.R) =======
OPT     <- "lcl"
SPATIAL <- c(rhoS = TRUE, Hw = TRUE, U = TRUE, V = TRUE, gamma = FALSE, zeta = TRUE) #sélection du cube lcl_no_gamma

# --------------------------------------------------- reglages propres au CALAGE
# PARAMS <- list(
#   tau_c = list(mode = "class", fix = NULL),
#   tau_f = list(link = "tau_c"),
#   fdry  = list(mode = "global", fix = NULL),
#   beta  = list(mode = "global",  fix = list("5" = 0)) # fix = c(NA,NA,NA,NA,0) possible aussi pour éteindre le lpm en classe 5
# )
PARAMS <- list(
  tau_c = list(mode = "class", fix = NULL), #ou "global"
  tau_f = list(link = "tau_c", fix = NULL),
  fdry  = list(mode = "global", fix = NULL),
  beta  = list(mode = "class",  fix = c(NA,NA,NA,NA,0)) 
)
C_SPEC    <- list(mode = "daily")   #coefficient c 
BETA_GRID <- seq(0.0, 1.0, by = 0.1)   #plage de valeurs balayée pour beta 
FDRY_GRID <- seq(0.0, 1.0, by = 0.1)   #plage de valeurs balayée pour fdry
CLAMP_C   <- c(-1, 1.5) #bornes de  c

# ============================ REGLAGES DE VALIDATION CROISEE =================
CV_SEED       <- 17        # graine du tirage (reproductibilite)
CV_FRAC_DAYS  <- 0.7      # fraction de JOURS  en apprentissage (par classe)
CV_FRAC_STA   <- 0.7      # fraction de POSTES en apprentissage
CV_SPLIT_STA  <- TRUE     # TRUE = split aussi les postes ; FALSE = tous en train
# ============================================================================

OUT_DIR <- path.expand(sprintf("%s/calib_cubes_ortho_%s_%s",    #charge les fichiers du cube : un fichier par année
                               CUBE_ROOT, OPT, spatial_tag(SPATIAL)))
files <- list.files(OUT_DIR, pattern = "^calib_cube_\\d+\\.rds$", full.names = TRUE)
stopifnot(length(files) > 0)
cat(sprintf("%d cubes annuels dans %s\n", length(files), OUT_DIR))

# ordre des postes : lu dans le 1er cube
sta_id <- readRDS(files[1])$sta_id
postes <- sf::st_read(POSTES, quiet = TRUE)
postes <- postes[match(sta_id, postes$idMF), ]

# --------- fournisseurs obs / ERA5 par annee (charges a la demande) ----------
RR <- read.table(RR_CSV, header = TRUE, sep = ";")
RR[RR < 0] <- NA
col_ord <- match(sta_id, names(RR)[-1]); stopifnot(!anyNA(col_ord))
RR_dates <- as.Date(RR$Date)

P_obs_f <- function(tsy) {       #permet de récupérer les données observées pour les jours et les postes demandés
  jc <- as.Date(format(tsy, "%Y-%m-%d"))
  m  <- match(jc, RR_dates); stopifnot(!anyNA(m))
  t(as.matrix(RR[, -1])[m, col_ord, drop = FALSE])
}

#extrait les précipitations journalières ERA5 aux coordonnées des postes :
P_era5_f <- function(tsy) read_era5_daily_at_stations(ERA5_DAILY, postes, tsy)  

# classification journaliere : TOUJOURS chargee (reporting par classe systematique ;
# la calibration/validation stratifie de toute facon par classe).
classe_all <- NULL
ts_all <- do.call(c, lapply(files, function(f) readRDS(f)$t_start))   #permet d'assigner le jour à sa classe correspondante
classes_obj <- readRDS(CLASSES_RDS)
jc <- as.Date(format(ts_all, "%Y-%m-%d"))
classe_all <- classes_obj$Classe[match(jc, as.Date(classes_obj$Date))]
if (anyNA(classe_all))
  stop(sprintf("%d jour(s) sans classe.", sum(is.na(classe_all))))   #arrêt si jour a pas de classe assignée
cat("Jours par classe : "); print(table(classe_all))

# ============================================================================
#  TIRAGE DU SPLIT train/test, stratifie par classe (jours) + tirage postes.
# ============================================================================
set.seed(CV_SEED)
n_days  <- length(ts_all)
grp_day <- if (is.null(classe_all)) rep(1L, n_days) else classe_all
date_all <- as.Date(format(ts_all, "%Y-%m-%d"))

# jours d'apprentissage : tires PAR CLASSE (meme proportion dans chaque classe)
is_train_day <- logical(n_days)
for (g in unique(grp_day)) {
  idx <- which(grp_day == g)
  is_train_day[sample(idx, floor(CV_FRAC_DAYS * length(idx)))] <- TRUE
}
cat(sprintf("Jours  : %d train / %d test\n", sum(is_train_day), sum(!is_train_day)))

# postes d'apprentissage (optionnel)
nS <- length(sta_id)
is_train_sta <- rep(TRUE, nS)
if (CV_SPLIT_STA) {
  is_train_sta[sample(seq_len(nS), round((1 - CV_FRAC_STA) * nS))] <- FALSE
  cat(sprintf("Postes : %d train / %d test\n", sum(is_train_sta), sum(!is_train_sta)))
}

# fournisseur obs MASQUE : ne montre au calage que le bloc TRAIN x TRAIN, remplace les postes et jours tests par NA
mask_train <- function(M, tsy) {
  jc <- as.Date(format(tsy, "%Y-%m-%d"))
  kd <- is_train_day[match(jc, date_all)]        # jours train de cette annee
  M[!is_train_sta, ] <- NA                        # postes test  -> NA
  M[, !kd]           <- NA                        # jours  test  -> NA
  M
}
P_obs_train_f <- function(tsy) mask_train(P_obs_f(tsy), tsy) #c'est ce qui passe au calage, ne contient que les données train

# ------------------------------------------------------------------- calage
# NB : classe_all reste l'ensemble des jours ; le masquage des obs suffit a
# retirer les jours-test du KGE (kge_daily_median ignore les colonnes tout-NA).
opt <- optimize_tau_ortho_streaming(files, P_obs_train_f, P_era5_f, classe = classe_all,
                                    params = PARAMS, c_spec = C_SPEC,
                                    beta_grid = BETA_GRID, fdry_grid = FDRY_GRID,
                                    clamp_c = CLAMP_C)

# reference ERA5 seul, evaluee sur le MEME bloc TRAIN que le modele (obs masquees
# par mask_train) -> comparaison a perimetre egal. KGE journalier + par classe.
kj_era5 <- c()   #calcul kge journalier
for (f in files) {
  tsy <- readRDS(f)$t_start
  Po_tr <- mask_train(P_obs_f(tsy), tsy)             # obs du train uniquement
  cc <- kge_daily_median(pmax(P_era5_f(tsy), 0), Po_tr, sd_min = 0.5, n_min = 20)
  kj_era5 <- c(kj_era5, cc$kge_j)
}
kge_era5_train <- median(kj_era5, na.rm = TRUE)    #calcul kge par classe
era5_by_class <- NULL
if (!is.null(classe_all)) {
  era5_by_class <- setNames(lapply(sort(unique(classe_all)), function(k)
    median(kj_era5[classe_all == k], na.rm = TRUE)), as.character(sort(unique(classe_all))))
}

cat(sprintf("\nOPTIMUM CV (c_mode=%s) cale sur le TRAIN :\n", opt$c_mode))
for (k in names(opt$per_class)) {
  p <- opt$per_class[[k]]
  cat(sprintf("  classe %s (n=%d) : tau_c=%.0f tau_f=%.0f fdry=%.2f beta=%.3f\n",
              k, p$n, p$tau_c, p$tau_f, p$fdry, p$beta))
}

# --- resume tabulaire des parametres cales sur le TRAIN (modele vs ERA5) ---
param_table <- do.call(rbind, lapply(names(opt$per_class), function(k) {
  p <- opt$per_class[[k]]
  ke <- if (!is.null(era5_by_class)) era5_by_class[[k]] else NA_real_
  data.frame(classe = k, n = p$n,
             tau_c = round(p$tau_c), tau_f = round(p$tau_f),
             fdry = round(p$fdry, 2), beta = round(p$beta, 3),
             kge_mod = round(p$kge, 4), kge_era5 = round(ke, 4),
             gain = round(p$kge - ke, 4))
}))
cat("\n=== RESUME DES PARAMETRES CALES (TRAIN, modele vs ERA5) ===\n")
print(param_table, row.names = FALSE)
cat(sprintf("KGE global TRAIN : modele=%.4f  ERA5=%.4f\n", opt$comp$kge, kge_era5_train))
opt$param_table <- param_table
cat("(scores sur le TRAIN ; evaluation train/test complete via diagnostics_cv.R)\n")

# ---- production de P_mod sur TOUS les postes/jours (params du train) ----
# P_obs COMPLET ici : on veut pouvoir evaluer sur train ET test ensuite.
pmod <- produce_pmod_stations(opt, files, P_obs_f, P_era5_f)

# ------------------------------------------------------------------- sauvegardes
tag <- make_calib_tag(PARAMS, c_mode = C_SPEC$mode, seed = CV_SEED, clamp_c = CLAMP_C)

saveRDS(opt,  file.path(OUT_DIR, sprintf("optimum_ortho_cv_%s.rds", tag)))
saveRDS(pmod, file.path(OUT_DIR, sprintf("pmod_stations_cv_%s.rds", tag)))

# masques et reglages du split : tout ce qu'il faut pour reproduire et evaluer train/test dans diagnostic_cv.
cv_split <- list(
  is_train_day = is_train_day,      # [jour] dans l'ordre de ts_all
  is_train_sta = is_train_sta,      # [poste] dans l'ordre de sta_id
  date_all = date_all, ts_all = ts_all, sta_id = sta_id,
  classe_all = classe_all,
  frac_days = CV_FRAC_DAYS, frac_sta = CV_FRAC_STA,
  split_sta = CV_SPLIT_STA, seed = CV_SEED)
saveRDS(cv_split, file.path(OUT_DIR, sprintf("cv_split_%s.rds", tag)))

cat(sprintf("\nSauvegarde : optimum_ortho_cv_%s.rds, pmod_stations_cv_%s.rds, cv_split_%s.rds\n",
            tag, tag, tag))
cat(sprintf("=> lancer diagnostics_cv.R pour l'evaluation train/test et les figures,\n avec TAG <- \"%s\"",tag))
