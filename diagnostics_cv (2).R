# ============================================================================
#  diagnostics_cv.R  -  EVALUATION train/test et FIGURES (validation croisee).
#
#  Charge les objets ecrits par calib_fit_cv.R (opt, pmod, cv_split) et evalue
#  les performances SEPAREMENT sur les 4 blocs du split jours x postes :
#     train      : jours train  & postes train   (a servi au calage)
#     test_day   : jours test   & postes train   (generalisation temporelle)
#     test_sta   : jours train  & postes test    (generalisation spatiale)
#     test_both  : jours test   & postes test    (validation stricte)
#  Puis figures : KGE journalier par bloc, KGE spatial par echelle x bloc,
#  composantes par echelle x bloc.
#
#  Ne recale rien : purement lecture + traces (rejouable a volonte).
# ============================================================================

source("C:/Users/cmartin/ownCloud/partage_LPM_20260824/calib_ortho/calib_cube_auxDEF.R")
library(ggplot2); library(patchwork); library(dplyr); library(tidyr); library(ggh4x)

CUBE_ROOT <- "C:/Users/cmartin/ownCloud/partage_LPM_20260824/calib_ortho"
OPT       <- "lcl"
SPATIAL   <- c(rhoS = TRUE, Hw = TRUE, U = TRUE, V = TRUE, gamma = FALSE, zeta = TRUE)
TAG       <- "tau-cls-link_fdry-glob_beta-cls-fix5-0_c-daily_clamp-1.1.5_seed1"   # suffixe des fichiers _cv

OUT_DIR <- path.expand(sprintf("%s/calib_cubes_ortho_%s_%s",
                               CUBE_ROOT, OPT, spatial_tag(SPATIAL)))
opt      <- readRDS(file.path(OUT_DIR, sprintf("optimum_ortho_cv_%s.rds", TAG)))
pmod     <- readRDS(file.path(OUT_DIR, sprintf("pmod_stations_cv_%s.rds", TAG)))
cv       <- readRDS(file.path(OUT_DIR, sprintf("cv_split_%s.rds", TAG)))

P_mod <- pmod$P_mod; P_obs <- pmod$P_obs; P_era5 <- pmod$P_era5
ts    <- pmod$t_start

era5color <- "lightsalmon"
modelcolor <- "steelblue"

# ----------------------------------------------------------------------------
#  Masques de bloc [poste, jour], alignes sur pmod.
# ----------------------------------------------------------------------------
jc_pmod <- as.Date(format(ts, "%Y-%m-%d"))
kd <- cv$is_train_day[match(jc_pmod, cv$date_all)]     # jour train ? (aligne pmod)
st <- cv$is_train_sta                                  # poste train ? (ordre sta_id)
stopifnot(length(st) == nrow(P_mod), length(kd) == ncol(P_mod))

cell <- list(
  train     = outer( st,  kd),
  test_day  = outer( st, !kd),
  test_sta  = outer(!st,  kd),
  test_both = outer(!st, !kd))

# ============================================================================
#  1. KGE MEDIAN JOURNALIER par bloc (modele vs ERA5).
# ============================================================================
kge_block <- function(Pm, Po, keep_cell, sd_min = 0.5, n_min = 20) {
  Po2 <- Po; Po2[!keep_cell] <- NA
  kge_daily_median(pmax(Pm, 0), Po2, sd_min = sd_min, n_min = n_min)
}

cat("\n--- KGE median journalier par bloc (modele / ERA5) ---\n")
tab_block <- bind_rows(lapply(names(cell), function(nm) {
  cm <- kge_block(P_mod,  P_obs, cell[[nm]])
  ce <- kge_block(P_era5, P_obs, cell[[nm]])
  cat(sprintf("  %-9s : modele KGE=%.4f (r=%.3f a=%.3f b=%.3f) | ERA5 KGE=%.4f  [n_j=%d]\n",
              nm, cm$kge, cm$r, cm$alpha, cm$beta, ce$kge, cm$n_used))
  data.frame(bloc = nm, source = c("modele","ERA5"),
             kge = c(cm$kge, ce$kge), r = c(cm$r, ce$r),
             alpha = c(cm$alpha, ce$alpha), beta = c(cm$beta, ce$beta))
}))
tab_block$bloc <- factor(tab_block$bloc, levels = names(cell))

# barres : KGE par bloc, modele vs ERA5 (train a gauche, tests a droite)
gA <- ggplot(tab_block, aes(bloc, kge, fill = source)) +
  geom_col(position = position_dodge(0.8), width = 0.7, alpha = 0.85) +
  scale_fill_manual(values = c(ERA5 = era5color, modele = modelcolor)) +
  labs(title = sprintf("KGE median journalier par bloc de validation croisee (%s)", TAG),
       subtitle = "train = calage ; test_day/test_sta/test_both = jamais vus au calage",
       x = NULL, y = "KGE median journalier") +
  theme_minimal(base_size = 10)
print(gA)

# ============================================================================
#  2. KGE SPATIAL par ECHELLE d'agregation x BLOC.
#     Pour chaque bloc, on masque les cellules hors-bloc puis on agrege.
# ============================================================================
keys <- list(
  journalier = format(ts, "%Y-%m-%d"),
  mensuel    = format(ts, "%Y-%m"),
  saisonnier = paste0(format(ts, "%Y"), "-S", (as.integer(format(ts,"%m"))-1)%/%3 + 1),
  annuel     = format(ts, "%Y"))

# agrege [poste,jour]->[poste,periode] en sommant les jours ou l'obs est presente
# ET la cellule appartient au bloc (masque combine).
agg_masked_block <- function(M, obs, key, keep_cell) {
  M2 <- M; M2[!(is.finite(obs) & keep_cell)] <- NA
  t(apply(M2, 1, function(x) tapply(x, key, function(v)
    if (all(is.na(v))) NA_real_ else sum(v, na.rm = TRUE))))
}

# KGE spatial + composantes apparie mod/ERA5, par periode (comme workflow initial)
kge_paired_by_period <- function(Am, Ae, Ao, sd_min = 0.5, n_min = 20) {
  np <- ncol(Ao)
  out <- data.frame(kge_mod=NA_real_, r_mod=NA_real_, alpha_mod=NA_real_, beta_mod=NA_real_,
                    kge_era5=NA_real_, r_era5=NA_real_, alpha_era5=NA_real_, beta_era5=NA_real_)[rep(1,np),]
  comp <- function(sim, o, ok) {
    m <- sim[ok]; oo <- o[ok]; so <- sd(oo)
    if (!is.finite(so) || so < sd_min) return(rep(NA_real_,4))
    r <- suppressWarnings(cor(m, oo)); a <- sd(m)/so; b <- mean(m)/mean(oo)
    if (!is.finite(r) || !is.finite(a) || !is.finite(b)) return(rep(NA_real_,4))
    c(1 - sqrt((r-1)^2 + (a-1)^2 + (b-1)^2), r, a, b)
  }
  for (t in seq_len(np)) {
    o <- Ao[, t]; ok <- is.finite(o) & is.finite(Am[, t]) & is.finite(Ae[, t])
    if (sum(ok) < n_min) next
    out[t, ] <- c(comp(Am[, t], o, ok), comp(Ae[, t], o, ok))
  }
  out[is.finite(out$kge_mod) & is.finite(out$kge_era5), ]
}

# boucle blocs x echelles -> data.frame apparie avec composantes
paired <- bind_rows(lapply(names(cell), function(nm) {
  kc <- cell[[nm]]
  bind_rows(lapply(names(keys), function(sc) {
    key <- keys[[sc]]
    d <- kge_paired_by_period(agg_masked_block(P_mod,  P_obs, key, kc),
                              agg_masked_block(P_era5, P_obs, key, kc),
                              agg_masked_block(P_obs,  P_obs, key, kc))
    if (nrow(d) == 0) return(NULL)
    data.frame(bloc = nm, echelle = sc, d)
  }))
}))
paired$bloc    <- factor(paired$bloc,    levels = names(cell))
paired$echelle <- factor(paired$echelle, levels = names(keys))

# --- KGE spatial : boxplots, une facette par bloc, ERA5 vs modele par echelle ---
df_long <- paired |>
  pivot_longer(c(kge_mod, kge_era5), names_to = "source", values_to = "kge") |>
  mutate(source = factor(ifelse(source == "kge_mod", "modele", "ERA5"),
                         levels = c("ERA5", "modele")))

gB <- ggplot(df_long, aes(echelle, kge, fill = source)) +
  geom_boxplot(outlier.shape = NA, position = position_dodge(0.8), alpha = 0.7) +
  facet_wrap(~bloc, nrow = 1) +
  scale_fill_manual(values = c(ERA5 = era5color, modele = modelcolor)) +
  labs(title = "KGE spatial par echelle d'agregation, par bloc de validation",
       x = NULL, y = "KGE spatial") +
  theme_minimal(base_size = 10) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
  coord_cartesian(ylim = c(-0.5, 1))
print(gB)

# --- difference appariee (modele - ERA5) par echelle, une facette par bloc ---
paired <- paired |> mutate(diff = kge_mod - kge_era5)
frac_win <- paired |> group_by(bloc, echelle) |>
  summarise(f = mean(diff > 0), med = median(diff), .groups = "drop")
cat("\n--- fraction de periodes ou le modele bat ERA5 (par bloc x echelle) ---\n")
print(as.data.frame(frac_win))

gC <- ggplot(paired, aes(echelle, diff)) +
  geom_hline(yintercept = 0, linetype = 2, colour = "grey50") +
  geom_boxplot(fill = "#9467bd", outlier.shape = NA, alpha = 0.6, width = 0.5) +
  facet_wrap(~bloc, nrow = 1) +
  labs(title = "Difference appariee (modele - ERA5) par bloc",
       subtitle = "au-dessus de 0 = le LPM ameliore", x = NULL, y = "delta KGE") +
  theme_minimal(base_size = 10) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
  coord_cartesian(ylim = c(-0.5, 0.5))
print(gC)

# ============================================================================
#  3. COMPOSANTES (kge/r/alpha/beta) par echelle, focalise sur test_both.
#     (le bloc le plus severe ; changer 'bloc_focus' pour un autre)
# ============================================================================
bloc_focus <- "test_both"
comp_long <- paired |>
  filter(bloc == bloc_focus) |>
  select(echelle, kge_mod, r_mod, alpha_mod, beta_mod,
         kge_era5, r_era5, alpha_era5, beta_era5) |>
  pivot_longer(-echelle, names_to = "var", values_to = "val") |>
  separate(var, into = c("comp", "source"), sep = "_", extra = "merge") |>
  mutate(source = factor(ifelse(source == "mod", "modele", "ERA5"),
                         levels = c("ERA5", "modele")),
         comp = factor(comp, levels = c("kge", "r", "alpha", "beta")))

gD <- ggplot(comp_long, aes(echelle, val, fill = source)) +
  geom_boxplot(outlier.shape = NA, position = position_dodge(0.8), alpha = 0.7) +
  geom_hline(data = data.frame(comp = factor(c("r","alpha","beta"),
                                             levels = c("kge","r","alpha","beta"))),
             aes(yintercept = 1), linetype = 3, colour = "grey60") +
  facet_wrap(~comp, scales = "free_y") +
  facetted_pos_scales(y = list(
    scale_y_continuous(limits = c(-0.5, 1)),     # kge
    scale_y_continuous(limits = c(0, 1)),        # r
    scale_y_continuous(limits = c(0, 1.5)),      # alpha
    scale_y_continuous(limits = c(0.5, 1.5)))) + # beta
  scale_fill_manual(values = c(ERA5 = era5color, modele = modelcolor)) +
  labs(title = sprintf("KGE et composantes par echelle - bloc %s", bloc_focus),
       x = NULL, y = NULL) +
  theme_minimal(base_size = 10) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
print(gD)

# ============================================================================
#  4. KGE *TEMPOREL* AUX STATIONS, par echelle d'agregation, sur un bloc CV.
#     Pour chaque station et chaque echelle, on calcule le KGE de sa SERIE
#     temporelle (a travers les periodes), restreinte aux cellules du bloc.
#     Distribution sur l'echantillon de stations. Defaut : bloc test_both.
# ============================================================================
bloc_temp <- "test_sta"
kc_t <- cell[[bloc_temp]]

# construit, pour un bloc et une echelle, les series [poste, periode] agregees
# sur les jours du bloc ; njours = nb de jours agreges par cellule [poste,periode].
series_block <- function(M, obs, key, keep_cell) {
  M2 <- M; keep <- is.finite(obs) & keep_cell
  M2[!keep] <- NA
  A <- t(apply(M2, 1, function(x) tapply(x, key, function(v)
    if (all(is.na(v))) NA_real_ else sum(v, na.rm = TRUE))))
  N <- t(apply(keep, 1, function(x) tapply(x, key, sum)))   # nb jours par cellule
  list(A = A, N = N)
}

# KGE temporel d'une station : sa serie sim vs obs a travers les periodes.
kge_temporal <- function(sim, obs, n_min = 6, sd_min = 0.5) {
  ok <- is.finite(sim) & is.finite(obs)
  if (sum(ok) < n_min) return(rep(NA_real_, 4))
  s <- sim[ok]; o <- obs[ok]; so <- sd(o)
  if (!is.finite(so) || so < sd_min) return(rep(NA_real_, 4))
  r <- suppressWarnings(cor(s, o)); a <- sd(s)/so; b <- mean(s)/mean(o)
  if (!is.finite(r) || !is.finite(a) || !is.finite(b)) return(rep(NA_real_, 4))
  c(1 - sqrt((r-1)^2 + (a-1)^2 + (b-1)^2), r, a, b)
}

# nb minimal de jours pour qu'un cumul de periode compte (evite cumuls trop
# partiels : au bloc test, une periode n'agrege que ~30% des jours) ; et nb
# minimal de periodes pour un KGE temporel stable.
NJ_MIN  <- c(journalier = 1, mensuel = 8, saisonnier = 20, annuel = 60)
NP_MIN  <- c(journalier = 60, mensuel = 6, saisonnier = 4, annuel = 3)

temp <- bind_rows(lapply(names(keys), function(sc) {
  key <- keys[[sc]]
  sm <- series_block(P_mod,  P_obs, key, kc_t)
  se <- series_block(P_era5, P_obs, key, kc_t)
  so <- series_block(P_obs,  P_obs, key, kc_t)
  # invalide les cellules a trop peu de jours (cumul de periode partiel)
  bad <- so$N < NJ_MIN[[sc]]
  Ao <- so$A; Am <- sm$A; Ae <- se$A
  Ao[bad] <- NA; Am[bad] <- NA; Ae[bad] <- NA
  np <- nrow(Ao)
  rows <- lapply(seq_len(np), function(i) {
    cm <- kge_temporal(Am[i, ], Ao[i, ], n_min = NP_MIN[[sc]])
    ce <- kge_temporal(Ae[i, ], Ao[i, ], n_min = NP_MIN[[sc]])
    c(cm, ce)
  })
  M <- do.call(rbind, rows)
  d <- data.frame(echelle = sc, poste = seq_len(np),
                  kge_mod = M[,1], r_mod = M[,2], alpha_mod = M[,3], beta_mod = M[,4],
                  kge_era5 = M[,5], r_era5 = M[,6], alpha_era5 = M[,7], beta_era5 = M[,8])
  d[is.finite(d$kge_mod) & is.finite(d$kge_era5), ]
}))
temp$echelle <- factor(temp$echelle, levels = names(keys))

cat(sprintf("\n--- KGE temporel median aux stations, bloc %s (modele / ERA5) ---\n", bloc_temp))
temp_synth <- temp |> group_by(echelle) |>
  summarise(kge_mod = median(kge_mod), kge_era5 = median(kge_era5),
            alpha_mod = median(alpha_mod), alpha_era5 = median(alpha_era5),
            r_mod = median(r_mod), r_era5 = median(r_era5),
            beta_mod = median(beta_mod), beta_era5 = median(beta_era5),
            n = n(), .groups = "drop")
print(as.data.frame(temp_synth))

# --- KGE temporel : boxplots ERA5 vs modele, par echelle ---
temp_kge_long <- temp |>
  pivot_longer(c(kge_mod, kge_era5), names_to = "source", values_to = "kge") |>
  mutate(source = factor(ifelse(source == "kge_mod", "modele", "ERA5"),
                         levels = c("ERA5", "modele")))

gE <- ggplot(temp_kge_long, aes(echelle, kge, fill = source)) +
  geom_boxplot(outlier.shape = NA, position = position_dodge(0.8), alpha = 0.7) +
  scale_fill_manual(values = c(ERA5 = era5color, modele = modelcolor)) +
  labs(title = sprintf("KGE temporel aux stations par echelle - bloc %s", bloc_temp),
       subtitle = "distribution sur les stations ; KGE de la serie temporelle de chaque poste",
       x = NULL, y = "KGE temporel") +
  theme_minimal(base_size = 10) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
  coord_cartesian(ylim = c(-0.5, 1))
print(gE)

# --- composantes temporelles (kge/r/alpha/beta) par echelle ---
temp_comp_long <- temp |>
  select(echelle, kge_mod, r_mod, alpha_mod, beta_mod,
         kge_era5, r_era5, alpha_era5, beta_era5) |>
  pivot_longer(-echelle, names_to = "var", values_to = "val") |>
  separate(var, into = c("comp", "source"), sep = "_", extra = "merge") |>
  mutate(source = factor(ifelse(source == "mod", "modele", "ERA5"),
                         levels = c("ERA5", "modele")),
         comp = factor(comp, levels = c("kge", "r", "alpha", "beta")))

gF <- ggplot(temp_comp_long, aes(echelle, val, fill = source)) +
  geom_boxplot(outlier.shape = NA, position = position_dodge(0.8), alpha = 0.7) +
  geom_hline(data = data.frame(comp = factor(c("r","alpha","beta"),
                                             levels = c("kge","r","alpha","beta"))),
             aes(yintercept = 1), linetype = 3, colour = "grey60") +
  facet_wrap(~comp, scales = "free_y") +
  facetted_pos_scales(y = list(
    scale_y_continuous(limits = c(-0.5, 1)),     # kge
    scale_y_continuous(limits = c(0, 1)),        # r
    scale_y_continuous(limits = c(0, 2)),        # alpha
    scale_y_continuous(limits = c(0.5, 1.5)))) + # beta
  scale_fill_manual(values = c(ERA5 = era5color, modele = modelcolor)) +
  labs(title = sprintf("KGE temporel et composantes par echelle - bloc %s", bloc_temp),
       x = NULL, y = NULL) +
  theme_minimal(base_size = 10) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
print(gF)