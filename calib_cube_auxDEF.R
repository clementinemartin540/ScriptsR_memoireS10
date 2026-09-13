# ============================================================================
#  calib_cube_aux.R  -  Fonctions auxiliaires du calage GWALARN.
#  Charge par source("calib_cube_aux.R") depuis le driver.
# ============================================================================
library(terra)
library(ncdf4)
library(sf)
library(moistR)
# NB : make_coarse_buffer et locate_timestep sont fournis par le package moistR
# (avec l'argument opt pour la variante sfc/lcl). Ne PAS les redefinir ici : une
# redefinition sourcee masquerait la version du package dans .GlobalEnv.

# ---------------------------------------------------------------------------
#  Masque spectral PASSE-HAUT pour l'increment orographique.
#  Phi = exp(-rho^2 / k_c^2) est un passe-bas gaussien doux (k_c = 2*pi/Lc) ;
#  le passe-haut est hp = 1 - Phi. Applique a FFT_elev, il retire la topographie
#  de grande echelle (celle qu'ERA5 resout deja), ne gardant que le relief fin
#  sub-Lc -> LPM(h_fin - h_coarse). Lc voisin de la maille ERA5 (~28 km).
#  Retourne un vecteur de la longueur de pre$rho2 (grille spectrale aplatie).
# ---------------------------------------------------------------------------
make_hp_mask <- function(pre, Lc) {
  k_c <- 2 * pi / Lc
  phi <- exp(-(pre$rho2) / k_c^2)     # passe-bas gaussien
  1 - phi                              # passe-haut
}

# grille log de n longueurs de coupure autour d'une valeur centrale (defaut ERA5 ~28 km)
logspace_Lc <- function(Lc_min, Lc_max, n) {
  exp(seq(log(Lc_min), log(Lc_max), length.out = n))
}

# ---------------------------------------------------------------------------
#  KGE (Kling-Gupta Efficiency, Gupta et al. 2009).
#    KGE = 1 - sqrt( (r-1)^2 + (alpha-1)^2 + (beta-1)^2 )
#  r      = correlation de Pearson (mod, obs)
#  alpha  = sd(mod) / sd(obs)      (rapport de variabilite)
#  beta   = mean(mod) / mean(obs)  (rapport de biais)
#  KGE=1 parfait ; on MAXIMISE. NA-robuste (paires completes seulement).
#  `transform` : "identity" (defaut, sur les valeurs brutes en mm), "sqrt" ou
#  "log1p" pour attenuer le poids des forts cumuls.
#  Renvoie une liste (kge, r, alpha, beta) pour permettre le diagnostic des
#  composantes ; kge$kge est le scalaire a optimiser.
# ---------------------------------------------------------------------------
kge <- function(mod, obs, transform = c("identity", "sqrt", "log1p")) {
  transform <- match.arg(transform)
  tf <- switch(transform, identity = function(x) x,
               sqrt = function(x) sqrt(pmax(x, 0)),
               log1p = function(x) log1p(pmax(x, 0)))
  m <- as.vector(mod); o <- as.vector(obs)
  ok <- is.finite(m) & is.finite(o)
  m <- tf(m[ok]); o <- tf(o[ok])
  if (length(m) < 3 || sd(o) == 0) return(list(kge = NA_real_, r = NA, alpha = NA, beta = NA))
  r     <- suppressWarnings(cor(m, o))
  alpha <- sd(m) / sd(o)
  beta  <- mean(m) / mean(o)
  k <- 1 - sqrt((r - 1)^2 + (alpha - 1)^2 + (beta - 1)^2)
  list(kge = k, r = r, alpha = alpha, beta = beta)
}

# ---------------------------------------------------------------------------
#  KGE SPATIAL JOURNALIER puis MEDIANE sur les jours.
#  Pour chaque jour, on calcule le KGE spatial (sur les postes) du patron de
#  pluie de CE jour (mod vs obs), puis on prend la mediane des KGE journaliers.
#  Cette metrique evalue la reproduction de la STRUCTURE SPATIALE jour apres
#  jour, sans que les jours de forte pluie (grande variance) ecrasent les autres.
#
#  mod, obs : [poste, jour]. Pluie BRUTE (pas de transformation : on veut le
#  patron spatial reel en mm). Un jour est exclu si l'ecart-type spatial des obs
#  est < sd_min (pluie quasi uniforme -> KGE indefini) ou si < n_min postes
#  valides. Renvoie la mediane du KGE, les medianes des composantes (r, alpha,
#  beta) et le vecteur des KGE journaliers (diagnostic de distribution).
# ---------------------------------------------------------------------------
kge_daily_median <- function(mod, obs, sd_min = 0.5, n_min = 20) {
  nD <- ncol(obs)
  kge_j <- rep(NA_real_, nD); r_j <- kge_j; a_j <- kge_j; b_j <- kge_j
  for (t in seq_len(nD)) {
    m <- mod[, t]; o <- obs[, t]
    ok <- is.finite(m) & is.finite(o)
    if (sum(ok) < n_min) next
    m <- m[ok]; o <- o[ok]
    so <- sd(o)
    if (!is.finite(so) || so < sd_min) next        # jour spatialement quasi uniforme
    r  <- suppressWarnings(cor(m, o))
    a  <- sd(m) / so
    b  <- mean(m) / mean(o)
    if (!is.finite(r) || !is.finite(a) || !is.finite(b)) next
    r_j[t] <- r; a_j[t] <- a; b_j[t] <- b
    kge_j[t] <- 1 - sqrt((r - 1)^2 + (a - 1)^2 + (b - 1)^2)
  }
  list(kge = median(kge_j, na.rm = TRUE),
       r = median(r_j, na.rm = TRUE),
       alpha = median(a_j, na.rm = TRUE),
       beta = median(b_j, na.rm = TRUE),
       kge_j = kge_j, r_j = r_j, alpha_j = a_j, beta_j = b_j,
       n_used = sum(is.finite(kge_j)))
}

# ---------------------------------------------------------------------------
#  OPTIMISATION (tau_c, tau_f) au KGE median, FOND ERA5 + LPM ORTHOGONALISE.
#  Le LPM complet est orthogonalise vis-a-vis d'ERA5 SUR LE CHAMP (resolution
#  ERA5), pour retirer le double comptage grande echelle qui tirait les tau hors
#  des valeurs canoniques :
#     c(tau) = sum_j s_LB[j,tau] / sum_j s_BB[j]      (global, ou par classe)
#     P_LPM_perp = P_LPM - c * P_era5_postes          (aux postes)
#     P_mod = max(P_era5 + beta * P_LPM_perp, 0)
#  beta cale par balayage 1D au KGE median (beta=1 possible si beta_grid=1).
#
#  cube    : doit contenir s_LB, s_BB (calibrate_cube avec bg_daily).
#  P_obs   : [poste, jour] mm. P_era5 : [poste, jour] mm (pluie ERA5 aux postes).
#  classe  : NULL (c global) ou vecteur [jour] (c par classe).
#  beta_grid : grille de beta (defaut log 1e-2..3) ; mettre 1 pour beta fige.
# ---------------------------------------------------------------------------
#  c_mode : "global" (un seul c sur tous les jours), "class" (un c par classe,
#           via `classe`), ou "daily" (c_j par jour, adaptatif). En mode daily,
#           c_j n'est PAS un parametre cale mais la projection du LPM du jour sur
#           ERA5 du jour -> pas de ddl supplementaire, transferable (formule
#           applicable a tout jour futur). Garde-fou : c_j = 0 si s_BB[j] est
#           sous `see_frac` * median(s_BB) (jour trop sec -> rien a orthogonaliser).
#  fix_tau_c / fix_tau_f : si fournis (en s), le balayage est restreint a
#  l'indice de grille le plus proche de la valeur demandee (le tau est alors
#  FIXE, pas optimise). NULL = tau balaye normalement. Permet, sur un cube deja
#  calcule, de fixer l'un et/ou l'autre sans recalcul (test de sensibilite).
#  La valeur reellement utilisee (point de grille le plus proche) est renvoyee.
# ---------------------------------------------------------------------------
#  fix_tau_c / fix_tau_f / fix_fdry : si fournis, le balayage de l'axe concerne
#  est restreint a l'indice de grille le plus proche (parametre FIXE, pas
#  optimise). NULL = axe balaye normalement. Permet, sur un cube deja calcule,
#  de fixer un ou plusieurs axes sans recalcul. Valeur reellement utilisee (point
#  de grille le plus proche) renvoyee. fdry attenue l'assechement du LPM (axe du
#  cube car il change s_LB donc c).
# ---------------------------------------------------------------------------
optimize_tau_kge_ortho_era5 <- function(cube, P_obs, P_era5, classe = NULL,
                                        c_mode = c("global", "class", "daily"),
                                        beta_grid = exp(seq(log(1e-2), log(3), length.out = 40)),
                                        sd_min = 0.5, n_min = 20, see_frac = 0.05,
                                        fix_tau_c = NULL, fix_tau_f = NULL, fix_fdry = NULL,
                                        clamp_c = c(0, 1),
                                        fdry_grid = seq(0.4, 1.0, by = 0.1)) {
  c_mode <- match.arg(c_mode)
  clampc <- function(x) if (is.null(clamp_c)) x else pmin(pmax(x, clamp_c[1]), clamp_c[2])
  stopifnot(!is.null(cube$s_LB_pos), !is.null(cube$s_BB))
  d <- dim(cube$P_LPM_pos)                          # [poste, jour, tc, tf] (design pos/neg)
  nC <- d[3]; nF <- d[4]
  FDRY <- fdry_grid                                 # grille fdry BALAYEE (reconstruction lineaire)
  nX <- length(FDRY)
  nDj <- length(cube$t_start)
  if (is.null(classe)) classe <- rep(1L, nDj)
  classes <- sort(unique(classe))
  see_floor <- see_frac * median(cube$s_BB, na.rm = TRUE)

  # indices a explorer : tous, ou le plus proche si l'axe est fixe
  ic_set <- if (is.null(fix_tau_c)) seq_len(nC) else which.min(abs(cube$TAU_C - fix_tau_c))
  if_set <- if (is.null(fix_tau_f)) seq_len(nF) else which.min(abs(cube$TAU_F - fix_tau_f))
  ix_set <- if (is.null(fix_fdry))  seq_len(nX) else which.min(abs(FDRY - fix_fdry))
  if (!is.null(fix_tau_c))
    message(sprintf("tau_c fixe a %.0f s (demande %.0f s)", cube$TAU_C[ic_set], fix_tau_c))
  if (!is.null(fix_tau_f))
    message(sprintf("tau_f fixe a %.0f s (demande %.0f s)", cube$TAU_F[if_set], fix_tau_f))
  if (!is.null(fix_fdry))
    message(sprintf("fdry fixe a %.2f (demande %.2f)", FDRY[ix_set], fix_fdry))

  # reconstruction lineaire fdry : s_LB(fdry) et P_LPM(fdry) depuis pos/neg
  sLE_at <- function(j, ic, iff, ix)
    cube$s_LB_pos[j, ic, iff] + FDRY[ix] * cube$s_LB_neg[j, ic, iff]
  PLPM_at <- function(ic, iff, ix)
    cube$P_LPM_pos[, , ic, iff] + FDRY[ix] * cube$P_LPM_neg[, , ic, iff]

  # construit Pperp[poste,jour] pour (ic,iff,ix) selon le mode ; PL = P_LPM attenue
  make_perp <- function(ic, iff, ix, PL) {
    Pperp <- PL
    if (c_mode == "daily") {
      c_j <- numeric(nDj)
      for (j in seq_len(nDj)) {
        sBB <- cube$s_BB[j]
        c_j[j] <- if (is.finite(sBB) && sBB > see_floor) clampc(sLE_at(j, ic, iff, ix) / sBB) else 0
        Pperp[, j] <- PL[, j] - c_j[j] * P_era5[, j]
      }
      attr(Pperp, "c") <- c_j
    } else {                                        # global ou class
      grp <- if (c_mode == "global") rep(1L, nDj) else classe
      gv  <- sort(unique(grp))
      c_g <- setNames(numeric(length(gv)), gv)
      for (g in gv) {
        cols <- which(grp == g)
        sLE <- sum(sLE_at(cols, ic, iff, ix), na.rm = TRUE)
        sBB <- sum(cube$s_BB[cols], na.rm = TRUE)
        c_g[as.character(g)] <- if (sBB > 0) clampc(sLE / sBB) else 0
        Pperp[, cols] <- PL[, cols] - c_g[as.character(g)] * P_era5[, cols, drop = FALSE]
      }
      attr(Pperp, "c") <- c_g
    }
    Pperp
  }

  dn <- list(round(cube$TAU_C), round(cube$TAU_F), round(FDRY, 2))
  kge_m    <- array(NA_real_, dim = c(nC, nF, nX), dimnames = dn)
  beta_opt <- array(NA_real_, dim = c(nC, nF, nX))

  for (ic in ic_set) for (iff in if_set) for (ix in ix_set) {
    PL <- PLPM_at(ic, iff, ix)
    Pperp <- make_perp(ic, iff, ix, PL)
    best_k <- -Inf; best_b <- NA_real_
    for (b in beta_grid) {
      Pmod <- pmax(P_era5 + b * Pperp, 0)
      k <- kge_daily_median(Pmod, P_obs, sd_min = sd_min, n_min = n_min)$kge
      if (is.finite(k) && k > best_k) { best_k <- k; best_b <- b }
    }
    kge_m[ic, iff, ix] <- best_k; beta_opt[ic, iff, ix] <- best_b
  }

  best <- which(kge_m == max(kge_m, na.rm = TRUE), arr.ind = TRUE)[1, ]
  PLb   <- PLPM_at(best[1], best[2], best[3])
  Pperp <- make_perp(best[1], best[2], best[3], PLb)
  c_at  <- attr(Pperp, "c")
  comp <- kge_daily_median(pmax(P_era5 + beta_opt[best[1],best[2],best[3]] * Pperp, 0),
                           P_obs, sd_min = sd_min, n_min = n_min)
  # c_j : vecteur c par jour au point optimal (necessaire pour reproduire P_mod).
  if (c_mode == "daily") { c_used <- NA; c_j <- c_at }
  else {                                            # global/class : etaler par groupe
    grp <- if (c_mode == "global") rep(1L, nDj) else classe
    gv  <- sort(unique(grp)); c_j <- numeric(nDj)
    for (g in gv) c_j[grp == g] <- c_at[as.character(g)]
    c_used <- c_at
  }
  names(c_j) <- format(cube$t_start, "%Y-%m-%d")
  # triplet_j / beta_j : par jour (ici commun a tous, tau_mode="global"), pour
  # compatibilite avec produce_pmod_stations.
  triplet_j <- matrix(best, nDj, 3, byrow = TRUE)
  beta_j <- rep(beta_opt[best[1], best[2], best[3]], nDj)
  list(kge = kge_m, beta_opt = beta_opt, c = c_used, c_j = c_j, c_mode = c_mode,
       tau_mode = "global", triplet_j = triplet_j, beta_j = beta_j,
       best = best, t_start = cube$t_start,
       tau_c = cube$TAU_C[best[1]], tau_f = cube$TAU_F[best[2]], fdry = FDRY[best[3]],
       beta = beta_opt[best[1], best[2], best[3]], comp = comp,
       TAU_C = cube$TAU_C, TAU_F = cube$TAU_F, FDRY = FDRY)
}

# ---------------------------------------------------------------------------
#  VERSION STREAMING de optimize_tau_kge_ortho_era5 : ne combine JAMAIS P_LPM
#  en RAM (necessaire quand le cube global est trop gros, > qq Go). Traite les
#  cubes annuels UN PAR UN.
#
#  Principe : le KGE median journalier evalue chaque jour independamment, puis
#  prend la mediane. On peut donc calculer les composantes (r, alpha, beta) de
#  chaque jour au fil de la lecture des annees, sans tout charger. Les produits
#  scalaires s_LB/s_BB (petits, sans axe poste) sont combines en RAM pour former
#  le coefficient c (modes global/class) ; en mode daily c_j est deja local.
#
#  files    : vecteur de chemins .rds (cubes annuels, cf. calibrate_cube_par).
#  P_obs_f, P_era5_f : FONCTIONS (t_start_annee) -> matrice [poste, jour] pour
#             l'annee, alignee sur les jours du cube annuel. Permet de ne charger
#             les obs / ERA5 que d'une annee a la fois. (sta_id suppose homogene.)
#  Retour : meme structure que optimize_tau_kge_ortho_era5.
# ---------------------------------------------------------------------------
optimize_tau_ortho_streaming <- function(files, P_obs_f, P_era5_f, classe = NULL,
                                         params,
                                         c_spec = list(mode = "daily"),
                                         beta_grid = exp(seq(log(1e-2), log(3), length.out = 40)),
                                         sd_min = 0.5, n_min = 20, see_frac = 0.05,
                                         clamp_c = c(0, 1),
                                         fdry_grid = seq(0.4, 1.0, by = 0.1),
                                         classe_f = NULL) {
  # -------------------------------------------------------------------------
  #  SPEC PARAMS : liste par parametre (tau_c, tau_f, fdry, beta), chacun avec
  #    $mode "global"|"class" ; $fix NULL|scalaire(global)|liste-nommee(par classe)
  #    $link (tau_c) "tau_f" pour imposer tau_c = tau_f.
  #  c_spec$mode "global"|"class"|"daily" (c reste a part : projection, pas grille).
  # -------------------------------------------------------------------------
  need <- c("tau_c", "tau_f", "fdry", "beta")
  if (missing(params) || !all(need %in% names(params)))
    stop("`params` doit lister tau_c, tau_f, fdry, beta (chacun $mode, $fix).")
  for (nm in need) if (is.null(params[[nm]]$mode) && is.null(params[[nm]]$link)) params[[nm]]$mode <- "global"
  c_mode <- match.arg(c_spec$mode, c("global", "class", "daily"))

  # tau_link : vrai SEULEMENT si un lien EXPLICITE relie tau_c et tau_f (l'un
  # declare suivre l'autre via $link). A capturer AVANT la resolution, qui efface
  # le $link en copiant le meneur. NE PAS deriver de identical() : deux specs
  # libres identiques (ex. tau_c et tau_f tous deux mode=class,fix=NULL) ne sont
  # PAS liees -- elles doivent pouvoir diverger.
  tau_link <- identical(params$tau_c$link, "tau_f") || identical(params$tau_f$link, "tau_c")

  # --- resolution des LIENS : un parametre { link = "X" } devient une COPIE de X
  #     (herite mode, fix, exclusions). Le suiveur ne specifie que $link. Interdit
  #     les chaines (X suit Y qui suit Z) et les cibles inconnues. ---
  for (nm in need) {
    lk <- params[[nm]]$link
    if (!is.null(lk)) {
      if (!lk %in% need) stop(sprintf("params$%s$link pointe vers '%s' inconnu.", nm, lk))
      if (!is.null(params[[lk]]$link))
        stop(sprintf("params$%s suit '%s' qui suit lui-meme un autre (chaine de liens interdite).", nm, lk))
      params[[nm]] <- params[[lk]]                     # copie totale du meneur
    }
  }

  # --- validation de coherence de la spec ---
  is_class_fix <- function(f) !is.null(f) && (length(f) > 1 || !is.null(names(f)))
  for (nm in need) {
    p <- params[[nm]]
    if (!p$mode %in% c("global", "class"))
      stop(sprintf("params$%s$mode doit etre 'global' ou 'class'.", nm))
    # mode global + fix par classe = "global avec exclusions" : VALIDE.
    if (p$mode == "class" && !is.null(p$fix) && !is_class_fix(p$fix))
      stop(sprintf("params$%s : fix scalaire sur un mode 'class' (declarez 'global', ou donnez un fix par classe).", nm))
  }
  any_class <- any(vapply(need, function(nm) params[[nm]]$mode == "class", logical(1)))
  if (any_class && is.null(classe))
    stop("au moins un parametre est en mode 'class' : le vecteur `classe` est requis.")
  if (any_class && c_mode == "global")
    stop("un parametre est par classe : c_mode doit etre 'class' ou 'daily' (c au moins aussi fin que le triplet).")
  clampc <- function(x) if (is.null(clamp_c)) x else pmin(pmax(x, clamp_c[1]), clamp_c[2])
  nB <- length(beta_grid)

  # --- passe 1 : lire les PETITS tableaux (s_LB_pos/neg, s_BB, t_start) ---
  # DESIGN pos/neg : s_LB(fdry) = s_LB_pos + fdry*s_LB_neg reconstruit a la volee.
  # fdry n'est PAS stocke dans le cube : c'est la grille de calage `fdry_grid`
  # qu'on balaie (indice ix), reconstruite lineairement.
  meta <- lapply(files, function(f) {
    cu <- readRDS(f)
    list(s_LB_pos = cu$s_LB_pos, s_LB_neg = cu$s_LB_neg, s_BB = cu$s_BB,
         t_start = cu$t_start, TAU_C = cu$TAU_C, TAU_F = cu$TAU_F)
  })
  TAU_C <- meta[[1]]$TAU_C; TAU_F <- meta[[1]]$TAU_F; FDRY <- fdry_grid
  nC <- length(TAU_C); nF <- length(TAU_F); nX <- length(FDRY)
  # helper : s_LB d'une annee y au fdry d'indice ix (reconstruction lineaire).
  # s_LB_pos/neg sont [jour, ic, iff] -> s_LB(fdry) meme forme.
  sLB_y <- function(y, ix) meta[[y]]$s_LB_pos + FDRY[ix] * meta[[y]]$s_LB_neg

  # --- resolution des fix_* eventuellement PAR CLASSE ---------------------
  # fix_* accepte : NULL (libre pour tous), un scalaire (fixe pour tous), un
  # vecteur nomme/positionnel avec NA=libre (une entree par classe), ou une liste
  # nommee par classe {"5"=0}. Renvoie, pour une classe donnee, la valeur fixee
  # ou NULL (libre). Les classes non citees dans un vecteur/liste sont libres.
  fix_for <- function(fix, cl) {                     # cl : etiquette de classe (chr)
    if (is.null(fix)) return(NULL)
    if (length(fix) == 1 && is.null(names(fix))) {    # scalaire -> tous
      return(if (is.na(fix)) NULL else as.numeric(fix))
    }
    v <- if (!is.null(names(fix))) fix[[cl]] else fix[suppressWarnings(as.integer(cl))]
    if (is.null(v) || (length(v) == 1 && is.na(v))) NULL else as.numeric(v)
  }
  # indices d'un axe pour une valeur fixee (ou tous si NULL)
  idx_axis <- function(grid, val) if (is.null(val)) seq_along(grid) else which.min(abs(grid - val))

  gv0 <- if (is.null(classe)) 1L else sort(unique(classe))
  # indices a EXPLORER en passe 2 pour un parametre, selon sa spec :
  #  - mode global : la valeur fixee (fix scalaire) ou toute la grille (libre) ;
  #  - mode class  : l'union sur les classes des indices autorises par le fix
  #    de chaque classe (une classe libre -> toute la grille).
  union_idx_param <- function(nm, grid) {
    p <- params[[nm]]
    if (p$mode == "global") return(idx_axis(grid, fix_for(p$fix, "1")))
    sort(Reduce(union, lapply(as.character(gv0),
                              function(cl) idx_axis(grid, fix_for(p$fix, cl)))))
  }
  ic_set <- sort(union_idx_param("tau_c", TAU_C))
  if_set <- sort(union_idx_param("tau_f", TAU_F))
  ix_set <- sort(union_idx_param("fdry",  FDRY))

  # coefficient c global/class : sommes de s_LB/s_BB sur TOUTES les annees.
  # En daily, c_j est calcule par annee dans la passe 2 (local).
  see_all <- do.call(c, lapply(meta, `[[`, "s_BB"))
  see_floor <- see_frac * median(see_all, na.rm = TRUE)
  # somme s_BB par groupe (global : 1 groupe ; class : par classe) sur tous les jours
  ts_all <- do.call(c, lapply(meta, `[[`, "t_start"))
  nDtot <- length(ts_all)
  if (is.null(classe)) classe <- rep(1L, nDtot)
  # grp_all / gv : VRAIES classes, toujours -> ventilation du raffinement et du
  # per_class par classe, quel que soit le mode de calage (option A : reporting
  # par classe systematique).
  grp_all <- classe
  gv <- sort(unique(grp_all))
  # grp_c / gv_c : groupement pour le calcul du coefficient c. Suit c_mode :
  # "global" -> un seul groupe ; "class"/"daily" -> vraies classes (daily n'utilise
  # pas cgrid mais on garde la coherence).
  grp_c <- if (c_mode == "global") rep(1L, nDtot) else classe
  gv_c  <- sort(unique(grp_c))

  # pour global/class : accumuler sum(s_LB) et sum(s_BB) par groupe de c et (ic,iff,ix)
  cgrid <- NULL
  if (c_mode != "daily") {
    sLE_sum <- array(0, dim = c(nC, nF, nX, length(gv_c)))
    sEE_sum <- setNames(numeric(length(gv_c)), gv_c)
    off <- 0
    for (y in seq_along(meta)) {
      ny <- length(meta[[y]]$t_start); idx <- (off + 1):(off + ny); off <- off + ny
      g_y <- grp_c[idx]
      for (kk in seq_along(gv_c)) {
        cols <- which(g_y == gv_c[kk])
        if (!length(cols)) next
        sEE_sum[kk] <- sEE_sum[kk] + sum(meta[[y]]$s_BB[cols], na.rm = TRUE)
        # s_LB(fdry) reconstruit pos + fdry*neg ; [jour, ic, iff] par ix
        for (ix in ix_set) {
          sLB_full <- meta[[y]]$s_LB_pos + FDRY[ix] * meta[[y]]$s_LB_neg   # [jour,ic,iff]
          for (ic in ic_set) for (iff in if_set)
            sLE_sum[ic, iff, ix, kk] <- sLE_sum[ic, iff, ix, kk] +
              sum(sLB_full[cols, ic, iff], na.rm = TRUE)
        }
      }
    }
    cgrid <- array(0, dim = c(nC, nF, nX, length(gv_c)))
    for (kk in seq_along(gv_c))
      cgrid[, , , kk] <- if (sEE_sum[kk] > 0) clampc(sLE_sum[, , , kk] / sEE_sum[kk]) else 0
  }

  # --- passe 2 : accumuler les composantes journalieres, annee par annee ---
  # stockage des composantes journalieres. On indexe sur les seuls indices
  # EXPLORES (ic_set/if_set/ix_set) : si tau/fdry sont fixes, les hypercubes sont
  # bien plus petits (gain memoire et temps majeur). map_* : indice global -> position locale.
  nic <- length(ic_set); nif <- length(if_set); nix <- length(ix_set)
  map_ic <- match(seq_len(nC), ic_set)               # NA hors exploration
  map_if <- match(seq_len(nF), if_set)
  map_ix <- match(seq_len(nX), ix_set)
  dimK <- c(nDtot, nic, nif, nix, nB)
  R_j <- array(NA_real_, dimK); A_j <- array(NA_real_, dimK); B_j <- array(NA_real_, dimK)

  off <- 0
  for (y in seq_along(meta)) {
    cu <- readRDS(files[y])                         # charge UNE annee (~1 Go)
    tsy <- cu$t_start; ny <- length(tsy)
    idx_glob <- (off + 1):(off + ny); off <- off + ny
    Po <- P_obs_f(tsy)                              # [poste, jour] obs de l'annee
    Pe <- P_era5_f(tsy)                             # [poste, jour] ERA5 de l'annee
    g_y_c <- grp_c[idx_glob]                          # groupe de c (suit c_mode)

    for (ic in ic_set) for (iff in if_set) for (ix in ix_set) {
      # tau_link : ne calculer QUE la diagonale tau_c = tau_f. Ce saut est ici
      # (passe 2, point chaud) et non seulement dans agg_over, sinon on paierait
      # le calcul des ~90 couples hors-diagonale inutiles.
      if (isTRUE(tau_link) && iff != which.min(abs(TAU_F - TAU_C[ic]))) next
      # reconstruction lineaire fdry : P_LPM(fdry) et s_LB(fdry)
      PL <- cu$P_LPM_pos[, , ic, iff] + FDRY[ix] * cu$P_LPM_neg[, , ic, iff]   # [poste, jour]
      # coefficient c applique a cette annee
      if (c_mode == "daily") {
        sLB_ijk <- cu$s_LB_pos[, ic, iff] + FDRY[ix] * cu$s_LB_neg[, ic, iff]  # [jour]
        cvec <- numeric(ny)
        for (j in seq_len(ny)) {
          sBB <- cu$s_BB[j]
          cvec[j] <- if (is.finite(sBB) && sBB > see_floor) clampc(sLB_ijk[j] / sBB) else 0
        }
      } else {
        cvec <- numeric(ny)
        for (kk in seq_along(gv_c)) cvec[g_y_c == gv_c[kk]] <- cgrid[ic, iff, ix, kk]
      }
      Pperp <- PL - sweep(Pe, 2, cvec, `*`)         # [poste, jour] : LPM - c*ERA5
      for (ib in seq_len(nB)) {
        Pmod <- pmax(Pe + beta_grid[ib] * Pperp, 0)
        # composantes journalieres (meme regles que kge_daily_median)
        for (j in seq_len(ny)) {
          m <- Pmod[, j]; o <- Po[, j]
          ok <- is.finite(m) & is.finite(o)
          if (sum(ok) < n_min) next
          m <- m[ok]; o <- o[ok]; so <- sd(o)
          if (!is.finite(so) || so < sd_min) next
          r <- suppressWarnings(cor(m, o)); a <- sd(m) / so; b <- mean(m) / mean(o)
          if (!is.finite(r) || !is.finite(a) || !is.finite(b)) next
          jg <- idx_glob[j]
          R_j[jg, map_ic[ic], map_if[iff], map_ix[ix], ib] <- r
          A_j[jg, map_ic[ic], map_if[iff], map_ix[ix], ib] <- a
          B_j[jg, map_ic[ic], map_if[iff], map_ix[ix], ib] <- b
        }
      }
    }
    rm(cu, PL, Pperp, Po, Pe); gc(FALSE)
  }

  # --- agregation : recherche du triplet optimal + HYPERCUBES medians complets ---
  # ics/ifs/ixs/ibs : indices A EXPLORER (defaut = tous ceux calcules en passe 2).
  # Permet de restreindre PAR CLASSE (fix_* par classe). Les hypercubes de sortie
  # restent pleins [nC,nF,nX,nB] ; seuls les noeuds explores sont remplis.
  agg_over <- function(days, ics = ic_set, ifs = if_set, ixs = ix_set, ibs = seq_len(nB)) {
    dh <- c(nC, nF, nX, nB)
    KGE <- array(NA_real_, dh); Rm <- array(NA_real_, dh)
    Am  <- array(NA_real_, dh); Bm <- array(NA_real_, dh); Nn <- array(0L, dh)
    bk <- -Inf; bbest <- NULL; bb <- NA_real_; bcomp <- NULL
    for (ic in ics) for (iff in ifs) for (ix in ixs) for (ib in ibs) {
      # tau_link : ne garder que la diagonale tau_c = tau_f (crete tau_c+tau_f).
      # On apparie iff a l'indice de TAU_F le plus proche de TAU_C[ic] ; les
      # autres couples sont sautes -> exploration 1D en tau.
      if (isTRUE(tau_link) && iff != which.min(abs(TAU_F - TAU_C[ic]))) next
      lic <- map_ic[ic]; lif <- map_if[iff]; lix <- map_ix[ix]   # positions locales
      rj <- R_j[days, lic, lif, lix, ib]
      aj <- A_j[days, lic, lif, lix, ib]
      bj <- B_j[days, lic, lif, lix, ib]
      kj <- 1 - sqrt((rj - 1)^2 + (aj - 1)^2 + (bj - 1)^2)
      km <- median(kj, na.rm = TRUE)
      KGE[ic,iff,ix,ib] <- km
      Rm[ic,iff,ix,ib]  <- median(rj, na.rm = TRUE)
      Am[ic,iff,ix,ib]  <- median(aj, na.rm = TRUE)
      Bm[ic,iff,ix,ib]  <- median(bj, na.rm = TRUE)
      Nn[ic,iff,ix,ib]  <- sum(is.finite(kj))
      if (is.finite(km) && km > bk) {
        bk <- km; bbest <- c(ic, iff, ix); bb <- beta_grid[ib]
        bcomp <- list(kge = km, r = Rm[ic,iff,ix,ib], alpha = Am[ic,iff,ix,ib],
                      beta = Bm[ic,iff,ix,ib], n_used = Nn[ic,iff,ix,ib])
      }
    }
    dn <- list(round(TAU_C), round(TAU_F), round(FDRY, 2), signif(beta_grid, 3))
    dimnames(KGE) <- dn; dimnames(Rm) <- dn; dimnames(Am) <- dn; dimnames(Bm) <- dn
    list(best = bbest, beta = bb, kge = bk, comp = bcomp,
         kge_cube = KGE, r_cube = Rm, alpha_cube = Am, beta_cube = Bm, n_cube = Nn)
  }
  # indices beta a explorer pour une valeur fixee (ou tous)
  ib_for <- function(val) if (is.null(val)) seq_len(nB) else which.min(abs(beta_grid - val))

  # ==========================================================================
  #  SELECTION UNIFORME DES OPTIMA (2 niveaux, KGE global exact - option 3).
  #  Chaque parametre (tau_c, tau_f, fdry, beta) a un mode global|class et un fix.
  #   - EXTERNE : boucle sur la grille produit des axes des params GLOBAL.
  #   - INTERNE : pour chaque point externe, raffiner par classe les params CLASS.
  #   - SCORE   : KGE global exact (mediane sur tous les jours, chaque jour a
  #     l'indice retenu de sa classe), recompose depuis R_j/A_j/B_j.
  # ==========================================================================
  triplet_j <- matrix(NA_integer_, nDtot, 3); beta_j <- numeric(nDtot)
  gvc <- as.character(gv)
  is_glob <- function(nm) params[[nm]]$mode == "global"
  # un param GLOBAL peut avoir des EXCLUSIONS : un fix par classe (liste/vecteur)
  # signifie "valeur commune calee sur les classes NON fixees, valeur imposee aux
  # classes fixees". has_excl(nm) : ce param global a-t-il des exclusions ?
  has_excl <- function(nm) is_glob(nm) &&
    !is.null(params[[nm]]$fix) && (length(params[[nm]]$fix) > 1 || !is.null(names(params[[nm]]$fix)))
  # valeur fixee d'un param pour une classe (NULL si non fixee / non exclue)
  excl_val <- function(nm, cl) if (has_excl(nm)) fix_for(params[[nm]]$fix, cl) else NULL
  grids_p <- list(tau_c = TAU_C, tau_f = TAU_F, fdry = FDRY, beta = beta_grid)

  # KGE journalier d'un ensemble de jours a un point (indices globaux ic,iff,ix,ib)
  kj_at <- function(days, ic, iff, ix, ib) {
    lic <- map_ic[ic]; lif <- map_if[iff]; lix <- map_ix[ix]
    rj <- R_j[days, lic, lif, lix, ib]; aj <- A_j[days, lic, lif, lix, ib]
    bj <- B_j[days, lic, lif, lix, ib]
    1 - sqrt((rj - 1)^2 + (aj - 1)^2 + (bj - 1)^2)
  }
  # raffine les params CLASS d'une classe, params globaux figes a `gfix` (indices).
  # Pour un param GLOBAL : la classe prend la valeur commune `gfix`, SAUF si elle
  # est exclue (fix par classe) -> elle prend alors sa valeur fixee.
  refine_class <- function(cl, gfix) {
    days <- which(grp_all == gv[match(cl, gvc)])
    ax_glob <- function(nm, grid, gval) {              # indice pour un param global
      ev <- excl_val(nm, cl)
      if (!is.null(ev)) which.min(abs(grid - ev)) else gval   # exclue -> fixee ; sinon commune
    }
    icc <- if (is_glob("tau_c")) ax_glob("tau_c", TAU_C, gfix$tau_c) else idx_axis(TAU_C, fix_for(params$tau_c$fix, cl))
    iff <- if (is_glob("tau_f")) ax_glob("tau_f", TAU_F, gfix$tau_f) else idx_axis(TAU_F, fix_for(params$tau_f$fix, cl))
    ixx <- if (is_glob("fdry"))  ax_glob("fdry",  FDRY,  gfix$fdry)  else idx_axis(FDRY,  fix_for(params$fdry$fix,  cl))
    ibb <- if (is_glob("beta"))  ax_glob("beta",  beta_grid, gfix$beta) else ib_for(fix_for(params$beta$fix, cl))
    res <- agg_over(days, ics = icc, ifs = iff, ixs = ixx, ibs = ibb)
    list(days = days, res = res, ib = which.min(abs(beta_grid - res$beta)))
  }
  # axes GLOBAUX pour la boucle externe. Un param global (avec ou sans exclusions)
  # balaie sa grille pour la VALEUR COMMUNE (fix_for(...,"1") ne s'applique qu'au
  # fix scalaire ; un fix par classe = exclusions, n'affecte pas la grille commune).
  gax <- function(nm) {
    if (!is_glob(nm)) return(NA_integer_)
    sc <- if (has_excl(nm)) NULL else fix_for(params[[nm]]$fix, "1")  # scalaire eventuel
    idx_axis(grids_p[[nm]], sc)
  }
  ext_grid <- expand.grid(tau_c = gax("tau_c"), tau_f = gax("tau_f"),
                          fdry = gax("fdry"), beta = gax("beta"), KEEP.OUT.ATTRS = FALSE)

  best_ext <- list(kge = -Inf)
  for (r in seq_len(nrow(ext_grid))) {
    gfix <- list(tau_c = ext_grid$tau_c[r], tau_f = ext_grid$tau_f[r],
                 fdry = ext_grid$fdry[r], beta = ext_grid$beta[r])
    # tau_link au niveau externe : si tau_c et tau_f globaux, ne garder que la
    # diagonale (agg_over gere le cas class ; ici on filtre le couple global).
    if (tau_link && is_glob("tau_c") && is_glob("tau_f") &&
        gfix$tau_f != which.min(abs(TAU_F - TAU_C[gfix$tau_c]))) next
    per_cl <- lapply(gvc, function(cl) refine_class(cl, gfix)); names(per_cl) <- gvc
    kj_all <- rep(NA_real_, nDtot)
    for (cl in gvc) { x <- per_cl[[cl]]; b <- x$res$best
      kj_all[x$days] <- kj_at(x$days, b[1], b[2], b[3], x$ib) }
    score <- median(kj_all, na.rm = TRUE)
    if (is.finite(score) && score > best_ext$kge)
      best_ext <- list(kge = score, gfix = gfix, per_cl = per_cl, kj_all = kj_all)
  }

  # assemblage final depuis le meilleur point externe (deja calcule)
  per_class <- list()
  for (cl in gvc) {
    x <- best_ext$per_cl[[cl]]; bb <- x$res$best
    triplet_j[x$days, ] <- matrix(bb, length(x$days), 3, byrow = TRUE)
    beta_j[x$days] <- x$res$beta
    per_class[[cl]] <- list(tau_c = TAU_C[bb[1]], tau_f = TAU_F[bb[2]], fdry = FDRY[bb[3]],
                            beta = x$res$beta, kge = x$res$kge, comp = x$res$comp,
                            n = length(x$days), kge_cube = x$res$kge_cube,
                            r_cube = x$res$r_cube, alpha_cube = x$res$alpha_cube,
                            beta_cube = x$res$beta_cube, n_cube = x$res$n_cube)
  }
  gf <- best_ext$gfix
  best <- c(gf$tau_c, gf$tau_f, gf$fdry)               # indices globaux du triplet
  beta_star <- beta_grid[which.min(abs(beta_grid - median(beta_j, na.rm = TRUE)))]
  kge_global_exact <- best_ext$kge; kj_all <- best_ext$kj_all

  # coefficient c_j au triplet retenu pour chaque jour (necessaire pour P_mod).
  # c depend du triplet du jour : on lit s_LB au (ic,iff,ix) propre a ce jour.
  c_j <- numeric(nDtot); off <- 0
  for (y in seq_along(meta)) {
    ny <- length(meta[[y]]$t_start); idx <- (off + 1):(off + ny); off <- off + ny
    for (jj in seq_len(ny)) {
      jg <- idx[jj]; tr <- triplet_j[jg, ]
      if (anyNA(tr)) next
      if (c_mode == "daily") {
        sBB <- meta[[y]]$s_BB[jj]
        # s_LB(fdry=tr[3]) = pos + FDRY[tr[3]]*neg au (jj, tr[1], tr[2])
        sLB_v <- meta[[y]]$s_LB_pos[jj, tr[1], tr[2]] + FDRY[tr[3]] * meta[[y]]$s_LB_neg[jj, tr[1], tr[2]]
        c_j[jg] <- if (is.finite(sBB) && sBB > see_floor) clampc(sLB_v / sBB) else 0
      } else {                                        # class/global : valeur du groupe de c
        kk <- match(grp_c[jg], gv_c)
        c_j[jg] <- cgrid[tr[1], tr[2], tr[3], kk]
      }
    }
  }
  names(c_j) <- format(ts_all, "%Y-%m-%d")

  # KGE global : mediane sur TOUS les jours, chacun evalue avec SON triplet.
  kj_all <- numeric(nDtot) * NA_real_
  for (jg in seq_len(nDtot)) {
    tr <- triplet_j[jg, ]; if (anyNA(tr)) next
    ib <- which.min(abs(beta_grid - beta_j[jg]))
    lic <- map_ic[tr[1]]; lif <- map_if[tr[2]]; lix <- map_ix[tr[3]]
    r <- R_j[jg, lic, lif, lix, ib]; a <- A_j[jg, lic, lif, lix, ib]; b <- B_j[jg, lic, lif, lix, ib]
    if (is.finite(r) && is.finite(a) && is.finite(b))
      kj_all[jg] <- 1 - sqrt((r-1)^2 + (a-1)^2 + (b-1)^2)
  }
  comp_global <- list(kge = median(kj_all, na.rm = TRUE), n_used = sum(is.finite(kj_all)))
  kge_global_exact <- comp_global$kge               # coherence : etat final assemble

  # sortie uniforme : champs par jour + per_class (toujours) + valeurs globales
  # du triplet (best) + scalaires de compat. per_class contient les 4 params par
  # classe (ceux en mode global ont la meme valeur repetee sur les classes).
  out <- list(c_mode = c_mode, c_j = c_j, beta_j = beta_j,
              triplet_j = triplet_j, t_start = ts_all, params = params,
              per_class = per_class, classes = gv,
              kge_global = kge_global_exact, kj_all = kj_all,
              TAU_C = TAU_C, TAU_F = TAU_F, FDRY = FDRY)
  # valeurs globales du triplet (indices best) : utiles pour production gridded et
  # affichage. En presence de params class, ce sont les valeurs GLOBALES retenues
  # (les params class variant, voir per_class). tau_c/tau_f/fdry scalaires = compat.
  out$best <- best
  out$tau_c <- TAU_C[best[1]]; out$tau_f <- TAU_F[best[2]]; out$fdry <- FDRY[best[3]]
  out$beta  <- beta_star
  out$comp <- comp_global                           # KGE global exact (etat final)
  # hypercubes de diagnostic : dans per_class[[cl]]. Pour un acces "global"
  # (cas tout-global, une seule classe), exposer le cube de la 1re classe.
  if (length(per_class) >= 1) {
    p1 <- per_class[[1]]
    out$kge_cube <- p1$kge_cube; out$r_cube <- p1$r_cube
    out$alpha_cube <- p1$alpha_cube; out$beta_cube <- p1$beta_cube; out$n_cube <- p1$n_cube
  }
  if (c_mode != "daily") { cc <- cgrid[best[1], best[2], best[3], ]; names(cc) <- as.character(gv_c); out$c <- cc }
  else out$c <- NA
  out
}

# ---------------------------------------------------------------------------
#  make_calib_tag : construit un tag de calage SANS AMBIGUITE depuis la spec.
#  Convention (portable Windows : [A-Za-z0-9._-] uniquement) :
#    - un segment par parametre : param-granularite[-statut], _ entre params.
#    - granularite toujours presente : glob | cls.
#    - global cale : `fdry-glob` ; global fixe : `fdry-glob-1` (valeur = fixation).
#    - global AVEC EXCLUSIONS (valeur commune calee + classes fixees) :
#      `beta-glob-fix5-0` (glob = commune calee ; fix5-0 = classe 5 exclue a 0).
#    - class cale : `beta-cls` ; class fixe : `beta-cls-fix5-0` (groupes de classes
#      par valeur commune, chiffres colles trie, groupes lies par -).
#    - lien : `-link` sur tau (tau_c$link == "tau_f").
#    - suffixes : _seed<n> (CV). clamp optionnel si non-standard.
#  Ex. : tau-cls-link_fdry-glob_beta-cls-fix5-0_c-daily_seed1
# ---------------------------------------------------------------------------
make_calib_tag <- function(params, c_mode = "daily", seed = NULL, clamp_c = NULL) {
  fmt <- function(v) sub("\\.?0+$", "", sprintf("%.3f", v))   # 0.700 -> 0.7 ; 1.000 -> 1
  # encode un fix PAR CLASSE en groupes "fix<classes>-<valeur>" (classes par
  # valeur commune, chiffres colles, groupes tries et lies par -).
  fix_by_class_str <- function(f) {
    cls <- if (!is.null(names(f))) names(f) else as.character(seq_along(f))
    val <- if (!is.null(names(f))) unlist(f) else f
    keep <- !is.na(val); cls <- cls[keep]; val <- val[keep]
    if (!length(val)) return("")
    uv <- sort(unique(val))
    grps <- vapply(uv, function(vv) {
      ks <- sort(as.integer(cls[val == vv]))
      paste0("fix", paste(ks, collapse = ""), "-", fmt(vv))
    }, character(1))
    paste0("-", paste(grps, collapse = "-"))
  }
  seg_param <- function(nm, short) {
    p <- params[[nm]]; g <- if (p$mode == "global") "glob" else "cls"
    stat <- ""
    if (!is.null(p$fix)) {
      is_cls_fix <- length(p$fix) > 1 || !is.null(names(p$fix))
      if (!is_cls_fix) {                               # fix scalaire (global) : valeur
        v <- p$fix; if (!is.na(v)) stat <- paste0("-", fmt(v))
      } else {                                         # fix par classe : groupes
        # en mode global -> exclusions (glob = valeur commune calee) ;
        # en mode class  -> classes fixees. Meme encodage des groupes.
        stat <- fix_by_class_str(p$fix)
      }
    }
    paste0(short, "-", g, stat)
  }
  # tau_c / tau_f : lies si l'un declare suivre l'autre (link). Un seul segment
  # `tau-...-link` alors (le suiveur est identique au meneur, resolu ici).
  linked <- identical(params$tau_c$link, "tau_f") || identical(params$tau_f$link, "tau_c")
  tau_seg <- if (linked) {
    lead <- if (identical(params$tau_f$link, "tau_c")) "tau_c" else "tau_f"  # le meneur
    base <- seg_param(lead, "tau")            # "tau-<glob|cls>[-stat]"
    sub("^tau-([a-z]+)", "tau-\\1-link", base)          # insere -link apres la granularite
  } else if (identical(params$tau_c$mode, params$tau_f$mode) &&
             is.null(params$tau_c$fix) && is.null(params$tau_f$fix)) {
    g <- if (params$tau_c$mode == "global") "glob" else "cls"
    paste0("tau-", g)
  } else {
    paste(seg_param("tau_c", "tauc"), seg_param("tau_f", "tauf"), sep = "_")
  }
  parts <- c(tau_seg, seg_param("fdry", "fdry"), seg_param("beta", "beta"),
             paste0("c-", c_mode))
  if (!is.null(clamp_c) && !identical(clamp_c, c(0, 1)))
    parts <- c(parts, sprintf("clamp%s.%s", fmt(clamp_c[1]), fmt(clamp_c[2])))
  if (!is.null(seed)) parts <- c(parts, sprintf("seed%d", seed))
  paste(parts, collapse = "_")
}

# ---------------------------------------------------------------------------
#  PRODUCTION de P_mod AUX POSTES au point optimal d'un calage (mode "stations").
#  Reconstruit, en streaming (une annee a la fois), le champ downscale evalue
#  aux postes :
#     P_mod = max(P_era5 + beta * (P_LPM_opt - c_j * P_era5), 0)
#  a partir de la tranche optimale (tau_c, tau_f, fdry) des cubes annuels.
#  Ne recalcule PAS le LPM (il est deja dans le cube aux postes) -> pour la
#  version GRIDDED, une fonction distincte rejouera le LPM sur la grille fine.
#
#  opt      : sortie d'un calage streaming. Utilise opt$triplet_j (indice
#             (ic,iff,ix) par jour), opt$beta_j (beta par jour) et opt$c_j
#             (c par jour) -> gere aussi bien tau_mode global que class (chaque
#             jour applique le triplet de sa classe).
#  files    : cubes annuels (memes que le calage).
#  P_obs_f, P_era5_f : fonctions (t_start_annee) -> [poste, jour] (comme le calage).
#  Retour   : liste P_mod / P_obs / P_era5 [poste, jour_total], t_start, sta_id.
# ---------------------------------------------------------------------------
produce_pmod_stations <- function(opt, files, P_obs_f, P_era5_f) {
  m1 <- readRDS(files[1]); sta_id <- m1$sta_id

  # index par date : triplet (ic,iff,ix), beta, c pour chaque jour
  dts <- if (!is.null(names(opt$c_j))) names(opt$c_j) else format(opt$t_start, "%Y-%m-%d")
  tri_by_date  <- opt$triplet_j; rownames(tri_by_date) <- dts
  beta_by_date <- setNames(opt$beta_j, dts)
  cj_by_date   <- setNames(as.numeric(opt$c_j), dts)

  Pmod <- vector("list", length(files)); Pobs <- Pmod; Pera <- Pmod; ts <- Pmod
  FDRY <- opt$FDRY                                   # grille fdry du calage
  for (y in seq_along(files)) {
    cu  <- readRDS(files[y]); tsy <- cu$t_start; ny <- length(tsy)
    dy  <- format(tsy, "%Y-%m-%d")
    Po  <- P_obs_f(tsy); Pe <- P_era5_f(tsy)
    Pm  <- matrix(NA_real_, nrow(Pe), ny)
    for (j in seq_len(ny)) {
      tr <- tri_by_date[dy[j], ]                     # (ic,iff,ix) du jour
      if (anyNA(tr)) next
      fdry <- FDRY[tr[3]]                             # fdry du triplet (grille de calage)
      PLj <- cu$P_LPM_pos[, j, tr[1], tr[2]] + fdry * cu$P_LPM_neg[, j, tr[1], tr[2]]  # [poste]
      cj  <- cj_by_date[dy[j]]; if (is.na(cj)) cj <- 0
      bj  <- beta_by_date[dy[j]]
      Pm[, j] <- pmax(Pe[, j] + bj * (PLj - cj * Pe[, j]), 0)
    }
    Pmod[[y]] <- Pm; Pobs[[y]] <- Po; Pera[[y]] <- Pe; ts[[y]] <- tsy
    rm(cu, Po, Pe, Pm); gc(FALSE)
  }
  list(P_mod  = do.call(cbind, Pmod),
       P_obs  = do.call(cbind, Pobs),
       P_era5 = do.call(cbind, Pera),
       t_start = do.call(c, ts), sta_id = sta_id,
       tau_mode = opt$tau_mode, c_mode = opt$c_mode)
}

# ---------------------------------------------------------------------------
#  Generation ROBUSTE des journees d'observation 06-06 UTC sur une plage.
#  seq(..., by="1 day") gere automatiquement les annees bissextiles.
#  Retourne un data.frame $t_start (POSIXct 06 UTC, 1er jan year0 -> 31 dec year1).
# ---------------------------------------------------------------------------
gen_jours <- function(year0, year1) {
  d0 <- as.POSIXct(sprintf("%d-01-01 06:00:00", year0), tz = "UTC")
  d1 <- as.POSIXct(sprintf("%d-12-31 06:00:00", year1), tz = "UTC")
  data.frame(t_start = seq(d0, d1, by = "1 day"))
}

# ---------------------------------------------------------------------------
#  Disponibilite d'une journee : les 8 pas 06->03(j+1) doivent tous exister
#  dans les fichiers (annee j et, pour les 2 derniers pas, annee j+1).
#  Teste l'existence des fichiers requis SANS les charger. Renvoie TRUE/FALSE.
#  path_common(year), path_sfc(year) : memes fonctions que make_coarse_buffer.
# ---------------------------------------------------------------------------
# ============================================================================

# ---------------------------------------------------------------------------
#  RR      : data.frame [jour x (Date + postes)] ; RR$Date en "YYYY-MM-DD".
#            NB : l'ordre des colonnes postes DOIT correspondre a idx_sta.
#  jours   : data.frame $t_start (POSIXct 06 UTC) des journees candidates.
#  seuil   : seuil de pluie moyenne (mm) pour retenir un jour (defaut 1).
#  n_min   : nombre minimal de postes disponibles (non NA) ce jour-la (defaut 20).
#  Renvoie : liste $jours_wet (sous-df de jours retenus), $mean_rr (pluie moyenne
#            par jour retenu), $idx (indices dans jours), $rr_obs (matrice
#            [jour_retenu x poste] des obs alignees, pour la fonction de cout).
# ---------------------------------------------------------------------------
select_wet_days <- function(RR, jours, seuil = 1, n_min = 20) {
  rr_dates <- as.Date(RR$Date)
  obs_mat  <- as.matrix(RR[, -1, drop = FALSE])       # [jour x poste]

  # date calendaire du jour j = debut de fenetre 06-06
  jour_dates <- as.Date(format(jours$t_start, "%Y-%m-%d"))

  # ligne de RR correspondant a chaque journee de calage
  row_idx <- match(jour_dates, rr_dates)              # NA si date absente de RR

  keep <- logical(nrow(jours))
  mean_rr <- rep(NA_real_, nrow(jours))
  navail  <- rep(0L, nrow(jours))

  for (d in seq_len(nrow(jours))) {
    r <- row_idx[d]
    if (is.na(r)) next                                # date non couverte par RR
    v <- obs_mat[r, ]
    na  <- sum(is.finite(v))
    navail[d] <- na
    if (na < n_min) next                              # trop peu de postes -> non retenu
    mr <- mean(v, na.rm = TRUE)
    mean_rr[d] <- mr
    if (mr > seuil) keep[d] <- TRUE
  }

  idx <- which(keep)
  list(
    jours_wet = jours[idx, , drop = FALSE],
    mean_rr   = mean_rr[idx],
    idx       = idx,
    n_avail   = navail[idx],
    rr_obs    = obs_mat[row_idx[idx], , drop = FALSE], # [jour_retenu x poste]
    n_total   = nrow(jours),
    n_kept    = length(idx)
  )
}

# ---------------------------------------------------------------------------
#  Lit la pluie ERA5 journaliere 06-06 EN CHAMP (grille LPM = celle de A_agg),
#  pour l'orthogonalisation. Renvoie [nera5, nD] aligne sur t_start (mm), ordre
#  [lon,lat] col-major (comme A_agg attend). nc_daily : era5_daily_0606 (tp mm).
#  Appariement par DATE CALENDAIRE (robuste a un decalage horaire eventuel).
# ---------------------------------------------------------------------------
read_era5_daily_field <- function(nc_daily, t_start) {
  nc <- ncdf4::nc_open(nc_daily)
  tim <- ncdf4::ncvar_get(nc, "time")
  nlon <- nc$dim$longitude$len; nlat <- nc$dim$latitude$len
  date_file <- as.Date(as.POSIXct(tim, origin = "1970-01-01", tz = "UTC"))
  date_want <- as.Date(as.POSIXct(as.numeric(t_start), origin = "1970-01-01", tz = "UTC"))
  it <- match(date_want, date_file)
  nD <- length(t_start)
  out <- matrix(NA_real_, nlon * nlat, nD)
  for (d in seq_len(nD)) {
    if (is.na(it[d])) next
    tp <- ncdf4::ncvar_get(nc, "tp", start = c(1,1,it[d]), count = c(-1,-1,1))
    tp[tp == -9999] <- NA_real_
    out[, d] <- as.vector(tp)                            # deja en mm, ordre [lon,lat]
  }
  ncdf4::nc_close(nc)
  out
}

# ---------------------------------------------------------------------------
#  Lit la pluie journaliere ERA5 (fichier tp[time,lat,lon], unite m, _Fill=-9999)
#  et l'interpole aux postes. Grille propre (differente des lpm_common).
#  Renvoie [poste, jour] en mm, aligne sur les timestamps demandes (t_start).
# ---------------------------------------------------------------------------
read_era5_daily_at_stations <- function(nc_path, pts_sf, t_start) {
  nc <- ncdf4::nc_open(nc_path)
  lon <- ncdf4::ncvar_get(nc, "longitude"); lat <- ncdf4::ncvar_get(nc, "latitude")
  tim <- ncdf4::ncvar_get(nc, "time")                        # secondes depuis 1970
  # matrice d'interpolation propre a CETTE grille
  bl <- build_bilinear_to_stations(lon, lat, pts_sf)
  I  <- bl$I
  # appariement par DATE CALENDAIRE : le fichier date les cumuls 06-06 a 00:00 UTC
  # (debut du jour calendaire) alors que t_start est a 06:00. On compare donc les
  # dates, pas les timestamps exacts (evite le decalage systematique de 6 h).
  date_file <- as.Date(as.POSIXct(tim, origin = "1970-01-01", tz = "UTC"))
  date_want <- as.Date(as.POSIXct(as.numeric(t_start), origin = "1970-01-01", tz = "UTC"))
  it <- match(date_want, date_file)
  if (anyNA(it)) warning(sprintf("%d jour(s) absent(s) du fichier ERA5 daily.", sum(is.na(it))))
  npt <- nrow(I); nD <- length(t_start)
  out <- matrix(NA_real_, npt, nD)
  for (d in seq_len(nD)) {
    if (is.na(it[d])) next
    tp <- ncdf4::ncvar_get(nc, "tp", start = c(1, 1, it[d]), count = c(-1, -1, 1))
    # tp lu en [lon,lat] (R inverse l'ordre C [time,lat,lon]) ; _Fill -> NA ; m -> mm
    tp[tp == -9999] <- NA_real_
    v <- as.vector(tp) * 1000                               # [lon*lat] col-major, mm
    out[, d] <- as.vector(I %*% v)
  }
  ncdf4::nc_close(nc)
  out                                                        # [poste, jour] mm
}

# ---------------------------------------------------------------------------
#  EXEMPLE

day_available <- function(t_start, path_common, path_sfc) {
  ts0    <- as.numeric(t_start)
  stamps <- ts0 + 3 * 3600 * (0:7)
  years  <- unique(as.integer(format(
    as.POSIXct(stamps, origin = "1970-01-01", tz = "UTC"), "%Y")))
  all(vapply(years, function(y)
    file.exists(path_common(y)) && file.exists(path_sfc(y)),
    logical(1)))
}

# ---------------------------------------------------------------------------
#  Masque polygone : vecteur logique [ncell] des mailles fines dans le polygone.
#  poly_sf : objet sf POLYGON (meme CRS que pre$template, sinon reprojete).
# ---------------------------------------------------------------------------
build_domain_mask <- function(pre, poly_sf) {
  poly_v <- terra::vect(poly_sf)
  if (!terra::same.crs(poly_v, pre$template))
    poly_v <- terra::project(poly_v, pre$template)
  r <- terra::rasterize(poly_v, pre$template, background = NA, touches = TRUE)
  mask_ras <- !is.na(terra::values(r))                # [ncell] logique, ordre raster
  as.logical(mask_ras)
}

# ---------------------------------------------------------------------------
#  Masque CONTINENTAL a la resolution ERA5 : mailles ERA5 dont une fraction
#  >= `frac_min` de la surface fine est dans le polygone terrestre (mask_vec).
#  Sert a restreindre l'orthogonalisation ERA5 au continent : en mer le LPM est
#  quasi nul mais ERA5 pleut, ce qui gonfle <P_era5,P_era5> et ecrase le coef c.
#  A_agg : agregation fine->ERA5 (lignes = mailles ERA5). mask_vec : [ncell] terre.
#  Retourne un vecteur logique [nera5].
# ---------------------------------------------------------------------------
build_continental_mask <- function(A_agg, mask_vec, frac_min = 0.5) {
  # fraction terrestre de chaque maille ERA5 = (surface fine terrestre) / (surface fine totale)
  w_land  <- as.vector(A_agg %*% as.numeric(mask_vec))   # somme des poids fins terrestres
  w_tot   <- as.vector(Matrix::rowSums(A_agg))           # somme totale des poids fins
  frac    <- ifelse(w_tot > 0, w_land / w_tot, 0)
  frac >= frac_min                                        # [nera5] logique
}

# ---------------------------------------------------------------------------
#  Tag court des ddl spatiaux DESACTIVES, pour nommer le dossier des cubes de
#  facon reproductible : "full" si tout TRUE, sinon "no-<liste des FALSE>".
#  Centralise ICI (pas dans les drivers) pour que la construction et le calage
#  reconstruisent le MEME OUT_DIR : ne modifier que la variable SPATIAL, jamais
#  le vecteur de defaut ci-dessous.
#    ex. c(rhoS=F,Hw=F,U=T,V=T,gamma=F,zeta=T) -> "no-rhoS-Hw-gamma"
# ---------------------------------------------------------------------------
spatial_tag <- function(spatial) {
  def <- c(rhoS = TRUE, Hw = TRUE, U = TRUE, V = TRUE, gamma = TRUE, zeta = TRUE)
  s <- def; if (!is.null(names(spatial))) s[names(spatial)] <- spatial else s[] <- spatial
  off <- names(s)[!s]
  if (length(off) == 0) "full" else paste0("no-", paste(off, collapse = "-"))
}

# ---------------------------------------------------------------------------
#  Vecteur d'indices postes : idx[p] = cellule fine contenant le poste p.
#  pts_sf : objet sf POINTS (postes pluvio).
# ---------------------------------------------------------------------------
build_station_index <- function(pre, pts_sf) {
  pts_v <- terra::vect(pts_sf)
  if (!terra::same.crs(pts_v, pre$template))
    pts_v <- terra::project(pts_v, pre$template)
  xy  <- terra::crds(pts_v)
  idx <- terra::cellFromXY(pre$template, xy)          # indice lineaire raster
  idx                                                 # [n_postes]
}

# ---------------------------------------------------------------------------
#  SETUP partage par pas de temps : tout ce qui NE depend PAS de (tau_c,tau_f).
#  mask_vec : masque polygone [ncell] pour l'etat uniforme (build_domain_mask).
#  Retourne une liste d'objets spectraux reutilisables par lpm_apply_tau.
# ---------------------------------------------------------------------------
lpm_setup_timestep <- function(pre, fields, mask_vec,
                               Lc_wind = 3000, coef_sigma_eps = 0.25,
                               opt = c("sfc", "lcl"),
                               spatial = c(rhoS = TRUE, Hw = TRUE, U = TRUE,
                                           V = TRUE, gamma = TRUE, zeta = TRUE)) {
  opt <- match.arg(opt)
  # normalise le vecteur `spatial` : les ddl non mentionnes gardent le defaut TRUE.
  # TRUE = le parametre varie spatialement (correction de Taylor active) ;
  # FALSE = parametre maintenu uniforme a sa valeur de reference.
  spat_def <- c(rhoS = TRUE, Hw = TRUE, U = TRUE, V = TRUE, gamma = TRUE, zeta = TRUE)
  if (!is.null(names(spatial))) spat_def[names(spatial)] <- spatial else spat_def[] <- spatial
  spatial <- spat_def
  if (opt != "lcl") spatial["zeta"] <- FALSE     # pas de zeta hors mode lcl
  cst <- get_constants(); g <- cst$g

  Hw   <- fields$Hw_sfc;   gam <- fields$gamma_sfc; Gm <- fields$Gamma_m_sfc
  Tb   <- fields$Tbar_sfc; rhoS <- fields$rhoSref_sfc
  U    <- fields$U_L93_sfc; V   <- fields$V_L93_sfc
  zeta <- if (opt == "lcl") fields$zeta_LCL else NULL

  Nm2 <- (g / Tb) * (gam - Gm)
  Cw  <- rhoS * Gm / gam

  # --- etat uniforme moyenne sur le POLYGONE (mask_vec en ordre raster [row,col]).
  #     Les champs sont en [ny,nx] ; on aplati en ordre raster pour appliquer le
  #     masque, coherent avec build_domain_mask (terra::values, row-major).
  vmask <- function(x) {
    xr <- as.vector(t(x))          # [ny,nx] -> ordre raster row-major [ncell]
    mean(xr[mask_vec], na.rm = TRUE)
  }
  U0 <- vmask(U); V0 <- vmask(V); Hw0 <- vmask(Hw); rhoS0 <- vmask(rhoS)
  gam0 <- vmask(gam); Gm0 <- vmask(Gm); Tb0 <- vmask(Tb)
  zeta0 <- if (opt == "lcl") vmask(zeta) else 0
  Nm2_0 <- (g / Tb0) * (gam0 - Gm0)
  Cw0   <- rhoS0 * Gm0 / gam0

  sigma_eps <- coef_sigma_eps * (2 * pi / Lc_wind) * sqrt(U0^2 + V0^2)

  k <- pre$k; l <- pre$l; rho2 <- pre$rho2
  sigma0 <- U0 * k + V0 * l
  s2 <- sigma0^2
  m0 <- complex(length.out = length(sigma0)); dim(m0) <- dim(sigma0)
  prop <- (s2 > 0) & (s2 <= Nm2_0); evan <- (s2 > Nm2_0)
  rho <- sqrt(rho2)
  m0[prop] <- sqrt(Nm2_0 / s2[prop] - 1) * rho[prop] * sign(sigma0[prop])
  m0[evan] <- 1i * sqrt(1 - Nm2_0 / s2[evan]) * rho[evan]

  one_m <- 1 - 1i * m0 * Hw0
  dPm   <- (1i * Hw0) / one_m

  dm_dsig <- complex(length.out = length(m0)); dim(dm_dsig) <- dim(m0)
  dm_dNm2 <- complex(length.out = length(m0)); dim(dm_dNm2) <- dim(m0)
  ok <- prop & (Mod(m0) > 0)
  dm_dsig[ok] <- -Nm2_0 * rho2[ok] / (sigma0[ok]^3 * m0[ok])
  dm_dNm2[ok] <-  rho2[ok] / (2 * s2[ok] * m0[ok])
  if (sigma_eps > 0) {
    dm_dsig[abs(sigma0) < sigma_eps] <- 0+0i
    dm_dNm2[abs(sigma0) < sigma_eps] <- 0+0i
  }

  d_rhoS <- rhoS - rhoS0; d_Hw <- Hw - Hw0
  d_U <- U - U0; d_V <- V - V0; d_gam <- gam - gam0

  # --- facteur de troncature LCL a l'etat de reference (=1 en sfc) ---
  #     T0 = exp(zeta0*(i*m0 - 1/Hw0)) ; coefT = (i*m0 - 1/Hw0) = dT/dzeta / T.
  #     T0 multipliera Phat0 (herite par toutes les derivees) ; d_zeta = ecart
  #     spatial du ddl zeta_LCL (correction de Taylor dans lpm_apply_tau).
  coefT <- 1i * m0 - 1 / Hw0
  if (opt == "lcl") {
    T0     <- exp(zeta0 * coefT); dim(T0) <- dim(m0)
    d_zeta <- zeta - zeta0
  } else {
    T0     <- rep(1+0i, length(m0)); dim(T0) <- dim(m0)
    d_zeta <- NULL
  }

  lp <- 1
  if (!is.null(Lc_wind)) { rho_c <- 2 * pi / Lc_wind; lp <- exp(-(rho2) / rho_c^2) }

  coef_gam_Cw <- (-1/gam0)
  coef_gam_m  <- dPm * dm_dNm2 * (g / Tb0)

  list(pre = pre, sigma0 = sigma0, s2 = s2, m0 = m0, one_m = one_m, dPm = dPm,
       dm_dsig = dm_dsig, Cw0 = Cw0, Hw0 = Hw0, rhoS0 = rhoS0,
       k = k, l = l, lp = lp,
       d_rhoS = d_rhoS, d_Hw = d_Hw, d_U = d_U, d_V = d_V, d_gam = d_gam,
       coef_gam_Cw = coef_gam_Cw, coef_gam_m = coef_gam_m,
       opt = opt, T0 = T0, coefT = coefT, d_zeta = d_zeta,
       spatial = spatial,
       ncell = pre$ncell, FFT_elev = pre$FFT_elev)
}

# ---------------------------------------------------------------------------
#  APPLIQUE un jeu (tau_c, tau_f) au setup : les 6 iFFT + assemblage.
#  Retourne la carte de pluie P en mm/jour (matrice [ny,nx]), NON seuillee ici
#  (le seuillage/Pmin se fait a la fin, apres background).
#
#  hp_mask (optionnel) : masque spectral passe-haut applique a la topographie
#  (meme longueur que FFT_elev). Si fourni, on calcule l'INCREMENT orographique
#  P_LPM(h_fin - h_coarse) = LPM applique a (1-Phi)*h_fin, ou hp_mask = 1-Phi.
#  Par linearite du LPM en la topo, il suffit de multiplier FFT_elev par hp_mask.
#  hp_mask = NULL -> topo complete (comportement standard).
# ---------------------------------------------------------------------------
lpm_apply_tau <- function(S, tau_c, tau_f, hp_mask = NULL) {
  sigma0 <- S$sigma0; s2 <- S$s2; m0 <- S$m0

  FFT_elev <- if (is.null(hp_mask)) S$FFT_elev else S$FFT_elev * hp_mask

  denom0 <- (1 - 1i * m0 * S$Hw0) *
            (1 + 1i * sigma0 * tau_c) * (1 + 1i * sigma0 * tau_f)
  Phat0  <- (S$Cw0 * 1i) * sigma0 * FFT_elev / denom0 * S$T0   # T0=1 en sfc
  Phat0[s2 == 0] <- 0+0i

  dPsig_direct <- 1/sigma0 -
                  (1i * tau_c) / (1 + 1i * sigma0 * tau_c) -
                  (1i * tau_f) / (1 + 1i * sigma0 * tau_f)
  dPsig_direct[s2 == 0] <- 0+0i
  dPsig_tot <- dPsig_direct + S$dPm * S$dm_dsig

  inv <- function(Ph) {
    Ph[!is.finite(Re(Ph)) | !is.finite(Im(Ph))] <- 0+0i
    Ph[s2 == 0] <- 0+0i
    Re(fftwtools::fftw_c2c_2d(Ph, inverse = 1) / S$ncell)
  }

  P0r <- inv(Phat0)                             # terme de reference (toujours calcule)
  sp  <- S$spatial
  P   <- P0r
  # chaque correction de Taylor n'est calculee (iFFT) et ajoutee que si son ddl
  # est active dans `spatial`. Un ddl desactive = parametre maintenu uniforme.
  if (isTRUE(sp["rhoS"]))  P <- P + (P0r / S$rhoS0) * S$d_rhoS      # reutilise P0r (pas d'iFFT)
  if (isTRUE(sp["Hw"]))    P <- P + inv(Phat0 * (1i * m0) / S$one_m) * S$d_Hw
  if (isTRUE(sp["U"]))     P <- P + inv(S$k * Phat0 * dPsig_tot * S$lp) * S$d_U
  if (isTRUE(sp["V"]))     P <- P + inv(S$l * Phat0 * dPsig_tot * S$lp) * S$d_V
  if (isTRUE(sp["gamma"])) P <- P + inv(Phat0 * (S$coef_gam_Cw + S$coef_gam_m * S$lp)) * S$d_gam
  if (identical(S$opt, "lcl") && isTRUE(sp["zeta"]))
    P <- P + inv(S$coefT * Phat0 * S$lp) * S$d_zeta
  86400 * P                                   # mm/jour (NON seuille)
}

# ---------------------------------------------------------------------------
#  PILOTE de calage (design GENERIQUE). Balaye (TAU_C x TAU_F) pour une liste de
#  journees 06-06. Pour chaque couple et jour, stocke aux postes le LPM fin nu,
#  DECOMPOSE en parties POSITIVE et NEGATIVE (P_LPM_pos/neg) -> fdry reconstruit
#  au calage par P_LPM(fdry)=P_LPM_pos+fdry*P_LPM_neg (lineaire, exact).
#
#  Orthogonalisation vs un BACKGROUND generique X (bg_daily), suppose >= 0 partout
#  (ERA5, k*IVT, 20CRv3...). On stocke les produits scalaires du LPM AGREGE avec X,
#  eux aussi decomposes par signe du LPM :
#    s_LB_pos = <A.P_LPM^+, X>,  s_LB_neg = <A.P_LPM^-, X>,  s_BB = <X, X>
#  -> au calage : s_LX(fdry) = s_LB_pos + fdry*s_LB_neg, et c = s_LX/s_BB.
#  Un facteur d'echelle k sur le background (cas k*IVT) se cale EN AVAL sans
#  recalcul : c(k) = (1/k)*s_LX/s_BB (linearite du produit scalaire).
#
#  Plus de voie IVT dediee (s_BP/s_PP/P_basse/B supprimes) : une seule voie,
#  le background est un champ fourni en entree.
#
#  A_agg   : matrice agregation fine->grille background (build_aggregation_matrix)
#  idx_sta : indices postes (build_station_index) pour P_LPM nu
#  mask_vec: masque polygone pour l'etat uniforme
#  bg_daily: background X [maille_bg, jour] (>= 0), sur la grille de A_agg
#  ortho_mask: masque continental sur les mailles couvertes (frac_min)
# ---------------------------------------------------------------------------
calibrate_cube <- function(pre, buf, jours, TAU_C, TAU_F,
                           idx_sta, mask_vec, A_agg, Pmin = 0,
                           Lc_wind = 3000, coef_sigma_eps = 0.25,
                           progress = TRUE, sta_id = NULL,
                           spatial = c(rhoS = TRUE, Hw = TRUE, U = TRUE,
                                       V = TRUE, gamma = TRUE, zeta = TRUE),
                           bg_daily = NULL, ortho_mask = NULL) {
  opt <- if (!is.null(buf$opt)) buf$opt else "sfc"   # variante lue du buffer
  nS <- length(idx_sta); nD <- nrow(jours)
  nC <- length(TAU_C); nF <- length(TAU_F)
  frac <- 3 / 24
  step <- 3 * 3600

  # sorties : LPM aux postes decompose pos/neg [poste, jour, tau_c, tau_f]
  P_LPM_pos <- array(0, dim = c(nS, nD, nC, nF)); P_LPM_neg <- array(0, dim = c(nS, nD, nC, nF))
  covered <- as.vector(Matrix::rowSums(A_agg) > 0)   # mailles background couvertes par le fin

  do_bg <- !is.null(bg_daily)
  if (do_bg) {
    s_LB_pos <- array(0, dim = c(nD, nC, nF)); s_LB_neg <- array(0, dim = c(nD, nC, nF))
    s_BB <- numeric(nD)                            # <X,X> ; independant de tau et fdry
    keep_cov <- if (is.null(ortho_mask)) rep(TRUE, sum(covered)) else ortho_mask[covered]
  }

  avail <- vapply(jours$t_start, day_available,
                  logical(1), buf$path_common, buf$path_sfc)
  n_na <- sum(!avail)
  if (n_na > 0)
    cat(sprintf("ATTENTION : %d journee(s) incomplete(s) marquee(s) NA.\n", n_na))

  pb <- if (progress) utils::txtProgressBar(min = 0, max = nD, style = 3) else NULL
  ncov <- sum(covered)
  for (d in seq_len(nD)) {
    if (!avail[d]) {
      P_LPM_pos[, d, , ] <- NA_real_; P_LPM_neg[, d, , ] <- NA_real_
      if (do_bg) { s_LB_pos[d, , ] <- NA_real_; s_LB_neg[d, , ] <- NA_real_; s_BB[d] <- NA_real_ }
      if (progress) utils::setTxtProgressBar(pb, d); next
    }
    ts0    <- as.numeric(jours$t_start[d])
    stamps <- ts0 + step * (0:7)

    # cumuls journaliers du LPM agrege, parties pos/neg [ncov, nC, nF]
    Pagg_pos <- array(0, dim = c(ncov, nC, nF)); Pagg_neg <- array(0, dim = c(ncov, nC, nF))

    for (ts in stamps) {
      loc    <- locate_timestep(buf, ts)
      fields <- lpm_refine(pre, loc$block, loc$it)
      S <- lpm_setup_timestep(pre, fields, mask_vec,
                              Lc_wind = Lc_wind, coef_sigma_eps = coef_sigma_eps,
                              opt = opt, spatial = spatial)

      for (ic in seq_len(nC)) for (iff in seq_len(nF)) {
        Pv <- as.vector(t(lpm_apply_tau(S, TAU_C[ic], TAU_F[iff])))   # LPM fin nu [ncell]
        Pv_pos <- pmax(Pv, 0); Pv_neg <- pmin(Pv, 0)   # separation AU NIVEAU FIN
        # aux postes (pos/neg)
        P_LPM_pos[, d, ic, iff] <- P_LPM_pos[, d, ic, iff] + frac * Pv_pos[idx_sta]
        P_LPM_neg[, d, ic, iff] <- P_LPM_neg[, d, ic, iff] + frac * Pv_neg[idx_sta]
        # agregation a la grille background, pos/neg
        if (do_bg) {
          Pagg_pos[, ic, iff] <- Pagg_pos[, ic, iff] + frac * as.vector(A_agg %*% Pv_pos)[covered]
          Pagg_neg[, ic, iff] <- Pagg_neg[, ic, iff] + frac * as.vector(A_agg %*% Pv_neg)[covered]
        }
      }
    }

    # produits scalaires avec le background X : lineaires, decomposes par signe LPM
    if (do_bg) {
      x  <- bg_daily[covered, d]
      ok <- is.finite(x) & keep_cov
      s_BB[d] <- sum(x[ok] * x[ok])
      for (ic in seq_len(nC)) for (iff in seq_len(nF)) {
        s_LB_pos[d, ic, iff] <- sum(Pagg_pos[ok, ic, iff] * x[ok])
        s_LB_neg[d, ic, iff] <- sum(Pagg_neg[ok, ic, iff] * x[ok])
      }
    }
    if (progress) utils::setTxtProgressBar(pb, d)
  }
  if (progress) close(pb)
  out <- list(P_LPM_pos = P_LPM_pos, P_LPM_neg = P_LPM_neg,
              TAU_C = TAU_C, TAU_F = TAU_F,
              t_start = jours$t_start, sta_id = sta_id)
  if (do_bg) { out$s_LB_pos <- s_LB_pos; out$s_LB_neg <- s_LB_neg; out$s_BB <- s_BB }
  out
}

# ---------------------------------------------------------------------------
#  CALAGE INCREMENT : explore (tau_c, tau_f, Lc) et stocke, aux postes,
#  l'INCREMENT orographique P_LPM(h_fin - h_coarse) obtenu en filtrant la topo
#  au passe-haut de longueur de coupure Lc (masque spectral, cf. make_hp_mask).
#  Destine au cadre "fond ERA5 + increment" : P_mod = P_era5 + beta*incr.
#
#  Sortie : incr[poste, jour, tau_c, tau_f, Lc] (mm/jour, cumul 06-06, NON seuille),
#  plus B[poste, jour] (module IVT, garde au cas ou) et les axes.
#  hp_mask precalcule une fois par Lc (independant de tau et de l'etat).
# ---------------------------------------------------------------------------
calibrate_cube_incr <- function(pre, buf, jours, TAU_C, TAU_F, LC_HP,
                                idx_sta, mask_vec, I_sta, Pmin = 0,
                                Lc_wind = 3000, coef_sigma_eps = 0.25,
                                progress = TRUE, sta_id = NULL,
                                spatial = c(rhoS = TRUE, Hw = TRUE, U = TRUE,
                                            V = TRUE, gamma = TRUE, zeta = TRUE)) {
  opt <- if (!is.null(buf$opt)) buf$opt else "sfc"
  nS <- length(idx_sta); nD <- nrow(jours)
  nC <- length(TAU_C); nF <- length(TAU_F); nL <- length(LC_HP)
  frac <- 3 / 24; step <- 3 * 3600

  # masques passe-haut precalcules (un par Lc) : independants de l'etat et de tau
  hp_masks <- lapply(LC_HP, function(Lc) make_hp_mask(pre, Lc))

  incr <- array(0, dim = c(nS, nD, nC, nF, nL))
  B    <- matrix(0, nS, nD)

  avail <- vapply(jours$t_start, day_available,
                  logical(1), buf$path_common, buf$path_sfc)
  n_na <- sum(!avail)
  if (n_na > 0)
    cat(sprintf("ATTENTION : %d journee(s) incomplete(s) marquee(s) NA.\n", n_na))

  pb <- if (progress) utils::txtProgressBar(min = 0, max = nD, style = 3) else NULL
  for (d in seq_len(nD)) {
    if (!avail[d]) {
      incr[, d, , , ] <- NA_real_; B[, d] <- NA_real_
      if (progress) utils::setTxtProgressBar(pb, d); next
    }
    ts0 <- as.numeric(jours$t_start[d]); stamps <- ts0 + step * (0:7)
    for (ts in stamps) {
      loc    <- locate_timestep(buf, ts)
      fields <- lpm_refine(pre, loc$block, loc$it)
      S <- lpm_setup_timestep(pre, fields, mask_vec,
                              Lc_wind = Lc_wind, coef_sigma_eps = coef_sigma_eps,
                              opt = opt, spatial = spatial)
      ivu <- as.vector(loc$block$vars$IVTu[, , loc$it])
      ivv <- as.vector(loc$block$vars$IVTv[, , loc$it])
      B[, d] <- B[, d] + frac * as.vector(I_sta %*% sqrt(ivu^2 + ivv^2))

      for (ic in seq_len(nC)) for (iff in seq_len(nF)) for (il in seq_len(nL)) {
        P  <- lpm_apply_tau(S, TAU_C[ic], TAU_F[iff], hp_mask = hp_masks[[il]])
        Pv <- as.vector(t(P))
        incr[, d, ic, iff, il] <- incr[, d, ic, iff, il] + frac * Pv[idx_sta]
      }
    }
    if (progress) utils::setTxtProgressBar(pb, d)
  }
  if (progress) close(pb)
  list(incr = incr, B = B, TAU_C = TAU_C, TAU_F = TAU_F, LC_HP = LC_HP,
       t_start = jours$t_start, sta_id = sta_id)
}

# ---------------------------------------------------------------------------
#  PILOTE PARALLELE : un worker par ANNEE (mclapply). Chaque worker cree son
#  PROPRE buffer (isolation), calcule le sous-cube de son annee, l'ecrit sur
#  disque, et NE RETOURNE QUE le chemin (pas le cube) -> pic memoire maitrise.
#
#  path_common/path_sfc : fonctions chemin (comme make_coarse_buffer).
#  years    : vecteur d'annees a traiter (ex. 2019:2025).
#  out_dir  : repertoire de sortie des sous-cubes (un .rds par annee).
#  mc.cores : nombre de workers simultanes (defaut 6 : pic ~6*3Go sur 32Go).
#  Renvoie la liste des fichiers ecrits. Reprise : les annees deja ecrites
#  (fichier present) sont sautees.
# ---------------------------------------------------------------------------
calibrate_cube_par <- function(pre, path_common, path_sfc, years,
                               TAU_C, TAU_F, idx_sta, mask_vec, A_agg,
                               out_dir, Pmin = 0, Lc_wind = 3000,
                               coef_sigma_eps = 0.25, mc.cores = 6,
                               skip_existing = TRUE, jours_sel = NULL,
                               sta_id = NULL, opt = c("sfc", "lcl"),
                               spatial = c(rhoS = TRUE, Hw = TRUE, U = TRUE,
                                           V = TRUE, gamma = TRUE, zeta = TRUE),
                               bg_daily_pattern = NULL, ortho_mask = NULL) {
  opt <- match.arg(opt)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  nY <- length(years); t0 <- Sys.time()
  message(sprintf("Calage parallele : %d annees, %d workers (variante %s%s).",
                  nY, mc.cores, opt,
                  if (!is.null(bg_daily_pattern)) ", orthog. background" else ""))

  one_year <- function(year) {
    out_f <- file.path(out_dir, sprintf("calib_cube_%d.rds", year))
    if (skip_existing && file.exists(out_f)) {
      message(sprintf("  [%d] deja calcule, saute.", year)); return(out_f)
    }
    # buffer PROPRE a ce worker (isolation entre forks), variante opt
    buf <- make_coarse_buffer(path_common, path_sfc, opt = opt)
    jours <- gen_jours(year, year)                             # 06-06 de l'annee

    # si une pre-selection de jours est fournie (ex. jours pluvieux), ne garder
    # que ceux de l'annee courante. jours_sel accepte un VECTEUR de temps
    # (POSIXct ou numerique) ou un data.frame a colonne t_start (compat).
    if (!is.null(jours_sel)) {
      keep_ts <- as.numeric(if (is.data.frame(jours_sel)) jours_sel$t_start else jours_sel)
      sel <- as.numeric(jours$t_start) %in% keep_ts
      jours <- jours[sel, , drop = FALSE]
      if (nrow(jours) == 0) {
        message(sprintf("  [%d] aucun jour selectionne, saute.", year)); return(NA_character_)
      }
    }

    # background daily EN CHAMP (pour orthogonalisation), si demande (>= 0)
    bg_field <- NULL
    if (!is.null(bg_daily_pattern)) {
      f_bg <- bg_daily_pattern(year)
      if (file.exists(f_bg)) bg_field <- read_era5_daily_field(f_bg, jours$t_start)
      else message(sprintf("  [%d] fichier background daily absent, orthog. ignoree.", year))
    }

    # barre interne coupee (progress=FALSE) : illisible avec forks concurrents.
    res <- calibrate_cube(pre, buf, jours, TAU_C, TAU_F,
                          idx_sta, mask_vec, A_agg, Pmin = Pmin,
                          Lc_wind = Lc_wind, coef_sigma_eps = coef_sigma_eps,
                          progress = FALSE, sta_id = sta_id, spatial = spatial,
                          bg_daily = bg_field, ortho_mask = ortho_mask)
    saveRDS(res, out_f)
    rm(res, buf); gc()                                         # libere avant de rendre
    el <- as.numeric(Sys.time() - t0, units = "mins")
    message(sprintf("  [%d] termine  (%.1f min ecoulees)", year, el))
    out_f                                                      # retourne le CHEMIN, pas le cube
  }

  # mc.preschedule = FALSE : distribution DYNAMIQUE des annees (une tache prise
  # des qu'un worker se libere) -> meilleur equilibrage quand les annees ont des
  # charges inegales (nb de jours pluvieux variable). Cout : un fork par annee,
  # negligeable devant le temps de calcul d'une annee.
  files <- parallel::mclapply(years, one_year, mc.cores = mc.cores,
                              mc.preschedule = FALSE)
  message(sprintf("Calage termine en %.1f min.",
                  as.numeric(Sys.time() - t0, units = "mins")))
  unlist(files)
}

# ---------------------------------------------------------------------------
#  Relecture / recombinaison des sous-cubes annuels en un seul objet.
#  Colle P_LPM_pos/neg (poste-jour) et s_LB_pos/neg (jour) ; concatene t_start.
# ---------------------------------------------------------------------------
combine_year_cubes <- function(files) {
  files <- files[!is.na(files)]
  parts <- lapply(files, readRDS)
  ab2 <- function(nm) do.call(function(...) abind::abind(..., along = 2),
                              lapply(parts, `[[`, nm))   # recolle sur l'axe jour (poste-jour)
  ab1 <- function(nm) do.call(function(...) abind::abind(..., along = 1),
                              lapply(parts, `[[`, nm))   # recolle sur l'axe jour (jour-...)
  out <- list(P_LPM_pos = ab2("P_LPM_pos"), P_LPM_neg = ab2("P_LPM_neg"),
              TAU_C = parts[[1]]$TAU_C, TAU_F = parts[[1]]$TAU_F,
              t_start = do.call(c, lapply(parts, `[[`, "t_start")),
              sta_id = parts[[1]]$sta_id)
  # orthogonalisation background (si les sous-cubes contiennent s_LB_pos/neg/s_BB)
  if (!is.null(parts[[1]]$s_LB_pos)) {
    out$s_LB_pos <- ab1("s_LB_pos"); out$s_LB_neg <- ab1("s_LB_neg")
    out$s_BB <- do.call(c, lapply(parts, `[[`, "s_BB"))
  }
  out
}

# ---------------------------------------------------------------------------
#  PILOTE PARALLELE pour l'INCREMENT (3e axe Lc). Un worker par annee ; ecrit
#  incr_cube_YYYY.rds. Meme logique que calibrate_cube_par.
# ---------------------------------------------------------------------------
calibrate_cube_incr_par <- function(pre, path_common, path_sfc, years,
                                    TAU_C, TAU_F, LC_HP, idx_sta, mask_vec, I_sta,
                                    out_dir, Pmin = 0, Lc_wind = 3000,
                                    coef_sigma_eps = 0.25, mc.cores = 6,
                                    skip_existing = TRUE, jours_sel = NULL,
                                    sta_id = NULL, opt = c("sfc", "lcl"),
                                    spatial = c(rhoS = TRUE, Hw = TRUE, U = TRUE,
                                                V = TRUE, gamma = TRUE, zeta = TRUE)) {
  opt <- match.arg(opt)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  nY <- length(years); t0 <- Sys.time()
  message(sprintf("Calage increment : %d annees, %d workers (variante %s, %d Lc).",
                  nY, mc.cores, opt, length(LC_HP)))

  one_year <- function(year) {
    out_f <- file.path(out_dir, sprintf("incr_cube_%d.rds", year))
    if (skip_existing && file.exists(out_f)) {
      message(sprintf("  [%d] deja calcule, saute.", year)); return(out_f)
    }
    buf <- make_coarse_buffer(path_common, path_sfc, opt = opt)
    jours <- gen_jours(year, year)
    if (!is.null(jours_sel)) {
      keep_ts <- as.numeric(if (is.data.frame(jours_sel)) jours_sel$t_start else jours_sel)
      jours <- jours[as.numeric(jours$t_start) %in% keep_ts, , drop = FALSE]
      if (nrow(jours) == 0) {
        message(sprintf("  [%d] aucun jour selectionne, saute.", year)); return(NA_character_)
      }
    }
    res <- calibrate_cube_incr(pre, buf, jours, TAU_C, TAU_F, LC_HP,
                               idx_sta, mask_vec, I_sta, Pmin = Pmin,
                               Lc_wind = Lc_wind, coef_sigma_eps = coef_sigma_eps,
                               progress = FALSE, sta_id = sta_id, spatial = spatial)
    saveRDS(res, out_f)
    rm(res, buf); gc()
    message(sprintf("  [%d] termine  (%.1f min)", year,
                    as.numeric(Sys.time() - t0, units = "mins")))
    out_f
  }
  files <- parallel::mclapply(years, one_year, mc.cores = mc.cores,
                              mc.preschedule = FALSE)
  message(sprintf("Increment termine en %.1f min.",
                  as.numeric(Sys.time() - t0, units = "mins")))
  unlist(files)
}

# recombine les sous-cubes d'increment sur l'axe JOUR (dim 2 de incr, cols de B)
combine_incr_cubes <- function(files) {
  files <- files[!is.na(files)]
  parts <- lapply(files, readRDS)
  incr <- do.call(function(...) abind::abind(..., along = 2),
                  lapply(parts, `[[`, "incr"))
  B    <- do.call(cbind, lapply(parts, `[[`, "B"))
  t_start <- do.call(c, lapply(parts, `[[`, "t_start"))
  list(incr = incr, B = B, TAU_C = parts[[1]]$TAU_C, TAU_F = parts[[1]]$TAU_F,
       LC_HP = parts[[1]]$LC_HP, t_start = t_start, sta_id = parts[[1]]$sta_id)
}


# ===========================================================================
#  BACKGROUND IVT : matrices A (agregation), I_sta (ERA5->postes), calage beta
# ===========================================================================

# ---------------------------------------------------------------------------
#  A : matrice creuse [n_era5, ncell_fine] de MOYENNAGE par boite (fine -> ERA5).
#  Chaque maille ERA5 = moyenne des mailles fines dont le centre tombe dedans.
#  (Le domaine fin est inclus dans ERA5, donc toute maille fine a une mere ERA5.)
# ---------------------------------------------------------------------------
build_aggregation_matrix <- function(pre) {
  lon <- pre$lon_ref; lat <- pre$lat_ref
  nlon <- length(lon); nlat <- length(lat); nera5 <- nlon * nlat
  ncell <- pre$ncell

  lon_inc <- lon[2] > lon[1]; lat_inc <- lat[2] > lat[1]
  lon_s <- if (lon_inc) lon else rev(lon)
  lat_s <- if (lat_inc) lat else rev(lat)
  perm_lon <- if (lon_inc) seq_len(nlon) else rev(seq_len(nlon))
  perm_lat <- if (lat_inc) seq_len(nlat) else rev(seq_len(nlat))

  # centres des mailles fines (L93) -> WGS84
  xy <- terra::xyFromCell(pre$template, seq_len(ncell))     # ordre raster row-major
  ll <- terra::project(xy, from = "EPSG:2154", to = "EPSG:4326")
  lon_c <- ll[, 1]; lat_c <- ll[, 2]

  # bords de mailles ERA5 (grille reguliere 0.25) : maille CONTENANTE par arrondi.
  dlon <- median(diff(lon_s)); dlat <- median(diff(lat_s))
  ii <- round((lon_c - lon_s[1]) / dlon) + 1L               # indice trie 1..nlon
  jj <- round((lat_c - lat_s[1]) / dlat) + 1L
  ii <- pmin(pmax(ii, 1L), nlon); jj <- pmin(pmax(jj, 1L), nlat)

  # colonne ERA5 aplatie [lon,lat] (as.vector col-major), avec remap du sens d'axe
  era5_col <- (perm_lat[jj] - 1L) * nlon + perm_lon[ii]

  # moyennage : poids 1/n_fine par maille ERA5 (compte des mailles fines par mere)
  cnt <- tabulate(era5_col, nbins = nera5)
  w   <- 1 / cnt[era5_col]
  A <- Matrix::sparseMatrix(i = era5_col, j = seq_len(ncell), x = w,
                            dims = c(nera5, ncell))
  A                                                          # [nera5, ncell]
}

# ---------------------------------------------------------------------------
#  I : matrice creuse [n_postes, n_era5] d'interpolation BILINEAIRE (ERA5 -> postes).
#  Meme logique que build_weight_matrix, mais cibles = postes (sf POINTS).
# ---------------------------------------------------------------------------
# ---------------------------------------------------------------------------
#  Interpolation bilineaire d'une grille reguliere lon/lat (WGS84) vers des
#  postes (sf POINTS). Generique : marche pour la grille des lpm_common comme
#  pour toute autre grille ERA5 (ex. precip journaliere, emprise differente).
#  Renvoie la matrice creuse [n_postes, nlon*nlat] (ordre ERA5 [lon,lat] aplati).
# ---------------------------------------------------------------------------
build_bilinear_to_stations <- function(lon, lat, pts_sf) {
  nlon <- length(lon); nlat <- length(lat)
  lon_inc <- lon[2] > lon[1]; lat_inc <- lat[2] > lat[1]
  lon_s <- if (lon_inc) lon else rev(lon)
  lat_s <- if (lat_inc) lat else rev(lat)
  perm_lon <- if (lon_inc) seq_len(nlon) else rev(seq_len(nlon))
  perm_lat <- if (lat_inc) seq_len(nlat) else rev(seq_len(nlat))

  pts_v <- terra::vect(pts_sf)
  ll <- terra::crds(terra::project(pts_v, "EPSG:4326"))
  lon_c <- ll[, 1]; lat_c <- ll[, 2]; npt <- length(lon_c)

  i <- findInterval(lon_c, lon_s); j <- findInterval(lat_c, lat_s)
  inside <- i >= 1 & i < nlon & j >= 1 & j < nlat
  i <- pmin(pmax(i, 1), nlon - 1L); j <- pmin(pmax(j, 1), nlat - 1L)
  tx <- (lon_c - lon_s[i]) / (lon_s[i + 1] - lon_s[i])
  ty <- (lat_c - lat_s[j]) / (lat_s[j + 1] - lat_s[j])
  tx <- pmin(pmax(tx, 0), 1); ty <- pmin(pmax(ty, 0), 1)

  flat <- function(ii, jj) (perm_lat[jj] - 1L) * nlon + perm_lon[ii]
  c00 <- flat(i,     j    ); w00 <- (1 - tx) * (1 - ty)
  c10 <- flat(i + 1, j    ); w10 <- tx * (1 - ty)
  c01 <- flat(i,     j + 1); w01 <- (1 - tx) * ty
  c11 <- flat(i + 1, j + 1); w11 <- tx * ty

  rows <- rep(seq_len(npt), 4); cols <- c(c00, c10, c01, c11)
  vals <- c(w00, w10, w01, w11); ok <- rep(inside, 4)
  list(I = Matrix::sparseMatrix(i = rows[ok], j = cols[ok], x = vals[ok],
                                dims = c(npt, nlon * nlat)),
       inside = inside)                                       # inside : postes dans l'emprise
}

# ERA5 (grille des lpm_common) -> postes : cas particulier de la fonction generique.
build_era5_to_stations <- function(pre, pts_sf) {
  build_bilinear_to_stations(pre$lon_ref, pre$lat_ref, pts_sf)$I
}

# ---------------------------------------------------------------------------
#  build_production_plan : croise un calage `opt` avec les classes des jours de
#  production pour produire un PLAN journalier [t_start, tau_c, tau_f, fdry, beta, c].
#  c = NA -> recalculer chaque jour (c_mode "daily") ; sinon valeur fixe (global/class).
#  Isole toute la logique "classe -> parametres" : la production, elle, est agnostique.
# ---------------------------------------------------------------------------
build_production_plan <- function(opt, dates, classe_prod = NULL, c_mode = NULL) {
  if (is.null(c_mode)) c_mode <- if (!is.null(opt$c_mode)) opt$c_mode else "daily"
  pc <- opt$per_class
  if (is.null(pc))
    pc <- list("1" = list(tau_c = opt$tau_c, tau_f = opt$tau_f, fdry = opt$fdry, beta = opt$beta))
  nD <- length(dates)
  cl <- if (!is.null(classe_prod)) as.character(classe_prod) else rep(names(pc)[1], nD)
  # garde-fou : si un parametre varie entre classes, les classes sont requises
  varie <- function(nm) length(unique(vapply(pc, function(p) p[[nm]], numeric(1)))) > 1
  if (is.null(classe_prod) && (varie("tau_c")||varie("tau_f")||varie("fdry")||varie("beta")))
    stop("parametres par classe : fournir classe_prod (classe de chaque jour de production).")
  getp <- function(nm) vapply(seq_len(nD), function(d) {
    p <- pc[[cl[d]]]; if (is.null(p)) stop(sprintf("classe '%s' (jour %d) absente du calage.", cl[d], d))
    p[[nm]] }, numeric(1))
  cval <- if (c_mode == "daily") rep(NA_real_, nD)
          else if (c_mode == "global") rep(as.numeric(opt$c)[1], nD)
          else vapply(seq_len(nD), function(d) {                  # class
            v <- as.numeric(opt$c[cl[d]]); if (length(v) == 0 || is.na(v)) 0 else v }, numeric(1))
  data.frame(t_start = dates, tau_c = getp("tau_c"), tau_f = getp("tau_f"),
             fdry = getp("fdry"), beta = getp("beta"), c = cval)
}

# cle d'agregation temporelle d'un vecteur de dates (Date ou POSIXct).
.aggreg_key <- function(dates, agreg) {
  d <- as.Date(dates)
  switch(agreg,
    daily    = format(d, "%Y-%m-%d"),
    monthly  = format(d, "%Y-%m"),
    seasonal = paste0(format(d, "%Y"), "-S", (as.integer(format(d, "%m")) - 1) %/% 3 + 1),
    annual   = format(d, "%Y"),
    stop("agreg doit etre daily, monthly, seasonal ou annual."))
}

# ---------------------------------------------------------------------------
#  produce_gridded_year : applique un PLAN journalier sur la grille fine et ecrit
#  un NetCDF. `agreg` controle le STOCKAGE : "daily" ecrit chaque jour ; "monthly"
#  /"seasonal"/"annual" ecrivent le CUMUL par periode. Le calcul reste journalier
#  (seuillage max(.,0) par jour) et n'est accumule qu'ENSUITE -> agregation exacte
#  (le max ne commute pas avec la somme : on somme des P_mod deja seuilles).
#  Agnostique aux classes : ne lit que le plan (parametres par jour).
# ---------------------------------------------------------------------------
produce_gridded_year <- function(plan, pre, buf, jours, mask_vec, A_agg, I_fine,
                                 bg_daily, ortho_mask, out_nc,
                                 agreg = "daily", clamp_c = c(0, 1),
                                 Lc_wind = 3000, coef_sigma_eps = 0.25,
                                 spatial = NULL, opt_variant = NULL, progress = TRUE) {
  clampc <- function(x) if (is.null(clamp_c)) x else pmin(pmax(x, clamp_c[1]), clamp_c[2])
  if (is.null(spatial)) stop("`spatial` requis (vecteur des drapeaux LPM).")
  if (is.null(opt_variant)) opt_variant <- if (!is.null(buf$opt)) buf$opt else "sfc"
  stopifnot(nrow(plan) == nrow(jours))

  ncell <- length(mask_vec)
  nD    <- nrow(jours)
  frac  <- 3 / 24; step <- 3 * 3600
  covered  <- as.vector(Matrix::rowSums(A_agg) > 0)
  keep_cov <- if (is.null(ortho_mask)) rep(TRUE, sum(covered)) else ortho_mask[covered]

  # periodes de sortie : cle d'agregation par jour, periodes uniques ordonnees.
  key    <- .aggreg_key(jours$t_start, agreg)
  periods <- unique(key)                                  # dans l'ordre chronologique
  nP     <- length(periods)
  p_of_d <- match(key, periods)                           # jour -> index periode
  # timestamp representatif de chaque periode : 1er jour de la periode
  p_time <- as.numeric(jours$t_start[match(periods, key)])
  is_daily <- (agreg == "daily")
  units <- if (is_daily) "mm/day" else sprintf("mm/%s", sub("ly$","", agreg))  # month/season/annual

  # buffer de sortie : [ncell_fine, nP]. En daily nP=nD (un champ par jour) ;
  # sinon on accumule les cumuls par periode. Modeste en RAM (nP petit : <=365).
  acc <- matrix(0, nrow = ncell, ncol = nP)
  filled <- logical(nP)                                   # periode a-t-elle un jour valide ?

  pb <- if (progress) utils::txtProgressBar(min = 0, max = nD, style = 3) else NULL
  for (d in seq_len(nD)) {
    tau_c <- plan$tau_c[d]; tau_f <- plan$tau_f[d]; fdry <- plan$fdry[d]; beta <- plan$beta[d]
    ts0 <- as.numeric(jours$t_start[d]); stamps <- ts0 + step * (0:7)
    PL_fine <- numeric(ncell); Pagg <- numeric(sum(covered)); ok_day <- TRUE
    for (ts in stamps) {
      # ne capture QUE le cas legitime "timestamp introuvable" (donnee manquante
      # -> on saute le jour). Toute autre erreur (buffer malforme, fonction
      # cassee) DOIT remonter : sinon des fichiers vides seraient produits sans
      # signal (piege du tryCatch aveugle).
      loc <- tryCatch(locate_timestep(buf, ts),
                      error = function(e) {
                        if (grepl("introuvable", conditionMessage(e))) NULL
                        else stop(e)
                      })
      if (is.null(loc)) { ok_day <- FALSE; break }
      fields <- lpm_refine(pre, loc$block, loc$it)
      S <- lpm_setup_timestep(pre, fields, mask_vec, Lc_wind = Lc_wind,
                              coef_sigma_eps = coef_sigma_eps, opt = opt_variant, spatial = spatial)
      Pv <- as.vector(t(lpm_apply_tau(S, tau_c, tau_f)))
      Pv[Pv < 0] <- fdry * Pv[Pv < 0]
      PL_fine <- PL_fine + frac * Pv
      Pagg    <- Pagg    + frac * as.vector(A_agg %*% Pv)[covered]
    }
    if (!ok_day) { if (progress) utils::setTxtProgressBar(pb, d); next }

    # coefficient c du jour : recalcule si NA (daily), sinon valeur du plan.
    if (is.na(plan$c[d])) {
      x  <- bg_daily[covered, d]; okc <- is.finite(x) & keep_cov
      sBB <- sum(x[okc] * x[okc])
      cj  <- if (sBB > 0) clampc(sum(Pagg[okc] * x[okc]) / sBB) else 0
    } else cj <- plan$c[d]

    bg_fin <- as.vector(I_fine %*% bg_daily[, d])
    Pmod   <- pmax(bg_fin + beta * (PL_fine - cj * bg_fin), 0)   # seuillage JOURNALIER
    Pmod[!is.finite(Pmod)] <- 0                                  # jour manquant -> 0 au cumul
    acc[, p_of_d[d]] <- acc[, p_of_d[d]] + Pmod                  # accumulation dans la periode
    filled[p_of_d[d]] <- TRUE
    if (progress) utils::setTxtProgressBar(pb, d)
  }
  if (progress) close(pb)

  # --- ecriture NetCDF sur la GRILLE FINE (pre$template = MNT), pas la grille
  # grossiere pre$lon_ref/lat_ref. Les champs acc/Pmod sont en ordre raster
  # (row-major : cell 1 = coin haut-gauche, lignes de haut en bas). On construit
  # les axes depuis le template et on reordonne en [x, y] (y croissant) pour un
  # NetCDF standard. ---
  tmpl  <- pre$template
  xfin  <- terra::xFromCol(tmpl, seq_len(terra::ncol(tmpl)))   # x uniques (gauche->droite)
  yfin  <- terra::yFromRow(tmpl, seq_len(terra::nrow(tmpl)))   # y uniques (haut->bas, decroissant)
  nx <- length(xfin); ny <- length(yfin)
  stopifnot(nx * ny == ncell)
  # acc est en ordre RASTER row-major : cellule 1 = coin haut-gauche (nord), on
  # parcourt chaque ligne (x croissant) en descendant (y decroissant). On garde
  # cet ordre : axe y DECROISSANT (nord->sud, convention terra) et matrice [x, y]
  # sans inversion de lignes. Ainsi field[,j] correspond a la ligne j (du nord).
  to_xy <- function(v) t(matrix(v, nrow = ny, ncol = nx, byrow = TRUE))  # [x, y] y decroissant
  dimx <- ncdf4::ncdim_def("x", "m", xfin, longname = "projection_x_coordinate")
  dimy <- ncdf4::ncdim_def("y", "m", yfin, longname = "projection_y_coordinate")  # decroissant
  dimt <- ncdf4::ncdim_def("time", "seconds since 1970-01-01", p_time, unlim = TRUE)
  var  <- ncdf4::ncvar_def("pr", units, list(dimx, dimy, dimt),
                           missval = -9999, prec = "float", compression = 4)
  # CRS en convention CF : variable scalaire 'crs' portant le WKT/EPSG, referencee
  # par 'pr' via grid_mapping -> terra et QGIS lisent le georeferencement seuls.
  crsvar <- ncdf4::ncvar_def("crs", "", list(), prec = "integer")
  nc <- ncdf4::nc_create(out_nc, list(var, crsvar))
  ncdf4::ncatt_put(nc, "pr", "long_name",
                   if (is_daily) "daily precipitation" else sprintf("%s precipitation total", sub("ly$","", agreg)))
  ncdf4::ncatt_put(nc, "pr", "aggregation", agreg)
  ncdf4::ncatt_put(nc, "pr", "grid_mapping", "crs")
  ncdf4::ncatt_put(nc, "crs", "grid_mapping_name", "lambert_conformal_conic")
  ncdf4::ncatt_put(nc, "crs", "spatial_ref", "EPSG:2154")   # terra lit cet attribut
  ncdf4::ncatt_put(nc, "crs", "crs_wkt", sf::st_crs(2154)$wkt)
  ncdf4::ncatt_put(nc, 0, "crs", "EPSG:2154 (RGF93 / Lambert-93)")
  for (p in seq_len(nP)) {
    field <- if (filled[p]) to_xy(acc[, p]) else matrix(-9999, nx, ny)
    ncdf4::ncvar_put(nc, var, field, start = c(1, 1, p), count = c(nx, ny, 1))
  }
  ncdf4::nc_close(nc)
  out_nc
}
