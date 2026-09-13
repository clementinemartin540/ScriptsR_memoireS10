# ============================================================================
#  calib_build_cube.R  -  CONSTRUCTION des cubes de calage (partie lourde).
#
#  Calcule, par annee et en parallele, le cube [poste, jour, tau_c, tau_f] en
#  DESIGN pos/neg : P_LPM_pos/neg, P_bas_pos/neg (parties positive/negative du
#  LPM fin), s_LE_pos/neg + s_EE pour l'orthogonalisation continentale. fdry
#  n'est plus un axe stocke : il est reconstruit au calage par X = pos + fdry*neg
#  (lineaire, exact), ce qui rend le cube ~3.5x plus petit et fdry continu.
#  s_LB_pos/neg + s_BB pour l'orthogonalisation vs un background generique.
#  Ecrit un .rds par annee dans OUT_DIR. Le CALAGE (choix tau/fdry/beta sur obs)
#  se fait ensuite avec calib_fit.R, sans recalcul.
#
#  A relancer seulement si un reglage DEFINISSANT le cube change (OPT, SPATIAL,
#  TAU_C/F, FDRY_GRID, jours, FRAC_MIN). skip_existing saute les annees deja
#  calculees.
# ============================================================================

source("C:/Users/cmartin/ownCloud/scriptsR/calib_ortho/calib_cube_auxDEF.R")
library(ncdf4); library(sf); library(terra); library(lubridate); library(dplyr)

# ------------------------------------------------------------------- chemins
ERA5_DIR       <- "C:/Users/cmartin/donnees/LPM"  #données d'entrée tri-horaires, intégrales des différentes variables
ERA5_DAILY_DIR <- "C:/Users/cmartin/ownCloud/daily"  #précipitations journalières pour le background
MNT         <- "C:/Users/cmartin/ownCloud/scriptsR/calib/bdalti2000m_LPM_smooth.tif" 
DOMAIN      <- "C:/Users/cmartin/ownCloud/scriptsR/shape/DEP_GrandOuest.gpkg"
POSTES      <- "C:/Users/cmartin/ownCloud/scriptsR/calib/stations_MF.gpkg"
RR_CSV      <- "C:/Users/cmartin/ownCloud/scriptsR/QUOT/RR_data_light.csv"   #données pluviométriques de MétéoFrance à partir de 1940
CUBE_ROOT   <- "C:/Users/cmartin/ownCloud/scriptsR/calib_ortho"

# ================= REGLAGES DU CUBE (a garder IDENTIQUES dans calib_fit.R) ===
YEAR0 <- 1993
YEAR1 <- 2008
OPT      <- "sfc"
SPATIAL  <- c(rhoS = TRUE, Hw = TRUE, U = TRUE, V = TRUE, gamma = FALSE, zeta = FALSE)
TAU_C <- exp(seq(log(1000), log(6000), length.out = 8))
TAU_F <- exp(seq(log(1000), log(6000), length.out = 8))
FRAC_MIN  <- 0.5                        # seuil de continentalite (masque orthog. background)
# NB : fdry n'est PLUS un axe du cube (design pos/neg) -> la grille fdry est un
# reglage de CALAGE (FDRY_GRID dans calib_fit.R), reconstruite lineairement.
# ============================================================================
MC_CORES <- 1

path_common <- function(y) sprintf("%s/inputs/lpm_common_%d.nc", ERA5_DIR, y) #accès aux données 
path_sfc    <- function(y) sprintf("%s/inputs/lpm_%s_%d.nc", ERA5_DIR, OPT, y)
# background daily EN CHAMP (>= 0), grille de A_agg. Ici ERA5 ; remplacable par
# 20CRv3, etc. (meme format [maille, jour]).
bg_daily_pattern <- function(y) sprintf("%s/total_precip_daily_06utc_ouest_%d.nc",
                                        ERA5_DAILY_DIR, y)

# dossier des cubes : depend de OPT et SPATIAL (ce qui definit le cube).
# spatial_tag() vient de calib_cube_auxDEF.R. NB : le meme chemin doit etre
# reconstruit a l'identique dans calib_fit.R (memes OPT et SPATIAL).
OUT_DIR <- path.expand(sprintf("%s/calib_cubes_ortho_%s_%s",
                               CUBE_ROOT, OPT, spatial_tag(SPATIAL)))

# ------------------------------------------------------------- preparation
Elev <- terra::rast(MNT)  #prépare tout ce qui est lié au relief, FFT
poly <- sf::st_read(DOMAIN, quiet = TRUE) |> sf::st_union() |> sf::st_simplify(dTolerance = 2000)
pre  <- lpm_precompute(Elev, path_common(YEAR0))
postes <- sf::st_read(POSTES, quiet = TRUE)

mask_vec <- build_domain_mask(pre, poly)
A_agg    <- build_aggregation_matrix(pre)
idx_sta  <- build_station_index(pre, postes)
ortho_mask <- build_continental_mask(A_agg, mask_vec, frac_min = FRAC_MIN)
cat(sprintf("Mailles ERA5 continentales pour orthog. : %d\n", sum(ortho_mask)))
cat(sprintf("Cubes -> %s\n", OUT_DIR))

RR <- read.table(RR_CSV, header = TRUE, sep = ";")
RR[RR < 0] <- NA
RR <- RR[substr(RR$Date, 1, 4) %in% as.character(seq(YEAR0, YEAR1)), ]

# jours a caler : VECTEUR de temps (POSIXct 06 UTC). Ici via select_wet_days sur
# un evenement ; peut aussi venir d'un code externe (classif) qui fournit
# directement un vecteur de dates -> calibrate_cube_par accepte un vecteur nu.

jours_s <- readRDS("C:/Users/cmartin/ownCloud/scriptsR/jours_calib_random.RDS")
jours_sel0 <- jours_s[format(jours_s, "%Y") %in% as.character(seq(YEAR0, YEAR1))]
jours_sel <- data.frame(t_start = as.POSIXct(paste(jours_sel0, "06:00:00"),format = "%Y-%m-%d %H:%M:%S", tz = "UTC"))

cat(sprintf("%d jours retenus.\n", length(jours_sel$t_start)))

# ------------------------------------------------------- calcul des cubes
years <- seq(YEAR0, YEAR1)
files <- calibrate_cube_par(pre, path_common, path_sfc, years,
                            TAU_C, TAU_F, idx_sta, mask_vec, A_agg,
                            out_dir = OUT_DIR, jours_sel = jours_sel,
                            Pmin = 0, Lc_wind = 3000, coef_sigma_eps = 0.25,
                            mc.cores = MC_CORES, sta_id = postes$idMF,
                            opt = OPT, spatial = SPATIAL,
                            bg_daily_pattern = bg_daily_pattern,
                            ortho_mask = ortho_mask)
cat(sprintf("Sous-cubes ecrits (%d annees).\n", length(files)))
