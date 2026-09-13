# ============================================================================
#  produce_gridded_driver.R  -  PRODUCTION du forcage downscalé sur grille FINE.
#
#  Rejoue le modèle sur la grille du MNT au triplet optimal d'un calage,
#  pour une plage d'annees, et ecrit un NetCDF par annee [lon, lat, jour].
#  Le background ERA5 est interpole bilineairement a la grille fine ; c_j est
#  recalcule chaque jour (c_mode="daily") ou lu du calage (global/class).
#
#  Prerequis : inputs LPM (lpm_common + lpm_<opt>) et background ERA5 daily
#  disponibles pour les annees demandees. Reglages OPT/SPATIAL IDENTIQUES au
#  calage qui a produit `opt` (sinon le LPM rejoue differe du LPM cale).
# ============================================================================

library(terra); library(sf); library(ncdf4); library(Matrix); library(moistR)
source("C:/Users/cmartin/ownCloud/partage_LPM_20260824/calib_ortho/calib_cube_auxDEF.R")

# ----------------------------------------------------------------- chemins
ERA5_DIR   <- path.expand("C:/Users/cmartin/donnees/LPM")
ERA5_DAILY_DIR <- path.expand("C:/Users/cmartin/ownCloud/daily")
MNT        <- path.expand("C:/Users/cmartin/ownCloud/scriptsR/calib/bdalti2000m_LPM_smooth.tif")
DOMAIN     <- path.expand("C:/Users/cmartin/ownCloud/scriptsR/shape/DEP_GrandOuest.gpkg")

CUBE_ROOT <- "C:/Users/cmartin/ownCloud/partage_LPM_20260824/calib_ortho"
OPT       <- "lcl"
SPATIAL   <- c(rhoS = TRUE, Hw = TRUE, U = TRUE, V = TRUE, gamma = FALSE, zeta = TRUE)
TAG       <-  "tau-cls_fdry-glob_beta-cls-fix5-0_c-daily_clamp-1.1.5_seed1"
CV        <- TRUE   # TRUE = produire depuis un calage CV (optimum_ortho_cv_) ;
# FALSE = calage complet (optimum_ortho_).

# ------------------------------------------------------------------ arborescence
# MEME modele que les cubes/calage : <CUBE_ROOT>/calib_cubes_<OPT>_<spatial>_ortho.
# L'optimum du calage y est lu. Les sorties de production vont dans une
# arborescence PARALLELE indexee par OPT+SPATIAL+TAG, pour qu'une production soit
# tracable jusqu'a son calage exact (pas de melange entre versions).
CALIB_DIR <- path.expand(sprintf("%s/calib_cubes_ortho_%s_%s",
                                 CUBE_ROOT, OPT, spatial_tag(SPATIAL)))
opt_prefix <- if (isTRUE(CV)) "optimum_ortho_cv_" else "optimum_ortho_"
OPT_RDS   <- file.path(CALIB_DIR, sprintf("%s%s.rds", opt_prefix, TAG))
if (!file.exists(OPT_RDS))
  stop(sprintf("optimum introuvable :\n  %s\n(verifier OPT, SPATIAL, TAG, CV)", OPT_RDS))

# sorties : <PROD_ROOT>/prod_<OPT>_<spatial>/<TAG>/  -> un dossier par calage.
PROD_ROOT <- path.expand("C:/Users/cmartin/ownCloud/partage_LPM_20260824/production_2km")
PROD_DIR  <- file.path(PROD_ROOT, sprintf("prod_%s_%s", OPT, spatial_tag(SPATIAL)), TAG)
dir.create(PROD_DIR, showWarnings = FALSE, recursive = TRUE)
cat(sprintf("Calage : %s\nSorties : %s\n", OPT_RDS, PROD_DIR))

# =================== REGLAGES (a garder COHERENTS avec le calage) ===========
YEARS   <- 1991:2020
FRAC_MIN <- 0.5                         # masque continental pour c_j (idem calage)
CLAMP_C  <- c(-1, 1.5)                  # bornage c_j : DOIT matcher le calage (cf. TAG)
C_MODE   <- NULL                        # NULL = lit opt$c_mode ; sinon force
AGREG    <- "daily"                     # "daily" | "monthly" | "seasonal" | "annual"
# stockage agrege (cumul par periode) ; le
# calcul reste journalier (seuillage exact).
# CLASSES_PROD requis des qu'un parametre varie par classe (tau/fdry/beta en mode
# class) OU si c_mode="class". Classe de CHAQUE jour de production (classification
# TW des geopotentiels de la periode contre les centroides du calage).
CLASSES_PROD_RDS <- "C:/Users/cmartin/ownCloud/scriptsR/classif_jours_60_25.RDS"
# ============================================================================

opt <- readRDS(OPT_RDS)
cat(sprintf("Calage charge : c_mode=%s\n", opt$c_mode))
cat("  parametres par classe :\n")
for (k in names(opt$per_class)) {
  p <- opt$per_class[[k]]
  cat(sprintf("    classe %s : tau_c=%.0f tau_f=%.0f fdry=%.2f beta=%.3f\n",
              k, p$tau_c, p$tau_f, p$fdry, p$beta))
}

# ------------------------------------------------------------- preparation
path_common <- function(y) sprintf("%s/inputs/lpm_common_%d.nc", ERA5_DIR, y)
path_sfc    <- function(y) sprintf("%s/inputs/lpm_%s_%d.nc", ERA5_DIR, OPT, y)
bg_daily_pattern <- function(y) sprintf("%s/total_precip_daily_06utc_ouest_%d.nc",
                                        ERA5_DAILY_DIR, y)

Elev <- terra::rast(MNT)
poly <- sf::st_read(DOMAIN, quiet = TRUE) |> sf::st_union() |> sf::st_simplify(dTolerance = 2000)
pre  <- lpm_precompute(Elev, path_common(YEARS[1]))

mask_vec   <- build_domain_mask(pre, poly)
A_agg      <- build_aggregation_matrix(pre)
ortho_mask <- build_continental_mask(A_agg, mask_vec, frac_min = FRAC_MIN)
cat(sprintf("Grille fine : %d x %d = %d mailles ; %d mailles ERA5 continentales.\n",
            length(pre$lon_ref), length(pre$lat_ref),
            length(pre$lon_ref) * length(pre$lat_ref), sum(ortho_mask)))

# interpolation BACKGROUND (grille ERA5 daily) -> grille FINE (centres MNT).
# Les centres des cellules du MNT, en points sf, servent de cibles bilineaires.
cells   <- terra::xyFromCell(Elev, seq_len(terra::ncell(Elev)))
pts_fin <- sf::st_as_sf(data.frame(x = cells[, 1], y = cells[, 2]),
                        coords = c("x", "y"), crs = sf::st_crs(Elev))
# la grille du background daily est lue depuis un fichier (lon/lat propres)
nc0  <- ncdf4::nc_open(bg_daily_pattern(YEARS[1]))
lon_bg <- ncdf4::ncvar_get(nc0, "longitude"); lat_bg <- ncdf4::ncvar_get(nc0, "latitude")
ncdf4::nc_close(nc0)
I_fine <- build_bilinear_to_stations(lon_bg, lat_bg, pts_fin)$I   # [ncell_fine, n_bg]
cat(sprintf("Matrice background->fin : %d x %d\n", nrow(I_fine), ncol(I_fine)))

# classes de production. REQUISES dès qu'un paramètre varie par classe (tau, fdry
# ou beta en mode class), OU si c_mode="class". La production applique à chaque
# jour le quadruplet de SA classe ; sans classes, produce_gridded_year s'arrete.
# En calage tout-global, elles sont facultatives (params identiques entre classes).
classe_prod_all <- NULL
if (!is.null(CLASSES_PROD_RDS)) {
  cp <- readRDS(CLASSES_PROD_RDS)                       # $date (Date), $Numero (1..K)
  classe_prod_all <- setNames(cp$Classe, as.character(cp$Date))
}

# --------------------------------------------------------- boucle annees
for (year in YEARS) {
  out_nc <- file.path(PROD_DIR, sprintf("pmod_2km_%s_%d.nc", AGREG, year))
  if (file.exists(out_nc)) { cat(sprintf("[%d] existe, saute.\n", year)); next }
  
  buf   <- make_coarse_buffer(path_common, path_sfc, opt = OPT)
  jours <- gen_jours(year, year)
  
  # background daily grossier de l'annee [maille_bg, jour], aligne sur jours
  f_bg <- bg_daily_pattern(year)
  if (!file.exists(f_bg)) { cat(sprintf("[%d] background absent, saute.\n", year)); next }
  bg_daily <- read_era5_daily_field(f_bg, jours$t_start)    # [n_bg, nD]
  
  # classe de chaque jour (si mode class)
  classe_prod <- NULL
  if (!is.null(classe_prod_all)) {
    dts <- format(jours$t_start, "%Y-%m-%d")
    classe_prod <- as.integer(classe_prod_all[dts])
  }
  
  cat(sprintf("[%d] production %d jours (agreg=%s) -> %s\n", year, nrow(jours), AGREG, basename(out_nc)))
  plan <- build_production_plan(opt, jours$t_start, classe_prod = classe_prod, c_mode = C_MODE)
  produce_gridded_year(plan, pre, buf, jours, mask_vec, A_agg, I_fine,
                       bg_daily = bg_daily, ortho_mask = ortho_mask,
                       out_nc = out_nc, agreg = AGREG,
                       clamp_c = CLAMP_C, spatial = SPATIAL, opt_variant = OPT,
                       progress = TRUE)
  rm(buf, bg_daily); gc()
}
cat("Production terminee.\n")

library(sf)
library(ggplot2)
library(dplyr)
library(tidyr)
library(lubridate)
library(terra)
library(tidyterra)



poly <- sf::st_read("C:/Users/cmartin/ownCloud/scriptsR/shape/DEP_GrandOuest.gpkg")
class(poly)

p2000 <- terra::rast("C:/Users/cmartin/ownCloud/partage_LPM_20260824/production_2km/prod_lcl_no-gamma/tau-cls_fdry-glob_beta-cls-fix5-0_c-daily_clamp-1.1.5_seed1/pmod_2km_daily_2000.nc")
p2000_era5 <- terra::rast("C:/Users/cmartin/ownCloud/daily/total_precip_daily_06utc_ouest_2000.nc")
crs(p2000) <- "EPSG:2154"

j <- 363

plot_mod <- ggplot() + 
  geom_spatraster(data=sum(p2000)) + 
  geom_sf(data=poly,fill=NA) +
  scale_fill_viridis_b(breaks=seq(0,2000,100))

plot_era5 <- ggplot() + 
  geom_spatraster(data=sum(p2000_era5)) + 
  geom_sf(data=poly,fill=NA) +
  scale_fill_viridis_b(breaks=seq(0,2000,100))

library(patchwork)
plot_era5 | plot_mod

#pour faire moyenne : sum(p2000) / nlayers(p2000)



library(terra)
library(sf)
library(ggplot2)
library(tidyterra)
library(patchwork)
library(viridis)

poly <- sf::st_read("C:/Users/cmartin/ownCloud/scriptsR/shape/DEP_GrandOuest.gpkg")
class(poly)

p2000 <- terra::rast("C:/Users/cmartin/ownCloud/partage_LPM_20260824/production_2km/prod_lcl_no-gamma/tau-cls_fdry-glob_beta-cls-fix5-0_c-daily_clamp-1.1.5_seed1/pmod_2km_daily_2000.nc")
p2000_era5 <- terra::rast("C:/Users/cmartin/ownCloud/daily/total_precip_daily_06utc_ouest_2000.nc")
crs(p2000) <- "EPSG:2154"

j <- 363

# On sélectionne la couche j (au lieu de faire sum())
plot_mod <- ggplot() + 
  geom_spatraster(data = p2000[[j]]) + 
  geom_sf(data = poly, fill = NA) +
  scale_fill_viridis_c(breaks = seq(0, 50, 5)) # Échelle adaptée à un jour (0-50mm)

plot_era5 <- ggplot() + 
  geom_spatraster(data = p2000_era5[[j]]) + 
  geom_sf(data = poly, fill = NA) +
  scale_fill_viridis_c(breaks = seq(0, 50, 5)) # Échelle adaptée à un jour

library(patchwork)
plot_era5 | plot_mod










library(terra)
library(sf)
library(ggplot2)
library(tidyterra)
library(patchwork)
library(viridis)

# ---------------------------------------------------------------------------
# 1. Chargement et Configuration
# ---------------------------------------------------------------------------

poly <- sf::st_read("C:/Users/cmartin/ownCloud/scriptsR/shape/DEP_GrandOuest.gpkg", quiet = TRUE)

p2000 <- terra::rast("C:/Users/cmartin/ownCloud/partage_LPM_20260824/production_2km/prod_lcl_no-gamma/tau-cls_fdry-glob_beta-cls-fix5-0_c-daily_clamp-1.1.5_seed1/pmod_2km_daily_2000.nc")
p2000_era5 <- terra::rast("C:/Users/cmartin/ownCloud/daily/total_precip_daily_06utc_ouest_2000.nc")

# Définition des CRS
crs(p2000) <- "EPSG:2154"
# ERA5 est souvent en WGS84, on le laissera tel quel, ggplot gérera la projection si on utilise coord_sf

# Jour cible
j <- 363
message(sprintf("Sélection du jour : %d", j))

# ---------------------------------------------------------------------------
# 2. Extraction du jour unique (SANS faire la somme)
# ---------------------------------------------------------------------------

# On extrait la couche j pour le modèle
r_mod_day <- p2000[[j]]

# On extrait la couche j pour ERA5
r_era5_day <- p2000_era5[[j]]

# ---------------------------------------------------------------------------
# 3. Application du Masque (Pour ne colorier que les départements)
# ---------------------------------------------------------------------------

# Conversion du polygone en objet terra
v_poly <- terra::vect(poly)

# Masque pour le Modèle (le polygone est projeté automatiquement dans le CRS du raster par mask)
r_mod_masked <- terra::mask(r_mod_day, v_poly)

# Masque pour ERA5 
# Note: Si ERA5 est en WGS84 et poly en WGS84, ça marche direct. 
# Si besoin de projeter le polygone vers ERA5 avant : v_poly_era5 <- terra::project(v_poly, crs(r_era5_day))
r_era5_masked <- terra::mask(r_era5_day, v_poly)

# On transforme aussi le polygone pour l'affichage (pour qu'il soit dans le même CRS que le raster affiché)
# Pour simplifier, on va forcer l'affichage en Lambert-93 (2154) pour les deux cartes
poly_l93 <- sf::st_transform(poly, 2154)

# ---------------------------------------------------------------------------
# 4. Création des graphiques (Jour unique)
# ---------------------------------------------------------------------------

# Échelle adaptée à un jour (0 à 50mm, ajustable si gros épisode)
# Si vous voulez voir les détails, gardez des breaks fins.
breaks_seq <- seq(0, 50, by = 5) 

plot_mod <- ggplot() + 
  geom_spatraster(data = r_mod_masked) + 
  geom_sf(data = poly_l93, fill = NA, colour = "white", linewidth = 0.3) +
  scale_fill_viridis_c(option = "D", name = "mm/jour", breaks = breaks_seq, limits = c(0, 50), na.value = "transparent") +
  labs(title = sprintf("Modèle LPM - Jour %d (2000)", j)) +
  coord_sf(crs = sf::st_crs(2154)) + # Force la projection Lambert-93
  theme_minimal() +
  theme(panel.background = element_rect(fill = "white"),
        plot.title = element_text(hjust = 0.5))

plot_era5 <- ggplot() + 
  geom_spatraster(data = r_era5_masked) + 
  geom_sf(data = poly_l93, fill = NA, colour = "white", linewidth = 0.3) +
  scale_fill_viridis_c(option = "D", name = "mm/jour", breaks = breaks_seq, limits = c(0, 50), na.value = "transparent") +
  labs(title = sprintf("ERA5 - Jour %d (2000)", j)) +
  coord_sf(crs = sf::st_crs(2154)) + # Force la projection Lambert-93
  theme_minimal() +
  theme(panel.background = element_rect(fill = "white"),
        plot.title = element_text(hjust = 0.5))

# Affichage côte à côte
plot_era5 | plot_mod

# Sauvegarde optionnelle
# ggsave(sprintf("comparaison_jour_%d_2000.png", j), plot_era5 | plot_mod, width = 14, height = 6, dpi = 300)
