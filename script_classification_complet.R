library(sf)
library(ggplot2)
library(dplyr)
library(tidyr)
library(lubridate)
library(terra)

#############################--------------------Chargement des données-------------------------------------------------

root_dir <- "C:/Users/cmartin/ownCloud/scriptsR/"

RR <- read.table( paste0(root_dir,"QUOT/RR_data.csv"),     #données des précipitations quotidiennes de MétéoFrance sur 2047 postes
                  header = TRUE, sep = ";")
RR[RR<0] <- NA

stations <- read.table( paste0(root_dir,"QUOT/metadonnees_stations_v2.csv"),   #données des postes (lat,lon,idMF)
                        header = TRUE, sep = ";")

stations_sf <- st_as_sf(stations, coords = c("lon","lat"), crs = 4326 ) # WGS84

DEP <- read_sf(paste0(root_dir, "shape/DEP_GrandOuest.gpkg"))   

Europe <- read_sf("C:/Users/cmartin/ownCloud/r/shapefile_Europe/Europe.shp")


###########################------------------Sélection des jours pluvieux----------------------------------------------

RR$Date <- lubridate::ymd(RR$Date)
RR_sub <- RR%>%
  mutate(annee=lubridate::year(Date))%>%
  filter(annee>=1993 & annee<=2008)

seuil_pluie <- RR_sub%>%  #calcul de la moyenne et application du seuil de 1mm
  dplyr::mutate(moyenne = rowMeans(across(-c(Date,annee)), na.rm = TRUE))%>%
  dplyr::filter(moyenne >1)

#sélection des stations avec -10% de NA

nb_na <- colSums(is.na(seuil_pluie)) #calcul du nombre de NA par postes
pourcentage_na <- nb_na/nrow(seuil_pluie)

stations10 <- names(pourcentage_na[pourcentage_na < 0.10])
postes <- seuil_pluie[,stations10]  #dataframe trié avec seulement les 393 postes pluvio ayant moins de 10% de nA

postes_na <- na.omit(postes) #suppression des jours contenant des NA sur les 393 postes sélectionnés

Dates_pluvieuses <- postes_na$Date  #correspond au dates à 0h
Dates_pluvieuses1 <- Dates_pluvieuses+1 #correspond aux dates à 24h


########################-------------------------Fusion fichiers raster-------------------------------------------------

Dates_int <- sort(unique(c(Dates_pluvieuses, Dates_pluvieuses+1)))

valid_time_JP <- (Dates_pluvieuses - as.Date("1970-01-01")) * 86400  #dates en seconde écoulées de puis le 01/01/1970
valid_time_J1 <- valid_time_JP +86400
valid_time_utile <- sort(unique(c(valid_time_JP, valid_time_J1)))


DIR_IN      <- "C:/Users/cmartin/donnees/era5/6_hourly/geopotential/crop/crop_net/joursP_bons1/" 
#DIR_IN = chemin vers les fichiers netcdf ERA5 avec les valeurs de géopotentiel à 500 et 850 hPa pour les dates à 0h et 24h
YEARS       <- 1993:2008


fusion_raster <- function(y) {       #un fichier netcdf par an, fusionne en un fichier 
  ncfile <- paste0(DIR_IN, "geopotential_500_850hPa_", y, "_crops.nc")
  message(ncfile)
  r <- terra::rast(ncfile)
  return(r)
}

r_tot <- lapply(YEARS, fusion_raster)

r_tot1 <- do.call(c, r_tot)  #spatraster avec valeurs géopotentiels des jours pluvieux et j+1

dates_utiles <- unique(time(r_tot1))

########################---------------------------Calculs gradients--------------------------------------------------

terra_gradient <- function(r)  
{
  dx <- terra::res(r)[1] 
  
  mx <- matrix(nrow = 3, ncol = 3,
               c(0, 0, 0,-1, 0, 1,0, 0, 0)/dx/2) 
  
  my <- matrix(nrow = 3, ncol = 3,
               c( 0, 1, 0,0, 0, 0,0,-1, 0)/dx/2) 
  
  grad_x <- terra::focal(r, w = mx, fun = "sum")
  grad_y <- terra::focal(r, w = my, fun = "sum")
  
  return(list(grad_x, grad_y))
  
}

grad <- terra_gradient(r_tot1)

ind_grad850 <- which(terra::depth(grad[[1]]) == 850) #séparation selon les valeurs de géopotentiel
ind_grad500 <- which(terra::depth(grad[[1]]) == 500)

gradx_850 <- grad[[1]][[ind_grad850]]  #spatraster des gradients ouest-est 
grady_850 <- grad[[2]][[ind_grad850]]  #spatraster des gradients nord-sud
grad_850 <- list(gradx_850,grady_850)

gradx_500 <- grad[[1]][[ind_grad500]]  #pareil à 500 hPa
grady_500 <- grad[[2]][[ind_grad500]]
grad_500 <- list(gradx_500,grady_500)

##sélection dates à 0h + la journée suivante pour que la matrice puisse être complète 
Dates_pluvieuses_alt <- c(Dates_pluvieuses,Dates_pluvieuses[1691]+1) 
ind_grad_alt <- which(time(grad_850[[1]]) %in% as.POSIXct(Dates_pluvieuses_alt))

gradx_850_JP_alt <- subset(grad_850[[1]], ind_grad_alt)
grady_850_JP_alt <- subset(grad_850[[2]], ind_grad_alt)

gradx_500_JP_alt <- subset(grad_500[[1]], ind_grad_alt) 
grady_500_JP_alt <- subset(grad_500[[2]], ind_grad_alt)

fenetre <- c(-823871.4, 1282339,5746220, 7895016) #réduction à la fenêtre spatiale utilisée
gradx_850_JP_alt <- terra::crop(gradx_850_JP_alt, fenetre)
grady_850_JP_alt <- terra::crop(grady_850_JP_alt, fenetre)
gradx_500_JP_alt <- terra::crop(gradx_500_JP_alt, fenetre)
grady_500_JP_alt <- terra::crop(grady_500_JP_alt, fenetre)

###################--------------------------Calculs scores de Teweles-Wobus----------------------------------------

indices_jpluvieux <- which(as.Date(terra::time(gradx_850_JP)) %in% as.character(Dates_pluvieuses))

precompute_norms <- function(gradx_list, grady_list) {  #calcul des normes pour le dénominateur en amont afin de limiter le temps de calcul
  n <- terra::nlyr(gradx_list)
  sx <- numeric(n)
  sy <- numeric(n)
  for (i in seq_len(n)) {
    sx[i] <- terra::global(abs(gradx_list[[i]]), "sum", na.rm = TRUE)[1,1]
    sy[i] <- terra::global(abs(grady_list[[i]]), "sum", na.rm = TRUE)[1,1]
  }
  list(sx = sx, sy = sy)
}

norms_850 <- precompute_norms(gradx_850_JP_alt, grady_850_JP_alt)
norms_500 <- precompute_norms(gradx_500_JP_alt, grady_500_JP_alt)

# Retourne une matrice n×n des score de Teweles-Wobus pour un niveau de pression donné

tw_matrix_level <- function(gradx_list, grady_list, norms, indices_jpluvieux, offset = 0L){
  
  n <- length(indices_jpluvieux)
  idx <- indices_jpluvieux + offset
  sx  <- norms$sx[idx]
  sy  <- norms$sy[idx]
  
  mat_num_x <- matrix(0, n, n)
  mat_num_y <- matrix(0, n, n)
  
  for (j in seq_len(n - 1)) {   
    
    k_idx <- (j + 1):n   
    
    # empilement de (gradx[k] - gradx[j]) pour tous les k en une seule soustraction
    diff_x <- abs(gradx_list[[idx[k_idx]]] - gradx_list[[idx[j]]])  
    diff_y <- abs(grady_list[[idx[k_idx]]] - grady_list[[idx[j]]])
    
    sums_x <- terra::global(diff_x, "sum", na.rm = TRUE)[, 1] #calcule séparément termes numérateur
    sums_y <- terra::global(diff_y, "sum", na.rm = TRUE)[, 1] 
    
    mat_num_x[j, k_idx] <- sums_x
    mat_num_y[j, k_idx] <- sums_y
  }
  
  # symétrisation
  mat_num_x <- mat_num_x + t(mat_num_x)
  mat_num_y <- mat_num_y + t(mat_num_y)
  
  # dénominateur scalaire
  sx_mat  <- outer(sx, sx, `+`)
  sy_mat  <- outer(sy, sy, `+`)
  mat_den <- (sx_mat + mat_num_x) / 2 + (sy_mat + mat_num_y) / 2
  
  100 * (mat_num_x + mat_num_y) / mat_den
}

#calcul final
compute_tw_matrix <- function(DatesP) {
  
  tw1 <- tw_matrix_level(gradx_850_JP_alt, grady_850_JP_alt, norms_850, indices_jpluvieux, offset = 0L)
  tw2 <- tw_matrix_level(gradx_850_JP_alt, grady_850_JP_alt, norms_850, indices_jpluvieux, offset = 1L)
  tw3 <- tw_matrix_level(gradx_500_JP_alt, grady_500_JP_alt, norms_500, indices_jpluvieux, offset = 0L)
  tw4 <- tw_matrix_level(gradx_500_JP_alt, grady_500_JP_alt, norms_500, indices_jpluvieux, offset = 1L)
  
  dist_mat <- tw1 + tw2 + tw3 + tw4
  
  # NA sur la diagonale et le triangle inférieur pour plus de rapidité
  dist_mat[lower.tri(dist_mat, diag = TRUE)] <- NA
  
  rownames(dist_mat) <- colnames(dist_mat) <- as.character(DatesP)
  dist_mat
}

library(tictoc)
tic()
dist_matrix <- compute_tw_matrix(as.character(Dates_pluvieuses)) 
toc()

##############--------------------------------------Classification-------------------------------------------------

diag(dist_matrix) <- 0 #diagonale à zéro

#copier triangle sup pour remplir triangle inf et enlever NA
dist_matrix[lower.tri(dist_matrix)] <- t(dist_matrix)[lower.tri(dist_matrix)]
all.equal(dist_matrix, t(dist_matrix)) # vérifier que matrice symétrique

dist_mat <- as.dist(dist_matrix)  #matrice de distance finale

library(cluster)
classif <- hclust(dist_mat, method="ward.D2") #CAH en utilisant l'algorithme de Ward
dendro <- plot(classif)    #obtention du dendrogramme
abline(h = 1900, col = "red")

classes_syn <- cutree(classif,k=4)  #découpage en quatre classes

##############-------------------------Calcul des centroïdes----------------------------------------

classes <- as.data.frame(classes_syn)
colnames(classes) <- c("Numero")
Classe <- classes%>%
  mutate(Dates = rownames(classes), .before = Numero)

library(tidyterra) 

centroidej1 <- list() #icic alcul des centro¨des à 24h mais même script à 0h

for(wp in seq(1,max(Classe$Numero))){      
  DatesClasse <- as.Date(Classe$Dates[Classe$Numero==wp]) 
  d <- DatesClasse+1
  Annee <- year(d)                #prendre d ou Datesclasse pour calculer respectivement à 24h ou 0h
  
  sum_r <- NULL
  n_r <- 0 
  
  for(a in unique(Annee))
  {
    DatesClasseAn <- d[which(Annee == a)]
    valid_time_cible <- (DatesClasseAn - as.Date("1970-01-01"))*86400
    
    ncfile_z <- sprintf("C:/Users/cmartin/donnees/era5/6_hourly/geopotential/crop/crop_net/joursP_bons1/geopotential_500_850hPa_%4d_crops.nc",a)
    nc_z <- ncdf4::nc_open( ncfile_z )
    
    valid_time <- ncdf4::ncvar_get(nc_z,"valid_time")
    plevels <- ncdf4::ncvar_get(nc_z,"pressure_level")
    ix_t <- which(valid_time %in% valid_time_cible)
    
    varsize <- nc_z$var[["z"]]$varsize
    ndims <- nc_z$var[["z"]]$ndims
    nlon <- varsize[1]
    nlat <- varsize[2]
    nplev <- varsize[3]
    ntime <- varsize[4]
    
    for(i in ix_t) 
    {
      # récupération valeurs géopotentiel
      
      start <- c(1,1,1,i)              
      count <- c(nlon,nlat,nplev,1)    #remplacer nplev par 1 pour n'avoir qu'un niveau de pression (graphiques)
                                      #et start <- c(1,1,1,i) pour 850 hPa ou c(1,1,2,i) pour 500 hPa
      r <- ncdf4::ncvar_get( nc_z, "z", start = start, count = count )
      
      if (is.null(sum_r)) {
        sum_r <- array(0,dim = dim(r))}
      
      sum_r <- sum_r + r
      n_r <- n_r + 1
      
    }     
    ncdf4::nc_close(nc_z)
    #ncdf4::nc_close(nc_t)
  }
  
  moyenne <- sum_r / n_r
  centroidej1[[as.character(wp)]] <- moyenne
}

######---------------------------------Classification des journées de 1960 à 2025--------------------------

#Les journées ont été récupérées commme précédemment en utilisant la fonction fusion_raster
#Les gradients ont également été calculés avec la fonction terra_gradient et les normes avec precompute_norms

#Le même procédé a été réalisé pour les centroïdes des classes, le script étant très similaire, il n'est pas présenté ici

#fonction pour calculer score de Teweles-Wobus entre un jour et un centroïde
tw_score_pair <- function(gradx_a, grady_a, gradx_b, grady_b,  
                          sx_a, sy_a, sx_b, sy_b) {
  
  diff_x <- abs(gradx_a - gradx_b)
  diff_y <- abs(grady_a - grady_b)
  
  num_x <- terra::global(diff_x, "sum", na.rm = TRUE)[1, 1]
  num_y <- terra::global(diff_y, "sum", na.rm = TRUE)[1, 1]
  
  den <- (sx_a + sx_b + num_x) / 2 + (sy_a + sy_b + num_y) / 2
  
  100 * (num_x + num_y) / den
}

#obtention de la matrice de distance joursxcentroïdes
library(tictoc)
tic()

n_days    <- length(ind)
classes   <- names(centroide_gradj)
dist_mat_days_classes <- matrix(NA, nrow = n_days, ncol = length(classes),
                                dimnames = list(as.character(dates_ind), classes))

for (d in seq_len(n_days)) {
  
  iJ  <- ind[d]
  iJ1 <- ind1[d]
  
  gx850_J  <- gradx_850_j[[iJ]];  gy850_J  <- grady_850_j[[iJ]]
  gx850_J1 <- gradx_850_j[[iJ1]]; gy850_J1 <- grady_850_j[[iJ1]]
  gx500_J  <- gradx_500_j[[iJ]];  gy500_J  <- grady_500_j[[iJ]]
  gx500_J1 <- gradx_500_j[[iJ1]]; gy500_J1 <- grady_500_j[[iJ1]]
  
  sx850_J  <- norms_850$sx[iJ];  sy850_J  <- norms_850$sy[iJ]
  sx850_J1 <- norms_850$sx[iJ1]; sy850_J1 <- norms_850$sy[iJ1]
  sx500_J  <- norms_500$sx[iJ];  sy500_J  <- norms_500$sy[iJ]
  sx500_J1 <- norms_500$sx[iJ1]; sy500_J1 <- norms_500$sy[iJ1]
  
  for (wp in classes) {
    
    cj  <- centroide_gradj[[wp]]
    cj1  <- centroide_gradj1[[wp]]
    nj <- norms_centroidej[[wp]]
    nj1 <- norms_centroidej1[[wp]]
    
    tw_850_J  <- tw_score_pair(gx850_J,  gy850_J,  cj$gradx_850,  cj$grady_850, 
                               sx850_J,  sy850_J,  nj$sx850, nj$sy850)           
    tw_850_J1 <- tw_score_pair(gx850_J1, gy850_J1, cj1$gradx_850, cj1$grady_850,
                               sx850_J1, sy850_J1, nj1$sx850, nj1$sy850)
    tw_500_J  <- tw_score_pair(gx500_J,  gy500_J,  cj$gradx_500,  cj$grady_500,
                               sx500_J,  sy500_J,  nj$sx500, nj$sy500)
    tw_500_J1 <- tw_score_pair(gx500_J1, gy500_J1, cj1$gradx_500, cj1$grady_500,
                               sx500_J1, sy500_J1, nj1$sx500, nj1$sy500)
    
    dist_mat_days_classes[d, wp] <- tw_850_J + tw_850_J1 + tw_500_J + tw_500_J1
  }
}
toc()

#répartition des jours dans la bonne classe en choisissant la classe pour laquelle le score est minimal
classes_assignees <- apply(dist_mat_days_classes, 1, function(x) classes[which.min(x)])

resultat <- data.frame(
  Date   = dates_ind,
  Classe = classes_assignees)



















