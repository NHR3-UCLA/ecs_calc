#This file contains all the R functions that are used to compute the event 
#coordinate system (ECS)

#libraries
#spatial data libraries
library(sp)
#parallel processing libraries
library(parallel)
#spline libraries
library(mgcv)

#function to compute u&t coordinates from lat&lon coordinates
latlon2reftrace <- function(fault_disp, ref_trace){
  
  #compute utm coorinates and add at the end of dataframe
  out <- longlat2UTM(long = ref_trace$Longitude, lat = ref_trace$Latitude)
  ref_trace_utm <- out[[1]]
  utm_zone <- out[[2]]
  rm(out)
  
  
  #center x&y UTM coordinates of ECS
  utm_offset_xy = c(mean(ref_trace_utm$x), mean(ref_trace_utm$y))
  #offset Cartesian coordinates to zero mean
  ref_trace_utm$x <- ref_trace_utm$x - utm_offset_xy[1]
  ref_trace_utm$y <- ref_trace_utm$y - utm_offset_xy[2]
  #compute reference trace u & t coordinates
  ref_trace_ut <- gc2ext_ut(ref_trace_utm, ref_trace_utm)
  ref_trace_utm$u <- ref_trace_ut[,1]
  ref_trace_utm$t <- ref_trace_ut[,2]
    
  #compute utm coordinates and add at the end of dataframe
  out <- longlat2UTM(long = fault_disp$Longitude, lat = fault_disp$Latitude, utm_zone=utm_zone)
  fault_disp_utm <- out[[1]]
  rm(out)
  
  #center x&y UTM coordinates of fault displacements
  fault_disp_utm$x <- fault_disp_utm$x - utm_offset_xy[1]
  fault_disp_utm$y <- fault_disp_utm$y - utm_offset_xy[2]
  
  #compute reference trace for disp. points
  disp_pt_ut <- gc2ext_ut(fault_disp_utm, ref_trace_utm)
  fault_gc2 <- data.frame(u = disp_pt_ut[,1], t = disp_pt_ut[,2]) 
  
  #compute curvature of reference trace
  ref_trace_utm$curv <- compute_curvature(ref_trace_utm$u, ref_trace_utm$x, ref_trace_utm$y)
  #compute curvature at fault disp locations
  fault_gc2$curv     <- approx(ref_trace_utm$u, y = ref_trace_utm$curv, fault_gc2$u, method="linear", rule = 2)$y
  
  return(fault_gc2)
}

#function to determine the UTM zone
longlat2UTMzone <- function(long, lat) {
  utm_zone <- (floor((long + 180)/6) %% 60) + 1

  return(as.numeric(utm_zone))  
}

#function to convert longitude and latitude coordinates to UTM
longlat2UTM <- function(long, lat, utm_zone = NA){
  
  #create long lat dataframe
  df_longlat <- data.frame(Longitude = long ,Latitude = lat)
  
  #mean longitude and latitude
  mean_longlat <- colMeans(df_longlat[c("Longitude","Latitude")])  
  
  #utm zone
  if (is.na(utm_zone)){
    utm_zone <- longlat2UTMzone(long = mean_longlat[1], lat = mean_longlat[2])    
  }
  
  #define current datum 
  coordinates(df_longlat) <- c("Longitude","Latitude")
  proj4string(df_longlat) <- CRS("+proj=longlat +datum=WGS84")
  
  #transform to utm coordinates
  df_utmxy <- spTransform(df_longlat, CRS(paste("+proj=utm +zone=",utm_zone," +datum=WGS84",sep='')))
  df_utmxy <- as.data.frame(df_utmxy)
  colnames(df_utmxy) <- c('x', 'y')
  
  return(list(df_utmxy, utm_zone))
}

#function to convert UTM coordinates to longitude and latitude
UTM2longlat <- function(x_utm, y_utm, utm_zone){
  
  #create long lat dataframe
  df_utm_xy<- data.frame(x_utm = x_utm ,y_utm = y_utm)
  

  #define current datum 
  coordinates(df_utm_xy) <- c("x_utm","y_utm")
  proj4string(df_utm_xy) <- CRS(paste("+proj=utm +zone=",utm_zone," +datum=WGS84",sep=''))

  #transform to utm coordinates
  df_longlat <- spTransform(df_utm_xy, CRS("+proj=longlat +datum=WGS84"))
  df_longlat <- as.data.frame(df_longlat)
  colnames(df_longlat) <- c("Longitude","Latitude")
  
  return(df_longlat)
}

#R wrapper function to compute the GC2 coordinates
gc2_ut <- function(fault_disp, ecs_df){
  #number of disp points in each batch
  n_batch <- 50000
  
  #check if distance metrics FORTRAN library is loaded
  if (!is.loaded('r_wrapper_gc2')) {
    #set path FORTRAN library
    #file path for FORTRAN subroutine
    filepath_lib <- file.path('R_lib','R_dist_metrics.so')
    dyn.load(filepath_lib)
  }
  
  #coordinates disp points
  disp_pt_xy <- cbind(fault_disp$x,fault_disp$y)
  #number of displacement points
  n_disp_pt <- nrow(disp_pt_xy)
  
  #ECS coordinates
  flt_xy <- cbind(ecs_df$x, ecs_df$y)
  #number of ECS points
  n_pt_flt <- nrow(flt_xy)
  
  #compute distance metrics
  if (n_disp_pt <= n_batch){
    disp_pt_ut <- .Fortran('r_wrapper_gc2', flt_xy = flt_xy, sta_xy = disp_pt_xy, 
                           n_pt_flt = n_pt_flt, n_sta = n_disp_pt, sta_ut = disp_pt_xy)$sta_ut    
  }else{
    #number of cores for parallelization
    n_cores <- max(detectCores()-1,1)
    #indices for frist and last batch point
    i_b_s <- seq(1,n_disp_pt-1,n_batch)
    i_b_e <- c(i_b_s[2:length(i_b_s)]-1, n_disp_pt)
    #batches of displacement points 
    disp_pt_b_xy <- lapply( seq(length(i_b_s)), function(k) disp_pt_xy[i_b_s[k]:i_b_e[k], ] )
    n_disp_pt_b <-  lapply( disp_pt_b_xy, function(data) nrow(data) )
    #compute distance metrics of each disp batch
    disp_pt_ut <-   mcMap( function(d_pt_xy, n_d_pt) .Fortran('r_wrapper_gc2', flt_xy = flt_xy, sta_xy = d_pt_xy, 
                                                            n_pt_flt = n_pt_flt, n_sta = n_d_pt, sta_ut = d_pt_xy)$sta_ut,
                           disp_pt_b_xy, n_disp_pt_b, mc.cores = getOption("mc.cores", n_cores) )
    disp_pt_ut <- do.call(rbind, disp_pt_ut)
  }

  return(disp_pt_ut)
}

#GC2 projection extended
gc2ext_ut <-  function(fault_disp, ecs_df, ext_len = 50000){

  #ecs unit vectors
  ecs_n_df <- data.frame(x=diff(ecs_df$x), y=diff(ecs_df$y))
  ecs_n_df <- ecs_n_df/sqrt(ecs_n_df$x^2 + ecs_n_df$y^2)

  #extended starting point
  ecs_extstr_df <- ecs_df[1,c('x','y')]            - ext_len * ecs_n_df[1,] 
  #extended end point
  ecs_extend_df <- ecs_df[nrow(ecs_df),c('x','y')] + ext_len * ecs_n_df[nrow(ecs_n_df),] 
  #extended ecs
  ecs_df_ext <- rbind.fill(ecs_extstr_df, ecs_df, ecs_extend_df)
  rownames(ecs_df_ext) <- NULL
  
  #compute u&t coordinates with original gc2 system
  disp_pt_ut <- gc2_ut(fault_disp, ecs_df_ext)
  #correct for extended points
  ecs_str_ut <-  gc2_ut(ecs_df[1,], ecs_df_ext) #u&t coordinates of first original ecs point
  assert_that(abs(ecs_str_ut[1,1] - ext_len)<3e-2)
  disp_pt_ut[,1] <- disp_pt_ut[,1] - ecs_str_ut[1,1] 
  
  return(disp_pt_ut)
}

#Reference trace update
update_ecs <- function(fault_disp, nom_trace_orig, lamb_spline = 1, ecs_du = 250, wt = NULL, alpha = 1 ){
  
  #compute GC2 coordinates of displacement points based on the original ECS
  disp_pt_ut <- gc2ext_ut(fault_disp, nom_trace_orig)
  #store the u and t coordiantes in the fault_disp data frame
  fault_disp$u <- disp_pt_ut[,1] #along strike distance of displacement points
  fault_disp$t <- disp_pt_ut[,2] #perpendicular distance of displacement points
  
  #fault lenght 
  fault_len <- max(disp_pt_ut[,1]) - min(disp_pt_ut[,1])

  #fit the ECS spline on the x and y coordinates as a function of the along strike dist.
  #along strike points to evaluate the ECS
  # ecs_u <- seq(floor(min(fault_disp$u)/ecs_du)*ecs_du, ceiling(max(fault_disp$u)/ecs_du)*ecs_du, ecs_du)
  ecs_u <- seq(round(min(fault_disp$u)/ecs_du)*ecs_du, round(max(fault_disp$u)/ecs_du)*ecs_du, ecs_du)
  
  # #smooth spline function
  # #estimate coefficients
  # fit_ss_x <- smooth.spline(fault_disp$u, fault_disp$x , lambda = lamb_spline/fault_len, w = weights)
  # fit_ss_y <- smooth.spline(fault_disp$u, fault_disp$y , lambda = lamb_spline/fault_len, w = weights)
  # #compute updated ECS
  # ecs_upd = data.frame(x = predict(fit_ss_x,ecs_u)$y, y = predict(fit_ss_y,ecs_u)$y)
  
  #GAM estimate x & y simultaniously
  #estimate coefficients
  # fit_gam_xy <- gam(list(x ~ s(u, k = 5) -1, y ~  s(u, k = 5) -1), data = fault_disp, family = mvn(d = 2), sp = c(1,1)*lamb_spline/fault_len )
  fit_gam_xy <- mgcv::gam(list(x ~ s(u) -1, y ~  s(u) -1), data = fault_disp, family = mvn(d = 2), sp = c(1,1)*lamb_spline/fault_len )
  # fit_gam_xy <- gam(list(x ~ s(u) -1, y ~  s(u) -1), data = fault_disp, family = mvn(d = 2), sp = c(1,1)*lamb_spline/fault_len, weights = wt )
  #compute updated ECS
  ecs_upd = data.frame(predict(fit_gam_xy,data.frame(u=ecs_u)))
  colnames(ecs_upd) <- c("x","y")
  nom_trace_ut <- gc2ext_ut(ecs_upd, ecs_upd)
  ecs_upd$u <- nom_trace_ut[,1]
  ecs_upd$t <- nom_trace_ut[,2]
  
  #update GC2 coordinates of displacement points for new ECS
  disp_pt_ut <- gc2ext_ut(fault_disp, ecs_upd)
  fault_disp$u <- disp_pt_ut[,1] #along strike distance of displacement points
  fault_disp$t <- disp_pt_ut[,2] #perpendicular distance of displacement points
  return(list(ecs_upd, fault_disp, fit_gam_xy))
}


#Reference trace update
update_ecs_wolambda <- function(fault_disp, nom_trace_orig, ecs_du = 250, wt = NULL, alpha = 1 ){
  
  #compute GC2 coordinates of displacement points based on the original ECS
  disp_pt_ut <- gc2ext_ut(fault_disp, nom_trace_orig)
  #store the u and t coordiantes in the fault_disp data frame
  fault_disp$u <- disp_pt_ut[,1] #along strike distance of displacement points
  fault_disp$t <- disp_pt_ut[,2] #perpendicular distance of displacement points
  
  #fault lenght 
  fault_len = max(disp_pt_ut[,1]) - min(disp_pt_ut[,2])
  
  #fit the ECS spline on the x and y coordinates as a function of the along strike dist.
  #along strike points to evaluate the ECS
  ecs_u <- seq(floor(min(fault_disp$u)/ecs_du)*ecs_du, ceiling(max(fault_disp$u)/ecs_du)*ecs_du, ecs_du)
  
  # #smooth spline function
  # #estimate coefficients
  # fit_ss_x <- smooth.spline(fault_disp$u, fault_disp$x , lambda = lamb_spline/fault_len, w = weights)
  # fit_ss_y <- smooth.spline(fault_disp$u, fault_disp$y , lambda = lamb_spline/fault_len, w = weights)
  # #compute updated ECS
  # ecs_upd = data.frame(x = predict(fit_ss_x,ecs_u)$y, y = predict(fit_ss_y,ecs_u)$y)
  
  #GAM estimate x & y simultaniously
  #estimate coefficients
  # fit_gam_xy <- gam(list(x ~ s(u, k = 5) -1, y ~  s(u, k = 5) -1), data = fault_disp, family = mvn(d = 2) )
  fit_gam_xy <- mgcv::gam(list(x ~ s(u) -1, y ~  s(u) -1), data = fault_disp, family = mvn(d = 2)  )
  # fit_gam_xy <- gam(list(x ~ s(u) -1, y ~  s(u) -1), data = fault_disp, family = mvn(d = 2) )
  #compute updated ECS
  ecs_upd = data.frame(predict(fit_gam_xy,data.frame(u=ecs_u)))
  colnames(ecs_upd) <- c("x","y")
  nom_trace_ut <- gc2ext_ut(ecs_upd, ecs_upd)
  ecs_upd$u <- nom_trace_ut[,1]
  ecs_upd$t <- nom_trace_ut[,2]
  
  
  #update GC2 coordinates of displacement points for new ECS
  disp_pt_ut <- gc2ext_ut(fault_disp, ecs_upd)
  fault_disp$u <- disp_pt_ut[,1] #along strike distance of displacement points
  fault_disp$t <- disp_pt_ut[,2] #perpendicular distance of displacement points
  
  return(list(ecs_upd, fault_disp, fit_gam_xy))
}

#Main function of definition of reference trace
ecs_main <- function(fault_disp, lamb_spline = NaN, ecs_du = 100, flt_max_ds = 50, calc_xy = FALSE, start_sol = 'PCA', ...){
  #ecs_main computes the ECS based on the location and magnitude of displacement measurements
  #Input Arguments:
  #   fault_disp [data-frame]: contains information about the location and magnitude of displacements
  #   lamb_spline [real number]: penalty coefficient controling the smoothness of the ECS
  #   ecs_du [real number]: along strike interval of ECS
  #   flt_max_ds [real number]: average thresshold offset for ECS location between iterations
  #   calc_xy [logical var.]: if true compute caritcian coordinates based on latitude and longitude 
  #   start_sol [string]: option for coordinate system starting solution
  #                         'PCA': pricipal component
  #                         'MRS': mean rupture strike
  
  #parse extra arguments
  xargs <- list(...)
  
  if (calc_xy){
    #compute utm coordinates and add at the end of dataframe
    out <- longlat2UTM(long = fault_disp$Longitude, lat = fault_disp$Latitude)
    df_utmxy <- out[[1]]
    utm_zone <- out[[2]]
    fault_disp <- cbind(fault_disp, df_utmxy)
    rm(df_utmxy, out)
  }

  #center x and y UTM coordinates
  utm_offset_xy = c(mean(fault_disp$x), mean(fault_disp$y))
  #offset Cartesian coordinates to zero mean
  fault_disp$x <- fault_disp$x - utm_offset_xy[1]
  fault_disp$y <- fault_disp$y - utm_offset_xy[2]
  
    
  #compute the starting ECS
  if (start_sol == 'PCA'){ #principal component
    #total least squares coefficients
    xy_pca <- fault_disp[c('x','y')] * fault_disp$wt
    xy_pca <- t(t(xy_pca) - colMeans(xy_pca))
    #compute principal directions
    v <- prcomp(xy_pca)$rotation
    rot_th <- -1 * atan2(v[2,1], v[1,1])  #PCA angle
    rm(v, xy_pca)     
  }else if(start_sol == 'MRS'){ #mean rupture strike
    #browser()
    rup_data <- xargs$rup_data
    if (calc_xy){
      #compute utm coordinates and add at the end of dataframe
      out <- longlat2UTM(long = rup_data$Longitude, lat = rup_data$Latitude, utm_zone = utm_zone)
      df_utmxy <- out[[1]]
      rup_data <- cbind(rup_data, df_utmxy)
      rm(df_utmxy, out)
    }
    #offset Cartesian coordinates based on utm offset of disp data
    rup_data$x <- rup_data$x - utm_offset_xy[1]
    rup_data$y <- rup_data$y - utm_offset_xy[2]
    
    #unique ruptures
    rup_ids <- unique(rup_data$RUP_ID)
    
    #compute average fualt_rupture strike
    rup_strike_ang <- c()
    rup_len <- c()
    for (r_id in rup_ids){
      rup_seg <- subset(rup_data, RUP_ID == r_id)
      #compute rupture length and strike
      rup_strike_ang <- c( rup_strike_ang, rup_avg_strike(rup_seg) )
      rup_len <- c( rup_len, rup_length(rup_seg) )
    }
    
    #weighted average fault strike angle
    fault_stike_ang <- weighted.mean(rup_strike_ang, rup_len^2)
    #rotation angle
    rot_th <- -1 * fault_stike_ang
    rm(rup_data,rup_ids, rup_strike_ang, rup_len, rup_seg, fault_stike_ang)
  }else{
    stop('Unavailable starting solution')
  }
  
  #number of points, stating solution
  n_pt = 10
  #define rotation matrix
  rot_mat <- matrix(data = c(cos(rot_th), sin(rot_th), -sin(rot_th), cos(rot_th)), nrow = 2, ncol = 2)
  #rotate displacement points to find the extents of the nominal trace
  disp_pt_ut <- t( rot_mat %*% t( as.matrix(fault_disp[c('x','y')]) ) ) 
  #ECS in rotated space
  ecs_ut <- cbind(seq(min(disp_pt_ut[,1]), max(disp_pt_ut[,1]), length.out = n_pt), rep(0, n_pt))
  #rotate ECS to original axes
  ecs_xy <- t( solve(rot_mat) %*% t(ecs_ut) )
  #ECS, starting solution
  ecs_df = data.frame(x = ecs_xy[,1], y = ecs_xy[,2], u = ecs_ut[,1], t = ecs_ut[,2]) 
  ecs_ut <- gc2ext_ut(ecs_df, ecs_df)
  ecs_df$u <- ecs_ut[,1]
  ecs_df$t <- ecs_ut[,2]
  
  #save all ecs solutions 
  ecs_df2save <- list()
  ecs_df2save[[1]] <- ecs_df
  
  #compute GC2 coordinates of displacement points with respect to the ECS, calculate it again with the GC2 strike definition
  disp_pt_ut <- gc2ext_ut(fault_disp, ecs_df)
  fault_disp$u <- disp_pt_ut[,1] #along strike distance of displacement points
  fault_disp$t <- disp_pt_ut[,2] #perpendicular distance of displacement points
  rm(n_pt, rot_th, rot_mat, ecs_xy, disp_pt_ut, ecs_ut)  
  
  #update ECS
  wt  =  fault_disp$wt / sum(fault_disp$wt)
  # weights = NULL
  
  #iterate until ECS is stable
  flag_stop = FALSE

  k <- 1
  while (!flag_stop) {
    k <- k + 1
    #update ECS
    if (!is.na(lamb_spline)){
      out <- update_ecs(fault_disp, ecs_df, lamb_spline = lamb_spline, ecs_du = ecs_du, wt = wt)
    } else {
      #out <- update_ecs(fault_disp, ecs_df, lamb_spline = lamb_spline, ecs_du = ecs_du, wt = wt)
      out <- update_ecs_wolambda(fault_disp, ecs_df, ecs_du = ecs_du, wt = wt)
    }

    ecs_df <- out[[1]]
    ecs_df2save[[k]] <- ecs_df
    fault_disp_new <- out[[2]]
    fit_nom_trace <- out[[3]]

    #average offset of ECS with respect to disp points
    flt_ds <- mean(sqrt((fault_disp_new$u - fault_disp$u)^2 + (fault_disp_new$t - fault_disp$t)^2))
    if( flt_ds < flt_max_ds){
      flag_stop = TRUE
    }
    
    fault_disp <- fault_disp_new
  }
  
  #convert back to utm coordinates 
  ecs_df$x <- ecs_df$x + utm_offset_xy[1]
  ecs_df$y <- ecs_df$y + utm_offset_xy[2]
  for ( k in seq(length(ecs_df2save)) ){
    ecs_df2save[[k]]$x <- ecs_df2save[[k]]$x + utm_offset_xy[1]
    ecs_df2save[[k]]$y <- ecs_df2save[[k]]$y + utm_offset_xy[2]
  }
  fault_disp$x <- fault_disp$x + utm_offset_xy[1]
  fault_disp$y <- fault_disp$y + utm_offset_xy[2]
  
  #add reference line node ID
  ecs_df$REF_ID <- seq(nrow(ecs_df))
  
  #compute latitude and longitude coordinates of the ECS
  if (calc_xy){
    nom_trace_latlon <- UTM2longlat(x_utm = ecs_df$x, y_utm = ecs_df$y, utm_zone = utm_zone)
  
    #ECS latitude and longitude (final verison)
    ecs_df$Longitude <- nom_trace_latlon[,1] 
    ecs_df$Latitude <- nom_trace_latlon[,2]
    
    for ( k in seq(length(ecs_df2save)) ){
      nom_trace_latlon <- UTM2longlat(x_utm = ecs_df2save[[k]]$x, y_utm = ecs_df2save[[k]]$y, utm_zone = utm_zone)
      ecs_df2save[[k]]$Longitude <- nom_trace_latlon[,1] 
      ecs_df2save[[k]]$Latitude  <- nom_trace_latlon[,2]
      ecs_df2save[[k]]$REF_ID    <- seq(nrow(ecs_df2save[[k]]))
    }
  }
  
  #function to return the ECS coordinates as a function of along strike distance
  fun_nom_trace <- function(u_array){
    
    #compute cartesian coordinates of ECS
    ecs_df = data.frame(predict(fit_nom_trace, data.frame(u = u_array)))
    colnames(ecs_df) <- c("x","y")
    #convert cartesian to utm coordinates 
    ecs_df$x <- ecs_df$x +  utm_offset_xy[1] 
    ecs_df$y <- ecs_df$y +  utm_offset_xy[2]
    
    #compute latitude and longitude coordinates
    if (calc_xy){
      nom_trace_latlon <- UTM2longlat(x_utm = ecs_df$x, y_utm = ecs_df$y, utm_zone = utm_zone)
      ecs_df$Longitude <- nom_trace_latlon[,1] #ECS latitude and longitude
      ecs_df$Latitude <- nom_trace_latlon[,2]
    }

    #along strike and perpendicular distance
    ecs_df$u <- u_array
    ecs_df$t <- 0
    
    #add reference line node ID
    ecs_df$REF_ID <- seq(nrow(ecs_df))
    return(ecs_df)
  }
  
  #compute curvature
  ecs_df$curv <- compute_curvature(ecs_df$u, ecs_df$x, ecs_df$y)

  #return fault displacement and ECS data-frame
  return(list(ecs_df = ecs_df, fault_disp = fault_disp, fit_nom_trace = fit_nom_trace, fun_nom_trace = fun_nom_trace, ecs_iter = ecs_df2save))
}


## Auxiliary Functions
##-----------------------

weight_ecs_data <- function(data_ecs, ratio_thres = 100){
  #weight_ecs_data returns a data-frame in which each observation is repeated
  #based on each weight  

  #number of cores for parallelization
  n_cores <- max(detectCores()-1,1)

  #number of points
  n_pt <- nrow(data_ecs)
  
  #weights array
  wt_array <- data_ecs$wt
  #normalize weights
  wt_array <- wt_array / min(wt_array)
  #minimum, maximum and median normalized weights
  min_wt <- min(wt_array)
  max_wt <- max(wt_array)
  med_wt <- median(wt_array)
  
  #ratio maximum to minimum weight
  ratio_min2max <- max_wt/min_wt
  
  #compute reference weight
  #check if max2min ratio is greater than threshold and cap it
  if(ratio_min2max > ratio_thres){
    print("Exceeding weight threshold ratio")
    wt_array <- pmin(wt_array, ratio_thres)
  }
  
  #repeat every observation based on the weight ratio
  n_rep_array <- pmin( pmax(round(wt_array), 1 ), ratio_thres)
  data_ecs_wt <- mcMap(function(k, n_rep) data_ecs[rep(k, each=n_rep), ] , seq(n_pt), n_rep_array,
                       mc.cores = getOption("mc.cores", n_cores) )
  data_ecs_wt <- do.call(rbind, data_ecs_wt)
  rownames(data_ecs_wt) <- NULL
  #reset weights
  data_ecs_wt$wt <- 1
  
  return(data_ecs_wt)  
}

rup_avg_strike <- function(rup_data){
  
  #ensure that all rupture points belong to a single rupture
  assert_that( length(unique(rup_data$RUP_ID)) == 1 )
  
  #sort all rupture points based on the node id
  rup_data <- rup_data[order(rup_data$NODE_ID),]

  #compute rupture segment dx and dy intervals
  seg_dx <- diff(rup_data$x)
  seg_dy <- diff(rup_data$y)
  #compute segment lengths
  seg_len <- sqrt(seg_dx^2 + seg_dy^2)
  #compute segment strike vectors
  seg_n <- cbind(dx=seg_dx, dy=seg_dy) / seg_len
  
  #compute average segment strike angles
  seg_s_ang <- as.double(atan2(seg_n[,2],seg_n[,1]))
  seg_s_ang <- seg_s_ang %% pi

  #average rupture strike angle
  avg_strike_ang <- as.double( weighted.mean(seg_s_ang, seg_len) )
  
  
  return(avg_strike_ang)
}

rup_length <- function(rup_data){
  
  #ensure that all rupture points belong to a single rupture
  assert_that( length(unique(rup_data$RUP_ID)) == 1 )

  #compute cartesian coordinates
  if ( !("x"%in% names(rup_data) & "y"%in% names(rup_data)) ){
    df_utmxy <- longlat2UTM(long = rup_data$Longitude, lat = rup_data$Latitude)[[1]]
    rup_data <- cbind(rup_data, df_utmxy)
  }
  
  #sort all rupture points based on the node id
  rup_data <- rup_data[order(rup_data$NODE_ID),]
  
  #compute rupture segment dx and dy intervals
  seg_dx <- diff(rup_data$x)
  seg_dy <- diff(rup_data$y)
  #compute segment lengths
  seg_len <- sqrt(seg_dx^2 + seg_dy^2)
  #compute rupture length
  rup_len <- sum(seg_len)
  
  return(rup_len)
}

compute_deriv_forward <- function(x_array, y_array){
  
  #number of points
  n_pt <- length(x_array)
  assert_that(length(x_array) == length(y_array))
  
  #compute differences
  dx_array <- c(diff(x_array, lag=1), x_array[n_pt]-x_array[n_pt-1])
  dy_array <- c(diff(y_array, lag=1), y_array[n_pt]-y_array[n_pt-1])
  
  #compute derivatives
  dydx_array <- dy_array/dx_array
  
  return(dydx_array)
}

compute_deriv_backward <- function(x_array, y_array){
  
  #number of points
  n_pt <- length(x_array)
  assert_that(length(x_array) == length(y_array))
  
  #compute differences
  dx_array <- c(x_array[2]-x_array[1], diff(x_array, lag=1))
  dy_array <- c(y_array[2]-y_array[1], diff(y_array, lag=1))
  
  #compute derivatives
  dydx_array <- dy_array/dx_array
  
  return(dydx_array)
}

compute_deriv_average <- function(x_array, y_array){
  
  #compute forward and backward derivatives
  dydx_fwrd <- compute_deriv_forward(x_array, y_array)
  dydx_bwrd <- compute_deriv_backward(x_array, y_array)

  #average derivative
  dydx_avg <- 0.5*(dydx_fwrd + dydx_bwrd)
  
  return(dydx_avg)
}

compute_deriv_central <- function(x_array, y_array){
  
  #number of points
  n_pt <- length(x_array)
  assert_that(length(x_array) == length(y_array))
  
  #compute differences
  dx_array <- c(x_array[2]-x_array[1], diff(x_array, lag=2), x_array[n_pt]-x_array[n_pt-1])
  dy_array <- c(y_array[2]-y_array[1], diff(y_array, lag=2), y_array[n_pt]-y_array[n_pt-1])
  
  #compute derivatives
  dydx_array <- dy_array/dx_array
  
  return(dydx_array)
}

compute_2deriv_central <- function(x_array, y_array){
  
  #number of points
  n_pt <- length(x_array)
  assert_that(length(x_array) == length(y_array))
  
  #compute differences
  d2y_array <-  diff(y_array[-1], lag=1) - diff(y_array[-n_pt], lag=1)
  dx2_array <- (diff(x_array, lag=2)/2)^2

  #compute derivatives
  d2ydx2_array <- d2y_array/dx2_array
  d2ydx2_array <- c(d2ydx2_array[1], d2ydx2_array, d2ydx2_array[n_pt])
  
  return(d2ydx2_array)
}


compute_curvature <- function(u_array, x_array, y_array){
  
  #compute derivatives
  d2xdu2_array <- compute_2deriv_central(u_array, x_array)
  d2ydu2_array <- compute_2deriv_central(u_array, y_array)
  
  #compute curvature
  curv_array <- sqrt(d2xdu2_array^2 + d2ydu2_array^2)
  
  return(curv_array)
}
