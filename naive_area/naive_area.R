library(raster)
library(rgdal)

# For reading filenames
landsatRead <- function(filenames){
  landsatImgs = list()
  k = 1
  for(i in filenames){
    temp = stack(i)
    landsatImgs[k] = temp
    k  = k + 1
  }
  return(landsatImgs)
}

ndsi_calc <- function(green, swir){
  numer = green - swir
  denom = green + swir
  return (numer / denom)
}

ndwi_calc <- function(nir, swir){
  numer = nir - swir
  denom = nir + swir
  return (numer / denom)
}

# For adding NDSI/NDWI to the rasterStack
add_layers <- function(filenames){
  rasterStack=landsatRead(filenames)
  for (i in 1:(length(rasterStack) - 1)){
    green = rasterStack[[i]][[2]]
    nir = rasterStack[[i]][[4]]
    swir = rasterStack[[i]][[5]]
    ndsi = ndsi_calc(green, swir)
    ndwi = ndwi_calc(nir, swir)
    rasterStack[[i]] <- addLayer(rasterStack[[i]], ndsi, ndwi)
  }
  return(rasterStack)
}

# For getting a single threshold
single_thresh <- function(rasterStack, r, row){
  dict = list("blue" = 1, 
              "green" = 2, 
              "red" = 3,
              "nir" = 4,
              "swir" = 5,
              "ti" = 6,
              "ndsi" = 7,
              "ndwi" = 8)
  if (row[["greater"]]){
    return(rasterStack[[r]][[dict[[row[["band"]]]]]] > row[["value"]])
  }
  return(rasterStack[[r]][[dict[[row[["band"]]]]]] < row[["value"]])
}

# For combining multiple thresholds
multi_thresh <- function(rasterStack, r, thresholds, plot=FALSE) {
  # start with first layer
  curr_rast = single_thresh(rasterStack, r, thresholds[1,])
  # `&` the rest of the layers, TODO: vectorize
  if (nrow(thresholds) > 1) {
    for (i in 2:(length(thresholds) - 1)) {
      temp = single_thresh(rasterStack, r, thresholds[i,])
      curr_rast = curr_rast & temp
    }
  }
  # create matrix for pixel area calculation
  combinedMat = as.matrix(curr_rast)
  combinedMat[is.na(combinedMat)] <- 0
  if (plot) {
    return (list("raster" = curr_rast, 
                 "pixelArea" = sum(combinedMat)))
  }
  return (sum(combinedMat))
}

decimaldate <- function(date){
  
  year = as.Date(date)
  year = format(year,"%Y")
  year = as.numeric(as.character(year))
  first = as.Date(paste(year, "-01-01", sep = ''))
  
  diff = as.numeric(as.Date(date) - first)
  
  number_of_days = ifelse(year %% 4 == 0, 366, 365)
  decimal=diff/number_of_days
  
  decimal.date = year + decimal
  
  return(decimal.date)
}

timeSeries <- function(filenames, rasterStack, threshs) {
  dates = substr(filenames,1,nchar(filenames)-4)
  dates = dates[-length(dates)]
  dates = decimaldate(dates)
  
  areas = c()
  for (i in 1:(length(rasterStack) - 1)) {
    areas <- c(areas, multi_thresh(rasterStack, i, threshs))  
  }
  return (data.frame(dates, areas))
} 

# Thresholds
threshs = data.frame(band=NaN, value=NaN, greater=NaN)
threshs[1, ] = list("band"="red", "value"=0.1, "greater"=TRUE)
threshs[2, ] = list("band"="nir", "value"=0.2, "greater"=TRUE)
threshs[3, ] = list("band"="ndsi", "value"=0.4, "greater"=TRUE)
threshs[4, ] = list("band"="ndwi", "value"=0.2, "greater"=TRUE)

df = data.frame(id=NaN, start_date=NaN, end_date=NaN, slope=NaN, slope_error=NaN, intercept=NaN, intercept_error=NaN, adj_r2=NaN)
print('created empty df')

# ids = list.files("/Users/DarrenPC/Desktop/Landsat/new")
ids = c("G014028E66722N")

print('loaded IDs')

for (i in 1:length(ids)) {
  # load in rasterStack
  id = ids[i]
  setwd(paste("C:/Users/DarrenPC/Desktop/Landsat/data/", id, sep=""))
  filenames=list.files(getwd())
  rasterStack=add_layers(filenames)
  print(paste("Loaded ID", id))  

  # create time series
  ts = timeSeries(filenames, rasterStack, threshs)
  ts$time = ts$dates - min(ts$dates)
  
  # linear model
  lm = lm(areas ~ time + time*cos(2*pi*time) + time*sin(2*pi*time) + time*cos(4*pi*time) + time*sin(4*pi*time), 
          data=ts)
  sumlm = summary(lm)
  
  df[i, ] = list("glac_id"=id,
                 "start_date"=floor(min(ts$dates)), 
                 "end_date"=floor(max(ts$dates)), 
                 "slope"=sumlm$coefficients[2, 1],
                 "slope_error"=sumlm$coefficients[2, 2],
                 "intercept"=sumlm$coefficients[1, 1],
                 "intercept_error"=sumlm$coefficients[1, 2],
                 "adj_r2"=sumlm$adj.r.squared)
  
  plotlm = lm(areas ~ dates + dates*cos(2*pi*dates) + dates*sin(2*pi*dates) + dates*cos(4*pi*dates) + dates*sin(4*pi*dates), 
              data=ts)
  sumplotlm = summary(plotlm)
  plot(ts$dates, ts$areas, main=id, xlab='Year', ylab='Area (km^2)', type='l')
  abline(sumplotlm$coefficients[1:2], col="blue")
  print(paste("Finished ID", id))
}
#setwd("/datasets/home/home-03/66/966/dsl030/Landsat")
#write.csv(df,'r_df.csv')
