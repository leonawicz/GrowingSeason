setwd("C:/github/GrowingSeason_BACKUP/data")

f <- function(x, proj4=NA){
	require(raster)
	r <- raster(x)
	if(is.na(projection(r)) & !is.na(proj4)) projection(r) <- proj4
	con <- file(x, "rb")
	on.exit(close(con))
	v <- readBin(con, integer(), size=2, n=ncell(r), endian="little")
	r <- setValues(r, v)
}

proj4 <- "+proj=lcc +lat_1=50 +lat_0=50 +lon_0=-107 +k_0=1 +x_0=0 +y_0=0 +a=6370997 +b=6370997 +units=m +no_defs"

files <- list.files(pattern=".bil$", full=T)
#b <- brick(lapply(files, f, proj4))
r <- f(files, proj4)

#writeRaster(b, file.path(outDir,"pct20_tdd_spring_1979_2010.tif"))
writeRaster(r, "clim_tdd_1979_2010.tif")
