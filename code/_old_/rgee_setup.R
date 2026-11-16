# Following installation instructions from here: https://r-spatial.github.io/rgee/articles/rgee01.html

# more info here: https://r-spatial.github.io/rgee/

#install.packages(c("remotes", "googledrive", "rgee"))
#remotes::install_github("r-spatial/rgee")

library(rgee)

# Get the username
HOME <- Sys.getenv("HOME")

# ee install
ee_install(py_env = "rgee") # It is just necessary once!

library(rgee)

# authorize GEE and drive
ee_Initialize(user = 'samuelczipper@gmail.com', drive = TRUE)


# test
library(rgee)
ee_Initialize(user = 'samuelczipper@gmail.com', drive = TRUE)
srtm <- ee$Image("USGS/SRTMGL1_003")

viz <- list(
  max = 4000,
  min = 0,
  palette = c("#000000","#5AAD5A","#A9AD84","#FFFFFF")
)

rgee::Map$addLayer(
  eeObject = srtm,
  visParams =  viz,
  name = 'SRTM'
)

# check
ee_check()

rgee::ee_install_upgrade()
