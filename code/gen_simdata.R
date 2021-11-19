#################################
#### Generate simulated data ####
#################################

## Generate simulated data based on [World Input-Output Tables](http://wiod.org/database/wiots16).
## The raw data covers 44 economies (28 EU countries, 15 other major countries, and Rest of 
## the World) each of which has 35 sectors.
## The goal is to generate the input-output table covering 30 economies and 30 sectors.

## number of sectors per province
sec <- 30

## number of province
prov <- 30

## total number of nodes
tot <- sec*prov

province_name <- c("Beijing", "Tianjin", "Hebei", "Shanxi", "InnerMongolia", 
                   "Liaoning", "Jilin", "Heilongjiang", "Shanghai", "Jiangsu", 
                   "Zhejiang", "Anhui", "Fujian", "Jiangxi", "Shandong", 
                   "Henan", "Hubei", "Hunan", "Guangdong", "Guangxi", 
                   "Hainan", "Chongqing", "Sichuan", "Guizhou", "Yunnan", 
                   "Shaanxi", "Gansu", "Qinghai", "Ningxia", "Xinjiang")

tot_name <- paste0(rep(province_name, each = sec), 1:sec)

gen_sim <- function(wiot) {
  rownames(wiot) <- paste0(wiot$Country, wiot$RNr)
  rownames(wiot)[nrow(wiot)] <- "TOT"
  wiot <- wiot[-which(wiot$IndustryDescription == "Total intermediate consumption"), -c(1:5)]
  flow_row <- which(gsub("[A-z]", "", rownames(wiot)) %in% 1:sec)[1:tot]
  flow_col <- which(gsub("[A-z]", "", colnames(wiot)) %in% 1:sec)[1:tot]
  
  ## intermediate flow matrix
  flow <- wiot[flow_row, flow_col]
  wiot[flow_row, flow_col] <- 0
  ## final use
  fnl <- rowSums(wiot[, 1:(ncol(wiot)-1)])[flow_row]
  ## imports
  im <- colSums(wiot[1:which(rownames(wiot) == "ROW56"),])[flow_col]
  ## value added
  tva <- colSums(wiot[(which(rownames(wiot) == "ROW56")+1):(nrow(wiot)-1), ])[flow_col]
  
  result <- cbind(rbind(flow, im, tva, wiot[nrow(wiot), flow_col]), 
                  rbind(cbind(fnl, wiot$TOT[flow_row]), array(NA, dim = c(3, 2))))
  rownames(result) <- c(tot_name, "Imports", "Total value added", "Total input")
  colnames(result) <- c(tot_name, "Total.final.use", "Total.output")
  return(result)
}

if (! dir.exists("../data")) dir.create("../data")

for (year in c(2007L, 2012L)) {
  fname <- paste0("WIOT", year, "_October16_ROW.RData")
  load(file.path("wiot", fname))
  mriot <- gen_sim(wiot)
  Fname <- file.path("..", "data", paste0("data_", year, ".RData"))
  save(mriot, file = Fname)
  rm(mriot, wiot)
}
