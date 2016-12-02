rm(list=ls())
library(jsonlite)
library(magrittr)
library(ggplot2)
library(stringr)
library(assertthat)
library(dbscan)
library(kernlab)
setwd("~/Documents/sanger/paninisim")

#### how to import #####
# target directory contains a bunch of json files from PANINI
# the file names are read in, and then the paths are stored in a data.frame 
# with fields: 
# for each file (simulation):
#     the json is parsed and k is calculated
#     k-means clustering is done
#     for n <= k, what is max number of members of n in a single cluster
#     calculate percentage mis-clustered and store in data.frame
# visualise

genome.groups.from.id <- function(df,id="id") {
  # returns a vector the same length as the incoming data.frame
  return(apply(df[id],1,str_match,pattern="l(\\d+)g")[2,])
}

get.success.rate <- function(d, t, num.groups) {
  ### algorythm to determine success rate (loop over groups)
  ### t: data.frame with columns group & cluster
  ### returns float
  num.correct <- 0
  for (i in 0:(num.groups-1)) {
    t <- d %>% subset(group==i) %>% use_series(cluster) %>% table()
    print(t)
    num.correct <- num.correct + max(t)
    print(paste("i:",i,"correct:",max(t),"/",sum(t),sep=" "))
  }
  frac.correct <- num.correct / dim(d)[[1]]
  return(frac.correct)
}

do.kmeans.clustering <- function(d, num.groups) {
  # the clustering seems variable, perhaps do 10 times and pick top one?
  clustering <- kmeans(d[,c(2,3)], centers=num.groups, iter.max=1000, nstart=5)
  return(as.factor(clustering$cluster))
}

# do.dbscan.clustering <- function(d, num.groups) {
#   dbscan(as.matrix(d[,2:3]), eps = 0.1, minPts = num.groups)
# }

do.spectral.clustering <- function(d, num.groups) {
  x <- specc(x=as.matrix(d[,2:3]), centers=num.groups, kernel = "rbfdot", na.action = na.omit)
  return(as.factor(x@.Data))
}

analyse.json <- function(json.filename, spectral=FALSE) {
  d <- fromJSON(txt=json.filename)
  # d is a data.frame with 3 columns: id / x / y
  d$group <- genome.groups.from.id(d)
  num.groups <- length(unique(d$group))
  if (spectral==TRUE) {
    d$cluster <- do.spectral.clustering(d, num.groups)
  } else {
    d$cluster <- do.kmeans.clustering(d, num.groups)
  }
  success.rate <- get.success.rate(d, t, num.groups)
  print(success.rate)
  gg <- ggplot(data=d, aes(x=x,y=y,colour=group,shape=cluster)) + geom_point()
  return(list("success.rate" = success.rate, "gg" = gg))
}






folder.name <- "json/beta"
list.of.json.files <- list.files(path = folder.name, pattern = ".*\\.json")
master <- data.frame(n=integer(),
                     k=integer(),
                     s=integer(),
                     tau=integer(),
                     m=integer(),
                     mu_gain=integer(),
                     mu_loss=integer(),
                     r=integer(),
                     file.path=character(),
                     f=numeric(),
                     stringsAsFactors=FALSE) 
graphs <- list()


for (json.file in list.of.json.files) {
  # json.file <- list.of.json.files[[6]]
  json.path <- paste(folder.name,json.file,sep="/")
  params <- as.integer(str_extract_all(json.file,"(\\d+)")[[1]])
  print(params)
  assert_that(length(params) == 8)
  results <- analyse.json(json.path, spectral=FALSE)
  master[dim(master)[[1]]+1 , ] <- c(params, json.path, results$success.rate)
  graphs[[dim(master)[[1]]]] <- results$gg
}

master$f <- as.numeric(master$f)

### types in R are a pain
for (i in 1:8) {
  master[,i] <- as.factor(as.integer(master[,i]))
}
master[,10] <- as.numeric(master[,10])

ggplot(data=master, aes(x=k, y=f, colour=n, group=n)) + geom_line()







