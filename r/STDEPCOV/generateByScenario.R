# scp ~/shexport/PCORI/generateByScenario.R kikapp@sherlock:~/PCORI/lib
require(msm)
require(survival)
require(mvtnorm)
require(plyr)
require(MBESS)
require(data.table)

args <- commandArgs(trailingOnly = T)
print(args)

output_location <- args[1]
n_cores <- as.numeric(args[2])
range_low <- as.numeric(args[3])
range_high <- as.numeric(args[4])
source_location <- args[5]
scenario_location <- args[6]
scenario_row <- as.numeric(args[7])
# 
# output_location <- "~/shexport/PCORI/results/lim1_corr/localtest/data"
# n_cores <- 2
# range_low <- 1
# range_high <- 5
# source_location <- "~/shexport/PCORI/lib/"
# scenario_location <- "~/shexport/PCORI/lib/scenarios_corr_rv_2014_06_27.csv"
# scenario_row <- 2


if(is.na(output_location)) {output_location <- "/scratch/users/kikapp/dump"}
if(is.na(n_cores)) {n_cores <- 16}
if(is.na(range_low)) {range_low <- 0}
if(is.na(range_high)) {range_high <- 5}
if(is.na(source_location)) {source_location <- "~/PCORI/lib"}
if(is.na(scenario_location)) {scenario_location <- "~/PCORI/lib"}
if(is.na(scenario_row)) {scenario_row <- 1}


#register parallel backend
require(doSNOW)
cl<-makeCluster(n_cores)
registerDoSNOW(cl)

#generate worker ids in each worker
clusterApply(cl, seq(along=cl), function(id, scenario) WORKER.ID <<- paste0("job", scenario, "_core", id), scenario_row)
     
#### clean data
#rm(list=ls())

###function call and work directory
# setwd ("~/shexport/PCORI")
# source_location <- paste0(getwd(),"/shexport/PCORI/lib/")
# scenario_location <- paste0("~/shexport/PCORI/lib/")
# scenario_row <- 1
# output_location <- paste0(getwd(),"/shexport/PCORI/data/")
# source_location <- paste0(getwd(),"/lib/")
# output_location <- paste0("/scratch/users/kikapp/PCORI/data/lim1_large_betas_rv2")
print(list.files(source_location))

hextextToRaw <- function(text) {
  vals <- matrix(as.integer(as.hexmode(strsplit(text, "")[[1]])), ncol=2, byrow=TRUE)
  vals <- vals %*% c(16, 1)
  as.raw(vals)
}



scenarios <- read.table(scenario_location, sep = ",", header = T, stringsAsFactors = F)


#NUMBER OF OBSERVATIONS
n_subjects <- 1000
simrange <- range_low:range_high
nrow(unique(scenarios[, c("gm", "lim", "md", "w", "rv", "ns")]))
#combinations
scenario <- unique(scenarios[ , c("gm", "lim", "md", "w", "rv", "ns")])[scenario_row, ]
gm <- scenario$gm
lim <- scenario$lim
md <- scenario$md
w <- scenario$w
rv <- scenario$rv
ns <- scenario$ns
rho <- c(-0.8, -0.4, 0, 0.4, 0.8)[ns]
mediangoal <- c(35,75,150)[md]

trad_cens <- merge(x = scenario, y = scenarios, by = c("gm", "lim", "md", "w", "rv", "ns"))[ ,c("pc", "trad_nu", "trad_mdg")]
admin_cens <- merge(x = scenario, y = scenarios, by = c("gm", "lim", "md", "w", "rv", "ns"))[ ,c("pc", "admin_time")]

seed.init<-switch(lim,3182014,264064,856492)
comb <- as.numeric(paste0(scenario[, c("gm", "lim", "md", "w", "rv", "ns")], collapse = ""))

#file keeping run times
# sink("runtimes gm=1 lim=3.txt")
print(output_location)
source(paste0(source_location, "mysimRangeScenarioSB.R"))
mysim(seed.init, comb, n_subjects, gm, lim, w, md, rv, output_location, source_location, simrange, trad_cens, admin_cens, ns, mediangoal, rho)
