#!/usr/bin/env Rscript
source( '/usr/local/share/R/lib.R' )

#setup parallel backend to use many processors
cores=detectCores()
cl <- makeCluster(cores[1]-4) #not to overload your computer
registerDoParallel(cl)

# script to get coverage data for all baits in the PTMs folder

dir_path <- '/mnt/driveB/jdemeter/msrepo/fractionFiles/PTMs/'
setwd( dir_path )

files    <- list.files(path = list.dirs(), pattern = '*_mods_summary.tsv', full.names = T)
sumfiles <- data.table(dir = gsub('\\.\\/(.+)\\/.*', '\\1', files), sumfile = gsub('\\.\\/.+\\/(.+)$', '\\1', files))
sumfiles[, sumfile := paste0(dir, '/', sumfile)]
sumfiles[, outfile := paste0(dir, '/', dir, '_cov_summary.tsv')]

for( i in 1:nrow(sumfiles)){
    print( sumfiles[i, dir])
    detCoverage(sumfiles[i, outfile], sumfiles[i, sumfile])
}

