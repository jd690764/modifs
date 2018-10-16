#!/usr/bin/env Rscript
source( '/usr/local/share/R/lib.R' )
#source('/mnt/driveB/jdemeter/usr/R/lib.R')
options(warn=-1)
# create parser object
parser <- ArgumentParser()
# specify our desired options
# by default ArgumentParser will add an help option

parser$add_argument("-p", "--path", default = 'Ras',
                    help = 'fractionFiles folder name for the dataset' )
parser$add_argument("-k", '--keep', default = 'FALSE',
                    help = 'keep peptides that have no modifications in them.' )
parser$add_argument("-b", '--bait', default = 'KRAS',
                    help = 'gene symbol for the bait')
parser$add_argument("-o", '--offset', default = 65L,
                    help = "length of N-terminal tag to subtract from bait's coordinates" )
parser$add_argument("-w", '--overwrite', default = FALSE,
                    help = "overwrite *_all_data intermediate file?: TRUE/FALSE" )
parser$add_argument("-t", "--taxid", default = 9606,
                    help = 'taxid of the organism' )
parser$add_argument("-u", "--uid",
                    help = 'dproc.id of the experiment. If this is used, no other param is required.' )
parser$add_argument("-f", "--file",
                    help = 'csv list of dproc.id values. If this is used, no other param is required and a summary file is also generated.' )
parser$add_argument("-r", "--threshold", default = 0L,
                    help = "Threshold score to filter out bad quality peptides.")
parser$add_argument("-c", "--cuteff", default = FALSE,
                    help = "Calculate missed cuts?")


args          <- parser$parse_args()

datadir       <- "/srv/msrepo/fractionFiles/"
#datadir       <- "/mnt/driveB/jdemeter/msrepo/fractionFiles/"

files         <- character( )
groups        <- hash( ) # expt (dir) -> files

# set working directory
setwd(datadir)
makeDupMap()
ftsv        <-'/tmp/out.tsv'
outdir      <- paste0(datadir, 'PTMs')
matrixFile  <- '/mnt/driveB/jdemeter/ms_fraction_data/result.txt'

#args$path <- 'Ras'
#args$bait <- 'KRAS'
#args$keep <- FALSE
#args$offset <- '65'
#args$taxid <- '9606'
#args$overwrite <- TRUE
#args$uid <- '2038' # just use the id column value
#args$file <- '/usr/local/share/R/t.txt'
print( str(args))
if( class( args$keep ) == 'character' & nchar(args$keep) > 0 ){
    args$keep <- as.logical(args$keep)
}
if( class( args$offset ) == 'character' & nchar(args$offset) > 0 ){
    args$offset <- as.integer(args$offset)
}
if( class( args$threshold ) == 'character' & nchar(args$threshold) > 0){
    args$threshold <- as.integer(args$threshold)
}
if( class( args$overwrite ) == 'character' & nchar(args$overwrite) > 0 ){
    args$overwrite <- as.logical(args$overwrite)
}
if( class( args$taxid ) == 'character' & nchar(args$taxid) > 0 ){
    args$taxid <- as.integer(args$taxid)
}
if( 'uid' %in% names(args) & !is.null(args[['uid']])){

    if(grepl(',', args$uid)){
        uids_string <- gsub(' ', '', trimws(args$uid, which = 'b'))
        uids <- as.integer(str_split(uids_string, ',', simplify = T))
        for(uid in uids){
            args$uid <- uid
            # get params from the database table
            args <- retrieveArgs( args )

            l           <- findExpts( args, files, groups)
            files       <- l$files
            groups      <- l$groups
            makeSummary3( groups, outdir, args )
        }
    }
} else if( 'file' %in% names(args) & !is.null(args[['file']])){
    cwd  <- getwd()
    parfile <- args$file
    print(parfile)
    if( file.exists(parfile)){
        con  <- file(parfile, open = "r" )
        line <- readLines(con, n = 1)
        fname <- readLines(con, n = 1)
        map_desc <- as.data.table(str_split_fixed(readLines(con), '\\\t', n = 3))
        if(nrow(map_desc)> 0){
            colnames(map_desc) <- c('type', 'from', 'to')
            map_desc <- map_desc[type == 'map']
            if(nrow(map_desc) > 0){
                args$map_descr <- map_desc
            }
        }
        close(con)
        uids <- as.integer(str_split(line, ',', simplify = T))
        print(paste('uids: ', uids))
    }
    paths <- character()
    for( uid in uids ){
        args$uid <- uid
        args     <- retrieveArgs(args)
        paths    <- c(paths, args$path)
        l        <- findExpts( args, files, groups)
        files    <- l$files
        groups   <- l$groups
    }
    makeSummary3( groups, outdir, args )
    print( paste( 'paths: ', paths ))
    # combine datasets into a single report
    dset <- assembleDataset( paths, outdir, args )
    if(length(fname) > 0){
        fname <- paste0(outdir, '/', fname )
    }
    else{
        fname <- paste0(outdir, '/combined_modifications_report.tsv')
    }
    fwrite( dset, fname, sep = '\t')
} else {
    l           <- findExpts( args, files, groups)
    files       <- l$files
    groups      <- l$groups
    makeSummary3( groups, outdir, args )
}
#collectAllMods( files )

#makeMatrixFile( groups )
