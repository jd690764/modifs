options(java.parameters = "-Xmx8182m")
options( warn = -1)
# suppressPackageStartupMessages( library('XLConnect' ))
suppressPackageStartupMessages( library('data.table'))
suppressPackageStartupMessages( library('sqldf'))
suppressPackageStartupMessages( library('argparse'))
suppressPackageStartupMessages( library( "hash" ))
suppressPackageStartupMessages( library( "R.utils" ))
suppressPackageStartupMessages( library(dplyr))
suppressPackageStartupMessages( library(DBI))
suppressPackageStartupMessages( library('parallel'))
suppressPackageStartupMessages( library('foreach'))
suppressPackageStartupMessages( library('doParallel'))
suppressPackageStartupMessages( library('stringr'))
suppressPackageStartupMessages( require(readxl))
suppressPackageStartupMessages( require(rio))
suppressPackageStartupMessages( library(tidyr))

pogodir       <- '/usr/local/share/bin/'
sheetToRead   <- 'Spectra'
matchLength   <- 120

spec_mods     <- hash('3HyAsn'   = c('Oxidation', 'N'),
                      '4HyPro'   = c('Oxidation', 'P') )
modifications <- hash('+71.037'  = 'Propionamide',
                      '+15.995'  = 'Oxidation',
                      '+31.990'  = 'Double Oxidation',
                      '+0.984'   = 'Deamidated',
                      '-17.027'  = 'Gln->pyro-Glu',
                      '-18.011'  = 'Glu->pyro-Glu',
                      '+114.043' = 'GlyGly(Ubiquitylation)',
                      '+42.011'  = 'Acetyl',
                      '+43.006'  = 'Carbamoylation',
                      '+14.016'  = 'Methyl',
                      '+28.032'  = 'Dimethyl',
                      '+28.031'  = 'Dimethyl',
                      '+42.011'  = 'Acetyl',
                      '+79.966'  = 'Phosphorylation',
                      '+226.078' = 'Biotinylation',
                      '+361.146' = 'APEX2-Biotinylation')

psi_mods    <- hash('+71.037'  = 'Propionamide',
                    '+15.995'  = 'Oxidation',
                    '+31.990'  = 'Double Oxidation',
                    '+0.984'   = 'Deamidated',
                    '-17.027'  = 'Gln->pyro-Glu',
                    '-18.011'  = 'Glu->pyro-Glu',
                    '+114.043' = 'GlyGly',
                    '+42.011'  = 'Acetyl',
                    '+43.006'  = 'Carbamoylation',
                    '+14.016'  = 'Methyl',
                    '+28.032'  = 'Dimethyl',
                    '+42.011'  = 'Acetyl',
                    '+79.966'  = 'Phospho',
                    '+226.078' = 'Biotin')

gygi_mods    <- hash('*'  = '+15.995',
                     '#'  = '+28.032',
                     '@'  = '+14.016')

options( warn = 0 )
dup     <- data.table( )
newOnly <- FALSE

refseqre  <- '^.*([NXY]P_\\d+)\\.?\\d?.*$'
uniprotre <- '^.*(sp|tr)\\|([^\\|]+)\\|.*$'
gnre      <- '^.+GN=([^ ]+) ?.*$'

ref_files  <- data.table( taxid = c(9606, 10090, 9606, 10090),
                          ref = c('fasta', 'fasta', 'gtf', 'gtf'),
                          path = c('input/gencode.v30.pc_translations.fa',
                                   'input/gencode.vM21.pc_translations.fa',
                                   'input/gencode.v30.annotation.gtf',
                                   'input/gencode.vM21.annotation.gtf' ))
pogo       <- data.table( taxid = c(10090, 9606), cmd = c('./PoGo_mm', './PoGo'))

dupFile  <- ifelse(Sys.info()['sysname'] == 'Windows',
                   'Z:/djscripts/data/pickles/dup_latest.txt',
                   ifelse(Sys.info()['sysname'] == 'Darwin', '
                          ~/mnt/driveB/jdemeter/djscripts/data/pickles/dup_latest.txt',
                          '/usr/local/share/py/djscripts/data/pickles/dup_latest.txt' ))

getDBcon <- function(){
    if(Sys.info()['sysname'] == 'Windows'){
        x <- suppressWarnings(fread('D:/jdemeter/my.cnf', sep = '=', skip = 1, header = F))
        y <- x$V2
        names(y) <- x$V1
        con_sql  <- dbConnect( RMySQL::MySQL(), user = y[['user']], password = y[['password']],
                               dbname = y[['database']], host = y[['host']] )
    } else {
        con_sql  <- dbConnect( RMySQL::MySQL(), group = "tcga" )
    }
}

retrieveArgs <- function( args ){

    con_sql  <- getDBcon()
    sql      <- paste0("select bait_symbol, tag_length, ff_folder, taxid
                       from network_dproc p, network_psample s
                       where p.expid_id = s.expid_id
                       and p.id = ", trim(args$uid) )
    #print(paste('sql=', sql))
    query    <- dbSendQuery( con_sql, sql )
    rs       <- as.data.table(fetch( query, n = -1 ))
    print(rs)
    dbClearResult( query )
    dbDisconnect(con_sql)
    if( nrow(rs) > 1){
        print( paste('More than one datasets were found in network_sample for uid = ', args$uid))
        exit()
    } else if( nrow(rs) == 0 ){
        print( paste(args$uid, 'was not found in the database.'))
        exit()
    } else if( rs$ff_folder == '' ){
        print( paste('no ff_folder was found for uid =', args$uid ))
        exit()
    }
    args$bait   <- rs$bait_symbol
    args$path   <- rs$ff_folder
    args$offset <- rs$tag_length
    args$taxid  <- rs$taxid
    return( args )
}

assembleDataset <- function( paths, outdir, args ){

    mods_all   <- data.table()
    allPos_all <- data.table()

    for( path in paths ){
        print(paste('assembleDataset: ', path))
        path   <- gsub('\\/?$', '', path)
        # assumed that the file is one of '_mods_summary.tsv" files
        file   <- paste0( outdir, '/', path, '/', path, '_mods_summary.tsv')
        pept   <- fread(file)[!is.na(pept), .(dset, pept, symbol, descr, start, end, nseq)]
        mods   <- fread(file)[ is.na(pept), c()]
        mods   <- fread(file)[ is.na(pept), -c('pept', 'start', 'end', 'nseq', 'score')]

        mods_all <- bind_rows(mods_all, mods)
        print('start countAllPositions')
        start.time <- Sys.time()
        allPos <- countAllPositions( pept )
        end.time <- Sys.time()
        print(paste('It took', end.time - start.time, 'sec'))
        allPos_all <- bind_rows(allPos_all, allPos)
    }

    # this table contains all observed modification:
    modifs <- mods_all[, .N, c('symbol', 'descr', 'modif', 'pos', 'residue')]
    modifs <- modifs[, N := NULL]

    #  get cartesian join with all experiment names
    dsets  <- data.table(dset = unique(mods_all$dset))
    modifs <- setkey(dsets[, c(k=1,.SD)],k)[modifs[,c(k=1,.SD)], allow.cartesian = T]
    modifs <- modifs[, k:= NULL]

    # now, add back the mod counts
    modifs <- as.data.table(left_join(modifs, mods_all, by = c('dset' = 'dset',
                                                               'symbol' = 'symbol',
                                                               'descr' = 'descr',
                                                               'modif' = 'modif',
                                                               'pos' = 'pos',
                                                               'residue' = 'residue'
                                                               )))

    # now join in allPos_all counts and add percentages
    modifs <- as.data.table(left_join(modifs, allPos_all, by = c('dset' = 'dset',
                                                                 'symbol' = 'symbol',
                                                                 'descr' = 'descr',
                                                                 'pos' = 'coords')))
    modifs[, percent_mod := round(100*mod_count/count,2)]

    # remove the 'all_count' column
    modifs[ , all_count := NULL ]

    # last, convert table from long to wide
    if( args$cuteff == TRUE ){
        cols <- c('mod_count', 'count', 'percent_mod', 'missed', 'unmod_missed')
    } else {
        cols <- c('mod_count', 'count', 'percent_mod')
    }

    modifs <- dcast(modifs, symbol + descr + pos + residue + modif ~ dset, value.var = cols, fun.aggregate = sum)

    # sort the columns by dataset
    cnames <- data.table(expts = dsets$dset, cnames = colnames(modifs)[6:ncol(modifs)])[order(expts)]
    setcolorder(modifs, c(colnames(modifs)[1:5], cnames$cnames))

    return( modifs )
}

assembleDataset_ori <- function( paths, outdir ){

    dataset <- data.table()
    for( path in paths ){
        print(paste('assembleDataset: ', path))
        path   <- gsub('\\/?$', '', path)
        # assumed that the file is one of '_mods_summary.tsv" files
        file   <- paste0( outdir, '/', path, '/', path, '_mods_summary.tsv')
        pept   <- fread(file)[!is.na(pept), .(dset, pept, symbol, descr, start, end, nseq)]
        mods   <- fread(file)[ is.na(pept), .(dset, symbol, descr, modif,
                                              pos, residue, mod_count, all_count)]
        print('start countAllPositions')
        start.time <- Sys.time()
        allPos <- countAllPositions( pept )
        end.time <- Sys.time()
        print(paste('It took', end.time - start.time, 'sec'))
        comb   <- unique(as.data.table(left_join(allPos, mods,
                                                 by = c('dset' = 'dset',
                                                        'symbol' = 'symbol',
                                                        'descr' = 'descr',
                                                        'coords' = 'pos'))))

        setcolorder(comb, c('dset', 'symbol', 'descr', 'coords', 'modif',
                            'residue', 'count', 'mod_count', 'all_count'))
        colnames(comb)[7:9] <- paste(colnames(comb)[7:9], comb[1,dset], sep = '.')
        comb[, dset := NULL ]
        if( nrow(dataset) == 0 ){
            dataset <- comb
        }
        else {
            dataset <- as.data.table(left_join(dataset, comb,
                                               by = c('symbol' = 'symbol',
                                                      'descr' = 'descr',
                                                      'coords' = 'coords',
                                                      'modif' = 'modif',
                                                      'residue' = 'residue')))
        }
    }
    return( dataset )
}

countAllPositions <- function( pept ){

    coords <- data.table()
    cores  <- min(detectCores()-1, 7)
    cl     <- makeCluster(cores[1], type = "FORK", outfile = '') #not to overload your computer
    registerDoParallel(cl)
    print(paste0('cores: ', cores))
    coords <- foreach(i=1:nrow(pept), .combine=rbind) %dopar% {
        #print(i)
        st <- pept[i, start]
        fi <- pept[i, end]
        seq_coords <- seq(st, fi, length.out = (fi - st + 1))
        tempTable = as.data.table(
            cbind(
                pept[i, .(dset, symbol, descr)],
                data.table(coords = seq_coords )
            )
        )

        tempTable
    }

    stopCluster(cl)

    coords <- coords[, .(count = .N), c('dset', 'symbol', 'descr', 'coords')]
    return( coords )
}

collectAllMods <- function( fs ){
    # read in the 'Spectra' sheets from each file
    # and keep column5 that contains the modifications
    # the mods are in column5 "Modification Type(s)"
    modifs <- character()
    fh <- file( "log.txt", 'wt' )
    for( f in fs ){
        #cont <- readWorksheetFromFile( f, sheetName = sheetToRead )
        cont <- read_excel(f, sheet = sheetToRead) %>% as.data.table
        mods <- as.character(levels(cont[,5]))
        for( m in mods ){
            if( grepl( '\\[', m, perl=TRUE )){
                kws <- strsplit( m, ', ' )
                for( kw in kws ){
                    kw <- gsub( '\\*[0-9]+$', '', kw, perl = TRUE )
                    if( grepl( '\\]\\[|\\.\\[|\\[\\d', kw )){
                        writeLines( paste0(f, '=', kw), fh )
                    }
                    modifs <- c( modifs, kw )
                }
            }
        }
    }
    writeLines( paste( as.character(levels(as.factor(modifs))), collapse="\n"), fh)
    close( fh )
}

compare_datasets <- function( filea, fileb, namea, nameb, outfile ){
    a <- fread( paste0('PTMs/', filea))
    b <- fread( paste0('PTMs/', fileb))

    all <- as.data.table(full_join(
        a[is.na(pept), -c('pept', 'dset', 'start', 'end', 'nseq')],
        b[is.na(pept), -c('pept', 'dset', 'start', 'end', 'nseq')],
        by = c('symbol', 'descr', 'modif', 'pos'))
    )
    moda <- paste0('mod_count.', namea)
    modb <- paste0('mod_count.', nameb)
    alla <- paste0('all_count.', namea)
    allb <- paste0('all_count.', nameb)
    colnames(all) <- c('symbol', 'descr', 'modif', 'pos',
                       moda, alla, modb, allb)

    all[complete.cases(all), pval_fisher := apply(
        all[complete.cases(all), c(moda, alla, modb, allb), with = FALSE], 1, function(x){
            y <- matrix(c(
                as.integer(x[1]), as.integer(x[2])-as.integer(x[1]),
                as.integer(x[3]), as.integer(x[4])-as.integer(x[3])),
                nrow = 2, byrow = T)
            p <- fisher.test(y)$p.value
            return(p)
        })]
    all[complete.cases(all[, -c('pval_fisher')]), pval_chisq := apply(
        all[complete.cases(all[, -c('pval_fisher')]), c(moda, alla, modb, allb), with = FALSE], 1, function(x){
            y <- matrix(c(
                as.integer(x[1]), as.integer(x[2])-as.integer(x[1]),
                as.integer(x[3]), as.integer(x[4])-as.integer(x[3])),
                nrow = 2, byrow = T)
            p <- chisq.test(y, simulate.p.value = TRUE)$p.value
            return(p)
        })]
    write.table(all[order(pval_fisher)], paste0('PTMs/', outfile), sep = '\t',
                quote = FALSE, row.names = FALSE)
}

currversion <- function(table = 'entrez'){
    con_sql <- getDBcon()
    q <- dbSendQuery( con_sql,
                      paste0('select id from network_version ',
                             "where tablename = '",
                             table, "' and inuse = 1"))
    vid <- as.integer(fetch(q, n = -1))
    dbClearResult( q )
    dbDisconnect( con_sql )
    return(vid)
}

dbAnnots <- function( ){

    evid    <- currversion()
    con_sql <- getDBcon()
    query   <- dbSendQuery( con_sql,
                            paste("select name ref, symbol_ref from (
                            select SUBSTRING_INDEX(SUBSTRING_INDEX(entrez.peptide, ';', numbers.n), ';', -1) name,
                            entrez.symbol symbol_ref, entrez.taxid, entrez.vid
                            from
                            numbers inner join entrez
                            on CHAR_LENGTH(entrez.peptide)
                            -CHAR_LENGTH(REPLACE(entrez.peptide, ';', ''))>=numbers.n-1
                            order by symbol
                            ) x
                            where name <> ''
                            and vid =", evid ))
    rs    <- as.data.table(fetch( query, n = -1 ))
    dbClearResult( query )
    query <- dbSendQuery( con_sql,
                          paste("select swissprot sp, symbol symbol_sp",
                                "from entrez where swissprot <> ''",
                                'and vid =', evid,
                                'union',
                                "select trembl, symbol from entrez",
                                "where trembl <> ''",
                                'and vid =', evid))
    sp    <- as.data.table(fetch( query, n = -1 ))
    dbClearResult( query )
    dbDisconnect(con_sql)
    return(list(sp=sp, rs=rs))
}

detCoverage1 <- function( outfile, sumfile ){
    print( 'Calculating coverage...' )
    pept <- fread( sumfile ) %>%
        .[!symbol %in% c('CONTAM'),
          acc := gsub('^>(.P_\\d+)(\\.\\d)? .*', '\\1', descr)] %>%
        .[!is.na(acc) & start != '' & end != '']
    pepuniq <- pept[, .(acc, symbol, start, end)] %>% unique %>%
        .[order(acc, symbol, start, end)]

    # get max isoform length for each protein accession
    if(sum(grepl('\\|', unique(pept$descr))) < length(unique(pept$descr))*0.9){
        # most likely an older annotation, sep = ' '
        # protein length is in ind = -3
        pept[, len := as.integer(gsub('^.* (\\d+) [^ ]+ [^ ]+$', '\\1', descr))]
    } else {
        # standard, sep = ' | '
        pept[, len := as.integer(gsub('^.* \\| (\\d+) \\| [^ ]+ \\| [^ ]+$', '\\1', descr))]
    }

    # this is a matrix to keep track of covered positions
    covm <- matrix(0,
                   nrow = length(unique(pept$acc)),
                   ncol = max(pept$len, na.rm = T))
    rownames(covm) <- pepuniq$acc %>% unique()

    sapply(unique(pepuniq$acc), function(xacc){
        # print(xacc)
        d <- pepuniq[acc == xacc, .(start = min(start)), end] %>%
            .[, .(end = max(end)), start] %>%
            .[order(start)]
        pos <- d[, apply(.SD, 1, function(x){
            seq(x[1], x[2])
        })] %>% unlist %>% as.integer %>% sort %>% unique %>%
            .[.>0] # when dealing with the bait, we could get negative indeces
        covm[xacc, pos] <<- 1
    })

    # count covered aas
    covm <- covm %>%
        as.data.table(keep.rownames = 'acc') %>%
        .[, cov := rowSums(.SD[, -c('acc')], na.rm = T)] %>%
        .[, .(acc, cov)]

    # coverage per gene symbol
    cov <- pept %>%
        left_join(covm, by = 'acc') %>%
        as.data.table %>%
        .[, cov_per := round(100*cov/len, 2)] %>%
        .[,
          .(uniq_pept = uniqueN(nseq), len = max(cov),
            len.full = mean(unique(len)), cov = max(cov_per)),
          .(symbol)] %>%
        .[, .(symbol, len, len.full, cov, uniq_pept)]

    fwrite(cov, outfile, sep = '\t')
}
detCoverage <- function( outfile, sumfile ){
    print( 'Calculating coverage...' )
    pept <- fread( sumfile )
    pept_uniq <- unique(pept[!is.na(pept), c(3,5:7), with=FALSE])
    pept_uniq <- pept_uniq[order(symbol, start, end)]

    # get max isoform length for each gene symbol
    con_sql  <- getDBcon()
    sql      <- paste("select e.symbol, max(len) len",
                      'from ncbiprots n, entrez e',
                      'where n.eid = e.eid',
                      'and e.vid =', currversion('entrez'),
                      'and n.vid =', currversion('ncbiprot'),
                      "group by e.eid")
    query    <- dbSendQuery( con_sql, sql )
    entr_len <- as.data.table(fetch( query, n = -1 ))
    dbClearResult( query )
    dbDisconnect(con_sql)

    # an attempt to deal with overlapping peptides
    uniq_pepts   <- data.table()
    discard_next <- FALSE
    ref_line     <- NA
    for( i in 1:(nrow(pept_uniq)-1)){
        if( discard_next == TRUE ){
            # don't save this one, but we should still carry
            # another row from earlier
            curr_line    <- ref_line
            ref_line     <- NA
            discard_next <- FALSE
        } else {
            curr_line <- pept_uniq[i]
        }
        if(pept_uniq[i+1,symbol] == curr_line$symbol){

            if(pept_uniq[i+1,start] == curr_line$start){
                # do nothing -> discard this line
            }
            else if(pept_uniq[i+1,start] < curr_line$end ){
                if( pept_uniq[i+1,end] <= curr_line$end){
                    # the next item is contained in this one -> throw it away
                    discard_next <- TRUE
                    # but we don't know yet what to do with the current line
                    # that depends on the row in i+2
                    ref_line     <- curr_line
                }
                else if(pept_uniq[i+1,end] > curr_line$end){
                    curr_line$end <- pept_uniq[i+1, start]
                    uniq_pepts    <- bind_rows(uniq_pepts, curr_line)
                }
            }
            else if(pept_uniq[i+1,start] >= curr_line$end){
                uniq_pepts <- bind_rows(uniq_pepts, curr_line)
            }
            else{
                print('should not be here...')
            }
        } else {
            uniq_pepts <- bind_rows(uniq_pepts, curr_line)
        }
    }

    coverage <- uniq_pepts[, .(len = sum(end-start+1)), symbol]

    coverage <- as.data.table(left_join(coverage, entr_len, by = 'symbol', suffix = c('', '.full')))

    # calculate coverage
    coverage[, cov := 100*len/len.full]

    # include number of unique peptides
    coverage <- as.data.table(left_join(coverage, pept_uniq[, .(uniq_pept = .N), symbol][order(uniq_pept)],
                                        by = 'symbol', suffix = c('', '.')))
    fwrite(coverage, outfile, sep = '\t')
}

findExpts <- function( args, files, groups) {

    if( length( args$path ) > 0 ){
        if( args$path %in% c( 'all', 'new' ) ){
            dirs    <- list.dirs( recursive = FALSE )
            newOnly <<- ifelse( args$path == 'new', TRUE, FALSE )
        } else {
            dirs    <- strsplit(args$path, ',')[[1]]
        }
    } else{
        stop( '\n    Please, provide a dataset/all to process!', call. = F )
    }
    for (d in dirs) {
        if( dir.exists( d ) && d != '.'){
            if ( grepl('^\\./', d, perl = TRUE) ) {
                expt <- sub( './', '', d )
            } else {
                expt <- d
            }
            groups[ expt ] <- character()
            df <- list.files( d )
            for (f in df) {
                if (grepl('^[0-9a-zA-Z].*.xlsx$', f)) {
                    if( grepl( 'HeatMap', f )){
                        next()
                    } else if( grepl('.*Blank.*', f)){
                        next()
                    }
                    f <- paste0(d, '/', f)
                    files <- c(files, f)
                    groups[expt] <- c(groups[[expt]], f)
                }
            }
        }
    }

    return( list( "groups" = groups, "files" = files, 'expt' = expt ))
}

findPos <- function(gd, dtype = 'sums_byonic'){
    print( 'collect modifications...' )
    if(dtype != 'sums_byonic'){
        print(system.time(
            gdata <- gd %>% left_join(
                bind_rows( apply(gd, 1, function(x){
                    # get the seqn without the modifications
                    findPosition( as.character(x['seqn_f']), as.character(x['start']),
                                  as.integer(x['pept']), 'data.table', TRUE, datatype )
                })), by = 'pept'
            ) %>% as.data.table))
    } else {
        gdata <- findPosition2(gd)
    }

    return(gdata)
}

findPosition2 <- function( dtab ){

    print('findPosition2 ...')
    #getmod <- function()
    search_str <- '\\[[^a-zA-Z]+\\]'
    dt <- data.table::copy(dtab)

    # sometimes there are modifications that are actually not there, e.g.:
    #     HIIVAC[+71.037][-71.037]EGNPYVPVHFDASV
    # remove these
    dt[grepl('\\[(\\+|\\-)([^a-zA-Z]+)\\]\\[(\\+|\\-)\\2\\]', seqn_f),
       gsub('\\[(\\+|\\-)([^a-zA-Z]+)\\]\\[(\\+|\\-)\\2\\]', '', seqn_f)]

    dt[, nseq := seqn_f]

    # starting aa is modified
    dt[grepl('^\\[', seqn_f),
       `:=`(residue = 'N-term',
            pos     = as.character(start),
            modif   = sapply(seqn_f,
                             function(x){
                                 x <- gsub( '[\\[\\]]',
                                            '',
                                            str_extract(x, search_str),
                                            perl = T)
                                 ifelse(x %in% keys(modifications),
                                        modifications[[x]],
                                        paste('unknown:', x))}),
            nseq    = str_replace(seqn_f, search_str, ''))]

    # capture all the others
    while(nrow(dt[grepl(search_str, nseq)]) > 0){
        dt[grepl('\\[', nseq),
           `:=`(residue = paste(residue,
                                str_replace(nseq, '^[^\\[]*(.)\\[.*$', '\\1'),
                                sep = ';'),
                modif   = paste(modif,
                                sapply(gsub( '[\\[\\]]',
                                             '',
                                             str_extract(nseq, search_str),
                                             perl = T),
                                       function(x){
                                           ifelse( x %in% keys(modifications),
                                                   modifications[[x]],
                                                   paste('unknown:', x))
                                       }),
                                sep = ';'),
                pos     = paste(pos,
                                as.character(
                                    nchar(str_replace(nseq,
                                                      '^([^\\[]+)\\[.*$',
                                                      '\\1')) - 1 + start),
                                sep = ';'),
                # update this column with the modification removed
                nseq    = str_replace(nseq, search_str, ''))]
    }

    dt <- dt %>%
        separate_rows(residue, pos, modif, sep = ';') %>%
        as.data.table %>%
        .[!(residue == 'NA' & grepl('\\[', seqn_f))] %>%
        .[, pos := as.integer(pos)]

    # update certain modifications:
    for( k in keys(spec_mods)){
        dt[residue == spec_mods[[k]][2] & modif == spec_mods[[k]][1],
           modif := k ]
    }
    print('back from findPosition2')
    return(dt)
}

findPosition <- function( seqn, strt, index, str = 'hash',
                          nseq = FALSE, dtype = 'sums_byonic' ){
    # this function determines about the modification its:
    #    residue (what amino acid)
    #    position ( where in the protein sequence)
    #    modification type (acetylation, ...)

    # it returns a hash, with the key(s) position
    # and values a vector containing residue and mtype
    #print(paste(seqn, strt, str, dtype))
    h <- list(a = hash(),
              b = data.table(residue = character(), pos = numeric(),
                             modif = character(), pept = integer() ))
    names(h) <- c('hash', str)

    search_str   <- '\\[[^a-zA-Z]+\\]'
    if( dtype == 'gygi' ){
        search_str <- '[^a-zA-Z]'
    }

    fragment     <- seqn
    #print(paste('ind=',index))
    while( any(grepl(search_str, fragment)) ){
        #print(fragment)
        # index of the first '[' -
        ind      <- str_locate( fragment, search_str )[1]
        if( ind == 1 ){
            residue  <- 'N-term'
            pos      <- as.character(as.integer( strt ))
        } else {
            residue  <- substr(fragment, ind-1, ind-1 )
            # the position of the modified residue in the protein
            pos      <- as.character(as.integer( strt ) + ind - 2 )
        }
        # this is the matched modificaton (e.g.: +42.011)
        modif    <- gsub( '[\\[\\]]', '', str_extract(fragment, search_str), perl = T)
        if( dtype == 'gygi' ){
            modif <- gygi_mods[[modif]]
        }
        if(! modif %in% keys(modifications)){
            modifications[modif] <<- paste0('unknown: ', modif)
        }

        # remove the modification (e.g.: [+42.011])
        fragment <- sub( search_str, '', fragment )

        # save the values
        h[['hash']][ pos ] <- c( residue, pos, modifications[[modif]] )
        h[[ str ]]         <- bind_rows(h[[str]],
                                        data.table(residue = as.character(residue),
                                                   pos = as.numeric(pos),
                                                   modif = as.character(modifications[[modif]]),
                                                   pept = index))
    }

    if( str != 'hash' && nrow(h[[str]]) == 0 ){
        h[[str]] <- data.table(residue = as.character(NA), pos = as.numeric(NaN),
                               modif = 'unmodified', pept = index)
    }

    if( nseq == TRUE ){
        h[[str]][, nseq := rep(fragment, nrow(h[[str]])) ]
    }

    # HO-proline and HO-Asn have the same diff as oxidation - they have to be fixed
    # as special case
    for( k in keys(spec_mods)){
        h[[str]][ residue == spec_mods[[k]][2] & modif == spec_mods[[k]][1], modif := k ]
    }

    return( h[[str]] )
}

findSymbol <- function( desc ){
    # the symbol can be got from rbase structure 'dup'
    # truncate desc to 'matchLength'
    d <- substr( desc, 1, matchLength )
    if( nrow(dup[desc == d])>0 ){
        symbol <- dup[ desc == d, symbol ][1]
    } else if( grepl( '.*contam.*', desc, ignore.case=TRUE, perl=TRUE ) ){
        symbol <- 'CONTAM'
    } else {
        symbol <- 'UNKNOWN'
        # make a last effort to find the correct symbol
        r   <- gsub(refseqre, '\\1', d)
        uni <- gsub(uniprotre, '\\2', d)

        # some of the uniprot accession have a -n extension
        uni <- gsub('-\\d+$', '', uni)

        if( nrow(dup[ ref == r,])>0 ){
            symbol <- dup[ ref == r, symbol ][1]
        } else if( nrow(dup[ sp == uni,])>0){
            symbol <- dup[ tolower(sp) == tolower(uni), symbol ][1]
        }
    }
    #    print( paste(d, symbol) )
    return( symbol )
}

fixAnnots <- function( d ){
    # get current symbols
    db <- dbAnnots()

    # parse out refseq and swissprot ids
    d[symbol == 'UNKNOWN' & grepl(uniprotre, descr), sp := gsub('-\\d+$', '', gsub(uniprotre, '\\2', descr))]
    d[symbol == 'UNKNOWN' & grepl(refseqre, descr), ref := gsub(refseqre, '\\1', descr)]

    # find the corresponding symbols
    d <- as.data.table( left_join( d, db$rs, by = 'ref'))
    d <- as.data.table( left_join( d, db$sp, by = 'sp'))

    # update to current symbols
    d[ symbol == 'UNKNOWN' & !is.na(symbol_ref), symbol := symbol_ref]
    d[ symbol == 'UNKNOWN' & !is.na(symbol_sp), symbol := symbol_sp]

    # if there are still UNKNOWNs and they have GN
    d[ symbol == 'UNKNOWN' & grep('GN=', descr), symbol := gsub(gnre, '\\1', descr)]

    # drop unneeded columns
    d[, c('ref', 'sp', 'symbol_ref', 'symbol_sp') := NULL]
    return( d )
}

fixSeq <- function( x ){
    # remove 'context' residues at the
    # begining and end of detected peptides
    x <- sub( '^[A-Z-]\\.', '', x, perl=TRUE )
    x <- sub( '\\.[A-Z-]$', '',  x, perl=TRUE )
    return( x )
}

increaseTotals <- function( symbol, seqn, strt, sumData ){
    # increase total counts for positions with
    # modified residues
    fragment <- seqn
    end      <- as.integer( strt ) + nchar( fragment ) - 1

    for( i in seq( as.integer( strt ), end )){
        key <- makeKey( as.character( i ), symbol )
        if( has.key( key, sumData )){
            sumData[[ key ]][ 5 ] <- as.integer(sumData[[ key ]][ 5 ]) + 1
        }
    }
}

makeDupMap  <- function(){

    # the dup.txt file contains the symbols
    # it has to be updated before looking
    # at modifications in a new dataset
    #cwd         <- getwd()
    #setwd( '/usr/local/share/py/' )
    #    setwd( '/mnt/driveB/jdemeter/usr/py/' )

    # run 'convertDup.py' to update 'dup.txt'
    #system( 'pwd' )
    #print( 'updating dup.txt file' )
    #system( 'python3 convertDup.py' )
    #setwd( cwd )

    # a text file is writen when dup is updated,
    # so the file needs to be loaded only
    readDupMap( )
}

makeKey <- function( pos, symbol ){
    # create a key from symbol and position
    return( paste( symbol, pos, sep=':' ))
}

makeMatrixFile <- function( gs ){
    sheet   <- 'Proteins'
    allData <- data.table(  expt  = character(0),
                            genes = character(0),
                            spnum = integer(0) )
    failed  <- character(0)
    for( g in hash::keys(gs) ){
        # if( g != '20160120_Rabl2B'){
        #        next()
        #    }
        print( paste0( "expt=", g ) )
        gData <- data.table(  expt = character(0),
                              desc = character(0),
                              spnum = integer(0) )
        fs    <- gs[[ g ]] # files of the expt

        if( grepl( '^.+\\/.+$', g )){
            # don't look at internal directories
            next()
        }

        for( f in fs ){
            print( f )

            # system( paste( python, pyParser, paste0( '"', f, '"' ), ftsv, sheet ))
            # fraction <- fread( ftsv, skip = 1 )
            # fraction <- data.table( readWorksheetFromFile( f, sheet = sheet ))
            fraction <- data.table( read_excel( f, sheet = sheet ))
            fraction <- fraction[, c(2,7), with = FALSE ]
            if( length(fraction) > 0 ){
                fraction           <- cbind( g, fraction )
                colnames(fraction) <- c('expt', 'desc', 'spnum')
                fraction           <- fraction[ !is.na( spnum )]
                gData              <- rbindlist( list(gData, fraction ), use.names = TRUE )
            }
            fraction <- NULL
            unlink( ftsv )
        }
        # get gene symbols from dup file
        gData[ , genes := sapply( desc, findSymbol )]
        # add up counts within same experiment by gene -
        # for UNKNOWNs and CONTAMs, collapse them based
        # on the desc field
        gData   <- data.table( sqldf('select expt, genes, sum(spnum) spnum
                                      from gData
                                      where genes not in ("CONTAM", "UNKNOWN")
                                      group by expt, genes
                                      union
                                      select expt, genes || "|" || desc genes, sum(spnum) spnum
                                      from gData
                                      where genes in ("CONTAM", "UNKNOWN")
                                      group by expt, desc, genes'))
        # append to all previous data
        allData <- rbindlist( list( allData, gData ), use.names = TRUE )
    }

    # reorganize data into a matrix with columns: experiments and rows: genes
    expts <- unique( allData[ , expt ])
    genes <- data.table( genes = unique( allData[ , genes ]))

    for( e in expts ){
        res <- data.table( sqldf( paste0('select d.genes, d.spnum
                                         from allData d, genes g
                                         where d.expt = "', e, '"
                                         and d.genes = g.genes')))
        #res  <- allData[ genes %in% genes$genes & expt == e, .(genes, spnum) ]
        genes[ genes %in% res$genes, c(e) := res$spnum, with = FALSE ]
    }
    write.table( genes, matrixFile, quote = FALSE, sep = "\t", row.names = FALSE )
}

makeSummary <- function( gs ){ # input is a hash of arrays
    # the data files are sorted by protein and peptide
    # read all files for a given experiment group
    for( g in hash::keys(gs) ){
        summaryFile <- paste0(outdir, '/', g, "_mods_summary.tsv")
        print( summaryFile )
        if( (newOnly && !file.exists( summaryFile) || !newOnly )){
            print( paste0( "expt=", g ) )
            fs    <- gs[[ g ]] # files of the expt
            gData <- data.frame(seqn = character(0),
                                mod = character(0),
                                start = numeric(0),
                                desc = character(0) )
            if( length( fs ) > 0 ){
                for( f in fs ){
                    print( f )
                    #system( paste( python, pyParser, paste0( '"', f, '"' ), ftsv, sheetToRead ))
                    #system( paste( 'perl -p -i -e "s/\n/ /g"', ftsv))
                    #system( paste( 'perl -p -i -e "s/\r /\n/g"', ftsv))

                    #fraction           <- as.data.frame(fread( ftsv ))
                    #if( nrow( fraction ) == 0 ){
                    #    Sys.sleep( 5 )
                    #    fraction       <- as.data.frame(fread( ftsv ))
                    #}

                    # fraction           <- readWorksheetFromFile( f, sheet = sheetToRead )
                    fraction           <- read_excel( f, sheet = sheetToRead ) %>% as.data.frame
                    # we only need columns # 3, 5, 12, 19
                    fraction           <- fraction[, c(3,5,12,19)]
                    colnames(fraction) <- c('seqn', 'mods', 'start', 'descr')
                    gData              <- rbind( gData, fraction )
                    rm( fraction )
                    unlink( ftsv )
                }

                # sort columns ( the N-term is counted as position 1, and the 1st amino acid
                # is counted as 2, so we decrease start by 1)
                gData          <- sqldf("select seqn, mods, start - 1 start, descr
                                        from gData
                                        where descr not like '%Reverse%'
                                        order by 4 desc, 3 asc, 2 desc")
                sumData        <- hash() # sumData['symbol:position'][res, pos, mtype, desc, totalcount, modcount]
                totalModSites  <- 0
                for( i in seq(nrow(gData))){
                    seqn   <- fixSeq( gData[ i, 1] )
                    mod    <- gData[ i, 2]
                    strt   <- gData[ i, 3]
                    desc   <- gData[ i, 4]
                    #cat( paste( seqn, mod, strt, desc, sep="\t"), "\n" )
                    symbol <- findSymbol( desc )
                    totalModSites <- totalModSites + numModSites( seqn )
                    if( nchar(mod) > 0 && grepl('\\[', seqn, perl=TRUE)){
                        info <- findPosition( seqn, strt )
                        # add this info to summary data
                        recordModifs( symbol, desc, info, sumData )
                    } else {
                        # increase totalcounts for relevant positions
                        increaseTotals( symbol, seqn, strt, sumData )
                    }
                }
                printResult( g, summaryFile, sumData, totalModSites, totalPeptCount )
            }
        }
    }
}

makeSummary2 <- function( gs, keep_unmodified_peptides = TRUE ){ # input is a hash of arrays
    # the data files are sorted by protein and peptide
    # read all files for a given experiment group
    for( g in hash::keys(gs) ){
        summaryFile   <- paste0(outdir, '/', g, "_mods_summary.tsv")
        all_data_file <- paste0(outdir, '/', g, '_all_data.txt')
        print( summaryFile )
        if( (newOnly && !file.exists( summaryFile) || !newOnly )){
            print( paste0( "expt=", g ) )
            fs    <- gs[[ g ]] # files of the expt
            gData <- data.frame(seqn = character(0), mod = character(0), start = numeric(0), desc = character(0) )

            if( length( fs ) > 0 ){
                if( ! file.exists(all_data_file)){
                    for( f in fs ){
                        print( f )

                        # fraction           <- readWorksheetFromFile( f, sheet = sheetToRead )
                        fraction           <- read_excel( f, sheet = sheetToRead ) %>% as.data.frame
                        # we only need columns # 3, 5, 12, 19
                        fraction           <- fraction[, c(3,5,12,19)]
                        colnames(fraction) <- c('seqn', 'mods', 'start', 'descr')
                        gData              <- rbind( gData, fraction )
                        rm( fraction )
                    }

                    gData                  <- data.table( gData )
                    print( 'lookup symbols...' )
                    # sort columns ( the N-term is counted as position 1, and the 1st amino acid
                    # is counted as 2, so we decrease start by 1)
                    gData                 <- gData[ !grepl( '.*Reverse.*', descr), .(seqn, mods, start, descr )]
                    gData[, start         := start -1 ][order(-descr, start, -mods)]
                    gData[, symbol        := sapply(descr, function(x){findSymbol(x)})]
                    gData                 <- fixAnnots( gData )
                    gData[, numModSites   := sapply(seqn, function(x){numModSites(x)})]
                    gData[, seqn_f        := sapply(seqn, function(x){fixSeq(x)})]
                    gData[, pept          := .I]
                    print( 'write table with all data...' )
                    write.table( gData, all_data_file, quote = FALSE, row.names = FALSE, sep = '\t' )
                } else {
                    gData                 <- fread( all_data_file )
                }

                print( 'collect modifications...' )
                gdata                <- bind_rows( apply(gData, 1, function(x){
                    #print(x)
                    y  <- findPosition( as.character(x[7]), as.character(x[3]), 'data.table' )
                    x  <- data.table(matrix(x, nrow = 1, dimnames = list(1, names(x))))
                    while( nrow(y) > nrow(x)){
                        x <- bind_rows( x, x[1])
                    }
                    y  <- bind_cols(as.data.table(x), as.data.table(y))
                }))

                print( 'do counts and save...' )
                totalModSites        <- length(unique(gdata[ !mods == '', pept]))
                totalPeptCount       <- max(gdata$pept)
                counts_by_symb_start <- gdata[, length(unique(pept)), by = .(symbol, descr, start)][order(symbol, start)]
                mod_counts           <- gdata[!is.na(modif), length(unique(pept)),
                                              by = .(symbol, descr, start, pos, modif, res)][order(symbol, start)]
                setnames(mod_counts, 'V1', 'counts')
                setnames(counts_by_symb_start, 'V1', 'counts')

                counts_by_symb_start[, pos := NaN ]
                counts_by_symb_start[, modif := 'all']
                counts_by_symb_start[, res := NA ]
                counts_by_symb_start[, percent := '100.00' ]

                mod_counts[, percent := as.numeric(apply(mod_counts, 1, function(x){
                    y <- (as.numeric(x['counts']) /
                              counts_by_symb_start[descr == x['descr'] & start == x['start'],
                                                   as.numeric(counts)] ) * 100
                    return(y)}))]
                mod_counts[, percent := format(round(percent, 2), nsmall = 2)]
                setcolorder(counts_by_symb_start, colnames(mod_counts))

                res    <- bind_rows( mod_counts, counts_by_symb_start)[order(symbol, descr, start, pos, modif)]
                if( keep_unmodified_peptides == TRUE ){
                    keep   <- unique(res[ , .(symbol, descr, start) ])
                } else {
                    keep   <- unique(res[ !is.na( pos ), .(symbol, descr, start) ])
                }
                # remove unmodified peptides
                res    <- data.table( sqldf( 'select r.symbol, r.descr, r.start, r.pos position,
                                                     r.modif modification, r.res residue, r.counts,
                                                     r.percent from res r, keep k
                                             where r.symbol = k.symbol
                                             and r.start = k.start
                                             and r.descr = k.descr
                                             order by r.symbol, r.start, r.pos, r.modif'))
                # move peptide count in new column instead of separate row
                res    <- data.table( sqldf( 'select r1.symbol, r1.descr, r1.start,
                                                     r1.position, r1.modification, r1.residue,
                                                     r1.counts, r2.counts, r1.percent
                                             from res r1, res r2
                                             where r1.symbol = r2.symbol and r1.descr = r2.descr
                                             and r1.start = r2.start
                                             and r2.modification = "all"
                                             and r1.modification not in ("all")
                                             order by r1.symbol, r1.descr, r1.start,
                                                      r1.position, r1.modification'))
                colnames(res)[7]  <- paste0( 'mod_counts (all mod = ', totalModSites, ')')
                colnames(res)[8]  <- paste0( 'total_counts (all pept = ', totalPeptCount, ')')
                if( nrow(res[ symbol == args$bait]) > 0 & args$offset > 0 ){
                    res[ symbol == args$bait, start := as.integer(start) - args$offset ]
                    res[ symbol == args$bait, position := as.integer(position) - as.integer(args$offset) ]
                    res <- res[ position > -1, ]
                }
                # add dset name as the first column
                as.data.table(bind_cols(data.table(dset = rep(g, nrow(res))), res))

                write.table( res[ !(modification == 'unmodified' & percent < 100) & symbol != 'CONTAM',],
                             summaryFile, quote = FALSE, row.names = FALSE, sep = '\t' )
            }
        }
    }
    }

makeSummary3 <- function( gs, outdir, args ){ # input is a hash of arrays

    # create a report by position of modification (not by start of reported peptide as in makeSummary2)

    for( g in hash::keys(gs) ){
        subdir        <- paste0(outdir, '/', g)
        if( !dir.exists(subdir)){
            dir.create(subdir)
        }
        summaryFile   <- paste0(subdir, '/', g, "_mods_summary.tsv")
        summaryDFile  <- paste0(subdir, '/', g, "_mods_summary_data.tsv")
        all_data_file <- paste0(subdir, '/', g, '_all_data.txt')
        pogo_file     <- paste0(subdir, '/', g, '_all_data_pogo.txt')
        cov_file      <- paste0(subdir, '/', g, "_cov_summary.tsv")
        print( paste('summaryFile: ', summaryFile ))
        datatype      <- 'sums_byonic'

        if( (!file.exists( summaryFile) || args$overwrite )){
            print( paste0( "expt=", g ) )
            fs    <- gs[[ g ]] # files of the expt
            gData <- data.frame(seqn = character(0), mod = character(0),
                                start = numeric(0), desc = character(0) )

            if( length( fs ) > 0 ){
                if( ! file.exists(all_data_file) | args$overwrite == TRUE){
                    for( f in fs ){
                        print( f )
                        fraction           <- as.data.table(
                            read_excel( f, sheet = sheetToRead ))
                        #print(paste('colnames: ', colnames(fraction)))
                        if('Trimmed.Peptide' %in% colnames(fraction)){
                            # the data are from the Gygi lab
                            print( 'from gygi lab')
                            fraction           <- fraction[, c('Reference', 'Trimmed.Peptide',
                                                               'Start.Position', 'LDA.Score', 'MissedCleav')]
                            datatype           <- 'gygi'
                            colnames(fraction) <- c('descr', 'seqn', 'start', 'score', 'missed')
                            fraction[, mods := NA ]
                            setcolorder(fraction, c('seqn', 'mods', 'start', 'score', 'descr', 'missed'))
                            # remove designated contaminants
                            fraction           <- fraction[! grepl('_contaminant$', descr)]
                        }
                        else {
                            # the data are from sums
                            # we only need columns # 3, 5, 12, 14, 19
                            fraction           <- fraction[, c(3,5,12,14,19)]
                            colnames(fraction) <- c('seqn', 'mods', 'start', 'score', 'descr')
                            # find the number of missed cleavages
                            fraction[, missed := gsub('[RK]?$', '', gsub('^.?\\.(.+)\\..?$', '\\1', seqn))]
                            fraction[, missed := str_count(missed, '[RK]')]
                        }
                        # remove peptides that are below a threshold score
                        fraction           <- fraction[score >= args$threshold,]

                        fraction           <- bind_cols( data.frame(fraction = rep(f, nrow(fraction))),
                                                         as.data.frame(fraction) )
                        gData              <- rbind( gData, fraction )
                        rm( fraction )
                    }

                    gData                  <- data.table( gData )
                    # remove reverse hits
                    gData                 <- gData[ !grepl( '.*Reverse.*', descr),
                                                    .(fraction, seqn, mods, start, score, missed, descr )]
                    # remove contaminants
                    gData                 <- gData[ !grepl('contaminant', descr)]

                    # change descr if we have a need for that
                    if( 'map_descr' %in% names(args) ){
                        md <- args$map_descr
                        for( i in 1:nrow(md)){
                            print(paste('replacing', md[i,from], 'with', md[i,to]))
                            gData[descr == md[ i, from ], descr := md[ i, to ]]
                        }
                    }

                    print( 'lookup symbols...' )
                    print('findSymbol ...')
                    gData[, symbol        := sapply(descr, findSymbol)]
                    print('fixAnnots...')
                    gData                 <- fixAnnots( gData )
                    print('numModSites...')
                    gData[, numModSites   := sapply(seqn, numModSites, dtype = datatype)]
                    print('fixSeq...')
                    gData[, seqn_f        := sapply(seqn, fixSeq)]
                    gData[, pept          := .I]
                    print( 'write table with all data...' )
                    fwrite( gData, all_data_file, sep = '\t' )
                    if( 'usepogo' %in% names(args) & args$usepogo == TRUE ){
                        outputForPogo( gData, g,  pogo_file )
                        print( 'run pogo...')
                        runPogo( pogo_file, args$taxid )
                    }
                } else {
                    gData                 <- fread( all_data_file )
                }

                gData[, fraction := NULL ]
                gdata <- findPos( gData )

                print( 'do counts and save...' )

                # end coordinate of peptide
                gdata[, end := as.integer(start) + nchar(nseq) - 1]
                # group data by descr, modification and position - count uniqe modifications
                res <- gdata[modif != 'unmodified',
                             .(mod_count = uniqueN(pept), score = mean(as.numeric(score), na.rm = T),
                               missed = mean(as.numeric(missed), na.rm = T)),
                             by = .(symbol, descr, modif, pos, residue)]
                # count the number of unique peptides that cover each modification position
                res[, all_count := apply(res, 1, function(x){
                    gdata[as.integer(x['pos']) >= as.integer(start) &
                              as.integer(x['pos']) <= as.integer(end) &
                              x['descr'] == descr,
                          uniqueN(pept)]
                })]

                # mean 'missed' for unmodified peptides
                res[, unmod_missed := apply(res, 1, function(x){
                    gdata[modif == 'unmodified' & as.integer(x['pos']) >= as.integer(start) &
                              as.integer(x['pos']) <= as.integer(end) &
                              x['descr'] == descr,
                          mean(as.integer(missed))]
                })]

                res <- bind_rows(gdata[, uniqueN(pept),
                                       .(pept, symbol, descr, start, end, nseq)][,-'V1'],
                                 res)

                # correct coordinates of the baits for the N-terminal fragment of the tag
                # the endogeneous proteins might be expressed, or there might be
                # contamination from N/C-term tagged version run just before ...
                sym  <- args$baits[[g]]
                oset <- as.integer(args$offsets[[g]])
                desc <- res[grepl(sym, descr) & !grepl('^>.P_', descr), # assume refseq library
                            .N, descr][order(-N)][1, descr] # the desc field for the actual bait
                if( args$offsets[[g]] > 0 ){
                    res[ descr == desc, start := as.integer(start) - oset ]
                    res[ descr == desc, end := as.integer(end) - oset ]
                    res[ descr == desc, pos := as.integer(pos) - oset ]
                    res <- res[ pos > -1 | is.na(pos), ]
                }
                # add dset name as the first column
                res <- as.data.table(bind_cols(data.table(dset = rep(g, nrow(res))), res))

                fwrite( res, summaryFile, sep = '\t' )
                dt <- res[is.na(pept), .(dset, symbol, descr, modif, pos, residue, score, missed,
                                         unmod_missed, mod_count, all_count)]
                fwrite( dt, summaryDFile, sep = '\t' )
            }
        }
        if( (!file.exists( cov_file) || args$overwrite )){
            print(system.time(
            detCoverage1( cov_file, summaryFile )
            ))
        }
    }
}

numModSites <- function( seqn, dtype = 'sums_byonic' ){
    # determine the number of modified sites in
    # a sequence
    s <-  seqn
    if( dtype == 'sums_byonic'){
        s <- gsub( '\\[', '', s )
    }
    else if( dtype == 'gygi' ){
        s <- gsub( '[^a-zA-Z]', '', s )
    }
    numMods <- nchar( seqn ) - nchar( s )
    #    print(paste(seqn, s, numMods))
    return( numMods )
}

outputForPogo <- function( g, dset, file ){
    gc <- g
    gc <- gc[, .('Peptide' = seqn_f)]
    for( k in keys(psi_mods) ){
        v = psi_mods[[k]]
        gc[grepl(paste0('\\', k), Peptide), Peptide := gsub( paste0('\\[\\', k, '\\]'), paste0('\\(', v, '\\)'), Peptide)]
    }
    res <- gc[, .('PSMs' = .N), Peptide]
    res[ , Experiment := rep(dset, nrow(res))]
    #res[ , Quant := rep(1, nrow(res))]
    res[ , Quant := PSMs ]
    res <- res[, c('Experiment', 'Peptide', 'PSMs', 'Quant')]
    setcolorder(res, c('Experiment', 'Peptide', 'PSMs', 'Quant'))
    fwrite(res, file, sep = '\t' )
}

printResult <- function( expt, summaryFile, sumData, totalModSites, totalPeptCount ){
    # print the summary file
    fh <- file( summaryFile, open='wt' )
    writeLines( paste( "symbol", "description", "position", "residue",
                       "modification", "total count", paste0("modified count(total=", totalModSites, ")"),
                       "% modified", sep="\t"), fh, sep="\n" )
    for( k in sort(hash::keys( sumData ))){
        symb <- splitKey( k )
        writeLines( paste( symb,
                           sumData[[k]][4], #desc
                           sumData[[k]][2], #position
                           sumData[[k]][1], #residue
                           sumData[[k]][3], #modification
                           sumData[[k]][5], #total cound
                           sumData[[k]][6], #mod count
                           100*as.integer(sumData[[k]][6])/as.integer(sumData[[k]][5]), # % mod
                           sep="\t"
        ), fh, sep="\n")
    }
    close( fh )
    system( paste( 'chgrp', 'mscurators', summaryFile ))
}

readDupMap  <- function(){
    dup      <<- fread( dupFile )
    colnames(dup) <<- c('desc', 'eid', 'symbol', 'taxid')
    dup[, ref := ifelse(grepl(refseqre, desc), gsub(refseqre, '\\1', desc), NA)]
    dup[, sp := ifelse(grepl(uniprotre, tolower(desc)), gsub(uniprotre, '\\2', tolower(desc)), NA)]
}

recordModifs <- function( symbol, desc, info, sumData ){
    # create a new modified site in the protein,
    # increase it's count if it already exists
    for( pos in hash::keys(info) ){
        key <- makeKey( pos, symbol )
        if( has.key( key, sumData )){
            sumData[[ key ]][5] <- as.integer(sumData[[ key ]][5]) + 1 # total peptide count
            sumData[[ key ]][6] <- as.integer(sumData[[ key ]][6]) + 1 # modified peptide count
        } else {
            sumData[ key ] <- c( info[[ pos ]][seq(3)], desc, 1, 1 )
        }
    }
}

runPogo <- function( pfile, txid ){
    cwd <- getwd()
    setwd( pogodir )
    command <- paste( pogo[taxid == txid, cmd], '-fasta', ref_files[taxid == txid & ref == 'fasta', path ],
                     '-gtf', ref_files[taxid == txid & ref == 'gtf', path],
                     '-in', pfile)
    print( paste('pogo cmd: ', command) )
    system( command, wait = FALSE)
    setwd( cwd )
}

splitKey <- function( key ){
    # split a key and return the symbol
    ks <- unlist( strsplit( key, ':' ))
    return( ks[1] )
}

trim <- function (x){
    # remove leading and trailing spaces
    gsub("^\\s+|\\s+$", "", x)
    gsub("^'*|'*$", "", x)
}

create_summary_file <- function( dir, outfile ){
    # creates a SUMS-like summary file from fraction files in a folder
    library(openxlsx)
    if(!dir.exists(dir)){
        print( paste0('Please, check the directory: ', dir))
        return()
    }
    if(length(list.files(path = dir, pattern = '^[a-zA-Z0-9]+.*\\.xlsx$')) == 0){
        print( paste( dir, 'contains no excel files'))
        return()
    }
    ffiles <- list.files(path = dir, pattern = '^[a-zA-Z0-9]+.*\\.xlsx')
    count  <- 1
    sum_data <- data.table()
    cnames <- c('Rank Number', 'Protein Name')
    for(ff in ffiles){
        print( paste('read', ff))
        cnames <- c(cnames, ff)
        #wb <- loadWorkbook(paste0(dir, '/', ff))
        #ff_data <- as.data.table(readWorksheetFromFile(paste0(dir, '/', ff), sheet = 'Proteins'))
        ff_data <- read_excel(paste0(dir, '/', ff), sheet = 'Proteins') %>% as.data.table()
        #ff_data <- as.data.table(read.xlsx(paste0(dir, '/', ff), sheet = 'Proteins', cols = c(2,7)))
        # only care about columns 'Description' and '# of spectra' (2, 7)
        ff_data <- ff_data[, c(2,7),with = FALSE]
        colnames(ff_data) <- c('Protein Name', ff)
        # discard Reverse rows
        ff_data <- ff_data[!grepl('^>Reverse.*', `Protein Name`)]
        # discard rows with no values
        ff_data <- ff_data[!is.na(eval(as.name(ff)))]
        if( count == 1 ){
            colnames(ff_data)[2] <- ff
            count <- count + 1
            sum_data <- ff_data
        } else {
            sum_data <- as.data.table(full_join(sum_data, ff_data, by = 'Protein Name', suffix = c('', '.')))
        }
    }
    sum_data[, MAX := apply(sum_data[,2:ncol(sum_data)], 1, max, na.rm = T)]
    sum_data[, SUM := apply(sum_data[,2:(ncol(sum_data)-1)], 1, sum, na.rm = T)]
    cnames <- c(cnames, 'MAX', 'SUM')
    sum_data[, 'Rank Number' := rank(-MAX, ties.method = 'first')]
    sum_data <- sum_data[order(`Rank Number`)]
    setcolorder(sum_data, cnames)
    #wb <- XLConnect::loadWorkbook(outfile, create = TRUE)
    #createSheet( wb, name = 'GLOBAL')
    #writeWorksheet(wb, sum_data, sheet = 'GLOBAL' )
    #saveWorkbook(wb)
    #writeWorksheetToFile(outfile, sum_data, 'GLOBAL' )
    write.xlsx(sum_data, outfile, sheetName = 'GLOBAL', row.names = FALSE)
    return( sum_data )
}

#sf <- create_summary_file( '/mnt/driveB/jdemeter/msrepo/fractionFiles/Fuz_10090_redo/',
#                            '/mnt/driveB/jdemeter/msrepo/rawfiles/Fuz_10090_redo.xlsx')
#sf <- create_summary_file( '/mnt/driveB/jdemeter/msrepo/fractionFiles/Intu_10090_redo/',
#                            '/mnt/driveB/jdemeter/msrepo/rawfiles/Intu_10090_redo.xlsx')
#sf <- create_summary_file( '/mnt/driveB/jdemeter/msrepo/fractionFiles/Wdpcp_10090_redo/',
#                            '/mnt/driveB/jdemeter/msrepo/rawfiles/Wdpcp_10090_redo.xlsx')



#outdir <- '/mnt/driveB/jdemeter/msrepo/fractionFiles/PTMs/'
#paths  <- c('Anks3_10090_redo_oh', 'Anks6_10090_NT_redo_oh', 'Anks6_10090_R823W_redo_oh',
#              'Anks6_10090_CT_redo_oh', 'Nek8_10090_G442V_redo_oh', 'Nek8_10090_redo_oh',
#              'Nek8_10090_H425Y_redo_oh', 'Nek8_10090_K33M_redo_oh')
#all_data <- assembleDataset(paths, outdir)
#setwd('/mnt/driveB/jdemeter/msrepo/fractionFiles/')
#args <- list(bait = 'Trp53', file = NA, keep = FALSE, offset = 0, overwrite = TRUE,
#              path = 'Trp53_qm_v3', taxid = 9606, uid = 2054, threshold = 0)
#groups <- hash()
#groups$Trp53_qm_v3 <- dir(args$path, full.names = TRUE)
#args$map_descr <- data.table(type = 'map', from = '>pLAP7/Puro_p53_QM', to = '>pLAP7/Puro_p53')
#makeDupMap()
#makeSummary3(groups, outdir, args)
#makeSummary3(zgroups, zoutdir, zargs)
#pth <- '~/mnt/driveB/jdemeter/usr/R/'
#makeSummary3(readRDS(paste0(pth, 'groups.rds')), pth, readRDS(paste0(pth, 'args.rds')))