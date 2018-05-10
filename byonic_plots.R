require(data.table)
require(readxl)
require(dplyr)
require(ggplot2)
require(scales)
require(stringr)

byonic_plots <- function( folder , rexp ){

    print(folder)
    file_pattern <- '^[A-Za-z0-1].*xlsx$'
    ds <- gsub('^.*\\/([^\\/]+)\\/?$', '\\1', folder)
    setwd(folder)

    # peptide file
    pept <- data.table()
    prot <- data.table()
    for(f in dir(pattern = file_pattern)){
        print(f)
        pr <- as.data.table(read_xlsx(f, sheet = 'Spectra'))
        pr[, fraction := f]
        pept <- bind_rows(pept, pr)
        pr <- as.data.table(read_xlsx(f, sheet = 'Proteins'))
        pr[, fraction := f]
        prot <- bind_rows(prot, pr)
    }

    colnames(pept)[11] <- 'merr' # mass error
    colnames(pept)[5] <- 'mods'
    colnames(pept)[19] <- 'prot'
    colnames(pept)[3] <- 'seq'

    pept[, fraction := gsub(rexp, '\\1', fraction)]
    pept <- pept[! grepl( '^>Reverse', prot)] # remove reverse hits

    # do peptides with modification have lower scores?
    png(filename = paste0(folder, ds, '_mod_vs_unmod_scores.png'), width = 1000, height = 1000)
    print(ggplot(pept, aes(x = is.na(mods), y = Score)) +
        geom_boxplot() +
        ggtitle(paste(ds, 'scores for modified vs unmodified peptides')))
    dev.off()

    t.test(pept$Score ~ is.na(pept$mods))
    # yes! clearly

    # scores along the experiments
    png(filename = paste0(folder, ds, '_scores_by_scan.png'), width = 2000, height = 2000)
    print(ggplot(pept, aes(x = `Scan Time`, y = Score, color = fraction)) +
        geom_point(alpha = 0.1) +
        stat_smooth(method = 'lm', se = F) +
        facet_wrap( ~ fraction) +
        ggtitle('Scores along scans by fraction'))
    dev.off()

    # mass errors along the experiments
    png(filename = paste0(folder, ds, '_merr_by_scan.png'), width = 2000, height = 2000)
    print(ggplot(pept, aes(x = `Scan Time`, y = merr, color = fraction)) +
        geom_point(alpha = 0.1) +
        stat_smooth(method = 'lm', se = F) +
        facet_wrap( ~ fraction) +
        ggtitle('Mass errors along scans by fraction'))
    dev.off()

    # distribution of incomplete cuts
    pept[, seqn := gsub('^[A-Z\\-]\\.(.+)\\.[A-Z\\-]$', '\\1', seq)]
    pept[, seqn := gsub('[^A-Z]', '', seqn)]
    pept[, krs := str_count(seqn, '[RK]')]
    cuts_long <- pept[, .N, .(krs, fraction)][order(krs)]
    cuts_long <- cuts_long[, .(cuts = N/sum(N), krs), fraction]

    png(filename = paste0(folder, ds, 'pepts_krs_by_fr.png'), width = 1000, height = 1000)
    ggplot(cuts_long, aes(x = krs, y = cuts, color = fraction)) + geom_line() +
        labs(title = 'Density of KRs in detected peptides per fraction', x = 'Number of KRs',
             y = 'density')
    dev.off()

    # number of peptides per fraction
    png(filename = paste0(folder, ds, '_pepts_by_fr.png'), width = 1000, height = 1000)
    print(ggplot(pept, aes(x = fraction)) + geom_bar() +
        ggtitle('Peptide counts by fraction')) +
        theme(text = element_text(size=20))
    dev.off()

    # protein sheet
    colnames(prot)[6] <- 'totint'
    colnames(prot)[7] <- 'spc'
    colnames(prot)[11] <- 'aa'

    prot[, fraction := as.numeric(gsub(rexp, '\\1', fraction))]
    prot <- prot[! grepl('^>Reverse', Description) & !is.na(totint)]

    # boxplot of protein length along fractions
    png(filename = paste0(folder, ds, '_aalen_by_fraction.png'), width = 2000, height = 1000)
    print(ggplot(prot, aes(x = as.factor(fraction), y = log10(aa))) +
              geom_boxplot() + geom_jitter(aes(color = log(spc)), alpha = 0.3, width =  0.25) +
              scale_color_gradient(low = muted('blue'), high = 'red') +
              theme(text = element_text(size=20)) +
              labs(title = 'Protein length in fractions', x = 'fraction'))
    dev.off()

    genes <- prot[,
                  sum(spc),
                  .(Description)][
                      !grepl('keratin|contam',
                             Description, ignore.case = TRUE), ][order(V1)] %>%
        tail(6) %>%
        select(Description)

    png(filename = paste0(folder, ds, '_gene_profiles_by_fraction.png'), width = 2000, height = 1000)
    print(
        ggplot(prot, aes(x = as.factor(fraction), y = log10(aa)-min(log10(aa)))) +
            geom_boxplot() +
            geom_line(data = prot[Description %in% genes$Description, ],
                      aes(x = fraction,
                          y = spc/(ceiling(max(spc)/(max(log10(prot$aa))-min(log10(prot$aa))))),
                          color = Description), size = 2) +
            theme(text = element_text(size=20)) +
            labs(title = 'Protein profiles across fractions (with boxplot)',
                 x = 'Fraction'))
    dev.off()


    # totint vs spectral count
    print(paste0(folder, ds, '_totint_vs_spc.png'))
    png(filename = paste0(folder, ds, '_totint_vs_spc.png'), width = 1000, height = 1000)
    print(ggplot(prot, aes(x = log2(totint), y = log2(spc))) +
        geom_point() + geom_smooth() +
        ggtitle('Sectral count vs total intensity'))
    dev.off()
}

byonic_plots('/mnt/driveB/jdemeter/msrepo/fractionFiles/Cep164_10090_redo/', '^.+cep164(\\d\\d).+$')
#byonic_plots('/mnt/driveB/jdemeter/msrepo/fractionFiles/RAB28_wt/', '^.+nm(\\d).+$')

