################################################
##### DATA IMPORT, CLEANING, AND WRANGLING #####
################################################

#' Tidy encoded patient + sample metadata and normalizes the time variable for measurements
#' 
#' @param seq.type string; whether to return measurements corresponding to 16S or WGS features
#' @param mapfiles named list; kicked up to `clean.sample.map()` to unify across data types
#' @param patient.data data.frame; raw .tsv of static patient metrics and characteristics
#' @param sample.data data.frame; raw .tsv of logitudinal clinical patient metrics
#' 
#' @returns a tidy data frame with header cols `DataType, Site, PatientID, Visit, Time, Days` + measurements
#' where `Time` has been normalized to the baseline/days since first positive COVID test.
clean.clinical.data <- function(seq.type, mapfiles, patient.data, sample.data) {
    
    # browser()
    # 45 x 2 - check
    # TODO checks for patient data
    baseline.norm <- select(patient.data, PatientID, Positive) 
    
    # checks for sample map(s)
    if (all(exists('sample.map.16s', mapfiles),
            exists('sample.map.wgs', mapfiles),
            exists('sample.map.met', mapfiles))) {
        raw.sample.map <- clean.mapfiles(ret.type = 'all',
                                         sample.id = FALSE,
                                         mapfiles) } 
    else { stop('Missing one or more mapfiles, check input data.') }
    
    # combine all the sample info from data types
    # n=563 from belen's powerpoint, n=561 here
    # distinct() brings it down to n=559... 
    # C19-CB-164 has two TBS samples on one day/visit 
    sample.map <- raw.sample.map %>%
        # will be NA for controls... 
        left_join(baseline.norm) 
    
    # 562 (should be 561 after metabolomics data fix) samples!
    tmp <- sample.map %>%
        nest(data = -PatientID) %>%
        mutate(Baseline = map2(data, PatientID, function(.x, .y) {
            # V1 missing; taken from Belen's .xlxs date of symptom onset
            if (.y == 'C19-CB-197') '2020-05-29' 
            else {
                .x %>% 
                    dplyr::filter(Visit == 'V1') %>% 
                    arrange(Date) %>%
                    pull(Date) %>% 
                    extract2(1) %>%
                    as.character() }})) %>%
        unnest(Baseline) %>%
        distinct(across(matches('PatientID','Baseline')), .keep_all = TRUE) %>%
        mutate(Baseline = lubridate::ymd(Baseline)) %>%
        unnest(data) %>%
        mutate(Time = Date-Baseline-Positive) %>%
        dplyr::rename(Site = 'Measurement') %>%
        select(DataType, Site, PatientID, Visit, Time) %>%
        arrange(DataType, Site, PatientID, Visit)
    
    # load valid clinical data -- 366 measurements for 45 unique patient-visit 
    clinical.data.tmp <- sample.data %>%
        tidyr::gather('PatientID','Values', -c(Measurement, Days)) %>%
        dplyr::filter(!is.na(Values)) %>%
        spread('Measurement','Values') %>%
        select(-WGS) %>% # this X-VX means nothing/just confuses
        pivot_longer(cols = c(SDNA,OP,TBS,Plasma,Urine), 
                     names_to = 'Treatment', 
                     values_to = 'Visit') %>% 
        dplyr::filter(!is.na(Visit)) %>%   # here's where we have 366 rows
        dplyr::rename(Site = 'Treatment') %>%
        select(Site, PatientID, Visit, everything()) %>%
        nest(data = -Site) %>%
        mutate_at(vars(Site), ~ fct_recode(., P = 'Plasma', U = 'Urine')) %>%
        unnest(data) %>%
        mutate(Site = as.character(Site)) # %>%
    # results in 45 unique PatientID-Visits -- correct
    # distinct(across(matches('PatientID','Visit')), .keep_all=TRUE)
    
    clinical.data <- tmp %>%
        left_join(clinical.data.tmp, by = c('Site' = 'Site',
                                            'PatientID' = 'PatientID',
                                            'Visit' = 'Visit')) %>%
        dplyr::filter(DataType %in% c(seq.type, 'MET'))
}


#' Filters clinical data for WGS site of relevance; prep for statistical analysis
#' 
#' @param df.clin data.frame; fully expanded clinical data for all WGS
#' @param site string; which WGS-associated body site to filter for
#' @param osci.disc logical; coarse-grained OSCI categories T/F
#' @param days.disc logical; coarse-grained Days (early or late) T/F
#' 
#' @returns a tidy data frame 
filter.clinical.data <- function(df.clin, site, days.disc = TRUE) {
    # browser()
    tmp <- df.clin %>%
        filter(Site == site) %>%
        select(-c(Site, PatientID, Visit, Time)) %>%
        mutate(PatientID = str_split(UniqueID, '-') %>% 
                   map_chr(~ head(., 1))) %>%
        mutate(UniqueID = str_remove_all(UniqueID, paste0('-', site))) %>%
        # operationalize controls
        mutate_at(vars('OSCI_Class_Samp'), ~ replace_na(., 'Control')) %>%
        mutate_at(vars('OSCI_Score_Samp'), ~ replace_na(., '0')) %>%
        # remove VL and scRNA variables later updated in separate files.. 
        select(-c(ISG, IFNa, IL28, INFa, Leukocytes, PCR), 
                   # -starts_with('VL'),
                   -contains('Abs'),
                   -contains('Pct')) 
    
    ## Early or late timepoints only
    if (days.disc) {
        tmp <- tmp %>%
            mutate_at(vars('Days'), 
                      ~ fct_collapse(.,
                                     `Early` = c('5-7','8-10'),
                                     `Late` = c('14-16','11-13','17-19','26-31',
                                                '32-37','20-22','38-43','44-49',
                                                '23-25'))) }

        tmp %>%
            select(UniqueID, GenID, PatientID, 
                   Days, contains('OSCI'), everything()) %>%
            mutate_at(vars(c(7:21, 'OSCI_Score_Samp')), ~ as.numeric(.))
}


#' Reads and cleans up sample-donor mapfiles; utility fn for `load.clinical.data()`
#' 
#' @param ret.type string; whether to `bind_rows` on the different data types 
#' (default = 'all') or return just one type (options: '16S', 'WGS', 'MET')
#' @param sample.id logical; whether to return SampleID (not available for MET)
#' @param mapfiles named list; raw .tsv files about sampling information and Patient/SampleIDs
#' 
#' @returns a tidy data frame 
clean.mapfiles <- function(ret.type = 'all', sample.id = TRUE, mapfiles = list()) {
    
    # browser()
    if (length(mapfiles) == 0) 
        stop('List of sample-donor mapfiles is empty! Please check input.')
    
    # 16S sample data: n=194
    sample.map.16s <- mapfiles$sample.map.16s %>%
        select(c('SampleID','Donor_No','Treatment','Taking_Date','Visite')) %>%
        dplyr::rename(Date = 'Taking_Date') %>%
        dplyr::rename(PatientID = 'Donor_No') %>%
        dplyr::rename(Visit = 'Visite') %>%
        mutate(Date = lubridate::mdy(Date)) %>%
        filter(!is.na(Visit)) %>% # gets rid of blanks/water 
        add_column(DataType = '16S') 
    # shotgun sample data: n=194
    sample.map.wgs <- mapfiles$sample.map.wgs %>%
        # mutate_at(vars(SampleID), ~ paste0('S', .)) %>%
        mutate(Date = lubridate::dmy(Date)) %>%
        mutate(SampleID = paste0('S', SampleID)) %>%
        add_column(DataType = 'WGS') %>%
        # gets rid of blanks/water
        filter(!is.na(Visit))
    # urine & plasma sample data: n=173, NO SAMPLE IDS YET
    sample.map.met <- mapfiles$sample.map.met %>%
        select(-Original) %>%
        mutate(Date = lubridate::dmy(Date)) %>%
        dplyr::filter(Treatment %in% c('Plasma', 'Urine')) %>%
        nest(data = -Treatment) %>%
        mutate_at(vars(Treatment), ~ fct_recode(., P = 'Plasma', U = 'Urine')) %>%
        unnest(data) %>%
        mutate(Treatment = as.character(Treatment)) %>%
        add_column(DataType = 'MET')
    
    # return specified mapfile(s)
    switch(ret.type, 
           'all' = {
               if (all(exists('sample.map.16s'),
                       exists('sample.map.wgs'),
                       exists('sample.map.met'))) {
                   tmp <- sample.map.16s %>%
                       bind_rows(sample.map.wgs) %>%
                       bind_rows(sample.map.met) %>%
                       dplyr::rename(Measurement = 'Treatment') 
                   
               } else { stop('Could not parse all mapfiles, please check input.') }},
           '16S' = { 
               if (exists('sample.map.16s')) {
                   tmp <- dplyr::rename(sample.map.16s,
                                        Measurement = 'Treatment') }
               else { stop('Could not parse 16S mapfile, please check input.') }},
           'WGS' = { 
               if (exists('sample.map.wgs')) {
                   tmp <- dplyr::rename(sample.map.wgs,
                                        Measurement = 'Treatment') }
               else { stop('Could not parse WGS mapfile, please check input.') }},
           'MET' = { 
               if (exists('sample.map.met')) {
                   tmp <- dplyr::rename(sample.map.met,
                                        Measurement = 'Treatment') }
               else { stop('Could not parse all mapfiles, please check input.') }})
    
    # return SampleID (y/n for 16S & WGS, not available for metabolite data)
    if (sample.id) {
        if (ret.type == 'MET') 
            select(tmp, PatientID, Measurement, Date, Visit, DataType)
        else
            select(tmp, SampleID, PatientID, Measurement, Date, Visit, DataType) }
    else { select(tmp, PatientID, Measurement, Date, Visit, DataType) }
}


clean.master.mapfile <- function(list.in, clean = TRUE) {
    
    if (clean) {
        list.in %>%
            map(~ mutate(., Date = lubridate::ymd(Date))) %>%
            map_dfr(~ mutate(., Note = as.character(Note)), .id = 'Site') %>%
            mutate(GenID = paste0(PatientID,'-',Visit)) %>%
            mutate(DataType = ifelse(Site %in% c('P','U'),'MET','WGS')) %>%
            mutate(AltID = case_when(str_detect(Note,'x2') ~ paste0(GenID,'.2'),
                                     str_detect(Note,'V') ~ paste0(PatientID,'-',Note))) %>%
            pivot_longer(cols = c(GenID, AltID)) %>%
            dplyr::filter(!is.na(value)) %>%
            select(DataType, Site, PatientID, value, Date) %>%
            rename(GenID = 'value')
    } else { # find discrepancies and potential mislabels
        list.in %>%
            map(~ mutate(., Date = lubridate::ymd(Date))) %>%
            map_dfr(~ mutate(., Note = as.character(Note)), .id = 'Site') %>%
            dplyr::filter(!is.na(Note)) %>%
            mutate(GenID = paste0(PatientID,'-',Visit)) %>%
            mutate(DataType = ifelse(Site %in% c('P','U'),'MET','WGS')) %>%
            mutate(AltID = case_when(str_detect(Note,'x2') ~ paste0(GenID,'.2'),
                                     str_detect(Note,'V') ~ paste0(PatientID,'-',Note))) %>%
            pivot_longer(cols = c(GenID, AltID)) %>%
            # ADD these samples later
            dplyr::filter(name == 'AltID') %>%
            select(DataType, Site, PatientID, value, Date) %>%
            rename(GenID = 'value')
    }
}


combine.sample.maps <- function(clin.map, visit.map, discreps, short.ids) {
    
    # browser()
    
    # contains P11-V2.2 and P17-V4.2
    visit.map <- visit.map %>% #alt.map
        left_join(short.ids, by = c('PatientID'='ShortID')) %>%
        rename(LongID = 'PatientID.y') %>%
        select(DataType,Site,PatientID,GenID,LongID,Date) %>%
        mutate(GenID = ifelse(str_detect(PatientID,'K') & str_detect(GenID,'V2'),
                              paste0(PatientID,'-V7'),
                              GenID))
    
    map.tmp <- clin.map %>%
        tidyr::gather('PatientID','Values', -c(Measurement, Days)) %>%
        left_join(short.ids) %>%
        rename(LongID = 'PatientID', PatientID = 'ShortID') %>%
        dplyr::filter(!is.na(Values)) %>%
        spread('Measurement','Values') %>%
        select(LongID, PatientID,Days,OP,SDNA,TBS,Urine,Plasma,WGS) %>%
        pivot_longer(cols= -c(LongID,PatientID,Days)) %>%
        rename(Site = 'name', Visit = 'value') %>%
        mutate(Site = str_replace_all(Site, 'Urine','U')) %>%
        mutate(Site = str_replace_all(Site, 'Plasma','P')) %>%
        mutate(DataType = ifelse(Site %in% c('P','U'),'MET','WGS')) %>%
        filter(!is.na(Visit)) %>%
        mutate(GenID = paste0(PatientID,'-',Visit)) %>%
        select(DataType,Site,PatientID,GenID,LongID,Days) %>%
        filter(Site != 'WGS')
    
    tmp1 <- select(visit.map, -Date) %>% arrange(GenID)
    tmp2 <- select(map.tmp, -Days) %>% arrange(GenID)
    waldo::compare(tmp1$GenID, tmp2$GenID)
    # raw.otu.feat%>%select(1:10) %>% filter(str_detect(PatientID,'CB-119'))
    # raw.metabolites %>% filter(str_detect(Patient.ID,'CB-119')) %>% distinct(across(c(Patient.ID,Visit)))
    
    # take TBS dates from OP
    # keep x2 samples
    # remove discrepancies (keep what is saved in discreps df)
    full <- map.tmp %>%
        full_join(visit.map) %>%
        mutate(Unique = paste0(DataType,'-',Site,'-',GenID)) %>%
        filter(!Unique %in% c('MET-P-P05-V2', # V1 and V3 only
                              'WGS-OP-P05-V2', # V1 and V3 only
                              'WGS-SDNA-P15-V4', # V2 and V3 only
                              'MET-P-P16-V3', 'WGS-OP-P16-V3',
                              'MET-P-P22-V3','WGS-OP-P22-V3'))
    probs <- full %>%
        arrange(Unique) %>% 
        filter(is.na(Days) | is.na(Date)) %>% 
        filter(Site!='TBS')
    
    full %>% filter(PatientID %in% probs$PatientID)
    
    tmp.otu <- raw.otu.feat %>% 
        select(1:5) %>%
        left_join(short.ids) %>%
        mutate(GenID = paste0(ShortID,'-',Visit))
    tmp.met <- raw.metabolites %>%
        distinct(across(c(Patient.ID,Visit,Matrix))) %>%
        rename(PatientID='Patient.ID') %>%
        left_join(short.names) %>%
        mutate(ShortID = ifelse(str_detect(PatientID,'K'),
                                PatientID, ShortID)) %>%
        mutate(GenID = paste0(ShortID,'-',Visit))
    select(full,Site,GenID,Days,DataType) %>% 
        anti_join(tmp.otu) %>% 
        filter(DataType=='WGS')
    select(full,Site,GenID,Days,DataType) %>% 
        anti_join(tmp.met) %>% 
        filter(DataType=='MET')

}


#' Attaches appropriate sample data to raw OTU feature data frame
#' 
#' @param seq.type string; either '16S' or 'WGS' for amplicon or shotgun data
#' @param raw.file data.frame; raw .tsv data corresponding to `seq.type`
#' @param mapfile data.frame; single raw mapfile to pass one to `clean.sample.names()`
#' @param unmapped logical; iff WGS, whether unmapped or mapped reads should be returned
#' 
#' @returns a tidy data frame with header cols `Site, PatientID, Date, Visit, SampleID` + OTUs
#' corresponding to `seq.type` (either mOTUs or SILVA 138 OTUs)
clean.otu.data <- function(seq.type, raw.file, mapfile, unmapped = FALSE) {
    
    # browser()
    switch(seq.type,
           '16S' = {
               raw.file %>%
                   tidyr::gather('SampleID','abundance', -OTU) %>%
                   spread('OTU','abundance') %>%
                   clean.sample.names(mapfile, seq.type) %>%
                   select(-SampleID) },
           'WGS' = { 
               if (unmapped) {
                   raw.file %>%
                       tidyr::gather('SampleID','abundance', -X1) %>%
                       spread('X1','abundance') %>%
                       clean.sample.names(mapfile, seq.type) }
               else {
                   raw.file %>%
                       dplyr::filter(X1 != '-1') %>%
                       tidyr::gather('SampleID','abundance', -X1) %>%
                       spread('X1','abundance') %>%
                       clean.sample.names(mapfile, seq.type) }})
}


#' Reads and cleans up sample-donor mapfiles; utility fn for `load.otu.data()`
#' 
#' @param seq.type string; either '16S' or 'WGS' for amplicon or shotgun data
#' @param raw.file data.frame; raw .tsv data corresponding to `seq.type`
#' @param mapfile data.frame; single raw mapfile to pass one to `clean.sample.names()`
#' @param unmapped logical; iff WGS, whether unmapped or mapped reads should be returned
#' 
#' @returns a tidy data.frame
clean.sample.names <- function(df, mapfile, seq.type) {
    
    # browser()
    switch(seq.type, 
           '16S' = {
               map <- mapfile %>%
                   select(c('SampleID','Donor_No','Treatment','Taking_Date','Visite')) %>%
                   rename(Date = 'Taking_Date') %>%
                   rename(PatientID = 'Donor_No') %>%
                   rename(Visit = 'Visite') %>%
                   rename(Site = 'Treatment') %>%
                   mutate(Date = lubridate::mdy(Date)) 
               
               df %>%
                   left_join(map, by = 'SampleID') %>%
                   filter(!str_detect(PatientID, 'BLANK')) %>%
                   select(Site, PatientID, Date, Visit, everything()) },
           'WGS' = {
               map <- mapfile %>%
                   mutate(SampleID = paste0('S', SampleID)) %>%
                   rename(Site = 'Treatment') %>%
                   mutate(Date = lubridate::dmy(Date))
               
               df %>%
                   left_join(map, by = 'SampleID') %>%
                   dplyr::filter(!str_detect(PatientID, 'BLANK')) %>%
                   select(Site, PatientID, Date, Visit, everything()) 
               })
}


#' Joins static patient data with short IDs and factorizes appropriately
#' 
#' @param raw.data data.frame; raw .tsv of patient data
#' @param short.ids data.frame; map of long and short PatientIDs
#' 
#' @returns a tidy data.frame
clean.patient.data <- function(raw.data, short.ids) {
    
    raw.data %>%
        mutate_at(vars(Sex, Smoker, Immunosuppressives, Antivirals, 
                       Abx3Mo, AbxSamp, AbxFreeHosp), ~ as.factor(.)) %>%
        select(-Antifungals) %>%
        left_join(short.ids) %>%
        select(PatientID, ShortID, everything())
    
}


#' Futher cleans patient data to prepare for metadeconfoundR; OSCI variables
#' 
#' @param df.pat data.frame; cleaned static patient.data
#' 
#' @returns a tidy data.frame
filter.patient.data <- function(df.pat) {
    
    # browser()
    df.pat %>%
        select(-c(Positive, Pregnancy)) %>%
        mutate_at(vars('OSCI_Class'), ~ replace_na(., 'Control')) %>%
        mutate_at(vars('OSCI_Class'), ~ fct_collapse(.,
                                                     `Control` = 'Control',
                                                     `Mild_Moderate` = c('Asymptomatic','Mild','Moderate'),
                                                     `Severe_Critical` = c('Severe','Critical'))) %>%
        mutate_at(vars(DaysHospitalized), ~ ifelse(. == 'deceased', NA, .)) %>%
        mutate(DaysHospitalized = as.numeric(DaysHospitalized))
}


#' Rectifies the rtk bug that inserts an colshift of NAs in the first column
#' 
#' @param raw.file data.frame; raw .tsv data corresponding to `seq.type`
#' @param preloaded logical; whether or not this step has been completed
#' 
#' @returns a tidy data.frame
clean.functional.data <- function(raw.file, preloaded = FALSE) {
    
    if (preloaded) 
        return(vroom::vroom(here('data','proc-data','rarefied-gene-table.tsv')))
    else if (file.exists(here('data','proc-data','rarefied-gene-table.tsv')))
        return(vroom::vroom(here('data','proc-data','rarefied-gene-table.tsv'))) 
    else {
        cn <- colnames(raw.file)
        tmp <- raw.file %>%
            select(-S97) %>%
            separate(S198, into = c('tmp1','tmp2'), sep = '\t') %>%
            set_names(cn)
        vroom::vroom_write(tmp, 
                           here('data','proc-data','rarefied-gene-table.tsv'))
        return(tmp)
    }
}


#' Joins the gmgc map to the abundance table for later filtering/renaming/binning
#' 
#' @param df.func data.frame; raw .tsv data corresponding to `seq.type`
#' @param preloaded logical; whether or not this step has been completed
#' 
#' @returns a tidy data.frame
annotate.functional.data <- function(df.func = tibble(), preloaded = FALSE) {
    
    # stored on /fast/AG_Forslund/Morgan/projects/covid-opitz/data/proc-data-lg
    if (preloaded) 
        return(vroom::vroom(here('data','proc-data','annotated-gene-table.tsv')))
    else if (file.exists(here('data','proc-data','annotated-gene-table.tsv')))
        return(vroom::vroom(here('data','proc-data','annotated-gene-table.tsv')))
    else {
        map <- vroom::vroom(here('data','raw-data','metagenomic',
                                 'gmgc-filtered-map-kegg.tsv'))
        tmp <- df.func %>%
            left_join(map, by = c('Rarefied'='#query_name'))
        vroom::vroom_write(tmp, 
                           here('data','proc-data',
                                'annotated-gene-table.tsv'))
        return(tmp)
    }
}


#' Reads in and partitions original abundance table and writes it to temp/
#' Reads in chunked tables and joins with the gmgc mapfile
#' 
#' @param raw.gene.cts string; filepath to raw abundance table
#' @param raw.chunked.path string: filepath to chunked abundance tables
#' @param raw.mapfile string; filepath to raw gmgc mapfile
#' 
#' @returns a tidy data.frame
functional.annotation.from.scratch <- function(raw.gene.cts, # raw.chunked.path,
                                               raw.mapfile) {
    
    # make a cluster using half the cores on the VM
    message('+ creating cluster ...')
    worker.chunks <- future::availableCores()*(3/4)
    cluster <- multidplyr::new_cluster(worker.chunks)
    multidplyr::cluster_library(cluster, c('tidyverse','here'))
    message('... done !')
    
    # read raw file and prepare temp directory
    tmpdir <- here('data','proc-data','temp')
    # if (!dir.exists(tmpdir))
    #     dir.create(tmpdir)
    # if (length(dir(tmpdir)) > 0)
    #     unlink(tmpdir, recursive = TRUE)
    # no matter what, quite memory intensive (20GB RAM?)
    # raw.gene.cts <- vroom::vroom(raw.gene.cts)
    
    # write partitioned files - do once and comment out
    # idx <- rep_len(1:worker.chunks, 
    #                length.out = nrow(raw.gene.cts))
    # also very very slow step
    # raw.gene.cts %>%
    #     add_column(group = as_factor(idx), .before = 1) %>%
    #     group_by(group) %>%
    #     group_walk(~ vroom::vroom_write(.x,
    #                                     here(tmpdir,
    #                                          paste0(.y$group, 
    #                                                 '.tsv'))))
    # remove from memory and read in piecewise next step
    # rm(raw.gene.cts)
    
    raw.mapfile <- vroom::vroom(raw.mapfile)
    
    # send each chunk to a worker
    files <- dir(tmpdir, full.names = TRUE)
    multidplyr::cluster_assign_partition(cluster, files = files)
    multidplyr::cluster_send(cluster, 
                             part.tmp <- vroom::vroom(files))
    cleaned.df <- cluster %>%
        multidplyr::party_df('part.tmp', 
                             auto_rm = FALSE) %>%
        left_join(raw.mapfile, 
                  by = c('Rarefied'='#query_name'),
                  # send mapfile to workers, slow but necessary for now
                  copy = TRUE) %>%
        collect() %>%
        fix.rtk.bug() %>%
        select(-group)
    
    # unfinished .. in reality, do stepwise from plan
    # after this, clean.functional.map to get eg gene-KEGGko pairs
}


#' Colshift and NA remove for rtk-rarefied abundance tables
#' 
#' @param df data.frame; raw data frame
#' 
#' @returns a tidy data.frame
fix.rtk.bug <- function(df) {
    
    cn <- colnames(df)
    tmp <- df %>%
        select(-S97) %>%
        separate(S198, into = c('tmp1','tmp2'), sep = '\t') %>%
        set_names(cn)
    return(tmp)
}


#' Reduces space of annotated gene table to produce a map of genes found in
#' the samples and a more relevant functional category, e.g. KEGG modules
#' 
#' @param df.func data.frame; annotated gene abundance table -- not needed and 
#'                defaults to an empty tibble once annotated table has been
#'                chunked (will access from raw.path)
#' @param raw.path string; tmpdir for the raw chunked annotated table
#' @param clean.path string; tmpdir for the 'by' filtered chunked annotated table
#' @param by string; functional mapping variable
#' @param preloaded logical; whether or not this step has been completed
#' 
#' @returns a tidy data.frame
clean.functional.map <- function(df.func = tibble(), 
                                 raw.path, clean.path, by,
                                 preloaded = FALSE) {
    
    if (preloaded) 
        return(vroom::vroom(here('data','proc-data',
                                 paste0('gene-', by, '-map.tsv'))))
    else {
        # DON'T FORGET TO MOVE LARGE FILES/DIRS TO /FAST/ WHEN DONE
        # AND SOFTLINK FROM VM DIRS !
        stopifnot(length(by)>0)
        # make a cluster using half the cores on the VM
        message('+ creating cluster ...')
        worker.chunks <- future::availableCores()*(3/4)
        cluster <- multidplyr::new_cluster(worker.chunks)
        multidplyr::cluster_library(cluster, c('tidyverse','here'))
        message('... done !')
        
        # STEP 1: partition raw data for cleanup (makes smaller for next step)
        # if (!dir.exists(raw.path))
        #     dir.create(raw.path)
        # write partitioned files - do once and comment out
        # idx <- rep_len(1:worker.chunks, length.out = nrow(df.func))
        # df.func %>%
        #     add_column(Sgroup = as_factor(idx), .before = 1) %>%
        #     group_by(Sgroup) %>%
        # group_walk(~ vroom::vroom_write(.x,
        #                                 here(raw.path,
        #                                      paste0(.y$Sgroup, '.tsv'))))
        
        # START UNCOMMENTING HERE IF 'BY' PARAMETER CHANGED
        # browser()
        # send each annotated chunk to a worker
        # files <- dir(raw.path, full.names = TRUE)
        # multidplyr::cluster_assign_partition(cluster, files = files)
        # multidplyr::cluster_send(cluster,
        #                          # read in each chunk on respective worker
        #                          # name of variable needed for party_df
        #                          part.tmp <- vroom::vroom(files))
        # # perform filtering step and rename 'by' variable to generic 'temp'
        # cleaned.df <- cluster %>%
        #     # name must match what was saved on each worker^
        #     multidplyr::party_df('part.tmp',
        #                          auto_rm = FALSE) %>%
        #     clean.gene.map.pair(by = by) %>%
        #     # bring back to this worker
        #     collect()
        
        # STEP 2: partition filtered data and restructure
        # send each gene-'by' pair to a worker
        # if (!dir.exists(clean.path))
        #     dir.create(clean.path)
        # # write cleaned, partitioned files - do once and comment out
        # idx <- rep_len(1:worker.chunks, length.out = nrow(cleaned.df))
        # cleaned.df %>%
        #     add_column(Sgroup = as_factor(idx), .before = 1) %>%
        #     group_by(Sgroup) %>%
        #     group_walk(~ vroom::vroom_write(.x,
        #                                 here(clean.path,
        #                                      paste0(.y$Sgroup, '.tsv'))))
        
        # computation to create a clean (long) df of gene-'by' pairs
        # files <- dir(clean.path, full.names = TRUE)
        # multidplyr::cluster_assign_partition(cluster, files = files)
        # # can't unnest with a partitioned data frame, so 
        # # already do that when reading in files on workers
        # multidplyr::cluster_send(cluster, 
        #                          part <- vroom::vroom(files) %>%
        #                              # only needed for KO, not modules
        #                              mutate(temp = str_remove_all(temp, 'ko:')) %>%
        #                              mutate(temp = map(temp, ~ str_split(., ','))) %>%
        #                              unnest_longer(temp) %>%
        #                              unnest_longer(temp))
        # # save off filtered partitioned data frames
        # final.df <- cluster %>%
        #     multidplyr::party_df('part', 
        #                          auto_rm = FALSE) %>%
        #     # bring back to this worker
        #     collect() %>%
        #     rename(!!by := 'temp')
        # vroom::vroom_write(final.df, here('data','proc-data', 
        #                                   paste0('gene-', by, '-map.tsv')))
        return(final.df) }
}

#' Does basic filtering of NAs and generically renames for later processing
#' 
#' @param df.kegg.anno data.frame; annotated gene abundance table
#' @param by string; mapping variable of interest for next step
#' 
#' @returns a tidy data.frame
clean.gene.map.pair <- function(df.kegg.anno, by) {
    
    if (str_detect(all_of(by), 'KEGG')) {
            # browser()
            tmp <- df.kegg.anno %>%
                select(Rarefied, !!rlang::sym(by)) %>%
                filter(!is.na(!!rlang::sym(by))) %>%
                rename(temp = by) 
            }
        else {
            tmp <- df.kegg.anno %>%
                select(Rarefied, !!rlang::sym(by)) %>%
                filter(!is.na(!!rlang::sym(by))) }
        return(tmp)
}


#' Sums mOTUs into genera for higher power in modeling steps
#' 
#' @param df.abun data.frame; annotated gene abundance table
#' @param df.map data.frame; mapping variable of interest for next step
#' @param by string; mapping variable of interest for next step
#' @param preloaded logical; whether or not this step has been completed
#' 
#' @returns a tidy data.frame
bin.gene.abundances <- function(df.abun, df.map, by, preloaded = FALSE) {
    
    # browser() # can't access vroom-read tbls with this...
    if (preloaded) {
        if (by == 'KEGG_Module')
            return(vroom::vroom(here('data','proc-data',
                                     'binned-kegg-modules.tsv')))
        else if (by == 'KEGG_ko')
            return(vroom::vroom(here('data','proc-data',
                                     'binned-kegg-kos.tsv')))
        else browser()
    }
    else {
        df.abun %<>% select(Rarefied, starts_with('S'))
        tmp <- df.map %>%
            left_join(df.abun) %>%
            select(!!rlang::sym(by), everything()) %>%
            group_by(!!rlang::sym(by)) %>%
            summarize_if(., is.double, sum)

        if (by == 'KEGG_Module')
            vroom::vroom_write(tmp, here('data','proc-data',
                                     'binned-kegg-modules.tsv'))
        else if (by == 'KEGG_ko')
            vroom::vroom_write(tmp, here('data','proc-data',
                                         'binned-kegg-kos.tsv'))
        else
            browser()
        return(tmp)
    }
}


#' Sums mOTUs into genera for higher power in modeling steps
#' 
#' @param df.otu data.frame; motu abundance table
#' @param rel.ab logical; use unmapped to calc relative abundances of feat
#' 
#' @returns a tidy data.frame
bin.motus <- function(df.otu, df.clin, rel.ab = TRUE, time = 'Days') {
 
    # browser()
    times <- select(df.otu, UniqueID, rlang::sym(time))
    # tmp <- df.otu %>%
    #     # select(-c(Site,PatientID,Visit, rlang::sym(time))) %>%
    #     select(-Site, -GenID, -rlang::sym(time)) %>%
    #     tidyr::gather('mOTU','count',-UniqueID) %>%
    #     spread('UniqueID','count') %>%
    #     mutate(taxonomy = str_split(mOTU, ' ')) %>%
    #     mutate(bin = map_chr(taxonomy, ~ head(., 1))) %>%
    #     select(mOTU, taxonomy, bin, everything()) %>%
    #     mutate(bin = str_remove_all(bin, fixed('['))) %>%
    #     mutate(bin = str_remove_all(bin, fixed(']'))) %>%
    #     # remove odd species (~3)
    #     filter(bin != 'bacterium') %>%
    #     filter(bin != 'bacteria') %>%
    #     filter(bin != 'uncultured') %>%
    #     group_by(bin) %>% 
    #     summarize_if(is.double, sum) %>%
    #     distinct() %>%
    #     column_to_rownames('bin')
    
    tmp <- df.otu %>%
        select(-Site, -GenID, -rlang::sym(time), -UniqueID) %>%
        t() %>%
        as.data.frame() %>% rownames_to_column('mOTU') %>%
        mutate(taxonomy = str_split(mOTU, ' ')) %>%
        mutate(bin = map_chr(taxonomy, ~ head(., 1))) %>%
        select(mOTU, taxonomy, bin, everything()) %>%
        mutate(bin = str_remove_all(bin, fixed('['))) %>%
        mutate(bin = str_remove_all(bin, fixed(']'))) %>%
        # remove odd species (~3)
        filter(bin != 'bacterium') %>%
        filter(bin != 'bacteria') %>%
        filter(bin != 'uncultured') %>%
        group_by(bin) %>% 
        summarize_if(is.double, sum) %>%
        distinct() %>%
        column_to_rownames('bin') %>%
        magrittr::set_colnames(df.otu$UniqueID)
    
    if (rel.ab) {
        tmp <- tmp %>%
            # relative abundances
            remove.zero.cols() %>%
            as.matrix() %>%
            prop.table(2) %>%
            as.data.frame() }
    
    df.out <- tmp %>%
        # remove unmapped after rel.ab calculation
        rownames_to_column('mOTU') %>%
        filter(!str_detect(mOTU, 'Unmapped_')) %>%
        column_to_rownames('mOTU') %>%
        mutate(variance = matrixStats::rowVars(as.matrix(.))) %>%
        filter(variance != 0) %>%
        select(-variance) %>%
        rownames_to_column('mOTU') %>%
        as_tibble()
    
    prev.filtered <- df.out %>%
        tidyr::gather('UniqueID','rel.ab', -mOTU) %>%
        spread('mOTU','rel.ab') %>%
        summarize_if(is.double, ~ calc.single.prevalence(.)) %>%
        tidyr::gather('mOTU','prevalence') %>%
        filter(prevalence >= 0.1) %>% 
        pull(mOTU)
    
    removed <- df.out$mOTU[!(df.out$mOTU %in% prev.filtered)]
    
    df.out %>%
        filter(mOTU %in% prev.filtered) %>%
        tidyr::gather('UniqueID','rel.ab', -mOTU) %>%
        spread('mOTU','rel.ab') %>%
        left_join(times) %>%
        select(UniqueID, rlang::sym(time), everything())
}


bin.motus.simp <- function(df.otu) {
    
    # simply sum up existing relative abundances -- should after all be additive
    # now total rel. ab. of OSCI-assoc. feat is a sample/patient characteristic
    # ranges in stool from 0.01 to 0.9
    tmp <- df.otu %>%
        gather('mOTU','ra',-UniqueID) %>%
        spread('UniqueID','ra') %>%
        mutate(taxonomy = str_split(mOTU, ' ')) %>%
        mutate(bin = map_chr(taxonomy, ~ head(., 1))) %>%
        select(mOTU, taxonomy, bin, everything()) %>%
        mutate(bin = str_remove_all(bin, fixed('['))) %>%
        mutate(bin = str_remove_all(bin, fixed(']'))) %>%
        # remove odd species (~3)
        filter(bin != 'bacterium') %>%
        filter(bin != 'bacteria') %>%
        filter(bin != 'uncultured') %>%
        group_by(bin) %>% 
        summarize_if(is.double, sum) %>%
        distinct() %>%
        gather('UniqueID','ra',-bin) %>%
        spread('bin','ra') %>%
        arrange(UniqueID)
    
    # proportion of each sample which is osci-associated microbiota features
    # interesting, but not useful to re-normalize
    # props <- tmp %>% 
    #     column_to_rownames('UniqueID') %>% 
    #     rowSums() %>%
    #     enframe('UniqueID','RA.prop') %>%
    #     arrange(UniqueID)
    
    # mean abundance and variance thresholds
    keep <- tmp %>% 
        summarize_if(is.double, mean) %>%
        gather('genus','mean') %>% 
        mutate(variance = matrixStats::colVars(as.matrix(tmp[,-1]))) %>%
        filter(mean >= 0.001 & variance > 1e-6) %>% 
        pull(genus)
    
    select(tmp, UniqueID, all_of(keep))
}

bin.metabs.simp <- function(df.metab, df.annot) {
    
    df.metab %>%
        gather('R_Compound','val',-UniqueID) %>%
        mutate(Sample = str_split(R_Compound,'_') %>%
                   map_chr(~ head(., 1))) %>%
        mutate(R_Compound = str_remove_all(R_Compound,'P_')) %>%
        mutate(R_Compound = str_remove_all(R_Compound,'U_')) %>%
        left_join(select(df.annot, R_Compound, Compound, 
                         Measurement, ChemClass)) %>%
        # ratio of pro to cit,, tryptophan_13C
        filter(!is.na(ChemClass)) %>%
        # check and remove duplicate measurements
        # add_count(Compound) %>% 
        # filter(n==192)  %>% 
        # add_count(ChemClass) %>% 
        # arrange(desc(nn))%>% 
        # distinct(across(c(
        #     'R_Compound','Compound','ChemClass'))) %>% 
        # arrange(Compound)
        filter(!R_Compound %in% c('Gly','X3.IAA','X3.IPA',
                                  'Kynurenine','Trp')) %>%
        nest(data = -c(UniqueID, ChemClass)) %>%
        mutate(data = map(data, ~ pull(., val) %>%
                              sum(na.rm = TRUE))) %>%
        unnest(data) %>%
        mutate(ChemClass = str_replace_all(ChemClass,' ','_')) %>%
        spread('ChemClass','data')
    
    # R_Compound              ChemClass          
    # 1 Gly                     Amino acid         
    # 2 Glycine..3TMS.          Amino acid         
    # 3 indole_3_propionic_acid Indole derivatives 
    # 4 indole_acetic_acid      Indole derivatives 
    # 5 Trp                     Amino acid         
    # 6 tryptophan              Amino acid         
    # 7 X3.IAA                  Indole derivatives 
    # 8 X3.IPA                  Indole derivatives 
    # 9 Ornithine..3TMS.        Urea cycle         
    # 10 Ornithine..4TMS.        Urea cycle         
    # 11 kynurenine              Aromatic amino acid
    # 12 Kynurenine              AA related    
    
}

chem.class.stats <- function(df.feat, df.meta, group,
                                  stat = 'cliffs') {
    
    meta <- df.meta %>%
        as.data.frame() %>%
        rownames_to_column('GenID') %>%
        select(GenID, OSCI_Score_Samp) %>%
        mutate(OSCI = case_when(
            OSCI_Score_Samp==0~'Control',
            OSCI_Score_Samp<=4~'Mild',
            OSCI_Score_Samp>4~'Severe')) %>%
        fastDummies::dummy_cols('OSCI') %>%
        select(-OSCI_Score_Samp)
    
    switch(group,
           'mild' = {
               meta.sub <- filter(meta, OSCI != 'Severe') },
           'severe' = {
               meta.sub <- filter(meta, OSCI != 'Mild') })
    
    ctrl.idx <- which(meta.sub$OSCI=='Control')
    case.idx <- which(meta.sub$OSCI!='Control')
    
    if (stat == 'cliffs') {
        df.feat[meta.sub$GenID,] %>%
            map_dfc(function(col) {
                effsize::cliff.delta(
                    col[case.idx],
                    col[ctrl.idx])$estimate }) %>%
            gather('feature','cliffs')
    } else if (stat == 'wilcox') {
        df.feat[meta.sub$GenID,] %>%
            map_dfc(function(col) {
                wilcox.test(
                    col[case.idx],
                    col[ctrl.idx])$p.value }) %>%
            gather('feature','p.val') %>%
            mutate(p.val = p.adjust(p.val, method = 'fdr'))
    } else stop('do not recognize statistical metric')
}

bin.metab.results <- function(df.res, df.metab, by='ChemClass') {
    
    # prev params were df.res, df.metab, by='ChemClass'
    tmp <- df.res %>%
        # remove once updated annotation file available
        filter(!is.na(!!rlang::sym(by))) %>%
        mutate(Site = case_when(
            str_detect(feature, 'P_') ~ 'P',
            str_detect(feature, 'U_') ~ 'U',
            TRUE ~ 'NA')) %>%
        mutate(ChemClass = str_replace_all(ChemClass,' ','_')) %>%
        nest(data = -c(Site, {{by}})) %>%
        mutate(data = pmap(list(data, ChemClass, Site),
                           function(.x, .y, .z) {
            df.metab %>%
                select(UniqueID, all_of(.x$feature)) %>%
                gather('compound','val',-UniqueID) %>%
                group_by(UniqueID) %>%
                # sum OSCI-associated metabolites within this chemclass/pathway
                summarize(val = sum(val)) %>%
                magrittr::set_colnames(c('UniqueID', paste0(.z,'_',.y))) }))

    joined <- tmp$data %>%
        reduce(full_join)
}

#' Transpose to row-samples and create UniqueID for later metadeconfoundR runs
#' 
#' @param df.abun data.frame; annotated gene abundance table
#' @param clinical.feat data.frame; uniqueIDs and days
#' @param by string; functional mapping variable
#' @param sample.map data.frame; SXX sample nomenclature
#' @param short.names data.frame; shortIDs, long PatientIDs
#' 
#' @returns a tidy data.frame
structure.kegg.table <- function(df.abun, clinical.feat, by,
                                 sample.map, short.names) {

    # browser()
    sIDs <- sample.map %>%
        left_join(short.names) %>%
        mutate(SampleID = map_chr(SampleID, ~ paste0('S', .)))
    tmp1 <- df.abun %>%
        tidyr::gather('SampleID','abundance', -!!rlang::sym(by)) %>%
        left_join(sIDs) %>%
        select(colnames(sIDs), everything()) %>%
        mutate(UniqueID = paste0(ShortID,'-',Treatment,'-',Visit)) %>%
        spread(by, 'abundance') %>%
        left_join(select(clinical.feat, UniqueID, Days)) %>%
        select(SampleID, # PatientID, 
               Treatment, UniqueID, Days, 
               starts_with('M'), starts_with('K')) %>%
        rename(Site = 'Treatment') %>%
        distinct() %>%
        discretize.days(keep.numeric = TRUE)
    tmp2 <- tmp1 %>%
        filter(str_detect(UniqueID, 'P11-SDNA-V2|P17-SDNA-V4')) %>% 
        arrange(UniqueID) %>%
        mutate(rows = row_number()) %>%
        mutate(UniqueID = ifelse(rows%%2==0, 
                                 paste0(UniqueID,'.2'), UniqueID)) %>%
        select(-rows) 
    tmp1 %>%
        filter(!str_detect(UniqueID, 'P11-SDNA-V2|P17-SDNA-V4')) %>%
        bind_rows(tmp2)
}


ko.omixer.prep <- function(df.abun) {
    
    tmp <- df.abun %>%
        mutate(UniqueID = str_remove_all(UniqueID, '-SDNA')) %>%
        select(-SampleID, -Site, -Days) %>%
        tidyr::gather('KO','ab', -UniqueID) %>%
        spread('UniqueID','ab')
    fn <- here('data','raw-data','metagenomic','gmms-omixer')
    if (!file.exists(fn))
        vroom::vroom_write(tmp, here(fn, 'ko-table-formatted.tsv'))
    return(tmp)
}


load.omixer.table <- function(mods.dir, wide = TRUE) {

    # browser()
    samples <- mods.dir %>%
        dir() %>%
        map_chr(~ str_remove_all(., '.modules')) %>%
        set_names(vector('list', length(.)), .) # pre-allocate named list
    
    modules <- samples %>%
        imap_dfr(function(.x, .y) {
            read_tsv(here(mods.dir, 
                          paste0(.y, '.modules'))) }, 
            .id = 'SampleID') 
    
    if (wide) {
        modules <- modules %>%
            select(-Coverage) %>%
            spread('Module','Value') %>%
            map_df(~ replace_na(., 0)) }
    
    return(modules)
    
}


structure.gmm.table <- function(df.abun, clinical.feat) {
    
    # browser()
    tmp <- clinical.feat %>%
        filter(Site == 'SDNA') %>%
        select(UniqueID, Visit, Days) %>%
        mutate(UniqueID = str_split(UniqueID, '-') %>% 
                   map_chr(~ head(., 1))) %>%
        mutate(UniqueID = paste0(UniqueID, '-', Visit)) %>%
        select(-Visit)
    df.abun %>%
        rename(UniqueID = 'SampleID') %>%
        left_join(tmp, by = 'UniqueID') %>%
        select(UniqueID, Days, everything()) %>%
        distinct()
}


#' Prevalence and variance filters
#' 
#' @param df.otu data.frame; full motu abundance table
#' @param site string; SDNA or OP
#' 
#' @returns a tidy data.frame
clean.motus <- function(df.otu, site, rel.ab = FALSE,
                        h.cols = 5) {

    # browser()
    df.otu %<>% filter(Site == site)
    tmp <- df.otu %>%
        select(-c(1:h.cols)) %>%
        t() %>%
        as.data.frame() %>% 
        rownames_to_column('motu') %>%
        mutate(libsize = rowSums(.[,2:ncol(.)])) %>%
        filter(libsize != 0) %>%
        select(-libsize) %>%
        column_to_rownames('motu')
    
    if (rel.ab) {
        tmp <- tmp %>%
            as.matrix() %>%
            prop.table(2) %>%
            as.data.frame() %>%
            # need to keep to do accurate rel.abs
            rownames_to_column('mOTU') %>%
            filter(!str_detect(mOTU, 'Unmapped_')) %>%
            column_to_rownames('mOTU') }
    
    df.keep <- tmp %>%
        mutate(variance = matrixStats::rowVars(as.matrix(.)[,2:ncol(.)])) %>%
        filter(variance != 0) %>%
        select(-variance) %>%
        t() %>%
        as.data.frame()
    df.otu %>%
        select(c(1:5)) %>%
        as.data.frame() %>%
        bind_cols(df.keep) %>%
        as_tibble() 
}

near.zero.var.filter <- function(df.otu, n.header.cols = 5) {
    
    n.zeros <- df.otu %>%
        select((n.header.cols+1):ncol(df.otu)) %>%
        summarize_all(function(v) { length(v[v==0])/nrow(df.otu) }) %>%
        gather('otu','pct.zeros') %>%
        arrange(desc(pct.zeros)) %>%
        filter(pct.zeros < 0.9) %>%
        pull(otu)
    df.otu %>%
        select(1:n.header.cols, all_of(n.zeros))
}

log.TSS.transform <- function(df.count, pc = 1) {
    labels <- select(df.count, 1:2)
    tss <- df.count %>%
        select(-c(1:2)) %>%
        as.matrix() %>%
        prop.table(1)
    tss %>%
        as.data.frame() %>%
        summarize_all(mean) %>%
        gather('genus','mean.ab') %>%
        pull(mean.ab) %>%
        summary()
    lt <- vapply(colnames(tmp), function(c) { 
        t <-tmp[,c]; 
        t[t==0] <- 1; 
        return(t) }, FUN.VALUE = numeric(nrow(tmp)))
}


#' Attaches correct time variable for downstream analysis and add unique IDs to
#' ensure alignment between clinical and OTU data
#' Handles the duplicate samples; renames and maps correct clinical data to both
#' 
#' @param df.otu 
#' @param df.clin
#' @param time.type string; either 'continuous' or 'discrete'
#' @param seq.type string; either '16S' or 'WGS'
#' @param short.ids data.frame; map of long and short PatientIDs
#' @param df.out string; desired features to output, either 'OTU' or 'clinical' currently
#' @param feat.id string; fed to `contains()` select helper; common identifier among features
#' 
#' @returns a tidy data.frame
combine.features <- function(df.otu, df.clin, time.type, seq.type,
                             short.ids, df.out = 'OTU', feat.id = 'OTU') {
    
    # browser()
    # C19-CB-120 has 2x SDNA on V4 (df.clin will repeat)
    # C19-CB-87 has 2x SDNA on V2 (df.clin will repeat)
    switch(time.type,
           'continuous' = {
               dup.times <- df.clin %>%
                   select(1:5) %>%
                   dplyr::filter(DataType == seq.type) %>%
                   select(-DataType) %>%
                   mutate(tmp = paste0(Site, PatientID, Visit)) %>%
                   filter(tmp == 'SDNAC19-CB-120V4' | tmp == 'SDNAC19-CB-87V2') %>%
                   select(-tmp) %>%
                   mutate(rows = row_number()) %>%
                   mutate(Visit = ifelse(rows%%2==0, paste0(Visit,'.2'), Visit)) %>%
                   select(-rows)
               times <- df.clin %>%
                   select(1:5) %>%
                   dplyr::filter(DataType == seq.type) %>%
                   select(-DataType) %>%
                   mutate(tmp = paste0(Site, PatientID, Visit)) %>%
                   filter(tmp != 'SDNAC19-CB-120V4') %>%
                   filter(tmp != 'SDNAC19-CB-87V2') %>%
                   select(-tmp) },
           'discrete' = {
               dup.times <- df.clin %>%
                   select(1:4,6) %>%
                   dplyr::filter(DataType == seq.type) %>%
                   select(-DataType) %>%
                   mutate(tmp = paste0(Site, PatientID, Visit)) %>%
                   filter(tmp == 'SDNAC19-CB-120V4' | tmp == 'SDNAC19-CB-87V2') %>%
                   select(-tmp) %>%
                   mutate(rows = row_number()) %>%
                   mutate(Visit = ifelse(rows%%2==0, paste0(Visit,'.2'), Visit)) %>%
                   select(-rows) %>%
                   rename(Time = 'Days')
               times <- df.clin %>%
                   select(1:4,6) %>%
                   dplyr::filter(DataType == seq.type) %>%
                   select(-DataType) %>%
                   mutate(tmp = paste0(Site, PatientID, Visit)) %>%
                   filter(tmp != 'SDNAC19-CB-120V4') %>%
                   filter(tmp != 'SDNAC19-CB-87V2') %>%
                   select(-tmp) %>%
                   rename(Time = 'Days') })
    
    
    switch(df.out,
           'OTU' = {
               # this isn't a problem for WGS which has correct, unique SampleIDs...
               dup.otus <- df.otu %>%
                   mutate(tmp = paste0(Site, PatientID, Visit)) %>%
                   filter(tmp == 'SDNAC19-CB-120V4' | tmp == 'SDNAC19-CB-87V2') %>%
                   select(-c(Date, tmp)) %>%
                   mutate(rows = row_number()) %>%
                   mutate(Visit = ifelse(rows%%2==0, paste0(Visit,'.2'), Visit)) %>%
                   select(-rows) %>%
                   arrange(PatientID, Visit) 
               otus <- times %>%
                   left_join(select(df.otu, -Date)) %>%
                   select(Site, PatientID, Visit, Time, contains(feat.id)) %>%
                   arrange(Site, PatientID, Visit) 
               tmp <- dup.times %>%
                   left_join(dup.otus) %>%
                   distinct(across(c(Site, PatientID, Visit, 6)), 
                            .keep_all = TRUE) %>%
                   select(Site, PatientID, Visit, Time, contains(feat.id)) %>%
                   bind_rows(otus) %>%
                   arrange(Site, PatientID, Visit) },
           'metabolites' = { },
           'metab.otu' = { },
           'clinical' = { # match clinical data to otu output
               dup.clin <- df.clin %>%
                   mutate(tmp = paste0(Site, PatientID, Visit)) %>%
                   filter(tmp == 'SDNAC19-CB-120V4' | tmp == 'SDNAC19-CB-87V2') %>%
                   select(-c(tmp)) %>%
                   mutate(rows = row_number()) %>%
                   mutate(Visit = ifelse(rows%%2==0, paste0(Visit,'.2'), Visit)) %>%
                   select(-rows) %>%
                   arrange(PatientID, Visit) 
               clin <- df.clin %>%
                   # select(-Days) %>%
                   select(Site, PatientID, Visit, Time, DataType, everything())
               tmp <- dup.times %>%
                   left_join(dup.clin) %>%
                   distinct(across(c(Site, PatientID, Visit, 10)), 
                            .keep_all = TRUE) %>%
                   # select(-Days) %>%
                   bind_rows(clin) %>%
                   arrange(Site, PatientID, Visit)
           })
    
    tmp %>%
        left_join(short.ids) %>%
        mutate(UniqueID = paste(ShortID, Site, Visit, sep = '-')) %>%
        select(UniqueID, everything(), -ShortID)
    
}


#' Creates UniqueID and renames columns from Mathias's raw data
#' 
#' @param df.met data.frame; raw metabolites (long format)
#' @param short.ids data.frame; of the P13 / C19-CB-105 pairs
#' 
#' @returns a tidy data.frame
clean.metabolites <- function(df.met, short.ids) {
    
    df.met %>%
        rename(PatientID = 'Patient.ID') %>%
        rename(Site = 'Matrix') %>%
        mutate_at(vars(PatientID), ~ str_replace_all(., 'K','K-')) %>%
        left_join(short.ids) %>%
        mutate(UniqueID = paste(ShortID, Site, Visit, sep = '-')) %>%
        select(UniqueID, everything(), -ShortID)
}


#' Filters metabolites by their originating body fluid and updates UniqueID
#' to be of the PXX-VX type, without site intervening
#' 
#' @param df.met data.frame; metabolites in long format
#' @param df.days data.frame; UniqueID PatientID (discretized) Days
#' @param site string; data origin: P for plasma or U for urine
#' 
#' @returns a tidy data.frame
filter.metabolites <- function(df.met, df.days, site) {
    
    # browser()
    df.met %>%
        filter(Site == site) %>%
        select(-c(Site, Visit)) %>%
        mutate(UniqueID = str_remove_all(UniqueID, 
                                         paste0('-', site))) %>%
        mutate(PatientID = str_split(UniqueID, '-') %>% 
                   map_chr(~ head(., 1))) %>%
        # mutate_at(vars(UniqueID), ~ ifelse(. == 'P17-V2', 'P17-V1', .)) %>%
        # mutate_at(vars(UniqueID), ~ ifelse(. == 'P27-V1', 'P27-V2', .)) %>%
        left_join(df.days) %>%
        mutate(Days = replace_na(Days, 'Early')) 
}


#' Removes weird characters and leading numbers, spreads/makes wide
#' 
#' @param df.met data.frame; metabolites in long format
#' @param df.days data.frame; UniqueID PatientID (discretized) Days
#' @param site string; data origin: P for plasma or U for urine
#' 
#' @returns a tidy data.frame
r.friendly.metabolite.names <- function(df.met, site, incl.site = FALSE) {
    
    tmp <- df.met %>%
        # previous manual effort, for recordkeeping ... 
        # mutate(Compound = str_replace_all(Compound, ' ', '_')) %>%
        # mutate(Compound = str_replace_all(Compound, '-', '_')) %>%
        # mutate(Compound = str_replace_all(Compound, ':', '_')) %>%
        # mutate(Compound = str_replace_all(Compound, '/', '_')) %>%
        # mutate(Compound = str_replace_all(Compound, '^[:digit:]', 'x')) %>%
        mutate(Compound = make.names(Compound)) %>%
        select(-c(PatientID, Method, Batch, Days))
    
    if (incl.site) {
        tmp %>%
            mutate(Compound = paste0(site,'_',Compound)) %>%
            spread('Compound','Value') }
    else spread(tmp, 'Compound', 'Value')
}


# TODO! awaiting fix of CCM names
clean.metabolite.map <- function(df.met, df.map) {
    
    # browser()
    tmp <- df.met %>%
        mutate(R_Compound = make.names(Compound))
    
    # Data_Compound = how they're called in mathias's input file
    data.names <- select(tmp, Compound, R_Compound) %>%
        distinct() %>%
        rename(Data_Compound = 'Compound')
    
    map <- df.map %>%
        select(1:6) %>%
        set_names(c('Data_Compound','Compound','ChemClass',
                  'Measurement','Method','Pathway')) %>%
        mutate(Method = str_replace_all(Method,'tryptophan','TRYP')) %>%
        mutate(Method = str_replace_all(Method,'Quant500','Q500')) %>%
        mutate(Method = ifelse(Method == 'Q500', 
                               paste0(Method, Measurement),
                               Method)) %>%
        mutate(Pathway = str_replace_all(Pathway, 
                                         'Brach chain amino acids',
                                         'BCAA')) %>%
        # 1st map col for CCM does not match data -- make them match, here
        mutate(Data_Compound = ifelse(Method=='CCM',
                                      make.names(Data_Compound),
                                      Data_Compound))
    data.names %>%
        left_join(map, by='Data_Compound')
    
    # previous efforts -- 06.22 fixed input xlsx with Ulrike & Ela
    # ccm.metabs <- map %>%
    #     filter(Method == 'CCM') 
    
    # out <- data.names %>% 
    #     left_join(filter(map, Method != 'CCM'), 
    #               by = c('Data_Compound'='Abb_Compound')) %>%
    #     left_join(ccm.metabs, 
    #               by = c('Data_Compound'='CCM_Data_Compound')) %>%
    #     mutate(Method = coalesce(Method.x, Method.y)) %>%
    #     mutate(Measurement = coalesce(Measurement.x, Measurement.y)) %>%
    #     mutate(ChemClass = coalesce(ChemClass.x, ChemClass.y)) %>%
    #     mutate(Pathway = coalesce(Pathway.x, Pathway.y)) %>%
    #     #mutate(Compound = coalesce(Compound.x, Abb_Compound)) %>%
    #     mutate(Compound = ifelse(Method == 'TRYP', 
    #                              coalesce(Data_Compound, 
    #                                       Compound.x, 
    #                                       Abb_Compound),
    #                              coalesce(Compound.x, Abb_Compound) )) %>%
    #     mutate(R_Compound = coalesce(R_Compound, CCM_Data_Compound)) %>%
    #     select(Data_Compound, R_Compound, Compound, ChemClass, Pathway,
    #            Measurement, Method) %>%
    #     mutate(Compound = ifelse(is.na(Compound), Data_Compound, Compound))
    
}


#' Combines cytokine data, cleans, and assigns UniqueID (not yet)
#' 
#' @param plasma.ck data.frame; 7 plasma cytokine values (pcr?)
#' @param airway.ck data.frame; 2 airway cytokine values
#' @param short.ids data.frame; of the P13 / C19-CB-105 pairs
#' 
#' @returns a tidy data.frame
clean.cytokines <- function(plasma.ck, airway.ck, short.ids) {

    # browser()
    # 92 samples
    tmp.p <- plasma.ck %>%
        rename(Visit = 'Visiten ID') %>%
        # remove rows where all values are NA
        filter(if_any(c(3:9), ~ !is.na(.))) %>%
        tidyr::gather('Cytokine', 'Value', -c(PatientID, Visit))
    
    # 94 samples -- down to 82 with NAs removed
    tmp.a <- airway.ck %>% 
        # remove if both values NA
        filter(if_any(c(`IFNß-GAPDH`, `IFNL2-GAPDH`), ~ !is.na(.))) %>%
        tidyr::gather('Cytokine', 'Value', -c(PatientID, Visit)) 
    
    # more than 2 differences -- previous for wide formatted table
    # compare <- waldo::compare(plasma.ck$UniqueID, airway.ck$UniqueID)
    # plasma.ck only: C19-CB-120V2, C19-CB-135V14, C19-CB-164V7, 
        # C19-CB-197V2, C19-CB-200V2, C19-CB-86V2, K-07V7, K-09V7, K-10V7
    # airway.ck only: C19-CB-40V1, C19-CB-114V4, C19-CB-115V1, 19-CB-75V2
    # print(compare, n = Inf)
    
    # Site indicates 16S site data will be correlated with, **not** site of origin
    # plasma (systemic) cytokines correlated with SDNA, airway (OP) with OP, etc
    tmp.a %>%
        add_column(Site = 'OP', .before = 1) %>%
        bind_rows(add_column(tmp.p, Site = 'SDNA', .before = 1)) %>%
        left_join(short.ids) %>%
        mutate(UniqueID = paste(ShortID, Site, Visit, sep = '-')) %>%
        select(UniqueID, everything(), -ShortID) %>%
        arrange(UniqueID)
}


#' Filters and returns specific cytokine data; prep for statistical analysis
#' 
#' @param df.ck data.frame; all cytokine data with PXX-VXX UniqueIDs
#' @param df.days data.frame; data-independent UniqueIDs with discretized days
#' @param site string; either SDNA (ELISA) or OP (qRT-PCR)
#' 
#' @returns a tidy data.frame
filter.cytokines <- function(df.ck, df.days, site) {
    
    # browser()
    tmp <- df.ck %>%
        filter(Site == site) %>%
        select(-Site) %>%
        left_join(df.days)
    
    if (site == 'SDNA') {
        tmp <- tmp %>%
            select(UniqueID, everything(), -contains('-')) }
    else if (site == 'OP') {
        tmp <- tmp %>%
            select(UniqueID, PatientID, Days, contains('-')) %>%
            set_names(c('UniqueID','PatientID','Days',
                        'IFNL2_GAPDH','IFNB_GAPDH')) }
    else { message('+ cytokine-associated body site not recognized!') }
    
    select(tmp, UniqueID, PatientID, Days, everything()) 
}


#' Filters and cleans monocyte frequencies and ISG scores
#' Replaces European ##,## notation with ##.## (native numeric)
#' 
#' @param df.raw data.frame; raw data from .xlsx
#' 
#' @returns a tidy data.frame
clean.freqs.isg <- function(df.raw) {
    
    # browser()
    df.raw %>%
        dplyr::filter(!is.na(PatientID)) %>%
        magrittr::set_names(c('PatientID','ShortID','Visit',
                    'cMono','ncMono','mDC','pDC','pNK','NK','MK','ISG')) %>%
        mutate_at(c(4:11), ~ str_replace_all(., ',', '.')) %>%
        mutate_at(c(4:11), ~ as.numeric(.)) %>%
        mutate(UniqueID = paste0(ShortID, '-', Visit)) %>%
        select(-PatientID, -ShortID, -Visit) %>%
        # for some reason, controls named with C in sc data
        mutate(UniqueID = str_replace_all(UniqueID, 'C','K')) 
}


#' CLeans up and harmonizes with other single cell files for later combination
#' 
#' @param df.raw data.frame; raw CYP/AHR expression data from .xlsx
#' @param df.pbmc data.frame; cleaned up ISG scores and monocyte frequencies
#' 
#' @returns a tidy data.frame
clean.cyps <- function(df.raw, df.pbmc) {
    
    # browser()
    # get measurement times/visits from other sc data
    pv.map <- df.pbmc %>%
        tidyr::gather('Celltype','value',-UniqueID) %>%
        select(-value) %>%
        distinct() %>%
        separate(UniqueID, into = c('PatientID','Visit'), sep = '-') %>%
        distinct(PatientID, Visit)
    # end up with 25 patients across 5 cell types
    tmp <- df.raw %>%
        rename(PatientID = 'Sample_ID') %>%
        mutate(PatientID = str_replace_all(PatientID,'C','K')) %>%
        mutate(Celltype = str_replace_all(Celltype,
                                          'Non-classic Mono.',
                                          'ncMono')) %>%
        mutate(Celltype = str_replace_all(Celltype,
                                          'cMono.',
                                          'cMono')) %>%
        nest(data = -Celltype) %>%
        mutate(data = map(data, ~ left_join(., pv.map, by = 'PatientID') %>%
                              mutate(UniqueID = paste0(PatientID,'-',Visit)) %>%
                              select(-PatientID, -Visit) %>%
                              tidyr::gather('var','val',-UniqueID))) %>%
        mutate(data = map2(data, Celltype, ~ mutate(.x,
                                                    var = paste0(var,'_',.y)) %>%
                               spread('var','val') %>%
                               arrange(UniqueID))) %>%
        pull(data) %>%
        reduce(~ left_join(.x, .y, by = 'UniqueID')) %>%
        mutate_at(c(2:ncol(.)), ~ as.numeric(.))
        
}

################################################
############## MODELING FUNCTIONS ##############
################################################
create.master.table <- function(df.feat, df.clin, df.meta, metagenomic = TRUE,
                                patients.only = FALSE, site = NULL) {
    if (metagenomic) {
        tmp <- df.feat %>%
            left_join(select(df.clin, 
                             UniqueID, PatientID, contains('OSCI'),
                             VL_OP, Antibiotics)) %>% 
            mutate(Time = Days) %>%
            discretize.days(keep.numeric = FALSE) %>%
            discretize.days(var = 'Time', keep.numeric=TRUE) %>%
            mutate(Time.int = as.numeric(Time)) %>%
            mutate(Abx = factor(Antibiotics, levels = c(0,1))) %>%
            mutate(COVID = ifelse(str_detect(PatientID,'K'), 0, 1)) %>%
            mutate(COVID = factor(COVID, levels = c(0,1))) %>%
            mutate(OSCI_Fct_Samp = factor(OSCI_Score_Samp, levels=c(0,1,2,3,4,5,6,7,8))) %>%
            mutate(OSCI_Class_Samp = fct_collapse(as.character(OSCI_Fct_Samp),
                                             Control = '0',
                                             Mild = c('1','3','4'),
                                             Severe = c('5','6','7'))) %>%
            left_join(select(df.meta,
                             ShortID, DaysHospitalized, Comorbidity, Medication, 
                             Pneumonia, GISymptoms, contains('OSCI')),
                      by = c('PatientID'='ShortID')) %>%
            rename(nMeds = 'Medication', CCI = 'Comorbidity') %>%
            mutate(HAP = factor(Pneumonia, levels=c(0,1))) %>%
            mutate(GI = factor(GISymptoms, levels=c(0,1))) %>%
            select(-Pneumonia, -Antibiotics, -GISymptoms)
        if (patients.only) return(filter(tmp, !str_detect(PatientID,'K')))
        else return(tmp)
    } else {
        df.feat %>%
            left_join(filter(df.clin, Site == site) %>%
                          select(GenID, Days, Antibiotics, contains('OSCI')),
                      by = c('UniqueID'='GenID')) %>%
            mutate(PatientID = str_split(UniqueID,'-') %>%
                       map_chr(~ head(.,1))) %>%
            mutate(Time = Days) %>%
            discretize.days(keep.numeric = FALSE) %>%
            discretize.days(var = 'Time', keep.numeric=TRUE) %>%
            mutate(Time.int = as.numeric(Time)) %>%
            mutate(Abx = factor(Antibiotics, levels = c(0,1))) %>%
            mutate(COVID = ifelse(str_detect(PatientID,'K'), 0, 1)) %>%
            mutate(COVID = factor(COVID, levels = c(0,1))) %>%
            mutate(OSCI_Fct_Samp = factor(OSCI_Score_Samp, 
                                          levels=c(0,1,2,3,4,5,6,7,8))) %>%
            mutate(OSCI_Class_Samp = fct_collapse(as.character(OSCI_Fct_Samp),
                                             Control = '0',
                                             Mild = c('1','3','4'),
                                             Severe = c('5','6','7'))) %>%
            left_join(select(df.meta,
                             ShortID, DaysHospitalized, Comorbidity, Medication, 
                             Pneumonia, GISymptoms, contains('OSCI')),
                      by = c('PatientID'='ShortID')) %>%
            rename(nMeds = 'Medication', CCI = 'Comorbidity') %>%
            mutate(HAP = factor(Pneumonia, levels=c(0,1))) %>%
            mutate(GI = factor(GISymptoms, levels=c(0,1))) %>%
            select(-Pneumonia, -Antibiotics, -GISymptoms) %>%
            # left_join(select(s.cytokines, UniqueID, IP10, IL10, IFNa, IFNg),
            #           by = 'UniqueID') %>%
            filter(!is.na(Days))
    }

}

model.diagnostics <- function(model.list) {
    model.list %>%
        mutate(stats = map(mods, ~ broom.mixed::tidy(.))) %>%
        # mutate(tmp = map(mods, ~ effects::predictorEffects(.))) %>%
        # mutate(pred.diag = map(tmp, function(x) {
        #     plot(x)
        #     recordPlot() })) %>%
        # select(-tmp) %>%
        mutate(model.diag = map(mods, ~ check_model(.))) %>%
        # problem if this is < 0.05 (residuals aren't normally distributed)
        mutate(resid.diag = map_dbl(mods, ~ as.numeric(check_normality(.)))) %>%
        mutate(collin = map(mods, ~ as.data.frame(check_collinearity(.)))) %>%
        mutate(singular = map_lgl(mods, ~ check_singularity(.)))
}

clean.model.coefs <- function(df.stats) {
    df.stats %>% 
        select(feat, stats) %>% 
        unnest(stats) %>% 
        filter(effect=='fixed' & term!='(Intercept)') %>% 
        select(-group,-statistic,-df) %>% 
        mutate(err.pct = std.error/abs(estimate)*100) %>% 
        arrange(p.value, err.pct)
}

basic.volcano.plots <- function(df.clean.stats, df.prev, df.abun,
                                metagenomic = TRUE, labels = NULL) {
    
    if (metagenomic) {
        tmp <- df.clean.stats %>% 
            left_join(df.prev, by=c('feat'='feature')) %>%
            left_join(df.abun, by=c('feat'='feature')) %>%
            ggplot(aes(x=estimate,y=-log10(p.value))) + 
            geom_point(aes(alpha=prevalence, size=mean.ra))
    } else {
        tmp <- df.clean.stats %>%
            ggplot(aes(x=estimate,y=-log10(p.value))) + 
            geom_point()
        # later join some metabolite summary metrics
    }
    tmp <- tmp +
        geom_hline(yintercept = -log10(0.05)) + 
        facet_wrap(~ term) + 
        ggembl::theme_publication() 
    
    if (is.null(labels)) {
        tmp + 
            ggrepel::geom_label_repel(
                dplyr::filter(df.clean.stats, p.value<0.05), size = 7/.pt, 
                mapping = aes(label = feat), 
                min.segment.length = unit(0, 'lines'), force = 5)
    } else {
        tmp +
            ggrepel::geom_label_repel(
                dplyr::filter(df.clean.stats, feat %in% labels), 
                size = 7/.pt, 
                mapping = aes(label = feat), 
                min.segment.length = unit(0, 'lines'), force = 5)
    }
}

################################################
##### LIBRARY SIZE, ALPHA & BETA DIVERSITY #####
################################################


#' Sums up reads per sample when cols are features and rows are samples
#' 
#' @param df.ab data.frame; raw abundance table for a *single* body site
#' @param feat.id string; either 'OTU' or 'N' to look at bacterial or human reads
#' 
#' @returns a tidy data.frame with 2 columns: UniqueID and LibSize
calc.lib.size <- function(df.ab, feat.id = 'OTU') {
    
    if (feat.id == 'OTU') {
        df.ab %>%
            select(UniqueID, contains(feat.id)) %>%
            mutate(LibSize = rowSums(.[,2:ncol(.)])) %>%
            select(UniqueID, LibSize)

    } else if (feat.id == 'N') {
        df.ab %>%
            select(UniqueID, starts_with(feat.id)) %>%
            mutate(LibSize = rowSums(.[,2:ncol(.)])) %>%
            select(UniqueID, LibSize)
            #select(UniqueID, `-1`) 
        
    } else stop('feature identifier not recognized, please check abundance table.') 
}


#' Wrapper for `vegan::diversity(index = 'shannon')`
#' 
#' @param df.ab data.frame; raw abundance table for a *single* body site
#' @param index string; either 'shannon' 'simpson' or 'invsimpson'
#' 
#' @returns a tidy data.frame with 2 columns: UniqueID and index 
calc.diversity <- function(df.ab, index = 'shannon') {
    
    df.ab %>%
        arrange(UniqueID) %>%
        select(UniqueID, c(5:ncol(.))) %>%
        column_to_rownames('UniqueID') %>%
        vegan::diversity(index = index) %>%
        as.data.frame() %>% 
        rownames_to_column('UniqueID') %>%
        as_tibble() %>%
        set_names(c('UniqueID', index))
}


#' Plots the (raw) sample library sizes by body site and group
#' 
#' @param df.plot data.frame; *NOT* separated by body site, will make composite plot
#' @param stratify logical; whether or not to stratify by group as well as body site
#' defaults to FALSE
#' 
#' @returns a tidy data.frame with 2 columns: UniqueID and Shannon 
plot.lib.sizes <- function(df.plot, seq.type, stratify = FALSE, Site = ...) {
    
    # browser()
    plot <- df.plot %>%
        ggplot(aes(x = Site, y = LibSize, fill = Group)) +
        geom_boxplot(alpha = 0.7, outlier.shape = NA) + 
        geom_point(aes(shape = Group), 
                   position = position_jitterdodge(), size = 2) + 
        scale_fill_lancet() +
        theme_bw() +
        theme(legend.position = 'top', legend.justification = 'left')
    
    if ((seq.type == '16S') & stratify) {
        plot + labs(title = 'Sample Library Sizes by Body Site and COVID Status',
                 y = '16S Amplicon Reads', x = '', fill = '', shape = '') } 
    else if ((seq.type == 'WGS') & stratify) { 
        plot + labs(title = 'Sample Library Sizes by Body Site and COVID Status',
                    y = 'Raw Metagenomic Reads', x = '', fill = '', shape = '') }
    else if ((seq.type == '16S') & !stratify) {
        plot + labs(title = 'Sample Library Sizes by Body Site',
                 y = '16S Amplicon Reads', x = '', fill = '', shape = '') }
    else if ((seq.type == 'WGS') & !stratify) {
        plot + labs(title = 'Sample Library Sizes by Body Site',
                    y = 'Raw Metagenomic Reads', x = '', fill = '', shape = '') }
    else { stop('Issue with library size plotting function')}
}


pcoa <- function(df.abun, df.meta, svar = 'OSCI_Class_Worst', incl.abx = TRUE) {
    
    # browser()
    meta <- df.meta %>%
        filter(UniqueID %in% df.abun$UniqueID) %>%
        select(UniqueID, PatientID, all_of(svar), 
               contains('Abx'), contains('OSCI')) %>%
        arrange(UniqueID)
    if (!incl.abx) meta <- filter(meta, AbxCurr == 'No')
    pco.data <- df.abun %>%
        filter(UniqueID %in% meta$UniqueID) %>%
        arrange(UniqueID) %>%
        column_to_rownames('UniqueID') %>%
        as.matrix() %>%
        vegan::vegdist(method = 'bray', diag = TRUE, upper = TRUE) %>%
        stats::cmdscale(k = 2, eig = TRUE)
    plot.points <- pco.data$points %>%
        magrittr::set_colnames(c('V1','V2')) %>%
        cbind(column_to_rownames(meta, 'UniqueID')) %>%
        add.centroids(svar) %>%
        as.data.frame() %>% rownames_to_column('UniqueID')
    axes <- sprintf(fmt='%.2f',
                    (pco.data$eig[1:2]/sum(pco.data$eig[pco.data$eig > 0])) * 100)
    
    return(list('points' = plot.points,
                'axes' = axes))
}


add.centroids <- function(df, group.by) {
    
    # rlang might not work...
    pcoa.centroids <- df %>%
        group_by(!!rlang::sym(group.by)) %>%
        nest() %>%
        ungroup() %>%
        mutate(V1 = map_dbl(data, ~ mean(.$V1))) %>%
        mutate(V2 = map_dbl(data, ~ mean(.$V2))) %>%
        select(-data) %>%
        add_column(Pt.Type = 'Group Mean')
    
    df %>%
        add_column(Pt.Type = 'Individual') %>%
        bind_rows(pcoa.centroids)
}


alpha.coin.tests <- function(df.meta) {
    
    rich <- coin::spearman_test(
        Richness ~ OSCI_Score | AbxSamp,
        data = df.meta) %>% coin::pvalue()
    shan <- coin::spearman_test(
        Shannon ~ OSCI_Score | AbxSamp,
        data = df.meta) %>% coin::pvalue()
    simp <- coin::spearman_test(
        InvSimp ~ OSCI_Score | AbxSamp,
        data = df.meta) %>% coin::pvalue()
    
    return(list('Richness' = rich, 'Shannon' = shan,
                'InvSimp' = simp))
}

################################################
############## UTILITY FUNCTIONS ###############
################################################

filter.by.pct.nas <- function(df.met, cov.thresh) {
    
    sel <- df.met %>%
        select_if(is.numeric) %>%
        summarize_all(~ sum(is.na(.))) %>%
        gather('compound','n.nas') %>%
        arrange(desc(n.nas)) %>%
        mutate(pct.cov = 1-n.nas/nrow(df.met)) %>%
        filter(pct.cov >= cov.thresh) %>%
        pull(compound)
    df.met %>%
        select(contains('ID', ignore.case = FALSE), all_of(sel))
}

calc.single.prevalence <- function(vec) {
    if (all(vec == 0)) return(0)
    else return(length(vec[vec != 0])/length(vec))
}


discretize.days <- function(df, var = 'Days', keep.numeric = FALSE) {
    
    if (keep.numeric) {
        df %>% 
            mutate_at(vars(var),
                      ~ fct_relevel(., c('5-7','8-10','11-13','14-16',
                                         '17-19','20-22','23-25',
                                         '26-31','32-37','38-43','44-49'))) }
    else {
        df %>% 
            mutate_at(vars(var), 
                      ~ fct_collapse(.,
                                     `Early` = c('5-7','8-10'),
                                     `Late` = c('14-16','11-13','17-19','26-31',
                                                '32-37','20-22','38-43','44-49',
                                                '23-25'))) %>%
            mutate_at(vars(var), ~ fct_relevel(., c('Early','Late'))) }
}


remove.zero.cols <- function(df) {
    
    char <- df %>%
        select(where(is.character)) %>%
        colnames()
    keep <- df %>%
        summarize_if(is.double, ~ var(., na.rm = TRUE)) %>%
        tidyr::gather('feat','var') %>%
        filter(var != 0) %>%
        pull(feat)
    select(df, c(all_of(char), all_of(keep)))
}


rename.otus <- function(df.plot, phy.obj) {
    
    # browser()
    df.tax <- phy.obj@tax_table %>%
        as.data.frame() %>%
        rownames_to_column('OTU') %>% as_tibble() %>%
        tidyr::gather('Rank','Value',-OTU) %>%
        dplyr::filter(Rank != 'Species') %>%
        nest(data = -OTU) %>%
        mutate(TaxaID = map_chr(data, function(x) {
            x %>% 
                dplyr::filter(Value != '?') %>%
                use_series(Value) %>% 
                tail(1) })) %>% 
        mutate(TaxaUp = map_chr(data, function(x) {
            x %>% 
                dplyr::filter(Value != '?') %>%
                use_series(Value) %>% 
                tail(2) %>% head(1)
        })) %>%
        mutate(Duped = duplicated(TaxaID)) %>% 
        mutate(TaxaID = ifelse(Duped, paste(TaxaUp, TaxaID), TaxaID)) %>%
        select(-c(data, TaxaUp, Duped))
    
    df.plot %>%
        left_join(df.tax, by = c('feature' = 'OTU')) %>%
        rename(OTU = 'feature') %>% 
        rename(feature = 'TaxaID')
}


# split by body site, do some basic stats & plots
# write_tsv for each body site for potential GMM/omixer run
eval.write.functional.data <- function(dir, raw.ko.feat) {
    
    # browser()
    sample.map <- dir %>%
        load.sample.maps('WGS', TRUE) %>%
        # select(SampleID, PatientID, Measurement)
        select(SampleID, Measurement)
    
    tmp <- raw.ko.feat %>%
        tidyr::gather('SampleID', 'Value', -KO) %>%
        left_join(sample.map, by = 'SampleID') %>%
        nest(data = -Measurement) %>%
        mutate(data = map(data, ~ spread(., 'SampleID', 'Value')))
    
    ko.op <- tmp$data[[1]]
    ko.sdna <- tmp$data[[2]]
    ko.tbs <- tmp$data[[3]]
    
    # do some tests to check "library sizes" etc (what is func. equiv??)
    ko.op.func.ab <- ko.op %>%
        summarize_if(is.double, sum) %>%
        tidyr::gather('SampleID', 'FuncAb')
    ko.sdna.func.ab <- ko.sdna %>%
        summarize_if(is.double, sum) %>%
        tidyr::gather('SampleID', 'FuncAb')
    ko.tbs.func.ab <- ko.tbs %>%
        summarize_if(is.double, sum) %>%
        tidyr::gather('SampleID', 'FuncAb')
    tab.data <- summary(ko.op.func.ab$FuncAb) %>%
        cbind(summary(ko.sdna.func.ab$FuncAb), 
              summary(ko.tbs.func.ab$FuncAb)) %>%
        as_tibble() %>% set_names(c('OP','SDNA','TBS')) %>%
        add_column(Stat = c('Min.','1st Qu.','Median','Mean','3rd Qu.','Max.'),
                   .before = 1)
    plot.data <- ko.op.func.ab %>%
        bind_rows(ko.sdna.func.ab, 
                  ko.tbs.func.ab,
                  .id = 'DataType') %>%
        mutate_at(vars(DataType), ~ fct_recode(., OP = '1', 
                                               SDNA = '2',
                                               TBS = '3'))
    plot.data %>%
        ggplot(aes(x = DataType, y = FuncAb)) +
        geom_boxplot(aes(fill = DataType), alpha = 0.7, 
                     outlier.shape = NA) +
        geom_jitter() +
        scale_y_log10() +
        labs(x = '')
    
    proc.dir <- here(dir, 'proc-data')
    write_tsv(ko.op, paste0(proc.dir, '/ko-op.tsv'))
    write_tsv(ko.sdna, paste0(proc.dir, '/ko-sdna.tsv'))
    write_tsv(ko.tbs, paste0(proc.dir, '/ko-tbs.tsv'))
}


################################################
###### SANITY CHECKS - NO LONGER UPDATED #######
################################################