## SETTINGS ----
library(targets)
library(tarchetypes)
library(patchwork)

source(here::here('code','functions.R'))
source(here::here('code','mdc-util-functions.R'))
source(here::here('code','plotting-functions.R'))

tar_option_set(packages = c('magrittr','tidyverse','here',#'ggfortify',
                            'performance','phyloseq', 'ggsankey'))

## MAIN PIPELINE ----
list(
    ###################################################################
    tar_target(params, list(prevalence = 0.1,
                            fdr.cutoff = 0.05,
                            num.feat = 'auto',
                            nnodes = future::availableCores()*(3/4),
                            na.thresh.glob = 0.25)),
    ###################################################################
    ## INPUT FILES ----
    tar_target(short.names.fn,
               here('data','raw-metadata','patient-id-map.tsv'),
               format = 'file'),
    tar_target(map.all.samples.fn,
               here('data','raw-metadata','all-sampling-dates.xlsx'),
               format = 'file'),
    tar_target(map.16s.fn, 
               here('data','raw-metadata','sample-donor-map-16S.tsv'),
               format = 'file'), 
    tar_target(map.wgs.fn, 
               here('data','raw-metadata','sample-donor-map-wgs.tsv'),
               format = 'file'),
    tar_target(map.met.fn,
               here('data','raw-metadata','sampling-dates.tsv'),
               format = 'file'),
    # single measurement (static) patient characteristics
    tar_target(raw.patient.data.fn,    
               here('data','raw-metadata','me-encoded-full','patient-data-clean.tsv'), 
               format = 'file'),
    # longitudinal (repeated) clinical measurements
    tar_target(raw.sample.data.fn,   
               here('data','raw-metadata','me-encoded-full','sample-data-clean.tsv'), 
               format = 'file'),
    tar_target(raw.16s.tax.table.fn, 
               here('data','raw-metadata','all-clean-tax-table.biom'),
               format = 'file'),
    tar_target(raw.16s.fn,
               here('data','raw-data','metagenomic','all-otu.txt'),
               format = 'file'),
    tar_target(raw.motus.fn, 
               here('data','raw-data','metagenomic','all-motus-2.6.tsv'),
               format = 'file'),
    tar_target(raw.motus.order.fn,
               here('data','raw-data','metagenomic','all-motus-order.tsv'),
               format = 'file'),
    tar_target(raw.motus.family.fn,
               here('data','raw-data','metagenomic','all-motus-family.tsv'),
               format = 'file'),
    tar_target(gmms.fn, 
               here('data','raw-data','metagenomic','gmms-omixer','module-prediction')),
    tar_target(raw.isg.fn, 
               here('data','raw-data','gene-expression','cell-frequencies-isg.tsv'),
               format = 'file'),
    tar_target(raw.plasma.ck.fn,
               here('data','raw-data','cytokines','plasma.tsv'),
               format = 'file'),
    tar_target(raw.airway.ck.fn, 
               here('data','raw-data','cytokines','airway.tsv'),
               format = 'file'),
    tar_target(raw.metabolites.fn,
               here('data','raw-data','metabolites','metabolites-v3.csv'),
               format = 'file'),
    tar_target(raw.metabolites.scale.fn,
               here('data','raw-data','metabolites','metabolites-v3-median-scaled.csv'),
               format = 'file'),
    tar_target(metabolite.categ.fn, 
               here('data','raw-metadata','metabolite-classifications.xlsx'),
               format = 'file'),
    tar_target(raw.cyps.fn, 
               here('data','raw-data','gene-expression','cyps-ahr.xlsx'),
               format = 'file'),
    tar_target(raw.gene.cts.fn, 
               here('data','raw-data','metagenomic','gmgc-profiles-rare-oct.tsv'),
               format = 'file'),
    tar_target(raw.mapfile.fn, 
               here('data','raw-data','metagenomic','gmgc-filtered-map-kegg.tsv'),
               format = 'file'),
    tar_target(module.defs.fn,
               here('data','raw-metadata','module-defs.tsv'),
               format = 'file'),
    
    ## EXPERIMENT METADATA (MAP)FILES, SHORT IDS ----
    tar_target(sample.map.plasma, 
               readxl::read_xlsx(map.all.samples.fn, sheet = 1)),
    tar_target(sample.map.op, 
               readxl::read_xlsx(map.all.samples.fn, sheet = 2)),
    tar_target(sample.map.urine, 
               readxl::read_xlsx(map.all.samples.fn, sheet = 3)),
    tar_target(sample.map.stool, 
               readxl::read_xlsx(map.all.samples.fn, sheet = 4)),
    tar_target(sample.map.list, list('P' = sample.map.plasma,
                                     'U' = sample.map.urine,
                                     'OP' = sample.map.op,
                                     'SDNA' = sample.map.stool)),
    tar_target(sample.discrepancies, 
               clean.master.mapfile(sample.map.list, clean = FALSE)),
    tar_target(sample.map, 
               clean.master.mapfile(sample.map.list, clean = TRUE)),
    # tar_target(master.map, combine.sample.maps(sample.data, # clinical
    #                                            sample.map,  # visits
    #                                            sample.discrepancies, # visits
    #                                            short.names)),
    tar_target(sample.map.16s, read_tsv(map.16s.fn)),
    tar_target(sample.map.wgs, read_tsv(map.wgs.fn)),
    # same information as sample.map.plasma above...
    tar_target(sample.map.met, read_tsv(map.met.fn)),
    tar_target(sample.donor.mapfiles,
               list('sample.map.16s' = sample.map.16s, 
                    'sample.map.wgs' = sample.map.wgs, 
                    'sample.map.met' = sample.map.met)),
    tar_target(short.names, read_tsv(short.names.fn)),
    
    ## RAW FEATURE DATA TABLES ----
    tar_target(raw.patient.data, read_tsv(raw.patient.data.fn)),
    tar_target(raw.sample.data, read_tsv(raw.sample.data.fn)),
    tar_target(raw.motus, read_tsv(raw.motus.fn)),
    tar_target(raw.motus.order, read_tsv(raw.motus.order.fn)),
    tar_target(raw.motus.family, read_tsv(raw.motus.family.fn)),
    tar_target(raw.isg, read_tsv(raw.isg.fn)),
    tar_target(raw.plasma.ck, read_tsv(raw.plasma.ck.fn)),
    tar_target(raw.airway.ck, read_tsv(raw.airway.ck.fn)),
    tar_target(raw.cyps, readxl::read_xlsx(raw.cyps.fn)),
    tar_target(raw.metabolites, read_csv(raw.metabolites.fn)),
    tar_target(raw.metabolites.med, read_csv(raw.metabolites.scale.fn)),
    tar_target(raw.clinical.feat, clean.clinical.data(seq.type = 'WGS', 
                                                      mapfiles = sample.donor.mapfiles, 
                                                      raw.patient.data,
                                                      raw.sample.data)),
    tar_target(raw.motu.feat, clean.otu.data(seq.type = 'WGS', 
                                            raw.file = raw.motus,
                                            mapfile = sample.map.wgs, 
                                            unmapped = FALSE) %>%
                   rename(Unmapped_mOTUs = 'unassigned')),
    tar_target(raw.motu.order.feat, clean.otu.data(seq.type = 'WGS', 
                                             raw.file = raw.motus.order,
                                             mapfile = sample.map.wgs, 
                                             unmapped = FALSE) %>%
                   rename(Unmapped_mOTUs = 'unassigned')),
    tar_target(raw.motu.fam.feat, clean.otu.data(seq.type = 'WGS', 
                                                   raw.file = raw.motus.family,
                                                   mapfile = sample.map.wgs, 
                                                   unmapped = FALSE) %>%
                   rename(Unmapped_mOTUs = 'unassigned')),
    
    ###################################################################
    ## CLEAN PATIENT & CLINICAL DATA TABLES ----
    tar_target(patient.data, clean.patient.data(raw.patient.data, short.names) %>%
                   filter.patient.data() %>%
                   mutate(OSCI_Score = ifelse(ShortID=='P23',8,OSCI_Score)) %>%
                   mutate(Death = ifelse(ShortID=='P23',1,Death)) %>%
                   mutate(OSCI_Class = case_when(OSCI_Score==0 ~ 'Control',
                                                 OSCI_Score<=4 ~ 'Mild',
                                                 OSCI_Score>4~'Severe')) %>%
                   mutate(AbxHosp = fct_recode(AbxFreeHosp,`0` = '1', `1` = '0')) %>%
                   select(-AbxFreeHosp)),
    tar_target(clinical.feat.tmp, raw.clinical.feat %>%
                   # fix misnamed P05-V3 metabolite measurement
                   add_row(tibble(DataType = 'MET', Site = 'P', 
                                  raw.clinical.feat %>%
                                      dplyr::filter(Site=='OP' & PatientID=='C19-CB-50') %>%
                                      select(everything(), -DataType, -Site))) %>%
                   # fix missing P16-V4 metabolite measurement
                   add_row(tibble(DataType = 'MET', Site = 'U', 
                                  raw.clinical.feat %>%
                                      dplyr::filter(Site=='SDNA' & PatientID=='C19-CB-119') %>%
                                      select(everything(), -DataType, -Site))) %>%
                   left_join(short.names) %>%
                   mutate(GenID = paste(ShortID, Visit, sep = '-')) %>%
                   mutate(UniqueID = paste(ShortID, Site, Visit, sep = '-')) %>%
                   select(UniqueID, GenID, everything(), -ShortID, -DataType)),
    tar_target(clinical.feat, clinical.feat.tmp %>%
                   # add duplicate samplings
                   add_row(tibble(UniqueID = 'P17-OP-V4.2', clinical.feat.tmp %>% 
                                      dplyr::filter(UniqueID == 'P17-OP-V4') %>% 
                                      select(everything(), -UniqueID))) %>%
                   add_row(tibble(UniqueID = 'P11-OP-V2.2', clinical.feat.tmp %>% 
                                      dplyr::filter(UniqueID == 'P11-OP-V2') %>% 
                                      select(everything(), -UniqueID))) %>% 
                   add_row(tibble(UniqueID = 'P17-SDNA-V4.2', clinical.feat.tmp %>% 
                                      dplyr::filter(UniqueID == 'P17-SDNA-V4') %>% 
                                      select(everything(), -UniqueID))) %>%
                   add_row(tibble(UniqueID = 'P11-SDNA-V2.2', clinical.feat.tmp %>% 
                                      dplyr::filter(UniqueID == 'P11-SDNA-V2') %>% 
                                      select(everything(), -UniqueID))) %>%
                   add_row(tibble(UniqueID = 'P17-TBS-V4.2', clinical.feat.tmp %>% 
                                      dplyr::filter(UniqueID == 'P17-TBS-V4') %>% 
                                      select(everything(), -UniqueID))) %>%
                   dplyr::filter(!UniqueID %in% c('P05-P-V2')) %>%
                   mutate(OSCI_Score_Samp = ifelse(str_detect(UniqueID,'K'),
                                              0, as.numeric(OSCI_Score))) %>%
                   mutate(OSCI_Class_Samp = case_when(OSCI_Score_Samp == 0 ~ 'Control',
                                                 OSCI_Score_Samp <= 4 ~ 'Mild',
                                                 OSCI_Score_Samp > 4 ~ 'Severe')) %>%
                   select(-c(OSCI_Score, OSCI_Class)) %>%
                   mutate(PatientID = str_split(GenID,'-') %>% 
                              map_chr(~head(.,1)))),
    tar_target(vl, clinical.feat %>%
                   select(UniqueID, contains('VL'), Antibiotics) %>%
                   mutate(UniqueID = str_remove_all(UniqueID, '-OP')) %>%
                   mutate(UniqueID = str_remove_all(UniqueID, '-SDNA')) %>%
                   distinct(UniqueID, .keep_all = TRUE) %>%
                   mutate_at(vars(c(2:ncol(.))), ~ as.double(.)) ),
    tar_target(sdna.clinical.feat, filter.clinical.data(clinical.feat, 
                                                        site = 'SDNA', 
                                                        days.disc = FALSE) %>%
                   mutate(GenID = UniqueID) %>%
                   mutate(UniqueID = paste0(PatientID,'-SDNA-',
                                            str_split(UniqueID,'-') %>%
                                                map_chr(~ tail(.,1)))) %>%
                   distinct(UniqueID, .keep_all = TRUE)),
    tar_target(op.clinical.feat, filter.clinical.data(clinical.feat, 
                                                      site = 'OP', 
                                                      days.disc = FALSE) %>%
                   mutate(GenID = UniqueID) %>%
                   mutate(UniqueID = paste0(PatientID,'-OP-',
                                            str_split(UniqueID,'-') %>%
                                                map_chr(~ tail(.,1)))) %>%
                   distinct(UniqueID, .keep_all = TRUE)),
    tar_target(tbs.clinical.feat, filter.clinical.data(clinical.feat, 
                                                       site = 'TBS',
                                                       days.disc = FALSE) %>%
                   mutate(GenID = UniqueID) %>%
                   mutate(UniqueID = paste0(PatientID,'-TBS-',
                                            str_split(UniqueID,'-') %>%
                                                map_chr(~ tail(.,1)))) %>%
                   distinct(UniqueID, .keep_all = TRUE)),
    tar_target(gen.days, clinical.feat %>%
                   select(UniqueID, PatientID, Days) %>%
                   mutate(UniqueID = str_remove_all(UniqueID, '-SDNA')) %>%
                   mutate(UniqueID = str_remove_all(UniqueID, '-OP')) %>%
                   mutate(UniqueID = str_remove_all(UniqueID, '-TBS')) %>%
                   mutate(PatientID = str_split(UniqueID, '-') %>% 
                              map_chr(~ head(., 1))) %>%
                   mutate_at(vars('Days'), 
                             ~ fct_collapse(.,
                                            `Early` = c('5-7','8-10'),
                                            `Late` = c('14-16','11-13',
                                                       '17-19','26-31',
                                                       '32-37','20-22',
                                                       '38-43','44-49',
                                                       '23-25'))) %>%
                   distinct()),
    tar_target(sdna.meta, motu.feat %>%
                   filter(Site == 'SDNA') %>%
                   select(UniqueID, GenID) %>%
                   mutate(PatientID = str_split(UniqueID, '-') %>%
                              map_chr(~ head(., 1))) %>%
                   left_join(select(patient.data, -PatientID),
                             by=c('PatientID'='ShortID')) %>%
                   left_join(motu.lib.sizes) %>%
                   left_join(select(sdna.clinical.feat,
                                    UniqueID, contains('OSCI'), Days, CRP,
                                    Antibiotics, contains('VL')))),
    tar_target(op.meta, motu.feat %>%
                   filter(Site == 'OP') %>%
                   select(UniqueID, GenID) %>%
                   mutate(PatientID = str_split(UniqueID, '-') %>%
                              map_chr(~ head(., 1))) %>%
                   left_join(select(patient.data, -PatientID),
                             by=c('PatientID'='ShortID')) %>%
                   left_join(motu.lib.sizes) %>%
                   left_join(select(op.clinical.feat,
                                    UniqueID, contains('OSCI'), Days, CRP,
                                    Antibiotics, contains('VL')))),
    tar_target(tbs.meta, motu.feat %>%
                   filter(Site == 'TBS') %>%
                   select(UniqueID, GenID) %>%
                   mutate(PatientID = str_split(UniqueID, '-') %>%
                              map_chr(~ head(., 1))) %>%
                   left_join(select(patient.data, -PatientID),
                             by=c('PatientID'='ShortID')) %>%
                   left_join(motu.lib.sizes) %>%
                   left_join(select(op.clinical.feat,
                                    UniqueID, contains('OSCI'), Days, CRP,
                                    Antibiotics, contains('VL')))),
    
    ## CLEAN SINGLE CELL (PBMC) AND CYTOKINE DATA TABLES ----
    tar_target(isg.pbmc, clean.freqs.isg(raw.isg)),
    tar_target(cyps.pbmc, clean.cyps(raw.cyps, isg.pbmc)),
    tar_target(combined.pbmc, left_join(cyps.pbmc, 
                                        isg.pbmc,
                                        by = 'UniqueID')),
    tar_target(cytokines, clean.cytokines(raw.plasma.ck, 
                                          raw.airway.ck, 
                                          short.names) %>%
                   select(-Site, -Visit) %>%
                   mutate(UniqueID = str_remove_all(UniqueID, '-SDNA')) %>%
                   mutate(UniqueID = str_remove_all(UniqueID, '-OP')) %>%
                   spread('Cytokine','Value') %>%
                   mutate(PatientID = str_split(UniqueID,'-') %>%
                              map_chr(~ head(., 1))) %>%
                   mutate_at(vars(UniqueID), ~ ifelse(. == 'P15-V4',
                                                      'P15-V3', .)) %>%
                   left_join(gen.days) %>%
                   select(UniqueID, PatientID, Days, everything())),
    
    ## CLEAN METABOLITE DATA TABLES ----
    tar_target(metabolites, clean.metabolites(raw.metabolites, short.names) %>%
                   filter(!is.na(Site))),
    tar_target(metabolites.med, clean.metabolites(raw.metabolites.med,
               short.names) %>% filter(!is.na(Site))),
    tar_target(metab.annotated, metabolites %>%
                   dplyr::rename(LongID = 'PatientID') %>%
                   mutate(PatientID = str_split(UniqueID, '-') %>%
                              map_chr(~ head(., 1))) %>%
                   mutate(GenID = paste0(PatientID,'-',Visit)) %>%
                   left_join(select(clinical.feat, GenID, Site, Days, 
                                    contains('OSCI'))) %>%
                   left_join(select(patient.data, ShortID, contains('OSCI')),
                             by = c('PatientID'='ShortID')) %>%
                   rename(OSCI_Score_Worst = 'OSCI_Score',
                          OSCI_Class_Worst = 'OSCI_Class') %>%
                   select(UniqueID, LongID, PatientID, GenID,
                          everything(), -Visit) %>% 
                   discretize.days(keep.numeric = TRUE) %>% 
                   mutate(Days_Disc = Days) %>% 
                   discretize.days(var='Days_Disc') %>% 
                   select(1:10, 15, 11:14)),
    tar_target(metab.distinct, metab.annotated %>%
                   distinct(across(UniqueID), .keep_all = TRUE) %>%
                   select(-c(Method, Batch, Compound,Value)) %>% 
                   discretize.days(keep.numeric = TRUE) %>% 
                   mutate(Days_Disc = Days) %>% 
                   discretize.days(var = 'Days_Disc')), 
    # CCM metabolite names are messed up/not matchable!
    tar_target(metabolite.categ, readxl::read_xlsx(metabolite.categ.fn) %>%
                   clean.metabolite.map(metabolites, .)),
    tar_target(met.of.interest, metabolite.categ %>% 
                   filter(
                       str_detect(ChemClass,
                                  'Bile acids')|str_detect(
                                      Pathway,'Tryp|ynur|indole'))),
    tar_target(p.metabolites, filter.metabolites(metabolites, 
                                                 gen.days, 
                                                 site = 'P')),
    tar_target(p.metabolites.wide, p.metabolites %>%
                   r.friendly.metabolite.names(site = 'P',
                                               incl.site = TRUE)),
    # matched to microbiota 2x samplings for metadeconfoundR
    tar_target(p.metabolites.wide.full, p.metabolites.wide %>%
                   add_row(tibble(UniqueID = 'P17-V4.2',
                                  p.metabolites.wide %>%
                                      filter(UniqueID == 'P17-V4') %>%
                                      select(everything(), -UniqueID))) %>%
                   add_row(tibble(UniqueID = 'P11-V2.2',
                                  p.metabolites.wide %>%
                                      filter(UniqueID == 'P11-V2') %>%
                                      select(everything(), -UniqueID)))),
    tar_target(u.metabolites, filter.metabolites(metabolites, 
                                                gen.days,
                                                site = 'U')),
    tar_target(u.metabolites.wide, u.metabolites %>%
                   r.friendly.metabolite.names(site = 'U',
                                               incl.site = TRUE)),
    tar_target(u.metabolites.wide.full, u.metabolites.wide %>%
                   add_row(tibble(UniqueID = 'P17-V4.2', 
                                  u.metabolites.wide %>%
                                      filter(UniqueID == 'P17-V4') %>%
                                      select(everything(), -UniqueID))) %>%
                   add_row(tibble(UniqueID = 'P11-V2.2', 
                                  u.metabolites.wide %>%
                                      filter(UniqueID == 'P11-V2') %>%
                                      select(everything(), -UniqueID)))),
    # median-scaled
    tar_target(p.medscale.metabolites.wide, metabolites.med %>%
                   filter.metabolites(gen.days, site = 'P') %>%
                   r.friendly.metabolite.names(site = 'P', incl.site = TRUE)),
    tar_target(u.medscale.metabolites.wide, metabolites.med %>%
                   filter.metabolites(gen.days, site = 'U') %>%
                   r.friendly.metabolite.names(site = 'U', incl.site = TRUE)),
    # autoscaled
    tar_target(p.metabolites.as, p.metabolites.wide %>%
                   mutate_if(is.double, ~ scale(.))),
    tar_target(u.metabolites.as, u.metabolites.wide %>%
                   mutate_if(is.double, ~ scale(.))),
    
    ## CLEAN 16S PHYLOSEQ OBJECTS ----
    tar_target(raw.16s, read_delim(raw.16s.fn, delim='\t') %>%
                   gather('SampleID','Count',-OTU) %>%
                   left_join(select(sample.map.16s, 
                                    SampleID,Treatment,Donor_No,Visite)) %>%
                   filter(!str_detect(Donor_No, 'BLANK')) %>%
                   left_join(short.names, by=c('Donor_No'='PatientID')) %>%
                   mutate(UniqueID = paste0(ShortID,'-',Treatment,'-',Visite)) %>%
                   mutate(GenID = paste0(ShortID,'-',Visite)) %>%
                   select(Treatment, UniqueID, GenID, SampleID, OTU, Count) %>%
                   arrange(SampleID) %>%
                   spread(OTU, Count) %>%
                   nest(data = -Treatment) %>%
                   mutate(data = map(data, ~ mutate(., dup = duplicated(UniqueID)) %>%
                                         mutate(UniqueID = ifelse(
                                             dup,paste0(UniqueID,'.2'),UniqueID)) %>%
                                         mutate(GenID = ifelse(dup, paste0(GenID,'.2'), GenID)) %>%
                                         select(-dup))) %>%
                   unnest(data)),
    tar_target(raw.tax.table, read_delim(raw.16s.tax.table.fn, delim=',') %>%
                   mutate_all(~ trimws(.)) %>%
                   magrittr::set_colnames(c('OTU','Kingdom','Phylum','Class','Order',
                                  'Family','Genus','Species'))),
    tar_target(sdna.ps.raw, raw.16s %>% 
                   filter(Treatment == 'SDNA') %>% 
                   select(UniqueID, contains('OTU')) %>%
                   arrange(UniqueID) %>%
                   column_to_rownames('UniqueID') %>%
                   phyloseq::otu_table(taxa_are_rows = FALSE) %>%
                   phyloseq::phyloseq(tax_table(raw.tax.table %>%
                                                    column_to_rownames('OTU') %>%
                                                    as.matrix()),
                                      sample_data(raw.16s %>%
                                                      filter(Treatment == 'SDNA') %>%
                                                      select(contains('ID')) %>%
                                                      arrange(UniqueID) %>%
                                                      left_join(sdna.meta) %>%
                                                      column_to_rownames('UniqueID')))),
    tar_target(op.ps.raw, raw.16s %>% 
                   filter(Treatment == 'OP') %>% 
                   select(UniqueID, contains('OTU')) %>%
                   arrange(UniqueID) %>%
                   column_to_rownames('UniqueID') %>%
                   phyloseq::otu_table(taxa_are_rows = FALSE) %>%
                   phyloseq::phyloseq(tax_table(raw.tax.table %>%
                                                    column_to_rownames('OTU') %>%
                                                    as.matrix()),
                                      sample_data(raw.16s %>%
                                                      filter(Treatment == 'OP') %>%
                                                      select(contains('ID')) %>%
                                                      arrange(UniqueID) %>%
                                                      left_join(op.meta) %>%
                                                      column_to_rownames('UniqueID')))),
    tar_target(tbs.ps.raw, raw.16s %>% 
                   filter(Treatment == 'TBS') %>% 
                   select(UniqueID, contains('OTU')) %>%
                   arrange(UniqueID) %>%
                   column_to_rownames('UniqueID') %>%
                   phyloseq::otu_table(taxa_are_rows = FALSE) %>%
                   phyloseq::phyloseq(tax_table(raw.tax.table %>%
                                                    column_to_rownames('OTU') %>%
                                                    as.matrix()),
                                      sample_data(raw.16s %>%
                                                      filter(Treatment == 'TBS') %>%
                                                      select(contains('ID')) %>%
                                                      arrange(UniqueID) %>%
                                                      left_join(tbs.meta) %>%
                                                      column_to_rownames('UniqueID')))),
    tar_target(sdna.ps.phy, sdna.ps.raw %>%
                   transform_sample_counts(function(OTU) OTU/sum(OTU)) %>%
                   tax_glom('Phylum') %>%
                   subset_taxa(Phylum != '?') %>%
                   filter_taxa(function(x) calc.single.prevalence(x) > params$prevalence, 
                               TRUE)),
    tar_target(op.ps.phy, op.ps.raw %>%
                   transform_sample_counts(function(OTU) OTU/sum(OTU)) %>%
                   tax_glom('Phylum') %>%
                   subset_taxa(Phylum != '?') %>%
                   filter_taxa(function(x) calc.single.prevalence(x) > params$prevalence, 
                               TRUE)),
    tar_target(tbs.ps.phy, tbs.ps.raw %>%
                   transform_sample_counts(function(OTU) OTU/sum(OTU)) %>%
                   tax_glom('Phylum') %>%
                   subset_taxa(Phylum != '?') %>%
                   filter_taxa(function(x) calc.single.prevalence(x) > params$prevalence, 
                               TRUE)),
    
    ## CLEAN TAXONOMIC METAGENOMIC DATA TABLES ----
        # should be 81 SDNA samples and 94 OP samples ... 
        # species/motus level counts
    tar_target(motu.feat, raw.motu.feat %>%
                   left_join(short.names) %>%
                   mutate(GenID = paste0(ShortID,'-',Visit)) %>%
                   mutate(UniqueID = paste0(ShortID,'-',Site,'-',Visit)) %>%
                   mutate(dup = duplicated(UniqueID)) %>%
                   mutate(UniqueID = ifelse(
                       dup,paste0(UniqueID,'.2'),UniqueID)) %>%
                   mutate(GenID = ifelse(dup, paste0(GenID,'.2'), GenID)) %>%
                   select(-dup) %>%
                   left_join(select(clinical.feat, UniqueID, Days)) %>%
                   select(Site, GenID, UniqueID, Days, everything(), 
                          -PatientID, -SampleID, -ShortID, -Date, -Visit) %>%
                   distinct(UniqueID, .keep_all = TRUE)),
    tar_target(motu.order.feat, raw.motu.order.feat %>%
                   left_join(short.names) %>%
                   mutate(GenID = paste0(ShortID,'-',Visit)) %>%
                   mutate(UniqueID = paste0(ShortID,'-',Site,'-',Visit)) %>%
                   mutate(dup = duplicated(UniqueID)) %>%
                   mutate(UniqueID = ifelse(
                       dup,paste0(UniqueID,'.2'),UniqueID)) %>%
                   mutate(GenID = ifelse(dup, paste0(GenID,'.2'), GenID)) %>%
                   select(-dup) %>%
                   left_join(select(clinical.feat, UniqueID, Days)) %>%
                   select(Site, GenID, UniqueID, Days, everything(), 
                          -PatientID, -SampleID, -ShortID, -Date, -Visit) %>%
                   distinct(UniqueID, .keep_all = TRUE)),
    tar_target(order.sdna.feat, motu.order.feat %>%
                   clean.motus('SDNA', rel.ab = TRUE, 
                               h.cols = 4)),
    tar_target(order.op.feat, motu.order.feat %>%
                   clean.motus('OP', rel.ab = TRUE, 
                               h.cols = 4)),
    # not returning correctly
    tar_target(motu.fam.feat, raw.motu.fam.feat %>%
                   left_join(short.names) %>%
                   mutate(GenID = paste0(ShortID,'-',Visit)) %>%
                   mutate(UniqueID = paste0(ShortID,'-',Site,'-',Visit)) %>%
                   mutate(dup = duplicated(UniqueID)) %>%
                   mutate(UniqueID = ifelse(
                       dup,paste0(UniqueID,'.2'),UniqueID)) %>%
                   mutate(GenID = ifelse(dup, paste0(GenID,'.2'), GenID)) %>%
                   select(-dup) %>%
                   left_join(select(clinical.feat, UniqueID, Days)) %>%
                   select(Site, GenID, UniqueID, Days, everything(), 
                          -PatientID, -SampleID, -ShortID, -Date, -Visit) %>%
                   distinct(UniqueID, .keep_all = TRUE)),
    tar_target(fam.sdna.feat, motu.fam.feat %>%
                   clean.motus('SDNA', rel.ab = TRUE, 
                               h.cols = 4)),
    tar_target(fam.op.feat, motu.fam.feat %>%
                   clean.motus('OP', rel.ab = TRUE, 
                               h.cols = 4)),
    
    # filtered species abundances -- counts or relative abundances
    tar_target(sdna.feat, clean.motus(motu.feat, 
                                      site = 'SDNA', 
                                      rel.ab = FALSE)),
    tar_target(sdna.feat.ra, clean.motus(motu.feat,
                                         site = 'SDNA',
                                         rel.ab = TRUE)),
    tar_target(op.feat, clean.motus(motu.feat,
                                    site = 'OP',
                                    rel.ab = FALSE)),
    tar_target(op.feat.ra, clean.motus(motu.feat,
                                       site = 'OP',
                                       rel.ab = TRUE)),
    tar_target(tbs.feat, clean.motus(motu.feat,
                                     site = 'TBS',
                                     rel.ab = FALSE)),
    
    # filtered genus RELATIVE abundances
    tar_target(sdna.genus.feat, bin.motus(sdna.feat, 
                                          sdna.clinical.feat, 
                                          rel.ab = TRUE)),
    tar_target(op.genus.feat, bin.motus(op.feat,
                                        op.clinical.feat,
                                        rel.ab = TRUE)),
    tar_target(tbs.genus.feat, bin.motus(tbs.feat,
                                         tbs.clinical.feat,
                                         rel.ab = TRUE)),
    
    ## CLEAN FUNCTIONAL METAGENOMIC DATA TABLES (INCOMPLETE) ----
    tar_target(module.defs, vroom::vroom(module.defs.fn,
                                         delim = '\t',
                                         col_names = FALSE,
                                         col_select = c(1:2)) %>%
                   set_colnames(c('KEGG_Module', 'Description')) %>%
                   filter(!str_detect(KEGG_Module, '_'))), 
    tar_target(gmm.defs, omixerRpm::loadDB('GMMs.v1.07')@module.names %>%
                   rownames_to_column('Module') %>%
                   rename(Description = 'V2')),
    tar_target(mod.defs, module.defs %>%
                   rename(feature = 'KEGG_Module') %>%
                   bind_rows(rename(gmm.defs, feature = 'Module'))),
    # tar_target(raw.gene.cts, vroom::vroom(raw.gene.cts.fn,
    #                                       delim = '\t')),
    tar_target(ko.binned, bin.gene.abundances(tibble(), # annotated.genes,
                                              tibble(), # gene.ko.map,
                                              by = 'KEGG_ko',
                                              preloaded = TRUE)),
    tar_target(sdna.ko, structure.kegg.table(ko.binned, 
                                             clinical.feat,
                                             by = 'KEGG_ko',
                                             sample.map.wgs, 
                                             short.names)),
    tar_target(gmm.binned, load.omixer.table(gmms.fn, wide = TRUE)),
    tar_target(sdna.gmm, structure.gmm.table(gmm.binned, clinical.feat)),
    tar_target(kegg.binned, bin.gene.abundances(tibble(), # annotated.genes,
                                                tibble(), # gene.module.map,
                                                by = 'KEGG_Module',
                                                preloaded = TRUE)),
    tar_target(sdna.kegg, structure.kegg.table(kegg.binned, 
                                               clinical.feat,
                                               by = 'KEGG_Module',
                                               sample.map.wgs, 
                                               short.names)),
    
    ###################################################################
    ## ALPHA DIVERSITY, PREVALENCE, MEAN RELATIVE ABUNDANCES ----
    tar_target(motu.lib.sizes, motu.feat %>%
                   select(-GenID, -Days) %>%
                   nest(data = -Site) %>%
                   mutate(data = map(data, 
                                     ~ calc.lib.size(., feat.id = 'OTU'))) %>%
                   unnest(data) %>%
                   select(-Site)),
    tar_target(sdna.richness, sdna.feat %>%
                   select(-c(1,2,4)) %>%
                   mutate_if(is.double, ~ ifelse(. == 0, 0, 1)) %>%
                   mutate(Richness = rowSums(.[,2:ncol(.)])) %>%
                   select(UniqueID, Richness)),
    tar_target(op.richness, op.feat %>%
                   select(-c(1,2,4)) %>%
                   mutate_if(is.double, ~ ifelse(. == 0, 0, 1)) %>%
                   mutate(Richness = rowSums(.[,2:ncol(.)])) %>%
                   select(UniqueID, Richness)),
    tar_target(sdna.shannon, sdna.feat %>%
                   select(-c(1,2,4)) %>%
                   column_to_rownames('UniqueID') %>% 
                   vegan::diversity('shannon') %>% 
                   enframe() %>% 
                   set_names(c('UniqueID','Shannon'))),
    tar_target(op.shannon, op.feat %>%
                   select(-c(1,2,4)) %>%
                   column_to_rownames('UniqueID') %>% 
                   vegan::diversity('shannon') %>%
                   enframe() %>% 
                   set_names(c('UniqueID','Shannon'))),
    tar_target(sdna.invsimp, sdna.feat %>%
                   select(-c(1,2,4)) %>%
                   column_to_rownames('UniqueID') %>% 
                   vegan::diversity('invsimp') %>% 
                   enframe() %>% 
                   set_names(c('UniqueID','InvSimp'))),
    tar_target(op.invsimp, op.feat %>%
                   select(-c(1,2,4)) %>%
                   column_to_rownames('UniqueID') %>% 
                   vegan::diversity('invsimp') %>%
                   enframe() %>% 
                   set_names(c('UniqueID','InvSimp'))),
    tar_target(sdna.alpha, sdna.richness %>%
                   left_join(sdna.shannon) %>%
                   left_join(sdna.invsimp)),
    tar_target(op.alpha, op.richness %>%
                   left_join(op.shannon) %>%
                   left_join(op.invsimp)),
    # significance from coin tests
    tar_target(p.sdna.alpha, sdna.alpha %>% 
                   left_join(select(sdna.meta, 
                                    UniqueID,contains(c('OSCI','Abx')))) %>%
                   alpha.coin.tests()),
    tar_target(p.op.alpha, op.alpha %>% 
                   left_join(select(op.meta, 
                                    UniqueID,contains(c('OSCI','Abx')))) %>%
                   alpha.coin.tests()),
    
    # SDNA prevalences
    tar_target(sdna.motus.prev, sdna.feat.ra %>%
                   summarize_if(is.double, ~ calc.single.prevalence(.)) %>%
                   tidyr::gather('feature','prevalence')),
    tar_target(sdna.motus.pat.prev, sdna.feat.ra %>%
                   filter(!str_detect(UniqueID,'K')) %>%
                   summarize_if(is.double, ~ calc.single.prevalence(.)) %>%
                   tidyr::gather('feature','prevalence')),
    tar_target(sdna.genus.prev, sdna.genus.feat %>%
                   summarize_if(is.double, ~ calc.single.prevalence(.)) %>%
                   tidyr::gather('feature','prevalence')),
    tar_target(sdna.genus.pat.prev, sdna.genus.feat %>%
                   filter(!str_detect(UniqueID,'K')) %>%
                   summarize_if(is.double, ~ calc.single.prevalence(.)) %>%
                   tidyr::gather('feature','prevalence')),
    # OP prevalences
    tar_target(op.motus.prev, op.feat.ra %>%
                   summarize_if(is.double, ~ calc.single.prevalence(.)) %>%
                   tidyr::gather('feature','prevalence')),
    tar_target(op.motus.pat.prev, op.feat.ra %>%
                   filter(!str_detect(UniqueID,'K')) %>%
                   summarize_if(is.double, ~ calc.single.prevalence(.)) %>%
                   tidyr::gather('feature','prevalence')),
    tar_target(op.genus.prev, op.genus.feat %>%
                   summarize_if(is.double, ~ calc.single.prevalence(.)) %>%
                   tidyr::gather('feature','prevalence')),
    tar_target(op.genus.pat.prev, op.genus.feat %>%
                   filter(!str_detect(UniqueID,'K')) %>%
                   summarize_if(is.double, ~ calc.single.prevalence(.)) %>%
                   tidyr::gather('feature','prevalence')),
    # SDNA mean relative abundances
    tar_target(sdna.motus.mean.ra, sdna.feat.ra %>%
                   summarize_if(is.double, ~ mean(.)) %>%
                   tidyr::gather('feature','mean.ra')),
    tar_target(sdna.motus.pat.mean.ra, sdna.feat.ra %>%
                   filter(!str_detect(UniqueID,'K')) %>%
                   summarize_if(is.double, ~ mean(.)) %>%
                   tidyr::gather('feature','mean.ra')),
    tar_target(sdna.genus.mean.ra, sdna.genus.feat %>%
                   summarize_if(is.double, ~ mean(.)) %>%
                   tidyr::gather('feature','mean.ra')),
    tar_target(sdna.genus.pat.mean.ra, sdna.genus.feat %>%
                   filter(!str_detect(UniqueID,'K')) %>%
                   summarize_if(is.double, ~ mean(.)) %>%
                   tidyr::gather('feature','mean.ra')),
    # OP mean relative abundances
    tar_target(op.motus.mean.ra, op.feat.ra %>%
                   summarize_if(is.double, ~ mean(.)) %>%
                   tidyr::gather('feature','mean.ra')),
    tar_target(op.motus.pat.mean.ra, op.feat.ra %>%
                   filter(!str_detect(UniqueID,'K')) %>%
                   summarize_if(is.double, ~ mean(.)) %>%
                   tidyr::gather('feature','mean.ra')),
    tar_target(op.genus.mean.ra, op.genus.feat %>%
                   summarize_if(is.double, ~ mean(.)) %>%
                   tidyr::gather('feature','mean.ra')),
    tar_target(op.genus.pat.mean.ra, op.genus.feat %>%
                   filter(!str_detect(UniqueID,'K')) %>%
                   summarize_if(is.double, ~ mean(.)) %>%
                   tidyr::gather('feature','mean.ra')),
    # combined stats
    tar_target(sdna.motus.stats, sdna.motus.prev %>%
                   left_join(sdna.motus.mean.ra, by='feature') %>%
                   left_join(sdna.motus.pat.prev, by='feature') %>%
                   left_join(sdna.motus.pat.mean.ra, by='feature') %>%
                   set_names(c('feature','prev.all','mean.ra.all',
                               'prev.pats','mean.ra.pats'))),
    tar_target(sdna.genus.stats, sdna.genus.prev %>%
                   left_join(sdna.genus.mean.ra, by='feature') %>%
                   left_join(sdna.genus.pat.prev, by='feature') %>%
                   left_join(sdna.genus.pat.mean.ra, by='feature') %>%
                   set_names(c('feature','prev.all','mean.ra.all',
                               'prev.pats','mean.ra.pats'))),
    tar_target(op.motus.stats, op.motus.prev %>%
                   left_join(op.motus.mean.ra, by='feature') %>%
                   left_join(op.motus.pat.prev, by='feature') %>%
                   left_join(op.motus.pat.mean.ra, by='feature') %>%
                   set_names(c('feature','prev.all','mean.ra.all',
                               'prev.pats','mean.ra.pats'))),
    tar_target(op.genus.stats, op.genus.prev %>%
                   left_join(op.genus.mean.ra, by='feature') %>%
                   left_join(op.genus.pat.prev, by='feature') %>%
                   left_join(op.genus.pat.mean.ra, by='feature') %>%
                   set_names(c('feature','prev.all','mean.ra.all',
                               'prev.pats','mean.ra.pats'))),
    # lists of features at predetermined cutoffs
    tar_target(sdna.motus.cutoff.feat,
               list('prev_10' = sdna.motus.stats %>%
                        filter(prev.all >= 0.1) %>% pull(feature),
                    'prev_pat_10' = sdna.motus.stats %>%
                        filter(prev.pats >= 0.1) %>% pull(feature),
                    'prev_20' = sdna.motus.stats %>%
                        filter(prev.all >= 0.2) %>% pull(feature),
                    'prev_pat_20' = sdna.motus.stats %>%
                        filter(prev.pats >= 0.2) %>% pull(feature),
                    'prev_50' = sdna.motus.stats %>%
                        filter(prev.all >= 0.5) %>% pull(feature),
                    'prev_pat_50' = sdna.motus.stats %>%
                        filter(prev.pats >= 0.5) %>% pull(feature),
                    'ra_e4' = sdna.motus.stats %>%
                        filter(mean.ra.all >= 1e-4) %>% pull(feature),
                    'ra_pat_e4' = sdna.motus.stats %>%
                        filter(mean.ra.pats >= 1e-4) %>% pull(feature))),
    tar_target(op.motus.cutoff.feat,
               list('prev_10' = op.motus.stats %>%
                        filter(prev.all >= 0.1) %>% pull(feature),
                    'prev_pat_10' = op.motus.stats %>%
                        filter(prev.pats >= 0.1) %>% pull(feature),
                    'prev_20' = op.motus.stats %>%
                        filter(prev.all >= 0.2) %>% pull(feature),
                    'prev_pat_20' = op.motus.stats %>%
                        filter(prev.pats >= 0.2) %>% pull(feature),
                    'prev_50' = op.motus.stats %>%
                        filter(prev.all >= 0.5) %>% pull(feature),
                    'prev_pat_50' = op.motus.stats %>%
                        filter(prev.pats >= 0.5) %>% pull(feature),
                    'ra_e4' = op.motus.stats %>%
                        filter(mean.ra.all >= 1e-4) %>% pull(feature),
                    'ra_pat_e4' = op.motus.stats %>%
                        filter(mean.ra.pats >= 1e-4) %>% pull(feature))),
    tar_target(sdna.genus.cutoff.feat,
               list('prev_10' = sdna.genus.stats %>%
                        filter(prev.all >= 0.1) %>% pull(feature),
                    'prev_pat_10' = sdna.genus.stats %>%
                        filter(prev.pats >= 0.1) %>% pull(feature),
                    'prev_20' = sdna.genus.stats %>%
                        filter(prev.all >= 0.2) %>% pull(feature),
                    'prev_pat_20' = sdna.genus.stats %>%
                        filter(prev.pats >= 0.2) %>% pull(feature),
                    'prev_50' = sdna.genus.stats %>%
                        filter(prev.all >= 0.5) %>% pull(feature),
                    'prev_pat_50' = sdna.genus.stats %>%
                        filter(prev.pats >= 0.5) %>% pull(feature),
                    'ra_e4' = sdna.genus.stats %>%
                        filter(mean.ra.all >= 1e-4) %>% pull(feature),
                    'ra_pat_e4' = sdna.genus.stats %>%
                        filter(mean.ra.pats >= 1e-4) %>% pull(feature))),
    tar_target(op.genus.cutoff.feat,
               list('prev_10' = op.genus.stats %>%
                        filter(prev.all >= 0.1) %>% pull(feature),
                    'prev_pat_10' = op.genus.stats %>%
                        filter(prev.pats >= 0.1) %>% pull(feature),
                    'prev_20' = op.genus.stats %>%
                        filter(prev.all >= 0.2) %>% pull(feature),
                    'prev_pat_20' = op.genus.stats %>%
                        filter(prev.pats >= 0.2) %>% pull(feature),
                    'prev_50' = op.genus.stats %>%
                        filter(prev.all >= 0.5) %>% pull(feature),
                    'prev_pat_50' = op.genus.stats %>%
                        filter(prev.pats >= 0.5) %>% pull(feature),
                    'ra_e4' = op.genus.stats %>%
                        filter(mean.ra.all >= 1e-4) %>% pull(feature),
                    'ra_pat_e4' = op.genus.stats %>%
                        filter(mean.ra.pats >= 1e-4) %>% pull(feature))),
    
    ## PRINCIPAL COMPONENTS ANALYSIS: METABOLITES ----
    # https://cran.r-project.org/web/packages/ggfortify/vignettes/plot_pca.html
    # plasma metabolites only
    tar_target(p.met.imp, p.metabolites.wide %>%
                   filter.by.pct.nas(cov.thresh = 0.7) %>%
                   mutate(PatientID = str_split(UniqueID,'-') %>%
                              map_chr(~ head(.,1))) %>%
                   # impute w/ median of other samples from a given patient
                   nest(data = -PatientID) %>%
                   mutate(
                       data = map(
                           data, ~ tidyimpute::impute_median_if(., is.double))) %>%
                   unnest(data) %>%
                   filter.by.pct.nas(cov.thresh = 0.9) %>%
                   # then try to randomly sample from a metab. vector
                   tidyimpute::impute_all(na.tools::na.bootstrap) %>%
                   # finally, impute final NAs (WHY?) from med. of metab.vector
                   mutate_if(is.double, 
                             ~ gtools::na.replace(., median, na.rm=TRUE)) %>%
                   select(-PatientID) %>%
                   arrange(UniqueID)),
    tar_target(p.met.imp.annot, p.met.imp %>%
                   gather('R_Compound','Val',-UniqueID) %>%
                   mutate(R_Compound = str_remove_all(
                       R_Compound,'P_')) %>%
                   left_join(select(metabolite.categ, 
                                    R_Compound, ChemClass)) %>%
                   mutate(PatientID = str_split(UniqueID,'-') %>%
                              map_chr(~ head(.,1)))),
    # every metabolite w/ eigenvector -- no visual value
    # tar_target(p.met.pca, p.met.imp %>%
    #                column_to_rownames('UniqueID') %>%
    #                prcomp(center = TRUE, scale = TRUE) ),
    # group metabolites by chem class then perform PCA
    ## average across samples (per patient), sum across compounds (per class)
    tar_target(p.met.imp.bin, p.met.imp.annot %>%
                   group_by(PatientID, ChemClass, R_Compound) %>%
                   # average across samples from same individual
                   summarize(Val = mean(Val)) %>%
                   ungroup() %>%
                   group_by(PatientID,ChemClass) %>%
                   # sum/bin compounds of same group
                   summarize(Val = sum(Val))),
    tar_target(p.met.ccbin.pca, p.met.imp.bin %>%
                   spread('ChemClass','Val') %>%
                   column_to_rownames('PatientID') %>%
                   prcomp(center = TRUE, scale = TRUE)),
    tar_target(p.met.ccbin.biplot, plot.osci.biplot(p.met.ccbin.pca,
                                                    patient.data)),
    ## average across samples (per patient), average across compounds (per class)
    tar_target(p.met.imp.avg, p.met.imp.annot %>% 
                   group_by(PatientID, ChemClass,R_Compound) %>%
                   # average across samples from same individual - fm
                   mutate(Val = mean(Val)) %>%
                   ungroup() %>%
                   group_by(PatientID,ChemClass) %>% 
                   # sum/bin compounds of same group
                   summarize(Val = mean(Val))), 
    tar_target(p.met.ccavg.pca, p.met.imp.avg %>%
                   spread('ChemClass','Val') %>%
                   column_to_rownames('PatientID') %>%
                   prcomp(center = TRUE, scale = TRUE)), 
    tar_target(p.met.ccavg.biplot, plot.osci.biplot(p.met.ccavg.pca,
                                                    patient.data)),
    
    # urine metabolites only
    tar_target(u.met.imp, u.metabolites.wide %>%
                   filter.by.pct.nas(cov.thresh = 0.7) %>%
                   mutate(PatientID = str_split(UniqueID,'-') %>%
                              map_chr(~ head(.,1))) %>%
                   # impute w/ median of other samples from a given patient
                   nest(data = -PatientID) %>%
                   mutate(
                       data = map(
                           data, ~ tidyimpute::impute_median_if(., is.double))) %>%
                   unnest(data) %>%
                   filter.by.pct.nas(cov.thresh = 0.9) %>%
                   # then try to randomly sample from a metab. vector
                   tidyimpute::impute_all(na.tools::na.bootstrap) %>%
                   # finally, impute final NAs (WHY?) from med. of metab.vector
                   mutate_if(is.double, 
                             ~ gtools::na.replace(., median, na.rm=TRUE)) %>%
                   select(-PatientID) %>%
                   arrange(UniqueID)),
    tar_target(u.met.imp.annot, u.met.imp %>%
                   gather('R_Compound','Val',-UniqueID) %>%
                   mutate(R_Compound = str_remove_all(
                       R_Compound,'U_')) %>%
                   left_join(select(metabolite.categ, 
                                    R_Compound, ChemClass)) %>%
                   mutate(PatientID = str_split(UniqueID,'-') %>%
                              map_chr(~ head(.,1)))),
    tar_target(u.met.imp.bin, u.met.imp.annot %>%
                   group_by(PatientID, ChemClass, R_Compound) %>%
                   # average across samples from same individual
                   summarize(Val = mean(Val)) %>%
                   ungroup() %>%
                   group_by(PatientID,ChemClass) %>%
                   # sum/bin compounds of same group
                   summarize(Val = sum(Val))),
    tar_target(u.met.ccbin.pca, u.met.imp.bin %>%
                   spread('ChemClass','Val') %>%
                   column_to_rownames('PatientID') %>%
                   prcomp(center = TRUE, scale = TRUE)),
    tar_target(u.met.ccbin.biplot, plot.osci.biplot(u.met.ccbin.pca,
                                                    patient.data,
                                                    frame = FALSE)),
    
    ## CO-ASSOCIATION: MICROBIOME FEATURES ----
    tar_target(genus.only, sdna.genus.feat %>%
                   rename_if(is.double, ~ paste0('S_',.)) %>%
                   mutate(UniqueID = str_remove_all(UniqueID,'-SDNA')) %>%
                   full_join(op.genus.feat %>%
                                 rename_if(is.double,
                                           ~ paste0('O_',.)) %>%
                                 mutate(UniqueID = str_remove_all(
                                     UniqueID,'-OP'))) %>%
                   select_if(is.double)),
    tar_target(genus.corr.raw, genus.only %>%
                   corrr::correlate(method='spearman')),
    tar_target(genus.qs.raw, genus.only %>%
                   psych::corr.test(method='spearman',
                                    adjust='fdr',
                                    ci=FALSE) %>%
                   get('p',.) %>%
                   corrr::as_cordf()), 
    tar_target(genus.corr.filt, genus.qs.raw %>%
                   gather('term2','FDR',-term) %>%
                   filter(!is.na(FDR)) %>%
                   left_join(genus.corr.raw %>%
                                 gather('term2','rho',-term) %>%
                                 filter(!is.na(rho)),
                             by=c('term','term2'))),
    
    ## CO-ASSOCIATION: METABOLOME FEATURES ----
    tar_target(metabs.only, p.metabolites.wide %>%
                   left_join(u.metabolites.wide) %>%
                   select(-UniqueID)),
    tar_target(metab.corr.raw, metabs.only %>%
                   corrr::correlate(method='spearman')),
    tar_target(metab.qs.raw, metabs.only %>%
                   psych::corr.test(method='spearman',
                                    adjust='fdr',
                                    ci=FALSE) %>%
                   get('p',.) %>%
                   corrr::as_cordf()),
    tar_target(metab.corr.filt, metab.corr.raw %>% 
                   gather('term2','rho',-term) %>%
                   left_join(metab.qs.raw %>%
                                 gather('term2','FDR',-term),
                             by=c('term','term2')) %>%
                   filter(!is.na(rho) & !is.na(FDR)) %>%
                   filter(!str_detect(term,'_13C|_D3|_D5')) %>%
                   filter(!str_detect(term2,'_13C|_D5|_D3')) %>%
                   mutate(term1.name = str_remove_all(term,'P_')) %>%
                   mutate(term1.name = str_remove_all(term1.name,'U_')) %>%
                   mutate(term2.name = str_remove_all(term2,'P_')) %>%
                   mutate(term2.name = str_remove_all(term2.name,'U_')) %>%
                   left_join(select(metabolite.categ,R_Compound,ChemClass) %>%
                                 magrittr::set_names(
                                     c('term1.name','term1.class'))) %>%
                   left_join(select(metabolite.categ,R_Compound,ChemClass) %>%
                                 magrittr::set_names(
                                     c('term2.name','term2.class')))),
    tar_target(tryp.metab.corrs, metab.corr.filt %>%
                   filter(str_detect(term, 'rypt|ynur|onin|ndol'))),
    
    ## CONNECTING TRYP PATHWAYS: METABOLOME ----
    tar_target(p.tryp.select, metabolite.categ %>% 
                   filter(str_detect(Compound,'ryp|ynur|ndol|tonin')) %>%
                   # remove duplicate measurements
                   filter(Measurement == 'LC-MS') %>%
                   select(R_Compound, Compound) %>%
                   pull(R_Compound) %>%
                   paste0('P_', .)), 
    tar_target(p.tryp.metabs, p.metabolites.wide %>%
                   select(UniqueID, 
                          p.tryp.select[which(p.tryp.select%in%colnames(
                              p.metabolites.wide))] )),
    tar_target(u.tryp.select, metabolite.categ %>% 
                   filter(str_detect(Compound,'ryp|ynur|ndol|tonin')) %>%
                   # remove duplicate measurements
                   filter(Measurement == 'LC-MS') %>%
                   select(R_Compound, Compound) %>%
                   pull(R_Compound) %>%
                   paste0('U_', .)), 
    tar_target(u.tryp.metabs, u.metabolites.wide %>%
                   select(UniqueID, 
                          u.tryp.select[which(u.tryp.select%in%colnames(
                              u.metabolites.wide))] )),
    # annotated and long format
    tar_target(tryp.metabs, p.tryp.metabs %>%
                   full_join(u.tryp.metabs) %>%
                   gather('Metabolite','Val',-UniqueID) %>%
                   mutate(Sample = str_split(Metabolite, '_') %>%
                              map_chr(~ head(., 1))) %>%
                   mutate(Metabolite = str_remove_all(Metabolite, 'P_')) %>%
                   mutate(Metabolite = str_remove_all(Metabolite, 'U_')) %>%
                   left_join(select(metabolite.categ, 
                                    R_Compound, Compound), 
                             by = c('Metabolite'='R_Compound')) %>%
                   left_join(clinical.feat %>% 
                                 filter(Site%in%c('P','U')) %>% 
                                 select(GenID, PatientID, Days, 
                                        contains('OSCI')) %>% 
                                 distinct(), by=c('UniqueID'='GenID')) %>%
                   mutate(Sample = ifelse(Sample=='P','Plasma','Urine'))),
    
    ###################################################################
    ### METADECONFOUNDR RUNS: TAXONOMIC TABLES ----
    # motus -- all SDNA samples -- status = COVID ----
    tar_target(mdc.sdna.motus.all.feat, sdna.feat.ra %>%
                   select(c(1:4), all_of(
                       intersect(sdna.motus.cutoff.feat$prev_20,
                                 sdna.motus.cutoff.feat$ra_e4))) %>%
                   left_join(sdna.alpha) %>%
                   select(all_of(colnames(sdna.alpha)),
                          everything(),
                          -c(Site, UniqueID, Days)) %>%
                   arrange(GenID) %>%
                   remove.zero.cols() %>%
                   column_to_rownames('GenID')),
    tar_target(mdc.sdna.motus.all.meta, sdna.clinical.feat %>%
                   select(PatientID, GenID, Days, OSCI_Score_Samp,
                          Antibiotics, contains('VL')) %>%
                   left_join(select(patient.data, 
                                    -PatientID, -OSCI_Class, -GISymptoms),
                             by=c('PatientID'='ShortID')) %>%
                   dplyr::rename(OSCI_Score_Worst = 'OSCI_Score') %>%
                   mutate(Status = ifelse(str_detect(PatientID,'K'),0,1)) %>%
                   remove.sparse.vars(keep = c('Antibiotics','DaysHospitalized',
                                               'Pneumonia','Bacteremia'), 
                                      na.thresh = params$na.thresh.glob) %>%
                   discretize.days() %>%
                   prep.metadata() %>%
                   mutate(Sex = ifelse(Sex==0,1,0)) %>%
                   filter(GenID %in% rownames(mdc.sdna.motus.all.feat)) %>%
                   select(Status, everything()) %>%
                   arrange(GenID) %>%
                   column_to_rownames('GenID')),
    tar_target(mdc.sdna.motus.all.res,
               metadeconfoundR::MetaDeconfound(mdc.sdna.motus.all.feat,
                                               mdc.sdna.motus.all.meta,
                                               deconfT = c('Antibiotics'),
                                               randomVar = list('+ (1|PatientID)',
                                                                c('PatientID')),
                                               nnodes = params$nnodes,
                                               robustCutoff = 0,
                                               logfile = here('logfiles','log.txt')) %>%
                   map(~ rownames_to_column(as.data.frame(.), 'feature')) %>%
                   map(~ as_tibble(.)) %>%
                   tibble.results(., 'tbl.all', 'COVID') %>%
                   mutate(Confounding = str_replace_all(
                       Confounding,'Status','COVID')) %>%
                   filter(Covariate != 'PatientID')),
    
    # motus -- all OP samples -- status = COVID ----
    tar_target(mdc.op.motus.all.feat, op.feat.ra %>%
                   select(c(1:4), all_of(
                       intersect(op.motus.cutoff.feat$prev_20,
                                 op.motus.cutoff.feat$ra_e4))) %>%
                   left_join(op.alpha) %>%
                   select(all_of(colnames(op.alpha)),
                          everything(),
                          -c(Site, UniqueID, Days)) %>%
                   arrange(GenID) %>%
                   remove.zero.cols() %>%
                   column_to_rownames('GenID')),
    tar_target(mdc.op.motus.all.meta, op.clinical.feat %>%
                   select(PatientID, GenID, Days, OSCI_Score_Samp,
                          Antibiotics, contains('VL')) %>%
                   left_join(select(patient.data, 
                                    -PatientID, -OSCI_Class, -GISymptoms),
                             by=c('PatientID'='ShortID')) %>%
                   dplyr::rename(OSCI_Score_Worst = 'OSCI_Score') %>%
                   mutate(Status = ifelse(str_detect(PatientID,'K'),0,1)) %>%
                   remove.sparse.vars(keep = c('Antibiotics','DaysHospitalized',
                                               'Pneumonia','Bacteremia'), 
                                      na.thresh = params$na.thresh.glob) %>%
                   discretize.days() %>%
                   prep.metadata() %>%
                   mutate(Sex = ifelse(Sex==0,1,0)) %>%
                   filter(GenID %in% rownames(mdc.op.motus.all.feat)) %>%
                   select(Status, everything()) %>%
                   arrange(GenID) %>%
                   column_to_rownames('GenID')),
    tar_target(mdc.op.motus.all.res,
               metadeconfoundR::MetaDeconfound(mdc.op.motus.all.feat,
                                               mdc.op.motus.all.meta,
                                               deconfT = c('Antibiotics'),
                                               randomVar = list('+ (1|PatientID)',
                                                                c('PatientID')),
                                               nnodes = params$nnodes,
                                               robustCutoff = 0,
                                               logfile = here('logfiles','log.txt')) %>%
                   map(~ rownames_to_column(as.data.frame(.), 'feature')) %>%
                   map(~ as_tibble(.)) %>%
                   tibble.results(., 'tbl.all', 'COVID') %>%
                   mutate(Confounding = str_replace_all(
                       Confounding,'Status','COVID')) %>%
                   filter(Covariate != 'PatientID')),
    
    
    # motus -- patients only -- SDNA -- status = AbxCurr ----
    tar_target(mdc.sdna.motus.pat.feat,
               mdc.sdna.motus.all.feat[which(!str_detect(
                   rownames(mdc.sdna.motus.all.feat),'K')),
                   which(colnames(mdc.sdna.motus.all.feat) %in% c(
                       'Richness','Shannon','InvSimp', intersect(
                           sdna.motus.cutoff.feat$prev_pat_20,
                           sdna.motus.cutoff.feat$ra_pat_e4)) )]),
    tar_target(mdc.sdna.motus.pat.meta, sdna.clinical.feat %>%
                   left_join(select(patient.data, 
                                    -PatientID, -OSCI_Class,),
                             by=c('PatientID'='ShortID')) %>%
                   dplyr::rename(OSCI_Score_Worst = 'OSCI_Score') %>%
                   mutate(Status = Antibiotics) %>%
                   remove.sparse.vars(keep = c('DaysHospitalized','Pneumonia',
                                               'Bacteremia','VAP','Sepsis','Status'), 
                                      na.thresh = params$na.thresh.glob) %>%
                   discretize.days() %>%
                   prep.metadata() %>%
                   mutate(Sex = ifelse(Sex==0,1,0)) %>%
                   filter(GenID %in% rownames(mdc.sdna.motus.pat.feat)) %>%
                   select(Status, everything(), -UniqueID, -contains('OSCI_Class')) %>%
                   arrange(GenID) %>%
                   column_to_rownames('GenID')),
    tar_target(mdc.sdna.motus.pat.res,
               metadeconfoundR::MetaDeconfound(mdc.sdna.motus.pat.feat,
                                               mdc.sdna.motus.pat.meta,
                                               deconfT = c('Pneumonia'),
                                               randomVar = list('+ (1|PatientID)',
                                                                c('PatientID')),
                                               nnodes = params$nnodes,
                                               robustCutoff = 0,
                                               logfile = here('logfiles','log.txt')) %>%
                   map(~ rownames_to_column(as.data.frame(.), 'feature')) %>%
                   map(~ as_tibble(.)) %>%
                   tibble.results(., 'tbl.all', 'AbxCurr') %>%
                   mutate(Confounding = str_replace_all(
                       Confounding,'Status','AbxCurr')) %>%
                   filter(Covariate != 'PatientID')),
    
    # motus -- patients only -- OP -- status = AbxCurr ----
    tar_target(mdc.op.motus.pat.feat,
               mdc.op.motus.all.feat[which(!str_detect(
                   rownames(mdc.op.motus.all.feat),'K')),
                   which(colnames(mdc.op.motus.all.feat) %in% c(
                       'Richness','Shannon','InvSimp', intersect(
                           op.motus.cutoff.feat$prev_pat_20,
                           op.motus.cutoff.feat$ra_pat_e4)) )]),
    tar_target(mdc.op.motus.pat.meta, op.clinical.feat %>%
                   left_join(select(patient.data, 
                                    -PatientID, -OSCI_Class,),
                             by=c('PatientID'='ShortID')) %>%
                   dplyr::rename(OSCI_Score_Worst = 'OSCI_Score') %>%
                   mutate(Status = Antibiotics) %>%
                   remove.sparse.vars(keep = c('DaysHospitalized','Pneumonia',
                                               'Bacteremia','VAP','Sepsis','Status'), 
                                      na.thresh = params$na.thresh.glob) %>%
                   discretize.days() %>%
                   prep.metadata() %>%
                   mutate(Sex = ifelse(Sex==0,1,0)) %>%
                   filter(GenID %in% rownames(mdc.op.motus.pat.feat)) %>%
                   select(Status, everything(), -UniqueID, -contains('OSCI_Class')) %>%
                   arrange(GenID) %>%
                   column_to_rownames('GenID')),
    tar_target(mdc.op.motus.pat.res,
               metadeconfoundR::MetaDeconfound(mdc.op.motus.pat.feat,
                                               mdc.op.motus.pat.meta,
                                               deconfT = c('Pneumonia'),
                                               randomVar = list('+ (1|PatientID)',
                                                                c('PatientID')),
                                               nnodes = params$nnodes,
                                               robustCutoff = 0,
                                               logfile = here('logfiles','log.txt')) %>%
                   map(~ rownames_to_column(as.data.frame(.), 'feature')) %>%
                   map(~ as_tibble(.)) %>%
                   tibble.results(., 'tbl.all', 'AbxCurr') %>%
                   mutate(Confounding = str_replace_all(
                       Confounding,'Status','AbxCurr')) %>%
                   filter(Covariate != 'PatientID')),
    
    # genus -- all SDNA samples -- status = COVID ----
    tar_target(mdc.sdna.genus.all.feat, sdna.genus.feat %>%
                   select(c(1:2), all_of(
                       intersect(sdna.genus.cutoff.feat$prev_10,
                                 sdna.genus.cutoff.feat$ra_e4))) %>%
                   left_join(sdna.alpha) %>%
                   select(all_of(colnames(sdna.alpha)),
                          everything(),
                          -Days) %>%
                   arrange(UniqueID) %>%
                   remove.zero.cols() %>%
                   column_to_rownames('UniqueID')),
    tar_target(mdc.sdna.genus.all.meta, sdna.clinical.feat %>%
                   select(PatientID, UniqueID, Days, OSCI_Score_Samp,
                          Antibiotics, contains('VL')) %>%
                   left_join(select(patient.data, 
                                    -PatientID, -OSCI_Class, -GISymptoms),
                             by = c('PatientID'='ShortID')) %>%
                   dplyr::rename(OSCI_Score_Worst = 'OSCI_Score') %>%
                   mutate(Status = ifelse(str_detect(PatientID,'K'),0,1)) %>%
                   remove.sparse.vars(keep = c('Antibiotics','DaysHospitalized',
                                               'Pneumonia','Bacteremia'), 
                                      na.thresh = params$na.thresh.glob) %>%
                   discretize.days() %>%
                   prep.metadata() %>%
                   mutate(Sex = ifelse(Sex==0,1,0)) %>%
                   filter(UniqueID %in% rownames(mdc.sdna.genus.all.feat)) %>%
                   select(Status, everything()) %>%
                   arrange(UniqueID) %>%
                   column_to_rownames('UniqueID')),
    tar_target(mdc.sdna.genus.all.res,
               metadeconfoundR::MetaDeconfound(mdc.sdna.genus.all.feat,
                                               mdc.sdna.genus.all.meta,
                                               deconfT = c('Antibiotics'),
                                               randomVar = list('+ (1|PatientID)',
                                                                c('PatientID')),
                                               nnodes = params$nnodes,
                                               robustCutoff = 0,
                                               logfile = here('logfiles','log.txt')) %>%
                   map(~ rownames_to_column(as.data.frame(.), 'feature')) %>%
                   map(~ as_tibble(.)) %>%
                   tibble.results(., 'tbl.all', 'COVID') %>%
                   mutate(Confounding = str_replace_all(
                       Confounding,'Status','COVID')) %>%
                   filter(Covariate != 'PatientID')),
    
    # genus -- all SDNA samples w/ plasma metabolites -- status = COVID ----
    tar_target(mdc.sdna.genus.all.metab.meta, sdna.clinical.feat %>%
                   select(PatientID, UniqueID, Days, OSCI_Score_Samp,
                          Antibiotics, contains('VL')) %>%
                   left_join(select(patient.data, 
                                    -PatientID, -OSCI_Class, -GISymptoms),
                             by = c('PatientID'='ShortID')) %>%
                   mutate(Visit = str_split(UniqueID, '-') %>%
                              map_chr(~ tail(., 1))) %>%
                   mutate(TmpID = paste0(PatientID, '-', Visit)) %>%
                   select(TmpID, UniqueID, Visit, everything()) %>%
                   left_join(p.metabolites.wide.full, 
                             by = c('TmpID'='UniqueID')) %>%
                   select(-c(TmpID, Visit)) %>%
                   dplyr::rename(OSCI_Score_Worst = 'OSCI_Score') %>%
                   mutate(Status = ifelse(str_detect(PatientID,'K'),0,1)) %>%
                   remove.sparse.vars(keep = c('Antibiotics','DaysHospitalized',
                                               'Pneumonia','Bacteremia'), 
                                      na.thresh = params$na.thresh.glob) %>%
                   discretize.days() %>%
                   prep.metadata() %>%
                   mutate(Sex = ifelse(Sex==0,1,0)) %>%
                   dplyr::filter(UniqueID %in% rownames(mdc.sdna.genus.all.feat)) %>%
                   select(Status, everything()) %>%
                   arrange(UniqueID) %>%
                   column_to_rownames('UniqueID')),
    tar_target(mdc.sdna.genus.all.metab.res,
               metadeconfoundR::MetaDeconfound(mdc.sdna.genus.all.feat,
                                               mdc.sdna.genus.all.metab.meta,
                                               deconfT = c('Antibiotics'),
                                               randomVar = list('+ (1|PatientID)',
                                                                c('PatientID')),
                                               nnodes = params$nnodes,
                                               robustCutoff = 0,
                                               logfile = here('logfiles','log.txt')) %>%
                   map(~ rownames_to_column(as.data.frame(.), 'feature')) %>%
                   map(~ as_tibble(.)) %>%
                   tibble.results(., 'tbl.all', 'COVID') %>%
                   mutate(Confounding = str_replace_all(
                       Confounding,'Status','COVID')) %>%
                   dplyr::filter(Covariate != 'PatientID')),
    
    # genus -- all OP samples -- status = COVID ----
    tar_target(mdc.op.genus.all.feat, op.genus.feat %>%
                   select(c(1:2), all_of(
                       intersect(op.genus.cutoff.feat$prev_10,
                                 op.genus.cutoff.feat$ra_e4))) %>%
                   left_join(op.alpha) %>%
                   select(all_of(colnames(op.alpha)),
                          everything(),
                          -Days) %>%
                   arrange(UniqueID) %>%
                   remove.zero.cols() %>%
                   column_to_rownames('UniqueID')),
    tar_target(mdc.op.genus.all.meta, op.clinical.feat %>%
                   select(PatientID, UniqueID, Days, OSCI_Score_Samp,
                          Antibiotics, contains('VL')) %>%
                   left_join(select(patient.data, 
                                    -PatientID, -OSCI_Class, -GISymptoms),
                             by = c('PatientID'='ShortID')) %>%
                   dplyr::rename(OSCI_Score_Worst = 'OSCI_Score') %>%
                   mutate(Status = ifelse(str_detect(PatientID,'K'),0,1)) %>%
                   remove.sparse.vars(keep = c('Antibiotics','DaysHospitalized',
                                               'Pneumonia','Bacteremia'), 
                                      na.thresh = params$na.thresh.glob) %>%
                   discretize.days() %>%
                   prep.metadata() %>%
                   mutate(Sex = ifelse(Sex==0,1,0)) %>%
                   filter(UniqueID %in% rownames(mdc.op.genus.all.feat)) %>%
                   select(Status, everything()) %>%
                   arrange(UniqueID) %>%
                   column_to_rownames('UniqueID')),
    tar_target(mdc.op.genus.all.res,
               metadeconfoundR::MetaDeconfound(mdc.op.genus.all.feat,
                                               mdc.op.genus.all.meta,
                                               deconfT = c('Antibiotics'),
                                               randomVar = list('+ (1|PatientID)',
                                                                c('PatientID')),
                                               nnodes = params$nnodes,
                                               robustCutoff = 0,
                                               logfile = here('logfiles','log.txt')) %>%
                   map(~ rownames_to_column(as.data.frame(.), 'feature')) %>%
                   map(~ as_tibble(.)) %>%
                   tibble.results(., 'tbl.all', 'COVID') %>%
                   mutate(Confounding = str_replace_all(
                       Confounding,'Status','COVID')) %>%
                   filter(Covariate != 'PatientID')),
    
    # genus -- all OP samples w/ plasma metabolites -- status = COVID ----
    tar_target(mdc.op.genus.all.metab.meta, op.clinical.feat %>%
                   select(PatientID, UniqueID, Days, OSCI_Score_Samp,
                          Antibiotics, contains('VL')) %>%
                   left_join(select(patient.data, 
                                    -PatientID, -OSCI_Class, -GISymptoms),
                             by = c('PatientID'='ShortID')) %>%
                   mutate(Visit = str_split(UniqueID, '-') %>%
                              map_chr(~ tail(., 1))) %>%
                   mutate(TmpID = paste0(PatientID, '-', Visit)) %>%
                   select(TmpID, UniqueID, Visit, everything()) %>%
                   left_join(p.metabolites.wide.full, 
                             by = c('TmpID'='UniqueID')) %>%
                   select(-c(TmpID, Visit)) %>%
                   dplyr::rename(OSCI_Score_Worst = 'OSCI_Score') %>%
                   mutate(Status = ifelse(str_detect(PatientID,'K'),0,1)) %>%
                   remove.sparse.vars(keep = c('Antibiotics','DaysHospitalized',
                                               'Pneumonia','Bacteremia'), 
                                      na.thresh = params$na.thresh.glob) %>%
                   discretize.days() %>%
                   prep.metadata() %>%
                   mutate(Sex = ifelse(Sex==0,1,0)) %>%
                   dplyr::filter(UniqueID %in% rownames(mdc.op.genus.all.feat)) %>%
                   select(Status, everything()) %>%
                   arrange(UniqueID) %>%
                   column_to_rownames('UniqueID')),
    tar_target(mdc.op.genus.all.metab.res,
               metadeconfoundR::MetaDeconfound(mdc.op.genus.all.feat,
                                               mdc.op.genus.all.metab.meta,
                                               deconfT = c('Antibiotics'),
                                               randomVar = list('+ (1|PatientID)',
                                                                c('PatientID')),
                                               nnodes = params$nnodes,
                                               robustCutoff = 0,
                                               logfile = here('logfiles','log.txt')) %>%
                   map(~ rownames_to_column(as.data.frame(.), 'feature')) %>%
                   map(~ as_tibble(.)) %>%
                   tibble.results(., 'tbl.all', 'COVID') %>%
                   mutate(Confounding = str_replace_all(
                       Confounding,'Status','COVID')) %>%
                   dplyr::filter(Covariate != 'PatientID')),
    
    # genus -- patients only -- SDNA -- status = AbxCurr ----
    tar_target(mdc.sdna.genus.pat.feat,
               mdc.sdna.genus.all.feat[which(!str_detect(
                   rownames(mdc.sdna.genus.all.feat),'K')),
                   which(colnames(mdc.sdna.genus.all.feat) %in% c(
                       'Richness','Shannon','InvSimp', intersect(
                           sdna.genus.cutoff.feat$prev_pat_10,
                           sdna.genus.cutoff.feat$ra_pat_e4)) )]),
    tar_target(mdc.sdna.genus.pat.meta, sdna.clinical.feat %>%
                   filter(!str_detect(PatientID,'K')) %>%
                   left_join(select(patient.data, -PatientID, -OSCI_Class),
                             by = c('PatientID'='ShortID')) %>%
                   dplyr::rename(OSCI_Score_Worst = 'OSCI_Score') %>%
                   mutate(Status = Antibiotics) %>%
                   select(-contains('OSCI_Class'), -GenID, -Antibiotics) %>%
                   remove.sparse.vars(keep = c('Status','DaysHospitalized','VAP',
                                               'Pneumonia','Bacteremia','Sepsis'), 
                                      na.thresh = params$na.thresh.glob) %>%
                   discretize.days() %>%
                   prep.metadata() %>%
                   mutate(Sex = ifelse(Sex==0,1,0)) %>%
                   filter(UniqueID %in% rownames(mdc.sdna.genus.pat.feat)) %>%
                   select(Status, everything()) %>%
                   arrange(UniqueID) %>%
                   column_to_rownames('UniqueID')),
    tar_target(mdc.sdna.genus.pat.res,
               metadeconfoundR::MetaDeconfound(mdc.sdna.genus.pat.feat,
                                               mdc.sdna.genus.pat.meta,
                                               randomVar = list('+ (1|PatientID)',
                                                                c('PatientID')),
                                               nnodes = params$nnodes,
                                               robustCutoff = 0,
                                               logfile = here('logfiles','log.txt')) %>%
                   map(~ rownames_to_column(as.data.frame(.), 'feature')) %>%
                   map(~ as_tibble(.)) %>%
                   tibble.results(., 'tbl.all', 'AbxCurr') %>%
                   mutate(Confounding = str_replace_all(
                       Confounding,'Status','AbxCurr')) %>%
                   filter(Covariate != 'PatientID')),
    
    # genus -- patients only -- OP -- status = AbxCurr ----
    tar_target(mdc.op.genus.pat.feat,
               mdc.op.genus.all.feat[which(!str_detect(
                   rownames(mdc.op.genus.all.feat),'K')),
                   which(colnames(mdc.op.genus.all.feat) %in% c(
                       'Richness','Shannon','InvSimp', intersect(
                           op.genus.cutoff.feat$prev_pat_10,
                           op.genus.cutoff.feat$ra_pat_e4)) )]),
    tar_target(mdc.op.genus.pat.meta, op.clinical.feat %>%
                   filter(!str_detect(PatientID,'K')) %>%
                   left_join(select(patient.data, -PatientID, -OSCI_Class),
                             by = c('PatientID'='ShortID')) %>%
                   dplyr::rename(OSCI_Score_Worst = 'OSCI_Score') %>%
                   mutate(Status = Antibiotics) %>%
                   select(-contains('OSCI_Class'), -GenID, -Antibiotics) %>%
                   remove.sparse.vars(keep = c('Status','DaysHospitalized','VAP',
                                               'Pneumonia','Bacteremia','Sepsis'), 
                                      na.thresh = params$na.thresh.glob) %>%
                   discretize.days() %>%
                   prep.metadata() %>%
                   mutate(Sex = ifelse(Sex==0,1,0)) %>%
                   filter(UniqueID %in% rownames(mdc.op.genus.pat.feat)) %>%
                   select(Status, everything()) %>%
                   arrange(UniqueID) %>%
                   column_to_rownames('UniqueID')),
    tar_target(mdc.op.genus.pat.res,
               metadeconfoundR::MetaDeconfound(mdc.op.genus.pat.feat,
                                               mdc.op.genus.pat.meta,
                                               randomVar = list('+ (1|PatientID)',
                                                                c('PatientID')),
                                               nnodes = params$nnodes,
                                               robustCutoff = 0,
                                               logfile = here('logfiles','log.txt')) %>%
                   map(~ rownames_to_column(as.data.frame(.), 'feature')) %>%
                   map(~ as_tibble(.)) %>%
                   tibble.results(., 'tbl.all', 'AbxCurr') %>%
                   mutate(Confounding = str_replace_all(Confounding,'Status','AbxCurr')) %>%
                   filter(Covariate != 'PatientID')),
    
    # phylum -- all SDNA samples -- status = COVID ----
    tar_target(mdc.sdna.phy.all.feat, sdna.ps.phy@otu_table %>%
                   as.data.frame() %>%
                   magrittr::set_colnames(
                       as.data.frame(sdna.ps.phy@tax_table) %>% 
                           pull(Phylum))),
    tar_target(mdc.sdna.phy.all.meta, meta.from.phyloseq(sdna.ps.phy) %>%
                   dplyr::rename(OSCI_Score_Worst = 'OSCI_Score',
                                 AbxCurr = 'Antibiotics', HAP = 'Pneumonia') %>%
                   mutate(Status = ifelse(str_detect(PatientID,'K'),0,1)) %>%
                   remove.sparse.vars(keep = c('AbxCurr','Status','OSCI_Score_Worst',
                                               'OSCI_Score_Samp','HAP','Days'), 
                                      na.thresh = params$na.thresh.glob) %>%
                   discretize.days() %>%
                   prep.metadata() %>%
                   mutate(Sex = ifelse(Sex==0,1,0)) %>%
                   select(Status, everything(), -contains('OSCI_Class')) %>%
                   column_to_rownames('UniqueID')),
    tar_target(mdc.sdna.phy.all.res,
               metadeconfoundR::MetaDeconfound(mdc.sdna.phy.all.feat,
                                               mdc.sdna.phy.all.meta,
                                               deconfT = c('AbxCurr'),
                                               randomVar = list('+ (1|PatientID)',
                                                                c('PatientID')),
                                               nnodes = params$nnodes,
                                               robustCutoff = 0,
                                               logfile = here('logfiles','log.txt')) %>%
                   map(~ rownames_to_column(as.data.frame(.), 'feature')) %>%
                   map(~ as_tibble(.)) %>%
                   tibble.results(., 'tbl.all', 'COVID') %>%
                   mutate(Confounding = str_replace_all(
                       Confounding,'Status','COVID')) %>%
                   filter(Covariate != 'PatientID')),
    
    # phylum -- all OP samples -- status = COVID ----
    tar_target(mdc.op.phy.all.feat, op.ps.phy@otu_table %>%
                   as.data.frame() %>%
                   magrittr::set_colnames(
                       as.data.frame(op.ps.phy@tax_table) %>% 
                           pull(Phylum))),
    tar_target(mdc.op.phy.all.meta, meta.from.phyloseq(op.ps.phy) %>%
                   dplyr::rename(OSCI_Score_Worst = 'OSCI_Score',
                                 AbxCurr = 'Antibiotics', HAP = 'Pneumonia') %>%
                   mutate(Status = ifelse(str_detect(PatientID,'K'),0,1)) %>%
                   remove.sparse.vars(keep = c('AbxCurr','Status','OSCI_Score_Worst',
                                               'OSCI_Score_Samp','HAP','Days'), 
                                      na.thresh = params$na.thresh.glob) %>%
                   discretize.days() %>%
                   prep.metadata() %>%
                   select(Status, everything(), -contains('OSCI_Class')) %>%
                   column_to_rownames('UniqueID')),
    tar_target(mdc.op.phy.all.res,
               metadeconfoundR::MetaDeconfound(mdc.op.phy.all.feat,
                                               mdc.op.phy.all.meta,
                                               deconfT = c('AbxCurr'),
                                               randomVar = list('+ (1|PatientID)',
                                                                c('PatientID')),
                                               nnodes = params$nnodes,
                                               robustCutoff = 0,
                                               logfile = here('logfiles','log.txt')) %>%
                   map(~ rownames_to_column(as.data.frame(.), 'feature')) %>%
                   map(~ as_tibble(.)) %>%
                   tibble.results(., 'tbl.all', 'COVID') %>%
                   mutate(Confounding = str_replace_all(
                       Confounding,'Status','COVID')) %>%
                   filter(Covariate != 'PatientID')),
    
    ### METADECONFOUNDR RUNS: FUNCTIONAL TABLES ----
    # KEGG + GMM modules -- all SDNA samples -- status = COVID ----
    tar_target(mdc.sdna.mods.all.feat, sdna.kegg %>%
                   select(-SampleID, -Site, -Days) %>%
                   mutate(UniqueID = str_remove_all(UniqueID,'-SDNA')) %>%
                   left_join(select(sdna.gmm, -Days)) %>%
                   near.zero.var.filter() %>%
                   remove.zero.cols() %>%
                   mutate(Status = ifelse(str_detect(UniqueID,'K'),0,1)) %>%
                   remove.sparse.vars(keep = c(),
                                      na.thresh = params$na.thresh.glob) %>%
                   select(-Status) %>%
                   arrange(UniqueID) %>%
                   column_to_rownames('UniqueID')),
    tar_target(mdc.sdna.mods.all.meta, sdna.clinical.feat %>%
                   select(PatientID, GenID, Days, OSCI_Score_Samp,
                          Antibiotics, contains('VL')) %>%
                   left_join(select(patient.data, 
                                    -PatientID, -OSCI_Class, -GISymptoms),
                             by=c('PatientID'='ShortID')) %>%
                   dplyr::rename(OSCI_Score_Worst = 'OSCI_Score') %>%
                   mutate(Status = ifelse(str_detect(PatientID,'K'),0,1)) %>%
                   remove.sparse.vars(keep = c('Antibiotics','DaysHospitalized',
                                               'Pneumonia','Bacteremia'), 
                                      na.thresh = params$na.thresh.glob) %>%
                   discretize.days() %>%
                   prep.metadata() %>%
                   mutate(Sex = ifelse(Sex==0,1,0)) %>%
                   filter(GenID %in% rownames(mdc.sdna.mods.all.feat)) %>%
                   select(Status, everything()) %>%
                   arrange(GenID) %>%
                   column_to_rownames('GenID')),
    tar_target(mdc.sdna.mods.all.res,
               metadeconfoundR::MetaDeconfound(mdc.sdna.mods.all.feat,
                                               mdc.sdna.mods.all.meta,
                                               deconfT = c('Antibiotics'),
                                               randomVar = list('+ (1|PatientID)',
                                                                c('PatientID')),
                                               nnodes = params$nnodes,
                                               robustCutoff = 0,
                                               logfile = here('logfiles','log.txt')) %>%
                   map(~ rownames_to_column(as.data.frame(.), 'feature')) %>%
                   map(~ as_tibble(.)) %>%
                   tibble.results(., 'tbl.all', 'COVID') %>%
                   mutate(Confounding = str_replace_all(
                       Confounding,'Status','COVID')) %>%
                   filter(Covariate != 'PatientID')),
    
    # KEGG + GMM modules -- patients only -- status = AbxCurr ----
    tar_target(mdc.sdna.mods.pat.feat, mdc.sdna.mods.all.feat %>%
                   filter(!str_detect(rownames(.), 'K'))),
    tar_target(mdc.sdna.mods.pat.meta, mdc.sdna.mods.all.meta %>%
                   filter(!str_detect(rownames(.), 'K'))),
    tar_target(mdc.sdna.mods.pat.res,
               metadeconfoundR::MetaDeconfound(mdc.sdna.mods.pat.feat,
                                               mdc.sdna.mods.pat.meta,
                                               randomVar = list('+ (1|PatientID)',
                                                                c('PatientID')),
                                               nnodes = params$nnodes,
                                               robustCutoff = 0,
                                               logfile = here('logfiles','log.txt')) %>%
                   map(~ rownames_to_column(as.data.frame(.), 'feature')) %>%
                   map(~ as_tibble(.)) %>%
                   tibble.results(., 'tbl.all', 'AbxCurr') %>%
                   mutate(Confounding = str_replace_all(
                       Confounding,'Status','AbxCurr')) %>%
                   filter(Covariate != 'PatientID')),
    # KEGG + GMM modules -- all SDNA samples w/ bin mOTUs -- status = COVID ----
    tar_target(mdc.sdna.mods.motus.all.meta, sdna.clinical.feat %>%
                   select(PatientID, GenID, Days, OSCI_Score_Samp,
                          Antibiotics, contains('VL')) %>%
                   left_join(select(patient.data, 
                                    -PatientID, -OSCI_Class, -GISymptoms),
                             by=c('PatientID'='ShortID')) %>%
                   dplyr::rename(OSCI_Score_Worst = 'OSCI_Score') %>%
                   mutate(Status = ifelse(str_detect(PatientID,'K'),0,1)) %>%
                   remove.sparse.vars(keep = c('Antibiotics','DaysHospitalized',
                                               'Pneumonia','Bacteremia'), 
                                      na.thresh = params$na.thresh.glob) %>%
                   discretize.days() %>%
                   prep.metadata() %>%
                   mutate(Sex = ifelse(Sex==0,1,0)) %>%
                   left_join(all.alphas) %>%
                   left_join(osci.all.bin.motus) %>%
                   filter(GenID %in% rownames(mdc.sdna.mods.all.feat)) %>%
                   select(Status, everything()) %>%
                   arrange(GenID) %>%
                   column_to_rownames('GenID')),
    tar_target(mdc.sdna.mods.motus.all.res,
               metadeconfoundR::MetaDeconfound(mdc.sdna.mods.all.feat,
                                               mdc.sdna.mods.motus.all.meta,
                                               deconfT = c('Antibiotics'),
                                               randomVar = list('+ (1|PatientID)',
                                                                c('PatientID')),
                                               nnodes = params$nnodes,
                                               robustCutoff = 0,
                                               logfile = here('logfiles','log.txt')) %>%
                   map(~ rownames_to_column(as.data.frame(.), 'feature')) %>%
                   map(~ as_tibble(.)) %>%
                   tibble.results(., 'tbl.all', 'COVID') %>%
                   mutate(Confounding = str_replace_all(
                       Confounding,'Status','COVID')) %>%
                   filter(Covariate != 'PatientID')),
    
    # plasma & urine -- all SDNA samples -- status = COVID -----
    
    ### METADECONFOUNDR RUNS: METABOLITES -----
    # plasma & urine -- all timepoints -- status = COVID ----
    tar_target(mdc.metabs.all.feat, p.metabolites.wide %>%
                   left_join(u.metabolites.wide) %>%
                   filter(!str_detect(UniqueID, fixed('.2'))) %>%
                   arrange(UniqueID) %>%
                   column_to_rownames('UniqueID') %>%
                   remove.zero.cols()),
    tar_target(mdc.metabs.all.meta, clinical.feat %>%
                   filter(Site %in% c('P','U')) %>%
                   select(GenID, PatientID, Days, Antibiotics, 
                          all_of(contains('OSCI'))) %>%
                   distinct() %>%
                   filter(GenID %in% rownames(mdc.metabs.all.feat)) %>%
                   discretize.days() %>%
                   mutate(Days = ifelse(Days=='Early',0,1)) %>%
                   left_join(select(patient.data, -GISymptoms, -PatientID),
                             by = c('PatientID'='ShortID')) %>%
                   mutate(Status = ifelse(str_detect(PatientID,'K'),0,1)) %>%
                   select(-contains('OSCI_Class')) %>%
                   dplyr::rename(OSCI_Score_Worst = 'OSCI_Score') %>%
                   arrange(GenID) %>%
                   remove.sparse.vars(keep = c('Antibiotics','DaysHospitalized',
                                               'Pneumonia','Bacteremia'), 
                                      na.thresh = params$na.thresh.glob) %>% 
                   prep.metadata() %>%
                   mutate(Sex = ifelse(Sex==0,1,0)) %>%
                   select(Status, everything()) %>%
                   column_to_rownames('GenID')),
    tar_target(mdc.metabs.all.res, mdc.metabs.all.feat %>%
                   metadeconfoundR::MetaDeconfound(
                       mdc.metabs.all.meta[,-which(colnames(mdc.metabs.all.meta)=='Days')],
                       deconfT = c('Antibiotics'),
                       randomVar = list('+ (1|PatientID)',
                                        c('PatientID')),
                       nnodes = params$nnodes,
                       robustCutoff = 0,
                       logfile = here('logfiles','log.txt')) %>%
                   map(~ rownames_to_column(as.data.frame(.), 'feature')) %>%
                   map(~ as_tibble(.)) %>%
                   tibble.results(., 'tbl.all', 'COVID') %>%
                   mutate(Confounding = str_replace_all(Confounding,'Status','COVID')) %>%
                   filter(Covariate != 'PatientID')),

    # plasma & urine -- early timepoints only -- status = COVID ----
    tar_target(mdc.metabs.early.meta, mdc.metabs.all.meta %>%
                   filter(str_detect(PatientID,'K') | Days == 0) %>%
                   select(-Days)),
    tar_target(mdc.metabs.early.feat, 
               mdc.metabs.all.feat[rownames(mdc.metabs.early.meta),]),
    tar_target(mdc.metabs.early.res, mdc.metabs.early.feat %>%
                   metadeconfoundR::MetaDeconfound(mdc.metabs.early.meta,
                                                   deconfT = c('Antibiotics'),
                                                   randomVar = list('+ (1|PatientID)',
                                                                    c('PatientID')),
                                                   nnodes = params$nnodes,
                                                   robustCutoff = 0,
                                                   logfile = here('logfiles','log.txt')) %>%
                   map(~ rownames_to_column(as.data.frame(.), 'feature')) %>%
                   map(~ as_tibble(.)) %>%
                   tibble.results(., 'tbl.all', 'COVID') %>%
                   mutate(Confounding = str_replace_all(Confounding,'Status','COVID')) %>%
                   filter(Covariate != 'PatientID')),
    
    # plasma & urine -- late timepoints only -- status = COVID ----
    tar_target(mdc.metabs.late.meta, mdc.metabs.all.meta %>%
                   filter(str_detect(PatientID,'K') | Days == 1) %>%
                   select(-Days)),
    tar_target(mdc.metabs.late.feat, 
               mdc.metabs.all.feat[rownames(mdc.metabs.late.meta),]),
    tar_target(mdc.metabs.late.res, mdc.metabs.late.feat %>%
                   metadeconfoundR::MetaDeconfound(mdc.metabs.late.meta,
                                                   deconfT = c('Antibiotics'),
                                                   randomVar = list('+ (1|PatientID)',
                                                                    c('PatientID')),
                                                   nnodes = params$nnodes,
                                                   robustCutoff = 0,
                                                   logfile = here('logfiles','log.txt')) %>%
                   map(~ rownames_to_column(as.data.frame(.), 'feature')) %>%
                   map(~ as_tibble(.)) %>%
                   tibble.results(., 'tbl.all', 'COVID') %>%
                   mutate(Confounding = str_replace_all(Confounding,'Status','COVID')) %>%
                   filter(Covariate != 'PatientID')),
    # plasma -- all timepoints, chemclass -- status = COVID ----
    tar_target(mdc.p.chemclass.all.feat, bin.metabs.simp(
        p.metabolites.wide, 
        metabolite.categ %>%
            mutate(ChemClass = fct_collapse(
                ChemClass, `Amino acids & related` = c('Urea cycle',
                                                       'AA related','Amino acid','Amino acid derivative'),
                `Aromatic AAs & derivatives` = c(
                    'Aromatic amino acid','Aromatic amino acid derivative'),
                `Indole & derivatives` = c(
                    'Indole','Indole derivatives'),
                `Biogenic amines` = c(
                    'Monoamine','Monoamine derivative', 'Biogenic amines'),
                `Glycolysis & TCA` = c('Glycolysis','TCA'),
                `Carboxylic acids` = c('Carboxylic acids', 
                                       'Quinoline carboxylic acids'),
                `Choline & derivatives` = c(
                    'Amine oxides','Vitamins & Cofactors'),
                `Sugars` = c('Pentoses','Sugars'),
                `Ceramides` = c('Ceramides','Dihexosylceramides','Hexosylceramides'),
                `Nucleoside & related` = c(
                    'Nucleobase-related','Nucleoside','Nucleobase'))) ) %>%
            arrange(UniqueID) %>%
            column_to_rownames('UniqueID')),
    # run standard metadeconfoundr just with binned metabolites
    tar_target(mdc.p.chemclass.all.res, mdc.p.chemclass.all.feat %>%
                   metadeconfoundR::MetaDeconfound(
                       mdc.metabs.all.meta,
                       randomVar = list('+ (1|PatientID)',
                                        c('PatientID')),
                       nnodes = params$nnodes,
                       robustCutoff = 0,
                       logfile = here('logfiles','log.txt')) %>%
                   map(~ rownames_to_column(as.data.frame(.), 'feature')) %>%
                   map(~ as_tibble(.)) %>%
                   tibble.results(., 'tbl.all', 'COVID') %>%
                   mutate(Confounding = str_replace_all(Confounding,'Status','COVID')) %>%
                   filter(Covariate != 'PatientID')),
    
    # calculate cliff's delta and wilcoxon for mild/severe subsets by hand
    tar_target(mdc.p.chemclass.mild, list(
        'Effect Size' = chem.class.stats(
            mdc.p.chemclass.all.feat,
            mdc.metabs.all.meta,
            'mild', 'cliffs'),
        'FDR' = chem.class.stats(
            mdc.p.chemclass.all.feat,
            mdc.metabs.all.meta,
            'mild','wilcox')) %>%
            reduce(left_join) %>%
            left_join(filter(mdc.p.chemclass.all.res,
                             Covariate=='OSCI_Score_Samp') %>%
                          select(feature, Confounding))),
    tar_target(mdc.p.chemclass.severe, list(
        'Effect Size' = chem.class.stats(
            mdc.p.chemclass.all.feat,
            mdc.metabs.all.meta,
            'severe', 'cliffs'),
        'FDR' = chem.class.stats(
            mdc.p.chemclass.all.feat,
            mdc.metabs.all.meta,
            'severe','wilcox')) %>%
            reduce(left_join) %>%
            left_join(filter(mdc.p.chemclass.all.res,
                             Covariate=='OSCI_Score_Samp') %>%
                          select(feature, Confounding))),
    
    # two more metadeconfoundR runs with these subsets
    tar_target(ctrl.idx, which(mdc.metabs.all.meta$OSCI_Score_Samp==0)),
    tar_target(mild.idx, intersect(which(
        mdc.metabs.all.meta$OSCI_Score_Samp>0), which(
            mdc.metabs.all.meta$OSCI_Score_Samp<5))),
    tar_target(severe.idx, which(mdc.metabs.all.meta$OSCI_Score_Samp>=5)),
    tar_target(early.idx, which(mdc.metabs.all.meta$Days==0)),
    tar_target(late.idx, which(mdc.metabs.all.meta$Days==1)),
    tar_target(mdc.p.chemclass.mild.early.res,
                   metadeconfoundR::MetaDeconfound(
                       mdc.p.chemclass.all.feat[c(ctrl.idx,
                                                intersect(mild.idx,
                                                          early.idx)),],
                       mdc.metabs.all.meta[c(ctrl.idx, 
                                             intersect(mild.idx,
                                                       early.idx)),] %>%
                           select(-contains('OSCI')),
                       randomVar = list('+ (1|PatientID)',
                                        c('PatientID')),
                       nnodes = params$nnodes,
                       robustCutoff = 0,
                       logfile = here('logfiles','log.txt')) %>%
                   map(~ rownames_to_column(as.data.frame(.), 'feature')) %>%
                   map(~ as_tibble(.)) %>%
                   tibble.results(., 'tbl.all', 'COVID') %>%
                   mutate(Confounding = str_replace_all(
                       Confounding,'Status','COVID')) %>%
                   filter(Covariate != 'PatientID')),
    tar_target(mdc.p.chemclass.mild.late.res,
               metadeconfoundR::MetaDeconfound(
                   mdc.p.chemclass.all.feat[c(ctrl.idx,
                                            intersect(mild.idx,
                                                      late.idx)),],
                   mdc.metabs.all.meta[c(ctrl.idx, 
                                         intersect(mild.idx,
                                                   late.idx)),] %>%
                       select(-contains('OSCI')),
                   randomVar = list('+ (1|PatientID)',
                                    c('PatientID')),
                   nnodes = params$nnodes,
                   robustCutoff = 0,
                   logfile = here('logfiles','log.txt')) %>%
                   map(~ rownames_to_column(as.data.frame(.), 'feature')) %>%
                   map(~ as_tibble(.)) %>%
                   tibble.results(., 'tbl.all', 'COVID') %>%
                   mutate(Confounding = str_replace_all(
                       Confounding,'Status','COVID')) %>%
                   filter(Covariate != 'PatientID')),
    tar_target(mdc.p.chemclass.severe.early.res,
               metadeconfoundR::MetaDeconfound(
                   mdc.p.chemclass.all.feat[c(ctrl.idx, 
                                            intersect(severe.idx,
                                                      early.idx)),],
                   mdc.metabs.all.meta[c(ctrl.idx, 
                                         intersect(severe.idx,
                                                   early.idx)),] %>%
                       select(-contains('OSCI')),
                   randomVar = list('+ (1|PatientID)',
                                    c('PatientID')),
                   nnodes = params$nnodes,
                   robustCutoff = 0,
                   logfile = here('logfiles','log.txt')) %>%
                   map(~ rownames_to_column(as.data.frame(.), 'feature')) %>%
                   map(~ as_tibble(.)) %>%
                   tibble.results(., 'tbl.all', 'COVID') %>%
                   mutate(Confounding = str_replace_all(
                       Confounding,'Status','COVID')) %>%
                   filter(Covariate != 'PatientID')),
    tar_target(mdc.p.chemclass.severe.late.res,
               metadeconfoundR::MetaDeconfound(
                   mdc.p.chemclass.all.feat[c(ctrl.idx, 
                                            intersect(severe.idx,
                                                      late.idx)),],
                   mdc.metabs.all.meta[c(ctrl.idx, 
                                         intersect(severe.idx,
                                                   late.idx)),] %>%
                       select(-contains('OSCI')),
                   randomVar = list('+ (1|PatientID)',
                                    c('PatientID')),
                   nnodes = params$nnodes,
                   robustCutoff = 0,
                   logfile = here('logfiles','log.txt')) %>%
                   map(~ rownames_to_column(as.data.frame(.), 'feature')) %>%
                   map(~ as_tibble(.)) %>%
                   tibble.results(., 'tbl.all', 'COVID') %>%
                   mutate(Confounding = str_replace_all(
                       Confounding,'Status','COVID')) %>%
                   filter(Covariate != 'PatientID')),
    
    # urine -- all timepoints, chemclass -- status = COVID ----
    tar_target(mdc.u.chemclass.all.feat, bin.metabs.simp(
        u.metabolites.wide, 
        metabolite.categ %>%
            mutate(ChemClass = fct_collapse(
                ChemClass, `Amino acids & related` = c('Urea cycle',
                                                       'AA related','Amino acid','Amino acid derivative'),
                `Aromatic AAs & derivatives` = c(
                    'Aromatic amino acid','Aromatic amino acid derivative'),
                `Indole & derivatives` = c(
                    'Indole','Indole derivatives'),
                `Biogenic amines` = c(
                    'Monoamine','Monoamine derivative', 'Biogenic amines'),
                `Glycolysis & TCA` = c('Glycolysis','TCA'),
                `Carboxylic acids` = c('Carboxylic acids', 
                                       'Quinoline carboxylic acids'),
                `Choline & derivatives` = c(
                    'Amine oxides','Vitamins & Cofactors'),
                `Sugars` = c('Pentoses','Sugars'),
                `Ceramides` = c('Ceramides','Dihexosylceramides','Hexosylceramides'),
                `Nucleoside & related` = c(
                    'Nucleobase-related','Nucleoside','Nucleobase'))) ) %>%
            arrange(UniqueID) %>%
            column_to_rownames('UniqueID')),
    # run standard metadeconfoundr just with binned metabolites
    tar_target(mdc.u.chemclass.all.res, mdc.metabs.all.meta %>%
                   filter(rownames(.) %in% rownames(
                       mdc.u.chemclass.all.feat)) %>%
                   metadeconfoundR::MetaDeconfound(
                       mdc.u.chemclass.all.feat %>%
                           filter(rownames(.) %in% rownames(
                               mdc.metabs.all.meta)),
                       .,
                       randomVar = list('+ (1|PatientID)',
                                        c('PatientID')),
                       nnodes = params$nnodes,
                       robustCutoff = 0,
                       logfile = here('logfiles','log.txt')) %>%
                   map(~ rownames_to_column(as.data.frame(.), 'feature')) %>%
                   map(~ as_tibble(.)) %>%
                   tibble.results(., 'tbl.all', 'COVID') %>%
                   mutate(Confounding = str_replace_all(Confounding,'Status','COVID')) %>%
                   filter(Covariate != 'PatientID')),
    
    # calculate cliff's delta and wilcoxon for mild/severe subsets by hand
    tar_target(mdc.u.chemclass.mild, list(
        'Effect Size' = chem.class.stats(
            mdc.u.chemclass.all.feat,
            mdc.metabs.all.meta,
            'mild', 'cliffs'),
        'FDR' = chem.class.stats(
            mdc.u.chemclass.all.feat,
            mdc.metabs.all.meta,
            'mild','wilcox')) %>%
            reduce(left_join) %>%
            left_join(filter(mdc.u.chemclass.all.res,
                             Covariate=='OSCI_Score_Samp') %>%
                          select(feature, Confounding))),
    tar_target(mdc.u.chemclass.severe, list(
        'Effect Size' = chem.class.stats(
            mdc.u.chemclass.all.feat,
            mdc.metabs.all.meta,
            'severe', 'cliffs'),
        'FDR' = chem.class.stats(
            mdc.u.chemclass.all.feat,
            mdc.metabs.all.meta,
            'severe','wilcox')) %>%
            reduce(left_join) %>%
            left_join(filter(mdc.u.chemclass.all.res,
                             Covariate=='OSCI_Score_Samp') %>%
                          select(feature, Confounding))),
    
    tar_target(mdc.u.chemclass.mild.early.res,
               metadeconfoundR::MetaDeconfound(
                   mdc.u.chemclass.all.feat[c(ctrl.idx,
                                              intersect(mild.idx,
                                                        early.idx)),],
                   mdc.metabs.all.meta[c(ctrl.idx, 
                                         intersect(mild.idx,
                                                   early.idx)),] %>%
                       select(-contains('OSCI')),
                   randomVar = list('+ (1|PatientID)',
                                    c('PatientID')),
                   nnodes = params$nnodes,
                   robustCutoff = 0,
                   logfile = here('logfiles','log.txt')) %>%
                   map(~ rownames_to_column(as.data.frame(.), 'feature')) %>%
                   map(~ as_tibble(.)) %>%
                   tibble.results(., 'tbl.all', 'COVID') %>%
                   mutate(Confounding = str_replace_all(
                       Confounding,'Status','COVID')) %>%
                   filter(Covariate != 'PatientID')),
    tar_target(mdc.u.chemclass.mild.late.res,
               metadeconfoundR::MetaDeconfound(
                   mdc.u.chemclass.all.feat[c(ctrl.idx,
                                              intersect(mild.idx,
                                                        late.idx)),],
                   mdc.metabs.all.meta[c(ctrl.idx, 
                                         intersect(mild.idx,
                                                   late.idx)),] %>%
                       select(-contains('OSCI')),
                   randomVar = list('+ (1|PatientID)',
                                    c('PatientID')),
                   nnodes = params$nnodes,
                   robustCutoff = 0,
                   logfile = here('logfiles','log.txt')) %>%
                   map(~ rownames_to_column(as.data.frame(.), 'feature')) %>%
                   map(~ as_tibble(.)) %>%
                   tibble.results(., 'tbl.all', 'COVID') %>%
                   mutate(Confounding = str_replace_all(
                       Confounding,'Status','COVID')) %>%
                   filter(Covariate != 'PatientID')),
    tar_target(mdc.u.chemclass.severe.early.res,
               metadeconfoundR::MetaDeconfound(
                   mdc.u.chemclass.all.feat[c(ctrl.idx, 
                                              intersect(severe.idx,
                                                        early.idx)),],
                   mdc.metabs.all.meta[c(ctrl.idx, 
                                         intersect(severe.idx,
                                                   early.idx)),] %>%
                       select(-contains('OSCI')),
                   randomVar = list('+ (1|PatientID)',
                                    c('PatientID')),
                   nnodes = params$nnodes,
                   robustCutoff = 0,
                   logfile = here('logfiles','log.txt')) %>%
                   map(~ rownames_to_column(as.data.frame(.), 'feature')) %>%
                   map(~ as_tibble(.)) %>%
                   tibble.results(., 'tbl.all', 'COVID') %>%
                   mutate(Confounding = str_replace_all(
                       Confounding,'Status','COVID')) %>%
                   filter(Covariate != 'PatientID')),
    tar_target(mdc.u.chemclass.severe.late.res,
               metadeconfoundR::MetaDeconfound(
                   mdc.u.chemclass.all.feat[c(ctrl.idx, 
                                              intersect(severe.idx,
                                                        late.idx)),],
                   mdc.metabs.all.meta[c(ctrl.idx, 
                                         intersect(severe.idx,
                                                   late.idx)),] %>%
                       select(-contains('OSCI')),
                   randomVar = list('+ (1|PatientID)',
                                    c('PatientID')),
                   nnodes = params$nnodes,
                   robustCutoff = 0,
                   logfile = here('logfiles','log.txt')) %>%
                   map(~ rownames_to_column(as.data.frame(.), 'feature')) %>%
                   map(~ as_tibble(.)) %>%
                   tibble.results(., 'tbl.all', 'COVID') %>%
                   mutate(Confounding = str_replace_all(
                       Confounding,'Status','COVID')) %>%
                   filter(Covariate != 'PatientID')),
    
    # cleanup to send to ela & ulrike ----
    tar_target(metabs.clean.res, mdc.metabs.early.res %>%
                   bind_rows(mdc.metabs.late.res, mdc.metabs.all.res,
                             .id = 'Analysis') %>%
                   mutate(Analysis = fct_recode(Analysis,
                                                Early = '1', Late = '2',
                                                All = '3')) %>%
                   mutate(Sample = str_split(feature,'_') %>%
                              map_chr(~ head(., 1))) %>%
                   mutate(Sample = fct_recode(Sample,
                                                  Plasma = 'P',
                                                  Urine = 'U')) %>%
                   mutate(Metabolite = str_remove_all(feature,'P_')) %>%
                   mutate(Metabolite = str_remove_all(Metabolite,'U_')) %>%
                   filter(Covariate == 'OSCI_Score_Samp') %>%
                   mutate(ConfStatus = str_split(Confounding,',')) %>%
                   unnest(ConfStatus) %>%
                   mutate(ConfStatus = trimws(ConfStatus)) %>%
                   mutate(ConfStatus = case_when(
                       str_detect(ConfStatus,'SD|LD|NC|Worst') ~ 'R',
                       ConfStatus == 'COVID' ~ 'COVID',
                       ConfStatus == 'NS' ~ 'NS',
                       str_detect(ConfStatus, 'AbxSamp|Antibiotics|AbxHosp') ~ 'AbxCurr',
                       str_detect(ConfStatus, 'Abx3Mo') ~ 'Abx3Mo', 
                       str_detect(ConfStatus,'Pneu|DaysHosp|Bacter') ~ 'Hosp/HAP',
                       str_detect(ConfStatus,'Comorb|Medic') ~ 'Comorbidity related',
                       str_detect(ConfStatus,
                                  'Renal|Diur|Hyper|Antihyper') ~ 'Kidney/BP related',
                       TRUE ~ 'Other')) %>%
                   left_join(select(metabolite.categ, R_Compound, Compound),
                             by=c('Metabolite'='R_Compound')) %>%
                   select(Analysis, Sample, Metabolite, Compound, Covariate,
                          `Effect Size`, FDR, ConfStatus) ),
    tar_target(p.metabs.early.confcts, metabs.clean.res %>% 
                   filter(Analysis=='Early') %>% 
                   count(ConfStatus) %>% 
                   arrange(desc(n))),
    tar_target(p.metabs.late.confcts, metabs.clean.res %>% 
                   filter(Analysis=='Late') %>% 
                   count(ConfStatus) %>% 
                   arrange(desc(n))),
    tar_target(p.metabs.all.confcts, metabs.clean.res %>% 
                   filter(Analysis=='All') %>% 
                   count(ConfStatus) %>% 
                   arrange(desc(n))),
    
    ### METADECONFOUNDR RUNS: METABOLITES + MICROBIOTA -----
    tar_target(all.alphas, sdna.alpha %>%
                   mutate(GenID = str_remove_all(UniqueID, '-SDNA')) %>%
                   select(GenID, everything(), -UniqueID) %>%
                   magrittr::set_colnames(c('GenID','SDNA_Richness',
                                            'SDNA_Shannon','SDNA_InvSimp')) %>%
                   full_join(op.alpha %>%
                                 mutate(
                                     GenID = str_remove_all(UniqueID, '-OP')) %>%
                                 select(GenID, everything(), -UniqueID) %>%
                                 magrittr::set_colnames(c('GenID','OP_Richness',
                                                          'OP_Shannon','OP_InvSimp')))),
    # taxonomic
    tar_target(osci.sdna.motu.names, mdc.sdna.motus.all.res %>%
                   filter(Covariate == 'OSCI_Score_Samp' & Confounding %in% c(
                       'OSCI_Score_Worst','SD','LD','NC')) %>%
                   filter(!feature %in% c('Richness','Shannon','InvSimp')) %>%
                   pull(feature)),
    tar_target(osci.sdna.bin.motus, sdna.feat.ra %>%
                   select(UniqueID, all_of(osci.sdna.motu.names)) %>%
                   bin.motus.simp() %>%
                   rename_if(is.double, ~ paste0('S_', .)) %>%
                   mutate(GenID = str_remove_all(UniqueID,'-SDNA')) %>%
                   select(-UniqueID)),
    tar_target(osci.op.motu.names, mdc.op.motus.all.res %>%
                   filter(Covariate == 'OSCI_Score_Samp' & Confounding %in% c(
                       'OSCI_Score_Worst','SD','LD','NC')) %>%
                   filter(!feature %in% c('Richness','Shannon','InvSimp')) %>%
                   pull(feature)),
    tar_target(osci.op.bin.motus, op.feat.ra %>%
                   select(UniqueID, all_of(osci.op.motu.names)) %>%
                   bin.motus.simp() %>%
                   rename_if(is.double, ~ paste0('O_', .)) %>%
                   mutate(GenID = str_remove_all(UniqueID,'-OP')) %>%
                   select(-UniqueID)),
    tar_target(osci.all.bin.motus, osci.sdna.bin.motus %>%
                   full_join(osci.op.bin.motus)),
    # phylum
    tar_target(osci.sdna.phy.names, mdc.sdna.phy.all.res %>% 
                   filter(Covariate=='OSCI_Score_Samp'& Confounding %in% c(
                       'OSCI_Score_Worst','LD','SD','NC')) %>%
                   pull(feature)),
    tar_target(osci.sdna.phy, mdc.sdna.phy.all.feat %>%
                   rownames_to_column('UniqueID') %>%
                   mutate(UniqueID = str_remove_all(UniqueID,'-SDNA')) %>%
                   select(UniqueID, all_of(osci.sdna.phy.names)) %>%
                   rename_if(is.double, ~ paste0('S_', .))),
    tar_target(osci.op.phy.names, mdc.op.phy.all.res %>% 
                   filter(Covariate=='OSCI_Score_Samp'& Confounding %in% c(
                       'OSCI_Score_Worst','LD','SD','NC')) %>%
                   pull(feature)),
    tar_target(osci.op.phy, mdc.op.phy.all.feat %>%
                   rownames_to_column('UniqueID') %>%
                   mutate(UniqueID = str_remove_all(UniqueID,'-OP')) %>%
                   select(UniqueID, all_of(osci.op.phy.names)) %>%
                   rename_if(is.double, ~ paste0('O_',.))),
    # functional? - not unless able to bin @ brite hierarchies
    
    # plasma & urine -- all timepoints w/ bin mOTUs, phylum -- status = COVID ----
    # can use mdc.metabs.all.feat -- no longer contain V.2 repeats
    tar_target(mdc.metabs.motus.all.meta, mdc.metabs.all.meta %>%
                   rownames_to_column('GenID') %>%
                   as_tibble() %>%
                   left_join(all.alphas) %>%
                   left_join(osci.all.bin.motus) %>%
                   column_to_rownames('GenID')),
    tar_target(mdc.metabs.motus.all.res, 
               metadeconfoundR::MetaDeconfound(
                   mdc.metabs.all.feat,
                   mdc.metabs.motus.all.meta,
                   deconfT = c('Antibiotics'),
                   randomVar = list('+ (1|PatientID)',
                                    c('PatientID')),
                   nnodes = params$nnodes,
                   robustCutoff = 0,
                   logfile = here('logfiles','log.txt')) %>%
                   map(~ rownames_to_column(as.data.frame(.), 'feature')) %>%
                   map(~ as_tibble(.)) %>%
                   tibble.results(., 'tbl.all', 'COVID') %>%
                   mutate(Confounding = str_replace_all(Confounding,'Status','COVID')) %>%
                   dplyr::filter(Covariate != 'PatientID')),
    
    # plasma & urine -- early timepoints w/ bin mOTUs, phylum -- status = COVID ----
    tar_target(mdc.metabs.motus.early.meta, mdc.metabs.motus.all.meta %>%
                   filter(Days == 0 | str_detect(PatientID,'K')) %>%
                   select(-Days)),
    tar_target(mdc.metabs.motus.early.feat, mdc.metabs.all.feat %>%
                   filter(rownames(mdc.metabs.all.feat) %in% rownames(
                       mdc.metabs.motus.early.meta))),
    tar_target(mdc.metabs.motus.early.res, 
               metadeconfoundR::MetaDeconfound(
                   mdc.metabs.motus.early.feat,
                   mdc.metabs.motus.early.meta,
                   deconfT = c('Antibiotics'),
                   randomVar = list('+ (1|PatientID)',
                                    c('PatientID')),
                   nnodes = params$nnodes,
                   robustCutoff = 0,
                   logfile = here('logfiles','log.txt')) %>%
                   map(~ rownames_to_column(as.data.frame(.), 'feature')) %>%
                   map(~ as_tibble(.)) %>%
                   tibble.results(., 'tbl.all', 'COVID') %>%
                   mutate(Confounding = str_replace_all(Confounding,'Status','COVID')) %>%
                   filter(Covariate != 'PatientID')),
    
    # plasma & urine -- late timepoints  w/ bin mOTUs, phylum -- status = COVID ----
    tar_target(mdc.metabs.motus.late.meta, mdc.metabs.motus.all.meta %>%
                   filter(Days == 1 | str_detect(PatientID,'K')) %>%
                   select(-Days)),
    tar_target(mdc.metabs.motus.late.feat, mdc.metabs.all.feat %>%
                   filter(rownames(mdc.metabs.all.feat) %in% rownames(
                       mdc.metabs.motus.late.meta))),
    tar_target(mdc.metabs.motus.late.res, 
               metadeconfoundR::MetaDeconfound(
                   mdc.metabs.motus.late.feat,
                   mdc.metabs.motus.late.meta,
                   deconfT = c('Antibiotics'),
                   randomVar = list('+ (1|PatientID)',
                                    c('PatientID')),
                   nnodes = params$nnodes,
                   robustCutoff = 0,
                   logfile = here('logfiles','log.txt')) %>%
                   map(~ rownames_to_column(as.data.frame(.), 'feature')) %>%
                   map(~ as_tibble(.)) %>%
                   tibble.results(., 'tbl.all', 'COVID') %>%
                   mutate(Confounding = str_replace_all(Confounding,'Status','COVID')) %>%
                   filter(Covariate != 'PatientID')),
    
    ### METADECONFOUNDR RUNS: HOST RESPONSE ----
    # cytokines -- all timepoints -- status = COVID ----
    tar_target(mdc.host.ck.feat, cytokines %>%
                   left_join(clinical.feat %>%
                                 select(GenID, IL6) %>%
                                 distinct() %>%
                                 mutate(IL6 = as.numeric(IL6)),
                             by = c('UniqueID'='GenID')) %>%
                   arrange(desc(UniqueID)) %>%
                   select(-Days, -PatientID) %>%
                   rename_all(~ str_replace_all(., '-','_')) %>%
                   column_to_rownames('UniqueID')),
    tar_target(mdc.host.ck.basic.meta, clinical.feat %>%
                   select(GenID, PatientID, Days, Antibiotics,
                          OSCI_Score_Samp ,contains('VL')) %>%
                   distinct() %>%
                   dplyr::filter(GenID %in% rownames(mdc.host.ck.feat)) %>%
                   discretize.days() %>%
                   mutate(Days = ifelse(Days=='Early',0,1)) %>%
                   left_join(select(patient.data, -GISymptoms, -PatientID),
                             by = c('PatientID'='ShortID')) %>%
                   mutate(Status = ifelse(str_detect(PatientID,'K'),0,1)) %>%
                   select(-contains('OSCI_Class')) %>%
                   dplyr::rename(OSCI_Score_Worst = 'OSCI_Score') %>%
                   arrange(desc(GenID)) %>%
                   remove.sparse.vars(keep = c('Antibiotics','DaysHospitalized',
                                               'Pneumonia','Bacteremia'), 
                                      na.thresh = params$na.thresh.glob) %>% 
                   prep.metadata() %>%
                   mutate(Sex = ifelse(Sex==0,1,0)) %>%
                   select(Status, everything()) %>%
                   distinct(across(GenID), .keep_all = TRUE) %>%
                   column_to_rownames('GenID')),
    tar_target(mdc.host.ck.basic.res, 
               metadeconfoundR::MetaDeconfound(
                   mdc.host.ck.feat,
                   mdc.host.ck.basic.meta,
                   deconfT = c('Antibiotics'),
                   randomVar = list('+ (1|PatientID)',
                                    c('PatientID')),
                   nnodes = params$nnodes,
                   robustCutoff = 0,
                   logfile = here('logfiles','log.txt')) %>%
                   map(~ rownames_to_column(as.data.frame(.), 'feature')) %>%
                   map(~ as_tibble(.)) %>%
                   tibble.results(., 'tbl.all', 'COVID') %>%
                   mutate(Confounding = str_replace_all(
                       Confounding,'Status','COVID')) %>%
                   dplyr::filter(Covariate != 'PatientID')),
    
    # PBMCs -- all timepoints -- status = COVID ----
    tar_target(mdc.host.pbmc.feat, combined.pbmc %>%
                   arrange(UniqueID) %>%
                   remove.zero.cols() %>%
                   near.zero.var.filter(n.header.cols = 1) %>%
                   column_to_rownames('UniqueID')),
    tar_target(mdc.host.pbmc.basic.meta, clinical.feat %>%
                   select(GenID, PatientID, Antibiotics, OSCI_Score_Samp,
                          contains('VL')) %>%
                   distinct() %>%
                   dplyr::filter(GenID %in% rownames(mdc.host.pbmc.feat)) %>%
                   left_join(select(patient.data, -GISymptoms, -PatientID),
                             by = c('PatientID'='ShortID')) %>%
                   mutate(Status = ifelse(str_detect(PatientID,'K'),0,1)) %>%
                   select(-contains('OSCI_Class')) %>%
                   dplyr::rename(OSCI_Score_Worst = 'OSCI_Score') %>%
                   arrange(desc(GenID)) %>%
                   # remove.sparse.vars(keep = c('Antibiotics','DaysHospitalized',
                   #                             'Pneumonia','Bacteremia'), 
                   #                    na.thresh = params$na.thresh.glob) %>% 
                   prep.metadata() %>%
                   mutate(Sex = ifelse(Sex==0,1,0)) %>%
                   select(Status, everything(), -PatientID) %>%
                   column_to_rownames('GenID')),
    tar_target(mdc.host.pbmc.basic.res, 
               metadeconfoundR::MetaDeconfound(
                   mdc.host.pbmc.feat,
                   mdc.host.pbmc.basic.meta,
                   deconfT = c('Antibiotics'),
                   nnodes = params$nnodes,
                   robustCutoff = 0,
                   logfile = here('logfiles','log.txt')) %>%
                   map(~ rownames_to_column(as.data.frame(.), 'feature')) %>%
                   map(~ as_tibble(.)) %>%
                   tibble.results(., 'tbl.all', 'COVID') %>%
                   mutate(Confounding = str_replace_all(
                       Confounding,'Status','COVID')) %>%
                   dplyr::filter(Covariate != 'PatientID')),
    
    ### METADECONFOUNDR RUNS: HOST RESPONSE + MICROBIOTA + METABOLITES -----
    tar_target(osci.metab.names, mdc.metabs.all.res %>%
                   filter(!str_detect(feature, '_13C|_D3|_D5')) %>%
                   filter(Covariate == 'OSCI_Score_Samp' & Confounding %in% c(
                       'OSCI_Score_Worst','SD','LD','NC')) %>%
                   mutate(R_Compound = str_remove_all(feature,'P_|U_')) %>%
                   pull(feature)),
    tar_target(osci.metabs, p.metabolites.wide %>%
                   left_join(u.metabolites.wide) %>%
                   select(UniqueID, all_of(osci.metab.names))),
    tar_target(osci.metab.class.names, mdc.metabs.all.res %>%
                   filter(!str_detect(feature, '_13C|_D3|_D5')) %>%
                   filter(Covariate == 'OSCI_Score_Samp' & Confounding %in% c(
                       'OSCI_Score_Worst','SD','LD','NC')) %>%
                   mutate(R_Compound = str_remove_all(feature,'P_|U_')) %>%
                   left_join(select(metabolite.categ, R_Compound,
                                    ChemClass, Pathway)) %>%
                   bin.metab.results(full_join(p.metabolites.wide,
                                             u.metabolites.wide),
                                   by = 'ChemClass')),
    
    # cytokines -- all timepoints -- status = COVID ----
    tar_target(mdc.host.ck.meta, clinical.feat %>%
                   select(GenID, PatientID, Days, Antibiotics,
                          OSCI_Score_Samp ,contains('VL')) %>%
                   distinct() %>%
                   filter(GenID %in% rownames(mdc.host.ck.feat)) %>%
                   discretize.days() %>%
                   mutate(Days = ifelse(Days=='Early',0,1)) %>%
                   left_join(select(patient.data, -GISymptoms, -PatientID),
                             by = c('PatientID'='ShortID')) %>%
                   # add microbiota features
                   left_join(all.alphas) %>%
                   left_join(osci.all.bin.motus) %>%
                   left_join(osci.sdna.phy, by=c('GenID'='UniqueID')) %>%
                   left_join(osci.op.phy, by=c('GenID'='UniqueID')) %>%
                   # add metabolome features
                   left_join(osci.metabs, by=c('GenID'='UniqueID')) %>%
                   mutate(Status = ifelse(str_detect(PatientID,'K'),0,1)) %>%
                   select(-contains('OSCI_Class')) %>%
                   dplyr::rename(OSCI_Score_Worst = 'OSCI_Score') %>%
                   arrange(desc(GenID)) %>%
                   remove.sparse.vars(keep = c('Antibiotics','DaysHospitalized',
                                               'Pneumonia','Bacteremia',
                                               'P_GDCA'), 
                                      na.thresh = params$na.thresh.glob) %>% 
                   prep.metadata() %>%
                   mutate(Sex = ifelse(Sex==0,1,0)) %>%
                   select(Status, everything()) %>%
                   distinct(across(GenID), .keep_all = TRUE) %>%
                   column_to_rownames('GenID')),
    tar_target(mdc.host.ck.res, 
               metadeconfoundR::MetaDeconfound(
                   mdc.host.ck.feat,
                   mdc.host.ck.meta,
                   deconfT = c('Antibiotics'),
                   randomVar = list('+ (1|PatientID)',
                                    c('PatientID')),
                   nnodes = params$nnodes,
                   robustCutoff = 0,
                   logfile = here('logfiles','log.txt')) %>%
                   map(~ rownames_to_column(as.data.frame(.), 'feature')) %>%
                   map(~ as_tibble(.)) %>%
                   tibble.results(., 'tbl.all', 'COVID') %>%
                   mutate(Confounding = str_replace_all(
                       Confounding,'Status','COVID')) %>%
                   filter(Covariate != 'PatientID')),
    
    # PBMCs -- n=25 all timepoints (which are early) -- status = COVID ----
    tar_target(mdc.host.pbmc.meta, clinical.feat %>%
                   select(GenID, PatientID, Antibiotics, OSCI_Score_Samp,
                          contains('VL')) %>%
                   distinct() %>%
                   filter(GenID %in% rownames(mdc.host.pbmc.feat)) %>%
                   left_join(select(patient.data, -GISymptoms, -PatientID),
                             by = c('PatientID'='ShortID')) %>%
                   # add microbiota features
                   left_join(all.alphas) %>%
                   left_join(osci.all.bin.motus) %>%
                   left_join(osci.sdna.phy, by=c('GenID'='UniqueID')) %>%
                   left_join(osci.op.phy, by=c('GenID'='UniqueID')) %>%
                   # add metabolome features
                   left_join(osci.metabs, by=c('GenID'='UniqueID')) %>%
                   mutate(Status = ifelse(str_detect(PatientID,'K'),0,1)) %>%
                   select(-contains('OSCI_Class')) %>%
                   dplyr::rename(OSCI_Score_Worst = 'OSCI_Score') %>%
                   arrange(desc(GenID)) %>%
                   # remove.sparse.vars(keep = c('Antibiotics','DaysHospitalized',
                   #                             'Pneumonia','Bacteremia'), 
                   #                    na.thresh = params$na.thresh.glob) %>% 
                   prep.metadata() %>%
                   mutate(Sex = ifelse(Sex==0,1,0)) %>%
                   select(Status, everything(), -PatientID) %>%
                   column_to_rownames('GenID')),
    # no random effect needed -- 1 sample per patient
    tar_target(mdc.host.pbmc.res, 
               metadeconfoundR::MetaDeconfound(
                   mdc.host.pbmc.feat,
                   mdc.host.pbmc.meta,
                   deconfT = c('Antibiotics'),
                   nnodes = params$nnodes,
                   robustCutoff = 0,
                   logfile = here('logfiles','log.txt')) %>%
                   map(~ rownames_to_column(as.data.frame(.), 'feature')) %>%
                   map(~ as_tibble(.)) %>%
                   tibble.results(., 'tbl.all', 'COVID') %>%
                   mutate(Confounding = str_replace_all(
                       Confounding,'Status','COVID')) %>%
                   filter(Covariate != 'PatientID')),
    
    ###################################################################
    ### METADECONFOUNDR INTEGRATION: MICROBIOME-MEDIATED METABOLITES ----
    # metab as feat, OSCI covariate, microbiota confounding
    tar_target(tmp.res.1, mdc.metabs.motus.all.res %>% 
                   filter(str_detect(
                       Covariate,'Score_Samp') & str_detect(
                           Confounding,'S_|P_|SDNA_|OP_')) %>%
                   filter(!str_detect(feature, '_13C|_D3|_D5')) %>%
                   mutate(Confounding = str_split(Confounding,',') %>%
                              map(~ trimws(.))) %>%
                   unnest(Confounding) %>%
                   group_by(feature) %>%
                   mutate(n.conf = n()) %>%
                   filter(n.conf <= 3) %>%
                   select(-n.conf) %>%
                   filter(str_detect(Confounding,'S_|P_|SDNA_|OP_')) %>%
                   ungroup() %>%
                   add_column(Data = 'met.osci.mb', .before = 1)), 
    tar_target(metab.motus.sankey.dat, tmp.res.1 %>%
                   mutate(Site = str_split(feature,'_') %>%
                              map_chr(~ head(., 1))) %>%
                   mutate(Site = case_when(
                       Site == 'U' ~ 'Urine metabolome',
                       Site == 'P' ~ 'Plasma metabolome'
                   )) %>%
                   mutate(Feat = str_split(feature,'_') %>%
                              map_chr(~ tail(., 1))) %>%
                   left_join(select(metabolite.categ,
                                    R_Compound, Compound),
                             by = c('Feat'='R_Compound')) %>%
                   mutate(Compound = ifelse(is.na(Compound),
                                            '5-hydroxytryptophan',
                                            Compound)) %>%
                   mutate(Conf_Site = str_split(Confounding,'_') %>%
                              map_chr(~ head(., 1))) %>%
                   mutate(Conf_Site = case_when(
                       Conf_Site %in% c('S','SDNA') ~ 'Gut microbiome',
                       Conf_Site %in% c('O','OP') ~ 'Oropharyngeal microbiome'
                   )) %>%
                   mutate(Conf = str_split(Confounding, '_') %>%
                              map_chr(~ tail(., 1))) %>%
                   mutate(ES = case_when(
                       `Effect Size` > 0 ~ '+OSCI Correlation',
                       `Effect Size` < 0 ~ '-OSCI Correlation',
                       TRUE ~ 'NA')) %>%
                   select(Site, Compound, ES, Conf_Site, Conf)),
    tar_target(metab.motus.sankey, metab.motus.sankey.dat %>%
                   #make_long(Site, ES, Compound, Conf, Conf_Site) %>%
                   make_long(ES, Site, Compound, Conf, Conf_Site) %>%
                   ggplot(aes(x = x, 
                              next_x = next_x, 
                              node = node, 
                              next_node = next_node,
                              fill = factor(node),
                              label = node)) +
                   geom_sankey(flow.alpha = 0.5, node.color = 1) +
                   geom_sankey_label(size = 3.5, color = 1, fill = "white") +
                   # scico::scale_fill_scico_d(palette = 'vikO') +
                   #scico::scale_fill_scico_d(palette = 'batlow') +
                    scico::scale_fill_scico_d(palette = 'roma') +
                   theme_sankey(base_size = 14) +
                   theme(legend.position = 'none') +
                   labs(x = '') ),
    # microbiota predictors from above - osci-associations
    tar_target(tmp.res.2, mdc.sdna.genus.all.res %>%
                   filter(str_detect(Covariate,'Score_Samp')) %>%
                   filter(Confounding %in% 
                              c('OSCI_Score_Worst',
                                'SD','LD','NC') | feature %in% c(
                                    tmp.res.1 %>%
                                        filter(str_detect(Confounding,'S_')) %>%
                                        mutate(Confounding = str_split(
                                            Confounding, '_') %>%
                                                   map_chr(~ tail(.,1))) %>%
                                        pull(Confounding) %>%
                                        unique())) %>%
                   mutate(feature = paste0('S_',feature)) %>%
                   bind_rows(mdc.op.genus.all.res  %>%
                                 filter(str_detect(Covariate,'Score_Samp')) %>%
                                 filter(Confounding %in% 
                                            c('OSCI_Score_Worst',
                                              'SD','LD','NC') | feature %in% c(
                                                  tmp.res.1 %>%
                                                      filter(str_detect(Confounding,'S_')) %>%
                                                      mutate(Confounding = str_split(
                                                          Confounding, '_') %>%
                                                              map_chr(~ tail(.,1))) %>%
                                                      pull(Confounding) %>%
                                                      unique())) %>%
                                 mutate(feature = paste0('O_',feature))) %>%
                   add_column(Data = 'mb.osci', .before = 1) %>%
                   mutate(Confounding = case_when(
                       str_detect(Confounding,
                           'OSCI_Score_Worst|LD|SD|NC') ~ 'Robust',
                       TRUE ~ 'Confounded'))),
    # metabolite correlations with microbiota features 
    tar_target(tmp.res.3, mdc.metabs.motus.all.res %>%
                   mutate(tmp = paste0(feature,'_',Covariate)) %>%
                   filter(tmp %in% c(
                       tmp.res.1 %>%
                              mutate(tmp = paste0(feature,'_',Confounding)) %>%                 
                              pull(tmp) %>%
                              unique())) %>%
                   select(-tmp) %>%
                   add_column(Data = 'met.mb', .before = 1) %>%
                   mutate(Confounding = case_when(
                       Confounding %in% c('LD','SD') ~ 'Robust',
                       TRUE ~ Confounding ))),
    tar_target(microb.med.metab.combine, tmp.res.1 %>%
                   bind_rows(tmp.res.2) %>%
                   bind_rows(tmp.res.3)),
    
    ### METADECONFOUNDR INTEGRATION: MEDIATED HOST FEATURES ----
    tar_target(tmp.res.4, mdc.host.ck.res %>%
                   bind_rows(mdc.host.pbmc.res, .id = 'run') %>%
                   mutate(run = ifelse(run==1, 'Cytokines','PBMCs')) %>%
                   prep.threeway.sankey(metabolite.categ, 
                                        fdr.cutoff = 0.1) %>%
                   select(run, ES_OSCI, ES, feature, 
                          Cov_Site, Covariate, 
                          Confounding, Conf_Site)),
    tar_target(host.all.sankey, tmp.res.4 %>%
                   distinct() %>%
                   # filter(run=='Cytokines' & str_detect(feature,'IL6')) %>%
                   filter(!str_detect(Conf_Site,'Load')) %>%
                   filter(str_detect(Covariate, 'rypto|propi|ynur|onin')) %>%
                   filter(Cov_Site != Conf_Site) %>%
                   filter(str_detect(Conf_Site, 'microbiome')) %>%
                   # filter(feature %in% c('IL6')) %>%
                   # filter(str_detect(ES, 'Pos')) %>%
                   ggsankey::make_long(ES_OSCI, feature, ES, Covariate, 
                                       Confounding, Conf_Site) %>%
                   plot.threeway.sankey()),
    ### METADECONFOUNDR INTEGRATION: ALL RESULTS FROM A GIVEN -OMICS SPACE ----
    tar_target(all.sdna.feats, 
               add_column(mdc.sdna.motus.all.res, 
                          config='motus.all', .before=1) %>%
                   bind_rows(add_column(mdc.sdna.motus.pat.res, 
                                        config='motus.pats', .before=1), 
                             add_column(mdc.sdna.genus.all.res,
                                        config='genus.all', .before=1),
                             add_column(mdc.sdna.genus.pat.res,
                                        config='genus.pats', .before=1),
                             add_column(mdc.sdna.phy.all.res,
                                        config='phy.all', .before=1),
                             add_column(mdc.sdna.genus.all.metab.res,
                                        config='genus.all.metabs', .before=1),
                             .id = 'n') %>%
                   # dplyr::filter(FDR <= 0.05 & Confounding != 'NS') %>%
                   select(config, everything())),
    
    tar_target(all.op.feats, 
               add_column(mdc.op.motus.all.res,
                          config='motus.all', .before=1) %>%
                   bind_rows(add_column(mdc.op.motus.pat.res,
                                        config='motus.pats', .before=1),
                             add_column(mdc.op.genus.all.res,
                                        config='genus.all', .before=1),
                             add_column(mdc.op.genus.pat.res,
                                        config='genus.pats', .before=1),
                             add_column(mdc.op.phy.all.res,
                                        config='phy.all', .before=1),
                             add_column(mdc.op.genus.all.metab.res,
                                        config='genus.all.metabs', .before=1),
                             .id = 'n') %>%
                   # dplyr::filter(FDR <= 0.05 & Confounding != 'NS') %>%
                   select(config, n, everything())),
    
    tar_target(all.metab.feats.clin, 
               add_column(mdc.metabs.late.res,
                          config='metabs.late', .before=1) %>%
                   bind_rows(add_column(mdc.metabs.early.res,
                                        config='metabs.early', .before=1),
                             add_column(mdc.metabs.all.res,
                                        config='metabs.all', .before=1),
                             add_column(mdc.p.chemclass.all.res,
                                        config='chemclass.all', .before=1),
                             .id = 'n') %>%
                   # dplyr::filter(FDR <= 0.05 & Confounding != 'NS') %>%
                   mutate(site = str_split(feature, '_') %>%
                              map_chr(~ head(., 1))) %>%
                   mutate(feature = str_remove_all(feature, 'P_|U_')) %>%
                   left_join(select(metabolite.categ,
                                    R_Compound, Compound),
                             by = c('feature'='R_Compound')) %>%
                   select(config, n, site, Compound, Covariate,
                          `Effect Size`, FDR, Confounding, feature) %>%
                   # only 13C etc
                   dplyr::filter((!str_detect(
                       config, 'chem') || !is.na(Compound)))),
    
    tar_target(all.host.feats,
               add_column(mdc.host.ck.basic.res,
                          config='ck.basic.all', .before=1) %>%
                   bind_rows(add_column(mdc.host.pbmc.basic.res,
                                        config='pbmc.basic.all', .before=1),
                             add_column(mdc.host.ck.res,
                                        config='ck.all', .before=1),
                             add_column(mdc.host.pbmc.res,
                                        config='pbmc.all', .before=1),
                             .id = 'n') %>%
                   # dplyr::filter(FDR <= 0.05 & Confounding != 'NS') %>%
                   select(config, n, everything())),
    
    tar_target(all.metab.feats,
               add_column(mdc.metabs.motus.all.res,
                          config='metabs.motus.all', .before=1) %>%
                   bind_rows(
                       add_column(mdc.metabs.motus.early.res,
                                  config='metabs.motus.early', .before=1),
                       add_column(mdc.metabs.motus.late.res,
                                  config='metabs.motus.late', .before=1),
                       .id = 'n') %>%
                   # dplyr::filter(FDR <= 0.05 & Confounding != 'NS') %>%
                   select(config, n, everything())),
    
    ### METADECONFOUNDR INTEGRATION: OSCI-ASSOCIATED FEATURES ----
    # all osci-associated features regardless of confounder status
    tar_target(osci.sdna.feats, 
               add_column(mdc.sdna.motus.all.res, 
                          config='motus.all', .before=1) %>%
                   bind_rows(add_column(mdc.sdna.motus.pat.res, 
                                        config='motus.pats', .before=1), 
                             add_column(mdc.sdna.genus.all.res,
                                        config='genus.all', .before=1),
                             add_column(mdc.sdna.genus.pat.res,
                                        config='genus.pats', .before=1),
                             add_column(mdc.sdna.phy.all.res,
                                        config='phy.all', .before=1),
                             add_column(mdc.sdna.genus.all.metab.res,
                                        config='genus.all.metabs', .before=1),
                             .id = 'n') %>%
                   # incorrect now with genus.all.metabs but oh well
                   mutate(n = case_when(as.numeric(n)%%2 == 0 ~ 'n=51',
                                        as.numeric(n)%%2 == 1 ~ 'n=81')) %>%
                   dplyr::filter(Covariate == 'OSCI_Score_Samp') %>%
                   dplyr::filter(FDR <= 0.05) %>%
                   select(config, everything())),
    
    # all osci-associated features regardless of confounder status
    tar_target(osci.op.feats, 
               add_column(mdc.op.motus.all.res,
                          config='motus.all', .before=1) %>%
                   bind_rows(add_column(mdc.op.motus.pat.res,
                                        config='motus.pats', .before=1),
                             add_column(mdc.op.genus.all.res,
                                        config='genus.all', .before=1),
                             add_column(mdc.op.genus.pat.res,
                                        config='genus.pats', .before=1),
                             add_column(mdc.op.phy.all.res,
                                        config='phy.all', .before=1),
                             add_column(mdc.op.genus.all.metab.res,
                                        config='genus.all.metabs', .before=1),
                             .id = 'n') %>%
                   mutate(n = case_when(as.numeric(n)%%2 == 0 ~ 'n=64',
                                        as.numeric(n)%%2 == 1 ~ 'n=94')) %>%
                   dplyr::filter(Covariate == 'OSCI_Score_Samp') %>%
                   dplyr::filter(FDR <= 0.05) %>%
                   select(config, n, everything())),
    
    
    # all osci-associated features regardless of confounder status
    tar_target(osci.metab.feats.clinpred, 
               add_column(mdc.metabs.late.res,
                          config='metabs.late', .before=1) %>%
                   bind_rows(add_column(mdc.metabs.early.res,
                                        config='metabs.early', .before=1),
                             add_column(mdc.metabs.all.res,
                                        config='metabs.all', .before=1),
                             add_column(mdc.p.chemclass.all.res,
                                        config='chemclass.all', .before=1),
                             .id = 'n') %>%
                   mutate(n = case_when(as.numeric(n)==1 ~ 'n=70',
                                        as.numeric(n)==2 ~ 'n=56',
                                        as.numeric(n)==3 ~ 'n=96',
                                        as.numeric(n)==4 ~ 'n=96')) %>%
                   dplyr::filter(Covariate == 'OSCI_Score_Samp') %>%
                   dplyr::filter(FDR <= 0.05) %>%
                   mutate(site = str_split(feature, '_') %>%
                              map_chr(~ head(., 1))) %>%
                   mutate(feature = str_remove_all(feature, 'P_|U_')) %>%
                   left_join(select(metabolite.categ,
                                    R_Compound, Compound),
                             by = c('feature'='R_Compound')) %>%
                   select(config, n, site, Compound, Covariate,
                          `Effect Size`, FDR, Confounding, feature) %>%
                   # only 13C etc
                   dplyr::filter((!str_detect(
                       config, 'chem') || !is.na(Compound)))),
    
    # all osci-associated features regardless of confounder status
    tar_target(osci.host.feats,
               add_column(mdc.host.ck.basic.res,
                          config='ck.basic.all', .before=1) %>%
                   bind_rows(add_column(mdc.host.pbmc.basic.res,
                                        config='pbmc.basic.all', .before=1),
                             add_column(mdc.host.ck.res,
                                        config='ck.all', .before=1),
                             add_column(mdc.host.pbmc.res,
                                        config='pbmc.all', .before=1),
                             .id = 'n') %>%
                   mutate(n = case_when(as.numeric(n)%in%c(1,3) ~ 'n=97',
                                        as.numeric(n)%in%c(2,4) ~ 'n=25')) %>%
                   dplyr::filter(Covariate == 'OSCI_Score_Samp' & str_detect(
                       Confounding,'LD|SD|P_|S_|O_|SDNA_|OP_|Worst')) %>%
                   dplyr::filter(FDR <= 0.05) %>%
                   select(config, n, everything())),
    
    # all osci-associated features regardless of confounder status
    # WITH microbiota features as covariates
    tar_target(osci.metab.feats,
               add_column(mdc.metabs.motus.all.res,
                          config='metabs.motus.all', .before=1) %>%
                   bind_rows(
                       add_column(mdc.metabs.motus.early.res,
                                  config='metabs.motus.early', .before=1),
                       add_column(mdc.metabs.motus.late.res,
                                  config='metabs.motus.late', .before=1),
                       .id = 'n') %>%
                   mutate(n = case_when(as.numeric(n)%in%c(1,3) ~ 'n=97',
                                        as.numeric(n)%in%c(2,4) ~ 'n=25')) %>%
                   dplyr::filter(Covariate == 'OSCI_Score_Samp' & str_detect(
                       Confounding,'LD|SD|S_|O_|SDNA_|OP_|Worst')) %>%
                   dplyr::filter(FDR <= 0.05) %>%
                   select(config, n, everything())),
    
    ###################################################################
    ## MANUSCRIPT PLOTS ----
    ## statistical comparison group definition ----
    tar_target(grp.time.comp, list(c('Control','EarlyMild'), 
                                 c('Control','EarlySevere'),
                                 c('Control','LateMild'), 
                                 c('Control','LateSevere'),
                                 c('EarlyMild','LateMild'), 
                                 c('EarlySevere','LateSevere'),
                                 c('EarlyMild','LateSevere'), 
                                 c('EarlySevere','LateMild'))),
    tar_target(grp.comp, list(c('Control','Mild'),
                              c('Control','Severe'),
                              c('Mild','Severe'))),
    tar_target(grp.metabs.comp, list(c('P_Control','P_Mild'),
                                     c('P_Control','P_Severe'),
                                     c('P_Mild','P_Severe'),
                                     c('U_Control','U_Mild'),
                                     c('U_Control','U_Severe'),
                                     c('U_Mild','U_Severe'))),
    # FIG 2AB -- alpha diversity metric scatterplots ----
    tar_target(sdna.alpha.plot, plot.alpha.diversity(sdna.alpha, sdna.meta)),
    tar_target(op.alpha.plot, plot.alpha.diversity(op.alpha, op.meta)),
    tar_target(alpha.panel.plot, (sdna.alpha.plot | op.alpha.plot) +
                   plot_layout(guides = 'collect') &
                   theme(legend.position = 'right')),
    
    # FIG 2CD -- beta diversity PCoAs ----
    tar_target(sdna.beta.plot, plot.beta.diversity(sdna.feat.ra, sdna.meta,
                                                   color = 'full.range')),
    tar_target(op.beta.plot, plot.beta.diversity(op.feat.ra, op.meta,
                                                 color = 'full.range')),
    tar_target(beta.panel.plot, (sdna.beta.plot | op.beta.plot) +
                   plot_layout(guides = 'collect') &
                   theme(legend.position = 'right')),
    
    # FIG 2EF -- comparative effect size scatterplots (microbiome) ----
            # patient samples only
    tar_target(sdna.genus.es.plot, plot.comparative.es(mdc.sdna.genus.pat.res,
                                                       svar = 'AbxCurr')),
    tar_target(op.genus.es.plot, plot.comparative.es(mdc.op.genus.pat.res,
                                                     svar = 'AbxCurr')),
    tar_target(genus.es.panel.plot, (sdna.genus.es.plot | op.genus.es.plot) +
                   plot_layout(guides = 'collect') &
                   theme(legend.position = 'right')),
    
            # all samples
    tar_target(sdna.genus.all.es.plot, plot.comparative.es(mdc.sdna.genus.all.res,
                                                       svar = 'Antibiotics')),
    tar_target(op.genus.all.es.plot, plot.comparative.es(mdc.op.genus.all.res,
                                                     svar = 'Antibiotics')),
    tar_target(genus.all.es.panel.plot, (sdna.genus.all.es.plot | op.genus.all.es.plot) +
                   plot_layout(guides = 'collect') &
                   theme(legend.position = 'right')),
    
    tar_target(sdna.motus.es.plot, plot.comparative.es(mdc.sdna.motus.pat.res,
                                                       svar = 'AbxCurr')),
    tar_target(op.motus.es.plot, plot.comparative.es(mdc.op.motus.pat.res,
                                                     svar = 'AbxCurr')),
    tar_target(motus.es.panel.plot, (sdna.motus.es.plot | op.motus.es.plot) +
                   plot_layout(guides = 'collect') &
                   theme(legend.position = 'right')),
    
    ## FIG 2 -- MICROBIOTA: combine 2AB 2CD 2EF ----
    tar_target(combined.microbiota.layout, 
               layout <- "
               AAAAAAA
               AAAAAAA
               BBBBBBB
               BBBBBBB
               CCCCCCC
               CCCCCCC
               CCCCCCC"),
    tar_target(combined.microbiota.fig, 
               (alpha.panel.plot / beta.panel.plot / 
                    (genus.all.es.panel.plot & theme(legend.position='bottom'))) +
                   # plot_layout(design = combined.microbiota.layout) +
                   plot_annotation(tag_levels = 'a')),
    tar_target(combined.microbiota.fig.motus, 
               alpha.panel.plot / beta.panel.plot / motus.es.panel.plot),
    
    tar_target(combined.microbiota.fig.alt, 
               ((sdna.alpha.plot | sdna.beta.plot | sdna.genus.es.plot) /
                    (op.alpha.plot | op.beta.plot | op.genus.es.plot)) +
                   # plot_layout(widths = c(3,2)) +
                   plot_annotation(tag_levels = 'a')),
    
    # FIG 3AH -- cytokine boxplots by OSCI and timepoint ----
    tar_target(s.cytokines.plots, cytokines %>%
                   select(-contains('GAPDH')) %>%
                   left_join(select(clinical.feat, GenID, IL6),
                             by = c('UniqueID'='GenID')) %>%
                   distinct() %>%
                   left_join(select(patient.data, ShortID, contains('OSCI')),
                             by = c('PatientID'='ShortID')) %>%
                   mutate(IL6 = as.double(IL6)) %>%
                   pivot_longer(cols = c(IFNa,IFNg,IL28,IL6,IL10,IP10,MCP1,TNFa),
                                names_to = 'name', values_to = 'val') %>%
                   filter(!is.na(val) & !is.na(Days)) %>%
                   mutate(print = fct_recode(name,
                                             IFNα = 'IFNa',
                                             IFNγ = 'IFNg',
                                             `IFNλ2` = 'IL28',
                                             TNFα = 'TNFa',
                                             `IL-6` = 'IL6',
                                             `IL-10` = 'IL10',
                                             CCL2 = 'MCP1',
                                             `IP-10` = 'IP10')) %>%
                   mutate(name = factor(name,
                                        levels = c('IFNa','IFNg','IL28',
                                                   'TNFa','IL6','IL10',
                                                   'MCP1','IP10'))) %>%
                   plot.ck(comparisons = grp.time.comp)),
    # FIG IJ -- not computed by me ----
    # FIG 3K -- PBMC boxplots by OSCI ----
    tar_target(pbmc.plots, combined.pbmc %>%
                   select(-contains('_'), -ISG) %>%
                   left_join(select(clinical.feat, GenID),
                             by=c('UniqueID'='GenID')) %>%
                   distinct() %>%
                   mutate(ShortID = str_split(UniqueID, '-') %>% 
                              map_chr(~ head(.,1))) %>%
                   left_join(select(patient.data, ShortID, contains('OSCI')),
                             by = 'ShortID') %>%
                   pivot_longer(cols = c(cMono, ncMono, mDC, pDC, NK, pNK, MK),
                                names_to = 'name', values_to = 'val') %>%
                   mutate(name = factor(name, levels = c('cMono', 'ncMono',
                                                         'mDC', 'pDC',
                                                         'NK', 'pNK', 'MK'))) %>%
                   filter(!is.na(val)) %>%
                   plot.pbmc(comparisons = grp.comp)),
    ## FIG 3 -- HOST IMMUNE RESPONSE: combine panels manually ----
    ## FIG 4: METABOLOME: tryptophan and bile metabolites for pathways ---- 
    tar_target(tryp.plot.dat, tryp.metabs %>%
                   discretize.days() %>%
                   mutate(Group = paste0(Sample,'_',OSCI_Class_Samp)) %>%
                   mutate(Group = fct_relevel(Group, c('P_Control','P_Mild','P_Severe',
                                                       'U_Control','U_Mild','U_Severe'))) %>%
                   mutate(OSCI_Class_Samp = factor(
                       OSCI_Class_Samp, levels = c('Control',
                                                   'Mild',
                                                   'Severe'))) %>%
                   mutate(Sample = factor(Sample, levels=c('Plasma','Urine')))),
    tar_target(tryp.plotlist, tryp.plot.dat %>%
                   nest(data = -Compound) %>%
                   mutate(plots = map2(
                       data, Compound, 
                       ~ plot.metab.boxplots(.x, .y, grp.metabs.comp))) ),
    tar_target(tryp.plots, wrap_plots(tryp.plotlist$plots,
                                      ncol=3)),
    tar_target(tryp.conf.early.stat, metabs.clean.res %>%
                   filter(Compound %in% tryp.plot.dat$Compound) %>%
                   filter(Analysis == 'Early') %>%
                   select(Sample, Compound, ConfStatus)),
    tar_target(tryp.conf.late.stat, metabs.clean.res %>%
                   filter(Compound %in% tryp.plot.dat$Compound) %>%
                   filter(Analysis == 'Late') %>%
                   select(Sample, Compound, ConfStatus)),
    # plasma only
    tar_target(p.tryp.plot.dat, tryp.metabs %>%
                   filter(Sample == 'Plasma') %>%
                   discretize.days() %>%
                   mutate(OSCI_Class_Samp = factor(
                       OSCI_Class_Samp, levels = c('Control',
                                                   'Mild',
                                                   'Severe'))) %>%
                   mutate(Group = OSCI_Class_Samp) %>%
                   filter(!(UniqueID=='P19-V2' & Compound=='Melatonin'))),
    tar_target(p.tryp.plotlist, p.tryp.plot.dat %>%
                   nest(data = -Compound) %>%
                   mutate(plots = map2(
                       data, Compound,
                       ~ plot.metab.boxplots(.x, .y, grp.comp,
                                             side = 'left')))),
    tar_target(p.tryp.plots, wrap_plots(p.tryp.plotlist$plots,
                                        ncol=4)),
    # urine
    tar_target(u.tryp.plot.dat, tryp.metabs %>%
                   filter(Sample == 'Urine') %>%
                   discretize.days() %>%
                   mutate(OSCI_Class_Samp = factor(
                       OSCI_Class_Samp, levels = c('Control',
                                                   'Mild',
                                                   'Severe'))) %>%
                   mutate(Group = OSCI_Class_Samp)),
    tar_target(u.tryp.plotlist, u.tryp.plot.dat %>%
                   nest(data = -Compound) %>%
                   mutate(plots = map2(
                       data, Compound,
                       ~ plot.metab.boxplots(.x, .y, grp.comp,
                                             side = 'right')))),
    tar_target(u.tryp.plots, wrap_plots(u.tryp.plotlist$plots,
                                        ncol=4)),
    
    # bile acids (plasma only)
    tar_target(bile.plot.dat, p.metabolites %>%
                   filter(str_detect(Compound, 'CA|Glycine|Taurine')) %>%
                   select(UniqueID, Compound, Value, Days) %>%
                   rename(Val = 'Value') %>%
                   left_join(select(clinical.feat,
                                    GenID, contains('OSCI')),
                             by=c('UniqueID'='GenID')) %>%
                   # mutate(Group = ifelse(OSCI_Class_Samp == 'Control', 'Control',
                   #                   paste0(Days, OSCI_Class_Samp))) %>%
                   # mutate(Group = fct_relevel(Group, c('Control','EarlyMild','EarlySevere',
                   #                             'LateMild','LateSevere'))) %>%
                   mutate(Group = OSCI_Class_Samp) %>%
                   distinct() ),
    tar_target(p.bile.plots, bile.plot.dat %>%
                   nest(data = -Compound) %>%
                   mutate(plots = map2(
                       data, Compound, 
                       ~ plot.metab.boxplots(.x, .y, grp.comp))) %>%
                   pull(plots) %>%
                   wrap_plots()),
    
    # FIG 5A -- confounder status breakdown, all spaces ----
    tar_target(radar.input.abs, list(
        'sdna.all' = prep.radar.stats(
            mdc.sdna.genus.all.res, 'Stool',
            pct = FALSE),
        'op.all' = prep.radar.stats(
            mdc.op.genus.all.res, 'Oropharyngeal',
            pct = FALSE),
        'plasma.all' = prep.radar.stats(
            dplyr::filter(mdc.metabs.all.res,
                          str_detect(feature,'P_')), 'Plasma',
            pct = FALSE),
        'urine.all' = prep.radar.stats(
            dplyr::filter(mdc.metabs.all.res,
                          str_detect(feature,'U_')), 'Urine',
            pct = FALSE),
        'ck.all' = prep.radar.stats(
            bind_rows(mdc.host.ck.basic.res,
                      filter(mdc.host.pbmc.basic.res,
                             !str_detect(feature, '_'))),
            'Host Immune System', pct = FALSE)) %>%
            reduce(full_join, by = 'assoc.type') %>%
            mutate_if(is.double, ~ replace_na(.,0))),
    tar_target(rel.status.totals, radar.input.abs %>%
                   mutate(assoc.status = case_when(
                       assoc.type == 'Deconfounded' ~ 'Deconfounded',
                       assoc.type == 'NS' ~ 'Not Significant',
                       assoc.type == 'total' ~ 'Total',
                       TRUE ~ 'Confounded')) %>%
                   group_by(assoc.status) %>%
                   summarize_if(is.double, sum) %>%
                   gather('space','total.n',-assoc.status) %>%
                   filter(assoc.status == 'Total') %>%
                   select(space, total.n)),
    tar_target(rel.status.props, radar.input.abs %>%
                   mutate(assoc.status = case_when(
                       assoc.type == 'Deconfounded' ~ 'Deconfounded',
                       assoc.type == 'NS' ~ 'Not Significant',
                       assoc.type == 'total' ~ 'Total',
                       TRUE ~ 'Confounded')) %>%
                   group_by(assoc.status) %>%
                   summarize_if(is.double, sum) %>%
                   filter(assoc.status != 'Total') %>%
                   gather('space','n',-assoc.status) %>%
                   left_join(rel.status.totals) %>%
                   # mutate(nsig = total.n - n) %>%
                   # mutate(nconf = (total.n - n - nsig))
                   mutate(n = ifelse(
                       assoc.status == 'Confounded',
                       NA, n)) %>%
                   nest(data = -space) %>%
                   mutate(data = map(data, function(.x) {
                       mutate(.x, n = ifelse(
                           assoc.status=='Confounded',
                           .x$total.n - sum(.x$n, na.rm = TRUE),
                           n)) })) %>%
                   unnest(data) %>%
                   mutate(prop = n/total.n)),
    tar_target(rel.status.plot, 
               plot.confounder.status.comp(rel.status.props)),
    
    # FIG 5B -- radar plots of significant associations, all spaces ----
    tar_target(radar.input.pct, list(
        'sdna.all' = prep.radar.stats(
            mdc.sdna.genus.all.res, 'Stool'),
        'op.all' = prep.radar.stats(
            mdc.op.genus.all.res, 'Oropharyngeal'),
        'plasma.all' = prep.radar.stats(
            dplyr::filter(mdc.metabs.all.res,
                          str_detect(feature,'P_')), 'Plasma'),
        'urine.all' = prep.radar.stats(
            dplyr::filter(mdc.metabs.all.res,
                          str_detect(feature,'U_')), 'Urine'),
        'ck.all' = prep.radar.stats(
            bind_rows(mdc.host.ck.basic.res,
                      filter(mdc.host.pbmc.basic.res,
                             !str_detect(feature, '_'))),
            'Host Immune System')) %>%
            reduce(full_join, by = 'assoc.type') %>%
            mutate_if(is.double, ~ replace_na(.,0))),
    tar_target(radar.plots.pct.ns, radar.input.pct %>%
                   gather('Timepoint','Pct',-assoc.type) %>%
                   spread('assoc.type','Pct') %>%
                   plot.radar(slice.type = 'spaces.combo',
                              incl.ns = TRUE,
                              plot.title = '',
                              mid = 35,
                              max = 70)),
    tar_target(radar.plots.pct, radar.input.pct %>%
                   gather('Timepoint','Pct',-assoc.type) %>%
                   spread('assoc.type','Pct') %>%
                   plot.radar(slice.type = 'spaces.combo',
                              incl.ns = FALSE,
                              plot.title = '',
                              mid = 20,
                              max = 40)),
    # FIG 5C -- arc diagram, key spaces ----
    tar_target(arcdiag.input, prep.network(
        mdc.sdna.genus.all.metab.res,
        mdc.metabs.motus.all.res,
        mdc.host.ck.res,
        mdc.host.pbmc.res,
        list('metabolite.categ' = metabolite.categ,
             'sdna.genus.prev' = sdna.genus.prev,
             'sdna.genus.mean.ra' = sdna.genus.mean.ra))),
    tar_target(arcdiag.osci.feats, osci.sdna.feats %>%
                   filter(config == 'genus.all') %>%
                   bind_rows(osci.host.feats %>%
                                 filter(config == 'ck.all')) %>%
                   bind_rows(osci.metab.feats %>%
                                 filter(config == 'metabs.motus.all')) %>%
                   filter(str_detect(Confounding, 'OSCI|SD|LD|NC')) %>%
                   filter(!str_detect(feature, 'U_')) %>%
                   mutate(feature = str_remove_all(feature, 'P_'))),
    tar_target(arcdiag.plot.dat, 
               get.network.dat(arcdiag.input)),
    tar_target(arcdiag.plot, plot.arcdiagram(
        arcdiag.plot.dat, arcdiag.osci.feats)),
    
    ## FIG 5: -OMICS INTEGRATION: combine 5A 5B 5C + manual 5D ----
    
    ## SFIG 1 -- metadeconfoundR heatmap for genus-level microbiota ----
    tar_target(mdc.genus.all.pat.plot, 
               plot.genus.heatmap(
                   mdc.sdna.genus.pat.res, mdc.op.genus.pat.res,
                   mdc.sdna.genus.pat.feat, mdc.op.genus.pat.feat,
                   sdna.genus.stats, op.genus.stats,
                   all = FALSE)),
    # used in actual supplement:
    tar_target(mdc.genus.all.plot, 
               plot.genus.heatmap(
                   mdc.sdna.genus.all.res, mdc.op.genus.all.res,
                   mdc.sdna.genus.all.feat, mdc.op.genus.all.feat,
                   sdna.genus.stats, op.genus.stats,
                   all = TRUE)),
    
    ## SFIG 2 -- OP and TBS main genus-level composition ----
    tar_target(vap.op.comp.plot, plot.op.composition(op.genus.feat)),
    tar_target(vap.tbs.comp.plot, plot.tbs.composition(tbs.genus.feat)),
    tar_target(airway.comp.plot, (vap.op.comp.plot + theme(
        legend.position = 'top')) / (vap.tbs.comp.plot + 
                                         theme(legend.position = 'top'))),
    
    ## SFIG 3 -- ISG & cytokine boxplots by OSCI and/or timepoint ----
    tar_target(isg.ck.pbmc.data, cytokines %>%
                   select(UniqueID, PatientID, Days,
                          contains('GAPDH'), IFNa, IFNg, IL28) %>%
                   # full_join(select(s.cytokines, UniqueID, IFNa, IFNg, IL28),
                   #           by='UniqueID') %>%
                   full_join(select(combined.pbmc, UniqueID, ISG), 
                             by='UniqueID') %>%
                   left_join(select(patient.data, ShortID, contains('OSCI')),
                             by = c('PatientID'='ShortID'))),
    tar_target(isg.plot, 
               plot.isg.ckpbmc(isg.ck.pbmc.data,
                               comparisons = grp.comp,
                               type = 'ISG')),
    tar_target(isg.ifn.plots, 
               plot.isg.ckpbmc(isg.ck.pbmc.data,
                               type = 'ISG.IFN')),
    tar_target(ifn.gapdh.plots, 
               plot.isg.ckpbmc(isg.ck.pbmc.data,
                               comparisons = grp.time.comp,
                               type = 'IFN.GAPDH')),
    tar_target(isg.combo.plot, (((isg.plot | ifn.gapdh.plots) +
                                    plot_layout(widths = c(1,3))) / isg.ifn.plots) +
                   plot_annotation(tag_levels = 'a')),
    
    # SFIG 4A -- plasma comparative effect size scatterplot + heatmap ----
    tar_target(p.metab.time.es.plot, 
               plot.metab.comparative.es(mdc.metabs.early.res, 
                                         mdc.metabs.late.res, 
                                         metabolite.categ,
                                         cov = 'OSCI_Score_Samp',
                                         site = 'P_')),
    tar_target(p.chem.plot.input.orig, mdc.p.chemclass.mild %>%
                   bind_rows(mdc.p.chemclass.severe, .id = 'Group') %>%
                   mutate(Group = fct_recode(Group,
                                             Mild = '1',
                                             Severe = '2'))),
    tar_target(p.chem.plot.tmp1, mdc.p.chemclass.mild.early.res %>%
                   filter(Covariate=='COVID') %>%
                   bind_rows(mdc.p.chemclass.mild.late.res %>%
                                 filter(Covariate=='COVID'),
                             .id='Group') %>%
                   mutate(Group = fct_recode(Group,
                                             Mild.Early='1',
                                             Mild.Late='2'))),
    tar_target(p.chem.plot.tmp2, mdc.p.chemclass.severe.early.res %>%
                   filter(Covariate=='COVID') %>%
                   bind_rows(mdc.p.chemclass.severe.late.res %>%
                                 filter(Covariate=='COVID'),
                             .id='Group') %>%
                   mutate(Group = fct_recode(Group,
                                             Sev.Early='1',
                                             Sev.Late='2'))),
    tar_target(p.chem.plot.input, p.chem.plot.tmp1 %>%
                   bind_rows(p.chem.plot.tmp2)),
    tar_target(p.chemclass.plot, 
               plot.chemclass.heatmap(p.chem.plot.input)),
    
    # SFIG 4B -- urine comparative effect size scatterplot + heatmap ----
    tar_target(u.metab.time.es.plot, 
               plot.metab.comparative.es(mdc.metabs.early.res, 
                                         mdc.metabs.late.res, 
                                         metabolite.categ,
                                         cov = 'OSCI_Score_Samp',
                                         site = 'U_')),
    tar_target(u.chem.plot.input.orig, mdc.u.chemclass.mild %>%
                   bind_rows(mdc.u.chemclass.severe, .id = 'Group') %>%
                   mutate(Group = fct_recode(Group,
                                             Mild = '1',
                                             Severe = '2'))),
    tar_target(u.chem.plot.tmp1, mdc.u.chemclass.mild.early.res %>%
                   filter(Covariate=='COVID') %>%
                   bind_rows(mdc.u.chemclass.mild.late.res %>%
                                 filter(Covariate=='COVID'),
                             .id='Group') %>%
                   mutate(Group = fct_recode(Group,
                                             Mild.Early='1',
                                             Mild.Late='2'))),
    tar_target(u.chem.plot.tmp2, mdc.u.chemclass.severe.early.res %>%
                   filter(Covariate=='COVID') %>%
                   bind_rows(mdc.u.chemclass.severe.late.res %>%
                                 filter(Covariate=='COVID'),
                             .id='Group') %>%
                   mutate(Group = fct_recode(Group,
                                             Sev.Early='1',
                                             Sev.Late='2'))),
    tar_target(u.chem.plot.input, u.chem.plot.tmp1 %>%
                   bind_rows(u.chem.plot.tmp2)),
    tar_target(u.chemclass.plot,
               plot.chemclass.heatmap(u.chem.plot.input)),
    
    ## SFIG 4: PLASMA METABOLOME: combine 4A 4B ----
    # tar_target(combined.plasma.fig, (
    #     p.metab.time.es.plot | p.metab.conf.barplot) +
    #         plot_annotation(tag_levels = 'a') +
    #         plot_layout(widths = c(3,1))),
    tar_target(combined.plasma.fig.alt, (
        p.metab.time.es.plot | p.chemclass.plot) +
            plot_annotation(tag_levels = 'a') +
            plot_layout(widths = c(3,1))),
    tar_target(combined.urine.fig.alt, (
        u.metab.time.es.plot | u.chemclass.plot) +
            plot_annotation(tag_levels = 'a') +
            plot_layout(widths = c(3,1))),
    tar_target(combined.metab.fig,
               combined.plasma.fig.alt / combined.urine.fig.alt),
    ## SFIG 5 ---
    ###################################################################
    ## MANUSCRIPT TABLES AND DATA ----
    
    ###################################################################
    ## PLOTS THAT DIDN'T MAKE THE CUT ----
    ## metabolite volcano plots ----
    tar_target(osci.metplot.early,
               plot.metabolite.volcanos(mdc.metabs.early.res,
                                        metabolite.categ,
                                        'OSCI_Score_Samp',
                                        'Early samples, N=56')),
    tar_target(osci.metplot.late,
               plot.metabolite.volcanos(mdc.metabs.late.res,
                                        metabolite.categ,
                                        'OSCI_Score_Samp',
                                        'Late samples, N=70')),
    tar_target(covid.metplot.early,
               plot.metabolite.volcanos(mdc.metabs.early.res,
                                        metabolite.categ,
                                        'COVID',
                                        'Early samples, N=56',
                                        FDR.lim = 3.5)),
    tar_target(covid.metplot.late,
               plot.metabolite.volcanos(mdc.metabs.late.res,
                                        metabolite.categ,
                                        'COVID',
                                        'Late samples, N=70',
                                        FDR.lim = 2.5)),
    ## all metabolites heatmap ----
    tar_target(mdc.metabs.all.plot, 
               plot.metabolite.heatmap(mdc.metabs.all.res,
                                       metabolite.categ)),
    tar_target(mdc.metabs.early.plot, 
               mdc.metabs.early.res %>%
                   filter(!str_detect(Covariate, 
                                      'Pneumonia|Bacteremia|DaysHosp')) %>%
                   plot.metabolite.heatmap(metabolite.categ)),
    
    ## metabolite correlations with one another ----
    tar_target(metab.heat.1, metab.corr.filt %>% 
                   filter(str_detect(
                       term,'P_PC|P_TG|P_Cer|P_SM|P_FA|P_DG')) %>%
                   filter(str_detect(
                       term2.class,'aromatic|indole')|str_detect(
                           term2,'rypt|ynur|ndole')) %>% 
                   filter(str_detect(term2,'P_')) %>% 
                   select(term2,term,rho,FDR) %>% 
                   mutate(Stars = as.character(
                       cut(.$FDR,
                           breaks = c(-Inf, 0.001, 0.01, 0.1, Inf),
                           label = c("***", "**", "*", "")))) %>%
                   ggplot(aes(x=term2,y=term, fill=rho)) +
                   geom_tile() +
                   scale_fill_gradient2(limits = c(-1,1)) +
                   geom_text(aes(label = Stars), nudge_y = -0.5, size = 9) + 
                   ggembl::theme_publication() +
                   theme(axis.text.x = element_text(
                       angle = -45, hjust = 0, size = 9),
                       axis.text.y = element_text(size = 9),
                       axis.ticks.y = element_blank()) +
                   labs(x='',y='') ),
    ## metadeconfoundR heatmap for motus ----
    tar_target(mdc.motus.all.plot, 
               plot.motus.heatmap(
                   mdc.sdna.motus.pat.res,
                   mdc.op.motus.pat.res)),
    ## main confounding factors (metabolites) barplot ----
    # tar_target(p.metabs.bar.input, list(
    #     'early' = prep.radar.plot(mdc.metabs.early.res %>%
    #                                   filter(str_detect(feature,'P_')),
    #                               'Early', incl.ns = TRUE, sum = TRUE),
    #     'late' = prep.radar.plot(mdc.metabs.late.res %>%
    #                                  filter(str_detect(feature,'P_')),
    #                              'Late', incl.ns = TRUE, sum = TRUE)) %>%
    #         reduce(full_join, by = 'Confounding') %>%
    #         mutate_if(is.double, ~ replace_na(.,0)) %>%
    #         gather('Timepoint','Pct',-Confounding) ), 
    # tar_target(p.metab.conf.barplot, 
    #            plot.metab.barplot(p.metabs.bar.input)),
    
    ## urine metabolome and confounding factors ----
    # tar_target(u.metab.time.es.plot, 
    #            plot.metab.comparative.es(mdc.metabs.early.res, 
    #                                      mdc.metabs.late.res, 
    #                                      metabolite.categ,
    #                                      cov = 'OSCI_Score_Samp',
    #                                      site = 'U_')),
    # tar_target(u.metabs.bar.input, list(
    #     'early' = prep.radar.plot(mdc.metabs.early.res %>%
    #                                   filter(str_detect(feature,'U_')),
    #                               'Early', incl.ns = TRUE, sum = TRUE),
    #     'late' = prep.radar.plot(mdc.metabs.late.res %>%
    #                                  filter(str_detect(feature,'U_')),
    #                              'Late', incl.ns = TRUE, sum = TRUE)) %>%
    #         reduce(full_join, by = 'Confounding') %>%
    #         mutate_if(is.double, ~ replace_na(.,0)) %>%
    #         gather('Timepoint','Pct',-Confounding) ), 
    # tar_target(u.metab.conf.barplot, plot.metab.barplot(u.metabs.bar.input)),
    # tar_target(combined.urine.fig, (
    #     u.metab.time.es.plot | u.metab.conf.barplot) +
    #         plot_annotation(tag_levels = 'a') +
    #         plot_layout(widths = c(3,1))),
    ## radar plots for metabolic metadeconfoundR ----
    # tar_target(metabs.radar.input.sig, list(
    #     'early' = prep.radar.plot(mdc.metabs.early.res, 'Early'),
    #     'late' = prep.radar.plot(mdc.metabs.late.res, 'Late'))),
    # tar_target(metab.radar.sigonly, metabs.radar.input.sig %>%
    #                reduce(full_join, by = 'Confounding') %>%
    #                mutate_if(is.double, ~ replace_na(.,0)) %>%
    #                gather('Timepoint','Pct',-Confounding) %>%
    #                spread('Confounding','Pct') %>%
    #                plot.radar(slice.type = 'time',
    #                                 plot.title = '% total OSCI-assoc. (P<0.05, n=114)')),
    # tar_target(metabs.radar.input.wns, list(
    #     'early' = prep.radar.plot(mdc.metabs.early.res, 'Early',
    #                               incl.ns = TRUE),
    #     'late' = prep.radar.plot(mdc.metabs.late.res, 'Late',
    #                              incl.ns = TRUE))),
    # tar_target(metab.radar.all, metabs.radar.input.wns %>%
    #                reduce(full_join, by = 'Confounding') %>%
    #                mutate_if(is.double, ~ replace_na(.,0)) %>%
    #                gather('Timepoint','Pct',-Confounding) %>%
    #                spread('Confounding','Pct') %>%
    #                plot.radar(incl.ns = TRUE,
    #                                 slice.type = 'time',
    #                                 mid = 30, max = 60, 
    #                                 plot.title = '% total metabolites (n=256)')),
    # # combine into panel
    # tar_target(metab.radar.combine, (metab.radar.all | metab.radar.sigonly) +
    #                plot_layout(guides = 'collect') &
    #                theme(legend.position = 'bottom')),
    
    ## radar plots for metagenomic metadeconfoundR ----
    # tar_target(microb.radar.sep, (sdna.radar | op.radar) +
    #                plot_layout(guides = 'collect') &
    #                theme(legend.position = 'bottom') &
    #                plot_annotation(title = '% OSCI-associated mOTUs (P<0.05)')),
    ## radar plots for individual spaces (previous supp. fig) ----
    # separating plasma and urine metabolites
    # tar_target(p.met.radar.input.time, list(
    #     'early' = prep.radar.plot(
    #         mdc.metabs.early.res %>%
    #             filter(str_detect(feature,'P_')), 'Early'),
    #     'late' = prep.radar.plot(
    #         mdc.metabs.late.res %>%
    #             filter(str_detect(feature,'P_')), 'Late'))),
    # tar_target(p.met.radar.time, p.met.radar.input.time %>%
    #                reduce(full_join, by = 'Confounding') %>%
    #                mutate_if(is.double, ~ replace_na(.,0)) %>%
    #                gather('Timepoint','Pct',-Confounding) %>%
    #                spread('Confounding','Pct') %>%
    #                plot.radar(slice.type = 'time',
    #                           mid = 10, max = 100,
    #                           plot.title = 'plasma (n~90)')),
    # 
    # tar_target(u.met.radar.input.time, list(
    #     'early' = prep.radar.plot(
    #         mdc.metabs.early.res %>%
    #             filter(str_detect(feature,'U_')), 'Early'),
    #     'late' = prep.radar.plot(
    #         mdc.metabs.late.res %>%
    #             filter(str_detect(feature,'U_')), 'Late'))),
    # tar_target(u.met.radar.time, u.met.radar.input.time %>%
    #                reduce(full_join, by = 'Confounding') %>%
    #                mutate_if(is.double, ~ replace_na(.,0)) %>%
    #                gather('Timepoint','Pct',-Confounding) %>%
    #                spread('Confounding','Pct') %>%
    #                plot.radar(slice.type = 'time',
    #                           mid = 15, max = 75,
    #                           plot.title = 'urine (n~25)')),
    # 
    # # stool and oro microbiota
    # tar_target(sdna.radar.input.sig, list(
    #     'patients' = prep.radar.plot(mdc.sdna.motus.pat.res, 
    #                                  'Patients'),
    #     'all' = prep.radar.plot(mdc.sdna.motus.all.res, 'All'))),
    # tar_target(sdna.radar, sdna.radar.input.sig %>%
    #                reduce(full_join, by = 'Confounding') %>%
    #                mutate_if(is.double, ~ replace_na(.,0)) %>%
    #                gather('Timepoint','Pct',-Confounding) %>%
    #                spread('Confounding','Pct') %>%
    #                plot.radar(slice.type = 'patient.subset',
    #                           plot.title = 'SDNA')),
    # tar_target(op.radar.input.sig, list(
    #     'patients' = prep.radar.plot(mdc.op.motus.pat.res, 'Patients'),
    #     'all' = prep.radar.plot(mdc.op.motus.all.res, 'All'))),
    # tar_target(op.radar, op.radar.input.sig %>%
    #                reduce(full_join, by = 'Confounding') %>%
    #                mutate_if(is.double, ~ replace_na(.,0)) %>%
    #                gather('Timepoint','Pct',-Confounding) %>%
    #                spread('Confounding','Pct') %>%
    #                plot.radar(slice.type = 'patient.subset',
    #                           plot.title = 'OP', 
    #                           max = 60)),
    # 
    # # combine into panel
    # tar_target(radar.plots, ((sdna.radar | op.radar) / 
    #                              (p.met.radar.time | u.met.radar.time)) +
    #                plot_layout(guides = 'collect') &
    #                theme(legend.position = 'bottom') &
    #                plot_annotation(tag_levels = 'a')),
    ## chord plots for data integration figure ----
    # tar_target(chord.plot.dat.full, prep.osci.chord(osci.sdna.feats,
    #                                                 osci.op.feats,
    #                                                 osci.host.feats,
    #                                                 osci.metab.feats,
    #                                                 restricted = FALSE)),
    # 
    # tar_target(chord.plot.full.dir, plot.osci.chord(
    #     chord.plot.dat.full$data,
    #     chord.plot.dat.full$group.labels,
    #     directed = TRUE,
    #     'chord-osci-full-directed.pdf')),
    # 
    # tar_target(chord.plot.full.nondir, plot.osci.chord(
    #     chord.plot.dat.full$data,
    #     chord.plot.dat.full$group.labels,
    #     directed = FALSE,
    #     'chord-osci-full-nondirected.pdf')),
    # 
    # tar_target(chord.plot.dat.min, prep.osci.chord(osci.sdna.feats,
    #                                                osci.op.feats,
    #                                                osci.host.feats,
    #                                                osci.metab.feats,
    #                                                restricted = TRUE)),
    # 
    # tar_target(chord.plot.min.dir, plot.osci.chord(
    #     chord.plot.dat.min$data,
    #     chord.plot.dat.min$group.labels,
    #     directed = TRUE,
    #     'chord-osci-sel-directed.pdf')),
    # 
    # tar_target(chord.plot.min.nondir, plot.osci.chord(
    #     chord.plot.dat.min$data,
    #     chord.plot.dat.min$group.labels,
    #     directed = FALSE,
    #     'chord-osci-sel-nondirected.pdf')),
    # 
    # tar_target(chord.plot.dat.osci.only, prep.osci.chord(osci.sdna.feats,
    #                                                      osci.op.feats,
    #                                                      osci.host.feats,
    #                                                      osci.metab.feats,
    #                                                      restricted = TRUE,
    #                                                      osci.only = TRUE)),
    # 
    # tar_target(chord.plot.osci.only, plot.osci.chord(
    #     chord.plot.dat.osci.only$data,
    #     chord.plot.dat.osci.only$group.labels,
    #     directed = TRUE,
    #     'chord-osci-only-directed.pdf')),
    # 
    # tar_target(chord.plot.osci.only.feedback, plot.osci.chord(
    #     chord.plot.dat.osci.only$data,
    #     chord.plot.dat.osci.only$group.labels,
    #     directed = TRUE,
    #     'chord-osci-feedback-directed.pdf')),
    ###################################################################
    ## METABOLITE MIXED MODELS ----
    # master data table w/ both plasma & urine metabolites ----
    tar_target(p.met.dat, p.metabolites.wide %>%
                   # tidyimpute::impute_median_if(is.double) %>%
                   create.master.table(clinical.feat,
                                       patient.data,
                                       metagenomic = FALSE,
                                       patients.only = FALSE,
                                       site = 'P')),
    tar_target(u.met.dat, u.metabolites.wide %>%
                   # tidyimpute::impute_median_if(is.double) %>%
                   create.master.table(clinical.feat,
                                       patient.data,
                                       metagenomic = FALSE,
                                       patients.only = FALSE,
                                       site = 'U')),
    tar_target(met.dat, left_join(p.met.dat, u.met.dat) %>%
                   select(-contains(c('13C','_D3','_D5')))),
    tar_target(p.met.dat.pat, filter(p.met.dat, !str_detect(PatientID,'K'))),
    tar_target(u.met.dat.pat, filter(u.met.dat, !str_detect(PatientID,'K'))),
    tar_target(met.dat.pat, filter(met.dat, !str_detect(PatientID,'K'))),
    tar_target(p.mets, p.metabolites.wide[,-1] %>% colnames()),
    tar_target(u.mets, u.metabolites.wide[,-1] %>% colnames()),
    tar_target(all.mets, met.dat %>% select(starts_with(c('U_','P_'))) %>%
                   colnames()),
    # plasma & urine -- all samples --  OSCI, CCI, Meds, Days ----
    tar_target(met.models.1, all.mets %>%
                   tibble(feat = .) %>%
                   mutate(mods = map(feat, function(feat) {
                       form <- as.formula(
                           paste0('rank(',feat,
                                  ') ~ OSCI_Class+CCI+nMeds+Days+(1|PatientID)'))
                       m <- lmerTest::lmer(form,
                                           data = met.dat,
                                           REML = FALSE) }))),
    tar_target(met.stats.1, model.diagnostics(met.models.1)),
    tar_target(met.reg.1, clean.model.coefs(met.stats.1) %>%
                   mutate(R_Compound = str_split(feat, '_') %>%
                              map_chr(~ tail(., 1))) %>%
                   left_join(select(metabolite.categ,
                                    R_Compound, Compound, 
                                    ChemClass, Pathway))),
    tar_target(met.volcanos.1, basic.volcano.plots(met.reg.1,NULL,NULL,
                                                   metagenomic = FALSE)),
    tar_target(p.met.volcanos.1, met.reg.1 %>%
                   filter(str_detect(feat,'P_')) %>%
                   select(-feat) %>% dplyr::rename(feat = 'R_Compound') %>%
                   basic.volcano.plots(NULL, NULL, metagenomic = FALSE,
                                       labels = met.of.interest$Compound)),
    tar_target(u.met.volcanos.1, met.reg.1 %>%
                   filter(str_detect(feat, 'U_')) %>%
                   select(-feat) %>% dplyr::rename(feat = 'R_Compound') %>%
                   basic.volcano.plots(NULL, NULL, metagenomic = FALSE,
                                       labels = met.of.interest$Compound)),
    
    ### METAGENOMIC MIXED MODELS ----
    tar_target(sdna.dat, create.master.table(sdna.genus.feat,
                                             sdna.clinical.feat,
                                             patient.data,
                                             metagenomic = TRUE,
                                             patients.only = FALSE)),
    tar_target(sdna.pat.dat, filter(sdna.dat, !str_detect(PatientID,'K'))),
    tar_target(genera, sdna.genus.feat %>%
                   arrange(UniqueID) %>%
                   select(-UniqueID, -Days) %>%
                   colnames()),
    
    tar_target(op.dat, create.master.table(op.genus.feat,
                                             op.clinical.feat,
                                             patient.data,
                                             metagenomic = TRUE,
                                             patients.only = FALSE)),
    tar_target(op.genera, op.genus.feat %>%
                   arrange(UniqueID) %>%
                   select(-UniqueID, -Days) %>%
                   colnames()),
    
    
    # SDNA -- patients only -- OSCI, CCI, Meds, Days, Abx, GI ----
    tar_target(sdna.models.1, genera %>%
                   tibble(feat = .) %>%
                   mutate(mods = map(feat, function(feat) {
                       form <- as.formula(
                           paste0('rank(',feat,
                                  ') ~ OSCI_Class+CCI+nMeds+Days+Abx+GI+(1|PatientID)'))
                       m <- lmerTest::lmer(form, 
                                           data = sdna.pat.dat,
                                           REML = FALSE) }))),
    tar_target(sdna.stats.1, model.diagnostics(sdna.models.1)),
    tar_target(sdna.reg.1, clean.model.coefs(sdna.stats.1)),
    tar_target(sdna.volcanos.1, basic.volcano.plots(sdna.reg.1,
                                                    sdna.genus.pat.prev,
                                                    sdna.genus.pat.mean.ra)),
    
    # SDNA -- all samples --  OSCI, CCI, Meds, Days, Abx ----
    tar_target(sdna.models.2, genera %>%
                   tibble(feat = .) %>%
                   mutate(mods = map(feat, function(feat) {
                       form <- as.formula(
                           paste0('rank(',feat,
                                  ') ~ OSCI_Class+CCI+nMeds+Days+Abx+(1|PatientID)'))
                       m <- lmerTest::lmer(form, 
                                           data = sdna.dat,
                                           REML = FALSE) }))),
    tar_target(sdna.stats.2, model.diagnostics(sdna.models.2)),
    tar_target(sdna.reg.2, clean.model.coefs(sdna.stats.2)),
    tar_target(sdna.volcanos.2, basic.volcano.plots(sdna.reg.2,
                                                    sdna.genus.prev,
                                                    sdna.genus.mean.ra)),
    # SDNA -- all samples --  OSCI, Abx ----
    tar_target(sdna.models.simp, genera %>%
                   tibble(feat = .) %>%
                   mutate(mods = map(feat, function(feat) {
                       form <- as.formula(
                           paste0('rank(',feat,
                                  ') ~ OSCI_Class+Abx+(1|PatientID)'))
                       m <- lmerTest::lmer(form, 
                                           data = sdna.dat,
                                           REML = FALSE) }))),
    tar_target(sdna.stats.simp, model.diagnostics(sdna.models.simp)),
    tar_target(sdna.reg.simp, clean.model.coefs(sdna.stats.simp)),
    tar_target(sdna.volcanos.simp, basic.volcano.plots(sdna.reg.simp,
                                                    sdna.genus.prev,
                                                    sdna.genus.mean.ra)),
    # OP -- all samples --  OSCI, Abx ----
    tar_target(op.models.simp, op.genera %>%
                   tibble(feat = .) %>%
                   mutate(mods = map(feat, function(feat) {
                       form <- as.formula(
                           paste0('rank(',feat,
                                  ') ~ OSCI_Class+Abx+(1|PatientID)'))
                       m <- lmerTest::lmer(form, 
                                           data = op.dat,
                                           REML = FALSE) }))),
    tar_target(op.stats.simp, model.diagnostics(op.models.simp)),
    tar_target(op.reg.simp, clean.model.coefs(op.stats.simp)),
    tar_target(op.volcanos.simp, basic.volcano.plots(op.reg.simp,
                                                     op.genus.prev,
                                                     op.genus.mean.ra)),
    
    ###################################################################
    ## RMD TO TRIGGER OUTPUT -----
    tar_render(save.plots, here::here('rmd-notebooks',
                                '2022-03-03-manuscript-plots',
                                'final-plots.Rmd'))
    
)

###################################################################
### INTERACTING WITH TARGETS FROM CONSOLE ----
# tar_outdated()
# tar_visnetwork()
# view(tar_meta(fields=warnings))

## debugging: https://books.ropensci.org/targets/debugging.html
# tar_make(callr_function = NULL, names = 'variable-names-here', shortcut = TRUE)