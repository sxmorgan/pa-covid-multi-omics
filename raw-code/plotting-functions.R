################################################
####### HOST IMMUNE RESPONSE BOXPLOTS ETC ######
################################################
plot.ck <- function(df, comparisons) {
    
    colors <- c('#FFFFFF','#FFCB64','#FF5500')
    
    df %>%
        mutate(x = ifelse(OSCI_Class == 'Control', 'Control',
                          paste0(Days, OSCI_Class))) %>%
        mutate(x = fct_relevel(x, c('Control','EarlyMild','EarlySevere',
                                    'LateMild','LateSevere'))) %>%
        ggplot(aes(x=x, y=val)) +
        geom_boxplot(aes(fill=OSCI_Class), alpha = 0.2, outlier.shape = NA) +
        geom_jitter(aes(fill=OSCI_Class), shape=21) +
        # scale_fill_manual(values = colors) +
        scico::scale_fill_scico_d(palette = 'lapaz', direction=-1) +
        facet_wrap(~ print, scales='free_y', 
                   ncol=4) +
        ggpubr::stat_compare_means(comparisons=comparisons, aes(label=..p.signif..),
                                   method='wilcox.test', hide.ns=TRUE,
                                   tip.length=0, vjust=0.5) +
        scale_y_log10() + labs(x='',y='') + guides(fill='none') +
        ggpubr::theme_pubr() +
        theme(axis.text.x = element_text(angle = -45, hjust = 0, size = 9),
              axis.text.y = element_text(size = 9),
              axis.ticks.y = element_blank())
    
    ## double check stat_compare_means output ----
    # df.stat <- df %>% 
    #     mutate(x = ifelse(OSCI_Class == 'Control', 'Control',
    #                       paste0(Days, OSCI_Class))) %>%
    #     mutate(x = fct_relevel(x, c('Control','EarlyMild','EarlySevere',
    #                                 'LateMild','LateSevere'))) %>%
    #     select(-OSCI_Score,-print,-OSCI_Class) %>%
    #     spread('name','val') %>% 
    #     nest(data = -x)
    # cvar <- 'IFNg' #toggle
    # tmp <- df.stat %>%
    #     mutate(Control = map2_dbl(data, x, function(d, l) {
    #         if (l == 'Control') NA
    #         else {v2 <- filter(df.stat, x == 'Control') %>% 
    #                 pull(data) %>% magrittr::extract2(1) %>% pull(cvar)
    #         wilcox.test(get(cvar, d), v2)$p.value }} )) %>%
    #     mutate(EarlyMild = map2_dbl(data, x, function(d, l) {
    #         if (l == 'EarlyMild') NA
    #         else {v2 <- filter(df.stat, x == 'EarlyMild') %>% 
    #                 pull(data) %>% magrittr::extract2(1) %>% pull(cvar)
    #         wilcox.test(get(cvar,d), v2)$p.value }} )) %>%
    #     mutate(LateMild = map2_dbl(data, x, function(d, l) {
    #         if (l == 'LateMild') NA
    #         else {v2 <- filter(df.stat, x == 'LateMild') %>% 
    #                 pull(data) %>% magrittr::extract2(1) %>% pull(cvar)
    #         wilcox.test(get(cvar,d), v2)$p.value }} )) %>%
    #     mutate(EarlySevere = map2_dbl(data, x, function(d, l) {
    #         if (l == 'EarlySevere') NA
    #         else {v2 <- filter(df.stat, x == 'EarlySevere') %>% 
    #                 pull(data) %>% magrittr::extract2(1) %>% pull(cvar)
    #         wilcox.test(get(cvar,d), v2)$p.value }} )) %>%
    #     mutate(LateSevere = map2_dbl(data, x, function(d, l) {
    #         if (l == 'LateSevere') NA
    #         else {v2 <- filter(df.stat, x == 'LateSevere') %>% 
    #                 pull(data) %>% magrittr::extract2(1) %>% pull(cvar)
    #         wilcox.test(get(cvar,d), v2)$p.value }} ))
}

plot.pbmc <- function(df, comparisons) {
    
    # expressed in ratios, no log_y_scale needed
    df %>%
        mutate(x = fct_relevel(OSCI_Class, c('Control','Mild','Severe'))) %>%
        mutate(val = val*100) %>%
        ggplot(aes(x=x, y=val)) +
        geom_boxplot(aes(fill=x), alpha = 0.4, outlier.shape = NA) +
        geom_jitter(aes(fill=x), shape=21) +
        scico::scale_fill_scico_d(palette = 'lapaz', direction=-1) +
        labs(x='',y='') + guides(fill='none') +
        facet_wrap(~ name, scales='free_y', ncol = 7) +
        ggpubr::stat_compare_means(comparisons=comparisons, aes(label=..p.signif..),
                                   method='wilcox.test', hide.ns=TRUE,
                                   tip.length=0) +
        ggpubr::theme_pubr()+
        theme(axis.text.x = element_text(angle = -45, hjust = 0, size = 9),
              axis.text.y = element_text(size = 9),
              axis.ticks.y = element_blank())
}

plot.isg.ckpbmc <- function(df, comparisons = NULL, type) {
    
    switch(type,
           'ISG' = {
               df %>%
                   filter(!is.na(ISG) & !is.na(OSCI_Class)) %>%
                   mutate(OSCI_Class = fct_relevel(
                       OSCI_Class, c('Control','Mild','Severe'))) %>%
                   ggplot(aes(x=OSCI_Class, y=ISG)) +
                   geom_boxplot(aes(fill=OSCI_Class), alpha = 0.4, outlier.shape = NA) +
                   geom_jitter(aes(fill=OSCI_Class), shape=21) +
                   scico::scale_fill_scico_d(palette = 'lapaz', direction=-1) +
                   labs(x='',y='Mean ISG Score') + guides(fill='none') +
                   ggpubr::stat_compare_means(comparisons=comparisons, aes(label=..p.signif..),
                                              method='wilcox.test', hide.ns=TRUE,
                                              tip.length=0) +
                   ggpubr::theme_pubr()+
                   theme(axis.text.x = element_text(angle = -45, hjust = 0, size = 9),
                         axis.text.y = element_text(size = 9),
                         axis.ticks.y = element_blank())
           },
           'IFN.GAPDH' = {
               df %>%
                   pivot_longer(cols = contains('GAPDH')) %>%
                   filter(!is.na(value)) %>%
                   mutate(x = ifelse(OSCI_Class == 'Control', 'Control',
                                     paste0(Days, OSCI_Class))) %>%
                   mutate(x = fct_relevel(
                       x, c('Control','EarlyMild','EarlySevere',
                                               'LateMild','LateSevere'))) %>%
                   ggplot(aes(x=x, y=value)) +
                   geom_boxplot(aes(fill=OSCI_Class), alpha = 0.2, outlier.shape = NA) +
                   geom_jitter(aes(fill=OSCI_Class), shape=21) +
                   scico::scale_fill_scico_d(palette = 'lapaz', direction=-1) +
                   scale_y_log10() +
                   facet_wrap(~ name, scales='free_y', ncol=3) +
                   # hiding to avoid excessive height when NS values are removed later
                   ggpubr::stat_compare_means(comparisons=comparisons, aes(label=..p.signif..),
                                              method='wilcox.test', hide.ns=TRUE,
                                              tip.length=0, vjust=0.5) +
                   labs(x='',y='[ratio mRNA copies]') + guides(fill='none') +
                   ggpubr::theme_pubr() +
                   theme(axis.text.x = element_text(angle = -45, hjust = 0, size = 9),
                         axis.text.y = element_text(size = 9),
                         axis.ticks.y = element_blank())
           },
           'ISG.IFN' = {
               df %>%
                   pivot_longer(cols = c(IFNa,IFNg,IL28), 'IFN','value') %>%
                   select(-contains('GAPDH'), -Days) %>%
                   filter(!is.na(ISG)) %>%
                   mutate(OSCI_Class = fct_relevel(
                       OSCI_Class, c('Control','Mild','Severe'))) %>%
                   ggplot(aes(x=value, y=ISG)) +
                   geom_smooth(method = 'lm', se = FALSE, color ='black') +
                   geom_point(fill = 'grey', shape=21, size=2.5) +
                   scico::scale_fill_scico_d(palette = 'lapaz', direction=-1) +
                   scale_x_log10() + 
                   ggpubr::stat_cor(method = 'spearman', p.accuracy = 0.001, 
                                    r.accuracy = 0.01, label.y = c(0,0.5,1)) +
                   labs(x='(pg/ml)', y='Mean ISG Score') +
                   facet_wrap(~ IFN, scales='free') +
                   ggpubr::theme_pubr()
           })
    
}

################################################
######### MICROBIOTA FIGURES, HEATMAPS #########
################################################
plot.alpha.diversity <- function(df.alpha, df.meta) {
    
    withr::local_options(.new = list(warn = -1))
    
    # previous draft colors: rose chartreuse lilac
    # colors <- c('#b9a0b4','#cecb76','#ff9888')
    
    alpha.meta <- df.alpha %>%
        left_join(select(df.meta, 
                         UniqueID, GenID, PatientID,
                         contains('OSCI'), contains('Abx'),
                         contains('VL'), Antibiotics, Pneumonia, LibSize)) %>%
        rename(OSCI_Score_Worst = 'OSCI_Score', OSCI_Class_Worst = 'OSCI_Class',
               HAP = 'Pneumonia', AbxCurr = 'Antibiotics') %>%
        mutate(HAP = fct_recode(as.character(HAP), 
                                Yes = '1', No = '0')) %>%
        pivot_longer(c('Richness','Shannon','InvSimp'),'alpha.metric','value') %>%
        mutate(alpha.metric = factor(alpha.metric, levels = c('Richness',
                                                              'InvSimp',
                                                              'Shannon')))
    
    alpha.meta %>%
        ggplot(aes(x = OSCI_Score_Worst, y = value)) + 
        geom_smooth(method = 'lm', aes(color = alpha.metric)) +
        ggpubr::stat_cor(aes(color = alpha.metric),
                         method = 'spearman', p.accuracy = 0.001, r.accuracy = 0.01,
                         label.y = c(0,0.5,1)) +
        # previously: outline HAP patients
        # geom_jitter(data = filter(alpha.meta, HAP == 'No'), size = 2.5, alpha=0.8) + 
        # geom_jitter(aes(fill = alpha.metric), data = filter(alpha.meta, HAP == 'Yes'),
        #             pch=21, size=2.5, colour='black', stroke=1, alpha = 0.8) +
        geom_jitter(aes(fill = alpha.metric), 
                    pch=21, color='black', size=3, alpha=0.8) +
        scale_y_log10() +
        scico::scale_fill_scico_d(palette = 'davos', begin = 0.4, end = 0.9,
                                  direction = -1) +
        scico::scale_color_scico_d(palette = 'davos', begin = 0.4, end = 0.9,
                                   direction = -1) +
        labs(x = 'Worst OSCI Score', 
             y = 'Alpha Diversity', 
             color = 'Metric') +
        guides(fill='none') +
        ggembl::theme_embl()
}

plot.beta.diversity <- function(df.feat, df.meta, 
                                svar = 'OSCI_Class_Worst',
                                color = 'simple') {
    # browser()
    meta.tmp <- df.meta %>%
        rename(OSCI_Score_Worst = 'OSCI_Score', OSCI_Class_Worst = 'OSCI_Class',
               HAP = 'Pneumonia', AbxCurr = 'Antibiotics') %>%
        select(UniqueID, GenID, PatientID,
               contains('OSCI'), contains('Abx'), all_of(svar),
               contains('VL'), HAP, LibSize) %>%
        mutate(HAP = fct_recode(as.character(HAP), 
                                Yes = '1', No = '0')) 

    pco <- df.feat %>%
        select(-c(GenID, Site, Days)) %>%
        pcoa(meta.tmp, svar, incl.abx = TRUE)
    
    # add in marginal boxplots ??
    
    df.plot <- pco$points %>%
        as_tibble() %>%
        mutate(Pt.Type = as_factor(Pt.Type)) %>%
        mutate(OSCI_Score_Worst = ifelse(
            is.na(OSCI_Score_Worst), 
            case_when(
                OSCI_Class_Worst=='Control' ~ 0,
                OSCI_Class_Worst=='Mild' ~ 4,
                OSCI_Class_Worst=='Severe' ~ 7),
            OSCI_Score_Worst)) %>%
        mutate(OSCI_Score_Worst = factor(OSCI_Score_Worst,
                                         levels = c(0,1,2,3,4,5,6,7,8))) %>%
        mutate(OSCI_Class_Worst = factor(OSCI_Class_Worst,
                                         levels = c('Control','Mild','Severe')))
    
    if (color == 'discrete') {
        df.plot %>%
            ggplot(aes(x = V1, y = V2)) +
            geom_path(aes(group = PatientID), color = 'darkgrey', alpha = 0.8) +
            geom_point(aes(shape = Pt.Type, size = Pt.Type, alpha = Pt.Type,
                           color = OSCI_Class_Worst, fill = OSCI_Class_Worst),
                       pch=21, size=2.5, colour='black') +
            scale_shape_manual(values = c(21,22)) +
            scale_size_manual(values = c(3,5)) +
            scale_alpha_manual(values = c(0.6, 1)) +
            ggsci::scale_color_jama() +
            ggsci::scale_fill_jama() +
            labs(
                shape = '', size = '',
                x = paste0('PCo 1 [', pco$axes[1], '%]'), 
                y = paste0('PCo 2 [', pco$axes[2], '%]')) + 
            ggrepel::geom_label_repel(dplyr::filter(df.plot, Pt.Type == 'Group Mean'),
                                      mapping = aes(label = OSCI_Class_Worst), size = 9/.pt,
                                      min.segment.length = unit(0, 'lines'), force = 5) +
            ggembl::theme_embl() +
            guides(shape = 'none', size = 'none', 
                   color = 'none', fill = 'none', alpha='none')
    } else {
        df.plot %>%
            ggplot(aes(x = V1, y = V2)) +
            geom_path(aes(group = PatientID), color = 'darkgrey', alpha = 0.8) +
            geom_point(data = filter(df.plot, Pt.Type=='Individual'),
                       aes(fill = OSCI_Score_Worst),
                       shape=21, size=3, alpha=0.8, color='black') +
            geom_point(data = filter(df.plot, Pt.Type=='Group Mean'),
                       aes(fill = OSCI_Score_Worst),
                       shape=22, size=5, alpha=1, color='black') +
            scico::scale_fill_scico_d(palette = 'hawaii', direction = 1) +
            labs(
                fill = '',
                x = paste0('PCo 1 [', pco$axes[1], '%]'), 
                y = paste0('PCo 2 [', pco$axes[2], '%]')) + 
            ggrepel::geom_label_repel(dplyr::filter(df.plot, Pt.Type == 'Group Mean'),
                                      mapping = aes(label = OSCI_Class_Worst), size = 9/.pt,
                                      min.segment.length = unit(0, 'lines'), force = 5) +
            ggembl::theme_embl() 
        
        # ggExtra::ggMarginal(plot,
        #                     type = 'histogram',
        #                     margins = 'x',
        #                     size = 5,
        #                     groupColour = TRUE)
    }

}

plot.comparative.es <- function(df.res, svar = 'AbxCurr') {
    
    # hand-picked from first version introduced into the manuscript
    # colors <- c('#C19CB2','#B24745','#79AF97',
    #             '#008DA1','#DF8F44','#374E55',
    #             '#BEBEBE')
    colors <- c(scico::scico(7, palette = 'lajolla'), '#BEBEBE')
    values <- c('OSCI',
                'OSCI & AbxCurr',
                'OSCI & HAP',
                'OSCI, HAP & AbxCurr',
                'HAP',
                'AbxCurr & HAP',
                'AbxCurr',
                'Other or none')
    df.col <- tibble(Group = values, Color = colors)

    con <- df.res %>%
        select(feature, Covariate, Confounding) %>% 
        spread('Covariate','Confounding') %>%
        select(feature, OSCI_Score_Samp, OSCI_Score_Worst, !!rlang::sym(svar), Pneumonia) %>%
        rename(OSCI_s.FDR = 'OSCI_Score_Samp') %>%
        rename(OSCI_w.FDR = 'OSCI_Score_Worst') %>%
        rename(OTH.FDR = svar) %>% # antibiotics 
        rename(HAP.FDR = Pneumonia)
    
    eff.size <- df.res %>% 
        select(feature, Covariate, `Effect Size`) %>% 
        spread('Covariate','Effect Size') %>%
        select(feature, OSCI_Score_Samp, OSCI_Score_Worst, !!rlang::sym(svar), Pneumonia) 
    
    plot.data <- eff.size %>%
        left_join(con) %>%
        mutate(OSCI_s.FDR = ifelse(
            OSCI_s.FDR %in% c('SD','LD','NC'), TRUE, FALSE)) %>%
        mutate(OSCI_w.FDR = ifelse(
            OSCI_w.FDR %in% c('SD','LD','NC'), TRUE, FALSE)) %>%
        mutate(OTH.FDR = ifelse(
            OTH.FDR %in% c('SD','LD','NC'), TRUE, FALSE)) %>%
        mutate(HAP.FDR = ifelse(
            HAP.FDR %in% c('SD','LD','NC'), TRUE, FALSE)) %>%
        mutate(Group = pmap_chr(
            list(OSCI_s.FDR, OTH.FDR, HAP.FDR, OSCI_w.FDR),
            ~ case_when((..1|..4) && ..2 && !..3 ~ 'OSCI & AbxCurr',
                        (..1|..4) && !..2 && !..3 ~ 'OSCI',
                        !(..1|..4) && ..2 && !..3 ~ 'AbxCurr',
                        !(..1|..4) && ..2 && ..3 ~ 'AbxCurr & HAP',
                        !(..1|..4) && !..2 && ..3 ~ 'HAP',
                        (..1|..4) && ..3 && !..2 ~ 'OSCI & HAP',
                        (..1|..4) && ..2 && ..3 ~ 'OSCI, HAP & AbxCurr',
                        !(..1|..4) && !..2 & !..3 ~ 'Other or none'))) %>%
        mutate(Group = fct_relevel(Group, 'Other or none', after = Inf)) %>%
        mutate(alpha = fct_collapse(Group, 
                                    color = setdiff(levels(Group), 'Other or none'),
                                    grey = 'Other or none')) %>%
        left_join(df.col) %>%
        mutate(Group = factor(Group, levels = values)) %>%
        mutate(Color = factor(Color, levels = colors))
    
    
    labels <- plot.data %>%
        filter(!str_detect(
            Group,'Other') | feature%in%c('Pseudomonas',
                                          'Lactococcus')) %>%
        arrange(Group) %>%
        slice(1:30)
        
    
    plot.data %>%
        ggplot(aes(x = OSCI_Score_Samp, y = !!rlang::sym(svar), text=feature)) +
        geom_hline(yintercept = 0, color = 'grey') +
        geom_vline(xintercept = 0, color = 'grey') +
        geom_point(data = filter(plot.data, str_detect(Group,'Other')),
                   aes(fill = Group, alpha=alpha), pch=21, size = 3) +
        geom_point(data = filter(plot.data, !str_detect(Group,'Other')),
                   aes(fill = Group, alpha=alpha), pch=21, color='black', size = 3) +
        scale_alpha_manual(values=c(1, 0.3)) + 
        # patients-only
        # lims(x = c(-0.7, 0.7), y = c(-0.8, 0.8)) +
        # all samples (current MS 2022-11)
        lims(x = c(-1, 0.8), y = c(-1, 0.8)) +
        labs( 
            x = 'Effect Size of OSCI Score',
            y = paste0('Effect Size of Antibiotics'),
            color = 'Significance', fill = 'Significance') +
        scale_color_manual(values = colors, drop = FALSE) +
        scale_fill_manual(values = colors, drop = FALSE) +
        ggrepel::geom_text_repel(data = labels,
                                 mapping = aes(label = feature), size = 8/.pt,
                                 min.segment.length = unit(0, 'lines'), force = 7,
                                 max.overlaps = 20) +
        ggembl::theme_embl() +
        guides(alpha = 'none')
}

plot.genus.heatmap <- function(df.sdna, df.op, sdna.raw, op.raw, 
                               sdna.stats, op.stats, all = TRUE) {
    
    # browser()
    # standardize variable nomenclature
    df.sdna <- df.sdna %>%
        mutate(Covariate = str_replace_all(
            Covariate,'Antibiotics','AbxCurr')) %>%
        mutate(Covariate = str_replace_all(
            Covariate,'Pneumonia','HAP'))
    df.op <- df.op %>%
        mutate(Covariate = str_replace_all(
            Covariate,'Antibiotics','AbxCurr')) %>%
        mutate(Covariate = str_replace_all(
            Covariate,'Pneumonia','HAP'))
    
    # unused currently
    sdna.most.ab <- sdna.stats %>%
        arrange(desc(mean.ra.all)) %>%
        slice(1:10) %>% pull(feature)
    op.most.ab <- op.stats %>%
        arrange(desc(mean.ra.all)) %>%
        slice(1:10) %>% pull(feature)
    sdna.most.prev <- sdna.stats %>%
        arrange(desc(prev.all)) %>%
        slice(1:10) %>% pull(feature)
    op.most.prev <- op.stats %>%
        arrange(desc(prev.all)) %>%
        slice(1:10) %>% pull(feature)
    
    # get all naively-associated OSCI taxa
    sdna.taxa <- df.sdna %>%
        filter(str_detect(Covariate, 'OSCI') & Confounding!='NS') %>%
        pull(feature) %>% unique()
    op.taxa <- df.op %>%
        filter(str_detect(Covariate, 'OSCI') & Confounding!='NS') %>%
        pull(feature) %>% unique()
    
    # filter for robust OSCI associations only if all samples
    # (more significant taxa --> more to visualize)
    if (all) { 
        sdna.taxa <- df.sdna %>%
            filter(str_detect(
                Covariate, 'OSCI') & Confounding %in% c('LD','SD')) %>%
            pull(feature) %>% unique()
        op.taxa <- df.op %>%
            filter(str_detect(
                Covariate, 'OSCI') & Confounding %in% c('LD','SD')) %>%
            pull(feature) %>% unique()
    }
    
    gi <- c('Rothia','Bifidobacterium','Actinomyces','Turicibacter',
            'Coriobact','Lactococcus','Intestinibacter','Parvimonas',
            'Faecalibacterium', #'Eubacterium', # no Eu signal...
            'Suterellaceae',#'Clostridiales',
            'Ruminococcaceae','Enterococcus','Citrobacter',
            'Hungatella','Pseudomonas','Richness','Shannon','InvSimp')
    op <- c('Enterococcus',#'Lactococcus',
            'Streptococcus','Citrobacter',
            'Pseudomonas','Klebsiella','Actinobacteria','Hungatella',
            'Rothia','Richness','Shannon','InvSimp')
    taxa.remove <- c('Bacteria','Firmicutes',
                     'Clostridia','Angelakisella')
    
    covars <- c('OSCI_Score_Samp', 'OSCI_Score_Worst',
                'DaysHospitalized','HAP',
                'Bacteremia',
                'AbxCurr','Abx3Mo',
                'Comorbidity','Medication',
                'CRP',
                'Hemoglobin','Procalcitonin','Age','Sex','BMI')
    
    x.facets <- list('Severity' = c('OSCI_Score_Worst','OSCI_Score_Samp',
                                    'Sepsis','Bacteremia',
                                    'DaysHospitalized','HAP'),
                     'MedComor' = c('Comorbidity','AbxCurr','Abx3Mo',
                                    'AbxSamp','AbxHosp','Medication'),
                     'Clinical' = c('CRP','aPTT','Procalcitonin','Hemoglobin'),
                     'Demographic' = c('Age','Sex','BMI'))
    
    .assign.facet <- function(df.plot, list.meta) {
        df.plot %>% 
            mutate(x.facet = map_chr(Covariate, function(cov) {
                if (cov %in% list.meta$Severity) 'Disease Course'
                else if (cov %in% list.meta$MedComor) 
                    'Medication, Comorbidity'
                else if (cov %in% list.meta$Clinical) 'Clinical'
                else if (cov %in% list.meta$Demographic) 'Demographic'
                else 'Other' })) %>%
            mutate(x.facet = fct_relevel(x.facet, 'Disease Course',
                                         'Medication, Comorbidity',
                                         'Clinical',
                                         'Demographic','Other'))
    }
    
    sdna.data <- df.sdna %>%
        filter(feature %in% unique(c(sdna.taxa, gi))) %>%
        # filter(feature %in% unique(c(sdna.taxa, gi,
        #                              sdna.most.ab,
        #                              sdna.most.prev))) %>%
        add_column(y.facet = 'Gut Taxa') %>%
        mutate(`Effect Size` = ifelse(Covariate == 'Hemoglobin',
                                      -1*`Effect Size`, `Effect Size`)) %>%
        left_join(sdna.stats) %>%
        mutate_at(vars(c(prev.all, prev.pats, 
                         mean.ra.all, mean.ra.pats)),
                  ~ case_when(str_detect(
                      feature, 'InvSimp|Shann|Richne') ~ 0,
                      TRUE ~ .))
    
    op.data <- df.op %>%
        filter(feature %in% unique(c(op.taxa, op))) %>%
        # filter(feature %in% unique(c(op.taxa, op,
        #                              op.most.ab,
        #                              op.most.prev))) %>%
        add_column(y.facet = 'Oropharyngeal Taxa')  %>%
        mutate(`Effect Size` = ifelse(Covariate == 'Hemoglobin',
                                      -1*`Effect Size`, `Effect Size`)) %>%
        left_join(op.stats) %>%
        mutate_at(vars(c(prev.all, prev.pats,
                         mean.ra.all, mean.ra.pats)),
                  ~ case_when(str_detect(
                      feature, 'InvSimp|Shann|Richne') ~ 0,
                      TRUE ~ .))
    
    sdna.all.abs <- select(sdna.raw, unique(sdna.data$feature)) %>%
        rownames_to_column('UniqueID') %>%
        as_tibble() %>%
        gather('feature','ab', -UniqueID) %>%
        mutate(ab = ifelse(
            str_detect(feature,'Rich|InvSimp|Shann'),
            0, ab)) %>%
        mutate(y.facet = 'Gut Taxa')
    op.all.abs <- select(op.raw, unique(op.data$feature)) %>%
        rownames_to_column('UniqueID') %>%
        as_tibble() %>%
        gather('feature','ab', -UniqueID) %>%
        mutate(ab = ifelse(
            str_detect(feature,'Rich|InvSimp|Shann'),
            0, ab)) %>%
        mutate(y.facet = 'Oropharyngeal Taxa')
    
    combined.data <- op.data %>%
        bind_rows(sdna.data) %>%
        filter(Covariate %in% covars) %>%
        mutate(Covariate = factor(Covariate, levels = covars)) %>%
        .assign.facet(x.facets) %>%
        filter(!(feature %in% taxa.remove)) 
    
    all.plot.data <- combined.data %>%
        mutate(Stars = as.character(cut(.$FDR,
                                        breaks = c(-Inf, 0.001, 0.01, 0.1, Inf),
                                        label = c("***", "**", "*", "")))) %>%
        mutate(Rob.Stars = ifelse(Confounding %in% c('SD','LD','NC'), Stars, '')) %>%
        mutate(Conf.Stars = ifelse(Confounding %in% c('SD','LD','NC'), '', Stars)) %>%
        nest(data = -feature) %>%
        mutate(rob.star = map_lgl(data, 
                                  ~ ifelse(all(.$Rob.Stars == ''),
                                           FALSE, TRUE))) %>%
        mutate(conf.star = map_lgl(data, 
                                   ~ ifelse(all(.$Conf.Stars == ''),
                                            FALSE, TRUE))) %>%
        filter(rob.star | conf.star) %>%
        select(-c(rob.star, conf.star)) %>%
        unnest(data) %>%
        filter(Covariate %in% covars) %>%
        mutate(y.facet = fct_relevel(y.facet, 'Airway Taxa', 'Gut Taxa'))
    
    feat.order <- #all.plot.data %>%
        df.sdna %>% bind_rows(df.op, .id = 'y.facet') %>% 
        mutate(y.facet = ifelse(y.facet == 1, 'Gut Taxa',
                                'Oropharyngeal Taxa')) %>%
        # filter(!str_detect(feature,'InvSimp|Shannon|Richness')) %>%
        filter(Covariate == 'OSCI_Score_Samp') %>%
        mutate(tmp = paste0(y.facet, feature)) %>%
        filter(tmp %in% c(all.plot.data %>%
                   mutate(tmp = paste0(y.facet, feature)) %>%
                   pull(tmp) %>% unique())) %>%
        arrange(`Effect Size`) %>%
        pull(feature) %>% unique()
    feat.order <- c('Richness','Shannon','InvSimp',
                    setdiff(feat.order, c('Richness','Shannon','InvSimp')))
    
    # 17x13
    main <- all.plot.data %>% 
       mutate(feature = fct_relevel(feature, feat.order)) %>%
        ggplot(aes(x = Covariate, y = feature, fill = `Effect Size`)) +
        facet_grid(cols = vars(x.facet), rows = vars(y.facet),
                   drop = TRUE, scales = 'free', space = 'free') +
        geom_tile() +
        scale_fill_gradient2(limits = c(-1,1)) +
        geom_text(aes(label = Rob.Stars), nudge_y = -0.5, size = 9) + 
        geom_text(aes(label = Conf.Stars), nudge_y = -0.5,
                  color = '#f5f5f5', size = 9) + 
        labs(x = '', y = '') +
        ggembl::theme_embl() +
        ggpubr::labs_pubr() +
        theme(axis.text.x = element_text(angle = -45, hjust = 0, size = 9),
              axis.text.y = element_text(size = 9),
              axis.ticks.y = element_blank()) +
        guides(fill=FALSE)
    
    # browser()
    # take appropriate prevalence and abundances (all v patient-only)
    if (all) {
        all.plot.data <- all.plot.data %>%
            rename(prevalence = 'prev.all') }
    else { 
        all.plot.data <- all.plot.data %>%
            rename(prevalence = 'prev.pats')
        sdna.all.abs <- sdna.all.abs %>%
            filter(!str_detect(UniqueID, 'K')) %>%
            select(-UniqueID) 
        op.all.abs <- op.all.abs %>%
            filter(!str_detect(UniqueID, 'K')) %>%
            select(-UniqueID) }
    
    prev.barplot <- all.plot.data %>%
        distinct(across(c(feature, prevalence, y.facet))) %>%
        mutate(feature = fct_relevel(feature, feat.order)) %>%
        ggplot(aes(x = prevalence, y = feature)) +
        facet_grid(rows = vars(y.facet), scales='free', space = 'free') +
        geom_col() +
        ggpubr::theme_pubr() +
        labs(x='', y='') + 
        theme(axis.text.y = element_blank(),
              axis.ticks.y = element_blank())
    
    # browser()
    rel.ab.barplot <- sdna.all.abs %>%
        bind_rows(op.all.abs) %>%
        mutate(feature = fct_relevel(feature, feat.order)) %>%
        ggplot(aes(x = ab, y = feature)) + 
        geom_vline(xintercept = 0.001, color = scales::muted('red')) +
        ggembl::geom_quantile_box() +
        scale_x_sqrt() +
        facet_grid(rows = vars(y.facet), scales='free', space = 'free') +
        ggpubr::theme_pubr() +
        labs(x='', y='') + 
        theme(axis.text.y = element_blank(), 
              axis.ticks.y = element_blank())
    
    (main | prev.barplot | rel.ab.barplot) +
        plot_layout(widths = c(5,1,1))
                  
}

plot.op.composition <- function(df.feat) {
    
    hap.pats <- df.feat %>%
        mutate(PID = str_split(UniqueID, '-') %>% map_chr(~ head(., 1))) %>%
        filter(PID %in% c('P03','P17','P18','P23','P28')) %>%
        filter(UniqueID != 'P23-OP-V10') %>%
        select(PID, Days, everything(), -UniqueID)
    p <- hap.pats %>%
        summarize_if(is.double, ~ calc.single.prevalence(.)) %>%
        gather('taxa','prev') %>%
        filter(prev > 0) %>%
        arrange(desc(prev)) %>%
        slice(1:4) %>%
        pull(taxa) %>%
        c('Lactococcus','Campylobacter','Lactobacillus','Prevotella')
    
    df.other <- hap.pats %>%
        select(-Days, -PID, contains(p)) %>%
        rowSums() %>%
        bind_cols(select(hap.pats, PID, Days), .) %>%
        rename(Other = '...3')
    df.main <- hap.pats %>%
        select(PID, Days, p) %>%
        left_join(df.other) %>%
        discretize.days(var = 'Days', keep.numeric = TRUE) %>%
        gather('taxa','rel.ab', -c(PID, Days)) %>%
        mutate(taxa = fct_relevel(taxa, 'Other', after = Inf))
    
    df.main %>%
        ggplot(aes(fill=taxa, y=rel.ab, x=Days)) + 
        geom_bar(position="fill", stat="identity") +
        ggthemes::scale_fill_tableau('Miller Stone') +
        facet_grid(~ PID, #ncol = 4, 
                   scales = 'free', space = 'free') +
        labs(title = 'Oropharyngeal samples in HAP patients, dominant taxa',
             x = '', y = 'Relative abundance', fill = '') +
        theme_bw()
    
    # check most abundant taxa per sample ----
    # hap.pats %>% 
    #     filter(str_detect(PID,'P28')) %>% 
    #     summarize_if(is.double,mean) %>% 
    #     gather('genus','pct') %>% 
    #     filter(pct!=0) %>%
    #     arrange(desc(pct)) %>%
    #     mutate(pct = 100*pct)
    
}

plot.tbs.composition <- function(df.feat) {
    
    df.comp <- df.feat %>%
        select(-c(1:2)) %>%
        as.matrix() %>%
        prop.table(1) %>% 
        as_tibble() %>%
        bind_cols(select(df.feat, UniqueID, Days), .) %>%
        filter(!is.na(!!rlang::sym(colnames(.)[3])))
    # for motus-level analysis
    # df.s.comp <- df.spec.feat %>%
    #     select(-c(1:5)) %>%
    #     as.matrix() %>%
    #     prop.table(1) %>% 
    #     as_tibble() %>%
    #     bind_cols(select(tbs.feat, UniqueID, Days), .) %>%
    #     filter(!is.na(!!rlang::sym(colnames(.)[3])))
    df.other <- df.comp %>%
        select(-contains(c('UniqueID','Days','Klebsiella','Achromobacter',
                           'Staphylococcus','Chlorobi','Prevotella',
                           'Citrobacter','Unmapped','Streptococcus','Neisseria')),
               contains('ceae'), Alloprevotella) %>%
        rowSums() %>%
        bind_cols(select(df.comp, UniqueID, Days), .) %>%
        rename(Other = '...3')
    df.main <- df.comp %>%
        select(contains(c('UniqueID','Days','Klebsiella','Achromobacter',
                          'Staphylococcus','Chlorobi','Prevotella',
                          'Citrobacter','Unmapped','Streptococcus','Neisseria')),
               -contains('ceae'),-Alloprevotella) %>%
        left_join(df.other) %>%
        discretize.days(var = 'Days', keep.numeric = TRUE) %>%
        mutate(PID = str_split(UniqueID, '-') %>% map_chr(~ head(., 1))) %>%
        select(PID, Days, everything(), -UniqueID) %>%
        gather('taxa','rel.ab', -c(PID, Days)) %>%
        mutate(taxa = fct_relevel(taxa, 'Other', after = Inf))
    
    colors <- ggthemes::tableau_color_pal(palette = 'Nuriel Stone')(9)
    
    df.main %>%
        ggplot(aes(fill=taxa, y=rel.ab, x=Days)) + 
        geom_bar(position="fill", stat="identity") +
        # ggthemes::scale_fill_tableau('Nuriel Stone') +
        scale_fill_manual(values = c(colors[1:6], colors[8:9], colors[7])) +
        facet_grid(~ PID, #ncol = 4, 
                   scales = 'free', space = 'free') +
        labs(title = 'Trachiobronchial sputum samples in HAP patients, dominant taxa',
             x = '', y = 'Relative abundance', fill = '') +
        theme_bw()
    
    ## which species present in each patients TBS sample? ----
    # df.s.comp %>% 
    #     filter(str_detect(UniqueID,'P17')) %>% 
    #     summarize_if(is.double,mean) %>% 
    #     gather('mOTU','pct') %>% 
    #     filter(pct!=0) %>%
    #     arrange(desc(pct)) %>%
    #     mutate(pct = 100*pct)
    # 
    # df.s.comp %>% 
    #     filter(str_detect(UniqueID,'P17')) %>% 
    #     gather('mOTU','percent.sample',-c(UniqueID, Days)) %>% 
    #     mutate(percent.sample = round(percent.sample*100, digits=2)) %>%
    #     filter(percent.sample!=0) %>%
    #     discretize.days(var='Days',keep.numeric=TRUE) %>% 
    #     arrange(Days) %>%
    #     rename(Days = 'Days') %>%
    #     mutate_all(~ as.character(.))
    
}

# moved to scripts/curated-heatmap.R for now
plot.motus.heatmap <- function(df.sdna, df.op, ...) {
    
    args <- list(...)
    
    df.sdna <- df.sdna %>%
        mutate(Covariate = str_replace_all(Covariate,'Antibiotics','AbxCurr')) %>%
        mutate(Covariate = str_replace_all(Covariate,'Pneumonia','HAP'))
    df.op <- df.op %>%
        mutate(Covariate = str_replace_all(Covariate,'Antibiotics','AbxCurr')) %>%
        mutate(Covariate = str_replace_all(Covariate,'Pneumonia','HAP'))
    
    sdna.taxa <- df.sdna %>%
        filter(Confounding %in% c('LD','SD','NC') & FDR <= 0.05) %>%
        filter(str_detect(Covariate,'OSCI|Comorb||AbxCurr')) %>%
        pull(feature) %>% unique()
    op.taxa <- df.op %>%
        filter(Confounding %in% c('LD','SD','NC') & FDR <= 0.05) %>%
        filter(str_detect(Covariate,'OSCI|HAP|AbxCurr')) %>%
        pull(feature) %>% unique()
    
    covars <- c('OSCI_Score_Samp', 'OSCI_Score_Worst',
                'DaysHospitalized','HAP',
                'Bacteremia',
                'AbxCurr','Abx3Mo',
                'Comorbidity','Medication',
                'CRP',
                'Hemoglobin','Procalcitonin','Age','Sex','BMI')
    
    x.facets <- list('Severity' = c('OSCI_Score_Worst','OSCI_Score_Samp',
                                    'Sepsis','Bacteremia',
                                    'DaysHospitalized','HAP'),
                     'MedComor' = c('Comorbidity','AbxCurr','Abx3Mo',
                                    'AbxSamp','AbxHosp','Medication'),
                     'Clinical' = c('CRP','aPTT','Procalcitonin','Hemoglobin'),
                     'Demographic' = c('Age','Sex','BMI'))
    
    .assign.facet <- function(df.plot, list.meta) {
        df.plot %>% 
            mutate(x.facet = map_chr(Covariate, function(cov) {
                if (cov %in% list.meta$Severity) 'Disease Course'
                else if (cov %in% list.meta$MedComor) 
                    'Medication, Comorbidity'
                else if (cov %in% list.meta$Clinical) 'Clinical'
                else if (cov %in% list.meta$Demographic) 'Demographic'
                else 'Other' })) %>%
            mutate(x.facet = fct_relevel(x.facet, 'Disease Course',
                                         'Medication, Comorbidity',
                                         'Clinical',
                                         'Demographic','Other'))
    }
    
    sdna.data <- df.sdna %>%
        filter(feature %in% sdna.taxa) %>%
        add_column(y.facet = 'Gut Taxa') %>%
        mutate(`Effect Size` = ifelse(Covariate == 'Hemoglobin',
                                      -1*`Effect Size`, `Effect Size`))
    
    op.data <- df.op %>%
        filter(feature %in% op.taxa) %>%
        add_column(y.facet = 'Airway Taxa')  %>%
        mutate(`Effect Size` = ifelse(Covariate == 'Hemoglobin',
                                      -1*`Effect Size`, `Effect Size`))
    
    combined.data <- op.data %>%
        bind_rows(sdna.data) %>%
        filter(Covariate %in% covars) %>%
        mutate(Covariate = factor(Covariate, levels = covars)) %>%
        .assign.facet(x.facets) 
    
    all.plot.data <- combined.data %>%
        mutate(Stars = as.character(cut(.$FDR,
                                        breaks = c(-Inf, 0.001, 0.01, 0.1, Inf),
                                        label = c("***", "**", "*", "")))) %>%
        mutate(Rob.Stars = ifelse(Confounding %in% c('SD','LD','NC'), Stars, '')) %>%
        mutate(Conf.Stars = ifelse(Confounding %in% c('SD','LD','NC'), '', Stars)) %>%
        nest(data = -feature) %>%
        mutate(rob.star = map_lgl(data, 
                                  ~ ifelse(all(.$Rob.Stars == ''),
                                           FALSE, TRUE))) %>%
        mutate(conf.star = map_lgl(data, 
                                   ~ ifelse(all(.$Conf.Stars == ''),
                                            FALSE, TRUE))) %>%
        filter(rob.star | conf.star) %>%
        select(-c(rob.star, conf.star)) %>%
        unnest(data) %>%
        filter(Covariate %in% covars) %>%
        mutate(y.facet = fct_relevel(y.facet, 'Airway Taxa', 'Gut Taxa'))
    
    feat.order <- all.plot.data %>%
        filter(Covariate == 'OSCI_Score_Samp') %>%
        arrange(`Effect Size`) %>%
        # mutate(feature = as_factor(feature)) %>%
        pull(feature) %>% unique()
    feat.order <- c('Richness','Shannon','InvSimp',
                    setdiff(feat.order, c('Richness','Shannon','InvSimp')))
    
    # 17x13
    all.plot.data %>% 
        #mutate(feature = fct_relevel(feature, feat.order)) %>%
        ggplot(aes(x = Covariate, y = feature, fill = `Effect Size`)) +
        facet_grid(cols = vars(x.facet), rows = vars(y.facet),
                   drop = TRUE, scales = 'free', space = 'free') +
        geom_tile() +
        scale_fill_gradient2(limits = c(-1,1)) +
        geom_text(aes(label = Rob.Stars), nudge_y = -0.5, size = 9) + 
        geom_text(aes(label = Conf.Stars), nudge_y = -0.5,
                  color = '#f5f5f5', size = 9) + 
        labs(x = '', y = '') +
        ggembl::theme_embl() +
        ggpubr::labs_pubr() +
        theme(axis.text.x = element_text(angle = -45, hjust = 0, size = 9),
              axis.text.y = element_text(size = 9),
              axis.ticks.y = element_blank())
    
}

################################################
######### METABOLOME FIGURES, HEATMAPS #########
################################################
plot.metabolite.volcanos <- function(df.res, df.annot,
                                     cov = 'OSCI_Score_Samp',
                                     title, FDR.lim = 5) {
    
    df.plot <- df.res %>% 
        filter(!str_detect(feature,'_13C|_D5|_D3')) %>%
        filter(Covariate==cov) %>% 
        mutate(R_Compound=str_remove_all(feature,'P_')) %>%
        mutate(R_Compound=str_remove_all(R_Compound,'U_')) %>% 
        mutate(Site=str_split(feature,'_') %>% 
                   map_chr(~ head(.,1))) %>% 
        mutate(Site=case_when(
            Site=='P'~'Plasma',
            Site=='U'~'Urine')) %>%
        mutate(Confounding=case_when(
            Confounding%in%c('LD','SD','NC') & cov!='OSCI_Score_Samp' ~ 'Robust',
            cov=='OSCI_Score_Samp'&str_detect(
                Confounding,'Worst|LD|SD|NC') ~ 'Robust',
            Confounding=='NS' ~'Not Significant',
            TRUE ~ 'Confounded') ) %>%
        mutate(Confounding=factor(Confounding,
                                  levels=c('Not Significant',
                                           'Confounded',
                                           'Robust'))) %>%
        left_join(select(df.annot, 
                         R_Compound,ChemClass,Compound,Pathway)) %>%
        mutate(FDR = -log10(FDR))
    
    df.plot %>% 
        ggplot(aes(x=`Effect Size`, y=FDR)) +
        geom_point(aes(color=ChemClass, shape = Site,
                       alpha = Confounding), size=3) +
        geom_hline(yintercept=-log10(0.05), color='red') + 
        labs(title = paste0(cov,': ', title),
             y='-log10(FDR)') +
        scale_alpha_manual(values = c(0.2,0.6,1)) +
        ggrepel::geom_label_repel(
            dplyr::filter(df.plot, FDR>FDR.lim|str_detect(ChemClass,'Bile')), 
            size = 7/.pt, mapping = aes(label = Compound), 
            min.segment.length = unit(0, 'lines'), force = 5) +
        ggembl::theme_publication() 
}

cleanup.metabo.results <- function(df.res, cov, site) {
    
    df.res %>%
        filter(!str_detect(feature,'_13C|_D5|_D3')) %>%
        filter(Covariate == cov) %>%
        mutate(R_Compound=str_remove_all(feature, 'P_')) %>%
        mutate(R_Compound=str_remove_all(R_Compound, 'U_')) %>%
        mutate(Confounding=case_when(
            Confounding%in%c('LD','SD','NC') ~ 'R',
            cov=='OSCI_Score_Samp'&Confounding=='OSCI_Score_Worst'&FDR<0.1 ~ 'R',
            cov=='OSCI_Score_Samp'&Confounding=='COVID'&FDR<0.1 ~ 'R',
            Confounding=='NS' ~'NS',
            TRUE ~ 'C') ) %>%
        mutate(Confounding=factor(Confounding,
                                  levels=c('NS',
                                           'C',
                                           'R'))) %>%
        select(feature, R_Compound, `Effect Size`, FDR, Confounding)

}

plot.metab.comparative.es <- function(early.res, late.res, df.annot,
                                      cov = 'OSCI_Score_Samp',
                                      site = 'P_') {
    
    df.early <- early.res %>% 
        cleanup.metabo.results(cov, site) %>%
        mutate(FDR = -log10(FDR)) %>%
        magrittr::set_colnames(c('feature','R_Compound','early.es',
                                 'early.fdr','early.conf'))
    df.late <- late.res %>% 
        cleanup.metabo.results(cov, site) %>%
        mutate(FDR = -log10(FDR)) %>%
        magrittr::set_colnames(c('feature','R_Compound','late.es',
                                 'late.fdr','late.conf'))
    df.plot <- df.early %>%
        full_join(df.late, by=c('feature','R_Compound')) %>%
        left_join(select(df.annot, 
                         R_Compound,ChemClass,Compound),
                  by='R_Compound') %>%
        select(-R_Compound) %>%
        mutate(status = map2_chr(
            early.conf, late.conf,
            ~ case_when(.x=='NS' & .y=='NS' ~ 'Background',
                        .x=='R' & .y=='R' ~ 'Robust',
                        .x=='C' & .y=='C' ~ 'Confounded',
                        .x=='NS' & .y=='C' ~ 'Background',
                        .x=='C' & .y=='NS' ~ 'Background',
                        .x=='R' & .y=='NS' ~ 'Early only',
                        .x=='NS' & .y=='R' ~ 'Late only',
                        .x=='R' & .y=='C' ~ 'Early only',
                        .x=='C' & .y=='R' ~ 'Late only'))) %>%
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
                         'Nucleobase-related','Nucleoside','Nucleobase')))  %>%
        mutate(status = factor(status, levels = c('Robust',
                                                  'Early only',
                                                  'Late only',
                                                  'Confounded',
                                                  'Background'))) %>%
        mutate(ChemClass = fct_relevel(
            ChemClass, c('Amino acids & related',
                         'Aromatic AAs & derivatives',
                         'Indole & derivatives',
                         'Biogenic amines',
                         'Nucleoside & related',
                         'Sugars','Glycolysis & TCA',
                         'Choline & derivatives',
                         'Alkaloids','Acylcarnitines',
                         'Hormones','Bile acids',
                         'Carboxylic acids',
                         'Glycerol','Diacylglycerols','Fatty acids',
                         'SCFAs',
                         'Phosphatidylcholines',
                         'Triglycerides','Ceramides',
                         'Sphingolipids',
                         'Lysophosphatidylcholines') )) 
    
    # pick labels
    if (site == 'P_') {
        labels <- df.plot %>%
            filter(str_detect(feature, site) & (
                !str_detect(status,'Conf|Backg') | str_detect(Compound, 'Taurodeoxychol')))
    } 
    else if (site == 'U_') {
        labels <- df.plot %>%
            filter(str_detect(feature, site)) %>%
            filter(!str_detect(status,'Background'))
    } 
    else { stop('Sample site not recognized, try again') }
    
    colors <- scico::scico(23, palette = 'roma', direction = -1)
    
    df.plot.filt <- df.plot %>%
        filter(!is.na(ChemClass) & str_detect(feature, site)) %>%
        mutate(ChemClass = fct_drop(ChemClass))
    df.plot.filt %>%
        ggplot(aes(x = early.es, y = late.es)) +
        geom_hline(yintercept = 0, color = 'grey') +
        geom_vline(xintercept = 0, color = 'grey') +
        geom_abline(intercept = 0, slope = 1, color = 'black') +
        geom_point(data = filter(df.plot.filt, str_detect(status,'Background')),
                   aes(shape = status, color = ChemClass, fill = ChemClass), 
                   size = 4, alpha = 0.5, color = 'grey') +
        geom_point(data = filter(df.plot.filt, !str_detect(
            status,'Background')), 
            aes(color = ChemClass, fill = ChemClass, shape = status), 
            alpha = 1, size = 4) +
        lims(x=c(-1,1), y=c(-1,1)) +
        #scale_shape_manual(values = c(21,13,22:24)) +
        scale_shape_manual(values = c(16,13,15,18,17)) +
        scale_color_manual(values = colors, drop = FALSE) +
        scale_fill_manual(values = colors, drop = FALSE) +
        labs( 
            x = "Spearman's ρ OSCI (Early, ≤10 days symptomatic)",
            y = "Spearman's ρ OSCI (Late, >10 days symptomatic)",
            shape = 'Significance',
            color = '',
            fill = '') +
        ggrepel::geom_text_repel(data = labels,
                                 mapping = aes(label = Compound),
                                 size = 8/.pt,
                                 min.segment.length = unit(0, 'lines'),
                                 force = 7,
                                 max.overlaps = 20) +
        ggembl::theme_embl() +
        guides(fill = 'none')

}

plot.metabolite.heatmap <- function(df.feat, ...) {

    args <- list(...)
    
    covars <- c('COVID','OSCI_Score_Samp','OSCI_Score_Worst',
                'DaysHospitalized',
                'HAP','Bacteremia',
                'AbxCurr','Abx3Mo','Antivirals',
                'Comorbidity','Medication',
                'Age','Sex','BMI')
    
    x.facets <- list('Severity' = c('COVID','OSCI_Score_Worst','OSCI_Score_Samp',
                                    'HAP','Bacteremia',
                                    'DaysHospitalized'),
                     'MedComor' = c('Comorbidity','AbxCurr','Abx3Mo',
                                    'AbxHosp','Medication','Antivirals'),
                     'Demographic' = c('Age','Sex','BMI'))
    
    .assign.facet <- function(df.plot, list.meta) {
        df.plot %>% 
            mutate(x.facet = map_chr(Covariate, function(cov) {
                if (cov %in% list.meta$Severity) 'Disease Course'
                else if (cov %in% list.meta$MedComor) 
                    'Medication, Comorbidity'
                else if (cov %in% list.meta$Clinical) 'Clinical'
                else if (cov %in% list.meta$Demographic) 'Demographic'
                else 'Other' })) %>%
            mutate(x.facet = fct_relevel(x.facet, 'Disease Course',
                                         'Medication, Comorbidity',
                                         'Demographic','Other'))
    }
    
    dat <- df.feat %>%
        mutate(R_Compound = feature) %>%
        mutate(R_Compound = str_remove_all(R_Compound,'P_')) %>%
        mutate(R_Compound = str_remove_all(R_Compound,'U_')) %>%
        left_join(args[[1]]) %>%
        filter(Covariate != 'PatientID') %>%
        mutate(Covariate = str_replace_all(Covariate, 'Pneumonia','HAP')) %>%
        mutate(Covariate = str_replace_all(Covariate, 'Antibiotics','AbxCurr')) %>%
        filter(Covariate %in% covars) %>%
        mutate(Covariate = fct_relevel(Covariate, covars)) %>%
        .assign.facet(x.facets)
    
    of.interest <- dat %>%
        filter((str_detect(Covariate,
                           'OSCI|HAP|Comorbidity|Medication') & Confounding %in% c('LD','SD','NC')) | str_detect(
                               Pathway,'rypt|ynur|indole') | ChemClass=='Bile acids')
    
    plot.dat <- dat %>%
        filter(feature %in% of.interest$feature) %>%
        mutate(Stars = as.character(cut(.$FDR,
                                        breaks = c(-Inf, 0.001, 0.01, 0.1, Inf),
                                        label = c("***", "**", "*", "")))) %>%
        mutate(Rob.Stars = ifelse(Confounding %in% c('SD','LD','NC'), Stars, '')) %>%
        mutate(Conf.Stars = ifelse(Confounding %in% c('SD','LD','NC'), '', Stars)) %>%
        nest(data = -feature) %>%
        mutate(rob.star = map_lgl(data, 
                                  ~ ifelse(all(.$Rob.Stars == ''),
                                           FALSE, TRUE))) %>%
        mutate(conf.star = map_lgl(data, 
                                   ~ ifelse(all(.$Conf.Stars == ''),
                                            FALSE, TRUE))) %>%
        filter(rob.star | conf.star) %>%
        select(-c(rob.star, conf.star)) %>%
        unnest(data) %>%
        filter(Covariate %in% covars) %>%
        mutate(Site = str_split(feature,'_') %>% map_chr(~ head(.,1))) %>%
        mutate(Site = str_replace_all(Site,'P','Plasma')) %>%
        mutate(Site = str_replace_all(Site,'U','Urine')) %>%
        mutate(Method = ifelse(Compound=='indole_3_propionic_acid','TRYP',Method)) %>%
        filter(!is.na(Method))
    
    feat.order <- plot.dat %>%
        filter(Covariate == 'OSCI_Score') %>%
        arrange(`Effect Size`) %>%
        # mutate(feature = as_factor(feature)) %>%
        pull(Compound) %>% unique()
    
    plot.dat %>% 
        mutate(Compound = fct_relevel(Compound, feat.order)) %>%
        ggplot(aes(x = Covariate, y = Compound, fill = `Effect Size`)) +
        facet_grid(cols = vars(x.facet), rows = vars(Site, Method),
                   drop = TRUE, scales = 'free', space = 'free') +
        geom_tile() +
        scale_fill_gradient2(limits = c(-1,1)) +
        geom_text(aes(label = Rob.Stars), nudge_y = -0.5, size = 9) + 
        geom_text(aes(label = Conf.Stars), nudge_y = -0.5,
                  color = '#f5f5f5', size = 9) + 
        labs(x = '', y = '') +
        ggembl::theme_embl() +
        ggpubr::labs_pubr() +
        theme(axis.text.x = element_text(angle = -45, hjust = 0, size = 9),
              axis.text.y = element_text(size = 9),
              axis.ticks.y = element_blank())
}

plot.chemclass.heatmap <- function(df.feat) {
    
    plot.dat <- df.feat %>%
        mutate(Stars = as.character(cut(.$FDR,
                                        breaks = c(-Inf, 0.001, 0.01, 0.1, Inf),
                                        label = c("***", "**", "*", "")))) %>%
        mutate(Rob.Stars = ifelse(
            Confounding %in% c('SD','LD','NC'), Stars, '')) %>%
        mutate(feature = str_replace_all(feature,'_',' ')) %>%
        arrange(`Effect Size`) %>%
        mutate(feature = as_factor(feature)) #%>%
        # mutate(Group = ifelse(Group=='Mild',
        #                         'Mild vs. Control',
        #                         'Severe vs. Control')) 
    
    plot.dat %>% 
        ggplot(aes(x=Group, y=feature, fill=`Effect Size`)) + 
        geom_tile() +
        # scale_fill_gradient2(limits = c(-1,1)) +
        scale_fill_gradient2(
            low = scales::muted('darkorange'),
            mid = 'white',
            high = scales::muted('green'),
            limits = c(-1,1)) +
        geom_text(aes(label=Rob.Stars), 
                  nudge_y = -0.3, size = 9, color = 'white') +
        geom_text(data = filter(plot.dat, !Confounding%in%c('LD','SD','NS')),
                  aes(label=Confounding), color = 'white') +
        ggembl::theme_embl() +
        ggpubr::labs_pubr() +
        labs(x='',y='') +
        theme(axis.text.y = element_text(size = 9),
              axis.ticks.y = element_blank(),
              legend.position = 'top')
}

plot.osci.biplot <- function(pca.obj, df.meta, frame = TRUE) {

    if (frame) {
        ggplot2::autoplot(pca.obj, data = df.meta %>%
                              select(ShortID, contains('OSCI')) %>%
                              mutate(OSCI_Class = factor(OSCI_Class,
                                                         levels=c('Control',
                                                                  'Mild',
                                                                  'Severe'))) %>%
                              filter(ShortID %in% rownames(pca.obj$x)) %>%
                              column_to_rownames('ShortID'), 
                          colour = 'OSCI_Class', label = TRUE, shape = FALSE,
                          loadings = TRUE, loadings.colour = 'black',
                          loadings.label = TRUE, loadings.label.size = 4.5,
                          loadings.label.colour = 'black', frame = TRUE) +
            ggthemes::scale_color_tableau('Superfishel Stone') +
            ggthemes::scale_fill_tableau('Superfishel Stone') +
            ggpubr::theme_pubr() 
    } else {
        ggplot2::autoplot(pca.obj, data = df.meta %>%
                              select(ShortID, contains('OSCI')) %>%
                              mutate(OSCI_Class = factor(OSCI_Class,
                                                         levels=c('Control',
                                                                  'Mild',
                                                                  'Severe'))) %>%
                              filter(ShortID %in% rownames(pca.obj$x)) %>%
                              column_to_rownames('ShortID'), 
                          colour = 'OSCI_Class', label = TRUE, shape = FALSE,
                          loadings = TRUE, loadings.colour = 'black',
                          loadings.label = TRUE, loadings.label.size = 4.5,
                          loadings.label.colour = 'black') +
            ggthemes::scale_color_tableau('Superfishel Stone') +
            ggthemes::scale_fill_tableau('Superfishel Stone') +
            ggpubr::theme_pubr() 
    }
}

prep.radar.stats <- function(df.res, dataset.nm, pct = TRUE) {
    
    tmp <- df.res %>%
        mutate(Confounding = str_replace_all(Confounding,
                                             'Antibiotics|AbxHosp',
                                             'AbxCurr')) %>%
        mutate(Confounding = str_replace_all(Confounding,
                                             'Pneumonia',
                                             'HAP'))
    
    osci <- tmp %>%
        filter(str_detect(Covariate, 'Score_Samp')) %>%
        mutate(Covariate = str_split(Covariate,'_') %>%
                   map_chr(~ head(.,1))) %>%
        mutate(Covariate = paste0(Covariate,'_Sample'))
    
    plot.data <- osci %>%
        mutate(Confounding = str_split(Confounding, ',')) %>%
        unnest_longer(Confounding) %>%
        mutate_at(vars(Confounding), ~ trimws(.)) %>%
        distinct(across(c(feature,Confounding)), .keep_all = TRUE) %>%
        mutate(Confounding = str_replace_all(Confounding, 
                                             'SD|LD|NC','Deconfounded')) %>%
        mutate(Confounding = str_replace_all(Confounding, 
                                             'OSCI_Score_Worst','OSCI_Worst')) %>%
        mutate(Stars = as.character(
            cut(.$FDR,
                breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
                label = c("***", "**", "*", "")))) %>%
        count(Covariate, Confounding) %>%
        spread('Covariate','n') %>%
        mutate_if(is.numeric, ~ replace_na(., 0)) %>%
        # new additions
        mutate(assoc.status = case_when(
            Confounding == 'Deconfounded' ~ 'Deconfounded',
            Confounding == 'NS' ~ 'Not Significant',
            TRUE ~ 'Confounded')) %>%
        rename(assoc.type = 'Confounding', n = 'OSCI_Sample') %>%
        select(assoc.type, assoc.status, n) %>%
        arrange(desc(n))
    
    n.nonconf <- plot.data %>%
        filter(assoc.status != 'Confounded') %>%
        pull(n) %>% sum()
    pct.nonconf <- n.nonconf/nrow(osci)*100
    n.conf <- nrow(osci)-n.nonconf
    pct.conf <- n.conf/nrow(osci)*100
    # confounded and robust -- all shown in radar!
    n.sig <- plot.data %>%
        filter(str_detect(assoc.status, 'onf')) %>%
        pull(n) %>% sum()
    
    # plot.dat.out <- plot.data %>%
    #     arrange(desc(n)) %>%
    #     select(assoc.type, n) %>%
    #     filter(assoc.type != 'NS')
    
    conf.stats.abs <- plot.data %>%
        # filter(assoc.type != 'NS') %>%
        group_by(assoc.status) %>%
        mutate(status.grp.n = sum(n)) %>%
        ungroup() %>%
        # with_groups(assoc.status, 
        #             ~ mutate(.x, status.grp.n = sum(n))) %>%
        mutate(status.prop = n/status.grp.n) %>%
        # best estimate ?
        mutate(n.hat = ifelse(
            assoc.status=='Confounded',
            n.conf*status.prop,
            n)) %>%
        add_row(assoc.type = 'total',
                assoc.status = 'total',
                n = nrow(osci), 
                status.grp.n = NA, status.prop = NA, 
                n.hat = nrow(osci)) %>%
        filter(n.hat >= 1)
        
    if (pct) {
        conf.stats.out <- conf.stats.abs %>%
            mutate(n.hat = n.hat/nrow(osci)*100) }
    else { conf.stats.out <- conf.stats.abs }
    
    conf.stats.out %>%
        dplyr::rename(!!rlang::sym(dataset.nm) := 'n.hat') %>%
        select(assoc.type, !!rlang::sym(dataset.nm))
}

plot.radar <- function(df.plot, plot.title,
                       incl.ns = FALSE, slice.type = 'time',
                       min = 0, mid = 15, max = 75) {
    
    if (incl.ns) {
        tmp <- df.plot %>%
            select(Timepoint, NS, Deconfounded, 
                   COVID, OSCI_Worst, Abx3Mo, AbxCurr, HAP,
                   # AbxHosp,
                   # Age, Supplements, 
                   # Bacteremia, Renal,
                   Medication, Comorbidity,
                   DaysHospitalized) }
    else {
    tmp <- df.plot %>%
        select(Timepoint, Deconfounded, COVID, OSCI_Worst, 
               Abx3Mo, AbxCurr, HAP,
               # AbxHosp, NS,
               #Age, Supplements, 
               # Bacteremia, Renal,
               Medication, Comorbidity,
               DaysHospitalized) }
    
    if (slice.type == 'time') {
        dat <- tmp %>%
            filter(Timepoint != 'All') %>%
            mutate(Timepoint = factor(Timepoint, 
                                      levels = c('Early','Late'))) 
        legend.title <- 'Timepoint'
        group.colors = c("#00AFBB", "#E7B800")
    } else if (slice.type == 'patient.subset') {
        dat <- tmp %>%
            mutate(Timepoint = factor(Timepoint, 
                                      levels = c('All','Patients'))) 
        legend.title <- 'Group'
        group.colors = c("#00AFBB", "#E7B800")
    } else if (slice.type == 'spaces.combo') {
        dat <- tmp %>%
            mutate(Timepoint = factor(Timepoint, 
                                      levels = c('Oropharyngeal',
                                                 'Plasma',
                                                 'Stool',
                                                 'Urine'))) 
        legend.title <- 'Space'
        group.colors <- c('#787405', '#737373','#8B9E95','#372234')
    } else { stop('need to determine datasets being compared')}
    
    tmp %>%
        rename(OSCI = 'Deconfounded') %>%
        ggradar::ggradar(
            values.radar = c(as.character(min),
                             as.character(mid),
                             as.character(max)),
            grid.min = min, 
            grid.mid = mid, 
            grid.max = max,
            # Polygons
            group.line.width = 1, 
            group.point.size = 3,
            # group.colours = group.colors,
            # Background and grid lines
            background.circle.colour = "white",
            gridline.mid.colour = "grey",
            gridline.min.linetype = 'solid',
            gridline.max.linetype = NA,
            legend.position = "right",
            plot.title = plot.title,
            legend.title = legend.title,
            #fill = TRUE,
            fill.alpha = 0.3
        )
}

plot.metab.barplot <- function(df.plot) {
    
    df.plot %>%
        mutate(Confounding = str_replace_all(Confounding,
                                             'Deconfounded','Specific/Robust')) %>%
        mutate(Confounding = fct_rev(factor(Confounding,
                                    levels = c('Not Significant','Specific/Robust',
                                               'OSCI_Worst','COVID','Hospitalization_HAP',
                                               'AbxCurr','Abx3Mo','Kidney_BP',
                                               'Med_Supps','Antivirals','Comorbidity',
                                               'Age_BMI')))) %>%
        mutate(Timepoint = fct_rev(factor(Timepoint, levels = c('Early','Late')))) %>%
        ggplot(aes(x=Confounding, y=Pct, group=Timepoint)) + 
        geom_bar(aes(fill=Timepoint), stat="identity", 
                 position="dodge") + 
        # geom_errorbar(aes(ymin= Pct - sds, ymax = Pct + sds, width=0.2), 
        #               position=position_dodge(width=0.90)) +
        coord_flip() +
        scico::scale_fill_scico_d(begin = 0.3, end = 0.8) +
        labs(x='', y='', fill='') +
        ggembl::theme_publication() +
        theme(panel.grid.major.x = element_line(color = "grey"),
              panel.grid.minor.x = element_line(color = 'grey', linetype = 'dotted'))
}

plot.metab.boxplots <- function(df.plot, met.name, comparisons,
                                together = FALSE, side = 'left') {
    tmp <- df.plot %>%
        ggplot(aes(x=Group, y=Val)) +
        geom_boxplot(aes(fill=OSCI_Class_Samp), 
                     alpha = 0.2, outlier.shape = NA) +
        geom_jitter(aes(fill=OSCI_Class_Samp), 
                    shape=21) +
        scico::scale_fill_scico_d(palette = 'lapaz', 
                                  direction=-1) +
        ggpubr::stat_compare_means(comparisons=comparisons, aes(label=..p.signif..),
                                   method='wilcox.test', hide.ns=TRUE,
                                   tip.length=0, vjust=0.5) +
        scale_y_log10(position = side) +
        labs(x='',y='', title=met.name) + 
        guides(fill='none') +
        ggpubr::theme_pubr() +
        theme(#axis.text.x = element_text(angle = -45, hjust = 0, size = 9),
              axis.text.x = element_blank(),
              axis.ticks.x = element_blank(),
              axis.text.y = element_text(size = 9),
              axis.ticks.y = element_blank())
    
    # plasma and urine 
    if (together) { 
        tmp + scale_x_discrete(drop=FALSE) }
    else tmp

}

################################################
############# INTEGRATION FIGURES ##############
################################################
plot.confounder.status.comp <- function(df.stat) {
    
    df.stat %>%
        mutate(space = factor(space, levels=c(
            'Plasma','Host Immune System','Stool',
            'Urine','Oropharyngeal'))) %>%
        mutate(assoc.status = factor(assoc.status, levels=c(
            'Not Significant','Confounded','Deconfounded'))) %>%
        ggplot(aes(fill=assoc.status, y=prop, x=space)) + 
        geom_bar(position="fill", stat="identity") +
        ggthemes::scale_fill_tableau('Miller Stone') +
        labs(title = '',
             x = '', y = 'Proportion of associations', fill = '') +
        theme_minimal() +
        theme(axis.text.x = element_text(
            angle = -45, hjust = 0, size = 9)) +
        coord_flip()
    
}

prep.network <- function(df.mb_met, df.met_mb,
                            df.ck_mb.met, df.pbmc_mb.met,
                            ...) {
    args <- list(...)[[1]]
    # browser()
    
    sdna.metab.feats <- df.mb_met %>%
        filter(str_detect(Covariate, 'P_')) %>%
        filter(Confounding %in% c('LD','SD')) %>%
        left_join(args$sdna.genus.prev) %>%
        left_join(args$sdna.genus.mean.ra) %>%
        filter(prevalence > 0.25 & mean.ra > 0.001) %>%
        mutate(X_space = str_split(Covariate, '_') %>%
                   map_chr(~ head(., 1))) %>%
        mutate(Covariate = str_remove_all(Covariate, 'P_')) %>%
        left_join(select(args$metabolite.categ, R_Compound, Compound),
                  by = c('Covariate' = 'R_Compound')) %>%
        select(-Covariate) %>% rename(Covariate = 'Compound') %>%
        select(feature, Covariate, everything(), -prevalence, -mean.ra)
    metab.micro.feats <- df.met_mb %>%
        filter(str_detect(feature, 'P_')) %>%
        filter(str_detect(Covariate, 'S_|SDNA_')) %>%
        filter(Confounding %in% c('LD','SD')) %>%
        mutate(feature = str_remove_all(feature, 'P_')) %>%
        left_join(select(args$metabolite.categ, R_Compound, Compound),
                  by = c('feature' = 'R_Compound')) %>%
        select(-feature) %>% rename(feature = 'Compound') %>%
        mutate(X_space = str_split(Covariate, '_') %>%
                   map_chr(~ head(., 1))) %>%
        mutate(Covariate = str_remove_all(Covariate, 'S_|SDNA_')) %>%
        left_join(args$sdna.genus.prev, by = c('Covariate' = 'feature')) %>%
        left_join(args$sdna.genus.mean.ra, by = c('Covariate' = 'feature')) %>%
        mutate_at(vars(prevalence, mean.ra), ~ ifelse(is.na(.), 1, .)) %>%
        # were filtered before regression; no change
        filter(prevalence > 0.25 & mean.ra > 0.001) %>%
        select(feature, everything(), -prevalence, -mean.ra)
    ck.metab.micro.feats <- df.ck_mb.met %>%
        filter(str_detect(Covariate, 'P_|S_')) %>%
        filter(Confounding %in% c('LD','SD') |
                   str_detect(feature,'IL6|TNFa|IP10') &
                   Confounding=='P_kynurenine, P_Kynurenine') %>%
        mutate(X_space = str_split(Covariate, '_') %>%
                   map_chr(~ head(., 1))) %>%
        mutate(Covariate = str_remove_all(Covariate, 'P_|S_')) %>%
        left_join(select(args$metabolite.categ, R_Compound, Compound),
                  by = c('Covariate' = 'R_Compound')) %>%
        mutate(Covariate = ifelse(is.na(Compound), Covariate, Compound)) %>%
        select(-Compound) 
    pbmc.metab.micro.feats <- df.pbmc_mb.met %>%
        filter(str_detect(Covariate, 'P_|S_')) %>%
        filter(Confounding %in% c('LD','SD')) %>%
        mutate(X_space = str_split(Covariate, '_') %>%
                   map_chr(~ head(., 1))) %>%
        mutate(Covariate = str_remove_all(Covariate, 'P_|S_')) %>%
        left_join(select(args$metabolite.categ, R_Compound, Compound),
                  by = c('Covariate' = 'R_Compound')) %>%
        mutate(Covariate = ifelse(is.na(Compound), Covariate, Compound)) %>%
        select(-Compound) %>%
        filter(feature == 'ISG')
    
    return(list('sdna.metab.feats' = sdna.metab.feats,
                'metab.micro.feats' = metab.micro.feats,
                'ck.metab.micro.feats' = ck.metab.micro.feats,
                'pbmc.metab.micro.feats' = pbmc.metab.micro.feats))
    
}

get.network.dat <- function(res.list) {
    
    combined <- res.list$sdna.metab.feats %>%
        bind_rows(list(# op.metab.feats, 
            res.list$metab.micro.feats, 
            res.list$ck.metab.micro.feats,
            res.list$pbmc.metab.micro.feats), 
            .id = 'Y_space') %>%
        mutate(Y_space = case_when(
            Y_space == '1' ~ 'Gut microbiome',
            Y_space == '2' ~ 'Plasma metabolome',
            TRUE ~ 'Host Immune Response')) %>%
        rename(Y_dep.var = 'feature', X_ind.var = 'Covariate',
               XY_effsize = 'Effect Size', XY_FDR = 'FDR') %>%
        mutate(X_space = case_when(
            str_detect(X_space, 'S|SDNA') ~ 'Gut microbiome',
            X_space == 'P' ~ 'Plasma metabolome')) %>%
        # removes 13C etc, ~6 features
        filter(!is.na(Y_dep.var) & !is.na(X_ind.var)) %>%
        select(Y_space, Y_dep.var, X_ind.var, X_space, contains('XY_'))
    
    plot.dat <- combined %>%
        # remove 5 features with feedback/duplicate associations
        mutate(tmp1 = paste0(Y_dep.var, X_ind.var)) %>%
        mutate(tmp2 = paste0(X_ind.var, Y_dep.var)) %>%
        # filter(!tmp1 %in% tmp2) %>%
        select(-tmp1, -tmp2) %>%
        select(X_ind.var, Y_dep.var, everything()) %>%
        rename(from = 'X_ind.var', to = 'Y_dep.var') %>%
        filter(XY_FDR <= 0.05) %>%
        mutate(XY_corr.dir = ifelse(XY_effsize < 0, 'neg','pos')) %>%
        distinct(across(c(from, to)), .keep_all = TRUE) %>%
        mutate(space = case_when(
            (str_detect(X_space,'Gut') & str_detect(Y_space,'Plasma')) |
                (str_detect(Y_space,'Gut') & str_detect(X_space,'Plasma')) ~
                'Microbiome-Metabolome',
            str_detect(X_space,'Gut') & str_detect(Y_space,'Host') ~
                'Microbiome-Host Immune Response',
            str_detect(X_space,'Plasma') & str_detect(Y_space,'Host') ~
                'Metabolome-Host Immune Response'))
    
    return(plot.dat)
}

plot.arcdiagram <- function(plot.dat, osci.feats) {
    
    nodes <- plot.dat %>%
        select(to, from) %>%
        gather('point','name') %>%
        select(name) %>%
        distinct()
    node.spaces <- plot.dat %>%
        select(from, to, Y_space, X_space) %>%
        pivot_longer(cols = c(from, to), names_to = 'space', values_to = 'name') %>%
        mutate(space = case_when(
            space == 'from' ~ X_space,
            space == 'to' ~ Y_space)) %>%
        select(name, space) %>%
        distinct() %>%
        mutate(OSCI = ifelse(name %in% osci.feats$feature, 'bold', 'plain'))
    edges <- plot.dat %>%
        select(to, from, XY_effsize, XY_corr.dir, space)
    mygraph <- tidygraph::tbl_graph(directed = FALSE, node_key = 'name', 
                                    nodes = node.spaces, edges = edges)
    
    plot <- ggraph::ggraph(mygraph, layout = 'linear') +
        ggraph::geom_edge_arc(aes(color = factor(XY_corr.dir), 
                                  edge_alpha = abs(XY_effsize)),
                              strength = .75, fold = TRUE) +
        ggraph::geom_node_text(aes(label = name, fontface = OSCI), # repel = T,
                               angle = -90, nudge_y = -1, hjust = 0,
                               size = 2, family = 'Helvetica') +
        ggraph::geom_node_point(aes(color = space), size = 3) +
        ggraph::scale_edge_color_manual(values = c('#BF6130', # neg. corr,
                                                   '#245797')) +
        scale_color_manual(values = c(
            '#A2A048', # gut
            '#D9C35A', # host
            '#89A78F')) + # plasma
        ggraph::theme_graph()
}

prep.threeway.sankey <- function(df.res, df.names, fdr.cutoff = 0.05) {
    
    osci <- df.res %>%
        filter(str_detect(Covariate, 'OSCI')) %>%
        filter(Confounding %in% c('LD','SD','NC')) %>%
        distinct(across(c(run, feature, Confounding)),
                 .keep_all = TRUE)
    
    tmp <- df.res %>%
        filter(str_detect(Covariate,
                          'P_|U_')) %>%
        filter(str_detect(Confounding,
                          'S_|SDNA_|O_|OP_')) %>%
        filter(Confounding != 'NS' & FDR <= fdr.cutoff) %>%
        mutate(Confounding = str_split(Confounding,',')) %>%
        unnest_longer(Confounding) %>%
        mutate(Confounding = trimws(Confounding)) %>%
        filter(!str_detect(Confounding,'COVID|OSCI')) %>%
        mutate(Cov_Site = str_split(Covariate,'_') %>%
                   map_chr(~ head(., 1))) %>%
        mutate(Cov_Site = case_when(
            Cov_Site == 'OSCI' ~ 'OSCI',
            Cov_Site == 'VL' ~ 'SARS-CoV-2 Load',
            Cov_Site == 'U' ~ 'Urine metabolites',
            Cov_Site == 'P' ~ 'Plasma metabolites',
            Cov_Site %in% c('S','SDNA') ~ 'Gut microbiome',
            Cov_Site %in% c('O','OP') ~ 'Oropharyngeal microbiome',
            TRUE ~ 'Clinical')) %>%
        mutate(Conf_Site = str_split(Confounding,'_') %>%
                   map_chr(~ head(., 1))) %>%
        mutate(Conf_Site = case_when(
            Conf_Site == 'OSCI' ~ 'OSCI',
            Conf_Site == 'VL' ~ 'SARS-CoV-2 Load',
            Conf_Site == 'U' ~ 'Urine metabolites',
            Conf_Site == 'P' ~ 'Plasma metabolites',
            Conf_Site %in% c('S','SDNA') ~ 'Gut microbiome',
            Conf_Site %in% c('O','OP') ~ 'Oropharyngeal microbiome',
            Conf_Site %in% c('LD','SD','NC') ~ 'Robust',
            TRUE ~ 'Clinical')) %>%
        mutate_at(vars(c('Covariate','Confounding')),
                  ~ str_remove_all(.,
                                   'P_|OP_|S_|O_|SDNA_|U_')) %>%
        left_join(select(df.names,
                         R_Compound, Compound),
                  by = c('Covariate'='R_Compound')) %>%
        mutate(Covariate = ifelse(
            is.na(Compound), Covariate, Compound)) %>%
        select(-Compound) %>%
        left_join(select(df.names,
                         R_Compound, Compound),
                  by = c('Confounding'='R_Compound')) %>%
        mutate(Confounding = ifelse(
            is.na(Compound), Confounding, Compound)) %>%
        mutate(ES = case_when(
            `Effect Size` > 0 ~ 'Positive correlation',
            `Effect Size` < 0 ~ 'Negative correlation',
            TRUE ~ 'NA')) %>%
        select(-`Effect Size`) %>%
        left_join(select(osci, run, feature, `Effect Size`),
                  by = c('run','feature')) %>%
        mutate(ES_OSCI = ifelse(is.na(`Effect Size`),
                                'Not correlated',
                                'Positive correlation')) %>%
        # mutate(ES_OSCI = case_when(
        #     is.na(`Effect Size`) ~ 'Not correlated',
        #     `Effect Size` > 0 ~ 'Positive correlation',
        #     `Effect Size` < 0 ~ 'Negative correlation',
        #     TRUE ~ NA)) %>%
        select(-`Effect Size`) %>%
        mutate(Confounding = str_replace_all(
            Confounding, 'InvSimp|Shannon|Richness',
            'Alpha diversity')) %>%
        mutate(Confounding = str_replace_all(
            Confounding, 'AbxHosp|AbxSamp',
            'Current abx')) %>%
        mutate(Confounding = str_replace_all(
            Confounding, 'Abx3Mo',
            'Previous abx')) %>%
        mutate(Confounding = str_replace_all(
            Confounding, 'Comorbidity|Medication',
            'Comorbidity and meds')) %>%
        mutate(Confounding = str_replace_all(
            Confounding, 'Renal|Diuretics',
            'Renal disease')) %>%
        mutate(Confounding = str_replace_all(
            Confounding, 'Bacteremia|Pneumonia|DaysHospitalized',
            'HAP and hospitalization'))
}

plot.threeway.sankey <- function(df.plot) {
    
    df.plot %>%
        ggplot(aes(x = x,
                   next_x = next_x,
                   node = node,
                   next_node = next_node,
                   fill = factor(node),
                   label = node)) +
        ggsankey::geom_sankey(flow.alpha = 0.5, node.color = 1) +
        ggsankey::geom_sankey_label(size = 3.5, color = 1, fill = "white") +
        # scico::scale_fill_scico_d(palette = 'vikO') +
        scico::scale_fill_scico_d(palette = 'batlow') +
        #scico::scale_fill_scico_d(palette = 'roma') +
        ggsankey::theme_sankey(base_size = 14) +
        theme(legend.position = 'none') +
        labs(x = '')
}

#### NOT USED ####
# prep.osci.chord <- function(df.sdna, df.op, 
#                             df.host, df.metab,
#                             restricted = TRUE,
#                             osci.only = FALSE) {
#     
#     # plasma metabolites and clinical factors describing variation
#     # in GI and OP microbial taxa
#     pre.plot.mb.clin.metab <- df.sdna %>%
#         mutate(feat.space = 'S') %>%
#         bind_rows(df.op %>%
#                       mutate(feat.space = 'OP')) %>%
#         mutate(Confounding = str_split(Confounding, ',')) %>%
#         dplyr::filter(Confounding != 'NS') %>%
#         unnest_longer(Confounding) %>%
#         mutate_at(vars(Confounding), ~ trimws(.)) %>%
#         mutate(Confounding = case_when(
#             str_detect(Confounding,'OSCI|SD|LD|Worst') ~ 'OSCI',
#             str_detect(Confounding, 'P_') ~ Confounding,
#             str_detect(Confounding,'Abx|Antibiotics') ~ 'Abx',
#             str_detect(Confounding, 'Hosp|HAP|VAP|Pneumonia|Seps|Bactere') ~ 'Hosp. and HAP',
#             str_detect(Confounding,'Medic|Comorb|BMI') ~ 'Med. and Como.',
#             TRUE ~ 'Other')) %>%
#         dplyr::filter(Confounding != 'Other' & config == 'genus.all.metabs') %>%
#         distinct(across(c(feature, Confounding)), .keep_all = TRUE) %>%
#         count(feat.space, feature, Covariate, Confounding) 
#     plot.data.0 <- pre.plot.mb.clin.metab %>%
#         mutate(conf.space = case_when(
#             Confounding == 'OSCI' ~ 'OSCI',
#             str_detect(Confounding, 'P_') ~ 'P',
#             TRUE ~ Confounding)) %>%
#         select(feat.space, feature, conf.space, Confounding, n)
#     
#     # microbial taxa and clinical factors describing variation 
#     # in plasma metabolites
#     pre.plot.metab.mb <- df.metab %>%
#         mutate(Confounding = str_split(Confounding, ',')) %>%
#         unnest_longer(Confounding) %>%
#         mutate_at(vars(Confounding), ~ trimws(.)) %>%
#         mutate(Confounding = case_when(
#             str_detect(Confounding, 'SDNA_') ~ str_replace_all(
#                 Confounding,'SDNA_','S_'),
#             str_detect(Confounding, 'OP_') ~ str_replace_all(
#                 Confounding, 'OP_', 'O_'),
#             TRUE ~ Confounding)) %>%
#         mutate(Confounding = case_when(
#             str_detect(Confounding,'OSCI|SD|LD|Worst') ~ 'OSCI',
#             str_detect(Confounding, 'S_|O_') ~ Confounding,
#             str_detect(Confounding,'Abx|Antibiotics') ~ 'Abx',
#             str_detect(Confounding, 'Hosp|HAP|VAP|Pneumonia|Seps|Bactere') ~ 'Hosp. and HAP',
#             str_detect(Confounding,'Medic|Comorb|BMI') ~ 'Med. and Como.',
#             TRUE ~ 'Other')) %>%
#         # dplyr::filter(str_detect(
#         #     Confounding, 'OSCI|S_|O_|OP_|SDNA_')) %>%
#         # these are binned at genus level, prefiltered for osci assoc.
#         dplyr::filter(Confounding != 'Other' & config == 'metabs.motus.all') %>%
#         count(feature, Covariate, Confounding) 
#     
#     plot.data.1 <- pre.plot.metab.mb %>%
#         mutate(feat.space = str_split(feature, '_') %>%
#                    map_chr(~ head(., 1))) %>%
#         mutate(conf.space = str_split(Confounding, '_') %>%
#                    map_chr(~ head(., 1))) %>%
#         select(feat.space, feature, conf.space, Confounding, n)
#     
#     pre.plot.host.metab <- df.host %>%
#         mutate(Confounding = str_split(Confounding, ',')) %>%
#         unnest_longer(Confounding) %>%
#         mutate_at(vars(Confounding), ~ trimws(.)) %>%
#         mutate(Confounding = str_replace_all(Confounding, 'SD|LD|NC','OSCI')) %>%
#         dplyr::filter(str_detect(Confounding, 'OSCI|P_|U_')) %>%
#         dplyr::filter(config == 'ck.all') %>%
#         count(feature, Covariate, Confounding) 
#     
#     plot.data.2 <- pre.plot.host.metab %>%
#         mutate(feat.space = 'host') %>%
#         mutate(conf.space = str_split(Confounding, '_') %>%
#                    map_chr(~ head(., 1))) %>%
#         select(feat.space, feature, conf.space, Confounding, n)
#     
#     plot.data.tmp <- plot.data.0 %>%
#         bind_rows(plot.data.1,
#                   plot.data.2)
#     
#     # plot.data.tmp %>%
#     #     mutate(Confounding=str_remove_all(
#     #         Confounding,'S_|O_|P_')) %>%
#     #     mutate(feature=str_remove_all(
#     #         feature,'S_|O_|P_')) %>%
#     #     mutate(feat.robust=ifelse(
#     #         Confounding=='OSCI',TRUE,FALSE)) %>%
#     #     mutate(feat.is.conf=ifelse(
#     #         feature %in% .$Confounding, TRUE, FALSE)) %>%
#     #     dplyr::filter(dual.space)
#     
#     if (restricted) {
#         plot.data.tmp <- plot.data.tmp %>%
#             dplyr::filter(!feat.space %in% c('U','OP')) %>%
#             dplyr::filter(conf.space != 'O') # %>%
#             # mutate(tmp = paste0(feat.space,'-',conf.space)) %>%
#             # dplyr::filter(tmp != 'S-P') %>%
#             # select(-tmp) 
#     } 
#     
#     if (osci.only) {
#         plot.data.tmp <- plot.data.tmp %>%
#             dplyr::filter(!(Confounding %in% c(
#                 'Abx','Hosp. and HAP','Med. and Como.')))
#     }
#       
#     plot.data.out <- plot.data.tmp %>%
#         group_by(feat.space, feature, conf.space, Confounding) %>%
#         mutate(n = sum(n)) %>%
#         ungroup() %>%
#         distinct(across(c(feat.space, feature, 
#                           conf.space, Confounding)),
#                  .keep_all = TRUE) %>%
#         count(feat.space, conf.space) %>%
#         spread(conf.space, 'n') %>%
#         mutate_if(is.numeric, ~ replace_na(., 0)) %>% 
#         gather('conf.space','n',-feat.space) %>% 
#         distinct(across(c(feat.space, conf.space)), .keep_all = TRUE) %>%
#         mutate(feat.space = case_when(
#             feat.space=='P'~'Plasma metabolites',
#             feat.space=='U'~'Urine metabolites',
#             feat.space=='S'~'Gut microbiome',
#             feat.space=='OP'~'Oropharyngeal microbiome',
#             feat.space=='host'~'Host immune system',
#             TRUE ~ 'Clinical factors')) %>%
#         mutate(conf.space = case_when(
#             conf.space=='P'~'Plasma metabolites',
#             conf.space=='S'~'Gut microbiome',
#             conf.space=='O'~'Oropharyngeal microbiome',
#             TRUE ~ conf.space)) %>%
#         spread(conf.space,'n') %>%
#         mutate_if(is.numeric, ~ replace_na(., 0)) %>%
#         column_to_rownames('feat.space') %>%
#         as.matrix()
#     
#     if (restricted) {
#         group.names <- c('Gut microbiome','Plasma metabolites',
#                          'Host immune system','Abx','OSCI','Hosp. and HAP',
#                          'Med. and Como.')
#         group.names <- c('Med. and Como.','Abx','OSCI','Hosp. and HAP',
#                          'Gut microbiome','Plasma metabolites',
#                          'Host immune system')
#         group <- c('Microbiome','Metabolome','Host Immune Response',
#                    'Clinical','Clinical','Clinical','Clinical')
#         group <- c('Clinical','Clinical','Clinical','Clinical',
#                    'Microbiome','Metabolome','Host Immune Response')
#         names(group) <- group.names
#     } else {
#         group.names <- c('Gut microbiome','Oropharyngeal microbiome',
#                          'Plasma metabolites','Urine metabolites',
#                          'Host immune system','OSCI','Abx',
#                          'Med. and Como.','Hosp. and HAP')
#         group <- c('Microbiome','Microbiome','Metabolome','Metabolome',
#                    'Host Immune Response','Clinical','Clinical',
#                    'Clinical','Clinical')
#         names(group) <- group.names
#     }
#     
#     return(list('data' = plot.data.out,
#                 'group.labels' = group))
# }
# 
# plot.osci.chord <- function(df.plot, group.names, fn,
#                             directed = TRUE) {
#     
#     pal <- ggthemes::tableau_color_pal('Miller Stone')(10)
#     # color by status
#     colors.df <- group.names %>%
#         enframe() %>%
#         mutate(color = case_when(str_detect(value, 'Microbiome') ~ pal[1],
#                                  str_detect(value, 'Metabolome') ~ pal[7],
#                                  str_detect(value, 'Immune') ~ pal[8],
#                                  str_detect(value, 'Clinical') ~ pal[5])) %>%
#         mutate(color2 = case_when(str_detect(name, 'microbiome') ~ pal[1],
#                                   str_detect(name, 'metabolites') ~ pal[7],
#                                   str_detect(name, 'immune') ~ pal[8],
#                                   TRUE ~ pal[5]))
#     colors <- colors.df$color
#     names(colors) <- colors.df$value
#     
#     colors <- colors.df$color2
#     names(colors) <- colors.df$name
#     
#     # remove microbial variation explained by metabolites
#     if (str_detect(fn, 'sel|only')) {
#         df.plot['Gut microbiome',
#                 'Plasma metabolites'] <- 0 }
# 
#     pdf(here('figures','tmp',fn))
#     circos.clear()
#     circos.par(start.degree = 270)
#     
#     if (directed) {
#         plot <- 
#             chordDiagram(df.plot, 
#                          group = group.names,
#                          # col = colors,
#                          transparency = 0.5,
#                          directional = -1,
#                          annotationTrack = c("name", "grid"),
#                          direction.type = c("diffHeight", "arrows"),
#                          link.arr.type = "big.arrow", scale = TRUE, 
#                          big.gap = 15)
#     } else {
#         plot <- 
#             chordDiagram(df.plot, 
#                          group = group.names,
#                          # col = colors,
#                          transparency = 0.5,
#                          directional = -1,
#                          annotationTrack = c("name", "grid"),
#                          link.arr.type = "big.arrow", scale = TRUE, 
#                          big.gap = 15)
#     }
#     
#     # browser()
#     # will only return df of coordinates/info...
#     dev.off()
#     return(plot)
# }
# host features -> metabolites -> microbiota/clinical factors

# 
# plot.mediation.network <- function(df.med) {
#     # https://www.data-imaginist.com/2017/ggraph-introduction-layouts/
#     library(ggraph)
#     library(igraph)
#     graph <- graph_from_data_frame(highschool)
#     V(graph)$friends <- degree(graph, mode = 'in')
#     V(graph)$friends <- ifelse(V(graph)$friends < 5, 'few', 
#                                ifelse(V(graph)$friends >= 15, 'many', 'medium'))
#     ggraph(graph, 'hive', axis = 'friends', sort.by = 'degree') + 
#         geom_edge_hive(aes(colour = factor(year), alpha = ..index..)) + 
#         geom_axis_hive(aes(colour = friends), size = 3, label = FALSE) + 
#         coord_fixed()
# }
# plot.mediation.hive <- function(met.res, mb.res) {
#     
#     microb.med.metab.combine %>%
#         select(-Data) %>%
#         mutate(Effect = case_when(
#             `Effect Size` > 0 ~ '+ correlation',
#             `Effect Size` < 0 ~ '- correlation',
#             TRUE ~ 'no correlation')) %>%
#         # select(feature,Covariate,`Effect`) %>%
#         # as_tbl_graph(nodes = c('OSCI_Score_Samp', unique(.$feature)),
#         #              edges = select(., feature, Covariate)) %>%
#         as_tbl_graph(nodes = c(unique(.$feature), unique(.$Confounding))) %>%
#         # mutate(rank = -1*centrality_degree()) %>%
#         mutate(axis = case_when(
#             str_detect(name,'S_|O_|SDNA_|OP_') ~ 'Microbiome',
#             str_detect(name,'P_|U_') ~ 'Metabolome',
#             name=='OSCI_Score_Samp' ~ 'OSCI')) %>%
#         ggraph::ggraph(layout = 'linear', circular = TRUE) +
#         geom_edge_arc(aes(color = factor(Effect), label = Confounding)) +
#         #geom_node_point(aes(color = axis, label = name))
#         geom_node_text(aes(label = name, color = axis)) +
#         theme_graph()
#     
#     ggraph::ggraph('hive', axis = axis) +
#         geom_edge_hive(aes(color = factor(Effect), label = Confounding,
#                            #width = abs(`Effect Size`)
#         )) +
#         geom_axis_hive(aes(color = axis), size = 3, label = FALSE) +
#         scale_color_manual(values = wesanderson::wes_palette('Darjeeling1', 3,
#                                                              'discrete')) +
#         coord_fixed()
# }