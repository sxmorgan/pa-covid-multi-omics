## S1-cohort-characteristics.tsv
* *All yes/no variables besides Sex were encoded as yes=1 and no=0; NA = no value recorded*
* Demographics
  * Age, Sex, BMI, (active) Smoker
* Comorbidities and health conditions
  * Diabetes, IBD, (chronic) Hypertension, Dyslipidemia, Arrythmia, Coronary (disease), (active) Neoplasia
  * Renal = chronic renal insufficiency
  * Liver = hronic liver disease
  * Thyroid = altered thyroid hormones
  * Lung = chronic lung disease
  * Pregnancy
  * Charlson Comorbidity Index (CCI)
* Medications
  * Immunosuppressives, Metformin, Corticosteroids, Statins, Antiaggregants, Antihypertensives, ß-Blockers, Uricostatics, Diuretics, (nutritional) Supplements
    * Medication = sum/10 above 
  * Antivirals (during hospitalization), Antifungals
  * Abx3Mo = antibiotics in the last 3 months
  * AbxSamp = antibiotics during sampling period
  * AbxFree = antibiotic-free (before and during hospitalization)
* Course of COVID-19
  * GI symptoms, OSCI (severity score), Death
  * DaysHospitalized = days of hospitalization due to COVID-19
* Secondary infections during hospitalization
  * HAP = hospital-acquired pneumonia
  * VAP = ventilator-assisted pneumonia
  * Sepsis, Bacteremia, Urinary tract infection (UTI), Infection (cause unclear)

## S2-clinical-lab-parameters.tsv
  * OSCI_Score_Samp = severity score at time of sampling based on clinical parameters
  * OSCI_Class = severity class (0=control, 1-4=mild, 5-8=severe)
  * CRP (mg/l), Procalcitonin (ug/L), Ferritin (ug/l), LDH (U/l), Hemoglobin (g/dl), Leucocytes (/nl), Platelets (/nl), Neutrophils absolute (/nl), Neutrophils %, Immature Granulocytes absolute (/nl), Immature Granulocytes %, Lymphocytes absolute (/nl), Lymphocytes %, Monocytes absolute (/nl), Monocytes %, D-dimers (mg/L), aPTT (sec.), Fibrinogen (g/L), eGFR, ALT (U/L), AST (U/L), gammaGT (U/L), NT-pro BNP (ng/l)

## S3-taxonomic-profiles.tsv
  * Raw genus and mOTUs level relative abundance profiles from microbiota WGS
  * OP = Oropharyngeal swab, SDNA = Stool, TBS = tracheobronchial secretion (from VAP patients)

## S4-immune-profiles.tsv
  * Host_PBMC features from scRNA-seq (post-processed), in form of gene_celltype
  * Host_Cytokine features from multiplex ELISA (post-processed)

## S5-metabolite-profiles.tsv
  * CCM = targeted central carbon metabolism
  * TRYP = targeted tryptophan metabolism, measured in µg/mL
  * Q500LC = , measured in µM
  * Q500FIA = , mesured in µM

## S6-clinical-confounder-analysis-results.tsv

## S7-integrated-analysis-results.tsv
