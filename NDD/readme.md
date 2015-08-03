## Predicting NDD risk by de novo mutations and extra phenotypes at birth

## Workflow (updated 07/28/2015)

## 1. Prepare risk genes list
# Subset developmental disorder genes list by imposing cutoff of FDR (p.adj) <= 0.4
python split_fdr.py external_DD.txt 0 0.4 > external_DD.cutoff.txt

# Take union of ASD risk genes (De Rubeis et al.) and subsetted genes
python union.py asd_derubeis.txt external_DD.cutoff.txt > riskgenes.txt

## 2. Run annotation script
python highrisk.py cases_metaSVM.csv phenotypecombined.csv phenotypes_extracardiac_mb_wkc.csv riskgenes.txt 0.6 > annos2.csv

## 3. Run classification on patients with certain diagnosis (NDD/non-NDD)
awk -F ',' '{if($4!="-1") print $0}' annos2.csv > annos2.certain.csv
--> run logitNDD5.R, svmNDD5.R, rfNDD4.R
