library(TwoSampleMR)
library(readr)
library(MRPRESSO)
library(ieugwasr)
library(MendelianRandomization)
library(MRlap)
library(meta)
library(MungeSumstats)
library(SNPlocs.Hsapiens.dbSNP144.GRCh37)
library(SNPlocs.Hsapiens.dbSNP144.GRCh38)
library(BSgenome.Hsapiens.1000genomes.hs37d5)
library(BSgenome.Hsapiens.NCBI.GRCh38)

### East Asian (EAS) and African American (AFA) GWAS of T2D (PMID:38374256)
### Download original EAS and AFA T2D GWAS data from the website https://diagram-consortium.org/downloads.html

### format EAS T2D data using MungeSumstats
### we aimed to add the rsid to EAS T2D summary-level data

T2D <- read.table('F:/diabetes_VTE/EAS_Metal_LDSC-CORR_Neff.v2.txt', sep='\t', header = T, stringsAsFactors = F)
T2D$N <- T2D$Ncases + T2D$Ncontrols
T2D <- T2D[,-c(9:11)]
colnames(T2D) <- c('CHR', 'POS', 'Effect_allele', 'Other_allele', 'Beta', 'SE', 'EAF',
         'P_value', 'N')
format_sumstats(path = T2D, ref_genome = "GRCh37", dbSNP=144, save_path = "F:/diabetes_VTE/T2D.EAS.tsv.gz")

### format AFA T2D data using MungeSumstats
### we aimed to add the rsid to AFA T2D summary-level data

T2D <- read.table('F:/diabetes_VTE/AFA_Metal_LDSC-CORR_Neff.v2.txt', sep='\t', header = T, stringsAsFactors = F)
T2D$N <- T2D$Ncases + T2D$Ncontrols
T2D <- T2D[,-c(9:11)]
colnames(T2D) <- c('CHR', 'POS', 'Effect_allele', 'Other_allele', 'Beta', 'SE', 'EAF',
         'P_value', 'N')
format_sumstats(path = T2D, ref_genome = "GRCh37", dbSNP=144, save_path = "F:/diabetes_VTE/T2D.AFA.tsv.gz")

#### select instrumental variables for T2D in EAS and AFA ####

### EAS ###
T2D <- read.table('F:/diabetes_VTE/T2D.EAS.tsv.gz', sep = '\t', stringsAsFactors = F, header = T)
### genome-wide significant snps
T2D <- T2D[T2D$P<5e-8, ]
### clumping ###
dat <- data.frame(rsid=T2D$SNP, pval=T2D$P)
retained_snps <- ld_clump(dat, clump_kb = 10000, clump_r2 = 0.001, bfile = 'D:/LD_reference/EAS',
                          plink_bin = 'D:/LD_reference/plink.exe')
T2D <- T2D[T2D$SNP%in%retained_snps$rsid, ]
### Exclude SNPs with minor allele frequency <0.01 ###
T2D <- T2D[(T2D$FRQ>0.01)&(T2D$FRQ<0.99), ]
### Estimate the F-statistics (PMID:30861319)
### Exclude SNPs with F-statistics<10
T2D$F_stat <- (T2D$BETA/T2D$SE)^2
T2D <- T2D[T2D$F_stat>10, ]

write.csv(T2D, 'F:/diabetes_VTE/T2D_IVs_EA.csv', quote = F, row.names = F)

### AFA ###
T2D <- read.table('F:/diabetes_VTE/T2D.AFA.tsv.gz', sep = '\t', stringsAsFactors = F, header = T)

### genome-wide significant snps
T2D <- T2D[T2D$P<5e-8, ]

### clumping ###
dat <- data.frame(rsid=T2D$SNP, pval=T2D$P)
retained_snps <- ld_clump(dat, clump_kb = 10000, clump_r2 = 0.001, bfile = 'D:/LD_reference/AFR',
                          plink_bin = 'D:/LD_reference/plink.exe')
T2D <- T2D[T2D$SNP%in%retained_snps$rsid, ]

### Exclude SNPs with minor allele frequency <0.01 ###
T2D <- T2D[(T2D$FRQ>0.01)&(T2D$FRQ<0.99), ]

### Estimate the F-statistics (PMID:30861319)
### Exclude SNPs with F-statistics<10
T2D$F_stat <- (T2D$BETA/T2D$SE)^2
T2D <- T2D[T2D$F_stat>10, ]

write.csv(T2D, 'F:/diabetes_VTE/T2D_IVs_AFA.csv', quote = F, row.names = F)


### East Asian (EAS) and African American (AFA) GWAS of VTE (PMID:36777996)
### Download original EAS and AFA VTE GWAS data from the website https://www.globalbiobankmeta.org/sharable-links

#### select instrumental variables for VTE in EAS and AFA ####

### EAS
VTE <- read.table(gzfile('F:/diabetes_VTE/VTE_Bothsex_eas_inv_var_meta_GBMI_052021_nbbkgt1.txt.gz'), sep='\t', header=F,
                  stringsAsFactors = F)
VTE <- VTE[,-c(14:17)]
colnames(VTE) <- c('CHR', 'POS', 'REF', 'ALT', 'rsid', 'all_meta_AF', 'inv_var_meta_beta', 
                   'inv_var_meta_sebeta', 'inv_var_meta_p', 'inv_var_het_p', 'direction', 'N_case', 'N_ctrl')
### genome-wide significant snps
VTE <- VTE[VTE$inv_var_meta_p<5e-6, ]
### only one snp, we do not need to clump ###
### Exclude SNPs with minor allele frequency <0.01 ###
VTE <- VTE[(VTE$all_meta_AF>0.01)&(VTE$all_meta_AF<0.99), ]
### Estimate the F-statistics (PMID:30861319)
### Exclude SNPs with F-statistics<10
VTE$F_stat <- (VTE$inv_var_meta_beta/VTE$inv_var_meta_sebeta)^2
VTE <- VTE[VTE$F_stat>10, ]

write.csv(VTE, 'F:/diabetes_VTE/VTE_IVs_EAS.csv', quote = F, row.names = F)

### AFA
VTE <- read.table(gzfile('F:/diabetes_VTE/VTE_Bothsex_afr_inv_var_meta_GBMI_052021_nbbkgt1.txt.gz'), sep='\t', header=F,
                  stringsAsFactors = F)
VTE <- VTE[,-c(14:17)]
colnames(VTE) <- c('CHR', 'POS', 'REF', 'ALT', 'rsid', 'all_meta_AF', 'inv_var_meta_beta', 
                   'inv_var_meta_sebeta', 'inv_var_meta_p', 'inv_var_het_p', 'direction', 'N_case', 'N_ctrl')
### genome-wide significant snps
VTE <- VTE[VTE$inv_var_meta_p<5e-6, ]
### clumping ###
dat <- data.frame(rsid=VTE$SNP, pval=VTE$P)
retained_snps <- ld_clump(dat, clump_kb = 10000, clump_r2 = 0.001, bfile = 'D:/LD_reference/AFR',
                          plink_bin = 'D:/LD_reference/plink.exe')
VTE <- VTE[VTE$SNP%in%retained_snps$rsid, ]
### Exclude SNPs with minor allele frequency <0.01 ###
VTE <- VTE[(VTE$all_meta_AF>0.01)&(VTE$all_meta_AF<0.99), ]
### Estimate the F-statistics (PMID:30861319)
### Exclude SNPs with F-statistics<10
VTE$F_stat <- (VTE$inv_var_meta_beta/VTE$inv_var_meta_sebeta)^2
VTE <- VTE[VTE$F_stat>10, ]

write.csv(VTE, 'F:/diabetes_VTE/VTE_IVs_AFA.csv', quote = F, row.names = F)



#### East Asian T2D ---> VTE ####
T2D <- read.csv('F:/diabetes_VTE/T2D_IVs_EA.csv', sep=',', stringsAsFactors = F, header = T)
T2D <- format_data(T2D, chr_col = 'CHR', pos_col = 'BP', effect_allele_col = 'A2',
                   other_allele_col = 'A1', eaf_col = 'FRQ', beta_col = 'BETA', se_col = 'SE',
                   pval_col = 'P', samplesize_col = 'N', snp_col = 'SNP')

VTE <- read.table(gzfile('VTE_Bothsex_eas_inv_var_meta_GBMI_052021_nbbkgt1.txt.gz'), sep='\t', header = T, stringsAsFactors = F)
VTE <- VTE[,-c(14:17)]
colnames(VTE) <- c('CHR', 'POS', 'REF', 'ALT', 'rsid', 'all_meta_AF', 'inv_var_meta_beta', 'inv_var_meta_sebeta',
                   'inv_var_meta_p', 'inv_var_het_p', 'direction', 'N_case', 'N_ctrl')
VTE <- VTE[VTE$rsid%in%T2D$SNP, ]
VTE <- format_data(VTE, type='outcome', chr_col = 'CHR', pos_col = 'POS', snp_col='rsid', effect_allele_col = 'ALT', other_allele_col = 'REF',
                   eaf_col='all_meta_AF', beta_col = 'inv_var_meta_beta', se_col='inv_var_meta_sebeta',
                   pval_col='inv_var_meta_p', samplesize_col = 'N')
T2D.VTE.EAS <- harmonise_data(T2D, VTE)
write.csv(T2D.VTE.EAS, 'F:/diabetes_VTE/T2D_to_VTE_EA.csv', quote=F, row.names = F)

T2D.VTE.EA <- read.csv('F:/diabetes_VTE/T2D_to_VTE_EA.csv', sep=',', stringsAsFactors = F, header = T)
ivw_res.VTE.EA <- mr(T2D.VTE.EA, method_list = c('mr_ivw', 'mr_egger_regression', 'mr_weighted_median'));ivw_res.VTE.EA
#### To better interpret the result, we multiplied the causal estimate by 0.693 to reflect 
#### the increase in the risk of VTE associated with each doubling of the odds of genetic predisposition to T2D
ivw_res.VTE.EA$b <- 0.693*ivw_res.VTE.EA$b
ivw_res.VTE.EA$se <- 0.693*ivw_res.VTE.EA$se
generate_odds_ratios(ivw_res.VTE.EA)



#### East Asian VTE ---> T2D ####
VTE <- read.csv('F:/diabetes_VTE/VTE_IVs_EAS.csv', sep=',', stringsAsFactors = F, header = T)
VTE <- format_data(VTE, chr_col = 'CHR', pos_col = 'POS', snp_col='rsid', effect_allele_col = 'ALT', other_allele_col = 'REF',
                   eaf_col='all_meta_AF', beta_col = 'inv_var_meta_beta', se_col='inv_var_meta_sebeta',
                   pval_col='inv_var_meta_p', samplesize_col = 'N')

T2D <- read.table('F:/diabetes_VTE/T2D.EAS.tsv.gz', sep='\t', header = T, stringsAsFactors = F)
T2D <- T2D[T2D$SNP%in%VTE$SNP, ]
T2D <- format_data(T2D, type='outcome', chr_col = 'CHR', pos_col = 'BP', effect_allele_col = 'A2',
                   other_allele_col = 'A1', eaf_col = 'FRQ', beta_col = 'BETA', se_col = 'SE',
                   pval_col = 'P', samplesize_col = 'N', snp_col = 'SNP')
VTE.T2D.EAS <- harmonise_data(VTE, T2D)
write.csv(VTE.T2D.EAS, 'F:/diabetes_VTE/VTE_to_T2D_EAS.csv', quote=F, row.names = F)

VTE.T2D.EAS <- read.csv('F:/diabetes_VTE/VTE_to_T2D_EAS.csv', sep=',', stringsAsFactors = F, header = T)
ivw_res.T2D.EAS <- mr(VTE.T2D.EAS, method_list = c('mr_wald_ratio'));ivw_res.T2D.EAS
#### To better interpret the result, we multiplied the causal estimate by 0.693 to reflect 
#### the increase in the risk of CHD associated with each doubling of the odds of genetic predisposition to T2D
ivw_res.T2D.EAS$b <- 0.693*ivw_res.T2D.EAS$b
ivw_res.T2D.EAS$se <- 0.693*ivw_res.T2D.EAS$se
generate_odds_ratios(ivw_res.T2D.EAS)

#### African T2D ---> VTE ####
T2D <- read.csv('F:/diabetes_VTE/T2D_IVs_AFA.csv', sep=',', stringsAsFactors = F, header = T)
T2D <- format_data(T2D, chr_col = 'CHR', pos_col = 'BP', effect_allele_col = 'A2',
                   other_allele_col = 'A1', eaf_col = 'FRQ', beta_col = 'BETA', se_col = 'SE',
                   pval_col = 'P', samplesize_col = 'N', snp_col = 'SNP')

VTE <- read.table(gzfile('VTE_Bothsex_afr_inv_var_meta_GBMI_052021_nbbkgt1.txt.gz'), sep='\t', header = T, stringsAsFactors = F)
VTE <- VTE[,-c(14:17)]
colnames(VTE) <- c('CHR', 'POS', 'REF', 'ALT', 'rsid', 'all_meta_AF', 'inv_var_meta_beta', 'inv_var_meta_sebeta',
                   'inv_var_meta_p', 'inv_var_het_p', 'direction', 'N_case', 'N_ctrl')
VTE <- VTE[VTE$rsid%in%T2D$SNP, ]
VTE <- format_data(VTE, type='outcome', chr_col = 'CHR', pos_col = 'POS', snp_col='rsid', effect_allele_col = 'ALT', other_allele_col = 'REF',
                   eaf_col='all_meta_AF', beta_col = 'inv_var_meta_beta', se_col='inv_var_meta_sebeta',
                   pval_col='inv_var_meta_p', samplesize_col = 'N')
T2D.VTE.EAS <- harmonise_data(T2D, VTE)
write.csv(T2D.VTE.EAS, 'F:/diabetes_VTE/T2D_to_VTE_AFA.csv', quote=F, row.names = F)

T2D.VTE.AFR <- read.csv('F:/diabetes_VTE/T2D_to_VTE_AFA.csv', sep=',', stringsAsFactors = F, header = T)
ivw_res.VTE.AFR <- mr(T2D.VTE.AFR, method_list = c('mr_ivw', 'mr_egger_regression', 'mr_weighted_median'));ivw_res.VTE.EA
#### To better interpret the result, we multiplied the causal estimate by 0.693 to reflect 
#### the increase in the risk of CHD associated with each doubling of the odds of genetic predisposition to T2D
ivw_res.VTE.AFR$b <- 0.693*ivw_res.VTE.AFR$b
ivw_res.VTE.AFR$se <- 0.693*ivw_res.VTE.AFR$se
generate_odds_ratios(ivw_res.VTE.AFR)

#### African VTE ---> T2D ####
VTE <- read.csv('F:/diabetes_VTE/VTE_IVs_AFA.csv', sep=',', stringsAsFactors = F, header = T)
VTE <- format_data(VTE, chr_col = 'CHR', pos_col = 'POS', snp_col='rsid', effect_allele_col = 'ALT', other_allele_col = 'REF',
                   eaf_col='all_meta_AF', beta_col = 'inv_var_meta_beta', se_col='inv_var_meta_sebeta',
                   pval_col='inv_var_meta_p', samplesize_col = 'N')

T2D <- read.table('F:/diabetes_VTE/T2D.AFA.tsv.gz', sep='\t', header = T, stringsAsFactors = F)
T2D <- T2D[T2D$SNP%in%VTE$SNP, ]
T2D <- format_data(T2D, type='outcome', chr_col = 'CHR', pos_col = 'BP', effect_allele_col = 'A2',
                   other_allele_col = 'A1', eaf_col = 'FRQ', beta_col = 'BETA', se_col = 'SE',
                   pval_col = 'P', samplesize_col = 'N', snp_col = 'SNP')
VTE.T2D.AFR <- harmonise_data(VTE, T2D)
write.csv(VTE.T2D.AFR, 'F:/diabetes_VTE/VTE_to_T2D_AFR.csv', quote=F, row.names = F)

VTE.T2D.AFR <- read.csv('F:/diabetes_VTE/VTE_to_T2D_AFR.csv', sep=',', stringsAsFactors = F, header = T)
ivw_res.T2D.AFR <- mr(VTE.T2D.AFR, method_list = c('mr_ivw', 'mr_egger_regression', 'mr_weighted_median'));ivw_res.T2D.AFR
#### To better interpret the result, we multiplied the causal estimate by 0.693 to reflect 
#### the increase in the risk of CHD associated with each doubling of the odds of genetic predisposition to T2D
ivw_res.T2D.AFR$b <- 0.693*ivw_res.T2D.AFR$b
ivw_res.T2D.AFR$se <- 0.693*ivw_res.T2D.AFR$se
generate_odds_ratios(ivw_res.T2D.AFR)


##### MR-lap method #####

##### East Asian VTE data processing ####
VTE <- read.table(gzfile("F:/diabetes_VTE/VTE_Bothsex_eas_inv_var_meta_GBMI_052021_nbbkgt1.txt.gz"), sep='\t', header = F)
VTE <- VTE[,-c(14:17)]
colnames(VTE) <- c('CHR', 'POS', 'Other_allele', 'effect_allele', 'SNP', 'EAF', 'Beta', 'SE',
                      'P_value', 'inv_var_het_p', 'direction', 'N_case', 'N_ctrl')
VTE <- VTE[(VTE$CHR)!=23, ]
VTE$N <- VTE$N_case + VTE$N_ctrl
VTE <- VTE[,-c(10,11,12,13)]
format_sumstats(path = VTE, ref_genome = "GRCh38", dbSNP=144, save_path = "VTE.EA.tsv.gz")

##### African VTE data processing ####
VTE <- read.table(gzfile("F:/diabetes_VTE/VTE_Bothsex_afr_inv_var_meta_GBMI_052021_nbbkgt1.txt.gz"), sep='\t', header = F)
VTE <- VTE[,-c(14:17)]
colnames(VTE) <- c('CHR', 'POS', 'Other_allele', 'effect_allele', 'SNP', 'EAF', 'Beta', 'SE',
                   'P_value', 'inv_var_het_p', 'direction', 'N_case', 'N_ctrl')
VTE <- VTE[(VTE$CHR)!=23, ]
VTE$N <- VTE$N_case + VTE$N_ctrl
VTE <- VTE[,-c(10,11,12,13)]
format_sumstats(path = VTE, ref_genome = "GRCh38", dbSNP=144, save_path = "VTE.AFA.tsv.gz")

#### because MRlap function exclude the SNPs located within the MHC region
#### so we modified the MRlap function to retain the MHC region
#### the new function called MRlap_EAS.r

source('F:/Asthma_CHD/MRlap_EAS.r')

### East Asian VTE results ###
T2D <- read.table('F:/diabetes_VTE/T2D.EAS.tsv.gz', header = T, stringsAsFactors = F)
VTE <- read.table('F:/diabetes_VTE/VTE.EA.tsv.gz', header = T, stringsAsFactors = F)
T2D <- T2D[T2D$SNP%in%VTE$SNP, ]
VTE <- VTE[VTE$SNP%in%T2D$SNP, ]

colnames(VTE) <- c('snp', 'chr', 'pos', 'a2', 'a1', 'freq', 'beta', 'se', 'P', 'N')

pos.hg38 <- data.frame(SNP=VTE$snp, Chr_hg38=VTE$chr, pos_hg38=VTE$pos)
T2D <- merge(T2D, pos.hg38, by='SNP')
T2D <- T2D[,-c(2,3)]
colnames(T2D) <- c('snp',  'a2', 'a1',  'beta', 'se', 'freq', 'P', 'N', 'chr', 'pos')

res1 <- MRlap(exposure = T2D,
              exposure_name = 'T2D',
              outcome = VTE,
              outcome_name = 'VTE',
              ld = 'D:/LD_reference/eas_w_ld_chr',
              hm3 = 'D:/LD_reference/w_hm3.noMHC.snplist',
              MR_pruning_LD = 0.001,
              MR_pruning_dist = 10000,
              MR_reverse = 0.001)
res1$MRcorrection$corrected_effect <- 0.693*res1$MRcorrection$corrected_effect
res1$MRcorrection$corrected_effect_se <- 0.693*res1$MRcorrection$corrected_effect_se
exp(res1$MRcorrection$corrected_effect)
lo_ci <- res1$MRcorrection$corrected_effect - qnorm(0.975)*res1$MRcorrection$corrected_effect_se
up_ci <- res1$MRcorrection$corrected_effect + qnorm(0.975)*res1$MRcorrection$corrected_effect_se
exp(c(lo_ci, up_ci))


#### MRlap for African ####

source('F:/Asthma_CHD/MRlap_AFR.r')

### African VTE results ###
T2D <- read.table('F:/diabetes_VTE/T2D.AFA.tsv.gz', header = T, stringsAsFactors = F)
VTE <- read.table('F:/diabetes_VTE/VTE.AFA.tsv.gz', header = T, stringsAsFactors = F)
T2D <- T2D[T2D$SNP%in%VTE$SNP, ]
VTE <- VTE[VTE$SNP%in%T2D$SNP, ]

colnames(VTE) <- c('snp', 'chr', 'pos', 'a2', 'a1', 'freq', 'beta', 'se', 'P', 'N')

pos.hg38 <- data.frame(SNP=VTE$snp, Chr_hg38=VTE$chr, pos_hg38=VTE$pos)
T2D <- merge(T2D, pos.hg38, by='SNP')
T2D <- T2D[,-c(2,3)]
colnames(T2D) <- c('snp',  'a2', 'a1',  'beta', 'se', 'freq', 'P', 'N', 'chr', 'pos')

res2 <- MRlap(exposure = T2D,
              exposure_name = 'T2D',
              outcome = VTE,
              outcome_name = 'VTE',
              ld = 'D:/LD_reference/eas_w_ld_chr',
              hm3 = 'D:/LD_reference/w_hm3.noMHC.snplist',
              MR_pruning_LD = 0.001,
              MR_pruning_dist = 10000,
              MR_reverse = 0.001)
res2$MRcorrection$corrected_effect <- 0.693*res2$MRcorrection$corrected_effect
res2$MRcorrection$corrected_effect_se <- 0.693*res2$MRcorrection$corrected_effect_se
exp(res2$MRcorrection$corrected_effect)
lo_ci <- res2$MRcorrection$corrected_effect - qnorm(0.975)*res2$MRcorrection$corrected_effect_se
up_ci <- res2$MRcorrection$corrected_effect + qnorm(0.975)*res2$MRcorrection$corrected_effect_se
exp(c(lo_ci, up_ci))

res3 <- MRlap(exposure = VTE,
              exposure_name = 'VTE',
              outcome = T2D,
              outcome_name = 'T2D',
              ld = 'D:/LD_reference/eas_w_ld_chr',
              hm3 = 'D:/LD_reference/w_hm3.noMHC.snplist',
              MR_pruning_LD = 0.001,
              MR_pruning_dist = 10000,
              MR_threshold = 5e-06)
res3$MRcorrection$corrected_effect <- 0.693*res3$MRcorrection$corrected_effect
res3$MRcorrection$corrected_effect_se <- 0.693*res3$MRcorrection$corrected_effect_se
exp(res3$MRcorrection$corrected_effect)
lo_ci <- res3$MRcorrection$corrected_effect - qnorm(0.975)*res3$MRcorrection$corrected_effect_se
up_ci <- res3$MRcorrection$corrected_effect + qnorm(0.975)*res3$MRcorrection$corrected_effect_se
exp(c(lo_ci, up_ci))