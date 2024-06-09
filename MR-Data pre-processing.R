


# tiredness ---------------------------------------------------------------


## tiredness
vcf_data <- vcfR::read.vcfR("data/ukb-b-929.vcf.gz")

meta <- data.frame(vcf_data@meta)
fix <- data.frame(vcf_data@fix) 
gt <- data.frame(vcf_data@gt) 
head(gt)

gt_df <- strsplit(gt[[names(gt)[2]]], split=":")
gt_df <- as.data.frame(do.call(rbind, gt_df))
names(gt_df) <- c("Beta","SE","log10Pval","EAF","SNP")

a <- fix[,c(1:5)]
names(a)[3] <- "SNP"

df_GWAS <- left_join(a,gt_df,by="SNP")
df_GWAS <- df_GWAS %>% 
        mutate(Pval=10^(-as.numeric(log10Pval))) %>% 
        rename(effect_allele=ALT,
               other_allele=REF)

write.table(df_GWAS,"GWAS_ukb-b-929_tired.txt",row.names = F,quote = F,sep = "\t")

###
df_GWAS <- fread("GWAS_ukb-b-929_tired.txt",
                 header = T,sep = "\t",nThread=8)
df_GWAS$phynotype <- "Frequency of tiredness in last 2 weeks"
df_GWAS$population <- 449019
head(df_GWAS)

exp_dat1 <- format_data(
        df_GWAS,
        phenotype_col = "phynotype",
        type = "exposure",
        snp_col = "SNP",
        beta_col = "Beta",
        se_col = "SE",
        samplesize_col = "population",
        effect_allele_col = "effect_allele",
        other_allele_col = "other_allele",
        eaf_col="EAF",
        pval_col = "Pval"
)

exp_dat1_sig <- exp_dat1 %>% filter(pval.exposure<5e-8)

ld_df <- exp_dat1_sig[,c("SNP","pval.exposure")]
names(ld_df) <- c("rsid","pval")

exp_clumped <- ieugwasr::ld_clump_local(ld_df,
                                        clump_kb=10000, clump_r2=0.001, 
                                        clump_p=1,
                                        bfile="F:/R/Project/MR-test/1kg.v3/EUR",
                                        plink_bin="C:/Users/Huangbo/Downloads/Compressed/plink_win64_20230116/plink.exe")

exp_clumped_dat1 <- exp_dat1_sig %>% filter(SNP %in% exp_clumped$rsid)

exp_dat1_t <- format_data(
        df_GWAS,
        phenotype_col = "phynotype",
        type = "outcome",
        snp_col = "SNP",
        beta_col = "Beta",
        se_col = "SE",
        samplesize_col = "population",
        effect_allele_col = "effect_allele",
        other_allele_col = "other_allele",
        eaf_col="EAF",
        pval_col = "Pval"
)

# stroke outcome ----------------------------------------------------------

out_dat1 <- fread("data/outcome/giscome.012vs3456.age-gender-5PC.meta1.txt",
                 header = T,sep = "\t",nThread=8)

out_dat1$phynotype <- "mRS 0-2 vs 3-6 at 3 months"
out_dat1 <- out_dat1 %>% mutate(Freq1=1-Freq1)


out_dat1_t <- format_data(
        out_dat1,
        phenotype_col = "phynotype",
        type = "exposure",
        snp_col = "MarkerName",
        beta_col = "Effect",
        se_col = "StdErr",
        effect_allele_col = "Allele1",
        other_allele_col = "Allele2",
        eaf_col="Freq1",
        pval_col = "P-value"
)


out_dat1 <- format_data(
        out_dat1,
        phenotype_col = "phynotype",
        type = "outcome",
        snp_col = "MarkerName",
        beta_col = "Effect",
        se_col = "StdErr",
        effect_allele_col = "Allele1",
        other_allele_col = "Allele2",
        eaf_col="Freq1",
        pval_col = "P-value"
)



out_dat1_t_sig <- out_dat1_t %>% filter(pval.exposure<1e-6)

ld_df <- out_dat1_t_sig[,c("SNP","pval.exposure")]
names(ld_df) <- c("rsid","pval")

exp_clumped <- ieugwasr::ld_clump_local(ld_df,
                                        clump_kb=10000, clump_r2=0.01, 
                                        clump_p=1,
                                        bfile="F:/R/Project/MR-test/1kg.v3/EUR",
                                        plink_bin="C:/Users/Huangbo/Downloads/Compressed/plink_win64_20230116/plink.exe")

out_clumped_dat1_t <- out_dat1_t_sig %>% filter(SNP %in% exp_clumped$rsid)



# 抑郁 失眠 焦虑 ----------------------------------------------------------------

## 抑郁
vcf_data <- vcfR::read.vcfR("data/ieu-b-102.vcf.gz")

meta <- data.frame(vcf_data@meta)
fix <- data.frame(vcf_data@fix) 
gt <- data.frame(vcf_data@gt) 
head(gt)

gt_df <- strsplit(gt[[names(gt)[2]]], split=":")
gt_df <- as.data.frame(do.call(rbind, gt_df))
names(gt_df) <- c("Beta","SE","log10Pval","EAF","SNP")

a <- fix[,c(1:5)]
names(a)[3] <- "SNP"

df_GWAS <- left_join(a,gt_df,by="SNP")
df_GWAS <- df_GWAS %>% 
        mutate(Pval=10^(-as.numeric(log10Pval))) %>% 
        rename(effect_allele=ALT,
               other_allele=REF)

write.table(df_GWAS,"GWAS_ieu-b-102_depression.txt",row.names = F,quote = F,sep = "\t")

###
df_GWAS <- fread("GWAS_ukb-b-102_depression.txt",
                 header = T,sep = "\t",nThread=8)
df_GWAS$phynotype <- "Major depression"
df_GWAS$population <- 500199
head(df_GWAS)

exp_dat01 <- format_data(
        df_GWAS,
        phenotype_col = "phynotype",
        type = "exposure",
        snp_col = "SNP",
        beta_col = "Beta",
        se_col = "SE",
        samplesize_col = "population",
        effect_allele_col = "effect_allele",
        other_allele_col = "other_allele",
        eaf_col="EAF",
        pval_col = "Pval"
)

exp_dat01_sig <- exp_dat01 %>% filter(pval.exposure<5e-8)

###
df_GWAS <- fread("data/MDD2018_ex23andMe.gz",
                 header = T,sep = "\t",nThread=8)
head(df_GWAS)

df_GWAS$phynotype <- "Depression"
df_GWAS$population <- 480359
df_GWAS <- df_GWAS %>% mutate(Beta=log(OR, base = exp(1)))
head(df_GWAS)

exp_dat01 <- format_data(
        df_GWAS,
        phenotype_col = "phynotype",
        type = "exposure",
        snp_col = "SNP",
        beta_col = "Beta",
        se_col = "SE",
        samplesize_col = "population",
        effect_allele_col = "A1",
        other_allele_col = "A2",
        eaf_col="FRQ_A_59851",
        pval_col = "P"
)

exp_dat01_sig <- exp_dat01 %>% filter(pval.exposure<5e-8)


exp_dat01 <- exp_dat01_1
exp_dat01_sig <- exp_dat01_1_sig

## 失眠
vcf_data <- vcfR::read.vcfR("data/ukb-b-3957.vcf.gz")

meta <- data.frame(vcf_data@meta)
fix <- data.frame(vcf_data@fix) 
gt <- data.frame(vcf_data@gt) 
head(gt)

gt_df <- strsplit(gt[[names(gt)[2]]], split=":")
gt_df <- as.data.frame(do.call(rbind, gt_df))
names(gt_df) <- c("Beta","SE","log10Pval","EAF","SNP")

a <- fix[,c(1:5)]
names(a)[3] <- "SNP"

df_GWAS <- left_join(a,gt_df,by="SNP")
df_GWAS <- df_GWAS %>% 
        mutate(Pval=10^(-as.numeric(log10Pval))) %>% 
        rename(effect_allele=ALT,
               other_allele=REF)

write.table(df_GWAS,"GWAS_ukb-b-3957_insomnia.txt",row.names = F,quote = F,sep = "\t")

###
df_GWAS <- fread("GWAS_ukb-b-3957_insomnia.txt",
                 header = T,sep = "\t",nThread=8)
df_GWAS$phynotype <- "Insomnia"
df_GWAS$population <- 462341
head(df_GWAS)

exp_dat02 <- format_data(
        df_GWAS,
        phenotype_col = "phynotype",
        type = "exposure",
        snp_col = "SNP",
        beta_col = "Beta",
        se_col = "SE",
        samplesize_col = "population",
        effect_allele_col = "effect_allele",
        other_allele_col = "other_allele",
        eaf_col="EAF",
        pval_col = "Pval"
)

exp_dat02_sig <- exp_dat02 %>% filter(pval.exposure<5e-8)



## Anxiety
df_GWAS <- fread("data/HillWD_30867560_Anxiety_Tension_Special_Factor.txt.gz",
                 header = T,sep = "\t",nThread=8)

df_GWAS$phynotype <- "Anxiety"
df_GWAS$population <- 270059

exp_dat03 <- format_data(
        df_GWAS,
        phenotype_col = "phynotype",
        type = "exposure",
        snp_col = "SNP",
        beta_col = "Beta",
        se_col = "SE_of_Beta",
        samplesize_col = "population",
        effect_allele_col = "Effect_Allele",
        other_allele_col = "NonEffect_Allele",
        eaf_col = "Effect_Allele_Frequency",
        pval_col = "P"
)

exp_dat03_sig <- exp_dat03 %>% filter(pval.exposure<5e-8)


# CSVD --------------------------------------------------------------------
## ALL
df_GWAS <- fread("data/MTAG_ICH.all_SVS_trait_1.txt",
                 header = T,sep = "\t",nThread=8)

df_GWAS$phynotype <- "ICH or SVS"
df_GWAS$population <- 239313

out_dat01 <- format_data(
        df_GWAS,
        phenotype_col = "phynotype",
        type = "outcome",
        snp_col = "SNP",
        beta_col = "mtag_beta",
        se_col = "mtag_se",
        samplesize_col = "population",
        effect_allele_col = "A1",
        other_allele_col = "A2",
        eaf_col = "freq",
        pval_col = "mtag_pval"
)

## lobar
df_GWAS <- fread("data/MTAG_ICH.lobar_SVS_trait_1.txt",
                 header = T,sep = "\t",nThread=8)

df_GWAS$phynotype <- "Lobar ICH or SVS"
df_GWAS$population <- 238266

out_dat02 <- format_data(
        df_GWAS,
        phenotype_col = "phynotype",
        type = "outcome",
        snp_col = "SNP",
        beta_col = "mtag_beta",
        se_col = "mtag_se",
        samplesize_col = "population",
        effect_allele_col = "A1",
        other_allele_col = "A2",
        eaf_col = "freq",
        pval_col = "mtag_pval"
)

## non-lobar
df_GWAS <- fread("data/MTAG_ICH.nonlobar_SVS_trait_1.txt",
                 header = T,sep = "\t",nThread=8)

df_GWAS$phynotype <- "Non-lobar ICH or SVS"
df_GWAS$population <- 238526

out_dat03 <- format_data(
        df_GWAS,
        phenotype_col = "phynotype",
        type = "outcome",
        snp_col = "SNP",
        beta_col = "mtag_beta",
        se_col = "mtag_se",
        samplesize_col = "population",
        effect_allele_col = "A1",
        other_allele_col = "A2",
        eaf_col = "freq",
        pval_col = "mtag_pval"
)

out_dat00_comb <- rbind(out_dat01,out_dat02) %>% 
        rbind(out_dat03)



# 多变量 ---------------------------------------------------------------------

## 高血压
df_GWAS <- fread("GWAS_ukb-b-12493_hypertension.txt",
                 header = T,sep = "\t",nThread=8)
df_GWAS$phynotype <- "Hypertension"
df_GWAS$population <- 463010
head(df_GWAS)

exp_dat04 <- format_data(
        df_GWAS,
        phenotype_col = "phynotype",
        type = "exposure",
        snp_col = "SNP",
        beta_col = "Beta",
        se_col = "SE",
        samplesize_col = "population",
        effect_allele_col = "effect_allele",
        other_allele_col = "other_allele",
        eaf_col="EAF",
        pval_col = "Pval"
)

exp_dat04_sig <- exp_dat04 %>% filter(pval.exposure<5e-8)

## 糖尿病
df_GWAS <- fread("data/30054458-GCST006867-EFO_0001360-build37.f.tsv.gz",
                 header = T,sep = "\t",nThread=8)
df_GWAS$phynotype <- "Type 2 diabetes"
df_GWAS$population <- 655666
head(df_GWAS)

exp_dat05 <- format_data(
        df_GWAS,
        phenotype_col = "phynotype",
        type = "exposure",
        snp_col = "variant_id",
        beta_col = "beta",
        se_col = "standard_error",
        samplesize_col = "population",
        effect_allele_col = "effect_allele",
        other_allele_col = "other_allele",
        eaf_col="effect_allele_frequency",
        pval_col = "p_value"
)

exp_dat05_sig <- exp_dat05 %>% filter(pval.exposure<5e-8)

## 高血脂
df_GWAS <- fread("data/34906840-GCST90104005-MONDO_0001336-Build37.f.tsv.gz",
                 header = T,sep = "\t",nThread=8)
df_GWAS$phynotype <- "Hyperlipidemia "
df_GWAS$population <- 349222
head(df_GWAS)

exp_dat06 <- format_data(
        df_GWAS,
        phenotype_col = "phynotype",
        type = "exposure",
        snp_col = "variant_id",
        beta_col = "beta",
        se_col = "standard_error",
        samplesize_col = "population",
        effect_allele_col = "effect_allele",
        other_allele_col = "other_allele",
        eaf_col="effect_allele_frequency",
        pval_col = "p_value"
)

exp_dat06_sig <- exp_dat06 %>% filter(pval.exposure<5e-8)


