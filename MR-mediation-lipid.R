



## 
ao <- available_outcomes()
ao_clean <- ao[grepl("met-d",ao$id),]

## MR tired to ALL
cl <- makeCluster(8)
parallel::clusterExport(cl= cl,varlist = c("exp_clumped_dat1"))
parallel::clusterEvalQ(cl= cl,library(tidyverse))
parallel::clusterEvalQ(cl= cl,library(data.table))
parallel::clusterEvalQ(cl= cl,library(TwoSampleMR))

res_list = pblapply(cl = cl,
                    X = ao_clean$id,
                    FUN = function(i){
                            tryCatch({
                                    
                                    chd_out_dat <- extract_outcome_data(snps = exp_clumped_dat1$SNP, outcomes = i)
                                    dat <- harmonise_data(exp_clumped_dat1, chd_out_dat)
                                    res <- mr(dat)
                                    
                                    return(res)
                                    
                            }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
                    })

res_list_com = dplyr::bind_rows(res_list)

stopCluster(cl)

write.csv(res_list_com,"res_tired2Met.csv")

res_ivw <- res_list_com %>% filter(method %in% c("Inverse variance weighted"))
res_ivw_sig1 <- res_ivw %>% filter(pval<0.05)


## MR ALL to strock outcome
cl <- makeCluster(8)
parallel::clusterExport(cl= cl,varlist = c("out_dat1"))
parallel::clusterEvalQ(cl= cl,library(tidyverse))
parallel::clusterEvalQ(cl= cl,library(data.table))
parallel::clusterEvalQ(cl= cl,library(TwoSampleMR))

res_list = pblapply(cl = cl,
                    X = res_ivw_sig1$id.outcome,
                    FUN = function(i){
                            tryCatch({
                                    
            bmi_exp_dat <- extract_instruments(outcomes = i)
            dat <- harmonise_data(bmi_exp_dat, out_dat1)
            res <- mr(dat)
                                    
            return(res)
                                    
                            }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
                    })

res_list_com = dplyr::bind_rows(res_list)

stopCluster(cl)

write.csv(res_list_com,"res_Met2strock.csv")

res_ivw <- res_list_com %>% filter(method %in% c("Inverse variance weighted"))
res_ivw_sig2 <- res_ivw %>% filter(pval<0.05)



res_list <- read.csv("res_Met2strock.csv")

res_ivw <- res_list %>% filter(method %in% c("Inverse variance weighted"))
res_ivw_sig <- res_ivw %>% filter(pval<0.05&nsnp>2)
res_ivw_sig <- res_ivw_sig %>% filter(id.exposure%in%ids)


write.csv(res_ivw_sig,"MRres_Met2Strocksig.csv")

# MVMR --------------------------------------------------------------------
cl <- makeCluster(8)
parallel::clusterExport(cl= cl,varlist = c("out_dat1"))
parallel::clusterEvalQ(cl= cl,library(tidyverse))
parallel::clusterEvalQ(cl= cl,library(data.table))
parallel::clusterEvalQ(cl= cl,library(TwoSampleMR))

res_list = pblapply(cl = cl,
                    X = res_ivw_sig2$id.exposure,
                    FUN = function(i){
                            tryCatch({
                                    
                                    id_exposure <- c("ukb-b-929",i) 
                                    exposure_dat <- mv_extract_exposures(id_exposure)
                                    mvdat <- mv_harmonise_data(exposure_dat,  out_dat1) 
                                    res <- mv_multiple(mvdat) 
                                    
                                    return(res[["result"]])
                                    
                            }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
                    })

res_list_com = dplyr::bind_rows(res_list)

stopCluster(cl)

write.csv(res_list_com,"MVMRres_tiredMetsig2Strockoutcom.csv")


# med ---------------------------------------------------------------------

res <- read.csv("res_tired2Met.csv")

res_ivw <- res %>% filter(method %in% c("Inverse variance weighted"))
res_ivw_sig <- res_ivw %>% filter(pval<0.05)
res_ivw_sig <- res_ivw_sig %>% filter(id.outcome%in%ids)

write.csv(res_ivw_sig,"MRres_tiredMetsig.csv")


res_mvmr <- read.csv("MVMRres_tiredMetsig2Strockoutcom.csv") %>% 
        filter(!id.exposure=="ukb-b-929")

ids <- intersect(res_ivw_sig$id.outcome,res_mvmr$id.exposure)

res_list <- list()
for(i in ids){
        tryCatch({
                res1 <- res_ivw_sig %>% filter(id.outcome==i)
                res2 <- res_mvmr %>% filter(id.exposure==i)
                
                b.a <- res1$b
                se.a <- res1$se
                
                b.b <- res2$b
                se.b <- res2$se
                
                b=b.a*b.b
                se=((b.a*se.b)^2+(b.b*se.a)^2+se.a^2*se.b^2)^0.5
                
                ef_m <- c(sprintf("%0.2f",b),sprintf("%0.2f",b-1.96*se),sprintf("%0.2f",b+1.96*se))       
                res_list[[i]] <- ef_m
                
        }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}

res <- t(as.data.frame(res_list, check.names = FALSE))
res <- as.data.frame(res)
res$var <- row.names(res)

write.csv(res,"med_tiredMet2Strockoutcom.csv")


