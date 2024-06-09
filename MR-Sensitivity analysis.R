


##
## F stat
i="Tiredness"

dat$EAF2 <- 1-dat$eaf.exposure
dat$MAF <- pmin(dat$eaf.exposure,dat$EAF2)
PVEfx <- function(BETA,MAF,SE,N){
        pve <- (2*(BETA^2)*MAF*(1-MAF))/((2*(BETA^2)*MAF*(1-MAF))+((SE^2)*2*N*MAF*(1-MAF)))
        return(pve)
}

dat$PVE <- mapply(PVEfx,BETA=dat$beta.exposure,
                  MAF=dat$MAF,
                  SE=dat$se.exposure,
                  N=dat$samplesize.exposure)
dat$FSTAT <- ((dat$samplesize.exposure-1-1)/1)*(dat$PVE/(1-dat$PVE))

write.csv(dat,paste0("mr_report/dat_",i,".csv"))

## egger
het <- mr_heterogeneity(dat)

write.csv(het,paste0("mr_report/het_",i,".txt"))

I2 <- (het[2,]$Q-het[2,]$Q_df)/het[2,]$Q

write.csv(I2,paste0("mr_report/I2_",i,".txt"))

## 散点图
p <- mr_scatter_plot(res, dat)[[1]]+
        theme_classic()+theme(legend.position="top")

ggsave(paste0("mr_report/plot-scatter_",i,".pdf"), plot = p, device = NULL, 
       width = 10, height = 10, units =c("cm"))

## 森林图
res_single <- mr_singlesnp(dat)
p <- mr_forest_plot(res_single)[[1]]+
        theme_classic()+theme(legend.position="none")

ggsave(paste0("mr_report/plot-forest_",i,".pdf"), plot = p, device = NULL, 
       width = 10, height = 14, units =c("cm"))

## 水平多效性
pleiotropy <- mr_pleiotropy_test(dat)

write.csv(pleiotropy,paste0("mr_report/pleiotropy_",i,".csv"))

## funnel plot
p <- mr_funnel_plot(res_single)[[1]]+
        theme_classic()+theme(legend.position="top")

ggsave(paste0("mr_report/plot-funnel_",i,".pdf"), plot = p, device = NULL, 
       width = 10, height = 10, units =c("cm"))

##
leaveone <- mr_leaveoneout(dat)
write.csv(leaveone,paste0("mr_report/leaveone_",i,".txt"))

p <- mr_leaveoneout_plot(leaveoneout_results = mr_leaveoneout(dat))[[1]]+
        theme_classic()+theme(legend.position="none")

ggsave(paste0("mr_report/plot-leaveoneout_",i,".pdf"), plot = p, device = NULL, 
       width = 10, height = 14, units =c("cm"))
