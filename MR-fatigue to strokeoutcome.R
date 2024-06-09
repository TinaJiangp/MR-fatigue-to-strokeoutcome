



# MR to ----------------------------------------------------------------------
## MR

dat <- harmonise_data(exp_clumped_dat1, out_dat1)
res1 <- mr(dat)
res1


dat <- harmonise_data(out_clumped_dat1_t, exp_dat1_t)
res2 <- mr(dat)
res2

res_com <- rbind(res1,res2)

## result
res_or <- generate_odds_ratios(res_com)
res_or <- res_or %>% 
        mutate(b_2=sprintf("%0.2f",b),
               se_2=sprintf("%0.2f",se),
               pval_2=sprintf("%0.3f",pval)) %>% 
        mutate(OR95CI=str_c(sprintf("%0.2f",or)," (",sprintf("%0.2f",or_lci95),", ",sprintf("%0.2f",or_uci95),")"))

res_or <- res_or %>% arrange(exposure,method)

library(forestplot)

df_plot <- res_or
df_plot$pval_2 <- as.character(df_plot$pval_2)

pdf(file='13forestplot_tired2outcome.pdf',width=48/2.54,height=10/2.54)
df_plot |>
        forestplot(labeltext = c(
                exposure,
                outcome,
                method,
                nsnp,b_2,se_2,
                OR95CI, pval_2),
                #clip = c(0.1, 2.5),
                graph.pos=5, #为Pvalue箱线图所在的位置
                mean=or,
                lower=or_lci95, 
                upper=or_uci95,
                xlog = TRUE,
                cex=1, lineheight = "auto",
                colgap=unit(8,"mm"),
                #箱子大小，线的宽度
                lwd.ci=1, boxsize=0.3) |>
        fp_add_lines(h_2 = gpar(lty = 2)) |> 
        fp_set_style(box = "#288DEF",
                     line = "#667DDC",
                     summary = "#288DEF",
                     align = "rrrrrr",
                     hrz_lines = "#999999") |> 
        fp_add_header(
                exposure = c("Exposure"),
                outcome = c("Outcome"),
                method = c("Method"),
                nsnp = c("Number of SNPs"),
                b_2 = c("Beta"),
                se_2 = c("Standard error"),
                OR95CI = c("OR (95%CI)"),
                pval_2 = c("P value")) |>
        fp_set_zebra_style("#EFEFEF")
dev.off()


# MR to 脑小血管 --------------------------------------------------------------

ids <- unique(out_dat00_comb$outcome)

## MR

res_list <- list()
for(i in ids){
        
        print(i)
        out_dat <- out_dat00_comb %>% filter(outcome==i)
        dat <- harmonise_data(exp_clumped_dat1, out_dat)
        res <- mr(dat)
        res_list[[i]] <- res
}

res_list_bind <- bind_rows(res_list)

## result
res_or <- generate_odds_ratios(res_list_bind)
res_or <- res_or %>% 
        mutate(b_2=sprintf("%0.2f",b),
               se_2=sprintf("%0.2f",se),
               pval_2=sprintf("%0.3f",pval)) %>% 
        mutate(OR95CI=str_c(sprintf("%0.2f",or)," (",sprintf("%0.2f",or_lci95),", ",sprintf("%0.2f",or_uci95),")"))

res_or1 <- res_or %>% arrange(outcome,method)

library(forestplot)

df_plot <- res_or
df_plot$pval_2 <- as.character(df_plot$pval_2)

pdf(file='13forestplot_tired2CSVD.pdf',width=35/2.54,height=15/2.54)
df_plot |>
        forestplot(labeltext = c(
                outcome,
                method,
                nsnp,b_2,se_2,
                OR95CI, pval_2),
                #clip = c(0.1, 2.5),
                graph.pos=6, #为Pvalue箱线图所在的位置
                mean=or,
                lower=or_lci95, 
                upper=or_uci95,
                xlog = TRUE,
                cex=1, lineheight = "auto",
                colgap=unit(8,"mm"),
                #箱子大小，线的宽度
                lwd.ci=1, boxsize=0.3) |>
        fp_add_lines(h_2 = gpar(lty = 2)) |> 
        fp_set_style(box = "#288DEF",
                     line = "#667DDC",
                     summary = "#288DEF",
                     align = "lrrrrrr",
                     hrz_lines = "#999999") |> 
        fp_add_header(
                #exposure = c("Exposure"),
                outcome = c("Outcome"),
                method = c("Method"),
                nsnp = c("Number of SNPs"),
                b_2 = c("Beta"),
                se_2 = c("Standard error"),
                OR95CI = c("OR (95%CI)"),
                pval_2 = c("P value")) |>
        fp_set_zebra_style("#EFEFEF")
dev.off()

# MR to 脑小血管 marker--------------------------------------------------------------

ids <- unique(out_dat_comb$outcome)

## MR

res_list <- list()
for(i in ids){
        
        print(i)
        out_dat <- out_dat_comb %>% filter(outcome==i)
        dat <- harmonise_data(exp_clumped_dat1, out_dat)
        res <- mr(dat)
        res_list[[i]] <- res
}

res_list_bind <- bind_rows(res_list)

## result
res_or <- generate_odds_ratios(res_list_bind)
res_or <- res_or %>% 
        mutate(b_2=sprintf("%0.2f",b),
               se_2=sprintf("%0.2f",se),
               pval_2=sprintf("%0.3f",pval)) %>% 
        mutate(OR95CI=str_c(sprintf("%0.2f",or)," (",sprintf("%0.2f",or_lci95),", ",sprintf("%0.2f",or_uci95),")"))

res_or2 <- res_or %>% arrange(outcome,method)

res_or_comb <- rbind(res_or1,res_or2) %>% filter(method == "Inverse variance weighted")

#res_or_comb <- res_or_comb %>% arrange(outcome,method)

df_plot <- res_or_comb
df_plot$pval_2 <- as.character(df_plot$pval_2)

pdf(file='13forestplot_tired2CSVD_ALL.pdf',width=35/2.54,height=10/2.54)
df_plot |>
        forestplot(labeltext = c(
                outcome,
                #method,
                nsnp,b_2,se_2,
                OR95CI, pval_2),
                #clip = c(0.1, 2.5),
                graph.pos=5, #为Pvalue箱线图所在的位置
                mean=or,
                lower=or_lci95, 
                upper=or_uci95,
                xlog = TRUE,
                cex=1, lineheight = "auto",
                colgap=unit(8,"mm"),
                #箱子大小，线的宽度
                lwd.ci=1, boxsize=0.3) |>
        fp_add_lines(h_2 = gpar(lty = 2)) |> 
        fp_set_style(box = "#288DEF",
                     line = "#667DDC",
                     summary = "#288DEF",
                     align = "lrrrrr",
                     hrz_lines = "#999999") |> 
        fp_add_header(
                #exposure = c("Exposure"),
                outcome = c("Outcome"),
                #method = c("Method"),
                nsnp = c("Number of SNPs"),
                b_2 = c("Beta"),
                se_2 = c("Standard error"),
                OR95CI = c("OR (95%CI)"),
                pval_2 = c("P value")) |>
        fp_set_zebra_style("#EFEFEF")
dev.off()





