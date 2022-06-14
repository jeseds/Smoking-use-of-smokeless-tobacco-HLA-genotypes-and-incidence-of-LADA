#data analysis
setwd("C:/D/PhD study/MR/results/211128") 
#######################################################
#1.plot
library(haven)
SmokingInitiation_gscan <- read_dta("C:/D/PhD study/MR/data/SmokingInitiation_gscan.dta")
View(SmokingInitiation_gscan)
SmokingInitiation_gscan_t2d=SmokingInitiation_gscan %>% drop_na(p_t2d)   #1 SNP drop


attach(SmokingInitiation_gscan_t2d)
#scatter plot
#1.plot
bx=SmokingInitiation_gscan_t2d[,"b_smok"]
by=SmokingInitiation_gscan_t2d[,"b_t2d"]
bxse=SmokingInitiation_gscan_t2d[,"se_smok"]
byse=SmokingInitiation_gscan_t2d[,"se_t2d"]

tiff(filename='scatter_smok_t2d.tiff',  units="in", width=10, height=5, res=500)
plot(b_smok, b_t2d, xlim=c(min(b_smok-2*se_smok, 0), max(b_smok+2*se_smok, 0)),
  ylim=c(min(b_t2d-2*se_t2d, 0), max(b_t2d+2*se_t2d, 0)), 
  xlab=expression(bold(paste("Genetic association with birth weight"))),
  ylab=expression(bold(paste("Genetic association with T2D"))),)

for (j in 1:length(b_smok)) {
 lines(c(b_smok[j],b_smok[j]), c(b_t2d[j]-1.96*se_t2d[j], b_t2d[j]+1.96*se_t2d[j]))
 lines(c(b_smok[j]-1.96*se_smok[j],b_smok[j]+1.96*se_smok[j]), c(b_t2d[j], b_t2d[j]))
          }
 abline(a=0, b=sum(bx*by*byse^-2)/sum(bx^2*byse^-2),col="red")
dev.off()

library(readxl)
fig1_smok_t2d = read_excel("fig2_blank.xlsx")
View(fig1_smok_t2d)
########################################################
install.packages("pacman")
## Loading required package: pacman
pacman::p_load(MendelianRandomization)
mrobject_t2d = mr_input(bx = b_smok, bxse = se_smok, by = b_t2d, byse = se_t2d,snps = SNP)
ivw_fixed_smok=mr_ivw(mrobject_t2d, model="fixed")
ivw_fixed_smok

ivw_random_smok=mr_ivw(mrobject_t2d, model="random")
ivw_random_smok
fig1_smok_t2d[fig1_smok_t2d$method=="IVW_normal","no_snp"]=toString(round(ivw_random_smok$SNPs,digit=0)) 
fig1_smok_t2d[fig1_smok_t2d$method=="IVW_normal","or"]=toString(round(exp(ivw_random_smok$Estimate),digit=2))
fig1_smok_t2d[fig1_smok_t2d$method=="IVW_normal","lci"]=toString(round(exp(ivw_random_smok$CILower),digit=2)) 
fig1_smok_t2d[fig1_smok_t2d$method=="IVW_normal","uci"]=toString(round(exp(ivw_random_smok$CIUpper),digit=2)) 
fig1_smok_t2d[fig1_smok_t2d$method=="IVW_normal","p_estimate"]=toString(round(ivw_random_smok$Pvalue,digit=3)) 
fig1_smok_t2d[fig1_smok_t2d$method=="IVW_normal","Cochran_Q"]=toString(round(ivw_random_smok$Heter.Stat[1],digit=0))   
fig1_smok_t2d[fig1_smok_t2d$method=="IVW_normal","p_heterogeneity"]=toString(round(ivw_random_smok$Heter.Stat[2],digit=3)) 

ivw_rob_smok=mr_ivw(mrobject_t2d, robust=TRUE)
ivw_rob_smok
fig1_smok_t2d[fig1_smok_t2d$method=="Robust IVW","no_snp"]=toString(round(ivw_rob_smok$SNPs,digit=0)) 
fig1_smok_t2d[fig1_smok_t2d$method=="Robust IVW","or"]=toString(round(exp(ivw_rob_smok$Estimate),digit=2))
fig1_smok_t2d[fig1_smok_t2d$method=="Robust IVW","lci"]=toString(round(exp(ivw_rob_smok$CILower),digit=2)) 
fig1_smok_t2d[fig1_smok_t2d$method=="Robust IVW","uci"]=toString(round(exp(ivw_rob_smok$CIUpper),digit=2)) 
fig1_smok_t2d[fig1_smok_t2d$method=="Robust IVW","p_estimate"]=toString(round(ivw_rob_smok$Pvalue,digit=3)) 
fig1_smok_t2d[fig1_smok_t2d$method=="Robust IVW","Cochran_Q"]=toString(round(ivw_rob_smok$Heter.Stat[1],digit=0))   
fig1_smok_t2d[fig1_smok_t2d$method=="Robust IVW","p_heterogeneity"]=toString(round(ivw_rob_smok$Heter.Stat[2],digit=3)) 


ivw_penal_smok=mr_ivw(mrobject_t2d, robust=FALSE, penalized=TRUE)
ivw_penal_smok


ivw_rob_penal_smok=mr_ivw(mrobject_t2d, robust=TRUE, penalized=TRUE)
ivw_rob_penal_smok

#median-based model
simple_median_smok=mr_median(mrobject_t2d, weighting="simple")
simple_median_smok

weighted_median_smok=mr_median(mrobject_t2d, weighting="weighted")
weighted_median_smok
fig1_smok_t2d[fig1_smok_t2d$method=="Weighted median","no_snp"]=toString(round(weighted_median_smok$SNPs,digit=0)) 
fig1_smok_t2d[fig1_smok_t2d$method=="Weighted median","or"]=toString(round(exp(weighted_median_smok$Estimate),digit=2))
fig1_smok_t2d[fig1_smok_t2d$method=="Weighted median","lci"]=toString(round(exp(weighted_median_smok$CILower),digit=2)) 
fig1_smok_t2d[fig1_smok_t2d$method=="Weighted median","uci"]=toString(round(exp(weighted_median_smok$CIUpper),digit=2)) 
fig1_smok_t2d[fig1_smok_t2d$method=="Weighted median","p_estimate"]=toString(round(weighted_median_smok$Pvalue,digit=3)) 

penal_median_smok=mr_median(mrobject_t2d, weighting="penalized")
penal_median_smok

# fit some egger models
egger_smok=mr_egger(mrobject_t2d,robust = FALSE, penalized = FALSE)
egger_smok
fig1_smok_t2d[fig1_smok_t2d$method=="MR-Egger","no_snp"]=toString(round(egger_smok$SNPs,digit=0)) 
fig1_smok_t2d[fig1_smok_t2d$method=="MR-Egger","or"]=toString(round(exp(egger_smok$Estimate),digit=2))
fig1_smok_t2d[fig1_smok_t2d$method=="MR-Egger","lci"]=toString(round(exp(egger_smok$CILower.Est),digit=2)) 
fig1_smok_t2d[fig1_smok_t2d$method=="MR-Egger","uci"]=toString(round(exp(egger_smok$CIUpper.Est),digit=2)) 
fig1_smok_t2d[fig1_smok_t2d$method=="MR-Egger","p_estimate"]=toString(round(egger_smok$Pvalue.Est,digit=3)) 
fig1_smok_t2d[fig1_smok_t2d$method=="MR-Egger","Cochran_Q"]=toString(round(egger_smok$Heter.Stat[1],digit=0))   
fig1_smok_t2d[fig1_smok_t2d$method=="MR-Egger","p_heterogeneity"]=toString(round(egger_smok$Heter.Stat[2],digit=3)) 
fig1_smok_t2d[fig1_smok_t2d$method=="MR-Egger","intercept"]=toString(round(egger_smok$Intercept,digit=3)) 
fig1_smok_t2d[fig1_smok_t2d$method=="MR-Egger","p_pleiotropy"]=toString(round(egger_smok$Pvalue.Int,digit=3)) 

#MR-Egger indicated directional pleiortopy for the 49 SNPs:  Intercept=-0.059 P_pleiotropy=0.023 I^2_GX statistic: 71.5%

egger_robust_smok=mr_egger(mrobject_t2d,robust = TRUE, penalized = FALSE)
egger_robust_smok


egger_penalized__smok=mr_egger(mrobject_t2d,robust = FALSE, penalized = TRUE)
egger_penalized__smok

egger_rob_penal_smok=mr_egger(mrobject_t2d,robust = TRUE, penalized = TRUE)
egger_rob_penal_smok


#MR presso: among the 49 SNPs, global test of outliers: P=0.007, rs72851023 was the outlier and P for distorion test=0.398 
library(MRPRESSO)

#The input data needs to be stored in a dat.frame object like the following.
beta_y = SmokingInitiation_gscan_t2d$b_t2d
se_y = SmokingInitiation_gscan_t2d$se_t2d
beta_x = SmokingInitiation_gscan_t2d$b_smok
se_x = SmokingInitiation_gscan_t2d$se_smok
smok_input = data.frame(betaY=beta_y, seY=se_y, betaX=beta_x, seX=se_x, row.names = NULL)

mrpresso_out_smok = mr_presso(BetaOutcome = "betaY", BetaExposure = c("betaX"), SdOutcome = "seY", SdExposure = c("seX"), OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = smok_input, NbDistribution = 10000,  SignifThreshold = 0.05)

 # global test: detection of pleiotropy (MR-PRESSO global test)
  mrpresso_out_smok$`MR-PRESSO results`$`Global Test`$RSSobs
  mrpresso_out_smok$`MR-PRESSO results`$`Global Test`$Pvalue
  #identify  outliers
  out_indices = mrpresso_out_smok$`MR-PRESSO results`$`Distortion Test`$`Outliers Indices`
SNP[out_indices]
#testing of significant distortion in the causal estimate before and after MR-PRESSO correction (MR-PRESSO distortion test).
mrpresso_out_smok$`MR-PRESSO results`$`Distortion Test`$`Distortion Coefficient`
mrpresso_out_smok$`MR-PRESSO results`$`Distortion Test`$Pvalue
mrpresso_out_smok

fig1_smok_t2d[fig1_smok_t2d$method=="MR-PRESSO","or"]=toString(round(exp(mrpresso_out_smok$`Main MR results`$'Causal Estimate'[2]),digit=2))
fig1_smok_t2d[fig1_smok_t2d$method=="MR-PRESSO","lci"]=toString(round(exp(mrpresso_out_smok$`Main MR results`$'Causal Estimate'[2]-1.96*mrpresso_out_smok$`Main MR results`$'Sd'[2]),digit=2)) 
fig1_smok_t2d[fig1_smok_t2d$method=="MR-PRESSO","uci"]=toString(round(exp(mrpresso_out_smok$`Main MR results`$'Causal Estimate'[2]+1.96*mrpresso_out_smok$`Main MR results`$'Sd'[2]),digit=2)) 
fig1_smok_t2d[fig1_smok_t2d$method=="MR-PRESSO","p_estimate"]=toString(round(mrpresso_out_smok$`Main MR results`$'P-value'[2],digit=3)) 
fig1_smok_t2d[fig1_smok_t2d$method=="MR-PRESSO","p_pleiotropy"]=toString(round(mrpresso_out_smok$`MR-PRESSO results`$`Global Test`$Pvalue,digit=3)) 
fig1_smok_t2d[fig1_smok_t2d$method=="MR-PRESSO","outliers"]=SNP[out_indices]
fig1_smok_t2d[fig1_smok_t2d$method=="MR-PRESSO","p_distort"]=toString(round(mrpresso_out_smok$`MR-PRESSO results`$`Distortion Test`$Pvalue,digit=3)) 

library(tidyr)
fig1_smok_t2d= unite(fig1_smok_t2d,"or_ci",or, sep_str1,lci,sep_str2, uci,sep_str3,sep="", remove=FALSE)

install.packages("xlsx")
library(openxlsx)
write.xlsx(fig1_smok_t2d,file="fig1_smok_t2d.xlsx",overwrite=TRUE)
getwd()
detach(SmokingInitiation_gscan_t2d)

#leave-one-out analysis
attach(SmokingInitiation_gscan_t2d)
mrobject_t2d = mr_input(bx = b_smok, bxse = se_smok, by = b_t2d, byse = se_t2d,snps = SNP)
tiff(filename='loo_smok_lada.tiff',  units="in", width=5, height=25, res=500)
leave_one_out_birwei=mr_loo(mrobject_t2d, alpha = 0.05)
dev.off()
detach(SmokingInitiation_gscan_t2d)

attach(SmokingInitiation_gscan_t2d)
mrobject_t2d = mr_input(bx = b_smok, bxse = se_smok, by = b_t2d, byse = se_t2d,snps = SNP)
tiff(filename='funnel_t2d.tiff',  units="in", width=9, height=6, res=500)
funnel=mr_funnel(mrobject_t2d, CI = FALSE)
dev.off()
detach(SmokingInitiation_gscan_t2d)
###########conservative analysis by excluding SNPs associated with any trait at p<5e-8######
library(dplyr)
library(tidyr)
SmokingInitiation_gscan_t2d_noany=SmokingInitiation_gscan_t2d%>% filter(noanytrait==1)%>% select(SNP,e_smok,r_smok,b_smok,se_smok, e_t2d,r_t2d,b_t2d,se_t2d,p_smok,p_t2d,noanytrait) 

attach(SmokingInitiation_gscan_t2d_noany)
pacman::p_load(MendelianRandomization)
mrobject_t2d = mr_input(bx = b_smok, bxse = se_smok, by = b_t2d, byse = se_t2d,snps = SNP)
ivw_fixed_smok=mr_ivw(mrobject_t2d, model="fixed")
ivw_fixed_smok

ivw_random_smok=mr_ivw(mrobject_t2d, model="random")
ivw_random_smok
detach(SmokingInitiation_gscan_t2d_noany)  #1.14 (0.99-1.32), p for estimate=0.060, p for heterogeneity<0.0001
