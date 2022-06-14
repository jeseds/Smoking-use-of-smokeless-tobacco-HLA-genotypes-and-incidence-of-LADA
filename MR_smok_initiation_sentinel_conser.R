######################conservative analyses########################################################
setwd("C:/D/PhD study/MR/results/211128") 
pacman::p_load(MendelianRandomization)
library(readxl)
conser_smok_lada = read_excel("fig2_blank.xlsx")
View(conser_smok_lada)
smok_lada_noany=subset(SmokingInitiation_gscan, noanytrait==1)  #155 SNPs remained
attach(smok_lada_noany) 

mrobject_noany = mr_input(bx = b_smok, bxse = se_smok, by = b_lada, byse = se_lada,snps = SNP)
conser_no_anydm=mr_ivw(mrobject_noany, model="fixed")
conser_no_anydm
conser_smok_lada[conser_smok_lada$method=="Conservative 1","no_snp"]=toString(round(conser_no_anydm$SNPs,digit=0)) 
conser_smok_lada[conser_smok_lada$method=="Conservative 1","or"]=toString(round(exp(conser_no_anydm$Estimate),digit=2))
conser_smok_lada[conser_smok_lada$method=="Conservative 1","lci"]=toString(round(exp(conser_no_anydm$CILower),digit=2)) 
conser_smok_lada[conser_smok_lada$method=="Conservative 1","uci"]=toString(round(exp(conser_no_anydm$CIUpper),digit=2)) 
conser_smok_lada[conser_smok_lada$method=="Conservative 1","p_estimate"]=toString(round(conser_no_anydm$Pvalue,digit=3)) 
conser_smok_lada[conser_smok_lada$method=="Conservative 1","Cochran_Q"]=toString(round(conser_no_anydm$Heter.Stat[1],digit=0))   
conser_smok_lada[conser_smok_lada$method=="Conservative 1","p_heterogeneity"]=toString(round(conser_no_anydm$Heter.Stat[2],digit=3)) 
detach(smok_lada_noany)

smok_lada_noanydm=subset(SmokingInitiation_gscan, noanytrait==1 & noanydm==1)  
attach(smok_lada_noanydm) 

mrobject_noanydm = mr_input(bx = b_smok, bxse = se_smok, by = b_lada, byse = se_lada,snps = SNP)
conser_no_anydmstyle=mr_ivw(mrobject_noanydm, model="fixed")
conser_no_anydmstyle
conser_smok_lada[conser_smok_lada$method=="Conservative 2","no_snp"]=toString(round(conser_no_anydmstyle$SNPs,digit=0)) 
conser_smok_lada[conser_smok_lada$method=="Conservative 2","or"]=toString(round(exp(conser_no_anydmstyle$Estimate),digit=2))
conser_smok_lada[conser_smok_lada$method=="Conservative 2","lci"]=toString(round(exp(conser_no_anydmstyle$CILower),digit=2)) 
conser_smok_lada[conser_smok_lada$method=="Conservative 2","uci"]=toString(round(exp(conser_no_anydmstyle$CIUpper),digit=2)) 
conser_smok_lada[conser_smok_lada$method=="Conservative 2","p_estimate"]=toString(round(conser_no_anydmstyle$Pvalue,digit=3)) 
conser_smok_lada[conser_smok_lada$method=="Conservative 2","Cochran_Q"]=toString(round(conser_no_anydmstyle$Heter.Stat[1],digit=0))   
conser_smok_lada[conser_smok_lada$method=="Conservative 2","p_heterogeneity"]=toString(round(conser_no_anydmstyle$Heter.Stat[2],digit=3)) 

detach(smok_lada_noanydm) 
smok_lada_noanydmalc=subset(SmokingInitiation_gscan, noanytrait==1 & noanydm==1 & noanyalc==1)  
attach(smok_lada_noanydmalc) 

mrobject_noanydmalc = mr_input(bx = b_smok, bxse = se_smok, by = b_lada, byse = se_lada,snps = SNP)
conser_no_anydmstylebw=mr_ivw(mrobject_noanydmalc, model="fixed")
conser_no_anydmstylebw
conser_smok_lada[conser_smok_lada$method=="Conservative 3","no_snp"]=toString(round(conser_no_anydmstylebw$SNPs,digit=0)) 
conser_smok_lada[conser_smok_lada$method=="Conservative 3","or"]=toString(round(exp(conser_no_anydmstylebw$Estimate),digit=2))
conser_smok_lada[conser_smok_lada$method=="Conservative 3","lci"]=toString(round(exp(conser_no_anydmstylebw$CILower),digit=2)) 
conser_smok_lada[conser_smok_lada$method=="Conservative 3","uci"]=toString(round(exp(conser_no_anydmstylebw$CIUpper),digit=2)) 
conser_smok_lada[conser_smok_lada$method=="Conservative 3","p_estimate"]=toString(round(conser_no_anydmstylebw$Pvalue,digit=3)) 
conser_smok_lada[conser_smok_lada$method=="Conservative 3","Cochran_Q"]=toString(round(conser_no_anydmstylebw$Heter.Stat[1],digit=0))   
conser_smok_lada[conser_smok_lada$method=="Conservative 3","p_heterogeneity"]=toString(round(conser_no_anydmstylebw$Heter.Stat[2],digit=3)) 
detach(smok_lada_noanydmalc) 

library(tidyr)
conser_smok_lada= unite(conser_smok_lada,"or_ci",or, sep_str1,lci,sep_str2, uci,sep_str3,sep="", remove=FALSE)

library(openxlsx)
write.xlsx(conser_smok_lada,file="conser_smok_lada.xlsx",overwrite=TRUE)
getwd()