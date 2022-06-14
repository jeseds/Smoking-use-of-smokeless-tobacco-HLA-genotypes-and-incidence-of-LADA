##CVD亚组分析森林图
####
install.packages('forestplot')
library(forestplot)
library(grid)
library(magrittr)
library(checkmate)
#############################################################################
setwd("C:/D/PhD study/MR/results/211128")
tiff(filename='figure1.tiff',  units="in", width=9, height=6, res=500)
figure1 <- read.csv("C:/D/PhD study/MR/results/211128/figure1.csv",header = FALSE)
forestplot(as.matrix(figure1[,1:5]),hrzl_lines =list ("3"=gpar(col="#111111")),mean=figure1$V6, lower=figure1$V7,upper= figure1$V8,
           graph.pos=4,zero=1,graphwidth=unit(50,'mm'),lineheight='auto',xlog = TRUE,is.summary=c(TRUE,TRUE,TRUE,rep(FALSE,5),TRUE,rep(FALSE,5)),
           xticks=(c(0.30, 1.00,2.00,4.00)),col=fpColors(all.elements='black'),
           xlab=expression(bold(paste("OR (95% CI)"))),
           colgap = unit(2,'mm'),
           clip = c(0.1,2.0), 
           txt_gp=fpTxtGp(label = gpar(fontfamily = "serif"), ticks=gpar(cex=1.0), xlab=gpar(cex = 1.0,adj=0.5),title=gpar(cex = 1.2)))
dev.off()
#fontfamily in R: sans=Arial, serif=Times New Roman, mono=Courier, symbol=Standard Symbols L
# clip = c(0.1,2.5): add arrows when 95% CI exceeds the limits

#Only for LADA
setwd("C:/D/PhD study/MR/results/211128")
tiff(filename='figure1.tiff',  units="in", width=9, height=4, res=500)
figure1 <- read.csv("C:/D/PhD study/MR/results/211128/figure1.csv",header = FALSE)
forestplot(as.matrix(figure1[,1:5]),hrzl_lines =list ("3"=gpar(col="#111111")),mean=figure1$V6, lower=figure1$V7,upper= figure1$V8,
           graph.pos=4,zero=1,graphwidth=unit(50,'mm'),lineheight='auto',xlog = TRUE,is.summary=c(TRUE,TRUE,rep(FALSE,5)),
           xticks=(c(0.50,1.00,2.00,3.00)),col=fpColors(all.elements='black'),
           xlab=expression(bold(paste("OR (95% CI)"))),
           colgap = unit(2,'mm'),
           clip = c(0.1,2.5), 
           txt_gp=fpTxtGp(ticks=gpar(cex=1.0), xlab=gpar(cex = 1.0,adj=0.5),title=gpar(cex = 1.2)))
dev.off()
