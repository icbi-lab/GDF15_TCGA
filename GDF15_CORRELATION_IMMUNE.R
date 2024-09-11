library(ggplot2)
library(grid)

# Correlation of GDF15 expression with various immune signature scores for LUAD, LUSC, BLCA_LUM, BLCA_SQUAM
# Heatmap of correlations (Fig 2c)
# Input: LUAD.txt, LUSC.txt, BLCA_LUM.txt, BLCA_SQUAM.txt
# Output: MAT_RHO.txt, MAT_P.txt, COR_HEATMAP.pdf

BR<-c("T cell CD8+","TCELL_DYSF_POS_NES","TCELL_INFLAM_NES","IFNG_NES","INFLAM_NES","IMPRES","CYT","PDCD1","CD274","CTL_IFNG","CD8","CTL","T_EXCL_NES")
LUAD<-read.table("LUAD.txt", header=TRUE, sep="\t", dec = ".", row.names=1, check.names=FALSE)
LUSC<-read.table("LUSC.txt", header=TRUE, sep="\t", dec = ".", row.names=1, check.names=FALSE)
BLCA_LUM<-read.table("BLCA_LUM.txt", header=TRUE, sep="\t", dec = ".", row.names=1, check.names=FALSE)
BLCA_SQUAM<-read.table("BLCA_SQUAM.txt", header=TRUE, sep="\t", dec = ".", row.names=1, check.names=FALSE)
DRHO<-data.frame(LUAD=rep(0.0,13),LUSC=rep(0.0,13),BLCA_LUM=rep(0.0,13),BLCA_SQUAM=rep(0.0,13))
row.names(DRHO)<-BR
DP<-data.frame(LUAD=rep(0.0,13),LUSC=rep(0.0,13),BLCA_LUM=rep(0.0,13),BLCA_SQUAM=rep(0.0,13))
row.names(DP)<-BR
for (i in 1:13) {
  corr1<-cor.test(LUAD$GDF15,LUAD[[BR[i]]],method="spearman")
  DRHO[i,1]<-corr1$estimate
  DP[i,1]<-corr1$p.value
  corr2<-cor.test(LUSC$GDF15,LUSC[[BR[i]]],method="spearman")
  DRHO[i,2]<-corr2$estimate
  DP[i,2]<-corr2$p.value
  corr3<-cor.test(BLCA_LUM$GDF15,BLCA_LUM[[BR[i]]],method="spearman")
  DRHO[i,3]<-corr3$estimate
  DP[i,3]<-corr3$p.value
  corr4<-cor.test(BLCA_SQUAM$GDF15,BLCA_SQUAM[[BR[i]]],method="spearman")
  DRHO[i,4]<-corr4$estimate
  DP[i,4]<-corr4$p.value
}
write.table(DRHO,"MAT_RHO.txt",sep="\t",row.names=TRUE,quote=FALSE,col.names=NA)
write.table(DP,"MAT_P.txt",sep="\t",row.names=TRUE,quote=FALSE,col.names=NA)

myBreaks <- seq(-0.3, 0.3,by=0.01)

pdf("COR_HEATMAP.pdf",  width=8, height=10)
names(DRHO)<-c("LUAD [n=515]","LUSC [n=501]", "BLCA LUM [n=246]", "BLCA SQUAM [n=142]")
heatmap<- pheatmap::pheatmap(DRHO,show_colnames = T, fontsize_row = 15, fontsize_col = 17,cluster_cols = FALSE,cluster_rows = FALSE,cellheight = 40,cellwidth = 40,
                             color=colorRampPalette(c("blue", "white", "red"))(60),breaks = myBreaks, fontsize=15)
print(heatmap)
y1<-1.013
sz<-1
x1<-0.1268
x2<-0.785
y2<-0.545
ks<-0.005
for (i in 1:4) {
  for (j in 1:13) {
    if (DP[j,i]<0.05) {
      sz<-abs(log10(DP[j,i]))
      if (sz>4) {
        sz=4
      }
      grid.circle(x=x1+i*0.07, y=y1-j*0.0558, r=ks*sz,gp=gpar(fill="white",col=NA))
    }
  }
}
print(grid.text("p<0.05",x=x2, y=y2+0.025,just="left",gp=gpar(fontsize=15,col="black")))
print(grid.text("p=0.01",x=x2, y=y2-0.031,just="left",gp=gpar(fontsize=15,col="black")))
print(grid.text("p=0.001",x=x2, y=y2-0.087,just="left",gp=gpar(fontsize=15,col="black")))
print(grid.text("p<0.0001",x=x2, y=y2-0.143,just="left",gp=gpar(fontsize=15,col="black")))
print(grid.text("Rho",x=0.82, y=0.88,just="left",gp=gpar(fontsize=15,col="black")))
grid.circle(x=0.75, y=y2+0.025, r=ks*abs(log10(0.05)),gp=gpar(fill="gray",col=NA))
grid.circle(x=0.75, y=y2-0.031, r=ks*abs(log10(0.01)),gp=gpar(fill="gray",col=NA))
grid.circle(x=0.75, y=y2-0.087, r=ks*abs(log10(0.001)),gp=gpar(fill="gray",col=NA))
grid.circle(x=0.75, y=y2-0.143, r=ks*abs(log10(0.0001)),gp=gpar(fill="gray",col=NA))
dev.off()

