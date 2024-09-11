library(ggplot2)

# Mean purity corrected cell fractions from quanTIseq for LUAD, LUSC, BLCA LUM, BLCA SQUAM for 3 GDF15 expression classes
# Stacked barplot (Fig 2d)
# Input: LUAD.txt, LUSC.txt, BLCA_LUM.txt, BLCA_SQUAM.txt, PANC_QUANTISEQ_PURITY_CORRECTED.txt
# Output: CELL_FRACTIONS.pdf, CELL_FRACTIONS_PURITY_CORRECTED.txt

terc<-function(x) {
  tertile_limits <- quantile(x, seq(0, 1, 1/3), na.rm = TRUE)
  x_terc <- cut(x, tertile_limits, c('low', 'mid', 'high'), include.lowest = TRUE)
  return (x_terc)
}
  
LU<-read.table("LUAD.txt", header=TRUE, sep="\t", dec = ".", check.names=FALSE,row.names=1)
LS<-read.table("LUSC.txt", header=TRUE, sep="\t", dec = ".", check.names=FALSE,row.names=1)
BL<-read.table("BLCA_LUM.txt", header=TRUE, sep="\t", dec = ".", check.names=FALSE,row.names=1)
BS<-read.table("BLCA_SQUAM.txt", header=TRUE, sep="\t", dec = ".", check.names=FALSE,row.names=1)
IM<-read.table("PANC_QUANTISEQ_PURITY_CORRECTED.txt", header=TRUE, sep="\t", dec = ".", check.names=FALSE)

CELLS<-c("B cell","Macrophage M1","Macrophage M2","Monocyte","Neutrophil","NK cell","T cell CD4+ (non-regulatory)","T cell CD8+","T cell regulatory (Tregs)","Dendritic cell","ICFractions_sum")

PROJ<-"LUAD"
NLU<- LU[,CELLS]
for (i in 1:nrow(NLU)) {
  ID<-row.names(NLU)[i]
  for (C in CELLS) {
    if (length(which(IM$cell_type==C & IM$patient==ID & IM$project==PROJ))>0) {
      NLU[i,C]<-IM$fraction_value_tpurity_corrected[which(IM$cell_type==C & IM$patient==ID & IM$project==PROJ)[1]]
    } else {
      NLU[i,C]<-NA
    }
  }
}
NLU$TERT<-terc(LU$GDF15)
LUL<-NLU[which(NLU$TERT=="low"),CELLS]
LUM<-NLU[which(NLU$TERT=="mid"),CELLS]
LUH<-NLU[which(NLU$TERT=="high"),CELLS]
LU1<-apply(LUL,2, mean,na.rm=TRUE)
LU2<-apply(LUM,2, mean,na.rm=TRUE)
LU3<-apply(LUH,2, mean,na.rm=TRUE)
LUDF<-data.frame(class=c(rep("low",length(LU1)),rep("mid",length(LU2)),rep("high",length(LU3))),cells=c(names(LU1),names(LU2),names(LU3)),fraction=c(LU1,LU2,LU3))

PROJ<-"LUSC"
NLS<- LS[,CELLS]
for (i in 1:nrow(NLS)) {
  ID<-row.names(NLS)[i]
  for (C in CELLS) {
    if (length(which(IM$cell_type==C & IM$patient==ID & IM$project==PROJ))>0) {
      NLS[i,C]<-IM$fraction_value_tpurity_corrected[which(IM$cell_type==C & IM$patient==ID & IM$project==PROJ)[1]]
    } else {
      NLS[i,C]<-NA
    }
  }
}
NLS$TERT<-terc(LS$GDF15)
LSL<-NLS[which(NLS$TERT=="low"),CELLS]
LSM<-NLS[which(NLS$TERT=="mid"),CELLS]
LSH<-NLS[which(NLS$TERT=="high"),CELLS]
LS1<-apply(LSL,2, mean,na.rm=TRUE)
LS2<-apply(LSM,2, mean,na.rm=TRUE)
LS3<-apply(LSH,2, mean,na.rm=TRUE)
LSDF<-data.frame(class=c(rep("low",length(LS1)),rep("mid",length(LS2)),rep("high",length(LS3))),cells=c(names(LS1),names(LS2),names(LS3)),fraction=c(LS1,LS2,LS3))

PROJ<-"BLCA"
NBL<- BL[,CELLS]
for (i in 1:nrow(NBL)) {
  ID<-row.names(NBL)[i]
  for (C in CELLS) {
    if (length(which(IM$cell_type==C & IM$patient==ID & IM$project==PROJ))>0) {
      NBL[i,C]<-IM$fraction_value_tpurity_corrected[which(IM$cell_type==C & IM$patient==ID & IM$project==PROJ)[1]]
    } else {
      NBL[i,C]<-NA
    }
  }
}
NBL$TERT<-terc(BL$GDF15)
BLL<-NBL[which(NBL$TERT=="low"),CELLS]
BLM<-NBL[which(NBL$TERT=="mid"),CELLS]
BLH<-NBL[which(NBL$TERT=="high"),CELLS]
BL1<-apply(BLL,2, mean,na.rm=TRUE)
BL2<-apply(BLM,2, mean,na.rm=TRUE)
BL3<-apply(BLH,2, mean,na.rm=TRUE)
BLDF<-data.frame(class=c(rep("low",length(BL1)),rep("mid",length(BL2)),rep("high",length(BL3))),cells=c(names(BL1),names(BL2),names(BL3)),fraction=c(BL1,BL2,BL3))

PROJ<-"BLCA"
NBS<- BS[,CELLS]
for (i in 1:nrow(NBS)) {
  ID<-row.names(NBS)[i]
  for (C in CELLS) {
    if (length(which(IM$cell_type==C & IM$patient==ID & IM$project==PROJ))>0) {
      NBS[i,C]<-IM$fraction_value_tpurity_corrected[which(IM$cell_type==C & IM$patient==ID & IM$project==PROJ)[1]]
    } else {
      NBS[i,C]<-NA
    }
  }
}
NBS$TERT<-terc(BS$GDF15)
BSL<-NBS[which(NBS$TERT=="low"),CELLS]
BSM<-NBS[which(NBS$TERT=="mid"),CELLS]
BSH<-NBS[which(NBS$TERT=="high"),CELLS]
BS1<-apply(BSL,2, mean,na.rm=TRUE)
BS2<-apply(BSM,2, mean,na.rm=TRUE)
BS3<-apply(BSH,2, mean,na.rm=TRUE)
BSDF<-data.frame(class=c(rep("low",length(BS1)),rep("mid",length(BS2)),rep("high",length(BS3))),cells=c(names(BS1),names(BS2),names(BS3)),fraction=c(BS1,BS2,BS3))

DF5<-rbind(LUDF,LSDF,BLDF,BSDF)
V5<-c(rep("LUAD",nrow(LUDF)),rep("LUSC",nrow(LSDF)),rep("BLCA LUM",nrow(BLDF)),rep("BLCA SQUAM",nrow(BSDF)) )
DF6<-data.frame(DF5,project=V5)
DF6$project<-factor(DF6$project,levels=c("LUAD","LUSC","BLCA LUM","BLCA SQUAM"))
DF6$class<-factor(DF6$class,levels=c("low","mid","high"))

for (proj in c("LUAD","LUSC","BLCA LUM","BLCA SQUAM")) {
    for (class in c("low","mid","high")){
      DF7<-DF6[which(DF6$class==class & DF6$project==proj),]
      ind<-which(DF6$class==class & DF6$project==proj & DF6$cells=="ICFractions_sum")
      DF6$fraction[ind]<-(1-sum(DF7$fraction[1:10]))
    }
}
VC<-rep(CELLS3,nrow(DF6)/length(CELLS3))
DF6$cells<-VC
DF6$cells<-factor(DF6$cells,levels=c("Other cell","B cell","Dendritic cell","Macrophage M1","Macrophage M2","Monocyte","Neutrophil","NK cell","T cell CD4+ (non-regulatory)","T cell CD8+","T cell regulatory (Tregs)"))
ccol<-c("#FF0000","#FFA602","#D6C468","#FFF81F","#FF3881","#C10092","#FF85E8","#B8E6FF","#3AD9FF","#0084FF","#0000FF")

g2<-ggplot(DF6, aes(x=class, y=fraction, fill = cells)) +
    geom_col(position = "stack") +
    scale_fill_manual(values=ccol) +
    geom_bar(position="fill", color="black",stat="identity", width=0.96) +
    facet_grid(~project, scales = "free_x", space = "free_x") +
    theme(panel.spacing = unit(0.5, "cm"),strip.text = element_text(size = 10, family = "arial")) +
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    labs(y= "", x = "GDF15 expression class")


pdf("CELL_FRACTIONS.pdf")
print (g2)
dev.off()
write.table(DF6,"CELL_FRACTIONS_PURITY_CORRECTED.txt",sep="\t",row.names=FALSE,quote=FALSE)









