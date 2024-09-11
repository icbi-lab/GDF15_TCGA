library(multtest)

# Association of GDF15 expression with 13 immune signatures in 30 cancer indications
# Matrix indicating significant associations  (Fig 2a)
# Input:<CANCER>_primary_TPM.tsv, PANC_IMMUNE_SIGNATURES_SSGSEA.txt, PANC_IMPRES.txt, PANC_QUANTISEQ.tsv
# Output: GDF15_RHO_INDICATIONS.txt, GDF15_PVAL_INDICATIONS.txt, GDF15_FDR_INDICATIONS.txt, GDF_RES_INDICATIONS.txt

CANCER<-c("CRC","BLCA","BRCA","PAAD","LUAD","KIRC","CESC","TGCT","ESCA","STAD","PRAD","HNSC","THYM","LUSC","LIHC","OV","SKCM","CHOL","DLBC","MESO","UCEC","ACC","KICH","UCS","THCA","SARC","UVM","KIRP","GBMLGG","PCPG")
IMMSIG<-c("CD8_infiltration","Tcell_dysfunction_pos","Tcell_inflammation","IFNG","Inflamed","IMPRES","Cytolytic_activity","PD1","PDL1","CTL_IFNG","CD8_genes","CTL_genes","Tcell_exclusion")
CTL<-c("CD8A","CD8B", "GZMA", "GZMB","PRF1")
CTL_IFNG<-c("IFNG","STAT1","IDO1","CXCL9","CXCL10","HLA-DRA")
CD274<-c("CD274")
PDCD1<-c("PDCD1")
CYT<-c("PRF1","GZMA")
CD8<-c("CD8A","CD8B")
GDF15<-c("GDF15")

RHO_MAT<-data.frame(matrix(0.0,nrow=13,ncol=30))
names(RHO_MAT)<-CANCER
row.names(RHO_MAT)<-IMMSIG
PVAL_MAT<-data.frame(matrix(0.0,nrow=13,ncol=30))
names(PVAL_MAT)<-CANCER
row.names(PVAL_MAT)<-IMMSIG

IMM<-read.table("PANC_IMMUNE_SIGNATURES_SSGSEA.txt", header=TRUE, sep="\t", dec = ".", check.names=FALSE)
IMP<-read.table("PANC_IMPRES.txt", header=TRUE, sep="\t", dec = ".", check.names=FALSE)
FR<-read.table("PANC_QUANTISEQ.tsv", header=TRUE, sep="\t", dec = ".", check.names=FALSE)


for (i in 1:length(CANCER)) {
  i<-9
    IMM1<-IMM[which(IMM$project==CANCER[i]),]
    IMP1<-IMP[which(IMP$project==CANCER[i]),]
    FR1<-FR[which(FR$cell_type=="T cell CD8+" & FR$project==CANCER[i]),]
    FILE<-paste(CANCER[i],"_primary_TPM.tsv",sep="")
    EX<-read.table(FILE, header=TRUE, sep="\t", dec = ".", check.names=FALSE,row.names=1)
    EX1<-log2(EX+1)
    
    CYT_SCORE<-apply(EX1[CYT,],2,mean,na.rm=TRUE)
    CD8_SCORE<-apply(EX1[CD8,],2,mean,na.rm=TRUE)
    CTL_SCORE<-apply(EX1[CTL,],2,mean,na.rm=TRUE)
    CTL_IFNG_SCORE<-apply(EX1[CTL_IFNG,],2,mean,na.rm=TRUE)
    IFNG_SCORE<-apply(EX1[IFNG_SIGNATURE,],2,mean,na.rm=TRUE)
    GDF15_SCORE<-EX1[GDF15,]
    PD1_SCORE<-EX1[PDCD1,]
    PDL1_SCORE<-EX1[CD274,]
    
    DF1<-data.frame(SAMPLE=names(EX1),GDF15=as.numeric(GDF15_SCORE),CD8=as.numeric(CD8_SCORE),PDL1=as.numeric(PDL1_SCORE),PD1=as.numeric(PD1_SCORE),CTL=as.numeric(CTL_SCORE),IFNG=as.numeric(IFNG_SCORE),CTL_IFNG=as.numeric(CTL_IFNG_SCORE),CYT=as.numeric(CYT_SCORE))
    row.names(DF1)<-DF1$SAMPLE
    
    CYT_COR<-cor.test(DF1$GDF15,DF1$CYT, method="spearman")
    CD8_COR<-cor.test(DF1$GDF15,DF1$CD8, method="spearman")
    CTL_COR<-cor.test(DF1$GDF15,DF1$CTL, method="spearman")
    CTL_IFNG_COR<-cor.test(DF1$GDF15,DF1$CTL_IFNG, method="spearman")
    IFNG_COR<-cor.test(DF1$GDF15,DF1$IFNG, method="spearman")
    PD1_COR<-cor.test(DF1$GDF15,DF1$PD1, method="spearman")
    PDL1_COR<-cor.test(DF1$GDF15,DF1$PDL1, method="spearman")
    
    T_EXCL_SCORE<-IMM1[which(IMM1$cell.type=="T_EXCL"),c(2,3)]
    row.names(T_EXCL_SCORE)<-T_EXCL_SCORE$patient
    DF2<- merge.data.frame(DF1,T_EXCL_SCORE,by = 'row.names')
    T_EXCL_COR<-cor.test(DF2$GDF15,DF2$NES, method="spearman")
    
    INFLAM<-IMM1[which(IMM1$cell.type=="INFLAM"),c(2,3)]
    row.names(INFLAM)<-INFLAM$patient
    DF2<- merge.data.frame(DF1,INFLAM,by = 'row.names')
    INFLAM_COR<-cor.test(DF2$GDF15,DF2$NES, method="spearman")
    
    IFNGS<-IMM1[which(IMM1$cell.type=="IFNG"),c(2,3)]
    row.names(IFNGS)<-IFNGS$patient
    DF2<- merge.data.frame(DF1,IFNGS,by = 'row.names')
    IFNGS_COR<-cor.test(DF2$GDF15,DF2$NES, method="spearman")
    
    EXPAND_IMM<-IMM1[which(IMM1$cell.type=="EXPAND_IMM"),c(2,3)]
    row.names(EXPAND_IMM)<-EXPAND_IMM$patient
    DF2<- merge.data.frame(DF1,EXPAND_IMM,by = 'row.names')
    EXPAND_IMM_COR<-cor.test(DF2$GDF15,DF2$NES, method="spearman")
    
    TCELL_INFLAM<-IMM1[which(IMM1$cell.type=="TCELL_INFLAM"),c(2,3)]
    row.names(TCELL_INFLAM)<-TCELL_INFLAM$patient
    DF2<- merge.data.frame(DF1,TCELL_INFLAM,by = 'row.names')
    TCELL_INFLAM_COR<-cor.test(DF2$GDF15,DF2$NES, method="spearman")
    
    TCELL_DYSF_POS<-IMM1[which(IMM1$cell.type=="TCELL_DYSF_POS"),c(2,3)]
    row.names(TCELL_DYSF_POS)<-TCELL_DYSF_POS$patient
    DF2<- merge.data.frame(DF1,TCELL_DYSF_POS,by = 'row.names')
    TCELL_DYSF_POS_COR<-cor.test(DF2$GDF15,DF2$NES, method="spearman")
    
    IMPRES<-IMP1[,c(2,3)]
    row.names(IMPRES)<-IMPRES$patient
    DF2<- merge.data.frame(DF1,IMPRES,by = 'row.names')
    IMPRES_COR<-cor.test(DF2$GDF15,DF2$impres.score, method="spearman")
    
    CD8_INFILT<-FR1[,c(2,3)]
    row.names(CD8_INFILT)<-CD8_INFILT$sample
    DF2<- merge.data.frame(DF1,CD8_INFILT,by = 'row.names')
    CD8_INFILT_COR<-cor.test(DF2$GDF15,DF2$fraction_value, method="spearman")
    
    RHO_MAT[1,i]<-CD8_INFILT_COR$estimate
    RHO_MAT[2,i]<-TCELL_DYSF_POS_COR$estimate
    RHO_MAT[3,i]<-TCELL_INFLAM_COR$estimate
    RHO_MAT[4,i]<-IFNGS_COR$estimate    
    RHO_MAT[5,i]<-INFLAM_COR$estimate    
    RHO_MAT[6,i]<-IMPRES_COR$estimate    
    RHO_MAT[7,i]<-CYT_COR$estimate   
    RHO_MAT[8,i]<-PD1_COR$estimate    
    RHO_MAT[9,i]<-PDL1_COR$estimate    
    RHO_MAT[10,i]<-CTL_IFNG_COR$estimate    
    RHO_MAT[11,i]<-CD8_COR$estimate    
    RHO_MAT[12,i]<-CTL_COR$estimate    
    r1<-INFLAM_COR$estimate
    z1<- 0.5*(log(1+r1) - log(1-r1))
    r2<-T_EXCL_COR$estimate 
    z2<- 0.5*(log(1+r2) - log(1-r2))
    RHO_MAT[13,i]<-z1-z2
    PVAL_MAT[1,i]<-CD8_INFILT_COR$p.value
    PVAL_MAT[2,i]<-TCELL_DYSF_POS_COR$p.value
    PVAL_MAT[3,i]<-TCELL_INFLAM_COR$p.value
    PVAL_MAT[4,i]<-IFNGS_COR$p.value    
    PVAL_MAT[5,i]<-INFLAM_COR$p.value    
    PVAL_MAT[6,i]<-IMPRES_COR$p.value    
    PVAL_MAT[7,i]<-CYT_COR$p.value   
    PVAL_MAT[8,i]<-PD1_COR$p.value    
    PVAL_MAT[9,i]<-PDL1_COR$p.value    
    PVAL_MAT[10,i]<-CTL_IFNG_COR$p.value    
    PVAL_MAT[11,i]<-CD8_COR$p.value    
    PVAL_MAT[12,i]<-CTL_COR$p.value    
    PVAL_MAT[13,i]<-T_EXCL_COR$p.value    
}
write.table(RHO_MAT,"GDF15_RHO_INDICATIONS.txt",row.names=TRUE,col.names=NA, sep="\t",quote=FALSE)
write.table(PVAL_MAT,"GDF15_PVAL_INDICATIONS.txt",row.names=TRUE,col.names=NA, sep="\t",quote=FALSE)
FDR_MAT<-PVAL_MAT
for (s in 1:ncol(PVAL_MAT)) {
  P<-PVAL_MAT[,s]
  adjpBH<-mt.rawp2adjp(P, proc="BH")
  FDR_MAT[,s]<-adjpBH$adjp[order(adjpBH$index),2]
}
write.table(FDR_MAT,"GDF15_FDR_INDICATIONS.txt",row.names=TRUE,col.names=NA, sep="\t",quote=FALSE)

RES_MAT<-RHO_MAT
for (k in 1:ncol(RHO_MAT)) {
  for (j in 1:12) {
    if (RHO_MAT[j,k]>0.1 & FDR_MAT[j,k]<0.1) {
        RES_MAT[j,k]<- 1
    } else {
      if (RHO_MAT[j,k]< -0.1 & FDR_MAT[j,k]<0.1) {
        RES_MAT[j,k]<- -1
      } else {
        RES_MAT[j,k]<- 0
      }
    }
  }
  if (RHO_MAT[13,k]>0.2) {
    RES_MAT[13,k]<- 1
  } else {
    if (RHO_MAT[13,k]< -0.2) {
      RES_MAT[13,k]<- -1
    } else {
      RES_MAT[13,k]<- 0
    }
  }
}  
write.table(RES_MAT,"GDF15_RES_INDICATIONS.txt",row.names=TRUE,col.names=NA, sep="\t",quote=FALSE)

