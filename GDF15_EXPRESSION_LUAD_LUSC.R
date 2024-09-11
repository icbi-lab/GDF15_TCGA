library(ggplot2)

# Distribution of GDF15 expression (log2 (TPM+1)) in LUAD and LUSC
# Violin plot with boxplot (Fig 2b)
# Input: LUAD_primary_TPM.tsv, LUSC_primary_TPM.tsv
# Output: LUAD_GDF15.txt, LUSC_GDF15.txt, GDF15_LUAD_LUSC.pdf

LUAD<-read.table("LUAD_primary_TPM.tsv", header=TRUE, sep="\t", dec = ".", row.names=1, check.names=FALSE)
LUAD_GDF15<-log2(as.numeric(LUAD["GDF15",])+1)
LUAD_GDF15_DF<-data.frame(PATID=names(LUAD), GDF15=LUAD_GDF15)
write.table(LUAD_GDF15_DF,"LUAD_GDF15.txt",row.names=TRUE,col.names=NA, sep="\t", quote=FALSE)

LUSC<-read.table("LUSC_primary_TPM.tsv", header=TRUE, sep="\t", dec = ".", row.names=1, check.names=FALSE)
LUSC_GDF15<-log2(as.numeric(LUSC["GDF15",])+1)
LUSC_GDF15_DF<-data.frame(PATID=names(LUSC), GDF15=LUSC_GDF15)
write.table(LUSC_GDF15_DF,"LUSC_GDF15.txt",row.names=TRUE,col.names=NA, sep="\t", quote=FALSE)

F3<-data.frame(GDF15=c(LUAD_GDF15,LUSC_GDF15), TUM=c(rep("LUAD",length(LUAD_GDF15)),rep("LUSC",length(LUSC_GDF15))))
W<-wilcox.test(GDF15~TUM,data=F3)
if (W$p.value<0.0001) {
  pvtext<-"p<0.0001"
} else {
  pvtext<-paste("p=",sprintf("%.4f",W$p.value),sep="")
}
e <- ggplot(F3, aes(x = TUM, y = GDF15))
e +  theme_classic() + 
  geom_violin(aes(fill = TUM), size=0.8, trim = FALSE) +
  geom_boxplot(width = 0.2,lwd=0.8) +
  labs(y= "GDF15 log2(TPM+1)") +
  theme(plot.margin=margin(1, 1, 1, 1,"cm")) +
  geom_segment(aes(x = 1, y = 14, xend = 2, yend = 14),size=0.8) +
  annotate("text", x=1.5, y=14.9, label= pvtext,size=7) + 
  scale_fill_manual(values = c("red", "blue"))+
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(size = 20,vjust=-2,color = "black")) +
  theme(axis.text.y = element_text(margin = margin(r = 10),size = 20,hjust=1,color = "black")) +
  theme(axis.title.y = element_text(margin = margin(r = 14),color = "black", size = 20)) +
  theme(axis.title.x = element_text(color = "white", size = 10)) +
  theme(axis.line=element_line(linewidth = 0.8)) +
  theme(axis.ticks.length=unit(-0.3, "cm")) +
  theme(axis.ticks=element_line(linewidth=0.8)) +
  ylim(0,15)
  ggsave("GDF15_LUAD_LUSC.pdf",width=5, height=6)

