library(ggplot2)

# Distribution of GDF15 expression (log2 (TPM+1)) in BLCA LUM and BLCA SQUAM
# Violin plot with boxplot (Fig 2b)
# Input: BLCA_LUM.txt, BLCA_SQUAM.txt
# Output: BLCA_LUM_GDF15.txt, BLCA_SQUAM_GDF15.txt, GDF15_BLCA_LUM_SQUAM.pdf

BLCA_LUM<-read.table("BLCA_LUM.txt", header=TRUE, sep="\t", dec = ".", row.names=1, check.names=FALSE)
BLCA_LUM_GDF15_DF<-data.frame(PATID=row.names(BLCA_LUM), GDF15=BLCA_LUM$GDF15)
write.table(BLCA_LUM_GDF15_DF,"BLCA_LUM_GDF15.txt",row.names=TRUE,col.names=NA, sep="\t", quote=FALSE)

BLCA_SQUAM<-read.table("BLCA_SQUAM.txt", header=TRUE, sep="\t", dec = ".", row.names=1, check.names=FALSE)
BLCA_SQUAM_GDF15_DF<-data.frame(PATID=row.names(BLCA_SQUAM), GDF15=BLCA_SQUAM$GDF15)
write.table(LUSC_GDF15_DF,"LUSC_GDF15.txt",row.names=TRUE,col.names=NA, sep="\t", quote=FALSE)

F3<-data.frame(GDF15=c(BLCA_LUM$GDF15,BLCA_SQUAM$GDF15), TUM=c(rep("BLCA_LUM",length(BLCA_LUM$GDF15)),rep("BLCA_SQUAM",length(BLCA_SQUAM$GDF15))))

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
  ggsave("GDF15_BLCA_LUM_SQUAM.pdf",width=5, height=6)

