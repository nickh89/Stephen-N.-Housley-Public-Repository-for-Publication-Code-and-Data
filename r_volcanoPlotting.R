
library(EnhancedVolcano)
  library(airway)
  library(magrittr)
library("DESeq2")

library(ggplot2)
library(ggpubr)





  
Res1<-read.csv(file = "Desktop/Projects /Papers/GT_2018_CancerChemoInteraction_Paper_1/GT_2018_CancerChemoInteraction_Paper_1/Data/Transcriptome/volcanoPlotting/EnhancedVolcano_Data.csv")
Res1<- as.data.frame(Res1)
channelstrings<-read.csv(file = "Desktop/Projects /Papers/GT_2018_CancerChemoInteraction_Paper_1/GT_2018_CancerChemoInteraction_Paper_1/Data/Transcriptome/volcanoPlotting/ChannelStrings.csv")
channelstrings$sodium<-as.character(channelstrings$sodium)
channelstrings<-as.data.frame(channelstrings)
  # create custom key-value pairs for high, low, mid expression
  keyvals <- rep("grey50", nrow(Res1))
  names(keyvals) <- rep("Mid / NA", nrow(Res1))
  
  keyvals[which(Res1$Log2FC_Con_POX > 1.0 & Res1$p.value.POX.control <0.01 )] <- "red"
  names(keyvals)[which(Res1$Log2FC_Con_POX > 1.0 & Res1$p.value.POX.control <0.01)] <- "High"
  
  
  keyvals[which(Res1$Log2FC_Con_POX < -1.0 & Res1$p.value.POX.control <0.01  )] <- "royalblue"
  names(keyvals)[which(Res1$Log2FC_Con_POX < -1.0 & Res1$p.value.POX.control <0.01  )] <- "Low"
  png(file="/Users/nickhousley/Desktop/Projects /Papers/GT_2018_CancerChemoInteraction_Paper_1/GT_2018_CancerChemoInteraction_Paper_1/Figures/Ia Neurons/transcriptome/NorthwesternVolcanoPlotFigure.png",res = 300, width = 3000, height =  2000 ,bg = "transparent")
  EnhancedVolcano(Res1,
                        lab = as.character(Res1$Gene.Symbol),
                        x = "Log2FC_Con_POX",
                        y = "p.value.POX.control",
                        pCutoff = 0.01,
                        FCcutoff = 1.0,
                        # selectLab = Res1$ID[which(names(keyvals) %in% c("High", "Low"))],
                        selectLab = c(""),
                       # selectLab = c("Scn8a","Scn9a","Scn1a","Scnn1a","Scnn1d","Scn7a","Scnn1g","Scnn1b","Kcnc3","Kcnn2","Asic2", "Asic1", "Piezo1","Piezo2","Hcn1","Hcn2","Hcn3","Slc17a7", "Cacna1c", "Cacna1d", "Cacna1s", "Cacna1f"),
                       #selectLab = as.character(channelstrings$pox_potassium),
                        
                        xlab = bquote(~Log[2]~ "fold change"),
                        ylab = bquote(~-Log[10]~adjusted~italic(P)),
                        xlim = c(-4,6),
                        ylim = c(0,12),
                        transcriptLabSize = 6.0,
                        colAlpha = 1,
                        legend = c(""),
                        legendPosition = "bottom",
                        transcriptPointSize = 1.2,
                        # legendLabSize = 15,
                        # legendIconSize = 5.0,
                        DrawConnectors = TRUE,
                        widthConnectors = .5,
                        colConnectors = "black",
                        border = "partial",
                        borderWidth = 1.5,
                        borderColour = "black",
                        gridlines.major = FALSE,
                        gridlines.minor = FALSE,
                        colOverride = keyvals)+
    theme_transparent()


  dev.off()
