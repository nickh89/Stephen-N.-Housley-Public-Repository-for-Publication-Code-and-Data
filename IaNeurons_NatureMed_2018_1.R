rm(list=ls())  # Careful! This clears all of R's memory!
graphics.off() # This closes all of R's graphics windows.


###load in required libraries
require(FactoMineR)
require(factoextra)
require(missMDA)
require(corrplot)
require(MASS)
require(scatterplot3d)
require(devtools)
require(ggplot2)
require(pca3d)
require(ggord)
require(parallel)
require(devtools)
require(ggpubr)
require(gridExtra)
require(dplyr)
require(tibble)
require(PerformanceAnalytics)
require(rgl)


source("/Users/nickhousley/Desktop/Projects /Papers/GT_2018_Cancer_Chemo_Interaction_Paper/GT_2018_CancerChemoInteraction_Paper/Code/Neurophysiology/multiplot.R")
paletteNature2018 = c("#808080","#3e5ce0","#db0415","#933bd3")


### read in data

IaData<-read.csv(file="/Users/nickhousley/Desktop/Projects /Papers/GT_2018_CancerChemoInteraction_Paper_1/GT_2018_CancerChemoInteraction_Paper_1/Data/Neurophysiology/Ia/Ia_Final_1.csv", header=TRUE, sep=",", check.names = FALSE) # read in data and filter
head(IaData)

### reformat to something useful 1. Raw Data is the physiological recorded numbers 2. Standardized is mean-centered/standardized 3. is #2 multiplied by the LDA coefficients 

DataIaRaw<-IaData[,c(2,5,9:25,29:43,117)]
DataIaStandardized<-IaData[,c(2,5,47:62,66:80,117)]
#DataAllScaled<-IaData[,c(2,5)]

chart.Correlation(IaData[,c(66:80)], histogram=TRUE, pch=19)
chart.Correlation(IaData[,c(47:62)], histogram=TRUE, pch=19)


### quick MANOVA check on the independent variables for 'significance'
res.man <- manova(cbind(Dyn.pfr,  Stat.mSfr ,  Stat.afr , Stat.HzSTD , Stat.HzEnd , Dyn.IB, Thr.T, Thr.F ,  Thr.L,  Stat.LSpkT ,  Stat.LSpkF ,  Dyn.spkNum ,  Stat.spkNum ,  Dyn.DI , Dyn.slp ,  Stat.slp , Dyn.F ,Dyn.pfr1 ,Dyn.pfr3 , Thr.T1 , Thr.F1 ,  Thr.L1 ,  Thr.T3 ,Thr.F3 , Thr.L3 , Dyn.slp1 , Dyn.slp3 ,Dyn.spkNum1, Dyn.spkNum3,Dyn.RDR,  Dyn.IFRdrop)~ treatment, data = DataIaStandardized)
summary(res.man, tol = 0)###overall
summary.aov(res.man)###individual tests


res.man <- aov(DataIaRaw$Thr.L ~ as.factor(treatment), data = DataIaRaw)
summary.aov(res.man)###individual tests
TukeyHSD(res.man, conf.level = 0.95)


### 
DataIaRaw %>%
  group_by(treatment)%>%
  summarise(meanForce=mean(Dyn.RDR,na.rm = TRUE))

df.sum<-DataIaRaw %>%
  group_by(treatment)%>%
  summarise(meanForce=mean(Stat.HzEnd))
            
            ,na.rm = TRUE), STDerr=stdErr(Stat.HzEnd),stDEV=StdDev(Stat.HzEnd))


ggplot(df.sum, aes(x=treatment, y=df.sum$meanForce, ymin = (df.sum$meanForce-df.sum$STDerr), ymax = (df.sum$meanForce+df.sum$STDerr)),fill=treatment) +
  geom_bar(stat = "identity")+
  geom_errorbar(
  width=.2,                    # Width of the error bars
position=position_dodge(.9))+
  scale_fill_manual(values=c(paletteNature2018),"Treatment")+
  theme_transparent()+
  
  theme(axis.text.x = element_text(face="bold", color="black", 
                                   size=15, angle=0),
        axis.title.x = element_blank(),
        axis.ticks = element_blank(),
        axis.text.y = element_text(face="bold", color="black", 
                                   size=15, angle=0),
        axis.title.y = element_blank(),
        plot.title = element_text(face="bold", color="black", 
                                  size=15, angle=0),
        axis.title = element_text(face="bold", color="black", 
                                  size=15, angle=0),
        legend.position="none")



############## ############## ############## Analysis with just WT and POX ############## ############## ############## 
paletteNature2018WT_POX = c("#808080","#933bd3")




DataIaWT_POX<-DataIaStandardized %>% filter(treatment == 4 | treatment == 1)
DataIaWT_POX_pca<-DataIaWT_POX[,c(-2)] #removes the autoClass column
DataIaWT_POX.pca<-PCA(DataIaWT_POX_pca[,-1],graph=T) #PCA
eig.val <- get_eigenvalue(DataIaWT_POX.pca)

fviz_eig(DataIaWT_POX.pca, addlabels = TRUE, ylim = c(0, 50)) # generate scree plot
var <- get_pca_var(DataIaWT_POX.pca)

# Contributions of variables to PC1
fviz_contrib(DataIaWT_POX.pca, choice = "var", axes = 1, top = 10)
# Contributions of variables to PC2
fviz_contrib(DataIaWT_POX.pca, choice = "var", axes = 2, top = 10) 

fviz_pca_biplot(DataIaWT_POX.pca, 
                # Fill individuals by groups
                geom.ind = "point",
                pointshape = 21,
                pointsize = 1.5,
                palette = paletteNature2018WT_POX,
                alpha.var = 1,
                alpha.ind = 1,
                addEllipses = T, ellipse.type = "norm", ellipse.level = 0.68, ellipse.alpha=.5,
                fill.ind = as.factor(DataIaWT_POX_pca$treatment),
                col.ind = "black",
                col.var = factor(c("Dynamic","Static","Static","Static","Static",
                                   "Dynamic","Threshold","Threshold","Threshold",
                                   "Static","Static","Dynamic","Static",
                                   "Dynamic","Dynamic","Static","Dynamic","Dynamic","Dynamic","Threshold","Threshold","Threshold","Threshold","Threshold","Threshold","Dynamic","Dynamic","Dynamic","Dynamic","Dynamic","Dynamic")),
                legend.title = list(fill = "treatment", color = "Clusters"),
                repel = TRUE        # Avoid label overplotting
)+
  #  ggpubr::fill_palette("paletteNature2018")+      # Indiviual fill color
  ggpubr::color_palette("npg")      # Variable colors


data_all<-t(DataIaWT_POX_pca)
colnames(data_all)<-(data_all[1,])
data_all = data_all[-1, ] 

colnames(data_all)   
colnames(data_all)[1:11]<-"WT"
colnames(data_all)[12:21]<-"Pirc+OX"
colnames(data_all)

X=data_all

X = t(scale(t(data_all),center=TRUE,scale=FALSE))
dim(X)
dim(data_all)
# we transpose X again for svd
sv = svd(t(X))
U = sv$u
V = sv$v
D = sv$d

plot(sv$d)

## in R calculate the rank of a matrix is by
qr(t(X))$rank


##SVD of PCA Space for PC1&2, unscaled
cols = as.numeric(as.factor(colnames(data_all)))
plot(U[,1],U[,2],type="n",xlab="PC1",ylab="PC2")
text(U[,1],U[,2],colnames(X),col=cols)

##SVD of PCA Space for PC1&2, scaled
par(mfrow=c(1,1))
Z = t(X)%*%V

#### 3D PCA
#### http://www.sthda.com/english/wiki/amazing-interactive-3d-scatter-plots-r-software-and-data-visualization
colors<-c("#933bd3","#808080")
colors<-colors[as.numeric(cols)]
png(file="/Users/nickhousley/Desktop/Projects /Papers/GT_2018_CancerChemoInteraction_Paper_1/GT_2018_CancerChemoInteraction_Paper_1/Figures/Ia Neurons/3dPCAPlotWT_POX_noTicks.png",res = 600, width = 6000, height =  6000 ,bg = "transparent")

S3D<- scatterplot3d(Z[,1],Z[,2],Z[,3], cex.symbols=5, color = colors, pch =16, angle = 75,xlab = "PC1", ylab = "PC2", zlab = "PC3", zlim = c(-10,10), ylim = c(-10,10), box = FALSE, grid = FALSE, tick.marks=FALSE)
dev.off()

DataIaWT_POX.lda<- lda(treatment ~ Dyn.pfr + 
                               Stat.mSfr + 
                               Stat.afr + 
                               Stat.HzSTD + 
                               Stat.HzEnd + 
                               Dyn.IB+ 
                               Thr.T + 
                               Thr.F + 
                               Thr.L + 
                               Stat.LSpkT + 
                               Stat.LSpkF + 
                               Dyn.spkNum + 
                               Stat.spkNum + 
                               Dyn.DI + 
                               Dyn.slp + 
                               Stat.slp +
                               Dyn.F +
                               Dyn.pfr1 +
                               Dyn.pfr3 +
                               Thr.T1 +
                               Thr.F1 +
                               Thr.L1 +
                               Thr.T3 +
                               Thr.F3 +
                               Thr.L3 +
                               Dyn.slp1 +
                               Dyn.slp3 +
                               Dyn.spkNum1 +
                               Dyn.spkNum3 +
                               Dyn.RDR +
                               Dyn.IFRdrop, 
                             data = DataIaWT_POX_pca, CV=TRUE)



discriminantFunctionsIa<-as.data.frame(DataIaStandardized.lda$scaling)
discriminantFunctionsIa<- discriminantFunctionsIa%>% rownames_to_column('features')
discriminantFunctionsIaLD1<-discriminantFunctionsIa  %>% arrange(desc(LD1)) 
discriminantFunctionsIaLD2<-discriminantFunctionsIa  %>% arrange(desc(LD2)) 
discriminantFunctionsIaLD3<-discriminantFunctionsIa  %>% arrange(desc(LD3)) 

### save the discriminant functions
AA<-discriminantFunctionsIaLD1[,c(1,2)]
AA[,c(3,4)]<-discriminantFunctionsIaLD2[,c(1,3)]
AA[,c(5,6)]<-discriminantFunctionsIaLD3[,c(1,4)]


############## ############## ############## ############## All Groups ############## ############## ############## ############## 
### perform the LDA on the entire dataset of interest, includes leave-one-out cross validation/jackknife validation (CV)
DataIaStandardized.lda<- lda(treatment ~ Dyn.pfr + 
                                Stat.mSfr + 
                                Stat.afr + 
                                Stat.HzSTD + 
                                Stat.HzEnd + 
                                Dyn.IB+ 
                                Thr.T + 
                                Thr.F + 
                                Thr.L + 
                                Stat.LSpkT + 
                                Stat.LSpkF + 
                                Dyn.spkNum + 
                                Stat.spkNum + 
                                Dyn.DI + 
                                Dyn.slp + 
                                Stat.slp +
                                Dyn.F +
                                Dyn.pfr1 +
                                Dyn.pfr3 +
                                Thr.T1 +
                                Thr.F1 +
                                Thr.L1 +
                                Thr.T3 +
                                Thr.F3 +
                                Thr.L3 +
                                Dyn.slp1 +
                                Dyn.slp3 +
                                Dyn.spkNum1 +
                                Dyn.spkNum3 +
                                Dyn.RDR +
                                Dyn.IFRdrop, 
                              data = DataIaStandardized)


discriminantFunctionsIa<-as.data.frame(DataIaStandardized.lda$scaling)
discriminantFunctionsIa<- discriminantFunctionsIa%>% rownames_to_column('features')
discriminantFunctionsIaLD1<-discriminantFunctionsIa  %>% arrange(desc(LD1)) 
discriminantFunctionsIaLD2<-discriminantFunctionsIa  %>% arrange(desc(LD2)) 
discriminantFunctionsIaLD3<-discriminantFunctionsIa  %>% arrange(desc(LD3)) 

### save the discriminant functions
AA<-discriminantFunctionsIaLD1[,c(1,2)]
AA[,c(3,4)]<-discriminantFunctionsIaLD2[,c(1,3)]
AA[,c(5,6)]<-discriminantFunctionsIaLD3[,c(1,4)]

png(file="/Users/nickhousley/Desktop/Projects /Papers/GT_2018_Cancer_Chemo_Interaction_Paper/GT_2018_CancerChemoInteraction_Paper/Figures/Ia Neurons/discriminantFunctionsIa.png",res = 600, width = 6000, height =  6000 ,bg = "transparent")
grid.table(AA)
dev.off()
write.csv(AA, file ="/Users/nickhousley/Desktop/Projects /Papers/GT_2018_CancerChemoInteraction_Paper_1/GT_2018_CancerChemoInteraction_Paper_1/Data/Neurophysiology/Ia/IaDiscriminantFunction.csv")



# # ## leave-one-out cross validation
# ct <- table(DataIaStandardized$treatment, DataIaStandardized.lda$class)
# diag(prop.table(ct, 1))
# sum(diag(prop.table(ct)))
# # 


table(DataIaStandardized$treatment, DataIaStandardized.lda$class, dnn = c('Actual Group','Predicted Group'))


### posterior prediction check on all data
training_sample <- sample(c(TRUE, FALSE), nrow(DataIaStandardized), replace = T, prob = c(0.6,0.4))
train <- DataIaStandardized[training_sample, ]
test <- DataIaStandardized[!training_sample, ]

lda.test <- predict(DataIaStandardized.lda,test)
test$lda <- lda.test$class

pred.train <- predict(DataIaStandardized.lda,train,method = c("plug-in", "predictive", "debiased"))$class
pred.test <- predict(DataIaStandardized.lda,test,method = c("plug-in", "predictive", "debiased"))$class
#accuracy on training data
TrainPerformance<-mean(pred.train == train$treatment)
TestPerformance<-mean(pred.test == test$treatment)

###barplot of model accuracy 
png(file="/Users/nickhousley/Desktop/Projects /Papers/GT_2018_Cancer_Chemo_Interaction_Paper/GT_2018_CancerChemoInteraction_Paper/Figures/Ia Neurons/LDAmodel accuracy.png",res = 600, width = 6000, height =  6000 ,bg = "transparent")
barplot(rbind(TestPerformance, TrainPerformance),beside = T, ylab = "model accuracy", legend=row.names(rbind(TestPerformance, TrainPerformance)),args.legend = list(x = "topright", bty = "n",inset=c(0,1)))
dev.off()


### save the confusion matrix of model performance
png(file="/Users/nickhousley/Desktop/Projects /Papers/GT_2018_Cancer_Chemo_Interaction_Paper/GT_2018_CancerChemoInteraction_Paper/Figures/Ia Neurons/confusion_matrixAll.png",res = 600, width = 6000, height =  4800 ,bg = "transparent")
confusion_matrix<-table(test$lda,test$treatment)
confusion_matrix<-tableGrob(confusion_matrix)
confusion_matrix<-grid.arrange(confusion_matrix)
dev.off()

lda.train1 <- predict(DataIaStandardized.lda)
train$lda <- lda.train$class
table(train$lda,train$Species)


###generate the canonical functions (basis vectors for the 3d plotting in Matlab)
CanonicalVariables<-as.matrix(DataIaStandardized[,c(3:33)]) %*% as.matrix(DataIaStandardized.lda$scaling)
CanonicalVariables<-as.data.frame(CanonicalVariables)
CanonicalVariables$treatment<-DataIaStandardized$treatment
write.csv(CanonicalVariables, file ="/Users/nickhousley/Desktop/Projects /Papers/GT_2018_Cancer_Chemo_Interaction_Paper/GT_2018_CancerChemoInteraction_Paper/Data/Ia/CanonicalVariablesAll.csv")


### save the 2D LDA plot
png(file="/Users/nickhousley/Desktop/Projects /Papers/GT_2018_Cancer_Chemo_Interaction_Paper/GT_2018_CancerChemoInteraction_Paper/Figures/Ia Neurons/DataIaStandardized.lda.png",res = 600, width = 6000, height =  4800 ,bg = "transparent")
ggord(DataIaStandardized.lda, as.factor(DataIaStandardized$treatment), ellipse_pro =0.68, cols = c(paletteNature2018), repel = TRUE, size=3, arrow=0.4, alpha_el=0.8)+
  theme(legend.justification=c(1,1), legend.position=c(1,1), legend.box.just = c("top"), legend.background = element_rect(fill=alpha(0.9)), legend.direction = "horizontal",legend.box="horizontal")+
  ggtitle("Ia Neurons")+
  labs(x="Canonical Variable 1", y= "Canonical Variable 2")+
  #  ylim(-2,4)+
  #theme(panel.background = element_blank())+
  theme_transparent()+
  theme(axis.text.x = element_text(face="bold", color="black", 
                                   size=15, angle=0),
        #axis.title.x = element_blank(),
        axis.text.y = element_text(face="bold", color="black", 
                                   size=15, angle=0),
        # axis.title.y = element_blank(),
        plot.title = element_text(face="bold", color="black", 
                                  size=15, angle=0),
        axis.title = element_text(face="bold", color="black", 
                                  size=15, angle=0),
        legend.position="none")
dev.off()





### boxplots for Canonical Variable Comparisons
my_comparisons <- list( c("1", "2"), c("1", "3"), c("1", "4"),c("2", "3"),c("2", "4"),c("3", "4") )

png(file="/Users/nickhousley/Desktop/Projects /Papers/GT_2018_Cancer_Chemo_Interaction_Paper/GT_2018_CancerChemoInteraction_Paper/Figures/Ia Neurons/CanonicaltreatmentgroupsIaCanVar1.png",res = 600, width = 3000, height =  6000, bg = "transparent")
ggplot(data = CanonicalVariables, 
       aes(x=as.factor(treatment), 
           y= LD1, 
           fill=as.factor(treatment)))+
  geom_boxplot(outlier.shape = NA, position = "dodge",varwidth = TRUE)+
  geom_point(alpha = 0.3, size = 3)+
  scale_fill_manual(values=c(paletteNature2018),"Treatment")+
  
  theme(legend.justification=c(1,1), legend.position=c(1,1), legend.box.just = c("top"), legend.background = element_rect(fill=alpha(0.4)), legend.direction = "horizontal",legend.box="horizontal")+
  xlab("Treatment")+
  
  #  ylim(-2,4)+
  #theme(panel.background = element_blank())+
  theme_transparent()+
  
  theme(axis.text.x = element_text(face="bold", color="black", 
                                   size=15, angle=0),
        axis.title.x = element_blank(),
        axis.ticks = element_blank(),
        axis.text.y = element_text(face="bold", color="black", 
                                   size=15, angle=0),
        axis.title.y = element_blank(),
        plot.title = element_text(face="bold", color="black", 
                                  size=15, angle=0),
        axis.title = element_text(face="bold", color="black", 
                                  size=15, angle=0),
        legend.position="none")
  #stat_compare_means(method = "anova", label.y = 10)+
  #stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "t.test",
   #                  ref.group = ".all.")

dev.off()


png(file="/Users/nickhousley/Desktop/Projects /Papers/GT_2018_Cancer_Chemo_Interaction_Paper/GT_2018_CancerChemoInteraction_Paper/Figures/Ia Neurons/CanonicaltreatmentgroupsIaCanVar2.png",res = 600, width = 3000, height =  6000, bg = "transparent")
ggplot(data = CanonicalVariables, 
       aes(x=as.factor(treatment), 
           y= LD2, 
           fill=as.factor(treatment)))+
  geom_boxplot(outlier.shape = NA, position = "dodge",varwidth = TRUE)+
  geom_point(alpha = 0.3, size = 3)+
  scale_fill_manual(values=c(paletteNature2018),"Treatment")+
  
  theme(legend.justification=c(1,1), legend.position=c(1,1), legend.box.just = c("top"), legend.background = element_rect(fill=alpha(0.4)), legend.direction = "horizontal",legend.box="horizontal")+
  xlab("Treatment")+
  
  #  ylim(-2,4)+
  #theme(panel.background = element_blank())+
  theme_transparent()+
  
  theme(axis.text.x = element_text(face="bold", color="black", 
                                   size=15, angle=0),
        axis.title.x = element_blank(),
        axis.ticks = element_blank(),
        axis.text.y = element_text(face="bold", color="black", 
                                   size=15, angle=0),
        axis.title.y = element_blank(),
        plot.title = element_text(face="bold", color="black", 
                                  size=15, angle=0),
        axis.title = element_text(face="bold", color="black", 
                                  size=15, angle=0),
        legend.position="none")+
  #stat_compare_means(method = "anova", label.y = 10)+
  #stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "t.test",
 #                    ref.group = ".all.")

dev.off()


png(file="/Users/nickhousley/Desktop/Projects /Papers/GT_2018_Cancer_Chemo_Interaction_Paper/GT_2018_CancerChemoInteraction_Paper/Figures/Ia Neurons/CanonicaltreatmentgroupsIaCanVar3.png",res = 600, width = 3000, height =  6000, bg = "transparent")
ggplot(data = CanonicalVariables, 
       aes(x=as.factor(treatment), 
           y= LD3, 
           fill=as.factor(treatment)))+
  geom_boxplot(outlier.shape = NA, position = "dodge",varwidth = TRUE)+
  geom_point(alpha = 0.3, size = 3)+
  scale_fill_manual(values=c(paletteNature2018),"Treatment")+
  
  theme(legend.justification=c(1,1), legend.position=c(1,1), legend.box.just = c("top"), legend.background = element_rect(fill=alpha(0.4)), legend.direction = "horizontal",legend.box="horizontal")+
  xlab("Treatment")+
  
  #  ylim(-2,4)+
  #theme(panel.background = element_blank())+
  theme_transparent()+
  
  theme(axis.text.x = element_text(face="bold", color="black", 
                                   size=15, angle=0),
        axis.title.x = element_blank(),
        axis.ticks = element_blank(),
        axis.text.y = element_text(face="bold", color="black", 
                                   size=15, angle=0),
        axis.title.y = element_blank(),
        plot.title = element_text(face="bold", color="black", 
                                  size=15, angle=0),
        axis.title = element_text(face="bold", color="black", 
                                  size=15, angle=0),
        legend.position="none")+
  # stat_compare_means(method = "anova", label.y = 10)+
  stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "t.test",
                     ref.group = ".all.")

dev.off()




#######PCA
DataIaStandardized_pca<-DataIaStandardized[,c(-2)] #removes the autoClass column
DataIaStandardized.pca<-PCA(DataIaStandardized_pca[,-1],graph=T) #PCA


png(file="/Users/nickhousley/Desktop/Projects /Papers/GT_2018_Cancer_Chemo_Interaction_Paper/GT_2018_CancerChemoInteraction_Paper/Figures/Ia Neurons/IaDataClusteredPCA.png",res = 600, width = 6000, height =  6000, bg = "transparent")
## @knitr All_Afferents_Treatments_Nature_2018_PCA_Graph_clustering
fviz_pca_biplot(DataIaStandardized.pca, 
                # Fill individuals by groups
                geom.ind = "point",
                pointshape = 21,
                pointsize = 1.5,
                palette = paletteNature2018,
                alpha.var = 1,
                alpha.ind = 1,
                addEllipses = T, ellipse.type = "norm", ellipse.level = 0.68, ellipse.alpha=.5,
                fill.ind = as.factor(DataIaStandardized$treatment),
                col.ind = "black",
                col.var = factor(c("Dynamic","Static","Static","Static","Static",
                                   "Dynamic","Threshold","Threshold","Threshold",
                                   "Static","Static","Dynamic","Static",
                                   "Dynamic","Dynamic","Static","Dynamic","Dynamic","Dynamic","Threshold","Threshold","Threshold","Threshold","Threshold","Threshold","Dynamic","Dynamic","Dynamic","Dynamic","Dynamic","Dynamic")),
                legend.title = list(fill = "treatment", color = "Clusters"),
                repel = TRUE        # Avoid label overplotting
)+
  #  ggpubr::fill_palette("paletteNature2018")+      # Indiviual fill color
  ggpubr::color_palette("npg")      # Variable colors

dev.off()







png(file="/Users/nickhousley/Desktop/Projects /Papers/GT_2018_Cancer_Chemo_Interaction_Paper/GT_2018_CancerChemoInteraction_Paper/Figures/Ia Neurons/FeatureGraphsDyn.png",res = 600, width = 8000, height =  8000, bg = "transparent")


Dyn.pfr<-ggplot(data = DataIaRaw, 
                aes(x=as.factor(treatment), 
                    y= Dyn.pfr, 
                    #colour=interaction(autoClass,as.factor(treatment))))+ 
                    fill=as.factor(treatment)))+
  geom_boxplot(outlier.shape = NA, position = "dodge",varwidth = TRUE)+
  geom_point(alpha = 0.3, size = 3 )+
  scale_fill_manual(values=c(paletteNature2018),"Treatment")+
  
  theme(legend.justification=c(1,1), legend.position=c(1,1), legend.box.just = c("top"), legend.background = element_rect(fill=alpha(0.4)), legend.direction = "horizontal",legend.box="horizontal")+
  xlab("Treatment")+
  ggtitle("Dyn.pfr")+
  #  ylim(-2,4)+
  #theme(panel.background = element_blank())+
  theme_transparent()+
  
  theme(axis.text.x = element_text(face="bold", color="black", 
                                   size=15, angle=0),
        axis.title.x = element_blank(),
        axis.ticks = element_blank(),
        axis.text.y = element_text(face="bold", color="black", 
                                   size=15, angle=0),
        axis.title.y = element_blank(),
        plot.title = element_text(face="bold", color="black", 
                                  size=15, angle=0),
        axis.title = element_text(face="bold", color="black", 
                                  size=15, angle=0),
        legend.position="none")+
  stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "t.test",
                     ref.group = ".all.") 


Dyn.IB<-ggplot(data = DataIaRaw, 
               aes(x=as.factor(treatment), 
                   y= Dyn.IB, 
                   #colour=interaction(autoClass,as.factor(treatment))))+ 
                   fill=as.factor(treatment)))+
  geom_boxplot(outlier.shape = NA, position = "dodge",varwidth = TRUE)+
  geom_point(alpha = 0.3, size = 3 )+
  scale_fill_manual(values=c(paletteNature2018),"Treatment")+
  
  theme(legend.justification=c(1,1), legend.position=c(1,1), legend.box.just = c("top"), legend.background = element_rect(fill=alpha(0.4)), legend.direction = "horizontal",legend.box="horizontal")+
  xlab("Treatment")+
  ggtitle("Dyn.IB")+
  #  ylim(-2,4)+
  #theme(panel.background = element_blank())+
  theme_transparent()+
  
  theme(axis.text.x = element_text(face="bold", color="black", 
                                   size=15, angle=0),
        axis.title.x = element_blank(),
        axis.ticks = element_blank(),
        axis.text.y = element_text(face="bold", color="black", 
                                   size=15, angle=0),
        axis.title.y = element_blank(),
        plot.title = element_text(face="bold", color="black", 
                                  size=15, angle=0),
        axis.title = element_text(face="bold", color="black", 
                                  size=15, angle=0),
        legend.position="none")+
  stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "t.test",
                     ref.group = ".all.") 




Dyn.spkNum<-ggplot(data = DataIaRaw, 
                   aes(x=as.factor(treatment), 
                       y= Dyn.spkNum, 
                       #colour=interaction(autoClass,as.factor(treatment))))+ 
                       fill=as.factor(treatment)))+
  geom_boxplot(outlier.shape = NA, position = "dodge",varwidth = TRUE)+
  geom_point(alpha = 0.3, size = 3 )+
  scale_fill_manual(values=c(paletteNature2018),"Treatment")+
  
  theme(legend.justification=c(1,1), legend.position=c(1,1), legend.box.just = c("top"), legend.background = element_rect(fill=alpha(0.4)), legend.direction = "horizontal",legend.box="horizontal")+
  xlab("Treatment")+
  ggtitle("Dyn.spkNum")+
  #  ylim(-2,4)+
  #theme(panel.background = element_blank())+
  theme_transparent()+
  
  theme(axis.text.x = element_text(face="bold", color="black", 
                                   size=15, angle=0),
        axis.title.x = element_blank(),
        axis.ticks = element_blank(),
        axis.text.y = element_text(face="bold", color="black", 
                                   size=15, angle=0),
        axis.title.y = element_blank(),
        plot.title = element_text(face="bold", color="black", 
                                  size=15, angle=0),
        axis.title = element_text(face="bold", color="black", 
                                  size=15, angle=0),
        legend.position="none")+
  stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "t.test",
                     ref.group = ".all.") 




Dyn.DI<-ggplot(data = DataIaRaw, 
               aes(x=as.factor(treatment), 
                   y= Dyn.DI, 
                   #colour=interaction(autoClass,as.factor(treatment))))+ 
                   fill=as.factor(treatment)))+
  geom_boxplot(outlier.shape = NA, position = "dodge",varwidth = TRUE)+
  geom_point(alpha = 0.3, size = 3 )+
  scale_fill_manual(values=c(paletteNature2018),"Treatment")+
  
  theme(legend.justification=c(1,1), legend.position=c(1,1), legend.box.just = c("top"), legend.background = element_rect(fill=alpha(0.4)), legend.direction = "horizontal",legend.box="horizontal")+
  xlab("Treatment")+
  ggtitle("Dyn.DI")+
  #  ylim(-2,4)+
  #theme(panel.background = element_blank())+
  theme_transparent()+
  
  theme(axis.text.x = element_text(face="bold", color="black", 
                                   size=15, angle=0),
        axis.title.x = element_blank(),
        axis.ticks = element_blank(),
        axis.text.y = element_text(face="bold", color="black", 
                                   size=15, angle=0),
        axis.title.y = element_blank(),
        plot.title = element_text(face="bold", color="black", 
                                  size=15, angle=0),
        axis.title = element_text(face="bold", color="black", 
                                  size=15, angle=0),
        legend.position="none")+
  stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "t.test",
                     ref.group = ".all.") 



Dyn.slp<-ggplot(data = DataIaRaw, 
                aes(x=as.factor(treatment), 
                    y= Dyn.slp, 
                    #colour=interaction(autoClass,as.factor(treatment))))+ 
                    fill=as.factor(treatment)))+
  geom_boxplot(outlier.shape = NA, position = "dodge",varwidth = TRUE)+
  geom_point(alpha = 0.3, size = 3 )+
  scale_fill_manual(values=c(paletteNature2018),"Treatment")+
  
  theme(legend.justification=c(1,1), legend.position=c(1,1), legend.box.just = c("top"), legend.background = element_rect(fill=alpha(0.4)), legend.direction = "horizontal",legend.box="horizontal")+
  xlab("Treatment")+
  ggtitle("Dyn.slp")+
  #  ylim(-2,4)+
  #theme(panel.background = element_blank())+
  theme_transparent()+
  
  theme(axis.text.x = element_text(face="bold", color="black", 
                                   size=15, angle=0),
        axis.title.x = element_blank(),
        axis.ticks = element_blank(),
        axis.text.y = element_text(face="bold", color="black", 
                                   size=15, angle=0),
        axis.title.y = element_blank(),
        plot.title = element_text(face="bold", color="black", 
                                  size=15, angle=0),
        axis.title = element_text(face="bold", color="black", 
                                  size=15, angle=0),
        legend.position="none")+
  stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "t.test",
                     ref.group = ".all.") 


Dyn.F<-ggplot(data = DataIaRaw, 
              aes(x=as.factor(treatment), 
                  y= Dyn.F, 
                  #colour=interaction(autoClass,as.factor(treatment))))+ 
                  fill=as.factor(treatment)))+
  geom_boxplot(outlier.shape = NA, position = "dodge",varwidth = TRUE)+
  geom_point(alpha = 0.3, size = 3 )+
  scale_fill_manual(values=c(paletteNature2018),"Treatment")+
  
  theme(legend.justification=c(1,1), legend.position=c(1,1), legend.box.just = c("top"), legend.background = element_rect(fill=alpha(0.4)), legend.direction = "horizontal",legend.box="horizontal")+
  xlab("Treatment")+
  ggtitle("Dyn.F")+
  #  ylim(-2,4)+
  #theme(panel.background = element_blank())+
  theme_transparent()+
  
  theme(axis.text.x = element_text(face="bold", color="black", 
                                   size=15, angle=0),
        axis.title.x = element_blank(),
        axis.ticks = element_blank(),
        axis.text.y = element_text(face="bold", color="black", 
                                   size=15, angle=0),
        axis.title.y = element_blank(),
        plot.title = element_text(face="bold", color="black", 
                                  size=15, angle=0),
        axis.title = element_text(face="bold", color="black", 
                                  size=15, angle=0),
        legend.position="none")+
  stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "t.test",
                     ref.group = ".all.") 


Dyn.pfr1<-ggplot(data = DataIaRaw, 
                 aes(x=as.factor(treatment), 
                     y= Dyn.pfr1, 
                     #colour=interaction(autoClass,as.factor(treatment))))+ 
                     fill=as.factor(treatment)))+
  geom_boxplot(outlier.shape = NA, position = "dodge",varwidth = TRUE)+
  geom_point(alpha = 0.3, size = 3 )+
  scale_fill_manual(values=c(paletteNature2018),"Treatment")+
  
  theme(legend.justification=c(1,1), legend.position=c(1,1), legend.box.just = c("top"), legend.background = element_rect(fill=alpha(0.4)), legend.direction = "horizontal",legend.box="horizontal")+
  xlab("Treatment")+
  ggtitle("Dyn.pfr1")+
  #  ylim(-2,4)+
  #theme(panel.background = element_blank())+
  theme_transparent()+
  
  theme(axis.text.x = element_text(face="bold", color="black", 
                                   size=15, angle=0),
        axis.title.x = element_blank(),
        axis.ticks = element_blank(),
        axis.text.y = element_text(face="bold", color="black", 
                                   size=15, angle=0),
        axis.title.y = element_blank(),
        plot.title = element_text(face="bold", color="black", 
                                  size=15, angle=0),
        axis.title = element_text(face="bold", color="black", 
                                  size=15, angle=0),
        legend.position="none")+
  stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "t.test",
                     ref.group = ".all.") 


Dyn.pfr3<-ggplot(data = DataIaRaw, 
                 aes(x=as.factor(treatment), 
                     y= Dyn.pfr3, 
                     #colour=interaction(autoClass,as.factor(treatment))))+ 
                     fill=as.factor(treatment)))+
  geom_boxplot(outlier.shape = NA, position = "dodge",varwidth = TRUE)+
  geom_point(alpha = 0.3, size = 3 )+
  scale_fill_manual(values=c(paletteNature2018),"Treatment")+
  
  theme(legend.justification=c(1,1), legend.position=c(1,1), legend.box.just = c("top"), legend.background = element_rect(fill=alpha(0.4)), legend.direction = "horizontal",legend.box="horizontal")+
  xlab("Treatment")+
  ggtitle("Dyn.pfr3")+
  #  ylim(-2,4)+
  #theme(panel.background = element_blank())+
  theme_transparent()+
  
  theme(axis.text.x = element_text(face="bold", color="black", 
                                   size=15, angle=0),
        axis.title.x = element_blank(),
        axis.ticks = element_blank(),
        axis.text.y = element_text(face="bold", color="black", 
                                   size=15, angle=0),
        axis.title.y = element_blank(),
        plot.title = element_text(face="bold", color="black", 
                                  size=15, angle=0),
        axis.title = element_text(face="bold", color="black", 
                                  size=15, angle=0),
        legend.position="none")+
  stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "t.test",
                     ref.group = ".all.") 



Dyn.slp1<-ggplot(data = DataIaRaw, 
                 aes(x=as.factor(treatment), 
                     y= Dyn.slp1, 
                     #colour=interaction(autoClass,as.factor(treatment))))+ 
                     fill=as.factor(treatment)))+
  geom_boxplot(outlier.shape = NA, position = "dodge",varwidth = TRUE)+
  geom_point(alpha = 0.3, size = 3 )+
  scale_fill_manual(values=c(paletteNature2018),"Treatment")+
  
  theme(legend.justification=c(1,1), legend.position=c(1,1), legend.box.just = c("top"), legend.background = element_rect(fill=alpha(0.4)), legend.direction = "horizontal",legend.box="horizontal")+
  xlab("Treatment")+
  ggtitle("Dyn.slp1")+
  #  ylim(-2,4)+
  #theme(panel.background = element_blank())+
  theme_transparent()+
  
  theme(axis.text.x = element_text(face="bold", color="black", 
                                   size=15, angle=0),
        axis.title.x = element_blank(),
        axis.ticks = element_blank(),
        axis.text.y = element_text(face="bold", color="black", 
                                   size=15, angle=0),
        axis.title.y = element_blank(),
        plot.title = element_text(face="bold", color="black", 
                                  size=15, angle=0),
        axis.title = element_text(face="bold", color="black", 
                                  size=15, angle=0),
        legend.position="none")+
  stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "t.test",
                     ref.group = ".all.") 


Dyn.slp3<-ggplot(data = DataIaRaw, 
                 aes(x=as.factor(treatment), 
                     y= Dyn.slp3, 
                     #colour=interaction(autoClass,as.factor(treatment))))+ 
                     fill=as.factor(treatment)))+
  geom_boxplot(outlier.shape = NA, position = "dodge",varwidth = TRUE)+
  geom_point(alpha = 0.3, size = 3 )+
  scale_fill_manual(values=c(paletteNature2018),"Treatment")+
  
  theme(legend.justification=c(1,1), legend.position=c(1,1), legend.box.just = c("top"), legend.background = element_rect(fill=alpha(0.4)), legend.direction = "horizontal",legend.box="horizontal")+
  xlab("Treatment")+
  ggtitle("Dyn.slp3")+
  #  ylim(-2,4)+
  #theme(panel.background = element_blank())+
  theme_transparent()+
  
  theme(axis.text.x = element_text(face="bold", color="black", 
                                   size=15, angle=0),
        axis.title.x = element_blank(),
        axis.ticks = element_blank(),
        axis.text.y = element_text(face="bold", color="black", 
                                   size=15, angle=0),
        axis.title.y = element_blank(),
        plot.title = element_text(face="bold", color="black", 
                                  size=15, angle=0),
        axis.title = element_text(face="bold", color="black", 
                                  size=15, angle=0),
        legend.position="none")+
  stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "t.test",
                     ref.group = ".all.") 


Dyn.spkNum1<-ggplot(data = DataIaRaw, 
                    aes(x=as.factor(treatment), 
                        y= Dyn.spkNum1, 
                        #colour=interaction(autoClass,as.factor(treatment))))+ 
                        fill=as.factor(treatment)))+
  geom_boxplot(outlier.shape = NA, position = "dodge",varwidth = TRUE)+
  geom_point(alpha = 0.3, size = 3 )+
  scale_fill_manual(values=c(paletteNature2018),"Treatment")+
  
  theme(legend.justification=c(1,1), legend.position=c(1,1), legend.box.just = c("top"), legend.background = element_rect(fill=alpha(0.4)), legend.direction = "horizontal",legend.box="horizontal")+
  xlab("Treatment")+
  ggtitle("Dyn.spkNum1")+
  #  ylim(-2,4)+
  #theme(panel.background = element_blank())+
  theme_transparent()+
  
  theme(axis.text.x = element_text(face="bold", color="black", 
                                   size=15, angle=0),
        axis.title.x = element_blank(),
        axis.ticks = element_blank(),
        axis.text.y = element_text(face="bold", color="black", 
                                   size=15, angle=0),
        axis.title.y = element_blank(),
        plot.title = element_text(face="bold", color="black", 
                                  size=15, angle=0),
        axis.title = element_text(face="bold", color="black", 
                                  size=15, angle=0),
        legend.position="none")+
  stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "t.test",
                     ref.group = ".all.") 


Dyn.spkNum3<-ggplot(data = DataIaRaw, 
                    aes(x=as.factor(treatment), 
                        y= Dyn.spkNum3, 
                        #colour=interaction(autoClass,as.factor(treatment))))+ 
                        fill=as.factor(treatment)))+
  geom_boxplot(outlier.shape = NA, position = "dodge",varwidth = TRUE)+
  geom_point(alpha = 0.3, size = 3 )+
  scale_fill_manual(values=c(paletteNature2018),"Treatment")+
  
  theme(legend.justification=c(1,1), legend.position=c(1,1), legend.box.just = c("top"), legend.background = element_rect(fill=alpha(0.4)), legend.direction = "horizontal",legend.box="horizontal")+
  xlab("Treatment")+
  ggtitle("Dyn.spkNum3")+
  #  ylim(-2,4)+
  #theme(panel.background = element_blank())+
  theme_transparent()+
  
  theme(axis.text.x = element_text(face="bold", color="black", 
                                   size=15, angle=0),
        axis.title.x = element_blank(),
        axis.ticks = element_blank(),
        axis.text.y = element_text(face="bold", color="black", 
                                   size=15, angle=0),
        axis.title.y = element_blank(),
        plot.title = element_text(face="bold", color="black", 
                                  size=15, angle=0),
        axis.title = element_text(face="bold", color="black", 
                                  size=15, angle=0),
        legend.position="none")+
  stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "t.test",
                     ref.group = ".all.") 


Dyn.RDR<-ggplot(data = DataIaRaw, 
                aes(x=as.factor(treatment), 
                    y= Dyn.RDR, 
                    #colour=interaction(autoClass,as.factor(treatment))))+ 
                    fill=as.factor(treatment)))+
  geom_boxplot(outlier.shape = NA, position = "dodge",varwidth = TRUE)+
  geom_point(alpha = 0.3, size = 3 )+
  scale_fill_manual(values=c(paletteNature2018),"Treatment")+
  
  theme(legend.justification=c(1,1), legend.position=c(1,1), legend.box.just = c("top"), legend.background = element_rect(fill=alpha(0.4)), legend.direction = "horizontal",legend.box="horizontal")+
  xlab("Treatment")+
  ggtitle("Dyn.RDR")+
  #  ylim(-2,4)+
  #theme(panel.background = element_blank())+
  theme_transparent()+
  
  theme(axis.text.x = element_text(face="bold", color="black", 
                                   size=15, angle=0),
        axis.title.x = element_blank(),
        axis.ticks = element_blank(),
        axis.text.y = element_text(face="bold", color="black", 
                                   size=15, angle=0),
        axis.title.y = element_blank(),
        plot.title = element_text(face="bold", color="black", 
                                  size=15, angle=0),
        axis.title = element_text(face="bold", color="black", 
                                  size=15, angle=0),
        legend.position="none")+
  stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "t.test",
                     ref.group = ".all.") 


Dyn.IFRdrop<-ggplot(data = DataIaRaw, 
                    aes(x=as.factor(treatment), 
                        y= Dyn.IFRdrop, 
                        #colour=interaction(autoClass,as.factor(treatment))))+ 
                        fill=as.factor(treatment)))+
  geom_boxplot(outlier.shape = NA, position = "dodge",varwidth = TRUE)+
  geom_point(alpha = 0.3, size = 3 )+
  scale_fill_manual(values=c(paletteNature2018),"Treatment")+
  
  theme(legend.justification=c(1,1), legend.position=c(1,1), legend.box.just = c("top"), legend.background = element_rect(fill=alpha(0.4)), legend.direction = "horizontal",legend.box="horizontal")+
  xlab("Treatment")+
  ggtitle("Dyn.IFRdrop")+
  #  ylim(-2,4)+
  #theme(panel.background = element_blank())+
  theme_transparent()+
  
  theme(axis.text.x = element_text(face="bold", color="black", 
                                   size=15, angle=0),
        axis.title.x = element_blank(),
        axis.ticks = element_blank(),
        axis.text.y = element_text(face="bold", color="black", 
                                   size=15, angle=0),
        axis.title.y = element_blank(),
        plot.title = element_text(face="bold", color="black", 
                                  size=15, angle=0),
        axis.title = element_text(face="bold", color="black", 
                                  size=15, angle=0),
        legend.position="none")+
  stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "t.test",
                     ref.group = ".all.") 


multiplot(Dyn.pfr,	Dyn.IB,		Dyn.spkNum,		Dyn.DI,	Dyn.slp,	 Dyn.F,	Dyn.pfr1,	Dyn.pfr3,	Dyn.slp1,	Dyn.slp3,	Dyn.spkNum1,	Dyn.spkNum3,	Dyn.RDR	,Dyn.IFRdrop, cols=4)
dev.off()





png(file="/Users/nickhousley/Desktop/Projects /Papers/GT_2018_Cancer_Chemo_Interaction_Paper/GT_2018_CancerChemoInteraction_Paper/Figures/Ia Neurons/FeatureGraphsStat.png",res = 600, width = 8000, height =  8000, bg = "transparent")

Stat.mSfr<-ggplot(data = DataIaRaw, 
                  aes(x=as.factor(treatment), 
                      y= Stat.mSfr, 
                      #colour=interaction(autoClass,as.factor(treatment))))+ 
                      fill=as.factor(treatment)))+
  geom_boxplot(outlier.shape = NA, position = "dodge",varwidth = TRUE)+
  geom_point(alpha = 0.3, size = 3 )+
  scale_fill_manual(values=c(paletteNature2018),"Treatment")+
  
  theme(legend.justification=c(1,1), legend.position=c(1,1), legend.box.just = c("top"), legend.background = element_rect(fill=alpha(0.4)), legend.direction = "horizontal",legend.box="horizontal")+
  xlab("Treatment")+
  ggtitle("Stat.mSfr")+
  #  ylim(-2,4)+
  #theme(panel.background = element_blank())+
  theme_transparent()+
  
  theme(axis.text.x = element_text(face="bold", color="black", 
                                   size=15, angle=0),
        axis.title.x = element_blank(),
        axis.ticks = element_blank(),
        axis.text.y = element_text(face="bold", color="black", 
                                   size=15, angle=0),
        axis.title.y = element_blank(),
        plot.title = element_text(face="bold", color="black", 
                                  size=15, angle=0),
        axis.title = element_text(face="bold", color="black", 
                                  size=15, angle=0),
        legend.position="none")+
  stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "t.test",
                     ref.group = ".all.") 


Stat.afr<-ggplot(data = DataIaRaw, 
                 aes(x=as.factor(treatment), 
                     y= Stat.afr, 
                     #colour=interaction(autoClass,as.factor(treatment))))+ 
                     fill=as.factor(treatment)))+
  geom_boxplot(outlier.shape = NA, position = "dodge",varwidth = TRUE)+
  geom_point(alpha = 0.3, size = 3 )+
  scale_fill_manual(values=c(paletteNature2018),"Treatment")+
  
  theme(legend.justification=c(1,1), legend.position=c(1,1), legend.box.just = c("top"), legend.background = element_rect(fill=alpha(0.4)), legend.direction = "horizontal",legend.box="horizontal")+
  xlab("Treatment")+
  ggtitle("Stat.afr")+
  #  ylim(-2,4)+
  #theme(panel.background = element_blank())+
  theme_transparent()+
  
  theme(axis.text.x = element_text(face="bold", color="black", 
                                   size=15, angle=0),
        axis.title.x = element_blank(),
        axis.ticks = element_blank(),
        axis.text.y = element_text(face="bold", color="black", 
                                   size=15, angle=0),
        axis.title.y = element_blank(),
        plot.title = element_text(face="bold", color="black", 
                                  size=15, angle=0),
        axis.title = element_text(face="bold", color="black", 
                                  size=15, angle=0),
        legend.position="none")+
  stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "t.test",
                     ref.group = ".all.") 


Stat.HzSTD<-ggplot(data = DataIaRaw, 
                   aes(x=as.factor(treatment), 
                       y= Stat.HzSTD, 
                       #colour=interaction(autoClass,as.factor(treatment))))+ 
                       fill=as.factor(treatment)))+
  geom_boxplot(outlier.shape = NA, position = "dodge",varwidth = TRUE)+
  geom_point(alpha = 0.3, size = 3 )+
  scale_fill_manual(values=c(paletteNature2018),"Treatment")+
  
  theme(legend.justification=c(1,1), legend.position=c(1,1), legend.box.just = c("top"), legend.background = element_rect(fill=alpha(0.4)), legend.direction = "horizontal",legend.box="horizontal")+
  xlab("Treatment")+
  ggtitle("Stat.HzSTD")+
  #  ylim(-2,4)+
  #theme(panel.background = element_blank())+
  theme_transparent()+
  
  theme(axis.text.x = element_text(face="bold", color="black", 
                                   size=15, angle=0),
        axis.title.x = element_blank(),
        axis.ticks = element_blank(),
        axis.text.y = element_text(face="bold", color="black", 
                                   size=15, angle=0),
        axis.title.y = element_blank(),
        plot.title = element_text(face="bold", color="black", 
                                  size=15, angle=0),
        axis.title = element_text(face="bold", color="black", 
                                  size=15, angle=0),
        legend.position="none")+
  stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "t.test",
                     ref.group = ".all.") 




Stat.HzEnd<-ggplot(data = DataIaRaw, 
                   aes(x=as.factor(treatment), 
                       y= Stat.HzEnd, 
                       #colour=interaction(autoClass,as.factor(treatment))))+ 
                       fill=as.factor(treatment)))+
  geom_boxplot(outlier.shape = NA, position = "dodge",varwidth = TRUE)+
  geom_point(alpha = 0.3, size = 3 )+
  scale_fill_manual(values=c(paletteNature2018),"Treatment")+
  
  theme(legend.justification=c(1,1), legend.position=c(1,1), legend.box.just = c("top"), legend.background = element_rect(fill=alpha(0.4)), legend.direction = "horizontal",legend.box="horizontal")+
  xlab("Treatment")+
  ggtitle("Stat.HzEnd")+
  #  ylim(-2,4)+
  #theme(panel.background = element_blank())+
  theme_transparent()+
  
  theme(axis.text.x = element_text(face="bold", color="black", 
                                   size=15, angle=0),
        axis.title.x = element_blank(),
        axis.ticks = element_blank(),
        axis.text.y = element_text(face="bold", color="black", 
                                   size=15, angle=0),
        axis.title.y = element_blank(),
        plot.title = element_text(face="bold", color="black", 
                                  size=15, angle=0),
        axis.title = element_text(face="bold", color="black", 
                                  size=15, angle=0),
        legend.position="none")+
  stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "t.test",
                     ref.group = ".all.") 



Stat.LSpkT<-ggplot(data = DataIaRaw, 
                   aes(x=as.factor(treatment), 
                       y= Stat.LSpkT, 
                       #colour=interaction(autoClass,as.factor(treatment))))+ 
                       fill=as.factor(treatment)))+
  geom_boxplot(outlier.shape = NA, position = "dodge",varwidth = TRUE)+
  geom_point(alpha = 0.3, size = 3 )+
  scale_fill_manual(values=c(paletteNature2018),"Treatment")+
  
  theme(legend.justification=c(1,1), legend.position=c(1,1), legend.box.just = c("top"), legend.background = element_rect(fill=alpha(0.4)), legend.direction = "horizontal",legend.box="horizontal")+
  xlab("Treatment")+
  ggtitle("Stat.LSpkT")+
  #  ylim(-2,4)+
  #theme(panel.background = element_blank())+
  theme_transparent()+
  
  theme(axis.text.x = element_text(face="bold", color="black", 
                                   size=15, angle=0),
        axis.title.x = element_blank(),
        axis.ticks = element_blank(),
        axis.text.y = element_text(face="bold", color="black", 
                                   size=15, angle=0),
        axis.title.y = element_blank(),
        plot.title = element_text(face="bold", color="black", 
                                  size=15, angle=0),
        axis.title = element_text(face="bold", color="black", 
                                  size=15, angle=0),
        legend.position="none")+
  stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "t.test",
                     ref.group = ".all.") 





Stat.LSpkF<-ggplot(data = DataIaRaw, 
                   aes(x=as.factor(treatment), 
                       y= Stat.LSpkF, 
                       #colour=interaction(autoClass,as.factor(treatment))))+ 
                       fill=as.factor(treatment)))+
  geom_boxplot(outlier.shape = NA, position = "dodge",varwidth = TRUE)+
  geom_point(alpha = 0.3, size = 3 )+
  scale_fill_manual(values=c(paletteNature2018),"Treatment")+
  
  theme(legend.justification=c(1,1), legend.position=c(1,1), legend.box.just = c("top"), legend.background = element_rect(fill=alpha(0.4)), legend.direction = "horizontal",legend.box="horizontal")+
  xlab("Treatment")+
  ggtitle("Stat.LSpkF")+
  #  ylim(-2,4)+
  #theme(panel.background = element_blank())+
  theme_transparent()+
  
  theme(axis.text.x = element_text(face="bold", color="black", 
                                   size=15, angle=0),
        axis.title.x = element_blank(),
        axis.ticks = element_blank(),
        axis.text.y = element_text(face="bold", color="black", 
                                   size=15, angle=0),
        axis.title.y = element_blank(),
        plot.title = element_text(face="bold", color="black", 
                                  size=15, angle=0),
        axis.title = element_text(face="bold", color="black", 
                                  size=15, angle=0),
        legend.position="none")+
  stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "t.test",
                     ref.group = ".all.") 



Stat.slp<-ggplot(data = DataIaRaw, 
                 aes(x=as.factor(treatment), 
                     y= Stat.slp, 
                     #colour=interaction(autoClass,as.factor(treatment))))+ 
                     fill=as.factor(treatment)))+
  geom_boxplot(outlier.shape = NA, position = "dodge",varwidth = TRUE)+
  geom_point(alpha = 0.3, size = 3 )+
  scale_fill_manual(values=c(paletteNature2018),"Treatment")+
  
  theme(legend.justification=c(1,1), legend.position=c(1,1), legend.box.just = c("top"), legend.background = element_rect(fill=alpha(0.4)), legend.direction = "horizontal",legend.box="horizontal")+
  xlab("Treatment")+
  ggtitle("Stat.slp")+
  #  ylim(-2,4)+
  #theme(panel.background = element_blank())+
  theme_transparent()+
  
  theme(axis.text.x = element_text(face="bold", color="black", 
                                   size=15, angle=0),
        axis.title.x = element_blank(),
        axis.ticks = element_blank(),
        axis.text.y = element_text(face="bold", color="black", 
                                   size=15, angle=0),
        axis.title.y = element_blank(),
        plot.title = element_text(face="bold", color="black", 
                                  size=15, angle=0),
        axis.title = element_text(face="bold", color="black", 
                                  size=15, angle=0),
        legend.position="none")+
  stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "t.test",
                     ref.group = ".all.") 





Stat.spkNum<-ggplot(data = DataIaRaw, 
                    aes(x=as.factor(treatment), 
                        y= Stat.spkNum, 
                        #colour=interaction(autoClass,as.factor(treatment))))+ 
                        fill=as.factor(treatment)))+
  geom_boxplot(outlier.shape = NA, position = "dodge",varwidth = TRUE)+
  geom_point(alpha = 0.3, size = 3 )+
  scale_fill_manual(values=c(paletteNature2018),"Treatment")+
  
  theme(legend.justification=c(1,1), legend.position=c(1,1), legend.box.just = c("top"), legend.background = element_rect(fill=alpha(0.4)), legend.direction = "horizontal",legend.box="horizontal")+
  xlab("Treatment")+
  ggtitle("Stat.spkNum")+
  #  ylim(-2,4)+
  #theme(panel.background = element_blank())+
  theme_transparent()+
  
  theme(axis.text.x = element_text(face="bold", color="black", 
                                   size=15, angle=0),
        axis.title.x = element_blank(),
        axis.ticks = element_blank(),
        axis.text.y = element_text(face="bold", color="black", 
                                   size=15, angle=0),
        axis.title.y = element_blank(),
        plot.title = element_text(face="bold", color="black", 
                                  size=15, angle=0),
        axis.title = element_text(face="bold", color="black", 
                                  size=15, angle=0),
        legend.position="none")+
  stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "t.test",
                     ref.group = ".all.") 

multiplot(	Stat.mSfr,	Stat.afr,	Stat.HzSTD,	Stat.HzEnd,	Stat.LSpkT,	Stat.LSpkF,	Stat.spkNum,		Stat.slp,  cols=4)
dev.off()







png(file="/Users/nickhousley/Desktop/Projects /Papers/GT_2018_Cancer_Chemo_Interaction_Paper/GT_2018_CancerChemoInteraction_Paper/Figures/Ia Neurons/FeatureGraphsThr.png",res = 600, width = 8000, height =  6000, bg = "transparent")


Thr.T<-ggplot(data = DataIaRaw, 
              aes(x=as.factor(treatment), 
                  y= Thr.T, 
                  #colour=interaction(autoClass,as.factor(treatment))))+ 
                  fill=as.factor(treatment)))+
  geom_boxplot(outlier.shape = NA, position = "dodge",varwidth = TRUE)+
  geom_point(alpha = 0.3, size = 3 )+
  scale_fill_manual(values=c(paletteNature2018),"Treatment")+
  
  theme(legend.justification=c(1,1), legend.position=c(1,1), legend.box.just = c("top"), legend.background = element_rect(fill=alpha(0.4)), legend.direction = "horizontal",legend.box="horizontal")+
  xlab("Treatment")+
  ggtitle("Thr.T")+
  #  ylim(-2,4)+
  #theme(panel.background = element_blank())+
  theme_transparent()+
  
  theme(axis.text.x = element_text(face="bold", color="black", 
                                   size=15, angle=0),
        axis.title.x = element_blank(),
        axis.ticks = element_blank(),
        axis.text.y = element_text(face="bold", color="black", 
                                   size=15, angle=0),
        axis.title.y = element_blank(),
        plot.title = element_text(face="bold", color="black", 
                                  size=15, angle=0),
        axis.title = element_text(face="bold", color="black", 
                                  size=15, angle=0),
        legend.position="none")+
  stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "t.test",
                     ref.group = ".all.") 




Thr.F<-ggplot(data = DataIaRaw, 
              aes(x=as.factor(treatment), 
                  y= Thr.F, 
                  #colour=interaction(autoClass,as.factor(treatment))))+ 
                  fill=as.factor(treatment)))+
  geom_boxplot(outlier.shape = NA, position = "dodge",varwidth = TRUE)+
  geom_point(alpha = 0.3, size = 3 )+
  scale_fill_manual(values=c(paletteNature2018),"Treatment")+
  
  theme(legend.justification=c(1,1), legend.position=c(1,1), legend.box.just = c("top"), legend.background = element_rect(fill=alpha(0.4)), legend.direction = "horizontal",legend.box="horizontal")+
  xlab("Treatment")+
  ggtitle("Thr.F")+
  #  ylim(-2,4)+
  #theme(panel.background = element_blank())+
  theme_transparent()+
  
  theme(axis.text.x = element_text(face="bold", color="black", 
                                   size=15, angle=0),
        axis.title.x = element_blank(),
        axis.ticks = element_blank(),
        axis.text.y = element_text(face="bold", color="black", 
                                   size=15, angle=0),
        axis.title.y = element_blank(),
        plot.title = element_text(face="bold", color="black", 
                                  size=15, angle=0),
        axis.title = element_text(face="bold", color="black", 
                                  size=15, angle=0),
        legend.position="none")+
  stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "t.test",
                     ref.group = ".all.") 



Thr.L<-ggplot(data = DataIaRaw, 
              aes(x=as.factor(treatment), 
                  y= Thr.L, 
                  #colour=interaction(autoClass,as.factor(treatment))))+ 
                  fill=as.factor(treatment)))+
  geom_boxplot(outlier.shape = NA, position = "dodge",varwidth = TRUE)+
  geom_point(alpha = 0.3, size = 3 )+
  scale_fill_manual(values=c(paletteNature2018),"Treatment")+
  
  theme(legend.justification=c(1,1), legend.position=c(1,1), legend.box.just = c("top"), legend.background = element_rect(fill=alpha(0.4)), legend.direction = "horizontal",legend.box="horizontal")+
  xlab("Treatment")+
  ggtitle("Thr.L")+
  #  ylim(-2,4)+
  #theme(panel.background = element_blank())+
  theme_transparent()+
  
  theme(axis.text.x = element_text(face="bold", color="black", 
                                   size=15, angle=0),
        axis.title.x = element_blank(),
        axis.ticks = element_blank(),
        axis.text.y = element_text(face="bold", color="black", 
                                   size=15, angle=0),
        axis.title.y = element_blank(),
        plot.title = element_text(face="bold", color="black", 
                                  size=15, angle=0),
        axis.title = element_text(face="bold", color="black", 
                                  size=15, angle=0),
        legend.position="none")+
  stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "t.test",
                     ref.group = ".all.") 




Thr.T1<-ggplot(data = DataIaRaw, 
               aes(x=as.factor(treatment), 
                   y= Thr.T1, 
                   #colour=interaction(autoClass,as.factor(treatment))))+ 
                   fill=as.factor(treatment)))+
  geom_boxplot(outlier.shape = NA, position = "dodge",varwidth = TRUE)+
  geom_point(alpha = 0.3, size = 3 )+
  scale_fill_manual(values=c(paletteNature2018),"Treatment")+
  
  theme(legend.justification=c(1,1), legend.position=c(1,1), legend.box.just = c("top"), legend.background = element_rect(fill=alpha(0.4)), legend.direction = "horizontal",legend.box="horizontal")+
  xlab("Treatment")+
  ggtitle("Thr.T1")+
  #  ylim(-2,4)+
  #theme(panel.background = element_blank())+
  theme_transparent()+
  
  theme(axis.text.x = element_text(face="bold", color="black", 
                                   size=15, angle=0),
        axis.title.x = element_blank(),
        axis.ticks = element_blank(),
        axis.text.y = element_text(face="bold", color="black", 
                                   size=15, angle=0),
        axis.title.y = element_blank(),
        plot.title = element_text(face="bold", color="black", 
                                  size=15, angle=0),
        axis.title = element_text(face="bold", color="black", 
                                  size=15, angle=0),
        legend.position="none")+
  stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "t.test",
                     ref.group = ".all.") 



Thr.F1<-ggplot(data = DataIaRaw, 
               aes(x=as.factor(treatment), 
                   y= Thr.F1, 
                   #colour=interaction(autoClass,as.factor(treatment))))+ 
                   fill=as.factor(treatment)))+
  geom_boxplot(outlier.shape = NA, position = "dodge",varwidth = TRUE)+
  geom_point(alpha = 0.3, size = 3 )+
  scale_fill_manual(values=c(paletteNature2018),"Treatment")+
  
  theme(legend.justification=c(1,1), legend.position=c(1,1), legend.box.just = c("top"), legend.background = element_rect(fill=alpha(0.4)), legend.direction = "horizontal",legend.box="horizontal")+
  xlab("Treatment")+
  ggtitle("Thr.F1")+
  #  ylim(-2,4)+
  #theme(panel.background = element_blank())+
  theme_transparent()+
  
  theme(axis.text.x = element_text(face="bold", color="black", 
                                   size=15, angle=0),
        axis.title.x = element_blank(),
        axis.ticks = element_blank(),
        axis.text.y = element_text(face="bold", color="black", 
                                   size=15, angle=0),
        axis.title.y = element_blank(),
        plot.title = element_text(face="bold", color="black", 
                                  size=15, angle=0),
        axis.title = element_text(face="bold", color="black", 
                                  size=15, angle=0),
        legend.position="none")+
  stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "t.test",
                     ref.group = ".all.") 



Thr.L1<-ggplot(data = DataIaRaw, 
               aes(x=as.factor(treatment), 
                   y= Thr.L1, 
                   #colour=interaction(autoClass,as.factor(treatment))))+ 
                   fill=as.factor(treatment)))+
  geom_boxplot(outlier.shape = NA, position = "dodge",varwidth = TRUE)+
  geom_point(alpha = 0.3, size = 3 )+
  scale_fill_manual(values=c(paletteNature2018),"Treatment")+
  
  theme(legend.justification=c(1,1), legend.position=c(1,1), legend.box.just = c("top"), legend.background = element_rect(fill=alpha(0.4)), legend.direction = "horizontal",legend.box="horizontal")+
  xlab("Treatment")+
  ggtitle("Thr.L1")+
  #  ylim(-2,4)+
  #theme(panel.background = element_blank())+
  theme_transparent()+
  
  theme(axis.text.x = element_text(face="bold", color="black", 
                                   size=15, angle=0),
        axis.title.x = element_blank(),
        axis.ticks = element_blank(),
        axis.text.y = element_text(face="bold", color="black", 
                                   size=15, angle=0),
        axis.title.y = element_blank(),
        plot.title = element_text(face="bold", color="black", 
                                  size=15, angle=0),
        axis.title = element_text(face="bold", color="black", 
                                  size=15, angle=0),
        legend.position="none")+
  stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "t.test",
                     ref.group = ".all.") 


Thr.T3<-ggplot(data = DataIaRaw, 
               aes(x=as.factor(treatment), 
                   y= Thr.T3, 
                   #colour=interaction(autoClass,as.factor(treatment))))+ 
                   fill=as.factor(treatment)))+
  geom_boxplot(outlier.shape = NA, position = "dodge",varwidth = TRUE)+
  geom_point(alpha = 0.3, size = 3 )+
  scale_fill_manual(values=c(paletteNature2018),"Treatment")+
  
  theme(legend.justification=c(1,1), legend.position=c(1,1), legend.box.just = c("top"), legend.background = element_rect(fill=alpha(0.4)), legend.direction = "horizontal",legend.box="horizontal")+
  xlab("Treatment")+
  ggtitle("Thr.T3")+
  #  ylim(-2,4)+
  #theme(panel.background = element_blank())+
  theme_transparent()+
  
  theme(axis.text.x = element_text(face="bold", color="black", 
                                   size=15, angle=0),
        axis.title.x = element_blank(),
        axis.ticks = element_blank(),
        axis.text.y = element_text(face="bold", color="black", 
                                   size=15, angle=0),
        axis.title.y = element_blank(),
        plot.title = element_text(face="bold", color="black", 
                                  size=15, angle=0),
        axis.title = element_text(face="bold", color="black", 
                                  size=15, angle=0),
        legend.position="none")+
  stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "t.test",
                     ref.group = ".all.") 


Thr.F3<-ggplot(data = DataIaRaw, 
               aes(x=as.factor(treatment), 
                   y= Thr.F3, 
                   #colour=interaction(autoClass,as.factor(treatment))))+ 
                   fill=as.factor(treatment)))+
  geom_boxplot(outlier.shape = NA, position = "dodge",varwidth = TRUE)+
  geom_point(alpha = 0.3, size = 3 )+
  scale_fill_manual(values=c(paletteNature2018),"Treatment")+
  
  theme(legend.justification=c(1,1), legend.position=c(1,1), legend.box.just = c("top"), legend.background = element_rect(fill=alpha(0.4)), legend.direction = "horizontal",legend.box="horizontal")+
  xlab("Treatment")+
  ggtitle("Thr.F3")+
  #  ylim(-2,4)+
  #theme(panel.background = element_blank())+
  theme_transparent()+
  
  theme(axis.text.x = element_text(face="bold", color="black", 
                                   size=15, angle=0),
        axis.title.x = element_blank(),
        axis.ticks = element_blank(),
        axis.text.y = element_text(face="bold", color="black", 
                                   size=15, angle=0),
        axis.title.y = element_blank(),
        plot.title = element_text(face="bold", color="black", 
                                  size=15, angle=0),
        axis.title = element_text(face="bold", color="black", 
                                  size=15, angle=0),
        legend.position="none")+
  stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "t.test",
                     ref.group = ".all.") 


Thr.L3<-ggplot(data = DataIaRaw, 
               aes(x=as.factor(treatment), 
                   y= Thr.L3, 
                   #colour=interaction(autoClass,as.factor(treatment))))+ 
                   fill=as.factor(treatment)))+
  geom_boxplot(outlier.shape = NA, position = "dodge",varwidth = TRUE)+
  geom_point(alpha = 0.3, size = 3 )+
  scale_fill_manual(values=c(paletteNature2018),"Treatment")+
  
  theme(legend.justification=c(1,1), legend.position=c(1,1), legend.box.just = c("top"), legend.background = element_rect(fill=alpha(0.4)), legend.direction = "horizontal",legend.box="horizontal")+
  xlab("Treatment")+
  ggtitle("Thr.L3")+
  #  ylim(-2,4)+
  #theme(panel.background = element_blank())+
  theme_transparent()+
  
  theme(axis.text.x = element_text(face="bold", color="black", 
                                   size=15, angle=0),
        axis.title.x = element_blank(),
        axis.ticks = element_blank(),
        axis.text.y = element_text(face="bold", color="black", 
                                   size=15, angle=0),
        axis.title.y = element_blank(),
        plot.title = element_text(face="bold", color="black", 
                                  size=15, angle=0),
        axis.title = element_text(face="bold", color="black", 
                                  size=15, angle=0),
        legend.position="none")+
  stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "t.test",
                     ref.group = ".all.") 


multiplot(Thr.T,	Thr.F,	Thr.L,	Thr.T1,	Thr.F1,	Thr.L1,	Thr.T3,	Thr.F3,	Thr.L3, cols=5)
dev.off()