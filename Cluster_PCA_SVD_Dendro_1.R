library(rafalib)
library(ISLR)
library(scatterplot3d)
library(ggplot2)
library(PMA)
library(genefilter)
library(RColorBrewer)
library(RColorBrewer)
library(gplots)
library(colorspace)



paletteNature2018 = c("#808080","#3e5ce0","#db0415","#933bd3")

### get colors hexs
## https://www.colorhexa.com/810081 

data_all<-read.csv(file = ("Desktop/ExpMatrix_1_for_PCA_Clustering.csv"), header=T, row.names=1)
data_all<-read.csv(file = ("ExpMatrix_1.csv"), header=T, row.names=1)
data_all<-data_all[,-1]
##rownames(data_all2)<-data_all[,1]
head(data_all)
data_all<-ExpMatrix_1

data_all<-ExpMatrix_1[,c(3:15)]

colnames(data_all)   
colnames(data_all)[1:3]<-"WT"
colnames(data_all)[4:6]<-"WT+OX"
colnames(data_all)[7:9]<-"Pirc"
colnames(data_all)[10:13]<-"Pirc+OX"
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

# plot PC1 vs PC2
plot(Z[,1], Z[,2], type ="n", xlab="PC1", ylab="PC2")
text(Z[,1], Z[,2], colnames(X), col=cols)

#### plot PC2 vs PC3
plot(Z[,2], Z[,3], type ="n", xlab = "PC2", ylab="PC3")
text(Z[,2], Z[,3], colnames(X), col=cols)

#### alternative pretty plot
pc_dat<- data.frame(type = rownames(Z), PC1 = Z[,1], PC2= Z[,2], PC3=Z[,3])
sp<-ggplot(pc_dat,aes(x=PC1, y=PC2, col=type)) + geom_point() + geom_text(aes(label = type), hjust=0, vjust=0)
sp + scale_color_manual(values=c("#db0415", "#933bd3", "#808080","#3e5ce0")) + theme(panel.background = element_rect(fill = 'white', colour = 'black'))

#### 3D
#### http://www.sthda.com/english/wiki/amazing-interactive-3d-scatter-plots-r-software-and-data-visualization
#png(file="/Users/nickhousley/Desktop/Projects /Papers/GT_2018_CancerChemoInteraction_Paper_1/GT_2018_CancerChemoInteraction_Paper_1/Figures/Ia Neurons/transcriptome/3dPCAPlotWT_POX_noTicks.png",res = 600, width = 6000, height =  6000 ,bg = "transparent")
colors<-c("#db0415", "#933bd3", "#808080","#3e5ce0")
colors<-colors[as.numeric(cols)]
S3D<- scatterplot3d(Z[,1],Z[,2],Z[,3], color = colors,pch =16, angle = 45,xlab = "PC1", ylab = "PC2", zlab = "PC3", grid = FALSE, box = FALSE, tick.marks = T, cex.symbols = 3)
#dev.off()
legend("bottomright", legend=c("Pirc", "Pirc+OX", "WT", "WT+OX"), 
       col=c("#db0415", "#933bd3", "#808080","#3e5ce0"), 
       pch =16)



## get a gradient of colors for grey, green, red.
aa<- grep("grey",colors())
bb<- grep("blue",colors())
cc<-  grep("red",colors())
gcol2<- colors()[c(aa[1:30],bb[1:20],rep(cc,2))]


## use the genes that drive the first PC1. This is the first major pattern in the data
## The matrix V contains the weigths for the features, and we can use V to select important features(genes) that contribute to the each PC.
k=1
ord1<- order(abs(V[,k]),decreasing=TRUE)
x1<- as.matrix(X[ord1[1:249],])
h<-heatmap(x1,col=gcol2,labRow = NULL, Rowv=TRUE)
###get rownames (genes/features)
rownames(x1)[h$rowInd]

# use the genes that drive the second PC (PC2)
j<- 1
ord<- order(abs(V[,j]),decreasing=TRUE)

## we just use the first 250 features(genes) to plot a heatmap, This is the second major pattern.
x<- as.matrix(X[ord[1:249],])
heatmap(x,col=gcol2)


###Variance Explained
varex = 0
cumvar = 0
denom = sum(D^2)
for(i in 1:15){
  varex[i] = D[i]^2/denom
  cumvar[i] = sum(D[1:i]^2)/denom
}

## variance explained by each PC cumulatively
cumvar

####screeplot
png(file="/Users/nickhousley/Desktop/Projects /Papers/GT_2018_CancerChemoInteraction_Paper_1/GT_2018_CancerChemoInteraction_Paper_1/Figures/Ia Neurons/transcriptome/CumVarExplainedPCA_trans.png",res = 600, width = 6000, height =  6000 ,bg = "transparent")

par(mfrow=c(1,2))
plot(1:15,varex,type="l",lwd=2,xlab="PC",ylab="% Variance Explained")
plot(1:15,cumvar,type="l",lwd=2,xlab="PC",ylab="Cummulative Variance Explained")

dev.off()

########### SParse PCA ########### ########### ########### ########### 
## we also look at the first 4 PCs
spc = SPC(t(X),sumabsv=10,K=4)

### how many genes do we need to consider
Genes_Need_to_Consider<-apply(spc$v!=0, 2, sum)
genes_need_to_consider_labs = c("PC1","PC2","PC3","PC4")


png(file="/Users/nickhousley/Desktop/Projects /Papers/GT_2018_CancerChemoInteraction_Paper_1/GT_2018_CancerChemoInteraction_Paper_1/Figures/Ia Neurons/transcriptome/Genes_Need_to_Consider.png",res = 600, width = 6000, height =  6000 ,bg = "transparent")
grid.table(Genes_Need_to_Consider,genes_need_to_consider_labs )
dev.off()

# sparce PC scatterplots
png(file="/Users/nickhousley/Desktop/Projects /Papers/GT_2018_CancerChemoInteraction_Paper_1/GT_2018_CancerChemoInteraction_Paper_1/Figures/Ia Neurons/transcriptome/Top_4D_sparcePCA_Plots_trans.png",res = 600, width = 6000, height =  6000 ,bg = "transparent")

cols = as.numeric(as.factor(colnames(data_all)))
K = 3
pclabs = c("SPC1","SPC2","SPC3","SPC4")
par(mfrow=c(1,K))
for(i in 1:K){
  j = i+1
  plot(spc$u[,i],spc$u[,j],type="n",xlab=pclabs[i], ylab=pclabs[j])
  text(spc$u[,i],spc$u[,j],colnames(X),col=cols)
}

dev.off()

#SPC loadings - visualize data by limiting to genes selected by the sparse PC loadings
aa = grep("grey",colors())
bb = grep("green",colors())
cc = grep("red",colors())
gcol2 = colors()[c(aa[1:30],bb[1:20],rep(cc,2))]

j = 1
ind = which(spc$v[,j]!=0)
x = as.matrix(X[ind,])
heatmap(x,col=gcol2)
length(ind)
## Variance Explained
spc$prop.var.explained

####### clustering

###COMPLETE LINKAGE
cols = as.numeric(as.factor(colnames(data_all)))
Dmat = dist(t(data_all))

png(file="/Users/nickhousley/Desktop/Projects /Papers/GT_2018_CancerChemoInteraction_Paper_1/GT_2018_CancerChemoInteraction_Paper_1/Figures/Ia Neurons/transcriptome/ClusterTree_Ward_trans.png",res = 600, width = 6000, height =  6000 ,bg = "transparent")
com.hclust = hclust(Dmat,method="ward.D")

myplclust(com.hclust, labels=colnames(data_all), lab.col=c("#808080","#808080","#808080","#3e5ce0","#3e5ce0","#3e5ce0","#db0415","#db0415","#db0415","#933bd3","#933bd3","#933bd3","#933bd3"), main = "ward linkage eucledian distance")
abline(h=60)
dev.off()
cl<- cutree(com.hclust, h= 60)
table(type=colnames(X), clusters=cl)



#PC loadings - visualize data by limiting to top genes in magnitude in the PC loadings
## get some colors

## use feature weigths for the first PC (PC1)
j = 1
ord = order(abs(V[,j]),decreasing=TRUE)
# x = as.matrix(X[ord[1:249],])
x = as.matrix(X[ord[1:2500],])


# the default is eucledian distance and complete linage for both rows and columns
hm<-heatmap(x,col=hmcols,hclustfun=function(x)hclust(x,method="ward.D"), labRow = FALSE)



rv<- rowVars(X)
idx<- order(-rv)[1:2500]
heatmap(X[idx,], col=hmcols,hclustfun=function(x)hclust(x,method="ward.D"))

hmcols<- colorRampPalette(rev(brewer.pal(9,"RdBu")))(100)
heatmap.2(x,col=hmcols,hclustfun=function(x)hclust(x,method="ward.D"),trace="none", scale = "row")


png(file="/Users/nickhousley/Desktop/Projects /Papers/GT_2018_CancerChemoInteraction_Paper_1/GT_2018_CancerChemoInteraction_Paper_1/Figures/Ia Neurons/transcriptome/HierarchicalCLusterPC1_POX_2019_redblue.png",res = 600, width = 6000, height =  6000 ,bg = "transparent")

hm<-heatmap.2(x,col=rev(redblue(250)),hclustfun=function(x)hclust(x,method="ward.D"),trace="none", scale = "row", labRow = FALSE)
dev.off()

## colors avaiable in the package 
display.brewer.all()
cols1<- palette(brewer.pal(8, "Dark2"))
cols2<- palette(brewer.pal(12, "Paired"))

cols<-c("#808080","#808080","#808080","#3e5ce0","#3e5ce0","#3e5ce0","#db0415","#db0415","#db0415","#933bd3","#933bd3","#933bd3","#933bd3")

cbind(colnames(x), cols)  # check which color maps to different cancer types

hm<- heatmap.2(x, distfun=function(x) as.dist(1-cor(t(x))), hclustfun=function(x)hclust(x,method="ward.D"),trace="column", scale = "row", ColSideColors=cols, col=hmcols, labCol=colnames(data_all), margins = c(6,6), main = "scaled gene and correlation distance",)
names(hm)

### returns a matrix after the clustering (should be in the correct vertical order but opposite, scaled and demeande)
m.afterclust<- x[rev(hm$rowInd),rev(hm$colInd)]



#Separating clusters
#convert the rowDendrogram to a hclust object
hc.rows<- as.hclust(hm$rowDendrogram)  ## needs to be hm derived from heatmap.2 object
hc.cols<- as.hclust(hm$colDendrogram)

table(type=colnames(X), clusters=cutree(hc.cols, k=5))
names(hc.rows)

plot(hc.rows)  # rotate the dendrogram 90 degree, it is the same as in the heatmap

rect.hclust(hc.rows,h=100)
ct<- cutree(hc.rows,h=100)


#get the members' names of each clusters


############# ct holds all the cluster identities ####### very important
write.csv(ct, file = "PC1 Clusters_final.csv")


# get the matrix after clustering in the order of the heatmap (up--->down)

tableclustn<-  data.frame(m.afterclust, cluster = rev(ct[hc.rows$order]))

# remake the heatmap adding the RowSide bar based on the subgroups

mycolhc<- palette(brewer.pal(8, "Dark2"))
mycolhc<-mycolhc[as.vector(ct)]


# remake the heatmap adding the RowSide bar based on the subgroups




heatmap.2(x, distfun=function(x) as.dist(1-cor(t(x))), hclustfun=function(x)hclust(x,method="ward.D"),trace="none", scale = "row", ColSideColors=cols, RowSideColors = mycolhc, col=hmcols, labCol=colnames(data_all), labRow = F, margins = c(6,6), density.info="none", main = "scaled gene and correlation distance")

hmcols<- colorRampPalette(rev(brewer.pal(9,"RdYlBu")))(100)

hmcols<- colorRampPalette<-pal(diverge_hcl(100))
hmcols <- colorRampPalette(brewer.pal(9, "GnBu"))(100)


png(file="/Users/nickhousley/Desktop/Projects /Papers/GT_2018_CancerChemoInteraction_Paper_1/GT_2018_CancerChemoInteraction_Paper_1/Figures/Ia Neurons/transcriptome/HierarchicalCLusterPC1_POX_2019.png",res = 600, width = 6000, height =  6000 ,bg = "transparent")


heatmap.2(x, 
          distfun=function(x) as.dist(1-cor(t(x))), 
          hclustfun=function(x)hclust(x,method="ward.D"),
          trace="none", 
          scale = "row", 
          ColSideColors=cols, 
          RowSideColors = mycolhc,
          col=diverge_hsv(100, s = 1, v = 1, power = 1, gamma = NULL, fixup = TRUE, alpha = 1),
          # col=hmcols,
          labCol=colnames(data_all), 
          margins = c(6,6), main = "", 
          labRow = F)
dev.off()


# assign the output of heatmap.2 to a variable hm
hm<- heatmap.2((as.matrix(data_all)), distfun=function(x) as.dist(1-cor(t(x))), hclustfun=function(x)hclust(x,method="ward.D"),trace="column", scale = "row", ColSideColors=cols, col=hmcols, labCol=colnames(data_all), margins = c(6,6), main = "scaled gene and correlation distance",)
names(hm)
