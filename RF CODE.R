#step 1 getGEO
rm(list = ls())
library(GEOquery)
gse ="GSE22255"
Sys.setenv("VROOM_CONNECTION_SIZE"= 2621440*2)
eSet <- getGEO(gse,
               destdir ='.',
               getGPL = F)
exp <- exprs(eSet[[1]])
exp[1:4,1:4]
pd <- pData(eSet[[1]])
p = identical(rownames(pd),colnames(exp))
p
if(!p) exp = exp[,match(rownames(pd),colnames(exp))]
gpl <- eSet[[1]]@annotation
gpl570<-getGEO(gpl,destdir =".")
colnames(Table(gpl570))
fd <- Table(gpl570)[,c(1,11,12)]
getwd()
setwd("C:/Users/Janice/Desktop/学位论文材料/github工作空间/RandomForest")
save(eSet,gse,pd,exp,gpl,fd,file ="getGEO.Rdata")

#step 2 group data
rm(list = ls())
load("getGEO.Rdata")
library(stringr)
library(dplyr)
gr = pd$characteristics_ch1
k1 = str_detect(gr,"IS patient")
table(k1)
k2 = str_detect(gr,"control")
table(k2)
group = ifelse(k1,"disease",'normal')
table(group)
group = factor(group,levels = c("normal","disease"))

#RF data
rm(list = ls())
setwd("C:/Users/Janice/Desktop/学位论文材料/github工作空间/ML&RNA")
load("probe_HMLE.Rdata")
load("probe_LMHE.Rdata")
setwd("C:/Users/Janice/Desktop/学位论文材料/github工作空间/RandomForest")
load("getGEO.Rdata")
load("group.Rdata")
ME<-rbind(HMLE,LMHE)
rownames(ME)<-ME$.id
exp <- exp[match(rownames(ME),rownames(exp)),]
result <- data.frame(ex=rowMeans(exp),
                     symbol=ME$symbol)
result<-data.frame(result,ID=rownames(result))
result <- arrange(result,symbol,desc(ex))
de <- result[!duplicated(result$symbol),]
exp <- exp[match(de$ID,rownames(exp)),]
rownames(exp)<-de$symbol
exp<-as.data.frame(t(exp))
exp<-cbind(group,exp)
save(exp,fd,ME,pd,group,file ="RFdata.Rdata")

#Random Forest
library(randomForest)
library(ggplot2)
library(cowplot)
str(exp)
set.seed(1423)
fit.forest <- randomForest(group~., data=exp,
                           ntree=3000,proximity=TRUE)
fit.forest
oob.error.data <- data.frame(
  Trees=rep(1:nrow(fit.forest$err.rate), times=3),
  Type=rep(c("OOB", "normal", "disease"), each=nrow(fit.forest$err.rate)),
  Error=c(fit.forest$err.rate[,"OOB"], 
          fit.forest$err.rate[,"normal"], 
          fit.forest$err.rate[,"disease"]))

ggplot(data=oob.error.data, aes(x=Trees, y=Error)) +
  geom_line(aes(color=Type))

set.seed(1423)
fit.forest_2 <- randomForest(group~., data=exp,
                             ntree=500,proximity=TRUE)
fit.forest_2


oob.values <- vector(length=37)
set.seed(1423)
for(i in 1:37) {
  fit.forest_3 <- randomForest(group~ ., data=exp, mtry=i, ntree=1000)
  oob.values[i] <- fit.forest_3$err.rate[nrow(fit.forest_3$err.rate),1]
}
oob.values
min(oob.values)
which(oob.values == min(oob.values))

set.seed(1423)
fit.forest_4 <- randomForest(group~ ., 
                             data=exp,
                             ntree=500, 
                             importance=T,
                             proximity=TRUE, 
                             mtry=34)
fit.forest_4
varImpPlot(fit.forest_4)

#MDS plot
distance.matrix <- as.dist(1-fit.forest_4$proximity)
mds.stuff <- cmdscale(distance.matrix, eig=TRUE, x.ret=TRUE)
mds.var.per <- round(mds.stuff$eig/sum(mds.stuff$eig)*100, 1)
mds.values <- mds.stuff$points
mds.data <- data.frame(Sample=rownames(mds.values),
                       X=mds.values[,1],
                       Y=mds.values[,2],
                       Status=group)
ggplot(data=mds.data, aes(x=X, y=Y, label=Sample)) + 
  geom_text(aes(color=Status)) +
  theme_bw() +
  xlab(paste("MDS1 - ", mds.var.per[1], "%", sep="")) +
  ylab(paste("MDS2 - ", mds.var.per[2], "%", sep="")) +
  ggtitle("MDS plot using (1 - Random Forest Proximities)")
ggsave(file="random_forest_mds_plot.pdf")
save(fit.forest,fit.forest_2,fit.forest_3,fit.forest_4,file="RF.Rdata")

