#!/usr/bin/env Rscript

library(docopt)
library(ggplot2)
library(broom)
library(randomForest)
library(cvTools)
library(gridExtra)


'usage: tool.R <input1> <input2> [--out=<output1>] [--out2=<output2>] [--out3=<output3>] [--out4=<output4>] [--out5=<output5>]
tool.R -h | --help

options:
<input1>        Data
<input2>        Vector of policices
--out=<output1>  linear regression selection table [Default: LR.pdf]
--out2=<output2>  random forest selection table [Default: RF.pdf]
--out3=<output3>  performance table  [Default: PT.pdf]
--out4=<output4>  table of quality  [Default: QT.statv]
--out5=<output5>  feature importance  [Default: FI.statv]
' -> doc


args<- docopt(doc)

#----------------reading aaks file----------------------------------------------------------
aaks_read <- function(f)
{
  df <- read.table(f)
  colnames(df) <- c('avgwait', 'maxwait', 'avgbsld', 'maxbsld', 'avgflow', 'maxflow',
                    'name', 'id1', 'id2', 'policy', 'sqff',
                    'AvgRun','MaxRun','MinRun','SumRun', 'TenRun', 'MedRun', 'goRun',
                    'AvgTReq','MaxTReq','MinTReq','SumTReq', 'TenTReq', 'MedTReq', 'goTReq',
                    'AvgSize','MaxSize','MinSize','SumSize', 'TenSize', 'MedSize', 'goSize',
                    'AvgArea','MaxArea','MinArea','SumArea', 'TenArea', 'MedArea', 'goArea',
                    'AvgRR','MaxRR','MinRR','SumRR', 'TenRR', 'MedRR', 'goRR',
                    'AvgRJ','MaxRJ','MinRJ','SumRJ', 'TenRJ', 'MedRJ', 'goRJ',
                    'AvgCos','MaxCos','MinCos','SumCos', 'TenCos', 'MedCos', 'goCos',
                    'AvgSin','MaxSin','MinSin','SumSin', 'TenSin', 'MedSin', 'goSin')
  return(df)
}
#-------------------------------------------------------------------------------------------

df <- aaks_read(file)

df <- aaks_read(args$input1)

inp2<-args$input2

inp2 = "fcfs, lcfs, lpf, spf, sqf, lqf, lexp, sexp, lrf, srf, laf, saf"
policies <- array(unlist(strsplit(inp2, ", ")))

K <- length(policies)


#absolutely useless features
drops <- c("MinRun","MaxTReq","MinTReq","MinSize","TenSize",
           "MinArea","MaxCos","MinCos","MaxSin", "MinSin")
data <- df[,!(names(df) %in% drops)]


#transformation of data (gathering of simulations for each week)----------------------------
for(i in 1:K){
  data[,paste(policies[i])] <- 0
}

maxid1<-max(data$id1)
maxid2<-max(data$id2)
l <- maxid1 * maxid2
w <- length(data)
L <- length(data[,1])
temp <- data
for(k in 1:L){
  if(k%%K!=0){
    i<-(k %/% K) + 1
    j<-w-K+(k%%K)
  } else {
    i<-(k %/% K)
    j<-w
    for (t in 1:(w-K)){
      temp[i,t]<-data[k,t]
    }
  }
  temp[i,j]<-data[k,1]
}
data <- temp[1:l,-(1:11)]
#-------------------------------------------------------------------------------------------

P <- ncol(data) - K




#--------------Performance-Table-----------------------------------
L <- length(data[,1])
BestTime<-apply(data[,(P+1):(P+K)],1,min)
BestTimeInd<-numeric(L)
NumTimeBest<-numeric(K)
for(i in 1:L){
  BestTimeInd[i]<-min(which(data[i,(P+1):(P+K)]== BestTime[i]))
}
for(i in 1:K){
  NumTimeBest[i] <- length(which(BestTimeInd == i))
}
names(NumTimeBest)<-NULL

AverOfPerf<-apply(data[,(P+1):(P+K)],2,mean)
names(AverOfPerf)<-NULL
tab <- data.frame(policies,NumTimeBest,AverOfPerf)
tab$AvgBest <- 0
for(i in 1:K){
  tab$AvgBest[i]<-mean(subset(BestTime, BestTimeInd == i))
}
#-------------------------------------------------------------------

#-----shifting relatively to average--------------------------------
datashift<-data
l<-length(datashift)
avg<-apply(datashift[,((l-K+1):l)],1,mean)
for (i in 0:(K-1)){
  datashift[,l-i]<- datashift[,l-i] - avg
}
data <- datashift
#-------------------------------------------------------------------



#-----------------linear-regression-t-test--------------------------
rss <- numeric(K)
for(i in 1:K){
  fo <- paste(names(data)[P+i], "~", paste(names(data)[(1:P)], collapse=" + "))
  fit <- lm(fo,data)
  rss[i] <- cvFit(fit, data = data, y = data[,(P+i)], K = 5, cost = mspe)$cv
  ti <- tidy(fit)
  ti$pcon <- ti$p.value<0.1
  ti$name <- policies[i]
  summar <- summary(fit)
  if(i==1){
    tibig<-ti[,c(1,5,6,7)]
    bigsummar <- c(summar$r.squared, unname(summar$fstatistic[1]),policies[i])
  } else {
    tibig<-rbind(tibig,ti[,c(1,5,6,7)]) #result table
    bigsummar <- rbind(bigsummar,c(summar$r.squared, unname(summar$fstatistic[1]),policies[i]))
  }
}
#rss #residual sum of squares by cross-validation

#plotting a feature selection table for t-test
tibig$name = factor(tibig$name,levels = policies)
tbl<-ggplot(tibig, aes(term, name))+
  geom_tile(aes(fill = pcon),colour = "black",size=2)+
  labs(title = NULL, x = "Feature", y = "Policy", colour="Is the cofficient\n significant?")+
  scale_fill_manual(values=c("red","green"))+
  guides(fill=guide_legend(title="The cofficient\n is significant"))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.y=element_text(size=25),
        axis.title.x=element_text(size=25),
        legend.text=element_text(size=17),
        legend.title=element_text(size=19),
        title=element_text(size=14),
        axis.text=element_text(size=19),
        plot.title=element_text(size=17))

#-----------------------------------------------------------------



#-----------Random Forest-----------------------------------------------------------------
loss<-numeric(K)
for (i in 1:K){
  fitRF <- randomForest(x=data[,(1:P)],y=data[,(P+i)],importance=TRUE)
  loss[i] <- sum(fitRF$mse)
  mm <- apply(fitRF$importance,2,sum)
  if(i==1){
    resrf <- data.frame(term = names(data)[1:P],
                        sig = ifelse((fitRF$importance[,1]/mm[1])>=0.05, "TRUE", "FALSE"),
                        prob = (fitRF$importance[,1]/mm[1]),
                        policy = policies[i])
    row.names(resrf)<-NULL
  } else {
    resrf <- rbind(resrf,data.frame(term = names(data)[1:P],
                                    sig = ifelse((fitRF$importance[,1]/mm[1])>=0.05, "TRUE", "FALSE"),
                                    prob = (fitRF$importance[,1]/mm[1]),
                                    policy = policies[i]))
    row.names(resrf)<-NULL
  }
}

#loss/nrow(data) # out-of-bad RSS


#plotting a feature selection table for random forest
resrf$policy = factor(resrf$policy,levels = policies)
tblrf<-ggplot(resrf, aes(term, policy))+
  geom_tile(aes(fill = sig),colour = "black",size=2)+
  labs(title =NULL, x = "Feature", y = "Policy", colour="Is the cofficient\n significant?")+
  scale_fill_manual(values=c("red","green"))+
  guides(fill=guide_legend(title="The cofficient\n is significant"))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.y=element_text(size=25),
        axis.title.x=element_text(size=25),
        legend.text=element_text(size=17),
        legend.title=element_text(size=19),
        title=element_text(size=14),
        axis.text=element_text(size=19),
        plot.title=element_text(size=17))

#-------------------------------------------------------------  





#----------------Output-of-results--------------------------------------------------------

#Output 1 - t-test feature selection
pdf(file=args$out,width=40,height=14)
tbl

#Output 2 - random forest feature selection
pdf(file=args$out2,width=40,height=14)
tblrf

#Output 3 - performance table from conclusion section
tab<-tab[K:1,]
row.names(tab)<-1:K
ggsave(file=args$out3,plot=tableGrob(tab),width=5,height=7)

#Output 4 - table of quality for linear regression 
write.table(bigsummar,file=args$out4)
for (i in 1:nrow(bigsummar)){
  cat(unlist(bigsummar[i,]))
  cat("\n")
}

#Output 5 - table of random forest feature importances
write.table(resrf,file=args$out5)
for (i in 1:nrow(resrf)){
  cat(unlist(resrf[i,]))
  cat("\n")
}
