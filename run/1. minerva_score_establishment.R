library(survival)
data<-read.csv("multigene-171_combined.csv",row.names=1)

#coxph(Surv(DFS,DFS_status) ~ Adj, data= data)
#coxph(Surv(OS,OS_status) ~ Adj, data= data)

#### gene-by-treatment interaction test
uni_cox <-list()
for (i in 12:ncol(data)){
  uni<-coxph(Surv(DFS,DFS_status)~data[,i]*data[,1],data=data, na.action=na.exclude)
  uni<-as.data.frame(as.matrix(coef(summary(uni))))
  uni$feature<-colnames(data)[[i]]
  uni_cox[[i]]<-uni
}
cox_result <- do.call(rbind,uni_cox)

idx=c()
for(i in 1:dim(cox_result)[1]){
  if(i%%3==0){idx=c(idx,i)} 
}
new_cox_result=cox_result[idx,]
mark <- new_cox_result[which(as.numeric(as.character(new_cox_result[,5]))<=0.05), ]
write.csv(mark, file="output/cox.interaction.filter0.05.csv",row.names = F,quote = F)

#### minerva score
mnvgenes = data[,names(data) %in% mark$feature]
score = c()
for(i in 1:171){
  score[i] = sum(as.numeric(mnvgenes[i, mark$feature])*mark$z)
}
mnvgenes$score = score
write.csv(mnvgenes, file="output/minerva_scores.csv", row.names = T,  quote = F)




#### prognostic biomarkers
prog_cox <-list()
for (i in 12:ncol(data)){
  prog_uni<-coxph(Surv(DFS,DFS_status)~data[,i]+data[,1],data=data, na.action=na.exclude)
  prog_uni<-as.data.frame(as.matrix(coef(summary(prog_uni))))
  prog_uni$feature<-colnames(data)[[i]]
  prog_cox[[i]]<-prog_uni
}
pcox_result <- do.call(rbind,prog_cox)

prog_idx=c()
for(i in 1:dim(pcox_result)[1]){
  if(i%%2==0){prog_idx=c(prog_idx,i)} 
}
new_pcox_result=pcox_result[-prog_idx,]
prog_mark <- new_pcox_result[which(as.numeric(as.character(new_pcox_result[,5]))<=0.05), ]
write.csv(prog_mark, file="cox.prognostic.filter0.05.csv", row.names = F, quote = F)


