library(survival)
library(survminer)
library(cvTools)

## cross validation
data<-read.csv("multigene-171_combined.csv",row.names=1)
two_year_results<-matrix(NA, nrow = 100, ncol = 29)
data$score <- NA
data$group <- NA

times=0
f=0
loop=0

set.seed(99)
for (b in 1:100)
{
  times=times+1
  set.seed(2*b)
  folds <- cvFolds(nrow(data),K=10)
  for(c in 1:10)
 {    
  x <- data[folds$subsets[folds$which != c], ]
  t<-x[,13:33]
  x<-cbind(x[,1:12],t[,colSums(t)>8])

  uni_cox<-list()
  for (i in 13:ncol(x)){
    uni<-coxph(Surv(DFS, DFS_status) ~x[,i]*x[,1],data=x,na.action=na.exclude)
    uni<-as.data.frame(as.matrix(coef(summary(uni))))[3,]
    uni$feature<-colnames(x)[[i]]
    uni_cox[[i]]<-uni
    }  
  
  cox_result<-do.call(rbind,uni_cox)
  cox_result<-na.omit(cox_result)
  tmp<-cox_result[ cox_result[,5]<0.05,c(1,4,5,6)]
  # cox_result<-cbind(cox_result,p.adjust(cox_result[,5],method="fdr"))
  # tmp<-cox_result[cox_result[,7]<0.1 & cox_result[,5]<0.01,c(1,4,5,6)]
  
  if(nrow(tmp)>0){
    loop=loop+1
    for (j in 1:nrow(tmp))
    {
      if(tmp[j,1]<0)
      {flag=T}
      
    }
    
    if(flag==T)
    {   
      f=f+1
      score=0
      for (j in 1:nrow(tmp))
      {
        score=score+data[folds$subsets[folds$which == c],tmp[j,4]]*tmp[j,2]
      }
      data$score[folds$subsets[folds$which == c]]<-score
      
      data$group[data$score < -0.5] <- "Chemo"
      data$group[data$score <=0.5 & data$score>= -0.5] <- "Target1"
      data$group[data$score > 0.5] <- "Target2"
    }
   }
  }
    
      cox_fit<-summary(coxph(Surv(DFS, DFS_status) ~Adj*score,data = data))$coefficients[3,]
      cor1<-summary(coxph(Surv(DFS, DFS_status) ~Adj*score,data = data))$concordance[1]
      cor2<-summary(coxph(Surv(DFS, DFS_status) ~Adj*group,data = data))$concordance[1]
      
      fit <- survfit(Surv(DFS, DFS_status) ~ Adj+group,data = data) 
      #ggsurvplot(fit, data = data, risk.table = TRUE)
      data1<-data[data$group=="Chemo",]
      data2<-data[data$group=="Target1",]
      data3<-data[data$group=="Target2",]
      outcome<-try(summary(fit,time=24),silent=T)
      
      
      if(sum(attr(outcome,"class") %in% "try-error")==0)
      {
        if(length(summary(fit,time=24)$surv)==6)
        {
          s1<-summary(coxph(Surv(DFS, DFS_status) ~Adj,data = data1))$coefficients
          s2<-summary(coxph(Surv(DFS, DFS_status) ~Adj,data = data2))$coefficients
          s3<-summary(coxph(Surv(DFS, DFS_status) ~Adj,data = data3))$coefficients
          
          
          two_year_survival<-c(fit$n,summary(fit,time=24)$surv,surv_median(fit)$median,
                               cox_fit[2],cox_fit[3],cox_fit[5],s1[2],s1[5],s2[2],s2[5],s3[2],s3[5], cor1,cor2)
        } 
        else
        {
          two_year_survival<-c(rep(NA,29))
        }
        
        #mpfs<-surv_median(fit)
        #uni<-as.data.frame(as.matrix(summary(fit,time=24))))[3,]
        #mpfs_results[[i]]<-mpfs
        
      } 
      else {
        two_year_survival<-c(rep(NA,29))
      }
  two_year_results[b,]<-two_year_survival
 }

colnames(two_year_results)<-c("cp_g","tp_g","htp_g","cp_c","tp_c","htp_c",
                              "cp_g_sur","tp_g_sur","htp_g_sur","cp_c_sur","tp_c_suv","htp_c_suv",
                              "cp_g_dfs","tp_g_dfs","htp_g_dfs","cp_c_dfs","tp_c_dfs","htp_c_dfs",
                              "ihr","ihr_se","ihr_p","cp_hr","cp_p","tp_hr","tp_p","htp_hr","htp_p","concordance_score","concordance_grop")

write.csv(file="cross_validation_score_prediction.csv",two_year_results,quote=F)

##
df = data.frame(two_year_results)
# 2-year surv ratio
mean(as.numeric(df$htp_g_sur)/as.numeric(df$htp_c_suv),na.rm=T)
mean(as.numeric(df$tp_g_sur)/as.numeric(df$tp_c_suv),na.rm=T) 
mean(as.numeric(df$cp_g_sur)/as.numeric(df$cp_c_sur),na.rm=T) 
# mDFS diff
mean(as.numeric(df$htp_g_dfs)-as.numeric(df$htp_c_dfs),na.rm=T) 
mean(as.numeric(df$tp_g_dfs)-as.numeric(df$tp_c_dfs),na.rm=T) 
mean(as.numeric(df$cp_g_dfs)-as.numeric(df$cp_c_dfs),na.rm=T) 


