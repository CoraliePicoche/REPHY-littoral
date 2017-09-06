#################################################
## 31/08 I want to compare the different time series in different places
#################################################

rm(list=ls())
graphics.off()
library('lubridate')

acex=2
alwd=2
colo=c("black","red","blue")

tab_sp=read.table('data/lieu_sp_post_reconstruct_pour_MAR.csv',header=TRUE,na.strings="",sep=";")
lieu=colnames(tab_sp)
lieu=gsub('.',' ',lieu,fixed=TRUE) #useful for Men er Roue
colnames(tab_sp)=lieu

#lieu=c("LEperon","Cornard","Auger")
#lieu=c("Men er Roue","Loscolo","Croisic")
#lieu=c("Antoine","Lazaret")
option_lieu=c("Men er Roue","Loscolo","Croisic","LEperon","Cornard","Auger","Antoine","Lazaret")
pdf("Rapport/graphe/covar_temp_sali.pdf",width=18)
par(mar=c(2,4,3,5))
for (l in 1:length(lieu)){
        tab_cov=read.table(paste("data/",lieu[l],'hydro.txt',sep=''),sep=";",na="NA",header=TRUE)
	a=cor(tab_cov[,'TEMP'],tab_cov[,'SALI'],method="spearman",use="complete.obs")
	plot(as.Date(tab_cov$Date),tab_cov[,'TEMP'],t="o",col="red",pch=16,ylab="TEMP",main=paste(lieu[l],format(a,digits=2),sep=" "))
	par(new=T)
	plot(as.Date(tab_cov$Date),tab_cov[,'SALI'],t="o",col="blue",pch=16,axes=F,ylab="")
	axis(4,col="blue")
	mtext(side=4,line=3,'SALI',col="blue")
}
dev.off()
