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
lieu=c("Antoine","Lazaret")
#option_lieu=c("Men er Roue","Loscolo","Croisic","LEperon","Cornard","Auger","Antoine","Lazaret")

a=table(as.matrix(tab_sp[,lieu]))
liste_sp=dimnames(a[a==length(lieu)])[[1]]

cov3_tot=c("TEMP","SALI","CHLOROA")

xmin=as.Date("1995-01-01")
xmax=as.Date("2018-01-01")

if(1==0){
pdf("Rapport/graphe/compare_plankton_species_SU.pdf",width=22)
        #Biotic variables
for (s in liste_sp){
	for (l in 1:length(lieu)){
        tab=read.table(paste("data/corres_hernandez_",lieu[l],'.txt',sep=''),sep=";",na="NA",header=TRUE)
        dates=as.Date(tab$Date)
        tab=tab[year(dates)>=1996,]#Using data from 1996
        dates=dates[year(dates)>=1996]
	if(l==1){
	plot(dates,log10(tab[,s]),col=colo[l],t="o",pch=16,cex=acex,lwd=alwd,main=s,ylab="",cex.main=acex,xlim=c(xmin,xmax))
	}else{
	lines(dates,log10(tab[,s]),col=colo[l],t="o",pch=16,cex=acex,lwd=alwd)
	}
        }
	legend("topright",lieu,col=colo,pch=16,lty=1,cex=2)

}
dev.off()
        #Hydro variables
}

pdf("Rapport/graphe/compare_hydrology_SU.pdf",width=22)
for (c in cov3_tot){
	for (l in 1:length(lieu)){
        tab_cov=read.table(paste("data/",lieu[l],'hydro.txt',sep=''),sep=";",na="NA",header=TRUE)
        dates_cov=as.Date(tab_cov$Date)
	if(l==1){
	if(c=="CHLOROA"){
	plot(dates_cov,tab_cov[,c],col=colo[l],t="o",pch=16,cex=acex,lwd=alwd,main=c,ylab="",xlim=c(xmin,xmax),ylim=c(min(tab_cov[,c]*0.5,na.rm=TRUE),max(tab_cov[,c]*2.,na.rm=TRUE)))
	}else{
	plot(dates_cov,tab_cov[,c],col=colo[l],t="o",pch=16,cex=acex,lwd=alwd,main=c,ylab="",xlim=c(xmin,xmax),ylim=c(min(tab_cov[,c]*0.5,na.rm=TRUE),max(tab_cov[,c]*1.25,na.rm=TRUE)))
	}
	}else{
	lines(dates_cov,tab_cov[,c],col=colo[l],t="o",pch=16,cex=acex,lwd=alwd)
	}
}
	pos="bottomright"
	if(c=="CHLOROA"){
		pos="topright"
	}
	legend(pos,lieu,col=colo,pch=16,lty=1,cex=2)
}

dev.off()
