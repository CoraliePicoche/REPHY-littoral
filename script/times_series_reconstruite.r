rm(list=ls())
graphics.off()
library('stringr')
library('zoo')
library("spectral.methods")
set.seed(42)

consecutif=2
timestep=14

acex=1.75
alwd=2
aft_size=2
amain=2

tab_sp=read.table('data/lieu_sp.csv',header=TRUE,na.strings="",sep=";")
lieu=colnames(tab_sp)
lieu=gsub('.',' ',lieu,fixed=TRUE)

for (l in 1:length(lieu)){
#for (l in 1:1){
	pdf(paste('Rapport/graphe/comparaison_SSA_random_',lieu[l],'.pdf',sep=""),width=15)
	liste_sp=tab_sp[!is.na(tab_sp[,l]),l]
	tab=read.table(paste("data/corres_hernandez_",lieu[l],'.txt',sep=''),sep=";",na="NA",header=TRUE)
        dates=as.Date(tab$Date)
        dates_bis=seq(dates[1],dates[length(dates)],timestep) #Regular time grid
        y1=0
        y2=log10(max(tab[,2:dim(tab)[2]],na.rm=TRUE))
	for (s in 1:length(liste_sp)){
#	for (s in 1:1){
		sure=!(grepl('?',liste_sp[s],fixed=TRUE))
		sp=str_trim(gsub('?','',liste_sp[s],fixed=TRUE))
		if(!sure){
			linestyle=2
		}else{
			linestyle=1
		}
	
		rec1=na.approx(tab[,sp],maxgap=consecutif,x=dates,xout=dates_bis,na.rm=FALSE)
		rec2=gapfillSSA(series=rec1,amnt.iters=c(20,20),fill.margins=TRUE,SSA.methods='svd')$filled.series #20 is enough for AST to converge for most methods
		rec2[rec2<=0]=runif(sum(rec2<=0),0,min(tab[,sp],na.rm=TRUE))
		if(s<3){
			rec3=gapfillSSA(series=rec1,amnt.iters=c(20,20),fill.margins=TRUE,SSA.methods='propack')$filled.series
			rec3[rec3<=0]=runif(sum(rec3<=0),0,min(tab[,sp],na.rm=TRUE))
			rec4=gapfillSSA(series=rec1,amnt.iters=c(20,20),fill.margins=TRUE,SSA.methods='nutrlan')$filled.series
			rec4[rec4<=0]=runif(sum(rec4<=0),0,min(tab[,sp],na.rm=TRUE))
			rec5=gapfillSSA(series=rec1,amnt.iters=c(20,20),fill.margins=TRUE,SSA.methods='eigen')$filled.series
			rec5[rec5<=0]=runif(sum(rec5<=0),0,min(tab[,sp],na.rm=TRUE))
			par(mfrow=c(2,1),mar=c(1,2,1,0.5))
		}else{
			par(mfrow=c(1,1),oma=c(1,2,1,0.5))
		}
		rec1[is.na(rec1)]=runif(sum(is.na(rec1)),0,min(tab[,sp],na.rm=TRUE))

		plot(dates,log10(tab[,sp]),col='red',t="p",ylim=c(y1,y2),pch=16,cex=acex,lwd=alwd,ylab="Log10(Abundance)",cex.main=amain,main=paste(lieu[l],sp))
		lines(dates_bis,log10(rec1),lty=linestyle,col="black",lwd=alwd)	
		lines(dates_bis,log10(rec2),lty=linestyle,col="blue",lwd=alwd)
		legend("topleft",c('Random','SSA svd'),col=c('black','blue'),lty=linestyle)	
		if(s<3){
			plot(dates,log10(tab[,sp]),col='red',t="p",ylim=c(y1,y2),pch=16,cex=acex,lwd=alwd,ylab="Log10(Abundance)")
			lines(dates_bis,log10(rec2),lty=linestyle,col="blue",lwd=alwd)	
			lines(dates_bis,log10(rec3),lty=linestyle,col="black",lwd=alwd)	
			lines(dates_bis,log10(rec4),lty=linestyle,col="orange",lwd=alwd)	
			lines(dates_bis,log10(rec5),lty=linestyle,col="green",lwd=alwd)
			legend("topleft",c('svd','propack','nutrlan','eigen'),col=c('blue','black','orange','green'),lty=linestyle)	
		}
	}
	dev.off()
}
