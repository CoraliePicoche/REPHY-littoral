#This script gets eigen values from SSA for different sites and compare them (this is completely exploratory) (and this completely doesn't work for now)

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

sp="CHA"

#for (l in 1:length(lieu)){
for (l in 1:1){
        liste_sp=tab_sp[!is.na(tab_sp[,l]),l]
        tab=read.table(paste("data/corres_hernandez_",lieu[l],'.txt',sep=''),sep=";",na="NA",header=TRUE)
        dates=as.Date(tab$Date)
        dates_bis=seq(dates[1],dates[length(dates)],timestep) #Regular time grid
        y1=0
        y2=log10(max(tab[,2:dim(tab)[2]],na.rm=TRUE))
	
                rec1=na.approx(tab[,sp],maxgap=consecutif,x=dates,xout=dates_bis,na.rm=FALSE)
                objet_SSA=gapfillSSA(series=rec1,amnt.iters=c(20,20),fill.margins=TRUE,SSA.methods='svd',plot.results=TRUE) #20 is enough for AST to converge for most methods
#                plotGapfillCube(gapfillNcdf(series=rec1,amnt.iters=c(20,20),fill.margins=TRUE,SSA.methods='svd')) #Doesn't work for now
                rec2=objet_SSA$filles.series
		rec2[rec2<=0]=runif(sum(rec2<=0),0,min(tab[,sp],na.rm=TRUE))
}

