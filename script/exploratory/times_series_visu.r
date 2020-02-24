#Plot time series for old files, using Arcachon classification and handling of missing values

graphics.off()
rm(list=ls())
library("zoo")

#sites=c("Eperon","Cornard","Men er Roue","At so","Antoine","Loscolo","Bois de","Croisic","Auger","Boyard","Lazaret","Dunkerque")
#fil=c()
#for (s in sites){
#	fil=c(fil,list.files(path=".",pattern=paste(s,"*.flortot_only.csv")))
#}
fil=list.files(path="data/",pattern="^[A-Z].*flortot_only.csv")
sp=c("AST","NIT","PSE","SKE","CHA","GUI","LEP","RHI","GYM","PRP","CRY","EUG","NEI")
pen=c("AST","NIT","PSE","SKE")
cen=c("CHA","GUI","LEP","RHI")
din=c("GYM","PRP")

consecutif=2 #Number of missing values above which we keep the NA
timestep=14 #Regular time lapses between two observations

set.seed(42)
for (f in 1:length(fil)){
	sisi=strsplit(strsplit(fil[f],split="_")[[1]][1],split="-")[[1]][1]
	tab_flore=read.table(paste("data/",fil[f],sep=""),header=TRUE,sep=";")
	dates=as.Date(tab_flore$Date)
	dates_bis=seq(dates[1],dates[length(dates)],timestep) #Regular time grid
	tab_flore_reconstruct=matrix(NA,nrow=length(dates_bis),ncol=length(sp))
	colnames(tab_flore_reconstruct)=sp
	pdf(paste("./graphe/",sisi,"_per_species.pdf",sep=""),width=15)
	par(mfrow=c(2,2))
	for (s in sp){
		if(sum(!is.na(tab_flore[,s]))>0){
			tab_flore_reconstruct[,s]=na.approx(tab_flore[,s],maxgap=consecutif,x=dates,xout=dates_bis,na.rm=FALSE) #Interpolation over regular time grid
	       		tab_flore_reconstruct[is.na(tab_flore_reconstruct[,s]),s]=runif(sum(is.na(tab_flore_reconstruct[,s])),0,min(tab_flore_reconstruct[,s],na.rm=TRUE))
		}
		plot(dates_bis,log10(tab_flore_reconstruct[,s]),t="l",col=c("black"),main=s,lwd=2,ylim=c(0,7))
		points(dates,log10(tab_flore[,s]),pch=16,col=c("red"))
		if(s=="AST"){
			legend("topleft",c("reconstruct","real"),lty=c(1,NA),pch=c(NA,16),col=c("black","red"),bty="n")
		}
	}

	tab_pen_reconstruct=apply(tab_flore_reconstruct[,pen],1,sum,na.rm=TRUE)
	tab_pen=apply(tab_flore[,pen],1,sum,na.rm=TRUE)
        plot(dates_bis,log10(tab_pen_reconstruct),t="l",col=c("black"),main="PEN",lwd=2,ylim=c(0,7))
        points(dates,log10(tab_pen),pch=16,col=c("red"))	

	tab_cen_reconstruct=apply(tab_flore_reconstruct[,cen],1,sum,na.rm=TRUE)
	tab_cen=apply(tab_flore[,cen],1,sum,na.rm=TRUE)
        plot(dates_bis,log10(tab_cen_reconstruct),t="l",col=c("black"),main="CEN",lwd=2,ylim=c(0,7))
        points(dates,log10(tab_cen),pch=16,col=c("red"))	
	
	tab_din_reconstruct=apply(tab_flore_reconstruct[,din],1,sum,na.rm=TRUE)
	tab_din=apply(tab_flore[,din],1,sum,na.rm=TRUE)
        plot(dates_bis,log10(tab_din_reconstruct),t="l",col=c("black"),main="DIN",lwd=2,ylim=c(0,7))
        points(dates,log10(tab_din),pch=16,col=c("red"))	

	dev.off()
}
