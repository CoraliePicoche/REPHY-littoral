graphics.off()
rm(list=ls())

library('lubridate')
library('zoo')
library('MARSS')
source("/home/cpicoche/Documents/Plankton/script/MARSS_clean.r")

set.seed(42)
timestep=14
consecutif=2
which_model="unconstrained" #for now, can be "null", "unconstrained" and "pencen"
which_NEI="cov" #NEI can be either a population or a covariate
which_timestep="regular" #for now, only regular is implemented, but I do think we should do monthly average in order to prepare an inter-site comparison
which_sp="common" #can be single (species depend on each site) or "common" (species are the ones in all sites)

tab_sp=read.table('data/lieu_sp_post_reconstruct_pour_MAR.csv',header=TRUE,na.strings="",sep=";")
lieu=colnames(tab_sp)
lieu=gsub('.',' ',lieu,fixed=TRUE) #useful for Men er Roue

cov3_tot=c("TEMP","SALI")

for (l in 1:length(lieu)){
	#Biotic variables
	tab=read.table(paste("data/corres_hernandez_",lieu[l],'.txt',sep=''),sep=";",na="NA",header=TRUE)
        dates=as.Date(tab$Date)
        tab=tab[year(dates)>=1996,]#Using data from 1996
	dates=dates[year(dates)>=1996]

	dates_bis=seq(dates[1],dates[length(dates)],timestep) #Regular time grid

	if(which_sp=="single"){
		liste_sp=as.character(tab_sp[!is.na(tab_sp[,l]),l])
	}else{
		a=table(as.matrix(tab_sp))
		liste_sp=dimnames(a[a==dim(tab_sp)[2]])[[1]]
	}
	
	if(which_NEI=="cov"){
			liste_sp=liste_sp[-which(liste_sp=="NEI")]
	}
	tab_plankton=na.approx(tab[,liste_sp],maxgap=consecutif,x=dates,xout=dates_bis,na.rm=FALSE) #Interpolation over regular time grid
	#Replace missing values
	for (s in liste_sp){
        	tab_plankton[is.na(tab_plankton[,s]),s]=runif(sum(is.na(tab_plankton[,s])),0,min(tab_plankton[,s],na.rm=TRUE))
	}

	#Hydro variables
	tab_cov=read.table(paste("data/",lieu[l],'hydro.txt',sep=''),sep=";",na="NA",header=TRUE)
        dates_cov=as.Date(tab_cov$Date)
	tab_cov_bis=matrix(NA,length(dates_bis),length(cov3_tot))
	colnames(tab_cov_bis)=cov3_tot
	for (c in cov3_tot){
        	tab_cov_bis[,c]=approx(tab_cov[,c],x=dates_cov,xout=dates_bis)$y
	}
	if(which_NEI=="cov"){
		cov3_tot_bis=c(cov3_tot,"NEI")
		tab_cov_bis=cbind(tab_cov_bis,approx(tab[,"NEI"],x=dates,xout=dates_bis)$y)
	}else{
		cov3_tot_bis=cov3_tot
	}

	#Log transfo for species abundance and scaling for all time series
	tab_plankton=log(tab_plankton)
	tab_plankton=t(scale(tab_plankton[2:(length(dates_bis)-1),]))
	tab_cov=t(scale(tab_cov_bis[2:(length(dates_bis)-1),]))
	rownames(tab_plankton)=liste_sp
	rownames(tab_cov)=cov3_tot_bis


	#Setting null model
	if(which_model=="null"){
		B1="diagonal and unequal"
	}else if(which_model=="unconstrained"){
		B1='unconstrained'
	}else if(which_model=="pencen"){
        B2=matrix(list(0),nrow=length(liste_sp),ncol=length(liste_sp),dimnames=list(liste_sp,liste_sp))
        corres=read.table(paste("corres_hernandez.csv",sep=''),sep=";",na="NA",header=TRUE)
        for (i in 1:length(liste_sp)){
                s=liste_sp[i]
                if(s=="NEI"){
                        for (j in 1:length(liste_sp)){
                                s2=liste_sp[j]
                                B2[i,j]=paste(s2,s,sep="")
                        }
                }else{
                t1=as.character(corres$Type[corres$Code==s])
                for (j in 1:length(liste_sp)){
                        s2=liste_sp[j]
                        if(s2=="NEI"){
                                B2[i,j]=paste(s2,s,sep="")
                        }else{
                        t2=as.character(corres$Type[corres$Code==s2])
                        if(t1==t2){
                                B2[i,j]=paste(s2,s,sep="")
                        }else{
                                if(s==s2){
                                        B2[i,j]=paste(s2,s,sep="")
                                }
                        }
                }
                }}
        }
		B1=B2
	}
	analyse_MARSS(tab_plankton,tab_cov,B1,paste(lieu[l],which_model,which_NEI,which_timestep,which_sp,"essai.RData",sep="_"),boot=FALSE)
	
}
