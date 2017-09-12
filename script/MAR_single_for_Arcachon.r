graphics.off()
rm(list=ls())

library('lubridate')
library('zoo')
library('MARSS')
source("./script/MARSS_clean.r")

which_NEI="null" #NEI can be either a population or a covariate, or not used
which_timestep="regular" #for now, only regular is implemented, but I do think we should do monthly average in order to prepare an inter-site comparison
which_sp="common" #can be single (species depend on each site) or "common" (species are the ones in all sites)
set.seed(42)
model_option=c("unconstrained","null","pencen","diatdin","inter") #for now, can be "null", "unconstrained" and "pencen"
corres=read.table(paste("corres_hernandez.csv",sep=''),sep=";",na="NA",header=TRUE)

sp=c("AST","NIT","PSE","SKE","CHA","GUI","LEP","RHI","GYM","PRP","CRY","EUG")
cov3_tot=c("TEMP","SAL")


#lieu="Teychan"
lieu="B7"
tabbis=read.csv(paste("./data/",lieu,"_base.csv",sep=""),na.strings="NA",header=TRUE,sep=";",dec=".")
dates=as.Date(tabbis$Date)
tab=tabbis[year(dates)>1996,]#Using data from 1997, because cryptophytes were not counted (or badly, for the first year, before
dates=dates[year(dates)>1996]
tab_sp=tab[,sp]

###Treating missing values for both species and covariates###
consecutif=2 #Number of missing values above which we keep the NA
timestep=14 #Regular time lapses between two observations
dates_bis=seq(dates[1],dates[length(dates)],timestep) #Regular time grid
tab_sp=na.approx(tab_sp,maxgap=consecutif,x=dates,xout=dates_bis,na.rm=FALSE) #Interpolation over regular time grid
#Replace missing values
for (s in sp){
        tab_sp[is.na(tab_sp[,s]),s]=runif(sum(is.na(tab_sp[,s])),0,min(tab_sp[,s],na.rm=TRUE))
}
tab_cov_bis=matrix(NA,length(dates_bis),length(cov3_tot))
colnames(tab_cov_bis)=cov3_tot
for (c in cov3_tot){
        tab_cov_bis[,c]=approx(tab[,c],x=dates,xout=dates_bis)$y
}

	#ON range les espèces en les regroupant
	pen=c()
	cen=c()
	dino=c()
	other=c()
	liste_sp=sp
	for (s in 1:length(liste_sp)){
                t1=as.character(corres$Type[corres$Code==liste_sp[s]])
		if(t1=="P"){
			pen=c(pen,liste_sp[s])
		}else if(t1=="C"){
                        cen=c(cen,liste_sp[s])
                }else if(t1=="D"){
                        dino=c(dino,liste_sp[s])
                }else if(t1=="O"){
                        other=c(other,liste_sp[s])
                }else{
			print(s)
			print(t1)
		}
	}
	liste_sp=c(cen,pen,dino,other)
	tab_plankton=tab_sp[,liste_sp]


	#Log transfo for species abundance and scaling for all time series
	tab_plankton=log(tab_plankton)
	tab_plankton=t(scale(tab_plankton[2:(length(dates_bis)-1),]))
	tab_cov=t(scale(tab_cov_bis[2:(length(dates_bis)-1),]))
	rownames(tab_plankton)=liste_sp
	rownames(tab_cov)=cov3_tot


	#Setting interaction matrix
	for(which_model in model_option){
	if(which_model=="null"){
		B1="diagonal and unequal"
	}else if(which_model=="unconstrained"){
		B1='unconstrained'
	}else if(which_model=="pencen"){
        B2=matrix(list(0),nrow=length(liste_sp),ncol=length(liste_sp),dimnames=list(liste_sp,liste_sp))
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
	}else if(which_model=="diatdin"){
        B2=matrix(list(0),nrow=length(liste_sp),ncol=length(liste_sp),dimnames=list(liste_sp,liste_sp))
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
                        if(t1==t2 || (t1 %in% c('P','C')&& t2 %in% c('P','C'))){
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

	}else if(which_model=="inter"){
        B2=matrix(list(0),nrow=length(liste_sp),ncol=length(liste_sp),dimnames=list(liste_sp,liste_sp))
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
                        if((t1 %in% c('P','C')&& t2 %in% c("D","O"))|(t1 %in% c('D','O')&& t2!=t1)){
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
	analyse_MARSS(tab_plankton,tab_cov,B1,paste(lieu,which_model,which_NEI,which_timestep,which_sp,"AR.RData",sep="_"),boot=TRUE)
	}
