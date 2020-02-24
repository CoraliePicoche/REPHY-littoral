#CP 02/10/2017 Trying to code the wavelet-version of synchrony according to Vasseur et al. 2014, using Keitt's (2008) mvcwt package

rm(list=ls())
graphics.off()
library("zoo")
library('lubridate')
library('mvcwt')
library('foreach')
#library('parallel')
#library('doParallel')
source("script/wmr.boot_perso.r")

###Loading data
#We have to standardize time and scales at which we sample the wavelet transformation
tab_global=read.table('data/lieu_sp_post_reconstruct_pour_MAR.csv',header=TRUE,na.strings="",sep=";")
lieu=colnames(tab_global)
lieu=gsub('.',' ',lieu,fixed=TRUE) #useful for Men er Roue

groupe=c("BZ","MO","SU","AR")
id_lieu=0
id_g=0
time_between_samples=c()
total_period_sampling=c()
annee_min=c()
annee_max=c()
for (g in groupe){
        if(g=="BZ"){
                option_lieu=c("Men er Roue","Loscolo","Croisic")
        }else if(g=="MO"){
                option_lieu=c("LEperon","Cornard","Auger")
        }else if(g=="SU"){
                option_lieu=c("Antoine","Lazaret")
        }else if(g=="AR"){
                option_lieu=c("Teychan","B7")
        }

	for (ll in 1:length(option_lieu)){
#for (ll in 2:2){
        #Hydro variables : Seasonality
	        if(g=="AR"){
		        tab=read.csv(paste("./data/",option_lieu[ll],"_base.csv",sep=""),na.strings="NA",header=TRUE,sep=";",dec=".")
	        	dates=as.Date(tab$Date)
		}else{
			tab=read.table(paste("data/corres_hernandez_",option_lieu[ll],'.txt',sep=''),sep=";",na="NA",header=TRUE)
			dates=as.Date(tab$Date)
		}
		time_between_samples=c(time_between_samples,3*mean(diff(dates),na.rm=TRUE))
		total_period_sampling=c(total_period_sampling,0.5*(dates[length(dates)]-dates[1]))
		annee_min=c(annee_min,year(min(dates,na.rm=TRUE)))
		annee_max=c(annee_max,year(max(dates,na.rm=TRUE)))
	}
}
time_between_samples_global=max(time_between_samples,na.rm=TRUE)
total_period_sampling_global=min(total_period_sampling,na.rm=TRUE)
annee_min_global=max(annee_min)
annee_min_global=1997
annee_max_global=min(annee_max)
annee_max_global=2015
dt=15
seq_time=seq(as.Date(paste(annee_min_global,"-01-01",sep="")),as.Date(paste(annee_max_global,"-12-31",sep="")),by="15 days")


ds=0.1
j_min_scale=ceiling(1/ds*log2(time_between_samples_global/365.25))
j_max_scale=floor(1/ds*log2(total_period_sampling_global/365.25))
seq_scale=365.25*2^(ds*seq(j_min_scale,j_max_scale))

#Creating regular table
tab_seq_ns=matrix(1,nrow=length(seq_time),length(seq_scale))

for (g in groupe){
        id_g=id_g+1
        if(g=="BZ"){
                option_lieu=c("Men er Roue","Loscolo","Croisic")
        }else if(g=="MO"){
                option_lieu=c("LEperon","Cornard","Auger")
        }else if(g=="SU"){
                option_lieu=c("Antoine","Lazaret")
        }else if(g=="AR"){
                option_lieu=c("Teychan")#,"B7")
        }
        
        for (ll in 1:length(option_lieu)){
		print(option_lieu[ll])
#                id_lieu=id_lieu+1
#for (ll in 1:1){
        #Hydro variables : Seasonality
                if(g=="AR"){
                        tab=read.csv(paste("./data/",option_lieu[ll],"_base.csv",sep=""),na.strings="NA",header=TRUE,sep=";",dec=".")
                        sp=c("AST","NIT","PSE","SKE","CHA","GUI","LEP","RHI","GYM","PRP","CRY","EUG")
                        dates=as.Date(tab$Date)
                }else{  
                        tab=read.table(paste("data/corres_hernandez_",option_lieu[ll],'.txt',sep=''),sep=";",na="NA",header=TRUE)
                        a=table(as.matrix(tab_global[, lieu %in% option_lieu]))
                        sp=dimnames(a[a==length(option_lieu)])[[1]]
                        dates=as.Date(tab$Date)
                }
	#	dmin=min(dates)
	#	dmax=max(dates)
	#	tab_tmp=tab_seq_ns
	#for (s in 1:length(seq_scale)){
	#	coi=c(dmin+sqrt(2)*seq_scale[s],dmax-sqrt(2)*seq_scale[s])
	#	[seq_time<coi[1]|seq_time>coi[2],s]=NA
	#}


	tab_sp=na.approx(tab[,sp],maxgap=2,x=dates,xout=dates,na.rm=FALSE) #We only interpolate when there are missing values and we think there should not be. We don't regularize the tieme grid as this will be done by the wavelet function itself

#Don't tolerate the Nan values
	for (s in sp){
#	tab_sp[is.na(tab_sp[,s]),s]=runif(sum(is.na(tab_sp[,s])),0,min(tab_sp[,s],na.rm=TRUE))
		tab_sp[is.na(tab_sp[,s]),s]=rep(0,sum(is.na(tab_sp[,s])))
	}
	tab_sp=log(tab_sp+1)

#	cl=makeCluster(2)
#	registerDoParallel(cl,cores=2)

	mm=mvcwt(dates-dates[1],tab_sp,scales=seq_scale,loc=as.numeric(seq_time-dates[1]))
	
	ww=wmr.boot_perso(mm,reps=1000,aquantile=c(0.025,0.5,0.975))
	filename=paste(option_lieu[ll],"_synch_vass_keitt.RData",sep="")
	save(mm,ww,seq_time,seq_scale,dates,file=filename)
#	stopCluster(cl)
#	readLines()
}
}




