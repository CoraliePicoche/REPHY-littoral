##CP 26/09/201 7 : Evaluating synchrony in the community using Vasseur and Gaedke (2007) method

rm(list=ls())
graphics.off()
library(lubridate)
library(lomb)
library(zoo)

print("WARNING: this script is not totally over: I've not completely tested it yet (I think the Hernandez dataset is not enough, we should go back to all genera in Quadrige) + I have not detrended for the main frequency, so I can't swear it works. Moving on to Vasseur & et. 2014")

M=6 #bartlett window
corres=read.table(paste("corres_hernandez.csv",sep=''),sep=";",na="NA",header=TRUE) #We want to know which species is pennate, centric, or dino.
tab_sp=read.table('data/lieu_sp_post_reconstruct_pour_MAR.csv',header=TRUE,na.strings="",sep=";")
lieu=gsub('.',' ',colnames(tab_sp),fixed=TRUE) #useful for Men er Roue
colnames(tab_sp)=lieu

groupe=c("BZ","MO","SU")#,"AR")
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
#for (ll in 1:1){
        if(g=="AR"){
		tab=read.csv(paste("./data/",option_lieu[ll],"_base.csv",sep=""),na.strings="NA",header=TRUE,sep=";",dec=".")
	}else{
        	tab=read.table(paste("data/corres_hernandez_",option_lieu[ll],'.txt',sep=''),sep=";",na="NA",header=TRUE)
	}
        dates=as.Date(tab$Date)
        tab=tab[year(dates)>=1996,]#Using data from 1996
        dates=dates[year(dates)>=1996]

	consecutif=2 #Number of missing values above which we keep the NA
	timestep=14 #Regular time lapses between two observations
	dates_bis=seq(dates[1],dates[length(dates)],timestep) #Regular time grid

	#Species identification
#	sp=colnames(tab)[2:(dim(tab)[2]-1)]
	liste_sp=as.character(tab_sp[!is.na(tab_sp[,option_lieu[ll]]),option_lieu[ll]])
	sp=liste_sp[-which(liste_sp=="NEI")]
	pen=c()
	cen=c()
	din=c()
	oth=c()
	t_all=c()
        for (s in 1:length(sp)){
                t1=as.character(corres$Type[corres$Code==sp[s]])
		t_all=c(t_all,t1)
		if(t1=="P"){
			pen=c(pen,sp[s])
		}else if(t1=="C"){
			cen=c(cen,sp[s])
		}else if(t1=="D"){
			din=c(din,sp[s])
		}else{
			oth=c(oth,sp[s])
		}
        }
	diat=c(pen,cen)


	subset=sp #could be pen, cen, etc.
	#This is dirty
	mat_freq=matrix(NA,nrow=250,ncol=length(sp))	
	mat_spec=matrix(NA,nrow=250,ncol=length(sp))	
	#After a first screening, the minimum set of value we can have is between 0.0009 and 0.010, with 83 points, that is between 1/3 cycle per year and more than 3 cycles per year
	for (s in 1:length(sp)){
	#First, we compute the Lomb-Scargle periodogram
	p=lsp(log(tab[,sp[s]]),as.numeric(dates-dates[1]),main=option_lieu[ll],plot=F)
	#Smoothing to get the spectrum, with Bartlett of breadth M
	aspectrum=rep(0,length(p$scanned))
	for (f in (M+1):(length(p$scanned)-M)){
		for (j in (-M):M){
			aspectrum[f]=aspectrum[f]+p$power[f+j]*(1-abs(j)/M)
		}
		aspectrum[f]=aspectrum[f]/(2*M+1)
		mat_freq[f,s]=p$scanned[f]
	}
	aspectrum[aspectrum==0]=NA
	mat_spec[1:length(p$scanned),s]=aspectrum
	#plot(p$scanned,aspectrum,t="l")
	}
	#I take the shortest of all series of frequency and interpolate each spectrum on this one. I am really decreasing the quality of the time series
	sh=apply(mat_freq,2,function(x) sum(!is.na(x))) 
	find_min=which(sh==min(sh,na.rm=TRUE))
	val_out=mat_freq[,find_min]
	val_out=val_out[!is.na(val_out)]
	
	mat_bis=matrix(NA,nrow=length(val_out),ncol=length(sp))
	mat_bis[,find_min]=mat_spec[!is.na(mat_spec[,find_min]),find_min]
	colnames(mat_bis)=sp
	for (s in 1:length(sp)){
		if(s!=find_min){
			val_in=mat_freq[!is.na(mat_freq[,s]),s]
			y_in=mat_spec[!is.na(mat_spec[,s]),s]
			y_out=approx(y=y_in,x=val_in,xout=val_out)$y
			mat_bis[,s]=y_out
		}
	}
	
	#I compute the community spectrum
	com=apply(tab[,subset],1,sum,na.rm=TRUE)
        p=lsp(log(com+1),as.numeric(dates-dates[1]),main=option_lieu[ll],plot=F) #I check, the addition of 1 only leaves to a power difference of maximum 10^-16
        #Smoothing to get the spectrum, with Bartlett of breadth M
        aspectrum_com=rep(0,length(p$scanned))
        freq_com=rep(0,length(p$scanned))
        for (f in (M+1):(length(p$scanned)-M)){
                for (j in (-M):M){
                        aspectrum_com[f]=aspectrum_com[f]+p$power[f+j]*(1-abs(j)/M)
                }
                aspectrum_com[f]=aspectrum_com[f]/(2*M+1)
                freq_com[f]=p$scanned[f]
        }
	val_in=freq_com[!is.na(freq_com)]
	y_in=aspectrum_com[!is.na(aspectrum_com)]
	y_com=approx(y=y_in,x=val_in,xout=val_out)$y

	#I take the shortest of all series of frequency and interpolate each spectrum on this one. I am really decreasing the quality of the time series
	
	#NOW, we can compute the mean population spectrum for each subcommunity, that is pen, cen, dino
	# I compute the weights
	weight=apply(tab[,subset],2,sum,na.rm=TRUE)/sum(tab[,sp],na.rm=TRUE)
	#I compute the mean population spectrum
	mat_pop=exp(apply(weight*log(mat_bis[,subset]),1,sum))

	#DELTA EV
	dev=y_com-mat_pop
	#EV ratio
	evr=y_com/mat_pop


	#Now, we compute the transition value
	#For that, we need the unevenness
	weight=sort(weight,decreasing=TRUE)
	plou=log(weight[2:length(weight)])-log(weight[1])
	x=2:length(weight)-1
        model_weight=lm(plou~0+x)
	a=-as.numeric(model_weight$coefficients)	
	n=length(weight)

	EVRT=(1-exp(-a))*(1+exp(-a*n))/(1+exp(-a))*(1-exp(-a*n))
	par(mfrow=c(1,2))
	plot(val_out*365,dev,t="l")
	abline(h=0,col="grey")

	plot(val_out*365,evr,t="l")
	abline(h=EVRT,col="grey")

#Here, I tried to check the difference between the Lomb-Scargle periodogram and the one you get from the spectrum function (which is the one you can use directly for synchrony computation). For now, I think I don't exactly get it... so I suggest we work with LS on Vasseur and Gaedke, and keep that good old spectrum function for coherence
#	sp=c('AST','CHA')
#	tab_sp=na.approx(tab[,sp],maxgap=consecutif,x=dates,xout=dates_bis,na.rm=FALSE) #Interpolation over regular time grid
#	for (s in sp){
#        	tab_sp[is.na(tab_sp[,s]),s]=runif(sum(is.na(tab_sp[,s])),0,min(tab_sp[,s],na.rm=TRUE))
#	}
#	pbis=spectrum(log(tab_sp[,'CHA']),method="pgram",main=option_lieu[ll],plot=T,taper=0.0,detrend=F)
#	aspectrumbis=rep(0,length(pbis$freq))
 #       for (f in (M+1):(length(pbis$freq)-M)){
 #               for (j in (-M):M){
 #                       aspectrumbis[f]=aspectrumbis[f]+pbis$spec[f+j]*(1-abs(j)/M)
 #               }
 #               aspectrumbis[f]=aspectrumbis[f]/(2*M+1)
 #       }
 #	lines(pbis$freq,aspectrumbis,t="l")


	#Detrending for the 1 cyle/year. We won't do the 2 cycles/yr, because it is not as obvious as in Vasseur and Gaedke
	print("WARNING! I am not detrending yet, see paragraph on p. 2062 on the subject") #Actually, I stopped there and should get to this to compute community and population-level synchrony


	
}#end option_lieu
}#end groupe
