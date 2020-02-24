##CP 26/09/2017. Wanted to test the correlation between seasonality and pairwise coefficients, as in Usinowicz 2017 in Nature

rm(list=ls())
graphics.off()
library('lubridate')

season=c()
groupe=c("BZ","MO","SU","AR")
pdf("Rapport/graphe/competition_vs_seasonality.pdf")
plot(0,0,t="n",xlim=c(-.8,-.3),ylim=c(-0.1,0.17),xlab="Seasonality",ylab="Pairwise competition") #Taking into account both inter and intragroup competition
#plot(0,0,t="n",xlim=c(-.8,-.3),ylim=c(-0.005,0.005),xlab="Seasonality",ylab="Pairwise competition") #Taking into account only intergroup competition
id_lieu=0
colog=c("blue","brown","red","green")
id_g=0
for (g in groupe){
	id_g=id_g+1
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
	id_lieu=id_lieu+1
#for (ll in 2:2){
        #Hydro variables : Seasonality
	if(g=="AR"){
	tab_cov=read.csv(paste("./data/",option_lieu[ll],"_base.csv",sep=""),na.strings="NA",header=TRUE,sep=";",dec=".")
	}else{
	tab_cov=read.table(paste("data/",option_lieu[ll],'hydro.txt',sep=''),sep=";",na="NA",header=TRUE)
	}
        dates_cov=as.Date(tab_cov$Date)
	length_d=length(dates_cov)
	date=seq(as.Date(paste(year(dates_cov[1]),month(dates_cov[1]),'01',sep='-')),as.Date(paste(year(dates_cov[length_d]),month(dates_cov[length_d]),'01',sep='-')),by="month")
	min_t=rep(NA,length(date))
	max_t=rep(NA,length(date))
	for (d in 2:length(date)){
		id=dates_cov>=date[d-1]&dates_cov<date[d]
		if(sum(!is.na(tab_cov$TEMP[id]))>0){
		min_t[d-1]=min(tab_cov$TEMP[id],na.rm=TRUE)
		max_t[d-1]=max(tab_cov$TEMP[id],na.rm=TRUE)
		}
	}
	season=c(season,log(sd(min_t,na.rm=TRUE)/mean(min_t,na.rm=TRUE)+sd(max_t,na.rm=TRUE)/mean(max_t,na.rm=TRUE)))

        f1=paste("data/analyse_MAR/",g,"/site_specific/",option_lieu[ll],"_unconstrained_null_regular_common_",g,".RData",sep="")
        load(f1)
        nom=dimnames(fit_log$par$B)[[1]] #Names of the interactions
	B=diag(-1,sqrt(length(nom)),sqrt(length(nom))) #I am heavily using the fact that we are using an unconstrained matrix
	for (n in 1:length(nom)){
        	if(grepl("^\\(",nom[n])){ #default when using unconstrained or diagonal, names are of the form (1,2) for (i,j)
                	i=as.numeric(strsplit(nom[n],split="[\\(,\\)]")[[1]][2])
                        j=as.numeric(strsplit(nom[n],split="[\\(,\\)]")[[1]][3])
			B[i,j]=B[i,j]+cis$par$B[n]
		}else{
			print(nom[n])
		}
	}
	coeff=c()
	for(i in 1:(dim(B)[1]-1)){
#		for(j in (i+1):dim(B)[2]){
		for(j in i:dim(B)[2]){
			coeff=c(coeff,B[i,j]*B[j,i])
		}
	}
	points(season[id_lieu],median(coeff),t="p",pch=16,col=colog[id_g])
	segments(season[id_lieu],median(coeff)-sd(coeff) ,season[id_lieu],median(coeff)+sd(coeff),col=colog[id_g])
	print(mean(coeff))
	abline(h=0,lty=2,col="grey")
}

}
legend("bottomleft",groupe,col=colog,lty=1)
dev.off()	
