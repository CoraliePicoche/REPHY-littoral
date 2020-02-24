#CP 29/09/2017 Trying to code the wavelet-version of synchrony according to Vasseur et al. 2014

#rm(list=ls())
graphics.off()
library("zoo")
library("WaveletComp") #Morlet wavenumber is 6 by default

#Parameter

#Here, we compute the wavelet transform for all
###Loading data
sp=c("AST","NIT","PSE","SKE","CHA","GUI","LEP","RHI","GYM","PRP","CRY","EUG")
filename=paste("/home/cpicoche/Documents/Plankton/data/treated/Teychan_base.csv",sep="")
tabbis=read.csv(filename,na.strings="NA",header=TRUE,sep=";",dec=".")
dates=as.Date(tabbis$Date)
tab_sp=tabbis[,sp]

consecutif=2 #Number of missing values above which we keep the NA
timestep=14 #Regular time lapses between two observations
dates_bis=seq(dates[1],dates[length(dates)],timestep) #Regular time grid
tab_sp=na.approx(tab_sp,maxgap=consecutif,x=dates,xout=dates,na.rm=FALSE) #Interpolation on a regular grid. 
#We COULD take into account the fact the the sampling is somewhat irregular. Results can be slightly different but Vasseur et al. indicate that this can bias the wavelet power
tab_sp=data.frame("date"=dates,tab_sp)

adt=14/(365.25)
tmp=tab_sp
        for (s in c('CHA','AST')){
#                if(s=="CRY"|s=="EUG"|s=="GYM"){
#                        tmp=tab_sp[year(dates_bis)>1996,c("date",s)]
#                }else{
#                        tmp=tab_sp[,c("date",s)]
#                }
                tmp[is.na(tmp[,s]),s]=runif(sum(is.na(tmp[,s])),0,min(tmp[,s],na.rm=TRUE))
}
#                tt=analyze.wavelet(my.data=data.frame(tmp),my.series=s,dt=adt,lowerPeriod=0.25,loess.span=0,method="white.noise",n.sim=1000,dj=1/20,make.pval=F) #loess.span does not help: no need to detrend (tried loess.span=3),default value of dj is 1/20, but 1/50 is smoother ad nicer and results do not change


#We have to standardize time and scales at which we sample the wavelet transformation
#tt$Wave[scale,t]
#Scale = tt$Scale
#Time=? I still need to find the time localization, which is not the sampling date. Vasseur et al. choose every 10 day. Do they interpolate ? I have to compare what happens when we interpolate the time series THEN use the wavelet ; or when we interpolate the wavelet


#Then we compute for each time and scale chosen the module of the sum and the sum of the module at each place #Mod(Wave(scale,time))
#The ratio is a first measure of synchrony: 0 is compensatory, 1 is synchrony


#Then, we test this ratio
#For each site, each scale and each time, we generate 1000 null model and compute the modulus ratio

#We compare what we have at one place to the median of the 1000 relation. We have 1 if its over, 0 if not. We can then know the fraction of lakes for which we are syncrhonysed
