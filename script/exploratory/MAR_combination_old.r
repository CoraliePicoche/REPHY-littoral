###This is clearly incomplete, because I'm not sure which models we should try (for now, they're on paper)

graphics.off()
rm(list=ls())

library('lubridate')
library('zoo')
library('MARSS')
source("/home/cpicoche/Documents/Plankton/script/MARSS_clean.r")

set.seed(42)
timestep=14
consecutif=2

tab_sp=read.table('data/lieu_sp_post_reconstruct_pour_MAR.csv',header=TRUE,na.strings="",sep=";")
lieu=colnames(tab_sp)
lieu=gsub('.',' ',lieu,fixed=TRUE) #useful for Men er Roue

a=table(as.matrix(tab_sp))
liste_sp=dimnames(a[a==dim(tab_sp)[2]])[[1]]

cov3_tot=c("TEMP","SALI")

tab_global
for (l in 1:length(lieu)){
        #Biotic variables
        tab=read.table(paste("data/corres_hernandez_",lieu[l],'.txt',sep=''),sep=";",na="NA",header=TRUE)
        dates=as.Date(tab$Date)
        tab=tab[year(dates)>=1996,]#Using data from 1996
        dates=dates[year(dates)>=1996]

        dates_bis=seq(dates[1],dates[length(dates)],timestep) #Regular time grid

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

        #Log transfo for species abundance and scaling for all time series
        tab_plankton=log(tab_plankton)
        tab_plankton=t(scale(tab_plankton[2:(length(dates_bis)-1),]))
        tab_cov=t(scale(tab_cov_bis[2:(length(dates_bis)-1),]))
        rownames(tab_plankton)=liste_sp
        rownames(tab_cov)=cov3_tot
}
