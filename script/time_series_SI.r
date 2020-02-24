#04/01/2019 CP: plotting time series for SI

rm(list=ls())
graphics.off()
library('lubridate')
library('zoo')

set.seed(42)
timestep=14
consecutif=2

#liste_sp=c("CHA","PSE","PRP")#,"SKE") #Wondering if we keep SKE; the first three are from different groups, the last one is another centric diatom
liste_sp=c("AST","CHA","DIT","GUI","LEP","NIT","PLE","PSE","RHI","SKE","THP","THL","GYM","PRO","PRP","SCR","CRY","EUG")

tab_sp=read.table('data/lieu_sp_post_reconstruct_pour_MAR.csv',header=TRUE,na.strings="",sep=";")
lieu=colnames(tab_sp)
lieu=gsub('.',' ',lieu,fixed=TRUE) #useful for Men er Roue

groupe=c("BZ","MO","AR","SU")

tab_tot=array(NA,dim=c(10,550,1+length(liste_sp)),dimnames=list(1:10,1:550,c("Date",liste_sp))) #10 sites, (max) sample points, up to 4 species and a date column

xlimi=c(as.Date("1996-01-01"),as.Date("2016-01-01"))
colo=c()
id_lieu=0
for (g in groupe){
        if(g=="BZ"){
                option_lieu=c("Men er Roue","Loscolo","Croisic")
                colo=c(colo,rep("green",length(option_lieu)))
        }else if(g=="MO"){
                option_lieu=c("LEperon","Cornard","Auger")
                colo=c(colo,rep("darkblue",length(option_lieu)))
        }else if(g=="SU"){
                option_lieu=c("Antoine","Lazaret")
                colo=c(colo,rep("darkred",length(option_lieu)))
        }else if(g=="AR"){
                option_lieu=c("Teychan","B7")
                colo=c(colo,rep("cyan",length(option_lieu)))
        }
        for (ll in 1:length(option_lieu)){
		id_lieu=id_lieu+1
		dimnames(tab_tot)[[1]][id_lieu]=option_lieu[ll]
		
                if(g=="AR"){
                        tab=read.csv(paste("./data/",option_lieu[ll],"_base.csv",sep=""),na.strings="NA",header=TRUE,sep=";",dec=".")
                }else{
                        tab=read.table(paste("data/corres_hernandez_",option_lieu[ll],'.txt',sep=''),sep=";",na="NA",header=TRUE)
                }

	        dates=as.Date(tab$Date)
        	tab=tab[year(dates)>=1996&year(dates)<2016,]#Using data from 1996
        	dates=dates[year(dates)>=1996&year(dates)<2016]

	        dates_bis=seq(dates[1],dates[length(dates)],timestep) #Regular time grid
		
		tab_tot[id_lieu,1:length(dates_bis),"Date"]=dates_bis
		
		if(g!="AR"){
		a=table(as.matrix(tab_sp[,lieu %in% option_lieu]))
                liste_sp_1=dimnames(a[a==length(option_lieu)])[[1]]
		liste_sp_1=liste_sp_1[liste_sp_1!="NEI"]
		}else{
			liste_sp_1=c("AST","NIT","PSE","SKE","CHA","GUI","LEP","RHI","GYM","PRP","CRY","EUG")
		}
		tab_plankton=na.approx(tab[,liste_sp_1],maxgap=consecutif,x=dates,xout=dates_bis,na.rm=FALSE) #Interpolation over regular time grid
#Replace missing values
                for (s in liste_sp_1){
                        tab_plankton[is.na(tab_plankton[,s]),s]=runif(sum(is.na(tab_plankton[,s])),0,min(tab_plankton[,s],na.rm=TRUE))
			tab_tot[id_lieu,1:length(dates_bis),s]=log10(tab_plankton[,s])
                }
	}
}


#All sites, plankton species separately
pdf("article/graphe/time_series_plankton_allsite.pdf",width=10)
liste_sp=c("CHA","PSE","PRP")#,"SKE") #Wondering if we keep SKE; the first three are from different groups, the last one is another centric diatom
maxi=max(c(tab_tot[,,liste_sp]),na.rm=T)
par(mfrow=c(length(liste_sp),1),mar=c(4,5,1,1))
for(s in liste_sp){
	plot(0,0,t="n",xlim=xlimi,ylim=c(0,maxi),ylab=paste(s),xaxt="n",xlab="",cex=2,cex.axis=2,cex.lab=2)
	for(l in 1:10){
		points(tab_tot[l,,"Date"],tab_tot[l,,s],col=colo[l],pch=16)
		lines(tab_tot[l,,"Date"],tab_tot[l,,s],col=colo[l],lty=1)
	}
	axis.Date(1,at=seq(xlimi[1], xlimi[2], by = "year"),format="%Y",cex.axis=2,cex.lab=2)
}
dev.off()

pdf("article/graphe/time_series_plankton_site_by_site_focus_on_CHA.pdf",width=20)
liste_sp=c("AST","DIT","GUI","LEP","NIT","PLE","PSE","RHI","SKE","THP","THL","GYM","PRO","PRP","SCR","CRY","EUG","CHA")
par(mfrow=c(3,1),mar=c(4,5,1,1))
for(l in 1:10){
	maxi=max(c(tab_tot[l,,liste_sp]),na.rm=T)
        plot(0,0,t="n",xlim=xlimi,ylim=c(0,maxi),ylab=dimnames(tab_tot)[[1]][l],xaxt="n",xlab="",cex=2,cex.axis=2,cex.lab=2)
	for(s in 1:length(liste_sp)){
		if(liste_sp[s]=="CHA"){
			colo1="red"
		}else{
			colo1="black"
		}
                points(tab_tot[l,,"Date"],tab_tot[l,,liste_sp[s]],col=colo1,pch=16)
                lines(tab_tot[l,,"Date"],tab_tot[l,,liste_sp[s]],col=colo1,lty=1)
        }
        axis.Date(1,at=seq(xlimi[1], xlimi[2], by = "year"),format="%Y",cex.axis=2,cex.lab=2)
}
dev.off()

pdf("article/graphe/time_series_plankton_site_by_site_rainbow_plankton.pdf",width=20)
liste_sp=c("AST","DIT","GUI","LEP","NIT","PLE","PSE","RHI","SKE","THP","THL","GYM","PRO","PRP","SCR","CRY","EUG","CHA")
par(mfrow=c(3,1),mar=c(4,5,1,1))
colo=rainbow(length(liste_sp))
for(l in 1:10){
	maxi=max(c(tab_tot[l,,liste_sp]),na.rm=T)
        plot(0,0,t="n",xlim=xlimi,ylim=c(0,maxi),ylab=dimnames(tab_tot)[[1]][l],xaxt="n",xlab="",cex=2,cex.axis=2,cex.lab=2)
	for(s in 1:length(liste_sp)){
                points(tab_tot[l,,"Date"],tab_tot[l,,liste_sp[s]],col=colo[s],pch=16)
                lines(tab_tot[l,,"Date"],tab_tot[l,,liste_sp[s]],col=colo[s],lty=1)
        }
        axis.Date(1,at=seq(xlimi[1], xlimi[2], by = "year"),format="%Y",cex.axis=2,cex.lab=2)
}
dev.off()

#Test on centric diatoms
pdf("article/graphe/time_series_plankton_site_by_site_centric.pdf",width=20)
liste_sp=c("CHA","DIT","GUI","LEP","RHI","SKE","THP")
par(mfrow=c(3,1),mar=c(4,5,1,1))
colo=rainbow(length(liste_sp))
for(l in 1:10){
        maxi=max(c(tab_tot[l,,liste_sp]),na.rm=T)
        plot(0,0,t="n",xlim=xlimi,ylim=c(0,maxi),ylab=dimnames(tab_tot)[[1]][l],xaxt="n",xlab="",cex=2,cex.axis=2,cex.lab=2)
        for(s in 1:length(liste_sp)){
                points(tab_tot[l,,"Date"],tab_tot[l,,liste_sp[s]],col=colo[s],pch=16)
                lines(tab_tot[l,,"Date"],tab_tot[l,,liste_sp[s]],col=colo[s],lty=1)
        }
        axis.Date(1,at=seq(xlimi[1], xlimi[2], by = "year"),format="%Y",cex.axis=2,cex.lab=2)
}
dev.off()

#5 most abundant
xlimi=c(as.Date("1996-01-01"),as.Date("2017-06-01"))
liste_sp=c("AST","CHA","DIT","GUI","LEP","NIT","PLE","PSE","RHI","SKE","THP","THL","GYM","PRO","PRP","SCR","CRY","EUG")
sizes=c(3,3,2,2)
colo=rainbow(5)
lili=list(1:3,4:6,7:8,9:10)
for(g in 1:4){
pdf(paste("article/graphe/time_series_plankton_site_by_site_most_abundant",groupe[g],".pdf",sep=""),width=20)
par(mfrow=c(sizes[g],1),mar=c(4,5,1,1))
for(l in lili[[g]]){
        maxi=max(c(tab_tot[l,,liste_sp]),na.rm=T)
	plou=apply(tab_tot[l,,2:dim(tab_tot)[[3]]],2,mean,na.rm=T)
	plou=plou[!is.na(plou)]
	or=order(plou,decreasing=T)
        liste_sp_bis=names(plou)[or][1:5]
	plot(0,0,t="n",xlim=xlimi,ylim=c(0,maxi),ylab=dimnames(tab_tot)[[1]][l],xaxt="n",xlab="",cex=2,cex.axis=2,cex.lab=2)
        for(s in 1:length(liste_sp_bis)){
                points(tab_tot[l,,"Date"],tab_tot[l,,liste_sp_bis[s]],col=colo[s],pch=16,cex=1.25)
                lines(tab_tot[l,,"Date"],tab_tot[l,,liste_sp_bis[s]],col=colo[s],lty=1,lwd=1.5)
        }
        axis.Date(1,at=seq(xlimi[1], xlimi[2], by = "year"),format="%Y",cex.axis=2,cex.lab=2)
	legend('topright',liste_sp_bis,col=colo,pch=16,cex=2,bty="n")
}
dev.off()
}

#5 most abundant
xlimi=c(as.Date("1996-01-01"),as.Date("2017-06-01"))
liste_sp=c("AST","CHA","DIT","GUI","LEP","NIT","PLE","PSE","RHI","SKE","THP","THL","GYM","PRO","PRP","SCR","CRY","EUG")
pdf(paste("article/graphe/time_series_plankton_site_by_site_most_abundant.pdf",sep=""),width=20,height=25)
colo=rainbow(5)
par(mfrow=c(10,1),mar=c(4,5,1,1))
for(l in 1:10){
        maxi=max(c(tab_tot[l,,liste_sp]),na.rm=T)
        plou=apply(tab_tot[l,,2:dim(tab_tot)[[3]]],2,mean,na.rm=T)
        plou=plou[!is.na(plou)]
        or=order(plou,decreasing=T)
        liste_sp_bis=names(plou)[or][1:5]
        plot(0,0,t="n",xlim=xlimi,ylim=c(0,maxi),ylab=dimnames(tab_tot)[[1]][l],xaxt="n",xlab="",cex=2,cex.axis=2,cex.lab=2)
        for(s in 1:length(liste_sp_bis)){
                points(tab_tot[l,,"Date"],tab_tot[l,,liste_sp_bis[s]],col=colo[s],pch=16,cex=1.25)
                lines(tab_tot[l,,"Date"],tab_tot[l,,liste_sp_bis[s]],col=colo[s],lty=1,lwd=1.5)
        }
        axis.Date(1,at=seq(xlimi[1], xlimi[2], by = "year"),format="%Y",cex.axis=2,cex.lab=2)
        legend('topright',liste_sp_bis,col=colo,pch=16,cex=2,bty="n")
}
dev.off()
