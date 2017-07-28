#Plot all files grouped in genera or subclasses for the sites we want to study

rm(list=ls())
graphics.off()

sites=c("Men er Roue","Loscolo","Croisic","LEperon","Cornard","Auger","Antoine","Lazaret")
corres=c("corres_proposition","corres_hernandez","corres_Arcachon")
corres=c("corres_hernandez")
#corres=c("subclass")


acex=3
acex=3
alwd=1.5
aft_size=2
amain=2

threshold_abondance=seq(0,0.2,0.005)
threshold_missing=seq(0.375,0.575,0.005)
mat_missing=array(0,dim=c(length(sites),length(threshold_missing),36))
mat_abondance=array(0,dim=c(length(sites),length(threshold_abondance),36))


for (t in 1:length(threshold_missing)){
tab_missing=rep(0,36)
tab_abondance=rep(0,36)
for (s in 1:length(sites)){
	for (c in 1:length(corres)){
#		print(corres[c])
#		file_title=read.table(paste(corres[c],'.csv',sep=''),sep=";",header=TRUE)
		pdf(paste(corres[c],sites[s],"visu.pdf",sep="_"),width=30,height=10)
        	tab=read.table(paste("data/",corres[c],'_',sites[s],'.txt',sep=''),sep=";",na="NA",header=TRUE)
		tab$Date=as.Date(tab$Date)
		name=colnames(tab[-1])
		#par(mfrow=c(4,4),mar=c(2,5,3,2))
		layout(matrix(1:2,nrow=1),widths=c(0.9,0.1))
		y1=log10(min(tab[,2:dim(tab)[2]],na.rm=TRUE))
		y2=log10(max(tab[,2:dim(tab)[2]],na.rm=TRUE))
		acol=rainbow(15)
		acol=rev(acol[1:11])
		for (n in 2:dim(tab)[2]){
		par(mar=c(2,5,3,2))
			prop=tab[,n]/rowSums(tab[2:dim(tab)[2]],na.rm=TRUE)*1.0
			moy_prop=mean(prop,na.rm=TRUE)
#			print(name[n-1])
#			print(mean(prop,na.rm=TRUE))
			#plot(tab$Date,log10(tab[,n]),col="black",xaxt="n",t="o",main=paste(name[n-1],file_title$Genus[n-1]),ylim=c(y1,y2),pch=21,bg=acol[1+round(prop*10)],ylab="Log10(abundance)",xlab="",cex=acex,lwd=alwd,cex.lab=aft_size,cex.axis=aft_size,cex.main=amain)
#			if((corres[c]=="subclass")&&(sum(!is.na(tab[,n]))>0)){
			plot(tab$Date,log10(tab[,n]),col="black",xaxt="n",t="o",main=paste(name[n-1]),ylim=c(y1,y2),pch=21,bg=acol[1+round(prop*10)],ylab="Log10(abundance)",xlab="",cex=acex,lwd=alwd,cex.lab=aft_size,cex.axis=aft_size,cex.main=amain)
			legend("topleft",c(paste(format(round(sum(is.na(tab[,n]))/length(tab[,n]),2),nsmall=2,ndigits=2),"% NA"),paste(format(round(moy_prop,2),nsmall=2,ndigits=2),"% abundance")),bty="n",cex=acex)
			if((sum(is.na(tab[,n]))/length(tab[,n]))<(threshold_missing[t])){tab_missing[n-1]=tab_missing[n-1]+1}
			if((!is.na(moy_prop))&(moy_prop>threshold_abondance[t])){tab_abondance[n-1]=tab_abondance[n-1]+1}
#			axis(2)
#			par(new=T)
#			plot(tab$Date,prop,col="red",t="o",ylim=c(0,1),xlab="",ylab="",axes=F,pch=16,lty=2)
			axis.Date(1,tab$Date,cex.axis=aft_size)
			par(mar=c(2,1,3,0.1))
			xl <- 1.0
yb <- 1
xr <- 1.5
yt <- 2
plot(NA,type="n",ann=FALSE,xlim=c(1,2),ylim=c(1,2),xaxt="n",yaxt="n",bty="n")
rect(
     xl,
     head(seq(yb,yt,(yt-yb)/11),-1),
     xr,
     tail(seq(yb,yt,(yt-yb)/11),-1),
     col=acol
    )

mtext(seq(0,1,0.1),side=2,at=tail(seq(yb,yt,(yt-yb)/11)-0.05,-1),las=2,cex=acex/1.5,line=0)
#			axis(4,at=seq(0.,1.0,0.2),label=seq(0.,1.0,0.2))
		}#}
	dev.off()
	}
}

for (si in 1:length(sites)){
for (n in 1:length(name)){
	if(tab_missing[n]>=si){mat_missing[si,t,n]=mat_missing[si,t,n]+1}
	if(tab_abondance[n]>=si){mat_abondance[si,t,n]=mat_abondance[si,t,n]+1}
}
}
}

acol=rainbow(length(sites)-1)
pdf("nb_groups.pdf",width=15)
par(mfrow=c(1,2),mar=c(5,5,1,1))
 plot(threshold_abondance,rowSums(mat_abondance[8,,]),ylab="Nb groups",xlab="Minimum proportion in the community",t="o",pch=21,cex=acex*0.9,lwd=alwd,cex.axis=aft_size,cex.lab=aft_size,bg="black",ylim=c(0,36))
for (si in 2:7){
	lines(threshold_abondance,rowSums(mat_abondance[si,,]),t="o",pch=21,cex=acex*0.9-si*0.2,lwd=alwd,cex.axis=aft_size,cex.lab=aft_size,bg=acol[si],col=acol[si])

}
legend("topright",as.character(paste(c(seq(2,7),8),"sites")),col=c(acol[2:7],"black"),pch=16,cex=acex*2/3,bty="n")
abline(h=10,col="grey",lty=2,lwd=alwd)
 plot(threshold_missing,rowSums(mat_missing[8,,]),ylab="",xlab="Maximum proportion of missing values",t="o",pch=21,cex=acex*0.9,lwd=alwd,cex.axis=aft_size,cex.lab=aft_size,bg="black",ylim=c(0,21))
for (si in 2:7){
	lines(threshold_missing,rowSums(mat_missing[si,,]),t="o",pch=21,cex=acex*0.9-si*0.2,lwd=alwd,cex.axis=aft_size,cex.lab=aft_size,bg=acol[si],col=acol[si])
}
abline(h=10,col="grey",lty=2,lwd=alwd)
dev.off()

choice_threshold_abondance=0.03
choice_threshold_missing=0.5
nb_sites=6
id1=which(threshold_abondance==choice_threshold_abondance)
id2=which(threshold_missing==choice_threshold_missing)
print("Abondance")
print(colnames(tab)[1+which(mat_abondance[nb_sites,id1,]==1)])
print("Missing")
print(colnames(tab)[1+which(mat_missing[nb_sites,id2,]==1)])
