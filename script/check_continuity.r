#This file to see if gaps are homogeneous in the times series, and to see their average duration per site

graphics.off()
rm(list=ls())
library('lubridate')

fil=list.files(path="data/",pattern="^[A-Z].*flortot_only.csv")

ccolfunc <- colorRampPalette(c("lightblue1", "darkblue"))
colo=ccolfunc(736-210)
pdf("graphe/boxplot_gap_all.pdf",width=20)
par(oma=c(5,2,1,1))
pdf("graphe/gap_per_season.pdf")
for (f in 1:length(fil)){
	missing_month=0
	no_missing_month=0
	missing_year=0
	missing_hiver=c()
	missing_ete=c()
	missing_printemps=c()
	missing_automne=c()
	duration_gap=c()
	duration_gap_under=c()
	tab=read.table(paste("data/",fil[f],sep=""),header=TRUE,sep=";")
        dates=as.Date(tab$Date)
	d1=dates[1]
	for (d in 2:length(dates)){
		d2=dates[d]
		nb_days=d2-d1
		if (nb_days<31){
			no_missing_month=no_missing_month+1
			duration_gap_under=c(duration_gap_under,nb_days)
		}
		if (nb_days>31&nb_days<365){
			if(month(d1)%in%c(1,2,3)){
				missing_hiver=c(missing_hiver,nb_days)
			}
			if(month(d1)%in%c(4,5,6)){
				missing_printemps=c(missing_printemps,nb_days)
			}
			if(month(d1)%in%c(10,11,12)){
				missing_automne=c(missing_automne,nb_days)
			}
			if(month(d1)%in%c(7,8,9)){
				missing_ete=c(missing_ete,nb_days)
			}
			missing_month=missing_month+as.period(interval(d1,d2))@month
			duration_gap=c(duration_gap,nb_days)
		}
		if(nb_days>=365){
			missing_year=missing_year+as.period(interval(d1,d2))@year
		}
		d1=d2
	}
	#The seasonal plot is active
	mini=min(c(missing_hiver,missing_ete,missing_printemps,missing_automne))
	maxi=max(c(missing_hiver,missing_ete,missing_printemps,missing_automne))
	plot(rep(1,length(missing_hiver)),missing_hiver,xlim=c(1,4),ylim=c(mini,maxi),main=strsplit(fil[f],split="_")[[1]][1],xlab="",ylab="Gap duration>1 month (days)",xaxt='n',pch=21,cex=2,col="black",bg=adjustcolor("blue",alpha=.6))
	axis(1,at=c(1,2,3,4),labels=c(paste("Winter",length(missing_hiver),sep="\n"),paste('Spring',length(missing_printemps),sep="\n"),paste('Summer',length(missing_ete),sep="\n"),paste('Autumn',length(missing_automne),sep="\n")))
	points(rep(2,length(missing_printemps)),missing_printemps,pch=21,cex=2,bg=adjustcolor("green",alpha=.6),col="black")
	points(rep(3,length(missing_ete)),missing_ete,pch=21,cex=2,bg=adjustcolor("orange",alpha=.6),col="black")
	points(rep(4,length(missing_automne)),missing_automne,pch=21,bg=adjustcolor("red",alpha=.6),cex=2,col="black")
	dev.set(dev.prev(which=dev.cur()))
	if (f==1||f==11||f==21){
		b=boxplot(at=f%%10,x=duration_gap,add=FALSE,xlab="",ylab="Gap duration>1 month (days)",xaxt="n",xlim=c(-0.5,10),ylim=c(0,180),col=colo[length(dates)-210+1],pch=16,cex=2,outcol="slategray3",medcol="cyan")
	}else{
		b=boxplot(at=f%%10,x=duration_gap,add=TRUE,yaxt="n",xaxt="n",col=colo[length(dates)-210+1],pch=16,cex=2,outcol="slategray3",medcol="cyan")
	}
	text(f%%10,y=b$stats[5,1]+25,labels=paste(year(min(dates)),year(max(dates)),sep="-"))
	text(f%%10,y=b$stats[5,1]+10,labels=paste(missing_year,"MY",missing_month,"MM",format(round(mean(duration_gap_under,na.rm=TRUE),0),nsmall=0),"/",format(round(mean(c(duration_gap_under,duration_gap),na.rm=TRUE),0),nsmall=0),sep=" "))
	text(f%%10,y=30,labels=paste(length(dates),"lines",sep=" "))
	text(f%%10,y=0,labels=paste(sum(b$out>178),"gap>1/2y",sep=" "))
	text(f%%10,y=par("usr")[3]-20,labels=strsplit(strsplit(fil[f],split="_")[[1]][1],split="-")[[1]][1],srt=45,xpd=TRUE)
	dev.set(dev.next(which=dev.cur()))
}
graphics.off()
