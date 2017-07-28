#Shows all variable monitoring for each site (which ones for which site ? For how long? Are there long gaps?)

graphics.off()
rm(list=ls())

ref=read.table("data/raw/Q2_REPHY_SRN_points.csv",sep=";",header=TRUE)
corres=read.table("data/lieu_corres.csv",sep=";",header=TRUE)

pdf("graphe/Gantt_REPHY.pdf",width=20)
#filename_to_search=c("Q2_170418_other_UTF8_only_good.csv","Q2_170418_site_Tania_UTF8_only_good.csv","Q2_170511_SRN_Chnord_UTF8_only_good.csv")
filename_tot=list.files(path="data/",pattern="^[A-Z].*flortot_only.csv")

for (f in 1:length(filename_tot)){
tab=read.table(paste("data/",filename_tot[f],sep=""),sep=";",header=TRUE)
date=as.Date(tab$Date)
name=strsplit(strsplit(filename_tot[f],split="_")[[1]][1],split="-")[[1]][1]
id=corres$Lieu_id[grep(strsplit(strsplit(filename_tot[f],split="_")[[1]][1],split="-")[[1]][1],corres$Lieu_mnemonique,fixed=TRUE)]
if(length(strsplit(strsplit(filename_tot[f],split="_")[[1]][1],split="-")[[1]])>1){
id=c(id,corres$Lieu_id[grep(strsplit(strsplit(filename_tot[f],split="_")[[1]][1],split="-")[[1]][2],corres$Lieu_mnemonique,fixed=TRUE)])
}
d_min=as.Date("1987-01-01")
d_max=as.Date("2017-01-01")

	if (f%%8==1){
		plot(d_min,0,t='n',xlab='',ylab='',xlim=c(d_min,d_max),ylim=c(1,8),xaxt="n",yaxt="n")
		axis.Date(side=1,unique(date))
	}
	text(y=f%%8+1,par("usr")[1],labels=paste(name,"\n (",round(ref$Latitude[ref$Lieu_id==id[1]]),",",round(ref$Longitude[ref$Lieu_id==id[1]]),")",sep=""),srt=45,xpd=TRUE)
	dd=unique(date)
	points(dd,rep(f%%8+1,length(dd)),pch=16,col="black")

#Now we check for abiotic variables
	file=paste("graphe/",name,"_var_time.pdf",sep="")
	pdf(file,width=20)
	tab=read.table("data/raw/Q2_170418_site_Tania_UTF8_only_good.csv",sep=";",header=TRUE)
	tab_lieu=subset(tab,is.element(Lieu_id,id))
	date_lieu=as.Date(tab_lieu$Date,"%d/%m/%Y")
	boo=FALSE
	v=0
	var_u=unique(tab_lieu$Res_code)
	var_u=var_u[!grepl('SOU',var_u)] #Removing all tests on mice, we can't use them
		
	tab2=read.table("data/raw/Q2_170511_SRN_Chnord_UTF8_only_good.csv",sep=";",header=TRUE)
	tab_lieu2=subset(tab2,is.element(Lieu_id,id))
	if(dim(tab_lieu2)[1]>0){
		boo=TRUE
		date_lieu2=as.Date(tab_lieu2$Date,"%d/%m/%Y")
		var_u2=unique(tab_lieu2$Res_code)
		var_u2=var_u2[!grepl('SOU',var_u2)] #Removing all tests on mice, we can't use them
	}
	
	if(length(var_u)>0){
	for (v in 1:length(var_u)){
		if(v%%12==1){
			plot(d_min,0,t='n',xlab='',ylab='',xlim=c(d_min,d_max),ylim=c(1,12),xaxt="n",yaxt="n")
			axis.Date(side=1,unique(date))
			legend("topleft",c("REPHY","SRN"),pch=16,col=c("black","grey"),cex=2,bty="n")
		}
			text(y=v%%12+1,par("usr")[1],labels=paste(tab_lieu$Res_code[tab_lieu$Res_code==var_u[v]][1],sep=""),srt=45,xpd=TRUE)
			dd=unique(date_lieu[grep(var_u[v],tab_lieu$Res_code,fixed=TRUE)])
			points(dd,rep(v%%12+1,length(dd)),pch=16,col="black",cex=1.3)
			if(boo){	
			dd2=unique(date_lieu2[grep(var_u[v],tab_lieu2$Res_code,fixed=TRUE)])
			points(dd2,rep(v%%12+1,length(dd2)),pch=16,col="grey",cex=0.75)
			idv=grep(var_u[v],var_u2,fixed=TRUE)
			var_u2=var_u2[-idv]}
	}
	}
	if(boo){
	if(length(var_u2)>0){
        for (vprime in 1:length(var_u2)){
                if((v+vprime)%%12==1){
                        plot(d_min,0,t='n',xlab='',ylab='',xlim=c(d_min,d_max),ylim=c(1,12),xaxt="n",yaxt="n")
                        axis.Date(side=1,unique(date))
                }       
                        text(y=(v+vprime)%%12+1,par("usr")[1],labels=paste(tab_lieu2$Res_code[tab_lieu2$Res_code==var_u2[vprime]][1],sep=""),srt=45,xpd=TRUE)
                        dd2=unique(date_lieu2[grep(var_u2[vprime],tab_lieu2$Res_code,fixed=TRUE)])
                        points(dd2,rep((v+vprime)%%12+1,length(dd2)),pch=16,col="grey")
        }
	}}
		dev.off()
	}
dev.off()
#}

