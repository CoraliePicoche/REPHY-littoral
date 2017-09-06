rm(list=ls())
graphics.off()

corres=read.table("data/lieu_corres.csv",sep=";",header=TRUE)

filename_tot=c("data/raw/Q2_170418_site_Tania_UTF8_only_good.csv")#,"data/raw/Q2_170511_SRN_Chnord_UTF8_only_good.csv")
#filename_tot=c("Q2_170511_SRN_Chnord_UTF8_only_good.csv")
#par(mfrow=c(1,2))
pdf("points_FLORPARIND.pdf",height=12)
par(mar=c(3,10,1,1))
#file.remove(paste("data/",list.files(pattern='*.flortot_only.csv'),sep=""))
for (filename in filename_tot){
tab=read.table(filename,sep=";",header=TRUE)
date=as.Date(tab$Date,"%d/%m/%Y")

l_u=unique(tab$Lieu_id)
d_min=min(date)
d_max=max(date)
libel=c()
for (l in 1:length(l_u)){
		print(l)
                tab_flore=subset(tab,((Lieu_id==l_u[l])&(Res_code %in% c("FLORPAR","FLORIND"))),select=c("Lieu_id","Lieu_libel","Date","Imm","Taxon","Val"))
     	        libel=c(libel,as.character(tab_flore$Lieu_libel[1]))
		if(l==1){
		plot(0,0,t="n",xlim=c(as.Date("1987-01-01"),as.Date("2016-01-01")),ylim=c(0,length(l_u)+1),ylab="",yaxt="n",xaxt="n")
		}
                if (dim(tab_flore)[1]!=0){
	                tab_flore$Date=as.Date(tab_flore$Date,"%d/%m/%Y")
                	date_flore=sort(unique(as.Date(tab_flore$Date,"%d/%m/%Y")))
			print(libel[l])
			print(c(min(date_flore, na.rm=TRUE),max(date_flore,na.rm=TRUE)))
			points(date_flore,rep(l,length(date_flore)),col="black",pch=16)
		}
	}	
axis(2,at=1:length(l_u),lab=libel,cex=.75,las=2)
axis.Date(1, at=seq(from=as.Date("1987-01-01"),to=as.Date("2016-01-01"),by="year"), format="%Y")
#axis(1,at=seq(from=as.Date("1987-01-01"),to=as.Date("2016-01-01"),by="year"))
}	
dev.off()

pdf("time_series_florparind.pdf")
for (filename in filename_tot){
tab=read.table(filename,sep=";",header=TRUE)
date=as.Date(tab$Date,"%d/%m/%Y")

l_u=unique(tab$Lieu_id)

for (l in 1:length(l_u)){
                tab_flore=subset(tab,((Lieu_id==l_u[l])&(Res_code %in% c("FLORPAR","FLORIND"))),select=c("Lieu_id","Lieu_libel","Date","Imm","Taxon","Val"))
                if (dim(tab_flore)[1]!=0){
                        tab_flore$Date=as.Date(tab_flore$Date,"%d/%m/%Y")
                        date_flore=sort(unique(as.Date(tab_flore$Date,"%d/%m/%Y")))
			plot(date_flore,log10(aggregate(tab_flore$Val,by=list(tab_flore$Date),FUN=sum)$x+1),main=as.character(tab_flore$Lieu_libel[1]),t="p",col="black",pch=16,ylab="Log10 Abondance",xlab="Date")
			axis.Date(1, at=seq(from=as.Date("1987-01-01"),to=as.Date("2016-01-01"),by="year"), format="%Y")
		}
}
}
dev.off()
