graphics.off()
rm(list=ls())


filename=c("data/raw/Q2_170418_site_Tania_UTF8_only_good.csv","data/raw/Q2_170511_SRN_Chnord_UTF8_only_good.csv")

tab_REPHY=read.table(filename[1],sep=";",header=TRUE)
tab_SRN=read.table(filename[2],sep=";",header=TRUE)

id=intersect(unique(tab_REPHY$Lieu_id),unique(tab_SRN$Lieu_id))

tab_REPHY=subset(tab_REPHY,((Lieu_id %in% id)&(Res_code=="FLORTOT")),select=c("Lieu_id","Lieu_libel","Date","Taxon","Val"))
tab_SRN=subset(tab_SRN,((Lieu_id %in% id)&(Res_code=="FLORTOT")),select=c("Lieu_id","Lieu_libel","Date","Taxon","Val"))

tab_REPHY$Date=as.Date(tab_REPHY$Date,'%d/%m/%Y')
tab_SRN$Date=as.Date(tab_SRN$Date,'%d/%m/%Y')

#Comparaison visuelle
pdf("graphe/comparaison_REPHY_SRN.pdf",width=15,height=15)
par(mfrow=c(2,2))
for (i in 1:length(id)){
	t1=tapply(tab_REPHY$Val[tab_REPHY$Lieu_id==id[i]],tab_REPHY$Date[tab_REPHY$Lieu_id==id[i]],sum,na.rm=TRUE)
	t2=tapply(tab_SRN$Val[tab_SRN$Lieu_id==id[i]],tab_SRN$Date[tab_SRN$Lieu_id==id[i]],sum,na.rm=TRUE)

	d1=min(as.Date(c(rownames(t1),rownames(t2))))
	d2=max(as.Date(c(rownames(t1),rownames(t2))))
	plot(as.Date(rownames(t1)),t1,pch=16,cex=2,col="red",main=tab_REPHY$Lieu_libel[tab_REPHY$Lieu_id==id[i]][1],xlab="",ylab="Abundance totale",xlim=c(d1,d2))
	points(as.Date(rownames(t2)),t2,pch=16,cex=1,col="black")
	legend("topleft",c("REPHY","SRN"),col=c("red","black"),pch=16)
}
dev.off()

#Comparaison de date uniquement
results_common=vector("list",length(id))
results_left_out=vector("list",length(id))
for (i in 1:length(id)){
	d1=unique(as.Date(tab_REPHY$Date[tab_REPHY$Lieu_id==id[i]]))
	d2=unique(as.Date(tab_SRN$Date[tab_SRN$Lieu_id==id[i]]))

	d3=as.Date(intersect(d1,d2),origin="1970-01-01")

	print(paste("Pour",tab_REPHY$Lieu_libel[tab_REPHY$Lieu_id==id[i]][1],'il y a',length(d1),'dans REPHY',length(d2),'dans SRN',length(d3),'dates communes',sep=" "))

	t1=subset(tab_REPHY,((Lieu_id==id[i])&(Date %in% d3)),select=c("Date","Taxon","Val"))
	t2=subset(tab_SRN,((Lieu_id==id[i])&(Date %in% d3)),select=c("Date","Taxon","Val"))

	or=order(as.Date(t1$Date),t1$Taxon)
	t1=t1[or,]
	or=order(as.Date(t2$Date),t2$Taxon)
	t2=t2[or,]
	
	diff1=c()
	diff2=c()
	if(dim(t1)[1]!=dim(t2)[1]){
		for(d in 1:length(d3)){
			taxa1=t1$Taxon[t1$Date==d3[d]]
			taxa2=t2$Taxon[t2$Date==d3[d]]
			diff1=c(diff1,which(t1$Date==d3[d]&t1$Taxon %in% setdiff(taxa1,taxa2)))
			diff2=c(diff2,which(t2$Date==d3[d]&t2$Taxon %in% setdiff(taxa2,taxa1)))
		}
		results_left_out[[i]]=list("REPHY_spec"=t1[diff1,],"SRN_spec"=t2[diff2,])
		t1=t1[-diff1,]
		t2=t2[-diff2,]
	}else{
		results_left_out[[i]]=NULL
	}
	t3=t1
	t3$Val=t1$Val-t2$Val
	t3=t3[t3$Val!=0,]
	results_common[[i]]=t3
}

