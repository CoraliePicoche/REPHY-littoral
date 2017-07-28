#Number of genus/classes/families... for each site

#WARNING: TO DO: check if there are some sites that changed names? (see exploitation_FLORTOT)

graphics.off()
rm(list=ls())

unique_element=function(x){
	require(stringr)
#A tuned version of unique in order to remove elements "containing" each other :
	x=unique(unique(str_replace(x,"incertae sedis","")))
	tmp=x
	for (f in 1:length(x)){
		if(sum(grepl(x[f],tmp))>1){
			id=grep(x[f],tmp)
			ll=max(str_count(x[id]),na.rm=TRUE)
			for (i in id){
				if(str_count(x[i])<ll){
					tmp=tmp[-i]
				}
			}
		}
	}
	return(tmp)
}

load("taxonomy_pencen.RData") #final_pencen =c("Name_REPHY","Phylum","Subphylum","Infraphylum","Class","Subclass","Order","Family","Genus")
xi=2:dim(final_pencen)[2]-1
xi_label=dimnames(final_pencen)[[2]][2:dim(final_pencen)[2]]

sites=c("Men er Roue","Loscolo","Croisic","LEperon","Cornard","Auger","Antoine","Lazaret")

filename_tot=c("data/raw/Q2_170418_site_Tania_UTF8_only_good.csv")#,"Q2_170511_SRN_Chnord_UTF8_only_good.csv")

#file.remove(paste("data/",list.files(pattern='*.flortot_only.csv'),sep=""))
for (filename in filename_tot){
tab=read.table(filename,sep=";",header=TRUE)
date=as.Date(tab$Date,"%d/%m/%Y")
d_min=min(date)
d_max=max(date)

#On peut regarder la répartition en terme de nombre de lignes utilisable...
pdf("graphe/repartition_taxonomy_precise.pdf",width=15)
par(mfrow=c(3,3),mar=c(2,2,2,0.5),oma=c(1,2,1,0.5))
for (l in 1:length(sites)){
                tab_flore=subset(tab,((grepl(sites[l],Lieu_libel))&(Res_code=="FLORTOT")),select=c("Lieu_id","Lieu_libel","Date","Imm","Taxon","Val"))
                if (dim(tab_flore)[1]!=0){
	                tab_flore$Date=as.Date(tab_flore$Date,"%d/%m/%Y")
			tax=unique(tab_flore$Taxon)
			print(sites[l])
			print(paste("Nb dates",length(unique(tab_flore$Date)),"Nb taxa",length(tax),"Nb observ",dim(tab_flore)[1],sep=" "))
			plot(0,0,xlim=c(0,dim(final_pencen)[2]),ylim=c(0.8,1.0),t="n",xaxt="n",xlab="",ylab="",yaxt="n",main=tab_flore$Lieu_libel[1],cex.main=2)
			txt_leg=c()
			for (l_t in 2:dim(final_pencen)[2]){
				count=0
				tmp=c()
				for (t in tax){
					in_db=!(is.na(final_pencen[final_pencen[,1]==t,l_t]))
					tmp=c(tmp,final_pencen[final_pencen[,1]==t,l_t])
					count=count+sum(tab_flore$Taxon==t)*in_db[1]
				}
				tmp=unique_element(tmp[!is.na(tmp)])
				text(x=l_t-1,y=count/length(tab_flore$Taxon)*1.0,labels=as.character(length(tmp)),cex=1.7)
				if(l_t==4||l_t==6){
					txt_leg=c(txt_leg,paste(format(round(count/length(tab_flore$Taxon)*1.0,2),nsmall=2),"for",length(tmp),xi_label[l_t-1],sep=" "))
				}
				
			}
			legend("bottomleft",txt_leg,pch=NA,bty="y",cex=1.4)
			if(l>6){
			axlab="Taxo level"
			axis(1,at=xi,label=c("Phylum","Sub","Infra","Class","Sub","Order","Fam.","Genus"),cex.axis=1.3)
			}else{
			axlab=""
			axis(1,at=xi,lab=NA,cex=1.5)
			}
			if(l%%3==1){
			aylab="Nb observations"
			axis(2,at=seq(0.2,1.0,0.05),lab=seq(0.2,1.0,0.05),cex.axis=1.5)
			}else{
			aylab=""
			axis(2,at=seq(0.2,1.0,0.05),lab=NA,cex=1.5)
			}
		}
	}
dev.off()

#Mais aussi en terme d'abondance
pdf("graphe/repartition_taxonomy_abondance_precise.pdf",width=15)
par(mfrow=c(3,3),mar=c(2,2,2,0.5),oma=c(1,2,1,0.5))
for (l in 1:length(sites)){
                tab_flore=subset(tab,((grepl(sites[l],Lieu_libel))&(Res_code=="FLORTOT")),select=c("Lieu_id","Lieu_libel","Date","Imm","Taxon","Val"))
                if (dim(tab_flore)[1]!=0){
			tot=sum(tab_flore$Val,na.rm=TRUE)
                        tab_flore$Date=as.Date(tab_flore$Date,"%d/%m/%Y")
                        tax=unique(tab_flore$Taxon)
                        plot(0,0,xlim=c(0,dim(final_pencen)[2]),ylim=c(0.8,1.0),t="n",xaxt="n",xlab="",ylab="",yaxt="n",main=tab_flore$Lieu_libel[1],cex.main=2)
			txt_leg=c()
                        for (l_t in 2:dim(final_pencen)[2]){
                                count=0
                                tmp=c()
                                for (t in tax){
                                        in_db=!(is.na(final_pencen[final_pencen[,1]==t,l_t]))
                                        tmp=c(tmp,final_pencen[final_pencen[,1]==t,l_t])
                                        count=count+sum(tab_flore$Val[tab_flore$Taxon==t])*in_db[1]
                                }
                                tmp=unique_element(tmp[!is.na(tmp)])
                                text(x=l_t-1,y=count/tot*1.0,labels=as.character(length(tmp)),cex=1.7)
				if(l_t==4||l_t==6){
					txt_leg=c(txt_leg,paste(format(round(count/tot*1.0,2),nsmall=2),"for",length(tmp),xi_label[l_t-1],sep=" "))
				}
                        }
			legend("bottomleft",txt_leg,pch=NA,bty="y",cex=1.4)
                        if(l>6){
                        axlab="Taxo level"
                        axis(1,at=xi,label=c("Phylum","Sub","Infra","Class","Sub","Order","Fami.","Genus"),cex.axis=1.3)
                        }else{
                        axlab=""
                        axis(1,at=xi,lab=NA)
                        }
                        if(l%%3==1){
                        aylab="Nb observations"
                        axis(2,at=seq(0.2,1.0,0.05),lab=seq(0.2,1.0,0.05),cex.axis=1.5)
                        }else{
                        aylab=""
                        axis(2,at=seq(0.2,1.0,0.05),lab=NA)
                        }
                }
        }
dev.off()
}
