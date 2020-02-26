#Builds the files per correspondance

#WARNING: TO DO: check if there are some sites that changed names? (see exploitation_FLORTOT)

options(warn=2)

graphics.off()
rm(list=ls())

matrix_classif=function(file_corres,which_taxo){
	require("stringr")
	load("./data/taxonomy/taxonomy_pencen.RData") #final_pencen =c("Name_REPHY","Phylum","Subphylum","Infraphylum","Class","Subclass","Order","Family","Genus")
	tab_corres=read.table(file_corres,header=TRUE,sep=";")
	matrix_tmp=matrix(0,nrow=dim(final_pencen)[1],ncol=1+dim(tab_corres)[1]+1)
	colnames(matrix_tmp)=c("Name_REPHY",as.character(tab_corres$Code),"NEI")
	matrix_tmp[,"NEI"]=rep(1,dim(matrix_tmp)[1])
	matrix_tmp[,"Name_REPHY"]=final_pencen[,"Name_REPHY"]
	for (i in 1:dim(tab_corres)[1]){
		group=unlist(strsplit(as.character(tab_corres[i,which_taxo]),split="+",fixed=TRUE))
		for (g in 1:length(group)){
			id=c()
			if(group[g]=="Gymnodiniaceae"){
				id=which(grepl(group[g],final_pencen[,"Family"])&is.na(final_pencen[,'Genus']))
			}else if(grepl("Cry",group[g])){
                                id=which(grepl("Cryptophyta",final_pencen[,"Phylum"])&is.na(final_pencen[,'Genus']))
                        }else if(grepl("Eug",group[g])){
                                id=which(grepl("Euglenozoa",final_pencen[,"Phylum"])&is.na(final_pencen[,'Genus']))
                        }else{
				id=grep(group[g],final_pencen[,which_taxo])
			}
			if(length(id)>0){
				matrix_tmp[id,1+i]=rep(1,length(id))
				matrix_tmp[id,"NEI"]=rep(0,length(id))
			}
		}
	}
	return(matrix_tmp)
}

ff=c('./data/taxonomy/corres_hernandez')
#ff=c('corres_Arcachon','corres_hernandez')
for (file_corres in ff){
matrix_corres=matrix_classif(paste(file_corres,'.csv',sep=""),"Genus")

sites=c("Men er Roue","Loscolo","Croisic","LEperon","Cornard","Auger","Antoine","Lazaret")
sites=c("Antoine")

filename_tot=c("data/raw/Q2_170418_site_Tania_UTF8_only_good.csv")#,"Q2_170511_SRN_Chnord_UTF8_only_good.csv")

#file.remove(paste("data/",list.files(pattern='*.flortot_only.csv'),sep=""))
for (filename in filename_tot){
	tab=read.table(filename,sep=";",header=TRUE)
	date=as.Date(tab$Date,"%d/%m/%Y")
	d_min=min(date)
	d_max=max(date)

	for (l in 1:length(sites)){
                tab_flore=subset(tab,((grepl(sites[l],Lieu_libel))&(Res_code=="FLORTOT")),select=c("Lieu_id","Lieu_libel","Date","Imm","Taxon","Val"))
                if (dim(tab_flore)[1]!=0){
                        tab_flore$Date=as.Date(tab_flore$Date,"%d/%m/%Y")
			if(sites[l]=="Antoine"){
                		tab_flore_bis=subset(tab,((grepl("Carteau",Lieu_libel))&(Res_code=="FLORTOT")),select=c("Lieu_id","Lieu_libel","Date","Imm","Taxon","Val"))
				d1=max(tab_flore$Date,na.rm=TRUE)
				tab_flore_bis$Date=as.Date(tab_flore_bis$Date, "%d/%m/%Y")
				tab_flore_bis=subset(tab_flore_bis,Date>d1)
				tab_flore=rbind(tab_flore,tab_flore_bis)
			}
			date_site=sort(unique(tab_flore$Date))
			tab_new=matrix(0,nrow=length(date_site),ncol=dim(matrix_corres)[2])
			colnames(tab_new)=c("Date",colnames(matrix_corres)[2:dim(matrix_corres)[2]])
			for (d in 1:length(date_site)){
				id_d=which(tab_flore$Date==date_site[d])
				if(length(unique(tab_flore$Imm[id_d]))>1){ #On prend en compte plusieurs immersions, c'est à dire surface et fond
					plop=which(tab_flore$Date==date_site[d]&tab_flore$Imm>min(tab_flore$Imm[id_d],na.rm=TRUE))
					if(length(plop)>0){
					tab_flore=tab_flore[-plop,]
					}
				}
				id_d=which(tab_flore$Date==date_site[d])
				taxon=tab_flore$Taxon[id_d]
				for(t in taxon){
					id=which(matrix_corres[,1]==t)
					tmp=tab_flore[(tab_flore$Date==date_site[d])&(tab_flore$Taxon==t),]
					if(length(tmp$Val)>1){
						if(length(tmp$Val)>1){
							tmp$Val[1]=mean(tmp$Val)
						}
					}
					tab_new[d,2:dim(tab_new)[2]]=tab_new[d,2:dim(tab_new)[2]]+as.numeric(matrix_corres[id,2:dim(matrix_corres)[2]])*tmp$Val[1]
				}
				tab_new[d,tab_new[d,]==0]=rep(NA,sum(tab_new[d,]==0))
			}
			tab_new[,'Date']=as.character(date_site)
		}
	write.table(tab_new,file=paste('data/',file_corres,'_',sites[l],'.txt',sep=''),sep=";",na="NA",col.names=TRUE,row.names=FALSE)
	}
}
}
