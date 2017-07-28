#Buils the subclass files

#WARNING: TO DO: check if there are some sites that changed names? (see exploitation_FLORTOT)

rm(list=ls())
graphics.off()

unique_element=function(x){
        require(stringr)
#A tuned version of unique in order to remove elements "containing" each other :
       # x=unique(unique(str_replace(x,"incertae sedis","")))
        x=unique(x)
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
file_corres="subclass"
subclass=unique_element(final_pencen[!is.na(final_pencen[,'Subclass']),'Subclass'])
id=grep("Euglenia",subclass)
subclass=subclass[-id]

sites=c("Men er Roue","Loscolo","Croisic","LEperon","Cornard","Auger","Antoine","Lazaret")

filename_tot=c("data/raw/Q2_170418_site_Tania_UTF8_only_good.csv")#,"Q2_170511_SRN_Chnord_UTF8_only_good.csv")

for (filename in filename_tot){
        tab=read.table(filename,sep=";",header=TRUE)
        date=as.Date(tab$Date,"%d/%m/%Y")
        d_min=min(date)
        d_max=max(date)

        for (l in 1:length(sites)){
                tab_flore=subset(tab,((grepl(sites[l],Lieu_libel))&(Res_code=="FLORTOT")),select=c("Lieu_id","Lieu_libel","Date","Imm","Taxon","Val"))
                if (dim(tab_flore)[1]!=0){
                        tab_flore$Date=as.Date(tab_flore$Date,"%d/%m/%Y")
                        date_site=sort(unique(tab_flore$Date))
                        tab_new=matrix(0,nrow=length(date_site),ncol=length(subclass)+1+4)
                        colnames(tab_new)=c("Date",subclass,"Dinophyceae","EUG","CRY","NEI")
                        for (d in 1:length(date_site)){
                                id_d=which(tab_flore$Date==date_site[d])
                                if(length(unique(tab_flore$Imm[id_d]))>1){ #On prend en compte plusieurs immersions, c'est à dire surface et fond
                                        plop=which(tab_flore$Date==date_site[d]&tab_flore$Imm>min(tab_flore$Imm[id_d],na.rm=TRUE))
                                        tab_flore=tab_flore[-plop,]
                                }
                                id_d=which(tab_flore$Date==date_site[d])
                                taxon=tab_flore$Taxon[id_d]
                                for(t in taxon){
                                        tmp=tab_flore[(tab_flore$Date==date_site[d])&(tab_flore$Taxon==t),]
                                        if(length(tmp$Val)>1){
                                                if(length(tmp$Val)>1){
                                                        tmp$Val[1]=mean(tmp$Val)
                                                }
                                        }
					if(final_pencen[final_pencen[,'Name_REPHY']==t,"Class"] %in% c("Dinophyceae")){
						tab_new[d,"Dinophyceae"]=tab_new[d,"Dinophyceae"]+tmp$Val[1]
					}else if(final_pencen[final_pencen[,'Name_REPHY']==t,"Phylum"] %in% c("Euglenozoa")){
						tab_new[d,"EUG"]=tab_new[d,"EUG"]+tmp$Val[1]
					}else if(final_pencen[final_pencen[,'Name_REPHY']==t,"Phylum"] %in% c("Cryptophyta")){
						tab_new[d,"CRY"]=tab_new[d,"CRY"]+tmp$Val[1]
					}else if(sum(!is.na(grep(final_pencen[final_pencen[,'Name_REPHY']==t,"Subclass"],subclass)))>0){
						id=grep(final_pencen[final_pencen[,'Name_REPHY']==t,"Subclass"],subclass)
                                       	 	tab_new[d,id+1]=tab_new[d,id+1]+tmp$Val[1]
					}else{
						tab_new[d,"NEI"]=tab_new[d,"NEI"]+tmp$Val[1]

					}
                                }
                                tab_new[d,tab_new[d,]==0]=rep(NA,sum(tab_new[d,]==0))
                        }
                        tab_new[,'Date']=as.character(date_site)
                }
        write.table(tab_new,file=paste(file_corres,'_',sites[l],'.txt',sep=''),sep=";",na="NA",col.names=TRUE,row.names=FALSE)
        }
}

