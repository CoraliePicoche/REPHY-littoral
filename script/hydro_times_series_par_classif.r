graphics.off()
rm(list=ls())

sites=c("Men er Roue","Loscolo","Croisic","LEperon","Cornard","Auger","Antoine","Lazaret")
sites="Antoine"

filename_tot=c("data/raw/Q2_170418_site_Tania_UTF8_only_good.csv")#,"Q2_170511_SRN_Chnord_UTF8_only_good.csv")

#file.remove(paste("data/",list.files(pattern='*.flortot_only.csv'),sep=""))
for (filename in filename_tot){
        tab=read.table(filename,sep=";",header=TRUE)
        date=as.Date(tab$Date,"%d/%m/%Y")
        d_min=min(date)
        d_max=max(date)

        for (l in 1:length(sites)){
		print(sites[l])
                tab_cov=subset(tab,grepl(sites[l],Lieu_libel),select=c("Lieu_id","Lieu_libel","Date","Imm","Res_code","Val"))
		liste_hydro=unique(tab_cov$Res_code)
		liste_hydro=as.character(liste_hydro[-which(grepl('FLOR',liste_hydro))])
                if (dim(tab_cov)[1]!=0){
                        tab_cov$Date=as.Date(tab_cov$Date,"%d/%m/%Y")
                        if(sites[l]=="Antoine"){
                                tab_cov_bis=subset(tab,((grepl("Carteau",Lieu_libel))),select=c("Lieu_id","Lieu_libel","Date","Imm","Res_code","Val"))
                                d1=max(tab_cov$Date,na.rm=TRUE)
                                tab_cov_bis$Date=as.Date(tab_cov_bis$Date, "%d/%m/%Y")
                                tab_cov_bis=subset(tab_cov_bis,Date>d1)
                                tab_cov=rbind(tab_cov,tab_cov_bis)
                        }

                        date_site=sort(unique(tab_cov$Date))
                        tab_new=matrix(0,nrow=length(date_site),ncol=1+length(liste_hydro))
                        colnames(tab_new)=c("Date",liste_hydro)
                        for (d in 1:length(date_site)){
                                id_d=which(tab_cov$Date==date_site[d])
                                if(length(unique(tab_cov$Imm[id_d]))>1){ #On prend en compte plusieurs immersions, c'est à dire surface et fond
                                        plop=which(tab_cov$Date==date_site[d]&tab_cov$Imm>min(tab_cov$Imm[id_d],na.rm=TRUE))
                                        if(length(plop)>0){
					tab_cov=tab_cov[-plop,]
					}
                                }
                                for(h in 1:length(liste_hydro)){
                                	id_d=which(tab_cov$Date==date_site[d]&tab_cov$Res_code==liste_hydro[h])
					if(length(id_d)==0){
						tab_new[d,1+h]=NA
					}else{
                                       		tab_new[d,1+h]=mean(tab_cov$Val[id_d],na.rm=TRUE)
					}
                                }
                        }
                        tab_new[,'Date']=as.character(date_site)
                }
        write.table(tab_new,file=paste('data/',sites[l],'hydro.txt',sep=''),sep=";",na="NA",col.names=TRUE,row.names=FALSE)
        }
}

