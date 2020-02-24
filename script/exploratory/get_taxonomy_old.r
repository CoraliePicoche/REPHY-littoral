rm(list=ls())
graphics.off()

library("taxize")
library("stringr")

#filename_tot=c("Q2_170418_other_UTF8_only_good.csv","Q2_170418_site_Tania_UTF8_only_good.csv","Q2_170511_SRN_Chnord_UTF8_only_good.csv")
filename_tot=c("Q2_170418_site_Tania_UTF8_only_good.csv")

for (filename in filename_tot){
tab=read.table(filename,sep=";",header=TRUE)
date=as.Date(tab$Date,"%d/%m/%Y")

l_u=unique(tab$Lieu_id)
d_min=min(date)
d_max=max(date)

phyla_ok=c("Cyanobacteria","Glaucophyta","Chlorophyta","Euglenozoa","Cryptophyta","Chrysophyta","Bacillariophyta","Haptophyta","Ciliophora","Ochrophyta","Myzozoa","Kingdom_protista")
#Not in WORMS but in Reynolds 2006 : prochlorobacteria, anoxyphotobacteria,prasinophyta are in chlorophyta phylum, raphidophyta does not exist, xanthophyta are for arthropoda, eustigmatophyta are in Ochrophyta. After Haptophyta, I added ochrophyta, myzozoa for dinoflagellates

for (l in 1:length(l_u)){
#for (l in 1:1){
		print('******************')
		print(l_u[l])	
                tab_flore=subset(tab,((Lieu_id==l_u[l])&(Res_code=="FLORTOT")),select=c("Lieu_id","Lieu_libel","Date","Imm","Taxon","Val"))
                if (dim(tab_flore)[1]!=0){
                        ta=as.character(unique(tab_flore$Taxon))
                        final_pencen=matrix(NA,nrow=length(ta),ncol=8+1)
			colnames(final_pencen)=c("Name_REPHY","Phylum","Subphylum","Infraphylum","Class","Subclass","Order","Family","Genus")
			final_pencen[,1]=ta
			#For now, I want to see if I can get clear families, even though I have multiple genus
                        #First, only take what's before ',' which are there only to specify species
                        ta_2=strsplit(ta,split=",")
                        #Then, take all the things between + and keep only those beginning with a upper case
                        ta_3= lapply(ta,function(x) strsplit(x[[1]][1],split="+",fixed=TRUE))
                        ta_4=list(0)
                        for (t in 1:length(ta_3)){
                                ta_4[[t]]=lapply(ta_3[[t]],function(x) x[grepl(str_trim(x),pattern="^[A-Z]")])
                                tmp=unlist(strsplit(unlist(ta_4[[t]][[1]])," "))
                                ta_4[[t]]=lapply(tmp[grep(tmp,pattern="^[A-Z]")],function(x) sub("[ ,]","",x))
                        }
                        
			#So, now we have a nice list of genera/class/anything above, and we want from phylum to family if available
			list_problem=c()
                        for (t in 1:length(ta_4)){
#                        for (t in 90:100){
                                #Warning, better to make sure that every genus correspond to the same Class
                                tryCatch({
				id=get_wormsid(unlist(ta_4[t]),accepted=TRUE,rows=1) #WARNING : we may have strange answers : Nitzschia is not a platelminth!
				classif=classification(id,db="worms")
				for (level in 2:dim(final_pencen)[2]){
					taxo=c()
					nn=names(classif)
					for (c in 1:length(nn)){
						a=classif[[nn[c]]]$name[grep(colnames(final_pencen)[level],classif[[nn[c]]]$rank)]
						if(colnames(final_pencen)[level]=="Phylum"){
							right_phylum=(a %in% phyla_ok)
							i=1
							if(!right_phylum){
								classif[[nn[c]]]=NULL
							}
							while(!right_phylum){
								i=i+1
                               					id=get_wormsid(unlist(ta_4[t][[1]][[c]]),accepted=TRUE,rows=i) #WARNING : we may have strange answers : Nitzschia is not a platelminth!
								tmp=classification(id,db="worms")
	                                                        a=tmp[[1]]$name[grep(colnames(final_pencen)[level],tmp[[1]]$rank)]
								print(ta_4[t])
								print("?")
								print(a)
        	                                                right_phylum=(a %in% phyla_ok)
								if(right_phylum){
									classif[[id[1]]]=tmp[[1]]
								}
							}
							
						}
						if(length(a)>0){
							if (!(a %in% taxo)){
								taxo=c(taxo,a)
							}
						}
					}
					taxo=paste(taxo,sep="+",collapse="+")
					if(length(taxo)>0){
						final_pencen[t,level]=taxo
					}
				}
				},error=function(e){print("Error")})
				#},error=function(e){stop()})
			}

			#Handling of values that are not present in WORMS
			for(t in 1:dim(final_pencen)[1]){
				#Values that might have problems : Protoctista, Pennées, Choanofila, Ebria, Centriques, Phytoflagellés excepté dinoflagellés,Chromista,Hermesinum,"Tous Dinophysis ronds avec épithèque bien visible",Ebriaceae 
				#We decide to ignore Staurastrum and Closterium which appear only once for 100 cells, in only one site. All others appear at least at two different sites. 
				#I don't do anything with Protoctista+Protista, they are a different kingdom ; same thing for Ebria, Ebria tripartita and Ebriaceae, Phytoflagellés ; Chromista ; and Hermesinum. TO be noticed: Protista, Ebria and Hermesinum are all proctista, maybe it would be intereseting

#,"Phylum","Subphylum","Infraphylum","Class","Subclass","Order","Family")#,"Genus")

				na=sum(is.na(final_pencen[t,]))
				if(na==(dim(final_pencen)[2]-1)){#There was a problem
					if(grepl("Pennées",final_pencen[t,1])){
						final_pencen[t,"Phylum"]="Ochrophyta"
						final_pencen[t,"Subphylum"]="Khakista"
						final_pencen[t,"Class"]="Bacillariophyceae" 
						final_pencen[t,"Subclass"]="Bacillariophycidae+Fragilariophycidae" #We have no information on the presence or not of raphs
					}else if(grepl("Centriques",final_pencen[t,1])){
						final_pencen[t,"Phylum"]="Ochrophyta"
						final_pencen[t,"Subphylum"]="Khakista"
						final_pencen[t,"Class"]="Bacillariophyceae" 
						final_pencen[t,"Subclass"]="Coscinodiscophycidae" #Wondering if coscino is not only polar centrique and not multipolar centrique -> at least Chaetoceros and Skeletonma (and Lepto and Rhizo) 
					}else if(grepl("Choanofila",final_pencen[t,1])){ #Bug from WORMS, it does exist
						final_pencen[t,"Phylum"]="Choanozoa"
						final_pencen[t,"Subphylum"]="Choanofila"
					}else if(grepl("Dinophysis",final_pencen[t,1])){
						final_pencen[t,"Phylum"]="Myzozoa"
						final_pencen[t,"Subphylum"]="Dinozoa"
						final_pencen[t,"Infraphylum"]="Dinoflagellata"
						final_pencen[t,"Class"]="Dinophyceae" 
						final_pencen[t,"Order"]="Dinophysiales"
						final_pencen[t,"Family"]="Dinophysiaceae"
						final_pencen[t,"Genus"]="Dinophysis"
					}else if(grepl("Hermesinum",final_pencen[t,1])|grepl("Ebria",final_pence[t,1])|grepl("Protoctista",final_pencen[t,1])|grepl("Protista",final_pencen[t,1])){
						final_pencen[t,"Phylum"]="Kingdom_protista"
					}
				}
			}
		}
	}
}

