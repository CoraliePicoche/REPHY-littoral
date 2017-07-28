#Build the database for the correspondance between WORMS and REPHY

rm(list=ls())
graphics.off()

library("taxize")
library("stringr")

phyla_ok=c("Cyanobacteria","Glaucophyta","Chlorophyta","Euglenozoa","Cryptophyta","Chrysophyta","Bacillariophyta","Haptophyta","Ciliophora","Ochrophyta","Myzozoa","Kingdom_protista")
#Not in WORMS but in Reynolds 2006 : prochlorobacteria, anoxyphotobacteria,prasinophyta are in chlorophyta phylum, raphidophyta does not exist, xanthophyta are for arthropoda, eustigmatophyta are in Ochrophyta. After Haptophyta, I added ochrophyta, myzozoa for dinoflagellates

filename_tot=c("data/raw/Q2_170418_other_UTF8_only_good.csv","data/raw/Q2_170418_site_Tania_UTF8_only_good.csv","data/raw/Q2_170511_SRN_Chnord_UTF8_only_good.csv")

total_taxa=c()
for (filename in filename_tot){
	tab=read.table(filename,sep=";",header=TRUE)
	
	total_taxa=c(total_taxa,unlist(as.character(unique(tab$Taxon[tab$Res_code=="FLORTOT"]))))
}
total_taxa=unique(total_taxa)

final_pencen=matrix(NA,nrow=length(total_taxa),ncol=8+1)
colnames(final_pencen)=c("Name_REPHY","Phylum","Subphylum","Infraphylum","Class","Subclass","Order","Family","Genus")

for (t in 1:length(total_taxa)){
#for (t in 1:10){
	final_pencen[t,1]=total_taxa[t]

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
        }else if(grepl("Hermesinum",final_pencen[t,1])|grepl("Ebria",final_pencen[t,1])|grepl("Protoctista",final_pencen[t,1])|grepl("Protista",final_pencen[t,1])){
        	final_pencen[t,"Phylum"]="Kingdom_protista"
        }else{
		ta_2=strsplit(total_taxa[t],split=",")
        #Then, take all the things between + and keep only those beginning with a upper case
 	   	ta_3= lapply(ta_2,function(x) strsplit(x[[1]][1],split="+",fixed=TRUE))
        	ta_4=list(0)
                for (i in 1:length(ta_3)){
                	ta_4[[i]]=lapply(ta_3[[i]],function(x) x[grepl(str_trim(x),pattern="^[A-Z]")])
                        tmp=unlist(strsplit(unlist(ta_4[[i]][[1]])," "))
                        ta_4[[i]]=lapply(tmp[grep(tmp,pattern="^[A-Z]")],function(x) sub("[ ,]","",x))
                }
                        
			#So, now we have a nice list of genera/class/anything above, and we want from phylum to family if available
		list_problem=c()
                for (j in 1:length(ta_4)){
                                #Warning, better to make sure that every genus correspond to the same Class
                                tryCatch({
				id=get_wormsid(unlist(ta_4[j]),accepted=TRUE,rows=1) #WARNING : we may have strange answers : Nitzschia is not a platelminth!
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
                               					id=get_wormsid(unlist(ta_4[j][[1]][[c]]),accepted=TRUE,rows=i) #WARNING : we may have strange answers : Nitzschia is not a platelminth!
								tmp=classification(id,db="worms")
	                                                        a=tmp[[1]]$name[grep(colnames(final_pencen)[level],tmp[[1]]$rank)]
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
					if(length(taxo)>0){
						taxo=paste(taxo,sep="+",collapse="+")
						final_pencen[t,level]=taxo
					}
				}
			},error=function(e){print("Error")})
				#},error=function(e){stop()})
		}

			#Handling of values that are not present in WORMS
	}
}

write.table(final_pencen,"db_taxonomy_REPHY.csv",sep=";",row.names=FALSE,col.names=TRUE,append=FALSE)
save(final_pencen,file="taxonomy_pencen.RData")

plot(0,0,t="n",xlim=c(0,dim(final_pencen)[2]),ylim=c(0,dim(final_pencen)[1]))
for (l in 2:dim(final_pencen)[2]){
	points(l,sum(is.na(final_pencen[,l])))
}
