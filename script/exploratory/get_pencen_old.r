rm(list=ls())
graphics.off()
library("stringr")

f_atlas=read.csv("/home/cpicoche/Documents/Biblio/Plankton/rephy_flore_totale_atlas.csv",sep=";",header=TRUE)

#filename_tot=c("Q2_170418_other_UTF8_only_good.csv","Q2_170418_site_Tania_UTF8_only_good.csv","Q2_170511_SRN_Chnord_UTF8_only_good.csv")
filename_tot=c("Q2_170418_site_Tania_UTF8_only_good.csv")

for (filename in filename_tot){
tab=read.table(filename,sep=";",header=TRUE)
date=as.Date(tab$Date,"%d/%m/%Y")

l_u=unique(tab$Lieu_id)
d_min=min(date)
d_max=max(date)

#for (l in 1:length(l_u)){
for (l in 1:1){
                tab_flore=subset(tab,((Lieu_id==l_u[l])&(Res_code=="FLORTOT")),select=c("Lieu_id","Lieu_libel","Date","Imm","Taxon","Val"))
                if (dim(tab_flore)[1]!=0){
			ta=as.character(unique(tab_flore$Taxon))
			final_pencen=rep('NA',length(ta))
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
			#Finally, search in Genera
			print("search in Genus")
			list_problem=c()
			for (t in 1:length(ta_4)){
				#Warning, better to make sure that every genus correspond to the same Class
				a=as.character(unique(unlist(lapply(ta_4[[t]],function(x) f_atlas$Classe[grep(x,f_atlas$Genre)]))))
				if(length(a)==1){
					final_pencen[t]=a
				}else if(length(a)==0){
					list_problem=c(list_problem,t)
				}else if(length(a)>1){
					print(ta_4[[t]])
					print(a)
				}
			}

			#If not found, search in Class
			print("search in Class")
                        list_problem_2=c()
                        for (t in list_problem){
				a=unlist(ta_4[[t]]) #Assuming there's only one element
                                if(length(a)==1){
					if(a %in% f_atlas$Classe){
                                        	final_pencen[t]=a
                                	}else{
 	                                     	list_problem_2=c(list_problem_2,t)
					}
				}else{
 	                                     	list_problem_2=c(list_problem_2,t)
				}
                        }

			#If not found, search in family or order
			print("search in Family or Order")
                        list_problem_3=c()
			for (t in list_problem_2){
				a=as.character(unique(unlist(lapply(ta_4[[t]],function(x) f_atlas$Classe[grep(x,f_atlas$Family)]))))
				if(length(a)==0){
					a=as.character(unique(unlist(lapply(ta_4[[t]],function(x) f_atlas$Classe[grep(x,f_atlas$Order)]))))
				}	
				if(length(a)==0){
                        		list_problem_3=c(list_problem_3,t)
				}else{
					if(length(a)==1){
						final_pencen[t]=a
					}else{
						print(ta_4[[t]])
						print(a)
					}
				}
			}

			#If not found, raise an alert
			print('still problems')
			for (t in list_problem_3){
				print(ta_4[[t]])
			}
			
		}

}
}
