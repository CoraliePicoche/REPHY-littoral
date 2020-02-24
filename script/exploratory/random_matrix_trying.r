#2018/07/09 CP
#This script tests the stability of random matrix in different scenarios, according to the relation between interaction competitions

rm(list=ls())
graphics.off()
library(MARSS)
source("script/matrix_MAR_clean.r")

option_model="unconstrained"
groupe=c("BZ","MO","AR","SU")

methods=c("intra","mean_abs","std")

name_species=c("AST","CHA","DIT","GUI","LEP","NIT","PLE","PSE","RHI","SKE","THP","THL","GYM","PRO","PRP","SCR","CRY","EUG")
id_lieu=0


#First, we define the relations we are going to use
tab_vulnerability=array(NA,dim=c(10,length(option_model),length(name_species),length(methods)),dimnames=list(rep("",10),option_model,name_species,methods))
tab_generality=array(NA,dim=c(10,length(option_model),length(name_species),length(methods)),dimnames=list(rep("",10),option_model,name_species,methods))

for (g in groupe){
        if(g=="BZ"){
                option_lieu=c("Men er Roue","Loscolo","Croisic")
        }else if(g=="MO"){
                option_lieu=c("LEperon","Cornard","Auger")
        }else if(g=="SU"){
                option_lieu=c("Antoine","Lazaret")
        }else if(g=="AR"){
                option_lieu=c("Teychan","B7")
        }
	for (ll in 1:length(option_lieu)){
		id_lieu=id_lieu+1
		dimnames(tab_vulnerability)[[1]][id_lieu]=option_lieu[ll]
		dimnames(tab_generality)[[1]][id_lieu]=option_lieu[ll]
        	for (m in 1:length(option_model)){
                        f1=paste("data/analyse_MAR/",g,"/site_specific/",option_lieu[ll],"_",option_model[m],"_null_regular_common_",g,".RData",sep="")
                        load(f1)
                        sp=dimnames(cis$model$data)[[1]] #Studied species
			B=clean_matrix(cis)
			for (i in 1:length(sp)){
             			tab_vulnerability[id_lieu,option_model[m],sp[i],"intra"]=B[i,i]
	             		tab_vulnerability[id_lieu,option_model[m],sp[i],"mean_abs"]=mean(abs(B[i,B[i,]!=0]))
	             		tab_vulnerability[id_lieu,option_model[m],sp[i],"std"]=sd(abs(B[i,B[i,]!=0]))
                	        tab_generality[id_lieu,option_model[m],sp[i],"mean_abs"]=mean(abs(B[B[,i]!=0,i]))
	             		tab_generality[id_lieu,option_model[m],sp[i],"std"]=sd(abs(B[B[,i]!=0,i]))
                                tab_generality[id_lieu,option_model[m],sp[i],"intra"]=B[i,i]
			}
		}
	}
}

for(m in 1:length(option_model)){
	for(i in 1:2){
	#Set-up
	if(i==1){
		tab_tmp=tab_vulnerability
		atitle="Vulnerability"
	}else{
		tab_tmp=tab_generality
		atitle="Generality"
	}

	xlimi=c(min(c(tab_tmp[,option_model[m],,"intra"]),na.rm=T),max(c(tab_tmp[,option_model[m],,"intra"]),na.rm=T))
        ylimi=c(min(c(tab_tmp[,option_model[m],,"mean_abs"]),na.rm=T),max(c(tab_tmp[,option_model[m],,"mean_abs"]),na.rm=T))
	x1=c(tab_tmp[,option_model[m],,"intra"])
	y1=c(tab_tmp[,option_model[m],,"mean_abs"])
	lm1=lm(y1~x1)
	}
}

#Second : we build random matrices according to this scenario
sp_random=10
random_B=matrix(NA,nrow=sp_random,ncol=sp_random)
diag(random_B)=sample(x1,sp_random)
