#2018/09/06 CP
#This script computes several metrics on the interaction matrices we study

rm(list=ls())
graphics.off()
library(MARSS)
library(stringr)

option_model=c("unconstrained","pencen")
#option_model="pencen"
option_NEI=c("null") #We don't take the "Not Elsewhere Identified" species into account
groupe=c("BZ","MO","AR","SU")
option_sp="common" #Species are the same 
take_0=TRUE #Do we consider the forced 0 in the interaction matrix when computing the mean values of inter-group competition?

methods=c("intra","mean_raw","mean_abs","summed_abs","summed_raw")

name_species=c("AST","CHA","DIT","GUI","LEP","NIT","PLE","PSE","RHI","SKE","THP","THL","GYM","PRO","PRP","SCR","CRY","EUG")

id_lieu=0

tab_vulnerability=array(NA,dim=c(10,length(option_model),length(name_species),length(methods)),dimnames=list(rep("",10),option_model,name_species,methods)) #10 subsites, 2 models, 21 species, 4 methods
tab_generality=array(NA,dim=c(10,length(option_model),length(name_species),length(methods)),dimnames=list(rep("",10),option_model,name_species,methods)) #10 subsites, 2 models, 21 species, 4 methods

colo=c()
for (g in groupe){
        if(g=="BZ"){
                option_lieu=c("Men er Roue","Loscolo","Croisic")
		colo=c(colo,rep("green",length(option_lieu)))
        }else if(g=="MO"){
                option_lieu=c("LEperon","Cornard","Auger")
		colo=c(colo,rep("darkblue",length(option_lieu)))
        }else if(g=="SU"){
                option_lieu=c("Antoine","Lazaret")
		colo=c(colo,rep("darkred",length(option_lieu)))
        }else if(g=="AR"){
                option_lieu=c("Teychan","B7")
		colo=c(colo,rep("cyan",length(option_lieu)))
        }
	for (ll in 1:length(option_lieu)){
		id_lieu=id_lieu+1
		dimnames(tab_vulnerability)[[1]][id_lieu]=option_lieu[ll]
		dimnames(tab_generality)[[1]][id_lieu]=option_lieu[ll]
        	for (m in 1:length(option_model)){
                        f1=paste("data/analyse_MAR/",g,"/site_specific/",option_lieu[ll],"_",option_model[m],"_",option_NEI,"_regular_common_",g,".RData",sep="")
                        load(f1)

                        sp=dimnames(cis$model$data)[[1]] #Studied species
                        var=dimnames(cis$call$model$c)[[1]] #Studied covariates
                        nom=dimnames(cis$par$B)[[1]] #Names of the interactions
                        #get the value from the MARSS object, according to the names of the coefficients
			B=diag(-1,length(sp),length(sp)) #I am heavily using the fact that we are using an unconstrained matrix
                        for (n in 1:length(nom)){
                                if(grepl("^\\(",nom[n])){ #default when using unconstrained or diagonal, names are of the form (1,2) for (i,j)
                                        i=as.numeric(strsplit(nom[n],split="[\\(,\\)]")[[1]][2])
                                        j=as.numeric(strsplit(nom[n],split="[\\(,\\)]")[[1]][3])
               				if(is.na(i)){ #i and j might not be numeric
			                        i=which(unlist(lapply(sp,function(x) grepl(x,strsplit(nom[n],split="[\\(,\\)]")[[1]][2]))))
                        			j=which(unlist(lapply(sp,function(x) grepl(x,strsplit(nom[n],split="[\\(,\\)]")[[1]][3]))))
                			}
      				  }else{
				 	a=str_split(nom[n],sp) #The coefficient is XY where X eats Y, X and Y being in the species list
                			for (abis in 1:length(a)){
                        			if(length(a[[abis]])==2){
                               				if((a[[abis]][1]=="")&&(a[[abis]][2]!="")){
                                        			j=abis
                                			}
                                			if((a[[abis]][2]=="")&&(a[[abis]][1]!="")){
                                        			i=abis
                                			}
                        			}else if(length(a[[abis]])==3){
                                			if((a[[abis]][2]=="")&&(a[[abis]][1]=="")&&(a[[abis]][3]=="")){
                                        			j=abis
                                        			i=abis
                                			}
                        			}
					}
				}
				B[i,j]=B[i,j]+cis$par$B[n]
			}
			B_nodiag=B
			diag(B_nodiag)=NA
			for (i in 1:length(sp)){
             			tab_vulnerability[id_lieu,option_model[m],sp[i],"intra"]=B[i,i]
				if(take_0){
             				tab_vulnerability[id_lieu,option_model[m],sp[i],"mean_raw"]=mean(B_nodiag[i,],na.rm=T)
	             			tab_vulnerability[id_lieu,option_model[m],sp[i],"mean_abs"]=mean(abs(B_nodiag[i,]),na.rm=T)
        	                        tab_generality[id_lieu,option_model[m],sp[i],"mean_raw"]=mean(B_nodiag[,i],na.rm=T)
                	                tab_generality[id_lieu,option_model[m],sp[i],"mean_abs"]=mean(abs(B_nodiag[,i]),na.rm=T)
				}else{
             				tab_vulnerability[id_lieu,option_model[m],sp[i],"mean_raw"]=mean(B_nodiag[i,B_nodiag[i,]!=0],na.rm=T)
	             			tab_vulnerability[id_lieu,option_model[m],sp[i],"mean_abs"]=mean(abs(B_nodiag[i,B_nodiag[i,]!=0]),na.rm=T)
        	                        tab_generality[id_lieu,option_model[m],sp[i],"mean_raw"]=mean(B_nodiag[B_nodiag[,i]!=0,i],na.rm=T)
                	                tab_generality[id_lieu,option_model[m],sp[i],"mean_abs"]=mean(abs(B_nodiag[B_nodiag[,i]!=0,i]),na.rm=T)
				}
             			tab_vulnerability[id_lieu,option_model[m],sp[i],"summed_abs"]=sum(abs(B_nodiag[i,]),na.rm=T)
             			tab_vulnerability[id_lieu,option_model[m],sp[i],"summed_raw"]=sum(B_nodiag[i,],na.rm=T)

                                tab_generality[id_lieu,option_model[m],sp[i],"intra"]=B[i,i]
                                tab_generality[id_lieu,option_model[m],sp[i],"summed_abs"]=sum(abs(B_nodiag[,i]),na.rm=T)
                                tab_generality[id_lieu,option_model[m],sp[i],"summed_raw"]=sum(B_nodiag[,i],na.rm=T)
			}
		}
	}
}

#Plotting now
for(m in 1:length(option_model)){
	pdf(paste("Rapport/graphe/MAR_estimates/",option_model[m],"metrics.pdf",sep=""))

	for(i in 1:2){
	#Set-up
	par(mfrow=c(2,2),mar=c(3,3,1,.5),oma=c(0.5,0.5,2,0.5))
	if(i==1){
		tab_tmp=tab_vulnerability
		atitle="Vulnerability"
	}else{
		tab_tmp=tab_generality
		atitle="Generality"
	}

	xlimi=c(min(c(tab_tmp[,option_model[m],,"intra"]),na.rm=T),max(c(tab_tmp[,option_model[m],,"intra"]),na.rm=T))
	ylimi=c(min(c(tab_tmp[,option_model[m],,"mean_raw"]),na.rm=T),max(c(tab_tmp[,option_model[m],,"mean_raw"]),na.rm=T))
	plot(0,0,t="n",xlab="",ylab="",xlim=xlimi,ylim=ylimi)
	mtext("Mean",side=2,font=2,line=2)
	mtext("Raw",side=3,font=2)
	for(l in 1:dim(tab_tmp)[1]){
		points(tab_tmp[l,option_model[m],,"intra"],tab_tmp[l,option_model[m],,"mean_raw"],xlab="",ylab="",pch=16,col=colo[l])
	}	
	x1=c(tab_tmp[,option_model[m],,"intra"])
	y1=c(tab_tmp[,option_model[m],,"mean_raw"])
	lm1=lm(y1~x1)
	sumlm1=summary(lm1)
#	abline(a=lm1$coefficients[1],b=lm1$coefficients[2])
#	legend("topleft",paste("R2=",format(sumlm1$r.squared,digit=2),sep=""),bty="n")

        ylimi=c(min(c(tab_tmp[,option_model[m],,"mean_abs"]),na.rm=T),max(c(tab_tmp[,option_model[m],,"mean_abs"]),na.rm=T))
        plot(0,0,t="n",xlab="",ylab="",xlim=xlimi,ylim=ylimi)
	legend("topright",c("Brittany","Oléron","Arcachon","Mediterranean"),col=c("green","darkblue","cyan","darkred"),bty="n",pch=16)
	mtext("Abs",side=3,font=2)
	for(l in 1:dim(tab_tmp)[1]){
		points(tab_tmp[l,option_model[m],,"intra"],tab_tmp[l,option_model[m],,"mean_abs"],xlab="",ylab="",pch=16,col=colo[l])
	}
	x1=c(tab_tmp[,option_model[m],,"intra"])
	y1=c(tab_tmp[,option_model[m],,"mean_abs"])
	lm1=lm(y1~x1)
	sumlm1=summary(lm1)
#	abline(a=lm1$coefficients[1],b=lm1$coefficients[2])
#	legend("topleft",paste("R2=",format(sumlm1$r.squared,digit=2),sep=""),bty="n")
	
        ylimi=c(min(c(tab_tmp[,option_model[m],,"summed_raw"]),na.rm=T),max(c(tab_tmp[,option_model[m],,"summed_raw"]),na.rm=T))
        plot(0,0,t="n",xlab="",ylab="",xlim=xlimi,ylim=ylimi)
	mtext("Summed",side=2,font=2,line=2)
	mtext("intra",side=1,line=2)
	for(l in 1:dim(tab_tmp)[1]){
		points(tab_tmp[l,option_model[m],,"intra"],tab_tmp[l,option_model[m],,"summed_raw"],xlab="",ylab="",pch=16,col=colo[l])
	}	
	x1=c(tab_tmp[,option_model[m],,"intra"])
	y1=c(tab_tmp[,option_model[m],,"summed_raw"])
	lm1=lm(y1~x1)
	sumlm1=summary(lm1)
#	abline(a=lm1$coefficients[1],b=lm1$coefficients[2])
#	legend("topleft",paste("R2=",format(sumlm1$r.squared,digit=2),sep=""),bty="n")

        ylimi=c(min(c(tab_tmp[,option_model[m],,"summed_abs"]),na.rm=T),max(c(tab_tmp[,option_model[m],,"summed_abs"]),na.rm=T))
        plot(0,0,t="n",xlab="",ylab="",xlim=xlimi,ylim=ylimi)
	title(atitle,outer=T)
	mtext("intra",side=1,line=2)
	for(l in 1:dim(tab_tmp)[1]){
		points(tab_tmp[l,option_model[m],,"intra"],tab_tmp[l,option_model[m],,"summed_abs"],xlab="",ylab="",pch=16,col=colo[l])
	}
	x1=c(tab_tmp[,option_model[m],,"intra"])
	y1=c(tab_tmp[,option_model[m],,"summed_abs"])
	lm1=lm(y1~x1)
	sumlm1=summary(lm1)
#	abline(a=lm1$coefficients[1],b=lm1$coefficients[2])
#	legend("topleft",paste("R2=",format(sumlm1$r.squared,digit=2),sep=""),bty="n")
	}
	dev.off()
}

#Second plot
for(m in 1:length(option_model)){
	pdf(paste("Rapport/graphe/MAR_estimates/",option_model[m],"_intra_for_each_species.pdf",sep=""))
	plot(0,0,xlim=c(0,dim(tab_vulnerability)[3]),ylim=c(-0.65,-0.15),t="n",xlab="",ylab="Intra",xaxt="n")
	axis(1,labels=dimnames(tab_vulnerability)[[3]],at=1:(dim(tab_vulnerability)[3]),las=2)
	for(s in 1:dim(tab_vulnerability)[3]){
		for(l in 1:dim(tab_vulnerability)[1]){
			points(s,tab_vulnerability[l,option_model[m],s,"intra"],pch=16,col=colo[l])
		}
	}
	dev.off()
}
