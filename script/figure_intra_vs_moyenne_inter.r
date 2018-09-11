#2018/09/10 CP
#This script plots the relation between intra- and inter-group interaction coefficients

rm(list=ls())
graphics.off()
library(MARSS)
source("script/matrix_MAR_clean.r")

option_model=c("unconstrained","pencen")
#option_model="pencen"
option_NEI=c("null") #We don't take the "Not Elsewhere Identified" species into account
groupe=c("BZ","MO","AR","SU")
option_sp="common" #Species are the same 

methods=c("intra","mean_raw","mean_abs")

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
                        B=clean_matrix(cis)
			B_nodiag=B
			diag(B_nodiag)=NA
                        for (i in 1:length(sp)){
                                tab_vulnerability[id_lieu,option_model[m],sp[i],"intra"]=B[i,i]
                                tab_generality[id_lieu,option_model[m],sp[i],"intra"]=B[i,i]
             			tab_vulnerability[id_lieu,option_model[m],sp[i],"mean_raw"]=mean(B_nodiag[i,B_nodiag[i,]!=0],na.rm=T)
	             		tab_vulnerability[id_lieu,option_model[m],sp[i],"mean_abs"]=mean(abs(B_nodiag[i,B_nodiag[i,]!=0]),na.rm=T)
        	                tab_generality[id_lieu,option_model[m],sp[i],"mean_raw"]=mean(B_nodiag[B_nodiag[,i]!=0,i],na.rm=T)
                	        tab_generality[id_lieu,option_model[m],sp[i],"mean_abs"]=mean(abs(B_nodiag[B_nodiag[,i]!=0,i]),na.rm=T)
                        }

		}
	}
}

#Plotting now
let_gg=c("a)","b)","c)","d)")
for(m in 1:length(option_model)){
	pdf(paste("Rapport/graphe/MAR_estimates/",option_model[m],"intra_vs_mean_inter_no0.pdf",sep=""))

	par(mfrow=c(2,2),mar=c(4,3,.5,.75),oma=c(0.5,0.5,2,0.5))
	id_gg=0
	for(i in 1:2){
	#Set-up
	if(i==1){
		tab_tmp=tab_vulnerability
		aylab=c(expression(tilde(b[i.])),expression(bgroup("|",tilde(b[i.]),"|")))
		axlab=""
	}else{
		tab_tmp=tab_generality
		aylab=c(expression(tilde(b[.i])),expression(bgroup("|",tilde(b[.i]),"|")))
		axlab=expression(b[ii])
	}

	xlimi=c(min(c(tab_tmp[,option_model[m],,"intra"]),na.rm=T),max(c(tab_tmp[,option_model[m],,"intra"]),na.rm=T))
	ylimi=c(min(c(tab_tmp[,option_model[m],,"mean_raw"]),na.rm=T),max(c(tab_tmp[,option_model[m],,"mean_raw"]),na.rm=T))
	id_gg=id_gg+1
	plot(0,0,t="n",ylab="",xlab=axlab,xlim=xlimi,ylim=ylimi)
        mtext(let_gg[id_gg],side=3,cex=1.,xpd=NA,font=2,line=0,adj=0)
        mtext(aylab[1],side=2,cex=1,line=1.75)
	for(l in 1:dim(tab_tmp)[1]){
		points(tab_tmp[l,option_model[m],,"intra"],tab_tmp[l,option_model[m],,"mean_raw"],pch=16,col=colo[l])
	}
	x1=c(tab_tmp[,option_model[m],,"intra"])
	y1=c(tab_tmp[,option_model[m],,"mean_raw"])
	lm1=lm(y1~x1)
	sumlm1=summary(lm1)
	abline(a=lm1$coefficients[1],b=lm1$coefficients[2])
	#legend("topleft",paste("R2=",format(sumlm1$r.squared,digit=2),sep=""),bty="n")
	legend("topleft",paste("cor=",format(cor(x1,y1,use="complete.obs"),digit=2),sep=""),bty="n")

        ylimi=c(min(c(tab_tmp[,option_model[m],,"mean_abs"]),na.rm=T),max(c(tab_tmp[,option_model[m],,"mean_abs"]),na.rm=T))
        plot(0,0,t="n",xlab=axlab,ylab="",xlim=xlimi,ylim=ylimi)
	legend("topright",c("Brittany","Oléron","Arcachon","Mediterranean"),col=c("green","darkblue","cyan","darkred"),bty="n",pch=16)
	id_gg=id_gg+1
        mtext(let_gg[id_gg],side=3,cex=1.,xpd=NA,font=2,line=0,adj=0)
        mtext(aylab[2],side=2,cex=1,line=1.75)
	for(l in 1:dim(tab_tmp)[1]){
		points(tab_tmp[l,option_model[m],,"intra"],tab_tmp[l,option_model[m],,"mean_abs"],pch=16,col=colo[l])
	}
	x1=c(tab_tmp[,option_model[m],,"intra"])
	y1=c(tab_tmp[,option_model[m],,"mean_abs"])
	lm1=lm(y1~x1)
	sumlm1=summary(lm1)
	abline(a=lm1$coefficients[1],b=lm1$coefficients[2])
	#legend("topleft",paste("R2=",format(sumlm1$r.squared,digit=2),sep=""),bty="n")
	legend("topleft",paste("cor=",format(cor(x1,y1,use="complete.obs"),digit=2),sep=""),bty="n")
	}	
	dev.off()
}
