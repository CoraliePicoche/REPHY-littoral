#2018/09/17 CP
#Compare the strength of interactions with the effect of covariates

rm(list=ls())
graphics.off()
library(MARSS)
source("script/matrix_MAR_clean.r")

option_model="pencen"
option_NEI="null" #We don't take the "Not Elsewhere Identified" species into account
groupe=c("BZ","MO","AR","SU")
option_sp="common" #Species are the same 
take0=TRUE

methods=c("intra_mean","intra_var","inter_mean_abs","inter_mean_raw","inter_var_raw","inter_var_abs","evt_abs","evt_raw")

id_lieu=0

if(take0){
        extant="with0"
}else{
        extant="no0"
}

tab_value=matrix(NA,nrow=10,ncol=length(methods))
colnames(tab_value)=methods
rownames(tab_value)=rep("",10)

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
                rownames(tab_value)[id_lieu]=option_lieu[ll]
                        f1=paste("data/analyse_MAR/",g,"/site_specific/",option_lieu[ll],"_",option_model,"_",option_NEI,"_regular_common_",g,".RData",sep="")
                        load(f1)

                        B=clean_matrix(cis)
                        B_nodiag=B
                        diag(B_nodiag)=NA

                        tab_value[id_lieu,"intra_mean"]=mean(diag(B))
                        tab_value[id_lieu,"intra_var"]=var(diag(B))
                        if(!take0){
                                B_nodiag[B_nodiag==0]=NA
                        }
                        tab_value[id_lieu,"inter_mean_abs"]=mean(c(abs(B_nodiag)),na.rm=T)
                        tab_value[id_lieu,"inter_mean_raw"]=mean(c(B_nodiag),na.rm=T)
                        tab_value[id_lieu,"inter_var_raw"]=var(c(B_nodiag),na.rm=T)
                        tab_value[id_lieu,"inter_var_abs"]=var(c(abs(B_nodiag)),na.rm=T)
			tab_value[id_lieu,"evt_abs"]=mean(abs(cis$par$U))
			tab_value[id_lieu,"evt_raw"]=mean(cis$par$U)
		}
	}
pdf(paste("Rapport/graphe/MAR_estimates/inter_vs_evt_",option_model,"_",extant,".pdf",sep=""),width=11)
par(mfrow=c(1,2),mar=c(4,4,1,1))
plot(tab_value[,"intra_mean"],tab_value[,"evt_abs"],xlab="Intra",ylab="Evt",pch=16,col=colo,cex=2)
plot(tab_value[,"inter_mean_abs"],tab_value[,"evt_abs"],xlab="Inter",ylab="",pch=16,col=colo,cex=2)
dev.off()
