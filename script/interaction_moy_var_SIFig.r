#2018/09/10 CP
#

rm(list=ls())
graphics.off()
library(MARSS)
library(bipartite)
source("script/matrix_MAR_clean.r")

option_model="pencen"
option_NEI="null" #We don't take the "Not Elsewhere Identified" species into account
groupe=c("BZ","MO","AR","SU")
option_sp="common" #Species are the same 
take0=FALSE
only_signif=FALSE

methods=c("intra_mean","intra_var","inter_mean_abs","inter_mean_raw","inter_var_raw","inter_var_abs")

id_lieu=0

if(take0){
	extant="with0"
}else{
	extant="no0"
}
if(only_signif){
	extant2="coeffsignif"
}else{
	extant2="allcoeff"
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

                        B=clean_matrix(cis,signif=only_signif)
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
		}
}

pdf(paste("article/graphe/moy_vs_var_for_all_interactions_",option_model,"_SI.pdf",sep=""),width=14)
par(mfrow=c(1,2),mar=c(4,4.5,3,.25))
x1=c(tab_value[,"intra_mean"],tab_value[,'inter_mean_raw'])
x2=c(tab_value[,"intra_mean"],tab_value[,'inter_mean_abs'])
y2=c(tab_value[,"intra_var"],tab_value[,'inter_var_abs'])
y1=c(tab_value[,"intra_var"],tab_value[,'inter_var_raw'])
plot(x1,y1,col=colo,pch=16,xlab="Mean coefficient",ylab="Variance of coefficients",cex=2.5,cex.axis=2,cex.lab=2)
mtext("a)",side=3,cex=1.75,xpd=NA,font=2,line=.8,adj=0)
text(mean(tab_value[,'intra_mean']),0.01,"Intra",xpd=NA,cex=1.5)
legend("topright",c("Brittany","Oléron","Arcachon","Mediterranean"),col=c("green","darkblue","cyan","darkred"),bty="n",pch=16,cex=2.0)
text(mean(tab_value[,'inter_mean_raw'])-0.025,0.006,"Inter raw",xpd=NA,cex=1.5)
plot(x2,y2,col=colo,pch=16,xlab="Mean coefficient",ylab="",cex=2.5,cex.axis=2,cex.lab=2,yaxt="n")
mtext("b)",side=3,cex=1.75,xpd=NA,font=2,line=.8,adj=0)
text(mean(tab_value[,'intra_mean']),0.01,"Intra",xpd=NA,cex=1.5)
text(mean(tab_value[,'inter_mean_abs'])-0.025,min(tab_value[,"intra_var"]),"Inter abs",xpd=NA,cex=1.5)
dev.off()
