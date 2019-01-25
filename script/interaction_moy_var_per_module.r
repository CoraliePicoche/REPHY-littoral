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

methods=c("intra_mean","intra_var","intra_cv","inter_mean_abs","inter_mean_raw","inter_var_raw","inter_var_abs","inter_cv_raw","inter_cv_abs")

name_species=c("AST","CHA","DIT","GUI","LEP","NIT","PLE","PSE","RHI","SKE","THP","THL","GYM","PRO","PRP","SCR","CRY","EUG")
what_sp=c("P","C","C","C","C","P","P","P","C","C","C","P","D","D","D","D","O","O")


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

tmp_v=c("C","P","D","O")
tab_value=array(NA,dim=c(10,length(methods),4),dimnames=list(rep("",10),methods,c("C","P","D","O")))

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
                dimnames(tab_value)[[1]][id_lieu]=option_lieu[ll]
                        f1=paste("data/analyse_MAR/",g,"/site_specific/",option_lieu[ll],"_",option_model,"_",option_NEI,"_regular_common_",g,".RData",sep="")
                        load(f1)
			sp=dimnames(cis$model$data)[[1]]
			print(sp)
                        B=clean_matrix(cis,signif=only_signif)

			for(gr in tmp_v){
				id=intersect(sp,name_species[what_sp==gr])
				print(gr)
				print(id)
				if(length(id)>0){
					print(sp %in% id)	
					sub_g=B[sp %in% id,sp %in% id]	
					print(sub_g)
                        		sub_g_no_diag=sub_g
                
			        	diag(sub_g_no_diag)=NA
			
			tab_value[id_lieu,"intra_mean",gr]=mean(diag(sub_g))
			tab_value[id_lieu,"intra_var",gr]=var(diag(sub_g))
			tab_value[id_lieu,"intra_cv",gr]=mean(diag(sub_g),na.rm=T)/sd(diag(sub_g),na.rm=T)
			
			tab_value[id_lieu,"inter_mean_abs",gr]=mean(c(abs(sub_g_no_diag)),na.rm=T)
			tab_value[id_lieu,"inter_mean_raw",gr]=mean(c(sub_g_no_diag),na.rm=T)
			tab_value[id_lieu,"inter_var_raw",gr]=var(c(sub_g_no_diag),na.rm=T)
			tab_value[id_lieu,"inter_var_abs",gr]=var(c(abs(sub_g_no_diag)),na.rm=T)
			tab_value[id_lieu,"inter_cv_abs",gr]=mean(c(abs(sub_g_no_diag)),na.rm=T)/sd(c(abs(sub_g_no_diag)),na.rm=T)
			tab_value[id_lieu,"inter_cv_raw",gr]=mean(c(sub_g_no_diag),na.rm=T)/sd(c(sub_g_no_diag),na.rm=T)
		}
		}
		}
}
pdf(paste("article/graphe/moy_vs_var_for_all_interactions_",option_model,"_per_cluster.pdf",sep=""),width=14)
par(mfrow=c(2,2),mar=c(4,4.5,3,.25))
xmini=min(c(tab_value[,"intra_mean",]),na.rm=T)
xmaxi=max(c(tab_value[,"inter_mean_abs",]),na.rm=T)
ymini=min(c(tab_value[,"intra_var",]),na.rm=T)
ymaxi=max(c(tab_value[,"intra_var",]),na.rm=T)
for(gr in tmp_v){
x2_intra=c(tab_value[,"intra_mean",gr])
x2_inter=c(tab_value[,'inter_mean_abs',gr])
y2_intra=c(tab_value[,"intra_var",gr])
y2_inter=c(tab_value[,"inter_var_abs",gr])
plot(0,0,xlab="Mean",ylab="Var",cex=2.5,cex.axis=2,cex.lab=2,main=gr,t="n",xlim=c(xmini,xmaxi),ylim=c(ymini,ymaxi))
points(x2_intra,y2_intra,col=colo,pch=16,cex=2)
points(x2_inter,y2_inter,col=colo,pch=17,cex=2)
legend("topright",c("Brittany","Oléron","Arcachon","Mediterranean"),col=c("green","darkblue","cyan","darkred"),bty="n",pch=16,cex=2.0)
}
xmini=min(c(tab_value[,"inter_mean_raw",]),na.rm=T)
xmaxi=max(c(tab_value[,"inter_mean_abs",]),na.rm=T)
ymini=min(c(tab_value[,"inter_var_raw",]),na.rm=T)
ymaxi=max(c(tab_value[,"inter_var_raw",]),na.rm=T)
for(gr in tmp_v){ 
x2_intra=c(tab_value[,"inter_mean_raw",gr])
x2_inter=c(tab_value[,'inter_mean_abs',gr])
y2_intra=c(tab_value[,"inter_var_raw",gr])
y2_inter=c(tab_value[,"inter_var_abs",gr])
plot(0,0,xlab="Mean",ylab="Var",cex=2.5,cex.axis=2,cex.lab=2,main=gr,t="n",xlim=c(xmini,xmaxi),ylim=c(ymini,ymaxi))
points(x2_intra,y2_intra,col=colo,pch=16,cex=2)
points(x2_inter,y2_inter,col=colo,pch=17,cex=2)
legend("topright",c("Brittany","Oléron","Arcachon","Mediterranean"),col=c("green","darkblue","cyan","darkred"),bty="n",pch=16,cex=2.0)
}
dev.off()
