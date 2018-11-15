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
only_signif=TRUE

methods=c("intra_mean","intra_var","inter_mean_abs","inter_mean_raw","inter_var_raw","inter_var_abs","stability","linkage density","weighted connectance","nb_pos","nb_neg","nb_negpos")

name_species=c("AST","CHA","DIT","GUI","LEP","NIT","PLE","PSE","RHI","SKE","THP","THL","GYM","PRO","PRP","SCR","CRY","EUG")

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
tab_value[,"nb_pos"]=rep(0,10)
tab_value[,"nb_neg"]=rep(0,10)
tab_value[,"nb_negpos"]=rep(0,10)

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
			tab_value[id_lieu,"stability"]=max(Mod(eigen(B)$values))
			metrics=networklevel(web=abs(B),index=c("linkage density","weighted connectance"),empty.web=F)
			tab_value[id_lieu,"linkage density"]=metrics["linkage density"]
			tab_value[id_lieu,"weighted connectance"]=metrics["weighted connectance"]
			if(!take0){
				B_nodiag[B_nodiag==0]=NA
			}
			tab_value[id_lieu,"inter_mean_abs"]=mean(c(abs(B_nodiag)),na.rm=T)
			tab_value[id_lieu,"inter_mean_raw"]=mean(c(B_nodiag),na.rm=T)
			tab_value[id_lieu,"inter_var_raw"]=var(c(B_nodiag),na.rm=T)
			tab_value[id_lieu,"inter_var_abs"]=var(c(abs(B_nodiag)),na.rm=T)

			for(i in 1:dim(B)[1]){	
				for(j in 1:dim(B)[1]){
					if(B[i,j]>0&B[j,i]>0){tab_value[id_lieu,"nb_pos"]=tab_value[id_lieu,"nb_pos"]+1}
					if(B[i,j]<0&B[j,i]<0){tab_value[id_lieu,"nb_neg"]=tab_value[id_lieu,"nb_neg"]+1}
					if((B[i,j]<0&B[j,i]>0)||(B[i,j]>0&B[j,i]<0)){tab_value[id_lieu,"nb_negpos"]=tab_value[id_lieu,"nb_negpos"]+1}
				}
			}
			tab_value[id_lieu,"nb_pos"]=tab_value[id_lieu,"nb_pos"]/sum(B!=0)
			tab_value[id_lieu,"nb_neg"]=tab_value[id_lieu,"nb_neg"]/sum(B!=0)
			tab_value[id_lieu,"nb_negpos"]=tab_value[id_lieu,"nb_negpos"]/sum(B!=0)
		}
}

pdf(paste("Rapport/graphe/MAR_estimates/moy_var_per_subsite_",option_model,"_",extant,"_",extant2,".pdf",sep=""))
par(mfrow=c(3,1),mar=c(2,4,0.5,0.5))
ylim1=min(tab_value[,"intra_mean"])-sqrt(max(tab_value[,"intra_var"]))/2
ylim2=max(tab_value[,"intra_mean"])+sqrt(max(tab_value[,"intra_var"]))/2
plot(0,0,t="n",xlab="",xaxt="n",ylab="Intra",xlim=c(1,10),ylim=c(ylim1,ylim2))
for(i in 1:dim(tab_value)[1]){
	points(i,tab_value[i,"intra_mean"],pch=16,col=colo[i])
	arrows(i,tab_value[i,"intra_mean"]-sqrt(tab_value[i,"intra_var"])/2,i,tab_value[i,"intra_mean"]+sqrt(tab_value[i,"intra_var"])/2,col=colo[i],length=0)
}
ylim1=min(c(tab_value[,"inter_mean_raw"]-sqrt(tab_value[,"inter_var_raw"])/2))
ylim2=max(c(tab_value[,"inter_mean_raw"]+sqrt(tab_value[,"inter_var_raw"])/2))
#ylim2=max(c(tab_value[,"inter_mean_raw"],tab_value[,"inter_mean_abs"]))+sqrt(max(tab_value[,"inter_var_raw"]))/2
plot(0,0,t="n",xlab="",xaxt="n",ylab="raw inter",xlim=c(1,10),ylim=c(ylim1,ylim2))
for(i in 1:dim(tab_value)[1]){
	points(i,tab_value[i,"inter_mean_raw"],pch=16,col=colo[i])
	arrows(i,tab_value[i,"inter_mean_raw"]-sqrt(tab_value[i,"inter_var_raw"])/2,i,tab_value[i,"inter_mean_raw"]+sqrt(tab_value[i,"inter_var_raw"])/2,col=colo[i],length=0)
}
ylim1=min(c(tab_value[,"inter_mean_abs"])-sqrt(tab_value[,"inter_var_abs"])/2)
ylim2=max(c(tab_value[,"inter_mean_abs"])+sqrt(tab_value[,"inter_var_abs"])/2)
plot(0,0,t="n",xlab="",xaxt="n",ylab="abs inter",xlim=c(1,10),ylim=c(ylim1,ylim2))
axis(1,at=1:10,labels=rownames(tab_value))
for(i in 1:dim(tab_value)[1]){
	points(i,tab_value[i,"inter_mean_abs"],pch=16,col=colo[i])
	arrows(i,tab_value[i,"inter_mean_abs"]-sqrt(tab_value[i,"inter_var_abs"])/2,i,tab_value[i,"inter_mean_abs"]+sqrt(tab_value[i,"inter_var_abs"])/2,col=colo[i],length=0)
}
dev.off()

pdf(paste("Rapport/graphe/MAR_estimates/moy_var_ratio_per_subsite_",option_model,"_",extant,"_",extant2,".pdf",sep=""),width=14)
par(mfrow=c(2,2),mar=c(4,4,3,0.025))
ylim1=min(tab_value[,"intra_mean"]/tab_value[,"inter_mean_raw"])
ylim2=max(tab_value[,"intra_mean"]/tab_value[,"inter_mean_raw"])
#ylim1=-20
#ylim2=25
plot(0,0,t="n",xlab="",xaxt="n",ylab="Mean(Intra)/Mean(inter raw)",xlim=c(1,10),ylim=c(ylim1,ylim2))
for(i in 1:dim(tab_value)[1]){
        points(i,tab_value[i,"intra_mean"]/tab_value[i,"inter_mean_raw"],pch=16,col=colo[i])
	if(tab_value[i,"intra_mean"]/tab_value[i,"inter_mean_raw"]>ylim2){
        	points(i,ylim2,pch=18,col=colo[i])
        	text(i,ylim2,format(tab_value[i,"intra_mean"]/tab_value[i,"inter_mean_raw"],digits=1),pch=18,col=colo[i],adj=0,xpd=NA)
	}
}
ylim1=min(tab_value[,"intra_mean"]/tab_value[,"inter_mean_abs"])
ylim2=max(tab_value[,"intra_mean"]/tab_value[,"inter_mean_abs"])
plot(0,0,t="n",xlab="",xaxt="n",ylab="Mean(Intra)/Mean(inter abs)",xlim=c(1,10),ylim=c(ylim1,ylim2))
for(i in 1:dim(tab_value)[1]){
        points(i,tab_value[i,"intra_mean"]/tab_value[i,"inter_mean_abs"],pch=16,col=colo[i])
}
ylim1=min(tab_value[,"intra_var"]/tab_value[,"inter_var_raw"])
ylim2=max(tab_value[,"intra_var"]/tab_value[,"inter_var_raw"])
plot(0,0,t="n",xlab="",xaxt="n",ylab="Var(Intra)/Var(inter raw)",xlim=c(1,10),ylim=c(ylim1,ylim2))
axis(1,at=1:10,labels=rownames(tab_value))
for(i in 1:dim(tab_value)[1]){
        points(i,tab_value[i,"intra_var"]/tab_value[i,"inter_var_raw"],pch=16,col=colo[i])
	print(tab_value[i,"intra_var"]/tab_value[i,"inter_var_raw"])
}
ylim1=min(tab_value[,"intra_var"]/tab_value[,"inter_var_abs"])
ylim2=max(tab_value[,"intra_var"]/tab_value[,"inter_var_abs"])
plot(0,0,t="n",xlab="",xaxt="n",ylab="Var(Intra)/Var(inter abs)",xlim=c(1,10),ylim=c(ylim1,ylim2))
axis(1,at=1:10,labels=rownames(tab_value))
for(i in 1:dim(tab_value)[1]){
        points(i,tab_value[i,"intra_var"]/tab_value[i,"inter_var_abs"],pch=16,col=colo[i])
	print(tab_value[i,"intra_var"]/tab_value[i,"inter_var_abs"])
}

dev.off()

pdf(paste("Rapport/graphe/MAR_estimates/moy_vs_var_",option_model,"_",extant,"_",extant2,".pdf",sep=""),width=12)
par(mfrow=c(1,3),mar=c(5,5,4,3))
xlimi=c(min(tab_value[,"intra_mean"]),max(tab_value[,"intra_mean"]))
ylimi=c(min(tab_value[,"intra_var"]),max(tab_value[,"intra_var"]))
plot(0,0,t="n",xlab="Mean",ylab="Variance",xlim=xlimi,ylim=ylimi,main="Intra",cex=3.0,cex.lab=3.0,cex.axis=2.25,cex.main=3.0)
for(i in 1:dim(tab_value)[1]){
	points(tab_value[i,"intra_mean"],tab_value[i,"intra_var"],pch=16,col=colo[i],cex=3.0)
        text(tab_value[i,"intra_mean"],tab_value[i,"intra_var"],rownames(tab_value)[i],col=colo[i],adj=0,xpd=NA,cex=2.5)
}
legend("topright",paste(format(mean(tab_value[,"intra_mean"]),digits=2),"(",format(mean(tab_value[,"intra_var"]),digits=2),")",sep=""),bty="n",cex=2.5)
xlimi=c(min(tab_value[,"inter_mean_raw"]),max(tab_value[,"inter_mean_raw"]))
ylimi=c(min(tab_value[,"inter_var_raw"]),max(tab_value[,"inter_var_raw"]))
plot(0,0,t="n",xlab="Mean",ylab="",xlim=xlimi,ylim=ylimi,main="Inter raw",cex=3.0,cex.lab=3.0,cex.axis=2.25,cex.main=3.0)
for(i in 1:dim(tab_value)[1]){
        points(tab_value[i,"inter_mean_raw"],tab_value[i,"inter_var_raw"],pch=16,col=colo[i],cex=3.0)
        text(tab_value[i,"inter_mean_raw"],tab_value[i,"inter_var_raw"],rownames(tab_value)[i],col=colo[i],adj=0,xpd=NA,cex=2.5)
}
legend("topleft",paste(format(mean(tab_value[,"inter_mean_raw"]),digits=2),"(",format(mean(tab_value[,"inter_var_raw"]),digits=2),")",sep=""),bty="n",cex=2.5)
xlimi=c(min(tab_value[,"inter_mean_abs"]),max(tab_value[,"inter_mean_abs"]))
ylimi=c(min(tab_value[,"inter_var_abs"]),max(tab_value[,"inter_var_abs"]))
plot(0,0,t="n",xlab="Mean",ylab="",xlim=xlimi,ylim=ylimi,main="Inter abs",cex=3.0,cex.lab=3.0,cex.axis=2.25,cex.main=3.0)
for(i in 1:dim(tab_value)[1]){
        points(tab_value[i,"inter_mean_abs"],tab_value[i,"inter_var_abs"],pch=16,col=colo[i],cex=3.0)
        text(tab_value[i,"inter_mean_abs"],tab_value[i,"inter_var_abs"],rownames(tab_value)[i],col=colo[i],adj=0,xpd=NA,cex=2.5)
}
legend("topleft",paste(format(mean(tab_value[,"inter_mean_abs"]),digits=2),"(",format(mean(tab_value[,"inter_var_abs"]),digits=2),")",sep=""),bty="n",cex=2.5)
dev.off()

pdf(paste("Rapport/graphe/MAR_estimates/stability_vs_var_",option_model,"_",extant,"_",extant2,".pdf",sep=""),width=12)
par(mfrow=c(1,3),mar=c(5,5,4,.5))
ylimi=c(min(tab_value[,"stability"]),max(tab_value[,"stability"]))
xlimi=c(min(tab_value[,"intra_var"]),max(tab_value[,"intra_var"]))
plot(0,0,t="n",xlab="Var",ylab="Stability",xlim=xlimi,ylim=ylimi,main="Intra",cex=3.0,cex.lab=3.0,cex.axis=2.25,cex.main=3.0)
for(i in 1:dim(tab_value)[1]){
        points(tab_value[i,"intra_var"],tab_value[i,"stability"],pch=16,col=colo[i],cex=3)
#        text(tab_value[i,"intra_var"],tab_value[i,"stability"],rownames(tab_value)[i],col=colo[i],adj=0,xpd=NA,cex=2.5)
}
xlimi=c(min(tab_value[,"inter_var_raw"]),max(tab_value[,"inter_var_raw"]))
plot(0,0,t="n",xlab="Var",ylab="",xlim=xlimi,ylim=ylimi,main="Inter raw",cex=3.0,cex.lab=3.0,cex.axis=2.25,cex.main=3.0)
for(i in 1:dim(tab_value)[1]){
        points(tab_value[i,"inter_var_raw"],tab_value[i,"stability"],pch=16,col=colo[i],cex=3.0)
#        text(tab_value[i,"inter_var_raw"],tab_value[i,"stability"],rownames(tab_value)[i],col=colo[i],adj=0,xpd=NA,cex=2.5)
}
xlimi=c(min(tab_value[,"inter_var_abs"]),max(tab_value[,"inter_var_abs"]))
plot(0,0,t="n",xlab="Var",ylab="",xlim=xlimi,ylim=ylimi,main="Inter abs",cex=3.0,cex.lab=3.0,cex.axis=2.25,cex.main=3.0)
for(i in 1:dim(tab_value)[1]){
        points(tab_value[i,"inter_var_abs"],tab_value[i,"stability"],pch=16,col=colo[i],cex=3.0)
#        text(tab_value[i,"inter_var_abs"],tab_value[i,"stability"],rownames(tab_value)[i],col=colo[i],adj=0,xpd=NA,cex=2.5)
}
legend(x=0.0019,y=0.475,c("Brittany","Oléron","Arcachon","Mediterr."),col=c("green","darkblue","cyan","darkred"),bty="n",pch=16,cex=2.5,xpd=NA)
dev.off()

pdf(paste("Rapport/graphe/MAR_estimates/metrics_vs_stability_",option_model,"_",extant,"_",extant2,".pdf",sep=""),width=12)
par(mfrow=c(1,2),mar=c(4,4.5,1,.5))
plot(tab_value[,'weighted connectance'],tab_value[,"stability"],col=colo,xlab="weighted connectance",ylab="Max_eig_B",pch=16,cex=2.5,cex.axis=2,cex.lab=2)
legend("topleft",c("Brittany","Oléron","Arcachon","Mediterranean"),col=c("green","darkblue","cyan","darkred"),bty="n",pch=16,cex=2.0)
plot(tab_value[,'linkage density'],tab_value[,"stability"],col=colo,xlab="linkage density",ylab="",pch=16,cex=2.5,cex.axis=2,cex.lab=2)
dev.off()

pdf(paste("Rapport/graphe/MAR_estimates/moy_vs_var_for_all_interactions_",option_model,"_",extant,"_",extant2,".pdf",sep=""),width=12)
par(mfrow=c(1,2),mar=c(4,4.5,.5,.5))
x1=c(tab_value[,"intra_mean"],tab_value[,'inter_mean_raw'])
x2=c(tab_value[,"intra_mean"],tab_value[,'inter_mean_abs'])
y2=c(tab_value[,"intra_var"],tab_value[,'inter_var_abs'])
y1=c(tab_value[,"intra_var"],tab_value[,'inter_var_raw'])
plot(x1,y1,col=colo,pch=16,xlab="Mean",ylab="Variance",cex=2,cex.axis=2,cex.lab=2)
text(mean(tab_value[,'intra_mean']),0.01,"Intra",xpd=NA,cex=1.5)
legend("topright",c("Brittany","Oléron","Arcachon","Mediterranean"),col=c("green","darkblue","cyan","darkred"),bty="n",pch=16,cex=2.0)
text(mean(tab_value[,'inter_mean_raw']),min(tab_value[,"intra_var"]),"Inter raw",xpd=NA,cex=1.5)
plot(x2,y2,col=colo,pch=16,xlab="Mean",ylab="",cex=2,cex.axis=2,cex.lab=2)
text(mean(tab_value[,'intra_mean']),0.01,"Intra",xpd=NA,cex=1.5)
text(mean(tab_value[,'inter_mean_abs']),min(tab_value[,"intra_var"]),"Inter abs",xpd=NA,cex=1.5)
dev.off()
