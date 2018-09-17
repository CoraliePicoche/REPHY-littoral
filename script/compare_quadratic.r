#18/09/14 CP
#Compare BIC and AIC from models with and without a polynomial form of the covariate (niche) + checks the shape of the curves

rm(list=ls())
graphics.off()
library(MARSS)
library(bipartite)
source("script/matrix_MAR_clean.r")

option_model=c("null","pencen","unconstrained")
option_NEI="null" #We don't take the "Not Elsewhere Identified" species into account
groupe=c("BZ","MO","AR","SU")
option_sp="common" #Species are the same 
criterium=c("BIC","AIC")

colo=c("darkblue","lightblue","darkorchid")
pdf("Rapport/graphe/MAR_estimates/compare_aic_bic_polynom.pdf")
par(mfrow=c(2,1),mar=c(5,5,.5,.5))
for (c in 1:length(criterium)){
if(c==1){
yli1=-150
yli2=150
}else{
yli1=-300
yli2=100
}
plot(0,0,xlim=c(0,length(option_model)*13),ylim=c(yli1,yli2),t="n",xlab="",ylab=paste(expression(Delta),criterium[c],"=",criterium[c],"(V²)-",criterium[c],sep=""),xaxt="n",lwd=2.5,cex.lab=1.5,cex.axis=1.5)
abline(h=0,lty=1,lwd=2.5)
id_lieu=0
id=-1
opp=c()
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
	opp=c(opp,option_lieu)
        for (ll in 1:length(option_lieu)){
		id_lieu=id_lieu+1
		min=10000000
		id_min=0
	for(m in 1:length(option_model)){
		if(m==1){
			id=id+2
		}else{
			id=id+1
		}
                f1=paste("data/analyse_MAR/",g,"/site_specific/",option_lieu[ll],"_",option_model[m],"_",option_NEI,"_regular_common_",g,".RData",sep="")
                load(f1)
		aic_1=fit_log$AICc
		bic_1=-2*fit_log$logLik+log(fit_log$samp.size/dim(fit_log$states)[1])*fit_log$num.params
                f2=paste("data/analyse_MAR/",g,"/site_specific/",option_lieu[ll],"_",option_model[m],"_",option_NEI,"_regular_common_ALL_squared.RData",sep="")
                load(f2)
		aic_2=fit_log$AICc
		bic_2=-2*fit_log$logLik+log(fit_log$samp.size/dim(fit_log$states)[1])*fit_log$num.params
		if(criterium[c]=="BIC"){
	        	rect(id,0,id+1,bic_2-bic_1,col=colo[m])
			if(bic_2<min){
				min=bic_2
				id_min=id+0.5
			}
		}else if(criterium[c]=="AIC"){
	        	rect(id,0,id+1,aic_2-aic_1,col=colo[m])
			if(aic_2<min){
				min=aic_2
				id_min=id+0.5
			}
		}
	}
		points(id_min,0,pch="*",cex=2,col="red")
	}
}
}
legend("bottomright",c("Null","Pencen","Unconstrained","Lowest crit."),fill=c(colo,NA),pch=c(NA,NA,NA,"*"),col=c(NA,NA,NA,"red"),bty="n",pt.cex=2.5,border=NA)
axis(1,at=seq(2,id,length.out=length(opp)),labels=opp,cex.axis=1.25,las=2)
dev.off()
