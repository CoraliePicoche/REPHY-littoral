graphics.off()
rm(list=ls())
library(MARSS)

#First, set the 0
which_model0="null" #for now, can be "null", "unconstrained" and "pencen"
which_NEI0="null" #NEI can be either a population or a covariate

option_model=c("null","unconstrained","pencen")
option_NEI=c("cov","sp")
option_NEI=c("null")
#option_lieu=c("LEperon","Cornard","Auger")
#option_lieu=c("Men er Roue","Loscolo","Croisic")
#option_lieu=c("Antoine","Lazaret")
option_lieu=c("LEperon","Cornard","Auger","Men er Roue","Loscolo","Croisic","Antoine","Lazaret")
option_sp='common'
criterium="AICc"

colo=rainbow(length(option_model)*length(option_NEI)-1)
colo=palette()
yli1=-250
yli2=30
pdf(paste("Rapport/graphe/comp_",criterium,"_all_sites_common.pdf",sep=""),width=16)
par(mar=c(3,6,1,1))
plot(0,0,xlim=c(0,length(option_lieu)*10),ylim=c(yli1,yli2),t="n",xlab="",ylab=paste(expression(Delta),criterium,"=",criterium,"-",criterium,"(null model)",sep=""),xaxt="n",lwd=2.5,cex.lab=1.75,cex.axis=1.75,cex=2)
axis(1,at=seq(2,length(option_lieu)*10,10),lab=option_lieu,cex=2.5,cex.axis=1.75)
abline(h=0,lty=1,lwd=2.5)
if(option_sp=="common"){
	abline(v=25,lty=1,lwd=2.5)
	abline(v=55,lty=1,lwd=2.5)
}

#Loop on site
for (l in 1:length(option_lieu)){
                                if(option_sp=="common"){
                                        if(option_lieu[l] %in% c("LEperon","Cornard","Auger")){
                                        SS="MO"
                                        }else if(option_lieu[l] %in% c("Men er Roue","Loscolo","Croisic")){
                                        SS="BZ"
                                        }else{
                                        SS="SU"
                                        }
                                }else{
                                        SS="essai"
                                }
#for (l in 1:1){
	f0=paste("data/analyse_MAR/",option_lieu[l],"_",which_model0,"_",which_NEI0,"_regular_",option_sp,"_",SS,".RData",sep="")
	load(f0)
	AICc0=fit_log$AICc
	BIC0=-2*fit_log$logLik+log(fit_log$samp.size/dim(fit_log$states)[1])*fit_log$num.params

	mn=0
	leg=c()
	for(n in 1:length(option_NEI)){
		for (m in 1:length(option_model)){
			if(!(option_model[m]==which_model0&&option_NEI[n]==which_NEI0)){
				mn=mn+1
				leg=c(leg,paste(option_model[m],option_NEI[n],sep="_"))
				if(option_sp=="common"){
					if(option_lieu[l] %in% c("LEperon","Cornard","Auger")){
					SS="MO"
					}else if(option_lieu[l] %in% c("Men er Roue","Loscolo","Croisic")){
					SS="BZ"
					}else{
					SS="SU"
					}
				}else{
					SS="essai"
				}
					
				f1=paste("data/analyse_MAR/",option_lieu[l],"_",option_model[m],"_",option_NEI[n],"_regular_",option_sp,"_",SS,".RData",sep="")
				load(f1)
				diff_AICc=fit_log$AICc-AICc0
				diff_BIC=-2*fit_log$logLik+log(fit_log$samp.size/dim(fit_log$states)[1])*fit_log$num.params-BIC0
				
				if(criterium=="AICc"){
					diff=diff_AICc
				}else if(criterium=="BIC"){
					diff=diff_BIC
				}
				rect((l-1)*10+mn,0,(l-1)*10+mn+1,diff,col=colo[mn])
			}
		}
	}
}
legend("topleft",leg,fill=colo,bty="n")
dev.off()
