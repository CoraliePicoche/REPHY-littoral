graphics.off()
rm(list=ls())
library(MARSS)

#First, set the 0
which_model0="null" #for now, can be "null", "unconstrained" and "pencen"
which_NEI0="null" #NEI can be either a population or a covariate

option_model=c("null","unconstrained","pencen","diatdin","inter")
option_NEI=c("null")
groupe=c("BZ","MO","SU") #and, later, "AR"
option_sp="common"
criterium="BIC"
yli1=-00
yli2=800
#option_lieu=c("Men er Roue","Loscolo","Croisic","LEperon","Cornard","Auger","Antoine","Lazaret")
for (g in groupe){
        if(g=="BZ"){
                option_lieu=c("Men er Roue","Loscolo","Croisic")
        }else if(g=="MO"){
                option_lieu=c("LEperon","Cornard","Auger")
        }else if(g=="SU"){
                option_lieu=c("Antoine","Lazaret")
        }


colo=rainbow(length(option_model)*length(option_NEI)-1)
colo=palette()
pdf(paste("Rapport/graphe/comp_",criterium,"_",g,"_common_with_bootstrap.pdf",sep=""),width=16)
par(mar=c(3,6,1,1))
plot(0,0,xlim=c(0,length(option_lieu)*10),ylim=c(yli1,yli2),t="n",xlab="",ylab=paste(expression(Delta),criterium,"=",criterium,"-",criterium,"(null model)",sep=""),xaxt="n",lwd=2.5,cex.lab=1.75,cex.axis=1.75,cex=2)
axis(1,at=seq(2,length(option_lieu)*10,10),lab=option_lieu,cex=2.5,cex.axis=1.75)
abline(h=0,lty=1,lwd=2.5)
#if(option_sp=="common"){
#	abline(v=25,lty=1,lwd=2.5)
#	abline(v=55,lty=1,lwd=2.5)
#}

#Loop on site
for (l in 1:length(option_lieu)){
#for (l in 1:1){
	f0=paste("data/analyse_MAR/",option_lieu[l],"_",which_model0,"_",which_NEI0,"_regular_",option_sp,"_",g,".RData",sep="")
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
					
				f1=paste("data/analyse_MAR/",option_lieu[l],"_",option_model[m],"_",option_NEI[n],"_regular_",option_sp,"_",g,".RData",sep="")
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
legend("topright",leg,fill=colo,bty="n")
dev.off()
}
