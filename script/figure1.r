#2018/09/04 CP: first figure of the paper, show all BIC for our 4 different sites

graphics.off()
rm(list=ls())
library(MARSS)

#First, set the 0
which_model0="null" #for now, can be "null", "unconstrained" and "pencen"
which_NEI0="null" #we do not take NEI into account

option_model=c("null","unconstrained","pencen","diatdin","inter")
option_NEI=c("null")
groupe=c("BZ","MO","AR","SU")
option_sp="common" #here, we will consider species that are specific to several subsites in one site
criterium="BIC"
yli1=0
yli2=800
option_lieu=c("Men er Roue","Loscolo","Croisic","LEperon","Cornard","Auger","Teychan","B7","Antoine","Lazaret")

colo=rainbow(length(option_model)*length(option_NEI)-1)
colo=palette()
pdf(paste("Rapport/graphe/MAR_estimates/COMMON/with_AR/comp_",criterium,"_per_site.pdf",sep=""),width=17,family="Helvetica")
par(mar=c(3,6,1,1))
plot(0,0,xlim=c(0,length(option_lieu)*10-5),ylim=c(yli1,yli2),t="n",xlab="",ylab=bquote(paste(Delta,.(criterium),"=",.(criterium),"-",.(criterium)[0],sep="")),xaxt="n",lwd=2.5,cex.lab=2.25,cex.axis=2.25,cex=2.5)
axis(1,at=seq(2,length(option_lieu)*10,10),lab=option_lieu,cex=2.25,cex.axis=1.65)
abline(h=0,lty=1,lwd=2.5)
       abline(v=27.5,lty=2,lwd=2.5)
       abline(v=57.5,lty=2,lwd=2.5)
       abline(v=77.5,lty=2,lwd=2.5)

id_lieu=0
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

#Loop on site
for (l in 1:length(option_lieu)){
        id_lieu=id_lieu+1
#for (l in 1:1){
        f0=paste("data/analyse_MAR/",g,"/site_specific/",option_lieu[l],"_",which_model0,"_",which_NEI0,"_regular_common_",g,".RData",sep="")
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

                                f1=paste("data/analyse_MAR/",g,"/site_specific/",option_lieu[l],"_",option_model[m],"_",option_NEI[n],"_regular_common_",g,".RData",sep="")
                                load(f1)
                                diff_AICc=fit_log$AICc-AICc0
                                diff_BIC=-2*fit_log$logLik+log(fit_log$samp.size/dim(fit_log$states)[1])*fit_log$num.params-BIC0

                                if(criterium=="AICc"){
                                        diff=diff_AICc
                                }else if(criterium=="BIC"){
                                        diff=diff_BIC
                                }
                                rect((id_lieu-1)*10+mn,0,(id_lieu-1)*10+mn+1,diff,col=colo[mn])
                        }
                }
        }
}
}
legend("topright",leg,fill=colo,bty="n",cex=1.75)
dev.off()
