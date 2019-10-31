#2018/09/14 Just following an idea from Barabas&Allesina 2015, looking at the spectrum of each matrix' eigenvalues

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

if(take0){
        extant="with0"
}else{
        extant="no0"
}

pdf(paste("report/graphe/MAR_estimates/eigen_spectra_",option_model,"_",extant,"_diag_unchanged.pdf",sep=""))
par(mfrow=c(2,2),mar=c(4,4,1.5,.25))
#colo=c(rgb(0, 56, 168, maxColorValue = 255),rgb(115, 79, 150, maxColorValue = 255),rgb(215, 2, 112, maxColorValue = 255))
colo=c("darkblue","lightblue","darkorchid")
for (g in groupe){
        if(g=="BZ"){
                option_lieu=c("Men er Roue","Loscolo","Croisic")
		axlab=""
		aylab="Im"
        }else if(g=="MO"){
                option_lieu=c("LEperon","Cornard","Auger")
		axlab=""
		aylab=""
        }else if(g=="SU"){
                option_lieu=c("Antoine","Lazaret")
		axlab="Re"
		aylab=""
        }else if(g=="AR"){
                option_lieu=c("Teychan","B7")
		axlab="Re"
		aylab="Im"
        }
	#plot(0,0,t="n",xlim=c(-0.65,-.15),ylim=c(-0.2,0.2),xlab=axlab,ylab=aylab,main=g)
	plot(0,0,t="n",xlim=c(-1,1),ylim=c(-1,1),xlab=axlab,ylab=aylab,main=g)
        for (ll in 1:length(option_lieu)){
                f1=paste("data/analyse_MAR/",g,"/site_specific/",option_lieu[ll],"_",option_model,"_",option_NEI,"_regular_common_",g,".RData",sep="")
                load(f1)
                B=clean_matrix(cis,diag_bool=F)
		eig_B=eigen(B)$values
		id=which(Mod(eig_B)==max(Mod(eig_B)))
		points(Re(eig_B),Im(eig_B),col=colo[ll],pch=16,cex=2)
		points(Re(eig_B[id]),Im(eig_B[id]),col="red",cex=2.2)
	}
}
dev.off()
