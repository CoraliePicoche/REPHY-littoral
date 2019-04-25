#2019/01/03: mutualistic interactions? Commensalist interactions?

rm(list=ls())
graphics.off()
source("script/matrix_MAR_clean.r")
library(plotrix)

option_model="pencen"
option_NEI=c("null") #We don't take the "Not Elsewhere Identified" species into account
groupe=c("BZ","MO","AR","SU")
option_sp="common" #Species are the same 

pdf("article/graphe/pair_coefficients.pdf")
plot(0,0,xlim=c(-.23,.23),ylim=c(-.23,.23),t="n",xlab=expression("b"["ij"]),ylab=expression('b'['ji']))
abline(h=0)
abline(v=0)
for (g in groupe){
        if(g=="BZ"){
                option_lieu=c("Men er Roue","Loscolo","Croisic")
		colo="green"
        }else if(g=="MO"){
                option_lieu=c("LEperon","Cornard","Auger")
		colo="darkblue"
        }else if(g=="SU"){
                option_lieu=c("Antoine","Lazaret")
		colo="darkred"
        }else if(g=="AR"){
                option_lieu=c("Teychan","B7")
		colo="cyan"
        }
        for (ll in 1:length(option_lieu)){
                f1=paste("data/analyse_MAR/",g,"/site_specific/",option_lieu[ll],"_",option_model,"_",option_NEI,"_regular_common_",g,".RData",sep="")
		load(f1)
                sp=dimnames(cis$model$data)[[1]] #Studied species
                B=clean_matrix(cis)
		B_nodiag=B
		diag(B_nodiag)=NA
		B_nodiag[B_nodiag==0]=NA
		print(summary(c(abs(B_nodiag))))
                for(i in 2:length(sp)){
                	for(j in 1:(i-1)){
			points(B[i,j],B[j,i],pch=16,col=colo)
			}
		}
	}
}

text(0.2,0.2,"Mutualistic")
text(-0.2,-0.1,"Competitive")
text(-0.2,0.2,"Mixed")
#text(0.2,-0.2,"Mixed")
text(0.2,-0.1,"Mixed")
legend("bottomright",c("Brittany","Oléron","Arcachon","Mediterranean"),col=c("green","darkblue","cyan","darkred"),pch=16,bty="n")

draw.arc(0,0,0.02,0,2*pi,lwd=2)
#draw.arc(0,0,0.02,3*pi/2,2*pi,lwd=2)
dev.off()
