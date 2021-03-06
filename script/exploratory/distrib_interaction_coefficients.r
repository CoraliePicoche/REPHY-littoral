#2018/09/07 CP
#We want to have a look at the distributions of interaction coefficients
#2019/12/26 CP: look at the distribution of coefficients when pooled together

rm(list=ls())
graphics.off()
library(MARSS)
source("script/matrix_MAR_clean.r")

#option_model=c("unconstrained","pencen")
option_model=c("pencen")
option_NEI=c("null") #We don't take the "Not Elsewhere Identified" species into account
groupe=c("BZ","MO","AR","SU")
option_sp="common" #Species are the same 
take_0=FALSE #Do we consider the forced 0 in the interaction matrix when computing the mean values of inter-group competition?

all_intra=c()
all_inter=c()
for (m in 1:length(option_model)){
pdf(paste("article/graphe/distrib_comp_",option_model[m],"n0.pdf",sep=""))
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
	par(mfrow=c(length(option_lieu),2))
        for (ll in 1:length(option_lieu)){
                        f1=paste("data/analyse_MAR/",g,"/site_specific/",option_lieu[ll],"_",option_model[m],"_",option_NEI,"_regular_common_",g,".RData",sep="")
                        load(f1)
			B=clean_matrix(cis)
			id_diag=seq(1,dim(B)[1]*dim(B)[2],dim(B)[1]+1)
			hist(diag(B),xlab="",ylab=option_lieu[ll],main="")
			all_intra=c(all_intra,diag(B))
			if(ll==1){
				mtext("Intra",side=3)
			}
			B_nodiag=c(B)[-id_diag]
			if(!take_0){
				B_nodiag=B_nodiag[B_nodiag!=0]
			}
			hist(B_nodiag,xlab="",ylab="",main="")
			if(ll==1){
				mtext("Inter",side=3)
			}
			all_inter=c(all_inter,B_nodiag)
			
		}
	}
dev.off()
}
			

pdf(paste("article/graphe/distrib_comp_",option_model[m],"n0_allcoefftogether.pdf",sep=""),width=10,height=6)
par(mfrow=c(1,2))
hist(all_intra,xlab="",ylab="",breaks=10,main="Intra")
hist(all_inter,xlab="",ylab="",breaks=10,main="Inter")
dev.off()

print(shapiro.test(all_intra))
print(shapiro.test(all_inter))
