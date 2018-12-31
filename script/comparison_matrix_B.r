#2018/12/31 CP Having a doubt on the results we have on network metrics with |B|

rm(list=ls())
graphics.off()
library('bipartite')
source("script/matrix_MAR_clean.r")
groupe=c("BZ","MO","SU","AR")
option_model=c("pencen")

results=array(NA,dim=c(10,3,3),dimnames=list(1:10,c("weighted connectance","linkage density","mean"),c("abs","pos","neg")))

id_lieu=0
id_g=0
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
                rownames(results)[id_lieu]=option_lieu[ll]

                f1=paste("data/analyse_MAR/",g,"/site_specific/",option_lieu[ll],"_",option_model,"_null_regular_common_",g,".RData",sep="")
                load(f1)
                B=clean_matrix(cis,diag_bool=F)
		B_nodiag=B
		diag(B_nodiag)=0
		B_pos=B_nodiag
		B_neg=B_nodiag
		B_pos[B_nodiag<0]=0 #Keep only value above 0
		B_neg[B_nodiag>0]=0 #Keep only value below 0

	
                ll_abs=networklevel(web=abs(B_nodiag),index=c("linkage density","weighted connectance"),empty.web=F) 
                ll_pos=networklevel(web=abs(B_pos),index=c("linkage density","weighted connectance"),empty.web=F) 
                ll_neg=networklevel(web=abs(B_neg),index=c("linkage density","weighted connectance"),empty.web=F) 


                results[id_lieu,"linkage density","abs"]=ll_abs["linkage density"]
                results[id_lieu,"linkage density","pos"]=ll_pos["linkage density"]
                results[id_lieu,"linkage density","neg"]=ll_neg["linkage density"]
                results[id_lieu,"weighted connectance","abs"]=ll_abs["weighted connectance"]
                results[id_lieu,"weighted connectance","pos"]=ll_pos["weighted connectance"]
                results[id_lieu,"weighted connectance","neg"]=ll_neg["weighted connectance"]
                

		B_nodiag=B
		diag(B_nodiag)=NA
                B_pos=B_nodiag
		B_neg=B_nodiag
                B_pos[B_nodiag<=0]=NA #Keep only value above 0
                B_neg[B_nodiag>=0]=NA #Keep only value below 0

		results[id_lieu,"mean","abs"]=mean(abs(B_nodiag),na.rm=T)
                results[id_lieu,"mean","pos"]=mean(abs(B_pos),na.rm=T)
                results[id_lieu,"mean","neg"]=mean(abs(B_neg),na.rm=T)
	}
}

pdf("article/graphe/test_on_coefficient_sign.pdf",width=12)
par(mfrow=c(2,2))

tmp=c(results[,'weighted connectance',c("pos","neg")])
ylimi=c(min(tmp,na.rm=T),max(tmp,na.rm=T))
plot(results[,"weighted connectance","abs"],results[,"weighted connectance","pos"],pch=16,col="blue",ylab="",xlab="abs B",main="Weighted conn",ylim=ylimi,cex=2)
points(results[,"weighted connectance","abs"],results[,"weighted connectance","neg"],pch=16,col="red",cex=2)
abline(b=1,a=0,lty=2)
legend("topleft",c("|B>0|","|B<0|"),col=c("blue","red"),pch=16,bty="n")

tmp=c(results[,'linkage density',c("pos","neg")])
ylimi=c(min(tmp,na.rm=T),max(tmp,na.rm=T))
plot(results[,"linkage density","abs"],results[,"linkage density","pos"],pch=16,col="blue",ylab="",xlab="abs B",main="Linkage density",ylim=ylimi,cex=2)
points(results[,"linkage density","abs"],results[,"linkage density","neg"],pch=16,col="red",cex=2)
abline(b=1,a=0,lty=2)

tmp=c(results[,'mean',c("pos","neg")])
ylimi=c(min(tmp,na.rm=T),max(tmp,na.rm=T))
plot(results[,"mean","abs"],results[,"mean","pos"],pch=16,col="blue",ylab="",xlab="abs B",main="Mean",ylim=ylimi,cex=2)
points(results[,"mean","abs"],results[,"mean","neg"],pch=16,col="red",cex=2)
abline(b=1,a=0,lty=2)

plot(results[,"mean","neg"],results[,"mean","pos"],pch=16,col="black",ylab="pos B",xlab="neg B",main="Mean",ylim=ylimi,cex=2)
abline(b=1,a=0,lty=2)
dev.off()
