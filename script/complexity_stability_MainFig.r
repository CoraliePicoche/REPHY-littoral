#2018/12/05 Link between stability (maximum eigen values) and complexity/metrics

rm(list=ls())
graphics.off()
library('bipartite')
source("script/matrix_MAR_clean.r")
groupe=c("BZ","MO","SU","AR")
option_model=c("unconstrained","pencen")

results=array(NA,dim=c(10,7,length(option_model)),dimnames=list(1:10,c("stability","positive","weighted connectance","linkage density","mutualism","all_neg","predation"),option_model)) #10 places, 4 indices (max eigen values, positive values, weighted connectance, linkage density), 2 models (pencen and unconstrained)
#There can be no commensalism, at least i we're looking at B[i,j]>0 and B[j,i]==0

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

	for (m in 1:length(option_model)){
	        f1=paste("data/analyse_MAR/",g,"/site_specific/",option_lieu[ll],"_",option_model[m],"_null_regular_common_",g,".RData",sep="")
        	load(f1)
	        B=clean_matrix(cis,diag_bool=F)
		results[id_lieu,"stability",option_model[m]]=max(Mod(eigen(B)$values))
		results[id_lieu,"positive",option_model[m]]=(sum(B>0)-dim(B)[1])/(sum(B!=0)-dim(B)[1]) #We only take into account interspecific interactions
		ll_abs=networklevel(web=abs(B),index=c("linkage density","weighted connectance"),empty.web=F) #At first, I wanted all but 'quantitative' takes a VERY long time to compute
		results[id_lieu,"linkage density",option_model[m]]=ll_abs["linkage density"]
		results[id_lieu,"weighted connectance",option_model[m]]=ll_abs["weighted connectance"]
		mm=0
		nn=0
		pp=0
		for(j in 2:dim(B)[1]){
			for(i in 1:(j-1)){
				if(B[i,j]>0&B[j,i]>0){
					mm=mm+2
				}else if((B[i,j]<0&B[j,i]>0)||(B[i,j]>0&B[j,i]<0)){
					pp=pp+2
				}else if(B[i,j]<0&B[j,i]<0){
					nn=nn+2
				}
			}
		}
		results[id_lieu,"mutualism",option_model[m]]=mm/(sum(B!=0)-dim(B)[1])
		results[id_lieu,"all_neg",option_model[m]]=nn/(sum(B!=0)-dim(B)[1])
		results[id_lieu,"predation",option_model[m]]=pp/(sum(B!=0)-dim(B)[1])
	}
	}
}


pdf(paste("./article/graphe/complexity_stability_MainFig_compare_unconstrained_pencen_justB_log2.pdf",sep=""),width=12,height=5)
par(mfrow=c(1,3),mar=c(2,2,0.5,0.5),oma=c(3,3,2,0.5),xpd=NA)

yli1=min(c(results[,"stability",]))
yli2=max(c(results[,"stability",]))

xli1=min(c(results[,"positive",]))*100
xli2=max(c(results[,"positive",]))*100
plot(results[,"positive","unconstrained"]*100,results[,"stability","unconstrained"],t="p",pch=16,cex=2,col="black",ylim=c(yli1,yli2),xlim=c(xli1,xli2),xlab="% positive values",ylab="stability",cex.axis=2,cex.lab=2)
points(results[,"positive","pencen"]*100,results[,"stability","pencen"],pch=16,cex=2,col="red")
legend("bottomleft",c("unconstrained","pencen"),col=c("black","red"),pch=16,cex=2,bty="n")

xli1=min(c(results[,"weighted connectance",]))
xli2=max(c(results[,"weighted connectance",]))
plot(results[,"weighted connectance","unconstrained"],results[,"stability","unconstrained"],t="p",pch=16,cex=2,col="black",ylim=c(yli1,yli2),xlim=c(xli1,xli2),xlab="weighted connectance",ylab="",yaxt="n",cex.axis=2,cex.lab=2)
points(results[,"weighted connectance","pencen"],results[,"stability","pencen"],pch=16,cex=2,col="red")

xli1=min(c(results[,"linkage density",]))
xli2=max(c(results[,"linkage density",]))
plot(results[,"linkage density","unconstrained"],results[,"stability","unconstrained"],t="p",pch=16,cex=2,col="black",ylim=c(yli1,yli2),xlim=c(xli1,xli2),xlab="linkage density",ylab="",yaxt="n",cex.axis=2,cex.lab=2)
points(results[,"linkage density","pencen"],results[,"stability","pencen"],pch=16,cex=2,col="red")

dev.off()

for (m in 1:length(option_model)){
	pdf(paste("./article/graphe/complexity_stability_MainFig_",option_model[m],"justB_log2.pdf",sep=""),width=12,height=5)
	par(mfrow=c(1,3),mar=c(2,2,0.5,0.5),oma=c(3,3,2,0.5),xpd=NA)

yli1=min(c(results[,"stability",option_model[m]]))
yli2=max(c(results[,"stability",option_model[m]]))

xli1=min(c(results[,"positive",option_model[m]]))*100
xli2=max(c(results[,"positive",option_model[m]]))*100
plot(results[,"positive",option_model[m]]*100,results[,"stability",option_model[m]],t="p",pch=16,cex=3,col=colo,ylim=c(yli1,yli2),xlim=c(xli1,xli2),xlab="% positive values",ylab="maximum eigenvalue",cex.axis=2,cex.lab=2,tck=-0.0075)
mtext("a)",side=3,cex=1.5,xpd=NA,font=2,line=1,adj=0)
legend('topleft',c("Brittany","Oléron","Arcachon","Mediterranean"),pch=16,col=c("green","darkblue","cyan","darkred"),bty="n",cex=2)

xli1=min(c(results[,"weighted connectance",option_model[m]]))
xli2=max(c(results[,"weighted connectance",option_model[m]]))
plot(results[,"weighted connectance",option_model[m]],results[,"stability",option_model[m]],t="p",pch=16,cex=3,col=colo,ylim=c(yli1,yli2),xlim=c(xli1,xli2),xlab="weighted connectance",ylab="",yaxt="n",cex.axis=2,cex.lab=2,tck=-0.0075)
mtext("b)",side=3,cex=1.5,xpd=NA,font=2,line=1,adj=0)

xli1=min(c(results[,"linkage density",option_model[m]]))
xli2=max(c(results[,"linkage density",option_model[m]]))
plot(results[,"linkage density",option_model[m]],results[,"stability",option_model[m]],t="p",pch=16,cex=3,col=colo,ylim=c(yli1,yli2),xlim=c(xli1,xli2),xlab="linkage density",ylab="",yaxt="n",cex.axis=2,cex.lab=2,tck=-0.0075)
mtext("c)",side=3,cex=1.5,xpd=NA,font=2,line=1,adj=0)

dev.off()
}

print(results[,'linkage density','pencen'])
