#2018/12/05 Link between stability (maximum eigen values) and complexity/metrics
#2019/10/31 Added mean value and variance of coefficients outside of the diagonal to answer R2

rm(list=ls())
graphics.off()
library('bipartite')
source("script/matrix_MAR_clean.r")
groupe=c("BZ","MO","SU","AR")
#option_model=c("unconstrained","pencen")
option_model=c("pencen")

results=array(NA,dim=c(10,12,length(option_model)),dimnames=list(1:10,c("stability","positive","weighted connectance","linkage density","vulnerability.LL","generality.HL","mutualism","all_neg","predation","E","V","SV"),option_model)) #10 places, all indices, 2 models (pencen and unconstrained)
#There can be no commensalism, at least i we're looking at B[i,j]>0 and B[j,i]==0

id_lieu=0
id_g=0
colo=c()
pch_sty=c()
for (g in groupe){
        if(g=="BZ"){
                option_lieu=c("Men er Roue","Loscolo","Croisic")
                colo=c(colo,rep("green",length(option_lieu)))
		pch_sty=c(pch_sty,rep(15,length(option_lieu)))
        }else if(g=="MO"){
                option_lieu=c("LEperon","Cornard","Auger")
                colo=c(colo,rep("darkblue",length(option_lieu)))
		pch_sty=c(pch_sty,rep(16,length(option_lieu)))
        }else if(g=="SU"){
                option_lieu=c("Antoine","Lazaret")
                colo=c(colo,rep("darkred",length(option_lieu)))
		pch_sty=c(pch_sty,rep(18,length(option_lieu)))
        }else if(g=="AR"){
                option_lieu=c("Teychan","B7")
                colo=c(colo,rep("cyan",length(option_lieu)))
		pch_sty=c(pch_sty,rep(17,length(option_lieu)))
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
		ll_abs=networklevel(web=abs(B),index=c("linkage density","weighted connectance","generality"),empty.web=F,logbase='2',level="both") #At first, I wanted all but 'quantitative' takes a VERY long time to compute
		results[id_lieu,"linkage density",option_model[m]]=ll_abs["linkage density"]
		results[id_lieu,"weighted connectance",option_model[m]]=ll_abs["weighted connectance"]
		results[id_lieu,"generality.HL",option_model[m]]=ll_abs["generality.HL"]
		results[id_lieu,"vulnerability.LL",option_model[m]]=ll_abs["vulnerability.LL"]
		mm=0
		nn=0
		pp=0
		upp_vec=c()
		low_vec=c()
		for(j in 2:dim(B)[1]){
			for(i in 1:(j-1)){
				if(B[i,j]!=0&B[j,i]!=0){
				upp_vec=c(upp_vec,B[i,j])
				low_vec=c(low_vec,B[j,i])
				}
				if(B[i,j]>0&B[j,i]>0){
					mm=mm+2
				}else if((B[i,j]<0&B[j,i]>0)||(B[i,j]>0&B[j,i]<0)){
					pp=pp+2
				}else if(B[i,j]<0&B[j,i]<0){
					nn=nn+2
				}
			
			}
		}
		rho=mean(upp_vec*low_vec)
		B_nodiag=B
		diag(B_nodiag)=NA
		B_nodiag_no0=B_nodiag
#		B_nodiag_no0[B_nodiag==0]=NA Add this line if you want the same numbers without the interactions set to 0
		results[id_lieu,"E",option_model[m]]=mean(c(B_nodiag_no0),na.rm=T)
		results[id_lieu,"V",option_model[m]]=var(c(B_nodiag_no0),na.rm=T)
		results[id_lieu,"SV",option_model[m]]=var(c(B_nodiag_no0),na.rm=T)*dim(B)[1]
		results[id_lieu,"mutualism",option_model[m]]=mm/(sum(B!=0)-dim(B)[1])
		results[id_lieu,"all_neg",option_model[m]]=nn/(sum(B!=0)-dim(B)[1])
		results[id_lieu,"predation",option_model[m]]=pp/(sum(B!=0)-dim(B)[1])

		

		S=dim(B)[1]
		tmp=max(sqrt(S*results[id_lieu,"V",option_model[m]])*(1+rho)-results[id_lieu,"E",option_model[m]],(S-1)*results[id_lieu,"E",option_model[m]])
		print(f1)
		print(rho)
		print(length(upp_vec))
		#print(cor.test(upp_vec,low_vec))
		print(tmp)
		print(mean(diag(B)))

	}
	}
}

if(1==0){
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
	pdf(paste("./article/graphe/complexity_stability_MainFig_",option_model[m],"justB_log2_differentsty.pdf",sep=""),width=12,height=5)
	par(mfrow=c(1,3),mar=c(2,2,0.75,0.5),oma=c(3,3,2,0.5),xpd=NA)

yli1=min(c(results[,"stability",option_model[m]]))
yli2=max(c(results[,"stability",option_model[m]]))

xli1=min(c(results[,"positive",option_model[m]]))*100
xli2=max(c(results[,"positive",option_model[m]]))*100
plot(results[,"positive",option_model[m]]*100,results[,"stability",option_model[m]],t="p",pch=pch_sty,cex=3,col=colo,ylim=c(yli1,yli2),xlim=c(xli1,xli2),xlab="% positive values",ylab="maximum eigenvalue",cex.axis=2,cex.lab=2,tck=-0.0075)
mtext("a)",side=3,cex=1.5,xpd=NA,font=2,line=1,adj=0)
legend('topleft',c("Brittany","Oléron","Arcachon","Mediterranean"),pch=c(15,16,17,18),col=c("green","darkblue","cyan","darkred"),bty="n",cex=2)

xli1=min(c(results[,"weighted connectance",option_model[m]]))
xli2=max(c(results[,"weighted connectance",option_model[m]]))
plot(results[,"weighted connectance",option_model[m]],results[,"stability",option_model[m]],t="p",pch=pch_sty,cex=3,col=colo,ylim=c(yli1,yli2),xlim=c(xli1,xli2),xlab="weighted connectance",ylab="",yaxt="n",cex.axis=2,cex.lab=2,tck=-0.0075)
mtext("b)",side=3,cex=1.5,xpd=NA,font=2,line=1,adj=0)

xli1=min(c(results[,"linkage density",option_model[m]]))
xli2=max(c(results[,"linkage density",option_model[m]]))
plot(results[,"linkage density",option_model[m]],results[,"stability",option_model[m]],t="p",pch=pch_sty,cex=3,col=colo,ylim=c(yli1,yli2),xlim=c(xli1,xli2),xlab="linkage density",ylab="",yaxt="n",cex.axis=2,cex.lab=2,tck=-0.0075)
mtext("c)",side=3,cex=1.5,xpd=NA,font=2,line=1,adj=0)

dev.off()
}

for (m in 1:length(option_model)){
        pdf(paste("./article/graphe/genvul_stability_FYI_",option_model[m],"justB_log2_differentsty.pdf",sep=""),width=12,height=5)
        par(mfrow=c(1,2),mar=c(2,2,0.5,0.5),oma=c(3,3,2,0.5),xpd=NA)

yli1=min(c(results[,"stability",option_model[m]]))
yli2=max(c(results[,"stability",option_model[m]]))

xli1=min(c(results[,"generality.HL",option_model[m]]))
xli2=max(c(results[,"generality.HL",option_model[m]]))
plot(results[,"generality.HL",option_model[m]],results[,"stability",option_model[m]],t="p",pch=pch_sty,cex=3,col=colo,ylim=c(yli1,yli2),xlim=c(xli1,xli2),xlab="generality",ylab="",cex.axis=2,cex.lab=2,tck=-0.0075)
mtext("b)",side=3,cex=1.5,xpd=NA,font=2,line=1,adj=0)

xli1=min(c(results[,"vulnerability.LL",option_model[m]]))
xli2=max(c(results[,"vulnerability.LL",option_model[m]]))
plot(results[,"vulnerability.LL",option_model[m]],results[,"stability",option_model[m]],t="p",pch=pch_sty,cex=3,col=colo,ylim=c(yli1,yli2),xlim=c(xli1,xli2),xlab="vulnerability",ylab="",yaxt="n",cex.axis=2,cex.lab=2,tck=-0.0075)
mtext("c)",side=3,cex=1.5,xpd=NA,font=2,line=1,adj=0)

dev.off()
}



print(results[,'linkage density','pencen'])
print(summary(results[,'stability','pencen']))
print(results[,'vulnerability.LL','pencen'])
print(results[,'generality.HL','pencen'])
}
#End 1==0


for (m in 1:length(option_model)){
	filename=paste("complexity_stability_MainFig_",option_model[m],"justB_with_E_and_SV_with0.pdf",sep="")
        pdf(paste("./article/graphe/",filename,sep=""),width=12,height=12)
        #par(mfrow=c(1,3),mar=c(2,2,0.75,0.5),oma=c(3,3,2,0.5),xpd=NA)
        par(mfrow=c(2,2),mar=c(4.5,2,1.5,0.5),oma=c(3,3,2,0.5),xpd=NA)

yli1=min(c(results[,"stability",option_model[m]]))
yli2=max(c(results[,"stability",option_model[m]]))

xli1=min(c(results[,"positive",option_model[m]]))*100
xli2=max(c(results[,"positive",option_model[m]]))*100
plot(results[,"positive",option_model[m]]*100,results[,"stability",option_model[m]],t="p",pch=pch_sty,cex=3,col=colo,ylim=c(yli1,yli2),xlim=c(xli1,xli2),xlab="% positive values",ylab=expression(paste("max(|",lambda,"|)",sep="")),cex.axis=2,cex.lab=2,tck=-0.0075)
mtext("a)",side=3,cex=1.5,xpd=NA,font=2,line=1,adj=0)
legend('topleft',c("Brittany","Oléron","Arcachon","Mediterranean"),pch=c(15,16,17,18),col=c("green","darkblue","cyan","darkred"),bty="n",cex=2)

xli1=min(c(results[,"weighted connectance",option_model[m]]))
xli2=max(c(results[,"weighted connectance",option_model[m]]))
plot(results[,"weighted connectance",option_model[m]],results[,"stability",option_model[m]],t="p",pch=pch_sty,cex=3,col=colo,ylim=c(yli1,yli2),xlim=c(xli1,xli2),xlab="weighted connectance",ylab="",yaxt="n",cex.axis=2,cex.lab=2,tck=-0.0075)
mtext("b)",side=3,cex=1.5,xpd=NA,font=2,line=1,adj=0)

xli1=min(c(results[,"E",option_model[m]]))
xli2=max(c(results[,"E",option_model[m]]))
plot(results[,"E",option_model[m]],results[,"stability",option_model[m]],t="p",pch=pch_sty,cex=3,col=colo,ylim=c(yli1,yli2),xlim=c(xli1,xli2),xlab="mean off-diagonal coefficients",ylab=expression(paste("max(|",lambda,"|)",sep="")),cex.axis=2,cex.lab=2,tck=-0.0075)
mtext("c)",side=3,cex=1.5,xpd=NA,font=2,line=1,adj=0)

xli1=min(c(results[,"SV",option_model[m]]))
xli2=max(c(results[,"SV",option_model[m]]))
plot(results[,"SV",option_model[m]],results[,"stability",option_model[m]],t="p",pch=pch_sty,cex=3,col=colo,ylim=c(yli1,yli2),xlim=c(xli1,xli2),xlab="S x off-diagonal coefficient variance ",ylab="",yaxt="n",cex.axis=2,cex.lab=2,tck=-0.0075)
mtext("d)",side=3,cex=1.5,xpd=NA,font=2,line=1,adj=0)

dev.off()
}
system(paste("cp article/graphe/",filename," article/submit_JEcol/response/",filename,sep=""))

