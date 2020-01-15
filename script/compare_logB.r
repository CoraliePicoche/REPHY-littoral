#######
#06/12/2020 CP: compare B-I eigenvalue spectrum to log(B) eigenvalue spectrum
######

source("script/matrix_MAR_clean.r")
library("expm")
groupe=c("BZ","MO","SU","AR")

#name_species=c("AST","CHA","DIT","GUI","LEP","NIT","PLE","PSE","RHI","SKE","THP","THL","GYM","PRO","PRP","SCR","CRY","EUG")

mat_Bprim=matrix(NA,10,2)
mat_Blog=matrix(NA,10,2)
mat_B=matrix(NA,10,2)
SV=rep(NA,10)

id_lieu=0

colo=c()
pch_sty=c()
pdf("article/submit_JEcol/response_R2/compare_spectra.pdf")
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
                pch_sty=c(pch_sty,rep(17,length(option_lieu)))
        }else if(g=="AR"){
                option_lieu=c("Teychan","B7")
                colo=c(colo,rep("cyan",length(option_lieu)))
                pch_sty=c(pch_sty,rep(18,length(option_lieu)))
        }
        for (ll in 1:length(option_lieu)){
                id_lieu=id_lieu+1
                f1=paste("data/analyse_MAR/",g,"/site_specific/",option_lieu[ll],"_pencen_null_regular_common_",g,".RData",sep="")
                load(f1)

                sp=dimnames(cis$model$data)[[1]] #Studied species
                B=clean_matrix(cis,diag_bool=F) #We do not remove the identity matrix of B, we can do it after

		B_prim=B-diag(nrow(B)) #B-I
		B_log=logm(B)

		eig_Bprim=eigen(B_prim)$values
		eig_Blog=eigen(B_log)$values
		eig_B=eigen(B)$values
	
		#Max modulus
                id=which(Mod(eig_B)==max(Mod(eig_B)))
                mat_B[id_lieu,1]=Mod(eig_B[id][1])
                #Max real part
                id=which(Re(eig_B)==max(Re(eig_B)))
                mat_B[id_lieu,2]=Re(eig_B[id][1])
	
		plot(0,0,t="n",xlim=c(-1,-0.10),ylim=c(-0.25,0.25),xlab="Re",ylab="Im",main=option_lieu[ll])
		points(Re(eig_Bprim),Im(eig_Bprim),col="black",pch=16,cex=2)
                #Max modulus
		id=which(Mod(eig_Bprim)==max(Mod(eig_Bprim)))
                points(Re(eig_Bprim[id]),Im(eig_Bprim[id]),col="red",cex=2.2,lwd=2)
		mat_Bprim[id_lieu,1]=Mod(eig_Bprim[id][1])
                #Max real part
		id=which(Re(eig_Bprim)==max(Re(eig_Bprim)))
                points(Re(eig_Bprim[id]),Im(eig_Bprim[id]),col="red",cex=2.2,pch="+")
		mat_Bprim[id_lieu,2]=Re(eig_Bprim[id][1])
                

		points(Re(eig_Blog),Im(eig_Blog),col="grey",pch=16,cex=2)
                #Max modulus
		id=which(Mod(eig_Blog)==max(Mod(eig_Blog)))
                points(Re(eig_Blog[id]),Im(eig_Blog[id]),col="orange",cex=2.2,lwd=2)
		mat_Blog[id_lieu,1]=Mod(eig_Blog[id][1])
                #Max real part
		id=which(Re(eig_Blog)==max(Re(eig_Blog)))
                points(Re(eig_Blog[id]),Im(eig_Blog[id]),col="orange",cex=2.2,pch="+")
		mat_Blog[id_lieu,2]=Re(eig_Blog[id][1])

		legend("topleft",c("B-I","log(B)"),col=c("black","grey"),bty="n",pch=16,cex=2)

	
		B_nodiag=B
		diag(B_nodiag)=NA
		SV[id_lieu]=var(c(B_nodiag),na.rm=T)*dim(B)[1]


	}
}
dev.off()

pdf("article/submit_JEcol/response_R2/maxRe_f_maxMod.pdf")
plot(0,0,xlim=c(-0.6,-0.1),ylim=c(0.3,1),t="n",xlab="Re(eig)",ylab="Mod(eig)")
points(mat_Bprim[,2],mat_Bprim[,1],col="black",pch=16)
points(mat_Blog[,2],mat_Blog[,1],col="grey",pch=16)
legend("topleft",c("B-I","log(B)"),col=c("black","grey"),bty="n",pch=16,cex=2)
dev.off()

pdf("article/submit_JEcol/response_R2/eigB_f_eigBprim.pdf")
par(mfrow=c(1,2))
plot(0,0,xlim=c(-1,0),ylim=c(0,1),t="n",xlab="B-I or Blog",ylab="B",main="Max Re")
points(mat_Bprim[,2],mat_B[,2],col="black",pch=16)
points(mat_Blog[,2],mat_B[,2],col="grey",pch=16)
legend("topleft",c("B-I","log(B)"),col=c("black","grey"),bty="n",pch=16,cex=2)
plot(0,0,xlim=c(0,1),ylim=c(0,1),t="n",xlab="B-I or Blog",ylab="B",main="Max Mod")
points(mat_Bprim[,1],mat_B[,1],col="black",pch=16)
points(mat_Blog[,1],mat_B[,1],col="grey",pch=16)
dev.off()

pdf("article/submit_JEcol/response_R2/modB_f_maxReA_singlepanel.pdf")
par(mfrow=c(1,1))
plot(0,0,xlim=c(-0.5,-0.1),ylim=c(0.6,0.85),t="n",xlab=expression(paste("max(Re(",lambda[A],"))",sep="")),ylab=expression(paste("max(|",lambda[B],"|)",sep="")))
points(mat_Bprim[,2],mat_B[,1],col="black",pch=16)
points(mat_Blog[,2],mat_B[,1],col="grey",pch=16)
legend("bottomright",c("A=B-I","A=log(B)"),col=c("black","black"),bty="n",pch=c(16,1))
dev.off()


pdf("article/submit_JEcol/response_R2/modB_f_maxReA.pdf",width=8,height=4)
par(mfrow=c(1,2),mar=c(4,4.5,2,1.5))
plot(0,0,xlim=c(-0.5,-0.1),ylim=c(0.6,0.85),t="n",xlab=expression(paste("max(Re(",lambda[A],"))",sep="")),ylab=expression(paste("max(|",lambda[B],"|)",sep="")))
print(mean(abs(mat_Bprim[,2]-mat_Blog[,2]))/mean(abs(mat_Bprim[,2])))
print(mean(abs(mat_Bprim[,2]-mat_Blog[,2]))/mean(abs(mat_Blog[,2])))
points(mat_Bprim[,2],mat_B[,1],col="black",pch=16)
points(mat_Blog[,2],mat_B[,1],col="black",pch=1,lwd=1.5)
print(paste("B-I",cor(mat_Bprim[,2],mat_B[,1])))
print(paste("log(B)",cor(mat_Blog[,2],mat_B[,1])))
legend("bottomright",c("A=B-I","A=log(B)"),col=c("black","black"),bty="n",pch=c(16,1))
mtext(" a)",side=3,line=-1,adj=0)

plot(SV,mat_Bprim[,2],t="p",pch=pch_sty,col=colo,xlab="S x off-diagonal coefficient variance",ylab=expression(paste("max(Re(",lambda[A],"))",sep="")),ylim=range(c(mat_Bprim[,2],mat_Blog[,2])))
pch_sty_2=pch_sty-15
pch_sty_2[pch_sty_2==3]=5
points(SV,mat_Blog[,2],col=colo,pch=pch_sty_2,lwd=1.5)
mtext(" b)",side=3,line=-1,adj=0)

dev.off()
