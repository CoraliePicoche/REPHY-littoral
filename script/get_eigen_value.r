#2018/09/05 CP 
#Computes eigen values of the interaction matrices

graphics.off()
rm(list=ls())
library(MARSS)
library("stringr")

option_model=c("unconstrained","pencen")
option_NEI=c("null")
groupe=c("BZ","MO","AR","SU")

id_g=0
pdf("Rapport/graphe/MAR_estimates/stability_vs_posi.pdf")
plot(0,0,t="n",xlim=c(0.2,0.6),ylim=c(0.4,0.65),xlab="%positive",ylab="max(eig)")
sink("data/analyse_MAR/eigenvalue_B.txt")
for (g in groupe){
        id_g=id_g+1
        if(g=="BZ"){
                option_lieu=c("Men er Roue","Loscolo","Croisic")
        }else if(g=="MO"){
                option_lieu=c("LEperon","Cornard","Auger")
        }else if(g=="SU"){
                option_lieu=c("Antoine","Lazaret")
        }else if(g=="AR"){
                option_lieu=c("Teychan","B7")
        }



for (l in 1:length(option_lieu)){
	print(option_lieu[l])
	for (m in 1:length(option_model)){
		print(option_model[m])
        	f1=paste("data/analyse_MAR/",g,"/site_specific/",option_lieu[l],"_",option_model[m],"_null_regular_common_",g,".RData",sep="")
                load(f1)
sp=dimnames(fit_log$model$data)[[1]]
B=matrix(0,length(sp),length(sp))
nom=dimnames(fit_log$par$B)[[1]]

for (n in 1:length(nom)){
        if(grepl("^\\(",nom[n])){ #default when using unconstrained or diagonal
                i=as.numeric(strsplit(nom[n],split="[\\(,\\)]")[[1]][2])
                j=as.numeric(strsplit(nom[n],split="[\\(,\\)]")[[1]][3])
                if(is.na(i)){ #i and j might not be numeric
                        i=which(unlist(lapply(sp,function(x) grepl(x,strsplit(nom[n],split="[\\(,\\)]")[[1]][2]))))
                        j=which(unlist(lapply(sp,function(x) grepl(x,strsplit(nom[n],split="[\\(,\\)]")[[1]][3]))))
                }
        }else{
                a=str_split(nom[n],sp)
                for (abis in 1:length(a)){
                        if(length(a[[abis]])==2){
                                if((a[[abis]][1]=="")&&(a[[abis]][2]!="")){
                                        j=abis
                                }
                                if((a[[abis]][2]=="")&&(a[[abis]][1]!="")){
                                        i=abis
                                }
                        }else if(length(a[[abis]])==3){
                                if((a[[abis]][2]=="")&&(a[[abis]][1]=="")&&(a[[abis]][3]=="")){
                                        j=abis
                                        i=abis
                                }
                        }
                }
        }
        if(i==j){
                B[i,j]=fit_log$par$B[n]-1
        }else{
                B[i,j]=fit_log$par$B[n]

        }
}


                eig_B=eigen(B)$values
		sum_pos=sum(B>0)/sum(B!=0)
                print(paste(max(Mod(eig_B)),sum_pos))
		
		if(m==1){
		cc="black"
		}else{
		cc="red"
		}
		points(sum_pos,max(Mod(eig_B)),pch=16,col=cc)
        }
}
}
sink()
legend("bottomleft",c("Unconstrained","Pencen"),col=c("black","red"),bty="n",pch=16)
dev.off()
