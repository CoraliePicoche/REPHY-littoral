#2018/12/31 CP: is there a link between self-regulation and abundance ?
#2019/12/26 CP: is there a link between interaction type (pos, neg?) with the others and abundance

rm(list=ls())
graphics.off()
library(MARSS)
source("script/matrix_MAR_clean.r")

#option_model=c("unconstrained","pencen")
option_model="pencen"
option_NEI=c("null") #We don't take the "Not Elsewhere Identified" species into account
groupe=c("BZ","MO","AR","SU")
option_sp="common" #Species are the same 

name_species=c("AST","CHA","DIT","GUI","LEP","NIT","PLE","PSE","RHI","SKE","THP","THL","GYM","PRO","PRP","SCR","CRY","EUG")

id_lieu=0

tab_tmp=array(NA,dim=c(10,length(name_species),8),dimnames=list(rep("",10),name_species,c("intra","abundance","mean_effect","mean_effect_V","mean_effect_abs","mean_effect_V_abs","sum_pos","sum_pos_V"))) #10 subsites, 21 species, 
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
                dimnames(tab_tmp)[[1]][id_lieu]=option_lieu[ll]
                f1=paste("data/analyse_MAR/",g,"/site_specific/",option_lieu[ll],"_",option_model,"_",option_NEI,"_regular_common_",g,".RData",sep="")
		if(g=="AR"){
			tab=read.csv(paste("./data/",option_lieu[ll],"_base.csv",sep=""),na.strings="NA",header=TRUE,sep=";",dec=".")
		}else{
			tab=read.table(paste("data/corres_hernandez_",option_lieu[ll],'.txt',sep=''),sep=";",na="NA",header=TRUE)
		}
                
	        load(f1)

                sp=dimnames(cis$model$data)[[1]] #Studied species
                B=clean_matrix(cis)
		B_nodiag=B
		diag(B_nodiag)=NA
		for(i in 1:length(sp)){
	        	tab_tmp[id_lieu,sp[i],"intra"]=B[i,i]
	        	tab_tmp[id_lieu,sp[i],"mean_effect"]=mean(B_nodiag[B_nodiag[,i]!=0,i],na.rm=T)
	        	tab_tmp[id_lieu,sp[i],"mean_effect_V"]=mean(B_nodiag[i,B_nodiag[i,]!=0],na.rm=T)
	        	tab_tmp[id_lieu,sp[i],"sum_pos"]=sum(B_nodiag[B_nodiag[,i]!=0,i]>0,na.rm=T)
	        	tab_tmp[id_lieu,sp[i],"sum_pos_V"]=sum(B_nodiag[i,B_nodiag[i,]!=0]>0,na.rm=T)
	        	tab_tmp[id_lieu,sp[i],"mean_effect_abs"]=mean(abs(B_nodiag[B_nodiag[,i]!=0,i]),na.rm=T)
	        	tab_tmp[id_lieu,sp[i],"mean_effect_V_abs"]=mean(abs(B_nodiag[i,B_nodiag[i,]!=0]),na.rm=T)
	                tab_tmp[id_lieu,sp[i],"abundance"]=mean(tab[,sp[i]],na.rm=T)
		}
	}
}

pdf("article/graphe/abundance_vs_inter.pdf",height=10)
par(mfrow=c(3,1),mar=c(4,4,1,1))
size_c=tab_tmp[,,"sum_pos"]/apply(tab_tmp[,,"sum_pos"],1,max,na.rm=T)
plot(log10(tab_tmp[,,"abundance"]),tab_tmp[,,"mean_effect"],pch=16,xlab="",ylab="Mean impact",col=colo,xaxt="n",cex=size_c+1)
legend("topleft",c("Brittany","Oléron","Arcachon","Mediterranean"),col=c("green","darkblue","cyan","darkred"),bty="n",pch=16,cex=1.5)
abline(h=0,lty=2)
y1=c(tab_tmp[,,"mean_effect"])
x1=c(log10(tab_tmp[,,"abundance"]))
lm1=lm(y1~x1)
sumlm1=summary(lm1)
abline(a=lm1$coefficients[1],b=lm1$coefficients[2])

size_c=tab_tmp[,,"sum_pos_V"]/apply(tab_tmp[,,"sum_pos_V"],1,max,na.rm=T)
plot(log10(tab_tmp[,,"abundance"]),tab_tmp[,,"mean_effect_V"],pch=16,xlab="",ylab="Mean vulnerability",col=colo,xaxt="n",cex=size_c+1)
y1=c(tab_tmp[,,"mean_effect_V"])
abline(h=0,lty=2)
lm1=lm(y1~x1)
sumlm1=summary(lm1)
abline(a=lm1$coefficients[1],b=lm1$coefficients[2])


plot(log10(tab_tmp[,,"abundance"]),tab_tmp[,,"intra"],pch=16,ylab="Self regulation",xlab="Log10 abundance",col=colo)
y1=c(tab_tmp[,,"intra"])
lm1=lm(y1~x1)
sumlm1=summary(lm1)
abline(a=lm1$coefficients[1],b=lm1$coefficients[2])
dev.off()

pdf("article/graphe/abundance_vs_inter_abs_value.pdf",height=10)
par(mfrow=c(3,1),mar=c(4,4,1,1))
size_c=tab_tmp[,,"sum_pos"]/apply(tab_tmp[,,"sum_pos"],1,max,na.rm=T)
plot(log10(tab_tmp[,,"abundance"]),tab_tmp[,,"mean_effect_abs"],pch=16,xlab="",ylab="Mean abs impact",col=colo,xaxt="n",cex=size_c+1)
y1=c(tab_tmp[,,"mean_effect_abs"])
x1=c(log10(tab_tmp[,,"abundance"]))
lm1=lm(y1~x1)
sumlm1=summary(lm1)
abline(a=lm1$coefficients[1],b=lm1$coefficients[2])

size_c=tab_tmp[,,"sum_pos_V"]/apply(tab_tmp[,,"sum_pos_V"],1,max,na.rm=T)
plot(log10(tab_tmp[,,"abundance"]),tab_tmp[,,"mean_effect_V_abs"],pch=16,xlab="",ylab="Mean abs vulnerability",col=colo,xaxt="n",cex=size_c+1)
y1=c(tab_tmp[,,"mean_effect_V_abs"])
lm1=lm(y1~x1)
sumlm1=summary(lm1)
abline(a=lm1$coefficients[1],b=lm1$coefficients[2])


plot(log10(tab_tmp[,,"abundance"]),tab_tmp[,,"intra"],pch=16,ylab="Self regulation",xlab="Log10 abundance",col=colo)
y1=c(tab_tmp[,,"intra"])
lm1=lm(y1~x1)
sumlm1=summary(lm1)
abline(a=lm1$coefficients[1],b=lm1$coefficients[2])
dev.off()


pdf("article/graphe/abundance_vs_regulation.pdf")
par(mfrow=c(1,1),mar=c(4.5,4.5,1,1))
plot(log10(tab_tmp[,,"abundance"]),tab_tmp[,,"intra"],pch=16,ylab="Self regulation",xlab="Log10 abundance",col=colo,cex=2,cex.lab=2,cex.axis=2)
abline(a=lm1$coefficients[1],b=lm1$coefficients[2])
legend("topleft",c("Brittany","Oléron","Arcachon","Mediterranean"),col=c("green","darkblue","cyan","darkred"),bty="n",pch=16,cex=1.5)
dev.off()
