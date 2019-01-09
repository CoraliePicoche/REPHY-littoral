#2018/09/12 Abiotic effects, for each species, for each site
rm(list=ls())
graphics.off()
library(MARSS)
library(Hmisc)

option_model="pencen"
option_NEI="null" #We don't take the "Not Elsewhere Identified" species into account
groupe=c("BZ","MO","AR","SU")
option_sp="common" #Species are the same 

name_species=c("AST","CHA","DIT","GUI","LEP","NIT","PLE","PSE","RHI","SKE","THP","THL","GYM","PRO","PRP","SCR","CRY","EUG")
covariate=c("TEMP","SAL")

tab_value=array(NA,dim=c(10,length(name_species),length(covariate)),dimnames=list(1:10,name_species,covariate)) #10 sites, 18 species, 2 covariates (TEMP and SALI)
tab_value_min=array(NA,dim=c(10,length(name_species),length(covariate)),dimnames=list(1:10,name_species,covariate)) #10 sites, 18 species, 2 covariates (TEMP and SALI)
tab_value_max=array(NA,dim=c(10,length(name_species),length(covariate)),dimnames=list(1:10,name_species,covariate)) #10 sites, 18 species, 2 covariates (TEMP and SALI)

id_lieu=0
colo=c()
a_pch=c()
name_places=c()
for (g in groupe){
        if(g=="BZ"){
                option_lieu=c("Men er Roue","Loscolo","Croisic")
                colo=c(colo,rep("green",length(option_lieu)))
		a_pch=c(a_pch,16,17,18)
		name_places=c(name_places,option_lieu)
        }else if(g=="MO"){
                option_lieu=c("LEperon","Cornard","Auger")
                colo=c(colo,rep("darkblue",length(option_lieu)))
		a_pch=c(a_pch,16,17,18)
		name_places=c(name_places,option_lieu)
        }else if(g=="SU"){
                option_lieu=c("Antoine","Lazaret")
                colo=c(colo,rep("darkred",length(option_lieu)))
		a_pch=c(a_pch,16,17)
		name_places=c(name_places,option_lieu)
        }else if(g=="AR"){
                option_lieu=c("Teychan","B7")
                colo=c(colo,rep("cyan",length(option_lieu)))
		a_pch=c(a_pch,16,17)
		name_places=c(name_places,option_lieu)
        }
        for (ll in 1:length(option_lieu)){
                id_lieu=id_lieu+1
                f1=paste("data/analyse_MAR/",g,"/site_specific/",option_lieu[ll],"_",option_model,"_",option_NEI,"_regular_common_",g,".RData",sep="")
                load(f1)
		name_all=rownames(cis$par$U)

		for (n in 1:length(name_all)){
                	i=which(unlist(lapply(name_species,function(x) grepl(x,strsplit(name_all[n],split="[\\(,\\)]")[[1]][2]))))
                        j=which(unlist(lapply(covariate,function(x) grepl(x,strsplit(name_all[n],split="[\\(,\\)]")[[1]][3]))))
			tab_value[id_lieu,i,j]=cis$par$U[n]
			tab_value_min[id_lieu,i,j]=cis$par.lowCI$U[n]
			tab_value_max[id_lieu,i,j]=cis$par.upCI$U[n]
		}
	}
}

#pdf("article/graphe/abiotic_first_try.pdf",width=15,height=10)
for(var in 1:length(covariate)){
mini=min(c(tab_value_min[,,var]),na.rm=T)
maxi=max(c(tab_value_max[,,var]),na.rm=T)
plot(0,0,t="n",ylim=c(mini,maxi),xlim=c(1.,18),xaxt="n",xlab="",ylab=covariate[var])
axis(1,at=1:18,lab=name_species)
for(i in 1:length(name_species)){
	for(j in 1:10){
		#points(i-0.2+(j-1)/10*0.4,tab_value[j,i,var],col=colo[j],pch=a_pch[j])
                errbar(i-0.3+(j-1)/10*0.6,tab_value[j,i,var],tab_value_min[j,i,var],tab_value_max[j,i,var],add=TRUE,col=colo[j],pch=a_pch[j],errbar.col=colo[j],cex=2)
	}
}
abline(h=0,lty=2)
abline(v=seq(0.5,18.5,by=1),lty=2)
}
#dev.off()

#pdf("article/graphe/abiotic_second_try.pdf",width=13,height=16)
par(mfcol=c(4,2),mar=c(6,5,4,2))
colo=c("blue","black","red")
pos="bottomright"
let=c('a)','c)','e)','g)','b)','d)','f)','h)')
let_id=0
mini=min(c(tab_value_min[,,]),na.rm=T)
maxi=max(c(tab_value_max[,,]),na.rm=T)
for(var in 1:length(covariate)){
site=1
for(j in 1:10){
	if(name_places[j] %in% c("Antoine","Teychan","LEperon","Men er Roue")){
		if(j>1&var==1){
			legend(pos,name_places[(j-site):(j-1)],col=colo,pch=16,bty="n",cex=2)
		}
		let_id=let_id+1
		sp=which(!is.na(tab_value[j,,1]))
		site=1
		plot(0,0,t="n",ylim=c(mini,maxi),xlim=c(1.,length(sp)+0.5),xaxt="n",xlab="",ylab=covariate[var],cex.axis=2,cex.lab=2)
		mtext(let[let_id],side=3,cex=2,xpd=NA,font=2,line=.75,adj=0)
		axis(1,at=1:length(sp),lab=name_species[sp],cex.axis=2,las=2)
	}else{
		site=site+1
	}
	for(i in 1:length(sp)){
                #points(i-0.2+(j-1)/10*0.4,tab_value[j,i,var],col=colo[j],pch=a_pch[j])
                errbar(i-0.2+site/3*0.6,tab_value[j,sp[i],var],tab_value_min[j,sp[i],var],tab_value_max[j,sp[i],var],add=TRUE,col=colo[site],pch=16,errbar.col=colo[site],cex=2)
        }
abline(h=0,lty=2)
}
if(var==1){
legend(pos,name_places[(j-site+1):(j)],col=colo,pch=16,bty="n",cex=2)
}
}
#dev.off()

#Count positive and significant parameters
var="SAL"
signif_val=rep(0,10)
posi_signif_val=rep(0,10) #POsitive AND significant values
posi_val=rep(0,10) #Positive values, significant or not
val=rep(0,10)
for(j in 1:10){
	sp=which(!is.na(tab_value[j,,1]))
        for(i in 1:length(sp)){
		if(tab_value_min[j,sp[i],var]*tab_value_max[j,sp[i],var]>0){
			signif_val[j]=signif_val[j]+1
			if(tab_value[j,sp[i],var]>0){
				posi_signif_val[j]=posi_signif_val[j]+1
			}
		}
		if(tab_value[j,sp[i],var]>0){
			posi_val[j]=posi_val[j]+1
		}
	}
	val[j]=mean(abs(tab_value[j,,var]),na.rm=T)
	posi_signif_val[j]=posi_signif_val[j]/signif_val[j]
	signif_val[j]=signif_val[j]/(length(sp))
	posi_val[j]=posi_val[j]/(length(sp))
}
