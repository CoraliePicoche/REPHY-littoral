#2018/09/14 CP
#Had a look at the responses to polynomial model of covariates
#TBC!!

rm(list=ls())
graphics.off()

option_model=c("pencen")
option_NEI="null" #We don't take the "Not Elsewhere Identified" species into account
groupe=c("BZ","MO","AR","SU")
option_sp="common" #Species are the same 
name_species=c("AST","CHA","DIT","GUI","LEP","NIT","PLE","PSE","RHI","SKE","THP","THL","GYM","PRO","PRP","SCR","CRY","EUG")

coeff=c("TEMP","TEMP2","SALI","SALI2")
tab_coeff=array(NA,dim=c(length(name_species),length(coeff),10))
tab_signif=array(NA,dim=c(length(name_species),length(coeff),10))
rownames(tab_coeff)=name_species
colnames(tab_coeff)=coeff
rownames(tab_signif)=name_species
colnames(tab_signif)=coeff

id_lieu=0
lieu=c()
colo=c()
for (g in groupe){
        if(g=="BZ"){
                option_lieu=c("Men er Roue","Loscolo","Croisic")
		lieu=c(lieu,option_lieu)
		colo=c(colo,"green","darkgreen","mediumseagreen")
		coeff=c("TEMP","TEMP2","SALI","SALI2")
        }else if(g=="MO"){
                option_lieu=c("LEperon","Cornard","Auger")
		lieu=c(lieu,option_lieu)
		colo=c(colo,"darkblue","navyblue","mediumblue")
		coeff=c("TEMP","TEMP2","SALI","SALI2")
        }else if(g=="SU"){
                option_lieu=c("Antoine","Lazaret")
		lieu=c(lieu,option_lieu)
		colo=c(colo,"red","red4")
		coeff=c("TEMP","TEMP2","SALI","SALI2")
        }else if(g=="AR"){
                option_lieu=c("Teychan","B7")
		lieu=c(lieu,option_lieu)
		colo=c(colo,"cyan","lightskyblue")
		coeff=c("TEMP","TEMP2","SAL","SAL2")
        }
        for (ll in 1:length(option_lieu)){
		id_lieu=id_lieu+1
                f2=paste("data/analyse_MAR/",g,"/site_specific/",option_lieu[ll],"_",option_model,"_",option_NEI,"_regular_common_ALL_squared.RData",sep="")
                load(f2)
		name_C=rownames(cis$par$U)
		for(n in 1:length(name_species)){
			for(c in 1:length(coeff)){
				id=which(grepl(name_species[n],name_C)&grepl(coeff[c],name_C))
				if(length(id)>0){
				tab_coeff[n,c,id_lieu]=cis$par$U[id[1]]
				if(cis$par.upCI$U[id[1]]*cis$par.lowCI$U[id[1]]>0){
					tab_signif[n,c,id_lieu]=cis$par$U[id[1]]
					}
				}
			}
		}
	}
}

xmax=3
xmin=-3
pdf("Rapport/graphe/MAR_estimates/niche_temperature.pdf")
#colo=rainbow(10)
for(n in 1:length(name_species)){
	if(n%%9==1){
		par(mfrow=c(3,3),mar=c(2,2,2,2))
	}
	coef_abs=max(abs(tab_coeff[n,"TEMP",]),na.rm=T)
	coef_abs2=max(abs(tab_coeff[n,"TEMP2",]),na.rm=T)
	yli2=xmax*coef_abs+xmax^2*coef_abs2
	plot(0,0,xlim=c(xmin,xmax),ylim=c(-yli2,yli2),main=name_species[n],t="n",xlab="",ylab="")
	for(l in 1:10){
		if(lieu[l] %in% c("Teychan","B7")){
                	tab_cov=read.table(paste("data/",lieu[l],'_base.csv',sep=''),sep=";",na="NA",header=TRUE)
		}else{
                	tab_cov=read.table(paste("data/",lieu[l],'hydro.txt',sep=''),sep=";",na="NA",header=TRUE)
		}
		mm=mean(tab_cov[,'TEMP'],na.rm=T)
		ss=sd(tab_cov[,'TEMP'],na.rm=T)
		xx=sort((tab_cov[,'TEMP']-mm)/ss)
		yy=tab_coeff[n,"TEMP",l]*xx+tab_coeff[n,"TEMP2",l]*xx^2
#		points(temp,yy,col=colo[l]) #We'll use dots if we decide to go for real temperature values #Or not
#		points(xx,yy,col=colo[l],pch="+")
                lines(xx,yy,col=colo[l],lty=2)
		if(!(is.na(tab_signif[n,"TEMP2",l]))){
#			points(xx,yy,col=colo[l],pch=16)
                	lines(xx,yy,col=colo[l],lty=1)
		}
	}
}
dev.off()


pdf("Rapport/graphe/MAR_estimates/niche_salinity.pdf")
#colo=rainbow(10)
for(n in 1:length(name_species)){
        if(n%%9==1){
                par(mfrow=c(3,3),mar=c(2,2,2,2))
        }
        coef_abs=max(abs(tab_coeff[n,3,]),na.rm=T)
        coef_abs2=max(abs(tab_coeff[n,4,]),na.rm=T)
	yli2=xmax*coef_abs+xmax^2*coef_abs2
        plot(0,0,xlim=c(xmin,xmax),ylim=c(-yli2,yli2),main=name_species[n],t="n",xlab="",ylab="")
        for(l in 1:10){
		if(lieu[l] %in% c("Teychan","B7")){
			name_sali="SAL"
                	tab_cov=read.table(paste("data/",lieu[l],'_base.csv',sep=''),sep=";",na="NA",header=TRUE)
		}else{
			name_sali="SALI"
                	tab_cov=read.table(paste("data/",lieu[l],'hydro.txt',sep=''),sep=";",na="NA",header=TRUE)
		}
                mm=mean(tab_cov[,name_sali],na.rm=T)
                ss=sd(tab_cov[,name_sali],na.rm=T)
                xx=sort((tab_cov[,name_sali]-mm)/ss)
                yy=tab_coeff[n,3,l]*xx+tab_coeff[n,4,l]*xx^2
#               points(temp,yy,col=colo[l]) #We'll use dots if we decide to go for real temperature values #Nah, dots are not easy to read, and they're much heavier
#                points(xx,yy,col=colo[l],pch="+")
                lines(xx,yy,col=colo[l],lty=2)
                if(!(is.na(tab_signif[n,4,l]))){
#                        points(xx,yy,col=colo[l],pch=16)
                	lines(xx,yy,col=colo[l],lty=1)
                }
        }
}
dev.off()
