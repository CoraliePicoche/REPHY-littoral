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

id_lieu=0
for (g in groupe){
        if(g=="BZ"){
                option_lieu=c("Men er Roue","Loscolo")#,"Croisic") Croisic not available yet
        }else if(g=="MO"){
                option_lieu=c("LEperon","Cornard","Auger")
        }else if(g=="SU"){
                option_lieu=c("Antoine","Lazaret")
        }else if(g=="AR"){
                option_lieu=c("Teychan","B7")
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

temp=seq(-2,2,by=0.01)

colo=rainbow(10)
for(n in 1:length(name_species)){
	if(n%%9==1){
		par(mfrow=c(3,3))
	}
	a=max(tab_coeff[n,"TEMP",],na.rm=T)
	b=max(tab_coeff[n,"TEMP2",],na.rm=T)
#	yli1=min(a*temp+b*temp^2),
	yli2=max(a*temp+b*temp^2)
	plot(0,0,xlim=c(min(temp),max(temp)),ylim=c(-yli2,yli2))
	for(l in 1:10){
		yy=tab_coeff[n,"TEMP",l]*temp+tab_coeff[n,"TEMP2",l]*temp^2
		points(temp,yy,col=colo[l])
	}
}
