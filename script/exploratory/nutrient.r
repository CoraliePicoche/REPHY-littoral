#25/01/2019 Checking nutrient concentrations when we can (see limitations in Burson et al. 2018 Ecology)

rm(list=ls())
graphics.off()

#option_model=c("unconstrained","pencen")
groupe=c("BZ","MO","AR")

pdf("article/graphe/nutrients.pdf",width=15)
for (g in groupe){
        if(g=="BZ"){
                option_lieu=c("Men er Roue","Loscolo")
        }else if(g=="MO"){
                option_lieu=c("Cornard","Auger")
        }else if(g=="AR"){
                option_lieu=c("Teychan","B7")
        }
        for (ll in 1:length(option_lieu)){
		plot(0,0,t="n",xlim=c(as.Date("1996-01-01"),as.Date("2015-01-01")),ylim=c(0,150),xlab="",ylab="",main=option_lieu[ll])
                if(g=="AR"){
                        tab=read.csv(paste("./data/",option_lieu[ll],"_base.csv",sep=""),na.strings="NA",header=TRUE,sep=";",dec=".")
			n_var=tab$NOx+tab$NH4
			p_var=tab$PHOS
			s_var=tab$SI
                }else{
                        tab=read.table(paste("data/",option_lieu[ll],'hydro.txt',sep=''),sep=";",na="NA",header=TRUE)
			if(g=="BZ"){
				tmp_n=cbind(tab$NO2,tab$NO3,tab$NO3.NO2)
				whichn=which(is.na(tab$NO2)&is.na(tab$NO3)&is.na(tab$NO3.NO2))
				n_var=apply(tmp_n,1,sum,na.rm=T)
				n_var[whichn]=NA
				n_var=n_var+tab$NH4
			}else{
				n_var=tab$NO3.NO2+tab$NH4
			}
			p_var=tab$PO4
			s_var=tab$SIOH
                }
		points(as.Date(tab$Date),n_var,col="red",pch=15+ll,lty=1,t="o")
		points(as.Date(tab$Date),s_var,col="grey",pch=15+ll,lty=1,t="o")
		par(new=T)
		plot(as.Date(tab$Date),p_var,col="blue",pch=15+ll,t="o",xaxt="n",yaxt="n",lty=1,xlab="",ylab="")
		legend("topleft",as.character(format(c(mean(n_var,na.rm=T),mean(p_var,na.rm=T),mean(s_var,na.rm=T)),digits=2)),bty="n",pch=c("N","P","S"),col=c("red","blue","grey"))
	}
}

dev.off()

