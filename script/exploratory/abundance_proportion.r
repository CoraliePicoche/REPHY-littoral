#05/02/2019 Checking the proportion of cells we take into account per site
#24/10/2019 Computing abundance of dinoflagellates vs abundance of diatoms

library("lubridate")
rm(list=ls())
graphics.off()

groupe=c("BZ","MO","SU")

tab_sp=read.table('data/lieu_sp_post_reconstruct_pour_MAR.csv',header=TRUE,na.strings="",sep=";")
lieu=colnames(tab_sp)
lieu=gsub('.',' ',lieu,fixed=TRUE) #useful for Men er Roue

pdf("article/graphe/abundance.pdf",width=15)
for (g in groupe){
        if(g=="BZ"){
                option_lieu=c("Men er Roue","Loscolo","Croisic")
        }else if(g=="MO"){
                option_lieu=c("Cornard","Auger","LEperon")
        }else if(g=="SU"){
                option_lieu=c("Antoine","Lazaret")
        }
        a=table(as.matrix(tab_sp[,lieu %in% option_lieu]))
        liste_sp=dimnames(a[a==length(option_lieu)])[[1]]
        if("NEI"%in%liste_sp){
        	liste_sp=liste_sp[-which(liste_sp=="NEI")]
        }

        for (ll in 1:length(option_lieu)){
        	tab=read.table(paste("data/corres_hernandez_",option_lieu[ll],'.txt',sep=''),sep=";",na="NA",header=TRUE)
		dd=as.Date(tab$Date)
        	tab=tab[year(dd)>=1996,]#Using data from 1996
		dd=dd[year(dd)>=1996]
		
		biom_in=apply(tab[,liste_sp],1,sum,na.rm=T)
		biom_tot=apply(tab[,2:dim(tab)[2]],1,sum,na.rm=T)
		biom_prop=biom_in/biom_tot

		plot(dd,log10(biom_tot+1),pch=16,xlim=c(as.Date("1995-01-01"),as.Date("2015-01-01")),xlab="",ylab="",main=option_lieu[ll],col="red",t="o")
                par(new=T)
                plot(dd,biom_prop,col="blue",pch=16,t="p",xaxt="n",yaxt="n",lty=1,xlab="",ylab="",xlim=c(as.Date("1995-01-01"),as.Date("2015-01-01")))
		axis(4)
                legend("bottomleft",c("Log10 Total abundance","Proportion used"),pch=16,bty="n",col=c("red","blue"))
		legend("topleft",as.character(format(c(min(biom_prop,na.rm=T),mean(biom_prop,na.rm=T),max(biom_prop,na.rm=T)),digits=2)),bty="n",col=c("blue"))
	}
}
dev.off()

mm=c()
for (g in groupe){
        if(g=="BZ"){
                option_lieu=c("Men er Roue","Loscolo","Croisic")
        }else if(g=="MO"){
                option_lieu=c("Cornard","Auger","LEperon")
        }else if(g=="SU"){
                option_lieu=c("Antoine","Lazaret")
        }
	list_diat=c("AST","CHA","DIT","GUI","LEP","RHI","SKE","THP","NIT","PLE","PSE","THL")
	list_dino=c("GYM","PRO","PRP","SCR")
	
        a=table(as.matrix(tab_sp[,lieu %in% option_lieu]))
        liste_sp=dimnames(a[a==length(option_lieu)])[[1]]
        if("NEI"%in%liste_sp){
                liste_sp=liste_sp[-which(liste_sp=="NEI")]
        }

        for (ll in 1:length(option_lieu)){
                tab=read.table(paste("data/corres_hernandez_",option_lieu[ll],'.txt',sep=''),sep=";",na="NA",header=TRUE)
                dd=as.Date(tab$Date)
                tab=tab[year(dd)>=1996,]#Using data from 1996
                dd=dd[year(dd)>=1996]

                biom_diat=mean(apply(tab[,list_diat],1,mean,na.rm=T),na.rm=T)
                biom_dino=mean(apply(tab[,list_dino],1,mean,na.rm=T),na.rm=T)
		mm=c(mm,biom_diat/biom_dino)
		print(option_lieu[ll])
		print(biom_diat)
		print(biom_dino)
		print(biom_diat/biom_dino)
	}
}
