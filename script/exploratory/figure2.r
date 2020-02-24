#2018/09/05 CP
#Figure 2 showing the coefficients of B and C for the 10 subsites we study

graphics.off()
rm(list=ls())
library(MARSS)
library(stringr)

#Graphical parameters
pm=0.1
ab=0.1
fact=1.5
acex=3.5

fac_axis=2.3
alwd=2.5
apc=2

corres=read.table(paste("corres_hernandez.csv",sep=''),sep=";",na="NA",header=TRUE)

option_model=c("pencen")
option_NEI=c("null")
groupe=c("BZ","MO","AR","SU")
groupe_nice=c("Brittany","Oléron","Arcachon","Mediterranean")
let=c("a)","b)","c)","d)")
option_sp="common"

pdf("Rapport/graphe/MAR_estimates/figure2.pdf",height=26,width=26.75)
par(mar=c(6.5,4,5,.5))
layout(matrix(c(1,1,1,2,3,3,3,4,5,5,5,6,7,7,7,8),2,8,byrow=TRUE))
color=c("black","blue","red")
for (gg in 1:length(groupe)){
	g=groupe[gg]
        if(g=="BZ"){
                option_lieu=c("Men er Roue","Loscolo","Croisic")
        }else if(g=="MO"){
                option_lieu=c("LEperon","Cornard","Auger")
        }else if(g=="SU"){
                option_lieu=c("Antoine","Lazaret")
        }else if(g=="AR"){
                option_lieu=c("Teychan","B7")
        }

                for (ll in 1:length(option_lieu)){
                        f1=paste("data/analyse_MAR/",g,"/site_specific/",option_lieu[ll],"_",option_model,"_",option_NEI,"_regular_common_",g,".RData",sep="")
                        load(f1)

                        sp=dimnames(cis$model$data)[[1]] #Studied species
                        var=dimnames(cis$call$model$c)[[1]] #Studied covariates
                        nom=dimnames(cis$par$B)[[1]] #Names of the interactions

        #ON range les espèces en les regroupant
        delim=c()
        t1=as.character(corres$Type[corres$Code==sp[1]])
        for (s in 2:length(sp)){
                t2=as.character(corres$Type[corres$Code==sp[s]])
                if(t1!=t2){
                        delim=c(delim,s-0.5)
                }
                t1=t2
        }
                        #Draw the graph frame
                        nb_pos_inter=0
                        nb_neg_inter=0
                        nb_pos_cov=0
                        nb_neg_cov=0
                        nb_par_tot_inter=0
                        nb_par_tot_cov=0
                        if((ll==1)||(ll>1&&option_sp!="common")){
                        consensus_inter=matrix(FALSE,nrow=length(nom)-length(sp),ncol=length(option_lieu))

                        plot(0,0,t='n',xlim=c(0.75,(length(var)+length(sp)+0.5)),ylim=c(0.5,length(sp)+0.25),yaxt="n",xaxt="n",xlab="",ylab="",main=groupe_nice[gg],cex.main=3.5)
               		mtext(let[gg],side=3,cex=2.,xpd=NA,font=2,line=1,adj=0)
	                axis(2,at=1:length(sp),lab=rev(sp),cex.axis=fac_axis)
                        axis(1,at=1:(length(sp)+length(var)),lab=c(sp,var),cex.axis=fac_axis,las=2)
                        for (i in 1:length(sp)){
                                abline(h=i,lty=2)
                        }
                        abline(v=length(sp)+0.75,lty=2)
                        abline(v=delim+0.05,lty=4,col="grey",lwd=4)
                        abline(h=length(sp)-delim+1,lty=4,col="grey",lwd=4)
                        }

                        #get the value from the MARSS object, according to the names of the coefficients
                        for (n in 1:length(nom)){
                                if(grepl("^\\(",nom[n])){ #default when using unconstrained or diagonal, names are of the form (1,2) for (i,j)
                                        i=as.numeric(strsplit(nom[n],split="[\\(,\\)]")[[1]][2])
                                        j=as.numeric(strsplit(nom[n],split="[\\(,\\)]")[[1]][3])
                if(is.na(i)){ #i and j might not be numeric
                        i=which(unlist(lapply(sp,function(x) grepl(x,strsplit(nom[n],split="[\\(,\\)]")[[1]][2]))))
                        j=which(unlist(lapply(sp,function(x) grepl(x,strsplit(nom[n],split="[\\(,\\)]")[[1]][3]))))
                }
        }else{
                a=str_split(nom[n],sp) #The coefficient is XY where X eats Y, X and Y being in the species list
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
        baseline=length(sp)-i+1 #Compute the position of the species on the y-axis
        if(i==j){ #For intragroup interactions, withdraw 1 to the value, to make them comparable to intergroup interactions
                mini=min(0,(cis$par$B[n]-1)*fact)
                maxi=max(0,(cis$par$B[n]-1)*fact)
                if(option_sp=="common"){
                        colo=color[ll]
                        xshift=ll*pm+0.05
                }else{
                xshift=0
                if(mini==0){
                        colo="blue"
                }else{
                        colo="red"
                }
                }
                rect(j-pm+xshift,baseline+mini,j+pm+xshift,baseline+maxi,col=colo) #Draw the bar corresponding to the coefficient value
        }else{
                mini=min(0,(cis$par$B[n])*fact)
                maxi=max(0,(cis$par$B[n])*fact)


                nb_pos_inter=nb_pos_inter+as.numeric(cis$par$B[n]>0)
                nb_neg_inter=nb_neg_inter+as.numeric(cis$par$B[n]<0)
                nb_par_tot_inter=nb_par_tot_inter+as.numeric(cis$par$B[n]!=0)

                consensus_inter[nb_par_tot_inter,ll]=cis$par$B[n]>0

                if(option_sp=="common"){
                        colo=color[ll]
                        xshift=ll*pm+0.05
                }else{
                        xshift=0
                        if(mini==0){
                                colo="blue"
                        }else{
                                colo="red"
                        }
                }
                rect(j-pm+xshift,baseline+mini,j+pm+xshift,baseline+maxi,col=colo)
        }
        if(!is.na(cis$par.se$B[n])){
                if(cis$par.upCI$B[n]*cis$par.lowCI$B[n]>0){ #If upper and lower values of the confidence intervals have the same sign, the coefficient is deemed significant
                        points(j+xshift,baseline+maxi+ab,pch='*',cex=acex,col=colo)
                }
        }
         }

#Covariate effects
l=0
for (j in 1:length(var)){
        for (i in 1:length(sp)){
                l=l+1
                ibis=(i-1)*2+1 #Position on the x-axis
                baseline=length(sp)-i+1 #Position on the y-axis
                mini=min(0,(cis$par$U[l])*fact)
                maxi=max(0,(cis$par$U[l])*fact)

                if(option_sp=="common"){
                        colo=color[ll]
                        xshift=ll*pm
                }else{
                        xshift=0
                if(mini==0){
                        colo="blue"
                }else{
                        colo="red"
                }
                }
                rect(length(sp)+j-pm+xshift,baseline+mini,length(sp)+j+pm+xshift,baseline+maxi,col=colo) #Draw a line
                if(!is.na(cis$par.se$U[l])){
                        if(cis$par.upCI$U[l]*cis$par.lowCI$U[l]>0){
                                points(j+length(sp)+xshift,baseline+maxi+ab,pch='*',col=colo,cex=acex)
                        }
                }
        }
}
}
        a=rowSums(consensus_inter)
        val_pos=sum(a==length(option_lieu))
        val_neg=sum(a==0)
#        mtext(paste(format(val_pos/nb_par_tot_inter,digits=2,nsmall=2,trim=TRUE),'com. inter. >0 /',format(val_neg/nb_par_tot_inter,digits=2,nsmall=2,trim=TRUE),'com. inter. <0',sep=" "),side=1,line=-2,outer=TRUE,cex=2)
        plot(0,0,t="n",bty="n",xaxt="n",yaxt="n",ylab="",xlab="")
        legend(-1.4,.5,c(option_lieu,paste(format(val_pos/nb_par_tot_inter,digits=2,nsmall=2,trim=TRUE),' +',sep=""),paste(format(val_neg/nb_par_tot_inter,digits=2,nsmall=2,trim=TRUE),' -',sep="")),fill=c(color[1:ll],NA,NA),bty="n",cex=3.25,border=NA)

}
dev.off()
