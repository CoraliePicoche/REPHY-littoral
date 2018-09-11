##CP 15/10/2017. Compute some network metrics on our networks. At first, I wanted a 'simple' weighted version of connectance, realized it did not exist. Using bipartite library which implements the Bersier et al. 2002 (slightly different from Altena et al. 2016) version of weighted linkage density. Still wondering if using a 'bipartite' library on a 'regular' network is a good idea, but could not find the same index in other libraries. 

rm(list=ls())
graphics.off()
library('bipartite')
source("script/matrix_MAR_clean.r")

groupe=c("BZ","MO","SU","AR")
option_model="unconstrained"
#metric='linkage density'
#metric='connectance'
metric='weighted connectance'
pdf(paste("./Rapport/graphe/weighted_connectance_littoral_",option_model,".pdf",sep=""))
#pdf(paste("./Rapport/graphe/linkage_density_littoral_",option_model,".pdf",sep=""))
par(mar=c(8,5,0.5,0.5))

#plot(0,0,t="n",xlim=c(1,11),ylim=c(2.0,8.0),xlab='',ylab=metric,xaxt="n",cex.lab=2,cex.axis=2) #Taking into account both inter and intragroup competition
#plot(0,0,t="n",xlim=c(1,11),ylim=c(0.35,0.7),xlab='',ylab=metric,xaxt="n") #Taking into account both inter and intragroup competition
plot(0,0,t="n",xlim=c(1,11),ylim=c(0,0.3),xlab='',ylab=metric,xaxt="n",cex.lab=2,cex.axis=2) #Taking into account both inter and intragroup competition
axis(1,at=seq(1,10),lab=c("Men er.","Loscolo","Croisic","LEperon","Cornard","Auger","Antoine","Lazaret","Teychan","B7"),las=2,cex.axis=2,cex.lab=2)
id_lieu=0
id_g=0
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

for (ll in 1:length(option_lieu)){
	id_lieu=id_lieu+1
#for (ll in 2:2){
        #Hydro variables : Seasonality
#	if(g=="AR"){
#	tab_cov=read.csv(paste("./data/",option_lieu[ll],"_base.csv",sep=""),na.strings="NA",header=TRUE,sep=";",dec=".")
#	}else{
#	tab_cov=read.table(paste("data/",option_lieu[ll],'hydro.txt',sep=''),sep=";",na="NA",header=TRUE)
#	}

        f1=paste("data/analyse_MAR/",g,"/site_specific/",option_lieu[ll],"_",option_model,"_null_regular_common_",g,".RData",sep="")
        load(f1)
#        nom=dimnames(fit_log$par$B)[[1]] #Names of the interactions
#	B=diag(-1,sqrt(length(nom)),sqrt(length(nom))) #I am heavily using the fact that we are using an unconstrained matrix
#	for (n in 1:length(nom)){
#        	if(grepl("^\\(",nom[n])){ #default when using unconstrained or diagonal, names are of the form (1,2) for (i,j)
#                	i=as.numeric(strsplit(nom[n],split="[\\(,\\)]")[[1]][2])
#                        j=as.numeric(strsplit(nom[n],split="[\\(,\\)]")[[1]][3])
#			B[i,j]=B[i,j]+cis$par$B[n]
#		}else{
#			print(nom[n])
#		}
#	}

	B=clean_matrix(cis)
	prop_signif=sum(sign(cis$par.lowCI$B*cis$par.upCI$B)>0)/sum(B!=0) #Varies between 0.2 and 0.4 for unconstrained ; 0.13 and 0.21 for 
	if (option_model=="unconstrained"){
		mini=0.21
		maxi=0.39
	}else if(option_model=="pencen"){
		mini=0.35
		maxi=0.59
	}
	size_point=1.5*(prop_signif-mini)/(maxi-mini)


	ll_abs=networklevel(web=abs(B),index=c("connectance","linkage density","weighted connectance")) #At first, I wanted all but 'quantitative' takes a VERY long time to compute
	print(ll_abs["connectance"])
	B_pos=matrix(0,nrow=dim(B)[1],ncol=dim(B)[2])
	B_pos[B>0]=B[B>0]
	ll_pos=networklevel(web=B_pos,index=c("connectance","linkage density","weighted connectance")) #At first, I wanted all but 'quantitative' takes a VERY long time to compute
	B_neg=matrix(0,nrow=dim(B)[1],ncol=dim(B)[2])
	B_neg[B<0]=B[B<0]
	ll_neg=networklevel(web=abs(B_neg),index=c("connectance","linkage density","weighted connectance")) #At first, I wanted all but 'quantitative' takes a VERY long time to compute

	points(id_lieu,ll_abs[metric],pch=16,col="black",cex=1+exp(size_point))
	points(id_lieu+0.1,ll_pos[metric],pch=16,col="blue",cex=1+exp(size_point))
	points(id_lieu+0.2,ll_neg[metric],pch=16,col="red",cex=1+exp(size_point))
	abline(v=c(3.6,6.6,8.6),lty=2,col="grey")

	if(1==0){
	B_abs=abs(B)
	diag(B_abs)=0
	points(id_lieu,mean(abs(B)),pch=16,col="black",cex=1+exp(size_point))
	diag(B_pos)=0
	points(id_lieu+0.1,mean(B_pos),pch=16,col="blue",cex=1+exp(size_point))
	diag(B_neg)=0
	points(id_lieu+0.2,mean(abs(B_neg)),pch=16,col="red",cex=1+exp(size_point))
	}

#	points(season[id_lieu],median(coeff),t="p",pch=16,col=colog[id_g])
#	segments(season[id_lieu],median(coeff)-sd(coeff) ,season[id_lieu],median(coeff)+sd(coeff),col=colog[id_g])
#	print(mean(coeff))
#	abline(h=0,lty=2,col="grey")

}

}
legend("bottomleft",c("All inter","Positive inter","Negative inter"),col=c("black","blue","red"),pch=16,bty="n")
#legend("bottomleft",c("Positive inter","Negative inter"),col=c("blue","red"),pch=16,bty="n",cex=2)
dev.off()	
