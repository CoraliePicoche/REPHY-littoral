#2019/10/19 CP draw the same figure as in biotic_interaction but with segments corresponding to the confidence interval instead of bars

graphics.off()
rm(list=ls())
library(MARSS)
library(stringr)

#Graphical parameters
pm=0.1
ab=0.175
fact=1.75
acex=3.5

fac_axis=2.3
alwd=2.5
apc=2

corres=read.table(paste("corres_hernandez.csv",sep=''),sep=";",na="NA",header=TRUE)

option_model=c("pencen")
#option_model=c("unconstrained")
option_NEI=c("null")
groupe=c("BZ","MO","AR","SU")
llieux=c()
groupe_nice=c("Brittany","Oléron","Arcachon","Mediterranean")
let=c("a)","b)","d)","e)")
option_sp="common"
max_sp=16

#pdf("article/graphe/biotic_interaction_matrices_MainFig_allin1_v4.pdf",height=24,width=18)
pdf("article/submit_JEcol/response/biotic_interaction_matrices_pencen_CI.pdf",height=24,width=24)
par(mar=c(7.5,6.,5,3),oma=c(0.1,2.5,.5,0.5),mfrow=c(2,2))
color=c("black","blue","red")
colo_end=c()
for (gg in 1:length(groupe)){
	g=groupe[gg]
        if(g=="BZ"){
                option_lieu=c("Men er Roue","Loscolo","Croisic")
                colo_end=c(colo_end,rep("green",length(option_lieu)))
        }else if(g=="MO"){
                option_lieu=c("LEperon","Cornard","Auger")
                colo_end=c(colo_end,rep("darkblue",length(option_lieu)))
        }else if(g=="SU"){
                option_lieu=c("Antoine","Lazaret")
                colo_end=c(colo_end,rep("darkred",length(option_lieu)))
        }else if(g=="AR"){
                option_lieu=c("Teychan","B7")
                colo_end=c(colo_end,rep("cyan",length(option_lieu)))
        }

                for (ll in 1:length(option_lieu)){
                        f1=paste("data/analyse_MAR/",g,"/site_specific/",option_lieu[ll],"_",option_model,"_",option_NEI,"_regular_common_",g,".RData",sep="")
                        load(f1)
			if(ll==1&gg>1){
				llieux=c(llieux,"")
			}
			llieux=c(llieux,option_lieu[ll])
                        sp=dimnames(cis$model$data)[[1]] #Studied species
                        var=dimnames(cis$call$model$c)[[1]] #Studied covariates
                        nom=dimnames(cis$par$B)[[1]] #Names of the interactions

        #ON range les espèces en les regroupant
        delim=c()
        t1=as.character(corres$Type[corres$Code==sp[1]])
        for (s in 2:length(sp)){
                t2=as.character(corres$Type[corres$Code==sp[s]])
                if(t1!=t2){
                        delim=c(delim,s)
                }
                t1=t2
        }
                        #Draw the graph frame
		if((ll==1)||(ll>1&&option_sp!="common")){
                        plot(0,0,t='n',xlim=c(0.95,(max_sp+0.25)),ylim=c(0.7,max_sp),yaxt="n",xaxt="n",xlab="",ylab="",main=groupe_nice[gg],cex.main=3.)
			grad_sp=seq(1,max_sp,length.out=length(sp))
			if(gg==1){
			segments(1,1,1.,1+fact,col="blue",lwd=5)
			text(1+5*pm,1+fact,"1",cex=3.0)
			text(1+5*pm,1,"0",cex=3.0)
			text(1.+7*pm,1+fact*.5,"0.5",cex=3.0)
			text(1.+17*pm,1+fact*.45,expression("b"["ij"]),cex=3.5)#,srt=90)
			}
			mtext(let[gg],side=3,cex=2.0,xpd=NA,font=2,line=2.0,adj=0)
	                axis(2,at=grad_sp,lab=rev(sp),cex.axis=fac_axis,las=2)
                        axis(1,at=grad_sp,lab=c(sp),cex.axis=fac_axis,las=2)
                        for (i in 1:length(sp)){
				if(grad_sp[i]<3){
					if(gg==1){
						lines(c(3,max_sp+1),rep(grad_sp[i],2),lty=2)
					}else{
                              			abline(h=grad_sp[i],lty=2)
					}
				}else{
				if(grad_sp[i]<12){
                                abline(h=grad_sp[i],lty=2)
				}else{
				lines(c(0,10.25),rep(grad_sp[i],2),lty=2)
				}
				}
                        }
        		if(option_model!="unconstrained"){
			legend(x=9.8,y=max_sp+1-0.75,option_lieu,fill=color[1:3],cex=3.5,border=NA,xpd=NA,bty="n")
			}
        		#legend("topright",option_lieu,fill=color[1:3],cex=3.5,border=NA,bg="white")#xpd=NA)#bty="n")
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
        baseline=max_sp-grad_sp[i]+1 #Compute the position of the species on the y-axis
        if(i==j){ #For intragroup interactions, withdraw 1 to the value, to make them comparable to intergroup interactions
                mini=cis$par.lowCI$B[n]*fact
                maxi=cis$par.upCI$B[n]*fact
                if(option_sp=="common"){
                        colo=color[ll]
                        xshift=(ll-1)*2.75*pm
                }else{
                xshift=0
                if(mini==0){
                        colo="blue"
                }else{
                        colo="red"
                }
                }
                arrows(grad_sp[j]+xshift,baseline+mini,grad_sp[j]+xshift,baseline+maxi,col=colo,lwd=5,angle = 90,code = 3, length = 0.05) #Draw the bar corresponding to the coefficient value
        }else{
                mini=cis$par.lowCI$B[n]*fact
                maxi=cis$par.upCI$B[n]*fact


                if(option_sp=="common"){
                        colo=color[ll]
                        xshift=(ll-1)*2.75*pm
                }else{
                        xshift=0
                        if(mini==0){
                                colo="blue"
                        }else{
                                colo="red"
                        }
                }
                arrows(grad_sp[j]+xshift,baseline+mini,grad_sp[j]+xshift,baseline+maxi,col=colo,lwd=5,angle = 90,code = 3, length = 0.05)
        }
         }

}
}


dev.off()
