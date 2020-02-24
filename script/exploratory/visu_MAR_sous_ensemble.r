graphics.off()
rm(list=ls())
library(stringr)

#Graphical parameters
pm=0.05
ab=0.1
fact=1.5
acex=3.5

fac_main=3.0
fac_axis=1.34
fac_lg=1.5
fac_lab=1.25
alwd=2.5
apc=2

corres=read.table(paste("corres_hernandez.csv",sep=''),sep=";",na="NA",header=TRUE)

option_model=c("null","unconstrained","pencen","diatdin","inter")
option_NEI=c("null")
groupe=c("BZ","MO","SU","AR")
option_sp="common"
liste_reduced=c("CHA","SKE","PSE","PRP")
        #ON range les espèces en les regroupant
        delim=c()
        t1=as.character(corres$Type[corres$Code==liste_reduced[1]])
	for (s in 2:length(liste_reduced)){
                t2=as.character(corres$Type[corres$Code==liste_reduced[s]])
                if(t1!=t2){
			delim=c(delim,s-0.5)
                }
		t1=t2
        }

lieu_tot=c("Men er Roue","Loscolo","Croisic","LEperon","Cornard","Auger","Antoine","Lazaret","Teychan","B7")
#option_lieu=c("Men er Roue","Loscolo","Croisic","LEperon","Cornard","Auger","Antoine","Lazaret")
for(ne in 1:length(option_NEI)){
        for (m in 1:length(option_model)){
		id_lieu=0
			pdf(paste("Rapport/graphe/MAR_estimates/COMMON/with_AR/analyse_MAR",option_model[m],"_common_regular_with_bootstrap_and_signif.pdf",sep=""),width=13)
			par(mar=c(5,6,3,1))
			layout(matrix(c(1,1,1,2),1,4,byrow=TRUE))
			color=c("black","red","blue","green","violet","cyan","orange","white","darkorchid","salmon")
for (g in groupe){
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
			f1=paste("data/analyse_MAR/",g,"/common/",option_lieu[ll],"_",option_model[m],"_",option_NEI[ne],"_regular_reduced_ALL.RData",sep="")
			#f1=paste(option_lieu[ll],"_",option_model[m],"_",option_NEI[ne],"_regular_",option_sp,"_",g,".RData",sep="")
			load(f1)

			#cis=fit_log
			sp=dimnames(cis$model$data)[[1]] #Studied species
			var=dimnames(cis$call$model$c)[[1]] #Studied covariates
			nom=dimnames(cis$par$B)[[1]] #Names of the interactions

			#Draw the graph frame
			nb_pos_inter=0
			nb_neg_inter=0
			nb_pos_cov=0
			nb_neg_cov=0
			nb_par_tot_inter=0
			nb_par_tot_cov=0
			if(id_lieu==1){
				consensus_inter=matrix(FALSE,nrow=length(liste_reduced)*(length(liste_reduced)-1),ncol=length(lieu_tot))

				plot(0,0,t='n',xlim=c(0.75,(length(var)+length(liste_reduced)+1)),ylim=c(0.5,length(liste_reduced)+0.4),yaxt="n",xaxt="n",xlab="",ylab="",main=paste(option_model[m],sep=""),cex.main=2)
				axis(2,at=1:length(liste_reduced),lab=rev(liste_reduced),cex.axis=fac_axis)
				axis(1,at=1:(length(liste_reduced)+length(var)),lab=c(liste_reduced,var),cex.axis=fac_axis)
				for (i in 1:length(liste_reduced)){
        				abline(h=i,lty=2)
				}
				abline(v=length(liste_reduced)+1,lty=2)
				abline(v=delim+0.1,lty=4,col="grey",lwd=4)
				abline(h=length(liste_reduced)-delim+1,lty=4,col="grey",lwd=4)
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
        			} #end options on different interaction names

				if(sp[i]%in%liste_reduced&&sp[j]%in%liste_reduced){
					ibis=grep(sp[i],liste_reduced)
					jbis=grep(sp[j],liste_reduced)
			        	baseline=length(liste_reduced)-ibis+1 #Compute the position of the species on the y-axis
				        if(i==j){ #For intragroup interactions, withdraw 1 to the value, to make them comparable to intergroup interactions
        	        			mini=min(0,(cis$par$B[n]-1)*fact)
			                	maxi=max(0,(cis$par$B[n]-1)*fact)
						colo=color[id_lieu]
						xshift=id_lieu*pm+0.025
                				rect(jbis-pm+xshift,baseline+mini,jbis+pm+xshift,baseline+maxi,col=colo) #Draw the bar corresponding to the coefficient value
				        }else{	
			                	mini=min(0,(cis$par$B[n])*fact)
	                			maxi=max(0,(cis$par$B[n])*fact)
						nb_pos_inter=nb_pos_inter+as.numeric(cis$par$B[n]>0)
       	        				nb_neg_inter=nb_neg_inter+as.numeric(cis$par$B[n]<0)
			        	        nb_par_tot_inter=nb_par_tot_inter+as.numeric(cis$par$B[n]!=0)

						consensus_inter[nb_par_tot_inter,id_lieu]=cis$par$B[n]>0
						colo=color[id_lieu]
						xshift=id_lieu*pm+0.025
			                	rect(jbis-pm+xshift,baseline+mini,jbis+pm+xshift,baseline+maxi,col=colo)
        				}
				        if(!is.na(cis$par.se$B[n])){
        	        			if(cis$par.upCI$B[n]*cis$par.lowCI$B[n]>0){ #If upper and lower values of the confidence intervals have the same sign, the coefficient is deemed significant
			                	        points(jbis+xshift,baseline+maxi+ab,pch='*',cex=acex,col=colo)
                				}
				        }
				 } #end test on being in liste_reduced
			}#end loops on names

			#Covariate effects
			l=0
			for (j in 1:length(var)){
        			for (i in 1:length(sp)){
					if(sp[i]%in%liste_reduced){
						iplou=grep(sp[i],liste_reduced)
                				l=l+1
				                ibis=(iplou-1)*2+1 #Position on the x-axis
                				baseline=length(liste_reduced)-iplou+1 #Position on the y-axis
				                mini=min(0,(cis$par$U[l])*fact)
                				maxi=max(0,(cis$par$U[l])*fact)
		
						colo=color[id_lieu]
						xshift=id_lieu*pm
                				rect(length(liste_reduced)+j-pm+xshift,baseline+mini,length(liste_reduced)+j+pm+xshift,baseline+maxi,col=colo) #Draw a line
		                		if(!is.na(cis$par.se$U[l])){
		        		                if(cis$par.upCI$U[l]*cis$par.lowCI$U[l]>0){
                				                points(j+length(liste_reduced)+xshift,baseline+maxi+ab,pch='*',col=colo,cex=acex)
                        				}
	                			}#end test on SE
					} #end test on liste_reduced
        			} #end liste sp
			}#end liste var
		} #end ll
	} #end groupe
	if(nb_par_tot_inter>0){
	consensus_inter=consensus_inter[1:nb_par_tot_inter,]
	a=rowSums(consensus_inter)
	val_pos=sum(a==length(lieu_tot))
	val_neg=sum(a==0)
	mtext(paste(format(val_pos/nb_par_tot_inter,digits=2,nsmall=2,trim=TRUE),'com. inter. >0 /',format(val_neg/nb_par_tot_inter,digits=2,nsmall=2,trim=TRUE),'com. inter. <0',sep=" "),side=1,line=-2,outer=TRUE,cex=2)
	}
	plot(0,0,t="n",bty="n",xaxt="n",yaxt="n",ylab="",xlab="")
	legend("topleft",lieu_tot,fill=color,bty="n",cex=3)
	dev.off()
} #end model
} #end NEI
