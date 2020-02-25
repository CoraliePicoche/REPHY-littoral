#CP 12/2017: this script extracts intra/inter ratios from other MAR analyses and other metrics on the matrices and compare them to our results on the REPHY dataset
#CP 10/2019: identify terrestrial/predation food webs with a different color
#CP 01/2020: removed Barraquand 2018 as it was redundant with our current data + changed paths to have figures in the response folder + added computation of the mean inter on the same number as intra
#CP 02/2020: cleaned up the code

graphics.off()
rm(list=ls())

source("./script/matrix_MAR_clean.r")

#This function extracts the ratios, absolute values and standard deviations of intraspecific and interspecific interaction strengths
whatwewant=function(x_intra,x_inter){
	dimension=length(x_intra)
	x_inter=x_inter[x_inter!=0]
	iter=100
	mean_ratio_with_sample=0
	for(i in 1:iter){
		random_inter=sample(x_inter,dimension,replace=T)
		mean_ratio_with_sample=mean_ratio_with_sample+(1/iter)*mean(abs(x_intra))/mean(abs(random_inter))
	}
	prop_signif=sum(c(x_inter,x_intra)!=0)/(dimension^2)
	val_signif=c(x_inter,rep(0,dimension*(dimension-1)-length(x_inter)))
	val=c(mean(abs(x_intra)),sd(abs(x_intra)),mean(abs(x_inter)),sd(abs(x_inter)),mean(abs(val_signif)),sd(abs(val_signif)),dimension,prop_signif,mean_ratio_with_sample)
}

#We begin with the data from other papers
tab_answer=matrix(NA,nrow=1,ncol=10)
colnames(tab_answer)=c("Code","MeanIntra","SdIntra","MeanInter","SdInter","MeanSignif","SdSignif","Dimension","Prop signif","SampledRatio")
code_name=c('Hampton2006a_biweekly','Hampton2006a_growing','Hampton2006b_full','Hampton2006b_simple','Griffiths2015_Phyto1','Griffiths2015_Phyto2','Huber2006_Phyto_signif','Huber2006_Ciliate_signif','Ives1999_Zoo1','Ives1999_Zoo2','Ives2003_Plank1','Ives2003_Plank2_signif','Ives2003_Plank3_signif','Klug2000_Phyto','Klug2000_Zoo','Klug2001_TaxoAlgae','Klug2001_MorphoAlgae','Lindegren2009_Fish_signif','Vik2008_LynxHare','Yamamura2006_Insects_signif')

for(c in 1:length(code_name)){
plou=as.matrix(read.csv(paste(code_name[c],'.csv',sep=""),sep=";",dec=".",header=T,row.names=1))
if(diag(plou)[1]>0){
	x_intra=diag(plou)-1
}else{
	x_intra=diag(plou)
}
nodiag=plou
diag(nodiag)=NA
nodiag=c(nodiag)
nodiag=nodiag[!is.na(nodiag)]
x_inter=nodiag
tab_answer=rbind(tab_answer,c(code_name[c],whatwewant(x_intra,x_inter)))
}

tab_answer=tab_answer[-which(is.na(tab_answer[,1])),]

#Alphabetical order
id=order(code_name)
tab_answer=tab_answer[id,]
code_name=code_name[id]

#Chronological order
years=unlist(regmatches(code_name, gregexpr("[[:digit:]]{4}", code_name)))
id=order(years)
tab_answer=tab_answer[id,]
code_name=code_name[id]

#This replaces names of the studies by their codes (1a, 2, etc.)
xlab_try=rep(NA,length(code_name))
names_1=rep(NA,length(code_name))
aletter=c('a','b','c')
xlab_try[1]=1
names_1[1]="1a"
al=1
for(i in 2:length(code_name)){
	beg1=strsplit(code_name[i-1],"_")[[1]][1]
	beg2=strsplit(code_name[i],"_")[[1]][1]
	if(beg1==beg2){
		al=al+1
		xlab_try[i]=xlab_try[i-1]	
	}else{
		al=1
		xlab_try[i]=xlab_try[i-1]+1
	}
	names_1[i]=paste(as.character(xlab_try[i]),aletter[al],sep="")
}

xlab_try=c(xlab_try,rep(xlab_try[length(xlab_try)]+1,10))

#This uses values from the present analysis on the REPHY Dataset
groupe=c("BZ","MO","AR","SU")
other=c('Br','Ol','Ar','Me')
app_T_this_study=c()
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
		c=c+1
		names_1=c(names_1,paste(other[gg],as.character(ll),sep=""))
		code_name=c(code_name,paste('ThisStudy','_',option_lieu[ll],sep=""))
        	f1=paste("~/Documents/Plankton/REPHY_littoral/data/analyse_MAR/",g,"/site_specific/",option_lieu[ll],"_pencen_null_regular_common_",g,".RData",sep="")
                load(f1)
		app_T_this_study=c(app_T_this_study,round(dim(fit_log$model$data)[2]/100)*100) #small trick so that all lengths are comparable
		plou=clean_matrix(cis,signif=T)
		x_intra=diag(plou)
		nodiag=plou
		diag(nodiag)=NA
		nodiag=c(nodiag)
		nodiag=nodiag[!is.na(nodiag)]
		x_inter=nodiag
		tab_answer=rbind(tab_answer,c(code_name[c],whatwewant(x_intra,x_inter)))
	}
}


tab_answer[,'Prop signif']=1-as.numeric(tab_answer[,'Prop signif']) #Sparsity index

#Comparing ratios when taking only the significant parameters. This produces Fig. S9 in the supplementary information
pdf("./comparaison_ratio_code_nolog_cleaner_ONLY_SIGNIF_withoutBarraquand_et_al_2018.pdf",width=18,height=8)
par(mfrow=c(1,2),xpd=NA,mar=c(4,4.5,3,0.5))
cod_cod=rep("cyan",dim(tab_answer)[1])
cod_cod[as.numeric(tab_answer[,"Prop signif"])>=0.65]="blue"
cod_cod[as.numeric(tab_answer[,"Prop signif"])>=0.75]="darkorchid"

symbol=rep(21,dim(tab_answer)[1])
symbol[as.numeric(tab_answer[,"Prop signif"])>=0.65]=22
symbol[as.numeric(tab_answer[,"Prop signif"])>=0.75]=23

ylab_try=as.numeric(tab_answer[,"MeanIntra"])/as.numeric(tab_answer[,"MeanInter"])

plot(tab_answer[,"Dimension"],ylab_try,t="p",pch=symbol,bg=cod_cod,col="black",cex=2,xlab="Number of taxa",ylab="|intra|/|inter|",cex.lab=2.0,cex.axis=1.5,xlim=c(1.5,14.5))
id_diff=16:18
points(tab_answer[id_diff,"Dimension"],ylab_try[id_diff],col="red",lwd=2,cex=1,pch=symbol[id_diff],bg="red")
mtext("a)",side=3,cex=1.5,xpd=NA,font=2,line=1.0,adj=0)
idx=c(1:3,5:29,22:30)
text(tab_answer[idx,'Dimension'],ylab_try[idx],names_1[idx],pos=c(4,2),cex=1.5)

idx=4
text(tab_answer[idx,'Dimension'],ylab_try[idx]+0.05,names_1[idx],pos=c(4),cex=1.5)

idx=20 
text(tab_answer[idx,'Dimension'],ylab_try[idx],names_1[idx],pos=c(4),cex=1.5)
#idx=23 #w Barraquand
idx=21 #wo Barraquand
text(tab_answer[idx,'Dimension'],ylab_try[idx]-0.15,names_1[idx],pos=c(4),cex=1.5)

legend("bottomright",pch=c(21,22,23),leg=c('sparsity<0.65',expression('0.65'<='sparsity<0.75'),expression('sparsity'>='0.75')),pt.bg=c("cyan","blue","darkorchid"),bty='n',cex=1.5,pt.cex=3)

#Second plot
xx=as.numeric(tab_answer[,"Prop signif"])
yy=as.numeric(tab_answer[,"MeanIntra"])/as.numeric(tab_answer[,"MeanInter"])
plot(xx,yy,pch=symbol,ylab="",xlab="",cex=2,cex.axis=1.5,xlim=c(-0.0,1),bg=cod_cod)
id_diff=16:18
points(xx[id_diff],yy[id_diff],col="red",lwd=2,cex=1,pch=symbol[id_diff],bg="red")
mtext("b)",side=3,cex=1.5,xpd=NA,font=2,line=1.0,adj=0)

names_1_bis=names_1
names_1_bis[c(2,8,19,26,22,23,29)]=""
text(xx,yy,names_1_bis,pos=c(4,2),cex=1.5)
idx=c(2,8)
text(xx[idx],yy[idx]-0.1,names_1[idx],pos=1,cex=1.5)

text(xx[19],yy[19],names_1[19],pos=2,cex=1.5)

text(xx[22]+0.01,yy[22]+0.1,names_1[22],pos=4,cex=1.5)
arrows(xx[22]+0.02,yy[22]+0.1,xx[22],yy[22],length=0)

text(xx[23]-0.02,yy[23]+0.25,names_1[23],pos=2,cex=1.5)
arrows(xx[23]-0.03,yy[23]+0.25,xx[23],yy[23],length=0)

text(xx[26]+0.04,yy[26]+0.175,names_1[26],pos=4,cex=1.5)
arrows(xx[26]+0.05,yy[26]+0.175,xx[26],yy[26],length=0)

text(xx[29]-0.04,yy[29]-0.025,names_1[29],pos=2,cex=1.5)
arrows(xx[29]-0.05,yy[29]-0.025,xx[29],yy[29],length=0)

mtext("Sparsity",side=1,line=2.5,cex=2.0)

dev.off()

#Comparing ratios using only significant parameters but replacing missing values by 0. This is not used anymore.
par(mfrow=c(1,1),xpd=NA,mar=c(3,4.5,1,1))

symbol=rep(21,dim(tab_answer)[1])
symbol[as.numeric(tab_answer[,"Prop signif"])>=0.65]=22
symbol[as.numeric(tab_answer[,"Prop signif"])>=0.75]=23


plot(xlab_try,log10(as.numeric(tab_answer[,"MeanIntra"])/as.numeric(tab_answer[,"MeanSignif"])),t="p",pch=symbol,bg=cod_cod,col="black",cex=2,xlab="",ylab="|intra|/|inter|",yaxt="n",ylim=c(-0.1,2),cex.lab=1.5,xaxt="n",xlim=c(0.8,max(xlab_try)+0.5))
lines(c(0.7,max(xlab_try)+0.6),rep(1,2),col="black",lty=2,lwd=2)
mtext("Chronological order",side=1,line=1.5,cex=1.5)
axis(2,at=c(0,log10(5),1,log10(50),2),lab=c("1","5","10","50","100"),cex.axis=1.5)
id=c(1,2,4:12,14:20)
text(xlab_try[id],log10(as.numeric(tab_answer[id,"MeanIntra"])/as.numeric(tab_answer[id,"MeanSignif"])),names_1[id],pos=c(2,4),cex=1.5)
id=3
text(xlab_try[id]-0.1,log10(as.numeric(tab_answer[id,"MeanIntra"])/as.numeric(tab_answer[id,"MeanSignif"])),names_1[id],pos=c(2,4),cex=1.5)
id=13
text(xlab_try[id],0.05+log10(as.numeric(tab_answer[id,"MeanIntra"])/as.numeric(tab_answer[id,"MeanSignif"])),names_1[id],pos=c(2,4),cex=1.5)

text(xlab_try[29:30],log10(as.numeric(tab_answer[29:30,"MeanIntra"])/as.numeric(tab_answer[29:30,"MeanSignif"])),names_1[29:30],pos=c(2,4),cex=1.5)
idx=c(25)
text(xlab_try[idx],log10(as.numeric(tab_answer[idx,"MeanIntra"])/as.numeric(tab_answer[idx,"MeanSignif"])),names_1[idx],pos=c(2,4),cex=1.5)
idx=28
text(xlab_try[idx],log10(as.numeric(tab_answer[idx,"MeanIntra"])/as.numeric(tab_answer[idx,"MeanSignif"])),names_1[idx],pos=c(4),cex=1.5)
id=order(log10(as.numeric(tab_answer[22:26,"MeanIntra"])/as.numeric(tab_answer[22:26,"MeanSignif"])),decreasing=T)
id=id+21
for(i in 1:length(id)){
	ylab_try=log10(as.numeric(tab_answer[id[i],"MeanIntra"])/as.numeric(tab_answer[id[i],"MeanSignif"]))
	text(xlab_try[id[i]]-0.2,ylab_try-i*0.05,names_1[id[i]],pos=2,cex=1.5)
	arrows(xlab_try[id[i]]-0.25,ylab_try-i*0.05,xlab_try[id[i]],ylab_try,length=0)
}

#This compares the ratios of intra and inter against the sparsity and length of the time-series of each studies (Fig. S10 in the supporting information)
pdf("./sparsity_vs_others_species2taxa_withoutBarraquand_et_al_2018.pdf",width=18,height=8)
par(mfrow=c(1,2),xpd=NA,mar=c(4,4.5,3,0.5))

app_T=c(100,100,100,50,300,300,100,100,100,300,200,400,400,300,300,50,100,30,1000,700)
col=rep('black',length(app_T))
sy=rep(16,length(app_T))
app_T=c(app_T,app_T_this_study)
col=c(col,rep('blue',length(app_T_this_study)))
sy=c(sy,rep(17,length(app_T_this_study)))

par(mfrow=c(1,2))
plot(app_T,tab_answer[,'Prop signif'],col=col,pch=sy,cex=2,cex.lab=2.0,cex.axis=1.5,ylab="Sparsity",xlab="Length of the time series")
id_diff=16:18
points(app_T[id_diff],tab_answer[id_diff,"Prop signif"],col="red",lwd=2,cex=1,pch=symbol[id_diff],bg="red")
mtext("a)",side=3,cex=1.5,xpd=NA,font=2,line=1.0,adj=0)
plot(tab_answer[,'Dimension'],tab_answer[,'Prop signif'],col=col,pch=sy,cex.lab=2.0,cex=2,cex.axis=1.5,ylab="",xlab="Number of taxa")
points(tab_answer[id_diff,'Dimension'],tab_answer[id_diff,"Prop signif"],col="red",lwd=2,cex=1,pch=symbol[id_diff],bg="red")
mtext("b)",side=3,cex=1.5,xpd=NA,font=2,line=1.0,adj=0)
legend("bottomright",c("Other studies","This study"),pch=c(16,17),col=c("black","blue"),bty="n",cex=2)

dev.off()

symbol=rep(21,dim(tab_answer)[1])
symbol[as.numeric(tab_answer[,"Prop signif"])>=0.65]=22
symbol[as.numeric(tab_answer[,"Prop signif"])>=0.75]=23

#This shows the ratio intra to inter against the number of taxa in the study (Fig. 4 in the main text)
pdf('~/Documents/Plankton/REPHY_littoral/article/submit_JEcol/response_R2/Ratio_function_dim_tmp_withoutBarraquand_et_al_2018.pdf',,width=10,height=8)
par(mfrow=c(1,1),mar=c(4,4.5,1,1))
plot(tab_answer[,"Dimension"],log10(as.numeric(tab_answer[,"MeanIntra"])/as.numeric(tab_answer[,"MeanSignif"])),t="p",pch=symbol,bg=cod_cod,col="black",cex=2,xlab="Number of taxa",ylab="|intra|/|inter|",yaxt="n",ylim=c(-0.1,2.1),cex.lab=1.5,xlim=c(1.9,14.5),cex.axis=1.5)
axis(2,at=c(0,log10(5),1,log10(50),2),lab=c("1","5","10","50","100"),cex.axis=1.5)
id_diff=16:18
points(tab_answer[id_diff,"Dimension"],log10(as.numeric(tab_answer[id_diff,"MeanIntra"])/as.numeric(tab_answer[id_diff,"MeanSignif"])),col="red",lwd=2,cex=1,pch=symbol[id_diff],bg="red")


idx=1:30
tmp=names_1[18]
names_1[18]=NA
text(tab_answer[idx,'Dimension'],log10(as.numeric(tab_answer[idx,"MeanIntra"])/as.numeric(tab_answer[idx,"MeanSignif"])),names_1[idx],pos=c(4,2),cex=1.5)
idx=18
text(tab_answer[idx,'Dimension'],log10(as.numeric(tab_answer[idx,"MeanIntra"])/as.numeric(tab_answer[idx,"MeanSignif"])),tmp,pos=c(4),cex=1.5)
legend("bottomright",leg=c('sparsity<0.65',expression("0.65"<="sparsity<0.75"),expression('sparsity'>='0.75')),pt.bg=c("cyan","blue","darkorchid"),pch=c(21,22,23),bty='n',cex=1.5,lty=NA,lwd=1)
dev.off()

##Comparing mean b_ij with all coefficients and only the same number of coefficients as b_ii. This is not used anymore.
pdf("Comparison_mean.pdf")
par(mfrow=c(1,1))
x=as.numeric(tab_answer[,"MeanIntra"])/as.numeric(tab_answer[,"MeanInter"])
y=as.numeric(tab_answer[,"SampledRatio"])
plot(x,y,xlab="All bij",ylab="Sampled bij")
abline(a=0,b=1)
dev.off()
