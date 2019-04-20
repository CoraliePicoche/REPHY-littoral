graphics.off()
rm(list=ls())

source("~/Documents/Plankton/REPHY_littoral/script/matrix_MAR_clean.r")

whatwewant=function(x_intra,x_inter){
	dimension=length(x_intra)
	x_inter=x_inter[x_inter!=0]
	prop_signif=sum(c(x_inter,x_intra)!=0)/(dimension^2)
	val_signif=c(x_inter,rep(0,dimension*(dimension-1)-length(x_inter)))
	val=c(mean(abs(x_intra)),sd(abs(x_intra)),mean(abs(x_inter)),sd(abs(x_inter)),mean(abs(val_signif)),sd(abs(val_signif)),dimension,prop_signif)
}


tab_answer=matrix(NA,nrow=1,ncol=9)
colnames(tab_answer)=c("Code","MeanIntra","SdIntra","MeanInter","SdInter","MeanSignif","SdSignif","Dimension","Prop signif")
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

#Barraquand2018
code_name=c(code_name,"Barraquand2018_Teychan","Barraquand2018_B7")
c=c+1
load("~/Documents/Plankton/script_propre/Fred_master/PhytoplanktonArcachon_MultivariateTimeSeriesAnalysis/MARSS_results/Teychan_physics_pencen.RData")
plou=clean_matrix(cis,signif=TRUE)
x_intra=diag(plou)
nodiag=plou
diag(nodiag)=NA
nodiag=c(nodiag)
nodiag=nodiag[!is.na(nodiag)]
x_inter=nodiag
tab_answer=rbind(tab_answer,c(code_name[c],whatwewant(x_intra,x_inter)))

c=c+1
load("~/Documents/Plankton/script_propre/Fred_master/PhytoplanktonArcachon_MultivariateTimeSeriesAnalysis/MARSS_results/B7_physics_pencen.RData")
plou=clean_matrix(cis,signif=TRUE)
x_intra=diag(plou)
nodiag=plou
diag(nodiag)=NA
nodiag=c(nodiag)
nodiag=nodiag[!is.na(nodiag)]
x_inter=nodiag
tab_answer=rbind(tab_answer,c(code_name[c],whatwewant(x_intra,x_inter)))

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

#names_1=paste('[',1:length(code_name),']',sep='')

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

#names_1=code_name

xlab_try=c(xlab_try,rep(xlab_try[length(xlab_try)]+1,10))

#This Study
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
	#	names_1=c(names_1,paste(substr(g,1,1),as.character(ll),sep=""))
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

#Comparing
pdf("~/Documents/Plankton/REPHY_littoral/article/graphe/comparaison_ratio_code_nolog_cleaner_ONLY_SIGNIF.pdf",width=18,height=8)
par(mfrow=c(1,2),xpd=NA,mar=c(4,4.5,3,0.5))
cod_cod=rep("cyan",dim(tab_answer)[1])
cod_cod[as.numeric(tab_answer[,"Prop signif"])>=0.65]="blue"
cod_cod[as.numeric(tab_answer[,"Prop signif"])>=0.75]="darkorchid"

symbol=rep(21,dim(tab_answer)[1])

ylab_try=as.numeric(tab_answer[,"MeanIntra"])/as.numeric(tab_answer[,"MeanInter"])

plot(xlab_try,ylab_try,t="p",pch=symbol,bg=cod_cod,col="black",cex=2,xlab="",ylab="|intra|/|inter|",cex.lab=2.0,xlim=c(0.85,max(xlab_try)+0.8),xaxt='n',cex.axis=1.5)
mtext("Chronological order",side=1,line=2.5,cex=2.0)
mtext("a)",side=3,cex=1.5,xpd=NA,font=2,line=1.0,adj=0)
idx=1:22
text(xlab_try[idx],ylab_try[idx],names_1[idx],pos=c(4,2),cex=1.5)

plou_id=order(ylab_try)
plou_id=plou_id[plou_id>22]
idx=plou_id[c(1:2,9:10)]
text(xlab_try[idx],ylab_try[idx],names_1[idx],pos=2,cex=1.5)
posi=c(2,4)
for(i in 8:3){
        text(xlab_try[plou_id[i]]+0.25,ylab_try[plou_id[i]]+(i-8)*0.25,names_1[plou_id[i]],pos=4,cex=1.5)
        arrows(xlab_try[plou_id[i]]+0.4,ylab_try[plou_id[i]]+(i-8)*0.25,xlab_try[plou_id[i]],ylab_try[plou_id[i]],length=0)
}

#Second plot
xx=as.numeric(tab_answer[,"Prop signif"])
yy=as.numeric(tab_answer[,"MeanIntra"])/as.numeric(tab_answer[,"MeanInter"])
plot(xx,yy,pch=symbol,ylab="",xlab="",cex=2,cex.axis=1.5,xlim=c(-0.0,1),bg=cod_cod)
mtext("b)",side=3,cex=1.5,xpd=NA,font=2,line=1.0,adj=0)
legend("topleft",pch=c(16),leg=c('sparsity<0.65','sparsity<0.75','sparsity>=0.75'),col=c("cyan","blue","darkorchid"),bty='n',cex=1.5,inset=c(0.025,0.065),pt.cex=3)

names_1_bis=names_1
names_1_bis[c(2,8,19,28,24,25,31)]=""
text(xx,yy,names_1_bis,pos=c(4,2),cex=1.5)
idx=c(2,8)
text(xx[idx],yy[idx]-0.1,names_1[idx],pos=1,cex=1.5)

text(xx[19],yy[19],names_1[19],pos=2,cex=1.5)

text(xx[24]+0.01,yy[24]+0.1,names_1[24],pos=4,cex=1.5)
arrows(xx[24]+0.02,yy[24]+0.1,xx[24],yy[24],length=0)

text(xx[25]-0.02,yy[25]+0.25,names_1[25],pos=2,cex=1.5)
arrows(xx[25]-0.03,yy[25]+0.25,xx[25],yy[25],length=0)

text(xx[28]+0.04,yy[28]+0.175,names_1[28],pos=4,cex=1.5)
arrows(xx[28]+0.05,yy[28]+0.175,xx[28],yy[28],length=0)

text(xx[31]-0.04,yy[31]-0.025,names_1[31],pos=2,cex=1.5)
arrows(xx[31]-0.05,yy[31]-0.025,xx[31],yy[31],length=0)

mtext("Sparsity",side=1,line=2.5,cex=2.0)

dev.off()

#Comparing
#pdf("~/Documents/Plankton/REPHY_littoral/article/graphe/comparaison_ratio_code_log_cleaner_allwith0_ONLY_SIGNIF.pdf",width=10,height=8)
par(mfrow=c(1,1),xpd=NA,mar=c(3,4.5,1,1))

symbol=rep(21,dim(tab_answer)[1])
symbol[as.numeric(tab_answer[,"Dimension"])>=6]=22
symbol[as.numeric(tab_answer[,"Dimension"])>=10]=23


plot(xlab_try,log10(as.numeric(tab_answer[,"MeanIntra"])/as.numeric(tab_answer[,"MeanSignif"])),t="p",pch=symbol,bg=cod_cod,col="black",cex=2,xlab="",ylab="|intra|/|inter|",yaxt="n",ylim=c(-0.1,2),cex.lab=1.5,xaxt="n",xlim=c(0.8,max(xlab_try)+0.5))
lines(c(0.7,max(xlab_try)+0.6),rep(1,2),col="black",lty=2,lwd=2)
mtext("Chronological order",side=1,line=1.5,cex=1.5)
axis(2,at=c(0,log10(5),1,log10(50),2),lab=c("1","5","10","50","100"),cex.axis=1.5)
#names_1=tab_answer[,1]
#text(1:dim(tab_answer)[1],log10(as.numeric(tab_answer[,"MeanIntra"])/as.numeric(tab_answer[,"MeanSignif"]))-.05,names_1)
id=c(1,2,4:12,14:21)
text(xlab_try[id],log10(as.numeric(tab_answer[id,"MeanIntra"])/as.numeric(tab_answer[id,"MeanSignif"])),names_1[id],pos=c(2,4),cex=1.5)
id=3
text(xlab_try[id]-0.1,log10(as.numeric(tab_answer[id,"MeanIntra"])/as.numeric(tab_answer[id,"MeanSignif"])),names_1[id],pos=c(2,4),cex=1.5)
id=13
text(xlab_try[id],0.05+log10(as.numeric(tab_answer[id,"MeanIntra"])/as.numeric(tab_answer[id,"MeanSignif"])),names_1[id],pos=c(2,4),cex=1.5)


text(xlab_try[31:32],log10(as.numeric(tab_answer[31:32,"MeanIntra"])/as.numeric(tab_answer[31:32,"MeanSignif"])),names_1[31:32],pos=c(2,4),cex=1.5)
idx=c(22,23,29)
text(xlab_try[idx],log10(as.numeric(tab_answer[idx,"MeanIntra"])/as.numeric(tab_answer[idx,"MeanSignif"])),names_1[idx],pos=c(2,4),cex=1.5)
idx=30
text(xlab_try[idx],log10(as.numeric(tab_answer[idx,"MeanIntra"])/as.numeric(tab_answer[idx,"MeanSignif"])),names_1[idx],pos=c(4),cex=1.5)
id=order(log10(as.numeric(tab_answer[24:28,"MeanIntra"])/as.numeric(tab_answer[24:28,"MeanSignif"])),decreasing=T)
id=id+23
for(i in 1:length(id)){
	ylab_try=log10(as.numeric(tab_answer[id[i],"MeanIntra"])/as.numeric(tab_answer[id[i],"MeanSignif"]))
	text(xlab_try[id[i]]-0.2,ylab_try-i*0.05,names_1[id[i]],pos=2,cex=1.5)
	arrows(xlab_try[id[i]]-0.25,ylab_try-i*0.05,xlab_try[id[i]],ylab_try,length=0)
}
legend("bottomleft",pch=c(16,16,16,21,22,23),leg=c('sparsity<0.65','sparsity<0.75','sparsity>=0.75','Dimension<6','Dimension<10','Dimension>=10'),col=c("cyan","blue","darkblue",NA,NA,NA),bty='n',cex=1.5,lty=NA,lwd=2,pt.bg=c(NA,NA,NA,"black","black","black"))
#dev.off()

symbol=rep(21,dim(tab_answer)[1])

#pdf("~/Documents/Plankton/REPHY_littoral/article/graphe/sparsity_vs_others.pdf",width=18,height=8)
par(mfrow=c(1,2),xpd=NA,mar=c(4,4.5,3,0.5))
app_T=c(100,100,100,50,300,300,100,100,100,300,200,400,400,300,300,50,100,30,1000,700,300,500)
col=rep('black',length(app_T))
app_T=c(app_T,app_T_this_study)
col=c(col,rep('blue',length(app_T_this_study)))
par(mfrow=c(1,2))
plot(app_T,tab_answer[,'Prop signif'],col=col,pch=16,cex=2,cex.lab=2.0,cex.axis=1.5,ylab="Sparsity",xlab="Length of the time series")
mtext("a)",side=3,cex=1.5,xpd=NA,font=2,line=1.0,adj=0)
plot(tab_answer[,'Dimension'],tab_answer[,'Prop signif'],col=col,pch=16,cex.lab=2.0,cex=2,cex.axis=1.5,ylab="",xlab="Number of species")
mtext("b)",side=3,cex=1.5,xpd=NA,font=2,line=1.0,adj=0)
legend("bottomright",c("Other studies","This study"),pch=16,col=c("black","blue"),bty="n",cex=2)
#dev.off()

pdf('~/Documents/Plankton/REPHY_littoral/article/graphe/Ratio_function_dim.pdf',,width=10,height=8)
par(mfrow=c(1,1),mar=c(4,4.5,1,1))
plot(tab_answer[,"Dimension"],log10(as.numeric(tab_answer[,"MeanIntra"])/as.numeric(tab_answer[,"MeanSignif"])),t="p",pch=symbol,bg=cod_cod,col="black",cex=2,xlab="Number of taxa",ylab="|intra|/|inter|",yaxt="n",ylim=c(-0.1,2.1),cex.lab=1.5,xlim=c(1.9,14.5),cex.axis=1.5)
axis(2,at=c(0,log10(5),1,log10(50),2),lab=c("1","5","10","50","100"),cex.axis=1.5)

idx=1:32
tmp=names_1[18]
names_1[18]=NA
text(tab_answer[idx,'Dimension'],log10(as.numeric(tab_answer[idx,"MeanIntra"])/as.numeric(tab_answer[idx,"MeanSignif"])),names_1[idx],pos=c(4,2),cex=1.5)
idx=18
text(tab_answer[idx,'Dimension'],log10(as.numeric(tab_answer[idx,"MeanIntra"])/as.numeric(tab_answer[idx,"MeanSignif"])),tmp,pos=c(4),cex=1.5)
legend("bottomright",pch=c(16),leg=c('sparsity<0.65','sparsity<0.75','sparsity>=0.75'),col=c("cyan","blue","darkorchid"),bty='n',cex=1.5,lty=NA,lwd=2)
dev.off()
