graphics.off()
rm(list=ls())
library('lubridate')

#We can look at all sites
#fil=list.files(path=".",pattern="^[A-Z].*flortot_only.csv")

#Or we can look only at sites that have previously been validated
sites=c("LEperon","Cornard","Men er Roue","Antoine","Loscolo","Croisic","Auger","Lazaret")
fil=c()
for (s in sites){
#	fil=c(fil,list.files(path="data/",pattern=paste(s,"*.flortot_only.csv")))
	fil=c(fil,list.files(path="data/",pattern=paste(s)))
}

maxi=20
name_place=c()
tab_abun=matrix(NA,nrow=maxi,ncol=length(fil))
tab_abun_quant=matrix(NA,nrow=58,ncol=length(fil))
tab_pres=matrix(NA,nrow=maxi,ncol=length(fil))
tab_pres_quant=matrix(NA,nrow=58,ncol=length(fil))
for (f in 1:length(fil)){
	tab=read.table(paste("data/",fil[f],sep=""),header=TRUE,sep=";")
	name_place=c(name_place,strsplit(strsplit(fil[f],split="_")[[1]][1],split="-")[[1]][1])
	
	#Two indices: most abundant groups over time ; and most present species over time

	t1=apply(tab[,2:dim(tab)[2]],2,sum,na.rm=TRUE)
	t1=t1/sum(t1)
	tab_abun_quant[,f]=t1
	t1=sort(t1,decreasing=TRUE)
	

	t2=apply(tab[,2:dim(tab)[2]],2,function(x) sum(is.na(x)|(x==0),na.rm=TRUE)/length(x))
	tab_pres_quant[,f]=t2
	t2=sort(t2)

	tab_abun[,f]=names(t1[1:maxi])
	tab_pres[,f]=names(t2[1:maxi])
}

rownames(tab_abun_quant)=colnames(tab[2:dim(tab)[2]])
rownames(tab_pres_quant)=colnames(tab[2:dim(tab)[2]])
print('Abundance')
#Warning: this is a strange mix of both ranking and averaging: i first count the number of times each groups is part of the top 20. It should be noted that 10 groups are ALWAYS in the top 20 in abundance. Then I rank them by abundance, which means that a species can appear only 6 times out of 8, but have a better rank than a species that appears 7 or 8 times out of 6 if it is really abundant in all samplings
print(sort(apply(tab_abun_quant[names(sort(table(tab_abun),decreasing=TRUE)),],1,mean),decreasing=TRUE))
print('Presence')
print(sort(apply(1-tab_pres_quant[names(sort(table(tab_pres),decreasing=TRUE)),],1,mean),decreasing=TRUE))
