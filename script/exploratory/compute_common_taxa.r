############################
## 21/06/2019 CP: This script computes the mean number of taxa each site has with each other
###########################

rm(list=ls())
graphics.off()

tab_sp=read.table('data/lieu_sp_post_reconstruct_pour_MAR.csv',header=TRUE,na.strings="",sep=";")
lieu=colnames(tab_sp)
lieu=gsub('.',' ',lieu,fixed=TRUE) #useful for Men er Roue
groupe=c("Men er Roue","Loscolo","Croisic","LEperon","Cornard","Auger","Antoine","Lazaret","Teychan","B7")
mean_val=c()

for (l in 1:length(groupe)){
	groupe1=groupe[-l]
	a=table(as.matrix(tab_sp[,lieu %in% groupe1]))
        liste_sp=dimnames(a[a==length(groupe1)])[[1]]
	mean_val=c(mean_val,length(liste_sp))
}

