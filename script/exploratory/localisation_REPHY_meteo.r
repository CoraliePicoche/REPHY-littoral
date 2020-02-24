#For each site, I want the closest Meteo France station

graphics.off()
rm(list=ls())
library(sp)

ref=read.table("data/local_sites.csv",sep=";",header=TRUE)
MF_db=read.table("data/meteofrance.csv",sep=";",header=TRUE)

sites=ref$Lieu_libel

mat=matrix(NA,nrow=length(sites),ncol=dim(MF_db)[1])
mat_sites=matrix(NA,nrow=length(sites),ncol=length(sites))
rownames(mat)=sites
colnames(mat)=MF_db$Name
rownames(mat_sites)=sites
colnames(mat_sites)=sites
for (s in 1:length(sites)){
	plou=matrix(c(MF_db$Longitude,MF_db$Latitude),ncol=2)
	plou_sites=matrix(c(ref$Longitude,ref$Latitude),ncol=2)
	mat[s,]=spDistsN1(plou,c(ref$Longitude[s],ref$Latitude[s]),longlat=TRUE)
	mat_sites[s,]=spDistsN1(plou_sites,c(ref$Longitude[s],ref$Latitude[s]),longlat=TRUE)
}
write.table(t(mat),file='data/rawdist_sites_meteofrance.csv',sep=';')
write.table(t(mat_sites),file='data/dist_intersites.csv',sep=';')

mat_bis=mat
for(s in 1:length(sites)){
	val=sort(mat[s,])[5]
	mat_bis[s,mat[s,]>val]=rep(NA,sum(mat[s,]>val))
}
write.table(t(mat_bis),file='data/smalldist_sites_meteofrance.csv',sep=';')


