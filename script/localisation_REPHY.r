#To locate all points on one map

graphics.off()
rm(list=ls())
library(ggplot2)
library(ggmap)

ref=read.table("data/raw/Q2_170419_Points_bis.csv",sep=";",header=TRUE)
corres=read.table("data/lieu_corres.csv",sep=";",header=TRUE)

filename_tot=c("data/raw/Q2_170418_site_Tania_UTF8.csv")

for (filename in filename_tot){
tab=read.table(filename,sep=";",header=TRUE)

l_u=unique(tab$Lieu_id)
lon=c()
lat=c()
name=c()
for (l in 1:length(l_u)){
	lon=c(lon,ref$Longitude[ref$Lieu_id==l_u[l]])	
	lat=c(lat,ref$Latitude[ref$Lieu_id==l_u[l]])
	name=c(name,as.character(tab$Lieu_libel[tab$Lieu_id==l_u[l]][1]))

}
}

pdf('graphe/localisation_specific_REPHY.pdf')
df=as.data.frame(cbind(lon,lat))
mapgilbert=get_map(location=c(lon=mean(df$lon),lat=mean(df$lat)),zoom=5,maptype="toner-lite",source="stamen",scale=2)
g=ggmap(mapgilbert)+geom_point(data=df,aes(x=lon,y=lat,fill="red"),size=5,shape=21)+guides(fill=FALSE,alpha=FALSE,size=FALSE)+annotate("text",x=jitter(lon,factor=5)-1.5,y=jitter(lat,factor=5)+rep(c(-0.5,0.5),16),label=name)#geom_text(aes(x=lon,y=lat,label=name),size=2,hjust=-0.5,vjust=0)
print(g)
dev.off()
