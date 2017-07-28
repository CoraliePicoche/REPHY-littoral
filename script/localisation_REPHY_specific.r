#To locate only the sites we want to follow

graphics.off()
rm(list=ls())
library(ggplot2)
library(ggmap)

ref=read.table("data/raw/Q2_170419_Points_bis.csv",sep=";",header=TRUE)
corres=read.table("data/lieu_corres.csv",sep=";",header=TRUE)

sites=c("LEperon","Cornard","Men er Roue","Antoine","Loscolo","Bois de la Chaise large","Large Croisic","Auger","Lazaret")

lon=c()
lat=c()
for (s in sites){
	id=grep(s,ref$Lieu_libel)
	print(s)
	print(ref$Lieu_libel[id])
	lon=c(lon,ref$Longitude[id])	
	lat=c(lat,ref$Latitude[id])
}

pdf('graphe/localisation_specific_REPHY.pdf')
df=as.data.frame(cbind(lon,lat))
#stye='feature:administrative.country|element:labels|visibility:off'
mapgilbert=get_map(location=c(lon=mean(df$lon),lat=mean(df$lat)),zoom=6,maptype="satellite",scale=2)
g=ggmap(mapgilbert)+geom_point(data=df,aes(x=lon,y=lat,fill="red"),size=5,shape=21)+guides(fill=FALSE,alpha=FALSE,size=FALSE)

###All three versions work, but we want to regulate more precisely the exact place where label should be written
#g=g+annotate("text",x=jitter(lon,factor=15),y=jitter(lat,factor=15),label=sites)
#for (i in 1:length(sites)){
#g=g+geom_text(x=lon[i],y=lat[i],label=sites[i],hjust=-0.5)
#}
# df2=data.frame("Lon"=lon,"Lat"=lat,"Nom"=sites)
# g=g+geom_text(data=df2,aes(label=df2$Nom,x=df2$Lon,y=df2$Lat))

#Sort by lat
id=order(lat)
plou=0
for (i in id){
plou=plou+1
g=g+geom_text(x=lon[i],y=lat[i],label=sites[i],hjust=c(-0.15,1.15)[plou%%2+1],colour=I("white"))
}

print(g)

dev.off()
