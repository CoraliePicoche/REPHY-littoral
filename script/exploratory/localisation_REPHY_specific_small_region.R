#To locate only the sites we want to follow

graphics.off()
rm(list=ls())
library(ggplot2)
library(ggmap)

ref=read.table("data/raw/Q2_170419_Points_bis.csv",sep=";",header=TRUE)
corres=read.table("data/lieu_corres.csv",sep=";",header=TRUE)

meteofrance=read.table("data/meteofrance.csv",sep=";",header=TRUE)
id=meteofrance$Name
lonMF=c()
latMF=c()
for(i in 1:length(id)){
	lonMF=c(lonMF,meteofrance$Longitude[i])
	latMF=c(latMF,meteofrance$Latitude[i])
}

buoy=read.table("data/postesMarine.csv",sep=";",header=TRUE)
id=buoy$Nom
lonBO=c()
latBO=c()
for(i in 1:length(id)){
        lonBO=c(lonBO,buoy$Longitude[i])
        latBO=c(latBO,buoy$Latitude[i])
}

sites=c("LEperon","Cornard","Men er Roue","Antoine","Loscolo","Large Croisic","Auger","Lazaret")
groupe1=c("Antoine","Lazaret")
#groupe1=c("LEperon","Cornard","Auger")
#groupe1=c("Men er Roue","Large Croisic","Loscolo")
sites=groupe1

lon=c()
lat=c()
#for (s in sites){
for (s in groupe1){
	id=grep(s,ref$Lieu_libel)
	print(s)
	print(ref$Lieu_libel[id])
	lon=c(lon,ref$Longitude[id])	
	lat=c(lat,ref$Latitude[id])
}

pdf('Rapport/graphe/localisation_regional_REPHY_SU.pdf')
df=as.data.frame(cbind(lon,lat))
dfMF=as.data.frame(cbind(lonMF,latMF))
dfBO=as.data.frame(cbind(lonBO,latBO))
#stye='feature:administrative.country|element:labels|visibility:off'
mapgilbert=get_map(location=c(lon=mean(df$lon),lat=mean(df$lat)),zoom=8,maptype="satellite",scale=2)
g=ggmap(mapgilbert)+geom_point(data=df,aes(x=lon,y=lat),fill="red",colour="black",size=5,shape=21)+guides(fill=FALSE,alpha=FALSE,size=FALSE)

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

#For MeteoFrance
plou=0
for(i in 1:length(latMF)){
plou=plou+1
g=g+geom_text(x=lonMF[i],y=latMF[i],label=meteofrance$Name[i],hjust=c(-0.15,1.15)[plou%%2+1],colour=I("grey"),size=1.75)+geom_point(data=dfMF,aes(x=lonMF,y=latMF),fill="green",colour="black",size=2,shape=21)
}

#For Buoys
plou=0
for(i in 1:length(latBO)){
plou=plou+1
g=g+geom_text(x=lonBO[i],y=latBO[i],label=buoy$ID[i],hjust=c(-0.15,1.15)[plou%%2+1],colour=I("white"),size=1.75)+geom_point(data=dfBO,aes(x=lonBO,y=latBO),fill="cyan",colour="black",size=2,shape=21)
}


print(g)

dev.off()
