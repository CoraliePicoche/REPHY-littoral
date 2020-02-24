library(date)

####################### FONCTIONS #################
#Correlation entre deux jeux de données. Par défaut, on suppose que l'indice 2 indique la plus longue série temporelle (celle sur laquelle on va interpoler les données)
correl=function(date_1,val_1,date_2,val_2){
#####METHODE 1 : on interpole sur l'ensemble des dates du jeu de données le plus grand
#plou=approx(date_1,val_1,date_2,method="linear")
##On trace l'interpolation pour verifier
#pdf("verif_interpolation.pdf")
#plot(date_2,plou$y,col="red",type="l")
#lines(date_2,val_2,type="l",lty=2,col="blue")
#dev.off()
#cc=cor(plou$y,val_2,use="pairwise.complete.obs",method="spearman")
#cct=cor.test(plou$y,val_2,method="spearman",na.rm=TRUE)

####METHODE 2 : on ne fait les corrélations que sur les données mesurées sur le terrain, en ne gardant que les dates communes
serie_correl1=rep(NA,length(date_1)) #Contient toutes les observations du jeu de données 1 pour lesquelles les mesures ont été faites en même temps pour 1 et 2
serie_correl2=rep(NA,length(date_1)) #Contient toutes les observations du jeu de données 2 pour lesquelles les mesures ont été faites en même temps pour 1 et 2
jmin=1
for (i in 1:length(date_1)){
        d1=date_1[i]
        j=jmin
        delai=d1-date_2[j]
        delai_av=d1-date_2[j]
        while (j<length(date_2)&&abs(delai)<=abs(delai_av)){
                delai_av=delai
                j=j+1
                delai=d1-date_2[j]
        }
        if (abs(delai)>abs(delai_av)){
                jmin=j-1
        }
        else{
                jmin=j
        }
        delai=d1-date_2[jmin]
        if (delai==0){
                serie_correl1[i]=val_1[i]
                serie_correl2[i]=val_2[jmin]
        }
        }
cc=cor(serie_correl1,serie_correl2,use="pairwise.complete.obs",method="spearman")
cct=cor.test(serie_correl1,serie_correl2,method="spearman",na.rm=TRUE) #Test si  la valeur est significative, on peut printer le résultat pour info
return(cc)
}

#Pour chaque date de mesure du jeu de données qui nous intéresse (flore, en l'occurrence), on regarde le délai minimum avec la mesure hydro faite un peu plus loin
delai=function(date_1,date_comp){
delaitab=rep(NA,length(date_1))
jmin=1
for (i in 1:length(date_1)){
        d1=date_1[i]
        j=jmin
        delai=d1-date_comp[j]
        delai_av=d1-date_comp[j]
        while (j<length(date_comp)&&abs(delai)<=abs(delai_av)){
                delai_av=delai
                j=j+1
                delai=d1-date_comp[j]
        }
        if (abs(delai)>abs(delai_av)){
                jmin=j-1
        }
        else{
                jmin=j
        }
        delaitab[i]=d1-date_comp[jmin]
}
        return(delaitab)
}

#Permet de faire le graphe des deux jeux de données 1 (considérées comme le plus long/complet) et 2, avec les légendes associées l1 et l2, ainsi que des paramètres graphiques (title est le titre du graphe et du fichier, unit est l'unite affichee dans les ordonnees. Le delai minimal entre date1 et date2 peut aussi être affichée 
plot_tout=function(date_1,obs_1,l1,date_2,obs_2,l2,delai,title,unit){
datee=as.Date(c("01/01/1986","01/01/1990","01/01/1995","01/01/2000","01/01/2003","01/01/2005","01/01/2008","01/01/2010","01/01/2013"),"%d/%m/%Y")
fileplot=paste(title,".pdf",sep="")
pdf(fileplot)
res=correl(date_2,obs_2,date_1,obs_1) #Calcul de la correlation entre les jeux de données d'après l'une des deux méthodes précédentes, pour l'afficher aussi dans la légende à titre indicatif
legtxt=c(l1,l2,format(res,digit=3))
colo=c("red","blue","white")

par(mar=c(5,4,4,5)+.1)
plot(date_1,obs_1,type="l",cex=0.5,col=colo[1],xlab="Date",ylab=unit,main=title,xaxt="n",xlim=c(date_2[1],date_2[length(date_2)]))
lines(date_2,obs_2,type="p",pch=16,col=colo[2],cex=0.5,xaxt="n")
axis.Date(1,at=datee,format(date_1[1],"%Y"),cex.axis=.7)

par(new=TRUE)
#Impression du délai entre deux observations proches des jeux 1 et 2
plot(date_2,abs(delai),type="p",cex=0.4,pch=17,xaxt="n",yaxt="n",xlab="",ylab="",col="green")
axis(4)
mtext("delay",side=4,line=3)

legend("topright",legtxt,pch=c(NA,16,NA),lty=c(1,NA,NA),col=colo)
dev.off()
}

#Graphe dynamique des especes de teychan_flore (matrice avec dates en ligne, espèce en colonnes) en fonction de date_teych. On représente nb_espece en tout dont st par graphe. La période temporelle est découpée en nb_frame. L'ordre des espèces est données dans comptage_sans, ti est le titre du graphe, long est la longueur du fichier pdf en sortie et relation est la fonction de lien (log10+1 par défaut,ctr pour centrée, crr pour centrée réduite, rel pour relatif, auc pour aucune transformation). 
graphe_dynamique=function(st,nb_espece,nb_frame,date_teych,teychan_flore,comptage_sans,ti,long=25,relation="log"){
spespe=seq(0,nb_espece,st) #Definition du nombre de graphes dont a besoin pour représenter toutes les espèces avec un maximum de st par graphe
if(spespe[length(spespe)]!=nb_espece){spespe=c(spespe,nb_espece)}
pdf(ti,width=long)
#Boucle sur les graphes par espèces
for (e in 1:(length(spespe)-1)){
liste_mini=comptage_sans[(spespe[e]+1):spespe[e+1]] #Liste des espèces à représenter sur la série de graphes considérées
liste_mini=liste_mini[!is.na(liste_mini)]
liste_mini_pour_legend=liste_mini
for (j in 1:length(liste_mini)){
	if(grepl("\\+",liste_mini_pour_legend[j])){
		liste_mini_pour_legend[j]=paste(strsplit(liste_mini_pour_legend[j],split="\\+")[[1]][1],"+...",sep="")
	}
}
colo=rainbow(length(liste_mini))
colo=sample(colo)
#Decoupage du temps en nb_frame
year_step=(date_teych[length(date_teych)]-date_teych[1])/nb_frame
plou=0:nb_frame
lim=date_teych[1]+year_step*plou
#Boucle temporelle sur les frames
for (j in 1:nb_frame){
eval(parse(text=paste("a=teychan_flore[,'",liste_mini[1],"']",sep="")))
a0=a
#On sépare les valeurs à 0 et les valeurs différentes de 0 qui seront plotées différemment
a0[a!=0]=NA
a1=a
a1[a==0]=NA
if (relation=="log"){
	if(sum(!is.na(a0))>0){
		a0=jitter(log10(a0+1))
	}
	a1=log10(a1+1)
	miniy=0
	maxiy=log10(max(teychan_flore+1,na.rm=TRUE))
	transfo="(Log10+1)"
}
else if (relation=="auc"){
	if(sum(!is.na(a0))>0){
		a0=jitter(a0)
	}
	a1=a1
	miniy=0
	maxiy=max(teychan_flore,na.rm=TRUE)+1
	transfo="réelle"
}
else if (relation=="rel"){
	TOTO=rowSums(teychan_flore)
	if(sum(!is.na(a0))>0){
		a0=jitter(a0/TOTO)
	}
	a1=a1/TOTO
	miniy=0
	maxiy=1.0
	transfo="relative"
}
else if (relation=="ctr"||relation=="crr"){
	stop("Le traitement pour centrées et centrées réduites n'est pas finalisé, puisque je ne suis pas sûre de devoir compter les 0 comme des 0 ou des nan")
}

#Le layout sert à tracer correctement la légende en dehors du graphe pour une meileure lisibilité
layout(matrix(c(1,2),nrow=1),width=c(0.85,0.15))
par(mar=c(2, 4, .5, .5))
#On trace la première espèce à part pour définir le cadre
plot(date_teych,a1,type='o',pch=16,cex=0.7,col=colo[1],xlim=c(lim[j],lim[j+1]),ylim=c(miniy,maxiy),ylab=paste("Abondance",transfo,sep=" "),xaxt='n')
lines(date_teych,a0,type='p',pch=16,cex=0.7,col=colo[1])
axis.Date(1,date_teych,at=seq(as.Date(paste("01/01/",format(lim[j],"%Y"),sep=""),"%d/%m/%Y"),as.Date(paste("01/01/",format(lim[j+1],"%Y"),sep=""),"%d/%m/%Y"),by="year"),"%Y")
axis.Date(1,date_teych,at=seq(as.Date(paste("01/",format(lim[j],"%m/%Y"),sep=""),"%d/%m/%Y"),as.Date(paste("01/",format(lim[j+1],"%m/%Y"),sep=""),"%d/%m/%Y"),by="month"),labels=FALSE,tck=-0.01)
#Boucle sur les autres espèces
if (length(liste_mini)>1){for (i in 2:length(liste_mini)){
if(!is.na(liste_mini[i])){
eval(parse(text=paste("a=teychan_flore[,'",liste_mini[i],"']",sep="")))
a0=a
#On sépare les valeurs à 0 et les valeurs différentes de 0 qui seront plotées différemment
a0[a!=0]=NA
a1=a
a1[a==0]=NA
if (relation=="log"){
	if(sum(!is.na(a0))>0){
 	       a0=jitter(log10(a0+1))
	}
        a1=log10(a1+1)
        miniy=0
        maxiy=log10(max(teychan_flore+1,na.rm=TRUE))
        transfo="(Log10+1)"
}
else if (relation=="auc"){
	if(sum(!is.na(a0))>0){
        	a0=jitter(a0)
	}
        a1=a1
        miniy=0
        maxiy=max(teychan_flore,na.rm=TRUE)+1
        transfo="réelle"
}
else if (relation=="rel"){
        TOTO=rowSums(teychan_flore)
	if(sum(!is.na(a0))>0){
        	a0=jitter(a0/TOTO)
	}
        a1=a1/TOTO
        miniy=0
        maxiy=1.0
        transfo="relative"
}
else if (relation=="ctr"||relation=="crr"){
        stop("Le traitement pour centrées et centrées réduites n'est pas finalisé, puisque je ne suis pas sûre de devoir compter les 0 comme des 0 ou des nan")
}

lines(date_teych,a1,type='o',pch=16,cex=0.7,col=colo[i])
lines(date_teych,a0,type='p',pch=16,cex=0.7,col=colo[i])
}}}
par(mar=c(2, 0, .5, .5))
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n",axes=FALSE,ann=FALSE)
legend("topright",liste_mini_pour_legend,lty=rep(1,length(liste_mini)),bty="n",col=colo)
}
}
dev.off()
}

#Differents types de moyenne
#Moyenne d'abondance sur nbmonth_span mois à partir de date_deb. Cette moyenne se fait sur l'ensemble de la période date_teych.
moyenne_temporelle=function(date_teych,abundance,date_deb,nbmonth_span){
require(lubridate)
nb_jours_tot=date_teych[length(date_teych)]-date_deb
nb_ligne=nb_jours_tot/(nbmonth_span*30.5)
#abundance[is.na(abundance)]=0.0
mama=matrix(0,nrow=nb_ligne,ncol=length(date_teych))
mamabis=matrix(0,nrow=nb_ligne,ncol=dim(abundance)[2])
dd=rep(as.Date("01/01/1970","%d/%m/%Y"),nb_ligne)
for (i in 1:nb_ligne){
        d0=date_deb
        month(d0)=month(d0)+nbmonth_span*(i-1)
        d1=date_deb
        month(d1)=month(d1)+nbmonth_span*i
        dd[i]=as.Date(d0+(d1-d0)/2,"%d/%m/%Y")
        for (j in 1:length(date_teych)){
                if(date_teych[j]>=d0&&date_teych[j]<=d1){
                        mama[i,j]=1
			mamabis[i,]=mamabis[i,]+abundance[j,]
                }
        }
}
#abundance_bis=(mama%*%abundance)/(rowSums(mama,na.rm=TRUE))
abundance_bis=mamabis/(rowSums(mama,na.rm=TRUE))
colnames(abundance_bis)=colnames(abundance)
rownames(abundance_bis)=dd
return(abundance_bis)
}

#Calcul de l'anomalie d'abundance par rapport à abundance_moyenne
anomalie=function(date_teych,abundance,abundance_moyenne){
dada=as.Date(as.numeric(rownames(abundance_moyenne)),origin="1970-01-01")
abundance_bis=matrix(0,nrow=dim(abundance)[1],ncol=dim(abundance)[2])
for (i in 1:dim(abundance)[2]){
                for (j in 1:length(date_teych)){
                        m=month(date_teych[j])
                        for (k in 1:(length(dada)-1)){
                                if(m>=month(dada[k])&&m<=month(dada[k+1])){
                                        abundance_bis[j,i]=abundance[j,i]-abundance_moyenne[k]
                                }
                        }
        }
}
colnames(abundance_bis)=colnames(abundance)
rownames(abundance_bis)=date_teych
return(abundance_bis)
}

#moyenne glissante sur nb_mesure
moyenne_glissante=function(date_teych,abundance,nb_mesure){
require(caTools)
abundance_bis=apply(abundance,2,function(x) runmean(x,nb_mesure))
dd=sapply(date_teych,function(x) runmean(x,nb_mesure))
colnames(abundance_bis)=colnames(abundance)
rownames(abundance_bis)=dd
return(abundance_bis)
}

#Calcul moyenne sur une année type de l'abondance par nb_mois (typiquement, nb_mois=1 pour avoir une idée de la moyenne mensuelle, mais on peut élarger ce pas de temps)
moyenne_mensuelle=function(date_teych,abundance,nb_mois){
require(lubridate)
nbspan=12/nb_mois
abundance_bis=matrix(0,nrow=nbspan,ncol=dim(abundance)[2])
nbmesure=matrix(0,nrow=nbspan,ncol=dim(abundance)[2])
dd=rep(as.Date("01/01/00","%d/%m/%y"),nbspan)
for (i in 1:nbspan){
        mm0=1+(i-1)*nb_mois
        mm1=1+i*nb_mois
        month(dd[i])=floor((mm0+mm1)/2)
        for (j in 1:length(date_teych)){
                if (month(date_teych[j])>=mm0&&month(date_teych[j])<=mm1){
                        for (k in 1:dim(abundance)[2]){
                                if(!is.na(abundance[j,k])){
                        		abundance_bis[i,k]=abundance_bis[i,k]+abundance[j,k]
                                        nbmesure[i,k]=nbmesure[i,k]+1
                                }
                        }
                }
        }
}
abundance_bis=abundance_bis/nbmesure
colnames(abundance_bis)=colnames(abundance)
rownames(abundance_bis)=dd
return(abundance_bis)
}

#On veut avoir toutes les espèces correspondant à un ensemble de critère, et pour une classification donnée. Par exemple, toutes les diatomées au niveau du genre ; ou encore toutes les diatomées formant des colonies au niveau des familles ; ou encore tous les dinoflagellés formant des familles, classées en régime alimentaire. On lui passera donc abundance raw, avec les sous-groupes que l'on retrouve dans atlas
sous_groupe_advanced=function(critere,classification,tableau){
	load("/home/cpicoche/Documents/Plankton/data/treated/atlas.RData")
#Vérification de cohérence entre les entrées et ce qui est disponible dans l'atlas
	li_sp1=colnames(tableau)
	li_sp2=atlas$ID
	if(sum(li_sp1==li_sp2)<length(li_sp1)||length(li_sp1)!=length(li_sp2)){
		stop("Les espèces du tableau données et les identifiants dans l'atlas ne sont pas les mêmes")
	}
	li_cr1=dimnames(critere)[[2]]
	li_cr2=dimnames(atlas)[[2]]
	for (j in li_cr1){
		 if (sum(grepl(j,li_cr2))==0){
			print(j)
			stop("L'un des critères définis n'est pas dans les headers de atlas")
		}
	}

#On met à 0 dans les abondances toutes les espèces qui ne nous intéressent pas	
	for (l in 1:length(li_sp1)){
		id=which(atlas$ID==li_sp1[l])
		corres=TRUE
		for (j in li_cr1){
			pouic=as.character(eval(parse(text=paste("atlas$",j,"[",id,"]",sep=""))))
			pouicbis=as.character(eval(parse(text=paste("critere$",j,"[1]",sep=""))))
			corres=corres&&pouic==pouicbis
		}
		if(!corres){    #Si on ne remplit pas toutes les conditions, on met à 0
			tableau[,l]=rep(0,length(tableau[,l]))
		}
	}

#Maintenant, on classe au niveau demandé
#Si la classification est taxonomique, on utilise les matrices déjà calculées (c'est plus pratique pour les regroupements au niveau du genre, par exemple)
	if (sum(grepl(classification,c("class","order","family","genus")))>0){
		load("/home/cpicoche/Documents/Plankton/data/treated/correspondance_taxonomy.RData")
		eval(parse(text=paste("ab=tableau%*%Mcorres_",classification,sep="")))
		TOTO=rowSums(ab)
		co=order(colSums(ab,na.rm=TRUE)/sum(TOTO,na.rm=TRUE),decreasing=TRUE)
		ab=ab[,co]
		eval(parse(text=paste("colnames(ab)=",classification,"name[co]",sep="")))
	}
	else{
		class=as.character(eval(parse(text=paste("atlas$",classification,sep=""))))
		class_type=unique(class)
		print(classification)
		Mcorres_class=matrix(data=NA,nrow=length(colnames(tableau)),ncol=length(class_type))
		for (i in 1:length(class_type)){
		        Mcorres_class[grep(class_type[i],class),]=0
		        Mcorres_class[grep(class_type[i],class),i]=1
		}
		ab=tableau%*%Mcorres_class
		colnames(ab)=class_type
	}

	tmp=colSums(ab)
	id=which(tmp==0)[1]
	ab=as.matrix(ab[,1:(id-1)])
	rownames(ab)=as.character(rownames(tableau),"%d/%m/%Y")
        return(ab)
}

#Fonction plus utilisé
sous_groupe=function(liste_inf,liste_sup,qui,abundance_inf){
        liste_nom=c()
        for (j in 1:length(liste_inf)){
                if((liste_sup[j]==qui)&&sum(grepl(liste_inf[j],liste_nom))==0){
                        liste_nom=c(liste_nom,liste_inf[j])
                }
        }
        ab=matrix(0,nrow=dim(abundance_inf)[1],ncol=length(liste_nom))
        for (i in 1:length(liste_nom)){
                x=grep(liste_nom[i],colnames(abundance_inf),fixed=TRUE)
		if(length(x)>0){
                ab[,i]=abundance_inf[,x]}
        }
        colnames(ab)=liste_nom
        rownames(ab)=rownames(abundance_inf)
        return(ab)
}

replace_na=function(mat,consecutif,val0){
#        dates=as.Date(as.numeric(rownames(mat)),origin="1970-01-01",tz="GMT")
	dates=as.Date(rownames(mat))
        species=colnames(mat)
        mat_bis=matrix(NA,nrow=dim(mat)[1],ncol=dim(mat)[2],dimnames=dimnames(mat))
        xli1=min(dates)
        xli2=max(dates)+500
        for (s in 1:length(species)){
                i=2
                med=median(mat[,s],na.rm=TRUE)
                matchange=mat[,s]
                mis_val0=c()
                while(i<(length(mat[,s])-consecutif+1)){
                        j=1
                        consecutif_NA=is.na(mat[i,s])
                        if((consecutif_NA)&(j<=(consecutif))){
                                consecutif_NA=consecutif_NA&is.na(mat[i+j,s])
                                j=j+1
                        }
                        if((!consecutif_NA)&(is.na(mat[i,s]))){
                        #Il y a entre 1 et consecutif-1 valeurs à NA
                                if(!is.na(mat[(i-1),s])){
                                        if((mat[(i-1),s]>med)|(mat[(i+j-1),s]>med)){
#Cas 1-b)                       
                                                matchange[(i-1):(i+j-1)]=c(approx(x=c(dates[(i-1)],dates[(i+j-1)]),y=c(mat[i-1,s],mat[i+j-1,s]),xout=dates[(i-1):(i+j-1)],method="linear"))$y
                                       }else{
#Cas 1-a)       
                                                matchange[(i):(i+j-2)]=rep(val0,length((i):(i+j-2)))
                                                mis_val0=c(mis_val0,i:(i+j-2))
                                        }
                                }else{
                                                matchange[(i-1):(i+j-1)]=c(approx(x=c(dates[(i-1)],dates[(i+j-1)]),y=c(val0,mat[i+j-1,s]),xout=dates[(i-1):(i+j-1)],method="linear"))$y
                                                mis_val0=c(mis_val0,i-1)
                                }
                                j=j-1
                        }else if(consecutif_NA){
                                if (!is.na(mat[i+j,s])){
#Il y a exactement consecutifs valeurs à NA
                                        if(!is.na(mat[(i-1),s])){
                                                if((mat[(i-1),s]>med)|(mat[(i+j),s]>med)){
#Cas 1-b)                               
                                                        matchange[(i-1):(i+j)]=c(approx(x=c(dates[(i-1)],dates[(i+j)]),y=c(mat[i-1,s],mat[i+j,s]),xout=dates[(i-1):(i+j)],method="linear"))$y
                                                }else{
#Cas 1-a)
                                                        matchange[(i):(i+j-1)]=rep(val0,length((i):(i+j-1)))
                                                        mis_val0=c(mis_val0,i:(i+j-1))
                                                }
                                        }else{
#Si i-1 est aussi égale à NA, on a plus de consecutif valeurs à NA      
                                                        matchange[(i-1):(i+j-1)]=rep(val0,length((i-1):(i+j-1)))
                                        }
                                }else{
#Il y a plus de consecutif valeurs à NA
                                        matchange[(i):(i+j)]=rep(val0,length((i):(i+j)))
                                        mis_val0=c(mis_val0,i:(i+j))
                                        if ((i+j)<length(dates)){
                                                j=j+1
                                                consecutif_NA=is.na(mat[i+j,s])
                                        }
                                        while((consecutif_NA)&((i+j)<length(dates))){
                                                matchange[i+j]=val0
                                                mis_val0=c(mis_val0,i+j)
                                                j=j+1
                                                consecutif_NA=is.na(mat[i+j,s])
                                        }
                                }
                        }
                        i=i+j
                }
                id=which(is.na(matchange))
                for (i in id){
                        if ((i==1)|(i==length(matchange))|(i==(length(matchange)-1))){matchange[i]=val0}
                        else{
                                print(species[s])
                                print(i)
                        }
                }
                mat_bis[,s]=matchange
        }
        return(mat_bis)
}

resampling_regulier=function(mat,timestep){
        dates=as.Date(as.numeric(rownames(mat)),origin="1970-01-01",tz="GMT")
#	dates=as.Date(rownames(mat))
        xli1=min(dates)
        xli2=max(dates)+500
        timereconstruct=seq(dates[1],dates[length(dates)],timestep)
        species=colnames(mat)
        mat_bis=matrix(NA,nrow=length(timereconstruct),ncol=dim(mat)[2],dimnames=list(timereconstruct,species))
        for (s in 1:length(species)){
                mat_bis[,s]=c(approx(x=dates,y=mat[,s],xout=timereconstruct,method="linear"))$y
        }
        return(mat_bis)
}

mat_taxon_pour_comparaison=function(my_tab){
	require("stringr")
	path_at="/home/cpicoche/Documents/Plankton/data/"
        regroupement=read.csv(paste(path_at,"groupe_daniele_REPHY_moi.csv",sep=""),header=TRUE,sep=";",na.strings="NA")
        groupe_rephy=strsplit(as.character(regroupement$CodeREPHY),"/")
        groupe_rephy=lapply(groupe_rephy,str_trim)
        regroupement$CodeDaniele=as.character(lapply(regroupement$CodeDaniele,str_trim))
        classement=order(regroupement$CodeDaniele)
        regroupement$CodeDaniele=regroupement$CodeDaniele[classement]
        groupe_rephy=groupe_rephy[classement]
	ay=my_tab
        a=sort(unique(ay$Date))
        if(length(ay$Taxon)!=0){
        mat_FLORTOT=matrix(NA,nrow=length(a),ncol=(length(regroupement$CodeDaniele)+1)) #On rajoute à la fin tous les inclassables
        for (i in 1:length(ay$Date)){
		found=FALSE
                if((length(grep("Pseudo-nitzschia",ay$Taxon[i]))==0)&(length(grep("Ceratium tripos",ay$Taxon[i]))==0)&(length(grep("Dinophysis",ay$Taxon[i]))==0)&(length(grep("Penn",ay$Taxon[i]))==0)){
                        id_tax=grep(ay$Taxon[i],groupe_rephy,fixed=TRUE)
			if(length(id_tax>1)){
				words=strsplit(as.character(ay$Taxon[i]),split=" ")
				z=1
				while((z<=length(id_tax))&(!found)){
					yy=length(groupe_rephy[id_tax[z]])
					uu=1
					while((uu<=yy)&(!found)){	
						found=(ay$Taxon[i] %in% groupe_rephy[id_tax[z]][[uu]])
						uu=uu+1
					 }
	                        	z=z+1
                                }
				if(found){
					id_tax=id_tax[z-1]
				}
		}
                }else{
                        if(length(grep("Pseudo-nitzschia",ay$Taxon[i]))>0){
				found=TRUE
                                id_tax=grep("Pseudo-nitzschia",groupe_rephy,fixed=TRUE)
                        }
                        if(length(grep("Dinophysis",ay$Taxon[i]))>0){
				found=TRUE
                                id_tax=grep("Dinophysis",groupe_rephy,fixed=TRUE)
                        }
                        if(length(grep("Ceratium trip",ay$Taxon[i]))>0){
				found=TRUE
                                id_tax=grep("Ceratium trip",groupe_rephy,fixed=TRUE)
                        }
                        if(length(grep("Penn",ay$Taxon[i]))>0){
				found=TRUE
                                id_tax=grep("Penn",groupe_rephy,fixed=TRUE)
                        }
                }
                id_date=grep(ay$Date[i],a)
		if(found){
				mat_FLORTOT[id_date,id_tax]=sum(ay$Val[i],mat_FLORTOT[id_date,id_tax],na.rm=TRUE)
		}else{
				mat_FLORTOT[id_date,dim(mat_FLORTOT)[2]]=sum(ay$Val[i],mat_FLORTOT[id_date,dim(mat_FLORTOT)[2]],na.rm=TRUE)
			}
        }
                rownames(mat_FLORTOT)=a
                colnames(mat_FLORTOT)=c(regroupement$CodeDaniele,"NEI")
        }else{
                mat_FLORTOT=NULL
        }
        mat_tax=mat_FLORTOT
        return(mat_tax)
}

mat_conv_genus=function(my_tab){
        path_at="/home/cpicoche/Documents/Plankton/data/"
        regroupement=read.csv(paste(path_at,"groupe_daniele_REPHY_moi.csv",sep=""),header=TRUE,sep=";",na.strings="NA")
	regroupement=regroupement[,-dim(regroupement)[2]] #to delete comments
	nb_group=sort(unique(regroupement$Groupement_Code))
	mat_tax=matrix(0,nrow=dim(my_tab)[1],ncol=length(nb_group))
	mat_tax_not_na=matrix(FALSE,nrow=dim(my_tab)[1],ncol=length(nb_group))
        colnames(mat_tax)=nb_group
	for (i in 1:dim(my_tab)[2]){
		id=which(regroupement$CodeDaniele==colnames(my_tab)[i])
		id=which(nb_group==regroupement$Groupement_Code[id])
		id_not_na=!is.na(my_tab[,i])
		mat_tax[id_not_na,id]=mat_tax[id_not_na,id]+my_tab[id_not_na,i]
		mat_tax_not_na[id_not_na,id]=mat_tax_not_na[id_not_na,id]+TRUE
	}
	mat_tax[!mat_tax_not_na]=NA
	mat_tax=cbind(mat_tax,my_tab[,"NEI"])
	colnames(mat_tax)[(length(nb_group)+1)]="NEI"
        rownames(mat_tax)=rownames(my_tab)
	return(mat_tax)
}

sh_expo=function(tab){
        H=rep(NA,dim(tab)[1])
        for (i in 1:dim(tab)[1]){
                total=sum(tab[i,],na.rm=TRUE)
                prop=rep(0,dim(tab)[2])
                if (total>0){
                        prop=tab[i,]/total
                        for (j in 1:length(prop)){
				if(!is.na(prop[j])){
                                if (prop[j]>0){
                                        prop[j]=-prop[j]*log(prop[j])
                                }
				}
                        }
                }
                H[i]=exp(sum(prop,na.rm=TRUE))
        }
	H=H/dim(tab)[2]
        return(H)
}
find_best_pacf=function(x,y=NULL,maxi=15){
        if(length(y)>0){
                lag_max=min(maxi,length(x)-1,floor(10*log10(length(x))),length(y)-1,floor(10*log10(length(y)))) #Taking the same threshold as ar
                x=c(x,rep(NA,lag_max+1),y) #concatenate with no dependence between x and y
        }else{
                lag_max=min(maxi,length(x)-1,floor(10*log10(length(x)))) #Taking the same threshold as ar
        }
        p=pacf(x,plot=FALSE,na.action=na.pass,lag.max=lag_max)
        significance_level <- qnorm((1 + 0.95)/2)/sqrt(sum(!is.na(x)))
        last=tail(which(abs(p$acf)>=significance_level),1)
        plot(p$acf)
        abline(h=significance_level,col="blue")
        abline(h=-significance_level,col="blue")
        points(p$lag[last],p$acf[last],col="red",pch=16)
        return(p$lag[last])
}

find_best_arma=function(x,y=NULL,type="p0",maxi=15){
        if(length(y)>0){
                lag_max=min(maxi,length(x)-1,floor(10*log10(length(x))),length(y)-1,floor(10*log10(length(y)))) #Taking the same threshold as ar
                x=c(x,rep(NA,lag_max+1),y) #concatenate with no dependence between x and y
        }else{
                lag_max=min(maxi,length(x)-1,floor(10*log10(length(x)))) #Taking the same threshold as ar
        }
        aic_mat=rep(NA,lag_max)
        for (l in 1:lag_max){
#               print(paste("Lag max",l))
		if(type=="p0"){
			q=0
		}else if (type=="pp"){
			q=l
		}else if (type=="pp-1"){
			q=l-1
		}else{
			stop("Precise the ARMA type")
		}
                tryCatch({
	                aic_mat[l]=arima(x,order=c(l,0,q),method="ML",optim.method="BFGS",optim.control=list(maxit=1000))$aic #Changed optim.method to avoid problems for non-finite difference ; changed the max iterations to try and converge: not always sufficient (but should only be a problem for the longer time lags which are not necessarily meaningful
        	        },error=function(err){
	                print(paste("For lag",l,err))
        	        })

        }
        plot(aic_mat,main=type)
        val=which(aic_mat==min(aic_mat,na.rm=TRUE))
        points(val,min(aic_mat,na.rm=TRUE),col="red",pch=16)
        return(aic_mat)
}

new_covar=function(tab,dates,accum_vent,start="t-1"){
        path_data_post="./data/treated/"
	load(paste(path_data_post,"vent.RData",sep=""))
	load(paste(path_data_post,"meteo.RData",sep=""))
	date_meteo=as.Date(meteo$Date)
	load(paste(path_data_post,"nao.RData",sep=""))
	
	cov_bis=c("CHL","CumDebit","MeanNAO","Ntot","PHEO","P/N","CumPrec","CumRg","SAL","SI/N","TEMP","MeanVent","PHOS","SI","MES","NH4","NOx","Rg","NAO_month","AMO_month")
	tab_cov=matrix(NA,nrow=dim(tab)[1],ncol=length(cov_bis),dimnames=list(rownames(tab),cov_bis))
	tab_cov[,"Rg"]=tab[,"Rg"]
	tab_cov[,"NAO_month"]=tab[,"NAO_month"]
	tab_cov[,"AMO_month"]=tab[,"AMO_month"]
	tab_cov[,"CHL"]=tab[,"CHL"]
	tab_cov[,"MES"]=tab[,"MES"]
	tab_cov[,"NOx"]=tab[,"NOx"]
	tab_cov[,"NH4"]=tab[,"NH4"]
	tab_cov[,"SI"]=tab[,"SI"]
	tab_cov[,"PHOS"]=tab[,"PHOS"]
	tab_cov[,"Ntot"]=tab[,"NH4"]+tab[,"NOx"]
	tab_cov[,"PHEO"]=tab[,"PHEO"]
	tab_cov[,"P/N"]=tab[,"PHOS"]/tab_cov[,"Ntot"]
	id=is.infinite(tab_cov[,"P/N"])
	tab_cov[id,"P/N"]=NA
	tab_cov[,"SAL"]=tab[,"SAL"]
	tab_cov[,"SI/N"]=tab[,"SI"]/tab_cov[,"Ntot"] #Chose this one instead of Ntot because less correlation with Ntot and PHOS
	id=is.infinite(tab_cov[,"SI/N"])
	tab_cov[id,"SI/N"]=NA
	tab_cov[,"TEMP"]=tab[,"TEMP"]

#Start defines the period over which meteorological variables are defined: it can be from t-1 to t in the case of correlation between variables ; or it can be from t to t+1 when doing an AR model

	if(start=="t-1"){
		sta=2
		sto=length(dates)
	}else if(start=="t"){
		sta=1
		sto=length(dates)-1
	}else{
		stop("You should define the period over which you wish to integrate your physical data")
	}

	for (d in sta:sto){
		if(start=="t-1"){
        		d1=dates[d-1]
	        	d2=dates[d]
		}else if(start=="t"){
        		d1=dates[d]
	        	d2=dates[d+1]
		}
	        tab_cov[d,"MeanNAO"]=mean(nao$Val[nao$Date>=d1&nao$Date<d2],na.rm=TRUE)
        	tab_cov[d,"CumPrec"]=tail(cumsum(meteo$Precipitation[date_meteo>=d1&date_meteo<d2]),1)
	        tab_cov[d,"CumRg"]=tail(cumsum(meteo$Rg[date_meteo>=d1&date_meteo<d2]),1)
	        tab_cov[d,"CumDebit"]=tail(cumsum(meteo$Debit.Eyre[date_meteo>=d1&date_meteo<d2]),1)
	        tab_vent=rep(NA,d2-d1)
	        for (dd in 1:(d2-d1)){
        	        d_ap=d1+dd
                	d_av=d_ap-accum_vent
	                id=which(vtot$date>=d_av&vtot$date<=d_ap)
        	        tab_vent[dd]=sum(as.numeric(vtot$vitesse[id])^2)
	        }
        	tab_cov[d,"MeanVent"]=mean(tab_vent,na.rm=TRUE)
}
	return(tab_cov)
}

#This function turns a raw concentration to a "stress" function, using the half-saturation value
stress_function=function(nit,threshold){
	lambda=log(2)/threshold
	vec=1-exp(-lambda*nit)
	return(vec)
	}

#This functions returns a random interaction matrix, with obligate diagonal interactions and the same probability of turning on and off the off-diagonal elements
initial_random_matrice=function(rowna,colna){
        tab1=matrix(list(0),nrow=length(rowna),ncol=length(colna))
        colnames(tab1)=colna
        rownames(tab1)=rowna
        for(i in 1:length(rowna)){
                for(j in 1:length(colna)){
                        if(i==j){
                                tab1[i,j]=paste(colnames(tab1)[j],rownames(tab1)[i],sep="")
                        }else{
                                tab1[i,j]=sample(list(0,paste(colnames(tab1)[j],rownames(tab1)[i],sep="")),1)[[1]]
                        }
                }
        }
        return(tab1)
}

#This function takes an interaction matrix and returns another one with only one element changed, with the same probability of removing or adding an interaction (once again, we don't touch the diagonal
followin_random_matrice=function(tab1){
        tab2=tab1
        val=sample(c("add","remove"),1)
        if(val=="add"){
                ind=which(tab1==0,arr.ind=TRUE)
                i=sample(ind[1],1)
                j=sample(ind[2],1)
                while(i==j){
                        i=sample(ind[1],1)
                        j=sample(ind[2],1)
                }
                tab2[i,j]=paste(colnames(tab1)[j],rownames(tab1)[i],sep="")
        }else if(val=="remove"){
                ind=which(is.na(tab1==0),arr.ind=TRUE)
                i=sample(ind[1],1)
                j=sample(ind[2],1)
                while(i==j){
                        i=sample(ind[1],1)
                        j=sample(ind[2],1)
                }
                tab2[i,j]=0
        }
        return(tab2)
}

#This fnctions takes and interaction matrix and the indices of the element of interest ; and it turns the estimation on/off, if it was off/on respectively
inverse_random_matrice=function(tab1,index1,index2){
        tab2=tab1
        if(index1!=index2){
                if(tab1[index1,index2]==0){
			tab2[index1,index2]=paste(colnames(tab1)[index2],rownames(tab1)[index1],sep="")
		}else{
			tab2[index1,index2]=0
		}
        }else{
		stop("I cannot change the diagonal")
        }
	return(tab2)
}

