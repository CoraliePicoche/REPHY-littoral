#13/01/2020 CP: compare abundances of species, because we use the ratio N_j/N_i when converting from linear MAR model to BH model

rm(list=ls())
graphics.off()
library('lubridate')
library('zoo')

set.seed(42)
timestep=14
consecutif=2

liste_sp=c("AST","CHA","DIT","GUI","LEP","NIT","PLE","PSE","RHI","SKE","THP","THL","GYM","PRO","PRP","SCR","CRY","EUG")

tab_sp=read.table('data/lieu_sp_post_reconstruct_pour_MAR.csv',header=TRUE,na.strings="",sep=";")
lieu=colnames(tab_sp)
lieu=gsub('.',' ',lieu,fixed=TRUE) #useful for Men er Roue

groupe=c("BZ","MO","AR","SU")

tab_tot=array(NA,dim=c(10,550,1+length(liste_sp)),dimnames=list(1:10,1:550,c("Date",liste_sp))) #10 sites, (max) sample points, up to 4 species and a date column

xlimi=c(as.Date("1996-01-01"),as.Date("2016-01-01"))
colo=c()
id_lieu=0
for (g in groupe){
        if(g=="BZ"){
                option_lieu=c("Men er Roue","Loscolo","Croisic")
                colo=c(colo,rep("green",length(option_lieu)))
        }else if(g=="MO"){
                option_lieu=c("LEperon","Cornard","Auger")
                colo=c(colo,rep("darkblue",length(option_lieu)))
        }else if(g=="SU"){
                option_lieu=c("Antoine","Lazaret")
                colo=c(colo,rep("darkred",length(option_lieu)))
        }else if(g=="AR"){
                option_lieu=c("Teychan","B7")
                colo=c(colo,rep("cyan",length(option_lieu)))
        }
        for (ll in 1:length(option_lieu)){
		id_lieu=id_lieu+1
		dimnames(tab_tot)[[1]][id_lieu]=option_lieu[ll]
		
                if(g=="AR"){
                        tab=read.csv(paste("./data/",option_lieu[ll],"_base.csv",sep=""),na.strings="NA",header=TRUE,sep=";",dec=".")
                }else{
                        tab=read.table(paste("data/corres_hernandez_",option_lieu[ll],'.txt',sep=''),sep=";",na="NA",header=TRUE)
                }

	        dates=as.Date(tab$Date)
        	tab=tab[year(dates)>=1996&year(dates)<2016,]#Using data from 1996
        	dates=dates[year(dates)>=1996&year(dates)<2016]

	        dates_bis=seq(dates[1],dates[length(dates)],timestep) #Regular time grid
		
		tab_tot[id_lieu,1:length(dates_bis),"Date"]=dates_bis
		
		if(g!="AR"){
		a=table(as.matrix(tab_sp[,lieu %in% option_lieu]))
                liste_sp_1=dimnames(a[a==length(option_lieu)])[[1]]
		liste_sp_1=liste_sp_1[liste_sp_1!="NEI"]
		}else{
			liste_sp_1=c("AST","NIT","PSE","SKE","CHA","GUI","LEP","RHI","GYM","PRP","CRY","EUG")
		}
		tab_plankton=na.approx(tab[,liste_sp_1],maxgap=consecutif,x=dates,xout=dates_bis,na.rm=FALSE) #Interpolation over regular time grid
#Replace missing values
                for (s in liste_sp_1){
                        tab_plankton[is.na(tab_plankton[,s]),s]=runif(sum(is.na(tab_plankton[,s])),0,min(tab_plankton[,s],na.rm=TRUE))
			tab_tot[id_lieu,1:length(dates_bis),s]=tab_plankton[,s]
                }
	}
}

#5 most abundant
xlimi=c(as.Date("1996-01-01"),as.Date("2017-06-01"))
liste_sp=c("AST","CHA","DIT","GUI","LEP","NIT","PLE","PSE","RHI","SKE","THP","THL","GYM","PRO","PRP","SCR","CRY","EUG")
list_centric=c("CHA","DIT","GUI","LEP","RHI","SKE","THP")
list_pennate=c("AST","NIT","PSE","PLE","THL")
list_dino=c("GYM","PRO","PRP","SCR")
colo=rainbow(5)
val=c()
for(l in 1:10){
	print("###########################")
        maxi=max(c(tab_tot[l,,liste_sp]),na.rm=T)
        plou=apply(tab_tot[l,,2:dim(tab_tot)[[3]]],2,mean,na.rm=T)
        plou=plou[!is.na(plou)]
        or=order(plou,decreasing=T)
        liste_sp_bis=names(plou)[or]
        for(s in 1:length(liste_sp_bis)){
		for(s_bis in (s+1):length(liste_sp_bis)){
			if(s!=s_bis){
			if((liste_sp_bis[s]%in%list_centric&liste_sp_bis[s_bis]%in%list_centric)|
				(liste_sp_bis[s]%in%list_pennate&liste_sp_bis[s_bis]%in%list_pennate)|
				(liste_sp_bis[s]%in%list_dino&liste_sp_bis[s_bis]%in%list_dino)){
				val=c(val,max(plou[s_bis]/plou[s],plou[s]/plou[s_bis]))
			}
			}
		}
        }
}

#pdf("article/submit_JEcol/response_R2/hist_ratio.pdf")
