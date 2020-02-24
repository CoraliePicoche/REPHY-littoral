#CP 04/10/17 : Using results from boostrapped wavelets to test for significant synchrony/compensation, trying to get the same kind of results as Fig. 2 in Vasseur et al. 2014

rm(list=ls())
graphics.off()

groupe=c("BZ","MO","SU","AR")
id_lieu=0
id_g=0
time_between_samples=c()
total_period_sampling=c()
annee_min=c()
annee_max=c()
values=c()
values_synchro=c()
values_comp=c()
init_matrix=matrix(0,nrow=463,ncol=55) #I'm clearly cheating
for (g in groupe){
        if(g=="BZ"){
                option_lieu=c("Men er Roue","Loscolo","Croisic")
        }else if(g=="MO"){
                option_lieu=c("LEperon","Cornard","Auger")
        }else if(g=="SU"){
                option_lieu=c("Antoine","Lazaret")
        }else if(g=="AR"){
                option_lieu=c("Teychan")#,"B7")
        }


        for (ll in 1:length(option_lieu)){

		filename=paste(option_lieu[ll],"_synch_vass_keitt.RData",sep="")
#		save(mm,ww,seq_time,seq_scale,dates,file=filename)
		load(filename)

		condition=1:length(seq_scale)
#		condition=which(seq_scale>=365.25/2&seq_scale<2*365.25)
#		condition=which(seq_scale<365.25/2)
#		condition=which(seq_scale>365.25*2)
		

		z_tmp=ww$z[,condition,]
		z_boot=ww$z.boot[,,condition]

		values=c(values,ww$z) #Marche pas
		find_synchro=which(z_tmp[,]>z_boot[2,,])
		find_comp=which(z_tmp[,]<=z_boot[2,,])
		values_synchro=c(values_synchro,unlist(z_tmp[find_synchro]))
		values_comp=c(values_comp,unlist(z_tmp[find_comp]))

		init_matrix[ww$z[,,1]>ww$z.boot[2,,]]=init_matrix[ww$z[,,1]>ww$z.boot[2,,]]+1
		
	}
}
new_tab=matrix(NA,nrow=length(values),ncol=2)
new_tab[1:length(values_synchro),1]=values_synchro
new_tab[1:length(values_comp),2]=values_comp

h=hist(new_tab,plot=F,breaks=seq(0,1,0.1))
h_synch=hist(new_tab[,1],plot=F,breaks=seq(0,1,0.05))
h_synch$counts=h_synch$counts/sum(h$counts)
h_comp=hist(new_tab[,2],plot=F,breaks=seq(0,1,0.05))
h_comp$counts=h_comp$counts/sum(h$counts)

plot(h_synch,col="red")
plot(h_comp,col="blue",add=T)


#Test of the binomial
library(fields)
proba_matrix=dbinom(init_matrix,9,0.5) 
proba_corrected=p.adjust(c(proba_matrix), method ="BH")
dim(proba_corrected)=c(463,55)
image.plot(seq_time,seq_scale,proba_corrected)
contour(seq_time,seq_scale,proba_corrected,add=T,levels=c(0.025,0.975))

