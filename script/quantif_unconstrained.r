##CP 13/10/19 : Analysis of unconstrained matrix (and differences to pencen), for review JEcol

graphics.off()
rm(list=ls())
library(xtable)
source("script/matrix_MAR_clean.r")
results=array(NA,dim=c(10,7),dimnames=list(1:10,c("signif","signif outside","positive","mean abs(val)","mean abs(val) signif","mean absolute difference","transfo sign"))) 

id_g=0
id_l=0
groupe=c("BZ","MO","SU","AR")
for (g in groupe){
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
		id_l=id_l+1
                f1=paste("data/analyse_MAR/",g,"/site_specific/",option_lieu[ll],"_unconstrained_null_regular_common_",g,".RData",sep="")
                load(f1)
                B_unconstrained_s=clean_matrix(cis,diag_bool=F,signif=T)

		results[id_l,"signif"]=sum(B_unconstrained_s!=0.0)/length(B_unconstrained_s)
                

		B_unconstrained_ns=clean_matrix(cis,diag_bool=F,signif=F)
		results[id_l,"positive"]=sum(B_unconstrained_ns>00.0)/length(B_unconstrained_ns)
		results[id_l,"mean abs(val)"]=mean(abs(B_unconstrained_ns))
		

                f2=paste("data/analyse_MAR/",g,"/site_specific/",option_lieu[ll],"_pencen_null_regular_common_",g,".RData",sep="")
                load(f2)
                B_pencen=clean_matrix(cis,diag_bool=F,signif=F)
		struct=B_pencen*0
		struct[which(abs(B_pencen)>0.0,arr.ind=TRUE)]=1
		
		results[id_l,"signif outside"]=sum(B_unconstrained_s[struct==1]!=0)/length(B_unconstrained_s)
		results[id_l,"mean abs(val) signif"]=mean(abs(B_unconstrained_s[struct==1]))
		results[id_l,"transfo sign"]=sum(sign(B_unconstrained_ns)==sign(B_pencen))/length(B_unconstrained_ns)
		results[id_l,"mean absolute difference"]=mean(abs((B_unconstrained_ns-B_pencen)*struct))

	}
}

print.xtable(xtable(results,digits=c(2),display=c(rep("f",6),"E","f")),"article/submit_JEcol/response/quantif_unconstrained.tex" ,type="latex")

