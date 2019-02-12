library(MARSS)
library(corrplot)

groupe=c("BZ","MO","SU","AR")
#groupe=c("SU","AR")
sink(paste("article/graphe/correlation_between_interaction.txt",sep=""))
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
print(option_lieu[ll])
f1=paste("data/analyse_MAR/",g,"/site_specific/",option_lieu[ll],"_pencen_null_regular_common_",g,".RData",sep="")
load(f1)

tab_sp=cis$model$data
cntl.list=cis$call$control
model.list=cis$call$model
amethod="kem"

#print(paste('MARSS estimates',option_lieu[ll]))
fit_log=MARSS(tab_sp, method=amethod,model=model.list,control=cntl.list,silent=T) #Checked on Croisic, same results as before

#hh=MARSShessian(fit_log) #GIves a fancy MARSS object
#Or Maybe 
#print(paste('MARSS hessian',option_lieu[ll]))
ff=MARSSFisherI(fit_log) #Only gives a matrix but these are basically the same stuff

nn=which(grepl("^B.",dimnames(ff)[[1]]))

cc=cov2cor(solve(ff[nn,nn]))

#Completely arbitrary
cc_tmp=cc
cc_tmp[lower.tri(cc,diag=T)]=NA
id=which(abs(cc_tmp)>0.1,arr.ind=T)
print("0.1 raw")
print(dim(id)[1])
id2=which(abs(cc_tmp)>0.25,arr.ind=T)
print("0.25 proportion")
print(dim(id2)[1]/(dim(cc_tmp)[1]*(dim(cc_tmp)[1]-1)/2))
id3=which(abs(cc_tmp)>0.5,arr.ind=T)
print("0.5")
print(dim(id3)[1])
if(length(id)>0){
for(i in 1:dim(id)[1]){
	print(paste(rownames(cc_tmp)[id[i,1]],colnames(cc_tmp)[id[i,2]],cc_tmp[id[i,1],id[i,2]]))
}
}
write.table(cc,paste("article/graphe/",option_lieu[ll],"_corcoef.txt",sep=""),sep=";")
print('\n')
pdf(paste("article/graphe/correlation_between_interaction",option_lieu[ll],".pdf",sep=""),width=15,height=15)
corrplot(cc)
dev.off()
}
}
sink()
