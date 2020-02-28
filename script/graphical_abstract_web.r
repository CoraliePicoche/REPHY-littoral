#############
# CP 26/02/2020: Draws a graphical abstract for JEcol, based on the interactions in the Auger site.
#############

rm(list=ls())
graphics.off()

library("igraph") #For web visualization
library("stringr") #For string handling
library("qgraph")
source("script/matrix_MAR_clean.r")



site="Auger"
g="MO"

f1=paste("data/analyse_MAR/",g,"/site_specific/",site,"_pencen_null_regular_common_",g,".RData",sep="")
load(f1)
sp=dimnames(cis$model$data)[[1]] #Studied species
B=clean_matrix(cis)
rownames(B)=sp
colnames(B)=sp

h=B
h=t(h) #Warning: it seems that for directed graph, igraph used A(i,j) where i acts on j
g=graph.adjacency(h,weighted=TRUE,mode="directed") #You can set diag to FALSE if you don't want to intraspecific competition

c_vertex=adjustcolor("seagreen3", alpha.f =.75)
c_plus="royalblue"

#Edge color depends on the sign of the interaction, width depends on the weight in the interaction matrix
ltype=rep(1,length(E(g)$weight))
lili=get.edgelist(g)
essai_width=rep(NA,length(E(g)$weight)) #Up to now, this one is the best
for (i in 1:length(E(g)$weight)){
        if (E(g)$weight[i]<0) {
                E(g)$color[i]='black'
                essai_width[i]=log10(1+abs(E(g)$weight[i]))
  #              essai_width[i]=sqrt(abs(E(g)$weight[i]))
                }
        if (E(g)$weight[i]>0) {
                E(g)$color[i]=c_plus
                ltype[i]=2
                essai_width[i]=log10(1+abs(E(g)$weight[i]))
 #               essai_width[i]=sqrt(abs(E(g)$weight[i]))
                }
        if(lili[i,1]==lili[i,2]){ #self_loop
                essai_width[i]=log10(1+abs(E(g)$weight[i]))
#                essai_width[i]=sqrt(abs(E(g)$weight[i]))
        }
}
essai_width=essai_width*60

l=layout_in_circle(g)

#Change the self loop positioning so that it doesn't interfere with other links. This is flexible enough to adapt to different number of species 
x=rep(0,length(E(g)))
x[1]=1.1
id=c()
li=as_edgelist(g, names = FALSE)
id=c()

for(i in 1:nrow(li)){
	if(li[i,1]==li[i,2]){
		id=c(id,i)
	}
}
val=seq(x[1],-(2*pi-x[1]-(2*pi)/length(id)),length.out=length(id))
x[id]=val


pdf("article/submit_JEcol/response_R3/comparison_graphical_abstract_prop_one_color_mix.pdf",width=10,height=10)
par(mar=c(4,4,4,4),mfrow=c(1,1))
qg=qgraph(get.adjacency(g,sparse=FALSE),layout=l,diag=TRUE,directed=TRUE,esize=1,edge.width=essai_width,curveDefault=.3,asize=1.9,vsize=7.5,lty=ltype,edge.color=E(g)$color,color=c_vertex,border.color=NA,arrowAngle=pi/5)
legend(x=-0.25,y=0.,c("+","-"),col=c(c_plus,"black"),lty=c(2,1),cex=2,lwd=3,bty="n",text.col=c(c_plus,"black"))

#Dummy graph for the line width in the legend
qgraph(cbind(1:2,2:1,(60*c(log10(1.5),log10(1.05)))),esize=1,DoNotPlot=TRUE)
lwd <- qg$graphAttributes$Edges$width 
legend(x=-0.2,y=-0.5,c("0.5","0.05"),col=c("black","black"),lty=c(1),cex=1.5,lwd=lwd,bty="n",text.col=c("black"))

dev.off()
