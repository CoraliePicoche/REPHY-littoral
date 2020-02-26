#############
# CP 26/02/2020: Draws a graphical abstract for JEcol, based on the interactions in the Auger site.
#############

rm(list=ls())
graphics.off()

library("igraph") #For web visualization
library("stringr") #For string handling
#library("png") #For png insertion
source("script/function_plot_igraph.r") #I had to change a bit the igraph plot in order to modify the self-loop
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

#Edge color depends on the sign of the interaction, width depends on the weight in the interaction matrix
ltype=rep(1,length(E(g)$weight))
lili=get.edgelist(g)
essai_width=abs(E(g)$weight)*30 #Up to now, this one is the best
##Tried a lot of idfferent widths
#essai_width=10*scale(x)
#essai_width=(43-order(essai_width))/42*10
#essai_width=log10(abs(E(g)$weight))-min(log10(abs(E(g)$weight)))
#essai_width=sqrt(abs(E(g)$weight))
for (i in 1:length(E(g)$weight)){
        if (E(g)$weight[i]<0) {
                E(g)$color[i]='black'
                essai_width[i]=log10(1+abs(E(g)$weight[i]))
  #              essai_width[i]=sqrt(abs(E(g)$weight[i]))
                }
        if (E(g)$weight[i]>0) {
                E(g)$color[i]='red'
                ltype[i]=2
                essai_width[i]=log10(1+abs(E(g)$weight[i]))
 #               essai_width[i]=sqrt(abs(E(g)$weight[i]))
                }
        if(lili[i,1]==lili[i,2]){ #self_loop
                essai_width[i]=essai_width[i]/10
                essai_width[i]=log10(1+abs(E(g)$weight[i]))
#                essai_width[i]=sqrt(abs(E(g)$weight[i]))
        }
}
essai_width=essai_width*50

#The plot changes according to the algorithm output: I chose the seed
set.seed(4) #I want to avoid the central positioning of one of the group in the vizualisation
#l=layout_with_dh(g,weight.edge.lengths=0.0001) #Tried all layouts, this one is the best in our case.
#l=layout_as_star(g) #Too much emphasis on only one species
#l=layout_as_tree(g) #Does not work
l=layout_in_circle(g) #I like this one
#l=layout_nicely(g,dim=3) #UGLY
#l=layout_on_grid(g) #Not nice
#l=layout_on_sphere(g) #Not readable
#l=layout_randomly(g) #Not readable
#l=layout_with_drl(g,dim=3)#Bug
#l=layout_with_fr(g) #Bug
#l=layout_with_gem(g) #Clusters are too small
#l=layout_with_graphopt(g) #Adapted to large graphs, so makes our clusters too small
#l=layout_with_kk(g) #Does not work with negative values
#l=layout_with_lgl(g) #Does not adapt to cluster
#l=layout_with_mds(g) #Overlays two species
#l=layout_with_sugiyama(g) #Not adapted to cycle

#Change the self loop positioning so that it doesn't interfere with other links. This one has to be changed every time you change the layout
x=rep(0,length(E(g)))
x[1]=1.1
id=c(seq(1,36,length.out=6),seq(37,45,length.out=3),seq(46,49,length.out=2))
val=seq(x[1],-(2*pi-x[1]-(2*pi)/length(id)),length.out=length(id))
x[id]=val


pdf("article/submit_JEcol/response_R3/comparison_graphical_abstract_prop.pdf",width=7.5,height=7.5)
par(mar=c(0.5,0,0,0),mfrow=c(1,1))
#function_plot_igraph(g,edge.curved=.2,edge.arrow.size=.5,layout=l,edge.width=essai_width,edge.loop.angle=x,vertex.size=20,edge.lty=ltype)
#l=layout_with_dh(g,weight.edge.lengths=0.0001) #Tried all layouts, this one is the best in our case.
#function_plot_igraph(g,edge.curved=.2,edge.arrow.size=.5,layout=l,edge.width=essai_width,vertex.size=20,edge.lty=ltype,label=sp)
l=layout_in_circle(g) #I like this one
function_plot_igraph(g,edge.curved=.2,edge.arrow.size=.5,layout=l,edge.width=essai_width,vertex.size=20,edge.lty=ltype,label=sp,edge.loop.angle=x,vertex.color=adjustcolor("seagreen1", alpha.f = 1))
legend(x=-0.25,y=0.,c("+","-"),col=c("blue","red"),lty=c(2,1),cex=2,lwd=3,bty="n",text.col=c("blue","red"))
dev.off()
