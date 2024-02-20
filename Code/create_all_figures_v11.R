setwd("~/Dropbox/Code/Inversity3/")
plist = c("doMC","sna","Matrix","igraph","xtable","plotrix")
lapply(plist, require, character.only=TRUE)
registerDoMC(detectCores()-2)

source("create_all_figures_functions.R")
affiliation = c("actor-movie","brunson_club-membership","ca-cit-HepPh","com-dblp")
human_inperson = c("moreno_innovation","moreno_health","contact","sociopatterns-infectious")
human_online = c("arenas-pgp","flickrEdges","advogato","youtube-u-growth","orkut-links","facebook-wosn-links","loc-brightkite_edges")
human_online = c("arenas-pgp","flickrEdges","advogato","loc-brightkite_edges","facebook-wosn-links","munmun_twitter_social")
human_online = c("arenas-pgp","flickrEdges","advogato","munmun_twitter_social")
computer = c("topology","web-Google","p2p-Gnutella31")
#computer = c("topology","p2p-Gnutella31")
infrastructure = c("opsahl-powergrid","opsahl-usairport","roadNet-CA")
biological = c("maayan-vidal","reactome", "moreno_propro", "arenas-meta")

n_affiliation = c("Actor-Movie","Club Mmebers","Citation (Physics)","Citation (CS)")
n_social = c("Physician","Adolescent","Contact","Conference")
n_online = c("PGP Users","Flickr","Advogato","Twitter")
n_computer = c("Internet Topology","WWW (Google)","Gnutella P2P")
n_infrastructure = c("Power Grid","US Airports","CA Roads")
n_biological = c("Human Protein 1","Human Protein 2","Yeast Protein","C. Elegans")
net.realnames = c(n_affiliation,n_social,n_online,n_computer,n_infrastructure,n_biological)

net_all = c(affiliation,human_inperson,human_online,computer,infrastructure,biological)

categories = c("Affiliation","In person Social","Onine Social","Computer","Infrastructure","Biological")

# setwd("~/Dropbox/Code/Inversity3")

#####################################################################################
# Uncomment to create and save image file
#####################################################################################

net_list = networks_to_image(net_all)
rsample = get_samples_from_networks(net_list,5000)
dsample = get_sample_degrees(net_list,rsample)
evsample = get_sample_eigencentrality(net_list,rsample)
# #betweenness_sample = get_sample_betweenness(net_list,rsample)
net_list$originalnames = net_all
net_list$categories = categories
net_list$names = net.realnames
crun = ceiling(as.numeric(Sys.time()))
ifile = paste("networks_koblenz_N=",length(net_all),"_",crun,".image",sep="")
save(list=c("net_list","rsample","dsample","evsample"),file=ifile)
# load(ifile)

#####################################################################################
###########################################################################################
# Get Basic Network Statistics (if not already DONE!)
###########################################################################################

setwd("~/Dropbox/Code/Inversity3")
nstats_all = get_all_network_stats(net_list)
crun = ceiling(as.numeric(Sys.time()))
ifile = paste("networks_koblenz_N=",length(net_all),"_",crun,".image",sep="")
save(list=c("net_list","rsample","dsample","evsample","nstats_all"),file=ifile)


#####################################################################################
# Load image file that's already been created
#####################################################################################
setwd("~/Dropbox/Code/Inversity3/")
#ifile = "networks_koblenz_N=23_1495038731.image"
# ifile = "networks_koblenz_N=22_1495136866.image"
ifile = "networks_koblenz_N=22_1587772354.image"
load(ifile)


#####################################################################################
# Print Network Details into LaTeX
#####################################################################################
# library(xtable)
categories = net_list$categories
ncategories = length(categories)
categories[2] = "F2F Social"
categories
num_in_category = c(length(affiliation),length(human_inperson),length(human_online),length(computer),length(infrastructure),length(biological))
net_list$names[18]="CA Roads"
catnames = NULL;
for (k in 1:ncategories)
{
  catnames = c(catnames,rep(categories[k],num_in_category[k]))
}
catshort = c("A","FS","OS","C","I","B")
netlabels = generate_network_labels(catshort,num_in_category)

get_edges = function(net_list)
{
  nnetworks = length(net_list$a)
  edges_all = rep(0,nnetworks)
  for (k in 1:nnetworks){
    xnet = net_list$a[[k]]
    edges_all[k] = sum(sign(xnet))/2
  }
  return(ceiling(edges_all))
}
options(digits=1)
netlabels2 = paste0("(",netlabels,")")
edges = round(get_edges(net_list)); nodes = unlist(nstats_all[,1]); dmax=unlist(nstats_all[,8]); dmin=unlist(nstats_all[,9])
net_df = data.frame(Category=catnames,Label=netlabels,"Network Name"=net_list$names,Nodes=nodes,Edges=edges,Min.Degree=dmin,Max.Degree=dmax)
net_df = data.frame(Label=netlabels2,"Network Name"=net_list$names,Nodes=nodes,Edges=edges,Min.Degree=0,Max.Degree=1)
print(xtable(net_df),include.rownames=FALSE)


#####################################################################################

nrow = length(net_list$categories)
num_in_category = c(length(affiliation),length(human_inperson),length(human_online),length(computer),length(infrastructure),length(biological))
ncol = max(num_in_category)


###########################################################################################
# Set the Figures Directory to Store Figures
###########################################################################################

setwd("~/Dropbox/Code/Inversity3/Figures")

###########################################################################################
# Figure F1: Network Motifs
###########################################################################################
F_SINGLECOMPONENT = 0
# Do the network diagrams separately
labelsize = 1.5/2
ecol=gray(0.75)
N=8
vcex=1+1.5*(1/(2*N)); vcol="darkgray"
pcex.names=0.65
pdf("F2_NetworkMotifs4.pdf",family="Times");
#par(mfrow=c(4,3))
# pnames = rev(c("Degree Mean","Local Mean","Global Mean"))
pnames = rev(c("Degree Mean","Ego-based Mean","Alter-based Mean"))
pcol = rev(c("gray","red","darkgreen"))
mlayout = make_layout_matrix(c(4,4,4,4,4),FALSE,TRUE)
nf=layout(mlayout,respect=TRUE)
#layout.show(nf)
par(mar=0.0*c(1,0,1,0))
catnames = c("Original", "Ego-based","Alter-based","")
startx_k = c(2,4,4)
for (k in 1:length(catnames))
{
  plot(0,xlim=c(0,20),type="n",bty="n",ylim=c(0,20),xaxt="n",yaxt="n",ann=FALSE)
  text(startx_k[k],2,catnames[k],font=2,cex=1.5,pos=4)
  #legend("bottom",bty="n",xjust=1,catnames[k],text.font=2,cex=2)
}

barspace = 3
par(mar=0.5*c(1,0,1,0))
# pdf("LowLocal_LowGlobal.pdf");
tgraph = watts.strogatz.game(1, 16, 2, 0, loops = FALSE, multiple = FALSE)
tnet = as.matrix(get.adjacency(tgraph))
g.ll=gplot(tnet,gmode="graph",vertex.col=vcol,vertex.cex=vcex,edge.col=ecol,label.cex=labelsize)
vplot_network_local(tnet,g.ll)
vplot_network_global(xnet=tnet,g.ll)
tnetstats = get_network_stats(tnet)
bdata = c(tnetstats$mu_g,tnetstats$mu_l,tnetstats$mu_d)
par(mar=c(3, 5, 2, 1))
barplot(bdata,horiz=TRUE,beside=TRUE,col=pcol,names.arg=pnames,cex.names=pcex.names,las=1,space=barspace)

par(mar=0.5*c(1,0,1,0))
tnet <- create.network(N, local="H", global="L");
g.hl=gplot(tnet, gmode="graph",vertex.col=vcol,vertex.cex=vcex,edge.col=ecol,label.cex=labelsize)
vplot_network_local(xnet=tnet,g.hl)
vplot_network_global(xnet=tnet,g.hl)
tnetstats = get_network_stats(tnet)
bdata = c(tnetstats$mu_g,tnetstats$mu_l,tnetstats$mu_d)
par(mar=c(3, 5, 2, 1))
barplot(bdata,horiz=TRUE,beside=TRUE,col=pcol,names.arg=pnames,cex.names=pcex.names,las=1,space=barspace)

par(mar=0.5*c(1,0,1,0))
tnet <- create.network(N, local="L", global="H");
g.lh = gplot(tnet, gmode="graph",vertex.col=vcol,vertex.cex=vcex,edge.col=ecol,label.cex=labelsize)
vplot_network_local(xnet=tnet,g.lh)
vplot_network_global(xnet=tnet,g.lh)
tnetstats = get_network_stats(tnet)
bdata = c(tnetstats$mu_g,tnetstats$mu_l,tnetstats$mu_d)
par(mar=c(3, 5, 2, 1))
barplot(bdata,horiz=TRUE,beside=TRUE,col=pcol,names.arg=pnames,cex.names=pcex.names,las=1,space=barspace)

par(mar=0.5*c(1,0,1,0))
tnet <- create.network(N, local="H", global="H");
g.hh=gplot(tnet, gmode="graph",vertex.col=vcol,vertex.cex=vcex,edge.col=ecol,label.cex=labelsize)
vplot_network_local(xnet=tnet,g.hh)
vplot_network_global(xnet=tnet,g.hh)
tnetstats = get_network_stats(tnet)
bdata = c(tnetstats$mu_g,tnetstats$mu_l,tnetstats$mu_d)
par(mar=c(3, 5, 2, 1))
barplot(bdata,horiz=TRUE,beside=TRUE,col=pcol,names.arg=pnames,cex.names=pcex.names,las=1,space=barspace)

#get_network_stats(tnet)
for (k in 1:5)
{
  plot(0,xlim=c(0,10),type="n",bty="n",ylim=c(0,10),xaxt="n",yaxt="n",ann=FALSE)
  #if (k>1) #text(9,5,LETTERS[k-1],font=2,cex=3)
  if (k>1) legend("center",bty="n",LETTERS[k-1],text.font=2,cex=2.4)
}
dev.off()
###########################################################################################
# Figure F2: Plot out eCDF for Random, Local and Global
###########################################################################################
# For each network or class, define a color or pch symbol.
crun = ceiling(as.numeric(Sys.time()))

pdf_fname = paste("eCDF_",length(net_all),"_Networks_",crun,".pdf",sep="")
pdf(pdf_fname,family="Times")
mlayout = make_layout_matrix(num_in_category,extra_row=TRUE)
nf = layout(mlayout,respect=TRUE)
#layout.show(nf)
nplot = sum(num_in_category)
xlab="degree"; ylab="ecdf";
#mtitle="Empirical CDF"
pcex=0.1; pcexaxis=0.5; pcexmain=0.75
par(mar=c(3,1,1,0))
pcol=c("darkgray","red","blue")
for (k in 1:nplot)
{
  mtitle=net_list$names[k]
  cy = dsample[[k]]
  plot(ecdf(cy$r),verticals=TRUE,cex.main=pcexmain,cex.axis=pcexaxis,cex=pcex,col=pcol[1],pch=4,main=mtitle,xlab=xlab,ylab=ylab)
  plot(ecdf(cy$l),verticals=TRUE,pch=2,cex=pcex,col=pcol[2],add=TRUE)
  plot(ecdf(cy$g),verticals=TRUE,pch=5,cex=pcex,col=pcol[3],add=TRUE)
  grid()
  # cymedian = c(median(cy$r),median(cy$l),median(cy$g))
  # abline(v=cymedian[1],col="black")
  # abline(v=cymedian[2],col="blue",lty=2,lwd=2)
  # abline(v=cymedian[3],col="red",lty=3,lwd=2)
  abline(h=0.5,col="black",lty=2,lwd=2)
}
for (k in 1:length(categories))
{
  #par(mar=c(1,5,1,0))
  plot(0, xlim=c(0,10),ylim=c(0,10),type="n",ann=FALSE,axes=FALSE)
  #legend("center",bty="n",categories[k],font=2)
  catname = paste(categories[k]," ")
  text(9,5,catname,font=2,pos=2,cex=0.8)
}
ltext = c("Random", "FoF Local","FoF Global")
ltext = c("Random", "Ego-based","Alter-based")
# plot(0, xlim=c(0,10),ylim=c(0,10),type="n",ann=FALSE,axes=FALSE)
plot.new()
legend("top",bty="n",ltext,horiz=TRUE,col=pcol,pch=15,cex=1.25)
dev.off()

###########################################################################################
# Figure F3: Heatmap Plot of Density (Random versus FoF)
###########################################################################################
crun = ceiling(as.numeric(Sys.time()))
pdf_fname = paste("SmoothScatter_",length(net_all),"_Networks_",crun,".pdf",sep="")
pdf(pdf_fname,family="Times")
par(mar=c(2,1,1,1))
mlayout = make_layout_matrix(num_in_category,TRUE,TRUE)
nf = layout(mlayout,respect=TRUE)
#layout.show(nf)
xlab="Degree"; ylab="FoF Local";
#Lab.palette <- colorRampPalette(c("white", "orange", "red","blue"), space = "Lab")
plot_smooth_scatter(net_list)
plot_categories(categories,fontsize=0.84)
color.bar(Lab.palette(100),0,1)
dev.off()

###########################################################################################
# Figure F4: Comparison Across Networks for Degree
# (a) Local mean versus Random Sample
# (b) Global mean versus Random Sample
# (c) Local versus Global Means
###########################################################################################


####################################################################################
# Generate a figure of the sampling process
####################################################################################

###########################################################################################
# Figure F5: Comparison Across Networks for Eigenvector Centrality
# (a) Local mean versus Random Sample
# (b) Global mean versus Random Sample
# (e) Local versus Global Means
###########################################################################################

###########################################################################################
# Figure F6: ER WS and BA Graphs
# Can you combine it with previous analysis in Figure F4 for Plotting?
###########################################################################################


###########################################################################################
# Figure F4: Comparison Across Networks for Degree
# (a) Local mean versus Random Sample
# (b) Global mean versus Random Sample
# (c) Local versus Global Means
###########################################################################################
########################################
# (a) Local mean versus Random Sample
########################################
networksize=unlist(nstats_all[,1])
networkedges=networksize*unlist(nstats_all[,2])/2
networkdensity=networkedges/((networksize-1)*(networksize))
degree.mean = unlist(nstats_all[,2])
local.mean = unlist(nstats_all[,4])
global.mean = unlist(nstats_all[,3])
lcex = log(networksize)/2
lpch = 16
#col_cat = c("#4682B4","red","darkgreen","darkorange","#9ACD32","#4B0082")
col_cat = c("#4682B4","#FF0000","#006400","#EE7600","#9ACD32","#4B0082")

ccol = NULL;
for (k in 1:length(num_in_category))
  ccol = c(ccol,rep(col_cat[k],num_in_category[k]))
lcol = ccol
lcol = paste0(ccol,"9f")

crun = ceiling(as.numeric(Sys.time()))
pdf_fname = paste("F4_NetworkComparison_",length(net_all),"_Networks_",crun,".pdf",sep="")
pdf(pdf_fname,family="Times")

mlayout = rbind(c(1,2),c(4,3))
par(oma=c(0,1,2,0))
layout(mlayout)

####################################################################################
# Generate a figure of the sampling process
####################################################################################
nnodes = 25
xnet = rgraph(nnodes,tprob=0.2,mode="graph")
#gplot(xnet,gmode="graph")
glabels=letters[1:nnodes]
pcol=c("black","red","darkgreen")
inodes = sample(1:nnodes,size=3,replace=FALSE)
vcol=rep("lightgray",nnodes)
ecol=matrix("lightgray",nnodes,nnodes)
vcol[inodes]=pcol[1]
rsample = sample(which(xnet[inodes[2],]==1),size=1)
vcol[rsample]=pcol[2]
gsample = which(xnet[inodes[3],]==1)
vcol[gsample]=pcol[3]
ecol[inodes[3],] = pcol[3]; ecol[,inodes[3]] = pcol[3]
#g1=gplot(xnet,gmode="graph",mode="kamadakawai",label=glabels,vertex.cex=2,vertex.col=vcol,edge.col=ecol)
#pdf("F10_Strategies.pdf")
#gplot(xnet,coord=g1,gmode="graph",mode="kamadakawai",label=glabels,vertex.cex=3,vertex.col=vcol,edge.col=ecol)
gplot(xnet,gmode="graph",mode="kamadakawai",label=glabels,vertex.cex=2,vertex.col=vcol,edge.col=ecol)



####################################################################################
# Figure F7: Global and Local Leverage in Real Networks
####################################################################################
cat_names = c("A","FS","OS","C","I","B")
num_in_category = c(4,4,4,3,3,4)
netlabels = generate_network_labels(cat_names,num_in_category)
length(netlabels)

# col_cat = c("#4682B4","red","darkgreen","darkorange","#9ACD32","#4B0082")
col_cat = c("#4682B4","#FF0000","#006400","#EE7600","#9ACD32","#4B0082")
ccol = NULL;
for (k in 1:length(num_in_category))
  ccol = c(ccol,rep(col_cat[k],num_in_category[k]))


plot(0,type="n",bty="n",ann=FALSE,axes=FALSE)
leg.txt = paste0("(",cat_names,") ",net_list$categories)
legend("bottomright",bty="n",leg.txt,fill=col_cat,cex=0.8)#,col=col_cat)
# mtext("A",outer=FALSE,cex=1.5,font=2,adj=-0.60)

####################################################################################

xylim = c(1,1.2*max(local.mean))
xlab="Mean Degree"
ylab= expression(frac((Local~Mean~+~Global~Mean),2))
ylab= expression(frac((Ego-based~Mean~+~Alter-based~Mean),2))
ldata = cbind(degree.mean,local.mean)
mtitle="Local Strategy \nLocal Mean and Mean Degree"
plotdata = cbind(degree.mean,(local.mean+global.mean)/2)
mtitle="Local and Global Mean"
mtitle="Ego-based and Alter-based Mean"
pdf("F3_LeverageinRealNetworksA.pdf",family="Times")
# par(mfrow=c(1,2))
par(oma=c(0,1,2,0))
par(mar=c(4,6,2,0))
p.xy = plot(plotdata,bty="n",log="xy",main=mtitle,cex.main=1.25,cex=lcex,col=lcol,pch=lpch,xlim=xylim,ylim=xylim,xlab=xlab,ylab=ylab)
points(plotdata,pch=3,col="white",cex=0.5)
text(plotdata,netlabels,cex=0.6,pos=1,font=2)
xvalues = c(1:500)
kset = c(1,2,5,10,100)
nkset = length(kset)
linecol= rep(NA,nkset)
for(ik in 1:length(kset)){
  k = kset[ik]
  linecol[ik] = gray(0.1+0.5*(1-k/max(kset)))
  yvalues = k*xvalues
  lines(xvalues,yvalues,lty=ik,col=linecol[ik])
}
# leg.txt = paste("leverage: ",kset)
leg.txt = paste(kset)
legend("bottomright",bty="n",title="Iso-Leverage Line",legend=leg.txt,lty=1:nkset,col=linecol,cex=1)#,col=col_cat)
# legend("bottomright",bty="n",title="Iso-Leverage Line",legend=leg.txt,lty=1:nkset,col=linecol,cex=0.75)#,col=col_cat)
grid(nx=10,ny=10,col=gray(0.80))
mtext("A",outer=FALSE,cex=1.5,font=2,adj=-0.0)
dev.off()
#legend("bottomright",bty="n",leg.txt,fill=col_cat,cex=0.75)#,col=col_cat)
#text(xylim[1],0.95*xylim[2],bty="n","A",cex=2,pos=3)#,col=col_cat)
# mtext("B",outer=FALSE,cex=1.5,font=2,adj=0)

########################################
# (b) Global Means versus Mean Degrees
########################################
# xylim = c(1,1.2*max(global.mean))
# xlab="Mean Degree - Random Sampling"
# ylab="Mean Degree - Global FoF"
# gdata = cbind(degree.mean,global.mean)
# mtitle="Global Strategy \nGlobal Mean and Mean Degree"
# #par(oma=c(0,1,2,0))
# par(mar=c(2,2,2,0))
# plot(gdata,bty="n",log="xy",main=mtitle,cex=lcex,col=lcol,pch=lpch,xlim=xylim,ylim=xylim,xlab=xlab,ylab=ylab)
# text(gdata,netlabels,cex=0.5,pos=1,font=2)
# points(gdata,pch=3,col="white",cex=0.5)
# abline(a=1e-100,b=1,lty=2,col="gray")
# legend("bottomright",bty="n",leg.txt,fill=col_cat,cex=0.75)#,col=col_cat)
# 
# grid()
# mtext("C",outer=FALSE,cex=1.5,font=2,adj=0)

########################################
# (c) Local versus Global Means
########################################
lratios = local.mean / global.mean;
nnetworks=length(lratios)
xlim=c(0,8); ylim=c(0,6.2);
xlab="";ylab=""
mx = mean(xlim); my = mean(ylim);
clist = get_all_coordinates(mx,my,lratios)
#csize = runif(nnetworks,1,4)
n_4 = ceiling(nnetworks/4)
# col_cat = c("#4682B4","red","darkgreen","darkorange","#9ACD32","#4B0082")
# ccol = NULL;
# for (k in 1:length(num_in_category))
#   ccol = c(ccol,rep(col_cat[k],num_in_category[k]))
reversename = rep(FALSE,nnetworks)
reversename[(ceiling(nnetworks/4)+1):floor(3*nnetworks/4)]=TRUE
positionvector=rep(4,nnetworks)
positionvector[reversename==TRUE]=2
# par(mar=1*c(1,1,1,0))

pdf("F3_LeverageinRealNetworksB.pdf",family="Times")
par(mar=c(2,0,2,0))
mtitle="Ratio of Local to Global"
mtitle="Ratio of Ego-based to Alter-based Mean"
#par(oma=c(0,1,2,0))
plot(0,type='n',bty="n",xaxt="n",yaxt="n",main=mtitle,cex.main=1.25,xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab,axes=FALSE,asp=1)
draw.circle(x=mx,y=my,radius=2,col="lightgray",border="darkgray")
draw.circle(x=mx,y=my,radius=1,lty=3,col=gray(0.5),lwd=1.5,border="black")
draw.circle(x=mx,y=my,radius=1/2,lty=3,col="black",border="darkgray")
for (k in 1:nnetworks)
{
  add_network_character(clist$i[k,],lcex[k]/2,ccol=lcol[k])
  add_network_character(clist$i[k,],csize=0.25,cpch=3,ccol="white")
  #text(clist$i[k,],netlabels[k],cex=0.5,pos=3)
  
  segments(clist$i[k,1],clist$i[k,2],clist$o[k,1],clist$o[k,2],lty=3,col="white")
  nname = paste("  (",netlabels[k],") ",net_list$names[k],"  ",sep="")
  anglek = 360*(k-1)/(nnetworks) + 180*reversename[k]
  text(clist$o[k,1],clist$o[k,2],nname,cex=0.75,pos=positionvector[k],col=ccol[k],srt=anglek)#,col="darkgray")
}
leg.txt = paste0("(",catshort,") ",categories)
legend("topleft",bty="n",leg.txt,fill=col_cat,cex=0.75)#,col=col_cat)
# text(xylim[1],0.95*xylim[2],bty="n","A",cex=2,pos=3)#,col=col_cat)
mtext("B",outer=FALSE,cex=1.5,font=2,adj=0.00)
dev.off()

###########################################################################################
# Figure F7: rho versus degree variance terms
###########################################################################################
ymax=10
xlim=c(-0.5,+0.5)
ymin=0.1
ylim=c(ymin,ymax)
xlab="Inversity"
ylab="Degree Distribution Multiplier"
xstats=data.frame(nstats_all)
nnetworks = length(xstats[[1]])
mtitle=""#"Source of Differences between Local and Global Means"
crun = ceiling(as.numeric(Sys.time()))
pdf_fname = paste("S12_Source_of_Differences_",length(net_all),"_Networks_",crun,".pdf",sep="")
pdf(pdf_fname,family="Times")
plot(1e-100,type='n',log="y",yaxt="n",main=mtitle,xlab=xlab,ylab=ylab,ylim=ylim,xlim=xlim)
abline(v=0,lty=2,lwd=2)

for(k in 1:nnetworks)
{
  pcex=log(xstats$nnodes[[k]])/2
  x = xstats$rho_d[[k]]
  #y = abs((xstats$mu_l[[k]]-xstats$mu_g[[k]])/xstats$rho_d[[k]])
  y = abs(((xstats$mu_l[[k]]/xstats$mu_g[[k]])-1)/xstats$rho_d[[k]])
  points(x=x,y=y,pch=16,col=lcol[k],cex=pcex)
  points(x=x,y=y,pch=3,col="white")
  text(x=x,y=y,netlabels[k],cex=0.5,pos=1)
}

lseq = function(from=1, to=100000, length.out=6) {
  # logarithmic spaced sequence
  # blatantly stolen from library("emdbook"), because need only this
  exp(seq(log(from), log(to), length.out = length.out))
}
n.out = 20
xrange = seq(-0.5,+0.5,length.out=100)
yrange = lseq(0.01,2000,n.out)#seq(0.01,40,length.out=10)

for (k in 1:n.out)
{
  yisocurve = abs(yrange[k]/xrange)
  lines(xrange,yisocurve,lty=2,lwd=0.75,col="darkgray")
}
#marks = seq(ymin,ymax,length.out=10)
#marks = c(0.1,1,2,5,10,100,500,1000,2000)#lseq(ymin,ymax,length.out=11)
marks = c(0.1,0.25,0.5,1,2,5,10)
axis(2,at=marks,labels=marks)
pexpr = expression("Increasing " ~ frac(mu[L],mu[G]))
text(-0.4,5,pexpr,pos=1,cex=0.8)
xseq = seq(0,-0.4,length.out=20)
seqfn = function(x){return(0.1+1000*(abs(x)^6))}
curve(seqfn,0,-0.38,lwd=0.5,lty=2,col="blue",add=TRUE)
curve(seqfn,0,0.38,lwd=0.5,lty=2,col="blue",add=TRUE)
text(0.4,5,pexpr,pos=1,cex=0.8)
dev.off()
###########################################################################################
# Figure F4: Comparison Across Networks
# (a) Local mean versus Random Sample
# (b) Global mean versus Random Sample
# (e) Local versus Global Means
###########################################################################################

###########################################################################################
# Figure F6: ER, WS, BA Graphs
###########################################################################################

###########################################################
# (B1) Generating Erdos Renyi Networks
###########################################################
N=500; Ngraphs=100;
p_set = c(0.0001,seq(0.05,1,0.05))

er_graph_generator = function(N,p)
{
  imat = matrix(runif(N*N),N,N)
  imat = lower.tri.remove(imat,1)
  imat[imat > p] = 0
  return(sign(imat))
}

graphs_all = list()
graph_all_counter=1
bern_graph = list(); bern_stats = data.frame()
for (ip in 1:length(p_set)){
  bern_graph[[ip]]=list();
  for (k in 1:Ngraphs){
    #dnetwork = as.matrix(rgraph(N, mode="graph", tprob=p_set[ip]))
    #dnetwork = as.matrix(get.adjacency(erdos.renyi.game(N,p_set[ip])))
    dnetwork = er_graph_generator(N,p_set[ip])
    #bern_graph[[ip]][[k]] = dnetwork
    graphs_all[[graph_all_counter]] = dnetwork
    graph_all_counter=graph_all_counter+1
    bern_stats = rbind(bern_stats,cbind(as.data.frame(get_network_stats(dnetwork)),density=p_set[ip]))
    print(paste(date(), " Bern completed: density= ",p_set[ip],", k:",k,sep=""))
  }
}
slist=c('bern_stats','p_set','Ngraphs')
# save(list=slist,file="FoF_Bern_graphs3.image")
# load("data/FoF_Bern_graphs3.image")

leverage_type="l"
# nf = layout(matrix(rbind(c(1,2),c(3,4)),2,2))
# layout.show(nf)


###########################################################
# (B2) Generating Barabasi Albert Scale Free Networks
###########################################################

N=2000; Ngraphs=100;
dmin=1;dmax=floor(N/2);
ba_graph = list(); ba_stats = data.frame()
dmin_set = c(1,2,5,10)
dmax_set = c(100,500,1000)
cgamma_set = seq(1,6,0.25)

dmin_set = 1:10
dmax_set = 500:1000

# Within gamma, choose different dmin and dmax
#source("inversity_functions.R")
for(icgamma in 1:length(cgamma_set))
{
  cgamma = cgamma_set[icgamma]
  ba_graph[[icgamma]]=list()
  for (k in 1:Ngraphs){
    dmin = sample(dmin_set,size=1,replace=FALSE)
    dmax = sample(dmax_set,size=1,replace=FALSE)
    dnetwork=get_ba_sample(dmin:dmax,cgamma,N)
    #ba_graph[[icgamma]][[k]] = dnetwork
    graphs_all[[graph_all_counter]] = dnetwork
    graph_all_counter=graph_all_counter+1
    ba_stats = rbind(ba_stats,cbind(as.data.frame(get_network_stats(dnetwork)),gamma=cgamma,dmin=dmin,dmax=dmax))
    print(paste(date(), " BA completed: gamma= ",cgamma,", dmin: ",dmin,", dmax:", dmax,", k: ",k,sep=""))
  }
}

slist=c('ba_stats','dmin_set','dmax_set','cgamma_set')
#save(list=slist,file="FoF_BAgraphs_DegSeq2.image")

# Do plot of distribution of leverage.
#xlim=c(1,max(ba_stats$mu_d));



###########################################################
# (B3) Generating Watts Strogatz Small World Networks
###########################################################
N=2000; Ngraphs=100;
ws_stats = data.frame()
ws_graph = list();
rprob = runif(Ngraphs)
neighbors = ceiling(runif(Ngraphs,1,20));

neighbor_set = c(1,2,5,10,20)
p_rewire_set = seq(0.05,0.95,0.1)

for(ineighbors in 1:length(neighbor_set))
{
  cneighbors = neighbor_set[ineighbors]
  #ws_graph[[ineighbors]]=list()
  for (ip in 1:length(p_rewire_set))
  {
    cp = p_rewire_set[ip]
    #ws_graph[[ineighbors]][[ip]]=list()
    for (k in 1:Ngraphs)
    {
      tnetwork = as_adjacency_matrix(sample_smallworld(dim=1,size=N,nei=cneighbors,p=cp))
      #ws_graph[[ineighbors]][[ip]] = tnetwork
      graphs_all[[graph_all_counter]] = tnetwork
      graph_all_counter=graph_all_counter+1
      ws_stats = rbind(ws_stats,cbind(as.data.frame(get_network_stats(tnetwork)),ineighbors=ineighbors,ip=ip,neighbors=cneighbors,p.rewire=cp))
      print(paste(date(), " WS completed: n:",cneighbors,", prob:",cp,", k:",k,sep=""))
    }
  }
}

slist=c("bern_stats","ba_stats","ws_stats",'neighbor_set','p_rewire_set','cgamma_set','p_set')
fname="data/FoF_ER_BA_WS_graphs.image"
save(list=slist,file=fname)
#rm(list=ls())
#load("data/FoF_WS_graphs3.image")
#ws_stats





########################################################################################
# Load up the networks
########################################################################################

N=500; Ngraphs=100;
p_set = c(0.0001,seq(0.05,1,0.05))
cgamma_set = seq(1,6,0.25)
dmin_set = 1:10
dmax_set = 500:1000
neighbor_set = c(1,2,5,10,20)
p_rewire_set = seq(0.05,0.95,0.1)

fname="data/FoF_ER_BA_WS_graphs.image"
load(fname)


########################################################################################
#(A) Plot of ER Networks
########################################################################################

crun = ceiling(as.numeric(Sys.time()))
pdfcurrent = "F5abc_"
pdf_fname = paste(pdfcurrent,"_",length(net_all),"_Networks_",crun,".pdf",sep="")
par(oma=c(0,5,2,0))
par(mar=0.1*c(1,1,1,0))
pdf(pdf_fname,family="Times")
nf = layout(matrix(rbind(c(1,2),c(3,4)),2,2))

select_iprob_set = seq(1,length(p_set),1); #1:length(cgamma_set)#c(1,5,7,9,11,13,17,21)
select_iprob_set = c(1,2,3,5,7,9,11,15,19)#,13,15,17,19)
n_size = length(select_iprob_set)
ylim=c(1,3000); xlim=c(1,1.15)#*max(bern_stats$mu_l/bern_stats$mu_d))
xlab="Leverage";ylab="Probability of Edge"; mtitle="Erdos-Renyi"
cex_axis = 0.8

plot(1,type="n",log="xy",main=mtitle,xlim=xlim,ylim=ylim,axes=FALSE,xlab=xlab,ylab=ylab,cex.lab=cex_axis)  
iprob_set = 1:length(p_set)
p_dist=0

cramp = colorRampPalette(colors=c("darkgreen","yellow","red","red","blue"),space = "Lab")
lplotcol = cramp(length(select_iprob_set))
lplotcol = paste0(lplotcol,"7f")
ctol=1e-10
leverage = matrix(0,length(p_set),2)
for(p_index in 1:length(select_iprob_set))
{
  iprob = select_iprob_set[p_index]
  cp = p_set[iprob]
  select = which(bern_stats$density==cp)
  xdata=cbind(bern_stats[select,]$mu_d, 
              bern_stats[select,]$mu_l/bern_stats[select,]$mu_d, 
              bern_stats[select,]$mu_g/bern_stats[select,]$mu_d)
  leverage[iprob,]=c(mean(xdata[,2]),mean(xdata[,3]))
  lpoints = density(xdata[,2],from=1)
  gpoints = density(xdata[,3],from=1)
  zpoints = lpoints
  if(leverage_type=="g"){zpoints = gpoints}
  z.x = zpoints$x
  z.y = zpoints$y+2*(p_index-1)
  lines(z.x,z.y,pch=1,cex=0.25,col=lplotcol[p_index],lwd=2)
  poly.xy = cbind(c(1,z.x),c(2*(p_index-1)+ctol,z.y))
  polygon(poly.xy,col=lplotcol[p_index],border=NA)
  segments(x0=xlim[1],y0=(p_index-1),x1=xlim[2],y1=(p_index),col="lightgray",lwd=0.5,lty=4)
}
lindex = ifelse(leverage_type=="l",1,2)
lines(x=leverage[select_iprob_set,lindex],y=p_set[select_iprob_set]+exp(select_iprob_set-1),type="b",lty=1,lwd=1.5,col="black")
segments(x0=xlim[1],y0=1,x1=xlim[2],y1=1,col="lightgray",lwd = 0.5,lty=4)
for (k in 1:5){segments(x0=1+(k/100),y0=ylim[1],x1=1+(k/100),y1=exp(p_index-1),lwd=0.5,lty=3)}
axis(side=2,at=exp(select_iprob_set-1),labels=p_set[select_iprob_set],cex.axis=cex_axis,las=2)
axis(side=1,at=seq(1,1.15,5/100),labels=seq(1,1.15,5/100),cex.axis=cex_axis,las=0)
mtext("A",outer=FALSE,cex=1.5,font=2,adj=0)



########################################################################################
#(B) Plot of Scale Free Networks
########################################################################################

select_igamma = seq(1,length(cgamma_set),2);
ngamma=length(cgamma_set)

NTYPE="F5b_ScaleFree"
leverage_type=""
pdfcurrent = paste(NTYPE,"_",leverage_type,"Leverage_Neighbors",sep="")
xlab=paste(leverage_type,"Leverage")
crun = ceiling(as.numeric(Sys.time()))
pdf_fname = paste(pdfcurrent,"_",length(net_all),"_Networks_",crun,".pdf",sep="")
# pdf(pdf_fname,family="Times")
# par(oma=c(5,3,2,0))
xlim=c(1,50);#max(ba_stats$mu_l/ba_stats$mu_d))
ylim = c(0,2*(ngamma))
leverage = matrix(0,ngamma,2)
main="Scale Free"; xlab=xlab;ylab=expression(gamma);
ltext = cgamma_set;

cramp = colorRampPalette(colors=c("darkgreen","yellow","orange","blue"),space = "Lab")
lplotcol = cramp(length(cgamma_set))
lplotcol = paste0(lplotcol,"7f")

plot(1,type="n",log="x",main=main,xlim=xlim,ylim=ylim,axes=FALSE,xlab=xlab,ylab=ylab,mgp=c(2,0,0),cex.lab=cex_axis) 
for (icgamma in select_igamma)
{
  p_index = icgamma
  cgamma = cgamma_set[icgamma]
  select = which(ba_stats$gamma==cgamma)
  xdata=cbind(ba_stats[select,]$mu_d, 
              ba_stats[select,]$mu_l/ba_stats[select,]$mu_d, 
              ba_stats[select,]$mu_g/ba_stats[select,]$mu_d)
  leverage[icgamma,]=c(mean(xdata[,2]),mean(xdata[,3]))
  lpoints = density(xdata[,2],from=1)
  zpoints=lpoints;
  lines(zpoints$x,zpoints$y+p_index,lwd=1.5,lty=1,cex=0.25,col=plotcol[p_index])
  polygon(c(1,zpoints$x),c(p_index,zpoints$y+p_index),col=lplotcol[p_index],border=NA)
  segments(x0=xlim[1],y0=p_index,x1=xlim[2],y1=p_index,col="gray",lwd = 0.5,lty=4)
}
grid(nx=20,ny=10,col=gray(0.9))
for (k in 1:40){segments(x0=k,y0=0,x1=k,y1=18,lwd=0.1,lty=3,col="lightgray")}
lines(x=leverage[select_igamma,1],y=select_igamma,type="b",lty=1,col=gray(0.25))
axis(side=2,at=select_igamma,labels=cgamma_set[select_igamma],cex.axis=cex_axis,las=2)
axis(side=1,at=c(1:10,50),labels=c(1:10,50),cex.axis=cex_axis,las=0)
mtext("B",outer=FALSE,cex=1.5,font=2,adj=0)
#dev.off()



########################################################################################
#(C) (1) Plot of Small World Networks with respect to neighbors
########################################################################################

leverage_local = matrix(0,length(neighbor_set),length(p_rewire_set));
leverage_global = leverage_local
for(ineighbors in 1:length(neighbor_set))
{
  cneighbors = neighbor_set[ineighbors]
  for (ip in 1:length(p_rewire_set))
  {
    select = which((ws_stats$neighbors==cneighbors)&(ws_stats$p.rewire==p_rewire_set[ip]))
    xdata=cbind(ws_stats[select,]$mu_d,ws_stats[select,]$mu_l/ws_stats[select,]$mu_d,ws_stats[select,]$mu_g/ws_stats[select,]$mu_d)
    xmeans = apply(xdata,2,mean)
    leverage_local[ineighbors,ip] = xmeans[2]; leverage_global[ineighbors,ip] = xmeans[3]
  }
}

NTYPE="F5c_SmallWorld"
leverage_type=""
xlab=paste(leverage_type,"Leverage")

pdfcurrent = paste(NTYPE,"_",leverage_type,"Leverage_Neighbors",sep="")
xlab=paste(leverage_type,"Leverage")
crun = ceiling(as.numeric(Sys.time()))
pdf_fname = paste(pdfcurrent,"_",length(net_all),"_Networks_",crun,".pdf",sep="")
nneighbors = length(neighbor_set)
leverage = matrix(0,nneighbors,2)
iselect_neighbors = 1:length(neighbor_set)
main="Small World"; xlab=xlab;ylab="Neighbors"
xlim=c(1,1.5);ylim = c(1,20*(nneighbors))
ltext = neighbor_set;
plot_color_transparency = 0.4
plotcol = adjustcolor((iselect_neighbors+1),alpha.f=plot_color_transparency)

cramp = colorRampPalette(colors=c("darkgreen","yellow","orange","blue"),space = "Lab")
lplotcol = cramp(length(neighbor_set))
lplotcol = paste0(lplotcol,"7f")

plot(1,type="n",log="xy",main=main,xlim=xlim,ylim=ylim,axes=FALSE,xlab=xlab,ylab=ylab,mgp=c(2,0,0),cex.lab=cex_axis)
for (ineighbors in iselect_neighbors)
{
  cneighbors = neighbor_set[ineighbors]
  select = which(ws_stats$neighbors==cneighbors)
  xdata=cbind(ws_stats[select,]$mu_d, 
              ws_stats[select,]$mu_l/ws_stats[select,]$mu_d, 
              ws_stats[select,]$mu_g/ws_stats[select,]$mu_d)
  leverage[ineighbors,]=c(mean(xdata[,2]),mean(xdata[,3]))
  lpoints = density(xdata[,2],from=1);gpoints = density(xdata[,3],from=1)
  zpoints=lpoints;
  lines(zpoints$x,zpoints$y+ineighbors,lwd=1.5,lty=1,cex=0.25,col=plotcol[ineighbors])
  polygon(c(1,zpoints$x),c(ineighbors,zpoints$y+ineighbors),col=lplotcol[ineighbors],border=NA)
  segments(x0=xlim[1],y0=ineighbors,x1=xlim[2],y1=ineighbors,col="black",lwd = 0.5,lty=4)
}
for (k in 1:nneighbors){segments(x0=k,y0=0,x1=k,y1=8,lwd=0.5,lty=3)}
lines(x=leverage[,1],y=iselect_neighbors,type="b",lty=1,col="black")

axis(side=2,at=iselect_neighbors,labels=neighbor_set[iselect_neighbors],cex.axis=cex_axis,las=2)
axis(side=1,at=seq(1,2,0.25),labels=seq(1,2,0.25),cex.axis=cex_axis,las=0)
mtext("C",outer=FALSE,cex=1.5,font=2,adj=0)


#########################################################################
# (C) (2): Plot Small World Leverage with respect to rewiring probabilities
#########################################################################

NTYPE="F5d_SmallWorld"; leverage_type=""
pdfcurrent = paste(NTYPE,"_",leverage_type,"Leverage_Rewire",sep="")
xlab=paste(leverage_type,"Leverage")
crun = ceiling(as.numeric(Sys.time()))
pdf_fname = paste(pdfcurrent,"_",length(net_all),"_Networks_",crun,".pdf",sep="")
p_rewire_set = sort(unique(ws_stats$p.rewire))
iselect_rewire = c(1,3,5,7,10)
n_rewire = length(iselect_rewire)
leverage = matrix(0,length(p_rewire_set),2)

main="Small World"; xlab=xlab;ylab="Rewire Probability"
xlim=c(1,1.5);ylim = c(1,30*(n_rewire))
ltext = p_rewire_set[iselect_rewire]

cramp = colorRampPalette(colors=c("darkgreen","yellow","orange","blue"),space = "Lab")
lplotcol = cramp(length(p_rewire_set))
lplotcol = paste0(lplotcol,"7f")

plot(1,type="n",log="xy",main=main,xlim=xlim,ylim=ylim,axes=FALSE,xlab=xlab,ylab=ylab,cex.lab=cex_axis)
for (irewire in iselect_rewire)
{
  cp_rewire = p_rewire_set[irewire]
  select = which(ws_stats$p.rewire==cp_rewire)
  xdata=cbind(ws_stats[select,]$mu_d, 
              ws_stats[select,]$mu_l/ws_stats[select,]$mu_d, 
              ws_stats[select,]$mu_g/ws_stats[select,]$mu_d)
  leverage[irewire,]=c(mean(xdata[,2]),mean(xdata[,3]))
  lpoints = density(xdata[,2],from=1);gpoints = density(xdata[,3],from=1)
  zpoints=lpoints;
  lines(zpoints$x,zpoints$y+irewire,lwd=1.5,lty=1,cex=0.25,col=lplotcol[irewire])
  polygon(c(1,zpoints$x),c(irewire,zpoints$y+irewire),col=lplotcol[irewire],border=NA)
  segments(x0=xlim[1],y0=irewire,x1=xlim[2],y1=irewire,col="black",lwd = 0.5,lty=4)
}

for (k in iselect_rewire){segments(x0=k,y0=0,x1=k,y1=8,lwd=0.5,lty=3)}
lines(x=leverage[iselect_rewire,1],y=iselect_rewire,type="b",lty=1,col="black")
axis(side=1,at=seq(1,2,0.25),labels=seq(1,2,0.25),cex.axis=cex_axis,las=0)
axis(side=2,at=iselect_rewire,labels=p_rewire_set[iselect_rewire],cex.axis=cex_axis,las=2)
mtext("D",outer=FALSE,cex=1.5,font=2,adj=0)
dev.off();



###########################################################################################
# Figure F8: Assortativity and Inversity Plots
###########################################################################################
n_sample = 2000;
rho_a_thresh = 0.4
ba_select = which(abs(ba_stats$rho_a)<rho_a_thresh)
ba_select = sample(ba_select,size=n_sample,replace=FALSE)
ws_select = sample(1:dim(ws_stats)[1],size=n_sample,replace=ifelse(dim(ws_stats)[1]>n_sample,FALSE,TRUE))
bern_select = sample(1:dim(bern_stats)[1],size=n_sample,replace=ifelse(dim(bern_stats)[1]>n_sample,FALSE,TRUE))

gstats = rbind(ba_stats[ba_select,c(1:9)],ws_stats[ws_select,c(1:9)],bern_stats[bern_select,c(1:9)])
plot_pch=c(rep(1,dim(ba_stats[ba_select,])[1]),rep(2,dim(ws_stats[ws_select,])[1]),rep(3,dim(bern_stats[bern_select,])[1]))

xlim=c(-0.2,+0.2); ylim=xlim;# ylim=c(-0.1,+0.1)
#plot_pch = rep(2,dim(gstats)[1]) #c(rep(3,dim(gstats)[1]),rep(2,dim(ba_stats)[1]))
#select = which((gstats$rho_a*gstats$rho_d)>0)
select_positive = which((gstats$rho_a > 0)&(gstats$rho_d>0))
select_negative = which((gstats$rho_a < 0)&(gstats$rho_d<0))

plot_col = rep("lightgray",dim(gstats)[1]); 
#plot_col[select]="red"
plot_col[select_positive]="blue"
plot_col[select_negative]="red"
xlab="Assortativity"
ylab="Inversity"
main=""#  "Inversity and Assortativity"

# c(gstats$rho_a[select[100]],gstats$rho_d[[select[100]]])

pdfmain=""
pdfcurrent="AssortativityInversity"
crun = ceiling(as.numeric(Sys.time()))
pdffile=paste(pdfmain,pdfcurrent,"_",crun,".pdf",sep="")
pdf(pdffile,family="Times")
pcol=c("darkgreen","darkgray","darkorange")
plot(gstats$rho_a,gstats$rho_d,add=TRUE,bty="n",cex=0.4,pch=plot_pch,col=plot_col,xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab,main=main)
Plot_ConvexHull(bern_stats$rho_a,bern_stats$rho_d,pcol[1])
Plot_ConvexHull(ba_stats$rho_a,ba_stats$rho_d,pcol[2])
Plot_ConvexHull(ws_stats$rho_a,ws_stats$rho_d,pcol[3])
abline(a=0,b=0,xlim=xlim,ylim=ylim,col="black",lty=3)
abline(a=0,b=1e100,xlim=xlim,ylim=ylim,col="black",lty=3)
legend(x="bottomright",inset=0,bty="n",xjust=0.0,#text.width=0.5,
       legend = c("Scale-Free (BA)","Small World (WS)","Bernoulli"),
       col="black",pch=c(1,2,3),lwd=1,cex=1,horiz=FALSE)
# legend(x="topright",inset=c(.1,.1),bty="n",xjust=0.5,#text.width=0.5,
#        legend = c("(+) Inversity and Assortativity"),
#        text.col="blue",cex=0.75)
# legend(x="bottomleft",inset=c(.1,.1),bty="n",xjust=0.5,#text.width=0.5,
#        legend = c("(-) Inversity and Assortativity"),
#        text.col="red",cex=0.75)
dev.off()

pdfcurrent="AssortativityInversity_ZoomView"
crun = ceiling(as.numeric(Sys.time()))
pdffile=paste(pdfmain,pdfcurrent,"_",crun,".pdf",sep="")
pdf(pdffile,family="Times")
pcol=c("darkgreen","darkgray","darkorange")
xlim=c(-0.02,+0.02); ylim=xlim;# ylim=c(-0.1,+0.1)
xlab="";ylab=""
plot(gstats$rho_a,gstats$rho_d,bty="n",cex=0.8,pch=plot_pch,col=plot_col,xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab,main=main)
abline(a=0,b=0,xlim=xlim,ylim=ylim,col="black",lty=3)
abline(a=0,b=1e100,xlim=xlim,ylim=ylim,col="black",lty=3)
dev.off()

# Pick a random network which is colored RED. Get its coordinates.
# From there, connect to outside the region, plot the structure.
# (x,y) => (x+5,y+5) # gplot(xnet,coord=(x+5,y+5))
#segment()
# (x,y) => (x+5,y+5) # gplot(xnet,coord=(x+5,y+5))
# Can also do the curve in powerpoint.

igraph = select_positive[1]
graphs_all[[igraph]]
xnet=graphs_all[[igraph]]
xnet=rgraph(10,tprob=.2)
gplot(xnet,gmode="graph")

get_same_inversity_assortativity = function(wsign,n.size)
{
  check = FALSE
  while (!check){
    tprob = runif(1)/2
    xg = rgraph(n.size,mode="graph",tprob=tprob)
    badgraph = (sd(rowSums(xg))==0) || (sum(rowSums(xg)==0)>0)
    if (!badgraph){
      xstats = get_network_stats(xg)
      signsum = (sign(xstats$rho_d) + sign(xstats$rho_a))/2
      check = (signsum == wsign)
    }
  }
  return(xg)
}
vcol=c("red","blue"); ecol="gray"
par(mfrow=c(1,2))
nsize=7

xg=get_same_inversity_assortativity(1,nsize)
pdf("example_positive_inversity_assortativity2.pdf",family="Times")
vplot_network_local(xg,gmode="graph",max.color="blue")
dev.off()

xg=get_same_inversity_assortativity(-1,nsize)
pdf("example_negative_inversity_assortativity.pdf",family="Times")
vplot_network_local(xg,gmode="graph",max.color="red")
dev.off()

g1=gplot(xg,gmode="graph",vertex.col=vcol[2],edge.col=ecol)
xg=get_same_inversity_assortativity(-1,nsize)
g2=gplot(xg,xlim=c(0,0.1),gmode="graph",vertex.col=vcol[1],edge.col=ecol)
gplot(xg,coord=g2*10,gmode="graph",vertex.col=vcol[1],edge.col=ecol)

gplot(xg,coord=g2/10,gmode="graph",vertex.col=vcol[1],edge.col=ecol)
xlim=c(-0.2,+0.2); ylim=xlim;# ylim=c(-0.1,+0.1)
plot(0,bty="n",cex=0.4,pch=plot_pch,col=plot_col,xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab,main=main)
points(g1+c(1,1),pch=3,col="black")


###########################################################################################
# Figure F10: Strategy illustration for RS, LS and GS strategies
###########################################################################################
nnodes = 20
xnet = rgraph(nnodes,tprob=0.15,mode="graph")
#gplot(xnet,gmode="graph")
glabels=letters[1:nnodes]
pcol=c("black","red","darkgreen")
inodes = sample(1:nnodes,size=3,replace=FALSE)
vcol=rep("lightgray",nnodes)
ecol=matrix("lightgray",nnodes,nnodes)
vcol[inodes]=pcol[1]

rsample = sample(which(xnet[inodes[2],]==1),size=1)
vcol[rsample]=pcol[2]
ecol[inodes[2],rsample]="red"
ecol[rsample,inodes[2]]="red"
gsample = which(xnet[inodes[3],]==1)
vcol[gsample]=pcol[3]
ecol[inodes[3],] = pcol[3]; ecol[,inodes[3]] = pcol[3]
#g1=gplot(xnet,gmode="graph",mode="kamadakawai",label=glabels,vertex.cex=2,vertex.col=vcol,edge.col=ecol)
pdf("F10_Strategies2.pdf")
gplot(xnet,coord=g1,gmode="graph",label.cex=1.5,label=glabels,vertex.cex=3,vertex.col=vcol,edge.col=ecol)
dev.off()


cat_names = c("A","FS","OS","C","I","B")
num_in_category = c(4,4,4,3,3,4)
netlabels = generate_network_labels(cat_names,num_in_category)
length(netlabels)


###########################################################################################
# Figure F11: Kendall Rau
###########################################################################################
xnet = rgraph(20,tprob=0.15,mode="graph")



###########################################################################################
###########################################################################################
# xnet = net_list$a[[2]]
# nnodes = dim(xnet)[1]
# xrank = as.matrix(get_rankdiff_matrix(xnet))
# ddegrees = rowSums(xnet)
# xrankdeciles = quantile(xrank[which(xrank>0,arr.ind=TRUE)],probs=seq(0,1,0.1))
# i10=2; i90=10;
# 
# edge.col= matrix("lightgray",dim(xnet)[1],dim(xnet)[1]);#gray(0.1)#gray(x2)
# edge.lwd= matrix(1,dim(xnet)[1],dim(xnet)[1]);#gray(0.1)#gray(x2)
# edge.lty= matrix(1,dim(xnet)[1],dim(xnet)[1]);#gray(0.1)#gray(x2)
# #edge.col[x2<=x2deciles[2]]="white";
# 
# xselect_high = which((xrank>0)&(xrank>=xrankdeciles[i90]),arr.ind=TRUE)
# edge.col[xselect_high]="red";
# edge.lwd[xselect_high]=4
# edge.lty[xselect_high]=1
# 
# xselect_low = which((xrank>0)&(xrank<=xrankdeciles[i10]),arr.ind=TRUE)
# edge.col[xselect_low]="blue";
# edge.lwd[xselect_low]=4
# edge.lty[xselect_low]=2
# v.cex=1.5*log(ddegrees)
# v.col=rep("red",nnodes)
# select_high_nodes = ddegrees > median(ddegrees)
# v.col[select_high_nodes]="black"
# gplot(as.matrix(xrank),gmode="graph",label=1:nnodes,edge.lwd=edge.lwd,edge.col=edge.col,edge.lty=edge.lty,vertex.cex=v.cex,vertex.col=v.col)
###########################################################################################
###########################################################################################


##################################################################
# Figure S13: India Villages Data
##################################################################
get_village_stats = function(flag)
{
  flagtext = c("_HH","")
  #fprefix = paste("./IndiaVillages/Data/1. Network Data/Adjacency Matrices/adj_allVillageRelationships",flagtext[flag+1],"_vilno_",sep="")
  fprefix = paste("~/Dropbox/Data/India_Villages_datav4-2.0/Data/1. Network Data/Adjacency Matrices/adj_allVillageRelationships",flagtext[flag+1],"_vilno_",sep="")
  villages = c(1:77)[-c(13,22)]
  fnames = list(); fdata = list()
  for (k in villages){
    fnames[[k]] = paste(fprefix,k,".csv",sep="")
    tnet = read.csv(fnames[[k]],header=FALSE)
    tnet = sign(tnet + t(tnet))
    fdata[[k]] = tnet
    print(paste("Village k:",k, "completed"))
  }
  villagestats = matrix(0,max(villages),7)
  for(k in villages){
    villagestats[k,] = unlist(get_network_stats(fdata[[k]]))
    print(paste("Village k:",k, "completed"))
  }
  return(villagestats)
}

hh_stats = get_village_stats(0)
ind_stats = get_village_stats(1)
vselect = -c(13,22)
nvillages = dim(hh_stats)[1] - length(vselect)
col1=rgb(0, 0, 1,0.5)
col2=rgb(1, 0, 0,0.5)
pcol = c(rep(col2,nvillages),rep(col1,nvillages))
ppch = c(rep(17,nvillages),rep(16,nvillages))
degree_mean = c(ind_stats[vselect,2],hh_stats[vselect,2])
degree_sd = c(ind_stats[vselect,5],hh_stats[vselect,5])
xlab="Mean Degree"; ylab="Standard Deviation of Degree"
pdf("IndiaVillage_Mean_SD_Degree.pdf",family="Times")
plot(degree_mean,degree_sd,col=pcol,pch=ppch,xlab=xlab,ylab=ylab)
ltext = c("Individual","Household")
lcol = c(col1,col2)
lpch = c(16,17)
legend("topleft",legend=ltext,bty="n",col=lcol,pch=lpch,cex=1.25)
dev.off()

inversities = cbind(ind_stats[,6],hh_stats[,6])
mtitle = "Inversity in India Villages at Individual and Household Level"
mtitle = ""
col1=rgb(0, 0, 1,0.25)
col2=rgb(1, 0, 0,0.25)
pdf("Inversity_India_Villages2.pdf",family="Times")
hist(inversities[,1],xlim=c(-0.5,0.25),breaks=31,col=col1,xlab="Inversity Value",main=mtitle)
hist(inversities[,2],breaks=31,col=col2,xlab="Inversity Value",add=TRUE)
points(density(inversities[,1]), pch=16, col=col1, cex=0.5)
points(density(inversities[,2]), pch=17, col=col2, cex=0.5)
col1=rgb(0, 0, 1,1)
col2=rgb(1, 0, 0,1)
abline(v=median((inversities[,1])),lwd=2,lty=2,col=col1)
abline(v=median((inversities[,2])),lwd=2,lty=2,col=col2)
abline(v=0,lwd=1,lty=1,col="black")
ltext = c("Individual","Household")
lcol = c(col1,col2)
lpch = c(16,17)
legend("top",legend=ltext,bty="n",col=lcol,pch=lpch,cex=1.5)
dev.off()


###########################################################################################
# Plot degree imbalance across edges
###########################################################################################

# generate a bunch of graphs with different degree assortativity
mlayout = rbind(c(1,2),c(3,4),c(5,5))
layout(mlayout)
par(mar=c(1,1,1,0))
for (m in 1:4) plot_ranked_matriximage2(xnet)

cpalettef = colorRampPalette(c("white","yellow","orange","red"))
color.bar(cpalettef(100),0,1,title="Degree Imbalance across Edges")



#gplot(mnet,gmode="graph")


dseq = sample(c(rep(1,4),2:50),size=200,replace=TRUE)
dpermgraph = list()
npermute = 500
for (k in 1:npermute) {
  tdpermgraph = degree.sequence.game(out.deg=dseq,method="vl")
  dpermgraph[[k]] = as.matrix(get.adjacency(tdpermgraph))
}
check_match_networks(dseq,dpermgraph)
sum(check_match_networks(dseq,dpermgraph))


#dseq = sample(c(rep(1,4),2:50),size=200,replace=TRUE)
dsplitpermgraph = list()
npermute = 500
for (k in 1:npermute) {
  tdpermgraph = split_degree_sequence_game(dseq,2)#,method="vl")
  dsplitpermgraph[[k]] = tdpermgraph
}
# check_match_networks(dseq,dsplitpermgraph)
sum(check_match_networks(dseq,dsplitpermgraph))


#get the range of inversities
npermstats = get_all_network_stats(dsplitpermgraph)
npermstats

mu_l = unlist(npermstats[,4])
mu_d = mean(dseq); var_d = var(dseq)
mu_g = mu_d + var_d/mu_d
xlim = c(0.9*mu_d,1.5*max(mu_l))

tpcol=c("darkgray","darkred","darkgreen")
pcol= tpcol
for (k in 1:length(tpcol)) pcol[k]=paste0(get_hex_from_colorname(tpcol[k]),"7f")

dev.off()
plot(density(mu_l),xlim=xlim,col=pcol[2],lwd=2,ann=FALSE)
abline(v=mu_d,col=pcol[1],lty=2,lwd=2)
abline(v=mu_g,col=pcol[3],lty=1,lwd=2)

dsplitpermgraph=list()
for (k in 1:npermute) {
  tdpermgraph = split_degree_sequence_game(dseq,1)#,method="vl")
  dsplitpermgraph[[k]] = tdpermgraph
}
sum(check_match_networks(dseq,dsplitpermgraph))
npermstats = get_all_network_stats(dsplitpermgraph)
mu_l = unlist(npermstats[,4])
lines(density(mu_l),col="orange",lty=2)





imin = which.min(mu_l)
imax = which.max(mu_l)

npermstats[imin,]
npermstats[imax,]

#par(mfrow=c(1,2))
mlayout = rbind(c(1,2),c(3,3))#,c(5,5))
layout(mlayout)
par(mar=c(2,1,0,0))
plot_ranked_matriximage2(dpermgraph[[imin]])
plot_ranked_matriximage2(dpermgraph[[imax]])
cpalettef = colorRampPalette(c("white","yellow","orange","red"))
color.bar(cpalettef(100),0,1,title="Degree Imbalance across Edges")

# Create 1 plot out of this to be placed in supplement.
