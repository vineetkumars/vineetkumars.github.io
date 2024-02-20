processnetwork <- function(net,excludeisolates=TRUE)
{
  diag(net) = 0
  mnet = sign(net + t(net))
  degrees = rowSums(mnet)
  selected = degrees>0
  if (excludeisolates==TRUE){
    rnet = mnet[selected,selected]
  } else {
    rnet = mnet
  }
  return(rnet)
}
var_pop = function(x)
{
  return(mean((x - mean(x))^2))
}
get_kappas = function(d)
{
  kappas = rep(NA,4)
  kappas[1] = mean(1/d)
  kappas[2] = mean(d)
  kappas[3] = mean(d^2)
  kappas[4] = mean(d^3)
  return(kappas)
}

get_network_stats <- function(xnet_)
{
  ddegrees_original = rowSums(processnetwork(xnet_,FALSE))
  xnet = processnetwork(xnet_,TRUE)
  ddegrees = rowSums(xnet)
  nedges = sum(sign(xnet)>0)
  nnodes = dim(xnet)[1]
  x_density = 2*nedges/(nnodes*(nnodes-1))
  d_left = matrix(0,nedges,1)
  d_inv_right = matrix(0,nedges,1)
  d_right = matrix(0,nedges,1)
  xnet_i = which(xnet>0,arr.ind = T)
  # ddegrees = rowSums(xnet)
  d_left = ddegrees[xnet_i[,1]]
  d_right = ddegrees[xnet_i[,2]]
  d_inv_right = 1/ddegrees[xnet_i[,2]]
  mu_l = sum((d_left / d_right) + (d_right / d_left))/(2*nnodes)
  mu_d = mean(ddegrees); sigma_d = sqrt(var_pop(ddegrees))
  mu_g = (1/mu_d)*(mu_d^2 + sigma_d^2)
  rho_a = cor(d_left,d_right)
  rho_d = cor(d_left,d_inv_right)

  k = get_kappas(ddegrees)
  sigmaO = sqrt(k[4]/k[2] - (k[3]/k[2])^2)
  sigmaID = (1/k[2])*(k[1] - 1/k[2])
  
  mu_l_formula = mu_g + rho_d * sqrt((k[1]-1/k[2])*(k[2]*k[4]-k[3]*k[3])/k[2])
  mu_l_proxy = mu_g + (-rho_a) * sqrt((k[1]-1/k[2])*(k[2]*k[4]-k[3]*k[3])/k[2])
  
  # mu_l_formula = sum((d_left / d_right) + (d_right / d_left))/(2*nnodes)
  rlist = list(nnodes=nnodes,density=x_density,mu_d=mean(ddegrees_original),kappa=k,sigmaO=sigmaO,sigmaID=sigmaID,mu_g=mu_g,mu_l=mu_l,mu_l_formula=mu_l_formula,mu_l_proxy=mu_l_proxy,sigma_d=sigma_d,rho_d=rho_d,rho_a=rho_a)
  return(rlist)
}

get_degrees_from_edges <- function(xnet_)
{
  xnet = processnetwork(xnet_,TRUE)
  ddegrees = rowSums(xnet)
  xnet_i = which(xnet>0,arr.ind = T)
  d_left = ddegrees[xnet_i[,1]]
  d_right = ddegrees[xnet_i[,2]]
  rlist = data.table(l=d_left,r=d_right)
  return(rlist)
}

generate_star_network = function(n)
{
  xnet = matrix(0,n+1,n+1)
  xnet[1,] = 1; xnet[,1]=1;
  xnet = sign(xnet + t(xnet))
  diag(xnet) = 0
  return(xnet)
}

generate_complete_network = function(n)
{
  xnet = matrix(1,n+1,n+1)
  diag(xnet) = 0
  return(xnet)
}

generate_giant_network = function(netlist)
{
  nnet = length(netlist)
  netsizes = rep(0,nnet)
  for (i in 1:nnet){ netsizes[i] = dim(netlist[[i]])[1]}
  n_total = sum(netsizes)
  net_all = matrix(0,n_total,n_total)
  istart = c(1,1+cumsum(netsizes))
  iend = istart[-1]-1
  for (i in 1:nnet){
    net_all[istart[i]:iend[i],istart[i]:iend[i]] = netlist[[i]]
  }
  return(net_all)
}

get_cor_ia = function(nstats,select=NULL)
{
  tstats = nstats
  if(sum(select)>0) tstats = nstats[select,]
  ivalues = as.vector(unlist(tstats[,7]))
  avalues = as.vector(unlist(tstats[,8]))
  cor_ia = round(cor(ivalues,avalues),digits=3)
  rlist = list(i=ivalues,a=avalues,cor_ia=cor_ia)
  return(rlist)
}

generate_clique = function(csize)
{
  cmat = matrix(1,csize,csize)
  diag(cmat) = 0
  return(cmat)
}

generate_many_cliques = function(csizes)
{
  ntotal = sum(csizes)
  nnets = length(csizes)
  net_all = list();
  for (i in 1:nnets){
    net_all[[i]] = generate_clique(csizes[i])
  }
  rgiantnet = generate_giant_network(net_all)
  return(rgiantnet)
}

get_node_degree = function(nodeid)
{
  tquery = "select * from twitter where o=:nodeid or d=:nodeid;"
  c_edges = dbGetQuery(con, tquery,params=list(nodeid=nodeid))
  return(dim(c_edges)[1])
}

get_node_records =function(cedge)
{
  d_o = get_node_degree(cedge[1])
  d_d = get_node_degree(cedge[2])
  tdlist = rbind(c(d_o,d_d),c(d_d,d_o))
  return(tdlist)
}

get_inversity_samples = function(nsamples=10)
{
  c_rows = sample(1:nrows,size=nsamples,replace=FALSE)
  dlist = NULL;
  for (i in 1:nsamples)
  {
    c_od = get_edge(c_rows[i])
    tdlist = get_node_records(c_od)
    dlist = rbind(dlist,tdlist)
  }
  return(dlist)
}
get_edge = function(c_rownum)
{
  tquery = "SELECT * FROM twitter WHERE rowid=:crow";
  c_edge = dbGetQuery(con, tquery,params=list(crow=c_rownum))
  c_o = c_edge$o; c_d = c_edge$d
  return(c(c_o,c_d))
}
set.seed(10) # For replicability.

library(igraph)
library(sna)
library(data.table)
library(knitr)

#############################################################################################
# Purpose: Replace Figure S6 from prior revision
#############################################################################################
figuresdir = "~/Documents/2024_PNAS_Submission_R4/Figures/"

networkstats = NULL;
cliquesizes = c(2:20)
starsizeset = seq(2,80,1)
gcliquenet = generate_many_cliques(cliquesizes)

mat_layout <- matrix(c(1,2,2,3,2,2),nrow=2,byrow=TRUE)
cliq_min = min(cliquesizes); cliq_max = max(cliquesizes)

for (cstarsize in starsizeset)
{
  tstarnet = generate_star_network(n=cstarsize-1)
  tnetall = generate_giant_network(list(tstarnet,gcliquenet))
  tnetstats_temp = get_network_stats(tnetall)
  tnetstats = cbind(tnetstats_temp,starsize=cstarsize)
  networkstats = rbind(networkstats,tnetstats)
  
  rho_d = round(tnetstats$rho_d,3); rho_a = round(tnetstats$rho_a,3)
  pname = paste0("Networks_","Starsize_",cstarsize,"_Cliques_",cliq_min,"_to_",cliq_max,".pdf")
  pname = paste0(figuresdir,pname)
  pdf(pname)#,width=3,height=3)
  layout(mat_layout)
  par(mar = c(0.1, 0.1,0.1,0.1), oma=c(1,1,1,1))
  gplot(tstarnet,gmode="graph",edge.col="gray90",vertex.col="blue")
  # replayPlot(gcliqueplot)
  set.seed(1234)
  gplot(gcliquenet,gmode="graph",edge.col="gray90")
  plot.new()
  ltext = paste0("Size of star = ",cstarsize,"\n","\n","Cliques with sizes \n ranging from 2 to 20","\n","\n","Inversity = ",rho_d,"\n","\n","Assortativity = ",rho_a)
  text(0.5,0.5,ltext,cex=2)
  dev.off()
  plot_crop(pname,quiet=TRUE) #Needs "knitr" packages
  print(paste("Completed: ",pname))
}

# for (cstarsize in starsizeset)
# {
#   pname = paste0("Networks_","Starsize_",cstarsize,"_Cliques_",cliq_min,"_to_",cliq_max)#,")#.pdf")
#   pdfjam = paste0("/usr/bin/pdfjam --papersize '{3in,3in}' ",pname,".pdf -o ",pname,"_3in.pdf")
#   system(pdfjam)
#   print(paste("Completed: ",pname))
# }

#############################################################################################
# Purpose: Get Correlation and R^2 Statistics
#############################################################################################

cor_ia = round(cor(networkstats$rho_d,networkstats$rho_a),3)
print(paste("Correlation between inversity and assortativity =",cor_ia))
r2_ia = round(summary(lm(rho_d~rho_a,data=networkstats))$r.squared,3)
print(paste("R2 =",r2_ia))
# Take only the subset of networks where both inversity and assortativity > 0.0
select_a = networkstats$rho_a>0.0
select_i = networkstats$rho_d>0.0
select_subset = select_a & select_i
r = get_cor_ia(networkstats,select_subset)
print(paste("Correlation between Assortativity and Inversity in Select Subset:",round(r$cor_ia,3)))
lm_ia = lm(r$i~r$a)
r_sq_ia = round(summary(lm_ia)$r.squared,3)
print(paste("R^2 for linear model between inversity and assortativity in Select Subset: ",r_sq_ia))


#############################################################################################
# Purpose: Replace Figure S6. Show how Inversity and Assortativity vary with star size
#############################################################################################

xlab="Size of star"
ylab = "Inversity and Assortativity"
ltext=c("Inversity","Assortativity")
atext1 = "Region of same sign \n for assortativity \n and inversity"
atext2 = paste0("Correlation between \n assortativity and inversity\n=\n",cor_ia)
atext2 = paste0("Correlation between \n assortativity and inversity = ",cor_ia)
pname = paste0(figuresdir,"InversityAssortativityExampleB12.pdf")
pdf(pname)
plot(networkstats$starsize,networkstats$rho_d,yaxt="n",type="b",col="blue",cex=0.5,cex.axis=1.1,cex.lab=1.5,xlim=c(0,60),ylim=c(-1,1.5),xlab=xlab,ylab=ylab)
lines(networkstats$starsize,networkstats$rho_a,type="b",col="darkred",cex=0.5,pch=3)
axis(side=2, at=c(-1,-0.5,0,0.5,1))
abline(h=0,lty=2)
abline(h=1,lty=2)
abline(h=-1,lty=2)
text(x=38.25,y=-0.5,atext1,cex=1)
text(x=15,y=1.30,atext2,cex=1)

legend("topright",legend=ltext,bty="n",pch=c(1,3),cex=1.25,col=c("blue","darkred"))
lgray1 = rgb(0.05,0.05,0.05,0.02)
lgray2 = rgb(0.15,0.15,0.15,0.08)
rect(29,-1,47,1,col=lgray2,lty=2)
# rect(0,-1,260,1,col=lgray1,lty=2,density=100)
# rect(260,-1,960,1,col=lgray2,lty=2)
# rect(960,-1,1000,1,col=lgray1,lty=2,density=100)
# rect(57,-1,100,1,col=lgray,lty=2)
dev.off()

#############################################################################################
# Get Twitter Network Statistics from SQLite
# Required local twitter database (large 60 GB+)
#############################################################################################
library(RSQLite)
sqlfile = "/home/uastra/Documents/2024_PNAS_Submission_R4/twitter.db"
con <- dbConnect(SQLite(),dbname=sqlfile)
dbListTables(con)
tquery = "select count(*) from twitter;"
nrows = as.numeric(dbGetQuery(con, tquery))
nrows=1468365182


date();
d = get_inversity_samples(5000)
date()


rho_a = cor(d[,1],d[,2])
rho_d = cor(d[,1],1/d[,2])
rholist[[i]] = c(rho_a,rho_d)
print(paste("Assortativity:",round(rho_a,3)," and Inversity:",round(rho_d,3)))

dbDisconnect(con)



#############################################################################################
# Purpose: Replicate reviewer 1's plots -- Part A (Cliques)
#############################################################################################
figuresdir = "~/Documents/2024_PNAS_Submission_R4/Figures/"
# cliquesizes = c(2:15)
# gcliquenet = generate_many_cliques(cliquesizes)
# netstats = get_network_stats(gcliquenet)
# gplot(gcliquenet,gmode="graph",edge.col="gray80")
# r = get_degrees_from_edges(gcliquenet)
# par(mar = c(2,2,1,4), oma=c(2,2,2,2))
# plot(r$l,r$r,col="blue",type="b",xlim=c(0,max(cliquesizes)),ylim=c(0,1.2*max(cliquesizes)),xlab=expression("d"[i]),ylab=expression("d"[j]),pch=1,cex.lab=1.5)
# atext1 = paste0("Assortativity: ",round(netstats$rho_a,3),"\nInversity: ",round(netstats$rho_d,3))
# text(x=04,y=15,atext1,cex=1.5)
# par(new=TRUE)
# plot(r$l,1/r$r,col="red",type="b",pch=2,xlim=c(0,max(cliquesizes)),ylim=c(0,1),axes=FALSE)#,xlab=expression("d"[i]),ylab=expression("1/d"[j]),pch=1,cex.lab=1.5)
# 
# # points(r$l,1/r$r,col="red",type="b",pch=2)
# mtext(expression("d"[i]),side=1,col="black",line=2,cex=1.5) 
# mtext(expression("d"[j]),side=2,col="blue",line=2,cex=1.5) 
# mtext(expression("1/d"[j]),side=4,col="red",line=3,cex=1.5) 
# axis(4, ylim=c(0,1), col="red",col.axis="red",las=1)

csize = list()
csize[[1]] = c(2:15)
csize[[2]] = c(2:6,8:15)
csize[[3]] = c(2:4,10:15)
csize[[4]] = c(2,12:15)
csize[[5]] = c(2,15)

# mat_layout <- matrix(c(1:10),nrow=5,byrow=TRUE)
# layout(mat_layout)

for(k in 1:5)
{
  cliquesizes = csize[[k]]
  fname = paste0("CliqueSizes_",paste0(csize[[k]],collapse="_"),".pdf")
  gcliquenet = generate_many_cliques(cliquesizes)
  netstats = get_network_stats(gcliquenet)
  round_places=3
  kappa = round(netstats$kappa,round_places)
  sigmaO = round(netstats$sigmaO,round_places)
  sigmaID = round(netstats$sigmaID,round_places)
  r = get_degrees_from_edges(gcliquenet)
  pdf(fname)
  par(oma=c(2,2,2,2))
  
  # par(mfrow=c(1,2))
  mat_layout <- matrix(c(1,1,2,2,2,3,3,2,2,2),nrow=2,byrow=TRUE)
  layout(mat_layout)
  par(mar = c(2,2,1,4))
  par(las=1)
  gplot(gcliquenet,gmode="graph",edge.col="gray70",vertex.col="gray95")
  plot(r$l,r$r,col="blue",type="b",xlim=c(0,max(cliquesizes)),ylim=c(0,1.2*max(cliquesizes)),xlab=expression("d"[i]),ylab=expression("d"[j]),pch=1,cex.lab=1.5)
  par(new=TRUE)
  plot(r$l,1/r$r,col="red",type="b",pch=2,xlim=c(0,max(cliquesizes)),ylim=c(0,1),axes=FALSE)#,xlab=expression("d"[i]),ylab=expression("1/d"[j]),pch=1,cex.lab=1.5)
  # points(r$l,1/r$r,col="red",type="b",pch=2)
  mtext(expression("d"[i]),side=1,col="black",line=2,cex=1.5) 
  mtext(expression("d"[j]),side=2,col="blue",line=2,cex=1.5) 
  mtext(expression("1/d"[j]),side=4,col="red",line=2,cex=1.5) 
  axis(4, ylim=c(0,1), col="red",col.axis="red",las=1)
  ltext=c(expression('Plot of (d'[i]*', d'[j]*')'),expression('Plot of (d'[i]*', 1/d'[j]*')'))
  legend("topright",legend=ltext,bty="n",col=c("blue","red"),pch=c(1,2),cex=1.5)
  atext1 = paste0("Assortativity: ",round(netstats$rho_a,round_places),"\nInversity: ",round(netstats$rho_d,round_places))
  # Get proxy local mean.
  text(x=5,y=12,atext1,cex=1.25,col="gray40")
  plot.new()
  

  atext2 = paste0("Local Mean: ",round(netstats$mu_l,round_places),"\nFormula Local Mean: ",round(netstats$mu_l_formula,round_places),"\nProxy Local Mean: ",round(netstats$mu_l_proxy,round_places),"\nGlobal Mean: ",round(netstats$mu_g,round_places))
  text(x=0.5,y=0.86,atext2,cex=1.25,col="gray40")
  atext3 = bquote(kappa[-1] == .(kappa[1])~" "~kappa[1] == .(kappa[2]))
  atext4 = bquote(kappa[2] == .(kappa[3])~" "~kappa[3] == .(kappa[4]))
  atext5 = bquote(sigma[O] == .(sigmaO)~" "~sigma[ID] == .(sigmaID))
  text(0.5,0.6,atext3,cex=1.25)
  text(0.5,0.5,atext4,cex=1.25)
  text(0.5,0.3,atext5,cex=1.25)
  
  dev.off()
}

#############################################################################################
# Purpose: How did R1 get those Inversity numbers?
# We confirm below that they just used the UNIQUE points on the plot to compute Inversity.
#############################################################################################
x_temp = list()
x_temp[[1]] = cbind(1:14,1:14)
x_temp[[2]] = cbind(1:14,1:14)[-6,]
x_temp[[3]] = cbind(1:14,1:14)[-c(4,5,6,7,8),]
x_temp[[4]] = cbind(1:14,1:14)[-c(2,3,4,5,6,7,8,9,10),]
x_temp[[5]] = cbind(c(1,14),c(1,14))

for (k in 1:5)
{
  x = rbind(x_temp[[k]],x_temp[[k]])
  rho_a = cor(x[,1],x[,2])
  rho_d = cor(x[,1],1/x[,2])
  print(paste("Assortativity: ",round(rho_a,3),"\nInversity:",round(rho_d,3)))
}
#############################################################################################
# Purpose: Replicate reviewer 1's plots -- Part B (Hubs and Pendants)
#############################################################################################
netlist = list()
round_places=3
starsizes = c(11:15)
starsizes = c(2:15)

for (i in starsizes){
  k = i - min(starsizes) +1
  netlist[[k]] = generate_star_network(n=i)
}

gstars_and_pendants = generate_giant_network(netlist)
r = get_degrees_from_edges(gstars_and_pendants)
netstats = get_network_stats(gstars_and_pendants)
kappa = round(netstats$kappa,round_places)
sigmaO = round(netstats$sigmaO,round_places)
sigmaID = round(netstats$sigmaID,round_places)

pname = paste0("Starsizes_",paste0(starsizes,collapse="_"),".pdf")
pdf(pname)
mat_layout <- matrix(c(1,1,2,2,2,3,3,2,2,2),nrow=2,byrow=TRUE)
layout(mat_layout)
par(oma=c(2,2,2,2),mar = c(2,2,1,4))
par(las=1)
gplot(gstars_and_pendants,gmode="graph")
xlim = c(0,max(starsizes)); ylim = c(0,1.2*max(starsizes))
# xlab = expression("d"[i]); ylab1 = expression("d"[j]); ylab2 = expression("1/d"[j])
plot(r$l,r$r,col="blue",type="p",xlim=xlim,ylim=ylim,pch=1,cex.lab=1.5,cex.axis=1.25)
atext1 = paste0("Assortativity: ",round(netstats$rho_a,3),"\nInversity: ",round(netstats$rho_d,3))
# Get proxy local mean.
text(x=07,y=15,atext1,cex=1.25,col="gray40")
par(new=TRUE)
plot(r$l,1/r$r,col="red",type="p",pch=2,xlim=c(0,max(starsizes)),ylim=c(0,1),axes=FALSE)#,xlab=expression("d"[i]),ylab=expression("1/d"[j]),pch=1,cex.lab=1.5)
mtext(expression("d"[i]),side=1,col="black",line=2,cex=1.5) 
mtext(expression("d"[j]),side=2,col="blue",line=2,cex=1.5) 
mtext(expression("1/d"[j]),side=4,col="red",line=1,cex=1.5) 
axis(4, ylim=c(0,1), col="red",col.axis="red",las=1,cex.axis=1.25)
ltext=c(expression('Plot of (d'[i]*', d'[j]*')'),expression('Plot of (d'[i]*', 1/d'[j]*')'))
legend("center",legend=ltext,bty="n",col=c("blue","red"),pch=c(1,2),cex=1.25)
plot.new()
atext2 = paste0("Local Mean: ",round(netstats$mu_l,round_places),"\nFormula Local Mean: ",round(netstats$mu_l_formula,round_places),"\nProxy Local Mean: ",round(netstats$mu_l_proxy,round_places),"\nGlobal Mean: ",round(netstats$mu_g,round_places))
text(x=0.5,y=0.86,atext2,cex=1.25,col="gray40")
atext3 = bquote(kappa[-1] == .(kappa[1])~" "~kappa[1] == .(kappa[2]))
atext4 = bquote(kappa[2] == .(kappa[3])~" "~kappa[3] == .(kappa[4]))
atext5 = bquote(sigma[O] == .(sigmaO)~" "~sigma[ID] == .(sigmaID))
text(0.5,0.6,atext3,cex=1.25)
text(0.5,0.5,atext4,cex=1.25)
text(0.5,0.3,atext5,cex=1.25)
dev.off()


#############################################################################################
# Get India Village Network Statistics
# Both household and individual networks
#############################################################################################
# datadir = "~/Documents/2024_PNAS_Submission_R4/"
# datafile = paste0(datadir,"FriendsBuzz2023.RData")
# load(datafile)

##################################################################
# Figure S13: India Villages Data
##################################################################
datadir = "~/Dropbox/Projects_Active/FoF_David_Vineet/InversityExamples/IndiaVillages/Data/1. Network Data/Adjacency Matrices/"
remove_isolates = function(xnet)
{
  iselect = which(rowSums(xnet)==0)
  rxnet = xnet[-iselect,-iselect]
  return(rxnet)
}
get_village_stats = function(villageset,flag)
{
  flagtext = c("_HH","")
  #fprefix = paste("./IndiaVillages/Data/1. Network Data/Adjacency Matrices/adj_allVillageRelationships",flagtext[flag+1],"_vilno_",sep="")
  fprefix = paste(datadir,"adj_allVillageRelationships",flagtext[flag+1],"_vilno_",sep="")
  fnames = list(); fdata = list()
  for (k in villageset){
    fnames[[k]] = paste(fprefix,k,".csv",sep="")
    tnet1 = read.csv(fnames[[k]],header=FALSE)
    tnet2 = sign(tnet1 + t(tnet1))
    tnet = remove_isolates(tnet2)
    fdata[[k]] = tnet
    print(paste("Village k:",k, "completed with N=",dim(tnet)[1]))
  }
  villagestats = NULL; #matrix(0,max(villages),7)
  for(k in villageset){
    villagestats = rbind(villagestats,unlist(get_network_stats(fdata[[k]])))
  }
  rlist = list(n=fdata,s=villagestats)
  # return(villagestats)
  return(rlist)
}


library(plotrix)
flagtext = c("HH","IND")
villageset = c(1:77)[-c(13,22)]
# villageset = villages_adoption
# v_net_stats_ind = get_village_stats(villageset,flag=1)

# hh_stats = v_net_stats_hh$s #get_village_stats(villageset,flag=0)
# ind_stats = v_net_stats_ind$s #get_village_stats(villageset,flag=1)

flag = 0;
v_net_stats = get_village_stats(villageset,flag=flag)

# vnetstats = hh_stats
# if (flag == 1) vnetstats = ind_stats

vstats = as.data.frame(v_net_stats$s)
# for (k in villageset)
# {
#   print(paste("Starting village: ",k))
#   tvillagestats = get_network_stats(village_HH_networks[[k]])
#   vnetstats = rbind(vnetstats,unlist(tvillagestats))
# }
# vnetstats$mu_l - vnetstats$mu_l_formula

pname=paste0("IndiaVillages_",flagtext[flag+1],"_InversityAssortativityDensity.pdf")
plot(density(vstats$rho_d),xlim=c(-1,1))
lines(density(vstats$rho_a))

pname=paste0("IndiaVillages_",flagtext[flag+1],"_InversityAssortativity.pdf")
pdf(pname)
par(las=1,mgp=c(2,1,0))
xlab = expression(rho[a]);
ylab = expression(rho)
title = "India Villages\n Inversity and Assortativity"
plot(vstats$rho_a,vstats$rho_d,col="blue",xlab=xlab,ylab=ylab,main=title)
select = vstats$rho_a > 0 & vstats$rho_d >0
if (sum(select)>0)
{
  points(vstats[select,]$rho_a,vstats[select,]$rho_d,col="black",pch=19)
  draw.circle(x=0.028,y=0.02,radius=0.025,lty=3)
}  
abline(h=0,lty=3)
abline(v=0,lty=3)
abline(a=0,b=-1,lty=2)
dev.off()

pname=paste0("IndiaVillages_",flagtext[flag+1],"_ProxyLocalMean.pdf")
pdf(pname)
par(las=1,mar=c(2,4,2,2))
xlab = expression(mu[L]);
ylab = expression(mu[L]^proxy)
title = "Proxy Local mean versus True Local mean"
plot(vstats$mu_l,vstats$mu_l_proxy,col="blue",xlab=xlab,ylab=ylab,main=title)
abline(a=0,b=1,lty=3)
dev.off()

villages_interesting = which(select)
if (sum(select)>0)
{
  for (k in villages_interesting)
  {
    points(vstats[select,]$rho_a,vstats[select,]$rho_d,col="black",pch=19)
    pname = paste0("VillageNetwork_",k,"_",flagtext[1+flag],"_Outlier.pdf")
    pdf(pname)
    title = paste0("Village ",k," Network")
    gplot(v_net_stats$n[[k]],gmode="graph",edge.col="gray70",main=title)
    dev.off()
  }
}





#############################################################################################
# Add the networks in Table 1 or 2 where proxy local mean computed using -assortativity 
# as a proxy for inversity
#############################################################################################



