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
  rlist = data.table(nnodes=nnodes,density=x_density,mu_d=mean(ddegrees_original),mu_g=mu_g,mu_l=mu_l,sigma_d=sigma_d,rho_d=rho_d,rho_a=rho_a)
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




set.seed(10) # For replicability.

library(igraph)
library(sna)
library(data.table)

#############################################################################################
# Purpose: Replace Figure S6 from prior revision
#############################################################################################
networkstats = NULL;
cliquesizes = c(2:20)
starsizeset = seq(2,80,1)
gcliquenet = generate_many_cliques(cliquesizes)

mat_layout <- matrix(c(1,2,2,3,2,2),nrow=2,byrow=TRUE)
cliq_min = min(cliquesizes); cliq_max = max(cliquesizes)

for (cstarsize in starsizeset)
{
  tstarnet = generate_star_network(n=cstarsize)
  tnetall = generate_giant_network(list(tstarnet,gcliquenet))
  tnetstats_temp = get_network_stats(tnetall)
  tnetstats = cbind(tnetstats_temp,starsize=cstarsize)
  networkstats = rbind(networkstats,tnetstats)
  
  rho_d = round(tnetstats$rho_d,3); rho_a = round(tnetstats$rho_a,3)
  pname = paste0("Networks_","Starsize_",cstarsize,"_Cliques_",cliq_min,"_to_",cliq_max,".pdf")
  pdf(pname)
  layout(mat_layout)
  par(mar = c(0.1, 0.1,0.1,0.1), oma=c(1,1,1,1))
  gplot(tstarnet,gmode="graph",edge.col="gray90",vertex.col="blue")
  gplot(gcliquenet,gmode="graph",edge.col="gray90")
  plot.new()
  ltext = paste0("Size of star = ",cstarsize,"\n","\n","Cliques with sizes \n ranging from 2 to 20","\n","\n","Inversity = ",rho_d,"\n","\n","Assortativity = ",rho_a)
  text(0.5,0.5,ltext,cex=2)
  dev.off()
  print(paste("Completed: ",pname))
}


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
r = round(get_cor_ia(networkstats,select_subset),3)
print(paste("Correlation between Assortativity and Inversity in Select Subset:",r$cor_ia))
lm_ia = lm(r$i~r$a)
r_sq_ia = round(summary(lm_ia)$r.squared,3)
print(paste("R^2 for linear model between inversity and assortativity in Select Subset: ",r_sq_ia))


#############################################################################################
# Purpose: Replace Figure S6. Show how Inversity and Assortativity vary with star size
#############################################################################################

xlab="Size of Star network"
ylab = "Inversity and Assortativity"
ltext=c("Inversity","Assortativity")
atext = "Region of same sign \n for assortativity \n and inversity"
pdf("InversityAssortativityExampleB10.pdf")
plot(networkstats$starsize,networkstats$rho_d,type="b",col="blue",cex=0.5,cex.axis=1.1,cex.lab=1.5,xlim=c(0,60),ylim=c(-1,1.5),xlab=xlab,ylab=ylab)
lines(networkstats$starsize,networkstats$rho_a,type="b",col="darkred",cex=0.5,pch=3)
abline(h=0,lty=2)
abline(h=1,lty=2)
abline(h=-1,lty=2)
text(x=37,y=-0.5,atext,cex=1)
legend("topright",legend=ltext,bty="n",pch=c(1,3),cex=1.25,col=c("blue","darkred"))
lgray1 = rgb(0.05,0.05,0.05,0.02)
lgray2 = rgb(0.15,0.15,0.15,0.08)
rect(28,-1,46,1,col=lgray2,lty=2)
# rect(0,-1,260,1,col=lgray1,lty=2,density=100)
# rect(260,-1,960,1,col=lgray2,lty=2)
# rect(960,-1,1000,1,col=lgray1,lty=2,density=100)
# rect(57,-1,100,1,col=lgray,lty=2)
dev.off()

#############################################################################################
s
