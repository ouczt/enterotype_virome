#load packages
library(cluster)
library(clusterSim)
library(ade4)
library("MASS")
library("vegan")
library("permute")
library("lattice")
library("ggplot2")
library("ggrepel")
library("compositions")
library("plyr")
library("phyloseq")


******************************************************
#begin enterotying
data=otu_GZ_RA
colnames(otu_GZ_RA)

data.dist2=vegdist(otu_GZ_RA_t,method="bray")

##基于bray的dis
pam.clustering=function(x,k) { # x is a distance matrix and k the number of clusters
  require(cluster)
  cluster = as.vector(pam(as.dist(x), k, diss=TRUE)$clustering)
  return(cluster)
}

require(clusterSim)

nclusters=NULL
for (k in 1:10) { 
  if (k==1) {
    nclusters[k]=NA 
  } else {
    data.cluster_temp1=pam.clustering(data.dist2, k)
    nclusters[k]=index.G1(t(data),data.cluster_temp1,  d = data.dist2,
                          centrotypes = "medoids")
  }
}

plot(nclusters, type="h", xlab="k clusters", ylab="CH index",main="Optimal number of clusters")
##t informs me that k=2 is optimal,however when plotting

data.cluster=pam.clustering(data.dist2, k=2)
nclusters = index.G1(t(data), data.cluster, d = data.dist2, centrotypes = "medoids")

colnames(data)
str(data.cluster)

subject_enterotype_GZ<-data.frame(sites=colnames(data),enterotype=data.cluster)
row.names(subject_enterotype_GZ)<-colnames(data)
##data即species_RPKM_RA
write.table(subject_enterotype_GZ,"subject_enterotype_GZ.txt",sep="\t")
head(subject_enterotype)

