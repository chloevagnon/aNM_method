# The following R codes are intended to provide step-by-step procedure to apply the allometric Niche model 
# to data from a simplified realistic lake foodweb.

library(ggplot2)
library(cheddar)
library(igraph)
library(gridExtra)

source("get_niche_attributes.R")
source("weights.R")
source("export_network3D2.R")

load("Param_regvert.Rdata")
load("Param_reginvert.Rdata")

  
fw <-read.delim("simplified_foodweb.txt", header=T, dec=",")
head(fw)

fw$log_mean_size<-log10(fw$mean_size)

p1<-ggplot(fw)+
geom_bar(aes( x=log_mean_size,fill = phylogenetic_category), position = "dodge")
p2<-ggplot(fw)+
  geom_bar(aes( x=lake_compartment,fill =lake_compartment))
grid.arrange(p1, p2, nrow = 1)


Param_regvert

Param_reginvert

consumers<-fw[fw$phylogenetic_category %in% c("vertebrate","invertebrate"),]

Niche_attributes<-data.frame(matrix(NA,nrow=nrow(consumers),ncol=4))
rownames(Niche_attributes)<-consumers$genus
colnames(Niche_attributes)<-c("n","c","low","upp")

for (i in 1:nrow(consumers)){
  Niche_attributes[i,]<-get_niche_attributes(consumer_size=consumers$log_mean_size[i], 
                                             consumer_category=consumers$phylogenetic_category[i])
  
}
Niche_attributes


binary_matrix<-matrix(0,ncol=nrow(fw),nrow=nrow(fw))
colnames(binary_matrix)<-fw$genus
rownames(binary_matrix)<-fw$genus

for (i in 1:nrow(consumers)){
 preys<- which(fw$log_mean_size>Niche_attributes[i,3] & fw$log_mean_size<Niche_attributes[i,4])
 binary_matrix[preys,i]<-1
}

  
for (i in 1:nrow(binary_matrix)){
  for (j in 1:nrow(consumers)){
  
   if(binary_matrix[i,j]==1 & 
       fw[fw$genus==rownames(binary_matrix)[i],]$lake_compartment == "pel" & 
     fw[fw$genus==colnames(binary_matrix)[j],]$lake_compartment == "ben" |
       binary_matrix[i,j]==1 & 
       fw[fw$genus==rownames(binary_matrix)[i],]$lake_compartment == "ben" &
      fw[fw$genus==colnames(binary_matrix)[j],]$lake_compartment == "pel")
    {binary_matrix[i,j]<-0}
  }
}
binary_matrix[1:10,1:10]


N_links<-nrow(binary_matrix)

inferred_predator_prey_links<-data.frame(matrix(NA,ncol=4,nrow=N_links*N_links))
colnames(inferred_predator_prey_links)<-c("predator","prey","predator_size","prey_size")

for (i in 1: nrow(consumers)){
  
  dietsp<-binary_matrix[which(binary_matrix[,i]>0),i]
  
  for (j in 1: length(dietsp)){
    inferred_predator_prey_links[j+(i-1)*N_links,1] <- colnames(binary_matrix)[i]
    inferred_predator_prey_links[j+(i-1)*N_links,2] <- names(dietsp)[j]
    inferred_predator_prey_links[j+(i-1)*N_links,3] <- fw[fw$genus==colnames(binary_matrix)[i],]$log_mean_size
    inferred_predator_prey_links[j+(i-1)*N_links,4] <- fw[fw$genus==names(dietsp)[j],]$log_mean_size
  }
}
inferred_predator_prey_links<-na.omit(inferred_predator_prey_links)
rownames(inferred_predator_prey_links)<-seq(from=1,to=nrow(inferred_predator_prey_links),by=1)
inferred_predator_prey_links$prednum<-as.numeric(as.factor(inferred_predator_prey_links$predator))
head(inferred_predator_prey_links)

weighted_links<-matrix(NA,ncol=1,nrow=nrow(inferred_predator_prey_links))
for(i in 1:max(inferred_predator_prey_links$prednum)){
  
  sp<-inferred_predator_prey_links[inferred_predator_prey_links$prednum==i,]
  
  niche_sp<-as.matrix(Niche_attributes[which(rownames(Niche_attributes)==sp[1,1]),c(2:4)])
 
  weighted_links[as.numeric(rownames(sp)),1]<-weights(prey_sizes = sp[,4],niche_center=niche_sp[1],niche_min=niche_sp[2],niche_max=niche_sp[3])
}
inferred_predator_prey_links$weighted_links<-weighted_links
head(inferred_predator_prey_links)



the pike niche attributes are (units = log10(µm))
pike_niche<-as.matrix(Niche_attributes[1,])
pike_niche


10^-4*10^pike_niche

10^-4*10^inferred_predator_prey_links[1:3,4]


prey_range_vector<-seq(from=pike_niche[3],to=pike_niche[4],length.out=100)

link_weight_density<-weights(prey_range_vector,niche_center=pike_niche[2],niche_min=pike_niche[3],niche_max=pike_niche[4])
                             

op <- par(mfrow=c(1,1),mar = c(2,2,1,1) + 0.5,mgp = c(1, 0.3, 0))
plot(link_weight_density~prey_range_vector,type="l",xlab="prey_size_range (r, log10(µm))")
abline(v=Niche_attributes[1,c(2:4)],col=c(1,2,2),lty=2)
points(weighted_links~prey_size,data=inferred_predator_prey_links[1:3,],col="green",pch=19,cex=2)


weighted_matrix<-matrix(0,ncol=nrow(fw),nrow=nrow(fw))
colnames(weighted_matrix)<-fw$genus
rownames(weighted_matrix)<-fw$genus

for(i in 1:nrow(inferred_predator_prey_links)){
  
  pred<-which(rownames(weighted_matrix)==inferred_predator_prey_links$predator[i])
  prey<-which(rownames(weighted_matrix)==inferred_predator_prey_links$prey[i])
  
  weighted_matrix[prey,pred]<-inferred_predator_prey_links$weighted_links[i]
}
weighted_matrix[1:10,1:10]


nodes<-data.frame(fw$genus,fw$log_mean_size);colnames(nodes)<-c("node","M")
trophic.links<-inferred_predator_prey_links[,c(2,1,6)];colnames(trophic.links)<-c("resource","consumer","weight")
property<-list(title=c("Simplified_foodweb"),M.units="Log10(µm)")

simplified_foodweb<-Community(nodes,properties=property,trophic.links=trophic.links)

topology<-as.matrix(c(DirectedConnectance(simplified_foodweb),
                              LinkageDensity(simplified_foodweb),
                              length(TopLevelNodes(simplified_foodweb))/nrow(simplified_foodweb$nodes),
                              length(IntermediateNodes(simplified_foodweb))/nrow(simplified_foodweb$nodes),
                              length(BasalNodes(simplified_foodweb))/nrow(simplified_foodweb$nodes)))

topological_metrics<-rbind(nrow(nodes), nrow(inferred_predator_prey_links), topology)
rownames(topological_metrics) <-c("Species richness","Links","DirectedConnectance", "Likage density", "Top %", "Inter %", "Basal %")
topological_metrics


graph <- graph.adjacency(weighted_matrix, mode="directed",diag=F,weighted=T)

metrics<-rbind(degree(graph),betweenness(graph,directed = TRUE, weights =E(graph)$weight),
                      transitivity(graph,type=c("weighted"), weights =E(graph)$weight),
                      TrophicVulnerability(simplified_foodweb),
                      TrophicGenerality(simplified_foodweb),
                      PreyAveragedTrophicLevel(simplified_foodweb))

rownames(metrics)<-c("Degree","Betweenness","Transitivity","Trophic_Vul","Trophic_Gen","Trophic_position")
species_metrics<-data.frame(fw[,c(2,4,6,7)],t(metrics))
species_metrics

ExportToNetwork3D2(simplified_foodweb,species_metrics,'D:/Desktop')
