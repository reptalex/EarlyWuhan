library(ape)
library(phangorn)
library(data.table)
library(magrittr)
library(ggplot2)
library(ggpubr)
library(parallel)
library(ggtree)
library(RRphylo)
library(igraph)
source('scripts/nb_outbreak.R')
theme_set(theme_bw(base_size=15))

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
# Universal parameters ----------------------------------------------------
### Serial interval edge lengths
## https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7128842/
## median 4d, mean 4.7d, sd 2.9d

## median = exp(u)
## mean = exp(u+sig^2/2)
## 
u=log(4)
sigma=sqrt((log(4.7)-u)*2)

q=0.02 ## probability of hospitalization
A=10    ## Alert threshold


# Early Outbreak Simulation -----------------------------------------------

### Outbreak: Galton-Watson process with NB offspring distribution
set.seed(2)
X <- GW_no_extinction()
tree=gw_tree(X,u,sigma) ##u, sigma define log-normal edge length distribution
# plot(tree)
# nodelabels(tree$node.label)


# One important comment: the resulting tree from this GW process is slighty different
# from typical sequence-based phylogenies in that nodes of our tree our patients
# that transmitted the virus a cluster of descendants. This corresponds to a mutation
# model in which every transmission event generates a mutation in a new pt, and the genome
# of the ancestral patient would align perfectly with the MRCA of its descendants.

## Consequently, our contact tracing will also draw the nodes of the tree.
## A similar analysis with a more complex mutation model is left for future research,
## and here we show the core, mathematical insight that contact tracing
## leads to biases in the sample phylogeny relative to the true phylogeny.



# Early cases + phylogeny ------------------------------------------------

### Hospitalizations + Medical Alert
pt_database <- hospitalize(X,q,tree)
id <- pt_database[cum_hosp>=A,id[1]]
index_pts <- pt_database[hospitalized==1 & cum_hosp<=A,id]

node_depths=node.depth.edgelength(tree)
names(node_depths) <- c(tree$tip.label,tree$node.label)
t0=node_depths[id]

tr <- trim_tree(tree,id)



# Contact Tracing ---------------------------------------------------------
### Sampling propensities under contact-tracing
half_life=0.1 # What is half-life of contact-tracing decay along xmsn/phylo distance?
lambda=log(2)/half_life

propensities <- ascertainment_propensities(tr,index_pts,t0,
                                              u,sigma,lambda)


#### Contact tracing vs. random sampling graph properties
DAT <- dists_contact_tracing(propensities,N=40,reps=100)

g_ix=DAT[,list(nearest_index_pt=nearest_index_pt[1]),by=c('rep','sampling_process')] %>%
  ggplot(aes(sampling_process,nearest_index_pt,color=sampling_process))+
  geom_boxplot(lwd=1.5)+
  geom_jitter(cex=4,alpha=0.5)+
  scale_y_continuous('Degrees of Separation')+
  ggtitle('Degrees of Separation to an Index Patient')+
  theme(legend.position='none')+
  geom_hline(yintercept = 2,lty=2)+
  annotate(geom='text',x=1.8,y=2.1,label='Patients from Same Cluster')

g_no=DAT[degrees_of_separation<=7] %>%
ggplot(aes(factor(degrees_of_separation),
       number_of_pts,color=sampling_process))+
  geom_boxplot(lwd=1.5)+
  geom_point(position=position_jitterdodge(),alpha=0.3)

g_perc=DAT[degrees_of_separation<=7] %>%
  ggplot(aes(factor(degrees_of_separation),
             number_of_pts*100,color=sampling_process))+
  geom_boxplot(lwd=1.5)+
  geom_point(position=position_jitterdodge(),alpha=0.5)+
  geom_smooth(aes(x=degrees_of_separation))+
  scale_x_discrete('Degrees of Separation')+
  scale_y_continuous('Percent of patient-pairs')+
  theme(legend.position=c(0.4,0.9))+
  ggtitle('Patient Degrees of Separation')+
  guides(color=guide_legend(ncol=2))
g_perc


# Phylogeny visualization -------------------------------------------------
set.seed(2)
focal_pt <- sample(index_pts,1)
remaining_pts <- setdiff(propensities$patient,focal_pt)
ct_5 <- ct(focal_pt,k=5)
ct_20 <- ct(ct_5,15)
ct_40 <- ct(ct_20,20)

r_5 <- c(focal_pt,sample(remaining_pts,5))
r_20 <- c(r_5,sample(setdiff(remaining_pts,r_5),15))
r_40 <- c(r_20,sample(setdiff(remaining_pts,r_20),20))




############## Plotting ######################
png('figures/contact_tracing_vs_random_sampling_trees.png',height=10,width=16,res = 400,units = 'in')
  par(mfrow=c(2,4))
  cls=gg_color_hue(2)
  
  ###### Random Patients
  plot(tr,show.tip.label = F,main='Index Patient')
  tiplabels(text=NA,tip = which(tr$tip.label==focal_pt),pch = 21,cex=4,bg=cls[2])
  
  plot(tr,show.tip.label = F,main='+5 Random Patients')
  tiplabels(text=NA,tip = which(tr$tip.label %in% r_5),pch = 21,cex=4,bg=cls[2])
  nodelabels(text=NA,node=Ntip(tr)+which(tr$node.label%in% r_5),pch = 21,cex=4,bg=cls[2])
  
  plot(tr,show.tip.label = F,main='+20 Random Patients')
  tiplabels(text=NA,tip = which(tr$tip.label %in% r_20),pch = 21,cex=4,bg=cls[2])
  nodelabels(text=NA,node=Ntip(tr)+which(tr$node.label%in% r_20),pch = 21,cex=4,bg=cls[2])
  
  plot(tr,show.tip.label = F,main='+40 Random Patients')
  tiplabels(text=NA,tip = which(tr$tip.label %in% r_40),pch = 21,cex=4,bg=cls[2])
  nodelabels(text=NA,node=Ntip(tr)+which(tr$node.label%in% r_40),pch = 21,cex=4,bg=cls[2])
  
  ########### Contact Tracing
  plot(tr,show.tip.label = F,main='Index Patient')
  tiplabels(text=NA,tip = which(tr$tip.label==focal_pt),pch = 21,cex=4,bg=cls[1])
  
  plot(tr,show.tip.label = F,main='+5 Contact-Traced Patients')
  tiplabels(text=NA,tip = which(tr$tip.label %in% ct_5),pch = 21,cex=4,bg=cls[1])
  nodelabels(text=NA,node=Ntip(tr)+which(tr$node.label%in% ct_5),pch = 21,cex=4,bg=cls[1])
  
  plot(tr,show.tip.label = F,main='+20 Contact-Traced Patients')
  tiplabels(text=NA,tip = which(tr$tip.label %in% ct_20),pch = 21,cex=4,bg=cls[1])
  nodelabels(text=NA,node=Ntip(tr)+which(tr$node.label%in% ct_20),pch = 21,cex=4,bg=cls[1])
  
  plot(tr,show.tip.label = F,main='+40 Contact-Traced Patients')
  tiplabels(text=NA,tip = which(tr$tip.label %in% ct_40),pch = 21,cex=4,bg=cls[1])
  nodelabels(text=NA,node=Ntip(tr)+which(tr$node.label%in% ct_40),pch = 21,cex=4,bg=cls[1])

dev.off()



DAT[,list(nearest_index_pt=nearest_index_pt[1]),by=c('rep','sampling_process')] %>%
  ggplot(aes(sampling_process,nearest_index_pt,color=sampling_process))+
  geom_boxplot(lwd=1.5)+
  geom_jitter(cex=4,alpha=0.5)+
  scale_y_continuous('Degrees of Separation')+
  ggtitle('Degrees of Separation to an Index Patient')+
  theme(legend.position='none')+
  geom_hline(yintercept = 2,lty=2)+
  annotate(geom='text',x=1.8,y=2.1,label='Patients from Same Cluster')
ggsave('figures/Degrees_of_separation_to_pt_0_by_sampling_process.png',height=10,width=6)


# etc ---------------------------------------------------------------------



#### Contact tracing vs. random sampling, tree visualization
##### This code uses the more complex contact tracing algorithm with
##### phylogenetic distance-based exponential decay in ascertainment propensity
#### It's fun & neat code, but too complex for the purposes of paper.
N=40
random_pts=c(pts,sample(propensities[ascertained==FALSE,patient],size=N,replace=F))
ct_pts=contact_trace(propensities,N=N,simultaneous=F,u=u,sigma=sigma,lambda=lambda)

cls=gg_color_hue(2)

par(mfrow=c(1,2))
  plot(tr,type='fan',rotate.tree = 90,show.tip.label = F,main='Random Sampling')
  tiplabels(text=NA,tip = which(tr$tip.label %in% index_pts),pch = 16,cex=2)
  nodelabels(text=NA,node=Ntip(tr)+which(tr$node.label%in% index_pts),pch=16,cex=2)
  tiplabels(text=NA,tip = which(tr$tip.label %in% setdiff(random_pts,index_pts)),
            pch = 16,cex=1.5,col=cls[2])
  nodelabels(text=NA,node=Ntip(tr)+which(tr$node.label%in% setdiff(random_pts,index_pts)),
             pch=16,cex=1.5,col=cls[2])
  
  plot(tr,type='fan',rotate.tree = 90,show.tip.label = F,main='Contact Tracing')
  tiplabels(text=NA,tip = which(tr$tip.label %in% index_pts),pch = 16,cex=2)
  nodelabels(text=NA,node=Ntip(tr)+which(tr$node.label%in% index_pts),pch=16,cex=2)
  tiplabels(text=NA,tip = which(tr$tip.label %in% setdiff(ct_pts,index_pts)),
            pch = 16,cex=1.5,col=cls[1])
  nodelabels(text=NA,node=Ntip(tr)+which(tr$node.label%in% setdiff(ct_pts,index_pts)),
             pch=16,cex=1.5,col=cls[1])
