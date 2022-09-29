library(ape)
library(data.table)
nb_gw <- function(R0=2.4,
                  Cases=1e3,
                  k=0.1){
  
    ### simulate first 1,000 cases in branching process
    R0           #=pr/(1-p)
    v=R0+k*R0^2  #=pr/(1-p)^2
    r=(R0^2)/(v^2-R0)
    
    N=0
    N=rnbinom(1,r,mu=R0)
    
    if (N>0){
    # while(N==0){
    #   N=rnbinom(1,r,mu=R0)
    # }
    X <- data.frame('patient'=1,
                    'descendants'=1+1:N,
                    'generation'=0)
    complete=0
    desc <- X$descendants
    pt=N+1
    tm=0
    while(max(desc)<Cases){
      
      dd=NULL  ## same as X, but for all descendants whose 
               ## transmission we haven't yet tracked
      
      tm=tm+1
      for (d in desc){  ## track the transmission of descendants
        N=rnbinom(1,size=r,mu=R0)
        if (N>0){
          dd <- rbind(dd,data.frame('patient'=d,
                                    'descendants'=pt+1:N,
                                    'generation'=tm))
          pt=max(dd$descendants)
        }
      }
      
      if (is.null(dd)){ ## all descendant lineages extinct
        break
      } else {
        X <- rbind(X,dd)
        desc <- dd$descendants
      }
      
    }
    } else {
      X <- data.frame('patient'=1,
                      'descendants'=1,
                      'generation'=0)
    }
    X <- as.data.table(X)
    return(X)
}


GW_no_extinction <- function(R0=2.4,
                             Cases=1e3,
                             k=0.1){
  complete=FALSE
  while(!complete){
    X=nb_gw(R0,Cases,k)
    if (max(X$descendants)>Cases){
      complete=TRUE
    }
  }
  return(X)
}

gw_tree <- function(X,u=1.39,sigma=0.568){
  tips <- setdiff(unique(X$descendants),unique(X$patient)) ## patients w/o descendants = tips
  parents <- setdiff(unique(X$patient),tips)
  
  
  node_map <- data.table('pt'=c(tips,parents),
                         'node'=1:(length(tips)+length(parents)))
  
  mat <- X[,c('patient','descendants')]
  mat[,patient:=node_map[match(patient,pt),node]]
  mat[,descendants:=node_map[match(descendants,pt),node]]
  
  tree <- NULL
  tree$tip.label <- as.character(tips)
  tree$Ntip=length(tips)
  tree$Nnode=length(parents)
  tree$node.label=paste0('n',parents)
  tree$edge <- as.matrix(mat)
  tree$edge.length <- rlnorm(nrow(tree$edge),meanlog=u,sdlog=sigma)
  class(tree) <- 'phylo'
  return(tree)
}

hospitalize <- function(X,q,tree){
  pt_database <- data.table('id'=c(tree$tip.label,tree$node.label))
  pt_database[,time:=node.depth.edgelength(tree)]
  pt_database <- pt_database[order(time)]
  pt_database[,hospitalized:=rbinom(.N,1,q)]
  pt_database[,cum_hosp:=cumsum(hospitalized)]
  return(pt_database)
}

trim_tree <- function(tree,id){
  #### Some manual preprocessing is needed for an edge case
  #### If we have a-->b-->c and we want to cut the branch b-->c
  #### the ape tools seem to 
  
  
  
  ### Subset tree to all patients at time of Alert
  node_depths=node.depth.edgelength(tree)
  names(node_depths) <- c(tree$tip.label,tree$node.label)
  t0=node_depths[id]
  tr=tree
  while(any(node_depths>t0)){
    
    excluded_tips <- tr$tip.label[intersect(which(node_depths>t0),1:Ntip(tr))]
    for (tp in excluded_tips){
    tr=drop.tip(tr,tp,
                trim.internal = F,
                collapse.singles = F)
    }
    node_depths=node.depth.edgelength(tr)
    names(node_depths) <- c(tr$tip.label,tr$node.label)
    
  }
  return(tr)
}

ascertainment_propensities <- function(tr,pts,t0,xmsn_distance=TRUE,
                                       u=1.30,sigma=0.568,lambda=1){
  ### Need pairwise distances frome very node & tip to pts
  ### Some pts are nodes
  tip_pts <- pts[pts %in% tr$tip.label]
  node_pts <- setdiff(pts,tip_pts)
  
  tip_ix=match(tip_pts,tr$tip.label)
  node_ix=match(node_pts,tr$node.label)+Ntip(tr)
  
  if (xmsn_distance){
    els=tr$edge.length
    tr$edge.length=rep(1,length(tr$edge.length))
    
    D <- dist.nodes(tr)
    phylo_distances=apply(D[,c(tip_ix,node_ix)],1,min)
    
    tr$edge.length <- els
  } else {
    D <- dist.nodes(tr)
    phylo_distances=rowSums(D[,c(tip_ix,node_ix)])
  }
  
  
  times_since_case=t0-node.depth.edgelength(tr)
  propensities=exp(-lambda*phylo_distances)*dlnorm(times_since_case,meanlog = u,sdlog = sigma)
  
  props=data.table('patient'=c(tr$tip.label,tr$node.label),
             'phylo_distances'=phylo_distances,
             'propensities'=propensities)
  props[,ascertained:=patient %in% pts]
  return(props)
}


contact_trace <- function(propensities,N=10,
                          simultaneous=TRUE,t0.=t0,tr.=tr,
                          u=1.30,sigma=0.568,lambda=1){
  pts=propensities[ascertained==TRUE,patient]
  if (simultaneous){
    uas=which(propensities$ascertained==FALSE)
    pts=c(pts,sample(propensities[uas,patient],size=N,
           replace=F,prob=propensities[uas,propensities]))
  } else {
    props=propensities
    for (i in 1:N){
      uas=which(props$ascertained==FALSE)
      pts=c(pts,sample(props[uas,patient],1,prob=props[uas,propensities]))
      props=ascertainment_propensities(tr,pts,t0,u,sigma,lambda)
    }
  }
  return(pts)
}

ct <- function(pts=NULL,k=5,tr.=tr,detection_prob=0.8,
               u=1.30,sigma=0.568,max_pcr_test_prob=0.9){
  
  ## All cases 1-degree away from focal case, or 2-degrees through shared ancestral node,
  ## are tested with prob detection_prob. Tests turn up positive with max_pcr_test_prob.
  
  
  md=exp(u-sigma^2)
  pk=dlnorm(md,u,sigma)
  scaling_factor=peak_prob/pk  # formula for pcr_test_prob : dlnorm(t0-t,u,sigma)*scaling_factor
  
  ### current patients can be tips or nodes of transmission tree
  tip_pts <- pts[pts %in% tr$tip.label]
  node_pts <- setdiff(pts,tip_pts)
  
  tip_ix=match(tip_pts,tr$tip.label)
  node_ix=match(node_pts,tr$node.label)+Ntip(tr)
  
  
  ### convert transmission tree to transmission distances
  els=tr$edge.length
  tr$edge.length=rep(1,length(tr$edge.length))
  D <- dist.nodes(tr)
  tr$edge.length <- els
  
  phylo_distances=apply(D[,c(tip_ix,node_ix),drop=F],1,min)
  
  times_since_case=t0-node.depth.edgelength(tr)
  
  ### pts must be <=2 dgrs away * not already a pt * may not test positive from PCR test.
  probs = (phylo_distances<=2)*(phylo_distances!=0)*(dlnorm(times_since_case,u,sigma)*scaling_factor)
  
  n_pts <- min(c(sum(probs!=0),k)) ## We'll draw n_pts total
  
  all_patients=c(tr$tip.label,tr$node.label)
  
  pts <- c(pts,sample(all_patients,size=n_pts,replace=F,prob=probs))
  return(pts)
}

#' ... optional input args for contact_trace
dists_contact_tracing <- function(propensities,N=10,reps=100,tr.=tr,...){
  els <- tr$edge.length
  
  ##phylogenetic distances
  Dphy <- dist.nodes(tr)
  rownames(Dphy) <- c(tr$tip.label,tr$node.label)
  colnames(Dphy) <- c(tr$tip.label,tr$node.label)
  
  ##Degrees of separation
  tr$edge.length <- rep(1,Nedge(tr))
  D <- dist.nodes(tr)
  rownames(D) <- c(tr$tip.label,tr$node.label)
  colnames(D) <- c(tr$tip.label,tr$node.label)
  
  DAT=NULL
  pts <- propensities[ascertained==TRUE,patient]
  for (i in 1:reps){
    random_pts=c(pts,sample(propensities[ascertained==FALSE,patient],size=N,replace=F))
    dr=D[random_pts,random_pts]
    dr_phy=Dphy[random_pts,random_pts]
    
    
    ct_pts=contact_trace(propensities,N=N,...)
    dc=D[ct_pts,ct_pts]
    dc_phy=Dphy[ct_pts,ct_pts]

    
    for (j in 1:max_dst){
      r_count=mean(rowSums(dr==j))
      ct_count=mean(rowSums(dc==j))
      
      r_cum=sum(dr[upper.tri(dr)]<=j)/sum(upper.tri(dr))
      ct_cum=sum(dc[upper.tri(dc)]<=j)/sum(upper.tri(dr))
      
      
      DAT=rbind(DAT,data.table('N'=N,'rep'=i,'degrees_of_separation'=j,
                               'sampling_process'=c('random','contact_tracing'),
                               'cum_fraction_of_pts'=c(r_cum,ct_cum),
                               'number_of_pts'=c(r_count,ct_count),
                               'nearest_index_pt'=-1))
    }
    r_nearest_index=mean(apply(dr[setdiff(random_pts,pts),pts],1,min))
    ct_nearest_index=mean(apply(dc[setdiff(ct_pts,pts),pts],1,min))
    
    DAT[rep==i & sampling_process=='random']$nearest_index_pt <- r_nearest_index
    DAT[rep==i & sampling_process=='contact_tracing']$nearest_index_pt <- ct_nearest_index
    
  }
  tr$edge.length <- els
  return(DAT)
}
