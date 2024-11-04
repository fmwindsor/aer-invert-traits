### APPENDIX F. R functions for estimating functional measures using fuzzy coded trait data.
## Written by Cayetano Gutiérrez-Cánovas. Aquatic Ecology. Universidad de Murcia (Spain). Email: cayeguti@um.es

### Input variables:
# ATTENTION: check that traits, taxa and stressor matrices have row and column names ordered in the same way
# traits: matrix of traits. Must be fuzzy coding or percentages of affinity for each trait category (taxa x trait categories)
# taxa: taxonomic matrix of abundances or presences (sites x taxa)
# col.blocks: vector with number of categories of each trait. Names are the trait names. 

### Input variables to run funcNiche() and funcNiche.pca()
# rdt: number of randomizations to generate the taxon niche space
# m: number of axes

# Output variables:
# fpc: functional space (functional principal components)
# tfric.data: area and centroid of each taxon in the functional space
# overlap.mat: triangular matrix of niche overlap between taxon pairs
# tFRic: mean taxon functional richness for each community
# FSim: mean functional similarity between taxon pairs
# FDis: functional dispersion of each community
# FRic: functional richness of each community

### Functions to prepare matrices, construct the functional space and to estimate functional diversity components
# prep.fuzzy.df (traits, col.blocks))  # Convert fuzzy code into percentages of affinity for each taxon and trait
# num.pca.axes (traits) ## Calculates the number of significant axes for the functional space (m)
# funcSpace (taxa,traits,m=2,rdt=9) # Calculates the functional space (fpc)

### Functions to estimate functional diversity components (not to run, included in funcNiche)
# pol.area (P) ## To calculate the area of a polygon
# int.area (A,B) ## To calculate the area of the overlapping area between two polygons
# taxon.fric (fpc,m,rdt) ## To estimate the area and centroid of each taxon in the functional space 
# mean.taxon.fric (taxa,tfric.data) ## To estimate mean area and centroid for each community
# taxon.overlap (taxa,fpc,m,rdt,tfric.data) ## To calculate the triangular matrix of niche overlap between taxon pairs
# fr.sim<-function(taxa,nic.overlap) ## To estimate Functional similarity and Functional Redundancy for each community
# fdis (taxa,fpc,m,rdt) ## To estimate Functional dispersion for each community
# fric<-function(taxa,fpc,m,rdt) ## To estimate Functional richness for each community

### Function to determine the optimal number of randomizations (rdt):
# det.rdt (taxa,traits,calc.pca=T, m=2,test.rdt,plots=T)

### Functions to estimate Functional Diversity components:
# funcNiche.pca (taxa,traits,m,rdt) # Also calculates functional space
# funcNiche (taxa,fpc,m,rdt) # It requires a previously calculated funcional space

### Functions to estimate null models and test for significant differences between observed and simulated model paramters:
# simul.func (taxa,fpc,m=2,rdt,runs=99) # Produces a number of simulated communities and the associated funcional components
# test.simul (taxa,real.res,sim.data,stressor) # Calculates model parameters for observed and simulated functional diversity components



#### Required libraries
library(ade4)
library(vegan)

### Functions to prepare matrices, construct the functional space and to estimate functional diversity components
prep.fuzzy.df<-function (traits, col.blocks) 
{
  if (!is.data.frame(traits)) 
    stop("Data.frame expected")
  if (sum(col.blocks) != ncol(traits)) {
    stop("Non convenient data in col.blocks")
  }
  if (is.null(names(col.blocks))) {
    names(col.blocks) <- paste("FV", as.character(1:length(col.blocks)), sep = "")
  }
  f1 <- function(x) {
    a <- sum(x)
    if (is.na(a)) 
      return(rep(0, length(x)))
    if (a == 0) 
      return(rep(0, length(x)))
    return(x/a)
  }
  k2 <- 0
  col.w <- rep(1, ncol(traits))
  
  for (k in 1:(length(col.blocks))) {
    k1 <- k2 + 1
    if (col.blocks[k]==1) k2<-k1 else k2 <- k2 + col.blocks[k]
    X <- as.matrix(traits[, k1:k2])
    if (col.blocks[k]==1) X[which(X[,1]>0),]<-1 else X <- t(apply(X, 1, f1))
    X.marge <- apply(X, 1, sum)
    X.marge <- X.marge
    X.marge <- X.marge/sum(X.marge)
    X.mean <- apply(X * X.marge, 2, sum)
    nr <- sum(X.marge == 0)
    cat(nr, "missing data found in block", k, "\n")
    traits[, k1:k2] <- X
    col.w[k1:k2] <- X.mean
  }
  attr(traits, "col.blocks") <- col.blocks
  attr(traits, "col.freq") <- col.w
  col.num <- factor(rep((1:length(col.blocks)), col.blocks))
  attr(traits, "col.num") <- col.num
  return(traits)
}


num.pca.axes<-function(traits){
  res<-list()
  par(mfrow=c(1,1))
  screeplot(princomp(traits)->pca.tr,bstick=T,type="lines")
  res$comp<-(pca.tr$sdev^2)
  res$m<-length(pca.tr$sdev[which(res$comp>mean(res$comp))])
  res$var<-sum(res$comp[1:res$m])/sum(res$comp)
  return(res)
}

funcSpace<-function(traits,col.blocks,m=2,rdt=9) {
  
  func.pca<-function(traits,m){
    trait.pca<-princomp(traits) # new PCA space is created based on the whole variability of genera. Covariance matrix (delete scale=F)
    res<-list()
    res$pca<-trait.pca
    res$pca.scores<-trait.pca$scores
    res$pca.var<-trait.pca$sd[1:m]^2/sum(trait.pca$sd^2)
    res$pca.cumvar<-sum(trait.pca$sd[1:m]^2/sum(trait.pca$sd^2))
    return(res)
  }
  
  func.pca(traits,m)->trait.pca.object
  cat("Explained variance by the first", m, "components:",trait.pca.object$pca.cumvar,"\n")
  
  prep.pca<-function (traits, col.blocks) {
    ## This function was adapted from prep.fuzzy.var() from ade4 package
    
    if (!is.data.frame(traits)) 
      stop("data.frame expected")
    
    if (sum(col.blocks) != ncol(traits)) {
      stop("non convenient data in col.blocks")
    }
    
    if (is.null(names(col.blocks))) {
      names(col.blocks) <- paste("FV", as.character(1:length(col.blocks)), sep = "")
    }
    
    f1 <- function(x) {
      v<-(rep(0,length(x)))
      prob=c(as.numeric(x[x>0]/sum(x)))
      if(sum(x)>0) (sample(which(x==x),1,prob=c(as.numeric(x/sum(x)))))->value
      if(sum(x)>0) v[value]<-1
      return(v)
    }
    
    k2 <- 0
    col.w <- rep(1, ncol(traits))
    for (k in 1:(length(col.blocks))) {
      k1 <- k2 + 1
      if (col.blocks[k]==1) k2<-k1 else k2 <- k2 + col.blocks[k]
      X <- as.matrix(traits[, k1:k2])
      if (col.blocks[k]==1) X[which(X[,1]>0),]<-1 else X <- t(apply(X, 1, f1))
      X.marge <- apply(X, 1, sum)
      X.marge <- X.marge
      X.marge <- X.marge/sum(X.marge)
      X.mean <- apply(X * X.marge, 2, sum)
      traits[, k1:k2] <- X
      col.w[k1:k2] <- X.mean
    }
    attr(traits, "col.blocks") <- col.blocks
    attr(traits, "col.freq") <- col.w
    col.num <- factor(rep((1:length(col.blocks)), col.blocks))
    attr(traits, "col.num") <- col.num
    return(traits)
  }
  
  pca.new<-1:(rdt*m*nrow(traits))
  dim(pca.new)<-c(nrow(traits),(rdt*m))
  colnames(pca.new)<-c(rep(c(1:m),rdt))
  rownames(pca.new)<-rownames(traits)
  for (trial in 1:rdt){
    prep.pca(traits,col.blocks)->random.traits
    predict(trait.pca.object$pca,random.traits)[,1:m]->pca.new[,(((trial-1)*(m)+1):(m*trial))]
  }
  return(pca.new[,order(colnames(pca.new))])
}


### Functions to estimate functional diversity components (not to run by the user, just to make work other functions)

pol.area <- function(P) { 
    P<-as.matrix(P)
    length(unique(as.vector(P)))->u.points
    if (ncol(P)<2)  (if(u.points<2) return(0) else return(abs(max(P)-min(P))))
    if (u.points<3) return(0) else 
P[chull(P),]->P ## to calculate bounded area
n<-nrow(P)
x<-P[,1]; y<-P[,2] 
p1<-sum(x[1:(n-1)]*y[2:n])+x[n]*y[1] 
p2<-sum(x[2:n]*y[1:(n-1)])+x[1]*y[n] 
## abs() is used to avoid negative areas
return(if (u.points<4) max((abs(max(P[,1])-min(P[,1]))),(abs(max(P[,2])-min(P[,2])))) else abs(0.5*(p1-p2)))
}


int.area<-function(A,B) {
    A<-as.matrix(A)
    B<-as.matrix(B)
    poly.merge<-rbind(A,B)
    return(if (pol.area(A)+pol.area(B)-pol.area(poly.merge)<=0) 0 else pol.area(A)+pol.area(B)-pol.area(poly.merge))   
}


taxon.fric<-function(fpc,m,rdt){
  res<-list()
  if (m>1) combn(m,2)->axis.pair else matrix(c(1,1))->axis.pair
  if (m>1) ncol(axis.pair)->m.ncol else m.ncol<-1
  if (m>1) ind.nic<-rep(0,m.ncol) else ind.nic<-0
  tFRicent<-tFRirange<-matrix(rep(0,m*nrow(fpc)),ncol=m)
  
  if (m>1) tfric.data<-matrix(rep(0,m.ncol*nrow(fpc)),ncol=m.ncol) else tfric.data<-matrix(rep(0,m*nrow(fpc)),ncol=m)# matrix of in values per axis pair
  
  for (sp in 1:nrow(fpc)){
    for (p.row in 1:m.ncol){ ### obtaining the coordinates for each axis pair
      fpc[sp,(((axis.pair[1,p.row]-1)*ncol(fpc)/m)+1):(ncol(fpc)/m*(axis.pair[1,p.row]-1)+rdt)]->a
      if (m>1) fpc[sp,(((axis.pair[2,p.row]-1)*ncol(fpc)/m)+1):(ncol(fpc)/m*(axis.pair[2,p.row]-1)+rdt)]->b
      a<-t(a)
      if (m>1) b<-t(b)
      if (m>1) ind.poly<-c(a,b) else ind.poly<-a
      if (m>1) dim(ind.poly)<-c(rdt,2) else dim(ind.poly)<-c(rdt,1)
      pol.area(ind.poly)->tfric.data[sp,p.row]
      tFRicent[sp,axis.pair[1,p.row]]<-mean(a)
      if (max(a)<0 & min(a)<0) abs(min(a)-max(a))->tFRirange[sp,axis.pair[1,p.row]] else max(a)-min(a)->tFRirange[sp,axis.pair[1,p.row]]
      if (m>1)  tFRicent[sp,axis.pair[2,p.row]]<-mean(b)
      if (m>1) if (max(b)<0 & min(b)<0) abs(min(b)-max(b))->tFRirange[sp,axis.pair[2,p.row]] else max(b)-min(b)->tFRirange[sp,axis.pair[2,p.row]]
    }
  }
  res$tn.dim<-tfric.data
  res$tn.sum<-rowSums(tfric.data)
  res$cent<-tFRicent
  res$range<-tFRirange
  return(res)
}

mean.taxon.fric<-function(taxa,tfric.data){  ##This function permit community mean estimation from individual niche data
  res<-list()
  if (is.vector(tfric.data)==T) tfric.data<-matrix(tfric.data,ncol=1)
  tFRicol<-ncol(tfric.data)
  rowSums(taxa)->tax.ric
  tFRic<-matrix(rep(0,tFRicol*nrow(taxa)),ncol=tFRicol)
  
  for (i in 1:tFRicol){
    apply(taxa,1,function(x){(x*tfric.data[,i])})->temp.matrix
    rowSums(t(temp.matrix))/tax.ric->tFRic[,i]
  }
  res$tFRic<-as.numeric(tFRic)
  return(res$tFRic)
}

### Funcional redundancy
taxon.overlap<-function(fpc,m,rdt,tfric.data){
  res<-list()
  nr<-nrow(fpc) # taxon number in this dataset, taxa are disposed in columns
  abs.overlap<-rel.overlap<-matrix(nrow=nr,ncol=nr) ## creating a square matrix of overlapping
  if (m>1) combn(m,2)->axis.pair else matrix(c(1,1))->axis.pair ## a matrix containing the different combination of axis pairs
  if (m>1) ncol(axis.pair)->m.ncol else m.ncol<-1
  combn(nr,2)->taxon.pair ## a matrix containing the different combination of species pairs
  r.overlap<-a.overlap<-rep(NA, m.ncol) ## creates a new variable to store the results of the absolute and relative species overlap per axis pair
  for (taxon in 1:ncol(taxon.pair)){
    for (p.row in 1:m.ncol){
      
      fpc[taxon.pair[1,taxon],(((axis.pair[1,p.row]-1)*ncol(fpc)/m)+1):(ncol(fpc)/m*(axis.pair[1,p.row]-1)+rdt)]->x1
      if (m>1) fpc[taxon.pair[1,taxon],(((axis.pair[2,p.row]-1)*ncol(fpc)/m)+1):(ncol(fpc)/m*(axis.pair[2,p.row]-1)+rdt)]->y1
      fpc[taxon.pair[2,taxon],(((axis.pair[1,p.row]-1)*ncol(fpc)/m)+1):(ncol(fpc)/m*(axis.pair[1,p.row]-1)+rdt)]->x2
      if (m>1) fpc[taxon.pair[2,taxon],(((axis.pair[2,p.row]-1)*ncol(fpc)/m)+1):(ncol(fpc)/m*(axis.pair[2,p.row]-1)+rdt)]->y2
      
      if (m>1) poly1<-cbind(x1,y1) else poly1<-cbind(x1)
      if (m>1) poly2<-cbind(x2,y2) else poly2<-cbind(x2)
      
      int.area(poly1,poly2)->a.overlap[p.row]
      pol.area(poly1)+pol.area(poly2)->area
      if (a.overlap[p.row]==0) r.overlap[p.row]<-0 else 2*a.overlap[p.row]/area->r.overlap[p.row]
      
    }
    sum(a.overlap)->abs.overlap[taxon.pair[2,taxon],taxon.pair[1,taxon]] # lower triangle, storing data in overlapping triangular matrix, note that diagonal should be 1
    sum(a.overlap)->abs.overlap[taxon.pair[1,taxon],taxon.pair[2,taxon]] # upper triangle,
    diag(abs.overlap)<-tfric.data$tn.sum ## filling matrix diagonal with individual niches
    mean(r.overlap)->rel.overlap[taxon.pair[2,taxon],taxon.pair[1,taxon]]
    mean(r.overlap)->rel.overlap[taxon.pair[1,taxon],taxon.pair[2,taxon]]
    diag(rel.overlap)<-rep(1,nr)
  }
  res$abs.overlap<-abs.overlap
  res$rel.overlap<-rel.overlap
  rownames(res$abs.overlap)<-rownames(res$rel.overlap)<-rownames(fpc)
  colnames(res$abs.overlap)<-colnames(res$rel.overlap)<-rownames(fpc)
  return(res)
}

fr.sim<-function(taxa,nic.overlap){
  res<-list()
  abs.fredun.com<-1:nrow(taxa)
  rel.fredun.com<-1:nrow(taxa)
  for (r in 1:nrow(taxa)) {
    abs.fredun<-NULL
    rel.fredun<-NULL
    fredun.taxa<-taxa[r,which(taxa[r,]>0)] ## select species with abundance > 0
    which(taxa[r,]>0)->taxa.pos ## in which position can we find each taxon
      for (i in taxa.pos){
        sum(nic.overlap$abs.overlap[i,taxa.pos]*fredun.taxa)->abs.fredun[which(taxa.pos==i)] ### store the total overlap of each taxon
        sum(nic.overlap$rel.overlap[i,taxa.pos]*fredun.taxa)/length(fredun.taxa)->rel.fredun[which(taxa.pos==i)] ### store the total overlap of each taxon
      
      }
    sum(abs.fredun)->res$FR[r] # Functional redundancy
    mean(rel.fredun)->res$FSim[r] # Functional similarity
  }
  return(res)
}


fdis<-function(taxa,fpc,m,rdt) { ## taxa: fauna, fpc: functional principal components
  if (m>1) combn(m,2)->axis.pair else matrix(c(1,1))->axis.pair ## a matrix containing the different combination of axis pairs
  if (m>1) ncol(axis.pair)->m.ncol else m.ncol<-1
  com.fdis<-rep(0,m.ncol)
  apply(taxa,1,function(com.row) 
  {
  taxon.ric<-sum(decostand(com.row,method="pa"))

    # selecting the coordinates for the community of the axis p.row, including all random trials
    com.row[which(com.row>0)]->com.abun ## abundances of the community
    for (p.row in 1:m.ncol){
      fpc[which(com.row>0),(((axis.pair[1,p.row]-1)*ncol(fpc)/m)+1):(ncol(fpc)/m*(axis.pair[1,p.row]-1)+rdt)]->a
      if (m>1) fpc[which(com.row>0),(((axis.pair[2,p.row]-1)*ncol(fpc)/m)+1):(ncol(fpc)/m*(axis.pair[2,p.row]-1)+rdt)]->b
      
      axis.a<-axis.pair[1,p.row]
      if (m>1) axis.b<-axis.pair[2,p.row]
      
      # weighted centroids
      if (taxon.ric<2) mean(a)->cent.a else sum(apply(a,2,function(x){x*t(com.abun/(sum(com.abun)*rdt))}))->cent.a ## weighted by 
      if (m>1) if (taxon.ric<2) mean(b)->cent.b else sum(apply(b,2,function(x){x*t(com.abun/(sum(com.abun)*rdt))}))->cent.b
      
      (a-cent.a)^2->a.temp
      if (m>1) (b-cent.b)^2->b.temp else b.temp<-0
      a.temp+b.temp->sq.dist
      if (taxon.ric<2) mean(sq.dist)->w.fdis else sum(apply(sq.dist,2,function(x){x*t(com.abun/(sum(com.abun)*rdt))}))->w.fdis
      
    
    w.fdis->com.fdis[p.row]
  }
  return(FDis<-sum(com.fdis))
}
  )
}


fric<-function(taxa,fpc,m,rdt) {## taxa: fauna, fpc: functional principal components
  if (m>1) combn(m,2)->axis.pair else matrix(c(1,1))->axis.pair ## a matrix containing the different combination of axis pairs
  if (m>1) ncol(axis.pair)->m.ncol else m.ncol<-1
  apply(taxa,1,function(com.row) 
  {
    com.nic<-rep(0,m.ncol)
    for (p.row in 1:m.ncol){
      fpc[which(com.row>0),(((axis.pair[1,p.row]-1)*ncol(fpc)/m)+1):(ncol(fpc)/m*(axis.pair[1,p.row]-1)+rdt)]->a
      if (m>1) fpc[which(com.row>0),(((axis.pair[2,p.row]-1)*ncol(fpc)/m)+1):(ncol(fpc)/m*(axis.pair[2,p.row]-1)+rdt)]->b
      
      a<-as.vector(t(a))
      if (m>1) b<-as.vector(t(b))
      if (m>1) com.poly<-c(a,b) else com.poly<-a
      if (m>1) dim(com.poly)<-c(length(a),2)
      pol.area(com.poly)->com.nic[p.row]
    }
    return(FRic<-sum(com.nic))
  }
  )
}

### Function to determine the optimal number of randomizations (rdt):
det.rdt<-function(taxa,traits,calc.pca=T, m=2,test.rdt,plots=T){
  if (min(test.rdt)<3) stop("'rdt' must be greater than 2")
  res<-list() # generating an empty list of results
  if (calc.pca==T) {funcSpace(taxa,traits,m,max(test.rdt))->fpc; fpc->res$fpc} else fpc<-traits ## storing pca axes

  ### creating variables to store funcNiche results
  tFRic<-rep(0,length(test.rdt)*nrow(taxa))
  dim(tFRic)<-c(nrow(taxa),length(test.rdt))
  tFRic->FDis->FRic->FR->FSim
  
  for (i in test.rdt){
  cat("Random trials:",i,fill=T)  
    funcNiche(taxa,fpc,m,i)->results
    results$tFRic->tFRic[,which(test.rdt==i)]
    results$FDis->FDis[,which(test.rdt==i)]
    results$FRic->FRic[,which(test.rdt==i)]
    results$FR->FR[,which(test.rdt==i)]
    results$FSim->FSim[,which(test.rdt==i)]
  }
  res$tFRic<-tFRic
  res$FDis<-FDis
  res$FRic<-FRic
  res$FR<-FR
  res$FSim<-FSim
  
  ### the last row is the dependent variable of the plot
  cor(tFRic)[length(test.rdt),]->res$tFRic.cor
  cor(FDis)[length(test.rdt),]->res$FDis.cor
  cor(FRic)[length(test.rdt),]->res$FRic.cor
  cor(FR)[length(test.rdt),]->res$FR.cor
  cor(FSim)[length(test.rdt),]->res$FSim.cor
  
  res$rdt.table<-data.frame(tRic=res$tFRic.cor, FSim=res$FSim.cor, FRic=res$FRic.cor, FDis=res$FDis.cor, FR=res$FR.cor)
  row.names(res$rdt.table)<-test
  
  if (plots==F) return(res) else
 
  ### Plotting results
  par(mar = c(4, 6, 4, 1))
  par(mfrow=c(3,2))
  plot(test.rdt,colMeans(res$tFRic),col="grey",xlab="Number of random trials",ylab="Correlation coefficient", cex=1.25,cex.axis=1.5,cex.lab=1.75,cex.main=2,main="Mean individual niche")
  plot(test.rdt,colMeans(res$FDis),col="grey",xlab="Number of random trials",ylab="", cex=1.25,cex.axis=1.5,cex.lab=1.75,cex.main=2,main="Functional dispersion")
  plot(test.rdt,colMeans(res$FRic),col="grey",xlab="Number of random trials",ylab="Correlation coefficient", cex=1.25,cex.axis=1.5,cex.lab=1.75,cex.main=2,main="Community niche")
  plot(test.rdt,colMeans(res$FR),col="grey",xlab="Number of random trials",ylab="", cex=1.25,cex.axis=1.5,cex.lab=1.75,cex.main=2,main="Absolute redundancy")
  plot(test.rdt,colMeans(res$FSim),col="grey",xlab="Number of random trials",ylab="Correlation coefficient", cex=1.25,cex.axis=1.5,cex.lab=1.75,cex.main=2,main="Relative redundancy")
  return(res)
  }


### Functions to estimate Functional Diversity components:

funcNiche.pca<-function(taxa,traits,m=2,rdt=9) {
  
  if(is.null(taxa)==T) stop("Taxonomic matrix is required")
  if(!exists("rdt")==T) stop("You must enter a number of randomizations") else if(is.null(rdt)==T) stop("You must enter a number of randomizations")
  if(ncol(taxa)!=nrow(traits)) stop("Taxonomic and trait matrices must have the same taxon number")
  if(!exists("m")==T) stop("You must enter a number of functional space dimensions") else if(is.null(m)==T) stop("You must enter a number of functional space dimensions")
  
  funcSpace(taxa,traits,m,rdt)->fpc
  tfric.data<-taxon.fric(fpc,m,rdt)
  overlap.mat<-taxon.overlap(fpc,m,rdt,tfric.data)
  
  res<-list()
  res$tFRic<-mean.taxon.fric(taxa,tfric.data$tn.sum)
  cat("Mean taxon functional richness calculated",fill=T)
  res$FDis<-fdis(taxa,fpc,m,rdt)
  cat("Functional disimilarity calculated",fill=T)
  res$FRic<-fric(taxa,fpc,m,rdt)
  cat("Functional richness calculated",fill=T)
  fredun<-fr.sim(taxa,overlap.mat)
  fredun$FR->res$FR
  cat("Functional redundancy calculated",fill=T)
  fredun$FSim->res$FSim
  cat("Functional similarity calculated",fill=T)
  return(res)
}

funcNiche<-function(taxa,fpc,m=2,rdt=9) {
  
  if(is.null(taxa)==T) stop("Taxonomic matrix is required")
  if(is.null(fpc)==T) stop("Functional space is required")
  if(ncol(taxa)!=nrow(fpc)) stop("Taxonomic and trait matrices must have the same taxon number")
  
  tfric.data<-taxon.fric(fpc,m,rdt)
  overlap.mat<-taxon.overlap(fpc,m,rdt,tfric.data)
  
  res<-list()
  res$tFRic<-mean.taxon.fric(taxa,tfric.data$tn.sum)
  cat("Mean taxon functional richness calculated",fill=T)
  res$FDis<-fdis(taxa,fpc,m,rdt)
  cat("Functional disimilarity calculated",fill=T)
  res$FRic<-fric(taxa,fpc,m,rdt)
  cat("Functional richness calculated",fill=T)
  fredun<-fr.sim(taxa,overlap.mat)
  fredun$FR->res$FR
  cat("Functional redundancy calculated",fill=T)
  fredun$FSim->res$FSim
  cat("Functional similarity calculated",fill=T)
  return(res)
}

### Functions to estimate null models and test for significant differences between observed and simulated model paramters:

simul.func<-function(taxa,fpc,m=2,rdt=9,runs=99){
  
  if(is.null(taxa)==T) stop("Taxonomic matrix is required")
  if(is.null(fpc)==T) stop("Functional space is required")
  if(ncol(taxa)!=nrow(fpc)) stop("Taxonomic and trait matrices must have the same taxon number")
  
  res<-list()
  tfric.data<-taxon.fric(fpc,m,rdt)
  overlap.mat<-taxon.overlap(fpc,m,rdt,tfric.data)
  
  f1<-function(x){ # This function performs funcNiche analyses
    
    res<-list()
    res$tFRic<-mean.taxon.fric(x,tfric.data$tn.sum)
    cat("Mean taxon functional richness calculated",fill=T)
    res$FDis<-fdis(x,fpc,m,rdt)
    cat("Functional disimilarity calculated",fill=T)
    res$FRic<-fric(x,fpc,m,rdt)
    cat("Functional richness calculated",fill=T)
    fredun<-fr.sim(x,overlap.mat)
    fredun$FR->res$FR
    cat("Functional redundancy calculated",fill=T)
    fredun$FSim->res$FSim
    cat("Functional similarity calculated",fill=T)
    return(res)
    
  } 
  # permatfull is configured to fix row and column totals in an absence/presence taxonomic matrix
  sapply(permatfull(taxa,fixedmar="both",mtype="prab",times=runs)$perm,f1)->temp.res 
  matrix(unlist(temp.res),ncol=5*runs)->res.mat
  colnames(res.mat)<-rep(c("tn.s","fd.s","cn.s","fr.s","fs.s"),runs)
  res$tn.s<-res.mat[,seq(1,ncol(res.mat),by=5)]
  res$fd.s<-res.mat[,seq(2,ncol(res.mat),by=5)]
  res$cn.s<-res.mat[,seq(3,ncol(res.mat),by=5)]
  res$fr.s<-res.mat[,seq(4,ncol(res.mat),by=5)]
  res$fs.s<-res.mat[,seq(5,ncol(res.mat),by=5)]
  res$runs<-runs
  return(res)
}

test.simul<-function (taxa,real.res,sim.data,stressor){
  
  if(is.null(taxa)==T) stop("Taxonomic matrix is required")
  if(is.null(real.res)==T) stop("Observed functional data are required")
  if(is.null(sim.data)==T) stop("Simulated functional data are required")
  if(is.null(stressor)==T) stop("Stress data are required")

  res.sim<-list()
  tn1.m<-tn1.b1<-fd1.m<-fd1.b1<-cn1.m<-cn1.b1<-fr1.m<-fr1.b1<-fs1.m<-fs1.b1<-rep(NA,sim.data$runs)
  ric<-apply(decostand(taxa,method="pa"),1,sum)
  tn.1<-res.sim$tn.1<-(glm(real.res$tFRic~scale(stressor)))
  fd.1<-res.sim$fd.1<-(glm((real.res$FDis)~scale(stressor)))
  cn.1<-res.sim$cn.1<-(glm((real.res$FRic)~scale(stressor)))
  fr.1<-res.sim$fr.1<-(glm(log(real.res$FR+1)~scale(stressor)))
  fs.1<-res.sim$fs.1<-(glm(log(real.res$FSim[which(ric>1)]+1)~scale(stressor[which(ric>1)])))
  
  for (a in 1:sim.data$runs){
    tn.sim1<-(glm(sim.data$tn.s[,a]~scale(stressor)))
    fd.sim1<-(glm((sim.data$fd.s[,a])~scale(stressor)))
    cn.sim1<-(glm((sim.data$cn.s[,a])~scale(stressor)))
    fr.sim1<-(glm(log(sim.data$fr.s[,a]+1)~scale(stressor)))
    fs.sim1<-(glm(log(sim.data$fs.s[which(ric>1),a]+1)~scale(stressor[which(ric>1)])))
    
      
    tn.sim1$coef[1]->res.sim$tn1.m[a]->tn1.m[a]
    tn.sim1$coef[2]->res.sim$tn1.b1[a]->tn1.b1[a]
    fd.sim1$coef[1]->res.sim$fd1.m[a]->fd1.m[a]
    fd.sim1$coef[2]->res.sim$fd1.b1[a]->fd1.b1[a]
    cn.sim1$coef[1]->res.sim$cn1.m[a]->cn1.m[a]
    cn.sim1$coef[2]->res.sim$cn1.b1[a]->cn1.b1[a]
    fr.sim1$coef[1]->res.sim$fr1.m[a]->fr1.m[a]
    fr.sim1$coef[2]->res.sim$fr1.b1[a]->fr1.b1[a]
    fs.sim1$coef[1]->res.sim$fs1.m[a]->fs1.m[a]
    fs.sim1$coef[2]->res.sim$fs1.b1[a]->fs1.b1[a]
  }
  
  par(mar = c(5, 4.25, 2.5, 2))
  par(mfrow=c(5,2))
  
  ### Individual niche
  hist(tn1.m,lwd=4,main="tFRic",ylab= "Frequency", xlab="", cex=1.25,cex.axis=1.5,cex.lab=1.5,cex.main=2,xlim=c(min(c(tn1.m,tn.1$coef[1])),max(c(tn1.m,tn.1$coef[1]))))
  abline(v=tn.1$coef[1], lty=3, lwd=7, col="gray")
  
  hist(tn1.b1,lwd=4,main="",ylab= "", xlab="", cex=1.25,cex.axis=1.5,cex.lab=1.5,cex.main=2,,xlim=c(min(c(tn1.b1,tn.1$coef[2])),max(c(tn1.b1,tn.1$coef[2]))))
  abline(v=tn.1$coef[2], lty=3, lwd=7, col="gray")
  
  ## Functional similarity
  hist(fr1.m,lwd=4,main="FSim",ylab= "Frequency", xlab="", cex=1.25,cex.axis=1.5,cex.lab=1.5,cex.main=2,,xlim=c(min(c(fr1.m,fr.1$coef[1])),max(c(fr1.m,fr.1$coef[1]))))
  abline(v=fr.1$coef[1], lty=3, lwd=7, col="gray")
  
  hist(fr1.b1,lwd=4,main="",ylab="", xlab="", cex=1.25,cex.axis=1.5,cex.lab=1.5,cex.main=2,,xlim=c(min(c(fr1.b1,fr.1$coef[2])),max(c(fr1.b1,fr.1$coef[2]))))
  abline(v=fr.1$coef[2], lty=3, lwd=7, col="gray")
  
  ## Community niche
  hist(cn1.m,lwd=4,main="FRic",ylab= "Frequency", xlab="", cex=1.25,cex.axis=1.5,cex.lab=1.5,cex.main=2,xlim=c(min(c(cn1.m,cn.1$coef[1])),max(c(cn1.m,cn.1$coef[1]))))
  abline(v=cn.1$coef[1], lty=3, lwd=7, col="gray")
  
  hist(cn1.b1,lwd=4,main="",ylab= "", xlab="", cex=1.25,cex.axis=1.5,cex.lab=1.5,cex.main=2,xlim=c(min(c(cn1.b1,cn.1$coef[2])),max(c(cn1.b1,cn.1$coef[2]))))
  abline(v=cn.1$coef[2], lty=3, lwd=7, col="gray")
  
  ## Functional dispersion
  hist(fd1.m,lwd=4,main="FDis",ylab= "Frequency", xlab="", cex=1.25,cex.axis=1.5,cex.lab=1.5,cex.main=2,xlim=c(min(c(fd1.m,fd.1$coef[1])),max(c(fd1.m,fd.1$coef[1]))))
  abline(v=fd.1$coef[1], lty=3, lwd=7, col="gray")
  
  hist(fd1.b1,lwd=4,main="",ylab= "", xlab="", cex=1.25,cex.axis=1.5,cex.lab=1.5,cex.main=2,,xlim=c(min(c(fd1.b1,fd.1$coef[2])),max(c(fd1.b1,fd.1$coef[2]))))
  abline(v=fd.1$coef[2], lty=3, lwd=7, col="gray")
  
  ## Functional redundancy
  hist(fs1.m,lwd=4,main="FR",ylab= "Frequency", xlab="Intercept", cex=1.25,cex.axis=1.5,cex.lab=1.5,cex.main=2,,xlim=c(min(c(fs1.m,fs.1$coef[1])),max(c(fs1.m,fs.1$coef[1]))))
  abline(v=fs.1$coef[1], lty=3, lwd=7, col="gray")
  
  hist(fs1.b1,lwd=4,main="",ylab="", xlab="Slope", cex=1.25,cex.axis=1.5,cex.lab=1.5,cex.main=2,xlim=c(min(c(fs1.b1,fs.1$coef[2])),max(c(fs1.b1,fs.1$coef[2]))))
  abline(v=fs.1$coef[2], lty=3, lwd=7, col="gray")

  p.val<-function (null.dist, obs, alternative = c("two.sided", "less", "greater"))
  { 
    alternative <- match.arg(alternative)
    nsimul<-length(null.dist)
    z <- (obs - mean(null.dist))/sd(null.dist) # z-score
    pless <- sum(obs <= null.dist, na.rm = TRUE) # null values less or equal than real observation
    pmore <- sum(obs >= null.dist, na.rm = TRUE) # null values greater or equal than real observation
    
    p <- switch(alternative, two.sided = 2 * pmin(pless, pmore), less = pless, greater = pmore)
    p <- pmin(1, (p + 1)/(nsimul + 1)) 
    res<-list(z.score=z,p.value=p,alternative=alternative)
    res.table<-data.frame(Real.data=obs,z=z,mean=mean(null.dist),q2.5=quantile(null.dist,0.025),q50=quantile(null.dist,0.50),q97.5=quantile(null.dist,0.975), P.value=p)
    return(res.table)
  }

p.val(tn1.m,tn.1$coef[1])$P.value->pv.tn.m
p.val(tn1.b1,tn.1$coef[2])$P.value->pv.tn.b1
p.val(fd1.m,fd.1$coef[1])$P.value->pv.fd.m
p.val(fd1.b1,fd.1$coef[2])$P.value->pv.fd.b1
p.val(cn1.m,cn.1$coef[1])$P.value->pv.cn.m
p.val(cn1.b1,cn.1$coef[2])$P.value->pv.cn.b1
p.val(fr1.m,fr.1$coef[1])$P.value->pv.fr.m
p.val(fr1.b1,fr.1$coef[2])$P.value->pv.fr.b1
p.val(fs1.m,fs.1$coef[1])$P.value->pv.fs.m
p.val(fs1.b1,fs.1$coef[2])$P.value->pv.fs.b1
  
p.val(tn1.m,tn.1$coef[1])$z->z.tn.m
p.val(tn1.b1,tn.1$coef[2])$z->z.tn.b1
p.val(fd1.m,fd.1$coef[1])$z->z.fd.m
p.val(fd1.b1,fd.1$coef[2])$z->z.fd.b1
p.val(cn1.m,cn.1$coef[1])$z->z.cn.m
p.val(cn1.b1,cn.1$coef[2])$z->z.cn.b1
p.val(fr1.m,fr.1$coef[1])$z->z.fr.m
p.val(fr1.b1,fr.1$coef[2])$z->z.fr.b1
p.val(fs1.m,fs.1$coef[1])$z->z.fs.m
p.val(fs1.b1,fs.1$coef[2])$z->z.fs.b1

res.sim$pv.test<-data.frame(tFRic=c(pv.tn.m,z.tn.m,pv.tn.b1,z.tn.b1),FSim=c(pv.fs.m,z.fs.m,pv.fs.b1,z.fs.b1),FDis=c(pv.fd.m,z.fd.m,pv.fd.b1,z.fd.b1),FRic=c(pv.cn.m,z.cn.m,pv.cn.b1,z.cn.b1),FR=c(pv.fr.m,z.fr.m,pv.fr.b1,z.fr.b1))  
row.names(res.sim$pv.test)<-c("Intercept","z-intercept","Slope","z-slope")
  
return(res.sim)
}
