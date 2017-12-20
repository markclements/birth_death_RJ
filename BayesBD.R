MCMC.BD<-function(phy,diversity,r=0.5,a=0,iter=1000,sample=1000,lambda=0.5){

rootnode <- length(phy$tip.label) + 1
edge<-phy$edge
stem.ages<-branching.times(phy)
edge<-cbind(edge,stem.ages[phy$edge[,1]-Ntip(phy)])
#diversity <- diversity# vector of terminal diversity
edge<-cbind(edge,phy$edge.length)
diversity[phy$edge[,2]] # assign diversity to termainals, NAs are for internals. 
edge<-cbind(edge,diversity[phy$edge[,2]])
edge<-cbind(edge,numeric(length(edge[,1])))#,edge<-cbind(edge,numeric(length(edge[,1]))))
edge<-edge[order(edge[,1]),] # reorder matrix descending by node 
#edge[,6]<-runif(phy$Nnode,0.0001,8)
#edge[,6]<-rgamma(phy$Nnode,r)
edge[,6]<-r
a=a

likelihood<-function(edge){
int <- edge[edge[, 2] > rootnode, ]
term <- edge[edge[, 2] < rootnode, ]
comp.term<-numeric(length(term[,1]))	
for (i in 1:length(comp.term)){
 		beta <-(expm1(term[i,6]*term[i,4]))/(exp(term[i,6]*term[i,4])-a) 
 		comp.term[i]<-(1-beta)*((beta)^(term[i,5]-1)) 
 	}
lik.term<-sum(log(comp.term))


comp.int<-numeric(length(int[,1]))
for (i  in 1:length(comp.int)){
   	comp.int[i]<-int[i,6]*((1-a*exp(-int[i,6]*int[i,3]))^(-1))*(exp(-int[i,6]*int[i,4]))  	 }
 lik.int<-sum(log(comp.int))

lik<-lik.int+lik.term
lik
}

logLik<-likelihood(edge)
#return(logLik)
#}


MCMC<-function(logLik,edge){
newedge<-edge
newnode<-sample(newedge[,1],1)
U<-runif(1,0,1)
#newedge[which(newedge[,1]==newnode),6]<-newedge[which(newedge[,1]==newnode),6]*exp(U-0.5)
newedge[which(newedge[,1]==newnode),6]<-newedge[which(newedge[,1]==newnode),6]+(U-0.5)*lambda
if (newedge[which(newedge[,1]==newnode),6][-1]<0)
	newedge[which(newedge[,1]==newnode),6]<-abs(newedge[which(newedge[,1]==newnode),6])
logLik<-logLik
newlogLik<-likelihood(newedge)
lr <- newlogLik - logLik
prior.ratio <- exp(edge[which(edge[,1]==newnode),6]-newedge[which(newedge[,1]==newnode),6])

ratio <- exp(lr) * prior.ratio[-1] * 1
#ratio <- exp(lr) * prior.ratio[-1] * U 

  if(runif(1,0,1) < ratio)
    return(list(newedge,newlogLik,1))
	 
  else
    return(list(edge,logLik,0))
	
	
}

#BEGIN CALCULATION
out.edge <- vector("list",iter)
out.logLik <- vector("list",iter)
out.accept <- vector("list",iter)
  for(i in (1:iter)) {
  	for (j in (1:sample)){
  step<-MCMC(logLik,edge)
  edge<-step[[1]]
  logLik<-step[[2]]	
  accept<-step[[3]]
	}
print(list(i*j,logLik,edge))
out.edge[[i]]<-edge[,6]
out.logLik[[i]]<-logLik
out.accept[[i]]<-c(accept,i*j)
}
#return(list(edge,logLik,edge))
list(iter=out.accept,logLik=out.logLik,edge=out.edge)
       
}






find.clade<-function(phy,nod){
	if (is.null(attr(phy, "order")) || attr(phy, "order") == "pruningwise") 
        phy <- reorder(phy,"c")
	edge<-phy$edge
	start<-which(edge[,1]==nod)
	tmp<-which(edge[,1]<edge[start[2],1])
	end<-tmp[tmp>start[2]]
	clade<-c(start[1],end[1]-1)
	if (is.na(clade[2])) clade[2]<-length(edge[,1])
	clade
	}








#BEGIN CALCULATION

  for(i in (1:nstep + 1)) {

  #progress bar
  if(i %% 100 == 0){
   z<-i/nstep
   zt<-(i-100)/(nstep)
   polygon(c(zt,zt,z,z), c(1,0,0,1), col="black")

    }

  # calculate jump probabilities without given lamda
  if(method.prior.changepoints=="hierarchical"){
    prior[1]<-rgamma(1,shape=gamma.shape,scale=gamma.scale)
    jump.prob <- matrix(ncol=4,nrow=prior[4]+1)
    p <- dpois(0:prior[4],prior[1])/ppois(prior[4]+1,prior[1])
    bk <- c(p[-1]/p[-length(p)],0)
    bk[bk > 1] <- 1
    dk <- c(0,p[-length(p)]/p[-1])
    dk[dk > 1] <- 1
    mx <- max(bk+dk)
    bk <- bk/mx*0.9
    dk <- dk/mx*0.9
    bk[is.na(bk)]<-0   # added
    dk[is.na(dk)]<-0   # added
    jump.prob[,3] <- bk
    jump.prob[,4] <- dk
    jump.prob[1,2] <- 0
    jump.prob[1,1] <- 1-bk[1]-dk[1]
    jump.prob[-1,1] <- jump.prob[-1,2] <-
    (1-jump.prob[-1,3]-jump.prob[-1,4])/2
  }

    # determine what type of jump to make
    wh <- sample(1:4,1,prob=jump.prob[length(h)-1,])

    if (i %% thinning == 0& i>burn.in) {save.steptype[[count.i]] <- wh}

    if(wh==1) {
      step <- ht.move(data,pos,h,curloglik,prior, b.lin, sk1, ci, prior.height.mean, prior.height.var)
      h <- step[[1]]
      curloglik <- step[[2]]
      if(i%%thinning==0 & i>burn.in){
         save.pos[[count.i]]<-pos
         save.h[[count.i]]<-h
         save.loglik[[count.i]]<-step[[2]]
         save.accept[[count.i]]<-step[[3]]
         }
    }
    else if(wh==2) {
      step <- pos.move(data,pos,h,curloglik, b.lin,sk1,ci)
      pos <- step[[1]]
      curloglik <- step[[2]]
      if(i%%thinning==0 & i>burn.in){
          save.pos[[count.i]]<-pos
          save.h[[count.i]]<-h
          save.loglik[[count.i]]<-step[[2]]
          save.accept[[count.i]]<-step[[3]]
          }
    }
    else if(wh==3) {
      step <- birth.step(data,pos,h,curloglik,prior,jump.prob, b.lin, sk1, ci, prior.height.mean, prior.height.var)
      pos <- step[[1]]
      h <- step[[2]]
      curloglik <- step[[3]]
      if(i%%thinning==0 & i>burn.in){
         save.pos[[count.i]]<-pos
         save.h[[count.i]]<-h
         save.loglik[[count.i]]<-step[[3]]
         save.accept[[count.i]]<-step[[4]]
         }
    }
    else {
      step <- death.step(data,pos,h,curloglik,prior,jump.prob, b.lin, sk1, ci, prior.height.mean, prior.height.var)
      pos <- step[[1]]
      h <- step[[2]]
      curloglik <- step[[3]]
      if(i%%thinning==0 & i>burn.in){
         save.pos[[count.i]]<-pos
         save.h[[count.i]]<-h
         save.loglik[[count.i]]<-step[[3]]
         save.accept[[count.i]]<-step[[4]]
         }
    }
    if (i %% thinning == 0& i>burn.in) {count.i<-count.i+1}
  }

fitNDR_1rate(phy,combined=FALSE)



# old stuff below!!
comp1<-numeric(length(edge[,1]))
comp2<-numeric(length(edge[,1]))
for (i in 1:length(comp1)){ # log like of taxanomic data
	if (is.na(edge[i,3])) {
		comp1[i]<-0
		comp2[i]<-0
		}
	else {
		beta <-(exp(edge[i,7]*edge[i,6])-1)/(exp(edge[i,7]*edge[i,6])-a) 
		comp1[i]<- log(1-beta) 
		comp2[i]<- (edge[i,3]-1)*log(beta)
		}
	comp<-sum(comp1)+sum(comp2)
	}
	
#simpler version of like of taxanomic data, maybe?

comp<-numeric(length(edge[,1]))	
for (i in 1:length(comp)){
if (is.na(edge[i,3])) comp[i]<-0
	else {
 		beta <-(exp(edge[i,7]*edge[i,6])-1)/(exp(edge[i,7]*edge[i,6])-a) 
 		comp[i]<-log(1-beta) + (edge[i,3]-1)*log(beta)
 		}
 	comp
 	}
 comp	
 
 Nint<-phy$Nnode
 for (i  in seq(from = 1, by = 2, length.out = length(comp1))){
 	comp[i]<-Nint*log(edge[i,7])-edge[i,7]*edge[i,6]-(log(1-a*(exp(-1*edge[i,7]*edge[i,5]))))
  	 }
  	 comp<-comp[comp!=0]
     lik<-sum(log(comp))
  	 
  	 
  	 
fitNDR_1rate(phy,combined=FALSE)