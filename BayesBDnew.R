MCMC.BD<-function(phy,diversity,root.rate=0.5,a=0,e=3,iter=1000,sample=1000,lambda=0.5){

require(ape)
if (is.null(attr(phy, "order")) || attr(tree, "order") == "cladewise")
phy<- reorder(phy,"c") 
rootnode <- length(phy$tip.label) + 1
edge<-phy$edge
stem.ages<-branching.times(phy)
edge<-cbind(edge,stem.ages[phy$edge[,1]-Ntip(phy)])
edge<-cbind(edge,phy$edge.length)
edge<-cbind(edge,diversity[phy$edge[,2]])# assign diversity to termainals, NAs are for internals. 
edge<-cbind(edge,NA) # rate column
edge<-cbind(edge,NA) # node split column
edge<-cbind(edge,0) # rate group column
dimnames(edge)<-NULL
rate.vect<-numeric(1)

#initiallize
edge<-rate.calc(edge=edge,root.rate=root.rate,rate.vect=rate.vect)
logLik<-likelihood(edge)

#BEGIN CALCULATION
out.edge <- vector("list",iter)
out.logLik <- vector("list",iter)
out.rates <- vector("list",iter)
  for(i in (1:iter)) {
  	for (j in (1:sample)){
  step<-MCMC(logLik=logLik,edge=edge,root.rate=root.rate, a=a, e=e, rate.vect=rate.vect)
  edge<-step[[1]]
  logLik<-step[[2]]	
  rate.vect<-step[[3]]
  root.rate<-step[[4]]	
	}
#cat(list(i*j,logLik,edge,rates))
out.edge[[i]]<-cbind(edge[,6],edge[,8])
out.logLik[[i]]<-logLik
out.rates[[i]]<-c(rate.vect,i*j)
}
#return(list(edge,logLik,edge))
list(iter=out.rates,logLik=out.logLik,edge=out.edge)
       
}

MCMC<-function(logLik,edge,root.rate,a,e,rate.vect){
N <-length(which(edge[,7]==1))
z<-sample(1:4,1,prob=c(5,5,60,30))
#z<-sample(1:2,1)

cat(root.rate)

if (N==0 && z==2){ #have to split
	x<-sample(1:length(edge[,1]),1)
	rate.vect.new<-rexp(1)
	newedge<-split(edge,x)
	newedge<-rate.calc(newedge,root.rate,rate.vect.new)
	newlogLik<-likelihood(newedge)
	HP<-e/2
	accept=1
	lr <- newlogLik/logLik
	ratio <- lr*HP
	cat(accept," ",log(logLik)," ", log(newlogLik)," ",N,"\n")
  if (runif(1,0,1) < ratio) return(list(newedge,newlogLik,rate.vect.new,root.rate))
  else return(list(edge,logLik,rate.vect,root.rate))
	}

else if (N==length(edge[,1]) && z==1){#must merge
	x<-sample(1:length(edge[,1]),1)
	rate.vect.new<-rate.vect[-x]
	newedge<-merge(edge,x)
	newedge<-rate.calc(newedge,root.rate,rate.vect.new)
	newlogLik<-likelihood(newedge)
	HP<-N/e
	accept=2
	lr <- newlogLik/logLik
	ratio <- lr*HP
	cat(accept," ",log(logLik)," ", log(newlogLik)," ",N,"\n")
  if (runif(1,0,1) < ratio) return(list(newedge,newlogLik,rate.vect.new,root.rate))
  else return(list(edge,logLik,rate.vect,root.rate))
  	}

else if (z==1 && N>0 && N<length(edge[,1])){# regular merge
	if (N==1){ 
		x<-which(edge[,7]==1)
		HP<-2/e
		}
	else { x<-sample(which(edge[,7]==1),1)
		HP<-N/e
		}
	s<-edge[x,8]
	rate.vect.new<-rate.vect[-s]
	newedge<-merge(edge,x)
	newedge<-rate.calc(newedge,root.rate,rate.vect.new)
	newlogLik<-likelihood(newedge)
	accept=3
	lr <- newlogLik/logLik
	ratio <- lr*HP
	cat(accept," ",log(logLik)," ", log(newlogLik)," ",N,"\n")
  if (runif(1,0,1) < ratio) return(list(newedge,newlogLik,rate.vect.new,root.rate))
  else return(list(edge,logLik,rate.vect,root.rate))
	}

else if (z==2 && N>0 && N<length(edge[,1])){# regular split
	x<-sample(which(is.na(edge[,7])),1)
	rate.vect.new<-c(rate.vect,rexp(1))
	newedge<-split(edge,x)	
	newedge<-rate.calc(newedge,root.rate,rate.vect.new)
	newlogLik<-likelihood(newedge)
	HP<-e/(N+1)
	accept=4
	lr <- newlogLik/logLik
	ratio <- lr*HP
	cat(accept," ",log(logLik)," ", log(newlogLik)," ",N,"\n")
  if(runif(1,0,1) < ratio) return(list(newedge,newlogLik,rate.vect.new,root.rate))
  else return(list(edge,logLik,rate.vect,root.rate))
	}

else if (z==4 && N>0){
	rm<-rexp(1)
	if (N==1) x<-which(edge[,7]==1)
	else x<-sample(which(edge[,7]==1),1)
	s<-edge[x,8]
	rate.vect.new<-rate.vect
	rate.vect.new[s]<-rm
	newedge<-rate.calc(edge,root.rate,rate.vect.new)
	newlogLik<-likelihood(newedge)
	accept=5
	lr <- newlogLik/logLik
	ratio <- lr
	cat(accept," ",log(logLik)," ", log(newlogLik)," ",N,"\n")
  if(runif(1,0,1) < ratio) return(list(newedge,newlogLik,rate.vect.new,root.rate))
  else return(list(edge,logLik,rate.vect,root.rate))
	}

else {
	rr<-root.rate + runif(1,-0.3,0.3)
	if (rr<=0) new.root.rate<-0+(0-rr)
	else if (rr>1) new.root.rate <-1-(rr-1)
	else new.root.rate<-rr
	newedge<-rate.calc(edge,new.root.rate,rate.vect)
	newlogLik<-likelihood(newedge)
	accept=6
	lr <- newlogLik/logLik
	ratio <- lr
	cat(accept," ",log(logLik)," ", log(newlogLik)," ",N,"\n")
  if(runif(1,0,1) < ratio) return(list(newedge,newlogLik,rate.vect,new.root.rate))
  else return(list(edge,logLik,rate.vect,root.rate))
	}
}


likelihood<-function(edge,a=0){
lik<-numeric(length(edge[,1]))
for (i in 1:length(edge[,1])){
 		
 	if (!is.na(edge[i,5])){
 		beta <-(expm1(edge[i,6]*edge[i,4]))/(exp(edge[i,6]*edge[i,4])-a) 
 		lik[i]<-(1-beta)*((beta)^(edge[i,5]-1)) 
 	}
	else {
   		lik[i]<-edge[i,6]*((1-a*exp(-edge[i,6]*edge[i,3]))^(-1))*(exp(-edge[i,6]*edge[i,4]))  	}
	#cat(lik[i],"\n")
}
#lik<-sum(log(lik))
lik<-prod(lik)
lik
}




#RJ code is below

find.clade<-function(edge,x){
	if (!is.na(edge[x,5])) clade<-x
	else{
	nod<-edge[x,2]	
	start<-which(edge[,1]==nod)
	tmp<-which(edge[,1]<edge[start[2],1])
	end<-tmp[tmp>start[2]]
	clade<-c(x[1],end[1]-1)
	if (is.na(clade[2])) clade[2]<-length(edge[,1])
	}
	clade
	}


split<-function(edge,x){
	if (!is.na(edge[x,7])) {
		stop("that lineage is already split, choose another")
		}
	edge[x,7]<-1
	tmp<-find.clade(edge,x)
	group<-max(unique(edge[,8]))+1
	oldgroup<-edge[x,8]
	if (length(tmp)==1) edge[x,8]<-group
	else for (i in tmp[1]:tmp[2]){
		if (edge[i,8]!=oldgroup) next
		edge[i,8]<-group
		}
	edge	
	}

merge<- function(edge,x){
	if (is.na(edge[x,7])) {
		stop("that lineage is not split, choose another")
		}
	edge[x,7]<-NA
	tmp<-find.clade(edge,x)
	oldgroup<-edge[x,8]
	d<-which(edge[,2]==edge[x,1])
	if (length(d)==0) anc.group<-0
	else anc.group<-edge[d,8]
	if (length(tmp)==1) edge[x,8]<-anc.group
	else for (i in tmp[1]:tmp[2]){
		if (edge[i,8]!=oldgroup) next
		edge[i,8]<-anc.group
		}
	edge[edge[,8]>oldgroup,8]<-edge[edge[,8]>oldgroup,8]-1
	edge	
	}



rate.calc<-function(edge,root.rate,rate.vect){
	for (i in 1:length(edge[,1])){
		d<-which(edge[,2]==edge[i,1])
		if (is.na(edge[i,7]) && length(d)==0) edge[i,6]<-root.rate
		else if (!is.na(edge[i,7]) && length(d)==0) edge[i,6]<-root.rate*rate.vect[edge[i,8]]
		else if (is.na(edge[i,7]) && length(d)==1) edge[i,6]<-edge[d,6]
		else if (!is.na(edge[i,7]) && length(d)==1) edge[i,6]<-edge[d,6]*rate.vect[edge[i,8]]
		}
	edge
	}





fitNDR_1rate(phy,combined=FALSE)



# Dirichlet process prior method below....did not work well

MCMC<-function(edge,iter=1000,sample=1000,r=0.5,alpha=5,k=5){
out.edge <- vector("list",iter)
out.logLik <- vector("list",iter)
out.accept <- vector("list",iter)
	for (j in 1:iter){
		for(i in 1:sample){
	lineage<-sample(1:length(edge[,1]),1)
	c<-cbind(Gibbs.k(edge,lineage,r),Gibbs.n(edge,lineage,r,alpha,k))
	edge[lineage,6]<-c[2,which(c[1,]==max(c[1,]))]	#cat(lineage," ",edge[lineage,6],"\n")
	}
	out.edge[[j]]<-edge[,6]
}
list(edge=out.edge)
}


Gibbs.k<-function(edge,lineage,r=r){
		edge[lineage,6]<-NA
		edge1<-edge
		k.rates<-unique(edge[,6])
		k.rates<-k.rates[!is.na(k.rates)]
		cond.lik<-length(k.rates)
		for (i in 1:length(k.rates)){
			edge[lineage,6]<-k.rates[i]
			cond.lik[i]<-likelihood(edge,r)
			cond.lik[i]<-cond.lik[i]*length(which(edge1[,6]==k.rates[i]))
		#cat(cond.lik[i],k.rates[i],length(which(edge1[,6]==k.rates[i])),"\n")	
	}	
#index<-which(cond.lik==max(cond.lik))	
#new.rate<-k.rates[index]
#new.rate<-c(new.rate, max(cond.lik))
new.rate<-rbind(cond.lik,k.rates)
return(new.rate)	
}


Gibbs.n<-function(edge,lineage,r=r,alpha=alpha,k=k){
		n.rates<-rexp(k)
		cond.lik<-length(n.rates)
		for (i in 1:length(n.rates)){
			edge[lineage,6]<-n.rates[i]
			cond.lik[i]<-likelihood(edge)
			cond.lik[i]<-cond.lik[i]*(alpha/k)
			#cat(cond.lik[i],n.rates[i],alpha/k,"\n")
			}
#index<-which(cond.lik==max(cond.lik))	
#new.rate<-n.rates[index]
#new.rate<-c(new.rate, max(cond.lik))
new.rate<-rbind(cond.lik,n.rates)
return(new.rate)	
}



# seperate likelihood calcs for int/ext pertions of tree
int.lik<-function(edge,rate,a=0){
	eps<-a
	r<-rate
	int<-edge[which(is.na(edge[,5])),]
	nint <- nrow(int)
	if (is.null(nint)){
	cond.lik<-(1*log(r)-r*int[4]-log(1-(eps*exp(-r*int[3]))))
	}
	else
	cond.lik<-(nint*log(r)-r*sum(int[1:nint, 4])-sum(log(1-(eps*exp(-r*int[1:nint, 3])))))
	return(cond.lik)
	}

ext.lik<-function(edge,rate,a=0){		
	eps<-a
	r<-rate
	term<-edge[which(!is.na(edge[,5])),]
	nterm<-nrow(term)
	betaF<-function(r,t1) {
        xf<-(exp(r*t1)-1)/(exp(r*t1)-eps)
        xf
    }
    if (is.null(nterm)){
    cond.lik<-(sum(log(1-betaF(r,term[4])))+sum((term[5]-1)*log(betaF(r,term[4]))))
	return(cond.lik)
	}
	
    else
	cond.lik<-(sum(log(1-betaF(r,term[1:nterm,4])))+sum((term[1:nterm,5]-1)*log(betaF(r,term[1:nterm, 4]))))
	return(cond.lik)
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

> death.step
function (data, pos, h, curloglik, prior, jump.prob, b.lin, sk1, 
    ci, prior.height.mean, prior.height.var) 
{
    if (length(pos) == 3) 
        j <- 2
    else j <- sample(2:(length(pos) - 1), 1)
    left <- pos[j - 1]
    right <- pos[j + 1]
    prior.mean <- prior.height.mean(pos[j])
    prior.var <- prior.height.var(pos[j])
    prior[3] <- prior.mean/prior.var
    prior[2] <- (prior.mean^2)/prior.var
    h.left <- h[j - 1]
    h.right <- h[j + 1]
    newheight <- (((pos[j] - left)/(right - left)) * (h.right - 
        h.left) + h.left)
    k <- length(pos) - 3
    L <- max(pos)
    prior.logratio <- log(k + 1) - log(prior[1]) - log(2 * (k + 
        1) * (2 * k + 3)) + 2 * log(L) - log(pos[j] - left) - 
        log(right - pos[j]) + log(right - left) - prior[2] * 
        log(prior[3]) + lgamma(prior[2]) - (prior[2] - 1) * log(newheight) - 
        prior[3] * (newheight)
    proposal.ratio <- (k + 1) * jump.prob[k + 1, 3]/jump.prob[k + 
        2, 4]/L
    jacobian <- ((pos[j] - left)/(right - left)) * (h[j + 1] - 
        h[j - 1]) + h[j - 1]
    newpos <- pos[-j]
    newh <- h[-j]
    newloglik <- loglik(data, newpos, newh, b.lin, sk1, ci)
    lr <- newloglik - curloglik
    ratio <- exp(lr + prior.logratio) * proposal.ratio * (jacobian^(-1))
    if (runif(1, 0, 1) < ratio) 
        return(list(newpos, newh, newloglik, 1))
    else return(list(pos, h, curloglik, 0))
}
<environment: namespace:ape>




fitNDR_1rate(phy,combined=FALSE)



# old stuff below!!

#edge<-sort(unique(edge[,3]),decreasing=TRUE) #branching time for internal nodes
	#cond.lik<-N*log(r)-r*(sum(abs(diff(edge))))-sum(log(1-a*exp(-r*edge)))
	#for (i in 1:length(edge)-1){		
	#		n<-i+1 # this works
	#		d<-n*r*exp(-n*r*(edge[i]-edge[i+1]))
	#		b<-(1-a*exp(-r*edge[i+1])^(n-1))
	#		c<-(1-a*exp(-r*edge[i]))^n
	#		cond.lik[i] <- d*(b/c)		
	#		}
	#return(sum(log(na.omit(cond.lik))))




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

getLambda.internal
function (zmat, rootnode, rbounds, para = 0.01, eps, combined = TRUE) 
{
    int <- zmat[zmat[, 2] > rootnode, ]
    term <- zmat[zmat[, 2] < rootnode, ]
    nint <- nrow(int)
    nterm <- nrow(term)
    betaF <- function(r, t1) {
        xf <- (exp(r * t1) - 1)/(exp(r * t1) - eps)
        xf
    }
    Lfunc_tax <- function(p) {
        r <- p
        (sum(log(1 - betaF(r, term[1:nterm, 4]))) + sum((term[1:nterm, 
            5] - 1) * log(betaF(r, term[1:nterm, 4]))))
    }
    Lfunc_phy <- function(p) {
        r <- p
        (nint * log(r) - r * sum(int[1:nint, 4]) - sum(log(1 - 
            (eps * exp(-r * int[1:nint, 3])))))
    }
    Lfunc_comb <- function(p) {
        r <- p
        (sum(log(1 - betaF(r, term[1:nterm, 4]))) + sum((term[1:nterm, 
            5] - 1) * log(betaF(r, term[1:nterm, 4]))) + nint * 
            log(r) - r * sum(int[1:nint, 4]) - sum(log(1 - (eps * 
            exp(-r * int[1:nint, 3])))))
    }
    res <- list()
    if (combined == TRUE) {
        if (nrow(int) == 0) 
            tempres <- optimize(Lfunc_tax, interval = rbounds, 
                maximum = TRUE)
        else tempres <- optimize(Lfunc_comb, interval = rbounds, 
            maximum = TRUE)
    }
    else {
        tempres <- optimize(Lfunc_tax, interval = rbounds, maximum = TRUE)
    }
    res$LH <- tempres$objective
    res$lambda <- tempres$maximum/(1 - eps)
    res$r <- tempres$maximum
    res$eps <- eps
    res <- as.data.frame(res)
    return(res)
}
> splitEdgeMatrix
function (phy, node) 
{
    x <- branching.times(phy)
    rootnode <- length(phy$tip.label) + 1
    phy$tag <- rep(1, nrow(phy$edge))
    if (node >= rootnode) {
        node.desc <- node
        pos <- 1
        phy$tag[phy$edge[, 1] == node.desc[1]] <- 2
        while (pos != (length(node.desc) + 1)) {
            temp <- node.sons(phy, node.desc[pos])
            temp <- temp[temp > rootnode]
            for (k in 1:length(temp)) {
                phy$tag[phy$edge[, 1] == temp[k]] <- 2
            }
            node.desc <- c(node.desc, temp)
            pos <- pos + 1
        }
    }
    else if (node > 0) 
        phy$tag[phy$edge[, 2] == node] <- 2
    z <- cbind(phy$edge, gsr(phy$edge[, 1], names(x), x), phy$edge.length, 
        phy$phenotype, phy$tag)
    z <- matrix(as.numeric(z), dim(z))
    z <- as.data.frame(z)
    return(z)
}
<environment: namespace:laser>


Descendants
function (x, node, type = c("tips", "children", "all")) 
{
    type <- match.arg(type)
    if (type == "children") 
        return(Children(x, node))
    desc = function(x, node, type) {
        isInternal = logical(max(x$edge))
        isInternal[unique(x$edge[, 1])] = TRUE
        if (!isInternal[node]) 
            return(node)
        ch = allChildren(x)
        res = NULL
        while (length(node) > 0) {
            tmp = unlist(ch[node])
            res = c(res, tmp)
            node = tmp[isInternal[tmp]]
        }
        if (type == "tips") 
            return(res[!isInternal[res]])
        res
    }
    if (length(node) > 1) 
        return(lapply(node, desc, x = x, type = type))
    desc(x, node, type)
}
<environment: namespace:phangorn>


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




find.clade<-function(edge,nod){
	start<-which(edge[,1]==nod)
	tmp<-which(edge[,1]<edge[start[2],1])
	end<-tmp[tmp>start[2]]
	clade<-c(start[1],end[1]-1)
	if (is.na(clade[2])) clade[2]<-length(edge[,1])
	clade
	}
	

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

