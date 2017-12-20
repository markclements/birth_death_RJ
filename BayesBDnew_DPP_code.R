MCMC.BD<-function(phy,diversity,root.rate=0.01,a=0,e=3,ngen=1000000,sfreq=100,lam=0.01,expl=0.01){
	require(ape)
	if (is.null(attr(phy, "order")) || attr(phy, "order") == "cladewise")
	phy<- reorder.phylo(phy,"c") 
	
	edge<-phy$edge
	stem.ages<-branching.times(phy)
	edge<-cbind(edge,stem.ages[phy$edge[,1]-Ntip(phy)]) # birth time xi =3
	edge<-cbind(edge,phy$edge.length) # branch length or waiting time, ti = 4
	edge<-cbind(edge,diversity[phy$edge[,2]])# assign diversity to termainals, NAs are for internals. =5
	edge<-cbind(edge,NA) # div rate column = 6
	edge<-cbind(edge,NA) # div node split column =7
	edge<-cbind(edge,0) # div rate group column = 8
	dimnames(edge)<-NULL
	rate.vect<-numeric(1)
#_______________________________________________________________________________



#______________________________________________________
# Main LOG likelihood code 
		likelihood<-function(edge,a){
		lik<-numeric(length(edge[,1]))
		for (i in 1:length(edge[,1])){
				
			if (!is.na(edge[i,5])){
				beta <-(expm1(edge[i,6]*edge[i,4]))/(exp(edge[i,6]*edge[i,4])-a) #taxonomic likelihood
				lik[i]<-log(1-beta) + ((edge[i,5]-1)* log (beta)) 
			}
			else {
				lik[i]<-log(edge[i,6])-(edge[i,6]*edge[i,4])-log(1-a*exp(-edge[i,6]*edge[i,3])) 
				#lik[i]<-log(edge[i,6]*((1-a*exp(-edge[i,6]*edge[i,3]))^(-1))*(exp(-edge[i,6]*edge[i,4]))) 
				}
			#cat(lik[i],"\n")
		}
		#cat(lik,"\n", edge[,6],"\n")
		lik<-sum(lik) ##
		}


#_________________________________________________________________

#_________________________________________________________________
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
			cond.lik[i]<-cond.lik[i]*length(which(edge1[,6]==k.rates[i])) ## b=normalizing constant so that values sum to 1
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
			cond.lik[i]<-cond.lik[i]*(alpha/k) ## b=normalizing constant so that values sum to 1
			#cat(cond.lik[i],n.rates[i],alpha/k,"\n")
			}
#index<-which(cond.lik==max(cond.lik))	
#new.rate<-n.rates[index]
#new.rate<-c(new.rate, max(cond.lik))
new.rate<-rbind(cond.lik,n.rates)
return(new.rate)	
}

unique(na.omit(edge[,6])) ## get the rate classes k from the edge matrix

Switch mode: uses Gibbs sampling to change the rate vector (Lart and Phil 2004)

stationary move: change a pre-existing rate (Lart and Phil 2004)

alpha: change the concentration parameter of the DPP (Lart and Phil 2004, Heulsenbeck 2011)
