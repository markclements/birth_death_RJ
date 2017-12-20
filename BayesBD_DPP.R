DPPMCMC.BD<-function(phy,div,div.rate=0.5,ext.rate=0.00,alpha=3,kay=5,ngen=1000000,sfreq=100,alpha.lam=0.01,expl=0.01){
	require(ape)
	
	build.edge<-function(phy,div){
			if (is.null(attr(phy, "order")) || attr(phy, "order") == "cladewise")
			phy<- reorder.phylo(phy,"c") 
			if (length(div)!=length(phy$tip.label)) stop("diversity does not equal tips")
			stem.ages<-branching.times(phy)
			edge<-phy$edge
			edge<-cbind(edge,stem.ages[phy$edge[,1]-Ntip(phy)],phy$edge.length,div[phy$edge[,2]],NA,NA,1,NA,NA,1) 
			# 1=anc node, 2=des node, 3=birth time xi, 4=branch length or waiting time, ti, 5=tip d1versity 6=div rate 7=div node split,8=div rate group, 9 rext rate,  
			dimnames(edge)<-NULL
			return(edge)
			}
	edge<-build.edge(phy,div)
	
	rate.vect<-div.rate
	#_______________________________________________________________________________
	# Main MCMC function
		MCMC<-function(logLik,edge,rate.vect,alpha,kay,alpha.lam,expl){
		
		z<-sample(c(1,2,3),1,prob=c(25,75,25)) #1=Gibbs Sampling, 2=change existing div.rate, 3=change alpha
		
		if (z==1){# Switch mode: uses Gibbs sampling to change the rate vector (Lart and Phil 2004)
			
			
			edge[index,6]<-NA ## remove the rate from rate vector			
			
			x<-c(rexp(kay),unique(na.omit(edge[,6]))) ##make new rate vector kay with new auxillary rates
			
			
			ind.lik<-numeric()
			
			
			for (i in 1:length(x)){ # calcualte ind probs
				edge[index,6]<-x[i]
					if (i<=kay){
						ind.lik[i]<-ind.likelihood(edge[index,])*(alpha/kay)					}
				
					else {
						ind.lik[i]<-ind.likelihood(edge[index,])*length(x)

					}
				}
			p<-(1/sum(ind.lik))*ind.lik #calcualte acceptance prob using normalizing constant
			edge[index,6]<-x[which(p==max(p))] # choose new rate based on ind.prob
			logLik<-likelihood(edge) # calcualte new logLikeliood
			accept=1
			cat(accept," ",logLik," ",x, "\n")
		  return(list(edge,logLik,alpha,accept,1))
			}
		
		else if (z==2){## stationary move: change a pre-existing rate (Lart and Phil 2004)
			
			lineage<-sample(1:length(edge[,1]),1)## select a lineage at random
			index<-which(edge[,6]==edge[lineage,6]) ## find all rates equal to rate of random lineage
			r<-exp(runif(1,0,1)-0.5)*expl
			newedge<-edge
			newedge[index,6]<-edge[index,6]*r
			newlogLik<-likelihood(newedge)
			lr <- newlogLik-logLik
			ratio <- exp(lr)*r
			accept=2
			cat(accept," ",logLik," ", newlogLik," ",r," ", ratio, "\n")
		  if (runif(1,0,1) < ratio) return(list(newedge,newlogLik,alpha,accept,1))
		  else return(list(edge,logLik,alpha,accept,0))
			}
		
		else if (z==3){#change alpha
			
			r<-(runif(1,0,1)-0.5)*lam
			if (alpha+r>=10) new.alpha<-10-((alpha+r)-10) 
			else if (alpha+r<=0) new.alpha<-abs(alpha+r) 
			else new.alpha<-alpha+r
			k<-length(unique(edge[,6]))
			N<-length(edge[,6])
			ratio<-((new.alpha^k)/(alpha^k))*prod((alpha + (1:N) -1)/(new.alpha + (1:N) -1)) 			
			accept=3
			
			cat(accept," ",r," ",ratio, "\n")
		  if (runif(1,0,1) < ratio) return(list(edge,logLik, new.alpha, accept,1))
		  else return(list(edge,logLik, alpha, accept,0))
			}
		
	}

#______________________________________________________
# Main LOG likelihood code 
	likelihood<-function(edge){
		lik<-numeric(length(edge[,1]))
		for (i in 1:length(edge[,1])){
			if (!is.na(edge[i,5])){
				beta <-(expm1(edge[i,6]*edge[i,4]))/(exp(edge[i,6]*edge[i,4])-edge[i,7]) #taxonomic likelihood
				lik[i]<-log(1-beta) + ((edge[i,5]-1)* log (beta)) 
			}
			else {
				lik[i]<-log(edge[i,6])-(edge[i,6]*edge[i,4])-log(1-edge[i,7]*exp(-edge[i,6]*edge[i,3])) 
				}
			#cat(lik[i],"\n")
		}
		#cat(lik,"\n", edge[,6],"\n")
	 lik<-sum(lik)
	}
#_______________________________________________________
#Main alt like code
	ind.likelihood<-function(edge){
		lik<-numeric()
			if (!is.na(edge[5])){
				beta <-(expm1(edge[6]*edge[4]))/(exp(edge[6]*edge[4])-edge[7]) #taxonomic likelihood
				lik<-log(1-beta) + ((edge[5]-1)* log (beta)) 
			}
			else {
				lik<-log(edge[6])-(edge[6]*edge[4])-log(1-edge[7]*exp(-edge[6]*edge[3])) 
				}
	 lik<-exp(lik)
	}
#_______________________________________________________

# initiallize chain
	logLik<-likelihood(edge)

# BEGIN CALCULATION
	count.it<-ngen/sfreq
	
	count.i<-1
	for(i in (1:ngen)) {
		step<-MCMC(logLik,edge,alpha,kay,lam,expl)
	  	edge<-step[[1]]
	  	logLik<-step[[2]]	
	  	alpha<-step[[3]]
	if (i==1 || i%%sfreq==0){
	
	#cat(e," ",logLik,"\n",edge[,6],"\n",edge[,8],"\n","\n")
	plot(phy,cex=0.5)
	edgelabels(round(edge[,6],3),frame="n",cex=0.6)
	cat(count.i,logLik,edge[,6],"\n", sep=" ", file="divrates", append=TRUE)
	#cat(edge[,7],"\n", file="divshift", append=TRUE)
	cat(alpha," ",step[[4]]," ", step[[5]], "\n", sep=" ", file="accept", append=TRUE)
	count.i<-count.i+1
		}
	}
}
