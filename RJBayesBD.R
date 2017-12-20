MCMC.BD<-function(phy,diversity,root.rate=0.5,a=0,e=3,iter=1000,sample=1000){

	require(ape)
	if (is.null(attr(phy, "order")) || attr(phy, "order") == "cladewise")
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
	
	# Main MCMC function
		MCMC<-function(logLik,edge,root.rate,a,e,rate.vect){
		N <-length(which(edge[,7]==1))
		z<-sample(1:4,1,prob=c(5,5,60,30))
		#z<-sample(1:2,1)
		#cat(root.rate)
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
		#	cat(accept," ",log(logLik)," ", log(newlogLik)," ",N,"\n")
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
		#	cat(accept," ",log(logLik)," ", log(newlogLik)," ",N,"\n")
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
		#	cat(accept," ",log(logLik)," ", log(newlogLik)," ",N,"\n")
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
		#	cat(accept," ",log(logLik)," ", log(newlogLik)," ",N,"\n")
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
		#	cat(accept," ",log(logLik)," ", log(newlogLik)," ",N,"\n")
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
		#	cat(accept," ",log(logLik)," ", log(newlogLik)," ",N,"\n")
		  if(runif(1,0,1) < ratio) return(list(newedge,newlogLik,rate.vect,new.root.rate))
		  else return(list(edge,logLik,rate.vect,root.rate))
			}
		}

# Main likelihood code
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


# RJ and tree functions
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


# initiallize and run chain
	edge<-rate.calc(edge=edge,root.rate=root.rate,rate.vect=rate.vect)
	logLik<-likelihood(edge)

# BEGIN CALCULATION
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
# cat(list(i*j,logLik,edge,rates))
	out.edge[[i]]<-cbind(edge[,6],edge[,8])
	out.logLik[[i]]<-logLik
	out.rates[[i]]<-c(rate.vect,i*j)
	cat(log(logLik),"\n",cbind(edge[,6],edge[,8]),"\n")
	}
#return(list(edge,logLik,edge))
list(iter=out.rates,logLik=out.logLik,edge=out.edge)
}


_______________________________________________________________________________
MCMC<-function(logLik,edge,root.rate,a,e,rate.vect){
N <-length(which(edge[,7]==1))
z<-sample(1:4,1,prob=c(5,5,60,30))
#z<-sample(1:2,1)

#cat(root.rate)

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
#	cat(accept," ",log(logLik)," ", log(newlogLik)," ",N,"\n")
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
#	cat(accept," ",log(logLik)," ", log(newlogLik)," ",N,"\n")
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
#	cat(accept," ",log(logLik)," ", log(newlogLik)," ",N,"\n")
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
#	cat(accept," ",log(logLik)," ", log(newlogLik)," ",N,"\n")
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
#	cat(accept," ",log(logLik)," ", log(newlogLik)," ",N,"\n")
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
#	cat(accept," ",log(logLik)," ", log(newlogLik)," ",N,"\n")
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






