MCMC.BD<-function(phy,diversity,root.rate=0.01,a=0,e=3,ngen=1000000,sfreq=100,lam=0.01,expl=0.01){
	require(ape)
	if (is.null(attr(phy, "order")) || attr(phy, "order") == "cladewise")
	phy<- reorder.phylo(phy,"c") 
	
	edge<-phy$edge
	stem.ages<-branching.times(phy)
	edge<-cbind(edge,stem.ages[phy$edge[,1]-Ntip(phy)]) # birth time xi =3
	edge<-cbind(edge,phy$edge.length) # branch length or waiting time, ti = 4
	edge<-cbind(edge,diversity[phy$edge[,2]])# assign diversity to termainals, NAs are for internals. = 5
	edge<-cbind(edge,NA) # div rate column = 6
	edge<-cbind(edge,NA) # div node split column = 7
	edge<-cbind(edge,0) # div rate group column = 8
	
	edge<-cbind(edge,NA) # rext rate column = 9
	edge<-cbind(edge,NA) # rext node split column = 10
	edge<-cbind(edge,0) # rext rate group column = 11

	dimnames(edge)<-NULL
	rate.vect<-numeric(1)
	rext.vect<-numeric(1)
	
#_______________________________________________________________________________
	# Main MCMC function
		MCMC<-function(logLik,edge,root.rate,a,e,rate.vect,lam,expl,rext.vect){
		
		divORext<-sample(c(1,2),1,prob=c(50,50))# 1 == div, 2 == ext
		
	if (divORext==1){	#### SPECIATION SIDE OF PHYLOGENY
		N <-length(which(edge[,7]==1))
		if (N==0) z<-sample(c(1,6),1,prob=c(30,60)) # 1=have to split, 6=change root rate
		else if (N==length(edge[,1])) z<-sample(c(2,5,6),1,prob=c(33,33,33)) #2=must merge, 5=change rate modifier
		else if (N>0 && N<length(edge[,1])) z<-sample(c(3:6),1,prob=c(10,10,30,50))  
		#cat(a,"\n")
		if (z==1){# must split diversity
			x<-sample(1:length(edge[,1]),1)
			rate.vect.new<-rexp(1)
			newedge<-split(edge,x)
			newedge<-rate.calc(newedge,root.rate,rate.vect.new)
			newlogLik<-likelihood(newedge)
			HP<-e/2
			accept=1
			lr <- newlogLik-logLik
			ratio <- exp(lr)*HP
			#cat(accept," ",logLik," ", newlogLik," ",N," ", ratio, "\n")
		  if (runif(1,0,1) < ratio) return(list(newedge,newlogLik,rate.vect.new,root.rate,accept,1,a,rext.vect))
		  else return(list(edge,logLik,rate.vect,root.rate,accept,0,a,rext.vect))
			}
		
		else if (z==2){# must merge diversity
			x<-sample(1:length(edge[,1]),1)
			rate.vect.new<-rate.vect[-x]
			newedge<-merge(edge,x)
			newedge<-rate.calc(newedge,root.rate,rate.vect.new)
			newlogLik<-likelihood(newedge)
			HP<-N/e
			accept=2
			lr <- newlogLik-logLik
			ratio <- exp(lr)*HP
			#cat(accept," ",logLik," ", newlogLik," ",N," ", ratio, "\n")
		  if (runif(1,0,1) < ratio) return(list(newedge,newlogLik,rate.vect.new,root.rate,accept,1,a,rext.vect))
		  else return(list(edge,logLik,rate.vect,root.rate,accept,0,a,rext.vect))
			}
		
		else if (z==3){# regular merge
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
			lr <- newlogLik-logLik
			ratio <- exp(lr)*HP
			#cat(accept," ",logLik," ", newlogLik," ",N," ", ratio, "\n")
		  if (runif(1,0,1) < ratio) return(list(newedge,newlogLik,rate.vect.new,root.rate,accept,1,a,rext.vect))
		  else return(list(edge,logLik,rate.vect,root.rate,accept,0,a,rext.vect))
			}
		
		else if (z==4){# regular split
			x<-which(is.na(edge[,7]))
			if (length(x)>1) x<-sample(x,1)
			rate.vect.new<-c(rate.vect,rexp(1))
			newedge<-split(edge,x)	
			newedge<-rate.calc(newedge,root.rate,rate.vect.new)
			newlogLik<-likelihood(newedge)
			HP<-e/(N+1)
			accept=4
			lr <- newlogLik-logLik
			ratio <- exp(lr)*HP
			#cat(accept," ",logLik," ", newlogLik," ",N," ", ratio, "\n")
		  if(runif(1,0,1) < ratio) return(list(newedge,newlogLik,rate.vect.new,root.rate,accept,1,a,rext.vect))
		  else return(list(edge,logLik,rate.vect,root.rate,accept,0,a,rext.vect))
			}
		
		else if (z==5){ #change exisiting rate multiplier
			#rm<-rexp(1)
			if (N==1) x<-which(edge[,7]==1)
			else x<-sample(which(edge[,7]==1),1)
			s<-edge[x,8]
			rate.vect.new<-rate.vect
			r<-exp(runif(1,0,1)-0.5)*expl
			rate.vect.new[s]<-rate.vect.new[s]*r
			newedge<-rate.calc(edge,root.rate,rate.vect.new)
			newlogLik<-likelihood(newedge)
			accept=5
			lr <- newlogLik-logLik
			ratio <- exp(lr)*r
			#cat(accept," ",logLik," ", newlogLik," ",N," ", ratio, "\n")
		  if(runif(1,0,1) < ratio) return(list(newedge,newlogLik,rate.vect.new,root.rate,accept,1,a,rext.vect))
		  else return(list(edge,logLik,rate.vect,root.rate,accept,0,a,rext.vect))
			}
		
		else if (z==6){ #change root rate
			#rr<-root.rate + runif(1,-0.05,0.05) # simply add a small value to rate
			rr<-runif(1,root.rate-(lam/2),root.rate+(lam/2)) # 
			if (rr<=0) new.root.rate<-0+(0-rr)	
			else if (rr>=10) new.root.rate <-10-(10-rr)
			else new.root.rate<-rr
			#r<-exp(runif(1,0,1)-0.5)*lam #proportionally grow/shrink
			#new.root.rate<-root.rate*r
			newedge<-rate.calc(edge,new.root.rate,rate.vect)
			newlogLik<-likelihood(newedge)
			accept=6
			lr <- newlogLik-logLik
			ratio <- exp(lr)
			#ratio <- exp(lr) * r
			#cat(accept," ",logLik," ", newlogLik," ",N," ", ratio, "\n")
		  if(runif(1,0,1) < ratio) return(list(newedge,newlogLik,rate.vect,new.root.rate,accept,1,a,rext.vect))
		  else return(list(edge,logLik,rate.vect,root.rate,accept,0,a,rext.vect))
			}
		
		}
		
	else if (divORext==2){ ##### EXTENCTION SIDE OF PHYLOGENY
		
		N <-length(which(edge[,10]==1))
		if (N==0) z<-sample(c(1,6),1,prob=c(30,60)) # 1=have to split, 6=change root ext
		else if (N==length(edge[,1])) z<-sample(c(2,5,6),1,prob=c(33,33,33)) #2=must merge, 5=change rate modifier
		else if (N>0 && N<length(edge[,1])) z<-sample(c(3:6),1,prob=c(10,10,30,50))  
		#cat(a,"\n")
		if (z==1){# must split ext
			x<-sample(1:length(edge[,1]),1)
			rext.vect.new<-runif(1)
			newedge<-split.e(edge,x)
			newedge<-rate.calc.e(newedge,a,rext.vect.new)##
			newlogLik<-likelihood(newedge)
			HP<-e/2
			accept=7
			lr <- newlogLik-logLik
			ratio <- exp(lr)*HP
			#cat(accept," ",logLik," ", newlogLik," ",N," ", ratio, "\n")
		  if (runif(1,0,1) < ratio) return(list(newedge,newlogLik,rate.vect,root.rate,accept,1,a,rext.vect.new))
		  else return(list(edge,logLik,rate.vect,root.rate,accept,0,a,rext.vect))
			}
		
		else if (z==2){# must merge diversity
			x<-sample(1:length(edge[,1]),1)
			rext.vect.new<-rext.vect[-x]
			newedge<-merge.e(edge,x)
			newedge<-rate.calc.e(newedge,a,rext.vect.new)
			newlogLik<-likelihood(newedge)
			HP<-N/e
			accept=8
			lr <- newlogLik-logLik
			ratio <- exp(lr)*HP
			#cat(accept," ",logLik," ", newlogLik," ",N," ", ratio, "\n")
		  if (runif(1,0,1) < ratio) return(list(newedge,newlogLik,rate.vect,root.rate,accept,1,a,rext.vect.new))
		  else return(list(edge,logLik,rate.vect,root.rate,accept,0,a,rext.vect))
			}
		
		else if (z==3){# regular merge
			if (N==1){ 
				x<-which(edge[,10]==1)
				HP<-2/e
				}
			else { x<-sample(which(edge[,10]==1),1)
				HP<-N/e
				}
			s<-edge[x,11]
			rext.vect.new<-rext.vect[-s]
			newedge<-merge.e(edge,x)
			newedge<-rate.calc.e(newedge,a,rext.vect.new)
			newlogLik<-likelihood(newedge)
			accept=9
			lr <- newlogLik-logLik
			ratio <- exp(lr)*HP
			#cat(accept," ",logLik," ", newlogLik," ",N," ", ratio, "\n")
		  if (runif(1,0,1) < ratio) return(list(newedge,newlogLik,rate.vect,root.rate,accept,1,a,rext.vect.new))
		  else return(list(edge,logLik,rate.vect,root.rate,accept,0,a,rext.vect))
			}
		
		else if (z==4){# regular split
			x<-which(is.na(edge[,10]))
			if (length(x)>1) x<-sample(x,1)
			rext.vect.new<-c(rext.vect,runif(1))
			newedge<-split.e(edge,x)	
			newedge<-rate.calc.e(newedge,a,rext.vect.new)
			newlogLik<-likelihood(newedge)
			HP<-e/(N+1)
			accept=10
			lr <- newlogLik-logLik
			ratio <- exp(lr)*HP
			#cat(accept," ",logLik," ", newlogLik," ",N," ", ratio, "\n")
		  if(runif(1,0,1) < ratio) return(list(newedge,newlogLik,rate.vect,root.rate,accept,1,a,rext.vect.new))
		  else return(list(edge,logLik,rate.vect,root.rate,accept,0,a,rext.vect))
			}
		
		else if (z==5){ #change exisiting rate multiplier
			if (N==1) x<-which(edge[,10]==1)
			else x<-sample(which(edge[,10]==1),1)
			s<-edge[x,11]
			rext.vect.new<-rext.vect
			r<-runif(1)
			rext.vect.new[s]<-r
			newedge<-rate.calc.e(edge,a,rext.vect.new)
			newlogLik<-likelihood(newedge)
			accept=11
			lr <- newlogLik-logLik
			ratio <- exp(lr)
			#cat(accept," ",logLik," ", newlogLik," ",N," ", ratio, "\n")
		  if(runif(1,0,1) < ratio) return(list(newedge,newlogLik,rate.vect,root.rate,accept,1,a,rext.vect.new))
		  else return(list(edge,logLik,rate.vect,root.rate,accept,0,a,rext.vect))
			}
		
		else if (z==6){ #change root rate
			a.new<-runif(1)
			newedge<-rate.calc.e(edge,a.new,rext.vect)
			newlogLik<-likelihood(newedge)
			accept=12
			lr <- newlogLik-logLik
			ratio <- exp(lr)
			#cat(accept," ",logLik," ", newlogLik," ",N," ", ratio, "\n")
		  if(runif(1,0,1) < ratio) return(list(newedge,newlogLik,rate.vect,root.rate,accept,1,a.new,rext.vect))
		  else return(list(edge,logLik,rate.vect,root.rate,accept,0,a,rext.vect))
			}
		}	

}

#______________________________________________________
# Main LOG likelihood code 
		likelihood<-function(edge){
		lik<-numeric(length(edge[,1]))
		for (i in 1:length(edge[,1])){
				
			if (!is.na(edge[i,5])){
				beta <-(expm1(edge[i,6]*edge[i,4]))/(exp(edge[i,6]*edge[i,4])-(edge[i,9])) #taxonomic likelihood
				lik[i]<-log(1-beta) + ((edge[i,5]-1)* log (beta)) 
			}
			else {
				lik[i]<-log(edge[i,6])-(edge[i,6]*edge[i,4])-log(1-(edge[i,9])*exp(-edge[i,6]*edge[i,3])) 
				#lik[i]<-log(edge[i,6]*((1-a*exp(-edge[i,6]*edge[i,3]))^(-1))*(exp(-edge[i,6]*edge[i,4]))) 
				}
			#cat(lik[i],"\n")
		}
		#cat(lik,"\n", edge[,6],"\n")
		lik<-sum(lik)
		}


#_________________________________________________________________
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
		
		split.e<-function(edge,x){
			if (!is.na(edge[x,10])) {
				stop("that ext lineage is already split, choose another")
				}
			edge[x,10]<-1
			tmp<-find.clade(edge,x)
			group<-max(unique(edge[,11]))+1
			oldgroup<-edge[x,11]
			if (length(tmp)==1) edge[x,11]<-group
			else for (i in tmp[1]:tmp[2]){
				if (edge[i,11]!=oldgroup) next
				edge[i,11]<-group
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
		
		
		merge.e<- function(edge,x){
			if (is.na(edge[x,10])) {
				stop("that ext lineage is not split, choose another")
				}
			edge[x,10]<-NA
			tmp<-find.clade(edge,x)
			oldgroup<-edge[x,11]
			d<-which(edge[,2]==edge[x,1])
			if (length(d)==0) anc.group<-0
			else anc.group<-edge[d,11]
			if (length(tmp)==1) edge[x,11]<-anc.group
			else for (i in tmp[1]:tmp[2]){
				if (edge[i,11]!=oldgroup) next
				edge[i,11]<-anc.group
				}
			edge[edge[,11]>oldgroup,11]<-edge[edge[,11]>oldgroup,11]-1
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


		rate.calc.e<-function(edge,a,rext.vect){
			for (i in 1:length(edge[,1])){
				d<-which(edge[,2]==edge[i,1])
				if (is.na(edge[i,10]) && length(d)==0) edge[i,9]<-a # non-split nodes that have background rate
				else if (!is.na(edge[i,10]) && length(d)==0) edge[i,9]<-add.rate(a,rext.vect[edge[i,11]]) # 
				else if (is.na(edge[i,10]) && length(d)==1) edge[i,9]<-edge[d,9]
				else if (!is.na(edge[i,10]) && length(d)==1) edge[i,9]<-add.rate(edge[d,9],rext.vect[edge[i,11]])
				}
			edge
			}
			
		add.rate<-function(ex,value,a,b){
			exprime<-ex+value	
			if (exprime<=a)	zz<-(2*a)-exprime
			else if (exprime>=b) zz<-(2*b)-exprime
			else zz<-exprime
			return (zz)	
			}
			

#_______________________________________________________

# initiallize chain
	edge<-rate.calc(edge,root.rate,rate.vect)
	edge<-rate.calc.e(edge,a,rext.vect)
	logLik<-likelihood(edge)

# BEGIN CALCULATION
	count.it<-ngen/sfreq
	#out.edge <- vector("list",count.it)
	#out.logLik <- 1:count.it
	#out.rates <- vector("list",count.it)
	#out.ratecat <- vector("list",count.it)
	#out.rateshift <- vector("list",count.it)
	#out.proposal<-1:count.it
	#out.accept<-1:count.it
	#out.root.rate<-1:count.it
	#out.a<-1:count.it
	#out.rext.rates <- vector("list",count.it)
	#out.rext.cat <- vector("list",count.it)
	#out.rext.shift <- vector("list",count.it)
	#out.edge.rext<- vector("list",count.it)
	count.i<-1
	for(i in (1:ngen)) {
		step<-MCMC(logLik,edge,root.rate,a,e,rate.vect,lam,expl,rext.vect)
	  	edge<-step[[1]]
	  	logLik<-step[[2]]	
	  	rate.vect<-step[[3]]
	  	root.rate<-step[[4]]
	  	a<-step[[7]] 
	  	rext.vect<-step[[8]]
		#e<-step[[7]]
	if (i==1 || i%%sfreq==0){
	#out.edge[[count.i]]<-edge[,6]
	#out.logLik[[count.i]]<-logLik
	#out.rates[[count.i]]<-rate.vect
	#out.ratecat[[count.i]]<-edge[,8]
	#out.rateshift[[count.i]] <-edge[,7]
	#out.proposal[[count.i]]<-step[[5]]
	#out.accept[[count.i]]<-step[[6]]
	#out.root.rate[[count.i]]<-root.rate
	#out.a[[count.i]]<-a
	#out.rext.rates[[count.i]] <- rext.vect
	#out.rext.cat[[count.i]] <- edge[,11]
	#out.rext.shift[[count.i]] <- edge[,10]	
	#out.edge.rext[[count.i]]<-edge[,9]
	plot(phy,cex=0.5)
	edgelabels(round(edge[,6],3),frame="n",cex=0.6)
	#cat(logLik,"\n",root.rate, rate.vect,"\n",edge[,8],"\n",a, rext.vect,"\n",edge[,11],"\n")
	cat(count.i,logLik,edge[,6],"\n", sep=" ", file="divrates", append=TRUE)
	cat(edge[,9],"\n", file="rextrates", append=TRUE)
	cat(edge[,7],"\n", file="divshift", append=TRUE)
	cat(edge[,10],"\n", file="rextshift", append=TRUE)
	cat(step[[5]], step[[6]],"\n", sep=" ", file="accept", append=TRUE)
	count.i<-count.i+1
		}
	}
#list(rate.vect=out.rates,logLik=out.logLik,edge.rates=out.edge,rate.cat=out.ratecat,rate.shift=out.rateshift,proposal=out.proposal, accept=out.accept, root.rate=out.root.rate,a=out.a, rext.rate=out.rext.rates, rext.cat=out.rext.cat, rext.shift=out.rext.shift, edge.rext<-out.edge.rext )
}
