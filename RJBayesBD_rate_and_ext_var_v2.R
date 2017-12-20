MCMC.BD<-function(phy,div,directory="~/Desktop/trials",root.rate=0.01,root.rate.e=0,e=3,ngen=1000000,sfreq=100,lam.root=0.01,lam.div=0.1,lamroot.e=0.1, lam.e=0.1,lamsplit.e=0.5){
	require(ape)
	edge<-build.edge(phy,div)
	rate.vect<-numeric(1)
	rext.vect<-numeric(1)
	setwd(directory)
	dir.create(paste("RJ_BD_run",date()))
	z<-file.info(dir())
	zq<-which(z[,5]==max(z[,5]))
	setwd(rownames(z[zq,]))
#_______________________________________________________________________________
	# Main MCMC function
	MCMC<-function(logLik,edge,root.rate,root.rate.e,e,rate.vect,lam.root,lam.div,rext.vect,lamroot.e,lam.e,lamsplit.e){
		
		divORext<-sample(c(1,2),1,prob=c(1,0))# 1 == div, 2 == ext
		
	if (divORext==1){	#### SPECIATION SIDE OF PHYLOGENY
		N <-length(which(edge[,7]==1))
		if (N==0) z<-sample(c(1,6),1,prob=c(.3,.7)) # 1=have to split, 6=change root rate
		else if (N==length(edge[,1])) z<-sample(c(2,5,6),1,prob=c(.33,.33,.33)) #2=must merge, 5=change rate modifier
		else if (N>0 && N<length(edge[,1])) z<-sample(c(3:6),1,prob=c(.10,.10,.30,.50))  
		#cat(root.rate.e,"\n")
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
			cat(accept," ",logLik," ", newlogLik," ",N," ", ratio, "\n")
		  if (runif(1,0,1) < ratio) return(list(newedge,newlogLik,rate.vect.new,root.rate,accept,1,root.rate.e,rext.vect))
		  else return(list(edge,logLik,rate.vect,root.rate,accept,0,root.rate.e,rext.vect))
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
			cat(accept," ",logLik," ", newlogLik," ",N," ", ratio, "\n")
		  if (runif(1,0,1) < ratio) return(list(newedge,newlogLik,rate.vect.new,root.rate,accept,1,root.rate.e,rext.vect))
		  else return(list(edge,logLik,rate.vect,root.rate,accept,0,root.rate.e,rext.vect))
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
			cat(accept," ",logLik," ", newlogLik," ",N," ", ratio, "\n")
		  if (runif(1,0,1) < ratio) return(list(newedge,newlogLik,rate.vect.new,root.rate,accept,1,root.rate.e,rext.vect))
		  else return(list(edge,logLik,rate.vect,root.rate,accept,0,root.rate.e,rext.vect))
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
			cat(accept," ",logLik," ", newlogLik," ",N," ", ratio, "\n")
		  if(runif(1,0,1) < ratio) return(list(newedge,newlogLik,rate.vect.new,root.rate,accept,1,root.rate.e,rext.vect))
		  else return(list(edge,logLik,rate.vect,root.rate,accept,0,root.rate.e,rext.vect))
			}
		
		else if (z==5){ #change exisiting rate multiplier
			if (N==1) x<-which(edge[,7]==1)
			else x<-sample(which(edge[,7]==1),1)
			s<-edge[x,8]
			rate.vect.new<-rate.vect
			r<-exp(runif(1,0,1)-0.5)*lam.div
			rate.vect.new[s]<-rate.vect.new[s]*r
			newedge<-rate.calc(edge,root.rate,rate.vect.new)
			newlogLik<-likelihood(newedge)
			accept=5
			lr <- newlogLik-logLik
			ratio <- exp(lr)*r
			cat(accept," ",logLik," ", newlogLik," ",N," ", ratio, "\n")
		  if(runif(1,0,1) < ratio) return(list(newedge,newlogLik,rate.vect.new,root.rate,accept,1,root.rate.e,rext.vect))
		  else return(list(edge,logLik,rate.vect,root.rate,accept,0,root.rate.e,rext.vect))
			}
		
		else if (z==6){ #change root rate
			new.root.rate<-add.rate(root.rate,lam.root*(runif(1)-0.5),0,10)		
			newedge<-rate.calc(edge,new.root.rate,rate.vect)
			newlogLik<-likelihood(newedge)
			accept=6
			lr <- newlogLik-logLik
			ratio <- exp(lr)
			#ratio <- exp(lr) * r
		  if(runif(1,0,1) < ratio) return(list(newedge,newlogLik,rate.vect,new.root.rate,accept,1,root.rate.e,rext.vect))
		  else return(list(edge,logLik,rate.vect,root.rate,accept,0,root.rate.e,rext.vect))
		  cat(accept," ",logLik," ", newlogLik," ",root.rate," ", new.root.rate, "\n")			
			}
		
		}
		
	else if (divORext==2){ ##### EXTENCTION SIDE OF PHYLOGENY
		
		N <-length(which(edge[,10]==1))
		if (N==0) z<-sample(c(1,6),1,prob=c(30,70)) # 1=have to split, 6=change root ext
		else if (N==length(edge[,1])) z<-sample(c(2,5,6),1,prob=c(33,33,33)) #2=must merge, 5=change rate modifier
		else if (N>0 && N<length(edge[,1])) z<-sample(c(3:6),1,prob=c(10,10,30,50))  
		#cat(root.rate.e,"\n")
		if (z==1){# must split ext
			x<-sample(1:length(edge[,1]),1)
			rext.vect.new<-lamsplit.e*runif(1)-0.5 
			newedge<-split.e(edge,x)
			newedge<-rate.calc.e(newedge,root.rate.e,rext.vect.new)##
			newlogLik<-likelihood(newedge)
			HP<-e/2
			accept=7
			lr <- newlogLik-logLik
			ratio <- exp(lr)*HP
			#cat(accept," ",logLik," ",ratio," ",root.rate.e," ",rext.vect.new,"\n")
		  if (runif(1,0,1) < ratio) return(list(newedge,newlogLik,rate.vect,root.rate,accept,1,root.rate.e,rext.vect.new))
		  else return(list(edge,logLik,rate.vect,root.rate,accept,0,root.rate.e,rext.vect))
			}
		
		else if (z==2){# must merge diversity
			x<-sample(1:length(edge[,1]),1)
			rext.vect.new<-rext.vect[-x]
			newedge<-merge.e(edge,x)
			newedge<-rate.calc.e(newedge,root.rate.e,rext.vect.new)
			newlogLik<-likelihood(newedge)
			HP<-N/e
			accept=8
			lr <- newlogLik-logLik
			ratio <- exp(lr)*HP
			#cat(accept," ",logLik," ", newlogLik," ",N," ", ratio, "\n")
		  if (runif(1,0,1) < ratio) return(list(newedge,newlogLik,rate.vect,root.rate,accept,1,root.rate.e,rext.vect.new))
		  else return(list(edge,logLik,rate.vect,root.rate,accept,0,root.rate.e,rext.vect))
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
			newedge<-rate.calc.e(newedge,root.rate.e,rext.vect.new)
			newlogLik<-likelihood(newedge)
			accept=9
			lr <- newlogLik-logLik
			ratio <- exp(lr)*HP
			#cat(accept," ",logLik," ", newlogLik," ",N," ", ratio, "\n")
		  if (runif(1,0,1) < ratio) return(list(newedge,newlogLik,rate.vect,root.rate,accept,1,root.rate.e,rext.vect.new))
		  else return(list(edge,logLik,rate.vect,root.rate,accept,0,root.rate.e,rext.vect))
			}
		
		else if (z==4){# regular split
			x<-which(is.na(edge[,10]))
			if (length(x)>1) x<-sample(x,1)
			rext.vect.new<-c(rext.vect,lamsplit.e*runif(1)-0.5) ### perhaps add a tuning parameter here 
			newedge<-split.e(edge,x)	
			newedge<-rate.calc.e(newedge,root.rate.e,rext.vect.new)
			newlogLik<-likelihood(newedge)
			HP<-e/(N+1)
			accept=10
			lr <- newlogLik-logLik
			ratio <- exp(lr)*HP
			#cat(accept," ",logLik," ",ratio," ",root.rate.e," ",rext.vect.new,"\n")
		  if(runif(1,0,1) < ratio) return(list(newedge,newlogLik,rate.vect,root.rate,accept,1,root.rate.e,rext.vect.new))
		  else return(list(edge,logLik,rate.vect,root.rate,accept,0,root.rate.e,rext.vect))
			}
		
		else if (z==5){ #change exisiting rate multiplier
			if (N==1) x<-which(edge[,10]==1)
			else x<-sample(which(edge[,10]==1),1)
			s<-edge[x,11]
			rext.vect.new<-rext.vect
			rext.vect.new[s]<-add.rate(rext.vect.new[s],lam.e*(runif(1)-0.5),0,1) 
			newedge<-rate.calc.e(edge,root.rate.e,rext.vect.new)
			newlogLik<-likelihood(newedge)
			accept=11
			lr <- newlogLik-logLik
			ratio <- exp(lr)
			#cat(accept," ",logLik," ",ratio," ",rext.vect.new," ", "\n")
		  if(runif(1,0,1) < ratio) return(list(newedge,newlogLik,rate.vect,root.rate,accept,1,root.rate.e,rext.vect.new))
		  else return(list(edge,logLik,rate.vect,root.rate,accept,0,root.rate.e,rext.vect))
			}
		
		else if (z==6){ #change root rate
			newroot.rate.e<-add.rate(root.rate.e,lamroot.e*(runif(1)-0.5),0,1)
			newedge<-rate.calc.e(edge,newroot.rate.e,rext.vect)
			newlogLik<-likelihood(newedge)
			accept=12
			lr <- newlogLik-logLik
			ratio <- exp(lr)
			#cat(accept," ",logLik," ", newlogLik," ",root.rate.e," ",newroot.rate.e," ", ratio, "\n")
		  if(runif(1,0,1) < ratio) return(list(newedge,newlogLik,rate.vect,root.rate,accept,1,newroot.rate.e,rext.vect))
		  else return(list(edge,logLik,rate.vect,root.rate,accept,0,root.rate.e,rext.vect))
			}
		}	
  
}

#______________________________________________________
# Main LOG likelihood code 
		likelihood<-function(edge, internal.only=FALSE){
		lik.e<-numeric()
		lik.i<-numeric()
		N<-0
		for (i in 1:length(edge[,1])){
				
			if (!is.na(edge[i,5])){
				if (internal.only) lik.e<-0
				else {
					beta <-(expm1(edge[i,6]*edge[i,4]))/(exp(edge[i,6]*edge[i,4])-(edge[i,9])) #taxonomic likelihood
					lik.e<-c(lik.e,log(1-beta) + ((edge[i,5]-1)* log (beta)))
					#lik[i]<-(1-beta)*(beta^(edge[i,5]-1))
				}
			}
			else {
				lik.i<-c(lik.i,log(edge[i,6])-(edge[i,6]*edge[i,4])-log(1-(edge[i,9])*exp(-edge[i,6]*edge[i,3])))
				#lik[i]<-(((edge[i,"d_rate"])*(1-edge[i,"e_rate"]*exp(-edge[i,"d_rate"]*edge[i,"bt"])))^(-1))*exp(-edge[i,"d_rate"]*edge[i,"brl"])
				N<-N+1
				}
			}
		
		lik<-sum(lik.e)+(sum(lik.i))
		#lik<-sum(log(lik))
		#lik<-prod(lik)
		return(lik)
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


		rate.calc.e<-function(edge,root.rate.e,rext.vect){
			for (i in 1:length(edge[,1])){
				d<-which(edge[,2]==edge[i,1])
				if (is.na(edge[i,10]) && length(d)==0) edge[i,9]<-root.rate.e # non-split nodes that have background rate
				else if (!is.na(edge[i,10]) && length(d)==0) edge[i,9]<-add.rate(root.rate.e,rext.vect[edge[i,11]],0,1) # 
				else if (is.na(edge[i,10]) && length(d)==1) edge[i,9]<-edge[d,9]
				else if (!is.na(edge[i,10]) && length(d)==1) edge[i,9]<-add.rate(edge[d,9],rext.vect[edge[i,11]],0,1)
				}
			edge
			}
			
		add.rate<-function(ex,value,lower,upper){
			exprime<-ex+value	
			if (exprime<=lower)	ans<-(2*lower)-exprime
			else if (exprime>=upper) ans<-(2*upper)-exprime
			else ans<-exprime
			return (ans)	
			}
		
		
		build.edge<-function(phy,div){
			phy<- reorder.phylo(phy,"c") 
			if (length(div)!=length(phy$tip.label)) stop("diversity does not equal tips")
			stem.ages<-branching.times(phy)
			edge<-phy$edge
			edge<-cbind(edge,stem.ages[phy$edge[,1]-Ntip(phy)],phy$edge.length,div[phy$edge[,2]],NA,NA,1,NA,NA,1) 
			colnames(edge)<-c("an","dn","bt","brl","tip_d","d_rate","d_nodesplit","d_rategroup","e_rate","e_nodesplit","e_rategroup")  
			rownames(edge)<-NULL
			return(edge)
			}

#_______________________________________________________

# initiallize chain
	edge<-rate.calc(edge,root.rate,rate.vect)
	edge<-rate.calc.e(edge,root.rate.e,rext.vect)
	logLik<-likelihood(edge)
	
# BEGIN CALCULATION
	count.it<-ngen/sfreq
	count.i<-1
	for(i in (1:ngen)) {
		step<-MCMC(logLik,edge,root.rate,root.rate.e,e,rate.vect,lam.root,lam.div,rext.vect,lamroot.e,lam.e,lamsplit.e)
	  	edge<-step[[1]]
	  	logLik<-step[[2]]	
	  	rate.vect<-step[[3]]
	  	root.rate<-step[[4]]
	  	root.rate.e<-step[[7]] 
	  	rext.vect<-step[[8]]
	  	#cat(logLik,"\n")
	if (i==1 || i%%sfreq==0){
	#plot(phy,cex=0.5,no.margin=TRUE)
	#edgelabels(text=round(edge[,6],3),adj=c(1,0),frame="n",col="black",cex=0.5)
	#edgelabels(text=round(edge[,9],3),adj=c(1,1),frame="n",col="red",cex=0.5)
	#edgelabels(round(edge[,9],3),frame="n",cex=0.6)
	#cat(logLik),"\n", edge,"\n")
	#cat(count.i,logLik,step[[5]], step[[6]], "\n", sep=" ", file="trace", append=TRUE)
	#cat(edge[,6],"\n", file="divrates", append=TRUE)
	#cat(edge[,9],"\n", file="rextrate", append=TRUE)
	#cat(edge[,7],"\n", file="divshift", append=TRUE)
	#cat(edge[,10],"\n",file="rextshift",append=TRUE)
	count.i<-count.i+1
		}
	}
}
