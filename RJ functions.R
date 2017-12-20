#______________________________________________________
#  Main LOG likelihood code updated 23 May 2012
#	Lnlikelihood<-function(edge, internal.only){
#		lik.e<-numeric()
#		lik.i<-numeric()
#		N<-0
#		for (i in 1:length(edge[,1])){
#				
#			if (!is.na(edge[i,5])){
#				if (internal.only) lik.e<-0
#				else {
#					beta <-(expm1(edge[i,6]*edge[i,4]))/(exp(edge[i,6]*edge[i,4])-(edge[i,9])) #taxonomic likelihood
#					lik.e<-c(lik.e,log(1-beta) + ((edge[i,5]-1)* log (beta)))
#					#lik[i]<-(1-beta)*(beta^(edge[i,5]-1))
#				}
#			}
#			else {
#				lik.i<-c(lik.i,log(edge[i,6])-(edge[i,6]*edge[i,4])-log(1-(edge[i,9])*exp(-edge[i,6]*edge[i,3])))
#				#lik[i]<-(((edge[i,"d_rate"])*(1-edge[i,"e_rate"]*exp(-edge[i,"d_rate"]*edge[i,"bt"])))^(-1))*exp(-edge[i,"d_rate"]*edge[i,"brl"])
#				N<-N+1
#				}
#			}
#		
#		lik<-sum(lik.e)+(sum(lik.i))
#		return(lik)
#		}
#_________________________________________________________________


# Vectorized LOG likelihood code updated 29 May 2012
	Lnlikelihood<-function(edge, internal.only){
		if (internal.only) lik.e<-0
			else {
			tips<-edge[which(!is.na(edge[,"tip_d"])),]
			beta <-(expm1(tips[,6]*tips[,4]))/(exp(tips[,6]*tips[,4])-(tips[,9])) #taxonomic likelihood
			lik.e<-log(1-beta) + ((tips[,5]-1)* log(beta))
				}	
		int<-edge[which(is.na(edge[,"tip_d"])),]
		lik.i<-log(int[,6])-(int[,6]*int[,4])-log(1-(int[,9])*exp(-int[,6]*int[,3]))
	
		lik<-sum(lik.e)+(sum(lik.i))
		return(lik)
		}



#_________________________________________________________________
		
		find.clade<-function(edge,x){
			if (!is.na(edge[x,"tip_d"])) clade<-x #tip clade
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
		
##_________________________			
		
		split<-function(edge,x){
			if (!is.na(edge[x,"d_nodesplit"])) {
				stop("that lineage is already split, choose another")
				}
			edge[x,"d_nodesplit"]<-1
			tmp<-find.clade(edge,x)
			group<-max(unique(edge[,"d_rategroup"]))+1
			oldgroup<-edge[x,"d_rategroup"]
			if (length(tmp)==1) edge[x,"d_rategroup"]<-group
			else for (i in tmp[1]:tmp[2]){
				if (edge[i,"d_rategroup"]!=oldgroup) next
				edge[i,"d_rategroup"]<-group
				}
			edge	
			}
##_________________________			
		
		split.e<-function(edge,x){
			if (!is.na(edge[x,"e_nodesplit"])) {
				stop("that ext lineage is already split, choose another")
				}
			edge[x,"e_nodesplit"]<-1
			tmp<-find.clade(edge,x)
			group<-max(unique(edge[,"e_rategroup"]))+1
			oldgroup<-edge[x,"e_rategroup"]
			if (length(tmp)==1) edge[x,"e_rategroup"]<-group
			else for (i in tmp[1]:tmp[2]){
				if (edge[i,"e_rategroup"]!=oldgroup) next
				edge[i,"e_rategroup"]<-group
				}
			edge	
			}
##_________________________			
		
		merge<- function(edge,x){
			if (is.na(edge[x,"d_nodesplit"])) {
				stop("that lineage is NOT split, choose another")
				}
			edge[x,"d_nodesplit"]<-NA
			tmp<-find.clade(edge,x)
			oldgroup<-edge[x,"d_rategroup"]
			d<-which(edge[,2]==edge[x,1])
			if (length(d)==0) anc.group<-1
			else anc.group<-edge[d,"d_rategroup"]
			if (length(tmp)==1) edge[x,"d_rategroup"]<-anc.group
			else for (i in tmp[1]:tmp[2]){
				if (edge[i,"d_rategroup"]!=oldgroup) next
				edge[i,"d_rategroup"]<-anc.group
				}
			edge[edge[,"d_rategroup"]>oldgroup,"d_rategroup"]<-edge[edge[,"d_rategroup"]>oldgroup,"d_rategroup"]-1
			edge	
			}
		
##_________________________			
		
		merge.e<- function(edge,x){
			if (is.na(edge[x,"e_nodesplit"])) {
				stop("that ext lineage is not split, choose another")
				}
			edge[x,"e_nodesplit"]<-NA
			tmp<-find.clade(edge,x)
			oldgroup<-edge[x,"e_rategroup"]
			d<-which(edge[,2]==edge[x,1])
			if (length(d)==0) anc.group<-1
			else anc.group<-edge[d,"e_rategroup"]
			if (length(tmp)==1) edge[x,"e_rategroup"]<-anc.group
			else for (i in tmp[1]:tmp[2]){
				if (edge[i,"e_rategroup"]!=oldgroup) next
				edge[i,"e_rategroup"]<-anc.group
				}
			edge[edge[,"e_rategroup"]>oldgroup,"e_rategroup"]<-edge[edge[,"e_rategroup"]>oldgroup,"e_rategroup"]-1
			edge	
			}
		
		
##_________________________			
		
		rate.calc<-function(edge,rate.vect,assign.d){
			if (assign.d) edge[,"d_rate"]<-rate.vect[edge[,"d_rategroup"]]
			else
			for (i in 1:length(edge[,1])){
				d<-which(edge[,2]==edge[i,1])
				if (is.na(edge[i,7]) && length(d)==0) edge[i,6]<-rate.vect[1]
				else if (!is.na(edge[i,7]) && length(d)==0) edge[i,6]<-rate.vect[1]*rate.vect[edge[i,8]]
				else if (is.na(edge[i,7]) && length(d)==1) edge[i,6]<-edge[d,6]
				else if (!is.na(edge[i,7]) && length(d)==1) edge[i,6]<-edge[d,6]*rate.vect[edge[i,8]]
				}
			edge
			}

##_________________________	
		
		rate.calc.e<-function(edge,rext.vect,assign.e){
			if (assign.e) edge[,"e_rate"]<-rext.vect[edge[,"e_rategroup"]]
			else
			for (i in 1:length(edge[,1])){
				d<-which(edge[,2]==edge[i,1])
				if (is.na(edge[i,10]) && length(d)==0) edge[i,9]<-rext.vect[1] # non-split nodes that have background rate
				else if (!is.na(edge[i,10]) && length(d)==0) edge[i,9]<-add.rate(rext.vect[1],rext.vect[edge[i,11]],0,1) # 
				else if (is.na(edge[i,10]) && length(d)==1) edge[i,9]<-edge[d,9]
				else if (!is.na(edge[i,10]) && length(d)==1) edge[i,9]<-add.rate(edge[d,9],rext.vect[edge[i,11]],0,1)
				}
			edge
			}
				
##_________________________				
		
		add.rate<-function(ex,value,lower,upper){ ####!!!! make sure value to add does not exceed interval
			exprime<-ex+value	
			if (exprime<=lower)	ans<-(2*lower)-exprime
			else if (exprime>=upper) ans<-(2*upper)-exprime
			else ans<-exprime
			return (ans)	
			}
		
##_________________________			
		
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
##_________________________		
		
		splitORmerge.d<-function(edge,rate.vect,e,internal.only,assign.d){
			N.R<-length(which(edge[,"d_nodesplit"]==1)) ##number of rate shifts 
			N.N<-length(edge[,1])						##number of lineages (nodes)
			if (internal.only) x<-sample(which(is.na(edge[,"tip_d"])),1) else x<-sample(1:length(edge[,1]),1)
			#split tree
			if (is.na(edge[x,"d_nodesplit"])){		
				prop="split_d"
				rate.vect.new<-c(rate.vect,rexp(1))	
				newedge<-rate.calc(split(edge,x),rate.vect.new,assign.d)
				H = (N.R+1)/(N.N-N.R)### from Drummond and Suchard 2010
				P = ptpois(N.R+1,e,N.N)/ptpois(N.R,e,N.N)
				HP=H*P
				#	if (N.R==0) 	HP<-e/2 #N=0 i.e. no split 
				#	else 		HP<-e/(N.R+1) #N>0 i.e. at least one split
				}
			#merge tree
			else {
				prop="merge_d"
				s<-edge[x,"d_rategroup"]
				rate.vect.new<-rate.vect[-s]
				newedge<-rate.calc(merge(edge,x),rate.vect.new,assign.d)
				H=(N.N-(N.R+1))/N.N
				P=ptpois(N.R-1,e,N.N)/ptpois(N.R,e,N.N)
				HP=H*P
				#	if (N==1) 	HP<-2/e #N=1 i.e. only one split to merge
				#	else 		HP<-N/e #N>1 i.e. more that one split to merge
				}
				#cat("N.R",N.R,"H=",H,"P=",P,"ratio=",H*P,"\n",sep=" ")
			return(list(newedge=newedge,rate.vect.new=rate.vect.new,HP=HP,prop=prop))
		}	

##_________________________		
		
		splitORmerge.e<-function(edge,rext.vect,e,internal.only,assign.e,lamsplit.e){
			N.R<-length(which(edge[,"e_nodesplit"]==1))
			N.N<-length(edge[,1])
			if (internal.only) x<-sample(which(is.na(edge[,"tip_d"])),1) else x<-sample(1:length(edge[,1]),1)
			#split tree
			if (is.na(edge[x,"e_nodesplit"])){		
				prop="split_e"
			if (assign.e) 	rext.vect.new<-c(rext.vect,runif(1))	
			else			rext.vect.new<-c(rext.vect,lamsplit.e*runif(1)-0.5)
				newedge<-rate.calc.e(split.e(edge,x),rext.vect.new,assign.e)
				H = (N.R+1)/(N.N-N.R)### from Drummond and Suchard 2010
				P = ptpois(N.R+1,e,N.N)/ptpois(N.R,e,N.N)
				HP=H*P
				#	if (N==0) 	HP<-e/2 #N=0 i.e. no split 
				#	else 		HP<-e/(N+1) #N>0 i.e. at least one split
				}
			#merge tree
			else {
				prop="merge_e"
				s<-edge[x,"e_rategroup"]
				rext.vect.new<-rext.vect[-s]
				newedge<-rate.calc.e(merge.e(edge,x),rext.vect.new,assign.e)
				H=(N.N-(N.R+1))/N.N
				P=ptpois(N.R-1,e,N.N)/ptpois(N.R,e,N.N)				
				HP=H*P
				#	if (N==1) 	HP<-2/e #N=1 i.e. only one split to merge
				#	else 		HP<-N/e #N>1 i.e. more that one split to merge
				}
			return(list(newedge=newedge,rext.vect.new=rext.vect.new,HP=HP,prop=prop))
		}

##_________________________		

#general statistical function for finding probability of a value within a truncated Poisson distribution
#author: JM EASTMAN 2010

ptpois <-function(x, lambda, k) {
	p.k=ppois(k, lambda, lower.tail=FALSE)
	p.x=ppois(x, lambda, lower.tail=FALSE)
	ptp=p.x/(1-p.k)
	return(ptp)
}		

#general statistical function for random samples from a truncated Poisson distribution
#author: JM EASTMAN 2010


##_________________________		
rtpois <-
function(N, lambda, k) {
	p=ppois(k, lambda, lower.tail=FALSE)
	out=qpois(runif(N, min=0, max=1-p), lambda)
	out
}


##_________________________		
		
		adjustrate.d<-function(edge,rate.vect,assign.d,rootprob.d,lam.d){
			if (length(rate.vect)==1) s=1 
			else if (length(rate.vect)>1 & runif(1)>(1-rootprob.d)) s=1
			else s=sample(2:length(rate.vect),1)
			
			r<-exp(runif(1,0,1)-0.5)*lam.d
			rate.vect.new<-rate.vect
			rate.vect.new[s]<-rate.vect[s]*r
			newedge<-rate.calc(edge,rate.vect.new,assign.d)
			if (s==1) prop="a_root.d" else prop="a_rate.d"
		return(list(newedge=newedge,rate.vect.new=rate.vect.new,HP=r,prop=prop))			
		
		}
##_________________________		
		
		adjustrate.e<-function(edge,rext.vect,assign.e,rootprob.e,lam.e){
			if (length(rext.vect)==1) s=1 
			else if (length(rext.vect)>1 & runif(1)>(1-rootprob.e)) s=1
			else s=sample(2:length(rext.vect),1)
			rext.vect.new<-rext.vect
			rext.vect.new[s]<-add.rate(rext.vect.new[s],lam.e*(runif(1)-0.5),0,1)		
			newedge<-rate.calc.e(edge,rext.vect.new,assign.e)			
			if (s==1) prop="a_root.e" else prop="a_rate.e"
			r=1
		return(list(newedge=newedge,rext.vect.new=rext.vect.new,HP=r,prop=prop))
		}
		
##--------------

add<-function(K,N){
	
	(K+1)/(N-K)	
}

subtract<-function(N,K){
	N-(K+1)/K
}