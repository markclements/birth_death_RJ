MCMC.BD<-function(phy,div,ngen=1000000,sfreq=100,root.rate=0.01,root.rate.e=0,e=3,lam.d=1,lam.e=0.5,lamsplit.e=0.5,msprob.d=0.25,msprob.e=0.25,rateprob.d=0.25,rateprob.e=0.25,rootprob.d=0.5,rootprob.e=0.5,internal.only=FALSE,assign.d=FALSE,assign.e=FALSE, prior.only=FALSE, directory="~/Desktop/trials",generate.log=TRUE){
	require(ape)
### Prepare data	
	edge<-build.edge(phy,div)
	if (generate.log){
	setwd(directory)
	dir.create(paste("RJ_BD_run",date()))
	z<-file.info(dir())
	zq<-which(z[,5]==max(z[,5]))
	setwd(rownames(z[zq,]))
	}
	if (sum(c(msprob.d,msprob.e,rateprob.d,rateprob.e))>1) stop("the probabilities must sum to 1")
	if (1-rootprob.d<0 || 1-rootprob.e<0) stop("the probabilies must be less than or equal to 1")
###--------------------------------------------------------------#####
# initiallize chain
	rate.vect<-root.rate
	rext.vect<-root.rate.e
	edge<-rate.calc(edge,rate.vect,assign.d)
	edge<-rate.calc.e(edge,rext.vect,assign.e)
	lnLik<-Lnlikelihood(edge,internal.only)
# Proposal tracking?? 

	
# BEGIN CALCULATION
	count.i<-1
	for(i in (1:ngen)) {  	
	  	while(1) {
	  	z<-sample(1:4,1,prob=c(msprob.d,msprob.e,rateprob.d,rateprob.e))  	
			 if (z==1){					#adjust number of d rates
				xx<-splitORmerge.d(edge,rate.vect,e,internal.only,assign.d)	
				newedge<-xx$newedge
				rate.vect.new<-xx$rate.vect.new
				HP<-xx$HP
				rext.vect.new=rext.vect
			  	prop<-xx$prop
			  	break()	
			  }	
		  	else if (z==2){					#adjust number of e rates
			  	xx<-splitORmerge.e(edge,rext.vect,e,internal.only,assign.e,lamsplit.e)
			  	newedge<-xx$newedge
			  	rext.vect.new<-xx$rext.vect.new
			  	HP<-xx$HP
			  	rate.vect.new<-rate.vect	
		  		prop<-xx$prop
		  		break()	
		  	}
		  	else if (z==3){		#adjust d rates, including root
		  		xx<-adjustrate.d(edge,rate.vect,assign.d,rootprob.d,lam.d)
		  		newedge<-xx$newedge
				rate.vect.new<-xx$rate.vect.new
				HP<-xx$HP
				rext.vect.new=rext.vect
		  		prop<-xx$prop
		  		break()	
		  	}
		  	else if (z==4){	#adjust e rates, including root 
				xx<-adjustrate.e(edge,rext.vect,assign.e,rootprob.e,lam.e)	  	
		  		newedge<-xx$newedge
			  	rext.vect.new<-xx$rext.vect.new
			  	HP<-xx$HP
			  	rate.vect.new<-rate.vect
				prop<-xx$prop
				break()			  	
		  	}
	  	}
	  	if (prior.only) newlnLik<-0
	  	else newlnLik<-Lnlikelihood(newedge,internal.only)
	  	lr <- newlnLik-lnLik
		ratio <- exp(lr)*HP
		
		if (runif(1) <= ratio) {			
			decision="accept"
			lnLik<-newlnLik
			edge<-newedge
			rate.vect<-rate.vect.new
			rext.vect<-rext.vect.new
		} else {							
			decision="reject" 
		}
		
		
	if (i==1 || i%%sfreq==0){
	#plot(phy,cex=0.5,no.margin=TRUE)
	#edgelabels(text=round(edge[,6],3),adj=c(1,0),frame="n",col="black",cex=0.5)
	#edgelabels(text=round(edge[,9],3),adj=c(1,1),frame="n",col="red",cex=0.5)
	#edgelabels(round(edge[,9],3),frame="n",cex=0.6)
	#print(edge)
	#cat(logLik,"\n", edge,"\n")
	if (generate.log){
	cat(count.i,lnLik,prop,decision,"\n", sep=" ", file="trace", append=TRUE)
	cat(edge[,6],"\n", file="divrates", append=TRUE)
	cat(edge[,9],"\n", file="rextrate", append=TRUE)
	cat(edge[,7],"\n", file="divshift", append=TRUE)
	cat(edge[,10],"\n",file="rextshift",append=TRUE)
	count.i<-count.i+1
			}
		}
	}
}
