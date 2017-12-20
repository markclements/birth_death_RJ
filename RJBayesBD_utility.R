
order.d<-function(phy,div){	
	x<-order(phy$tip.label)
	y<-order(div[,1])
	z<-numeric(length(x))
	z[x]<-div[y,2]
	names(z)<-phy$tip.label
	return(z)
}


trace.plot<-function(file="trace",burn.in=0){# produce a trace plot, check convergence
	x<-read.table(file=file)
	max=length(x[,1])
	if (burn.in==0) burn.in<-1	
	plot(burn.in:max,x[burn.in:max,2],"l",ylab="Likelihood",xlab="Generations")	
}

mixing<-function(file="trace"){# look at trace file to see acceptance rates of proposals
	x<-read.table(file=file)	
	zz<-matrix(NA,12,2)
	for (i in 1:12){
	x1<-x[which(x[,3]==i),]
	total<-length(x1[,4])
	perc.acc<-(sum(x1[,4])/total)*100	
	zz[i,]<-c(total,round(perc.acc,2))
	}
	rownames(zz)<-c("must split div","must merge div","regular merge div","regular split div","change div rate mod", "change div root rate","must split ext" ,"must merge ext" ,"reg merg ext", "reg split ext" , "change ext rate mod","change ext root rate")
	colnames(zz)<-c("# proposals", "% acceptance")
	return(zz)
	}

rs.calc<-function(file=file, burn.in=1000,cut.off=50){#calcuate posterior prob. of rate shift
	x<-read.table(file=file,skip=burn.in)	
	rs<-colSums(x,na.rm=TRUE)/length(x[,1])*100
	rs[rs<cut.off]<-NA
	return(rs)
	}
	

summary.rate<-function(file=file,fun=mean, burn.in=1000){#calc summary stats for each lineage
	x<-read.table(file=file,skip=burn.in)
	z<-apply(x,2,FUN=fun)
	return(z)
	}

rs.density<-function(file=file,burn.in=1000){
	x<-read.table(file=file,skip=burn.in)
	zz<-apply(x,1,function(x) length(which(x==1)))
	return(zz)
	
}



#general phylogenetic plotting utility, given a named vector or data.frame of values that can be associated with phy$edge[,2]
#author: JM EASTMAN 2010
#note: small values are given bluish hues, large values reddish hues; median values are given gray hues

branchcol.plot <- function(phy, cur.rates, color.length=8, digits=3, plot=TRUE, legend=TRUE, legend.title="", ...) {
	require(ape)
	warning("Rates assumed to be ordered as in 'phy$edge'")
	phy<-reorder(phy,"c")
	mm=cur.rates-median(cur.rates)
	colors<-colorRampPalette(c("blue","grey","red"))
	cce<-colors(2*color.length+1)
	e.seq=seq(-max(abs(mm)),max(abs(mm)),length=2*color.length+1)
	colors.branches=sapply(mm, function(x) cce[which(min(abs(e.seq-x))==abs(e.seq-x))])
	#lseq=e.seq+median(cur.rates)
	lseq=seq(min(cur.rates), max(cur.rates), length=color.length)
	lcce=cce[round(seq(1, length(cce), length=color.length))]
	lseq=rev(lseq)
	
	if(plot) {
		plot.phylo(phy, cex=0.5, edge.color=colors.branches, ...)
		if(legend) {
			legend("topright", title=legend.title, cex=0.5, pt.cex=1, text.col="darkgray", 
				   legend = sprintf(paste("%", 2*digits, paste(digits, "f", sep=""), sep="."), lseq), pch=21, ncol=1, col = "darkgray", 
				   pt.bg = rev(lcce), box.lty="blank", border="white")
		}
	} else {
		return(list(col=colors.branches,legend.seq=lseq,legend.col=lcce))
	}
}





#cut(x,breaks=??,labels=FALSE)

rcalc<-function(e,ne,N,NN){
	x<-((exp(-ne))*(ne^N))/((exp(-(e)))*(e^N))
	return(x)
	}

rcalc<-function(e,ne,avge,elam){
	x<-((ne/e)^avge)*exp(-elam*(ne-e))
	return(x)
	}

betacalc<-function(rate,ext,time,div) {
	beta<-expm1(rate*time)/(exp(rate*time)-ext)
	lik<-log(1-beta) + ((div-1) * log(beta))
	return(lik)	
	}		
	
#FUN<-function(x) length(x)
#	x<-lapply(XXXX$rate.vect,FUN)
#	x<-unlist(x)
#	hist(x)