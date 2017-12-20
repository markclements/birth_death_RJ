# set it up to do 1 time shift first, then work on random shifts per lineage birth
require(ape)

tree.sim<-function(b,d,time.switch,time.stop){
	edge<-initial.edge(time.switch)
		if (time.switch[1]==0) N<-2
		else N<-1
		repeat{ #main simulation loop
			edge<-evolve.lineage(edge,b,d,N,time.stop)
			if (sum(is.na(edge[,3]))==0) {
			break;}	
			if (N<=length(time.switch) & max(edge[,4]>=time.switch[N])) {
				x<-which(edge[,4]>=time.switch[N] & is.na(edge[,3]))
				edge[sample(x,1),5]<-N+1
				N<-N+1
				next;
				}
			}
		edge
	tr<-make.phylo(edge)
	tr$birth.p<-b
	tr$death.p<-d
	tr$time.swith<-time.switch
	tr$time.stop<-time.stop
	tr
}

evolve.lineage<-function(edge,b,d,N,time.stop){
	for (i in 1:N){ 
	alive<-which(edge[,5]==i & is.na(edge[,3]))
	if (length(alive)>0) {
		if (length(alive)>1){
			add<-sample(alive,1)
		}
		else add<-alive
		tn<-max(edge[alive,4])
		t<-tn+rexp(1,(length(alive)*(b[i]+d[i])))
		cat(i," ",t," ", tn," ", alive," ")
		if (t>=time.stop) {
			#t<-time.stop
			cat("end"," ","\n")
			edge<-end.edge(edge,alive,time.stop)
		}
		else if (runif(1)<=b[i]/(b[i]+d[i])) {# this creates a bifucation in the tree.
			cat("add"," ","\n")
			edge<-add.edge(edge,add,t)
		}
		else {# lineage goes extinct	
			cat("ext"," ","\n")
			edge<-end.edge(edge,add,t)
			edge[add,6]<-NA
		}
	}
	else next;
    }
edge
}

add.edge<-function(edge,add,t){
	new.lin<-rbind(c(max(edge[,1])+2,0),c(max(edge[,1])+2,0))
	new.lin<-cbind(new.lin,c(NA,NA),c(t,t),c(edge[add,5],edge[add,5]),c(0,0))
	new.edge<-rbind(cbind(edge[,1]+1,edge[,2:6]),new.lin)
	min<-min(new.edge[,1])
	ind.int<-which(new.edge[,2]>=min)
	new.edge[ind.int,2]<-new.edge[ind.int,2]+1
	max<-max(new.edge[,1])
	new.edge[add,2]<-max
	new.edge[add,3]<-t-new.edge[add,4]
	new.edge[which(new.edge[,2]<=min),2]<-1:(min-1)
  	new.edge
	}

end.edge<-function(edge,add,t){
	edge[add,3]<-t-edge[add,4]
	edge
	}

initial.edge<-function(time.switch){
	edge<-rbind(c(3,1),c(3,2)) # starter edge
	if (time.switch[1]==0) edge<-cbind(rbind(c(3,1),c(3,2)),c(NA,NA),c(0,0),c(1,2),c(0,0)) # edge,br-length,time of lineage birth,b-d switch column
	else edge<-cbind(rbind(c(3,1),c(3,2)),c(NA,NA),c(0,0),c(1,1),c(0,0)) # edge,br-length,time of lineage birth,b-d switch column
	edge
	}


make.phylo<-function(edge){#makes simulated edge matrix into phylo object.
		x=edge[,1:2]
		edge.length<-edge[,3]
		tip.label<-1:(min(edge[,1])-1)
		div.switch<-edge[,5]
		extinct<-edge[,6]
		Nnode<-length(unique(edge[,1]))
		mode(tip.label) <- "character"
    	mode(x)<-"integer"
    	mode(Nnode)<-"integer"
    	obj <- list(edge = x, edge.length = edge.length, tip.label=tip.label,Nnode=Nnode, div.switch=div.switch, extinct=extinct)
    	class(obj) <- "phylo"
    	#obj<-read.tree(text=write.tree(obj))	
    obj
}


exp.n<-function(b,d,t){ # expected # of species given clade age. 
	beta<-(b*(exp((b-d)*t)-1))/(b*(exp((b-d)*t)-d))
	n<-(1-beta)/((beta-1)^2)
	return(n)
	}	
	
# returns the heights of each node
# written by Liam J. Revell 2011

nodeHeights<-function(tree){
	# compute node heights
	root<-length(tree$tip)+1
	X<-matrix(NA,nrow(tree$edge),2)
	for(i in 1:nrow(tree$edge)){
		if(tree$edge[i,1]==root){
			X[i,1]<-0.0
			X[i,2]<-tree$edge.length[i]
		} else {
			X[i,1]<-X[match(tree$edge[i,1],tree$edge[,2]),2]
			X[i,2]<-X[i,1]+tree$edge.length[i]
		}
	}
	return(X)
}	

# get all the extant/extinct tip names
# written by Liam J. Revell 2012

getExtant<-function(tree,tol=1e-8){
	H<-nodeHeights(tree)
	tl<-max(H)
	x<-which(H[,2]>=(tl-tol))
	y<-tree$edge[x,2]
	y<-y[y<=length(tree$tip)]
	z<-tree$tip.label[y]
	return(z)
}

#getExtinct<-function(tree,tol=1e-8) setdiff(tree$tip.label,getExtant(tree,tol))
	
#t<-tree.sim(c(0.5,1,3),c(0,0,0),c(2,4),6)