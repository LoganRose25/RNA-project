#### Functions for use in analysing the structure profiles for RNA #####
### Functions ####
#dx=0.1
extend.by.zeros <- function(v,N){
  if(length(v)>N){
    w <- v
  } else {
    Z <- N-length(v)
    w <- c(v,integer(Z))
  }
  return(w)
}

extend.matrix.by.zeros <- function(Mat,N){
  l <- length(Mat[1,])
  if(l>N){
    w <- Mat
  }  else {
    Z <- t(apply(Mat,1,extend.by.zeros,N))
  }
  return(Z)
}

euclidean.length <- function(u){
  d <- sqrt(sum((u)^2))
}

euclidean.distance <- function(u,v){
  d <- sqrt(sum((u-v)^2))
}

files_to_matrix <- function(files){
  
}
max_length <- function(files,homology=1,thres=3.0,dx=0.1){
  M<-0
  hom <- homology+1
  pb <- txtProgressBar(min=1, max=length(files),style=3)
  i <- 0
  for (x in files){
    i <- i+1
    setTxtProgressBar(pb,i)
    mydata <- read.csv(x, header=FALSE) # load file
    ### Compute persistence diagram for mydata ###
    pd <- diagram(mydata, 'point-cloud', dim_max=homology, threshold=thres)
    if((length(pd$pairs[[hom]][,1]))>0){
      pl <- landscape(pd$pairs[[hom]], exact=FALSE, dx=dx, min_x=0, max_x=thres)
      v<-pl$getInternal()[,,2]
      pl_vec<-c(t(v))
      M<-max(M,length(pl_vec))
    } else {
      M<-M
    }
  }
  close(pb)
  return(M)
}

files_to_matrix <- function(files){
  M <- 0
  outdata <- vector(mode = "numeric")
  for(x in files){
    mydata <- read.csv(x, header=FALSE) # load file
    M <- max(M,length(mydata[1,]))
    mydata <- t(apply(mydata,1,extend.by.zeros,M))
    if(length(outdata)>0){
      outdata <- t(apply(outdata,1,extend.by.zeros,M))
    }
    outdata <- rbind(outdata,mydata)
  }
  return(outdata)
}

files_to_vector <- function(files){
  M <- 0
  outdata <- c()
  for(x in files){
    mydata <- read.csv(x, header=FALSE) # load file
    mydata <- c(t(mydata))
    M <- max(M,length(mydata))
    mydata <- extend.by.zeros(mydata,M)
    if(length(outdata)>0){
      outdata <- t(apply(outdata,1,extend.by.zeros,M))
    }
    outdata <- rbind(outdata,mydata)
  }
}

landscapes_to_matrix <- function(files,hom,M,dx=0.1){
  pl_mat <- vector(mode = "numeric")
  for (x in files){
    mydata <- read.csv(x, header=FALSE) # load file
    ### Compute persistence diagram for mydata ###
    pd <- diagram(mydata, 'point-cloud', dim_max=1, threshold=3.0)
    if((length(pd$pairs[[hom]][,1]))>0){
      pl <- landscape(pd$pairs[[hom]], exact=FALSE, dx=dx, min_x=0, max_x=3)
      v<-pl$getInternal()[,,2]
      pl_vec<-c(t(v))
      pl_vec<-extend.by.zeros(pl_vec,M)
      pl_mat<-rbind(pl_mat,pl_vec)
    } else {
      pl_vec<-integer(M)
      pl_mat<-rbind(pl_mat,pl_vec)
    }
  }  
  return(pl_mat)
}


get_list_landscapes_shifted_prof <- function(files,homology=1,thres=3.0,dx=0.1){
  ## This function loads structure profile data in as a list of dataframes
  ## containing the landscapes and average distances.  
  ##
  ## The intended data files should be csv files where each row corresponds to a 
  ## weighted structure profile and should be contained within a hypercube. The 
  ## last row of the csv file should be the native structure for the RNA strand. 
  ##
  ## The function calculates the distance for each structure profile to the 
  ## native one as well as the average distance. It then computes the persistence
  ## diagram followed by the persistence landscape.
  ##
  ## The output is a list where each object is a dataframe containing the landscape 
  ## and average distance for individual structure profiles
  pl_list <- list()
  hom <- homology + 1
  J<-0
  ave_dist <- c()
  label <- c()
  for (x in files){
    mat <- read.csv(x,header=FALSE,stringsAsFactors=FALSE) #read in data
    mat <- as.matrix(mat) #make sure data is interpreted as a matrix
    truth <- mat[nrow(mat),] #identify 
    mat.adjust <- t(apply(mat, 1, function(u) u-truth))
    distances <- apply(mat.adjust, 1, euclidean.length)
    mydata <- mat.adjust[-nrow(mat.adjust),]
    avg.dist <- mean(distances[-nrow(mat.adjust)])
    pd <- diagram(mydata, 'point-cloud', dim_max=homology, threshold=thres)
    if((length(pd$pairs[[hom]][,1]))>0){
      J <- J+1
      pl <- landscape(pd$pairs[[hom]], exact=FALSE, dx=dx, min_x=0, max_x=thres)
      pl_list[[J]] <- pl
      ave_dist <- c(ave_dist,avg.dist)
      label <- c(label,basename(x))
    } 
  }
  names(ave_dist)<-label
  names(pl_list)<-label
  pl_list[[J+1]] <- ave_dist
  return(pl_list)
}


get_matrix_landscapes_shifted <- function(files,homology=1,M,thres=3.0,dx=0.1,Norm=TRUE,Avg_distance=FALSE,reduced=TRUE){
  ## This function loads structure profile data in as a matrix.  
  ##
  ## The intended data files should be csv files where each row corresponds to a 
  ## weighted structure profile and should be contained within a hypercube. The 
  ## last row of the csv file should be the native structure for the RNA strand. 
  ##
  ## The function calculates the distance for each structure profile to the 
  ## native one as well as the average distance. It then computes the persistence
  ## diagram followed by the persistence landscape.
  ##
  ## The output is a matrix whose first column is the average distance repeated
  ## while the remaining columns correspond to the persistence landscape.
  pl_mat <- vector(mode = "numeric")
  label <- c()
  hom = homology+1
  pb = txtProgressBar(min = 1, max = length(files), style=3)
  i <- 0
  for (x in files){
    i <- i+1
    setTxtProgressBar(pb,i)
    label <- c(label,basename(x))
    mat <- read.csv(x,header=FALSE) #read in data
    mat <- as.matrix(mat) #make sure data is interpreted as a matrix
    d <- dist(mat)
    truth <- mat[nrow(mat),] #extract native structure of particular RNA
    mat.adjust <- t(apply(mat, 1, function(u) u-truth)) #shift data to be centered at native structure
    distances <- apply(mat.adjust, 1, euclidean.length) #compute norms for each point (i.e. each)
    mydata <- mat.adjust[-nrow(mat.adjust),] #remove native structure from data
    d <- dist(mydata)
    avg.dist <- mean(distances[-nrow(mat.adjust)]) #compute average distance
    pd <- diagram(mydata, 'point-cloud', dim_max=homology, threshold=thres) # compute persistence diagrams
    if((length(pd$pairs[[hom]][,1]))>0){ #compute persistence landscapes
      if(reduced&&hom==0){
        pdpairs <- pd$pairs[[hom]][-nrow(pd$pairs[[hom]]),]
      } else {
        pdpairs <- pd$pairs[[hom]]
      }
      pl <- landscape(pdpairs, exact=FALSE, dx=dx, min_x=0, max_x=thres)
      pl_vec <-pl$getInternal()[,,2] 
      #pl_vec <-c(avg.dist,t(pl_vec))
      pl_vec<-extend.by.zeros(pl_vec,M)
     # pl_mat<-rbind(pl_mat,pl_vec)
    } else {
      pl_vec<-integer(M)
    }
    if(Norm){
      norm <- euclidean.length(pl_vec)
      pl_vec <- c(norm,pl_vec)
    }
    if(Avg_distance){
      pl_vec <-c(avg.dist,pl_vec)
    }
    pl_mat<-rbind(pl_mat,pl_vec)
    pl_mat <- as(pl_mat,"sparseMatrix")
  }
  close(pb)
  rownames(pl_mat) <- label
  return(pl_mat)
}

landscapes_to_list <- function(files,hom,dx=0.1){
  pl_list <- list()
  for (x in files){
    mydata <- read.csv(x, header=FALSE) # load file
    ### Compute persistence diagram for mydata ###
    pd <- diagram(mydata, 'point-cloud', dim_max=1, threshold=3.0)
    if((length(pd$pairs[[hom]][,1]))>0){
      pl <- landscape(pd$pairs[[hom]], exact=FALSE, dx=dx, min_x=0, max_x=3)
      pl_list <- append(pl_list,pl)
    } else {
      pl_list <- append(pl_list,NA)
    }
  }  
  return(pl_mat)
}

sum_of_landscapes_list <- function(list,hom){
  N <- length(list)
  pl_sum <- list[[1]]
  for(i in 2:N){
    pl <- list[[i]]
    pl_sum <- PLsum(pl_sum,pl)
  }
  return(pl_sum)
}

subsetting_list <- function(list,range=c(0,5)){
  subsetted_list <- list()
  N<-length(list)
  dist <- list[[N]]
  i1 <- which(dist>range[1])
  i2 <- which(dist[i1]<range[2])
  index <- names(i2)
  distances <- dist[index]
  N <- length(index)
  L <- list[index]
  s <- sum_of_landscapes_list(L,2)
  s <- PLscale(1/N,s)
  subsetted_list[["names"]] <- index
  subsetted_list[["average_distances"]] <- distances
  subsetted_list[["average_landscape"]] <- s
  return(subsetted_list)
}


sum_of_landscapes <- function(files,hom,dx=0.1){
  row_numA0 <-vector(mode = "numeric")
  row_numA1 <-vector(mode = "numeric")
  
  mydata <- read.csv(files[1], header=FALSE) # load file
  pd <- diagram(mydata, 'point-cloud', dim_max=1, threshold=3.0)
  row_numA0<-append(row_numA0,nrow(pd$pairs[[1]]))
  row_numA1<-append(row_numA1,nrow(pd$pairs[[2]])) 
  pl <- landscape(pd$pairs[[hom]], exact=FALSE, dx=dx, min_x=0, max_x=3)
  PLplot(pl)
  pl_sum <- pl
  for (x in files[-1]){
    mydata <- read.csv(x, header=FALSE) # load file
    ### Compute persistence diagram for mydata ###
    pd <- diagram(mydata, 'point-cloud', dim_max=1, threshold=3.0)
    row_numA0<-append(row_numA0,nrow(pd$pairs[[1]]))
    row_numA1<-append(row_numA1,nrow(pd$pairs[[2]])) 
    if((length(pd$pairs[[hom]][,1]))>0){
      pl <- landscape(pd$pairs[[hom]], exact=FALSE, dx=dx, min_x=0, max_x=3)
      pl_sum <- PLsum(pl_sum,pl)
    }
  }  
  return(pl_sum)
}

random_distances <- function(mat,grp1size=0){
  L <- length(mat[,1])
  if(grp1size==0){
    L1 <- floor(L/2)
  }else{
    L1 <- grp1size
  }
  S <- sample(L)
  G1 <- mat[S[1:L1],]
  G2 <- mat[S[(L1+1):L],]
  d <- sum((colMeans(G1)-colMeans(G2))^2)
  return(d)
}

permutation_test_mat <- function(vec1,vec2,perm){
  a1 <- colMeans(vec1)
  a2 <- colMeans(vec2)
  d0 <- sum((a1-a2)^2)
  combined <- rbind(vec1,vec2)
  L1 <- length(vec1[,1])
  L2 <- length(vec2[,1])
  L <- length(combined[,1])
  d <- c()
  pb <- txtProgressBar(min=0, max=perm, style=3)
  j <- 0
  for(i in 1:perm){
    Sys.sleep(0.1)
    d1 <- random_distances(combined,grp1size=L1)
    d[i] <- d1-d0
    j<-j+1
    setTxtProgressBar(pb, j)
  }
  close(pb)
  
  p1 <- sum(d>0)/perm
  return(p1)
}



pers_lands <- function(persistence_diagram,hom,dx=0.1){
  pd <- persistence_diagram
  if((length(pd$pairs[[2]][,1]))>0){
    pl <- landscape(pd$pairs[[hom]],exact=FALSE, dx=dx, min_x=0, max_x=3)
    pl_vec <- pl$getInternal[,,2]
  } else {
    pl_vec <- 0
  }
  return(pl_vec)
}

landscape_length <- function(persistence_landscape){
  pl_vec <- persistence_landscape
  pl_vec <- c(t(pl_vec))
  length(pl_vec)
}

max_landscape_length <- function(persistence_landscape_list){
  lengths <- lapply(persistence_landscape_list,landscape_length)
  max_length <- max(lengths)
}