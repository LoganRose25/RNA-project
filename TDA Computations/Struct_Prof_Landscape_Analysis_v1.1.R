library(tdatools)
#setwd("/Users/matthewwheeler/Dropbox (UFL)/TDA_tools/SetA_no_shape/")
source("./struc_prof_tools_v1.0.R")

dx<-0.1


#### Files A #########
filesA <- list.files(path="./TDA_tools/SetA_no_shape", pattern="*no_shape.csv", full.names=TRUE, recursive=FALSE)

files<-filesA
N<-length(files)
M <- max_length(files)
pl_matA <- landscapes_to_matrix(files,2,M)
pl_sumA <- sum_of_landscapes(files,2)
pl_listA <- landscapes_to_list(files,2)
pl_averageA <- PLscale(1/N, pl_sumA)
PLplot(pl_averageA)

landscape.all <- pl_matA
column.means <- colMeans(pl_matA)
#landscape.short <- landscape.all[,column.means!=0]
#K <- length(landscape.short[1,])
#landscape.short <- landscape.short[,4:K]
print("Calculating PCA coordinates")
pca <-prcomp(landscape.all, retx=TRUE, center=FALSE)
print("Completed PCA calculation")
landscape.pca <- pca$x
landscape.eigenvec <- pca$rotation
vars <- apply(landscape.pca, 2, var)
props <- vars / sum(vars)
cumprop1 <- 100*cumsum(props)[1] # % variance explained by first 1 PCA coords
cumprop2 <- 100*cumsum(props)[2]
cumprop3 <- 100*cumsum(props)[3] # % variance explained by first 3 PCA coords
cumprop <- c(cumprop1,cumprop2,cumprop3)


landscape.center <- t(t(landscape.all)-column.means)
pca_cent <- prcomp(landscape.all, retx=TRUE, center=TRUE)
landscape.pca.center <- pca_cent$x
landscape.eigenvec.center <- pca_cent$rotation
vars <- apply(landscape.pca.center, 2, var)
props <- vars / sum(vars)
cumprop1 <- 100*cumsum(props)[1] # % variance explained by first 1 PCA coords
cumprop2 <- 100*cumsum(props)[2]
cumprop3 <- 100*cumsum(props)[3] # % variance explained by first 3 PCA coords
cumprop <- c(cumprop1,cumprop2,cumprop3)

plot(landscape.pca.center[,1:2])


filesAS <- list.files(path="./TDA_tools/SetA_with_shape", pattern="*with_shape.csv", full.names=TRUE, recursive=FALSE)

files<-filesAS[-473]
N<-length(files)
M <- max_length(files)
pl_matAS <- landscapes_to_matrix(files,2,M)
pl_sumAS <- sum_of_landscapes(files,2)
pl_listAS <- landscapes_to_list(files,2)
pl_averageAS <- PLscale(1/N, pl_sumAS)
PLplot(pl_averageAS)



#### Files B #########

filesB <- list.files(path="./TDA_tools/SetB_no_shape/", pattern="*no_shape.csv", full.names=TRUE, recursive=FALSE)

files<-filesB
N<-length(files)
M <- max_length(files)
pl_matB <- landscapes_to_matrix(files,2,M)
pl_sumB <- sum_of_landscapes(files,2)
pl_listB <- landscapes_to_list(files,2)
pl_averageB <- PLscale(1/N, pl_sumB)
PLplot(pl_averageB)

filesBS <- list.files(path="./TDA_tools/SetB_with_shape/", pattern="*with_shape.csv", full.names=TRUE, recursive=FALSE)

files<-filesBS
N<-length(files)
M <- max_length(files)
pl_matBS <- landscapes_to_matrix(files,2,M)
pl_sumBS <- sum_of_landscapes(files,2)
pl_listBS <- landscapes_to_list(files,2)
pl_averageBS <- PLscale(1/N, pl_sumBS)
PLplot(pl_averageBS)


#### Plotting Results #####

PLplot(pl_averageA)
title("Average Landscape of Group A")
PLplot(pl_averageAS)
title("Average Landscape of Group A with Shape")

PLplot(pl_averageB)
title("Average Landscape of Group B")
PLplot(pl_averageBS)
title("Average Landscape of Group B with Shape")


#### Statistical Analysis #########
M <- max(length(pl_matA[1,]),pl_matB[1,])

pl_matA <- t(apply(pl_matA,1,extend.by.zeros,M))
pl_matB <- t(apply(pl_matB,1,extend.by.zeros,M))

landscape.all <- rbind(pl_matA,pl_matB)
column.means <- colMeans(landscape.all)
#landscape.short <- landscape.all[,column.means!=0]
#K <- length(landscape.short[1,])
#landscape.short <- landscape.short[,4:K]
print("Calculating PCA coordinates")
pca <-prcomp(landscape.all, retx=TRUE, center=FALSE)
print("Completed PCA calculation")
landscape.pca <- pca$x
landscape.eigenvec <- pca$rotation
vars <- apply(landscape.pca, 2, var)
props <- vars / sum(vars)
cumprop1 <- 100*cumsum(props)[1] # % variance explained by first 1 PCA coords
cumprop2 <- 100*cumsum(props)[2]
cumprop3 <- 100*cumsum(props)[3] # % variance explained by first 3 PCA coords

cumprop <- c(cumprop1,cumprop2,cumprop3)

diff <- PLsum(pl_averageA,PLscale(-1,pl_averageB))
d0 <- PLinner(diff,diff)

#Using only first pca component to calculate significance
t_test <- t.test(landscape.pca[1:400,1],landscape.pca[401:800,1],alternative = "two.sided")
p1 <- t_test$p.value

#Kolmogorov-Smirnov Test
#KL Divergence
#Jensen (average)
#KL(x,y)+KL(y,x)-KL(x,y)*KL(y,x)

#Using first 7 pca components to calculate significance
# t <-(a1-a2/s^2(1/n1+1/n2))

a1 <- colMeans(landscape.pca[1:400,1:7])
a2 <- colMeans(landscape.pca[401:800,1:7])
a3 <- colMeans(landscape.pca[,1:7])
s <- sqrt(rowSums((t(landscape.pca[,1:7])-a3)^2)/799)
s2 <- sum(s^2)
t_score <- (a1-a2)/sqrt(s2*(1/200))
p_value <-pt(q=sqrt(sum(t_score^2)), df=398, lower.tail=TRUE)
1-p_value


## Permutation Test ##
# Calculate the average landscapes for the two groups. Then combine the 
# both groups into one collection and then split them randomly into two new equal
# groups A1 and A2 and calculate their average landscapes l1 and l2. Let d0 be
# the difference between the original averages and d1 the distance between l1
# and l2. Repeat this process 10,000 times. The p-value is the proportion of 
# cases in which d1 > d0.

L<-length(landscape.all[,1])
diff <- PLsum(pl_averageA,PLscale(-1,pl_averageB))
d0 <- PLinner(diff,diff)

pl_avgA <- c(t(pl_averageA$getInternal()[,,2]))
pl_avgB <- extend.by.zeros(c(t(pl_averageB$getInternal()[,,2])), N = length(pl_avgA))
d0a <-sum((pl_avgA-pl_avgB)^2)

ave1<-colMeans(landscape.all[1:400,])
ave2<-colMeans(landscape.all[401:800,])
d0b <- sum((ave1-ave2)^2)

M <- max(c(length(pl_matA[1,]),length(pl_matAS[1,]),length(pl_matB[1,]),length(pl_matBS[1,])))
pl_matA <- t(apply(pl_matA,1,extend.by.zeros,M))
pl_matAS <- t(apply(pl_matAS,1,extend.by.zeros,M))
pl_matB <- t(apply(pl_matB,1,extend.by.zeros,M))
pl_matBS <- t(apply(pl_matBS,1,extend.by.zeros,M))

pl_matAB <- rbind(pl_matA,pl_matB)
pl_matABS <- rbind(pl_matAS,pl_matBS)

p1 <- permutation_test_mat(pl_matA,pl_matAS,10000)
p2 <- permutation_test_mat(pl_matB,pl_matBS,10000)
p3 <- permutation_test_mat(pl_matA,pl_matB,10000)
p4 <- permutation_test_mat(pl_matAB,pl_matABS,10000)

d <- c()
perm <- 100 # Number of permutations
pb <- txtProgressBar(min=0, max=perm, style=3)
j <- 0
for(i in 1:perm){
  Sys.sleep(0.1)
  S<-sample(L)
  L1<-landscape.all[S[1:(L/2)],]
  L2<-landscape.all[S[((L/2)+1):L],]
  d1 <- sum((colMeans(L1)-colMeans(L2))^2)
  d[i] <- d1-d0b
  j<-j+1
  setTxtProgressBar(pb, j)
}
close(pb)

p1 <- sum(d>0)/perm

#### plot data ######
library(ggplot2)
x<-c(1:length(row_numA0))
dfA<-data.frame(x,row_numA0)
dfA1<-data.frame(x,row_numA1)

y<-c(1:length(row_numB0))
dfB<-data.frame(y,row_numB0)
dfB1<-data.frame(y,row_numB1)

ggplot(data=dfA,aes_q(x,row_numA0)) +
  geom_point() +
  geom_hline(yintercept = mean(row_numA0), color="blue")

ggplot(data=dfA1,aes_q(x,row_numA1)) +
  geom_point() +
  geom_hline(yintercept = mean(row_numA1), color="blue")

ggplot(dfA, aes(x=row_numA0)) + geom_histogram(binwidth=.5)
ggplot(dfA1, aes(x=row_numA1)) + geom_histogram(binwidth=10.0)

dat <- data.frame(set = factor(rep(c("A","B"), times=c(length(row_numA0),length(row_numB0)))), 
                  bars0 = c(row_numA0,row_numB0)
)

ggplot(data=dfB,aes_q(y,row_numB0)) +
  geom_point() +
  geom_hline(yintercept = mean(row_numB0), color="blue")

ggplot(data=dfB1,aes_q(y,row_numB1)) +
  geom_point() +
  geom_hline(yintercept = mean(row_numB1), color="blue")

ggplot(dfB, aes(x=row_numB0)) + geom_histogram(binwidth=.5)
ggplot(dfB1, aes(x=row_numB1)) + geom_histogram(binwidth=10.0)

#dat <- data.frame(set = factor(rep(c("A","B"), each=400)), 
 #                 bars0 = c(row_numA0,row_numB0)
#)
# View first few rows
head(dat)

## Basic histogram from the vector "row 0". Each bin is .5 wide.
## These both result in the same output:
ggplot(dat, aes(x=bars0, fill=set)) + geom_histogram(binwidth=.5)
# qplot(dat$rating, binwidth=.5)



library(plyr)
mu <- ddply(dat, "set", summarise, grp.mean=mean(bars0))
head(mu)

# Add mean lines
p<-ggplot(dat, aes(x=bars0, color=set)) +
  geom_histogram(fill="white", position="dodge")+
  geom_vline(data=mu, aes(xintercept=grp.mean, color=set),
             linetype="dashed")+
  theme(legend.position="top")
p

dat1 <- data.frame(set = factor(rep(c("A","B"), each=400)), 
                   bars1 = c(row_numA1,row_numB1)
)
# View first few rows
head(dat1)

## Basic histogram from the vector "rows 1". Each bin is .5 wide.
## These both result in the same output:
ggplot(dat1, aes(x=bars1, fill=set)) + geom_histogram(binwidth=.5)
# qplot(dat$rating, binwidth=.5)

library(plyr)
mu1 <- ddply(dat1, "set", summarise, grp.mean=mean(bars1))
head(mu1)

# Add mean lines
p1<-ggplot(dat1, aes(x=bars1, color=set)) +
  geom_histogram(fill="white", position="dodge")+
  geom_vline(data=mu1, aes(xintercept=grp.mean, color=set),
             linetype="dashed")+
  theme(legend.position="top")
p1
