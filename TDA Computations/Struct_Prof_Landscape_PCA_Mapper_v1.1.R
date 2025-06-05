library(tdatools)
library(tidyverse)
library(plotly)
library(umap)
library(Rtsne)

setwd('/Users/matthewwheeler/Dropbox (UFL)/Research/RNA Structure Profililng')
source("./struc_prof_tools_v1.0.R")

#setwd("/Users/matthewwheeler/Dropbox (UFL)/TDA_tools/SetA_no_shape/")

dx<-0.1

Load = "as_mat" 

#### Load the Data #########

files[3]

#### Files A #########

filesA <- list.files(path="./TDA_tools/SetA_no_shape", pattern="*no_shape.csv", full.names=TRUE, recursive=FALSE)
files<-filesA[file.size(filesA) > 0]
filesA <- files

N1<-length(files)
if(Load == "as_mat"){
  M <- max_length(files,thres=5)
  pl_matA <- get_matrix_landscapes_shifted(files,homology=1,M,thres=5)
}else{
  pl_listA <- get_list_landscapes_shifted_prof(files,thres=5)
}

acc <- read.csv("./TDA_tools/SetA_no_shape/distance_vs_accuracy_SetA_no_shape.csv",header=FALSE)
seq_names <- paste0(acc[,1],"_no_shape.csv")
rownames(acc) <- seq_names
pl_matacc <- pl_matA[seq_names,]
norms <- apply(pl_matacc[,2:ncol(pl_matacc)],1,euclidean.length)

A.df <-data.frame(
  row.names = rownames(pl_matacc),
  accuracy = acc[,4],
  min_distance = acc[,3],
  avg_distance = pl_matacc[,1],
  norms = norms,
  landscape = pl_matacc[,2:ncol(pl_matacc)]
)

pl_matAacc <- cbind(acc[,3:4],pl_matacc)



#pl_matA <- files_to_matrix(files)

#### Files A with Shape #########
filesAS <- list.files(path="./TDA_tools/SetA_with_shape", pattern="*with_shape.csv", full.names=TRUE, recursive=FALSE)
files<-filesAS[file.size(filesAS) > 0]
filesAS<-files
N2<-length(files)
if(Load == "as_mat"){
  M <- max_length(files,thres=5)
  pl_matAS <- get_matrix_landscapes_shifted(files,homology=1,M,thres=5)
}else{
  pl_listAS <- get_list_landscapes_shifted_prof(files)
}

acc <- read.csv("./TDA_tools/SetA_with_shape/distance_vs_accuracy_SetA_with_shape.csv",header=FALSE)
seq_names <- paste0(acc[,1],"_with_shape.csv")
rownames(acc) <- seq_names
pl_matacc <- pl_matAS[seq_names,]
pl_matASacc <- cbind(acc[,3:4],pl_matacc)

#### Files B #########
filesB <- list.files(path="./TDA_tools/SetB_no_shape", pattern="*no_shape.csv", full.names=TRUE, recursive=FALSE)
files<-filesB[file.size(filesB) > 0]
filesB<-files
N3<-length(files)
if(Load == "as_mat"){
  M <- max_length(files,thres=5)
  pl_matB <- get_matrix_landscapes_shifted(files,homology=1,M,thres=5)
}else{
  pl_listB <- get_list_landscapes_shifted_prof(files)
}

acc <- read.csv("./TDA_tools/SetB_no_shape/distance_vs_accuracy_SetB_no_shape.csv",header=FALSE)
seq_names <- paste0(acc[,1],"_no_shape.csv")
rownames(acc) <- seq_names
pl_matacc <- pl_matB[seq_names,]
pl_matBacc <- cbind(acc[,3:4],pl_matacc)

#### Files B with Shape #########
filesBS <- list.files(path="./TDA_tools/SetB_with_shape", pattern="*with_shape.csv", full.names=TRUE, recursive=FALSE)
files<-filesBS[file.size(filesBS) > 0]
filesBS<-files
if(Load == "as_mat"){
  M <- max_length(files,thres=5)
  pl_matBS <- get_matrix_landscapes_shifted(files,homology=1,M,thres=5)
}else{
  pl_listBS <- get_list_landscapes_shifted_prof(files)
}

acc <- read.csv("./TDA_tools/SetB_with_shape/distance_vs_accuracy_SetB_with_shape.csv",header=FALSE)
seq_names <- paste0(acc[,1],"_with_shape.csv")
rownames(acc) <- seq_names
pl_matacc <- pl_matBS[seq_names,]
pl_matBSacc <- cbind(acc[,3:4],pl_matacc)

### Lists ####

avg_distA <- pl_listA[[length(pl_listA)]]
avg_distAS <- pl_listAS[[length(pl_listAS)]]
avg_distB <- pl_listB[[length(pl_listB)]]
avg_distBS <- pl_listBS[[length(pl_listBS)]]

cutoffs <- c(0,2.5,5,7.5,10,100)

listA_all <- subsetting_list(pl_listA,range=c(cutoffs[1],cutoffs[6]))
listA_1<-subsetting_list(pl_listA,range=c(cutoffs[1],cutoffs[2]))
listA_2<-subsetting_list(pl_listA,range=c(cutoffs[2],cutoffs[3]))
listA_3<-subsetting_list(pl_listA,range=c(cutoffs[3],cutoffs[4]))
listA_4<-subsetting_list(pl_listA,range=c(cutoffs[4],cutoffs[5]))
listA_5<-subsetting_list(pl_listA,range=c(cutoffs[5],cutoffs[6]))

PLplot(listA_all$average_landscape)
PLplot(listA_1$average_landscape)
PLplot(listA_2$average_landscape)
PLplot(listA_3$average_landscape)
PLplot(listA_4$average_landscape)
PLplot(listA_5$average_landscape)





#### Imported as a matrix #####

M <- max(c(length(pl_matA[1,]),length(pl_matAS[1,]),length(pl_matB[1,]),length(pl_matBS[1,])))
pl_matA <- t(apply(pl_matA,1,extend.by.zeros,M))
pl_matAS <- t(apply(pl_matAS,1,extend.by.zeros,M))
pl_matB <- t(apply(pl_matB,1,extend.by.zeros,M))
pl_matBS <- t(apply(pl_matBS,1,extend.by.zeros,M))

# pl_matAdiff <- pl_matAS-pl_matA
# pl_matBdiff <- pl_matBS-pl_matB
#### Labeling #######
pl_matA.lab <- cbind(rep(0,nrow(pl_matA)),pl_matA)
pl_matAS.lab <- cbind(rep(1,nrow(pl_matAS)),pl_matAS)
pl_matB.lab <- cbind(rep(2,nrow(pl_matB)),pl_matB)
pl_matBS.lab <- cbind(rep(3,nrow(pl_matBS)),pl_matBS)

#rownames(pl_matA.lab) <- basename(filesA)
#rownames(pl_matAS.lab) <- basename(filesAS)
#rownames(pl_matB.lab) <- basename(filesB)
#rownames(pl_matBS.lab) <- basename(filesBS)

pl_matAdiff.lab <- cbind(rep(0,nrow(pl_matAdiff)),pl_matAdiff)
pl_matBdiff.lab <- cbind(rep(1,nrow(pl_matBdiff)),pl_matBdiff)
#### Landscape Averages #########

pl_sumA <- sum_of_landscapes(filesA,2)
pl_averageA <- PLscale(1/N1, pl_sumA)
PLplot(pl_averageA)

pl_sumAS <- sum_of_landscapes(filesAS,2)
pl_averageAS <- PLscale(1/N2, pl_sumAS)
PLplot(pl_averageAS)

pl_sumB <- sum_of_landscapes(filesB,2)
pl_averageB <- PLscale(1/N3, pl_sumB)
PLplot(pl_averageB)

pl_sumBS <- sum_of_landscapes(filesBS,2)
pl_averageBS <- PLscale(1/N4, pl_sumBS)
PLplot(pl_averageBS)

pl_aveA <- (1/nrow(pl_matA))*colSums(pl_matA[,3:M])
#### Combine #######

pl_matA_comb <- rbind(pl_matA.lab,pl_matAS.lab)
pl_matB_comb <- rbind(pl_matB.lab,pl_matBS.lab)
pl_mat_comb <- rbind(pl_matA.lab,pl_matAS.lab,pl_matB.lab,pl_matBS.lab)
pl_mat_comb <- pl_mat_comb[pl_mat_comb[,2]<10,]
pl_mat_comb <- pl_mat_comb[order(pl_mat_comb[,2],decreasing=TRUE),]

pl_matdiff <- rbind(pl_matAdiff.lab,pl_matBdiff.lab)

#### Analysis #######


##### Calculate PCA Coordinates for All Landscapes #####

landscape.all <- pl_mat_comb[,3:ncol(pl_mat_comb)]
column.means <- colMeans(landscape.all)
landscape.all <- t(t(landscape.all)-column.means)
#landscape.short <- landscape.all[,column.means!=0]
#K <- length(landscape.short[1,])
#landscape.short <- landscape.short[,4:K]
print("Calculating PCA coordinates")
pca <-prcomp(landscape.all, retx=TRUE, center=FALSE)
print("Completed PCA calculation")
landscape.pca.all <- pca$x

landscape.eigenvec <- pca$rotation
vars <- apply(landscape.pca.all, 2, var)
props <- vars / sum(vars)
cumprop1 <- 100*cumsum(props)[1] # % variance explained by first 1 PCA coords
cumprop2 <- 100*cumsum(props)[2]
cumprop3 <- 100*cumsum(props)[3]

pca.all.df <- data.frame(
  group = as.factor(pl_mat_comb[,1]),
  avg.dist = pl_mat_comb[,2],
  x = landscape.pca.all[,1],
  y = landscape.pca.all[,2],
  z = landscape.pca.all[,3],
  pca.coords = landscape.pca.all
)

##### Calculate PCA Coordinates for Group Pairs ####

##### Group A ####
landscape.all <- pl_matA_comb[,3:ncol(pl_matA_comb)]
column.means <- colMeans(landscape.all)
landscape.all <- t(t(landscape.all)-column.means)
print("Calculating PCA coordinates")
pca <-prcomp(landscape.all, retx=TRUE, center=FALSE)
print("Completed PCA calculation")
landscape.pca <- pca$x

landscape.eigenvec <- pca$rotation
vars <- apply(landscape.pca.all, 2, var)
props <- vars / sum(vars)
cumprop1 <- 100*cumsum(props)[1] # % variance explained by first 1 PCA coords
cumprop2 <- 100*cumsum(props)[2]
cumprop3 <- 100*cumsum(props)[3]

pca.A.df <- data.frame(
  group = as.factor(pl_matA_comb[,1]),
  avg.dist = pl_matA_comb[,2],
  x = landscape.pca[,1],
  y = landscape.pca[,2],
  z = landscape.pca[,3],
  pca.coords = landscape.pca
)
##### Group B ####
landscape.all <- pl_matB_comb[,3:ncol(pl_matB_comb)]
column.means <- colMeans(landscape.all)
landscape.all <- t(t(landscape.all)-column.means)
#landscape.short <- landscape.all[,column.means!=0]
#K <- length(landscape.short[1,])
#landscape.short <- landscape.short[,4:K]
print("Calculating PCA coordinates")
pca <-prcomp(landscape.all, retx=TRUE, center=FALSE)
print("Completed PCA calculation")
landscape.pca <- pca$x

landscape.eigenvec <- pca$rotation
vars <- apply(landscape.pca.all, 2, var)
props <- vars / sum(vars)
cumprop1 <- 100*cumsum(props)[1] # % variance explained by first 1 PCA coords
cumprop2 <- 100*cumsum(props)[2]
cumprop3 <- 100*cumsum(props)[3]

pca.B.df <- data.frame(
  group = as.factor(pl_matB_comb[,1]),
  avg.dist = pl_matB_comb[,2],
  x = landscape.pca[,1],
  y = landscape.pca[,2],
  z = landscape.pca[,3],
  pca.coords = landscape.pca
)
#### Calculate PCA Coordinates for Differences ####

landscape.all <- pl_matdiff[,3:ncol(pl_matdiff)]
column.means <- colMeans(landscape.all)
landscape.all <- t(t(landscape.all)-column.means)
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
cumprop3 <- 100*cumsum(props)[3]

pca.diff.df <- data.frame(
  group = as.factor(pl_matdiff[,1]),
  avg.dist = pl_matdiff[,2],
  x = landscape.pca[,1],
  y = landscape.pca[,2],
  z = landscape.pca[,3],
  pca.coords = landscape.pca
)
#### Visualization ######

### PCA coords of Landscapes ###

## All Coords 
ggplot(pca.all.df, aes(x=x,y=y,color=group)) +
  geom_point()

ggplot(pca.all.df, aes(x=x,y=y,color=avg.dist),) +
  geom_point()+
  scale_color_gradient(low="blue", high="red")

ggplot(pca.all.df, aes(x=x,y=y,color=group,alpha=avg.dist)) +
  geom_point()

## 3D Plot
fig <- plot_ly(pca.all.df, x = ~x, y = ~y, z = ~z, size =~avg.dist,
               marker = list(color = ~group, #colorscale = c('blue', 'red'), 
                             showscale = TRUE))
fig <- fig %>% add_markers()

fig



## Diff
ggplot(pca.diff.df, aes(x=x,y=y,color=group)) +
  geom_point()

ggplot(pca.diff.df, aes(x=x,y=y,color=group,alpha=avg.dist)) +
  geom_point()

## 3D Plot
fig <- plot_ly(pca.all.df, x = ~x, y = ~y, z = ~z, size =~avg.dist,
               marker = list(color = ~group, #colorscale = c('blue', 'red'), 
                             showscale = TRUE))
fig <- fig %>% add_markers()

fig

## Grp A
ggplot(pca.A.df, aes(x=x,y=y,color=group)) +
  geom_point()

ggplot(pca.A.df, aes(x=x,y=y,color=group,alpha=avg.dist)) +
  geom_point()

## 3D Plot
fig <- plot_ly(pca.A.df, x = ~x, y = ~y, z = ~z, size =~avg.dist,
               marker = list(color = ~group, #colorscale = c('blue', 'red'), 
                             showscale = TRUE))
fig <- fig %>% add_markers()

fig

fig <- plot_ly(pca.A.df, x = ~x, y = ~y, z = ~z, size =~avg.dist,
               marker = list(color = ~avg.dist, colorscale = c('blue', 'red'), 
                             showscale = TRUE))
fig <- fig %>% add_markers()

fig

##Grp B
ggplot(pca.B.df, aes(x=x,y=y,color=group)) +
  geom_point()

ggplot(pca.B.df, aes(x=x,y=y,color=group,alpha=avg.dist)) +
  geom_point()

fig <- plot_ly(pca.B.df, x = ~x, y = ~y, z = ~z, size =~avg.dist,
               marker = list(color = ~group, #colorscale = c('blue', 'red'), 
                             showscale = TRUE))
fig <- fig %>% add_markers()

fig

### Norms ####
norms <- apply(pl_mat_comb[,3:(M+1)],1,euclidean.length)
norms.df <- data.frame(
  group = as.factor(pl_mat_comb[,1]),
  ave.dist = pl_mat_comb[,2],
  norm = norms
)

ggplot(norms.df, aes(x=norm, color = group)) +
  geom_histogram(fill="white",alpha=0.5,position="identity")

ggplot(norms.df, aes(x=ave.dist, y=norm, color = group)) +
  geom_point()

 
landscape.all <- pl_mat_comb[,2:ncol(pl_mat_comb)]
column.means <- colMeans(landscape.all)
landscape.all <- t(t(landscape.all)-column.means)
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
cumprop3 <- 100*cumsum(props)[3]

df <- data.frame(
  pca <- landscape.pca,
  x <- landscape.pca[,1],
  y <- landscape.pca[,2],
  z <- landscape.pca[,3],
  group <- as.factor(pl_mat_comb[,1])
)
gg <- ggplot(df,aes(x=x,y=y,color=group)) +
  geom_point()
pl <- plot_ly(x=df$x, y=df$y, z=df$z, type="scatter3d", mode="markers", color=df$group)





landscape.umap <- umap(landscape.all[,3:ncol(landscape.all)])

data.umap <- landscape.umap$data
df <- data.frame(
  pca <- data.umap,
  x <- data.umap[,1],
  y <- data.umap[,2],
  z <- data.umap[,3],
  group <- as.factor(pl_mat_comb[,1])
)

average_distances <- landscape.all[,1]
layout <- landscape.umap[["layout"]] 
layout <- data.frame(layout) 
final <- cbind(layout, group, average_distances) 

fig1 <- plot_ly(final, x = ~X1, y = ~X2, color = ~average_distances, #colors = c('#636EFA','#EF553B','#00CC96','orange'), 
               type = 'scatter', mode = 'markers')%>%  
  layout(
    plot_bgcolor = "#e5ecf6",
    legend=list(title=list(text='avg dist')), 
    xaxis = list( 
      title = "0"),  
    yaxis = list( 
      title = "1")) 

fig2 <- plot_ly(final, x = ~X1, y = ~X2, color = ~group, #colors = c('#636EFA','#EF553B','#00CC96','orange'), 
                type = 'scatter', mode = 'markers')%>%  
  layout(
    plot_bgcolor = "#e5ecf6",
    legend=list(title=list(text='group')), 
    xaxis = list( 
      title = "0"),  
    yaxis = list( 
      title = "1")) 


iris.umap = umap(iris.data, n_components = 3, random_state = 15) 
layout <- iris.umap[["layout"]] 
layout <- data.frame(layout) 
final <- cbind(layout, iris$Species) 

fig2 <- plot_ly(final, x = ~X1, y = ~X2, z = ~X3, color = ~iris$Species, colors = c('#636EFA','#EF553B','#00CC96')) 
fig2 <- fig2 %>% add_markers() 
fig2 <- fig2 %>% layout(scene = list(xaxis = list(title = '0'), 
                                     yaxis = list(title = '1'), 
                                     zaxis = list(title = '2'))) 

fig 

ggplot(df,aes(x=x,y=y,color=group)) +
  geom_point()

