
### Simulating mtDNA gene trees
### Edward A. Myers; 07.19.21

###############
install.packages("rlang", dependencies = c("Depends", "Imports", "LinkingTo", "Suggests"))
install.packages("devtools")
library(devtools)
install_github("lliu1871/phybase")
library(phybase)
library(phytools)
library(beastio)
##############

#Below I have re-scaled branch lengths by fixing the root node at 37.3 mya.
ano <- "((((garmani:12.82,grahami:12.82):2.901,(opalinus:10.8953,opalinus_Blue_Mts:10.8953):4.8346):6.3691,(lineatopus:13.90776,reconditus:13.90776):8.1913):15.2009,valencienni:37.3);"

plottree(ano)
sp_name <- species.name(ano)
node_matrix <- read.tree.nodes(str = ano, name=sp_name)$nodes
node_matrix[,5]<-2

#Make sure it works before running loop
sim1 <- sim.coaltree.sp(rootnode=13, nodematrix=node_matrix,nspecies=7,seq=as.double(c("35", "20", "27", "15", "7", "4", "20")), name=sp_name)$gt
plottree(sim1)

# From Brian O'Meara (https://rdrr.io/github/bomeara/phybase/man/phybase2phylo.html)
phybase2phylo <- function(x) {
  if(grepl('#', x)) {
    x<-gsub("#\\.*\\d*\\.*\\d*","",x)
  }
  return(read.tree(text=x))
}

for(i in 1:10000){
  write.tree(phybase2phylo(sim.coaltree.sp(rootnode=13, nodematrix=node_matrix,nspecies=7,seq=as.double(c("35", "20", "27", "15", "7", "4", "20")), name=sp_name)$gt), "Anolis_sims.tre", append = T)
}
read.tree("Anolis_sims.tre", keep.multi = T)->Z
length(Z)
opalinus <- c("opalinuss1", "opalinuss2", "opalinuss3", "opalinuss4", "opalinuss5", "opalinuss6", "opalinuss7", "opalinuss8", "opalinuss9", "opalinuss10", "opalinuss11", "opalinuss12", "opalinuss13", "opalinuss14", "opalinuss15", "opalinus_Blue_Mtss1", "opalinus_Blue_Mtss2", "opalinus_Blue_Mtss3", "opalinus_Blue_Mtss4", "opalinus_Blue_Mtss5", "opalinus_Blue_Mtss6", "opalinus_Blue_Mtss7")
getCladeMonophyly(Z, tips= opalinus)
val_bl <- c("valenciennis1", "valenciennis2", "valenciennis3", "valenciennis4", "valenciennis5", "valenciennis6", "valenciennis7", "valenciennis8", "valenciennis9", "valenciennis10", "valenciennis11", "valenciennis12", "valenciennis13", "valenciennis14", "valenciennis15", "valenciennis16", "valenciennis17", "valenciennis18", "valenciennis19", "valenciennis20", "opalinus_Blue_Mtss1", "opalinus_Blue_Mtss2", "opalinus_Blue_Mtss3", "opalinus_Blue_Mtss4", "opalinus_Blue_Mtss5", "opalinus_Blue_Mtss6", "opalinus_Blue_Mtss7")
getCladeMonophyly(Z, tips= val_bl)
