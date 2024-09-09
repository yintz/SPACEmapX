



Dendrogram <- function(location){

clustering <- read.dendrogram(file=location)
clustering_phylo <- as.phylo(clustering)
# Visualizing Tree
#Next, we use the dendrogram file to visualize the dendrogram itself.
my.subtrees = subtrees(clustering_phylo)  # subtrees() to subset
png("Dendrogram.png",width=10000,height=10500, res = 300)
plot(clustering_phylo,show.tip.label = FALSE)
nodelabels(text=1:clustering_phylo$Nnode,node=1:clustering_phylo$Nnode+Ntip(clustering_phylo))
dev.off()
}
