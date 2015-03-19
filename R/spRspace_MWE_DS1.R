### This script is meant to be a first stab at using the sprspace pipeline (https://github.com/cwhidden/sprspace),
###  also see accompanying paper http://sysbio.oxfordjournals.org/content/early/2015/01/27/sysbio.syv006.abstract)
###  We will analyse a small posterior sample of 10,000 trees of the DS1 data set in Whidden & Matsen (2015) comprising 27 taxa.
library(ape)
source("read.nexusB.r")
source("unique_trees.r")
source("sort_trees.r")
source_url("https://raw.githubusercontent.com/maxbiostat/rspr/master/R/rspR.r")
source("spRspace_aux.R")
##
mcc.tree <- read.nexus("DATA/DS1_run_1.tree")
trees <- read.nexusB("DATA/DS1_run_1.trees")
###########
### Distances and all that
###########
seqdist <- sapply(seq_along(trees)[-1], function(p) distt2(pos = p, x = trees) ) # "sequential" distance, i.e, how different is a STATE from the previously kept STATE?
plot(seqdist, type = "l")

withmcc <- vector(length(trees) + 1, mode = "list")
withmcc[[1]] <- mcc.tree
for (i in 1:length(trees))  withmcc[[i+1]] <- trees[[i]] 
class(withmcc) <- "multiPhylo"
tomcc <- as.numeric(rspr.matrix(trees = withmcc, type = "o"))

plot(tomcc, type = "l")
png("RESULTS/DS1_spr_distances_to_mcc.png")
barplot(table(tomcc), main = "SPR distances to MCC tree")
dev.off()
############
#### Manipulating the tree posterior
utrees <- unique_trees(trees)# unique topologies
save(utrees, file = "DATA/dengue/DS1_10K_unique_trees.RData")
length(utrees[[1]]) # How many [U]nique trees?
uhits <- utrees[[2]] # How many times each?  
most.sampled <- which(uhits == max(uhits)) # most frequent topology. A rare commodity, though.
rspr(utrees$trees[[most.sampled]], mcc.tree)
summary(tomcc)
rutrees <- sort_trees(utrees) # sorting trees according to their frequency
## Now let's follow Whidden & Matsen (2015) and use the 'first' first trees in the list of [R]anked [U]nique trees
first <- min(length(rutrees), 256)
pp <- (uhits/sum(uhits))[1:first]
first.rutrees <- rutrees[1:first]
system.time(
  matrix.full <- rspr.matrix(first.rutrees, type = "full")
)
system.time(
  matrix.radius.one <- rspr.matrix(first.rutrees, type = "restricted", maxdist = 1)
)
write.csv(matrix.full, file = paste("RESULTS/DS1_spr_distances_credibleset_first_", first, ".csv", sep = ""), row.names = FALSE, col.names = FALSE)

# Now let's compute the distance of each tree in the ordered (ranked) posterior credible set to the mcc tree
plusmcc <- vector(first + 1, mode = "list")
plusmcc[[1]] <- mcc.tree
for (i in 1:length(first.rutrees))  plusmcc[[i+1]] <- first.rutrees [[i]] 
class(plusmcc) <- "multiPhylo"

ocs2mcc <- as.numeric(rspr.matrix(trees = plusmcc, type = "o"))[-1] # the first entry is just the MCC tree and should always be zero. OCS = ordered credible set.
plot(ocs2mcc, type = "l", ylab = "SPR distance", main = "Distance to MCC in the ordered posterior credible set")
write.table(data.frame(distance = matrix(ocs2mcc, ncol = 1)),
            file = paste("RESULTS/DS1_", first, "_dists2mcc_ordered_credible_set.txt", sep = ""),
            row.names = FALSE)

#####################
#### Graph analysis
#####################
library(igraph)
pp.col.pos <- match(pp, sort(unique(pp)))
dist.col.pos <- match(ocs2mcc, sort(unique(ocs2mcc)))
k.pp <-  length(unique(pp.col.pos))
k.dist <- length(unique(dist.col.pos))

### Radius graph
radius.graph <- graph.adjacency(adjmatrix = binarise.matrix(m = matrix.radius.one), mode = "undirected")
write.graph(radius.graph, file = paste("RESULTS/DS1_first_", first, "_graph_radius_1.gml", sep = ""), format = "gml")
degdist.ro <- degree.distribution(graph = radius.graph)
deg.ro <- degree(radius.graph)
( mean.degree.ro <- mean(deg.ro)) 
( g.diameter.ro <- igraph::diameter(radius.graph) )
modularity(x = radius.graph, membership = dist.col.pos)
modularity(x = radius.graph, membership = pp.col.pos)
assortativity.degree(graph  = radius.graph, directed = FALSE)
assortativity.nominal(radius.graph, types = as.integer(pp.col.pos), directed = FALSE)
assortativity.nominal(radius.graph, types = as.integer(dist.col.pos), directed = FALSE)

############
#### Plotting
############

# library(RColorBrewer)
# brewer.pal(k, "YlOrRd") ## Alternative pallete
library(fields)
## Heat map of the distances in the credible set
png(filename = paste("RESULTS/DS1_first_", first, "_posterior.png", sep = ""))
N <- ncol(matrix.full)
labs <- rep("", N)
par(mar = c(5.1, 2.1, 4.1, 4.1))
image(matrix.full, axes = FALSE, col = 'transparent')
axis(1, at = 1:N , labels = labs)
axis(2, at = 1:N , labels = labs)
image.plot(matrix.full, add = TRUE, legend.mar = 3.1)
dev.off()
# Now the graph. Vertex sizes are pp (normalised uhits) and colours are distance to MCC
png(paste("RESULTS/DS1_graph_first_", first, "_MCCdist.png", sep = ""))
plot(radius.graph, vertex.size = pp.col.pos, vertex.color = heat.colors(k.dist)[rev(dist.col.pos)],
     vertex.label = NA, layout = layout.fruchterman.reingold)
dev.off()
# Boxplot of degree versus distance to MCC
png("RESULTS/DS1_boxplots_degree_versus_distancetoMCC.png")
boxplot(degree(radius.graph)~ocs2mcc, col = rev(heat.colors(k.dist)), main = "Vertex (tree) degree versus distance to MCC tree",
        xlab = "SPR distance to MCC", ylab = "k")
dev.off()
