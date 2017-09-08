dummy.input.vector <- c(1.1, 0.25, 0, 3, 0.23, 0.23,  # what if 1.1 becomes 1.4
                        45, 45, #45, 45,              # what if 45 becomes 60
                        -0.5, 2.8, -0.2, -0.2, -2.5, -0.52, -0.05)# c(1000, 2, 3, 4)

source("/Users/delvaw/Documents/MiceABC/R/00-Functions.R")
seed <- 1
inputvector <- c(seed, dummy.input.vector)

library('igraph')
# We start with a transnetwork object (the 7th chain of the 10 chains)
transnetwork <- transm.ls[[7]] # 1 seed considered at a time
transnetw.attr <- transm.ls.attrib[[7]]


edges.df <- data.frame(from = transnetwork$parent,
                    to = transnetwork$id,
                    inf.time = 10 - transnetwork$itimes)
vertices.df <- data.frame(nodes = unique(c(transnetwork$parent, transnetwork$id)))

net <- graph_from_data_frame(d = edges.df, vertices = vertices.df, directed=T)
V(net)$size <- 3

#V(net)$frame.color <- "white"
#V(net)$color <- "orange"
#V(net)$label <- ""
E(net)$arrow.size <- .5


l <- layout_with_kk(net) # layout_with_lgl layout_with_kk
plot(net, layout = l, vertex.size = 3, vertex.label = NA)

epi.tree <- epi2tree2(transnetwork)


tree.dat.full <- simSeq(epi.tree,
                        l = 1000,
                        bf = freq,
                        #rootseq = hiv.seq.env.short,
                        type = "DNA",
                        rate = overall.rate)
fit.ini <- pml(epi.tree, tree.dat.full, k = 4)

fit <- optim.pml(fit.ini, optNni = FALSE, optBf = FALSE, optQ = TRUE, model = "GTR", optGamma = TRUE, optEdge = FALSE, optRate = FALSE, optRooted = FALSE)
# fit$tree stores the tree
plot(epi.tree)
plot(fit$tree)
all.equal(epi.tree, fit$tree) # We didn't change the tree, only the model for the viral evolution, and the simulated sequences.
tree.dat.full.GTR = simSeq(fit) # The phyDat object that stores the sequences
plot(fit$tree, root.edge = TRUE, type = "fan", open.angle = 20, rotate.tree = 10) #fan
add.scale.bar(x = 0, y = 0, length = 10, lwd = 2, lcol = "blue3", col = "blue3")
# This tree is the "hyperfull" tree: it's topology is given entirely by the transmission events and
# it does not check whether everybody is still alive at the time of the sequencing.

# Let's therefore now prune the tree, and then reconstruct it.
alive.selected.vector <- transnetw.attr$id.orig.recip %in% ID.vect.alive
sequences.keep <- transnetw.attr$id[alive.selected.vector] # Only the sequences of these tip labels (i.e. new IDs) can be used to build the phylo tree
numbered.tips <- 1:length(alive.selected.vector)

numbered.keep <- numbered.tips[as.numeric(names(tree.dat.full.GTR)) %in% sequences.keep]
#
tree.dat.pruned <- subset(tree.dat.full.GTR,
                          subset = numbered.keep)# Only keep sequences of people that were alive at datalist.agemix$itable$population.simtime[1]
  tree.ml.pruned <- dist.ml(tree.dat.pruned, model = "F81", bf = freq) # F81

  tree.sim.pruned <- nj(dist.dna(as.DNAbin(tree.dat.pruned), model = "TN93")) #nj(tree.ml.pruned) # phangorn::NJ(tree.ml.pruned) # changed upgma to neighbour-joining

  internal.node.times.element <- branching.times(tree.sim.pruned) # datalist.agemix$itable$hivseed.time[1] + branching.times(tree.sim.full) * (datalist.agemix$itable$population.simtime[1] - datalist.agemix$itable$hivseed.time[1])
  internal.node.times.unsorted <- c(internal.node.times, internal.node.times.element)
  internal.node.times.sorted <- sort(as.numeric(internal.node.times.unsorted))
  recent.transmissions <- sum(internal.node.times.sorted < 1 * overall.rate)


# tree4treedater <- root(fit$tree, outgroup = 1, resolve.root = TRUE)
# plot(tree4treedater)
# head(tree4treedater$tip.label)
# # Now the dates of the sequences in the tree
# tree4treedater.dates <- annot[, 2]
# names(tree4treedater.dates) <- annot[, 1]
# head(tree4treedater.dates)
#
# #help(dater)
# dated.tree <- dater(tre = fit$tree,
#                     sts = tree4treedater.dates,
#                     s = attr(dna4, "nr"))
# plot(dated.tree)
# axisPhylo(root.time = min(tree4treedater.dates),
#           backward = FALSE)
#
#
#
# data("H3N2")
# myPath <- system.file("files/usflu.fasta",package="adegenet")
# myPath2 <- system.file("data/usflu.annot.csv",package="adegenet")
#
# dna.test <- read.dna(file = myPath, format = "fasta") # See tutorial: http://adegenet.r-forge.r-project.org/files/MRC-session2-tuto.1.3.pdf
# annot <- read.csv(file = "https://raw.githubusercontent.com/reconhub/phylo-practical/master/data/usflu.annot.csv", header = TRUE, row.names = 1)
#
# D <- dist.dna(dna.test, model = "TN93")
# class(D)
# tre <- nj(D)
# plot(tre, cex = 0.6)
# tre2 <- root(tre, out = 1)
# dna2 <- as.phyDat(dna.test)
# # Initial fit
# pml(tre2, dna2, k = 4)
#
# na.posi <- which(apply(as.character(dna.test), 2, function(e) any(!e %in% c("a", "t", "g", "c"))))
# dna3 <- dna.test[, -na.posi]
# dna4 <- as.phyDat(dna3)
# tre.ini <- nj(dist.dna(dna3, model = "TN93"))
# tre.ini <- root(tre.ini, out = 1, resolve.root = TRUE)
#
# fit.ini <- pml(tre.ini, dna4, k = 4)
# fit.ini
#
# fit <- optim.pml(fit.ini, optNni = TRUE, optBf = TRUE, optQ = TRUE, optGamma = TRUE)
# anova(fit.ini, fit)
# BIC(fit.ini)
# BIC(fit)
#
# tree4treedater <- root(fit$tree, outgroup = 1, resolve.root = TRUE)
# plot(tree4treedater)
# head(tree4treedater$tip.label)
# # Now the dates of the sequences in the tree
# tree4treedater.dates <- annot[, 2]
# names(tree4treedater.dates) <- annot[, 1]
# head(tree4treedater.dates)
#
# #help(dater)
# dated.tree <- dater(tre = tree4treedater,
#                     sts = tree4treedater.dates,
#                     s = attr(dna4, "nr"))
# plot(dated.tree)
# axisPhylo(root.time = min(tree4treedater.dates),
#           backward = FALSE)
