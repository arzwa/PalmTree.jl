using PalmTree
using PhyloTrees

t = readtree("test/morris-9taxa.nw")
drawtree(t)
coltree(t, rand(length(t.tree.nodes)))
