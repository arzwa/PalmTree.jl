# PalmTree.jl

Visualization tools for phylogenetic trees in Julia as defined in my fork of PhyloTrees.jl. Figures are constructed in vector graphics using Luxor.jl. A bit of a mess currently.

Example

```julia
t = readtree("test/morris-9taxa.nw")
drawtree(t)
coltree(t, rand(length(t.tree.nodes)))
```
