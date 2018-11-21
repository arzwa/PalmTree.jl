__precompile__()

module PalmTree

    using Luxor
    using PhyloTrees
    using Random

    export
        LabeledTree,
        randomtree,
        drawtree

    include("types.jl")
    include("utils.jl")
    include("minimal.jl")

end # module
