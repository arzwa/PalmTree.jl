
module PalmTree

    using Luxor
    using PhyloTrees
    using Printf
    using ColorSchemes
    using Parameters
    import ColorSchemes: RGB

    export drawtree, TreeLayout #, coltree, draw_whaletree

    include("main.jl")
    include("core.jl")

end # module
