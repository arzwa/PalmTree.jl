# # PalmTree.jl

# Plot trees using Luxor. Works for abstract tree data structures consisting of
# nodes that implement a `children`, `distance` and `id` function. Allows flexible
# use of the Luxor library.

using NewickTree, PalmTree, Luxor

outdir  = mkpath(joinpath(@__DIR__, "assets/"))
outpath = joinpath(outdir, "tree.svg")
tree = readnw("((A:77.3,B:77.3):16.4,(C:46.9,(D:23.4,E:23.4):23.4):46.9);")

# Draw the tree
tl = TreeLayout(tree, dims=(280, 250))
@svg begin
    origin(Point(10,10))
    setline(4)
    drawtree(tl, color=n->Luxor.RGB(repeat(rand(1), 3)...))
    setfont("Noto sans italic", 16)
    sethue(Luxor.RGB(0,0,0))
    nodemap(tl, (k, p)->settext("  "*name(k), p, valign="center"))
    nodemap(tl, (k, p)->rand() < 0.5 ? box(p, 10, 10, :fill) : star(p, 8, 5, 0.5, 4, :fill))
end 350 300 outpath;

# ![](docs/assets/tree.svg)

using Literate  #src
Literate.markdown(joinpath(@__DIR__, "README.jl"), joinpath(@__DIR__, "../"), documenter=false, execute=true)  #src
