# PalmTree.jl

Plot trees using Luxor. Works for abstract tree data structures consisting
of nodes that implement a `children`, `distance` and `id` function. Allows
flexible use of the Luxor library.

```julia
using NewickTree, PalmTree, Luxor

tr = readnw("((A:77.3,B:77.3):16.4,(C:46.9,(D:23.4,E:23.4):23.4):46.9);")
tl = TreeLayout(tr[1])
@svg begin
    origin(Point(10,10))
    drawtree(tl)
    setfont("Noto mono", 16)
    nodemap(tl, (k, p)->
        settext(" "*tr.leaves[k], p, valign="center"), keys(tr.leaves))
    nodemap(tl, (k, p)->box(p - (5,5), p + (5,5), :fill))
end
```
