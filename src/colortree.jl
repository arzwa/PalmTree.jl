
"""
Colortree from a posterior distribution represented by a data frame. Columns
are assumed to be named as `parameter+index`, e.g. :λ1, ... , :λ17. If the
tree structure contains a `rindex` of `bindex`; that is used to related branches
to parameters; else the branch index is used.
"""
function colortree(tree::Arboreal, df::DataFrame, θ::Symbol;
        f::Function=mean, scale::Function=identity)
    cols = [x for x in names(df) if startswith(string(x), string(θ))]
    values = Dict(x=>f(df[x]) for x in cols)
end

# low-level function
function colortree(tree::Arboreal, x::Vector{Number})
    coords, paths = treecoords(tree.tree, width=width, height=height,
        margin=margin)

# colored tree drawing routine
function coltree(coords, paths, values; width::Int64=400,
        height::Int64=300)
    for p in paths
        c2 = get(ColorSchemes.viridis, values[p[2]])
        c1 = get(ColorSchemes.viridis, values[p[1]])
        p1, p2 = coords[p[1]], coords[p[2]]
        p3 = Luxor.Point(p1.x, p2.y)
        Luxor.setblend(Luxor.blend(p1, p3, c1, c2))
        drawhook(p1, p2)
        Luxor.strokepath()
    end
end
