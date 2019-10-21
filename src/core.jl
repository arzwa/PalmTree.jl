
mutable struct TreeLayout
    coords::Dict{Int64,Luxor.Point}     # coordinates of nodes
    colors::Dict{Int64,RGB}             # colors associated with nodes
    range::Tuple{Float64,Float64}       # values corresponding to color extremes
    paths::Array{Tuple{Int64,Int64},1}  # edge list
    bbox::Tuple                         # current calculated bounding box
end

function TreeLayout(coords::Dict, paths::Array, bbox::Tuple)
    colors = Dict(k=>RGB(0.,0.,0.) for (k,v) in coords)
    range = (0., 0.)
    TreeLayout(coords, colors, range, paths, bbox)
end

TreeLayout(tree::T, width, height) where T<:Arboreal =
    treecoords(tree.tree, width, height)

function TreeLayout(tree::T, width, height, nodevalues::Dict) where T<:Arboreal
    tl = treecoords(tree.tree, width, height)
    setcolors!(tl, nodevalues)
    tl
end

function TreeLayout(tree::T, width, height, nodevalues::Dict, mn, mx) where
        T<:Arboreal
    tl = treecoords(tree.tree, width, height)
    setcolors!(tl, nodevalues, mn, mx)
    tl
end

function setcolors!(tl::TreeLayout, nodevalues::Dict, mn, mx;
        cs=ColorSchemes.viridis, tf=identity)
    for (k, v) in nodevalues
        tl.colors[k] = get(cs, (v - mn)/(mx - mn))
    end
    tl.range = (mn, mx)
end

function treecoords(tree, width, height)
    leaves = findleaves(tree)
    root = findroots(tree)[1]
    coords = Dict{Int64,Luxor.Point}()
    yleaf = -1.
    xmax = 0.
    paths = []
    # postorder traversal
    function walk(node)
        if isleaf(tree, node)
            yleaf += 1.
            x = PhyloTrees.distance(tree, root, node)
            xmax = x > xmax ? x : xmax
            coords[node] = Luxor.Point(x, yleaf)
            return yleaf
        else
            ychildren = []
            for child in childnodes(tree,node)
                push!(ychildren, walk(child))
                push!(paths, (node, child))
            end
            y = sum(ychildren) / length(ychildren)
            x = PhyloTrees.distance(tree, root, node)
            coords[node] = Luxor.Point(x, y)
            return y
        end
    end
    walk(root)
    scalex, scaley = scale!(coords, xmax, yleaf, width, height)
    return TreeLayout(coords, paths, (xmax*scalex, yleaf*scaley))
end

function scale!(coords, xmx, ymx, w, h)
    scalex = w/xmx
    scaley = h/ymx
    for (k,v) in coords
        coords[k] = Point(v.x*scalex, v.y*scaley)
    end
    scalex, scaley
end

function drawtree(tl::TreeLayout)
    # quite general, for colored and normal trees
    # linewidth etc. are assumed to be set
    @unpack coords, paths, colors = tl
    for p in paths
        p1, p2 = coords[p[1]], coords[p[2]]
        p3 = Luxor.Point(p1.x, p2.y)
        setblend(blend(p1, p3, colors[p[1]], colors[p[2]]))
        poly([p1, p3, p2])
        Luxor.strokepath()
    end
end

function drawwgds(tl, tree, q; threshold=0.1)
    for (k,v) in tree.qindex
        setline(2)
        if q[k] > threshold
            sethue("black")
        else
            sethue("white")
        end
        Luxor.box(tl.coords[k], 12, 12, 0, :fill)
        sethue("black")
        Luxor.box(tl.coords[k], 12, 12, 1, :stroke)
    end
end

function leaflabels(tl::TreeLayout, labels)
    for (k, v) in labels
        p = tl.coords[k]
        q = Luxor.Point(0., p.y)
        settext(string(v), q, valign="center")
    end
end

function drawtree(tree::Arboreal;
        width=600, height=400, margin=40, fontsize=20, lw=3,
        font="Lato Italic", bg="random", leafcolumn=0.3)
    Drawing(width, height, :svg)
    sethue("black")
    setline(lw)
    t = Table([height], [width*(1. - leafcolumn) , width*leafcolumn])
    w = t.colwidths[1]-margin
    h = t.rowheights[1]-margin
    tl = treecoords(tree.tree, w, h)

    # tree
    bg == "random" ? randomhue() : sethue(bg)
    setopacity(0.2)
    box(BoundingBox(centered=false), :fill)
    translate(margin/2, margin/2)
    drawtree(tl)

    # leaf labels
    origin(t[2])
    translate(t[2].x, margin/2)
    setopacity(1)
    sethue("black")
    setfont(font, fontsize)
    leaflabels(tl, tree.leaves)

    finish()
    preview()
end

function colortree(tree::Arboreal, vals::NamedTuple;
        width=600, height=400, margin=40, fontsize=20, lw=6,
        font="Lato Italic", bg="random", leafcolumn=0.3, fname="")
    fname == "" ? Drawing(width, height, :svg) :
        Drawing(width, height, :svg, fname)
    sethue("black")
    setline(lw)
    t = Table([height], [width*(1. - leafcolumn) , width*leafcolumn])
    w = t.colwidths[1]-margin
    h = t.rowheights[1]-margin

    # tree
    x = [collect(values(vals.λ)) ; collect(values(vals.μ))]
    mn, mx = 0., 0.7

    tl = TreeLayout(tree, w, h, vals.λ, mn, mx)
    bg == "random" ? randomhue() : sethue(bg)
    setopacity(0.2)
    box(BoundingBox(centered=false), :fill)
    translate(margin/2, margin/2)
    drawtree(tl)

    tl = TreeLayout(tree, w, h, vals.μ, mn, mx)
    bg == "random" ? randomhue() : sethue(bg)
    setopacity(0.2)
    box(BoundingBox(centered=false), :fill)
    #translate(margin/2, margin/2)
    translate(lw+1,lw+1)
    drawtree(tl)

    setopacity(1)
    sethue("black")
    PalmTree.drawwgds(tl, tree, vals.q)

    # leaf labels
    origin(t[2])
    translate(t[2].x, margin/2)
    setopacity(1)
    sethue("black")
    setfont(font, fontsize)
    leaflabels(tl, tree.leaves)

    # colorbar
    translate(w*0.7, h*0.65)
    cb = get_colorbar(w, h, 50)
    draw_colorbar(cb, mn, mx)

    finish()
    preview()
end

function get_colorbar(width, height, pad; breaks=7)
    path = []
    c = 0.
    x = -width/2 + pad
    y = height/2 - pad
    y2 = height/2 - 3*pad
    Δy = (y2 - y) / breaks
    Δc = 1/breaks
    for i = 1:breaks
        push!(path, (Luxor.Point(x, y), Luxor.Point(x, y+Δy), c, c+Δc))
        y += Δy
        c += Δc
    end
    return path
end

# draw a colorbar
function draw_colorbar(cbar, minv, maxv; pad=2)
    for p in cbar
        Luxor.setline(12)
        ca = get(ColorSchemes.viridis, p[3])
        cb = get(ColorSchemes.viridis, p[4])
        bl = Luxor.blend(p[1], p[2], ca, cb)
        Luxor.line(p[1], p[2])
        Luxor.setblend(bl)
        Luxor.strokepath()
    end
    Luxor.setcolor("black")
    Luxor.setfont("Lato", 16)
    s1 = @sprintf " %3.2f" minv
    s2 = @sprintf " %3.2f" maxv
    p1 = Luxor.Point(cbar[1][1].x + pad, cbar[1][1].y)
    p2 = Luxor.Point(cbar[end][2].x + pad, cbar[end][2].y)
    Luxor.settext(s1, p1, valign="center", halign="left")
    Luxor.settext(s2, p2, valign="center", halign="left")
end
