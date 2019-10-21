# A new take on visuzalizations for some tree objects in Whale. In particular, (1) branch colored
# species trees (with WGD nodes) and (2) reconciled gene trees.

# draw a reconciled tree
function drawtree(rtree::RecTree; width::Int64=400, height::Int64=300,
        fname::String="", nonretained::Bool=true, fontsize::Int64=7,
        linewidth=1)
    d = fname == "" ? Luxor.Drawing(width, height, :svg) : Luxor.Drawing(width, height, :svg, fname)
    Luxor.background("white")
    coords, paths = treecoords(rtree.tree, width=width, height=height)
    Luxor.sethue("black")
    Luxor.origin()
    Luxor.setline(linewidth)
    drawtree(coords, paths, width=width, height=height)
    rectreenodes(rtree, coords, nonretained=nonretained)
    leaflabels(rtree.leaves, coords, fontfamily="monospace", fontsize=fontsize)
    Luxor.finish()
    Luxor.preview()
end

function drawtree(rtree::ConRecTree; width::Int64=400, height::Int64=300,
        fname::String="", nonretained::Bool=true, fontsize::Int64=7,
        linewidth=1)
    d = fname == "" ? Luxor.Drawing(width, height, :svg) :
        Luxor.Drawing(width, height, :svg, fname)
    Luxor.background("white")
    coords, paths = treecoords(rtree.tree, width=width, height=height)
    Luxor.sethue("black")
    Luxor.origin()
    Luxor.setline(linewidth)
    drawtree(coords, paths, width=width, height=height)
    rectreenodes(rtree, coords, nonretained=nonretained)
    leaflabels(rtree.leaves, coords, fontfamily="monospace", fontsize=fontsize)
    supportvalues(rtree.rsupport, coords, fontfamily="monospace", fontsize=fontsize)
    Luxor.finish()
    Luxor.preview()
end

# draw a species tree
function draw_whaletree(stree::Arboreal; width::Int64=400, height::Int64=300,
        linewidth=1, fname::String="", fontsize::Int64=8, markersize=10, star=5,
        starcolor="white", margin=0.8)
    d = fname == "" ? Luxor.Drawing(width, height, :svg) :
        Luxor.Drawing(width, height, :svg, fname)
    coords, paths = treecoords(stree.tree, width=width, height=height, margin=margin)
    Luxor.background("white")
    Luxor.sethue("black")
    Luxor.origin()
    Luxor.setline(linewidth)
    drawtree(coords, paths, width=width, height=height)
    wgdnodes(stree, coords, markersize, starcolor, star)
    leaflabels(stree.leaves, coords, fontfamily="Lato italic", fontsize=fontsize)
    Luxor.finish()
    Luxor.preview()
end

"""
    coltree_whale(stree::SpeciesTree, values::Array{Float64};
        q::Array{Float64}=[], [kwargs])
Draw a species tree with inferred rates and retention rates. `values`
should be an array with a value at index `i` for branch `i`.
"""
function coltree_whale(tree::Arboreal, values::Array{Float64};
        q::Array{Float64}=[], width::Int64=400, height::Int64=300,
        fname::String="", fontsize=9, linewidth=3, margin=0.8,
        fontsizeq=8, boxwidth=8)
    d = fname == "" ? Luxor.Drawing(width, height, :svg) :
        Luxor.Drawing(width, height, :svg, fname)
    coords, paths = treecoords(tree.tree, width=width, height=height,
        margin=margin)
    Luxor.sethue("black")
    Luxor.origin()
    Luxor.setline(linewidth)
    data, nval = process_values(tree.rindex, values)
    coltree_whale(coords, paths, data, tree, width=width, height=height)
    length(q) > 0 ? wgdnodes(tree, coords, q, fontsize=fontsizeq,
        boxwidth=boxwidth) : wgdnodes(tree, coords)
    leaflabels(tree.leaves, coords, fontfamily="Lato italic", fontsize=fontsize)
    cbar = get_colorbar(width, height, 25)
    draw_colorbar(cbar, nval[1], nval[2])
    Luxor.finish()
    Luxor.preview()
end

"""
    coltree(stree::SpeciesTree, values::Array{Float64}; [kwargs])
"""
function coltree(tree::Arboreal, values::Array{Float64}; width::Int64=400,
        height::Int64=300, fname::String="")
    d = fname == "" ? Luxor.Drawing(width, height, :svg) :
        Luxor.Drawing(width, height, :svg, fname)
    coords, paths = treecoords(tree.tree, width=width, height=height)
    Luxor.sethue("black")
    Luxor.origin()
    Luxor.setline(3)
    data, nval = process_values(values)
    coltree(coords, paths, data, width=width, height=height)
    Luxor.sethue("black")
    leaflabels(tree.leaves, coords, fontfamily="Lato italic", fontsize=9)
    cbar = get_colorbar(width, height, 25)
    draw_colorbar(cbar, nval[1], nval[2])
    Luxor.finish()
    Luxor.preview()
end

# get the tree paths, general (for a PhyloTrees Tree), should just assign a
# coordinate to every node actually
function treecoords(tree::Tree; width::Int64=400, height::Int64=300, margin=0.8)
    leaves = findleaves(tree)
    root = findroots(tree)[1]
    coords = Dict{Int64,Luxor.Point}()
    Δx = width * margin / 2  # offset x-coordinate
    Δy = height / (length(leaves) + 1)  # vertical space between leaves
    yleaf = -Δy * (length(leaves) / 2)  # initial leaf y coordinate
    xmax = maximum([PhyloTrees.distance(tree, root, x) for x in leaves])  # maximum distance from root
    a = width * 0.7 / xmax  # modifier
    paths = []

    # postorder traversal
    function walk(node)
        if isleaf(tree, node)
            yleaf += Δy
            x = PhyloTrees.distance(tree, root, node) * a - Δx
            coords[node] = Luxor.Point(x, yleaf)
            return yleaf
        else
            ychildren = []
            for child in childnodes(tree,node)
                push!(ychildren, walk(child))
                push!(paths, (node, child))
            end
            y = sum(ychildren) / length(ychildren)
            x = PhyloTrees.distance(tree, root, node) * a - Δx
            coords[node] = Luxor.Point(x, y)
            return y
        end
    end
    walk(root)
    return coords, paths
end

# minimal tree drawing routine
function drawtree(coords, paths; width::Int64=400, height::Int64=300,
        rect::Bool=true)
    for p in paths
        rect ? drawhook(coords[p[1]], coords[p[2]]) : Luxor.line(coords[p[1]], coords[p[2]])
        Luxor.strokepath()
    end
end

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

# colored tree drawing routine
function coltree_whale(coords, paths, values, stree; width::Int64=400,
        height::Int64=300)
    for p in paths
        c2 = haskey(values, p[2]) ? get(ColorSchemes.viridis, values[p[2]]) :
            get(ColorSchemes.viridis, values[non_wgd_child(stree, p[2])])
        c1 = haskey(values, p[1]) ? get(ColorSchemes.viridis, values[p[1]]) : c2
        p1, p2 = coords[p[1]], coords[p[2]]
        p3 = Luxor.Point(p1.x, p2.y)
        Luxor.setblend(Luxor.blend(p1, p3, c1, c2))
        drawhook(p1, p2)
        Luxor.strokepath()
    end
end

# connect two points with a square corner
function drawhook(p1, p2)
    p3 = Luxor.Point(p1.x, p2.y)
    Luxor.poly([p1, p3, p2])
end

# add node markers for rectree nodes
function rectreenodes(rtree, coords::Dict; nonretained::Bool=true)
    for (node, coord) in coords
        if rtree.labels[node] == "duplication"
            Luxor.box(coord, 5, 5, 0, :fill)
        elseif rtree.labels[node] == "wgd"
            if !(nonretained)
                if length(childnodes(rtree.tree, node)) == 2
                    Luxor.star(coord, 4, 8, 0.5, 0, :fill)
                end
            else
                Luxor.star(coord, 4, 8, 0.5, 0, :fill)
            end
        elseif rtree.labels[node] == "loss"
            Luxor.star(coord, 1, 4, 5, 0, :fill)
        end
    end
end

# add node markers for wgd nodes
function wgdnodes(stree::Arboreal, coords::Dict, markersize=4, color="black",
        star=5)
    for (node, i) in stree.qindex
        Luxor.sethue(color)
        Luxor.star(coords[node], markersize, star, 0.5, 0, :fill)
        Luxor.sethue("black")
        Luxor.star(coords[node], markersize, star, 0.5, 0, :stroke)
    end
end

# add node markers for wgd nodes
function wgdnodes(stree::Arboreal, coords::Dict, q::Array{Float64};
        fontfamily="Lato", fontsize=8, thresh=0.1, boxwidth=8)
    Luxor.setfont(fontfamily, fontsize)
    Luxor.sethue("black")
    Luxor.setline(1)
    for (node, i) in stree.qindex
        q[i] > thresh ? Luxor.box(coords[node], boxwidth, boxwidth, 0, :fill) :
            Luxor.box(coords[node], boxwidth, boxwidth, 0, :stroke)
        i % 2 == 1 ? Δy = -10 : Δy = 10
        tcoord = Luxor.Point(coords[node].x, coords[node].y + Δy)
        Luxor.settext(string(round(q[i], digits=2)), tcoord,
            valign="center", halign="center")
    end
end

# add node labels
function labelnodes(coords::Dict; fontfamily="monospace", fontsize=9)
    Luxor.setfont(fontfamily, fontsize)
    for (node, coord) in coords
        Luxor.box(coords[node], 5, 5, 0, :fill)
        tcoord = Luxor.Point(coords[node].x, coords[node].y + 10)
        Luxor.settext(" " * string(node), tcoord,
            valign="center", halign="left")
    end
end

# add leaf labels
function leaflabels(leaves::Dict, coords::Dict; fontfamily="monospace",
        fontsize=6)
    Luxor.setfont(fontfamily, fontsize)
    for (node, lab) in leaves
        Luxor.settext("  " * lab, coords[node], valign="center", halign="left")
    end
end

# add leaf labels
function supportvalues(support::Dict, coords::Dict; fontfamily="monospace",
        fontsize=6)
    Luxor.setfont(fontfamily, fontsize)
    for (node, lab) in support
        l = lab < 1. ? (@sprintf "  %.2f" lab) : ""
        Luxor.settext(l, coords[node], valign="center", halign="left")
    end
end

# preprocessing for colortree
function process_values(rindex, values::Array{Float64})
    nval = (minimum(values), maximum(values))
    nvalues = (values .- nval[1]) ./ (nval[2] - nval[1])
    data = Dict(k => nvalues[i] for (k,i) in rindex)
    return data, nval
end

# get a colorbar
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
        Luxor.setline(8)
        ca = get(ColorSchemes.viridis, p[3])
        cb = get(ColorSchemes.viridis, p[4])
        bl = Luxor.blend(p[1], p[2], ca, cb)
        Luxor.line(p[1], p[2])
        Luxor.setblend(bl)
        Luxor.strokepath()
    end
    Luxor.setcolor("black")
    Luxor.setfont("Lato", 8)
    s1 = @sprintf " %3.2f" minv
    s2 = @sprintf " %3.2f" maxv
    p1 = Luxor.Point(cbar[1][1].x + pad, cbar[1][1].y)
    p2 = Luxor.Point(cbar[end][2].x + pad, cbar[end][2].y)
    Luxor.settext(s1, p1, valign="center", halign="left")
    Luxor.settext(s2, p2, valign="center", halign="left")
end

function non_wgd_child(tree, n)
    while outdegree(tree.tree, n) == 1
        n = childnodes(tree.tree, n)[1]
    end
    return n
end
