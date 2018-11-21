# Minimal tree plotting (just black & white + leaf labels)

drawtree(tree::Tree) = drawtree(LabeledTree(tree))

function drawtree(tree::Arboreous; width::Int64=400, height::Int64=300, fname="")
    d = get_treepaths(tree)
    drawdict(d, width, height)
end

function get_treepaths(arboreal::Arboreous; width::Int64=400, height::Int64=300)
    # prepare
    tree = arboreal.tree
    pths = [] ; labs = []
    leaves = findleaves(tree)

    # (initial) settings
    Δy = height / (length(leaves) + 2)
    Δx = width * 0.8 / 2
    leaf_y = -Δy * (length(leaves) / 2)  # initial leaf y coordinate
    max_x = maximum([PhyloTrees.distance(tree, 1, x) for x  in leaves])
    a = width * 0.7 / max_x

    function walk(node)
        if isleaf(tree, node)
            leaf_y += Δy
            x = PhyloTrees.distance(tree, 1, node) * a - Δx
            push!(labs, [Luxor.Point(x, leaf_y), arboreal.leaflabels[node]])
            return Luxor.Point(x, leaf_y)
        else
            child = Luxor.Point[]
            children = childnodes(tree, node)
            for c in children
                child_point = walk(c)
                push!(child, child_point)
            end

            # construct path
            x = PhyloTrees.distance(tree, 1, node) * a - Δx
            y = sum([p.y for p in child])/length(child) # midpoint
            for i = 1:length(children)
                push!(pths, (Luxor.Point(x, y), child[i]))
            end
            return Luxor.Point(x, y)
        end
    end
    walk(1)
    return Dict("pths" => pths, "labs" => labs)
end


function drawdict(dict, width::Int64, height::Int64; fname::String="")
    fname == "" ? d = Luxor.Drawing(width, height, :svg) :
        d = Luxor.Drawing(width, height, :svg, fname)
    Luxor.sethue("black")
    Luxor.origin()
    Luxor.setline(2)
    Luxor.setfont("Lato Italic", 8)
    for p in dict["pths"]
        corner = Luxor.Point(p[1].x, p[2].y)
        Luxor.line(corner, p[2]); Luxor.line(p[1], corner)
        Luxor.strokepath()
    end
    for l in dict["labs"]
        Luxor.settext("  " * l[2], l[1], valign="center", halign="left")
    end
    Luxor.finish()
    Luxor.preview()
end
