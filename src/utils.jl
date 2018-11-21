"""
    randomtree(nleaves::Int)
"""
function randomtree(n::Int64)
    tree = Tree()
    addnode!(tree)
    leaves = Int64[1]
    while length(leaves) < n
        shuffle!(leaves)
        mom = pop!(leaves)
        addnode!(tree); kid = length(tree.nodes)
        addbranch!(tree, mom, kid, abs(randn()))
        push!(leaves, kid)
        addnode!(tree); kid = length(tree.nodes)
        addbranch!(tree, mom, kid, abs(randn()))
        push!(leaves, kid)
    end
    return tree
end
