abstract type Arboreous end

struct LabeledTree <: Arboreous
    tree::Tree
    leaflabels::Dict{Int64,AbstractString}

    LabeledTree(t) = new(t, Dict(x => string(x) for x in keys(t.nodes)))
    LabeledTree(t, leaflabels) = new(t, leaflabels)
end
