using PalmTree
using DataFrames, PhyloTrees, Whale, CSV, Distributions
base = "/home/arzwa/arthropods/"

# helper functions =============================================================
function leafrename!(tree, d, abbr=true)
    for (k, v) in tree.leaves
        if abbr
            s = split(d[v])
            tree.leaves[k] = s[1][1] * ". " * s[2]
        else
            tree.leaves[k] = d[v]
        end
    end
    tree
end

function posterior_means(df, tree)
    λ = [mean(df[c]) for c in names(df) if startswith(string(c), "λ")]
    μ = [mean(df[c]) for c in names(df) if startswith(string(c), "μ")]
    q = [mean(df[c]) for c in names(df) if startswith(string(c), "q")]
    λ = Dict(n => λ[i] for (n,i) in tree.rindex)
    μ = Dict(n => μ[i] for (n,i) in tree.rindex)
    q = Dict(n => q[i] for (n,i) in tree.qindex)
    return (λ=λ, μ=μ, q=q)
end

function colortrees(df, st, prefix)
    λ, μ, q = posterior_means(df)
    PalmTree.coltree_whale(st, λ, q=q, width=330, fname="$prefix-l.svg")
    PalmTree.coltree_whale(st, μ, q=q, width=330, fname="$prefix-m.svg")
end

# tree conf ====================================================================
conf1 = Dict(
    "cfe" => ("cfe"    , 1.21),
    "lep" => ("bmo,pra", 1.73),
    "atu" => ("atu"    , 0.83),
    "col" => ("atu,tca", 2.47),
    "hym" => ("ame,aro", 2.92))

conf2 = Dict(
    "foc" => ("foc"    , 1.70),
    "pol" => ("zne,bge", 1.95),
    "ins" => ("api,zne", 4.24),
    "fca" => ("fca"    , 0.99),
    "col" => ("fca,hdu", 3.38),
    "lp1" => ("lpo"    , 1.94),
    "lp2" => ("lpo"    , 3.89))

species = open("$base/nw/species.txt", "r") do f
    Dict(string(split(line, ",")[1]) => string(split(line, ",")[2])
        for line in readlines(f))
end

st2 = leafrename!(SlicedTree("$base/nw/species2.nw", conf1), species)
st3 = leafrename!(SlicedTree("$base/nw/species3.nw", conf2), species)
set_equalrootrates!(st3)

# ==============================================================================
df = CSV.read("$base/trace/ortho2-all-iid.trace")
vals = posterior_means(df, st2)

PalmTree.colortree(st2, vals, lw=6; bg="white",
    fname="../arthropods/img/doubletree1.svg")

df = CSV.read("$base/trace/ortho3-all-eqriid-eta07.trace")
vals = posterior_means(df, st3)

PalmTree.colortree(st3, vals, lw=6; bg="white",
    fname="../arthropods/img/doubletree2.svg")
