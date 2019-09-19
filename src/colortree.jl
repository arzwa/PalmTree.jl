
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
