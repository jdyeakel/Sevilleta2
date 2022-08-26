using DataFrames
using CSV
using UnicodePlots
using Distributed
using RCall
using Combinatorics
using MultivariateStats

@everywhere using ProgressMeter
@everywhere using Distributions
@everywhere using SharedArrays
@everywhere using LinearAlgebra
@everywhere using Arpack
@everywhere using Distances


if homedir() == "/home/z840"
    
    @everywhere include("$(homedir())/Sevilleta2/src/probaltcalc.jl")
    @everywhere include("$(homedir())/Sevilleta2/src/ksim.jl")
    # @everywhere include("$(homedir())/Sevilleta2/src/dailysim.jl")
    # @everywhere include("$(homedir())/Sevilleta2/src/dailysimcomb.jl")
    @everywhere include("$(homedir())/Sevilleta2/src/dailysimcombinatoric.jl")
    @everywhere include("$(homedir())/Sevilleta2/src/laplacian.jl")
    @everywhere include("$(homedir())/Sevilleta2/src/eigencluster.jl")
    @everywhere include("$(homedir())/Sevilleta2/src/smartpath.jl")
    @everywhere include("$(homedir())/Sevilleta2/src/vectoresource.jl")

else
    
    @everywhere include("$(homedir())/Dropbox/PostDoc/2020_Sevilleta2/src/probaltcalc.jl")
    @everywhere include("$(homedir())/Dropbox/PostDoc/2020_Sevilleta2/src/ksim.jl")
    # @everywhere include("$(homedir())/Dropbox/PostDoc/2020_Sevilleta2/src/dailysim.jl")
    # @everywhere include("$(homedir())/Dropbox/PostDoc/2020_Sevilleta2/src/dailysimcomb.jl")
    @everywhere include("$(homedir())/Dropbox/PostDoc/2020_Sevilleta2/src/dailysimcombinatoric.jl")
    @everywhere include("$(homedir())/Dropbox/PostDoc/2020_Sevilleta2/src/laplacian.jl")
    @everywhere include("$(homedir())/Dropbox/PostDoc/2020_Sevilleta2/src/eigencluster.jl")
    @everywhere include("$(homedir())/Dropbox/PostDoc/2020_Sevilleta2/src/smartpath.jl")
    @everywhere include("$(homedir())/Dropbox/PostDoc/2020_Sevilleta2/src/vectoresource.jl")

end
