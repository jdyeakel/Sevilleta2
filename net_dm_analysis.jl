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

@everywhere include("$(homedir())/Dropbox/PostDoc/2020_Sevilleta2/src/laplacian.jl")
@everywhere include("$(homedir())/Dropbox/PostDoc/2020_Sevilleta2/src/eigencluster.jl")



path = string("$(homedir())/Dropbox/PostDoc/2020_Sevilleta2/mean_data_combined.csv");
data = CSV.read(path,header=true,DataFrame);

adata = Array(data);
#Get rid of data columns that are invariant
sumdiff=Array{Float64}(undef,size(adata)[2]);
for i=1:(size(adata)[2])
    sumdiff[i] = sum(diff(data[:,i]));
end
zeromeasures = findall(iszero,sumdiff);
adata = adata[:,setdiff(collect(1:size(adata)[2]),zeromeasures)]

nc = size(adata)[1];
nmeasures = size(adata)[2];




scaled_adata = Array{Float64}(undef,size(adata));
for i = 1:nmeasures
    scaled_adata[:,i] = (maximum(abs.(adata[:,i])) .- abs.(adata[:,i])) ./ (maximum(abs.(adata[:,i])) - minimum(abs.(adata[:,i])));
end

# adata = adata[:,[collect(1:8);collect(10:46)]]


#test - make row 2 similar to row 1
# normdist = Normal(0,0.01);
# adata[2,:] = adata[1,:] .* (1 .+ rand(normdist,nmeasures))


# PC = Array{Float64}(undef,nc,nc);
PCalt = Array{Float64}(undef,nc,nc);
# measmeans = mean(pcdatatr[!isnothing(pcdatatr)],1);

#Build similarity matrix
@showprogress 1 "Computing..." for i = 0:(nc^2 - 1)
    a = mod(i,nc) + 1;
    b = Int64(floor(i/nc)) + 1;
    
    #Calculate dissimilarity
    #If a == b, then dissimilarity is zero
    # if a == b
    #     PC[a,b] = 0.0;
    #     continue
    # end
    #alternative - compare matrices
    # 1) create average individual
    m1 = adata[a,:];
    m2 = adata[b,:];
    
    #
    # dist_m1m2 = 1 .- evaluate(Euclidean(),m1,m2);
    dist_m1m2 = 1 .- Jaccard()(m1,m2);
    # Euclidean
    # Jaccard
    # Euclidean => Jaccard > CosineDist
    # Other interesting ones:
    # Chebyshev
    # TotalVariation

    # global ct = 0;
    # global ctones = 0;
    # for j = 1:nmeasures
        
    #     global ct += (minimum([adata[a,j],adata[b,j]])/maximum([adata[a,j],adata[b,j]]));
    #     # different needs to be 0; same needs to be 1
    #     # global ct += log(1 - (sqrt((adata[a,j] - adata[b,j])^2)));

    #     # global ct += (sqrt((cvec[a,j] - cvec[b,j])^2));
    #     # global ct += (sqrt((m_cvec_wks[a,j] - m_cvec_wks[b,j])^2))
    #     global ctones += 1;
        
    # end
    # # ctscaled = (ctones - ct)/ctones;
    # ctscaled = exp((ct)/ctones);
    # PC[a,b] = Float64(ctscaled); #/Float64(ctones);s

    PCalt[a,b] = mean(dist_m1m2);    
end

PCalt[diagind(PCalt)].=0.0;

S = laplacian(PCalt,10);
ev = eigs(S; nev=10,which=:SR);
evalues = (ev[1]);
evecs = (ev[2]);

#scale evecs by evalues
scaled_evecs = copy(evecs);
for i = 1:size(scaled_evecs)[2]
    scaled_evecs[:,i] = evecs[:,i] ./ evalues[i];
end

ranked = sortperm(evecs[:,2]);
ecluster = eigencluster(collect(1:nc),evecs,3);

# Pairwise distances
pdist = Distances.pairwise(Jaccard(),scaled_evecs[:,2:10],dims=1);
# Jaccard
# Chebyshev
# BrayCurtis
pcafit = fit(PCA,pdist; maxoutdim=2)
tpdist = MultivariateStats.transform(pcafit, pdist)
scatterplot(tpdist[1,:],tpdist[2,:])
#NOTE: we get different results for scaled_evecs vs. evecs


