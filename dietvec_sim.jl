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
@everywhere include("$(homedir())/Dropbox/PostDoc/2020_Sevilleta2/src/probaltcalc.jl")
@everywhere include("$(homedir())/Dropbox/PostDoc/2020_Sevilleta2/src/ksim.jl")
@everywhere include("$(homedir())/Dropbox/PostDoc/2020_Sevilleta2/src/dailysim.jl")
@everywhere include("$(homedir())/Dropbox/PostDoc/2020_Sevilleta2/src/dailysimcomb.jl")
@everywhere include("$(homedir())/Dropbox/PostDoc/2020_Sevilleta2/src/laplacian.jl")
@everywhere include("$(homedir())/Dropbox/PostDoc/2020_Sevilleta2/src/eigencluster.jl")


mass = 100.; #mass of rodent
#Right now, we are maxing out cache at 'infinite', so the value here isn't used


#EMPIRICAL DATA (from good months) in g/m^2
# m_ghc_fall = [4650,581,15500,1938,6500];
# m_ghc_spring = [9920,1240,620,78,3400]; #s

#Load in empirical data
rdata = CSV.read("$(homedir())/Dropbox/PostDoc/2020_Sevilleta2/rdata2.csv",header=true,DataFrame);
rdata[!,:Fall_Mean] .+= 1;
rdata[!,:Spring_Mean] .+= 1;
rdata[!,:Fall_SD] .+= 1;
rdata[!,:Spring_SD] .+= 1;

#Subgroups :: C3 Perennial, C3 Annual, C4 Perennial, C4 Annual
groups_photopath = unique(rdata.PhotoPath);
groups_lifehistory = unique(rdata.LifeHistory);
groups_functionalgroup = unique(rdata.FunctionalGroup);

nsubgroups = length(groups_photopath)*length(groups_lifehistory)*length(groups_functionalgroup);

df_array = [DataFrame() for _ in 1:nsubgroups];
df_length = Array{Int64}(undef,nsubgroups) .* 0;
global df_id = Array{String}(undef,0);
let tic = 0
    for i=1:length(groups_photopath)
        for j=1:length(groups_lifehistory)
            for k=1:length(groups_functionalgroup)
                sub = findall((rdata.PhotoPath .== groups_photopath[i]) .* (rdata.LifeHistory .== groups_lifehistory[j]) .* (rdata.FunctionalGroup .== groups_functionalgroup[k]));
                if length(sub) > 0
                    tic += 1;
                    df_array[tic] = rdata[sub,:];
                    df_length[tic] = length(sub);
                    global df_id = push!(df_id,string(groups_photopath[i],"_",groups_lifehistory[j],"_",groups_functionalgroup[k]));
                end
            end
        end
    end
end
#Final Array of nonzero categories
df_array = df_array[findall(!iszero,df_length)]
ngroups = length(df_array);


#Simple mean approach
m_ghc_fall = [mean(df_array[i].Fall_Mean) for i=1:ngroups];
m_ghc_spring = [mean(df_array[i].Spring_Mean) for i=1:ngroups];
pl = lineplot(sort(m_ghc_fall,rev=true));
lineplot!(pl,sort(m_ghc_spring,rev=true))


#Grab Nitrogen Concentration
nconc = [mean(df_array[i].Mean_N[findall(!ismissing,df_array[i].Mean_N)]) for i=1:ngroups];
nconc_sd = [std(df_array[i].Mean_N[findall(!ismissing,df_array[i].Mean_N)]) for i=1:ngroups];
# nconc = rdata[rdata.LifeHistory .== "Perennial",:Mean_N];

# #Subselect Annuals
# rdata2[rdata2[!,:LifeHistory] .== "Annuals",:Spring_Mean]

# m_ghc_fall = rdata[!,:Fall_Mean];
# m_ghc_spring = rdata[!,:Spring_Mean];


# mseasons = DataFrame([m_ghc_spring,m_ghc_fall],[:spring, :fall]);

#number of resources
nr = length(m_ghc_fall);


#Targeting range
# targetvalues = [0.5,0.75,1.0]; #0.5,0.75,0.9
targetvalues = collect(0.1:0.05:1.0);
tid = [0;repeat(collect(1:nr),inner=(length(targetvalues),1))];
tweight = [0;repeat(targetvalues,outer=(nr,1))];
tinfo = Tuple([tid,tweight]);

#Alternative targeting for sets
# targetvalues = collect(0.5:0.25:1.0);
# orders = collect(combinations(collect(1:nr)));
# tid = [0;repeat(orders,inner=(length(targetvalues),1))];
# tweight = [0;repeat(targetvalues,outer=(length(orders),1))];
# tinfo = Tuple([tid,tweight]);

# #kilojoules per gram
res_kjg = repeat([mean([15.0,21.0,15.0,21.0,25.0])],outer=nr);

# #Percent digestible
epsilon = repeat([mean([0.33,0.75,0.25,0.75,0.77])],outer=nr);


#Coefficient is for mass in KG, so convert to g/1000
#m/s
velocity = ((0.33/(1000^0.21)) * mass^0.21)/10;

#Handling time (seconds)
# ht = [50.0,100.0,50.0,100.0,150.0];
ht = repeat([0.0],outer=nr);

#Bout in seconds
activehours = 5;
tmax_bout = activehours*60*60;

# m_ghc_fall = Array(mseasons[!,:fall]);
# m_ghc_spring = Array(mseasons[!,:spring]);
# alpha_fall = Array(alphaseasons[!,:fall]);
# alpha_spring = Array(alphaseasons[!,:spring]);

#Standard homerange
areascaling = 0.75;
A0 = 100;
homerange = A0*mass^(areascaling); #square meters
perchect = homerange/10000; #homerange to percent of a hectare (=10,000 m^2)

#grams per homerange
# m_gs_fall = m_ghc_fall .* perchect;
# m_gs_spring = m_ghc_spring .* perchect;

#Scale m
mscale = 1;

#Total adjustment
scaleadj = mscale*perchect;

# m_fall = mscale*m_gs_fall;
# m_spring = mscale*m_gs_spring;

#Estimate mean and variance of subgroups
m_fall = scaleadj * m_ghc_fall;
m_spring = scaleadj * m_ghc_spring;

#Revising estimates of variance = 0 or variance = NaN
#Right now, I'm just saying variance = mean in these cases! 4/6/2021
var_fall = [var(scaleadj * df_array[i].Fall_Mean) for i=1:ngroups];
var_spring = [var(scaleadj * df_array[i].Spring_Mean) for i=1:ngroups];
for i=1:ngroups
    if var_fall[i] == 0 || isnan(var_fall[i])
        var_fall[i] = mean(scaleadj * df_array[i].Fall_Mean);
    end
    if var_spring[i] == 0 || isnan(var_spring[i])
        var_spring[i] = mean(scaleadj * df_array[i].Spring_Mean);
    end
end 


#Estimate alpha from variance of subgroups
alpha_fall = m_fall.^2 ./ var_fall;
alpha_spring = m_spring.^2 ./ var_spring;


#Compare data vs. dist
# falldraw = rand(Gamma(alpha_est_fall[1],m_fall[1]/alpha_est_fall[1]),100)
# pl = lineplot(sort(df_array[1].Fall_Mean,rev=true));



m = DataFrame([m_spring,m_fall],[:spring, :fall]);
alpha = DataFrame([alpha_spring,alpha_fall],[:spring, :fall])
# p = DataFrame(convert(Matrix,alpha)./(convert(Matrix,alpha) .+ convert(Matrix,m)),[:spring, :fall])
#convert alpha, m to c for Gamma Distribution
# c = DataFrame(convert(Matrix,alpha) ./ convert(Matrix,m), [:spring, :fall]);

#Calculate negative binomial for each food
#Configurations to build ksim
configurations = 100000;
max_enc = 50;
kmax = max_enc+1;

#Currently not utilized (if set to 1)
catchsuccess = 1.0;

#ksimulation distribution
# kdist_s,kinfo_s,tout_s,propres_s = ksim(nr,alpha[!,:spring],c[!,:spring],ht,catchsuccess,res_kjg,epsilon,velocity,kmax,tmax_bout,configurations,tid,tweight);
# lineplot(kdist_s[1,:])

# kdist_f,kinfo_f,tout_f,propres_f = ksim(nr,alpha[!,:fall],c[!,:fall],ht,catchsuccess,res_kjg,epsilon,velocity,kmax,tmax_bout,configurations,tid,tweight);

# #NOTE: trim last (empty) column (c++ holdover I think)
# kdist_s = kdist_s[:,1:kmax];
# kinfo_s = kinfo_s[:,1:kmax];
# kdist_f = kdist_f[:,1:kmax];
# kinfo_f = kinfo_f[:,1:kmax];


# kmeans_s = vec(sum(kinfo_s .* kdist_s,dims=2));
# kmeans_f = vec(sum(kinfo_f .* kdist_f,dims=2));

#Daily sim draws across different consumer types



#SEASONAL UNCERTAINTY
cycle = ["fall","spring","fall"];
tmax = 100;
sigma = 20;
smax = length(cycle);
probaltvec = probaltcalc(sigma,tmax,smax);
ltime = length(probaltvec);



# probline_f = cumsum(kdist_f[consumertype,:]);
# probline_f = probline_f/maximum(probline_f);
# probline_s = cumsum(kdist_s[consumertype,:]);
# probline_s = probline_s/maximum(probline_s);


#STILL TO DO - GET PROPORTIONAL CONTRIBUTION DATA
reps = 10;
nc = length(tid);
ns = SharedArray{Float64}(reps,nc,ltime);
cvec = SharedArray{Float64}(reps,nc,nr*ltime);
dailyreturn = SharedArray{Float64}(reps,nc,ltime);
dailynitrogen = SharedArray{Float64}(reps,nc,ltime);
nrt = nr*ltime;
nrt2 = Int64(floor(ltime/7)*nr);
cvec_wks = SharedArray{Float64}(reps,nc,nrt2);


@time @sync @distributed for i=1:nc
    consumertype = i;
    for r=1:reps
        let day = 0, week = 0
            for s = 1:smax
                for t=1:tmax
                    day += 1;
                    rindex = (day-1)*nr + 1;
                    sdraw = rand();
                    fdraw = rand();
                    if cycle[s] == "fall"
                        if sdraw < probaltvec[day];
                            # dailyreturndraw = findall(x->x>fdraw,probline_f);
                            # dailyreturn[day] = kinfo_f[consumertype,dailyreturndraw[1]];
                            drdraw,dndraw,propres,numsuccess = dailysim(nr,alpha[!,:fall],m[!,:fall],ht,catchsuccess,res_kjg,nconc,velocity,tmax_bout,configurations,tid,tweight,consumertype);
                            # drdraw,dndraw,propres,numsuccess = dailysimcomb(nr,alpha[!,:fall],m[!,:fall],ht,catchsuccess,res_kjg,nconc,velocity,tmax_bout,configurations,tid,tweight,consumertype);
                            dailyreturn[r,i,day] = drdraw;
                            dailynitrogen[r,i,day] = dndraw;
                            if sum(propres) > 0
                                pctd = propres./sum(propres);
                            else
                                pctd = repeat([0.0],outer=nr);
                            end
                            cvec[r,i,rindex:(rindex+nr-1)] .= pctd;
                            ns[r,i,day] = sum(numsuccess);
                        else 
                            # dailyreturndraw = findall(x->x>fdraw,probline_s);
                            # dailyreturn[day] = kinfo_s[consumertype,dailyreturndraw[1]];
                            drdraw,dndraw,propres,numsuccess = dailysim(nr,alpha[!,:spring],m[!,:spring],ht,catchsuccess,res_kjg,nconc,velocity,tmax_bout,configurations,tid,tweight,consumertype);
                            # drdraw,dndraw,propres,numsuccess = dailysimcomb(nr,alpha[!,:spring],m[!,:spring],ht,catchsuccess,res_kjg,nconc,velocity,tmax_bout,configurations,tid,tweight,consumertype);
                            dailyreturn[r,i,day] = drdraw;
                            dailynitrogen[r,i,day] = dndraw;
                            if sum(propres) > 0
                                pctd = propres./sum(propres);
                            else
                                pctd = repeat([0.0],outer=nr);
                            end
                            cvec[r,i,rindex:(rindex+nr-1)] .= pctd;
                            ns[r,i,day] = sum(numsuccess);
                        end
                    else
                        if sdraw < probaltvec[day];
                            # dailyreturndraw = findall(x->x>fdraw,probline_s);
                            # dailyreturn[day] = kinfo_s[consumertype,dailyreturndraw[1]];
                            drdraw,dndraw,propres,numsuccess = dailysim(nr,alpha[!,:spring],m[!,:spring],ht,catchsuccess,res_kjg,nconc,velocity,tmax_bout,configurations,tid,tweight,consumertype);
                            # drdraw,dndraw,propres,numsuccess = dailysimcomb(nr,alpha[!,:spring],m[!,:spring],ht,catchsuccess,res_kjg,nconc,velocity,tmax_bout,configurations,tid,tweight,consumertype);
                            dailyreturn[r,i,day] = drdraw;
                            dailynitrogen[r,i,day] = dndraw;
                            if sum(propres) > 0
                                pctd = propres./sum(propres);
                            else
                                pctd = repeat([0.0],outer=nr);
                            end
                            cvec[r,i,rindex:(rindex+nr-1)] .= pctd;
                            ns[r,i,day] = sum(numsuccess);
                        else 
                            # dailyreturndraw = findall(x->x>fdraw,probline_f);
                            # dailyreturn[day] = kinfo_f[consumertype,dailyreturndraw[1]];
                            drdraw,dndraw,propres,numsuccess = dailysim(nr,alpha[!,:fall],m[!,:fall],ht,catchsuccess,res_kjg,nconc,velocity,tmax_bout,configurations,tid,tweight,consumertype);
                            # drdraw,dndraw,propres,numsuccess = dailysimcomb(nr,alpha[!,:fall],m[!,:fall],ht,catchsuccess,res_kjg,nconc,velocity,tmax_bout,configurations,tid,tweight,consumertype);
                            dailyreturn[r,i,day] = drdraw;
                            dailynitrogen[r,i,day] = dndraw;
                            if sum(propres) > 0
                                pctd = propres./sum(propres);
                            else
                                pctd = repeat([0.0],outer=nr);
                            end
                            cvec[r,i,rindex:(rindex+nr-1)] .= pctd;
                            ns[r,i,day] = sum(numsuccess);
                        end
                    end
                    #Build weekly average cvec matrix
                    if mod(day,7) == 0
                        week += 1;
                        windex = (week-1)*nr + 1;
                        # because we are taking means, sum will not be 1
                        weekly_pctd = vec(mean(reshape(cvec[r,i,rindex-(nr*6):(rindex+nr-1)],nr,7)',dims=1));
                        # Normalize to sum to 1
                        # norm_weekly_pctd = weekly_pctd ./ sum(weekly_pctd);
                        cvec_wks[r,i,windex:(windex+nr-1)] = weekly_pctd;
                    end

                end
            end
        end
    end
end

# Either take mean of individuals or the single rep
if reps > 1
    mns = ns[1,:,:];
    mcvec = cvec[1,:,:];
    mcvec_wks = cvec_wks[1,:,:];
    mdailyreturn = dailyreturn[1,:,:];
    mdailynitrogen = dailynitrogen[1,:,:];
else 
    mns = mean(ns,dims=1)[1,:,:];
    mcvec = mean(cvec,dims=1)[1,:,:];
    mcvec_wks = mean(cvec_wks,dims=1)[1,:,:];
    mdailyreturn = mean(dailyreturn,dims=1)[1,:,:];
    mdailynitrogen = mean(dailynitrogen,dims=1)[1,:,:];
end




fitness = vec(std(mdailyreturn,dims=2) ./ mean(mdailyreturn,dims=2));
fitness_nitro = vec(std(mdailynitrogen,dims=2) ./ mean(mdailynitrogen,dims=2));
# nitrofitness = xxx

# cvec_wks = cvec_wks[:,338:(337*2)];

# # 337 + 337
# nrt2 = size(cvec_wks)[2]

# Diffusion Mapping
ntime = Int64(nrt2/nr);

PC = Array{Float64}(undef,nc,nc);
PCalt = Array{Float64}(undef,nc,nc);
# measmeans = mean(pcdatatr[!isnothing(pcdatatr)],1);
#Build similarity matrix
@showprogress 1 "Computing..." for i = 0:(nc^2 - 1)
    a = mod(i,nc) + 1;
    b = Int64(floor(i/nc)) + 1;
    if a == b
        PC[a,b] = 0.0;
        continue
    end
    #alternative - compare matrices
    # 1) create average individual
    m1 = reshape(mcvec_wks[a,:],(nr,ntime));
    m2 = reshape(mcvec_wks[b,:],(nr,ntime));
    dist_m1m2 = 1 .- Distances.colwise(CosineDist(),m1,m2);

    global ct = 0;
    global ctones = 0;
    for j = 1:nrt2
        
        # global ct += log(minimum([mcvec_wks[a,j],mcvec_wks[b,j]])/maximum([mcvec_wks[a,j],mcvec_wks[b,j]]));
        # different needs to be 0; same needs to be 1
        global ct += log(1 - (sqrt((mcvec_wks[a,j] - mcvec_wks[b,j])^2)));

        # global ct += (sqrt((cvec[a,j] - cvec[b,j])^2));
        # global ct += (sqrt((mcvec_wks[a,j] - mcvec_wks[b,j])^2))
        global ctones += 1;
        
    end
    # ctscaled = (ctones - ct)/ctones;
    ctscaled = exp((ct)/ctones);
    PC[a,b] = Float64(ctscaled); #/Float64(ctones);s

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

#Identifiers
mdf_id = [string(df_id[i],"_",Int64(floor(mean([m_fall[i],m_spring[i]])))) for i=1:nr]
resnames = ["GEN";mdf_id];
# TO compare with tinfo
# rn = rdata[!,:kartez];
rn = mdf_id;



# Visualize full eigenvector landscape with PHATE!
# From Ashkaan: 
# (1) Get eigenvectors from your diffusion map
# (2) Calculate pairwise distances in diffusion space
# (3) Perform principal coordinates analysis on diffusion distances

# Pairwise distances
pdist = Distances.pairwise(Euclidean(),scaled_evecs[:,2:10],dims=1);
pcafit = fit(PCA,pdist; maxoutdim=2)
tpdist = MultivariateStats.transform(pcafit, pdist)
scatterplot(tpdist[1,:],tpdist[2,:])
#NOTE: we get different results for scaled_evecs vs. evecs

dfout = DataFrame(tpdist',[:pca1,:pca2]);
insert!(dfout,3,fitness_nitro,:fitness)
namespace = "$(homedir())/Dropbox/PostDoc/2020_Sevilleta2/data/scaled_eigenvecs.csv";
CSV.write(namespace,  dfout, writeheader=false)


#And plot!
namespace = "$(homedir())/Dropbox/PostDoc/2020_Sevilleta2/figures2/fig_fitnessmanifold.pdf";
R"""
library(RColorBrewer)
fitness_nitro = floor($((fitness_nitro)) * 100)
fitleg = seq(min(fitness_nitro),max(fitness_nitro),length.out=5)
pal = colorRampPalette(brewer.pal(11,"Spectral"))(max(fitness_nitro))
pdf($namespace,width=6,height=6)
plot($(tpdist[1,:]),$(tpdist[2,:]),col=pal[fitness_nitro],pch=16,cex=1.5)
legend(-3,3,legend=fitleg,pch=21,pt.bg=pal[fitleg],col='black',bty='n',pt.cex=2,title='CV(nitrogen)')
dev.off()
"""

namespace = "$(homedir())/Dropbox/PostDoc/2020_Sevilleta2/figures2/fig_dietmanifold.pdf";
R"""
pal = c('black',colorRampPalette(brewer.pal(9,"Set1"))(max($tid)))
alphaw = $tweight*100
alphaw[1] = 100;
palalpha = numeric(length($tid))
legvec = numeric(0)
for (i in 1:length($tid)) {
    if (alphaw[i] < 100) {
        palalpha[i] = paste(pal[($tid+1)][i],alphaw[i],sep='')
    } else {
        palalpha[i] = pal[($tid+1)][i]
        legvec = c(legvec,i)
    }
}
pdf($namespace,width=6,height=6)
plot($(tpdist[1,:]),$(tpdist[2,:]),col=palalpha,pch=16,cex=1.5)
legend(-3,3,legend=$(resnames),pch=21,pt.bg=palalpha[legvec],col='black',bty='n',pt.cex=2)
dev.off()
"""













#what is the relationship between nitro fitness and mean availability?
pl = scatterplot(repeat(m_spring,inner=3),fitness_nitro[2:22])
scatterplot!(pl,repeat(m_fall,inner=3),fitness_nitro[2:22])

#### OLD PLOTS ####

# scatterplot(evecs[:,2],evecs[:,3])




namespace = "$(homedir())/Dropbox/PostDoc/2020_Sevilleta2/figures2/consumereigen_2_3_many.pdf";
R"""
library(RColorBrewer)
pal = c('black',colorRampPalette(brewer.pal(9,"Set1"))(max($tid)))
alphaw = $tweight*100
alphaw[1] = 100;
palalpha = numeric(length($tid))
legvec = numeric(0)
for (i in 1:length($tid)) {
    if (alphaw[i] < 100) {
        palalpha[i] = paste(pal[($tid+1)][i],alphaw[i],sep='')
    } else {
        palalpha[i] = pal[($tid+1)][i]
        legvec = c(legvec,i)
    }
}
pdf($namespace,width=6,height=6)
plot($(evecs[:,2]),$(evecs[:,3]),pch=21,bg=palalpha,col='black',xlab='Laplacian eigenvec 2',ylab='Laplacian eigenvec 3',cex=2) #,xlim=c(-0.2,0.2),ylim=c(-0.2,0.2)
points($(evecs[:,2][1]),$(evecs[:,3][1]),pch=21,bg=palalpha[1],col='black',cex=2)
# legend(0.6,0.5,legend=$(resnames),pch=16,col=palalpha[legvec],cex=0.45,bty='n',pt.cex=1.5) #pal[($tid+1)]
legend(-1,0.55,legend=$(resnames),pch=21,pt.bg=palalpha[legvec],col='black',bty='n',pt.cex=2)
dev.off()
"""

# namespace = "$(homedir())/Dropbox/PostDoc/2020_Sevilleta2/figures2/consumereigen_scaled_2_3.pdf";
# R"""
# library(RColorBrewer)
# pal = c('black',colorRampPalette(brewer.pal(9,"Set1"))(max($tid)))
# alphaw = $tweight*100
# alphaw[1] = 100;
# palalpha = numeric(length($tid))
# legvec = numeric(0)
# for (i in 1:length($tid)) {
#     if (alphaw[i] < 100) {
#         palalpha[i] = paste(pal[($tid+1)][i],alphaw[i],sep='')
#     } else {
#         palalpha[i] = pal[($tid+1)][i]
#         legvec = c(legvec,i)
#     }
# }
# pdf($namespace,width=6,height=6)
# plot($(scaled_evecs[:,2]),$(scaled_evecs[:,3]),pch=21,bg=palalpha,col=palalpha,xlab='Laplacian eigenvec 2',ylab='Laplacian eigenvec 3') #,xlim=c(-0.2,0.2),ylim=c(-0.2,0.2)
# points($(scaled_evecs[:,2][1]),$(scaled_evecs[:,3][1]),pch=21,bg=palalpha[1],col=palalpha[1])
# legend(0.6,0.5,legend=$(resnames),pch=16,col=palalpha[legvec],cex=0.45,bty='n',pt.cex=1.5) #pal[($tid+1)]
# dev.off()
# """


namespace = "$(homedir())/Dropbox/PostDoc/2020_Sevilleta2/figures2/consumereigen_3_4_many.pdf";
R"""
library(RColorBrewer)
pal = c('black',colorRampPalette(brewer.pal(9,"Set1"))(max($tid)))
alphaw = $tweight*100
alphaw[1] = 100;
palalpha = numeric(length($tid))
legvec = numeric(0)
for (i in 1:length($tid)) {
    if (alphaw[i] < 100) {
        palalpha[i] = paste(pal[($tid+1)][i],alphaw[i],sep='')
    } else {
        palalpha[i] = pal[($tid+1)][i]
        legvec = c(legvec,i)
    }
}
pdf($namespace,width=6,height=6)
plot($(evecs[:,3]),$(evecs[:,4]),pch=21,bg=palalpha,col='black',xlab='Laplacian eigenvec 3',ylab='Laplacian eigenvec 4',cex=2)
points($(evecs[:,3][1]),$(evecs[:,4][1]),pch=21,bg=palalpha[1],col='black',cex=2)
# legend(0.4,0.6,legend=$(resnames),pch=16,col=palalpha[legvec],cex=0.45,bty='n',pt.cex=1.5) #pal[($tid+1)]
legend(-0.7,-0.1,legend=$(resnames),pch=21,pt.bg=palalpha[legvec],col='black',bty='n',pt.cex=2)
dev.off()
"""

# namespace = "$(homedir())/Dropbox/PostDoc/2020_Sevilleta2/figures2/consumereigen_scaled_3_4.pdf";
# R"""
# library(RColorBrewer)
# pal = c('black',colorRampPalette(brewer.pal(9,"Set1"))(max($tid)))
# alphaw = $tweight*100
# alphaw[1] = 100;
# palalpha = numeric(length($tid))
# legvec = numeric(0)
# for (i in 1:length($tid)) {
#     if (alphaw[i] < 100) {
#         palalpha[i] = paste(pal[($tid+1)][i],alphaw[i],sep='')
#     } else {
#         palalpha[i] = pal[($tid+1)][i]
#         legvec = c(legvec,i)
#     }
# }
# pdf($namespace,width=6,height=6)
# plot($(scaled_evecs[:,3]),$(scaled_evecs[:,4]),pch=21,bg=palalpha,col=palalpha,xlab='Laplacian eigenvec 2',ylab='Laplacian eigenvec 3') #,xlim=c(-0.2,0.2),ylim=c(-0.2,0.2)
# points($(scaled_evecs[:,3][1]),$(scaled_evecs[:,4][1]),pch=21,bg=palalpha[1],col=palalpha[1])
# legend(0.3,0.5,legend=$(resnames),pch=16,col=palalpha[legvec],cex=0.45,bty='n',pt.cex=1.5) #pal[($tid+1)]
# dev.off()
# """



# namespace = "$(homedir())/Dropbox/PostDoc/2020_Sevilleta2/consumer_ensilica_DM3D2.pdf";
# R"""
# library(scatterplot3d) 
# library(RColorBrewer)
# pal = c('black',colorRampPalette(brewer.pal(9,"Set1"))(max($tid)))
# alphaw = $tweight*100
# alphaw[1] = 100;
# palalpha = numeric(length($tid))
# legvec = numeric(0)
# for (i in 1:length($tid)) {
#     if (alphaw[i] < 100) {
#         palalpha[i] = paste(pal[($tid+1)][i],alphaw[i],sep='')
#     } else {
#         palalpha[i] = pal[($tid+1)][i]
#         legvec = c(legvec,i)
#     }
# }
# pdf($namespace,height=6,width=6)
# s3d = scatterplot3d(x=cbind($(evecs[:,2]),$(evecs[:,3]),$(evecs[:,4])),pch=16,color=palalpha,xlab='Laplacian eigenvec 2',ylab='Laplacian eigenvec 3',zlab='Laplacian eigenvec 4',scale.y=0.9,angle=70,type='h')
# dev.off()
# """



namespace = "$(homedir())/Dropbox/PostDoc/2020_Sevilleta2/figures2/efitnesseigen_2_3.pdf";
R"""
library(RColorBrewer)
fitness = floor($((fitness)) * 100)
fitleg = seq(10,max(fitness),length.out=5)
pal = colorRampPalette(brewer.pal(11,"Spectral"))(max(fitness))
pdf($namespace,width=6,height=6)
plot($(evecs[:,2]),$(evecs[:,3]),pch=21,bg=pal[fitness],col=pal[fitness],xlab='Laplacian eigenvec 2',ylab='Laplacian eigenvec 3')
legend(0.6,0.6,legend=fitleg,pch=16,col=pal[fitleg],cex=1,bty='n',pt.cex=1,title='CV(returns)')
dev.off()
"""


namespace = "$(homedir())/Dropbox/PostDoc/2020_Sevilleta2/figures2/efitnesseigen_3_4.pdf";
R"""
library(RColorBrewer)
fitness = floor($((fitness)) * 100)
fitleg = seq(1,max(fitness),length.out=5)
pal = colorRampPalette(brewer.pal(11,"Spectral"))(max(fitness))
pdf($namespace,width=6,height=6)
plot($(evecs[:,3]),$(evecs[:,4]),pch=21,bg=pal[fitness],col=pal[fitness],xlab='Laplacian eigenvec 3',ylab='Laplacian eigenvec 4')
legend(0.4,0.6,legend=fitleg,pch=16,col=pal[fitleg],cex=1,bty='n',pt.cex=1,title='CV(returns)')
dev.off()
"""



namespace = "$(homedir())/Dropbox/PostDoc/2020_Sevilleta2/figures2/nfitnesseigen_2_3_many.pdf";
R"""
library(RColorBrewer)
fitness_nitro = floor($((fitness_nitro)) * 100)
fitleg = seq(min(fitness_nitro),max(fitness_nitro),length.out=5)
pal = colorRampPalette(brewer.pal(11,"Spectral"))(max(fitness_nitro))
pdf($namespace,width=6,height=6)
plot($(evecs[:,2]),$(evecs[:,3]),pch=21,bg=pal[fitness_nitro],col='black',xlab='Laplacian eigenvec 2',ylab='Laplacian eigenvec 3',cex=2)
legend(-1,0.5,legend=fitleg,pch=21,pt.bg=pal[fitleg],col='black',bty='n',pt.cex=2,title='CV(returns)')
dev.off()
"""


namespace = "$(homedir())/Dropbox/PostDoc/2020_Sevilleta2/figures2/nfitnesseigen_3_4_many.pdf";
R"""
library(RColorBrewer)
fitness_nitro = floor($((fitness_nitro)) * 100)
fitleg = seq(min(fitness_nitro),max(fitness_nitro),length.out=5)
pal = colorRampPalette(brewer.pal(11,"Spectral"))(max(fitness_nitro))
pdf($namespace,width=6,height=6)
plot($(evecs[:,3]),$(evecs[:,4]),pch=21,bg=pal[fitness_nitro],col='black',xlab='Laplacian eigenvec 3',ylab='Laplacian eigenvec 4',cex=2)
legend(-0.7,-0.1,legend=fitleg,pch=21,pt.bg=pal[fitleg],col='black',bty='n',pt.cex=2,title='CV(nitrogen)')
dev.off()
"""


namespace = "$(homedir())/Dropbox/PostDoc/2020_Sevilleta2/figures2/nfitnesseigen_4_5_many.pdf";
R"""
library(RColorBrewer)
fitness_nitro = floor($((fitness_nitro)) * 100)
fitleg = seq(min(fitness_nitro),max(fitness_nitro),length.out=5)
pal = colorRampPalette(brewer.pal(11,"Spectral"))(max(fitness_nitro))
pdf($namespace,width=6,height=6)
plot($(evecs[:,4]),$(evecs[:,5]),pch=21,bg=pal[fitness_nitro],col='black',xlab='Laplacian eigenvec 3',ylab='Laplacian eigenvec 4',cex=2)
legend(-0.7,-0.1,legend=fitleg,pch=21,pt.bg=pal[fitleg],col='black',bty='n',pt.cex=2,title='CV(nitrogen)')
dev.off()
"""




namespace = "$(homedir())/Dropbox/PostDoc/2020_Sevilleta2/figures2/efitnesseigen_4_5.pdf";
R"""
library(RColorBrewer)
fitness = floor($((fitness)) * 100)
fitleg = seq(1,max(fitness),10)
pal = colorRampPalette(brewer.pal(11,"Spectral"))(max(fitness))
pdf($namespace,width=6,height=6)
plot($(evecs[:,4]),$(evecs[:,5]),pch=21,bg=pal[fitness],col=pal[fitness],xlab='Laplacian eigenvec 4',ylab='Laplacian eigenvec 5')
legend(0.10,0.0,legend=fitleg,pch=16,col=pal[fitleg],cex=1,bty='n',pt.cex=1,title='CV(returns)')
dev.off()
"""



#3D plot with plotly
R"""
library(plotly)
library(RColorBrewer)
pal = c('black',colorRampPalette(brewer.pal(9,"Set1"))(max($tid)))
alphaw = $tweight*100
alphaw[1] = 100;
palalpha = numeric(length($tid))
legvec = numeric(0)
for (i in 1:length($tid)) {
    if (alphaw[i] < 100) {
        palalpha[i] = paste(pal[($tid+1)][i],alphaw[i],sep='')
    } else {
        palalpha[i] = pal[($tid+1)][i]
        legvec = c(legvec,i)
    }
}
t <- list(
  family = "sans serif",
  size = 14,
  color = toRGB("grey50"))
species = $sp;
df = data.frame(species,$(evecs[:,2]),$(evecs[:,3]),$(evecs[:,4]));
colnames(df) = c('sp','ev2','ev3','ev4');
p <- plot_ly(df, x = ~ev2, y = ~ev3, z = ~ev4,
        mode = 'text',
        text = ~species,
        textposition = 'middle right',
        marker = list(color = ~ev2, colorscale = c('#FFE1A1', '#683531'), showscale = TRUE)) %>%
        add_markers() %>%
        add_text(textfont = t, textposition = "top right") %>%
  layout(scene = list(xaxis = list(title = 'ev2'),
                     yaxis = list(title = 'ev3'),
                     zaxis = list(title = 'ev4')),
         annotations = F)
"""

