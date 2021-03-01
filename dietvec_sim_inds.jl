using DataFrames
using CSV
using UnicodePlots
using Distributed
using RCall

@everywhere using ProgressMeter
@everywhere using Distributions
@everywhere using SharedArrays
@everywhere using LinearAlgebra
@everywhere using Arpack
@everywhere include("$(homedir())/Dropbox/PostDoc/2020_Sevilleta2/src/probaltcalc.jl")
@everywhere include("$(homedir())/Dropbox/PostDoc/2020_Sevilleta2/src/ksim.jl")
@everywhere include("$(homedir())/Dropbox/PostDoc/2020_Sevilleta2/src/dailysim.jl")
@everywhere include("$(homedir())/Dropbox/PostDoc/2020_Sevilleta2/src/laplacian.jl")
@everywhere include("$(homedir())/Dropbox/PostDoc/2020_Sevilleta2/src/eigencluster.jl")


mass = 100.; #mass of rodent
#Right now, we are maxing out cache at 'infinite', so the value here isn't used


#EMPIRICAL DATA (from good months) in g/m^2
# m_ghc_fall = [4650,581,15500,1938,6500];
# m_ghc_spring = [9920,1240,620,78,3400]; #s

#Load in empirical data
rdata = CSV.read("$(homedir())/Dropbox/PostDoc/2020_Sevilleta2/rdata.csv",header=true);
m_ghc_fall = rdata[!,:Fall_Mean];
m_ghc_spring = rdata[!,:Spring_Mean];
mseasons = DataFrame([m_ghc_spring,m_ghc_fall],[:spring, :fall]);

#number of resources
nr = length(m_ghc_fall);
alpha_fall = repeat([10.],outer=nr);
alpha_spring = repeat([10.],outer=nr);
alphaseasons = DataFrame([alpha_spring,alpha_fall],[:spring, :fall]);


#poor - rich
# p_r_mscale = [0.00002,0.0001]; #./ 3;

#kilojoules per gram
res_kjg = repeat([mean([15.0,21.0,15.0,21.0,25.0])],outer=nr);

#Percent digestible
epsilon = repeat([mean([0.33,0.75,0.25,0.75,0.77])],outer=nr);

#Targeting range
# targetvalues = [0.5,0.75,1.0]; #0.5,0.75,0.9
targetvalues = collect(0.5:0.25:1.0);
tid = [0;repeat(collect(1:nr),inner=(length(targetvalues),1))];
tweight = [0;repeat(targetvalues,outer=(nr,1))];
tinfo = Tuple([tid,tweight]);


#Coefficient is for mass in KG, so convert to g/1000
#m/s
velocity = ((0.33/(1000^0.21)) * mass^0.21)/10;

#Handling time (seconds)
# ht = [50.0,100.0,50.0,100.0,150.0];
ht = repeat([0.0],outer=nr);

#Bout in seconds
activehours = 5;
tmax_bout = activehours*60*60;

#Standard homerange
areascaling = 0.75;
A0 = 100;
homerange = A0*mass^(areascaling); #square meters
perchect = homerange/10000; #homerange to percent of a hectare (=10,000 m^2)


m_ghc_fall = Array(mseasons[!,:fall]);
m_ghc_spring = Array(mseasons[!,:spring]);
alpha_fall = Array(alphaseasons[!,:fall]);
alpha_spring = Array(alphaseasons[!,:spring]);
  
#grams per homerange
m_gs_fall = m_ghc_fall .* perchect;
m_gs_spring = m_ghc_spring .* perchect;

#Scale m
mscale = .1;
m_fall = mscale*m_gs_fall;
m_spring = mscale*m_gs_spring;

m = DataFrame([m_spring,m_fall],[:spring, :fall]);
alpha = DataFrame([alpha_spring,alpha_fall],[:spring, :fall])
p = DataFrame(convert(Matrix,alpha)./(convert(Matrix,alpha) .+ convert(Matrix,m)),[:spring, :fall])
#convert alpha, m to c for Gamma Distribution
c = DataFrame(convert(Matrix,alpha) ./ convert(Matrix,m), [:spring, :fall]);

#Calculate negative binomial for each food
#Configurations to build ksim
configurations = 100000;
max_enc = 50;
kmax = max_enc+1;

catchsuccess = 0.1;

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
reps = 1;
nc = length(tid);
ns = SharedArray{Float64}(reps,nc,ltime);
cvec = SharedArray{Float64}(reps,nc,nr*ltime);
dailyreturn = SharedArray{Float64}(reps,nc,ltime);
nrt = nr*ltime;
nrt2 = Int64(floor(ltime/7)*nr);
cvec_wks = SharedArray{Float64}(reps,nc,nrt2);


@sync @distributed for i=1:nc
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
                            drdraw,propres,numsuccess = dailysim(nr,alpha[!,:fall],c[!,:fall],ht,catchsuccess,res_kjg,velocity,tmax_bout,configurations,tid,tweight,consumertype);
                            dailyreturn[r,i,day] = drdraw;
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
                            drdraw,propres,numsuccess = dailysim(nr,alpha[!,:spring],c[!,:spring],ht,catchsuccess,res_kjg,velocity,tmax_bout,configurations,tid,tweight,consumertype);
                            dailyreturn[r,i,day] = drdraw;
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
                            drdraw,propres,numsuccess = dailysim(nr,alpha[!,:spring],c[!,:spring],ht,catchsuccess,res_kjg,velocity,tmax_bout,configurations,tid,tweight,consumertype);
                            dailyreturn[r,i,day] = drdraw;
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
                            drdraw,propres,numsuccess = dailysim(nr,alpha[!,:fall],c[!,:fall],ht,catchsuccess,res_kjg,velocity,tmax_bout,configurations,tid,tweight,consumertype);
                            dailyreturn[r,i,day] = drdraw;
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

#Unroll individuals
cvec_wksinds = Array{Float64}(undef,nc*reps,nrt2);
let tic = 0
    for i = 1:nc
        for r = 1:reps
            tic += 1
            cvec_wksinds[tic,:] = cvec_wks[r,i,:];
        end
    end
end





mns = mean(ns,dims=1)[1,:,:];
mcvec = mean(cvec,dims=1)[1,:,:];
mcvec_wks = mean(cvec_wks,dims=1)[1,:,:];
mdailyreturn = mean(dailyreturn,dims=1)[1,:,:];

fitness = vec(std(mdailyreturn,dims=2) ./ mean(mdailyreturn,dims=2));


# cvec_wks = cvec_wks[:,338:(337*2)];

# # 337 + 337
# nrt2 = size(cvec_wks)[2]

# Diffusion Mapping
ntime = Int64(nrt2/nr);

# PC = Array{Float64}(undef,nc*reps,nc*reps);
PCalt = Array{Float64}(undef,nc*reps,nc*reps);
# measmeans = mean(pcdatatr[!isnothing(pcdatatr)],1);
#Build similarity matrix
@showprogress 1 "Computing..." for i = 0:((nc*reps)^2 - 1)
    a = mod(i,nc*reps) + 1;
    b = Int64(floor(i/(nc*reps))) + 1;
    if a == b
        # PC[a,b] = 0.0;
        PCalt[a,b] = 0.0;
        continue
    end
    #alternative - compare matrices
    # 1) create average individual
    # m1 = reshape(cvec_wksinds[a,:],(nr,ntime));
    # m2 = reshape(cvec_wksinds[b,:],(nr,ntime));
    # dist_m1m2 = 1 .- Distances.colwise(CosineDist(),m1,m2);

    # m1 = vec(cvec_wksinds[a,:]);
    # m2 = vec(cvec_wksinds[b,:]);
    # dist_m1m2 = 1 .- Distances.evaluate(CosineDist(),m1,m2);

    global ct = 0;
    global ctones = 0;
    for j = 1:nrt2
        
        # global ct += log(minimum([mcvec_wks[a,j],mcvec_wks[b,j]])/maximum([mcvec_wks[a,j],mcvec_wks[b,j]]));
        # global ct += (sqrt((cvec[a,j] - cvec[b,j])^2));
        global ct += (sqrt((cvec_wksinds[a,j] - cvec_wksinds[b,j])^2))
        global ctones += 1;
        
    end
    ctscaled = (ctones - ct)/ctones;
    # ctscaled = (ct)/ctones;
    PCalt[a,b] = Float64(ctscaled); #/Float64(ctones);s

    # PCalt[a,b] = mean(dist_m1m2);
end

S = laplacian(PCalt,10);
ev = eigs(S; nev=10,which=:SR);
evalues = ev[1];
evecs = ev[2];

ranked = sortperm(evecs[:,2]);

ecluster = eigencluster(collect(1:(nc*reps)),evecs,4);
resnames = ["GEN";rdata[!,:kartez]];

# scatterplot(evecs[:,2],evecs[:,3])


namespace = "$(homedir())/Dropbox/PostDoc/2020_Sevilleta2/consumer_ensilica_DM_inds.pdf";

R"""
library(RColorBrewer)
pal = c('black',colorRampPalette(brewer.pal(9,"Set1"))(max($tid)))
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
pdf($namespace,width=6,height=6)
plot($(evecs[:,2]),$(evecs[:,3]),pch='.',xlab='Laplacian eigenvec 2',ylab='Laplacian eigenvec 3')
dev.off()
"""


namespace = "$(homedir())/Dropbox/PostDoc/2020_Sevilleta2/consumer_ensilica_DM_inds2.pdf";

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
plot($(evecs[:,3]),$(evecs[:,4]),pch=16,xlab='Laplacian eigenvec 3',ylab='Laplacian eigenvec 4')
dev.off()
"""



namespace = "$(homedir())/Dropbox/PostDoc/2020_Sevilleta2/consumer_ensilica_DM3D2.pdf";
R"""
library(scatterplot3d) 
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
pdf($namespace,height=6,width=6)
s3d = scatterplot3d(x=cbind($(evecs[:,2]),$(evecs[:,3]),$(evecs[:,4])),pch=16,color=palalpha,xlab='Laplacian eigenvec 2',ylab='Laplacian eigenvec 3',zlab='Laplacian eigenvec 4',scale.y=0.9,angle=70,type='h')
dev.off()
"""



namespace = "$(homedir())/Dropbox/PostDoc/2020_Sevilleta2/consumer_fitness_DM2.pdf";
R"""
library(RColorBrewer)
fitness = floor($((fitness)) * 100)
fitleg = seq(10,100,10)
pal = colorRampPalette(brewer.pal(11,"Spectral"))(max(fitness))
pdf($namespace,width=6,height=6)
plot($(evecs[:,2]),$(evecs[:,3]),pch=21,bg=pal[fitness],col=pal[fitness],xlab='Laplacian eigenvec 2',ylab='Laplacian eigenvec 3')
legend(0.10,0.0,legend=fitleg,pch=16,col=pal[fitleg],cex=1,bty='n',pt.cex=1,title='CV(returns)')
dev.off()
"""


namespace = "$(homedir())/Dropbox/PostDoc/2020_Sevilleta2/consumer_fitness_DM22.pdf";
R"""
library(RColorBrewer)
fitness = floor($((fitness)) * 100)
fitleg = seq(1,max(fitness),10)
pal = colorRampPalette(brewer.pal(11,"Spectral"))(max(fitness))
pdf($namespace,width=6,height=6)
plot($(evecs[:,3]),$(evecs[:,4]),pch=21,bg=pal[fitness],col=pal[fitness],xlab='Laplacian eigenvec 3',ylab='Laplacian eigenvec 4')
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
