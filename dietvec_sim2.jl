if homedir() == "/home/z840"
    loadfunc = include("$(homedir())/Sevilleta2/src/loadfuncs.jl");
else
    loadfunc = include("$(homedir())/Dropbox/PostDoc/2020_Sevilleta2/src/loadfuncs.jl");
end

#NOTE: NAME FILE APPENDIX
fileappend = "sigma1_r1000b_4wks_resourceset";


mass = 100.; #mass of rodent
#Right now, we are maxing out cache at 'infinite', so the value here isn't used


#EMPIRICAL DATA (from good months) in g/m^2
# m_ghc_fall = [4650,581,15500,1938,6500];
# m_ghc_spring = [9920,1240,620,78,3400]; #s

#Load in empirical data
rdata = CSV.read(smartpath("rdata2.csv"),header=true,DataFrame);
rdata[!,:Fall_Mean] .+= 1;
rdata[!,:Spring_Mean] .+= 1;
rdata[!,:Fall_SD] .+= 1;
rdata[!,:Spring_SD] .+= 1;
todel = findall(x->x=="CAM",rdata.PhotoPath);
# NOTE: REMOVE CAM
rdata = filter(row -> !(row.PhotoPath == "CAM"),  rdata)

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
# targetvalues = collect(0.1:0.1:1.0); #collect(0.1:0.05:1.0);
# tid = [0;repeat(collect(1:nr),inner=(length(targetvalues),1))];
# tweight = [0;repeat(targetvalues,outer=(nr,1))];
# tinfo = Tuple([tid,tweight]);

#NOTE: This is to run a combinatoric version of the foraging model
rescomb = Array{Array{Int64}}(undef,0);
function frescomb()
local rescombpos = Array{Int64}(undef,0);
local lrescomb = 0;
    for i=1:nr
        xarray = collect(combinations(1:nr,i));
        push!(rescomb,vcat(transpose.(xarray)...));
        lrescomb += size(rescomb[i],1);
        rescombpos = [rescombpos;repeat([i],outer=size(rescomb[i],1))];
    end
    return rescomb, lrescomb, rescombpos
end
rescomb, lrescomb, rescombpos = frescomb();
targetvalues = collect(0.2:0.2:1.0);
tid = Array{Array{Int64,1},1}(undef,0);
push!(tid,[0]);
for i=1:length(rescomb)
    lcomb = size(rescomb[i],1);
    for j = 1:lcomb
        for k=1:length(targetvalues)
            push!(tid,rescomb[i][j,:]);
        end
    end
end
tweight = [0;repeat(targetvalues,outer=(lrescomb,1))];
tinfo = Tuple([tid,tweight]);


#Alternative targeting for sets
# targetvalues = collect(0.5:0.25:1.0);
# orders = collect(combinations(collect(1:nr)));
# tid = [0;repeat(orders,inner=(length(targetvalues),1))];
# tweight = [0;repeat(targetvalues,outer=(length(orders),1))];
# tinfo = Tuple([tid,tweight]);

# #kilojoules per gram for C3 leaves, seeds, C4 leaves, seeds, insects
# res_kjg = repeat([mean([15.0,21.0,15.0,21.0,25.0])],outer=nr);
# res_kjg = [21.0,21.0,21.0,15.0,15.0,15.0,21.0];
# epsilon = [0.75,0.75,0.75,0.75,0.75,0.75,0.75];
res_kjg = [repeat([21.],outer=3);repeat([0.9*0.21],outer=3)]

# #Percent digestible
epsilon = repeat([mean([0.33,0.75,0.25,0.75,0.77])],outer=nr);


#Coefficient is for mass in KG, so convert to g/1000
#m/s
velocity = ((0.33/(1000^0.21)) * mass^0.21)/10;

#Metabolic rate (1/s)
#Metabolic constants for the basal and field metabolic rate
b0_bmr = 0.018; #watts g^-0.75
b0_fmr = 0.047; #watts g^-0.75
#bout time in hours
# activehours = 5;
# nonactivehours = 24-activehours;
metrate_basal = (b0_bmr*(mass^0.75));
metrate_active = (b0_fmr*(mass^0.75));
metrate = tuple(metrate_basal,metrate_active);
#costs: f/df + sleeping over active hours
# cwh_df = (b0_bmr*(mass^0.75))*activehours + (b0_bmr*(mass^0.75))*nonactivehours; #watt*hour
# cwh_f = (b0_fmr*(mass^0.75))*activehours + (b0_bmr*(mass^0.75))*nonactivehours; #watt*hour
#Convert to kiloJoules
# whrkJ = 3.6;
# c_df = cwh_df*whrkJ;
# c_f = cwh_f*whrkJ;

#energy stores
fat_g = ((0.02*mass^1.19);
muscle_g = (0.1*0.38*mass^1.0));
fat_kJperg = (37700)/1000;
protein_kJperg = (17900)/1000;
muscle_kJperg = 0.4*protein_kJperg + 0.6*fat_kJperg
stores_kj = fat_g*fat_kJperg; #+ muscle_g*muscle_kJperg

#Handling time (seconds)
# ht = [50.0,100.0,50.0,100.0,150.0];
ht_C3 = 60. *10.; #60 * number of minutes
ht_C4 = ht_C3*1.5;
ht = [repeat([ht_C3],outer=3);repeat([ht_C4],outer=3)];
# ht = repeat([60.0*10],outer=nr);

#Bout in seconds
activehours = 2;
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
m_fall = m_ghc_fall .* .001; #scaleadj * m_ghc_fall;
m_spring = m_ghc_spring .* .001; #scaleadj * m_ghc_spring;

#Revising estimates of variance = 0 or variance = NaN
#Right now, I'm just saying variance = mean in these cases! 4/6/2021
# var_fall = [var(scaleadj * df_array[i].Fall_Mean) for i=1:ngroups];
# var_spring = [var(scaleadj * df_array[i].Spring_Mean) for i=1:ngroups];
# for i=1:ngroups
#     if var_fall[i] == 0 || isnan(var_fall[i])
#         var_fall[i] = mean(scaleadj * df_array[i].Fall_Mean);
#     end
#     if var_spring[i] == 0 || isnan(var_spring[i])
#         var_spring[i] = mean(scaleadj * df_array[i].Spring_Mean);
#     end
# end 

var_fall = [mean((df_array[i].Fall_SD).^2) for i=1:ngroups] .* 0.00000001;
var_spring = [mean((df_array[i].Spring_SD).^2) for i=1:ngroups] .* 0.00000001;
# for i=1:ngroups
#     if var_fall[i] == 0 || isnan(var_fall[i])
#         var_fall[i] = mean(scaleadj * df_array[i].Fall_Mean);
#     end
#     if var_spring[i] == 0 || isnan(var_spring[i])
#         var_spring[i] = mean(scaleadj * df_array[i].Spring_Mean);
#     end
# end 


#Estimate alpha from variance of subgroups
alpha_fall = m_fall.^2 ./ var_fall;
alpha_spring = m_spring.^2 ./ var_spring;

#mean distance_to_resource in meters
meandist_fall = alpha_fall ./ (m_fall .* (alpha_fall .- 1))
meandist_spring = alpha_spring ./ (m_spring .* (alpha_spring .- 1))

#test
# id=3;
# gammadist = Gamma(alpha_fall[id], m_fall[id] / alpha_fall[id]);
# distance_to_resource = mean(rand.(Exponential.(1.0 ./ rand(gammadist,1000))))


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
sigma = 10;
smax = length(cycle);
probaltvec = probaltcalc(sigma,tmax,smax);
ltime = length(probaltvec);



# probline_f = cumsum(kdist_f[consumertype,:]);
# probline_f = probline_f/maximum(probline_f);
# probline_s = cumsum(kdist_s[consumertype,:]);
# probline_s = probline_s/maximum(probline_s);


#STILL TO DO - GET PROPORTIONAL CONTRIBUTION DATA
reps = 50;
nc = length(tid);
ns = SharedArray{Float64}(reps,nc,ltime);
cvec = SharedArray{Float64}(reps,nc,nr*ltime);
dailyreturn = SharedArray{Float64}(reps,nc,ltime);
dailynitrogen = SharedArray{Float64}(reps,nc,ltime);
consumerstores = SharedArray{Float64}(reps,nc,ltime);
meanfatmass = SharedArray{Float64}(reps,nc);
survival = SharedArray{Bool}(reps,nc);
nrt = nr*ltime;
timespan = 7*4;
nrt2 = Int64(floor(ltime/timespan)*nr);
cvec_wks = SharedArray{Float64}(reps,nc,nrt2);

# NOTE: could speed up by combining nc and reps into one loop
ncr_paramvec = [repeat(collect(1:nc),inner=reps) repeat(collect(1:reps),outer=nc)];
its = reps*nc;

@time @sync @distributed for ii=1:its
    i = ncr_paramvec[ii,1];
    r = ncr_paramvec[ii,2];
    consumertype = i;
    # for r=1:reps
    cstores = Array{Float64}(undef,ltime+1);
    cstores[1] = stores_kj;
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
                        drdraw,dndraw,propres,numsuccess,cost_kj = dailysim(nr,alpha[!,:fall],m[!,:fall],ht,catchsuccess,res_kjg,nconc,velocity,tmax_bout,configurations,tid,tweight,consumertype,metrate);
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
                        drdraw,dndraw,propres,numsuccess,cost_kj = dailysim(nr,alpha[!,:spring],m[!,:spring],ht,catchsuccess,res_kjg,nconc,velocity,tmax_bout,configurations,tid,tweight,consumertype,metrate);
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
                        drdraw,dndraw,propres,numsuccess,cost_kj = dailysim(nr,alpha[!,:spring],m[!,:spring],ht,catchsuccess,res_kjg,nconc,velocity,tmax_bout,configurations,tid,tweight,consumertype,metrate);
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
                        drdraw,dndraw,propres,numsuccess,cost_kj = dailysim(nr,alpha[!,:fall],m[!,:fall],ht,catchsuccess,res_kjg,nconc,velocity,tmax_bout,configurations,tid,tweight,consumertype,metrate);
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

                #Energetic dynamics
                cstores[day+1] = maximum([minimum([cstores[day] + drdraw - cost_kj,stores_kj]),0]);

                #Build average cvec matrix
                if mod(day,timespan) == 0
                    week += 1;
                    windex = (week-1)*nr + 1;
                    # because we are taking means, sum will not be 1
                    weekly_pctd = vec(mean(reshape(cvec[r,i,rindex-(nr*(timespan-1)):(rindex+nr-1)],nr,timespan)',dims=1));
                    # Normalize to sum to 1
                    # norm_weekly_pctd = weekly_pctd ./ sum(weekly_pctd);
                    cvec_wks[r,i,windex:(windex+nr-1)] = weekly_pctd;
                end

            end
        end
    end
    # lineplot(cstores)

    #Save whole trajectory
    consumerstores[r,i,1:ltime] = cstores[1:ltime];

    meancstores = mean(cstores);
    meanfatmass[r,i] = meancstores;
    survival[r,i] = (minimum(cstores) > 0);
    # end
end

propyearsurvive = Array{Float64}(undef,reps,nc);
for r=1:reps
    for i=1:nc
        zerotime = findall(x->x==0,consumerstores[r,i,:]);
        if length(zerotime) == 0
            propyearsurvive[r,i] = 1;
        else
            propyearsurvive[r,i] = zerotime[1]/ltime;
        end
    end
end

m_meanfatstores = vec(mean(meanfatmass,dims=1));


scatterplot(m_meanfatstores,log.(vec(mean(propyearsurvive,dims=1))))

scatterplot(vec(meanfatmass),vec(log.(propyearsurvive)))


R"""
m = lm($(vec(log.(propyearsurvive))) ~ $(vec(meanfatmass)))
summary(m)
"""

namespace = smartpath(string("figures2/fig_fatmass_vs_survival_",fileappend,".pdf"));
# string("$(homedir())/Dropbox/PostDoc/2020_Sevilleta2/figures2/fig_dietmanifold_",fileappend,".pdf");
R"""
library(RColorBrewer)
pdf($namespace,width=8,height=5)
par(mfrow=c(1,2))
plot($(vec(meanfatmass)),$(vec(log.(propyearsurvive))),pch='.',xlab='Mean fat mass',ylab='log Prop. year surv.')
plot($(m_meanfatstores),$(log.(vec(mean(propyearsurvive,dims=1)))),pch=16,xlab='Mean fat mass (avg reps)',ylab='log Prop. year surv. (avg. reps)')
dev.off()
"""

namespace = smartpath(string("figures2/fig_fatmass_vs_survival_",fileappend,"_small.pdf"));
# string("$(homedir())/Dropbox/PostDoc/2020_Sevilleta2/figures2/fig_dietmanifold_",fileappend,".pdf");
R"""
library(RColorBrewer)
pdf($namespace,width=4,height=4)
par(mfrow=c(1,1))
plot($(vec(meanfatmass)),$(vec(log.(propyearsurvive))),pch='.',xlab='Mean fat mass (kJ)',ylab='log Proportion survive')
abline(m)
text(15,-0.1,paste("R2=",round(summary(m)$adj.r.squared,2),sep=""))
dev.off()
"""


namespace = smartpath(string("figures2/fig_fatmassinds_",fileappend,".pdf"));
R"""
library(RColorBrewer)
pdf($namespace,width=4,height=4)
par(mfrow=c(1,1))
pal=brewer.pal(11,"Set3")
plot($(consumerstores[1,12,1:50])+10,type='l',col=pal[1],ylim=c(0,150),lwd=2,ylab="Endogenous stores (kJ)",xlab="Time (days)")
"""
for i=1:5
    R"""
    ri = sample(seq(1,11),1)
    lines($(consumerstores[i,12,1:50])+10,type='l',col=pal[ri],lwd=2)
    """
end
R"""
dev.off()
"""




# Either take mean of individuals or the single rep
if reps == 1
    m_ns = ns[1,:,:];
    m_cvec = cvec[1,:,:];
    m_cvec_wks = cvec_wks[1,:,:];
    m_dailyreturn = dailyreturn[1,:,:];
    m_dailynitrogen = dailynitrogen[1,:,:];
else 
    m_ns = mean(ns,dims=1)[1,:,:];
    m_cvec = mean(cvec,dims=1)[1,:,:];
    m_cvec_wks = mean(cvec_wks,dims=1)[1,:,:];
    m_dailyreturn = mean(dailyreturn,dims=1)[1,:,:];
    m_dailynitrogen = mean(dailynitrogen,dims=1)[1,:,:];
end


probsurvival = vec(mean(survival,dims=1))

returns_cv = vec(std(m_dailyreturn,dims=2) ./ mean(m_dailyreturn,dims=2));
nitro_cv = vec(std(m_dailynitrogen,dims=2) ./ mean(m_dailynitrogen,dims=2));
# nitrofitness = xxx



#Import empirical consumer data
empdatanamespace = "data/empdietraw.csv";
emp_diet = empiricaldataimport(empdatanamespace,rn);



# cvec_wks = cvec_wks[:,338:(337*2)];

# # 337 + 337
# nrt2 = size(cvec_wks)[2]

# Diffusion Mapping
nec = size(emp_diet)[1];
# Turn off empirical assessment
# nec = 0;
ntime = Int64(nrt2/nr);
# PC = Array{Float64}(undef,(nc+nec),(nc+nec));
PCalt = Array{Float64}(undef,(nc+nec),(nc+nec));
# measmeans = mean(pcdatatr[!isnothing(pcdatatr)],1);
#Build similarity matrix
@showprogress 1 "Computing..." for i = 0:((nc+nec)^2 - 1)
    a = mod(i,(nc+nec)) + 1;
    b = Int64(floor(i/(nc+nec))) + 1;
    # if a == b
    #     PC[a,b] = 0.0;
    #     continue
    # end
    #alternative - compare matrices
    # 1) create average individual
    if a <= nc && b <= nc
        m1 = reshape(m_cvec_wks[a,:],(nr,ntime));
        m2 = reshape(m_cvec_wks[b,:],(nr,ntime));
    end
    if a > nc && b <= nc
        emp_pos = a - nc;
        m1 = emp_diet[emp_pos,:,:];
        m2 = reshape(m_cvec_wks[b,:],(nr,ntime));
    end
    if a <= nc && b > nc
        emp_pos = b - nc;
        m1 = reshape(m_cvec_wks[a,:],(nr,ntime));
        m2 = emp_diet[emp_pos,:,:];
    end
    if a > nc && b > nc
        emp_pos_a = a - nc;
        emp_pos_b = b - nc;
        m1 = emp_diet[emp_pos_a,:,:];
        m2 = emp_diet[emp_pos_b,:,:];
    end

    dist_m1m2 = 1 .- Distances.colwise(Jaccard(),m1,m2);
    # Euclidean
    # Jaccard
    # Euclidean => Jaccard > CosineDist
    # Other interesting ones:
    # Chebyshev
    # TotalVariation

    # global ct = 0;
    # global ctones = 0;
    # for j = 1:nrt2
        
    #     # global ct += log(minimum([m_cvec_wks[a,j],m_cvec_wks[b,j]])/maximum([m_cvec_wks[a,j],m_cvec_wks[b,j]]));
    #     # different needs to be 0; same needs to be 1
    #     global ct += log(1 - (sqrt((m_cvec_wks[a,j] - m_cvec_wks[b,j])^2)));

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
ecluster = eigencluster(collect(1:(nc+nec)),evecs,3);

#Identifiers
# mdf_id = [string(df_id[i],"_",Int64(floor(mean([m_fall[i],m_spring[i]])))) for i=1:nr]
mdf_id = [string(df_id[i]) for i=1:nr]
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
pdist = Distances.pairwise(Jaccard(),scaled_evecs[:,2:10],dims=1);
# Jaccard
# Chebyshev
# BrayCurtis
pcafit = fit(PCA,pdist; maxoutdim=2)
tpdist_all = MultivariateStats.transform(pcafit, pdist)
scatterplot(tpdist_all[1,:],tpdist_all[2,:])
#NOTE: we get different results for scaled_evecs vs. evecs

#Model coordinates
tpdist = tpdist_all[:,1:nc];
#Empirical coordinates
tpdist_emp = tpdist_all[:,(nc+1):(nc+nec)]


dfout = DataFrame(tpdist',[:pca1,:pca2]);
dfout[!,:fitness] = m_meanfatstores;
#nitro_cv;
# insert!(dfout,3,nitro_cv,:fitness)
namespace = smartpath(string("data/scaled_eigenvecs_",fileappend,".csv"));
CSV.write(namespace,  dfout, writeheader=false)

dfout_emp = DataFrame(tpdist_emp',[:pca1,:pca2]);
namespace = smartpath(string("data/scaled_eigenvecs_",fileappend,"_emp.csv"));
CSV.write(namespace,  dfout_emp, writeheader=false)


namespace = smartpath(string("data/tid_",fileappend,".csv"));
ltid = vec(length.(tid));
CSV.write(namespace,  DataFrame([ltid],[:ltid]), writeheader=false)

# deigen = CSV.read(namespace,header=false, DataFrame);
# tpdist = Array(deigen[:,1:2])';
# nitro_cv = Array(deigen[:,3]);

namespace = smartpath(string("figures2/fig_fitness_resourceset_",fileappend,".pdf"));
R"""
pdf($namespace,width=4,height=4)
plot($ltid,$(log.(m_meanfatstores)),xlab='Size of resource set', ylab='Log stores',pch=16,cex=0.5)
dev.off()
"""


# namespace = smartpath(string("figures2/fig_dietmanifold_",fileappend,".pdf"));
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
# pdf($namespace,width=5,height=5)
# plot($(tpdist[1,:]),$(tpdist[2,:]),col='#55C9E9',bg=palalpha,pch=21,cex=2,xlab='Dimension 1',ylab='Dimension 2',lwd=2)
# for (i in 1:$nr) {
#     if (i != 50) {
#         findtids = which($tid == i)
#         tpdist = $tpdist
#         xspline = c(tpdist[1,1],tpdist[1,findtids])
#         yspline = c(tpdist[2,1],tpdist[2,findtids])
#         # lines(spline(xspline,yspline,n=5))
#         lines(xspline,yspline,lwd=8,col=paste(pal[i+1],'30',sep=''))
#     }
    
# }
# points($(tpdist[1,:]),$(tpdist[2,:]),col='#55C9E9',bg=palalpha,pch=21,cex=2,xlab='Dimension 1',ylab='Dimension 2',lwd=2)
# # points($(tpdist[1,:]),$(tpdist[2,:]),col='black',pch=1,cex=1.5)
# # points($(tpdist[1,:]),$(tpdist[2,:]),col=palalpha,pch=16,cex=1.5)
# legend(0.8,2,legend=$(resnames),pch=21,pt.bg=palalpha[legvec],col='#55C9E9',bty='n',pt.cex=1.2,cex=0.8,pt.lwd=2)
# dev.off()
# """

minx = findmin(dfout[!,:pca1]);
maxx = findmax(dfout[!,:pca1]);
miny = findmin(dfout[!,:pca2]);
maxy = findmax(dfout[!,:pca2]);



# NOTE: 
# string("$(homedir())/Dropbox/PostDoc/2020_Sevilleta2/figures2/fig_dietmanifold_",fileappend,".pdf");
namespace = smartpath(string("figures2/fig_dietmanifold_",fileappend,"_emp.pdf"));
R"""
library(RColorBrewer)
pal = c("#000000",colorRampPalette(brewer.pal(9,"YlOrRd"))($nr))
alphaw = $tweight*100
alphaw[1] = 100;
palalpha = numeric(length($tid))
legvec = numeric(0)
palalpha[1] = "#000000"
for (i in 2:length($tid)) {
    #Count number of resources
    numres = length($(tid)[[i]])
    if (alphaw[i] < 100) {
        palalpha[i] = paste(pal[numres+1],alphaw[i],sep='')
    } else {
        palalpha[i] = pal[numres+1]
        legvec = c(legvec,i)
    }
}
pdf($namespace,width=5,height=5)
plot($(tpdist[1,:]),$(tpdist[2,:]),col='#55C9E9',bg=palalpha,pch=21,cex=1.5,xlab='Dimension 1',ylab='Dimension 2',lwd=2,xlim=c(-6,10))
# for (i in 1:$nr) {
#     if (i != 50) {
#         resnumvec = as.numeric(lapply($tid,length))
#         findtids = which(resnumvec == i)
#         tpdist = $tpdist
#         xspline = c(tpdist[1,1],tpdist[1,findtids])
#         yspline = c(tpdist[2,1],tpdist[2,findtids])
#         # lines(spline(xspline,yspline,n=5))
#         lines(xspline,yspline,lwd=8,col=paste(pal[i+1],'30',sep=''))
#     }
    
# }
# points($(tpdist[1,:]),$(tpdist[2,:]),col='#55C9E9',bg=palalpha,pch=21,cex=2,xlab='Dimension 1',ylab='Dimension 2',lwd=1)
points($(tpdist[1,1]),$(tpdist[2,1]),col='#55C9E9',bg=palalpha[1],pch=21,cex=1.5,xlab='Dimension 1',ylab='Dimension 2',lwd=2)
points($(tpdist_emp[1,:]),$(tpdist_emp[2,:]),col='#2C3F99',bg="#2C3F9900",pch=23,cex=2.5,lwd=3)
# points($(tpdist[1,:]),$(tpdist[2,:]),col='black',pch=1,cex=1.5)
# points($(tpdist[1,:]),$(tpdist[2,:]),col=palalpha,pch=16,cex=1.5)
legend(-6,-2,legend=c("1-Set","2-Set","3-Set","4-Set","5-Set","6-Set","Gen"),pch=21,pt.bg=c(pal[2:($nr+1)],pal[1]),col='#55C9E9',bty='n',pt.cex=1.2,cex=0.8,pt.lwd=2)
# text(tpdist[1,$(miny[2])]+1.2,tpdist[2,$(miny[2])]+0.0,$(vectoresource(tid[miny[2]])))
# text(tpdist[1,$(minx[2])]+5.5,tpdist[2,$(minx[2])]+0.8,$(vectoresource(tid[minx[2]])))
# text(tpdist[1,$(maxx[2])]-2.,tpdist[2,$(maxx[2])]-0.2,$(vectoresource(tid[maxx[2]])))
dev.off()
"""




#Analysis of Clusters
numclusters = length(ecluster);
resfreq_cluster = Array{Float64}(undef,numclusters,nr);
resfreqw_cluster = Array{Float64}(undef,numclusters,nr);
empclusterpos = zeros(Int64,numclusters);
for i=1:numclusters
    cluster = ecluster[i];
    lcluster = length(cluster);
    resfreq = zeros(Float64,nr);
    resfreqw = zeros(Float64,nr);
    wcount = zeros(Float64,nr);
    count = zeros(Float64,nr);
    for j = 1:lcluster
        if cluster[j] <= nc
            for k = 1:nr
                resfreq[k] += length(findall(x->x==k,tid[cluster[j]]))
                resfreqw[k] += length(findall(x->x==k,tid[cluster[j]]))*tweight[cluster[j]];
                wcount[k] += tweight[cluster[j]];
                count[k] += 1.;
            end
        else
            empclusterpos[i] += 1
        end
    end
    resfreq = resfreq ./ count;
    resfreq_cluster[i,:] = resfreq;
    resfreqw = resfreqw ./ wcount;
    resfreqw_cluster[i,:] = resfreqw;
end


namespace = smartpath(string("figures2/fig_resclusters_",fileappend,".pdf"));
R"""
pdf($namespace,width=8,height=6)
par(mfrow=c(2,2))
barplot($(resfreq_cluster[1,:]),ylim=c(0,1))
"""
for i=2:numclusters
    R"""barplot($(resfreq_cluster[i,:]),ylim=c(0,1))"""
end
R"""
dev.off()
"""

namespace = smartpath(string("figures2/fig_resclusters_weighted_",fileappend,".pdf"));
R"""
clusternames = c("Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4")
pdf($namespace,width=8,height=6)
par(mfrow=c(2,2))
barplot($(resfreqw_cluster[1,:]),ylim=c(0,1),names=c("R1","R2","R3","R4","R5","R6"),main=clusternames[1])
"""
for i=2:numclusters
    R"""barplot($(resfreqw_cluster[i,:]),ylim=c(0,1),names=c("R1","R2","R3","R4","R5","R6"),main=clusternames[$i])"""
end
R"""
dev.off()
"""

clusterid = Array{Int64}(undef,nc);
for i=1:nc
    for j=1:numclusters
        if length(findall(x->x==i,ecluster[j])) == 1
            clusterid[i] = j;
        end
    end
end


# NOTE: 
# string("$(homedir())/Dropbox/PostDoc/2020_Sevilleta2/figures2/fig_dietmanifold_",fileappend,".pdf");
namespace = smartpath(string("figures2/fig_cluster_dietmanifold_",fileappend,"_emp.pdf"));
R"""
library(RColorBrewer)
pal = brewer.pal(4, "Set3")
pdf($namespace,width=5,height=5)
clusterid = $clusterid
plot($(tpdist[1,:]),$(tpdist[2,:]),col='black',bg=pal[clusterid],pch=21,cex=1.5,xlab='Dimension 1',ylab='Dimension 2',lwd=2,xlim=c(-6,10))
# for (i in 1:$nr) {
#     if (i != 50) {
#         resnumvec = as.numeric(lapply($tid,length))
#         findtids = which(resnumvec == i)
#         tpdist = $tpdist
#         xspline = c(tpdist[1,1],tpdist[1,findtids])
#         yspline = c(tpdist[2,1],tpdist[2,findtids])
#         # lines(spline(xspline,yspline,n=5))
#         lines(xspline,yspline,lwd=8,col=paste(pal[i+1],'30',sep=''))
#     }
    
# }
# points($(tpdist[1,:]),$(tpdist[2,:]),col='#55C9E9',bg=palalpha,pch=21,cex=2,xlab='Dimension 1',ylab='Dimension 2',lwd=1)
# points($(tpdist[1,1]),$(tpdist[2,1]),col='#55C9E9',bg=palalpha[1],pch=21,cex=1.5,xlab='Dimension 1',ylab='Dimension 2',lwd=2)
points($(tpdist_emp[1,:]),$(tpdist_emp[2,:]),col='#2C3F99',bg="#2C3F9900",pch=23,cex=2.5,lwd=3)
# points($(tpdist[1,:]),$(tpdist[2,:]),col='black',pch=1,cex=1.5)
# points($(tpdist[1,:]),$(tpdist[2,:]),col=palalpha,pch=16,cex=1.5)
legend(-6,-4,legend=c("Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4"),pch=21,pt.bg=pal,col='black',bty='n',pt.cex=1.2,cex=0.8,pt.lwd=2)
# text(tpdist[1,$(miny[2])]+1.2,tpdist[2,$(miny[2])]+0.0,$(vectoresource(tid[miny[2]])))
# text(tpdist[1,$(minx[2])]+5.5,tpdist[2,$(minx[2])]+0.8,$(vectoresource(tid[minx[2]])))
# text(tpdist[1,$(maxx[2])]-2.,tpdist[2,$(maxx[2])]-0.2,$(vectoresource(tid[maxx[2]])))
dev.off()
"""





# text(tpdist[1,$(miny[2])]+1.2,tpdist[2,$(miny[2])]+0.0,$(string(tid[miny[2]])))
# text(tpdist[1,$(minx[2])]+4.5,tpdist[2,$(minx[2])]+0.8,$(string(tid[minx[2]])))
# text(tpdist[1,$(maxx[2])]-1.,tpdist[2,$(maxx[2])]-0.2,$(string(tid[maxx[2]])))


# namespace = smartpath(string("figures2/fig_legend_",fileappend,".pdf"));
# # string("$(homedir())/Dropbox/PostDoc/2020_Sevilleta2/figures2/fig_dietmanifold_",fileappend,".pdf");
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
# plot($(tpdist[1,:]),$(tpdist[2,:]),col='#55C9E9',bg=palalpha,pch=21,cex=2,xlab='Dimension 1',ylab='Dimension 2',lwd=2)
# legend(0.7,0.4,legend=$(resnames),pch=21,pt.bg=palalpha[legvec],col='#55C9E9',bty='n',pt.cex=1,cex=1,pt.lwd=2)
# dev.off()
# """


#NOTE: currrently colors wrap around on this version of the plot!

# #And plot!
# namespace = smartpath(string("figures2/fig_fitnessmanifold_",fileappend,".pdf"));
# R"""
# library(RColorBrewer)
# # nitro_cv = floor($((nitro_cv)) * 100)
# m_meanfatstores = $m_meanfatstores
# # m_meanfatstores = ((m_meanfatstores)/max(m_meanfatstores))
# fitleg = round(seq(min(m_meanfatstores),max(m_meanfatstores),length.out=8),0)
# pal = colorRampPalette(brewer.pal(11,"Spectral"))(max(m_meanfatstores))
# legend_image <- as.raster(matrix(rev(pal), ncol=1))
# pdf($namespace,width=6,height=6)
# plot($(tpdist[1,:]),$(tpdist[2,:]),col='black',bg=pal[m_meanfatstores],pch=21,cex=1.5,xlab='Embedding dimension 1',ylab='Embedding dimension 2')
# # legend(0.85,0.45,legend=fitleg,pch=21,pt.bg=pal[fitleg],col='black',bty='n',pt.cex=1,title='Mean stores',cex=0.8)
# rasterImage(legend_image, 0.85,-0.1, 1,0.4)
# text(x=1.01, y = seq(-0.1,0.4,l=5), labels = seq(0,100,l=5),adj=0)
# text(x=0.95, y = 0.45, labels = "Mean stores (%)")

# dev.off()
# """


















# #what is the relationship between nitro fitness and mean availability?
# pl = scatterplot(repeat(m_spring,inner=3),nitro_cv[2:22])
# scatterplot!(pl,repeat(m_fall,inner=3),nitro_cv[2:22])

# #### OLD PLOTS ####

# # scatterplot(evecs[:,2],evecs[:,3])




# namespace = "$(homedir())/Dropbox/PostDoc/2020_Sevilleta2/figures2/consumereigen_2_3_many.pdf";
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
# plot($(evecs[:,2]),$(evecs[:,3]),pch=21,bg=palalpha,col='black',xlab='Laplacian eigenvec 2',ylab='Laplacian eigenvec 3',cex=2) #,xlim=c(-0.2,0.2),ylim=c(-0.2,0.2)
# points($(evecs[:,2][1]),$(evecs[:,3][1]),pch=21,bg=palalpha[1],col='black',cex=2)
# # legend(0.6,0.5,legend=$(resnames),pch=16,col=palalpha[legvec],cex=0.45,bty='n',pt.cex=1.5) #pal[($tid+1)]
# legend(-1,0.55,legend=$(resnames),pch=21,pt.bg=palalpha[legvec],col='black',bty='n',pt.cex=2)
# dev.off()
# """

# # namespace = "$(homedir())/Dropbox/PostDoc/2020_Sevilleta2/figures2/consumereigen_scaled_2_3.pdf";
# # R"""
# # library(RColorBrewer)
# # pal = c('black',colorRampPalette(brewer.pal(9,"Set1"))(max($tid)))
# # alphaw = $tweight*100
# # alphaw[1] = 100;
# # palalpha = numeric(length($tid))
# # legvec = numeric(0)
# # for (i in 1:length($tid)) {
# #     if (alphaw[i] < 100) {
# #         palalpha[i] = paste(pal[($tid+1)][i],alphaw[i],sep='')
# #     } else {
# #         palalpha[i] = pal[($tid+1)][i]
# #         legvec = c(legvec,i)
# #     }
# # }
# # pdf($namespace,width=6,height=6)
# # plot($(scaled_evecs[:,2]),$(scaled_evecs[:,3]),pch=21,bg=palalpha,col=palalpha,xlab='Laplacian eigenvec 2',ylab='Laplacian eigenvec 3') #,xlim=c(-0.2,0.2),ylim=c(-0.2,0.2)
# # points($(scaled_evecs[:,2][1]),$(scaled_evecs[:,3][1]),pch=21,bg=palalpha[1],col=palalpha[1])
# # legend(0.6,0.5,legend=$(resnames),pch=16,col=palalpha[legvec],cex=0.45,bty='n',pt.cex=1.5) #pal[($tid+1)]
# # dev.off()
# # """


# namespace = "$(homedir())/Dropbox/PostDoc/2020_Sevilleta2/figures2/consumereigen_3_4_many.pdf";
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
# plot($(evecs[:,3]),$(evecs[:,4]),pch=21,bg=palalpha,col='black',xlab='Laplacian eigenvec 3',ylab='Laplacian eigenvec 4',cex=2)
# points($(evecs[:,3][1]),$(evecs[:,4][1]),pch=21,bg=palalpha[1],col='black',cex=2)
# # legend(0.4,0.6,legend=$(resnames),pch=16,col=palalpha[legvec],cex=0.45,bty='n',pt.cex=1.5) #pal[($tid+1)]
# legend(-0.7,-0.1,legend=$(resnames),pch=21,pt.bg=palalpha[legvec],col='black',bty='n',pt.cex=2)
# dev.off()
# """

# # namespace = "$(homedir())/Dropbox/PostDoc/2020_Sevilleta2/figures2/consumereigen_scaled_3_4.pdf";
# # R"""
# # library(RColorBrewer)
# # pal = c('black',colorRampPalette(brewer.pal(9,"Set1"))(max($tid)))
# # alphaw = $tweight*100
# # alphaw[1] = 100;
# # palalpha = numeric(length($tid))
# # legvec = numeric(0)
# # for (i in 1:length($tid)) {
# #     if (alphaw[i] < 100) {
# #         palalpha[i] = paste(pal[($tid+1)][i],alphaw[i],sep='')
# #     } else {
# #         palalpha[i] = pal[($tid+1)][i]
# #         legvec = c(legvec,i)
# #     }
# # }
# # pdf($namespace,width=6,height=6)
# # plot($(scaled_evecs[:,3]),$(scaled_evecs[:,4]),pch=21,bg=palalpha,col=palalpha,xlab='Laplacian eigenvec 2',ylab='Laplacian eigenvec 3') #,xlim=c(-0.2,0.2),ylim=c(-0.2,0.2)
# # points($(scaled_evecs[:,3][1]),$(scaled_evecs[:,4][1]),pch=21,bg=palalpha[1],col=palalpha[1])
# # legend(0.3,0.5,legend=$(resnames),pch=16,col=palalpha[legvec],cex=0.45,bty='n',pt.cex=1.5) #pal[($tid+1)]
# # dev.off()
# # """



# # namespace = "$(homedir())/Dropbox/PostDoc/2020_Sevilleta2/consumer_ensilica_DM3D2.pdf";
# # R"""
# # library(scatterplot3d) 
# # library(RColorBrewer)
# # pal = c('black',colorRampPalette(brewer.pal(9,"Set1"))(max($tid)))
# # alphaw = $tweight*100
# # alphaw[1] = 100;
# # palalpha = numeric(length($tid))
# # legvec = numeric(0)
# # for (i in 1:length($tid)) {
# #     if (alphaw[i] < 100) {
# #         palalpha[i] = paste(pal[($tid+1)][i],alphaw[i],sep='')
# #     } else {
# #         palalpha[i] = pal[($tid+1)][i]
# #         legvec = c(legvec,i)
# #     }
# # }
# # pdf($namespace,height=6,width=6)
# # s3d = scatterplot3d(x=cbind($(evecs[:,2]),$(evecs[:,3]),$(evecs[:,4])),pch=16,color=palalpha,xlab='Laplacian eigenvec 2',ylab='Laplacian eigenvec 3',zlab='Laplacian eigenvec 4',scale.y=0.9,angle=70,type='h')
# # dev.off()
# # """



# namespace = "$(homedir())/Dropbox/PostDoc/2020_Sevilleta2/figures2/efitnesseigen_2_3.pdf";
# R"""
# library(RColorBrewer)
# fitness = floor($((fitness)) * 100)
# fitleg = seq(10,max(fitness),length.out=5)
# pal = colorRampPalette(brewer.pal(11,"Spectral"))(max(fitness))
# pdf($namespace,width=6,height=6)
# plot($(evecs[:,2]),$(evecs[:,3]),pch=21,bg=pal[fitness],col=pal[fitness],xlab='Laplacian eigenvec 2',ylab='Laplacian eigenvec 3')
# legend(0.6,0.6,legend=fitleg,pch=16,col=pal[fitleg],cex=1,bty='n',pt.cex=1,title='CV(returns)')
# dev.off()
# """


# namespace = "$(homedir())/Dropbox/PostDoc/2020_Sevilleta2/figures2/efitnesseigen_3_4.pdf";
# R"""
# library(RColorBrewer)
# fitness = floor($((fitness)) * 100)
# fitleg = seq(1,max(fitness),length.out=5)
# pal = colorRampPalette(brewer.pal(11,"Spectral"))(max(fitness))
# pdf($namespace,width=6,height=6)
# plot($(evecs[:,3]),$(evecs[:,4]),pch=21,bg=pal[fitness],col=pal[fitness],xlab='Laplacian eigenvec 3',ylab='Laplacian eigenvec 4')
# legend(0.4,0.6,legend=fitleg,pch=16,col=pal[fitleg],cex=1,bty='n',pt.cex=1,title='CV(returns)')
# dev.off()
# """



# namespace = "$(homedir())/Dropbox/PostDoc/2020_Sevilleta2/figures2/nfitnesseigen_2_3_many.pdf";
# R"""
# library(RColorBrewer)
# nitro_cv = floor($((nitro_cv)) * 100)
# fitleg = seq(min(nitro_cv),max(nitro_cv),length.out=5)
# pal = colorRampPalette(brewer.pal(11,"Spectral"))(max(nitro_cv))
# pdf($namespace,width=6,height=6)
# plot($(evecs[:,2]),$(evecs[:,3]),pch=21,bg=pal[nitro_cv],col='black',xlab='Laplacian eigenvec 2',ylab='Laplacian eigenvec 3',cex=2)
# legend(-1,0.5,legend=fitleg,pch=21,pt.bg=pal[fitleg],col='black',bty='n',pt.cex=2,title='CV(returns)')
# dev.off()
# """


# namespace = "$(homedir())/Dropbox/PostDoc/2020_Sevilleta2/figures2/nfitnesseigen_3_4_many.pdf";
# R"""
# library(RColorBrewer)
# nitro_cv = floor($((nitro_cv)) * 100)
# fitleg = seq(min(nitro_cv),max(nitro_cv),length.out=5)
# pal = colorRampPalette(brewer.pal(11,"Spectral"))(max(nitro_cv))
# pdf($namespace,width=6,height=6)
# plot($(evecs[:,3]),$(evecs[:,4]),pch=21,bg=pal[nitro_cv],col='black',xlab='Laplacian eigenvec 3',ylab='Laplacian eigenvec 4',cex=2)
# legend(-0.7,-0.1,legend=fitleg,pch=21,pt.bg=pal[fitleg],col='black',bty='n',pt.cex=2,title='CV(nitrogen)')
# dev.off()
# """


# namespace = "$(homedir())/Dropbox/PostDoc/2020_Sevilleta2/figures2/nfitnesseigen_4_5_many.pdf";
# R"""
# library(RColorBrewer)
# nitro_cv = floor($((nitro_cv)) * 100)
# fitleg = seq(min(nitro_cv),max(nitro_cv),length.out=5)
# pal = colorRampPalette(brewer.pal(11,"Spectral"))(max(nitro_cv))
# pdf($namespace,width=6,height=6)
# plot($(evecs[:,4]),$(evecs[:,5]),pch=21,bg=pal[nitro_cv],col='black',xlab='Laplacian eigenvec 3',ylab='Laplacian eigenvec 4',cex=2)
# legend(-0.7,-0.1,legend=fitleg,pch=21,pt.bg=pal[fitleg],col='black',bty='n',pt.cex=2,title='CV(nitrogen)')
# dev.off()
# """




# namespace = "$(homedir())/Dropbox/PostDoc/2020_Sevilleta2/figures2/efitnesseigen_4_5.pdf";
# R"""
# library(RColorBrewer)
# fitness = floor($((fitness)) * 100)
# fitleg = seq(1,max(fitness),10)
# pal = colorRampPalette(brewer.pal(11,"Spectral"))(max(fitness))
# pdf($namespace,width=6,height=6)
# plot($(evecs[:,4]),$(evecs[:,5]),pch=21,bg=pal[fitness],col=pal[fitness],xlab='Laplacian eigenvec 4',ylab='Laplacian eigenvec 5')
# legend(0.10,0.0,legend=fitleg,pch=16,col=pal[fitleg],cex=1,bty='n',pt.cex=1,title='CV(returns)')
# dev.off()
# """



# #3D plot with plotly
# R"""
# library(plotly)
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
# t <- list(
#   family = "sans serif",
#   size = 14,
#   color = toRGB("grey50"))
# species = $sp;
# df = data.frame(species,$(evecs[:,2]),$(evecs[:,3]),$(evecs[:,4]));
# colnames(df) = c('sp','ev2','ev3','ev4');
# p <- plot_ly(df, x = ~ev2, y = ~ev3, z = ~ev4,
#         mode = 'text',
#         text = ~species,
#         textposition = 'middle right',
#         marker = list(color = ~ev2, colorscale = c('#FFE1A1', '#683531'), showscale = TRUE)) %>%
#         add_markers() %>%
#         add_text(textfont = t, textposition = "top right") %>%
#   layout(scene = list(xaxis = list(title = 'ev2'),
#                      yaxis = list(title = 'ev3'),
#                      zaxis = list(title = 'ev4')),
#          annotations = F)
# """

