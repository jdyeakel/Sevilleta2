function initializevar(mass,maxcacheunit,res_kjg,epsilon,mseasons,alphaseasons,mscale,kmax,tmax,cycle,targetvalues,configs,mortboost,cache_behavior,areascaling)
  
  #Rescale energy to 10 kj (rather than 1 kj)
  # xscale = 20.0;
  
  #how many seasons?
  smax = length(cycle);
  
  #how many resources?
  num_res = length(res_kjg);
  
  tid = [0;repeat(collect(1:num_res),inner=(length(targetvalues),1))];
  tweight = [0;repeat(targetvalues,outer=(num_res,1))];
  tinfo = Tuple([tid,tweight]);
  

  #define the boundary conditions for energetic state
  # mass_starve = round(mass - ((0.02*mass^1.19)+(0.1*0.38*mass^1.0))); #from NSM
  mass_starve = (mass - ((0.02*mass^1.19)+(0.1*0.38*mass^1.0))); #from NSM

  ############
  # ENERGETICS
  ############

  #Joules per gram
  jpg=20000; #varies between 7000 to 36000
  kjpg = jpg/1000;

  #how many kj units does this organism have?
  # xmaxkj = Int64(mass-mass_starve)*kjpg; #convert grams to kJ

  xmaxkj = (mass-mass_starve)*kjpg; #convert grams to kJ

  #Constraint the energetic values over which the individual can vary.
  #If x falls below xc, organism is dead
  xc = 1;
  xmax = 20;
  
  #set scale such that xmax is always 20
  # xmax_pre = round(Int64,xc+xmaxkj-1); #unchanged if xc = 1
  xscale = xmaxkj/Float64(xmax);
  
  # #Convert kjg to 10kjg
  # global xmax = Int64(floor(xmax_pre/xscale));
  # 
  # #adaptive scaling for smaller body sizes
  # while xmax < 10
  #     global xscale = xscale*0.1;
  #     xmax_pre = round(Int64,xc+xmaxkj-1);
  #     global xmax = Int64(floor(xmax_pre/xscale));
  # end
  # 
  if cache_behavior == "cache"
      maxcacheunit = Int64(floor(((xmax*tmax*smax)/epsilon[2])/3)); #set at 1/3 infinite
      # maxcacheunit = 500;
  else
      maxcacheunit = 0;
  end


  #theta now in 10*kJ
  thetaunit = collect(0:maxcacheunit);
  thetamax = length(thetaunit);
  # decay = (1/20)*last(thetaunit); #in theta units
  
  #Digestibility of the cache
  epcache = epsilon[2]; #the digestibility is equal to that of seeds (2,4)
  
  #Metabolic constants for the basal and field metabolic rate
  b0_bmr = 0.018; #watts g^-0.75
  b0_fmr = 0.047; #watts g^-0.75
  
  #bout time in hours
  activehours = 5;
  nonactivehours = 24-activehours;
  #costs: f/df + sleeping over active hours
  cwh_df = (b0_bmr*(mass^0.75))*activehours + (b0_bmr*(mass^0.75))*nonactivehours; #watt*hour
  cwh_f = (b0_fmr*(mass^0.75))*activehours + (b0_bmr*(mass^0.75))*nonactivehours; #watt*hour

  #Convert to kiloJoules
  whrkJ = 3.6;
  #Convert kjg to 10kjg
  whrkJ = whrkJ/xscale;

  c_df = cwh_df*whrkJ;
  c_f = cwh_f*whrkJ;
  
  #Lets assume that the max amount you can eat in a bout is gainlimit*xmax
  # gainlimit = 1/2;
  # xs = xmax*gainlimit; #min is 85 (for 100g rodent), because of digestibility + cost constraints
  
  #OR assume you stomach can hold 50% more than you daily costs...
  #i.e. you stomach can add a max of 50% per bout after expenses
  xs = 1.5*c_f

  if maximum(epsilon)*xs < c_f
    error("maximum gain is less than costs... won't run correctly")
  end
  
  #Coefficient is for mass in KG, so convert to g/1000
  #m/s
  velocity = ((0.33/(1000^0.21)) * mass^0.21)/10;
  
  #Handling time (seconds)
  # ht = [50.0,100.0,50.0,100.0,150.0];
  ht = [0.0,0.0,0.0,0.0,0.0];

  #Bout in seconds
  tmax_bout = activehours*60*60;
  
  #Standard homerange
  A0 = 100;
  homerange = A0*mass^(areascaling); #square meters
  perchect = homerange/10000; #homerange to percent of a hectare (=10,000 m^2)

  
  ##########################
  #LANDSCAPE ATTRIBUTES
    # if cycle[1] == "monsoon"
    #     m_ghc_monsoon = mseasons[1];
    #     alpha_monsoon = alphaseasons[1];
    #     m_ghc_winter = mseasons[2];
    #     alpha_winter = alphaseasons[2];
    # else
    #     m_ghc_monsoon = mseasons[2];
    #     alpha_monsoon = alphaseasons[2];
    #     m_ghc_winter = mseasons[1];
    #     alpha_winter = alphaseasons[1];
    # end
    # 
    
    m_ghc_monsoon = Array(mseasons[:monsoon]);
    m_ghc_winter = Array(mseasons[:winter]);
    alpha_monsoon = Array(alphaseasons[:monsoon]);
    alpha_winter = Array(alphaseasons[:winter]);
    
  #grams per homerange
  m_gs_monsoon = m_ghc_monsoon .* perchect;
  m_gs_winter = m_ghc_winter .* perchect;
  
  #Scale m
  m_monsoon = mscale*m_gs_monsoon;
  m_winter = mscale*m_gs_winter;
  
  m = DataFrame([m_winter,m_monsoon],[:winter, :monsoon]);

  alpha = DataFrame([alpha_winter,alpha_monsoon],[:winter, :monsoon])

  p = DataFrame(convert(Matrix,alpha)./(convert(Matrix,alpha) .+ convert(Matrix,m)),[:winter, :monsoon])

  #convert alpha, m to c for Gamma Distribution
  c = DataFrame(convert(Matrix,alpha) ./ convert(Matrix,m), [:winter, :monsoon]);
  
  

  #Convert kjg to 10kjg
  resgain = res_kjg/xscale;
  
  configurations = configs;
  
  #Probability of mortality when staying home
  pmort_df = 1/365; #0; #1/365;
  #Probability of mortality when foraging outside
  pmort_f = mortboost*pmort_df;
  pmvec= [pmort_df;repeat([pmort_f],inner=length(tid))];
  
  catchsuccess = 0.1;
  
  
  #ksimulation distribution
  kdist_w,kinfo_w,tout_w,propres_w = ksim(num_res,alpha[:winter],c[:winter],ht,catchsuccess,resgain,epsilon,velocity,kmax,tmax_bout,configurations,tid,tweight);
  kdist_m,kinfo_m,tout_m,propres_m = ksim(num_res,alpha[:monsoon],c[:monsoon],ht,catchsuccess,resgain,epsilon,velocity,kmax,tmax_bout,configurations,tid,tweight);
  
  #NOTE: trim last (empty) column (c++ holdover I think)
  kdist_w = kdist_w[:,1:kmax];
  kinfo_w = kinfo_w[:,1:kmax];
  kdist_m = kdist_m[:,1:kmax];
  kinfo_m = kinfo_m[:,1:kmax];
  
  
  kmeans_w = vec(sum(kinfo_w .* kdist_w,dims=2));
  kmeans_m = vec(sum(kinfo_m .* kdist_m,dims=2));
  
  # 
  # R"""
  # par(mfrow=c(1,2))
  # barplot(as.numeric($kmeans_w),names.arg=$tid,col='gray')
  # barplot(as.numeric($kmeans_m),names.arg=$tid,col='gray')
  # """
  # 
  
  # @time kdist, kinfo, tinfo = ksim(num_res,alpha,c,ht,resgain,velocity,kmax,tmax_bout,configurations);
  # 
  # namespace = string("$(homedir())/2017_Sevilleta/model/figures/kdist.pdf");
  # R"""
  # # pdf($namespace,height=15,width=18)
  # par(mfrow=c(5,5))
  # plot($(kinfo_w[1,:]),$(kdist_w[1,:]),type='h',xlim=c(min($(kinfo_w[:,1:201])),max($(kinfo_w[:,1:201]))))
  # """
  # [R"plot($(kinfo_w[i,:]),$(kdist_w[i,:]),type='h',xlim=c(min($(kinfo_w[:,1:201])),max($(kinfo_w[:,1:201]))))" for i=2:length(tinfo[1])]; #R"dev.off()"
  
  return(
  num_res,
  kmax,
  xc,
  xmax,
  xscale,
  thetaunit,
  thetamax,
  epcache,
  c_df,
  c_f,
  xs,
  pmvec,
  kdist_w,
  kinfo_w,
  tout_w,
  propres_w,
  kdist_m,
  kinfo_m,
  tout_m,
  propres_m,
  tinfo,
  ht,
  velocity,
  tmax_bout,
  resgain)
  
end
