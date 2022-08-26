function dailysim(nr,alpha,m,ht,catchsuccess,res_kjg,nconc,velocity,tmax_bout,configurations,tid,tweight,target,metrate)
        
    alpha = Array(alpha);
    # c = Array(c);
    m = Array(m);
    #Negative binomial parameters
    # r = alpha;
    # p = c ./ (c .+ 1);

    gammadist = Gamma.(alpha, m ./ alpha); #mean = alpha * (1/c)


    # probability = SharedArray{Float64}(length(tid),kmax+1);
    # kinfo = SharedArray{Float64}(length(tid),kmax+1);

    tdist = Array{Float64}(undef,1);
    thandle = Array{Float64}(undef,1);

    # propres = Array{Float64}(undef,nr);

    bernoulidist = Bernoulli(1-tweight[target]);
            
    data = zeros(Float64,configurations);

    #Slows down the consumer that is generalizing on foods (search)
    modvelocity = maximum([tweight[target],1/nr])*velocity;

    tdist[1] = 0.0;
    thandle[1] = 0.0;

    number_of_successes = zeros(Int64,nr);
    nearest_resource = 0;
    t=0.0;
    distance_to_resource = zeros(Float64,nr);
    nearest_distance = 0.0;    

    #WATTS
    basal_mr = metrate[1];
    field_mr = metrate[2];

    #Metabolic cost in watt*seconds
    cost_ws = 0;

    while t < tmax_bout
        
        for i=1:nr
            
            distance_to_resource[i] = rand(Exponential(1.0/rand(gammadist[i])));
            
        end
        
        #ALTERNATIVE
        distancetuple = findmin(distance_to_resource);
        nearest_distance = distancetuple[1];
        nearest_resource = distancetuple[2];
        
        if rand(bernoulidist) == 0
            #The rodent will move towards the Closest of targeted resources regardless if it's the closest as a whole
            specializedtarget = tid[target];
            sdistance_to_resource = distance_to_resource[specializedtarget];
            sdistance_traveled = sum(sdistance_to_resource);
            # specializedtargetdistance = findmin(distance_to_resource[specializedtarget]);
            # snearest_distance = specializedtargetdistance[1];
            # snearest_resource = specializedtargetdistance[2];

            deltat = sdistance_traveled/modvelocity;
            t += deltat;
            tdist[1] += deltat;

            cost_ws += field_mr*deltat;
            
            #Obtains the resource if there is time left in tmax_bout
            if tmax_bout > (sdistance_traveled/modvelocity + sum(ht[specializedtarget]))
                
                #If not an insect, success is gauranteed
                if in(7,tid[target]) == false
                    number_of_successes[specializedtarget] .+= 1;
                    deltat = sum(ht[specializedtarget]);
                    t += deltat;
                    thandle[1] += deltat;
                    cost_ws += field_mr*deltat;
                #If an insect, success is not gauranteed
                else 
                    catchinsect = rand();
                    if catchinsect < catchsuccess
                        number_of_successes[specializedtarget] .+= 1;
                        deltat = sum(ht[specializedtarget]);
                        t += deltat;
                        thandle[1] += deltat;
                        cost_ws += field_mr*deltat;
                    else
                        #No success; only time cost
                        # number_of_successes[tid[target]] += 0;
                        deltat = sum(ht[specializedtarget]);
                        t += deltat;
                        thandle[1] += deltat;
                        cost_ws += field_mr*deltat;
                    end
                end
                
            end
        
        else
            #The rodent will move towards the closest resource
            deltat = nearest_distance/modvelocity;
            t += deltat;
            tdist[1] += deltat;
            cost_ws += field_mr*deltat;
            
            if tmax_bout > (nearest_distance/modvelocity)
                number_of_successes[nearest_resource] += 1;
                deltat = ht[nearest_resource];
                t += deltat;
                thandle[1] += deltat;
                cost_ws += field_mr*deltat;
            end
                
            
        end
        
        
    end

    total_t = 60*60*24;
    rest_t = total_t - t;
    cost_ws += basal_mr*rest_t;

    #Convert cost to kJ
    cost_kj = cost_ws*0.001;

    total_kilojoules=dot((res_kjg),number_of_successes);
    total_nitrogen = dot(nconc,number_of_successes);
    

    propres = ((res_kjg).*number_of_successes);
    # propres = propres ./ sum(propres);

    return(total_kilojoules,total_nitrogen,propres,number_of_successes,cost_kj);

end
