function dailysimcomb(nr,alpha,m,ht,catchsuccess,res_kjg,nconc,velocity,tmax_bout,configurations,tid,tweight,target)
        
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

    modvelocity = maximum([tweight[target],1/nr])*velocity;

    tdist[1] = 0.0;
    thandle[1] = 0.0;

    number_of_successes = zeros(Int64,nr);
    nearest_resource = 0;
    t=0.0;
    distance_to_resource = zeros(Float64,nr);
    nearest_distance = 0.0;    

    while t < tmax_bout
        
        for i=1:nr
            
            distance_to_resource[i] = rand(Exponential(1.0/rand(gammadist[i])));
            
        end
        
        #ALTERNATIVE
        distancetuple = findmin(distance_to_resource);
        nearest_distance = distancetuple[1];
        nearest_resource = distancetuple[2];
        
        if rand(bernoulidist) == 0
            #The rodent will move towards the targeted resource regardless if it's the closest

            #randomly choose NEAREST resource in its target group:
            chosendistancetuple = findmin(distance_to_resource[tid[target]]);
            chosenresource = tid[target][chosendistancetuple[2]];
            # chosenresource = rand(tid[target],1)[1];
            
            t += distance_to_resource[chosenresource]/modvelocity;
            tdist[1] += distance_to_resource[chosenresource]/modvelocity;
            
            #Obtains the resource if there is time left in tmax_bout
            if tmax_bout > (distance_to_resource[chosenresource]/modvelocity + ht[chosenresource])
                
                catchfood = rand();
                if catchfood < catchsuccess
                    number_of_successes[chosenresource] += 1;
                    t += ht[chosenresource];
                    thandle[1] += ht[chosenresource];
                else
                    #No success; only time cost
                    # number_of_successes[chosenresource] += 0;
                    t += ht[chosenresource];
                    thandle[1] += ht[chosenresource];
                end
                
                
            end
        
        else
            #The rodent will move towards the closest resource
            t += nearest_distance/modvelocity;
            tdist[1] += nearest_distance/modvelocity;
            
            if tmax_bout > (nearest_distance/modvelocity)
                number_of_successes[nearest_resource] += 1;
                t += ht[nearest_resource];
                thandle[1] += ht[nearest_resource];
            end
                
            
        end
        
        
    end

    total_kilojoules=dot((res_kjg),number_of_successes);
    total_nitrogen = dot(nconc,number_of_successes);
    

    propres = ((res_kjg).*number_of_successes);
    # propres = propres ./ sum(propres);

    return(total_kilojoules,total_nitrogen,propres,number_of_successes);

end
