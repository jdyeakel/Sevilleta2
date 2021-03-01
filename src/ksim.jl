function ksim(nr,alpha,c,ht,catchsuccess,res_kjg,epsilon,velocity,kmax,tmax_bout,configurations,tid,tweight)
    alpha = Array(alpha);
    c = Array(c);
    #Negative binomial parameters
    r = alpha;
    p = c ./ (c .+ 1);
    
    #mean of the negative binomial: pr/(1-p)
    #Var of the negative binomial: pr/(1-p)^2
    # The above is a suspect!!!
    
    gammadist = Gamma.(alpha,1 ./ c); #mean = alpha * (1/c)
    
    
    # negbindist = NegativeBinomial.(r,p);
    
    # configurations = 1000;
    
    # targetvalues = collect(0.25:0.25:1);
    # tid = [0;repeat(collect(1:nr),inner=(length(targetvalues),1))];
    # tweight = [0;repeat(targetvalues,outer=(nr,1))];
    # tid = Array(tinfo[1]);
    # tweight = Array(tinfo[2]);
    # 
    probability = SharedArray{Float64}(length(tid),kmax+1);
    kinfo = SharedArray{Float64}(length(tid),kmax+1);
    
    tdist = SharedArray{Float64}(length(tid));
    thandle = SharedArray{Float64}(length(tid));
    
    propres = SharedArray{Float64}(length(tid),nr);
    
    @sync @distributed for target=1:length(tid)
        
        bernoulidist = Bernoulli(1-tweight[target]);
        
        data = zeros(Float64,configurations);
        
        modvelocity = maximum([tweight[target],1/nr])*velocity;
        
        tdist[target] = 0.0;
        thandle[target] = 0.0;
        
        for config = 1:configurations
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
                    
                    t += distance_to_resource[tid[target]]/modvelocity;
                    tdist[target] += distance_to_resource[tid[target]]/modvelocity;
                    
                    #Obtains the resource if there is time left in tmax_bout
                    if tmax_bout > (distance_to_resource[tid[target]]/modvelocity + ht[tid[target]])
                        
                        #If not an insect, success is gauranteed
                        if tid[target] != 5
                            number_of_successes[tid[target]] += 1;
                            t += ht[tid[target]];
                            thandle[target] += ht[tid[target]];
                        #If an insect, success is not gauranteed
                        else 
                            catchinsect = rand();
                            if catchinsect < catchsuccess
                                number_of_successes[tid[target]] += 1;
                                t += ht[tid[target]];
                                thandle[target] += ht[tid[target]];
                            else
                                #No success; only time cost
                                # number_of_successes[tid[target]] += 0;
                                t += ht[tid[target]];
                                thandle[target] += ht[tid[target]];
                            end
                        end
                        
                    end
                
                else
                    #The rodent will move towards the closest resource
                    t += nearest_distance/modvelocity;
                    tdist[target] += nearest_distance/modvelocity;
                    
                    if tmax_bout > (nearest_distance/modvelocity)
                        number_of_successes[nearest_resource] += 1;
                        t += ht[nearest_resource];
                        thandle[target] += ht[nearest_resource];
                    end
                        
                    
                end
                
                #ORIGINAL
                # distancetuple = findmin(distance_to_resource);
                # nearest_distance = distancetuple[1];
                # nearest_resource = distancetuple[2];
                # 
                # if nearest_resource == tid[target]
                #     number_of_successes[tid[target]] += 1;
                #     t += ht[tid[target]];
                #     thandle[target] += ht[tid[target]];
                # else
                #     if rand(bernoulidist) == 1
                #         number_of_successes[nearest_resource] += 1;
                #         t += ht[nearest_resource];
                #         thandle[target] += ht[nearest_resource];
                #     end
                # end
                # t += nearest_distance/modvelocity;
                # tdist[target] += nearest_distance/modvelocity;
            end
            
            total_kilojoules=dot((res_kjg),number_of_successes); #.*epsilon
            # avg_digestibility = epsilon .* (number_of_successes/sum(number_of_successes));
            data[config]=total_kilojoules;
            
            propres[target,:] += ((res_kjg).*number_of_successes); #.*epsilon
            
        end
        datamin = minimum(data);
        datamax = maximum(data);
        if datamin!=datamax
            for config = 1:configurations
		        probability[target,convert(Int64,trunc(kmax*(data[config]-datamin)/(datamax-datamin)))+1]+=1.0/configurations;
            end
	        k = collect(datamin:(datamax-datamin)/kmax:datamax);
            for i = 1:kmax
	            kinfo[target,i]=k[i];
	        end
        else
            probability[target,1]=1.;
            k=datamin;
            kinfo[target,1]=datamin;
        end
    
    
        
    end
    
    propres /= configurations;
    tdist /= configurations;
    thandle /= configurations;
    
    tout = tuple(tdist,thandle);
    
    
    #note: probability matrix should have dims (kmax,length(tid))
 
    return (probability, kinfo, tout, propres);
end
