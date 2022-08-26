function empiricaldataimport(datanamespace,rn)
    empdata = CSV.read(smartpath(datanamespace),header=true,DataFrame);
    num_ind = maximum(empdata[!,:ind]);
    start_t = findall(x->x=="t1",names(empdata))[1];
    end_t = length(names(empdata));
    tmax = end_t - start_t + 1;
    num_ind_res = Int64(size(empdata)[1]/num_ind);
    ind_diet = Array{Float64}(undef,num_ind,num_ind_res,tmax);
    for i=1:num_ind
        # indpos = findall(empdata.ind .== i);
        for j=1:num_ind_res
            indrespos = findall((empdata.ind .== i) .* (empdata.resource .== rn[j]))[1];
            ind_diet[i,j,:] = Array(empdata[indrespos,start_t:end_t]);
        end
    end
    return ind_diet
end

