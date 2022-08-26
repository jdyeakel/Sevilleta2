function vectoresource(nvec)
    lvec = length(nvec);
    res_vec = Array{String}(undef,lvec);
    for i = 1:lvec
        if i < lvec
            res_vec[i] = "R" * string(nvec[i]) * " ";
        else
            res_vec[i] = "R" * string(nvec[i])
        end
    end
    fullvec = join(res_vec);
    finalvec = "(" * fullvec * ")"
    return finalvec
end

