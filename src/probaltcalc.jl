function probaltcalc(sigma,tmax,smax)


    #First half of season
    s_start = floor(Int64,(tmax)/2);
    #second half of season
    s_end = s_start + Int64(mod((tmax)/2,1)*2);

    transitionvec1 = pdf.(Normal(0,sigma),collect(-s_end+1-s_start:1:s_start));

    transitionvec2 = pdf.(Normal(0,sigma),collect(-s_end:1:s_start-1+s_end));

    probaltvec = [transitionvec1;transitionvec2];
    
    probaltvec = probaltvec ./ (maximum(probaltvec)*2);

    return probaltvec
end
