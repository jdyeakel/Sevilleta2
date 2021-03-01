function laplacian(PC,branchnum)
    
    nsp = size(PC)[1];
    
    maxdiff = minimum([branchnum,nsp-1]);

    S = zeros(Float64,nsp,nsp);
    # S = copy(-PC);
    for i = 1:nsp

        val = zeros(Float64,maxdiff);
        lok = zeros(Int64,maxdiff);

        for j = 1:nsp
            if PC[i,j] > val[maxdiff]
                val[maxdiff] = PC[i,j];
                lok[maxdiff] = j;
            end
            for k=1:maxdiff-1
                if val[maxdiff+1-k] > val[maxdiff-k]
                    v = val[maxdiff+1-k];
                    val[maxdiff+1-k] = val[maxdiff-k];
                    val[maxdiff-k] = v;
                    l = lok[maxdiff+1-k];
                    lok[maxdiff+1-k] = lok[maxdiff-k];
                    lok[maxdiff-k] = l;
                end
            end
        end
        S[i,lok] = -PC[i,lok];
        S[lok,i] = -PC[lok,i];
    end
    rowsums = sum(S,dims=2);
    S[diagind(S)] = -rowsums;
    
    return S
end