function eigencluster(sp,evecs,n)
    carray = Array{Array}(undef);
    carray[1] = sortperm(evecs[:,2]);
    println("Ordered List ","=",sp[carray[1]])
    for i=2:n
        eigenvec = evecs[:,i];
        rarray = Array{Array}(undef,0);
        for j = 1:length(carray)
            set21 = findall(x->x>0,eigenvec);
            set22 = findall(x->x<0,eigenvec);
    
            set21 = intersect(carray[j],set21);
            set22 = intersect(carray[j],set22);
            if length(set21) > 0
                push!(rarray,set21);
            end
            if length(set22) > 0
                push!(rarray,set22)
            end
        end
        carray = copy(rarray);
    end
    for i=1:length(carray)
        println()
        println("Cluster ",i,"=",sp[carray[i]])
    end
    return(carray)
end
