function eigendistance(sp1,sp2,evecs,loc)
    
    dist = sqrt((evecs[loc,2][sp1] - evecs[loc,2][sp2])^2 + (evecs[loc,3][sp1] - evecs[loc,3][sp2])^2 + (evecs[loc,4][sp1] - evecs[loc,4][sp2])^2); #+ (evecs[loc,5][sp1] - evecs[loc,5][sp2])^2 
    
    return(dist)
    
end
    
