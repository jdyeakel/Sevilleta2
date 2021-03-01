function kmeans(data,cnum)
    R"""
    rk<-kmeans($data,$cnum);
    rcluster = rk$cluster;
    clusterid = numeric(length(rcluster));
    for (i in 1:length(rcluster)) {
        clusterid[i] = as.numeric(rcluster[i]);
    }
    """
    @rget clusterid;
    clusterid = convert(Array{Int64},clusterid);
    numclusters = maximum(clusterid);
    speciescluster = Array{Array}(undef,numclusters);
    for i=1:numclusters
        speciescluster[i] = vec(findall(x->x==i,clusterid);)
    end
    return(speciescluster)
end
