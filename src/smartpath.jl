function smartpath(filename::String,indices::Array{Int64}=[0])
    
    extension = filename[findall(x->x=='.',filename)[1]:length(filename)];
    filename_cut = replace(filename,extension=>"");
    
    if indices[1] == 0
        if homedir() == "/home/z840"
            namespace = string("$(homedir())/Sevilleta2/",filename_cut,extension);
        else
            # namespace = string("$(homedir())/Dropbox/Postdoc/2018_sharks/",filename_cut,extension);
            # namespace = string("$(homedir())/Dropbox/PostDoc/2018_sharks/",filename_cut,extension);
            namespace = string("$(homedir())/Dropbox/PostDoc/2020_Sevilleta2/",filename_cut,extension);
        end
    else
        
        sind = length(indices);
        
        indexstring = string();
        for i=1:sind
            indexstring = string(indexstring,"_",indices[i]);
        end
        
        if homedir() == "/home/z840"
            namespace = string("$(homedir())/Sevilleta2/",filename_cut,indexstring,extension);
        else
            # namespace = string("$(homedir())/Dropbox/Postdoc/2018_sharks/",filename_cut,indexstring,extension);
            # namespace = string("$(homedir())/Dropbox/PostDoc/2018_sharks/",filename_cut,indexstring,extension);
            namespace = string("$(homedir())/Dropbox/PostDoc/2020_Sevilleta2/",filename_cut,indexstring,extension);
        end
    end
    
    return namespace
    
end