function [S] = sim_spktrn(asdf)
%Calculates the similarity matrix for multiple spike
% trains. Input = = spike trains in asdf format. Bin
% == what bin size to use. Output == S; Similarity
% matrix with diagonals set to 0 and entries in only
% one side of the diagonal. 

[s1 s2] = size(asdf); 

S = zeros(max(s1,s2)-2,max(s1,s2)-2);
    for ii = 1:max(s1,s2)-2
        for jj =ii+1:max(s1,s2)-2
            inter = intersect(asdf{ii},asdf{jj});
            unio = union(asdf{ii},asdf{jj});
            S(ii,jj) = length(inter)/length(unio);
            clear inter unio
        end        
    end

end

