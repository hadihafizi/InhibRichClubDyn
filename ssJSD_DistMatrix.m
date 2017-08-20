function [ jsdMat ] = ssJSD_DistMatrix ( asdf, kernel )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

    duration = asdf{end}(2);
    N = length(asdf)-2;
    
    intspkints = cellfun(@(x) diff([reshape(x, numel(x),1) ; duration]), ...
        asdf(1:end-2), 'UniformOutput', 0);

    jsdMat = zeros(N);
    
    for ii=1:N
        disp(ii)
        for jj = (ii+1):N
            jsdMat(ii,jj) = findMaxDistSimilarity(intspkints{ii}, ...
                intspkints{jj}, 'Kernel', kernel);
        end
    end
    
end

