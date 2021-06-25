function [asdf_shuf] = shuf_isi(asdf)
% Input == spike train binned in asdf format. Generates
% shuffled spike timings keeping the distribution of 
% the ISI s same. The timing of the first spike is kept
% the same and all subsequent spike timings are changed

[s1 s2] = size(asdf);
isi_a = cellfun(@diff,(asdf),'uniformoutput',false);
  for ii = 1:max(s1,s2)-2
      N = length(isi_a{ii});
      isi_shuf = datasample(isi_a{ii},N,'replace',false);
      temp = horzcat(asdf{ii}(1,1),isi_shuf);
      asdf_shuf{ii,1} = cumsum(temp);
  end
  asdf_shuf{ii+1}= asdf{end-1};
  asdf_shuf{ii+2}(1,1) = asdf{end}(1,1);
  asdf_shuf{ii+2}(1,2) = asdf{end}(1,2);    
end

