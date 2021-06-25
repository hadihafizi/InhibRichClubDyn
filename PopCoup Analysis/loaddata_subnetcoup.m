%load sub network coupling files

% rich -rich
for ii = 1:89
    str = ['load popcoup_rr',num2str(ii)]; eval(str);
    Coup_r_r(ii,:) = x; clear x
end
save('Coup_r_r.mat','Coup_r_r')
clear ii

% nrich -nrich
for ii = 1:356
    str = ['load popcoup_nrnr',num2str(ii)]; eval(str);
    Coup_nr_nr(ii,:) = x; clear x
end
save('Coup_nr_nr.mat','Coup_nr_nr')

% rich -nrich
for ii = 1:89
    str = ['load popcoup_r_nr',num2str(ii)]; eval(str);
    Coup_r_nr(ii,:) = x; clear x
end
save('Coup_r_nr.mat','Coup_r_nr')

% nrich -rich
for ii = 1:356
    str = ['load popcoup_nr_r',num2str(ii)]; eval(str);
    Coup_nr_r(ii,:) = x; clear x
end
save('Coup_nr_r.mat','Coup_nr_r')