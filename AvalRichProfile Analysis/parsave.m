%
% parsave is capable of taking multiple variables and saving within a
% single file. The inputs are: 'fname' contains the desired filename; the
% input cell 'namesC' contains names of variables to be saved (e.g. namesC
% = {'x' 'y' 'z'}); the input cell 'valuesC' contains the variables to be
% saved (e.g. valuesC = {x y z}). I couldn't figure out an easy way of
% converting variable names into strings.
%
% Rashid Williams-Garcia, 2012

function parsave(fname,namesC,valuesC)
	varnames = genvarname(namesC);
    for i=1:numel(namesC)
        eval([varnames{i} ' = valuesC{i};'])
    end
    save(fname,varnames{:},'-v7.3')
% 	if numel(dir(fname))~=0
%         b = questdlg('Do you want to overwrite?','File already exist', ...
%             'Yes','No','No');
%         switch b
%             case 'Yes'
%                 save(fname,varnames{:})
%         end
% 	end

end

