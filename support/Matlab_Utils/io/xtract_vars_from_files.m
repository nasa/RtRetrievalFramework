function [dat] = xtract_vars_from_files(files,vars)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% extract data from a wild-carded list of files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[dum1,dum2]=unix(['ls ' files  ' | cat']);
fils=struct([]);
r=dum2;
while (length(r)~=0)
  [this_file,r]=strtok(r);
  fils{end+1}=this_file;
end
if length(fils{end})==0
  fils(end)=[];
end

% map h5 field to structure member name
labels=cellfun(@(c) regexp(c,'.*\/(\S+)$','tokens','once'), vars);

% extract values
for If=1:length(fils)
    dum=regexp(fils{If},'.*_(.*).h5','tokens','once'); % extract index
    dat(If).sounding_id=str2num(dum{1});
    for Iv=1:length(vars)
        eval(['dat(If).' labels{Iv} '=h5read(fils{If},vars{Iv});']);
    end
end
