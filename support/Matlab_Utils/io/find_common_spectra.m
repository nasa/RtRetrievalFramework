function [common_set,I,indices] = find_common_spectra(dat)

%%%  find indices of dat that correspond to l2/gfit data:
indices=struct('l2',[],'col',[],'ray',[],'grl',[],'array',[]);
[indices.l2] = find(cellfun(@(c) isfield(c,'fils'), dat));
[indices.col] = find(cellfun(@(c) isfield(c,'filcol'), dat));
[indices.ray] = find(cellfun(@(c) isfield(c,'Bend'), dat));
[indices.grl] = find(cellfun(@(c) isfield(c,'Hour'), dat));
[indices.array] = find(cellfun(@(c) isa(c,'double') | isa(c,'integer'), dat));

%%%  initialize common_set:
if (~isempty(indices.l2))
  common_set=dat{indices.l2(1)}.ispec;
elseif (~isempty(indices.col))
  common_set=dat{indices.col(1)}.col.spec_id;
elseif (~isempty(indices.ray))
  common_set=dat{indices.ray(1)}.spec_id;
elseif (~isempty(indices.grl))
  common_set=dat{indices.grl(1)}.spec_id;
elseif (~isempty(indices.array))
  common_set=dat{indices.array(1)};
else
  disp('error'); return;
end

%%%  intersect with all l2 and gfit spectrum IDs:

for il2=indices.l2
   common_set=intersect(dat{il2}.ispec, common_set);
end
for icol=indices.col
   common_set=intersect(dat{icol}.col.spec_id, common_set);
end
for iray=indices.ray
   common_set=intersect(dat{iray}.spec_id, common_set);
end
for igrl=indices.grl
   % first look for InGaAs spectrometers only:
   Ispec_InGaAs=find(cellfun(@(c) ~isempty(c),regexp(dat{igrl}.Spectrum_File_Name,'.*a[_\.].*')));
   common_set=intersect(dat{igrl}.spec_id(Ispec_InGaAs), common_set);
end
for iarray=indices.array
   common_set=intersect(dat{iarray}, common_set);
end

%%%  find mapping from dat indices to common_set

for il2=indices.l2
   I{il2}=find(ismember(dat{il2}.ispec,common_set));
end;
for icol=indices.col 
   I{icol}=find(ismember(dat{icol}.col.spec_id,common_set));
end;
for iray=indices.ray 
   I{iray}=find(ismember(dat{iray}.spec_id,common_set));
end;
for igrl=indices.grl 
   J = find(ismember(dat{igrl}.spec_id,common_set));
   % first look for InGaAs spectrometers only:
   Ispec_InGaAs=find(cellfun(@(c) ~isempty(c),regexp(dat{igrl}.Spectrum_File_Name,'.*a[_\.].*')));
   I{igrl}=intersect(J,Ispec_InGaAs);
end;
for iarray=indices.array 
   I{iarray}=find(ismember(dat{iarray},common_set));
end;
