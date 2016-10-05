function [obj] = gfit_read(fil_exp,varargin)
%
% return data in structure
%

    [dum,fils]=unix(['ls ',fil_exp]);
    [fil,remainder] = strtok(fils);

    if (~isempty(strtok(remainder)))
       disp(strcat('Warning: file expression is multivalued: ',remainder));
       disp(strcat('\n... taking first value: ', fil));
    end

    fid=fopen(fil,'r');
    if (fid==-1)
       disp(strcat('Unable to read file <',fil,'>')); 
       obj=[];
       return
    end

    % determine file type and position file ptr to header info
    dum=regexp(fil,'\.','split');
    ext = dum(end);

    if (length(dum) > 1)
        if (strcmp(ext,'mav'))
            for k=1:3; fgetl(fid); end; % discard header
        elseif (strcmp(ext,'ray'))
            for k=1:2; fgetl(fid); end; % discard header
        elseif (strcmp(ext,'grl'))
            for k=1:2; fgetl(fid); end; % discard header
        elseif (strcmp(ext,'col'))
            for k=1:22; fgetl(fid); end; % discard header
	elseif (~isempty(regexp(fil,'.*(asc_).*','once')))  % measured spectrum file
            for k=1:5; fgetl(fid); end; % discard header
	    dat=textscan(fid,'%f%f');
	    fclose(fid);
	    obj.Frequency = dat{1};
	    obj.Signal = dat{2};
	    return;
	elseif (~isempty(regexp(fil,'z.*','once')))      % generated spectrum
	    dum=textscan(fgetl(fid),'%d%d%f%f');
	    obj.CL_exact = dum{3};  % this is not in original gfit product
	    obj.CT_exact = dum{4};  % this is not in original gfit product
            for k=1:1; fgetl(fid); end; % discard header
        end
    else
        disp('unknown file type');
        obj=[];
        return;
    end
    
    % get header strings and remove illegal characters
    headers = textscan(fgetl(fid),'%s');
    Nhdr = length(headers{:});
    for k = 1:Nhdr
        headers{1}{k} = regexprep(headers{1}{k},'-','_minus_');
        headers{1}{k} = regexprep(headers{1}{k},'/','_over_');
        headers{1}{k} = regexprep(headers{1}{k},'_0','_zero');
        if (isstrprop(headers{1}{k}(1),'digit'))
            headers{1}{k} = [ 'X' headers{1}{k}];
        end
    end

    entries=cell(0);
    while 1
        this_line = fgetl(fid);
        if ~ischar(this_line), break, end

	dum=textscan(this_line,'%s');
        if ( length(dum{1}) == Nhdr + 1 & strcmp(ext,'grl'))
            % L2 runlog has an extra column; ignore it
            dum = {dum{1}(1:end-1)};
        elseif ( length(dum{1}) ~= Nhdr )
           disp(['WARNING: inconsistent number of columns; ignoring following entry:  ' this_line]);
           continue;
        end

	entries{end+1}=dum{1};
    end
    fclose(fid);

    % convert from cell array of cell array to 2D cell array
    datstr=[entries{:}]';

    % convert from 2D cell array to 2D double array
    datdbl=cellfun(@str2double,datstr);
    col_string=find(max(isnan(datdbl)));
    col_dbl=setdiff(1:Nhdr,col_string);

    % create cell array, each element of which is a column of values
    for j=col_string
        eval(['obj.' headers{1}{j} '= datstr(:,j);']);
    end
    for j=col_dbl
        eval(['obj.' headers{1}{j} '= datdbl(:,j);']);
    end

    % add some derived data
    if ( strcmp(ext,'col') & isfield(obj,'Spectrum') )
        obj.spec_id=zeros(size(obj.Spectrum));
        for k=1:length(obj.Spectrum)
            obj.spec_id(k)=str2num(regexprep(obj.Spectrum{k},'.*\.',''));
        end
    elseif ( strcmp(ext,'ray') & isfield(obj,'SpectrumName') )
        obj.spec_id=zeros(size(obj.SpectrumName));
        for k=1:length(obj.SpectrumName)
            obj.spec_id(k)=str2num(regexprep(obj.SpectrumName{k},'.*\.',''));
        end
    elseif ( strcmp(ext,'grl') & isfield(obj,'Spectrum_File_Name') )
        obj.spec_id=zeros(size(obj.Spectrum_File_Name));
        for k=1:length(obj.Spectrum_File_Name)
            obj.spec_id(k)=str2num(regexprep(obj.Spectrum_File_Name{k},'.*\.',''));
        end
    end
