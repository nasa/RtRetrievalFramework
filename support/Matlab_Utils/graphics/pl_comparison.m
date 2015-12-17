function [h] = pl_comparison(type, dat, labels, domain_filter, xf, yf, ...
                               colors, do_legend)

if (~exist('colors')) % default values
    colors={[1 0 0],[0 1 0],[0 0 1],[0 0 0]};
end

int_of_str=containers.Map({dat(:).desc},num2cell(uint8([1:length(dat)])));

hold on; box on; 

if ( strcmp(type,'all') )
    h=zeros(0); desc=cell(0);

    for j=1:length(labels)
        I=int_of_str(labels{j});
        x_all = xf(dat(I));
        y_all = yf(dat(I));
        domain=domain_filter(xf(dat(I)));
        h(end+1)=plot(x_all(domain), y_all(domain),'Color',colors{j});
        desc{end+1}=dat(j).desc;
    end

elseif ( strcmp(type,'diff') )
    h=zeros(0); desc=cell(0);

    % first element in cell array is the reference value
    I0=int_of_str(labels{1});
    x0 = xf(dat(I0));
    y0 = yf(dat(I0));
    domain0=domain_filter(xf(dat(I0)));

    for j=2:length(labels)
        I=int_of_str(labels{j});
        x_all = xf(dat(I));
        y_all = yf(dat(I));
        domain=domain_filter(xf(dat(I)));
        
        if (length(domain0) ~= length(domain))
            disp('Error: array dimensions are not commensurate');
        end

        h(end+1)=plot(x_all(domain), y_all(domain) - y0(domain0),'Color',colors{j-1});
        desc{end+1}=[dat(I).desc '-' dat(I0).desc];
    end

elseif ( strcmp(type,'reldiff') )
    h=zeros(0); desc=cell(0);

    % first element in cell array is the reference value
    I0=int_of_str(labels{1});
    x0 = xf(dat(I0));
    y0 = yf(dat(I0));
    domain0=domain_filter(xf(dat(I0)));

    for j=2:length(labels)
        I=int_of_str(labels{j});
        x_all = xf(dat(I));
        y_all = yf(dat(I));
        domain=domain_filter(xf(dat(I)));
        
        if (length(domain0) ~= length(domain))
            disp('Error: array dimensions are not commensurate');
        end

        h(end+1)=plot(x_all(domain), y_all(domain) ./ y0(domain0) - 1,'Color',colors{j-1});
        desc{end+1}=[dat(I).desc '-' dat(I0).desc];
    end
end


if ( exist('do_legend') & do_legend==1)
    [legend_h, object_h, plot_h, text_strings] = legend(h,desc);
    text_objects=object_h(strcmp(get(object_h,'Type'),'text'));
    for j=1:length(plot_h)
        set(text_objects(j),'Color',get(plot_h(j),'Color'));
    end
end
