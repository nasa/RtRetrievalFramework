function prettyPlot(usr)

%  prettyPlot:  [] <-- (string usr)
%
%==============================================================================
%
%  DESCRIPTION:
%
%    beautifies all objects in current figure.  Takes one argument
%    which is a string containing variables to be set.
%
%  INPUT:
%
%    string usr:  string evaluated in matlab; is used to define following
%                 variables:
%
%                 X   = xlimits (leave alone)
%                 Y   = ylimits (leave alone)
%                 Z   = zlimits (leave alone)
%                 sh  = shape {'square','rectangle'}  ('rectangle')
%                 sty = style {'journal','slide', 'notebook'}  ('journal')
%                 lw  = linewidth (1)
%                 ms  = markersize (6)
%                 fs  = fontsize (14 (slide), 9 (journal))
%                 lgs = legendsize (11 (slide), 6 (journal))
%                 fn  = fontname ('Helvetica')
%                 ps  = EPS file; if NULL don''t print to file ('')
%                 bx  = put box around figure? (leave alone)
%                 gr  = turn on grid? (0)
%                 nl  = no label? (0)
%
%  EXAMPLE:       prettyPlot('lw=0.2; fs=9; ps=''dum'';');
%
%  DEPENDENCIES:  none
%
%==============================================================================

if ( exist('usr') );  eval(usr);  end

% default values
if ( ~exist('sty') );  sty='journal';  end
if ( ~exist('sh') );  sh='rectangle';  end
if ( ~exist('nl') );  nl=0;  end

if ( strcmp(sty,'journal') )
   if ( ~exist('ms') );  ms=6;  end
   if ( ~exist('lw') );  lw=0.2;  end
   if ( ~exist('fs') );  fs=9;  end
   if ( ~exist('lgs') ); lgs=6;  end
elseif ( strcmp(sty,'slide') )
   if ( ~exist('ms') );  ms=10;  end
   if ( ~exist('lw') );  lw=2;  end
   if ( ~exist('fs') );  fs=24;  end
   if ( ~exist('lgs') ); lgs=9;  end
elseif ( strcmp(sty,'notebook') )
   if ( ~exist('ms') );  ms=6;  end
   if ( ~exist('lw') );  lw=2;  end
   if ( ~exist('fs') );  fs=14;  end
   if ( ~exist('lgs') ); lgs=9;  end
end

if ( ~exist('fn') );  fn='Helvetica';  end
if ( ~exist('ps') );  ps='';  end
if ( ~exist('bx') );  
   if ( strcmp(get(gca,'Box'),'on') ); bx=1;  else;  bx=0; end;
end
if ( ~exist('gr') );  gr=0;  end

if ( exist('X') );  set(gca,'XLim',X);  end
if ( exist('Y') );  set(gca,'YLim',Y);  end
if ( exist('Z') );  set(gca,'ZLim',Z);  end

if (nl==0)
   set(get(gca,'XLabel'),'FontSize', fs); 
   set(get(gca,'YLabel'),'FontSize', fs); 
   set(get(gca,'ZLabel'),'FontSize', fs); 
   set(get(gca,'XLabel'),'FontName', fn); 
   set(get(gca,'YLabel'),'FontName', fn); 
   set(get(gca,'ZLabel'),'FontName', fn); 
   set(get(gca,'Title'),'String', ''); 
   set(gca,'FontName',fn); 
   set(gca,'FontSize',fs);
else
   set(gca,'XTickLabel',''); xlabel('');
   set(gca,'YTickLabel',''); ylabel('');
   set(gca,'ZTickLabel',''); zlabel('');
end

ch=get(gcf,'Children'); 
for j=1:size(ch,1) 
  if (strcmp(get(ch(j),'Tag'),'legend')) 
     set(ch(j), 'FontSize', lgs);
  end; 
end; 

set(gca,'LineWidth',lw);

for h=get(gca,'Children'); 
  if (strcmp(get(h,'Type'),'line')) 
     set(h,'LineWidth', lw); 
     set(h,'MarkerSize', ms); 
  elseif (strcmp(get(h,'Type'),'surface')) 
     set(h,'LineWidth', lw); 
  elseif (strcmp(get(h,'Type'),'text')); 
     set(h,'FontSize', fs); 
  end; 
end; 

if (bx==0); set(gca,'Box','Off'); else; set(gca,'Box','On'); end

if (gr==0); grid off; else; grid on; end

if (~isempty(ps))
   if (strcmp(sty,'journal') & strcmp(sh,'rectangle'))
      set(gcf,'PaperPosition',[0.25 2.5 3.4 2.55])
      set(gca, 'Position', [0.2 0.22 .62 .652])
   elseif (strcmp(sty,'journal') & strcmp(sh,'square'))
      if (nl==0)
	 set(gcf,'PaperPosition',[0.25 2.5 2.25 2.25])
         set(gca, 'Position', [0.2 0.22 .62 .652])
      else
	 set(gcf,'PaperPosition',[0.25 2.5 1.75 1.75])
         set(gca, 'Position', [0.05 0.05 .9 .9])
      end
   end
   eval(['print -depsc ' ps '.eps']);
end
