function h=cplot(x,y,varargin)
%CPLOT 	colored linear plot
%	CPLOT(X,Y,C) plots a line (Y vs. X) with the color specified by C.
%	If X is a matrix, one line per column is created (Y&C should match its size then)
%	Colors are interpolated according to current colormap & caxis.
%	CPLOT(X,Y) is equivalent to CPLOT(X,Y,Y).
%  	CPLOT(...,linespec) specifies the style of the line used (cf. PLOT) 
%  	all other parameters are passed to the PATCH object (e.g. MarkerSize, etc.)
%	H=CPLOT(...) returns handle(s) of PATCH object.
%  
%	Example: cplot(-5:.1:5,cos(-5:.1:5));
%	
%	Note: you can use ('edgecolor','interp'), but neither line thickness, 
%			nor line style will work (at least with Matlab 5.3Win)
%	Note: PATCH objects are much slower than LINEs. Be patient with larger plots(~1e6 points) 

% 5/21/2001 ashcherbina@ucsd.edu
arg=varargin;

if isempty(arg) | isstr(arg{1})
   c=y;
else
   c=arg{1};
   arg={arg{2:end}};
end

% get the linestyle
line=[];
marker=[];
if ~isempty(arg)
   [line,marker]=lstyle(arg{1});
   if (~isempty(line) | ~isempty(marker)) 
      if length(arg)>1, 
         arg={arg{2:end}};
      else
         arg={};
      end
   end
end   
if (isempty(line) & ~isempty(marker)) line='none';
elseif (~isempty(line) & isempty(marker)) marker='none';
elseif (isempty(line) & isempty(marker)) line='-';marker='none';end

if any(size(x)==1) 
   
   x=[x(:);nan];
   y=[y(:);nan];
   c=[c(:);nan];
else
   x=[x;x(1,:)*nan];
   y=[y;y(1,:)*nan];
   c=[c;c(1,:)*nan];
end

cax = newplot;

hh=patch(x,y,0);
set(hh,'cdata',c,'LineStyle',line,'Marker',marker);
set(hh,'facecolor','none','edgecolor','flat');

if (~isempty(arg))
   set(hh,arg{:});
end

box on;

if nargout>0, h=hh;end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [line,marker]=lstyle(s)
% decode the linestyles

line=[];marker=[];
markers='+o*.xsd^v><ph';
if length(s)>3, return;end; % just in case a property was passed

for k=1:length(markers);
   if ~isempty(findstr(s,markers(k)))
      marker=markers(k);
      break
   end
end
s=[s,'   '];
if ~isempty(findstr(s,'--'))
   line='--';
elseif ~isempty(findstr(s,'-.'))
   line='-.';
   if (length(findstr(s,'.'))<2 & marker=='.'), marker=[];end;
elseif ~isempty(findstr(s,':'))
   line=':';
elseif ~isempty(findstr(s,'-'))
   line='-'; 
end
