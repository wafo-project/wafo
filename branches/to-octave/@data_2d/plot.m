function H1=plot(self,f,varargin)
%PLOT Plot WDATA objects of DATA_2D type.
%
% CALL: H = plot(self,f,plotflag,sym,shading)
%
%  H = handle to the created object
%
%  plot a WDATA object with the fields
%      
%      f.data   = matrix
%      f.args   = vector of values 
%  
% optional fields:   
%      f.dataCI
%      f.labels        = cellarray of label strings          (n=1:2) 
%      f.title         = title string
%
%      plotflag =  1 contour plot (default)
%                 2 mesh
%                 3 surf
%                 4 waterfall
%                 5 pcolor
%                 6 contour3
%      sym      = plot symbol (default '-')
%      shading  = controls the color shading of SURFACE and PATCH
%                 objects.  'faceted' (default), flat or 'interp'
%
% Example
%  x = linspace(0,10).';
%  [X1,X2] = meshgrid(x,x);
%  f = wdata(peaks(X1,X2),{x,x});
%  plot(f.type,f)
% 
% See also wdata/plot

% Note: is only able to handle 1D plot i.e. ndim=1

%Tested on: Matlab 5.3, 5.2
%History:
% by pab March 2007
% - based on old pdfplot





hold_state = ishold; % remember old hold state
Nff=length(f);
if Nff>1
 error('WAFO:DATA_2D:PLOT','Illegal size of WDATA object')
end

%cax  = newplot; % axes
%cfig = get(cax,'Parent'); %cfig=gcf;
if iscell(f.args)
  dim = length(f.args);
else
  dim = 1;
end
if dim>length(size(squeeze(f.data)))
  fsiz=size(f.data);
  dim=length(fsiz)-sum(fsiz==1); % Number of non-singleton dimensions
end
if dim~=2
  error('WAFO:DATA_2D:PLOT','Wrong dimension of WDATA object')
end


[plotflag,sym,shad] = plotchk(varargin);


if plotflag==0, 
  return,
end


switch plotflag
  case {1,6,7,8,9},
    isPL=0;
    if ~isempty(f.contourLevels) % check if contour levels are submitted
      CL=f.contourLevels;
      isPL=~isempty(f.percentLevels); % levels defines quantile levels? 0=no 1=yes
    else
      CL = max(f.data(:))-range(f.data(:))*(1-[0.01 0.025 0.05 0.1 0.2 0.4 0.5 0.75]);
      if 0 % automatic levels by using contours
        c=contours(f.args{:},f.data.'); % calculate 8 levels
        if isempty(c)
          c=contours(f.args{:},f.data);%,7); % calculate levels
        end
        CL = clevels(c);
      end
    end
    if isPL,
      clvec=sort(f.percentLevels);
    else
      clvec=sort(CL);
    end
    if any(plotflag==[1 8 9])
      [cs hcs] = contour(f);%.args{:},f.data,CL,sym);
    else
      [cs hcs] = contour3(f);%.args{:},f.data,CL,sym);
    end
    if any(plotflag==[1,6])
      ncl=length(clvec);
      if ncl>12, 
        ncl=12; 
        disp('   Only the first 12 levels will be listed in table.'),
      end
      %axcl = cltext(clvec(1:ncl),isPL);  % print contour level text
      cltext(clvec(1:ncl),isPL);  % print contour level text
    elseif any(plotflag==[7 9])
      clabel(cs);
    else
      clabel(cs,hcs);
    end
    H = hcs;
  case 2,	H = mesh(f); %mesh(f.args{:},f.data); % meshz
  case 3,	H = surf(f); %surf(f.args{:},f.data);  %shading interp % flat, faceted       % surfc
  case 4,	H = waterfall(f);%.args{:},f.data);
  case 5, H = pcolor(f);%.args{:},f.data); %shading interp % flat, faceted
  case 10,
    [cs,hcs] = contourf(f);%.args{:},f.data); 
    clabel(cs,hcs); 
    fcolorbar(cs);
    H = hcs;
  otherwise, error('unknown option for plotflag')
end
if any(plotflag==(2:5))
   shading(shad);
end
labelfig(f)
if ~hold_state
  hold('off')
end
if (nargout>=1)
  H1=H;
end
end




function [plotflag,sym,shad] = plotchk(P)
%plotCHK Helper function for plot.
%
% CALL  [plotflag,sym,shad]=plotchk(P,dim) 
%
%   P = the cell array P of input arguments (between 0 and 6 elements)


% initialize output to default values
plotflag = 1;
sym='k-'; % Black dots is default
shad ='faceted';


Np=length(P);
try
  strix = cellfun('isclass',P,'char');
catch
  strix=zeros(1,Np);
  for ix=1:Np, % finding symbol strings 
    strix(ix)=ischar(P{ix});
  end
end
k=find(strix);
if any(k) % remove strings
  Nk=length(k);
  if Nk>3
    warning('WAFO:WDATA:PLOT','More than 3 strings are not allowed in ')
  end
  for ix=1:length(k)
    switch lower(P{k(ix)})
      case {'flat','faceted','interp'}  , shad   = P{k(ix)};
      otherwise % plotsymbol is given
        sym = P{k(ix)};
    end
  end
  Np=Np-length(k);
  P={P{find(~strix)}}; % remove strings from input
end


if (Np>0) && ~isempty(P{1})
  plotflag=P{1};
end


end