function varargout = contourf(f,varargin)
%WDATA/CONTOURF  Filled Contour plot
%
% CALL [c,h] = contourf(f,varargin)
%
%         f = wdata object
%  varargin = list of additional arguments, see contourf for details.
%
% Example
%  x = linspace(-3,3); y=x;
%  [X,Y] = meshgrid(x,y);
%  wd = wdata(peaks(X,Y),{x,y});
%  clf()
%  contourf(wd)
%
%
% See also contourf, contour, mesh


Nf = numel(f);
if Nf>1 
  hold_state = ishold; % remember old hold state
  cfig = gcf;
  c = cell(1,Nf);
  h = zeros(1,Nf);
  for ix=1:Nf,
    if hold_state
      newplot
    else
      figure(cfig-1+ix)
    end
    [c{ix},h(ix)] = contourf(f(ix),varargin{:});
  end  
  if nargout>0
    varargout{1} = c;
  end
  if nargout>1
    varargout{2} = h;
  end
else
  if isempty(f.contourLevels) || nargin>1 && isnumeric(varargin{1})
    [varargout{1:nargout}] = contourf(f.args{:},f.data,varargin{:});
  else
    cl = f.contourLevels;
    if length(cl)==1
      cl(2) = cl(1);
    end
    zmin = min(f.data(:));
    if min(cl)>zmin
      cl(end+1) = zmin;
      cl = sort(cl);
    end
    [varargout{1:nargout}] = contourf(f.args{:},f.data,cl,varargin{:});
  end
  labelfig(f)
end