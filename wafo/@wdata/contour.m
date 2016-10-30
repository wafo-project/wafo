function varargout = contour(f,varargin)
%WDATA/CONTOUR   Contour plot
%
% CALL [c,h] = contour(f,varargin)
%
%         f = wdata object
%  varargin = list of additional arguments, see contour for details.
%
% Example
%  x = linspace(-3,3); y=x;
%  [X,Y] = meshgrid(x,y);
%  wd = wdata(peaks(X,Y),{x,y});
%  clf()
%  contour(wd)
%
%
% See also contour, mesh


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
    [c{ix},h(ix)] = contour(f(ix),varargin{:});
  end  
  if nargout>0
    varargout{1} = c;
  end
  if nargout>1
    varargout{2} = h;
  end
else
  if isempty(f.contourLevels) || nargin>1 && isnumeric(varargin{1})
    [varargout{1:nargout}] = contour(f.args{:},f.data,varargin{:});
  else
    cl = f.contourLevels;
    if length(cl)==1
      cl(2) = cl(1);
    end
    [varargout{1:nargout}] = contour(f.args{:},f.data,cl,varargin{:});
  end
  labelfig(f)
end