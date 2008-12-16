function varargout = surf(f,varargin)
%WDATA/SURF  3-D colored surface.
%
% CALL [c,h] = surf(f,varargin)
%
%         f = wdata object
%  varargin = list of additional arguments, see surf for details.
%
% Example
%  x = linspace(-3,3); y=x;
%  [X,Y] = meshgrid(x,y);
%  wd = wdata(peaks(X,Y),{x,y})
%  surf(wd)
%
% See also surf, contourf, contour, mesh

Nf = numel(f);
if Nf>1 
  hold_state = ishold; % remember old hold state
  cfig = gcf;
  h = zeros(1,Nf);
  for ix=1:Nf,
    if hold_state
      newplot
    else
      figure(cfig-1+ix)
    end
    h(ix) = surf(f(ix),varargin{:});
  end  
  if nargout>0
    varargout{1} = h;
  end
else
  [varargout{1:nargout}] = surf(f.args{:},f.data,varargin{:});
  labelfig(f)
end
