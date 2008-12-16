function varargout = pcolor(f,varargin)
%WDATA/PCOLOR   Pseudocolor (checkerboard) plot
%
% CALL  h = pcolor(f,varargin)
%
%         f = wdata object
%  varargin = list of additional arguments, see pcolor for details.
%
% Example
%  x = linspace(-3,3); y=x;
%  [X,Y] = meshgrid(x,y);
%  wd = wdata(peaks(X,Y),{x,y})
%  pcolor(wd)
%
%
% See also pcolor, contour, mesh

% History
% by pab 2007

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
    h(ix) = pcolor(f(ix),varargin{:});
  end  
  if nargout>0
    varargout{1} = h;
  end
else
  [varargout{1:nargout}] = pcolor(f.args{:},f.data,varargin{:});
  labelfig(f)
end