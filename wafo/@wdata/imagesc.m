function varargout = imagesc(f,varargin)
%WDATA/IMAGESC  Scale data and display as image
%
% CALL [c,h] = imagesc(f,varargin)
%
%         f = wdata object
%  varargin = list of additional arguments, see imagesc for details.
%
% Example
%  x = linspace(-3,3); y=x;
%  [X,Y] = meshgrid(x,y);
%  wd = wdata(peaks(X,Y),{x,y})
%  imagesc(wd)
%
% See also imagesc, surf, contourf, contour, mesh


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
    h(ix) = imagesc(f(ix),varargin{:});
  end  
  if nargout>0
    varargout{1} = h;
  end
else
  [varargout{1:nargout}] = imagesc(f.args{:},f.data,varargin{:});
  labelfig(f)
end
