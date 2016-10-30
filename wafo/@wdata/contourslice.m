function varargout = contourslice(f,varargin)
%WDATA/CONTOURSLICE   Contours in slice planes.
%
% CALL [c,h] = contourslice(f,varargin)
%
%         f = wdata object
%  varargin = list of additional arguments, see contourslice for details.
%
% Example
% %[x y z v] = flow;
%  x = -2:.2:2; y = -2:.25:2; z = -2:.16:2;
%  [X,Y,Z] = meshgrid(x,y,z);
%  v = X .* exp(-X.^2 - Y.^2 - Z.^2);
%  wd = wdata(v,{x,y,z});
%  h=contourslice(wd,[1:9],[],[0], linspace(-8,2,10));
%  axis([0 10 -3 3 -3 3]); daspect([1 1 1])
%  camva(24); camproj perspective;
%  campos([-3 -15 5])
%
% See also contourslice, isosurface


Nf = numel(f);
if Nf>1 
  hold_state = ishold; % remember old hold state
  cfig = gcf;
  %c = cell(1,Nf);
  h = zeros(1,Nf);
  for ix=1:Nf,
    if hold_state
      newplot
    else
      figure(cfig-1+ix)
    end
    [h(ix)] = contourslice(f(ix),varargin{:});
  end  
  if nargout>0
    varargout{1} = h;
  end
else
  [varargout{1:nargout}] = contourslice(f.args{:},f.data,varargin{:});
  labelfig(f)
end