function varargout = slice(f,varargin)
%WDATA/SLICE   Volumetric slice plot.
%
% CALL [c,h] = slice(f,varargin)
%
%         f = wdata object
%  varargin = list of additional arguments, see slice for details.
%
% Example
%  [x y z v] = flow;
%  wd = wdata(v,{x,y,z});
%  h = slice(wd,[1:2.5:9],[],[0]);
%  axis([0 10 -3 3 -3 3]); daspect([1 1 1])
%  camva(24); camproj perspective;
%  campos([-3 -15 5])
%
% See also slice, contourslice, isosurface


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
    [h(ix)] = slice(f(ix),varargin{:});
  end  
  if nargout>0
    varargout{1} = h;
  end
else
  [varargout{1:nargout}] = slice(f.args{:},f.data,varargin{:});
  labelfig(f)
end