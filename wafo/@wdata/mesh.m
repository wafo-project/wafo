function varargout = mesh(f,varargin)
%WDATA/MESH   3-D mesh surface
%
% CALL Hout = mesh(f,varargin)
%
%         f = wdata object
%  varargin = list of additional arguments, see mesh for details.
%
% Example
%  x = linspace(-3,3); y=x;
%  [X,Y] = meshgrid(x,y);
%  wd = wdata(peaks(X,Y),{x,y})
%  mesh(wd)
% %  Colorize peaks with clown image.
%  c1 = load('clown')
%  mesh(wd,flipud(c1.X),...
%        'FaceColor','texturemap',...
%        'EdgeColor','none',...
%        'CDataMapping','direct')
% colormap(c1.map)
%
% See also mesh


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
    h(ix) = mesh(f(ix),varargin{:});
  end  
  if nargout>0
    varargout{1} = h;
  end
else
  [varargout{1:nargout}] = mesh(f.args{:},f.data,varargin{:});
  labelfig(f)
end