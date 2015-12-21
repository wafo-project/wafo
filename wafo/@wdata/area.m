function varargout = area(f,varargin)
%WDATA/AREA  Discrete sequence or "stem" plot.
%
% CALL Hout = area(f,varargin)
%
%         f = wdata object
%  varargin = list of additional arguments, see area for details.
%
% Example
%  x = linspace(-3,3); 
%  y = humps(x);
%  wd = wdata(y,x);
%  area(wd)
%
% See also area wdata/plot wdata/bar wdata/stairs wdata/stem



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
    h(ix) = area(f(ix),varargin{:});
  end  
  if nargout>0
    varargout{1} = h;
  end
else
  [varargout{1:nargout}] = area(f.args,f.data,varargin{:});
  labelfig(f)
end