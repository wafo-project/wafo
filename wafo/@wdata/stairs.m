function varargout = stairs(f,varargin)
%WDATA/STAIRS   Stairstep plot
%
% CALL Hout = stairs(f,varargin)
%
%         f = wdata object
%  varargin = list of additional arguments, see stairs for details.
%
% Example
%  x = linspace(-3,3); 
%  y = humps(x);
%  wd = wdata(y,x);
%  stairs(wd)
%
% See also stairs wdata/plot wdata/bar wdata/stem wdata/area


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
    h(ix) = stairs(f(ix),varargin{:});
  end  
  if nargout>0
    varargout{1} = h;
  end
else
  [varargout{1:nargout}] = stairs(f.args,f.data,varargin{:});
  labelfig(f)
end