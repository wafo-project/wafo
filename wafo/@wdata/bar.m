function varargout = bar(f,varargin)
%WDATA/BAR  Bar graph.
%
% CALL Hout = bar(f,varargin)
%
%         f = wdata object
%  varargin = list of additional arguments, see bar for details.
%
% Example
%  x = linspace(-3,3); 
%  y = humps(x);
%  wd = wdata(y,x);
%  bar(wd)
%
% See also bar wdata/plot wdata/stairs wdata/stem wdata/area



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
    h(ix) = bar(f(ix),varargin{:});
  end  
  if nargout>0
    varargout{1} = h;
  end
else
  [varargout{1:nargout}] = bar(f.args,f.data,varargin{:});
  labelfig(f)
end