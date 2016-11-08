function plotcc(cc,ps)
%PLOTCC Plots a cycle count as a point process in the plane.
%
% CALL: plotcc(cc,psize)
%
%   cc    = a two column matrix with cycles.
%   psize = point size, (optional, default value is 12).
%
% See also  tp2mm, tp2rfc, tp2tc

% Tested  on Matlab  5.3
%
% History:
% Revised by PJ 26-Jul-2000
%   Now works when cc is 4 column matrix.
% Revised by PJ (P�r Johannesson) 18-May-2000
%   When input cc is a cell-array, 
%   then each cell i plotted in a subplot. 
% Revised by PJ (P�r Johannesson) 01-Nov-1999
%   updated for WAFO
% Copied from WAT Ver. 1.2

% Check input arguments

ni = nargin;
%no = nargout;
error(nargchk(1,2,ni));

if ni<2
  ps=12;
end

% If F is a cell-array, then plot each cell in a subplot
if iscell(cc)
  [N,M] = size(cc);
  for i=1:N
    for j=1:M
      subplot(N,M,(i-1)*M+j)
      plotcc(cc{i,j},ps);
    end
  end
  
elseif ~isempty(cc)
  
  m=max(max(cc(:,1:2)));
  n=min(min(cc(:,1:2)));
  
  border=max(abs(m),abs(n))*1.1;
  
  plot(cc(:,1),cc(:,2),'.','markersize',ps)
  axis([-border, border, -border, border]); axis('square')
  xlabel('min')
  ylabel('max')

end
