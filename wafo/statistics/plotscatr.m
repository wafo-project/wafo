function J = plotscatr(X,varargin)
%PLOTSCATR  Pairwise scatter plots.
%
% CALL:  plotscatr(X,label,sym)
%
%    X     = data
%    label = character array or a cellarray of strings to label the plot
%    sym   = plotsymbol    (default '.')
%
%	  The columns of X are plotted versus each other.  On the diagonal there is
%	  a normal probability plot for every variable.
%
%  Note: the order of label and sym is arbitrary
%
% Example
%   X = rand(30,3);
%   plotscatr(X,'*',{'wood','ply','none'}),shg
%
% See also plotmatrix

% Tested on matlab 5.x
% History
% revised pab 01.11.2000
%  - updated help header to WAFO style
%  - added N and correlation coefficient to the plot
%  - added label
%  - removed xtick and ytick from the axis
%       Copyright (c) Anders Holtsberg
%
% This program is free software; you can redistribute it and/or modify
%     it under the terms of the GNU Lesser General Public License as published by
%     the Free Software Foundation; either version 3 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU Lesser General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.



%error(nargchk(1,3,nargin))
narginchk(1,3)
plotsymbol = '.'; % default
labl =[];
[n,p] = size(X);

for ix=1:nargin-1,
   A = varargin{ix};
   if iscell(A),
      labl = char(A); % transform to a character array
   elseif size(A,1) == n 
      labl = A;
   else
      plotsymbol = A;
   end
end

clf
X = X - ones(n,1)*min(X);
X = X ./ (ones(n,1)*max(X));
X = X*0.8 + 0.1;


Z = zeros(p*p*n,2);
for i = 1:p
   for j = 1:p
      k = ((i-1)+(j-1)*p)*n+1;
      if i~=j,         
         Z(k:k+n-1, 1) = X(:,i) + i-1;
         Z(k:k+n-1, 2) = X(:,j) + p-j;
      else
         Z(k:k+n-1, 1) = sort(X(:,i)) + i-1;
         xx = ((1:n)-1/2)'/n;
         Z(k:k+n-1, 2) = erfinv(2*xx-1)/7 + 0.5 + p-j;
      end
   end
end

hold off;
H0 =plot(Z(:,1), Z(:,2), plotsymbol);
axis([0 p 0 p])
hold on;
for i = 0:p
   plot([0,p],[i,i],'-')
   plot([i,i],[0,p],'-')
end

rho = corr(X);% Put correlation coefficient in the plot
for i=1:p-1, 
   for j=i+1:p     
     figtext(i-1+0.01,p-j+0.01,['\rho = ' num2str(rho(i,j),2)],'data','data','left','bottom')
     figtext(j-0.01,p-i+0.01,['\rho = ' num2str(rho(i,j),2 )],'data','data','right','bottom')
  end
end
if ~isempty(labl)
   for j=1:p,
     set(gca,'DefaultTextRotation',90)
     figtext(-.3,p-j+0.5,labl(j,:),'data','data','left','middle')
     set(gca,'DefaultTextRotation',0)
     figtext(j-0.5,-.3,labl(j,:),'data','data','center','bottom')
   end
end

% Remove xtick and ytick
set(gca,'xtick',[],'ytick',[])
hold off
title(['N = ' num2str(n)])
if nargout>0
  J=H0;
end
