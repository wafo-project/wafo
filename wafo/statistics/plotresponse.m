function h = plotresponse(y,D)
%PLOTRESPONSE Cubic plot of responses
%
% CALL:  h = plotresponse(y,D);
%
%  y = responses 
%  D = Design matrix 
%
% Example
%   D = ffd(3);
%   y = [60 72 54 68 52 83 45 80]; % Responses to design D.
%   plotresponse(y,D)
%
% See also  ffd

%
%     This program is free software; you can redistribute it and/or modify
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


% History
% by pab 2001
% revised pab 2007
% renamed from cplot to plotresponse

%error(nargchk(1,2,nargin))
narginchk(1,2)
sz = size(y);

n = length(y);

if n == prod(sz),
  y = y(:);
  n = length(y);
else
  n = sz(1);
end

kmp = log2(n);
if nargin<2, D = ffd(kmp); end



[n1] = size(D,1);

if n~=n1, error('Something wrong'), end

%p = k-kmp;

Dmax = max(D(:));
Dmin = min(D(:));
range = Dmax-Dmin;
D = 2*(D-Dmin)/range-1;
%D(D<1)=-1;

r = 0.4;  % radius of circle
x2 = [-1+r 1-r ];
x1 = [-1 -1 ]; 
switch kmp
  case 2, % 2D square plot
      h1 = plot( x1,x2,'b',-x1,x2,'b',x2,x1,'b',x2,-x1,'b');
      axis([-1 1 -1 1]*1.5)
      set(gca,'xtick',-1:2:1,'ytick',-1:2:1,'Box','off')
      
      h2 = text(D(:,1),D(:,2),num2str(y(:)),...
	  'HorizontalAlignment','center','VerticalAlignment','middle');
      hold on,h3=circle2d(D(:,1),D(:,2),r,'b');hold off
     
  case 3, % 3D cubeplot
    h1 = plot3( x1,x1,x2,'b',x1,-x1,x2,'b',-x1,x1,x2,'b',-x1,-x1,x2,'b');
    hold on,
    h11 = plot3(x1,x2,x1,'b',x1,x2,-x1,'b',-x1,x2,x1,'b',-x1,x2,-x1,'b');
    h111 = plot3(x2,x1,x1,'b',x2,x1,-x1,'b',x2,-x1,x1,'b',x2,-x1,-x1,'b');
    h1 = [h1;h11;h111];
    axis([-1 1 -1 1 -1 1]*1.5)
    set(gca,'xtick',-1:2:1,'ytick',-1:2:1,'ztick',-1:2:1,'Box','off')
    
    h2 = text(D(:,1),D(:,2),D(:,3),num2str(y(:)),...
	'HorizontalAlignment','center','VerticalAlignment','middle');
    if 1,
      h3=circle3d(D(:,1),D(:,2),D(:,3),r,'b');
    else
      h3 =[];
    end
    hold off
    zlabel('C')  
  case 4, % 2D  cubeplot
     h1 = plot( x1,x2,'b',-x1,x2,'b',x2,x1,'b',x2,-x1,'b');
     x0 = 1; hold on
     h11 = plot(x1+x0,x2+x0,'b',-x1+x0,x2+x0,'b',x2+x0,x1+x0,'b',x2+x0,-x1+x0,'b');
     x3 = [-1+ r/sqrt(2)  -r/sqrt(2) ];
     h111 = plot(x3,x3,'b',2*x0+x3,x3,'b',x3,2*x0+x3,'b',2*x0+x3,2*x0+x3,'b');
     
      axis([-1 2 -1 2]*1.5)
      set(gca,'xtick',-1:2:1,'ytick',-1:2:1,'Box','off')
      ind1 = 1:4;
      ind2 = 5:8;
     
      h2 = text(D(ind1,1),D(ind1,2),num2str(y(ind1)),...
	  'HorizontalAlignment','center','VerticalAlignment','middle');
      h21 = text(D(ind2,1)+x0,D(ind2,2)+x0,num2str(y(ind2)),...
	  'HorizontalAlignment','center','VerticalAlignment','middle');
      h3=circle2d(D(ind1,1),D(ind1,2),r,'b');
      h31=circle2d(x0+D(ind2,1),x0+D(ind2,2),r,'b');
      hold off
      h3 = [h3(:);h31(:)];
      h2 = [h2(:);h21(:)];
      h1 = [h1(:);h11(:);h111(:)];
end
xlabel('A') 
ylabel('B')
axis square
if nargout>1,
  h=[h1(:);h2(:);h3(:)];
end

function h = circle2d(x,y,r,varargin)
  %error(nargchk(2,inf,nargin))
  narginchk(2,inf)
  if nargin<3,r=1;end
  [csize, x,y,r]=comnsize(x,y,r);
  
  n = length(x(:));
  h = zeros(n,1);
  th = linspace(0,2*pi,50)';
  for ix=1:n,
    h(ix) = plot(r(ix)*cos(th)+x(ix),r(ix)*sin(th)+y(ix),varargin{:});
  end
  
return



function h = circle3d(x,y,z,r,varargin)
  %error(nargchk(3,inf,nargin))
  narginchk(3,inf)
  if nargin<4,r=1;end
  [csize, x,y,z,r]=comnsize(x,y,z,r);
  
  n = length(x(:));
  h = zeros(n,1);
  m = 50;
  th = linspace(0,2*pi,m)';
  for ix=1:n,
    h(ix) = plot3(r(ix)*cos(th)+x(ix),r(ix)*sin(th)+y(ix),repmat(z(ix),m,1));
  end
return



