function [V,H] = rndmargcnd2d(N,phat,csm,lin)
%RNDMARGCND2D Random points from a MARGCND2D distribution 
%  
%  CALL:  [x1,x2] = rndmargcnd2d(N,phat,[csma,csmb,csmc],[lina,linb,linc]);
%
%    X1,X2   = N random points in R^2.
%    N       = number of points generated
%    phat    = parameter structure array (see fitmargcnd2d)
%    csma..c =  vector of internal smoothing parameters (default [1 1 1])
%               0 -> LS-straight line
%               1 -> cubic spline interpolant
%   lina..c  = vector defining the extrapolation of parameter A,B and C, respectively 
%              0 No linear extrapolation outside the range of data
%              1 Linear extrapolation outside the range of data (default)
%
% Example: Random points from a 2D Rayleigh distribution
%   x1=linspace(0,10)';
%   phat.x={[x1,exp(-0.1*x1)] 2 };
%   phat.dist={'rayl','rayl'};
%   [y1,y2] = rndmargcnd2d(1000,phat);
%   f = pdfmargcnd2d2(x1,x1,phat);
%   pdfplot(f), hold on
%   plot(y1,y2,'.'), hold off
%   
%
% See also  fitmargcnd2d , pdfmargcnd2d

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


% tested on: matlab 5.2
% history:
% by Per A. Brodtkorb 28.10.98

error(nargchk(2,4,nargin))

if (nargin< 3)||isempty(csm), 
  csm=[];
end
if (nargin< 4)||isempty(lin), 
  lin=[];
end
UDIST=phat.dist{2};
CDIST=phat.dist{1};

PH=phat.x{2};

  
% H is distributed
switch UDIST(1:2),
  case 'ra', H=rndray(PH(ones(N,1),:));
  case 'we', H=rndweib(PH(ones(N,1),1) , PH(ones(N,1),2));
  case 'tg', H=rndgumb(PH(ones(N,1),1) , PH(ones(N,1),2),1);
  case 'gu', H=rndgumb(PH(ones(N,1),1) , PH(ones(N,1),2),0);
  case 'lo', H=rndlognorm(PH(ones(N,1),1) , PH(ones(N,1),2));
  case 'ga', H=rndgam(PH(ones(N,1),1) , PH(ones(N,1),2));	
end
 

[Av , Bv, Cv]=margcnd2dsmfun(phat,H,csm,lin); %parameters of V given H 

% V conditioned on H  is distributed 
switch CDIST(1:2)
  case 'ra', V = rndray(Av)+Cv;
  case 'gu', V = rndgumb(Av,Bv,0)+Cv;% tGumbel
  case 'tg', V = rndgumb(Av,Bv,1)+Cv;% truncated  Gumbel
  case 'lo', V = rndlognorm(Av,Bv)+Cv;
  case 'ga', V = rndgam(Av,Bv)+Cv;	
  case 'we', V = rndweib(Av,Bv)+Cv;
  otherwise, error('Unknown distribution') 
end


