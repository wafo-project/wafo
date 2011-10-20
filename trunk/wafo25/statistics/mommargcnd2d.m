function [M,V] = mommargcnd2d(phat,condon,cvar,csm,lin)
%MOMMARGCND2D Mean and variance for the MARGCND2D distribution
%
%  CALL:  [M,V] = mommargcnd2d(phat,condon,cvar,csm,lin)
%    
%    M,V    = mean and variance, respectively
%    phat   = parameter structure array (see fitmargcnd2d)
%    condon = 0 returns marginal mean and variance for X1, X2 (default)
%             1 returns conditional mean and variance of X2 given X1 
%             2 returns conditional mean and variance of X1 given X2
%     cvar  = conditional variable, i.e.,x1 or x2 depending on condon.
%     csm   = smoothing vector (see margcnd2dsmfun) (default [1 1 1])
%     lin   = extrapolation vector (default [1 1 1])
%
% Example:
%  x1=linspace(0,10)';
%  phat.x={[x1,exp(-0.1*x1)] 2 };
%  phat.dist={'rayl','rayl'};
%  [M,V]=mommargcnd2d(phat,2,x1);
%  plot(x1,M,'r--',x1,sqrt(V),'k-')
%  title(' Conditional mean and standard deviation')
%  legend('E(x1|x2)','std(x1|x2)')
%  xlabel('x2')
%
% See also  fitmargcnd2d, margcnd2dsmfun

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


%tested on: matlab 5.2
% history:
%  by Per A. Brodtkorb 28.10.98

if (nargin< 4)||isempty(csm), 
  csm=[];
end
 
if (nargin< 5)||isempty(lin), 
  lin=[];
end
if nargin<3||isempty(cvar) ,
  cvar=[];
end
  
if (nargin <2) ||  isempty(condon), 
 condon=0;
end
if (condon~=0 )&&(nargin ==0), 
  error('Requires one input argument the levels to condition on.'); 
end

UDIST=lower(phat.dist{2});
CDIST=lower(phat.dist{1});
PH=phat.x{2};


switch UDIST(1:2)
  case 'ra',  [m2,v2]= momray(PH(1));
  case 'we' ,  [m2,v2]=momweib(PH(1),PH(2));
  case 'gu' ,  [m2,v2]=momgumb(PH(1),PH(2),0);
  case 'tg' ,  [m2,v2]=momgumb(PH(1),PH(2),1);
  case 'ga' ,  [m2,v2]=momgam(PH(1),PH(2));
  case 'lo' ,  [m2,v2]=momlognorm(PH(1),PH(2));
  otherwise, error('unknown distribution')
end 

switch condon, %marginal stats
  case 0,
    M=zeros(1,2); %initialize mean
    V=M;%initialize variance
    M(2)=m2; V(2,2)=v2;
    disp('Warning this option is not complete!')
    %    M(1)=m1; V(1,1)=v1;
    % v(1,2)=covar; M2,1)=covar;
    
 case 1 , % conditional stats given V
      %M=zeros(size(cvar)); %initialize mean
      %V=M;%initialize variance
      error('not implemented yet!')
%       switch UDIST(1:2) % this is not correct yet
%         case 'ra',  pdf1= momray(PH);
%         case 'we' ,  pdf1=momweib(PH(1),PH(2));
%         case 'gu' ,  pdf1=momgumb(PH(1),PH(2),0);
%         case 'tg' ,  pdf1=momgumb(PH(1),PH(2),1);
%         case 'ga' ,  pdf1=momgam(PH(1),PH(2));
%         case 'lo' ,  pdf1=momlognorm(PH(1),PH(2));
%         otherwise,
%           error('unknown distribution')
%       end
  case 2, % conditional stats given H
  
    if isempty(cvar),cvar=linspace(0 ,m2+3*sqrt(v2),30)'; end
    %M=zeros(size(cvar)); %initialize mean
    %V=M;%initialize variance
    %size(cvar)
    [Av , Bv, Cv]=margcnd2dsmfun(phat,cvar,csm,lin);
    switch CDIST(1:2)
      case 'ra', [M,V] =  momray(Av);
      case 'gu', [M,V] =  momgumb(Av,Bv,0);
      case 'tg', [M,V] =  momgumb(Av,Bv,1);
      case 'lo', [M,V] =  momlognorm(Av,Bv);
      case 'ga', [M,V] =  momgam(Av,Bv);
      case 'we', [M,V] =  momweib(Av,Bv);
      otherwise, error('Unknown distribution')
    end
    
  otherwise
    error('Unkown value for condon')
end
M=M+Cv;



