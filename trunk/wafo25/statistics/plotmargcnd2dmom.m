function [ax11, h11, h22  ]=plotmargcnd2dmom(V,H,varargin) 
%PLOTMARGCND2DMOM Plot conditional mean and standard deviation
% 
%   CALL: plotmargcnd2dmom(x1,x2,phat,condon,res,csm,lin,sym);
%
%      x1,x2 = data
%      phat  =  parameter structure array (see fitmargcnd2d)
%    condon  = 1 mean and standard deviation of X2 given X1
%              2 mean and standard deviation of X1 given X2 (default)
%      res   = resolution (default range(x2)/12)
%      csm   = smoothing vector     (default [1 1 1])
%      lin   = extrapolation vector (default [1 1 1])
%      sym   = {s1,s2,s3,s4} cell array of plot symbols for 
%              plotting empirical mean and standard deviation and 
%              theoretical mean and standard deviation, respectively.
%              (default {'bo','rx','b-','r--'})
%
% PLOTMARGCND2DMOM plots the  empirical conditional mean and standard
% deviation of X1 given X2 or X2 given X1 and optionally compares it with
% theoretical quantities.
%
% NOTE:  SYM can be given anywhere after X1 and X2.
%
% Example:
%  x1=linspace(0,10)';
%  phat.x={[x1,exp(-0.1*x1)] 2 };
%  phat.dist={'rayl','rayl'};
%  [y1,y2] = rndmargcnd2d(1000,phat);
%  plotmargcnd2dmom(y1,y2,phat,2);
%
% See also  mommargcnd2d
 
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
% revised pab 03.11.2000
% - replaced call to var with std(tmp)^2
% revised pab 22.06.2000
%  - added varargin
%  - fixed some bugs for condon ==1
%  - added sym 
%  - fixed a bug for res.
%  - updated header for csm and lin
% by Per A. Brodtkorb 23.11.98


if  nargin<2,  error('Two input arguments required'), end 

% default values
%~~~~~~~~~~~~~~~
sym = {'bo','rx','b-','r--'}; % default plot symbols for the empirical
                              % mean and std, and theoretical mean and std,respectively
phat =[];
lin    = [1 1 1]; 
csm    = [1 1 1];
res    = [];
condon = 2; 


if 1,%new call
  P  = varargin;
  Np = length(P);
  if Np>0
    strix = zeros(1,Np);
    cellix = strix;
    for ix=1:Np, % finding symbol strings or cell array of symbol strings
      strix(ix)  = ischar(P{ix});
      cellix(ix) = iscell(P{ix});
    end
    k  = find(strix);
    k1 = find(cellix);
    if any(k)
      Nk = length(k);
      if Nk>4,  warning('WAFO:PLOTMARGCND2DMOM','More than 4 strings are not allowed'),    end
      iy = 1;
      for ix=k      
	sym{iy} = P{ix};
	iy=iy+1;
      end
      Np = Np-Nk;
      P  = {P{find(~strix)}}; % remove strings from input
    elseif any(k1) % cell array of strings
      tmp = P{k1};
      Nk = length(tmp);
      if Nk>4,  warning('WAFO:PLOTMARGCND2DMOM','More than 4 strings are not allowed'),    end
      iy = 1;
      for ix=1:min(Nk,4)
	if ~isempty(tmp{ix}) && ischar(tmp{ix}), sym{ix}=tmp{ix};end
      end
      Np = Np-1;
      P  = {P{find(~cellix)}}; % remove cell array of strings from input
    end
    if Np>0 && ~isempty(P{1}), phat   = P{1};end
    if Np>1 && ~isempty(P{2}), condon = P{2};end
    if Np>2 && ~isempty(P{3}), res    = P{3};end
    if Np>3 && ~isempty(P{4}), csm    = P{4};end
    if Np>4 && ~isempty(P{5}), lin    = P{5};end
  end
else % old call
  %function [ax11, h11, h22  ]=mommargcnd2dplot(V,H,phat,condon ,res,csm,lin)
  if (nargin<7)||isempty(lin),    lin    = [1 1 1];  end
  if (nargin<6)||isempty(csm),    csm    = [1 1 1];  end
  if (nargin<5)||isempty(res),    res    = range(H(:))/12; end
  if (nargin<4)||isempty(condon), condon = 2; end
end

if isempty(res)
  if condon==1,    res = range(V(:))/12;  else    res = range(H(:))/12;  end
end


if condon==2,
  grp=floor(H/res)+1; % dividing the data into groups 
  xmax = max(H);
  Xc = V;         
else
  grp=floor(V/res)+1; % dividing the data into groups 
  xmax = max(V);
  Xc = H;
end

Ngrp  = max(grp);
Nmesh = 40;
h1    = linspace(res/2, (Ngrp-0.5)*res, Ngrp)';



%cvar=linspace(eps,max(V),Nmesh)';
if ~isempty(phat), % theoretical distribution not given
  cvar  = linspace(0, xmax, Nmesh)';
  [M1 V1] = mommargcnd2d(phat,condon,cvar,csm,lin);
end

m=h1;v=h1;
for ix=1:Ngrp,
  tmp = Xc(grp==ix);%find data in group number ix
  
  if length(tmp)>max(4,0),% if less than 4 observations in the group 
    m(ix)=mean(tmp); % mean of data in group ix
    v(ix)=std(tmp).^2;  % variance of data in group ix
  else 
    m(ix)=NaN;
    v(ix)=NaN;
  end
end

if 0,
  [ax1 h11 h22]=plotyy(h1,m,h1,sqrt(v));
  set(ax1(1),'nextplot','add' );
  set(h11, 'LineStyle' , 'x')
  axes(ax1(1));
  plot(cvar,M1,'-'),
  set(ax1(1),'nextplot','replace' );
  
  set(ax1(2),'nextplot','add' );
  set(h22, 'LineStyle' , 'o')
  axes(ax1(2));
  %axis([0 inf 0 max(v)])
  
  plot(cvar,sqrt(V1),'-'), 
  set(ax1(2),'nextplot','replace' );
else
  ih = ishold;
  plot(h1,m,sym{1},h1,sqrt(v),sym{2});  hold on
  
  if ~isempty(phat),    plot(cvar,M1,sym{3},cvar,sqrt(V1),sym{4}),  end
  
  if ~ih, hold off,end
end

xlabel(['x' num2str(condon)])

title('Conditional mean and standard deviation')
if nargout>0
  ax11=ax1;
end
