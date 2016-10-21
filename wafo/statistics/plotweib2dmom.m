function [ax11 ]=plotweib2dmom(V,H,varargin)
%PLOTWEIB2DMOM Plot conditional mean and standard deviation
% 
%   CALL: plotweib2dmom(x1,x2,phat,res,sym);
%
%      x1,x2 = data 
%      phat  = [A1 B1 A2 B2 C12], i.e., 2D weibull parameters.     
%      res   = resolution (default range(X2)/12)
%      sym   = {s1,s2,s3,s4} cell array of plot symbols for 
%              plotting empirical mean and standard deviation and 
%              theoretical mean and standard deviation, respectively.
%              (default {'bo','rx','b-','r--'})
%
% %PLOTWEIB2DMOM plots the  empirical conditional mean and standard
% deviation of X1 given X2 and optionally compares it with
% a 2D weibull distribution with parameters given in phat.
%
% Example
%  phat = {1 2 1.5 1 .8};
%  [y1,y2] = rndweib2d(phat{:},1000,1);
%  plotweib2dmom(y1,y2,phat{:});
%
% See also   momweib2d

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


% tested on matlab 5.2
% history:
% revised pab 03.11.2000
% changed var(tmp) to std(tmp)^2
% by Per A. Brodtkorb 23.11.98
  
%default values
sym = {'bo','rx','b-','r--'}; % default plot symbols for the empirical
                              % mean and std, and theoretical mean and
                              % std,respectivel
			      
error(nargchk(2,15,nargin))
Np = 5;
options = struct('condon',2,'resolution',[]); % default options
[params,options] = parsestatsinput(Np,options,varargin{:});
if numel(options)>1
  error('Multidimensional struct of distribution parameter not allowed!')
end

a1 =params{1};

res = options.resolution;

if isempty(res),   res=range(H(:))/12; end

grp=floor(H/res)+1; % dividing the data into groups 
Ngrp=max(grp);
Nmesh=40;
h1=linspace(res/2, (Ngrp-0.5)*res, Ngrp)';
%v1=linspace(eps,max(V),Nmesh)';
h2=linspace(0, max(H), Nmesh)';


m=h1;v=h1;
for ix=1:Ngrp,
  tmp=V(grp==ix);%find data in group number ix
  
  if length(tmp)>max(4,0),% if less than 4 observations in the group 
    m(ix)=mean(tmp);
    v(ix)=std(tmp).^2;
  else 
    m(ix)=NaN;
    v(ix)=NaN;
  end
end

ih = ishold;
plot(h1,m,sym{1},h1,sqrt(v),sym{2});  hold on
if ~isempty(a1), 
  [M1 V1]= momweib2d(params{:},'condon',2,'cvar',h2);
  plot(h2,M1,sym{3},h2,sqrt(V1),sym{4}),  
end
if ~ih, hold off,end

title('Conditional mean and standard deviation')
xlabel('x2')
if nargout>0
  ax11=gca;
end
