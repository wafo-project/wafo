function [phat, cov, pci]=fitweib3(data1, plotflag)
%FITWEIB3 Parameter estimates for 3 parameter Weibull data.
%
% CALL:  [phat, cov] = fitweib3(data, plotflag)
%
%     phat = [a,b,c] = the Least Squares estimates of the  
%            parameters of the 3 parameter  Weibull distribution
%            (see w3weibcdf) given the data.
%     cov  = asymptotic covariance matrix of estimates
%     data = data vector
% plotflag = 0, do not plot
%          > 0, plot the empiricial distribution function and the
%               estimated cdf (see edf for options)(default)
% 
% Example:
%   R=rndweib(1,3,1,200)+3;
%   [phat, cov] = fitweib3(R);
%   plotfitsumry(phat)
%
%   close all;
%
% See also  cdfweib

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


% Reference: Cohen & Whittle, (1988) "Parameter Estimation in Reliability
% and Life Span Models", p. 25 ff, Marcel Dekker.

%Tested on: matlab  5.3
% History:
% by pab januar 2006  

error(nargchk(1,2,nargin))
if nargin<2||isempty(plotflag),  plotflag=1; end

data  = data1(:);                            % make sure it is a vector
N     = length(data);

sd = sort(data);

F = (0.5:1:(N - 0.5))'/N;
if N>10000
  Fi = linspace(F(1),F(end-5),10000).';
  %Fi = fliplr(logspace(log10(F(end)),log10(F(1)),10000)).';
  sd =interp1(F,sd,Fi,'linear');
  F  = Fi;
end
c0 =  min(sd)-0.01*std(sd);

phat0 = fitweib(sd-c0,0);


%g0=inline('mean((-log(1-F)+((data+abs(x(3)))./x(1)).^x(2)-abs(x(3)/x(1)).^x(2)).^2 )','x','data','F');


monitor = false;

def=3; %  What is def? See 'w3weibfun' 

phat = fminsearch(@w3weibfun,[phat0,c0],optimset,sd,F,def,monitor);
%  phat =  fmins('w3weibfun',[phat0,c0],[],[],sd,F,def,monitor);



if nargout>1,
  [L,cov]= loglike(phat,data,'pdfweib');
end
if nargout>2,
  var=diag(cov)';
  alpha2=ones(1,3)*0.05/2;
  pci = invnorm([alpha2;1-alpha2],[phat;phat],[var;var]);
end


if plotflag 
  edf(sd,[sd, cdfweib(sd-phat(3),phat(1),phat(2))],plotflag)
  title(deblank('Empirical and 3 parameter Weibull estimated cdf'))
end



function y=w3weibfun(phat,x,F,def,monitor)
% W3WEIBFUN Is an internal routine for fitweib3
%

b =2; c1 = 0;

a = phat(1);
Np = length(phat);
if Np>1, b = phat(2);end
if Np>2, c1 = phat(3); end

c = -abs(c1);
N = length(F);
xmin = min(x(:));

penalty = abs(N*c*(xmin<c1));

%monitor = logical(1);
switch def
  case 1, % fit to sqrt(-log(1-F))
  if monitor
    plot(x,sqrt(-log(1-F)),....
	x,sqrt(((x+c)./a).^b)); drawnow
  end
  
  y=mean((-sqrt(-log(1-F))+...
      sqrt(((x+c)./a).^b).^2))*(1+penalty) + penalty;
case 2, % fit to (-log(1-F))
  if monitor
    plot(x,(-log(1-F)),...
	x,(((x+c)./a).^b)); drawnow
  end
  
  y=mean((-(-log(1-F))+...
      (abs((x+c)./a).^b)).^2)*(1+penalty) + penalty;

case   3, % fit to (-log(1-F)).^(1/b)
  if monitor
    plot(x,(-log(1-F)).^(1/b),x,...
	(((x+c)./a))); drawnow
  end
  y=mean((-(-log(1-F)).^(1/b)+...
      (((x+c)./a))).^2)*(1+penalty)+penalty;
case 4,  % fit to (-log(1-F)+(c/a).^b).^(1/b)
  if monitor
    plot(x,(-log(1-F)).^(1/b),x,(x+c)./a); drawnow
  end
  y=mean((-(-log(1-F)).^(1/b)+(x+c)./a).^2)*(1+penalty)+penalty;

end

if monitor
  disp(['err = ' num2str(y,10)   ' a b c = ' num2str([a,b,c1],4) ])
  shg
end

