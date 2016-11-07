function plotnorm(x,nr1,nr2,df)
%PLOTNORM Plot data on a Normal distribution paper
%
% CALL:  plotnorm(x,nr1,nr2,df)
%
%   x      = a vector with data, or a matrix with data from
%            several groups, one group per column.
%   nr1    = plot parameter;
%            if 1, estimated lines will be drawn (default).
%   nr2    = plot parameter;
%            if 1, percentiles will be drawn (default).
%   df     = degrees of freedom in the t-distribution. 
%            df = inf    : normal distribution paper will be used (default),
%            0 < df < inf: approximation of t-distribution will be used.
%
% Example:
%   R=rndnorm(0,1,100,2);
%   plotnorm(R),shg
%
% See also prbt, prbnorm, fitnorm, fitt

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


% Tested on: Matlab 5.3
% History: 
% revised pab 24.10.2000
%  - removed tdistr
%  - added checks on df and ishold as well as nargchk
%  - added example
% Revised by jr 22.12.1999

error(nargchk(1,4,nargin))
if nargin<4||isempty(df),  df=inf;end
if (df<=0 || (df~=round(df) && df<inf))
  error('df must be postive integer')
end
tdistr=(df<inf);

if (nargin<3)||isempty(nr2),  nr2=1;end
if (nargin<2)||isempty(nr1),  nr1=1;end

ih = ishold;

if size(x,1)==1
  x=x';
end

[n,m]=size(x);
x=sort(x);
X=((1:n)'-1/2)/n;
Y=sqrt(2)*erfinv(2*X-1);

if tdistr
  g1=1/4*(Y.^3+Y);
  g2=1/96*(5*Y.^5+16*Y.^3+3*Y);
  g3=1/384*(3*Y.^7+19*Y.^5+17*Y.^3-15*Y);
  g4=1/92160*(79*Y.^9+776*Y.^7+1482*Y.^5-1920*Y.^3-945*Y);
  Z=Y+g1./df+g2./df^2+g3./df^3+g4./df^4;
  Y=Z;
end

linregx=Y;
SXX=sum((linregx-mean(linregx)).^2);
Phi1=0.5*(1+erf(1/sqrt(2)));
lambda975=sqrt(2)*erfinv(2*(.975)-1);

levels=[ .5 .7 .9 .95 .98 .99 .995 .999 .9999];
levels=[1-fliplr(levels(2:9)) levels];
lev=sqrt(2)*erfinv(2*levels-1); 

data=zeros(2,m);
mx=zeros(1,m);
sx=zeros(1,m);
for i=1:m
  linregy=x(:,i);
  SXY=(linregx-mean(linregx))'*(linregy-mean(linregy));
  b=SXY/SXX;
  a=mean(linregy)-b*mean(linregx);
  mx(i)=a;
  sx(i)=b;
  data(1,i)=mx(i);
  data(2,i)=sx(i)^2;
  plot(x(:,i),Y,'b.','markersize',12);
  hold on
  if nr1
    plot([x(1,i) x(n,i)],[(x(1,i)-mx(i))/sx(i) (x(n,i)-mx(i))/sx(i)],'r--')
  end
end

span=max(max(x))-min(min(x));
xx1=min(min(x))-0.1*span;
xx2=max(max(x))+0.1*span;

axis([xx1 xx2 -4 4])
if ~tdistr
  title('Normal Probability Plot')
  ylabel('Quantiles of standard normal')
else
  title('Student probability t-plot')
  ylabel('Quantiles of Student t')
end

if nr2
  ax=axis;
  plot([ax(1) ax(2)],[lev; lev],'k');
  for l=1:length(levels)
    h=figtext(1.01,lev(l),[num2str(levels(l)*100) '%'],'normalized','data');
    set(h,'FontSize',10,'FontAngle','Italic')
  end
end

wafostamp;
if ~ih, hold off,end




