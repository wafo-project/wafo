function [y, trend] = detrendma(x,L)
%DETRENDMA Removes a trend from data using a moving average
%           of size 2*L+1.  If 2*L+1 > length(x) then the mean is removed
% 
% CALL:  [y, trend] = detrendma(x,L)
%
%    y     = x - trend
%    trend = moving average  which is removed  
%    x     = vector  of data or matrix of column vectors of data
%    L     = determines the size of the moving average window
%
% Example:
%  x = linspace(0,1,200)';
%  y = exp(x)+cos(5*2*pi*x)+1e-1*randn(size(x));
%  [y0, tr] = detrendma(y,20);
%  plot(x,y,x,y0,'r',x,exp(x),'k',x,tr,'m')
%
% See also  Reconstruct

% tested on : matlab 5.3
% revised pab 01.08.2001
% -added ; + nargchk + example + check on L
% - fixed a bug: y was always a column vector even if x was a row vector
% revised pab 13.01.2000
%  - made trend the same size as y
% By Per A. Brodtkorb  21.04.1999

error(nargchk(2,2,nargin))
if L<=0,error('L must be positive'),end
if L~=round(L), error('L must be an integer'),end

r=size(x,1);
if r==1,
  x=x(:);%make sure it is a column
end

[n, m]=size(x);
if n<2*L+1,% only able to remove the mean
  trend=mean(x);
  if m==1,
    y=x-trend;
  else
    y=x-trend(ones(n,1),:);
  end
  return
end

mn = mean(x(1:2*L+1,:));
y  = zeros(n,m);
y(1:L,:)=x(1:L,:)-mn(ones(L,1),:);

%if 1,%new call which is much faster
  ix      = (L+1):(n-L);
  trend   = cumsum([mn;(x(ix+L,:)-x(ix-L,:))/(2*L+1)],1);
  y(ix,:) = x(ix,:)-trend(2:end,:);
% else % old call slow
%   trend=zeros(n-2*L,m);
%   trend(1,:)=mn;
%   for ix=L+1:n-L,
%     mn=mn+ (x(ix+L,:)-x(ix-L,:))/(2*L+1);
%     y(ix,:)=x(ix,:)-mn;
%     trend(ix-L,:)=mn;
%   end
% end
mn2=trend(end,:);

if nargout>1
  trend=[mn(ones(L-1,1),:);trend;mn2(ones(L,1),:)];
  if r==1,  trend = trend.'; end
end

y(n-L+1:n,:)=x(n-L+1:n,:)-mn2(ones(L,1),:);

if r==1,  y = y.'; end
