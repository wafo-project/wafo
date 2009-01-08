function ind2=dat2midind(x,p,ind,h)
% DAT2MIDIND Finds indices to midpoints between a min and Max and Max and min.
% 
%   CALL: ind = dat2midind(x,p,TP_index/wdef,h);
%	
%	ind = indices to midpoints.
%         x = two column data matrix with sampled times and values.
%         p = level between min and max. The default value p=0.5
%		corresponds to midpoint. (0 < p < 1).    
%  TP_index = indices to turningpoints.        
%	wdef= defines the type of wave. Possible options are
%		'mw' 'Mw' or 'none'. Default is 'mw'.
%		If wdef='none' all rainflow filtered min and max 
%		will be used, otherwise only the rainflow filtered 
%		min and max  which define a wave according to the 
%		wave definition will be used.
% 	h  = a threshold
%             if  h<0, then  tp=x; 
%             if  h=0, then  tp  is a sequence of turning points (default); 
%             if  h>0, then all rainflow cycles with height smaller than
%                           h  are removed.

% last modified by Per A Brodtkorb 06.07.00 , 01.09.98


xn=x;

[n m]= size(xn);
if n<m
  b=m;m=n;n=b; 
  xn=xn';
end

if n<2, 
  error('The vector must have more than 2 elements!')
end

%istime=1;

switch m
  case 1, x2=xn; %istime=0;
  case 2, x2=xn(:,2);% dimension OK!
  otherwise, error('Wrong dimension of input! dim must be 2xN, 1xN, Nx2 or Nx1 ')  
end

if (nargin<4) || isempty(h)
  h=0;
end

if (nargin<3) || isempty(ind)
  wdef='mw';
  [tp ind]=dat2tp(xn,h,wdef);
else
  if ischar(ind),
    wdef=ind;
    [tp ind]=dat2tp(xn,h,wdef);
  else
    wdef=[];
    if (abs(diff(sign(diff(x2(ind)))))~=2),
      error('Wrong input! TP_index is not a sequence of turningpoints')
    end
  end
end

if ((nargin<2) || isempty(p)),
  p=0.5;
end

if (p<=0) || (p>=1),
  error('p must be between 0 and 1.')
end

Ntp=length(ind);
Nw=floor((Ntp-1)/2);


ind2=zeros(2*Nw,1);

for i=1:Nw,
  tmp=findcross(x2((ind(2*i-1)):ind(2*i)), ((1-p)*x2(ind(2*i-1))+p*x2(ind(2*i))));
  ind2(2*i-1)=tmp(end); % choosing the one closest to the Max (or min)
  tmp = findcross(x2((ind(2*i)):ind(2*i+1)), ...
      (p*x2(ind(2*i))+(1-p)*x2(ind(2*i+1))));
  ind2(2*i)=tmp(1); % choosing the one closest to the Max (or min)
end

%if 2*Nw+1<Ntp,
%  tmp = findcross(x2((ind(Ntp-1)+1):ind(Ntp)),((1-p)*x2(ind(Ntp-1))+ p*x2(ind(Ntp))) );
%  ind2(Ntp-1)=tmp(end);
%end

ind2=ind(1:(2*Nw))+ind2;