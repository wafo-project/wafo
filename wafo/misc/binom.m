function bi = binom(n,k,logb)
% BINOM Calculates the binomial coefficient n!/((n-k)!*k!)
%
% CALL:  b = binom(n,k,logb);
%
%  b    = n!/((n-k)!*k!)
%  logb = if TRUE, output, b, returned as log(b)
%
% The binomial coefficient BINOM(N,K} gives the number of ways that K objects
% can be chosen from a collection of N distinct objects, regardless of order.
%
% Example:%  
%   b52  = binom(5,2)       % Should be 10.
%   bmax = binom(realmax,1) % Should be realmax
%   assert(b52==10,'binom(5,2) should equal 10')
%   assert(bmax==realmax,'binom(realmax,1) should equal realmax')
% 
% See also pdfhyge, pdfbin, nchoosek

% tested on: matlab 5.x
% History:
% by pab 17.11.98

% Copyright (C) 2000, 2007  Per A. Brodtkorb
% 
%  This file, BINOM.M, is part of WAFO.
% 
%     BINOM is free software; you can redistribute it and/or modify
%     it under the terms of the GNU Lesser General Public License as published by
%     the Free Software Foundation; either version 3 of the License, or
%     (at your option) any later version.
% 
%     BINOM is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU Lesser General Public License for more details.
% 
%     You should have received a copy of the GNU Lesser General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.
    
if nargin<3
  logb = false;
end

%   bi = nchoosek(n,k);
%
% Old call kept just in case
%logbi = gammaln(n+1)-gammaln(n-k+1)-gammaln(k+1);
% New call handle more  cases correct
logbi = localbinom(n,k);



if logb
  bi = logbi;
else
  if all(round(n(:))==n(:)) && all(round(k(:))==k(:))
   bi = round(exp(logbi));
  else
    bi=exp(logbi);
  end
  spcase1 = n>0  & (k==1 | (n-k) == 1);
  if any(spcase1(:))
    if ~isscalar(n), n = n(spcase); end
    bi(spcase1) = n;
  end
end

function  logbi = localbinom(n,k)
if any(abs(n(:)*eps)>1) %any(abs(k(:)*eps)>2) 
  warning('WAFO:BINOM','Result may not be exact.')
end

nmk = n-k;
km1 = k-1;
nmkp1 = nmk+1;
% spcase0 = (k-km1)==1;
% if any(spcase0(:))
%   nmkp2 = n-km1;  
%   nmkp1(spcase0) = nmkp2(spcase0);
% end
%  swap = nmk>k;
% if any(swap(:))
%   tmp = nmk(swap);
%   nmk(swap) = k(swap);
%   k(swap) = tmp;
% end

%logbi = stirlerr(n)-stirlerr(nmk+1)-stirlerr(k+1)+(n+0.5)*log(n)-(k+0.5)*log(k+1)-(nmk+0.5)*log(nmk+1)+2-0.5*log(2*pi);
%logbi = stirlerr(n)-stirlerr(nmk+1)-stirlerr(k+1)-...
%  (n+0.5).*log1p(-(k-1)./n)-(k+0.5).*log(k+1)+k.*log(nmk+1)+2-0.5*log(2*pi);
%logbi = stirlerr(n)-stirlerr(nmk+1)-stirlerr(k+1)-(nmk+0.5).*log1p(-(k-1)./n)-(k+0.5).*log(k+1)+k.*log(n)+2-0.5*log(2*pi);

% Note the order of summation is important! Summing the smallest numbers
% first and positive and negative numbers separately before they are
% subtracted from eachother.
logbi = (((stirlerr(n)+2)-min(nmk+0.5,realmax).*log1p(-km1./n))-k.*log((k+1)./n)) -(((stirlerr(nmkp1)+stirlerr(k+1))+0.5*log(2*pi))+0.5*log1p(k));


% Handle special cases
logbi(nmk<0 | k<0) = -inf;
logbi((n==0 & k==0) | k==0) = 0;


%!test assert(binom(5,2), 10)
%!test assert(binom(realmax*1e-8,1), realmax*1e-8)
