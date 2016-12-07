function  X = rndbin(varargin)
%RNDBIN  Random numbers from the binomial distribution
%
% CALL X = rndbin(n,p,sz)
%      X = rndbin(phat,sz)
%
%        x = number of successes          (0<=x<=n)
%        n = number of independent trials
%        p = probability of succes in a given trial. (0<=p<=1)
%     phat = Distribution parameter struct
%            as returned from FITBIN. 
%
% Example:
%   par = {10,0.2}
%   X = rndbin(par{:},1000,1);
%   moments = [mean(X) var(X),skew(X),kurt(X)];   % Estimated moments
%   [mom{1:4}] = mombin(par{:}); % True mean and variance
%
% See also  pdfbin, cdfbin, invbin, fitbin, mombin


% Copyright (C) 2007 Per A. Brodtkorb
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


% Reference
% Luc Devroye (1986), 
% Non-Uniform Random Variate Generation, New York: Springer-Verlag, 1986. 
% Chapter X


Np = 2;
options = []; % default options
[params,options,rndsize] = parsestatsinput(Np,options,varargin{:});
% if numel(options)>1
%   error('Multidimensional struct of distribution parameter not allowed!')
% end

[n,p] = deal(params{:});

if isempty(rndsize)
  csize = comnsize(n,p);
else
  csize = comnsize(n,p,zeros(rndsize{:}));
end
if any(isnan(csize))
    error('n and p must be of common size or scalar.');
end



n(~isfinite(n) | (n ~= round (n))) = nan;
p(p<0 | 1<p) = nan;
[pmean,pvar] = mombin(n,p);

r1 = pmean>1000 | pvar==0 |~isfinite(pmean);
if isscalar(pmean)
  if any(r1)    
    % np> large
    X = round( pmean + sqrt(pvar)*randn( csize ) );
  else 
    X = zeros(csize);
    if p<=0.5
      pp =floor(n*p)/n;
      X(:) = rndbinreject(n,pp,prod(csize));
      Z = rndbinwait(n-X,(p-pp)/(1-pp),csize);
      X = X + Z;
    else
      q = 1-p;
      qq = floor(q*n)/n;
      X(:) = n-rndbinreject(n,qq,prod(csize));
      Z = rndbinwait(X,(q-qq)/(1-qq),(csize));
      X  = X-Z;
    end
  end
else
  if isscalar(p), p = p(ones(csize));end
  if isscalar(n), n = n(ones(csize));end
  X = zeros(csize);
  if any(r1(:))
    % np> large
    rsz = size(pmean(r1));
    X(r1) = round( pmean(r1) + sqrt(pvar(r1)).*randn( rsz) );
  end
  r2 = p<=0.5 & ~r1;
  if any(r2(:))
    nr2 = n(r2);
    pp =floor(pmean(r2))./nr2;
    
    xr2 = zeros(size(nr2));
    ppg0 = pp>0;
    xr2(ppg0) = rndbinrejectv(nr2(ppg0),pp(ppg0));
    Z = rndbinwait(nr2-xr2,(p(r2)-pp)./(1-pp),size(xr2));
    X(r2) = xr2 + Z;
  end
  r3 = 0.5<p & ~r1;
  if any(r3(:))
    q = 1-p(r3);
    nr3 = n(r3);
    qq = floor(q.*nr3)./nr3;
    qqg0 = qq>0;
    xr3 = (nr3);
    xr3(qqg0) = nr3(qqg0)-rndbinrejectv(nr3(qqg0),qq(qqg0));
    Z = rndbinwait(xr3,(q-qq)./(1-qq),size(xr3));
    X(r3)  = xr3-Z;
  end
end
return



function X = rndbinwait(n,p,csize)
%Random Binomial variables using the Waiting time algorithm
%n = n(:).';
X=-ones(csize);
Sum1  = zeros(csize);
left = Sum1;
left(:) = 1:numel(X);
Q = log1p(-p);
while any(left(:))
  G =  ceil(log(rand(size(left)))./Q); %Generate Geometric(p) random variate G.
  Sum1 = Sum1 + G;
  X(left)=X(left)+1;
  accept = Sum1>n;
  if any(accept)
    left(accept) = [];
    Sum1(accept) = [];
    if ~isscalar(Q)
      Q(accept) = [];
    end
    if ~isscalar(n)
      n(accept) = [];
    end
  end 
end
% 



function X = rndbinreject(n,p,N)
% Rejection sampling for the binomial distribution for p<=0.5
% Assumption n*p is an integer
if p==0;
  X = zeros(1,N);
  return
end
[pmean, pvar] = mombin(n,p);
[qmean] = mombin(n,1-p);

d1 = floor(max(1,sqrt(pvar.*log(128*pmean./(81*pi*(1-p))))));
d2 = floor(max(1,sqrt(pvar.*log(128*qmean./(pi*p)))));

s1 = sqrt(pvar).*(1+d1/(4.*pmean));
s2 = sqrt(pvar).*(1+d2/(4.*qmean));
c = 2.*d1./pmean;
a1 = 0.5.*exp(c).*s1*sqrt(2*pi);
a2 = 0.5.*s2*sqrt(2*pi);
a3 = exp(d1./qmean)*2*s1.^2./d1*exp(-0.5*(d1./s1).^2);
a4 = 2*s2.^2./d2*exp(-0.5*(d2./s2).^2);
s = a1+a2+a3+a4;

%logbnp  = pdfbin(pmean,n,p,'logp',true);
useLog = true;
logbnp = binom(n,pmean,useLog);

X = zeros(1,N);
V = X;
left = 1:N;
accept = zeros(1,N);
while any(left)
  u = s.*rand(size(left));
  r1 = find(u<=a1);
  if any(r1)
    N1 = randn(size(r1));
    Y = s1*abs(N1);
    accept(r1)  = Y<d1;
    if any(accept(r1))
      X(left(r1)) = floor(Y);
      E = rndexp(1,size(r1));
      V(r1) = -E-N1.^2/2+c;
    end;
  end
  r2 = find(a1<u & u<=a1+a2);
  if any(r2)
    N1 = randn(size(r2));
    Y = s2*abs(N1);
    accept(r2)  = Y<d2;
    if any(accept(r2))
      X(left(r2)) = floor(-Y);
      E = rndexp(1,size(r2));
      V(r2) = -E-N1.^2/2;
    end;
  end

  r3 = find(a1+a2<u & u<=a1+a2+a3);
  if any(r3)
    E1 = rndexp(1,size(r3));
    E2 = rndexp(1,size(r3));
    Y = d1+2*s1.^2.*E1/d1;
    
    accept(r3)  = true;

    X(left(r3)) = floor(Y);
    V(r3) = -E2-d1.*Y./(2*s1^2)+d1/qmean;

  end


  r4 = find(a1+a2+a3<u);
  if any(r4)
    E1 = rndexp(1,size(r4));
    E2 = rndexp(1,size(r4));
    Y = d2+2*s2^2*E1/d2;
  
    accept(r4)  = true;

    X(left(r4)) = floor(-Y);
    V(r4) = -E2-d2*Y./(2*s2.^2);

  end

  accept = accept & ~((X(left)< -pmean) | (X(left)>qmean));
  if any(accept)
    ok = find(accept);
    %logbnpx = pdfbin(pmean+X(left(ok)),n,p,'logp',true);
    logbnpx = binom(n,pmean+X(left(ok)),useLog) + X(left(ok)).*log(p./(1-p));
    accept(ok) =   V(ok)<=logbnpx-logbnp;
    left(accept) = [];
    V(accept) = [];
    
    accept(accept) = [];
  end
end % while
X = X +pmean;




function X = rndbinrejectv(n,p)
% Rejection sampling for the binomial distribution for p<=0.5
% Assumption n*p is an integer

X = zeros(size(n));

[pmean, pvar] = mombin(n,p);
[qmean] = mombin(n,1-p);

d1 = floor(max(1,sqrt(pvar.*log(128*pmean./(81*pi*(1-p))))));
d2 = floor(max(1,sqrt(pvar.*log(128*qmean./(pi*p)))));

s1 = sqrt(pvar).*(1+d1./(4.*pmean));
s2 = sqrt(pvar).*(1+d2./(4.*qmean));
c = 2.*d1./pmean;
a1 = 0.5.*exp(c).*s1*sqrt(2*pi);
a2 = 0.5.*s2.*sqrt(2*pi);
a3 = exp(d1./qmean)*2.*s1.^2./d1.*exp(-0.5*(d1./s1).^2);
a4 = 2*s2.^2./d2.*exp(-0.5*(d2./s2).^2);
s = a1+a2+a3+a4;

%logbnp  = pdfbin(pmean,n,p,'logp',true);

useLog = true;
logbnp = binom(n,pmean,useLog);

V = X;
left = X;
left(:) = 1:numel(X);
accept = X;
while any(left(:))
  u = s.*rand(size(left));
  r1 = find(u<=a1);
  if any(r1(:))
    N1 = randn(size(r1));
    Y = s1(r1).*abs(N1);
    accept(r1)  = Y<d1(r1);
    if any(accept(r1))
      X(left(r1)) = floor(Y);
      E = rndexp(1,size(r1));
      V(r1) = -E-N1.^2/2+c(r1);
    end;
  end
  r2 = find(a1<u & u<=a1+a2);
  if any(r2(:))
    N1 = randn(size(r2));
    Y = s2(r2).*abs(N1);
    accept(r2)  = Y<d2(r2);
    if any(accept(r2))
      X(left(r2)) = floor(-Y);
      E = rndexp(1,size(r2));
      V(r2) = -E-N1.^2/2;
    end;
  end

  r3 = find(a1+a2<u & u<=a1+a2+a3);
  if any(r3(:))
    E1 = rndexp(1,size(r3));
    E2 = rndexp(1,size(r3));
    Y = d1(r3)+2*s1(r3).^2.*E1./d1(r3);

    accept(r3)  = true;
    X(left(r3)) = floor(Y);
    V(r3) = -E2-d1(r3).*Y./(2*s1(r3).^2)+d1(r3)./qmean(r3);
  end


  r4 = find(a1+a2+a3<u);
  if any(r4(:))
    E1 = rndexp(1,size(r4));
    E2 = rndexp(1,size(r4));
    Y = d2(r4)+2*s2(r4).^2.*E1./d2(r4);
  
    accept(r4)  = true;
    X(left(r4)) = floor(-Y);
    V(r4) = -E2-d2(r4).*Y./(2*s2(r4).^2);
  end

  accept = accept & ~((X(left)< -pmean) | (X(left)>qmean));
 
  if any(accept)
    ok = find(accept);
    %logbnpx = pdfbin(pmean(ok)+X(left(ok)),n(ok),p(ok),'logp',true);
    logbnpx = binom(n(ok),pmean(ok)+X(left(ok)),useLog) + X(left(ok)).*log(p(ok)./(1-p(ok)));
    accept(ok) =   V(ok)<=logbnpx-logbnp(ok);
    X(left(accept)) = X(left(accept)) + pmean(accept);
    left(accept) = [];
    V(accept) = [];
    
    d1(accept) = [];
    d2(accept) = [];

    s1(accept) = [];
    s2(accept) = [];
    c(accept) = [];
    a1(accept) = [];
    a2(accept) = [];
    a3(accept) = [];
    a4(accept) = [];
    s(accept) = [];

    logbnp(accept) = [];
    p(accept) = [];
    n(accept) = [];
    pmean(accept) = [];
    qmean(accept) = [];
    
    accept(accept) = [];
  end
end % while


% Old test calls
% return
% X = invbin(rand(csize),n,p);
% return
% if isscalar(n) && isscalar(p)
%   X = zeros(csize);
%   X(:) = rndbinScalar(n,p,numel(X));
% else
%   X = rndbinVector(n,p);
% end
% function X = rndbinVector(n,p)
% [csize, n,p] = comnsize(n,p);
% 
% [pmean,pvar] = mombin(n(:),p(:));
% mode = floor((n+1).*p);
% M = 2./sqrt(2*pi.*pvar(:));
% % Universal rejection sampling algoritm
% 
% %s2 = pvar+ min(1,pmean);
% %rho = exp((log(3*s2)+2.*log(M))/3);
% rho = exp(log(6*(1+min(1,pmean)./pvar)/pi)/3);
% rhoSplit = rho./(3.*rho+M);
% X = zeros(csize);
% T = X(:);
% 
% left = (1:numel(X)).';
% while any(left)
%   csize = size(left);
%   U = rand(csize);
%   V  = rand(csize);
%   W  = 2*rand(csize)-1;
%   k1 = find(U<rhoSplit);
%   if any(k1)
%     Y = mode(k1) +sign(V(k1)).*(0.5+rho(k1)./(M(k1).*sqrt(abs(V(k1)))));
%     T(k1) = W(k1).*M(k1).*abs(V(k1)).^(2/3);
%     X(left(k1)) = round(Y);
%   end
%   k2 = find(U>=rhoSplit);
%   if any(k2)
%     Y = mode(k2)+(0.5+rho(k2)./M(k2)).*V(k2);
%     T(k2) = W(k2).*M(k2);
%     X(left(k2)) = round(Y);
%   end
%   accept  = (T>pdfbin(X(left),n(left),p(left)));
%   if any(accept)
%     left(accept) = [];
%     T(accept) = [];
%     mode(accept) = [];
%     M(accept) = [];
%     rho(accept) = [];
%     rhoSplit(accept) = [];
%   end
% end % while
% 
% 
% 
% function X = rndbinScalar(n,p,N)
% [pmean,pvar] = mombin(n,p);
% mode = floor((n+1).*p);
% M = 2/sqrt(2*pi.*pvar);
% % Universal rejection sampling algoritm
% 
% %s2 = pvar+ min(1,pmean);
% %rho = exp((log(3*s2)+2.*log(M))/3);
% rho = exp(log(6*(1+min(1,pmean)./pvar)/pi)/3);
% rhoSplit = rho/(3.*rho+M);
% X = zeros(1,N);
% T = X;
% left = 1:N;
% while any(left)
%   csize = size(left);
%   U = rand(csize);
%   V  = rand(csize);
%   W  = 2*rand(csize)-1;
%   k1 = find(U<rhoSplit);
%   if any(k1)
%     Y = mode +sign(V(k1)).*(0.5+rho./(M*sqrt(abs(V(k1)))));
%     T(k1) = W(k1).*M.*abs(V(k1)).^(2/3);
%     X(left(k1)) = round(Y);
%   end
%   k2 = find(U>=rhoSplit);
%   if any(k2)
%     Y = mode+(0.5+rho/M).*V(k2);
%     T(k2) = W(k2).*M;
%     X(left(k2)) = round(Y);
%   end
%   accept = (T>pdfbin(X(left),n,p));
%   if any(accept)
%     left(accept) = [];
%     T(accept) = [];
%   end
% end % while