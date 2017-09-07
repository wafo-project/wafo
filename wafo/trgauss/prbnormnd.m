function [value,err,inform,exTime] = prbnormnd(correl,A,B,abseps,releps,maxpts,method)
%PRBNORMND Multivariate Normal probability by Genz' algorithm.
%
%  CALL [value,error,inform]=prbnormnd(correl,A,B,abseps,releps,maxpts,method);
%
%     VALUE  REAL estimated value for the integral
%     ERROR  REAL estimated absolute error, with 99% confidence level.
%     INFORM INTEGER, termination status parameter:
%            if INFORM = 0, normal completion with ERROR < EPS;
%            if INFORM = 1, completion with ERROR > EPS and MAXPTS 
%                           function vaules used; increase MAXPTS to 
%                           decrease ERROR;
%            if INFORM = 2, N > NMAX or N < 1. where NMAX depends on the
%                           integration method
%
%     CORREL = Positive semidefinite correlation matrix
%     A	     = vector of lower integration limits.
%     B	     = vector of upper integration limits.
%     ABSEPS = absolute error tolerance.
%     RELEPS = relative error tolerance.
%     MAXPTS = maximum number of function values allowed. This 
%              parameter can be used to limit the time. A sensible 
%              strategy is to start with MAXPTS = 1000*N, and then
%              increase MAXPTS if ERROR is too large.
%     METHOD = integer defining the integration method
%             -1 KRBVRC randomized Korobov rules for the first 20
%                variables, randomized Richtmeyer rules for the rest, 
%                NMAX = 500 
%              0 KRBVRC, NMAX = 100 (default)
%              1 SADAPT Subregion Adaptive integration method, NMAX = 20 
%              2 KROBOV Randomized KOROBOV rules,              NMAX = 100
%              3 RCRUDE Crude Monte-Carlo Algorithm with simple
%                antithetic variates and weighted results on restart 
%              4 SPHMVN Monte-Carlo algorithm by Deak (1980),  NMAX = 100
%
% Example:% Compute the probability that X1<0,X2<0,X3<0,X4<0,X5<0,
%           % Xi are zero-mean Gaussian variables with variances one
%           % and correlations Cov(X(i),X(j))=0.3:
%           % indI=[0 5], and barriers B_lo=[-inf 0], B_lo=[0  inf]     
%           % gives H_lo = [-inf -inf -inf -inf -inf]  H_lo = [0 0 0 0 0] 
% 
%    N = 5; rho=0.3; NIT=3; Nt=N; indI=[0 N];
%    B_lo=-10; B_up=0; m=1.2*ones(N,1);
%    Sc=(ones(N)-eye(N))*rho+eye(N);
%    E = rind(Sc,m,B_lo,B_up,indI,[],Nt) % exact prob. 0.00195
%    A = [-inf -inf -inf -inf -inf],
%    B = [0 0 0 0 0]-m' 
%    [val,err,inform] = prbnormnd(Sc,A,B);  
%
% See also  prbnormndpc, rind

%History  
% revised pab April 2008
% -renamed from mvnormprb to prbnormnd
% revised pab 2007
% fixed a bug correlation is now oriented correct
% Revised pab April 2006
% -Fixed a bug: zero correlations are now allowed,  thanks to Daniel
% Lewandowski who made us aware of this.
% By pab 2002

%error(nargchk(3,7,nargin))
narginchk(3,7)
[m,n] = size(correl);
Na = length(A);
Nb = length(B);
if (m~=n || m~=Na || m~=Nb)
   error('Size of input is inconsistent!')
end

if nargin<4 || isempty(abseps), abseps = 1e-4; end
if nargin<5 || isempty(releps), releps = 1e-3; end
if nargin<6 || isempty(maxpts), maxpts = 1000*n; end
if nargin<7 || isempty(method), method = 0; end

maxpts = max(round(maxpts),10*n);

%            array of correlation coefficients; the correlation
%            coefficient in row I column J of the correlation matrix
%            should be stored in CORREL( J + ((I-2)*(I-1))/2 ), for J < I.
%            The correlation matrix must be positive semidefinite.

D = diag(correl);
if (any(D~=1))
   error('This is not a correlation matrix')
end

% Make sure integration limits are finite
A = min(max(A,-100),100);
B = max(min(B,100),-100);
%L = correl((tril(ones(m),-1)~=0));    % return only off diagonal elements
% New call:
L = correl((triu(ones(m),1)~=0));    % return only off diagonal elements

%CALL the mexroutine
t0 = clock;
if ((method==0) && (n<=100)),
  %NMAX = 100
  [value, err,inform] = mexmvnprb(L,A,B,abseps,releps,maxpts);
elseif ( (method<0) || ((method<=0) && (n>100)) ),
  % NMAX = 500
  [value, err,inform] = mexmvnprb2(L,A,B,abseps,releps,maxpts);
else
  [value, err,inform] = mexGenzMvnPrb(L,A,B,abseps,releps,maxpts,method);
end
exTime = etime(clock,t0);


return
% if m>100 | m<1
%    value = 0;
%    err = 1;
%    inform = 2;
%    return
% end