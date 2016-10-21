function [value,err,inform,exTime] = prbnormtnd(correl,A,B,delta,nu,abseps,releps,maxpts)
%PRBNORMTND Multivariate normal or T probability by Genz' algorithm.
%
%  CALL [value,error,inform]=prbnormtnd(cov,A,B,delta,nu,abseps,releps,maxpts);
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
%     COV = Positive semidefinite correlation matrix (size N x N)
%     A	     = vector of lower integration limits.  (length N)
%     B	     = vector of upper integration limits.
%     DELTA  = vector of non-centrality parameters. (default zeros(N,1)
%     NU     = the number of degrees of freedom.
%              If NU < 1, then an MVN probability is computed. (default -1)
%     ABSEPS = absolute error tolerance. (default 1e-4)
%     RELEPS = relative error tolerance. (default 1e-3)
%     MAXPTS = maximum number of function values allowed. This 
%              parameter can be used to limit the time. A sensible 
%              strategy is to start with MAXPTS = 1000*N, and then
%              increase MAXPTS if ERROR is too large.
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
%    E = rind(Sc,m,B_lo,B_up,indI) % exact prob. 0.00195
%    A = [-inf -inf -inf -inf -inf],
%    B = [0 0 0 0 0]-m' 
%    [val,err,inform] = prbnormtnd(Sc,A,B);  
%
% See also  prbnormndpc, rind

%History  
% revised pab April 2008
% -renamed from mvtnormprb to prbnormtnd
% By pab 2004
error(nargchk(3,8,nargin))
[m,n] = size(correl);
Na = length(A);
Nb = length(B);
if (m~=n || m~=Na || m~=Nb)
   error('Size of input is inconsistent!')
end
if nargin<4 || isempty(delta)
  delta = zeros(n,1);
end
if nargin<5 || isempty(nu)
  nu = -1;
end
if nargin<6 || isempty(abseps), abseps = 1e-4; end
if nargin<7 || isempty(releps), releps = 1e-3; end
if nargin<8 || isempty(maxpts), maxpts = 1000*n; end

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
constraint = eye(Na,n);

%CALL the mexroutine
t0 = clock;
%NMAX = 100
[value, err,inform] = mexmvtprb(correl,A,constraint,B,delta,nu,abseps, ...
				releps,maxpts);

exTime = etime(clock,t0);


return
