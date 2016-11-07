function [int, abserr,ix]= gaussq(fun,A,B,reltol,trace,varargin)
%GAUSSQ Numerically evaluate integral, Gauss quadrature.
%
% CALL:
%  [int, err] = gaussq('Fun',A,B,[reltol wfun],[trace,gn],p1,p2,....)
%  [int, err] = gaussq('Fun',A,B,[reltol 3],[trace,gn],alpha,p1,p2,....)
%  [int, err] = gaussq('Fun',A,B,[reltol 4],[trace,gn],alpha,beta,p1,p2,....)
%
%       int = evaluated integral
%       err = error estimate, absolute tolerance abs(int-intold)
%       Fun = inline object, function handle or a function string.  
%             The function may depend on the parameters alpha and beta.
%       A,B = lower and upper integration limits, respectively.
%    reltol = relative tolerance (default 1e-3).
%      wfun = integer defining the weight function, p(x). (default wfun = 1)     
%           1 : p(x) = 1                       a =-1,   b = 1   Gauss-Legendre 
%           2 : p(x) = exp(-x^2)               a =-inf, b = inf Hermite
%           3 : p(x) = x^alpha*exp(-x)         a = 0,   b = inf Laguerre
%           4 : p(x) = (x-a)^alpha*(b-x)^beta  a =-1,   b = 1 Jacobi 
%           5 : p(x) = 1/sqrt((x-a)*(b-x)),    a =-1,   b = 1 Chebyshev 1'st kind
%           6 : p(x) = sqrt((x-a)*(b-x)),      a =-1,   b = 1 Chebyshev 2'nd kind
%           7 : p(x) = sqrt((x-a)/(b-x)),      a = 0,   b = 1
%           8 : p(x) = 1/sqrt(b-x),            a = 0,   b = 1
%           9 : p(x) = sqrt(b-x),              a = 0,   b = 1
%   trace = for non-zero TRACE traces the function evaluations 
%              with a point plot of the integrand (default 0).
%      gn = number of base points and weight points to start the 
%            integration with (default 2).
% alpha, 
%    beta = Shape parameters of Laguerre or Jacobi weight function
%           (alpha,beta>-1) (default alpha=beta=0)
%p1,p2,...= coefficients to be passed directly to function Fun:   
%                  G = Fun(x,p1,p2,...).
%
% GAUSSQ numerically evaluate integral using a Gauss quadrature. 
% The Quadrature integrates a (2m-1)th order polynomial exactly and the 
% integral is of the form
%              b                              
%             Int (p(x)* Fun(x)) dx 
%              a 
% GAUSSQ is vectorized to accept integration limits A, B and 
% coefficients P1,P2,...Pn, as matrices or scalars and the 
% result INT is the common size of A, B and P1,P2,...,Pn.
%
% Examples :% a) integration of x.^2        from 0 to 2 and from 1 to 4
%           % b) integration of x^2*exp(-x) from zero to infinity 
%           % c) integrate humps            from 0 to 2 and from 1 to 4
%
%  A = [0 1]; B = [2,4];
%  [val1,err1] = gaussq('(x.^2)',A,B);                 % a)
%  [val2,err2] = gaussq('(1)',0,inf,[1e-3 3],[],2);    % b)
%  [val2,err2] = gaussq('(x.^2)',0,inf,[1e-3 3],[],0); % b)  
%  [val3,err3] = gaussq(@humps,A,B);                   % c)
%  assert(val1, [2.66666666666667, 21], 1e-10)
%  assert(val2, 2, 1e-10)
%  assert(val3, [34.92621394352225, 5.76237544601305], 1e-7)
%
% See also  qrule, gaussq2d


% References 
%
% [1]  Golub, G. H. and Welsch, J. H. (1969)
% 'Calculation of Gaussian Quadrature Rules'
%  Mathematics of Computation, vol 23,page 221-230,
%
% [2] Davis and Rabinowitz (1975) 'Methods of Numerical Integration', page 365,
%     Academic Press.
%
% [3]. Stroud and Secrest (1966), 'gaussian quadrature formulas', 
%      prentice-hall, Englewood cliffs, n.j.
% 
% [4] Abromowitz and Stegun (1954) 'Handbook of mathematical functions'




% tested on: Matlab  7.0
% history:
% Revised pab June 2007
% -fixed a bug for wfun==3: alpha1 was not transferred to lrule, now fixed
% Revised pab 5April2006
% -removed all eval statements.
% -new ordering of quadrature-rules.
% Revised pab 25nov2005
% -fixed a bug when size(A) = [1 1] ~= size(P1)= size(P2), ...
% Revised pab 2Nov2005
% -changed global variables to persistent variables.
% Revised pab 22Nov2004
% -Added the possibility of using a function handle.  
% Revised pab 09.09.2002
% -added the possibility of using a inline function.
% revised pab 27.03.2000
%  - fixed a bug for p1,p2,... changed to varargin in input
% revised pab 19.09.1999
%   documentation
%  by Per A. Brodtkorb 30.03.99, 17.02.99 :
%  -accept multiple integrationlimits, int is the common size of xlow and xhigh
%  -optimized by only computing the integrals which did not converge.
%  -enabled the integration of directly given functions enclosed in 
%     parenthesis. Example: integration from 0 to 2 and from 2 to 4 for x is done by:
%                        gaussq('(x.^2)',[0 2],[2 4])


persistent ALPHA1 BETA1 ALPHA2 cb cw
wfun    = 1;
maxIter = 11;

if nargin<4 || isempty(reltol),
  reltol=1e-3;
elseif length(reltol)==2,
  wfun=reltol(2);
   reltol=reltol(1);
end

P0 = varargin;
NP = length(P0);

istart      = 0; 
alpha1      = 0;
beta1       = 0;
FindWeights = 1;
wfun1       = mod(wfun-1,10)+1;

switch wfun1
  case 4,
    istart=2;
    if ((NP>=1) && (~isempty(P0{1}))),
      alpha1 = P0{1}; 
    end
    if ((NP>=2) && (~isempty(P0{2}))), 
      beta1  = P0{2}; 
    end
    if isempty(ALPHA1) || isempty(BETA1),
    elseif (ALPHA1==alpha1) && (BETA1==beta1),
      FindWeights = 0;
    end
    ALPHA1 = alpha1;
    BETA1  = beta1;
    %remember what type of weights are saved as global constants
  case 3,
    istart=1;
    if ((NP>=1) && (~isempty(P0{1}))),
      alpha1 = P0{1};  
    end
    if isempty(ALPHA2),
    elseif ALPHA2==alpha1,
      FindWeights=0;
    end
    ALPHA2=alpha1;
    %remember what type of weights are saved as global constants
  otherwise,
    FindWeights=0;
end

P0(1:istart)=[];

gn=2;
if nargin <5 || isempty(trace) ,
  trace = 0; 
elseif length(trace)==2,
  gn    = trace(2);
  trace = trace(1);
end

aSize = size(A); % remember the size of input
bSize = size(B);

if prod(aSize)==1,% make sure the integration limits have correct size
  A = A(ones(bSize));
  aSize = bSize;
elseif prod(bSize)==1,
  B = B(ones(aSize));
elseif any( aSize~=bSize)
  error('The integration limits must have equal size!')
end




isFunctionHandle = isa(fun,'function_handle');
if (not(isFunctionHandle) && isa(fun,'char') &&  any(fun=='(') ), %  & any(fun=='x'),
  %exec_string=['y=',fun ';']; %the call function is already setup
  fun = inline(fun);
end

num_parameters = NP-istart;
isvector1       = zeros(1,num_parameters);

nk    = prod(aSize); % # of integrals we have to compute
for ix=1:num_parameters,
  p0Size = size(P0{ix});
  Np0    = prod(p0Size);
  isvector1(ix) = isnumeric(P0{ix}) & (Np0 > 1);
  if isvector1(ix)
    if  nk==1,
      aSize = p0Size;
      nk    = Np0;
      A = A(ones(aSize));
      B = B(ones(aSize));
    elseif  nk~=Np0
      error('The input must have equal size!')
    end
    P0{ix} = P0{ix}(:); %make sure it is a column
  end
end


k       = (1:nk)';
int     = zeros(nk,1);
int_old = int;
abserr    = int;


if isempty(cb)
  % Initialize the size of cell array containing the basepoints and weights 
  cb = cell(maxIter,10);
  cw = cell(maxIter,10);
end
if gn~=length(cb{1,wfun})
  FindWeights = 1;
end

if FindWeights
  [cb{:,wfun1}] = deal([]);
  [cw{:,wfun1}] = deal([]);
end

%setup mapping parameters
A     = A(:);
jacob = (B(:)-A(:))/2;

shift = 1;
switch wfun1,
  case {1}, % Gauss-legendre
    dx = jacob;
  case {2,3}
    shift = 0;
    jacob = ones(nk,1);
    A     = zeros(nk,1);
    dx    = jacob;
  case {4}
     dx = jacob.^(alpha1+beta1+1);
  case 5,  
    dx = ones(nk,1);
  case 6, 
    dx = jacob.^2;
  case 7,  
    shift = 0;
    jacob = jacob*2;
    dx    = jacob;
  case 8,  
    shift = 0;
    jacob = jacob*2;
    dx    = sqrt(jacob);
  case 9,  
    shift = 0;
    jacob = jacob*2;
    dx    = sqrt(jacob).^3;
  otherwise
    error('unknown option')
end


if trace>0,
  x_trace = cell(maxIter,1);
  y_trace = cell(maxIter,1);
end

if num_parameters>0
  ixVector = find(isvector1);
  if any(ixVector)
    ixScalar = find(~isvector1);
    P1 = cell(1,num_parameters);
    if any(ixScalar)
      P1(ixScalar) = P0(ixScalar);
      %[P1{ixScalar}] = deal(P0{ixScalar});
    end
  end
end
% Break out of the iteration loop for three reasons:
%  1) the last update is very small (compared to int  and  compared to reltol)
%  2) There are more than 11 iterations. This should NEVER happen. 

converge='n';
for ix=1:maxIter,
  if isempty(cb{ix,wfun1}) || FindWeights ,  
    % calculate the weights and base points if necessary
     [cb{ix,wfun1},cw{ix,wfun1}]=qrule(gn,wfun,alpha1,beta1);
  end
  % calculate the x values
  x  = (cb{ix,wfun1}(ones(nk,1),:)+shift).*jacob(k,ones(1, gn )) + A(k,ones(1,gn ));
  
  
  % calculate function values  y=fun(x,p1,p2,....,pn)
  if num_parameters>0
    if any(ixVector)
      % Expand vector to the correct size
      for iy = ixVector,
        P1{iy} = P0{iy}(k,ones(1,gn));
      end
      y  = feval(fun,x,P1{:});
    else
      y  = feval(fun,x,P0{:});
    end
  else
    y = feval(fun,x);
  end
  
  int(k) = sum(cw{ix,wfun1}(ones(nk,1),:).*y,2).*dx(k); % do the integration sum(y.*w)
  
  
  if trace>0,
    x_trace{ix} = x(:).';
    y_trace{ix} = y(:).';
    
    hfig = plot(x,y,'r.'); hold on
    drawnow,shg
    if (trace>1)
      pause
    end
    set(hfig,'color','b')
  end

  abserr(k) = abs(int_old(k)-int(k)); %absolute tolerance
  if ix > 1
    k       = find(abserr > abs(reltol*int)); %| abserr > abs(reltol));%indices to integrals which did not converge
  end
  if any(k),% compute integrals again
    nk         = length(k);%# of integrals we have to compute again
    int_old(k) = int(k);
  else
    converge = 'y';
    break;
  end
  gn    = gn*2;% double the # of basepoints and weights
end

int = reshape(int,aSize); % make sure int is the same size as the integration  limits
if nargout>1,
  abserr=reshape(abserr,aSize);
end

if converge=='n'
  if nk>1
    if (nk==prod(aSize)),
      tmptxt = 'All integrals did not converge--singularities likely!';
    else
      tmptxt = sprintf('%d integrals did not converge--singularities likely!',nk);
    end
  else
    tmptxt = 'Integral did not converge--singularity likely!';
  end
  warning('GAUSSQ:SINGULARITY',tmptxt)
end

 if trace>0,
   clf
   plot([x_trace{:}],[y_trace{:}],'+')
 end
 
%!test
%! A = [0 1]; 
%! B = [2,4];
%! [val1,err1] = gaussq('(x.^2)',A,B);
%! assert(val1, [ 2.66666666666667, 21.00000000000000], 1e-13)

%!test
%! A = [0 1]; 
%! B = [2,4];
%! [val2,err2] = gaussq('(1)',0,inf,[1e-3 3],[],2);
%! assert(val2, 1, 1e-13)

%!test
%! A = [0 1]; 
%! B = [2,4];
%! [val3,err3] = gaussq('(x.^2)',0,inf,[1e-3 3],[],0);  
%! assert(val3, 2, 1e-13)
