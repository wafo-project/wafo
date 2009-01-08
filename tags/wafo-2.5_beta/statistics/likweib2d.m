function [LL, C] = likweib2d(param1,data1,data2,given,gparam)
% LIKWEIB2D 2D Weibull log-likelihood function.
%
% CALL:  [L, cov] = likweib2d(phat,data1,data2) 
%
%   L           = log-likelihood of the parameters given the data
%   cov         = Asymptotic covariance matrix of phat (if phat is estimated by
%                 a maximum likelihood method).
%   phat        = [A1 B1 A2 B2 C12] vector of distribution parameters
%   data1,data2 = data vectors
%
%   LIKWEIB2D is a utility function for maximum likelihood estimation. 
%   The PDF is defined by:
%
% f(X1,X2) = B1*B2*xn1^(B1-1)*xn2^(B2-1)/A1/B1/N*...
%            exp{-[xn1^B1 +xn2^B2 ]/N }*I0(2*C12*xn1^(B1/2)/N) 
%  where 
%    N=1-C12^2, xn1=X1/A1,  xn2=X2/A2 and 
%    I0 is the modified bessel function of zeroth order.
%
%   See also  pdfweib2d

%tested on: matlab 5.1
% history:
% revised pab 1.11.2000
% - improoved the calculation of cov.
%  by Per A. Brodtkorb 14.11.98 

% Secret options:
%   given       =  a vector with  Number to the given parameter: [ 1 3] means
%                  parameter 1 and 3 are fixed
%   gparam      = values of the given parameters which we consider fixed


error(nargchk(3,5,nargin))

data1 = data1(:);
data2 = data2(:);
n  = length(data1);
n2 = length(data2);
if n~=n2
  error('data1 and data2  must have equal size')
end

sparams=zeros(1,5);
sparams(given)=gparam;
iz=1:5;iz(given)=[];
sparams(iz)=param1;

x = pdfweib2d(data1,data2,sparams)+eps;
LL = -sum(log(x)); % log likelihood function

if nargout > 1
  params = sparams;
  delta = eps^.4;
  delta2=delta^2;
  np=length(param1);
  dist ='pdfweib2d';
  % Approximate 1/(nE( (d L(x|theta)/dtheta)^2)) with
  %             1/(d^2 L(theta|x)/dtheta^2) 
  %  using central differences
    
  H = zeros(np);             % Hessian matrix
  for ix=1:np,
    sparam = params;
    iw = iz(ix);
    sparam(iw)= params(iw)+delta;
    x  = feval(dist,data1,data2,sparam)+eps; 
    fp = sum(log(x));
    sparam(iw) = params(iw)-delta;
    x  = feval(dist,data1,data2,sparam)+eps; 
    fm = sum(log(x));
    H(ix,ix) = (fp+2*LL+fm)/delta2;
    for iy=ix+1:np,
      iu = iz(iy);
      sparam(iw) = params(iw)+delta;
      sparam(iu) = params(iu)+delta;
      x   = feval(dist,data1,data2,sparam)+eps; 
      fpp = sum(log(x));
      sparam(iu) = params(iu)-delta;
      x   = feval(dist,data1,data2,sparam)+eps; 
      fpm = sum(log(x));
      sparam(iw) = params(iw)-delta;
      x   = feval(dist,data1,data2,sparam)+eps; 
      fmm = sum(log(x));
      sparam(iu) = params(iu)+delta;
      x   = feval(dist,data1,data2,sparam)+eps; 
      fmp = sum(log(x));
      H(ix,iy) = (fpp-fmp-fpm+fmm)/(4*delta2);
      H(iy,ix) = H(ix,iy);
    end
  end
  % invert the Hessian matrix (i.e. invert the observed information number)
  C = -H\eye(np); 
end
