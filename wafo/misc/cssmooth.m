function [yy,p]= cssmooth(x,y,p,xx,LinExtrap,d2)
%SMOOTH Calculates a smoothing spline.
% 
%  CALL:  yy = cssmooth(x,y,p,xx,def,v)
%
%         x  = x-coordinates of data. (vector)
%         y  = y-coordinates of data. (vector or matrix)
%         p  = [0...1] is a smoothing parameter:
%              0 -> LS-straight line
%              1 -> cubic spline interpolant
%         xx = the x-coordinates in which to calculate the smoothed function.
%         yy = the calculated y-coordinates of the smoothed function if
%              xx is given. If xx is not given then yy contain the
%              pp-form of the spline, for use with PPVAL.
%        def = 0 regular smoothing spline (default)
%              1 a smoothing spline with a constraint on the ends to 
%                ensure linear extrapolation outside the range of the data
%          v = variance of each y(i) (default  ones(length(X),1)) 
%
%  Given the approximate values 
%  
%                y(i) = g(x(i))+e(i) 
%  
%  of some smooth function, g, where e(i) is the error. SMOOTH tries to 
%  recover g from y by constructing a function, f, which  minimizes
%
%              p * sum (Y(i) - f(X(i)))^2/d2(i)  +  (1-p) * int (f'')^2
%
%  The call  pp = cssmooth(x,y,p)  gives the pp-form of the spline, 
%  for use with PPVAL. 
%
% Example:%
%  x = linspace(0,1).';
%  y = exp(x)+1e-1*randn(size(x));
%  pp = cssmooth(x,y,.9); 
%  plot(x,y,x,cssmooth(x,y,.99,x,0,0.01),'g',x,ppval(pp,x),'k',x,exp(x),'r')
%
% See also  lc2tr, dat2tr, ppval


%References:
% Carl de Boor (1978)
% 'Practical Guide to Splines'
%  Springer Verlag
%  Uses EqXIV.6--9, pp 239

% tested on: Matlab 4.x 5.x 8.x
%History:
% revised by GL Feb 2011
% -changed name from smooth to cssmooth to avoid name collisions in MATLAB 
% revised pab Feb 2007
% -added extrapolate -> better extrapolation algorithm if few points
% -moved some code into computeU.
% revised pab  26.11.2000
% - added example
% - fixed a bug: n=2 is now possible
% revised by pab 21.09.99
%    - added d2
% revised by pab 29.08.99
%    -new extrapolation: ensuring that
%     the smoothed function has contionous derivatives
%     in the first and last knot
% modified by Per A. Brodtkorb  23.09.98 
%  secret option forcing linear extrapolation outside the ends when p>0
%  used in lc2tr

if nargin<3 || isempty(p)
  p = [];
else
  p = min(p,1);
end
if (nargin<5)||(isempty(LinExtrap)),
  LinExtrap=0; %do not force linear extrapolation in the ends (default)
end

dx       = diff(x(:));
mustSort = any(dx<0);
if mustSort
  [xi,ind] = sort(x(:));
  dx = diff(xi);
else
  xi = x(:);
end
n = length(xi);

if nargin<6||isempty(d2), 
  d2 = ones(n,1); 
 % d2 = ([dx;10]+[10;dx])./2;
elseif length(d2) == 1
  d2 = d2(ones(n,1),:); 
  %d2 = d2.*([dx;10]+[10;dx])/2;
else
  d2=d2(:);
end



ndy = ndims (y);
szy = size (y);
if (ndy == 2 && (szy(1) == 1 || szy(2) == 1))
  if (szy(1) == 1)
    yi = y.';
  else
    yi = y;
    szy = fliplr (szy);
  end
else
  yi = reshape (y, [prod(szy(1:end-1)), szy(end)]).';
end
nd = prod(szy(1:end-1));
ny = szy(end);


if n<2,
   error('There must be >=2 data points.')
elseif any(dx<=0),
   error('Two consecutive values in x can not be equal.')
elseif n~=ny,
   error('x and y must have the same length.')
end


if mustSort
  yi = yi(ind,:);
end

idx = ones(nd,1);

dydx = diff(yi)./dx(:,idx);

if (n==2)  % straight line
  coefs=[dydx(:) yi(1,:).'];
else
  if LinExtrap==2 && n==3,
    p = 0;  % Force LS-fit
  end
  dx1=1./dx;
  
  u = computeU();
  
  zrs = zeros(1,nd);
  if p<1
    ai = yi-6*(1-p)*D*diff([zrs;diff([zrs;u;zrs]).*dx1(:,idx);zrs]); % faster than yi-6*(1-p)*Q*u
  else
    ai = yi;
  end
  % The piecewise polynominals are written as
  % fi=ai+bi*(x-xi)+ci*(x-xi)^2+di*(x-xi)^3 
  % where the derivatives in the knots according to Carl de Boor are:
  %    ddfi  = 6*p*[0;u] = 2*ci;
  %    dddfi = 2*diff([ci;0])./dx = 6*di;
  %    dfi   = diff(ai)./dx-(ci+di.*dx).*dx = bi;
 
  ci = [zrs;3*p*u];
   

  if LinExtrap==2 && p~=0 && n>3, %Forcing linear extrapolation in the ends 
    ci([2,  end],:) = 0;
    % New call
    % fixing the coefficients so that we have continous
    % derivatives everywhere
    ai(1,:) = -(ai(3,:)-ai(2,:))*dx(1)/dx(2) +ai(2,:)+ ci(3,:)*dx(1)*dx(2)/3;
    ai(n,:) = (ai(n-1,:)-ai(n-2,:))*dx(n-1)/dx(n-2) +ai(n-1,:)+ ci(n-2,:)*dx(n-2)*dx(n-1)/3;  
  end
  
  
  di    = (diff([ci;zrs]).*dx1(:,idx)/3);
  bi    = (diff(ai).*dx1(:,idx)-(ci+di.*dx(:,idx)).*dx(:,idx));
  ai(n:end,:) = [];
  if nd>1
    di = di.';
    ci = ci.';
    ai = ai.';
  end
  if ~any(di)
    if ~any(ci)
      coefs = [bi(:) ai(:)]; 
    else
      coefs = [ci(:) bi(:) ai(:)]; 
    end
  else
    coefs = [di(:) ci(:) bi(:) ai(:)]; 
  end
end

pp = mkpp(xi,coefs,szy(1:end-1));

if ~any(LinExtrap==[0 2])
  % Extrapolation strategy
  pp = extrapolate(pp);
end


if (nargin<4)||(isempty(xx)),
  yy = pp;
else
  yy = ppval(pp,xx);
end

%% Nested functions
  function u = computeU()
    if isempty(p) || p~=0
      R = spdiags([dx(2:n-1) 2*(dx(1:n-2)+dx(2:n-1)) dx(1:n-2)],-1:1,n-2,n-2);
    end
    if isempty(p) || p<1
      Q = spdiags([dx1(1:n-2) -(dx1(1:n-2)+dx1(2:n-1)) dx1(2:n-1)],0:-1:-2,n,n-2);
      D = spdiags(d2,0,n,n);  % The variance
      QDQ = Q.'*D*Q;
      if isempty(p) || p<0
        % Crude estimate
        p = 1/(1+trace(QDQ)/(100*trace(R).^2));
      end
      if p==0
        QQ = (6*(1-p))*(QDQ);
      else
        QQ = (6*(1-p))*(QDQ)+p*R;
      end
      %clear QDQ
    else
      QQ = R;
    end
    % Make sure Matlab uses symmetric matrix solver
    u  = 2*((QQ+QQ.')\diff(dydx));  % faster than u=QQ\(Q'*yi);
  end
end % function cssmooth

%% Subfunctions
function g = extrapolate(pp)
%EXTRAPOLATE Extrapolate a 1D PP linearly outside its basic interval


maxOrder = 2;

[breaks,coefs,pieces,order,dim]=unmkpp(pp);
if order<=maxOrder
  g = pp;
  return
end



%Add new breaks beyond each end 
breaks2add = breaks([1,end]) + [-1,1];
newbreaks = [breaks2add(1),breaks, breaks2add(2)];


dx = newbreaks([1 end-1]) - breaks([1 end-1]);

ifl  = [1 pieces]; % index to first and last polynomial piece.
if dim>1 % repeat each point dim  times if necessary
  dx = repmat(dx,dim,1);
  ifl = repmat(dim*ifl,dim,1)+repmat((1-dim:0).',1,2);
end
dx = dx(:);
nx = length(dx);

% Get coefficients for the new last polynomial piece (aN)
% by just relocate the previous last polynomial and 
% then set all terms of order > maxOrder to zero
ix  = dim+1:nx;
aN  = coefs(ifl(ix),:); 
dxN = dx(ix);
% Relocate last polynomial using Horner's algorithm
aN = polyreloc(aN,dxN);

%set to zero all terms of order > maxOrder in first polynomial
aN(:,1:order-maxOrder) = 0;

%Get the coefficients for the new first piece (a1)
% by first setting all terms of order > maxOrder to zero and then
% relocate the polynomial.

ix = 1:dim;

%Set to zero all terms of order > maxOrder, i.e., not using them
a1 = coefs(ifl(ix),maxOrder+1:end);
% alternatively
% a1 = coefs(lu(ix),:);
% a1(:,1:order-maxOrder) = [];

dx1 = dx(ix);

% Relocate first polynomial using Horner's algorithm
a1 = polyreloc(a1,dx1);
a1 = [zeros(dim,order-maxOrder),a1];

% 
% % fixing the coefficients so that we have continous
% % derivatives everywhere
% a1=-(ai(2)-ai(1))*dx(1)/dx(2) +ai(1)+ ci(3)*dx(1)*dx(2)/3;
% an=(ai(n-2)-ai(n-3))*dx(n-1)/dx(n-2) +ai(n-2)+ ci(n-2)*dx(n-2)*dx(n-1)/3;
% ai=[a1;ai; an];


newcoefs = [ a1; coefs; aN];


g = mkpp(newbreaks,newcoefs,dim);


end % function extrapolate



function r = polyreloc(p,dx)
%POLYRELOC Relocate polynomial.
%
% CALL  R = POLYRELOC( P, dX) 
% 
% R  = vector/matrix of relocated polynomial coefficients.
% P  = vector/matrix of polynomial coefficients to relocate
% dx = distance to relocate P.
%
% POLYRELOC relocates the polynomial P by "moving" it dX
% units to the left along the x-axis. So R is
%   relative to the point (-dX,0) as P is relative to the point (0,0).
%
%   P is a matrix of row vectors of coefficients in decreasing order.

n = size( p ,2);
    

% Relocate polynomial using Horner's algorithm
r = p;
for ii=n:-1:2
  for i=2:ii
    r(:,i) = dx.*r(:,i-1)+r(:,i);
  end
end


end % function
