function [int, tol1,k] = gaussq2d(fun,Ax,Bx,Ay,By,tol,varargin)
%GAUSSQ2D Numerically evaluate2D integral, Gauss quadrature.
%
% CALL: [int, tol] = gaussq2d(fun,Ax,Bx,Ay,By,reltol,p1,p2,...)
%
%      int = evaluated integral
%      tol = estimated error, absolute tolerance abs(int-intold)
%      Fun = Fun(x,y,p1,p2,...), inline object, function handle, string containing
%            the name of the function or a directly given expression 
%            enclosed in parenthesis.
%            The function may depend on parameters pi.
%    Ax,Bx = lower and upper integration limits for x 
%    Ay,By = lower and upper integration limits for y 
%   reltol = relative tolerance (default 1e-3)
% p1,p2,...= coefficients to be passed directly to function Fun:   
%                  G = Fun(x,y,p1,p2,...).
%
% GAUSSQ numerically evaluate 2D integral using a Gauss quadrature. 
% GAUSSQ is vectorized to accept integration limits A, B and coefficients P1,P2,...
% as matrices or scalars and the result INT is the common size of A, B 
% and P1,P2,....
%
% Note that if there are discontinuities the region of integration 
% should be broken up into separate pieces.  And if there are singularities,
% a more appropriate integration quadrature should be used 
% (such as the Gauss-Chebyshev for a specific type of singularity).
%
%Example:%  integration for (x,y) in [0,1]x[-1,3] and [2,4]x[1,2]:
%
%   p1=2.; p2=0.5;
%   fun = '(p2*x.^2.*y+p1)';
%   assert(gaussq2d(fun,[0 2],[1 4],[-1 1],[3 2],[],p1,p2),[10, 12], 1e-10);
%
% See also  gaussq, qrule2d

% tested on: matlab 7
% history:
% Revised pab June 2007
% -removed all eval statements.
% -replace the use of globals with persistent variables
% revised pab Nov2004
% -wrong input to comnsize fixed  
% -enabled function handle as input  
% revised pab 12Nov2003
% replaced call to distchk with comnsize
%modified by Per A. Brodtkorb 17.11.98 : 
% -accept multiple integrations limits
% -optimized by saving the weights as global constants and only  
%   computing the integrals which did not converge 
% -enabled the integration of directly given functions enclosed in
%   parenthesis 
% -adopted from NIT toolbox changed name from quad2dg to gaussq2d  

persistent bpx bpy wfxy;

maxIter = 9;
if isempty(bpx)
  % Initialize the size of cell array containing 
  % the basepoints and weights 
  bpx = cell(maxIter,1);
  bpy = cell(maxIter,1);
  wfxy = cell(maxIter,1);
  [bpx{1},bpy{1},wfxy{1}] = qrule2d(2,2);
end



if nargin<6 || isempty(tol)
  tol=1e-3;
end
[csize ,Ax,Bx,Ay,By]=comnsize(Ax,Bx,Ay,By);
if any(isnan(csize))
  error('Requires non-scalar arguments to match in size.');
end

isFunctionHandle = isa(fun,'function_handle');
if (not(isFunctionHandle) && isa(fun,'char') &&  any(fun=='(') ), %  & any(fun=='x'),
  %exec_string=['y=',fun ';']; %the call function is already setup
  fun = inline(fun);
end

aSize = size(Ax);

P0 = varargin;
NP = length(P0);

isvector1       = zeros(1,NP);

nk    = prod(aSize); % # of integrals we have to compute
for ix=1:NP,
  p0Size = size(P0{ix});
  Np0    = prod(p0Size);
  isvector1(ix) = isnumeric(P0{ix}) & (Np0 > 1);
  if isvector1(ix)
    if  nk==1,
      aSize = p0Size;
      nk    = Np0;
      Ax = Ax(ones(aSize));
      Bx = Bx(ones(aSize));
      Ay = Ay(ones(aSize));
      By = By(ones(aSize));
    elseif  nk~=Np0
      error('The input must have equal size!')
    end
    P0{ix} = P0{ix}(:); %make sure it is a column
  end
end



%setup mapping parameters &  make sure they all are column vectors
Ax=Ax(:); 
Bx=Bx(:);
Ay=Ay(:);
By=By(:); 



%Map to x
qx=(Bx-Ax)/2;
px=(Bx+Ax)/2;

%Map to y
qy=(By-Ay)/2;
py=(By+Ay)/2;


k       = (1:nk).';
int     =zeros(nk,1);
int_old = int;
tol1    =int;


if NP>0
  ixVector = find(isvector1);
  if any(ixVector)
    ixScalar = find(~isvector1);
    P1 = cell(1,NP);
    if any(ixScalar)
      P1(ixScalar) = P0(ixScalar);
      %[P1{ixScalar}] = deal(P0{ixScalar});
    end
  end
end


% Break out of the iteration loop for three reasons:
%  1) the last update is very small (compared to int)
%  2) the last update is very small (compared to tol)
%  3) There are more than 8 iterations. This should NEVER happen. 

converge='n';
for ix=1:maxIter,
  gn = 2^ix;
  if isempty(bpx{ix})
    % calculate the weights and base points if necessary
     [bpx{ix},bpy{ix},wfxy{ix}] = qrule2d(gn,gn);
  end
  x=permute(qx(k,ones(1,gn), ones(1,gn) ),[ 2 3 1]).*bpx{ix}(:,:,ones(1,nk)) ...
    +permute(px(k,ones(1,gn), ones(1,gn) ),[2 3 1]);
  y=permute(qy(k,ones(1,gn), ones(1,gn) ),[ 2 3 1]).*bpy{ix}(:,:,ones(1,nk)) ...
    +permute(py(k,ones(1,gn), ones(1,gn) ),[2 3 1]);

  % calculate function values  y=fun(x,p1,p2,....,pn)
  if NP>0
    if any(ixVector)
      % Expand vector to the correct size
      for iy = ixVector,
        P1{iy} = permute(P0{iy}(k,ones(1,gn),ones(1,gn)),[2 3 1]);
      end
      fv  = feval(fun,x,y,P1{:});
    else
      fv  = feval(fun,x,y,P0{:});
    end
  else
    fv = feval(fun,x,y);
  end

  int(k) = squeeze(sum(sum((wfxy{ix}(:,:,ones(1,nk)).*fv)))).*qx(k).*qy(k);
    
  tol1(k)=abs(int_old(k)-int(k)); % absolute tolerance
  if ix>1
    k=find(tol1 > abs(tol*int)); %|tol1 > abs(tol));%indices to integrals which did not converge
  end  
  if any(k), % compute integrals again
    nk=length(k);%# of integrals we have to compute again
    int_old(k)=int(k);
  else
    converge='y';
    break;
  end
    
   
end
int=reshape(int,aSize); % reshape int to the same size as input arguments
tol1 = reshape(tol1,aSize);
if converge=='n',
  if nk>1
    if (nk==prod(aSize)),
      tmptxt = 'All integrals did not converge--singularities likely!';
    else
      tmptxt = sprintf('%d integrals did not converge--singularities likely!',nk);
    end
  else
    tmptxt = 'Integral did not converge--singularity likely!';
  end
  warning('GAUSSQ2D:SINGULARITY',tmptxt)
end



