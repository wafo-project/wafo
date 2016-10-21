function value = prbnorm2d(lower,upper,correl)
%PRBNORM2D Bivariate Normal probability   
%
%  CALL prb = prbnorm2d(a,b,r)
%
% prb = computed probability
% a   = vector with lower integration limits
% b   = vector with upper integration limits
% r   = correlation coefficient  
%
% PRBNORM2D computes bivariate normal probability, i.e.,
%  the probability Prob(A(1) <= X1 <= B(1) and A(2) <= X2 <= B(2)) 
%  with an absolute error less than 1e-15.
%
% Example
%  a = [-1 -2];
%  b = [1 1]; 
%  r = 0.3;
%  prbnorm2d(a,b,r)
%
% See also  cdfnorm2d, cdfnorm, prbnormndpc


% revised pab April 2008
% -renamed from binormprb to prbnorm2d
error(nargchk(3,3,nargin))
value = 1;
infinity = 37;
if all(lower <= -infinity & upper>= infinity),
  return, 
end  
if any(lower>=upper),
  value=0; 
  return,
end

infin =[2; 2]-(upper(:)>infinity)-2*(lower(:)<-infinity);

if all(infin==2)
  value =  bvd ( lower(1), lower(2), correl )...
      - bvd ( upper(1), lower(2), correl )...
      - bvd ( lower(1), upper(2), correl )...
      + bvd ( upper(1), upper(2), correl );
elseif ( infin(1) == 2  && infin(2) == 1 ) 
  value =  bvd ( lower(1), lower(2), correl )...
      - bvd ( upper(1), lower(2), correl );
elseif ( infin(1) == 1  && infin(2) == 2 ) 
  value =  bvd ( lower(1), lower(2), correl )....
      - bvd ( lower(1), upper(2), correl );
elseif ( infin(1) == 2  && infin(2) == 0 ) 
  value =  bvd ( -upper(1), -upper(2), correl )...
      - bvd ( -lower(1), -upper(2), correl );
elseif ( infin(1) == 0  && infin(2) == 2 ) 
  value =  bvd ( -upper(1), -upper(2), correl )....
      - bvd ( -upper(1), -lower(2), correl );
elseif ( infin(1) == 1  && infin(2) == 0 ) 
  value =  bvd ( lower(1), -upper(2), -correl );
elseif ( infin(1) == 0  && infin(2) == 1 ) 
  value =  bvd ( -upper(1), lower(2), -correl );
elseif ( infin(1) == 1  && infin(2) == 1 ) 
  value =  bvd ( lower(1), lower(2), correl );
elseif ( infin(1) == 0  && infin(2) == 0 ) 
  value =  bvd ( -upper(1), -upper(2), correl );
end 
return

function val = bvd(lo,up,r)
  val = cdfnorm2d(-lo,-up,r);
return