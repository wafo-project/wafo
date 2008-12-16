function [val, err,ier,extime]= prbnormndpc(rho,a,b,tol,useSimpson,useBreakPoints)
%PRBNORMNDPC Multivariate Normal probabilities with product correlation
% 
%  CALL:   [value, error] = prbnormndpc(rho,a,b,tol)
%  
%  value = value of integral
%  error = estimated absolute error
%  rho  = vector defining the correlation structure, i.e., 
%          corr(Xi,Xj) = rho(i)*rho(j) for i~=j
%                      = 1             for i==j  
%             -1 <= rho <= 1  
%  a,b   = lower and upper integration limits respectively.  
%  tol   = requested absolute tolerance
%
% PRBNORMNDPC calculates multivariate normal probability
% with product correlation structure for rectangular regions.
% The accuracy is up to almost double precision, i.e., about 1e-14.
%  
% Example:
%  rho2 = rand(1,2); 
%  a2   = zeros(1,2);
%  b2   = repmat(inf,1,2);
%  tol = 1e-4;  
%  [val2,err2] = prbnormndpc(rho2,a2,b2,tol);
%  g2 = inline('0.25+asin(x(1)*x(2))/(2*pi)');
%  E2 = g2(rho2)  % exact value
% 
%  rho3 = rand(1,3); 
%  a3   = zeros(1,3);
%  b3   = repmat(inf,1,3);
%  tol = 1e-4;  
%  [val3,err] = prbnormndpc(rho3,a3,b3,tol);  
%  g3 = inline('0.5-sum(sort(acos([x(1)*x(2),x(1)*x(3),x(2)*x(3)])))/(4*pi)');
%  E3 = g3(rho3)   %  Exact value
%
% See also  prbnormtndpc, prbnormnd, rind
  
%Reference
% P. A. Brodtkorb (2004), 
% Evaluating multinormal probabilites with product correlation structure.
% In Lund university report series
% and in the Dr.Ing thesis: 
% The probability of Occurrence of dangerous Wave Situations at Sea.
% Dr.Ing thesis, Norwegian University of Science and Technolgy, NTNU,
% Trondheim, Norway.
  
% History 
% revised pab April 2008
% -renamed from mvnormpcprb to prbnormndpc
% by pab 10.06.2003  
if nargin<4||isempty(tol)
  tol = 1e-4;
end
if nargin<5 || isempty(useSimpson)
  useSimpson = 1;
end
if nargin<5 || isempty(useBreakPoints)
  useBreakPoints = 0;
end

reltol = tol; %*1e-1;


  
a   = a(:);
b   = b(:);
rho = rho(:);



t1 = clock;
% Call fortran implementation
[val,err,ier] = mvnprodcorrprbmex(rho,a,b,tol,reltol,useBreakPoints,useSimpson);

extime = etime(clock,t1);
if (ier>0)
  
  warning('WAFO:prbnormndpc','Abnormal termination ier = %d\n\n%s',ier,getErrorMessage(ier))
end
return
% Old matlab calls kept just in case
trace = 0;
x_trace = [];
y_trace = []; 

if any(b<=a)
  val = 0;
  err = 0;
  return
end
%remove insignificant variables
ind = find(((a==-inf) &  (b == inf)));
if any(ind)
  rho(ind) = [];
  a(ind)   = [];
  b(ind)   = [];
 % base(ind) = [];
  if isempty(a)
    val = 1;
    err = 0;
    return
  end
end

val = min(cdfnorm(b)-cdfnorm(a));
if length(a)==1
  err = eps;
  return
end

den = sqrt(1-rho).*sqrt(1+rho);
%ind = find(base);
%if any(ind)
%  r0 = rho(ind);
%  r = abs(r0);
%  den(ind) = sqrt(r).*sqrt(2-r);
%  rho(ind) = (2*(r0>=0)-1).*(1-r); 
%end



val0 = 1;
ind = find(rho==0);
if any(ind)
  val0    = prod(cdfnorm(b(ind))-cdfnorm(a(ind)));
  a(ind)  = [];
  b(ind)  = [];
  rho(ind)= [];
  den(ind) = [];
  if (length(a)==0),
    val = val0;
    err = eps;
    return
  end
end
xCutOff    = abs(max(invnorm(0.05*reltol),-37));
xlo   = -xCutOff; %abs(max(invnorm(0.0001*reltol),-37))
xup   = -xlo;

[xlo, xup] = c1c2(xlo,xup,rho,a,b,xCutOff,den);

ind = find(abs(abs(rho)-1)<=eps);
if any(ind)
  a(ind) = [];
  b(ind) = [];
  rho(ind)= [];
  den(ind) = [];
  if length(a)==0
    val = (cdfnorm(xup)-cdfnorm(xlo))*val0;
    err = sqrt(eps);
    return
  end
end



limits = [xlo,xup];
isLimitsNarrowed = ((-7 < xlo)  |   (xup<7));
if (isLimitsNarrowed & useBreakPoints)
  
  %xCut = max(abs(xlo),abs(xup))
  xCut = 2*min(den);
  %[xlo1,xup1] = c1c2(xlo,xup,rho,a,b,xCut,den);
  %limits = unique([limits xlo1, xup1]);
  [xlo1,xup1] = c1c2(xlo,xup,rho,a,b,0,den);
  [xlo2,xup2] = c1c2(xlo,xup,rho,a,b,-xCut,den);
  [xlo3,xup3] = c1c2(xlo,xup,rho,a,b,xCut,den);
  limits = unique([limits xlo1, xup1 xlo2, xup2 xlo3,xup3])
  
  %x = linspace(xlo,xup,256);
  %f = intfun(x,rho,a,b,den);
  %[fMax,ixMax] = max(f);
% $$$   if 0 %(fMax>0),
% $$$     limits = unique([limits,x(ixMax)]);
% $$$     indLo = find(f(1:ixMax)<=0);
% $$$     if any(indLo)
% $$$       ind  = find(x(indLo(end))<limits);
% $$$       limits = unique([x(indLo(end)) limits(ind)]);
% $$$     end
% $$$     indUp = find(f(ixMax+1:end) <=0);
% $$$     if any(indUp)
% $$$       ind = find(limits<x(ixMax+indUp(1)));
% $$$       limits = unique([limits(ind) x(ixMax+indUp(1))]);
% $$$     end
% $$$   end
end

%limits
%limits = unique((limits+2)-2);
%ix = find(sqrt(eps)*10 < diff(limits));
%limits = [limits(1) limits(ix+1)]
%limits
if (useSimpson)
  [val,err,ier] = adaptivesimpson2('intfun', xlo,xup,reltol,rho,a,b,den);
else
  maxSubdivisions = 100;
  [val, err,neval,ier,alist,blist,rlist,elist,pts,iord, ...
   level,ndin, last] = ...
      dqagpe('intfun',xlo,xup,limits(2:end-1),reltol,0,...
	     maxSubdivisions,rho,a,b,den);
end

val = val*val0;
err = err*val0;

return
 

  
function [xMin,xMax ] = c1c2(xMin,xMax,rho,a,b,xCutOff,den)
%C1C2 uses the regression equation to limit the integration limits
%
  
if nargin<7
  den  =  sqrt(1-rho.^2);
end
%xCutOff = xMax;
k = find(rho>0);
if any(k)
  xMax = max(xMin, min(xMax,min((b(k)+den(k)*xCutOff)./rho(k))));
  xMin = min(xMax, max(xMin,max((a(k)-den(k)*xCutOff)./rho(k)))) ;
end
k1 = find(rho<0);
if any(k1)
  xMax = max(xMin,min(xMax,min((a(k1)-den(k1)*xCutOff)./rho(k1))));
  xMin = min(xMax,max(xMin,max((b(k1)+den(k1)*xCutOff)./rho(k1))));
end

return


function msg = getErrorMessage(ier)
  msg = '';
switch ier
  case 1,
 msg  = sprintf('%s\n',...
		'maximum number of subdivisions allowed',...
		'has been achieved. one can allow more',...
		'subdivisions by increasing the value of',...
		'limit (and taking the according dimension',...
		'adjustments into account). however, if',...
		'this yields no improvement it is advised',...
		'to analyze the integrand in order to',...
		'determine the integration difficulties. if',...
		'the position of a local difficulty can be',...
		'determined (i.e. singularity,',...
		'discontinuity within the interval), it',...
		'should be supplied to the routine as an',...
		'element of the vector points. If necessary',...
		'an appropriate special-purpose integrator',...
		'must be used, which is designed for',...
		'handling the type of difficulty involved.');
 case  2,
  msg = sprintf('%s\n',...
		'the occurrence of roundoff error is',...
		'detected, which prevents the requested',...
		'tolerance from being achieved.',...
		'the error may be under-estimated.');
 case 3,
  msg = sprintf('%s\n',...
		'extremely bad integrand behaviour occurs',...
		'at some points of the integration interval.');
  
 case 4,
  msg = sprintf('%s\n',...
		'the algorithm does not converge.',...
		'roundoff error is detected in the',...
		'extrapolation table. it is presumed that',...
		'the requested tolerance cannot be',...
		'achieved, and that the returned result is',...
		 'the best which can be obtained.');
 case 5,
  msg = sprintf('%s\n',...
		'the integral is probably divergent, or',...
		'slowly convergent. It must be noted that',...
		'divergence can occur with any other value of ier>0.');
 case 6,
  msg = sprintf('%s\n',...
		'the input is invalid because:',...
		 '1) npts2 < 2',...
		 '2) break points are specified outside the integration range',...
		 '3) (epsabs<=0 and epsrel<max(50*rel.mach.acc.,0.5d-28))',...
		 '4) limit < npts2.');
end
