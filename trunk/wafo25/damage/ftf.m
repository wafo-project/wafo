function [time,F]=ftf(e,d,sigma2,sigma_D,number)
%FTF  Calculates the fatigue failure time distribution.
%
%       F(t) = P[ T^f <= t ].
%
%  CALL: [t,F] = ftf(e,d,s2,s2D,number);
%
%  where
%
%        t      = an one column matrix with times  t,
%        F      = the distribution function F(t),
%        e      = epsilon, a constant,
%        d      = the damage intensity,
%        s2     = the residual variance,
%        s2D    = the variance of the total damage,
%        number = plot parameter (optionalinput argument); 
%                 if equal to 1 the distribution function will be plotted.
%  
% Example:  
%   T = 1;
%   tp = dat2tp(load('sea.dat'));
%   RFC = tp2rfc(tp);
%   [t,F] = ftf(5.5e-10,cc2dam(RFC,5)/T,0.06,0.5);
%
% See also  cc2dam, dat2tp, tp2cc
  
% Tested on: matlab 5.3
% History:
% Revised by PJ 10-Jan-2000
%   updated for WAFO
% Original version from FAT by Mats Frendahl 
%   Copyright 1993, Mats Frendahl, Dept. of Math. Stat., University of Lund.

timefailurecenter=1/d/e; number_of_t=99;
delta=timefailurecenter/number_of_t;
time=.5*timefailurecenter:delta:1.5*timefailurecenter;
F=.5+.5*erf(log(d*time.*e)/sqrt(sigma2));

number_of_x=99; x=-4:8/number_of_x:4; phi_x=phi(x,0,1);
I=0;
for i=1:length(time)
    t=log(d*e*time(i)+e*sigma_D*sqrt(time(i))*x)./sqrt(sigma2);
    y=(.5+.5*erf(t/sqrt(2))).*phi_x;
    I(i)=trapez(x,y);
end

if nargin==5
   if number==1
      plot(time,I)
      axis([min(time) max(time) -0.1 1.1])
      title('P[ T^f <= t ]'),xlabel('t')
      axis;
   end
end

function p=phi(x,m,v,nr)
%  Evalutes the phi-/Phi-function, density/distribution function 
%  for a Gaussian variable with mean  m  and variance  v.
%
%  CALL: f = phi(x,m,v,nr)
%
%  where
%
%        f  = the density/distribution function,
%        x  = a vector of x-values,
%        m  = the mean,
%        v  = the variance,
%        nr = plot parameter  (optional input argument)
%
%             0 => f = density function,
%             1 => f = distribution function.

%  Copyright 1993, Mats Frendahl, Dept. of Math. Stat., University of Lund.

if nargin==3, nr=0; end

p=1/sqrt(2*pi*v)*exp(-0.5*(x-m).^2/v);

if (nargin==4) && (nr==1)
  p=(1+erf((x-m)./sqrt(2*v)))./2;
end

function integral=trapez(x,y)
%  Calculates an integral according to the trapezodial rule given two 
%  vectors,  x  and  y,  with  x_k-  and  y_k-values.
%
%  CALL: I = trapez(x,y)
%
%  where
%
%        x = a vector with x_k-values,
%        y = a vector with y_k-values.

%  Copyright 1993, Mats Frendahl, Dept. of Math. Stat., University of Lund.

integral=.5*(y(2:length(y))+y(1:length(y)-1))*diff(x)';
