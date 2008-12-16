function [g ,t0]=hermitetr(x,data,def)
%HERMITETR Estimate transformation, g, from the first 4 moments.
%
%           Assumption: a Gaussian process, Y, is related to the
%                     non-Gaussian process, X, by Y = g(X). 
%
%  CALL:  [g,test] = hermitetr(x,data,def);
%         [g,test] = hermitetr(x,S,def);
%
%        g    = [x g(x)] a two column matrix with the transformation g(x).
%        test = int (g(x)-x)^2 dx  where int. limits is given by X. This
%               is a measure of departure from the Gaussian model.
%        x    = a row vector with x-values. 
%               (default linspace(-5*sigma,5*sigma,501)+mean)
%        data = [sigma skew kurt mean] is the  standard deviation,
%               skewness, kurtosis and mean of the process,
%               respectively. skew=kurt-3=0 for a Gaussian process.
%               This is fairly accurate if kurt>=0 and 
%               0<=skew^2 <= 8*kurt/9  (default  [1 0.16 3.04 0])
%        S    = spectral density struct from which 
%               [sigma skew kurt mean] is calculated using spec2skew.
%        def  = 1  Winterstein et. al. (1994) parametrization (default)
%               2  Winterstein (1988) parametrization
%
% HERMITETR is a hermite transformation model where the transformation is
% chosen to be monotonic cubic polynomial, calibrated such that the first 
% 4 moments of the transformed model G(y)=g^-1(y) match the moments of
% the true process. Information about the moments of the process can be
% obtained by site specific data, laboratory measurements or by resort to
% theoretical models (see spec2skew).
%
% If kurt<3 (hardening model)
%    g(x) =  xn - c3(xn^2-1) - c4*(xn^3-3*xn) 
% where 
%   xn = (x-mean)/sigma
%   c3 = skew/6
%   c4 = (kurt-3)/24
%
% If kurt>=3 (softening model)
%    G(y) = mean + K*sigma*[ y + c3(y^2-1) + c4*(y^3-3*y) ]
%  where
%    y  = g(x) = G^-1(x)
%    K  = 1/sqrt(1+2*c3^2+6*c4^2)
% If def = 2 :
%    c3 = skew/(6*(1+6*c4))
%    c4 = [sqrt(1+1.5*(kurt-3))-1]/18 
% If def = 1 :
%    c3  = skew/6*(1-0.015*abs(skew)+0.3*skew^2)/(1+0.2*(kurt-3))
%    c4  = 0.1*((1+1.25*(kurt-3))^(1/3)-1)*c41
%    c41 = (1-1.43*skew^2/(kurt-3))^(1-0.1*(kurt)^0.8)
%
%  NOTE: - by specifying NaN's in the data vector default values will be used.
%        - if length(data) is shorter than the parameters needed then the
%          default values are used for the parameters not specified. 
%        - The gaussian process in the transformed world is N(0,1)
%
% Example: Simulate a Transformed Gaussian process:
%  Hm0=7;Tp=11;
%  S = jonswap([],[Hm0 Tp]); g=hermitetr*Hm0/4; 
%  ys = spec2sdat(S,15000);   % Simulated in the Gaussian world
%  xs = gaus2dat(ys,g);      % Transformed to the real world
%
% See also  spec2skew, ochitr, lc2tr, dat2tr

% References:
% Winterstein, S.R. (1988)
% 'Nonlinear vibration models for extremes and fatigue.'
% J. Engng. Mech., ASCE, Vol 114, No 10, pp 1772-1790
%
% Winterstein, S.R, Ude, T.C. and Kleiven, G. (1994)
% "Springing and slow drift responses:
%  predicted extremes and fatigue vs. simulation"
% In Proc. 7th International behaviour of Offshore structures, (BOSS)
% Vol. 3, pp.1-15



% tested on matlab 5.2
% History:
% revised pab dec 2003  
% revised pab 02.04.2001
% -Added Spectrum as a possible input
%revised pab 02.01.2000
% - added t0, check on that g is monotone
% - moved some code to hermitefun
% by pab 21.02.2000


data2=[1 0.16 3.04 0]; % default values
if nargin>=2 && ~isempty(data)
  switch class(data)
    case 'double',
      if  any(~isnan(data))
        ind=find(~isnan(data(1:min(length(data),4))));
        data2(ind)=data(ind);
      end
    case 'struct', % data is a spctral density struct
      [skew, kurt, ma, sigma] = spec2skew(data);
      data2 =  [sigma skew, kurt, ma];
    otherwise
      warning('WAFO:HERMITETR','Wrong input data')
  end
end
sigma=data2(1); skew=data2(2); kurt=data2(3);  ma=data2(4);
if nargin<1||isempty(x),   x   = linspace(-5*sigma,5*sigma,501)+ma; end
if nargin<3||isempty(def), def = 1;end

%skew,ga2
ga2 = kurt-3;
if ga2<0
  c4 = ga2/24;
  c3 = skew/6;
else
  switch def
    case 2, % Winterstein 1988 parametrization
      if skew^2>8*(ga2+3)/9,
        disp('warning: kurtosis too low compared to the skewness')
      end
      c4 = (sqrt(1+1.5*ga2)-1)/18; 
      c3 = skew/(6*(1+6*c4));
    otherwise,   % Winterstein et. al. 1994 parametrization intended to
      % apply for the range:  0 <= ga2 < 12 and 0<= skew^2 < 2*ga2/3
      if skew^2>2*(ga2)/3,
        disp('Warning: kurtosis too low compared to the skewness')
      end
      if (ga2 < 0)|| (12 < ga2)
        disp('Warning: kurtosis must be between 0 and 12')
      end
      c3 = skew/6*(1-0.015*abs(skew)+0.3*skew^2)/(1+0.2*ga2);
      if ga2==0,
        c4=0;
      else
        c41= (1-1.43*skew.^2./ga2).^(1-0.1*(ga2+3).^0.8);
        c4 = 0.1*((1+1.25*ga2)^(1/3)-1)*c41;
      end
  end
end
if ~isreal(c4)||~isreal(c3)
  error('Unable to calculate the polynomial')
end

if nargout>1
  [g, t0]=hermitefun([c3 c4],x(:),ma,sigma,ga2);
else
  g = hermitefun([c3 c4],x(:),ma,sigma,ga2);
end
return





