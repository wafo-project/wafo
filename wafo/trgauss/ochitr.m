function [g,t0]=ochitr(y,data)
%OCHITR  Estimate transformation, g, from the first 3 moments.
%
%         Assumption: a Gaussian process, Y, is related to the
%                     non-Gaussian process, X, by Y = g(X). 
%
%  CALL:  [g, test]= ochitr(x,data);
%  
%     g    = [x g(x)] a two column matrix with the transformation g(x).  
%     test = int (g(x)-x)^2 dy  where int. limits is given by Y. This
%            is a measure of departure of the data from the Gaussian model.
%     x    = vector with x-values. 
%            (default linspace(-5*sigma,5*sigma,513)+mean)
%     data = the data vector [sigma skew mean]
%            is the  standard deviation, skewness and mean
%            of the process, respectively. skew=0 for a Gaussian process.
%            (abs(skew) <=2.82)(default  [1 0.16 0])
%
%  OCHITR is a transformation model where the transformation is chosen to
%  be a monotonic exponential function, calibrated such that the first 
%  3 moments of the transformed model G(y)=g^-1(y) match the moments of
%  the true  process. However, the skewness is limited by ABS(SKEW)<2.82.
%  Information about the moments of the process can be
%  obtained by site specific data, laboratory measurements or by resort
%  to theoretical models (see spec2skew). According to Ochi it is
%  appropriate for a process with very strong non-linear characteristics. 
%
%    g(x) = ((1-exp(-gamma*(x-mean)/sigma))/gamma-mean2)/sigma2
%  where
%    gamma  = 1.28*a  for x>=mean
%             3*a     otherwise
%    mean, 
%    sigma  = standard deviation and mean, respectively, of the process.
%    mean2,
%    sigma2 = normalizing parameters in the transformed world, i.e., to
%             make the gaussian process in the transformed world is
%             N(0,1).
%
% The unknown parameters a, mean2 and sigma2 are found by solving the
% following non-linear equations:
%
%        a*(sigma2^2+mean2^2)+mean2 = 0
%           sigma2^2-2*a^2*sigma2^4 = 1
%   2*a*sigma2^4*(3-8*a^2*sigma2^2) = skew
%
% NOTE: - by specifying NaN's in the data vector default values will be used.
%       - if length(data) is shorter than the parameters needed then the
%         default values are used for the parameters not specified. 
%       - The gaussian process in the transformed world is N(0,1)
%       - g does not have continous derivatives of 2'nd order or higher.
%
% Example: Simulate a Transformed Gaussian process:
%  Hm0=7;Tp=11;
%  S = jonswap([],[Hm0 Tp]); [sk ku ma]=spec2skew(S);
%  g = ochitr([],[Hm0/4,sk,ma]); g2=[g(:,1), g(:,2)*Hm0/4];
%  ys = spec2sdat(S,15000);   % Simulated in the Gaussian world
%  xs = gaus2dat(ys,g2);      % Transformed to the real world
%
% See also  spec2skew, hermitetr, lc2tr, dat2tr


% References:
% Ochi, M.K. and Ahn, K. (1994)
%  'Non-Gaussian probability distribution of coastal waves.'
%  In Proc. 24th Conf. Coastal Engng, Vol. 1, pp 482-496
%
% Michel K. Ochi (1998),
% "OCEAN WAVES, The stochastic approach",
%  OCEAN TECHNOLOGY series 6, Cambridge, pp 255-275.


% tested on matlab 5.1
% History:
% revised pab 02.01.2001
% - default x is now levels([-5 5 513])*sa+ma -> better to have the discretization
%  represented with exact numbers, especially when calculating
%  derivatives of the transformation numerically.
% - added fmins, fzero
% - moved some code into ochifun.
% revised pab 24.05.2000
%  - removed inline object with string object
% by pab 03.03.2000

data2=[1 0.16  0]; % default values
if nargin>=2 && any(~isnan(data))
  ind=find(~isnan(data(1:min(length(data),3))));
  data2(ind)=data(ind);
end
sig1=data2(1); skew=data2(2);  ma=data2(3);

if abs(skew)>2.82842712474619,
  error('Skewness must be less than 2.82842')
end

if nargin<1||isempty(y),  y=linspace(-5*sig1+ma,5*sig1+ma,513); end

%disp('1')
phat = ochitrfit(sig1,skew,ma);
%disp('2')


%disp('3')
if nargout>1,
  [g,t0] = ochifun(y,phat);
else
  g = ochifun(y,phat);
end
return




function phat=ochitrfit(sig1,skew,ma)

if skew==0
  phat = [0, 0, sig1, ma,1 0];
  return
end
% Solve the equations to obtain the gamma parameters:
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%          a*(sig2^2+ma2^2)+ma2 = 0
%           sig2^2-2*a^2*sig2^4 = E(y^2) % =1 
%   2*a*sig2^4*(3-8*a^2*sig2^2) = E(y^3) % = skew 

% Let x = [a sig2^2 ]
% Set up the 2D non-linear equations for a and sig2^2:
eqstr='[x(2)-2.*x(1).^2.*x(2).^2-P1, 2.*x(1).*x(2).^2.*(3-8.*x(1).^2.*x(2))-P2  ]';
% Or solve the following 1D non-linear equation for sig2^2:
eqstr2 = '-sqrt(abs(x-P1)*2).*(3.*x-4*abs(x-P1))+abs(P2)';
%g3 = inline(eqstr2,2);

g2 = eqstr2;
g1 = eqstr;

v = version;  ix = find(v=='.');
vnr= str2num(v(1:ix(min(2,length(ix)))-1));

% start value for: a   sig2^2  
X0=[0.01 sig1^2];
Xi =[1 2]; % Start interval where sig2^2 is located.
if vnr>5.2
  % sol = fsolve(g1,X0,[],1,skew);%sig1.^2,skew*sig1^3);
  opt = [];%optimset('disp','iter');
  if 1,
    sig22 = fzero(g2,Xi,opt,1,skew); % smallest solution for sig22
    a  =   sign(skew)*sqrt(abs(sig22-1)/2/sig22^2);
    sol = [a sig22];
  else
    % find the solution by least squares
    g1 = [ 'sum(' g1 '.^2)'];
    sol = fminsearch(g1,X0,opt,1,skew);%sig1.^2,skew*sig1^3);
  end  
else
  % sol = fsolve(g1,X0,[],[],1,skew);%sig1.^2,skew*sig1^3);
  trace1=[];
%  if 1,
    sig22 = fzero(g2,[1 2],[],trace1,1,skew);
    a  =   sign(skew)*sqrt(abs(sig22-1)/2/sig22^2);
    sol = [a sig22];
%   else
%     % find the solution by least squares
%     g1 = [ 'sum(' g1 '.^2)'];
%     %sol = fmins(g1,X0,[],1,sig1.^2,skew*sig1^3);
%     sol = fmins(g1,X0,[],trace1,1.^2,skew*1^3);
%   end
end

a     = sol(1);
sig22 = sol(2);
sig2 = sqrt(sig22);

% Solve the following 2nd order equation to obtain ma2
%        a*(sig2^2+ma2^2)+ma2 = 0
my2 =  (-1-sqrt(1-4*a^2*sig22))./a;  % Largest mean
ma2 = a*sig22/my2 ;                  % choose the smallest mean
gam_ab = [1.28 3]*a; % this is valid for processes with very strong
                      % nonlinear characteristics
phat = [gam_ab sig1 ma sig2 ma2];


return






