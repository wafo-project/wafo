function [g,t0]=ochifun(y,data)
%OCHIFUN Calculates the transformation g proposed by Ochi et al.
%         Assumption: a Gaussian process, Y, is related to the
%                     non-Gaussian process, X, by Y = g(X). 
%
%  CALL:  [g, test] = ochifun(y,data);
%
%     g    = [x g(x)], a two column matrix with the transformation g(x).
%    test  = int (g(x)-x)^2 dx  where int. limits is given by X. This
%            is a measure of departure of the data from the Gaussian model.
%     x    = a row vector with x-values. 
%            (default linspace(-5*sigma,5*sigma,513)+mean)
%     data = [gamma_a gamma_b sigma mean sigma2 mean2], 
%            transformation parameters, standard deviation, and mean,
%            respectively. (default [0.1 0.15 1 0 1 0])
%            Mean2 and sigma2 are normalizing parameters in the
%            transformed world, i.e., to make the gaussian process in
%            the transformed world is N(0,1).
%
%  This is a transformation model where the transformation is chosen to
%  be a monotonic exponential function:
%
%    g(x) = ((1-exp(-gamma*(x-mean)/sigma))/gamma-mean2)/sigma2
%  where
%    gamma  = gamma_a  for x>=mean
%             gamma_b     otherwise
%
%  According to Ochi it is appropriate for a process with very strong
%  non-linear characteristics.
%
%  NOTE: - g does not have continous derivatives of 2'nd order or higher
%          unless gamma_a==gamma_b.
%
% See also  ochitr, hermitetr, lc2tr, dat2tr


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
% - changed name to ochifun
% - fixed a bug: gamma =0 is now handled correctly
% - default x is now levels([-5 5 513])*sa+ma -> better to have the
% discretization
%  represented with exact numbers, especially when calculating
%  derivatives of the transformation numerically.
% revised pab 21.02.2000
%  - added ma ma2
%  - changed name to ochitr
%  - added references
%  - changed default value for y 
%  - fixed a normalization bug

%1./data2(1:2)

data2 = [0.1 0.15 1 0 1 0];
if nargin>=2 & any(~isnan(data))
  ind=find(~isnan(data(1:min(length(data),6))));
  data2(ind)=data(ind);
end
ga     =data2(1);   gb  = data2(2); 
sigma  = data2(3);  ma  = data2(4); 
sigma2 = data2(5);  ma2 = data2(6); 
if nargin<1|isempty(y);
  y = linspace(-5*sigma+ma,5*sigma+ma,513)';
else
  y=y(:);
end
g=zeros(length(y),2);
g(:,1)=y;
igp=(y>=ma);
igm=find(~igp);

yn = (y-ma)/sigma;
if ga==0,
  g(igp,2)=yn(igp);
else
  g(igp,2)=(1-exp(-ga*yn(igp)))/ga;
end
if gb==0,
  g(igm,2)=yn(igm);
else
  g(igm,2)=(1-exp(-gb*yn(igm)))/gb;
end

g(:,2)=(g(:,2)-ma2)/sigma2;

if nargout>1
  t0 = trapz(yn,(yn-g(:,2)).^2);
end
