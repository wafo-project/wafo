function [Pest,N,Ni] = estmc(z,r)
% ESTMC  estimates transition matrix P from a time series of a Markov chain.
%
% [Pest,N] = estmc(z,r);
%
% z    = time series of a Markov chain.    [1xT]
% r    = number of states in Markov chain.
%
% Pest = Estimated transition matrix.
% N    =


T = length(z);

N = zeros(r,r);
Pest = zeros(r,r);

for i = 1:r
  for j=1:r
    N(i,j) = sum( (z(1:T-1)==i) & (z(2:T)==j) );
%    if i==j
%      N(i,i) = sum(z(1:T-1)==i);
%    else
%      N(i,j) = sum((z(2:T)-z(1:T-1))==j-i);
%    end
  end
end

Ni = sum(N,2);

for i = 1:r
  if Ni(i) > 0
    Pest(i,:) = N(i,:)/Ni(i);
  else
    Pest(i,:) = ones(1,r)/r;
  end
  Pest(i,i) = 0;
  Pest(i,i) = 1-sum(Pest(i,:));
end



