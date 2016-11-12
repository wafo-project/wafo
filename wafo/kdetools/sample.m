function s=sample(A,m,r)
%SAMPLE  Random sampling of points from a data-set
%
% CALL: s = sample(data,m,R)
%  
%  s    = sampled selection from data,  size m x D
%  data = data matrix, size N x D (D = # dimensions)
%  m    = sampling size 
%  R    = 0 sampling without replacement 
%         1 sampling with replacement (default)
% 
%  SAMPLE(DATA,M,R) selects a random sample of M data points from the
%  multivariate data-set in the matrix DATA.
%
% Example:
%     data = rndnorm(0,1,500,3);
%     s    = sample(data,100,0);
%

% History:
% revised pab dec2003  
%  changed ind generation to avoid dependence on stats-toolbox
% revised pab 10.12.1999
%  - faster sampling
% by CB kdetools

if nargin<2,
  error('Incorrect number of function parameters');
end;
if nargin<3 ||isempty(r)
 r=1;
end

n =size(A,1);

if m>n && r==0,
  error('Requested sample size too large');
end;

if m==n && r==0,
  s=A;
  return;
end;

if r==0, % Sample without replacement.
 ind = randperm(n);
else  % sample with replacement
 ind = ceil(n*rand(m,1));
end
s=A(ind(1:m),:); 
