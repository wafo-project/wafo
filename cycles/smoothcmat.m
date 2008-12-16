function [Fsmooth,h] = smthcmat(F,method,h,NOsubzero,alpha)
%SMOOTHCMAT Smooth a cycle matrix using (adaptive) kernel smoothing
%
% CALL:  Fsmooth = smoothcmat(F,method);
%        Fsmooth = smoothcmat(F,method,[],NOsubzero);
%        Fsmooth = smoothcmat(F,2,h,NOsubzero,alpha);
%
% Input:
%        F       = Cycle matrix.           [nxn]
%        method  = 1: Kernel estimator (constant bandwidth). (Default)
%                  2: Adaptiv kernel estimator (local bandwidth). 
%        h       = Bandwidth (Optional, Default='automatic choice')
%      NOsubzero = Number of subdiagonals that are zero
%                  (Optional, Default = 0, only the diagonal is zero)
%        alpha   = Parameter for method (2) (Optional, Default=0.5).
%                  A number between 0 and 1.
%                  alpha=0 implies constant bandwidth (method 1).
%                  alpha=1 implies most varying bandwidth.
%
% Output:
% F       = Smoothed cycle matrix.   [nxn]
% h       = Selected bandwidth.
%
% See also  cc2cmat, tp2rfc, tp2mm, dat2tp

% Tested on Matlab 5.3
%
% History:
% Revised by PJ 18-May-2000
%   Updated help text.
% Revised by PJ  01-Nov-1999
%   updated for WAFO
% Created by PJ (Pär Johannesson) 1997
%   from 'Toolbox: Rainflow Cycles for Switching Processes V.1.0'

% Check input arguments

ni = nargin;
no = nargout;
error(nargchk(1,5,ni));

if ni<2, method=1; end
if ni<3, h=[]; end
if ni<4, NOsubzero=[]; end
if ni<5, alpha=[]; end

if method == 1 | method == 2  % Kernel estimator
  if isempty(h), aut_h=1; else aut_h=0; end
  if isempty(NOsubzero), NOsubzero=0; end
else
  error('Input argument "method" should be 1 or 2');
end

if method == 2   % Adaptive Kernel estimator
  if isempty(alpha), alpha=0.5; end
end

n = length(F);    % Size of matrix
N = sum(sum(F));  % Total number of cycles

Fsmooth = zeros(n,n);

if method == 1 | method == 2 % Kernel estimator

  d = 2;   % 2-dim
  [I,J] = meshgrid(1:n,1:n);

  % Choosing bandwidth
  % This choice is optimal if the sample is from a normal distr.
  % The normal bandwidth usualy oversmooths,
  % therefore we choose a slightly smaller bandwidth
  
  if aut_h == 1
    h_norm = smoothcmat_hnorm(F,NOsubzero);
    h = 0.7*h_norm;         % Don't oversmooth 
    
    %h0 = N^(-1/(d+4));
    %FF = F+F';
    %mean_F = sum(sum(FF).*(1:n))/N/2;
    %s2 = sum(sum(FF).*((1:n)-mean_F).^2)/N/2;
    %s = sqrt(s2);       % Mean of std in each direction
    %h_norm = s*h0;      % Optimal for Normal distr.
    %h = h_norm;         % Test
  end

  % Calculating kernel estimate
  % Kernel: 2-dim normal density function
  
  test=0;
  if test==0
    for i = 1:n-1
      for j = i+1:n
        if F(i,j) ~= 0
          F1 = exp(-1/(2*h^2)*((I-i).^2+(J-j).^2));  % Gaussian kernel
          F1 = F1+F1';                     % Mirror kernel in diagonal
          F1 = triu(F1,1+NOsubzero);       % Set to zero below and on diagonal
          F1 = F(i,j) * F1/sum(sum(F1));   % Normalize
          Fsmooth = Fsmooth+F1;
        end
      end
    end
  else
    
    [II,JJ]=find(F>0);
    for k = 1:length(II)
      i=II(k); j=JJ(k);
      F1 = exp(-1/(2*h^2)*((I-i).^2+(J-j).^2));  % Gaussian kernel
      F1 = F1+F1';                     % Mirror kernel in diagonal
      F1 = triu(F1,1+NOsubzero);       % Set to zero below and on diagonal
      F1 = F(i,j) * F1/sum(sum(F1));   % Normalize
      Fsmooth = Fsmooth+F1;
    end
    
  end
  
  
end

if method == 2
  Fpilot = Fsmooth/N;
  Fsmooth = zeros(n,n);
  [I1,I2] = find(F>0);
  logg = 0;
  for i =1:length(I1)
    logg = logg + F(I1(i),I2(i)) * log(Fpilot(I1(i),I2(i)));
  end
  g = exp(logg/N);
  lambda = (Fpilot/g).^(-alpha);

  for i = 1:n-1
    for j = i+1:n
      if F(i,j) ~= 0
        hi = h*lambda(i,j);
        F1 = exp(-1/(2*hi^2)*((I-i).^2+(J-j).^2));  % Gaussian kernel
        F1 = F1+F1';                     % Mirror kernel in diagonal
        F1 = triu(F1,1+NOsubzero);       % Set to zero below and on diagonal
        F1 = F(i,j) * F1/sum(sum(F1));   % Normalize
        Fsmooth = Fsmooth+F1;
      end
    end
  end

end

