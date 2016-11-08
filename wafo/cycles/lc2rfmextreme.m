function [Frfc,u,Nrfc,Nrfc0]=lc2nt(lc,vect)
%LC2RFMEXTREME Compute extreme RFM from level crossings.
%
% CALL: Frfc = lc2nt(lc)
%       [Frfc,u,Nrfc,Nrfc0] = lc2nt(lc)
%
%  lc     = Level upcrossing intensity.                    [nx2]
%
%  Frfc   = Extreme RFM.                                   [nxn]
%  u      = Discrete levels.                               [1xn]
%  Nrfc   = Cumulative extreme RFM.                        [nxn]
%  Nrfc0  = Cumulative extreme RFM, before correction.     [nxn]
%
% Computes the 'extreme RFM', an approximation of the RFM which 
% is good for large cycles (low minima and high maxima). 
% The approximation is based on the assumtion that level crossings
% of high and low levels are two independent Poisson processes. 
% This results in a simple expression for the intensity of interval 
% crossings, which is the same as the cumulative RFM.
%
% Example: % Gaussian process
%   u = (-5:0.2:5)'; lc = [u exp(-u.^2/2)];
%   [Frfc,u,Nrfc] = lc2rfmextreme(lc);
%   cmatplot(u,u,{Frfc Nrfc},3);
%
%   close all;
%
% See also  cmat2extralc, rfmextrapolate, tp2rfc, cc2cmat, dtp2rfm

% References:
%
%   Johannesson, P., and Thomas, J-.J. (2000): 
%   Extrapolation of Rainflow Matrices. 
%   Preprint 2000:82, Mathematical statistics, Chalmers, pp. 18. 

% Tested  on Matlab  5.3
%
% History:
% Created by PJ (P�r Johannesson) 14-Feb-2000
% Revised by PJ 20-Jul-2000
% Revised by PJ 03-Oct-2000
%   Vectorized some calculations. Now much faster.
%   vect = 0 : Not vectorized.
%          1 : Vectorized.

% Check input arguments
ni = nargin;
no = nargout;
error(nargchk(1,2,ni));

if ni<2, vect = []; end

% Default values, vectorized calculations
if isempty(vect), vect=1; end % Vectorized calculations

% Treat input data
n = size(lc,1); % Number of rows = number of levels
u = lc(:,1);    % Discrete levels
lc = lc(:,2);   % intensity of level crossings

% Compute the cumulative RFM
if vect == 0 % Not vectorized calculations
  
  Nrfc0 = zeros(n,n);
  for i = 2:n
    for j = i-1:n-1
      if (lc(i)+lc(j)) ~= 0
        Nrfc0(i,j) = lc(i)*lc(j)/(lc(i)+lc(j));
      end
    end
  end
  
else % Vectorized calculations

  [LCi,LCj] = meshgrid(lc,lc);
  LCi=triu(LCi);
  LCj=triu(LCj);
  I = (LCi+LCj~=0);
  Nrfc0 = zeros(n,n);
  Nrfc0(I) = LCi(I).*LCj(I)./(LCi(I)+LCj(I));
  
end

% Convert: Nrfc0  -->  Frfc
if vect == 0 % Not vectorized calculations
  
  Frfc = zeros(n,n);
  for i= 1:n-1
    for j= i+1:n
      Frfc(i,j) = Nrfc0(i+1,j-1) - Nrfc0(i,j-1) - Nrfc0(i+1,j) + Nrfc0(i,j);
    end
  end
  
else % Vectorized calculations
  
  Frfc=nt2cmat(Nrfc0);
  Frfc=triu(Frfc,1);
  
end

%keyboard

% Set negative elements to zero.
I = Frfc<0;
if any(I(:)) 
  %warning(['Negative elements in calculated rainflow matrix Frfc. Setting to zero!']);
  Frfc(I) = 0;
end

% Convert: Frfc  -->  Nrfc
Nrfc = cmat2nt(Frfc); % Calculate cumulative RFM (rainflow counting intensity)


