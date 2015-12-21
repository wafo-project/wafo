function [r]=armacov(C,A,s2,m)

% ARMACOV   Computes the covariance function for an AR- or ARMA-model.
%   The process is governed by the system equation
%     A(q) * x(t) = C(q) * sqrt(s2) * e(t) 
%
% r = armacov(C,A,s2,n)
%
% r   = Covariance function
%
% C   = Coefficients in C-polynomials. [1 c_1 ... c_nc]
% A   = Coefficients in A-polynomials. [1 a_1 ... a_na]
% s2  = Innovation variance.
% n+1 = Number of calculated values.
%
% Example: AR(2)-process.
%   r = armacov(1,[1 1 0.9],1,50);
%   plot(0:50,r)
% Example: ARMA(4,2)-process.
%   r = armacov([1 0.05 -0.88],[1 -2.06 1.64 -0.98 0.41],4.84e-6,150);
%   plot(0:150,r)

% Copyright (c) 1997 by Pär Johannesson
% Toolbox: Rainflow Cycles for Switching Processes V.1.0, 2-Oct-1997

p=length(A)-1;
q=length(C)-1;
n=max([p q]);
A1temp=[A zeros(1,n-p)]';
C1temp=[C zeros(1,n-q)]';
A2temp=A1temp;
A1=A1temp;
A2=zeros(n+1,1);
for i=2:n+1
   A1temp=[0;A1temp];
   A2temp=[A2temp;0];
   A1=[A1 A1temp(1:n+1,1)];
   A2=[A2 A2temp(i:n+i,1)];
end
A3=A1(1:q+1,1:q+1);
C1=C1temp;
for i=2:q+1
   C1temp=[C1temp;0];
   C1=[C1 C1temp(i:n+i,1)];
end
r=(inv(A1+A2)*C1*inv(A3)*C')';
if p>0
   for i=n+2:m+1
      r=[r -r(1,i-1:-1:i-p)*A(1,2:p+1)'];
   end
else
   r=[r zeros(1,m-length(r)+1)];
end
r=s2*r(1,1:m+1);


