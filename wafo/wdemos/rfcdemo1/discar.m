function [P,uu] = discar(u,a1,m,s)

% DISCAR   Discretizes an AR(1) process.
%     System equation: X(k) = -a1*X(k-1) + e(k),  e(k) = N(m,s^2)
%
% P = discar(u,a1,s,m)
%
% P = Transition matrix,
%     P(i,j) = P ( u(j)<X(k)<u(j+1) |  u(i)<X(k-1)<u(i+1) )
%
% u  = Discretization levels.
% a1 = System parameter.
% m  = Innovation mean.
% s  = Innovation std.

% Copyright (c) 1997 by Pär Johannesson
% Toolbox: Rainflow Cycles for Switching Processes V.1.0, 2-Oct-1997

n = length(u);
PP = zeros(n-1,n-1);
tol =2e-1;

for i = 1:n-1
  fprintf(1,'   i = %d\n',i);

  % P(i,j) = P ( u(j)<X(k)<u(j+1) ,  u(i)<X(k-1)<u(i+1) )
  for j = 1:n-1
    PP(i,j) = quad8('f_ar',u(i),u(i+1),tol,[],u(j),u(j+1),a1,m,s);
  end

  % Pi = P ( u(i)<X(k)<u(i+1) )
%  Pi = normcdf((u(i+1)-m)*sqrt(1-a1^2)/s)-normcdf((u(i)-m)*sqrt(1-a1^2)/s);

  % P(i,j) = P ( u(j)<X(k)<u(j+1) |  u(i)<X(k-1)<u(i+1) )
%  if Pi == 0  % Prob. outside working precision
%    PP(i,:) = PP(i,:)/sum(PP(i,:));
%  else
%    PP(i,:) = PP(i,:)/Pi;
%  end
end

P = PP;
uu = (u(1:n-1)+u(2:n))/2;
