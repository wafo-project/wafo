function [Nup,Ndown,u] = load2ud(param,L)
%  LOAD2UD Calculates the number of up/downcrossings of a load.
%  
%  CALL: [Nup,Ndown,u] = udcross([a b n],L);
%
%  where
%   
%        Nup   = the number of upcrossings of level  u_i,
%        Ndown = the number of downcrossings of level  u_i,
%        u     = a one column matrix with levels  u_i,
%        a,b   = limits for  u,  a <= u_i <= b,
%        n     = the number of levels  u_i,
%        L     = a two column matrix with times in the first column and 
%                load values in the second.
%
%  Note that the parameter matrix can be used as input.

a=param(1); b=param(2); n=param(3);

if a>b, disp('   a>b, program will terminate.'), end
if n<1,  disp('   n<1, program will terminate.'),  end

if n == 1
  u = a;
else
  u=levels(param);
end

Nup=zeros(1,n);
Ndown=zeros(1,n);

N = size(L,1);

for i = 1:n
  index = (L(1:N-1,2)<u(i)) & (L(2:N,2)>u(i));
  Nup(i) = sum(index);

  index = (L(1:N-1,2)>u(i)) & (L(2:N,2)<u(i));
  Ndown(i) = sum(index);
end
