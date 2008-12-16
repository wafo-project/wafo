function damM=damint(parm,beta,M)
% DAMINT Calculates damage intensity from counting distribution.
% 
%  CALL: D = damint(param,beta,M);
%
%
%        D     = an one column matrix with damages,
%        param = the parameter matrix,
%        beta  = an one column matrix with beta-values,
%        M     = the  nxn  matrix to be integrated.
%
% Model: 
%                  b_i
%     f(u,v) = (u-v)   ,  u > v.
% 
%
% 
%  Copyright 1993, Mats Frendahl, Dept. of Math. Stat., University of Lund.

u2=parm(2); u1=parm(1); c=(parm(3)-1)/2; delta=(u2-u1)/c/2;
[v,u]=meshgrid(u1:delta:u2,u1:delta:u2);

N=max(size(beta)); damM=zeros(1,length(beta));
z=u-v;

invbar=fliplr(triu(fliplr(z),1))+fliplr(tril(fliplr(ones(size(z))),-1));
n=max(size(invbar)); z=invbar; 

for a=1:N,
  tempsumma=0; f=invbar.^beta(a);
  tempsumma=tempsumma+f(1,1)*M(1,1);
  for i=2:n-2, tempsumma=tempsumma+M(1,i)*(f(1,i)-f(1,i-1)); end
  tempsumma=tempsumma-M(1,n-1)*f(1,n-2);
  for j=2:n-2, tempsumma=tempsumma+M(j,1)*(f(j,1)-f(j-1,1)); end
  tempsumma=tempsumma-M(n-1,1)*f(n-2,1);
  for i=2:n-1, tempsumma=tempsumma+M(i,n-i+1)*f(i-1,n-i); end
  for i=2:n-2, tempsumma=tempsumma-M(i,n-i)*f(i,n-i); end
  z=invbar;
  z=beta(a)*(beta(a)-1)*z.^(beta(a)-2);
  z=fliplr(triu(fliplr(z),1));
  N=fliplr(triu(fliplr(M(2:n-2,2:n-2))));
  damM(a)=tempsumma+sum(sum(N.*z(2:n-2,2:n-2)))*delta^2;
end

