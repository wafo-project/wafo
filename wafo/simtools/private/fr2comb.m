function comb=fr2comb(f,r)
%FR2COMB Combination matrix for placing out cycles given the frequency matrix

%  Given the frequency matrix of a cycle count the function 
%  calculates the combination matrix for placing out cycles

%  Copyright 1993, Mats Frendahl & Igor Rychlik,
%  Dept. of Math. Stat., University of Lund.

%nres=length(r);

N=length(f); nt=fr2nt(f); comb=zeros(N,N);

for i=1:N
  for j=1:N-i+1
    comb(i,j)=2*nt(i,j)+sum(f(1:i-1,j))+sum(f(i,1:j-1));
  end
end

for k=1:length(r)-1
  i=r(k); j=r(k+1);
  if ~isempty(j+1:i) %i>j+1 
    M=j+1:i; m=N+1-M;
    comb(M,m)=comb(M,m)+1;
  elseif ~isempty(i:j-1) %i<j-1
    M=i:j-1; m=N+1-M;
    comb(M,m)=comb(M,m)+1;
  end
end
 
comb=fliplr(triu(fliplr(comb),1));
