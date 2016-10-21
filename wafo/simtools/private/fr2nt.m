function NT = fr2nt(f)
%FR2NT  Calculates the counting distribution given the frequency matrix.
%   
%  CALL: NT = fr2nt(fr);
%
%  where
%
%        NT = a square counting distribution matrix for a cycle count,
%        fr = a square frequency matrix for a cycle count.

%  Copyright 1993, Mats Frendahl, Dept. of Math. Stat., University of Lund.

[n m]=size(f);
if (n==m) && (n>2)
   m1=cumsum(cumsum(f,2)-f);
   m2=zeros(n,n);
   m2(2:n-1,2:n-1)=m1(1:n-2,2:n-1);
   NT=fliplr(triu(fliplr(m2),0));
else
   disp('   The matrix is not square or dimension < 3.')
   disp('   Program will terminate.')
end
