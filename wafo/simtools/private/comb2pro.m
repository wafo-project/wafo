function prob=comb2pro(combinations,freq)
% COMB2PRO Probability vector of placing a cycle on a specific place given all combinations
% 
%  Calculates the probability vector of placing a cycle on a specific place 
%  given the number of possible places (combinations).

%  Copyright 1993, Mats Frendahl & Igor Rychlik, 
%  Dept. of Math. Stat., University of Lund.

denominator=combinations+freq-1;
add=(denominator==0);
%denominator=denominator+add;
%nominator=combinations+add-1;
%addcomb=(combinations==0);
%lambda=freq./(combinations+addcomb);
%pp=nominator./denominator;
%pp=1./(1+lambda);
pp=(1-ones(length(combinations+add),1)./(combinations+add)).^freq;
%pp=exp(-lambda);
pp=pp';

tmp=zeros(length(pp),1); 
%prob=zeros(length(pp),1);
tmp(1)=1;
for i=2:length(pp), tmp(i)=tmp(i-1)*pp(i-1); end
prob=[tmp(length(tmp))*pp(length(pp)) (1-pp).*tmp'];


