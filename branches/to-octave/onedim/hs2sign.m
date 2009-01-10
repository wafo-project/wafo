function hs=hs2sign(f,ratio)
% HS2SIGN Calculates a ratio-significant value of a histogram.
%
%  CALL: cs = hs2sign(f);
%  CALL: cs = hs2sign(f,ratio);
%
%        cs    = a three column matrix with ratio in the first column
%                the ratio-significant value of  H in the second and 
%                a ratio-quantile in the third column.
%        f     = a two column matrix containing the histogram.
%        ratio = a constant (0<ratio<1) defining the significant value,
%                (optional input, only one significant value is computed).

% Revised pab jan 2007
[y index]=sort(f(:,1));
f=f(index,:);


if (nargin<2)
  f  = flipud(f);
  Hs = cumsum(f(:,2)); 

  f=f(Hs>0,:);
  Hs=cumsum(f(:,2)); 
  Fs=flipud(cumsum(f(:,2).*f(:,1))./Hs);

  N=length(Hs);

  NN=sum(f(:,2));

  hs=ones(N,3);
  hs(:,1)=flipud(Hs)./NN;
  hs(:,2)=Fs;
  hs(:,3)=flipud(f(:,1));
else
  if ratio <= 0 
    error('ratio  must be >0.') 
  end
  if ratio >= 1
    error('ratio  must be <1.')
  end
  Hs=cumsum(f(:,2)); 
  hs=min(Hs./sum(f(:,2)),ones(length(Hs),1).*ratio);
  dhs=abs(diff([0; f(:,1)]'));
  hs=[ratio dhs*hs/ratio];
end
