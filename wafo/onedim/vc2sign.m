function hs=vc2sign(H,ratio)
% VC2SIGN Calculates ratio-significant value of the input vector.
%
%  CALL: hs=vc2sign(H);
%        hs=vc2sign(H,ratio);
%  where
%
%        hs    = a three column matrix with 'ratio' in the first column,
%                the ratio-significant value of  H  in the second, and 
%                a ratio-quantile in the third column.
%        H     = data vector.
%        ratio = a constant (0<ratio<1) defining the significant
%                value of the vector H. (optional input and only one
%                significant value is computed).
% 
% Example
%  x = load('sea.dat');
%  [vcf,h] = dat2steep(x); % extract Crest front velocity and wave height.  
%  hs = vc2sign(h);
%  plot(hs(:,1),hs(:,2));
%  xlabel('Ratio'); ylabel('Ratio-significant wave height.');
%  hs3 = vc2sign(h,1/3);
%  assert(hs3(2),  1.77039326808989, 1e-10);   % significant wave height  
%  S   = dat2spec(x);
%  assert(spec2char(S,'Hm0'), 1.89191906170386, 1e-10);
%
%  close all;
% 
% See also spec2char  
  
Hs=sort(H); 
test=size(Hs);

if test(1)>1
 Hs=Hs';
end

N=length(Hs);

if (nargin<2)
  hs=ones(N,3);

  hs(:,2)=flipud(cumsum(flipud(Hs')));

  hs(:,2)=hs(:,2)./flipud((1:N)');
  hs(:,1)=1-(0:N-1)'/N;
  hs(:,3)=Hs';
else
  if ratio <= 0 
    error('ratio  must be >0.')
  end

  if ratio > 1
    error('ratio  must be <=1.')
  end

  N3=ceil((1-ratio)*N); 
  hs=[ratio sum(Hs(N3:N))/(N-N3+1)];
end
