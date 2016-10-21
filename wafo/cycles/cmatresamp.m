function FF = cmatresamp(F,Method)
%CMATRESAMP Resamples a cycle matrix.
%
% CALL: FF = cmatresamp(F)
%
%   F  = Cycle matrix.                [n,n]
%   FF = Resampled cycle matrix.      [n,n]
%
% Resamles by picking the cycles at random with frequencies specified by F. 
% The result FF contains the same number of cyles as in F.
%
% Example:
%   F = round(5*triu(rand(4,4),1))
%   FF = cmatresamp(F)

% See also  cmatplot


% Check input arguments
ni = nargin;
%no = nargout;
error(nargchk(1,2,ni));

if ni<2, Method=[]; end

% Default values, vectorized calculations
if isempty(Method), Method=3; end % Vectorized calculations


[m,n]=size(F);
N=sum(sum(F));


switch Method
  % Method 1: Binomial sampling
case 1,
  F1=F(:);
  I=F1>0;
  F2=F1(I);
  FF2 = zeros(size(F2));
  NN=N;
  for i=1:length(F2)
    FF2(i) = rndbin(NN,F2(i)/NN);
    NN=NN-FF2(i);
  end
  FF1=zeros(size(F1));
  FF1(I)=FF2;
  FF=reshape(FF1,m,n);
  
  % Method 2: sum of Multinomial I
case 2
  
  F1=F(:);
  I=find(F1>0);
  F2=F1(I);
  tab = [[0;cumsum(F2)/sum(F2)] (0:length(F2))'];
  R=rand(N,1);
  x=ceil(interp1(tab,R));
  FF1=full(sparse(I(x),1,1,m*n,1));
  FF=reshape(FF1,size(F));
  
  % Method 3: sum of Multinomial II
case 3
  
  F1=F(:);
  I=find(F1>0);
  F2=F1(I);
  Fn = cumsum(F2)/sum(F2);
  FF2 = zeros(size(F2));
  for k = 1:N
    x=sum(rand>Fn)+1;
    FF2(x)=FF2(x)+1;
  end
  FF1=zeros(size(F1));
  FF1(I)=FF2;
  FF=reshape(FF1,m,n);
  
end

  
  
  