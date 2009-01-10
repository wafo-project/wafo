function [k,k2,ind]=w2k(w,th,h,g)
% W2K Translates from frequency to wave number
%     using the dispersion relation
%
% CALL:  [k,k2,ind]=w2k(w,theta,h,g)
%
%    k,k2 = wave numbers
%     ind = index list of frequencies for which the routine failed (rare)
%     w   = angular frequency   
%   theta = direction           (default 0)
%     h   = water depth         (default Inf)
%     g   = constant of gravity (default see gravity)
%
% Uses Newton Raphson method to find the wave number k
% in the dispersion relation w^2= g*k*tanh(k*h).
% The solution k(w) => k = k(w)*cos(theta)
%                      k2= k(w)*sin(theta)
% The size of k,k2 is the common size of w and theta if they are matrices
% OR length(theta) x length(w) if they are vectors. If w or theta is scalar
% it functions as a constant matrix of the same size as the other. 
%
% Example
%  w = linspace(0,3);
%  plot(w,w2k(w))
%
% See also  k2w

% Tested on: Matlab 5.3  
% History: 
%  revised by es 23.05.00, made 'If w or theta...' in help really true  
%  revised by es 21.01.2000  allow length(g)=2 for h=inf (3D normalization)
%  by es, pab 13.08.99

if nargin<1||isempty(w)
  k=[];k2=[];
  return
end
if nargin<2||isempty(th)
  th=0;
end
if nargin<3||isempty(h)
  h=inf;
end
if nargin<4||isempty(g)
  g=gravity;
end

% Control of dimension of w and theta
thtype=0;
sizeth = size(th); 
if any(sizeth>1) % non-scalar
  if sizeth(1)==1 || sizeth(2)==1 % theta vector
    th=th(:);
    thtype=1;
  end
  if size(w,1)==1||size(w,2)==1 % w vector or scalar
    if thtype==1 % theta also vector
      th = th(:,ones(1,length(w)));
%      w=ones(size(th,1),1)*w(:)';
    else
      error('Input dimensions do not match')
    end
  else % both matrices
    if any(size(w)~=size(th))
      error('Input dimensions do not match')
    end    
    thtype=2;
  end
end

k=sign(w).*w.^2;  % deep water
if h==inf,
  if thtype==1 % vector
    k=ones(size(th,1),1)*k(:)';
  end
  k2=k.*sin(th)/g(end); %size np x nf
  k=k.*cos(th)/g(1);
  return
end

if length(g)>1
  error('WAFO:W2K','Finite depth in combination with 3D normalization (=> length(g)=2) is not implemented yet.')
end
% Newton's Method
% Permit no more than count_limit interations.
count_limit = 100;
count = 0;
%disp('Finding the wave numbers')

%k = find(w);
hn=zeros(size(k));
ix=find(w>0 | w<0);

% Break out of the iteration loop for three reasons:
%  1) the last update is very small (compared to x)
%  2) the last update is very small (compared to sqrt(eps))
%  3) There are more than 100 iterations. This should NEVER happen. 
while(any(ix) && count < count_limit), 

  count = count + 1;
  ki=k(ix);
  hn(ix)  = (ki.*tanh(ki.*h)-w(ix).^2/g)./(tanh(ki.*h)+ki.*h./(cosh(ki.*h).^2));
  knew = ki - hn(ix);
  % Make sure that the current guess is not zero.
  % When Newton's Method suggests steps that lead to zero guesses
  % take a step 9/10ths of the way to zero:
  ksmall = find(abs(knew) == 0);
  if any(ksmall),
    knew(ksmall) = ki(ksmall) / 10;
    hn(ix(ksmall)) = ki(ksmall)-knew(ksmall);
  end
  
  k(ix) = knew;
%   disp(['Iteration ',num2str(count),'  Number of points left:  ' num2str(length(ix)) ]),
  
  ix=find((abs(hn) > sqrt(eps)*abs(k))  &  abs(hn) > sqrt(eps));   
end

if count == count_limit,   
  warning('WAFO:W2K','W2K did not converge. \nThe maximum error in the last step was: %13.8f', max(hn(ix)) )
  ind=ix;
else
  ind=[];
  %disp('All wave numbers found')
end

if thtype==1
  k=ones(size(th,1),1)*k(:)';
end
k2=k.*sin(th);
k =k.*cos(th);
