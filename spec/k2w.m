function [w,th]=k2w(k,k2,h,g,u,u2)
% K2W Translates from wave number to frequency
%     using the dispersion relation
%
% CALL:  [w,theta]=k2w(k,k2,h,g)
%
%     w   = angular frequency (rad/s) 
%   theta = direction         (rad)   
%     k   = wave numbers      (rad/m)          
%     k2  = second dimension wave number      (default 0)
%     h   = water depth       (m)             (default Inf)
%     g   = constant of gravity               (default see gravity)
%
% Dispersion relation:
%    w     = sqrt(g*K*tanh(K*h))   (  0 <   w   < inf)
%    theta = atan2(k2,k)           (-pi < theta <  pi)
% where
%     K = sqrt(k^2+k2^2)
%
% The size of w,theta is the common size of k and k2 if they are matrices,
% OR length(k2) x length(k) if they are vectors. If k or k2 is scalar
% it functions as a constant matrix of the same size as the other. 
% 
% See also  w2k

% secret options:
%     u   = current velocity                  (default 0)
%     u2  = second dimension current velocity (default 0)
% note: when u~=0 | u2~=0 then th is not calculated correctly


% Tested on: Matlab 5.3, 5.2
% History: 
%  revised by es 25.05.00, made 'If k or k2...' in help really true  
% revised by es 24.01.2000    
% revised pab 15.02.2000  - added current u,u2
%by es 13.08.99


if nargin<1||isempty(k)
  w=[];th=[];
  return
end
if nargin<2||isempty(k2)
  k2=0;
end
if nargin<3||isempty(h)
  h=inf;
end
if nargin<4||isempty(g)
  g=gravity;
end

if nargin<5||isempty(u)
  u=0;
end
if nargin<6||isempty(u2)
  u2=0;
end
% Control of dimension of k and k2
ktype=0;
if numel(k2)>1 % non-scalar
  if size(k2,1)==1||size(k2,2)==1 % k2 vector
    k2=k2(:);
    ktype=1;
  end
  if size(k,1)==1||size(k,2)==1 % k vector or scalar
    if ktype==1 % k2 also vector
      k2=k2(:,ones(1,length(k))); 
      k=k(:);
      k=k(:,ones(1,size(k2,1)))';
    else
      error('Input dimensions do not match')
    end
  else % both matrices
    if any(size(k)~=size(k2))
      error('Input dimensions do not match')
    end   
  end
end

ku=k.*u;
ku2=k2.*u2;

if nargout>1
  th=atan2(k2,k);
end
k=sqrt(k.^2+k2.^2);
w=zeros(size(k));
ix=find(k>0);
if length(ku2)==1 && length(k)>1
  ku2=ku2*ones(size(k));
end
if any(ix)
  w(ix)=ku(ix)+ku2(ix)+sqrt(g*k(ix).*tanh(k(ix)*h));
end
iy=find(w<0);
if any(iy)
  disp('Warning: waves and current are in opposite directions')
  disp('         making some of the frequencies negative.')
  disp('         Here we are forcing the negative frequencies to zero')
  w(iy)=0; % force w to zero 
end
