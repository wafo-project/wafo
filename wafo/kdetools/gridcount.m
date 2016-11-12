function c = gridcount(data,X)
%GRIDCOUNT D-dimensional histogram using linear binning.
%
% CALL  c = gridcount(data,X)
%
% c    = gridcount, size N x 1           if D = 1 
%                   size N x N x ... x N if D > 1  
% data = row vectors with D-dimensional data, size Nd x D
% X    = column vectors defining discretization, size N x D
%        Must include the range of the data.
% GRIDCOUNT obtains the grid counts using linear binning.
% There are 2 strategies: simple- or linear- binning.
% Suppose that an observation occurs at x and that the nearest point
% below and above is y and z, respectively. Then simple binning strategy
% assigns a unit weight to either y or z, whichever is closer. Linear
% binning, on the other hand, assigns the grid point at y with the weight
% of (z-x)/(z-y) and the gridpoint at z a weight of (y-x)/(z-y).
%  
% In terms of approximation error of using gridcounts as pdf-estimate,
% linear binning is significantly more accurate than simple binning.  
%
% NOTE: The interval [min(X);max(X)] must include the range of the data.
%       The order of C is permuted in the same order as 
%       meshgrid for D==2 or D==3.  
%   
% Example
%  N     = 500;
%  data  = rndray(1,N,1);
%  x = linspace(0,max(data)+1,50).';  
%  dx = x(2)-x(1);  
%  c = gridcount(data,x);
%  plot(x,c,'.')   % 1D histogram
%  plot(x,c/dx/N)  % 1D probability density plot
%  % trapz(x,c/dx/N)   
%
%  close all;
%
% See also  bincount, kdebin
  
%Reference
%  Wand,M.P. and Jones, M.C. (1995) 
% 'Kernel smoothing'
%  Chapman and Hall, pp 182-192  
  
%  History
% revised pab Jan2007
% -replaced old calls with accumarray many places.
% revised pab Feb2005
% -fixed a bug for d>2;  
% by pab Dec2003  
  
error(nargchk(2,2,nargin))

[n,d]    = size(data);  
[inc,d1] = size(X);

if d~=d1

    error('Dimension 2 of data and X do not match.')
end

dx  = diff(X(1:2,:),1);
xlo = X(1,  :);
xup = X(end,:);

if any(min(data,[],1)<xlo) || any(xup<max(data,[],1))
  error('WAFO:KDETOOLS:GRIDCOUNT','X does not include whole range of the data!')
end

if d>1
  csiz = inc(:,ones(1,d));
  n1   = n;
else
  csiz = [inc,1];
  n1   = 1;
end
  
binx = floor((data-xlo(ones(n1,1),:))./dx(ones(n1,1),:))+1;
w    = prod(dx);
if  d==1,
  
  try
    c = (accumarray(binx,(X(binx+1)-data),[inc,1]) + ...
      accumarray(binx+1,(data-X(binx)),[inc,1]))/w;

%     % binc is much faster than sparse for 1D data
%     % However, for N-Dimensional data sparse is sometimes faster.
%     c = zeros(csiz);
%     [len,bin,val] = bincount(binx,(X(binx+1)-data));
%     c(bin) = val;
%     [len,bin,val] = bincount(binx+1,(data-X(binx)));
%     c(bin) = c(bin)+val;
%     c = c/w;
   catch
    
    c = full(sparse(binx,1,(X(binx+1)-data),inc,1)+...
      sparse(binx+1,1,(data-X(binx)),inc,1))/w;
  end
elseif d==2
  b2 = binx(:,2);
  b1 = binx(:,1);
  try
    c = (accumarray([b1,b2] ,abs(prod(([X(b1+1,1) X(b2+1,2)]-data),2)),[inc,inc])+...
      accumarray([b1+1,b2]  ,abs(prod(([X(b1,1) X(b2+1,2)]-data),2)),[inc,inc])+...
      accumarray([b1  ,b2+1],abs(prod(([X(b1+1,1) X(b2,2)]-data),2)),[inc,inc])+...
      accumarray([b1+1,b2+1],abs(prod(([X(b1,1) X(b2,2)]-data),2)),[inc,inc]))/w;
  catch
    c = full(sparse(b1,b2 ,abs(prod(([X(b1+1,1) X(b2+1,2)]-data),2)),inc,inc)+...
      sparse(b1+1,b2  ,abs(prod(([X(b1,1) X(b2+1,2)]-data),2)),inc,inc)+...
      sparse(b1  ,b2+1,abs(prod(([X(b1+1,1) X(b2,2)]-data),2)),inc,inc)+...
      sparse(b1+1,b2+1,abs(prod(([X(b1,1) X(b2,2)]-data),2)),inc,inc))/w;
  end
  
else % d>2
 useSparse = 0;
  Nc = prod(csiz);
  c  = zeros(Nc,1);
 
  fact2 = [0 inc*(1:d-1)];
  fact1 = [1 cumprod(csiz(1:d-1))];
  fact1 = fact1(ones(n,1),:);
  for ir=0:2^(d-1)-1,
    bt0(:,:,1) = bitget(ir,1:d);
    bt0(:,:,2) = ~bt0(:,:,1);
    
    for ix = 0:1
      one = mod(ix,2)+1;
      two = mod(ix+1,2)+1;
      % Convert to linear index (faster than sub2ind)
      b1  = sum((binx + bt0(ones(n,1),:,one)-1).*fact1,2)+1; %linear index to c
      bt2 = bt0(:,:,two) + fact2;
      b2  = binx + bt2(ones(n,1),:);                     % linear index to X
      try
        if useSparse
          % Fast gridding using sparse
          c = c + sparse(b1,1,abs(prod(X(b2)-data,2)),Nc,1);
        else
          c = c + accumarray(b1,abs(prod(X(b2)-data,2)),[Nc,1]);
          %c = c + accumarray([b1,ones(n,1)],abs(prod(X(b2)-data,2)),[Nc,1]);
          %[len,bin,val] = bincount(b1,abs(prod(X(b2)-data,2)));
          %c(bin)        = c(bin)+val;
        end
      catch
         c = c + sparse(b1,1,abs(prod(X(b2)-data,2)),Nc,1);
      end
    end
  end
  c = reshape(c/w,csiz);
end
switch d % make sure c is stored in the same way as meshgrid
 case 2,  c = c.';
 case 3,  c = permute(c,[2 1 3]);
end
return

data =[ 0.75355792,  0.72779194,  0.94149169,  0.07841119,  2.32291887,...
    1.10419995,  0.77055114,  0.60288273,  1.36883635,  1.74754326,...
    1.09547561,  1.01671133,  0.73211143,  0.61891719,  0.75903487,...
    1.8919469 ,  0.72433808,  1.92973094,  0.44749838,  1.36508452]

data4 = [ 0.38103275,  0.35083136,  0.90024207,  1.88230239,  0.96815399,...
    0.57392873,  1.63367908,  1.20944125,  2.03887811,  0.81789145;...
    0.69302049,  1.40856592,  0.92156032,  2.14791432,  2.04373821,...
    0.69800708,  0.58428735,  1.59128776,  2.05771405,  0.87021964;...
    1.44080694,  0.39973751,  1.331243  ,  2.48895822,  1.18894158,...
    1.40526085,  1.01967897,  0.81196474,  1.37978932,  2.03334689;...
    0.870329  ,  1.25106862,  0.5346619 ,  0.47541236,  1.51930093,...
    0.58861519,  1.19780448,  0.81548296,  1.56859488,  1.60653533]

data3 = reshape([ 0.932896  ,  0.89522635,  0.80636346,  1.32283371,  0.27125435,...
    1.91666304,  2.30736635,  1.13662384,  1.73071287,  1.06061127,...
    0.99598512,  2.16396591,  1.23458213,  1.12406686,  1.16930431,...
    0.73700592,  1.21135139,  0.46671506,  1.3530304 ,  0.91419104,...
    0.62759088,  0.23988169,  2.04909823,  0.93766571,  1.19343762,...
    1.94954931,  0.84687514,  0.49284897,  1.05066204,  1.89088505,...
    0.840738  ,  1.02901457,  1.0758625 ,  1.76357967,  0.45792897,...
    1.54488066,  0.17644313,  1.6798871 ,  0.72583514,  2.22087245,...
    1.69496432,  0.81791905,  0.82534709,  0.71642389,  0.89294732,...
    1.66888649,  0.69036947,  0.99961448,  0.30657267,  0.98798713,...
    0.83298728,  1.83334948,  1.90144186,  1.25781913,  0.07122458,...
    2.42340852,  2.41342037,  0.87233305,  1.17537114,  1.69505988], 20,3)