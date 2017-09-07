function y = convlv(x,h,dim,flag)
%CONVLV Convolves real data set with a response function. 
%
%  CALL:  y = convlv(x,h,dim,flag);
%
%  y,x  = filtered data and raw data, respectively.
%  h    = vector with filter coefficients.
%  dim  = dimension to filter.
%  flag = 'periodic'    : if X is periodic
%         'nonperiodic' : otherwise (default)
%  
% CONVLV convolves X with H, i.e., :
%
%              nr
%       Y(i) = sum Cn(n)*X(i+n)
%              n=-nl
%  where
%    H(1) = Cn(0), H(2)=Cn(-1),...H(Nl+1)=Cn(-Nl), and H(Np) = Cn(1),
%    H(Np-1) = Cn(2),...H(Np-Nr) = Cn(Nr),
%
% Note: The filter, H, must be stored in wrap-around order, i.e., the
%       first half of the array H contains the impulse response
%       function at positive times, while the second half of the array
%       contains the impulse response function at negative times,
%       counting down from the highest element. 
%
% Example:
% dx = 2*pi/100; 
% x = linspace(0,2*pi-dx,100)';
% y = cos(x);
% c = savgol(2,2,4,1);          % differentiation filter
% yf = convlv(y,c)/dx;          % derivative
% yf2 = convlv(y,c,[],'p')/dx;  
% semilogy(x,abs(sin(x)+yf),x,abs(sin(x)+yf2),'g')
%
% See also  savgol  

%History  
% Revised pab 2003
% - fixed a bug when nr=nl=0 and when nh~=1  
% By Per A. Brodtkorb 2001
  
%error(nargchk(2,4,nargin))
narginchk(2,4)
if nargin<3,dim =[];end
if nargin<4||isempty(flag),flag ='nonperiodic';end

perm = []; nshifts = 0;
if ~isempty(dim)
  perm = [dim:max(ndims(x),dim) 1:dim-1];
  x = permute(x,perm);
else                      
  [x,nshifts] = shiftdim(x);
end
xsiz = size(x);

if isempty(perm) 
  h    = shiftdim(h);
  hsiz = size(h);
else
  hsiz = size(h);
  if prod(hsiz)==prod(xsiz)
    h = permute(h,perm);
    hsiz = size(h);
  end
end

m  = xsiz(1);
mh = hsiz(1);
nh = prod(hsiz(2:end));
p  = ceil(mh/2);
nhr = prod(xsiz(2:end)) - nh + 1;

if (p==mh)
  y  = x(:,:).*h(ones(m,1),:);
else
 % nl = max(find(h(1:p,1)));      % Number of non-zero elements to the left
 % nr = max(find(h(mh:-1:p+1,1)));% and to the right.
  nl = find(h(1:p,1),1,'last');      % Number of non-zero elements to the left
  nr = find(h(mh:-1:p+1,1),1,'last');% and to the right.
  nz = nr+nl;                    % Total number of non-zero elements
  
  
  if strncmpi(flag,'p',1)
    % periodic endconditions
    nfft = max(m,nz);
    hw   = fft([h(1:p,:); zeros(nfft-mh,nh); h(p+1:end,:)]);
    y    = real(ifft(fft(x(:,:),nfft).*repmat(hw,[1,nhr]))); % convolution
    y    = reshape(y(1:m,:),[ones(1,nshifts),xsiz]);
  else
    % nonperiodic end conditions
    % Extrapolate x linearly outside outside the range to lessen the end
    % effects
    y  = linear((1:m)',x(:,:),(-nl+1:m+nr)');
    %  y  = spline((1:m)',x(:,:),(-nl+1:m+nr)');
    if 0,
      %This is wrong
      y=(conv(y,-h(1:p,:))+flipud(conv(flipud(y),-[0; h(end:-1:p+1,:)]))); 
    else
      nfft = 2^nextpow2(m+nz);
      hw   = fft([h(1:p,:); zeros(nfft-mh,nh); h(p+1:end,:)]);
      y  = real(ifft(fft(y(:,:),nfft).*repmat(hw,[1,nhr]))); % convolution
    end
    y = reshape(y(nl+1:nl+m,:),[ones(1,nshifts),xsiz]);
    
  end
end
if ~isempty(perm), y = ipermute(y,perm); end  



function yi = linear(x,y,xi)
% LINEAR Interpolation and extrapolation
%
% CALL: yi = linear(x,y,xi);
%
%  x, xi   = column vectors
%  y, yi   = column vectors or matrices
%

siz = size(xi);
n = length(x);
ni = length(xi);
if ni~=1
    [xxi,k] = sort(xi);
    [dum,j] = sort([x;xxi]);
    r(j) = 1:n+ni;
    r    = r(n+1:end)-(1:ni);
    r(k) = r;
    r    = max(1,min(r,n-1));
    u    = (xi-x(r))./(x(r+1)-x(r));
    yi   = y(r,:)+(y(r+1,:)-y(r,:)).*u(:,ones(1,size(y,2)));
else
    % Special scalar xi case
    r = find(x <= xi,1,'last');
    if isempty(r),  % xi<x(1)
      r=1; 
    elseif r==n,    % x(n)<=xi
      r = n-1;
    end
    if (r>0) && (r<n),
        u = (xi-x(r))./(x(r+1)-x(r));
        yi=y(r,:)+(y(r+1,:)-y(r,:)).*u(:,ones(1,size(y,2)));
    else
        yi = NaN;
    end
end
if (min(size(yi))==1) && (prod(siz)>1), yi = reshape(yi,siz); end

return