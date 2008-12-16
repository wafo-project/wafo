function [zzgrid, xvec,options] = ffndgrid(x, f, varargin)
%FFNDGRID  Fast 'n' Furious N-D data gridding.
%
%  CALL:  [fgrid, xvec] = ffndgrid(x,f, options);
%
%  fgrid = Matrix of gridded data.
%  xvec  = cellarray of gridvectors 
%           xl1:dx1:xu1 or linspace(xl1,xu1,Nx1) depending on delta.
%  x     = [x1 x2,...xD] coordinate matrix for unevenly spaced data, f.
%          size NxD. 
%  f     = f(x), vector of function values length N.
%  options  = options defining the gridding ,i.e., structure with the
%             fields:
%      .delta   = [dx1, dx2 ,...,dxD] or [-Nx1, -Nx2,...,NxD], where 
%                dx1, dx2, ..., dxD and Nx1, Nx2, ..., NxD defines the
%                stepsize of grids and number of bins, respectively, in the 
%                D dimesional space.  If a single value of delta is given 
%                then this value is applied on all directions.(default -75) 
%      .xLimits = [xl1 xu1 xl2 ...xuN], defines the limits of the
%                D-dimensional x1-x2...xN-plane to grid.
%                (default [min(x1),max(x1),min(x2),.... max(xN)])
%      .fLimits = [fl fu] defines the limits of the
%                 function values to grid. (default [min(f),max(f)])
%      .method  = defines the method used on the function values falling  
%                 within each cell. Valid options are:
%            'sum'  : If each value of fgrid is the sum of all points                
%            'mean' : If each value of fgrid is the average of all points (default)
%            'min'  : If each value of fgrid is the minimum of all points
%            'max'  : If each value of fgrid is the maximum of all points
%      .Noutliers = OutlierThreshold. If there are grid points with Noutliers or
%                   less observations they will be considered as outliers 
%                   and removed outliers from the data set.
%                  When Noutliers is between 0 and 1 it is treated as the 
%                  percentage of the total number of data points.
%                  (default 0)
%  pad_value =  Value to pad empty gridpoints with.
%
% FFNDGRID grids unevenly spaced data in vector f into a matrix fgrid.
%
% NOTE: - The vector limits can be padded with NaNs if only
%         certain limits are desired, e g if xu1 and fu are wanted:
%
%            ffndgrid(x, f, 'xLimits',[nan,.5],'fLimits',[nan,45])
%
%       - If no output arguments are given, FFNDGRID will plot the gridded
%         function with the prescribed axes using PCOLOR.
%
% Examples:
% N = 500;D=2; sz = [N ,D ];
% x = randn(sz); z = ones(sz(1),1); 
% [nc, xv] = ffndgrid(x,z,'delta',-15,'method','sum');  % Histogram
% pcolor(xv{:},nc)     %
% [XV,YV]=meshgrid(xv{:});
% text(XV(:),YV(:),int2str(nc(:)))
% dx = [diff(xv{1}(1:2)) diff(xv{2}(1:2))];
% contourf(xv{:}, nc/(N*prod(dx))) % 2-D probability density plot.
% colorbar
% colormap jet
%
% See also  griddata

% Tested on: MatLab 4.2, 5.0, 5.1, 5.2, 5.3. 6.5, 7 and 7.01.
% History:
% revised pab 6Oct2005
% -added subfunction extract 
% -Options are now handled as a struct.
% -Included the methods min and max as suggestions by  R.BREGEON - Physicist, PHD - FRANCE
% revised pab 02.08.2001
% - made it general for D dimensions + changed name from ffgrid to ffndgrid
% -added nargchk + examples.
% -updated help header to wafo-style
% - moved dx and dy into delta =[dx,dy] 
% -removed call to bin
% modified by Per A. Brodtkorb 
% 05.10.98 secret option: aver
%          optionally do not take average of values for each point
% 12.06-98
% by
% 28.7.97, Oyvind.Breivik@gfi.uib.no.
%
% Oyvind Breivik
% Department of Geophysics
% University of Bergen
% NORWAY


% default values
Nx = 75;
defaultoptions = struct('delta',-Nx,'xLimits',[],'fLimits',[],'Noutliers',0,'method','mean','pad_value',0);

if (nargin==1) && strncmpi(x,'defaults',7)
  zzgrid = defaultoptions;
  return
end
error(nargchk(2,8,nargin))		

r = size(x,1);
if r==1,% Make sure x is a column vector.
  x = x(:);
end

[N,D]=size(x);
f = f(:);
if length(f)==1,
  f = f(ones(N,1),:) ;
elseif length(f)~=N,
  error('The length of f must equal size(x,1)!')
end

options = parseoptions(defaultoptions,varargin{:});

xLimits = zeros(1,2*D);
xLimits(1:2:2*D) = min(x,[],1);
xLimits(2:2:2*D) = max(x,[],1);

kx = find(isfinite(options.xLimits));
if any(kx<=2*D)
  kx(kx>2*D) = [];
  xLimits(kx) = options.xLimits(kx);
end
options.xLimits = xLimits;


fLimits = [min(f), max(f)];

kf = find(isfinite(options.fLimits));
if any(kf<=2)
  kf(kf>2) = [];
  fLimits(kf) = options.fLimits(kf);
end
options.fLimits = fLimits;


if isempty(options.pad_value)
  options.pad_value = defaultoptions.pad_value;
end

if isempty(options.method)
  options.method = defaultoptions.method;
elseif isnumeric(options.method)
  validMethods = {'sum','mean','min','max'};
  options.method = validMethods(option.method);
end

xL = options.xLimits(1:2:2*D);
xU = options.xLimits(2:2:2*D);

fL = options.fLimits(1);
fU = options.fLimits(2);

dx = repmat(defaultoptions.delta,1,D);

Nd = length(options.delta);
if Nd==0
  options.delta = defaultoptions.delta;
elseif Nd>D
  options.delta(D+1:end) = []; % remove superflous dimensions if any
end
  
if Nd ==1,
  options.delta = options.delta(ones(1,D));
end
ind = find(isfinite(options.delta) & options.delta~=0);
if any(ind),
  dx(ind)  = options.delta(ind);
end

ind = find(dx<0);
if any(ind), 
   if any(dx(ind)~=round(dx(ind))),
    error('Some of Nx1,...NxD in delta are not an integer!'),
  end
  dx(ind) = (xU(ind)-xL(ind))./(abs(dx(ind))-1);
end


% bin data in D-dimensional-space
binx = round((x - xL(ones(N,1),:))./dx(ones(N,1),:)) +1;

fgsiz = ones(1,max(D,2));
xvec  = cell(1,D);
for iy=1:D,
  xvec{iy} = xL(iy):dx(iy):xU(iy);
  fgsiz(iy) = length(xvec{iy});
end
if D>1
  in = all((binx >= 1) & (binx <= fgsiz(ones(N,1),:)),2) & (fL <= f) & (f <= fU);
else
  in = (binx >= 1) & (binx <= fgsiz(1)) & (fL <= f) & (f <= fU);
end
binx  = binx(in,:); 
f    = f(in);
N    = length(binx); % how many datapoints are left now?

Nf = prod(fgsiz);

if D>1,
  fact = [1 cumprod(fgsiz(1:D-1))];
  binx = sum((binx-1).*fact(ones(N,1),:),2)+1; % linear index to fgrid
end

if (strncmpi(options.method,'min',2) || strncmpi(options.method,'max',2))
  
    [sbinx, ind] = sort(binx); % general
    
    % Find indices to unique values
    i   = [ find(diff(sbinx)) ; length(sbinx) ];
    bin = sbinx(i);  % unique bins
    i   = [ 0 ; i ]+1; % indices marking the beginning and end of entries
    
    rho = length(bin)/Nf; % density of nonzero values in the grid
    
    if rho>0.4
      fgrid = zeros(Nf,1);
      ngrid      = zeros(Nf,1);
      ngrid(bin) = diff(i);
    else
      fgrid = spalloc(Nf,1,length(bin));
      ngrid = sparse(bin,1,diff(i),Nf,1);
    end
    
    % sort f and make a list
    sf = extract(f(ind), i(1:end-1),i(2:end),nan);
   
    
    if (strncmpi(options.method,'min',2) )
      fgrid(bin) = min(sf,[],2);
    else
      fgrid(bin) = max(sf,[],2);
      %else %mean or sum
      %sf(isnan(sf)) = 0;
      %fgrid(bin) = sum(sf,2)
      %fgrid(bin) = sum(sf,2)./ngrid(bin);
    end
   
else % sum or mean
    % Fast gridding
    fgrid = sparse(binx,1,f,Nf,1);% z-coordinate
    %fgrid = accumarray(binx,f,[Nf,1]);
end

removeOutliers     = (options.Noutliers~=0);
meanMethod         = strncmpi(options.method,'mean',2);
padEmptyGridPoints = (options.pad_value~=0);

if removeOutliers || meanMethod || padEmptyGridPoints
  ngrid = sparse(binx,1,ones(N,1),Nf, 1); % no. of obs per grid cell
  
  if removeOutliers,
    isPercentage =  (abs(options.Noutliers)<1);
    if (isPercentage),
      % no is given as percentage
      n0 = abs(options.Noutliers)*N;
    else
      n0 = options.Noutliers;
    end 
    fgrid(ngrid <= n0) = 0;
    ngrid(ngrid <= n0) = 0;
    N = full(sum(ngrid(:))); % how many datapoints are left now?
  end  
end


 ind = find(fgrid); % find nonzero values

if meanMethod
  fgrid(ind) = fgrid(ind)./ngrid(ind); % Make average of values for each point
end

if padEmptyGridPoints
  Nil = find(ngrid==0);
  if any(Nil)
    % Empty grid points are set to pad
    fgrid      = full(fgrid);
    fgrid(Nil) = options.pad_value;
  end
end

rho = length(ind)/Nf; % density of nonzero values in the grid
if ((rho>0.4) && issparse(fgrid))
  fgrid = full(fgrid); 
end
if r==1,
  fgrid = fgrid.';
else
  fgrid = reshape(fgrid,fgsiz);
  switch D % make sure fgrid is stored in the same way as meshgrid
    case 2,  fgrid=fgrid.';
    case 3,  fgrid=permute(fgrid,[2 1 3]);
  end
end




if (nargout > 0)
  zzgrid = fgrid;
elseif D==2,% no output, then plot
  colormap(flipud(hot)) %colormap jet
  if 1,
    %figure('Position', [100 100 size(fgrid)])    
    imagesc(xvec{:}, fgrid)
    set(gca,'YDir','normal')
  else
    pcolor(xvec{:}, fgrid)
    shading flat %interp
  end
  colorbar
  xlabel(inputname(1))
  ylabel(inputname(1))
  zstr=inputname(2);
  dum = full(size(fgrid'));
  if (~isempty(zstr)) % all this vital information ...
    str = sprintf('Color scale: %s, %d data points, grid: %dx%d, density: %4.2f', ...
	inputname(3), N, dum(1), dum(2), rho);
  else
    str = sprintf('%d data points, grid: %dx%d, density: %4.2f', ...
	N, dum(1), dum(2), rho);
 end
 title(str);
end




return

function  list = extract(v,no,nc,nfill,align)
%EXTRACT  Extract subentries of a VECTOR into a LIST (matrix)
%
%  CALL:  list = extract(V,NO,NC,NFILL,align)
% 
%   list  = extracted entries from V 
%   V     = character vector or numeric vector
%   NO,NC = vector of indices marking the beginning
%           and the end, respectively, of entries in V. 
%   Nfill = "filling" character or number in matrix LIST (default 0 for
%           numeric arrays and ' ' for character arrays)
%  align  = 'left'   if entries in list aligned to the left (default)
%           'right'  if entries in list aligned to the right. 
%  
% EXTRACT extracts the subentries of V into a LIST so that
%   LIST(i,1:N) = V(NO(i):NC(i)-1) where N = NC(i)-NO(i).
%   If NC(i) <= NO(i) then LIST(i,:) = NFILL 
%
% Example:
%  [names, p0, p1] = getnames('yy = foo(x1.^2,a_1*c,flag)');
%  names2  = extract('yy = foo(x1.^2,a_1*c,flag)',p0,p1+1);   % Same thing 
%  names3  = extract('yy = foo(x1.^2,a_1*c,flag)',p0,p1+1,[],'right');  
%
% See also  insert, getnames, cumsum


% History:
% revised pab 10oct2005
%  -updated help header.
% revised pab 02.11.2003
%  -added align  
% revised pab 18.07.2002
% -fixed a bug when e.g. NO = [1 6] and NC = [10 7] 
% -Removed the sorting
% revised pab 15.07.2002
% -fixed a bug when NC<=NO
% revised pab 04.10.2000
%  -updated documentation, changed isstr -> ischar
%  Kirill K. Pankratov,  kirill@plume.mit.edu
%  8-22-94

error(nargchk(3,5,nargin))
nfilldflt = [0 32]; % Defaults for filling (0 for numbers
                    % or blank for strings)
		    
 % Handle input ...........................................
if nargin<4||isempty(nfill), 
  nfill = nfilldflt(1+ischar(v));
end
if nargin<5||isempty(align)
  align = 'left';
end

v  = v(:);
nc = nc(:);
no = no(:);
Lv = length(v);

if any(Lv<no)
  warning('WAFO:ffndgrid','An index in NO can not be larger than L = %d',Lv-1)
end

ind   = find(nc<=no);
no    = min(max(no,1),Lv);      % Make sure 1<=no<Lv+1
nc    = min(max(nc,no+1),Lv+1); % Make sure no<nc<=Lv+1 pab 14.12.2000


d      = nc-no;       % length of every entry to extract.
lline  = max(d);      % maximum length of entry to extract
nlines = length(d);   % number of entrys to extract
  
  
num = cumsum(d);
L   = num(nlines); % number of letters to extract

% Make indices to the output list
nout                 = ones(L,1);
switch lower(align(1)),
 case 'r', %right
  nout(1+num(1:end)-d(1:end)) = lline-d(1:end)+1;
 otherwise % left adjust
  nout(1+num(1:end-1)) = lline-d(1:end-1)+1;
end
nout                 = cumsum(nout);


% Make index vector nin = [no(1):nc(1)-1, no(2):nc(2)-1,....]
nin                 = ones(L,1);
nin(1)              = no(1);
nin(1+num(1:end-1)) = no(2:nlines)-nc(1:nlines-1)+1;
nin                 = cumsum(nin);

list       = nfill(ones(lline,nlines)); % new call
list(nout) = v(nin);

list = list.';

if any(ind)
  list(ind,:) = nfill;
end
if ischar(v)
  list = char(list);
end


return




