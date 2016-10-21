function Slope=wav2slope(L,lev,dense,p)
%WAV2SLOPE Extracts slopes at up- and downcrossings after smoothing
%
%CALL: Slope=wav2slope(L,levels,dense,p)
%
%   Slope   = structure with observed up- and down-crossing slopes
%       
%   L       = data with fields L.t/L.u and L.Z
%   levels  = vector with absolute levels; (required field, no default)
%   dense   = interpolation rate,
%                 []: no interpolation or smoothing
%                  positive integer: data are interpolated at rate
%                  dense at equidistant points (default: dense=20)
%   p       = [0...1] is a smoothing parameter (default=1)
%                   0 -> LS-straight line
%                   1 -> cubic spline interpolation

%   Created by GL 2014 for general use, takes both time and space waves
%   Revised by GL 2015-02-16 to treat up- and downcrossings in a symmetric way

if isfield(L,'u'),
    Lx=L.u;
elseif isfield(L,'t'),
    Lx=L.t;
end
if min(diff(Lx))<=0,
    disp('Series argument must be in ascending order')
    return
end

Lz=L.Z;
[mL,nL]=size(Lx);
if mL>nL, Lx=Lx'; Lz=Lz'; end

if nargin<2 || isempty(lev),
    disp('No levels submitted to wav2slope')
    return
else
    levels=lev;
end

if nargin<4 || isempty(p),
    p=1;
end
if (nargin>2) && (~isempty(dense)),
    dens=dense;
else
    dens=10;
end

x=linspace(Lx(1),Lx(end),dens*(length(Lx)-1)+1);
Z=cssmooth(Lx,Lz,p,x);
dx=x(2)-x(1);

%gradU=gradient(Z,x);
%gradD=gradient(fliplr(Z),fliplr(x));

nlevel=length(levels);
Slope=struct('up',[],'down',[],'levels',[]);
Slope.levels=levels;
for j=1:nlevel,
    level=levels(j);
    ind=findcross(Z,level);
    indup=ind(Z(ind)<Z(ind+ones(size(ind))));
    inddo=ind(Z(ind)>Z(ind+ones(size(ind))));
    gradU=(Z(indup+ones(size(indup)))-Z(indup))/dx;
    gradD=(Z(inddo+ones(size(inddo)))-Z(inddo))/dx;
%    [cross_up,Ncu]=dat2crossind(Z,level,'u',true);
%    [cross_down,Ncd]=dat2crossind(fliplr(Z),level,'u',true);
    Slope.up{j}=gradU';
    Slope.down{j}=gradD';
end
