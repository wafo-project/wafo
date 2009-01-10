function fi = evalpdf(pdf,varargin) 
% EVALPDF evaluates a PDF struct by interpolation
%
% CALL: f = evalpdf(pdf,x1,x2,...,xd,method)
% 
%  f            = evaluated pdf at x1,x2,...,xd
%  pdf          =  PDF structure 
%  x1,x2,...,xd = vectors/matrices of points where the pdf is evaluated
%  method       = 'nearest' - nearest neighbor interpolation
%                 'linear'  - linear interpolation  (default)
%                 'spline'  - cubic spline interpolation
%                 'cubic'   - cubic interpolation
%
%  For faster interpolation when pdf.x{1},... pdf.x{d} are equally spaced and
%  monotonic, use the methods '*linear', '*cubic', or '*nearest'.
%  Out of range values are returned as NaN's.
%
% See also  datastructures
if isstruct(pdf)
 d=ndims(pdf.f);

 if d==2
  fsiz=size(pdf.f);
  if min(fsiz)==1
    d=1;
  end
 end
else
 error('pdf must be a structure') 
end
Xi=cell(d,1);
[Xi{1:d}]=deal(varargin{1:d});
nv=length(varargin);
if d<nv,
  method=varargin{nv};
else
  method='linear';
end
switch d
  case 1, 
    fi = interp1(pdf.x{1},pdf.f,Xi{1},method);
  case 2,  
    X=cell(d,1);
    [X{:}] = meshgrid(pdf.x{:});
    fi =interp2(X{:},pdf.f,Xi{:},method);
  case 3,
    X=cell(d,1);
    [X{:}] = meshgrid(pdf.x{:});
    fi =interp3(X{:},pdf.f,Xi{:},method);
  otherwise ,  
    disp('Dimension of data large, this will take a while.')
    [X{:}]    = ndgrid(pdf.x{:});
    fi =interpn(X{:},pdf.f,Xi{:},method);
end
