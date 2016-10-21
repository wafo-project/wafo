function f=createpdf(varargin)
% CREATEPDF PDF struct constructor
%
% CALL:  f=createpdf(ndim)
%        f=createpdf(list)
%
%  creates a struct with the following fields
%      
%      f.f      = pdf (ndim-dimensional matrix)
%      f.x      = cellarray of vectors of lags in ndim dimensions (ndim cells)
%      f.labx   = cellarray of strings of label strings (ndim cells)           
%      f.title  = title string
%      f.note   = note string
%      f.date   = creation date and time
%
%      ndim     = # dimensions (default 1)
%
% Examples: 
%  f = createpdf(2) %gives the structure
%  %      f: []
%  %      x: {2x1 cell}
%  %   labx: {2x1 cell}
%  %  title: []
%  %   note: []
%  %   date: '16-Oct-1999 16:56:54'
%  x = linspace(0,10);
%  f = createpdf('f',sin(x),'x',{x})
%  pdfplot(f)
%
% See also  datastructures, pdfplot, createcov

%Tested on: Matlab 5.3
%History:
% revised es 25.10.1999 help and cosmetics 
% revised pab 16.10.1999
% changed .xi and .labxi to .x and .labx cellarrays for 
%      faster computation and easier access
% by pab 16.09.99

ndim = 1;
names = {'f','x','labx','title','note','date'};
n    = length(names);
c    = cell(1,n);
c{n} = datestr(now);
f    = cell2struct(c,names,2);

if nargin<=1
  if nargin==1&&~isempty(varargin{1})
    ndim=varargin{1};
  end

  %f=struct('f',[]);
  f.x=cell(ndim,1);
  f.labx=cell(ndim,1);
  %f.title=[];
  %f.note=[];
  %f.date=datestr(now);
elseif nargin>1
  f = parseoptions(f,varargin{:});
end
