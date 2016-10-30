function R=createcov(nr,vari,ctype)
% CREATECOV Covariance struct constructor
%
% CALL:  R=createcov(nr,variables,ctype)
% 
%         nr = number of derivatives (default 0   )
%  variables = 'xyt'                 (default 't' )
%       ctype = 'enc', 'none'   (default none)
%
% Example: 
%  R = createcov(1,'tx'); 
%  R.date = '';
%  trueR = struct('R',[], 'x', [], 't', [], 'h', inf, 'tr', [], 'phi', 0, ...
%       'type', 'none', 'norm', [], 'Rx', [], 'Rt', [], 'note', [], 'date', '');
%  assert(fieldnames(R), fieldnames(trueR))
%  assert(struct2cell(R), struct2cell(trueR))
%
% See also  createspec, datastructures

%Tested on: Matlab 7.0
%History:
% revised pab feb2007
% -removed all eval statements
% by pab 12.08.99

if nargin<1||isempty(nr)
  nr = 0;
end
if nargin<2||isempty(vari)
  vari = 't';
  Nv   = 1;
else
  vari = sort(lower(vari));
  Nv   = length(vari);
  if ( (vari(1)=='t') && (Nv>1)),
    vari = [vari(2:Nv), 't'];
  end 
end

if nargin<3||isempty(ctype)
  ctype='none';
end

R=struct('R',[]);
for ix=1:Nv
  R.(vari(ix)) = [];
end
if Nv==1,
  R.stdev=[];
end
R.h=inf;
R.tr=[];
R.phi=0;
if strcmp(ctype(1:3),'enc')
   R.v=0;
end
R.type=ctype;
R.norm=[];

% initialize struct
switch Nv
 case 1, 
  fieldname = ['R' vari(ones(1,nr))]; 
  for ix=1:nr
    R.(fieldname(1:ix+1)) = [];
  end   
case 2,
  for ix=1:nr
    for iy=0:nr-ix
      fn=['R' vari(ones(1,ix),1)' vari(ones(1,iy),2)'];
      %eval(['R.R' tmp '=[];'  ])
      R.(fn) = [];
    end
  end   
  for ix=1:nr
    fn = ['R' vari(ones(1,ix),2)'];
    R.(fn)=[];
  end   
case 3,
 for ix=1:nr
    for iy=0:nr-ix
      for iz=0:nr-ix-iy
        fn=['R' vari(ones(1,ix),1)' vari(ones(1,iy),2)' vari(ones(1,iz),3)'];
        R.(fn) =[];
      end
    end
  end   
  for ix=1:nr
    for iy=0:nr-ix
      fn=['R' vari(ones(1,ix),2)' vari(ones(1,iy),3)'];
      R.(fn)=[];
    end
  end   
  for ix=1:nr
    fn = ['R' vari(ones(1,ix),3)'];
    R.(fn) =[];
  end   

otherwise, error('unknown # variables')
end

R.note=[];
R.date=datestr(now);
