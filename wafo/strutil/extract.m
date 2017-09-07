function  list = extract(v,no,nc,nfill,align)
%EXTRACT  Extract subentries of a VECTOR into a LIST (matrix or cell)
%
%  CALL:  list = extract(V,NO,NC,NFILL,align)
%         list = extract(V,NO,NC,Ltype)
% 
%   list  = matrix or cell array of extracted entries from V 
%   V     = character vector or numeric vector
%   NO,NC = vector of indices marking the beginning
%           and the end, respectively, of entries in V. 
%  Ltype  = Defines type of list. Options are
%           'cell' or 'matrix' (default).
%   Nfill = "filling" character or number in matrix LIST (default 0 if 
%           V is a numeric array and ' ' for a character arrays)
%  align  = 'left'   if entries in matrix list aligned to the left (default)
%           'right'  if entries in matrix list aligned to the right. 
%  
% EXTRACT extracts the subentries of V into a LIST so that
% For TYPE = 'matrix':
%   LIST(i,1:N) = V(NO(i):NC(i)-1) where N = NC(i)-NO(i).
%   If NC(i) <= NO(i) then LIST(i,:) = NFILL 
% For TYPE = 'cell'
%   LIST{i} = V(NO(i):NC(i)-1) where N = NC(i)-NO(i).
%   If NC(i) <= NO(i) then LIST{i} = [] 
%
% Example:
%  s = 'yy = foo(x1.^2,a_1*c,flag)';
%  [names, p0, p1] = getnames(s);
%  assert(names, strvcat('yy', 'foo', 'x1', 'a_1', 'c', 'flag'))
%  rnames = strvcat('  yy', ' foo', '  x1', ' a_1', '   c', 'flag');
%  cnames = {'yy', 'foo', 'x1', 'a_1', 'c', 'flag'};
%  assert(extract(s, p0, p1+1), names) % Same thing 
%  assert(extract(s, p0, p1+1, [], 'right'), rnames)  
%  assert(extract(s, p0, p1+1, 'cell'), cnames)  
%
% See also  insert, getnames, cumsum


% History:
% revised pab jan 2007
% revised pab aug 2006
% - added celloutput
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

%error(nargchk(3,5,nargin));
narginchk(3,5)
nfilldflt = [0 32]; % Defaults for filling (0 for numbers
                    % or blank for strings)
cellOutPut = (1==0);                    

 % Handle input ...........................................
if nargin<4||isempty(nfill) || (length(nfill)>1 && ~strncmpi(nfill,'cell',2)), 
  if ischar(v)
    nfill = char(nfilldflt(2));
  else
    nfill = nfilldflt(1);
  end
elseif strncmpi(nfill,'cell',2)
  cellOutPut = (1==1);
end

sz = size(no);

v  = v(:);
nc = nc(:);
no = no(:);
Lv = length(v);

if any(Lv<no)
  warning('STRUTIL:EXTRACT','An index in NO can not be larger than L = %d',Lv-1)
end
ind0  = find(no<nc);
ind   = find(nc<=no);
no    = min(max(no,1),Lv);      % Make sure 1<=no<Lv
nc    = min(max(nc,no+1),Lv+1); % Make sure no<nc<=Lv+1 pab 14.12.2000


if cellOutPut
  list = cell(sz);
  for ix = ind0(:).'
    list{ix} =  v(no(ix):nc(ix)-1).';
  end
else
  % array output
  if nargin<5||isempty(align)
    align = 'left';
  end

 

  d      = nc-no;       % length of every entry to extract.
  lline  = max(d);      % maximum length of entry to extract
  nlines = length(d);   % number of entrys to extract
  
  if nlines<=0
    list = '';
  else
  
    num = cumsum(d);
    L   = num(nlines); % number of letters to extract
    
    % Make indices to the output list
    nout = ones(L,1);
    switch lower(align(1)),
      case 'r', %right
        nout(1+num(1:end)-d(1:end)) = lline-d(1:end)+1;
      otherwise % left adjust
        nout(1+num(1:end-1)) = lline-d(1:end-1)+1;
    end
    nout = cumsum(nout);


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
    %if ischar(v)
    %  list = char(list);
    %end
  end
end

return
