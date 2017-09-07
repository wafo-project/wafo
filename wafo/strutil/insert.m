function  vout = insert(vin,list,no,nc,nfill)
%INSERT  Inserts entries from LIST into a VECTOR
%
%  CALL: vout = insert(vin,list,no,nc,nfill);
%
%  vout, vin = string or numeric vectors.
%  list      = cells of strings, character array or double matrix
%              to be inserted into vin.
%  no, nc    = vector of indices marking the beginnings and ends,
%              respectively, of where to insert the entries from the list.  
%              (default nc = no)
%  nfill     = filling number (default 32 if vin is a character
%              array, otherwise 0)
%
% INSERT inserts the rows of matrix LIST into vector VIN.
% If NO(i)<NC(i) then the row LIST(i,:) is inserted, replacing
% the characters from NO(i) to NC(i)-1 of VIN.
% If NO(i)==NC(i) then the row list(i,:) is inserted between
% the characters NO(i)-1 and NO(i) of VIN.
% NFILL specifies the "filling" number in matrix LIST:
% for each row of LIST only part up to the last number
% different from NFILL is inserted into VIN.
%
% Examples:
% s1 = 'How much wood would a Woodchuck chuck?';
% no = [5 15]; s0 = strvcat('test','j');
% assert(insert(s1,s0,no), 'How testmuch wood jwould a Woodchuck chuck?')
% assert(insert(s1,s0(1,:),no), 'How testmuch wood testwould a Woodchuck chuck?')  
% assert(insert(s1,s0,no,no+3), 'How testh wood jld a Woodchuck chuck?')  
% ix = findstr('wood',s1);
% assert(insert(s1,s0,ix), 'How much testwood would a Woodchuck chuck?')
% assert(insert(s1,s0,ix,ix+4), 'How much test would a Woodchuck chuck?')
%
% See also: strrepl, extract, getnames, findname

% Tested on Matlab 5.3
% History:
% revised pab 30.10.2003  
% revised pab 25.01.2001
% - updated help header
% - if no is empty do nothing and return vout=vin;
% revised pab 14.12.2000
% - changed isstr to ischar
% - updated header info + added example
% - fixed 2 bugs: 1) made sure 1<= no <= nc <=l  
%                 2) added the line a(nc) = a(nc)+1;
%  Kirill K. Pankratov,  kirill@plume.mit.edu
%  8/22/94


nfilldflt = [0 32]; % Defaults for filling (0 for numbers
                    % or blank for strings
% Handle input .......................................
%error(nargchk(2,5,nargin))
narginchk(2,5)
if nargin<3||isempty(no),
  vout = vin;
  return
end
if iscell(list)
  list = strvcat(list{:});
end
if nargin<5||isempty(nfill), 
  nfill = nfilldflt(1+ischar(list)); 
end
 % If inserts ends are not specified, make them equal
 % inserts beginnings
if nargin<4||isempty(nc), 
  nc = no; 
end

vin   = [vin(:).' nfill]; % added ' ' to vin pab 14.12.2000
no    = no(:).'; 
nc    = nc(:).';
lmin  = min(length(no),length(nc));

if size(list,1)==1,
  list = list(ones(lmin,1),:);
end

L     = length(vin);
szl   = size(list);

lmin  = min(lmin,szl(1));
list  = list(1:lmin,:);
%no    = no(1:lmin); % old call
%nc    = nc(1:lmin); % old call
no    = min(max(no(1:lmin),1),L);  % Make sure 1<=no<=L
nc    = min(max(nc(1:lmin),no),L); % Make sure no<=nc<=L pab 14.12.2000


a     = ones(1,L);
a(no) = a(no)+L;
a(nc) = a(nc)-L;

a     = cumsum(a);
a(nc) = a(nc)+1;           % Make sure insertion without replace works
                           % correctly. pab 14.12.2000 

aa    = flipud(cumprod(double(rot90(list==nfill))))*L;
aa    = aa+no(ones(szl(2),1),:);
[vout,a] = sort([a aa(:)']);
l0     = find(vout<=L, 1, 'last' );
list  = list.';
vout  = [vin list(:)'];
vout  = vout(a(1:l0));

if all(nc<L)
  vout(end) = []; % remove the added space pab 14.12.2000
end
return





