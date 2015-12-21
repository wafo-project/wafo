function [I00,I11] = alias(I0,order)
%ALIAS Alias structure of a fractional design.
%
% CALL:  str = alias(I0,n);
%
% str = string containing the alias structure.
% I0  = complete defining relation.
% n   = maximum order of alias structure.
% 
% Example
%   I = sudg(6,2);  % Design generator
%   I0 = cdr(I);    % Complete defining relation
%   alias(I0)       % The complete alias structure
%   alias(I0,3)     % Alias structure neglecting interactions larger than  3
%
% See also  cdr

%
%     This program is free software; you can redistribute it and/or modify
%     it under the terms of the GNU Lesser General Public License as published by
%     the Free Software Foundation; either version 3 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU Lesser General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.


%Known Bugs:   1) Number of variables must be less than 51
%              2) n<inf produse
error(nargchk(1,2,nargin))
if nargin<2||isempty(order),order = inf; else order = abs(order);end
if isempty(I0),
  str = 'No alias structure';
  disp(str)
  return
end


if isnumeric(I0)
  k = max(I0(:));         % Number of variables.
  I0 = cnr2cl(I0);
else
  k = max(cl2cnr(I0(:))); % Number of variables.
end
p0 = size(I0,1);

p = log2(p0+1);   % Number of generators.
if p~=round(p),   % May be it is a design generator
  warning('WAFO:ALIAS','This is not a Complete Defining Relation')
  I0 = cdr(I0);
  p0 = size(I0,1);
  p = log2(p0+1);   % Number of generators.
end

n = 2^(k-p);      % Number  of aliases

% Make a list of all possible main effects and interaction effects
%------------------------------------------------------------------
id = zeros(2^k-1,k);
iz = 0;
for ix = 1:k,
  iz       = iz+1;
  id(iz,1) = ix;
  iz0      = iz;
  for iy = 1:iz0-1,
    iz = iz+1;
    id(iz,:)   = id(iy,:);
    ind        = find(id(iy,:)==0, 1 );
    id(iz,ind) = ix;    
  end
end
% Make sure whitespace is at the end.
id = fliplr(sortrows(fliplr(cnr2cl(id))));

if 0, % This is not needed
  % Remove the defining relation effects from the effects list
  for ix=1:p0
    k1 = strmatch(I0(ix,:),id,'exact');
    if any(k1)
      id(k1,:)=[];
    end
  end
end

id = cellstr(id);

I0 = cellstr(I0);
I11 = cell(n,p0+1);


% Find all aliases
%------------------------------------------
[I11{end,1:p0}] = deal(I0{:});
I11{end,p0+1}= ' ';
wl = zeros(1,p0);    % word length
for ix =1:n-1
  I11(ix,1) = {id{ix}};
  for iy = 1:p0
    %disp([ix,iy])
    tmp          = setxor(id{ix},I0{iy});
    wl(iy) = length(tmp)-any(tmp=='-');  % Save the word length for later sorting.
    I11(ix,iy+1) = {tmp};
    k1 = strmatch(tmp,id,'exact');
    if any(k1),            % Remove found aliases from the effects list.
      ind     = ones(size(id));
      ind(k1) = 0;
      id      = id(logical(ind));
    end
  end
  % Sort by word length
  [wl,ind]=sort(wl);
  I11(ix,2:end) = I11(ix,ind+1);
  k1 = find(wl>order);
  if any(k1), % remove interactions of higher order than order.
    [I11{ix,k1+1}]= deal(' ');
  end
  if length(id)<=ix,
    warning('WAFO:ALIAS','Something is wrong')
    break,
  end
end

tmp = cell(2,p0+1);
[tmp{2,1:p0}]=deal( ' + ');
tmp{2,p0+1} = ' ';
I00= [];
for ix = 1:n
  [tmp{1,:}] = deal(I11{ix,:});
  % old call
  %tmp1 = char(tmp)';
 
  % New call
  tmp1 = cat(2,tmp{:});
  % Remove any multiple whitespaces.
  tmp1 = rstrrep(tmp1(:).','  ', ' '); % rstrrep from string utility toolbox
  k1 = find(tmp1 == '-');
  if any(k1),    
    tmp1(k1-2)='-';  % Change + to -
    tmp1(k1)='';     % Remove old - signs
  end
  I00 = strvcat(I00,tmp1);
end



  
  

