function addtag(tag,h);
% ADDTAG	adds string to an object's tag
%
% addtag(h,tag);
%
% tag	= string to add to object's tag
% h	= handle to object	(default = current object)
%
% See also SET

error(nargchk(1,2,nargin));
if nargin<1 || isempty(h),	h=gco; end

otag=get(h,'tag');

if any(findstr(otag,tag))	% already tagged with this tag
  return;
else				% tag is new
  if isempty(otag), sep=''; else sep='-'; end	% tag-separator
  ntag=strcat(tag,sep,otag);			% build whole tag
  set(h,'tag',ntag);				% set new tag
end
