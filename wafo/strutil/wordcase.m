function text=wordcase(text)
%WORDCASE Captializes first letters of words 
%         in a (cell)string. Trailing spaces are also removed.
%
% CALL:  text = wordcase(text);
% 
%Example
% s1 = 'How much wood would a Woodchuck chuck?';
% assert(wordcase(s1), 'How Much Wood Would A Woodchuck Chuck?')
%
% See also: upper, deblank

%Time-stamp:<Last updated on 00/11/23 at 14:40:55 by even@gfi.uib.no>
%File:<d:/home/matlab/wordcase.m>

cll=0;
%error(nargchk(1,1,nargin)); % INPUT-TEST
narginchk(1,1)
if isnumeric(text), error('String input only'); end

text = deblank(text);
if isempty(text), return; end

if iscell(text)
  cll=1;
  text=char(text);
end
text=lower(text);
[M]=size(text,1);
for i=1:M
  linje=deblank(text(i,:));
  spc=[findstr(' ',linje) findstr('-',linje)];
  caps=[0 spc]+1;
  text(i,caps)=upper(linje(caps));
end

if cll, text=cellstr(text); end