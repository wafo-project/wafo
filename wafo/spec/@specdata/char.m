function s = char(self)
%SPECDATA/CHAR Convert a SPECDATA object in a one line description string
%
%  CALL str = char(Hs) 
%  
%  CHAR is a class convertor from SPECDATA to a string, used
%  in online display.
% 
% Example
% d = specdata;
% char(d)
%
%  See also specdata/display

% History
% By pab 2007

if numel(self)<=1
  s = disp(struct(self));
else
    fn = fieldnames(self);
    fnTxt = sprintf('    %s\n',fn{:});
    szTxt = sprintf('%dx',size(self));
    szTxt(end) = [];
    
    s = sprintf('%s %s array with fields:\n%s',szTxt,class(self),fnTxt);
 end