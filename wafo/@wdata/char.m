function s = char(self)
%WDATA/CHAR Convert a WDATA object in a one line description string
%
%  CALL str = char(Hs) 
%  
%  CHAR is a class convertor from WDATA to a string, used
%  in online display.
%
% Example
% d = wdata;
% char(d)
%  
%  See also wdata/display

% History
% By pab 2007

if numel(self)<=1
  s = evalc('disp(struct(self))');
else
    fn = fieldnames(self);
    fnTxt = sprintf('    %s\n',fn{:});
    szTxt = sprintf('%dx',size(self));
    szTxt(end) = [];
    
    s = sprintf('%s %s array with fields:\n%s',szTxt,class(self),fnTxt);
 end