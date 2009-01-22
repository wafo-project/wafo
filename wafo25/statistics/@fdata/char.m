function s = char(self)
%CHAR Convert a class object in a one line description string
%
%  CALL str = char(thisClass) 
%  
%  CHAR is a class convertor from a class object to a string, used
%  in online display.
%
% Example
% d = pdffit;
% char(d)
%  
%  See also display

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