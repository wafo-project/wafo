function display(self,name)
%SPECDATA/DSIPLAY Display a SPECDATA object in Matlab window
% 
% Example
% d = specdata;
% display(d)
% 
% See also char

% Tested on matlab 7.3
% History:
% by pab 2007

if nargin>1 % secret call
  if isempty(name)
    name = inputname(1);
  end
  disp(' ');
  disp([name,' = ']);
  disp(' ');
end
disp(char(self))

