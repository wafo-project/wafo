function display(self,name)
%DISPLAY Display a class object in Matlab window
% 
% Example
% d = wdata;
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

