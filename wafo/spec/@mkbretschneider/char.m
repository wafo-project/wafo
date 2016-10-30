function s = char(self)
%MKBRETSCHNEIDER/CHAR Convert object in a one line description string
%
%  CALL str = char(Hs) 
%  
%  CHAR is a class convertor from object to a string, used
%  in online display.
% 
% Example
% d = mkjonswap;
% char(d)
%
%  See also mkbretschneider/display

% History
% By pab 2007


  s0 = struct(self);

  names = fieldnames(s0);
  n = length(names);
  tmp = cell(n,1);
  for i = [1:n-1]
    name = names{i};
    val = self.(name);
    tmp{i, 1} = strcat(name, '=', num2str(val), ', ');
  end
  name = names{n};
  val = self.(name);
  tmp{n, 1} = strcat(name, '=', num2str(val));
  s = strcat('mkbretschneider(', tmp{:}, ')');
end