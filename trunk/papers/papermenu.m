function papermenu(kind)
% PAPERMENU displays a user interface to the paper scripts of WAFO
%
% CALL:  papermenu(kind)
%
%    kind = 0, displays a point and click menu (default)
%           1, displays a text driven menu
%
% See also  wafomenu

%history
% by pab 20.05.2000

if nargin<1|isempty(kind)
  kind=0;
end

header = 'WAFO PAPER scripts';
labels =str2mat(...
    ' 1)  RECDEMO: A statistical procedure for reconstruction of 1D signals.', ...
    ' 2)  WAFODEMO: WAFO - a Matlab toolbox for analysis of random waves and loads.') ;
filenames = str2mat(...
    'recdemo', ...
    'wafodemo');

Noptions=size(filenames,1);
if kind==1, % text driven menu
  r=1;
  while ~isempty(r) & ~strcmpi(r,'q')
    clc
    disp(header)
    disp('  ')
    disp(labels)
    disp('  ')
    r=input(['Enter your choice 1,2,...,' num2str(Noptions-1),' or ' ...
	  num2str(Noptions) '.  ']);
    if ~isempty(r)
      clc
      r=round(r);
      if ((1<=r) & (r<=Noptions))
	eval(filenames(r,:))
	%pause
      end  
    end % if r
  end % while  
else % point and click menu
  callbacks = [filenames]; 
  choices('PAPER',header,labels,callbacks);
end % kind

return

