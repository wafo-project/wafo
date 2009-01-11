function URLout = wafohelp(ButtonName)
%WAFOHELP  Launch HTML help for the WAFO Toolbox.
%
%  CALL:  URL = wafohelp(BN)
% 
%   URL = string with the URL address of WAFO.
%   BN  = 'yes' If display URL site
%         'no'  If no display of URL site
%
% WAFOHELP displays and/or returns the HTML help
%   site for the WAFO Toolbox.  
%
% Examples
%  url = wafohelp('n'); web(url)
%  wafohelp 
%  wafohelp('yes')
%
% See also  web
 
% Tested on: Matlab 5.2
% History:
% by pab 19.11.2004

URL = fullfile(waforoot,'docs','wafodoc', 'index.html');

if nargin<1 || isempty(ButtonName),  
  if nargout>0, 
    ButtonName = 'No';
  else,    
    ButtonName ='Yes';  
  end
  %ButtonName =  questdlg('Go To WAFO Toolbox Home Page?', 'WWW', 'Yes', 'No', 'No'); 
end



if nargout > 0,
  URLout = URL;
else
  disp(['WAFO Toolbox Help Page:'])
  disp(URL)
end  
 
if strncmpi(ButtonName, 'y',1)
  
  Status = web(URL,'-browser');
  
  switch Status
    case 1,
      disp(' Web Browser not found.')
      disp(' See "help web".')
      help(web)	
    case 2,
      disp(' Web Browser found, but could not be launched.')
      disp(' See "help web".')
      help(web)
    otherwise
  end    
end
return
