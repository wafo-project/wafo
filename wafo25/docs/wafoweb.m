function URLout = wafoweb(ButtonName)
%WAFOWEB  World Wide Web site of the WAFO Toolbox.
%
%  CALL:  URL = wafoweb(BN)
% 
%   URL = string with the URL address of WAFO.
%   BN  = 'yes' If display URL site
%         'no'  If no display of URL site
%
% WAFOWEB displays and/or returns the WWW
%   site for the WAFO Toolbox.  If no input argument 
%   is givendisplayed, a dialog asks whether to go there.
%
% Examples
%  url = wafoweb('n'); web(url)
%  wafoweb 
%  wafoweb('yes')
%
% See also  web
 
% Tested on: Matlab 5.2
% History:
% by pab 19.06.2001

URL = 'http://www.maths.lth.se/matstat/wafo/';

if nargin<1 || isempty(ButtonName),  
  %if nargout>0, ButtonName = 'No';  else,    ButtonName ='Yes';  end
  ButtonName =  questdlg('Go To WAFO Toolbox Home Page?', 'WWW', 'Yes', 'No', 'No'); 
end



if nargout > 0,
  URLout = URL;
else
  disp(['WAFO Toolbox Home Page:'])
  disp(URL)
end  
 
if strncmpi(ButtonName, 'y',1)
  
  Status = web(URL);
  
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
