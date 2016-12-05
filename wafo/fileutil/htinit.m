function iniread()
% HTINIT  Initialize the HTMLTOOLS
% 
% CALL:  htinit
%
% HTINIT reads contents in the initialization file
% HTINI.DAT and saves it in the global structure HTINI
%
% The global struct HTINI contain the parameter fields:
% 
% name        = a string with user name indicated in ini-file
% email       = the e-mail address of user
% par         = a parameter vector
% color_sig1  = color of first line in signature file
% color_email = color of the string related to the e-mail address
% font        = type of font
% fontsize    = size of font
% html_home   = Home directory for the HTML toolbox
% cp          = Name of gif file beside user name (usually copyright.gif)
% bgfile      = background gif file for HTML
% gif_home    = Home directory of the Gif files
% color_idx   = Vector that defines the colour of the header in the 
%               index file, i.e. it will affect both the "INDEX FILE" 
%               string and the string containing the directory name.
%
% Example: View the global variables
%  htinit; global HTINI
%  HTINI
%
% See also: htini.dat

global HTINI

HTINI=struct('name',[]);

fid = fopen('htini.dat');
inp = fgetl(fid);
HTINI.name = deblank(inp);

inp = fgetl(fid);
HTINI.email = deblank(inp);

inp = fgetl(fid);
HTINI.par = eval(deblank(inp));

inp = fgetl(fid);
HTINI.color_sig1 = eval(deblank(inp));

inp = fgetl(fid);
HTINI.color_email = eval(deblank(inp));

inp = fgetl(fid);
HTINI.font = eval(deblank(inp));


inp = fgetl(fid);
HTINI.fontsize = eval(deblank(inp));

inp = fgetl(fid);
HTINI.html_home = eval(deblank(inp));

inp = fgetl(fid);
HTINI.cp = eval(deblank(inp));



inp = fgetl(fid);
HTINI.bg = eval(deblank(inp));

inp = fgetl(fid);
HTINI.gif_home = eval(deblank(inp));

inp = fgetl(fid);
HTINI.color_idx = eval(deblank(inp));

fclose(fid);


