function t = freadtxt(file)
% FREADTXT Opens file, reads the data into a string and close file.
% 
% CALL:  txt = freadtxt(file);
%
%   txt  = character vector containing the file data
%   file = string containing the name of the file
%
% Example
% txt = freadtxt(fullfile(waforoot, 'fileutil', 'freadtxt.m'));
% assert(txt(1:27), 'function t = freadtxt(file)');
%
% See also: fread

% history:
% revised pab 15.05.2001
% - changed fread(fid,inf,'char')  to fread(fid,inf,'uchar')
% revised pab 04.10.2000

if nargin == 0,
  error(' File name is not specified.');
end

[fid,msg] = fopen(file,'r');          % Open file to read
if fid==-1,
  disp(sprintf('Something wrong with opening file: %s',file));
  disp(msg);
  t = '';
 return, 
end
t = fread(fid,inf,'*char').';   % Read file
fclose(fid);                    % Close file
return
