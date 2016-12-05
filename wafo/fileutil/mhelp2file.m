function mhelp2file(rpwd,file)
%MHELP2FILE Writes all help texts to file
%
% CALL: mhelp2file(rpwd,file)
%
% rpwd  = path to root directory (Default current directory).
% file  = name of outfile
%
% Example
% mhelp2file(waforoot)
% 

error(nargchk(0,2,nargin))
if nargin<1||isempty(rpwd), rpwd = pwd;end
if nargin<2||isempty(file), file = 'helptexts.txt';end


edirs = findsubdir(rpwd,3);

[k2,m]=size(edirs);

[fid,msg] = fopen(file,'r');          % Open file to read
if fid==-1, disp(msg), return, end

   
for ix = 1:k2,
   DD = what(deblank(edirs(ix,:)));
   k3 = length(DD.m);
   for iy=1:k3
     fprintf(fid,'%s \n \n',DD.m{iy});
     str = help(DD.m{iy});
     fprintf(fid,str);
   end
end  
fclose(fid); 
