function sitemap(goUp)
% Web Site Map Generation
% Version 1.0
% Last modified 11.08.00
% For Matlab 5.x		 
% File sitemap.m						 
% Copyright (c) 2000 Igor Kagan					 
% kigor@tx.technion.ac.il
% http://igoresha.virtualave.net
%
% Usage: In Matlab Command Window, change to your site root directory and type:
% sitemap
%
% Description:
% This function uses <title>...</title> and
% <meta name="description" content="..."> html tags
% to fill entry for each html (htm) file.
% If these tags are not specified, sitemap wil just list the files
% and dirs in the current dir.
% Sitemap will create (and _overwrite_ if file with this name already
% exists!) the file Sitemap.html, which can be edited in any text
% editor later, if there is a need to do so.
% 
% Note that links in Sitemap.html will _not_ work while you are in your
% local hard drive! Once you upload this file in the root dir of your
% server, they will work.


 
global smc nfiles ndirs baseurl	% need access these variables during recursive calls 

if nargin < 1,

	% add to path
      	curdir = cd;
      	matlabpath = path;
      	if isempty(findstr(lower(matlabpath),lower(curdir))),
        	path(matlabpath,curdir);
      	end;

	% define startup parameters
	baseurl = pwd;
	goUp	= 0;
	smc	= []; % site map html contents
	nfiles	= 0;
	ndirs	= 1;  % first is root dir
	startentry = ['<html><head><title>Site Map</title><meta name="description" content="Generated site map"><LINK REL="stylesheet" TYPE="text/css" HREF="ss.css"></head><body><table><tr><td><b>Sitemap</b></td></tr>' sprintf('\n')];
	smc	= [smc startentry];
end

% finish prevoius dir
smc = [smc '<tr><td><HR ALIGN=LEFT NOSHADE size=1></td></tr>'];

% get current dir listing
d = dir;
dirname = [strrep(pwd,baseurl,'') '/']; % for root dir
dirname = strrep(dirname,'\','/');
	   
direntry = ['<tr><td>',dirname,'<HR ALIGN=LEFT NOSHADE size=1></td></tr>' sprintf('\n')];
ndirs = ndirs + 1;
smc = [smc direntry];
       
d = d(3:length(d));

[t,ind]=sort(cat(1,d.isdir));
d = d(ind);
for i=1:size(d,1),
	if d(i).isdir,
		eval(['cd ',d(i).name]);
		sitemap(1);
	else
		if findstr(d(i).name,'.htm'),
			disp(d(i).name);
			fidr=fopen(d(i).name,'r');
			F = fread(fidr);
			fclose(fidr);
			F = F((F~=10 & F~=13 & F~=9));
			S = char(F');
			s = lower(S);
			% find title
			i1 = findstr(s,'<title>'); % length 7
			i2 = findstr(s,'</title>');
			Title = S(i1+7:i2-1);
			% find description
			i1 = findstr(s,'<meta name="description" content="'); % length 34
			if ~isempty(i1),
				i2 = findstr(s(i1+1:length(s)),'">');
				Description = S(i1+34:i1+i2(1)-1);
			else
				Description = '(no content specified)';
			end

			% write file entry
			fileentry = ['<tr><td><a href="',dirname,d(i).name,'">',d(i).name,'</a><b> ',Title,'</b><br>',Description,'</td></tr>' sprintf('\n')];
			nfiles = nfiles + 1;
			smc = [smc fileentry];	     
		end
		
				       
	end

end
if (goUp),
	cd ..
else
	% finish html
	finishentry = ['<tr><td><br><br>',num2str(nfiles),' files in ',num2str(ndirs),' dirs<HR ALIGN=LEFT NOSHADE size=2>This site map was generated at ',date,' by <a href="Matlab/Matlab.html#sitemap">sitemap.m</a></td></tr></table></body></html>'];
	smc = [smc finishentry];
	% write Sitemap.html file
	if exist('Sitemap.html','file'),
		dos('del Sitemap.html')
	end
	fidw = fopen('Sitemap.html','wt');
	fwrite(fidw,smc,'char');
	fclose(fidw);
	       
end
return


