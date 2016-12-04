%FIXSEEALSO Remove all ':' from the See also line in all m-files

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% -*- Mode: Matlab -*- %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% fixSeeAlso.m --- Remove all ':' from the See also line in all m-files of WAFO
%% Author          : Per Andreas Brodtkorb
%% Created On      : Mon Sep 20 14:57:02 2004
%% Last Modified By: Per Andreas Brodtkorb
%% Last Modified On: Tue Sep 21 15:46:56 2004
%% Update Count    : 16
%% Status          : Unknown, Use with caution!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




proot = waforoot;
cd(proot)
w      = dir(proot); 
inddir = find(cat(1,w.isdir));
wdir   = {w(inddir).name};
for iy = 1:length(wdir)
  cd(wdir{iy})
  pwd
  s = dir('*.m');

  for ix=1:length(s)
    mfile = s(ix).name;
    [fid,msg] = fopen(mfile,'r'); 
    source = char(fread(fid,inf,'uchar'))';   % Read file
    fclose(fid);                         % Close file
    source = strrep(source,'See also:','See also ');
    source = strrep(source,'See also :','See also  ');
    fid = fopen(mfile,'w');
    fprintf(fid,'%s',source);
    fclose(fid);
  end
end