function  history2file(rpwd,outfile,files2include)
%HISTORY2FILE Writes all History notes in M-files to file
%
% CALL: history2file(rpwd,outfile,files2include)
%
% rpwd          = path to root directory (Default current directory).
% outfile       = name of outfile, or an integer file identifier obtained from FOPEN. 
%                It can also be 1 for standard output (the screen) or 2 for
%                standard error.           (default 'history.txt')
% files2include = list of files to include in history file.
%                 possibly with (Default '-all')
%
% Example
% history2file(fullfile(waforoot,'onedim'),[],'gaussq.m')
% 
error(nargchk(0,3,nargin))
if nargin<1||isempty(rpwd),
  rpwd = pwd;
end
if nargin<2||isempty(outfile),
  outfile = 'history.txt';
end

if nargin<2||isempty(files2include),
  %files2include = '-all';
  includeAll = true;
else
  includeAll =   any(strcmpi(files2include,'-all'));
  if ~includeAll
    if ischar(files2include)
      files2include = cellstr(files2include);
    end
    Nf = numel(files2include);
    pathstr = cell(1,Nf);
    for ix = 1:Nf
      [patstr{ix}, name, ext] = fileparts(files2include{ix});
      files2include{ix} = [name '.m'];
    end
  end
end


edirs = rpwd;

if ischar(edirs)
  edirs = cellstr(edirs);
end
%edirs = findsubdir(rpwd,0); % fileutil function



if isnumeric(outfile)
  fid = outfile;
else
  if ischar(outfile) && exist(fullfile(pwd,outfile),'file') ~= 0, 
    disp(['There allready exist a ', outfile , ' in this directory'])
    disp(['Copied this file to ', outfile  'old'])
    %Make backup copy before overwriting file
    copyfile(outfile,[outfile 'old']);
    delete(outfile)
  end

  [fid,msg] = fopen(outfile,'w');          % Open file to read
  if fid==-1,
    disp(msg),
    return
  end
end

k2= numel(edirs);

for ix = 1:k2,    %loop trough all  directories
   DD = what(edirs{ix});
   if any(findstr(edirs{ix},'spec'))
     disp('spec')
   end
   if includeAll || isempty(DD(1).m)
      k3 = 1:length(DD(1).m);
   else
     k3 = find(ismember(DD(1).m,files2include));
   end
  
   for iy=k3(:).',    % loop through selected files
     currentFile = fullfile(DD(1).path,DD(1).m{iy});
     str = gethistory(currentFile);          % Extract history
     fprintf(fid,'\n\n %% ----------------------------------------------\n');
     fprintf(fid,' %% HISTORY NOTES FOR %s : \n \n',currentFile);
     if isempty(str)
       fprintf(fid,'%s\n','History does not exist');
     else
       fprintf(fid,'%s',str);
     end
   end
end
if ~isnumeric(outfile)
  fclose(fid);
end

return  
function str = gethistory(mfile)  
%GETHISTORY Returns revision HISTORY notes from m-file.
%
% CALL: str = gethistory(mfile)
%    
error(nargchk(1,1,nargin))
filename = [strtok(mfile,'.') '.m'];
string = freadtxt(filename);

% initialise output
str = '';

LF          = char(10); % Linefeed character
space       = ' '; 

string = [space string(:)' space];  % Add blank edges
string = dm2unix(string,space);     % convert to unix line endings.
%lstr   = length(string);

 % Quotes and comments masks ...............
[mask_q,mask_c,i_nl,lineNumber] = findquot(string);

% Blank the quotes
string_q     = string;

string_q(mask_q) = space;%(ones(size(nn)));

% Find help header and comments
nn = find(mask_c);
if any(nn),
  nn2 = [0 find(diff(nn)>1), length(nn)];
  if length(nn2)>2, 
    [names,no,nc] = getnames(string_q);
    ind = [strmatch('history',lower(names));
	   strmatch('revised',lower(names))];
    
    if any(ind)
     for ix = 2:length(nn2)-1
       lo = nn(nn2(ix)+1);
       up = nn(nn2(ix+1));       
       if any( (lo<=  no(ind)) & (no(ind)<=up)),
         str = [str,LF,string(lo:up)];
       end
     end
    end
  end  
end