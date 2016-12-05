function mexamples2file(rpwd,file)
%MEXAMPLES2FILE Writes all examples in help headers to file
%
% CALL: mexamples2file(rpwd,file)
%
% rpwd  = path to root directory (Default current directory).
% file  = name of outfile
%
% SkipExample
% filename = tempname();
% mexamples2file(fullfile(waforoot, 'fileutil'), filename)
% 

error(nargchk(0,2,nargin))
if nargin<1||isempty(rpwd), rpwd = pwd;end
if nargin<2||isempty(file), file = 'example.m';end


edirs = findsubdir(rpwd,3); % fileutil function

[k2,m]=size(edirs);

if ischar(file) && exist(fullfile(pwd,file),'file') ~= 0, 
  disp(['There allready exist a ', file , ' in this directory'])
  disp(['Copied this file to ', file  'old'])
  %Make backup copy before overwriting file
  copyfile(file,[file 'old']);
  delete(file)
end



[fid,msg] = fopen(file,'w');          % Open file to read
if fid==-1, disp(msg), return, end

   
for ix = 1:k2,    %loop trough all  directories
   DD = what(deblank(edirs(ix,:)));
   k3 = length(DD.m);
   for iy=1:k3,    % loop through all files
     
     str = help(DD.m{iy});          % Extract help header
     if ~isempty(str)
       inl = findnl(str);             % Find newlines
       [names,p0 p1] = getnames(str); % Extract legal names
       names2 = lower(names);
       ind1 = strmatch('example',names2); % Find token string
       if ~isempty(ind1),
         ind1=ind1(end);               % Keeping only the last token
         ind2 = strmatch('see',names2(ind1:end,:));
         % Find the lines after the token
         ind0 =  find(inl-p0(ind1)>0);
         
       if ~isempty(ind0)
         fprintf(fid,'\n %% ----------------------------------------------\n');
         fprintf(fid,' %% HELP TEXT EXAMPLES FOR %s : \n \n',fullfile(DD.path,DD.m{iy}));
         ind0 = inl(ind0(1))+1; % First line after the token
         if isempty(ind2),
           str(ind0:end)
           fprintf(fid,str(ind0:end));
         else
           ind1 = p0(ind1+ind2-1)-1;
           fprintf(fid,str(ind0:ind1));
         end
       end
       
     end
   end
 end
end  
fclose(fid); 
