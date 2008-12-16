
% This file
str = freadtxt('nameconversionold.m');
names = cellstr(getnames(str));
names = reshape(names,2,numel(names)/2).';


oldn = names(:,1);
newn = names(:,2);

for ix=1:numel(oldn)  
  fulname = where(oldn{ix});
  lpath = fileparts(fulname{1});
  
  str = freadtxt(fulname{1});
  str2 = regexprep(str,oldn,newn,'preservecase');
  file = fullfile(lpath,[newn{ix} '.m']);
  [fid,msg] = fopen(file,'w');          % Open file to read
if fid==-1,
  disp(sprintf('Something wrong with opening file: %s',file))
  disp(msg),
  t = '';
else
  fprintf(fid,'%s',str2);   % write file
end
fclose(fid);       
end