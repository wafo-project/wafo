function  diary2m(fileout,filein)
% DIARY2M  Converts diary of matlab session into executable m-file.
%
% CALL: diary2m(fileout,filein)
%
%   filein  = diary file to clean up (default get(0,'diaryfile'))
%   fileout = the corresponding matlab executable of filein.
%             (default 'diaryex.m')
%
% DIARY2M processes the "diary" file FILEIN, 
% removes MATLAB prompt ">>", displayed results
% and messages so that only issued commands themselves
% are left. The result is written into a file FILEOUT.
%
% Warning: the output file may still contain errors and
% may need some additional editing before execution.
%
% SkipExample
%  diary_name = tempname();
%  diary(diary_name);
%  x = linspace(0,10);
%  diary off;
%  mfile = [tempname() '.m'];
%  diary2m(mfile, diary_name)
%  txt = freadtxt(mfile);
%  assert(ddeblank(txt)(1:19), ...
%         'x = linspace(0,10);')


% History
% revised pab 13.12.2000
% - updated documentation 
% - added nargchk
%  Kirill K. Pankratov,  kirill@plume.mit.edu
%  May 15 1994, August 13 1994.

 % Default file names ..........................................
fileindflt = get(0,'diaryfile');    % Default for input file
fileoutdflt = 'diaryex.m';            % Default for output m-file

 % Handle input ................................................
error(nargchk(0,2,nargin))
if nargin<2||isempty(filein),  filein  = fileindflt; end
if nargin<1||isempty(fileout), fileout = fileoutdflt; end

[fidi,msg] = fopen(filein);         % Open input file
if fidi<0, disp(msg), return, end

[fido,msg] = fopen(fileout,'wt+');  % Open output file
if fido<0, disp(msg), close(fidi); return, end

isbr = 0;            % Brackets
isnext = 0;          % Next line
morelines = 1;

 % Begin reading and processing lines from FILEIN ``````````````
while morelines      % While there are more lines in the file
  str = fgets(fidi);    % Get next line
  morelines = str(1)>0; % Is there more lines

  lstr = find(str(1:length(str)-1)~=' ', 1, 'last' );
  str = [str(1:lstr) '  '];
  lstr = length(str);
  shift = all(str(1:2)=='>>');
  isex = shift||isnext;
  if isex            % If executable, write into output
    fnd = find(str=='%');       % Find comments
    numbr = sort([1:lstr fnd]); % Copy comments character
    str = str(numbr);
    lstr = length(str);
    fprintf(fido,['\n' str(1+2*shift:lstr-2)]);
  end
  iscomm = cumsum(str=='%');
  numbr = ((str=='[')-(str==']')).*(iscomm==0);
  numbr = cumsum(numbr); % Open or close brackets
  isbr = isnext*isbr+numbr(lstr)*isex;

  % Is next line a continuation of executable line .....
  lstr = find(~iscomm&str~=' ', 1, 'last' );
  fnd = max(1,lstr-2:lstr);
  isnext = isbr || ~isempty(fnd) && all(str(fnd)=='...');
end                 % No more lines    '''''''''''''''''''''''
fprintf(fido,'\n');

fclose(fidi);       % Close input file
fclose(fido);       % Close output file

