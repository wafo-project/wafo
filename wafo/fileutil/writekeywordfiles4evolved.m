function writekeywordfiles4evolved(mdirs)
%writekeywordfiles4evolved
% 
% CALL: writekeywordfiles4evolved(mdirs)
%
% mdirs     = cellarray or character array with directories to search for
%             M-files. (default '.')
% 
% SkipExample % 
%   writekeywordfiles4evolved('signal')
%
% See also

% References
%

% Tested on: Matlab 7.1
% History:
% 

% Author: pab
% Date: 12.07.2006
% Time: 15:02:51
% Version: 1.0


if nargin<1
	mdirs = '.';	
end

filename = 'matlab';
fid = fopen(sprintf('%s.kw', filename),'w');
fprintf(fid,'Name: Matlab\n');
fprintf(fid,'Description: Matlab m-files\n');
fprintf(fid,'LineComment: \%% \n');
fprintf(fid,'Extensions: m\n');
fprintf(fid,'CaseSensitive: 1\n');
fprintf(fid,'StringChar: %s\n' ,'''');
fprintf(fid,'WholeWord: 1 \n');
fprintf(fid,'Keyword1:');

keywords = {...
    'break',
    'case',
    'catch',
    'continue',
    'else',
    'elseif',
    'end',
    'for',
    'function',
    'global',
    'if',
    'otherwise',
    'persistent',
    'return',
    'switch',
    'try',
    'while',
    };
fprintf(fid,'%s ',keywords{:});
recursive=1;
mfiles = getmfiles(mdirs, recursive);
funcNames = strrep(mfiles,'.m','');

fprintf(fid,'\nKeyword2:');
k = find(isfunction(funcNames) & ~ismember(funcNames,keywords));
if any(k)
	fprintf(fid,'%s ',funcNames{k});
end


fclose(fid)
fid = fopen(sprintf('%s.a',filename),'w');
for ix = k(:).'
[T,R] = strtok(geth1line(funcNames{ix}));
fprintf(fid,'%s %s\n',funcNames{ix},R);
end
fclose(fid)

