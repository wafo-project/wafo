function  [names,exst] = findref(filein,nex,varargin)
%FINDREF Find references to other functions in a file.
%
% CALL: [names,exst] = findref(file,nex,str1,str2,str3);
%    
%  names = character array containing function names used in 
%          the file.  
%  exst  = vector of the names exist status (from 1 to 5)
%  file  = string containing the name of the file to be searched.
%  nex   = specifies the exist status of functions to be returned.
%          (from 1 to 5), (default 2:5)
% str1,...
%  str3  = strings defining standard function names to be included in
%          output. Allowed options are: 
%         'trivia', 'math', 'controls' or 'exim' (see fnames)
%
%  FINDREF reads text from the file, parses it for all legal names
%  (excluding comments) and returns string matrix
%  Some simple functions and language commands
%  such as 'ones', 'rand', 'sin', 'for' , 'else',
%  etc. are normally excluded. To include them
%  the additional arguments such as 'trivia',
%  'math', 'controls' can be used.
%  See FNAMES function for details.
%
% Example: % Return only built-in and MEX functions
%   % (including all "math" functions) used in the
%   % findref function
%  [folder, root] = fileparts(waforoot);
%  if strcmpi(root,'wafo'),
%     ref = findref('findref',3:5,'math'); 
%     assert(ref, strvcat('error', 'exist', 'floor'));
%  end
% 
% See also: fnames, exist


% History:
% revised pab 16.10.2003
%  - streamlined some code and moved it into parsemfilestr    
% revised pab 05.10.2000
% - updated documentation
% - changed name from refer.m to findref.m due to the naming conflict
%   with the name of the directory refer
%  Kirill K. Pankratov, kirill@plume.mit.edu
%  02/06/95

 % Handle input ..................................

error(nargchk(1,6,nargin))
if nargin<2||isempty(nex)
  nex = 2:5; % Default exist status
end

% Read text from a file .........................
fn     = strtok(filein,'.'); % strip .m ending
string = freadtxt([fn '.m']);

[names,synopsis,subroutines] = parsemfilestr(string,varargin{:});


szn = size(names,1);

 % Exist status ...............................
msg = ' Extracted all names. Checking EXIST status...\n';
fprintf(msg)

exst = zeros(szn(1),1);
for jj = 1:szn(1)
  name = deblank(names(jj,:));  % Current name
  exst(jj) = exist(name);
  e = 2.1*(exist([name '.mat'])>0);
  exst(jj) = max(exst(jj),e);
end

 % Choose only function or file names .........
nn = zeros(size(exst));
for jj = nex;
  nn = (nn | floor(exst)==jj);
end
nn    = find(nn);
names = deblank(names(nn,:));
exst  = exst(nn);
return





