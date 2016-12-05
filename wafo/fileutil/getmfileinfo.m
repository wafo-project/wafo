function  info1 = getmfileinfo(mfile,varargin)
%GETMFILEINFO Return info about mfile such as subroutines, functions called etc.
%
% CALL: info = getmfileinfo(mfile,options);
% 
% info    = structure containing info of mfile with the field names:
%   .file_name            = name of file
%   .source               = source listing of m.file
%   .synopsis             = function synopsis line of the m-file. (empty if
%                          MFILE is a script)
%   .is_script            = 1 if MFILE is a script, 0 otherwise
%   .is_class_function    = 1 if MFILE is class function, 0 otherwise
%   .is_class_constructor = 1 if MFILE is class , 0 otherwise
%   .h1_line          = first comment line in help header of the m-file.
%   .help_header      = help header of m-file
%   .global           = cell array of global variables used in m-file.
%   .persistent       = cell array of persistent variables used in m-file.
%   .functions_called = character array of functions called by the m-file.
%   .todo             = structure containing information about todo tags
%   .bugs             = structure containing information about bugs tags                                                                  
%   .subroutine       = structure containing info about each subroutine. 
%        .names       = cellarray of names of the subroutines
%        .synopsis    = cellarray of function synopsis of each subroutine.
%        .h1_lines    = cellarray of first comment line in help header of
%                       each subroutine.
%        .hrefs       = boolean array with cross reference info for the 
%                       subroutines. hrefs(i+1,j+1) = 1 if 
%                         .subroutine.names{i}  calls  .subroutine.names{j} 
%                       (Column 1 / row 1 is the main function, 
%                        Column 2 / row 2 is sub function 1, 
%                        Column 3 / row 3 is sub function 2, ... etc).
%                       
% mfile       = name of m-file to load and dissect.  
% options     = strings defining standard function names to be included in
%                output. Allowed options are: 
%                'trivia', 'math', 'controls' or 'exim' (see fnames)
%
%  GETMFILEINFO loads the given file and dissects it into help header, H1
%  lines, subfunctions and functions called etc.  Some simple functions and
%  language commands such as 'ones', 'rand', 'sin', 'for' , 'else', etc. are 
%  normally excluded. To include them the additional arguments such as 
%  'trivia', 'math', 'controls' can be used.
%  See FNAMES function for details.
%
% Example: % Return functions (including all "math" functions) used in
%          % the getmfileinfo function
%
%  info1 = getmfileinfo('getmfileinfo.m', 'math','trivia');
%  assert(info1.todo.comment{1}(1:30), ' Unable to extract/distinguish')
%  assert(info1.synopsis{1}, 'info1 = getmfileinfo(mfile,varargin)');
%  assert(info1.subroutine.names{1}, 'removedoubleblanks');
%  assert(info1.subroutine.synopsis{1}, 'str = removedoubleblanks(str)')
%  assert(info1.functions_called(2:5)', {'all','any','brackets','cell'});
% 
% See also: fnames


% History:
% revised pab July 2006
% - changed name from parsemfilestr to getmfileinfo
% -added cross reference matrix for subfunctions.
% revised pab july 2005
% fixed a bug: code crashed when names is empty in the last setdiff(..)
% statement, now fixed
% revised pab 13.11.2003
% revised pab 11.11.2003
% by pab 16.10.2003
% - Based on refer.m by
%  Kirill K. Pankratov, kirill@plume.mit.edu
%  02/06/95

% TODO % Unable to extract/distinguish variable/function names in eval or feval statements
% BUGS % Mistakes a variable to be a function in scripts if it is initialized outside the script
% BUGS % Variables in nested functions or subfunctions are not handled properly if the name clashes with function names.
% BUGS % If Parts of code are never reached on execution, these parts may cause wrong names to be extracted as used in code.
% BUGS % functions called must be in path to be detected

% Handle input ..................................
error(nargchk(1,6,nargin))


info1 = struct('file_name',mfile,'source',freadtxt(mfile),...
  'synopsis','',...
  'is_script',1,...
  'is_class_function', 0,...
  'is_class_constructor', 0,...
  'h1_line','',...
  'help_header','',...
  'global','',...
  'persistent','',...
  'todo',[],...
  'bugs','',...
  'subroutine','',...
  'functions_called','');

[mfilepath,mfileName] = fileparts(which(mfile));
[dummy,mfileDir] = fileparts(mfilepath);
isclassfunction = any(mfilepath=='@');
if ~isempty(mfileDir) && mfileDir(1)=='@'
  isclassconstructor = strcmp(mfileDir(2:end),mfileName);
else
  isclassconstructor = 0;
end

info1.global = cell(1,0);
info1.persistent = cell(1,0);
todo        = struct('line',[],'comment',{{}});
info1.todo  = todo;
info1.bugs  = info1.todo;


% initialise output

subroutine = struct('names','','synopsis','','h1_lines','','hrefs', nan);
subroutine.names    = cell(1,0);
subroutine.synopsis = cell(1,0);
subroutine.h1_lines = cell(1,0);

%names       = '';
%helpHeader  = '';
%H1Line      = '';

%globalVariableNames = cell(1,0);

LF          = char(10); % Linefeed character
space       = ' '; 

startBlockKeywords = char('if','switch', 'try', 'while','for');
%- List of Matlab keywords (output from iskeyword)
matlabKeywords = char(startBlockKeywords,'break', 'case', 'catch', 'continue', 'elseif', 'else', ...
      'end',  'function', 'global',  'otherwise', ...
      'persistent', 'return' );


In_str = char('trivia','math','controls','exim');
In_n = []; % Default index to standard function names to include
for jj = 1:nargin-1,
  c_arg = varargin{jj};
  if ischar(c_arg),
    % String - names of matrices
    ind = strmatch(lower(c_arg),In_str);
    if any(ind);
      In_n = [In_n, ind.'];
    end
  end
end

string = dm2unix([space info1.source(:)' LF space],space);     % Add blank edges and convert to unix line endings.
lstr   = length(string);

 % Quotes and comments masks ...............
[mask_q,mask_c,i_LF,lineNumber] = findquot(string);

% Blank the quotes
string_q     = string;
%nn           = find(mask_q);
string_q(mask_q) = space;

% Find help header and comments
nn = find(mask_c);
if any(nn),
 % if nargout>3
    nn2 = [0,find(diff(nn)>1),length(nn)];
  
    info1.help_header = string(nn(nn2(1)+1):nn(nn2(2)));
  
    indLF = find(info1.help_header==LF);
    indP  = intersect(find(info1.help_header=='%'),[1, indLF+1]);
    
    info1.help_header(indP) = space; % blank the comment characters in HH
    if any(indLF),
      info1.h1_line = removedoubleblanks(info1.help_header(1:indLF(1)-1));
    end
  
   % if nargout>5
      no = findstr(string_q,'% TODO %');
      no = intersect(no,[1,i_LF+1]); % Dismiss all TODO comments not starting at beginning of line
      if any(no)
         info1.todo.line    = lineNumber(no);
         nc           = i_LF(info1.todo.line); % indices to end of line characters
         todoComments = extract(string,no+8,nc);
         Nc           = size(todoComments,1);
         info1.todo.comment(1:Nc) = cellstr(todoComments);
      end
      no = findstr(string_q,'% BUGS %');
      no = intersect(no,[1,i_LF+1]); % Dismiss all BUGS comments not starting at beginning of line
      if any(no)
         info1.bugs.line    = lineNumber(no);
         nc           = i_LF(info1.bugs.line); % indices to end of line characters
         bugsComments = extract(string,no+8,nc);
         Nc           = size(bugsComments,1);
         info1.bugs.comment(1:Nc) = cellstr(bugsComments);
      end
    %end
  %end % nargout > 3
  % Blank the comments
  string0 = string;
  string0(nn)  = space;  
  string_q(nn) = space;
  % but keep newlines
  string0(i_LF)  = LF;
  string_q(i_LF) = LF;
end

% Blank the openings of quotes (strings) ........
nn         = find(diff(mask_q)==1)+1;
string0(nn) = space;%(ones(size(nn)));

%string = char(string);

if ~any(string0>space); 
  info1.subroutine = subroutine;
  % If nothing left, return ......................
  return;
end

i_exclamation = find(string_q=='!');
if any(i_exclamation)
  % Call to the computer system System call
  % blank the statement
  for ix = 1:length(i_exclamation)
    iz = lineNumber(i_exclamation(ix));
    string_q(i_exclamation(ix):i_LF(iz)-1) = space;
  end
end
 % Determine bracketing .....
[brro,brrc,mask_brr,brrMisMatch] = brackets(string_q,'('); % parentheses
[brso,brsc,mask_brs,brsMisMatch] = brackets(string_q,'['); % Brackets
[brbo,brbc,mask_brb,brbMisMatch] = brackets(string_q,'{'); % Braces

if brrMisMatch || brsMisMatch || brbMisMatch
  warn_txt = sprintf('%s is probably not a valid m-file because one of the parantheses, brackets or braces do not match!\n Use mlint to check the m-file.',mfile);
  warning('fileutil:getmfileinfo:bracketsMisMatch',warn_txt)
end

 % Find end of statements, i.e., LF, semicolon or 
 % comma character which are not inside brackets.
i_end              = unique([i_LF,find(string_q==','|string_q==';'),lstr]);
notInsideBrackets  = (mask_brr(i_end)==0 & ...
                          mask_brs(i_end)==0 & ...
                          mask_brb(i_end)==0);
i_end              = i_end(notInsideBrackets);
i_start            = [1 i_end(1:end-1)+1];

statementList      = extract(string0,i_start,i_end+1);
nonEmptyStatements = (any(statementList>space,2));
clear statementList
% Keep only non-empty statements
i_end   = i_end(nonEmptyStatements);
i_start = i_start(nonEmptyStatements);

% merge empty statements to the preceding statement
i_end = max(i_end,[i_start(2:end)-1,lstr]);

% find continuation lines
i_continuationLines   = findstr(string_q,'...');
if any(i_continuationLines)
  notInsideBrackets   = (mask_brr(i_continuationLines)==0 & ...
                         mask_brs(i_continuationLines)==0 & ...
                         mask_brb(i_continuationLines)==0);
  i_continuationLines = i_continuationLines(notInsideBrackets);
  
  if any(i_continuationLines)
    % remove end statements from continuation lines.
    [dummy,ix] = intersect(lineNumber(i_end),lineNumber(i_continuationLines));
    i_end(ix)  = [];
  end
end

 % Statement numbering ...............
statementNumbers          = zeros(1,lstr);
statementNumbers(i_end+1) = ones(size(i_end));
statementNumbers(1)       = statementNumbers(1) + 1;
statementNumbers          = cumsum(statementNumbers);
numStatements             = statementNumbers(lstr);  % Total nmb. of statements


if 0
   disp('Numbering statements')
   for ix = 1:numStatements
      txt = string_q(statementNumbers==ix);
      if any(txt==LF)
         fprintf('%d %s',ix,txt)
      else
         disp(sprintf('%d %s',ix,txt))
      end
   end
end

% Find assignment operator for each statement.........
% (exclude those in comments and quotes)
i_eq = find(string_q == '=');  % All '=' characters

if ~isempty(i_eq) 
  % Auxillary ..................
  % Characters before '=' in non-assignment.
  relationalOperators = ('~=><')';  
  Nr     = length(relationalOperators); 
  onesRO = ones(Nr,1);
  ol     = ones(1,length(i_eq));
  % Separate logical from assignment operators
  logicalOperators = ...
     ((any(string_q(onesRO,i_eq-1)==relationalOperators(:,ol)) | ...
     any(string_q(onesRO,i_eq+1)==relationalOperators(:,ol))));
  % Keep indices to assignment operators
  i_eq(logicalOperators) = [];  
end

% Statement ## for assign. operators ................
i_assg = zeros(1,numStatements);
i_assg(statementNumbers(i_eq)) = i_eq;
%i_assg = sparse(1,statementNumbers(i_eq),i_eq,1,numStatements)

% Extract all variable and function names ...........
%Ltype = 'cell';
Ltype = 'matrix';
[names,no,nc] = getnames(string_q,Ltype);


 % Determine which side of the assignment operator ...
leftSide = (no<i_assg(statementNumbers(no))); 
% 1-left-hand side, 0-right-hand side

% Those on left hand side are variable names.
nn = find(leftSide);

if any(nn)
  variablenames = unique(names(nn,:),'rows');
else
  variablenames = '';
end
nn = find(string_q(no-1)=='.');
if any(nn)
  %names preceded by '.' is a struct name; blank them
  names(nn,:) = space;
end


i_globalWords = strmatch('global', lower(names),'exact');
if any(i_globalWords)
  % find GLOBAL variable names on the right hand side of 
  % GLOBAL statements if any.
  
  globalStatements = statementNumbers(no(i_globalWords));
  
  numNames = size(names,1);
  
  for ix = 1:length(globalStatements)
    currentStatement = globalStatements(ix);
    % index to next word
    iy = i_globalWords(ix) + 1; 
    while ((( iy <= numNames  ) && ...
	    ( currentStatement == statementNumbers(no(iy)) )))
      % parse global variable names
      info1.global{end+1} = deblank(names(iy,:));
      iy = iy + 1;
    end
  end 
end


i_persistentWords = strmatch('persistent', lower(names),'exact');
if any(i_persistentWords)
  % find PERSISTENT variable names on the right hand side of 
  % PERSISTENT statements if any.
  
  persistentStatements = statementNumbers(no(i_persistentWords));
  
  numNames = size(names,1);
  
  for ix = 1:length(persistentStatements)
    currentStatement = persistentStatements(ix);
    % index to next word
    iy = i_persistentWords(ix) + 1; 
    while ((( iy <= numNames  ) && ...
	    ( currentStatement == statementNumbers(no(iy)) )))
      % parse global variable names
      info1.persistent{end+1} = deblank(names(iy,:));
      iy = iy + 1;
    end
  end 
end





% Function statements
i_functionWords = strmatch('function', lower(names),'exact');
if any(i_functionWords)
  % find variable names on the right hand side of 
  % function statements if any.
  
  funStatements = statementNumbers(no(i_functionWords));
  index2functionAssignmentOps = i_assg(funStatements); 
  % index zero means no assignment operator in statement;
  
  numNames = size(names,1);
  
  i_functions = [];
  numFunctionStatements = length(funStatements);
  
  for ix = 1:numFunctionStatements
    currentStatement = funStatements(ix);
    % index to next word
    iy = i_functionWords(ix) + 1; 
    
    thisAssignmentOp = index2functionAssignmentOps(ix);
    
    % loop until we get to the function name      
    while ( (iy <= numNames) && ( no(iy) < thisAssignmentOp ) )
      % index to word < index to assignment operator
      iy = iy + 1;
    end
    if (iy <= numNames)
      
      % Extract synopsis for m-file and its subfunctions     
      iz = find(statementNumbers(i_end)==currentStatement);
      synopsistxt     = string_q(nc(i_functionWords(ix))+1:i_end(iz)-1);
      synopsistxt     = strrep(deblank(strrep(synopsistxt,'.','')),LF,'');
      
      i_functions(ix) = iy;
      if (ix>1)
        % Extract subroutine names.
        subroutine.names{ix-1}     = deblank(names(iy,:));
        
        %subroutine.synopsis{ix-1} = rstrrep(deblank(synopsistxt),'  ',' ');
        subroutine.synopsis{ix-1} = removedoubleblanks(synopsistxt);
        
        %subroutines = strvcat(subroutines,names(iy,:));
      else
        info1.is_script = 0;
        info1.is_class_function   = isclassfunction;
        info1.is_class_constructor = isclassconstructor;
        %info1.synopsis{ix} = rstrrep(synopsistxt,'  ',' ');
        info1.synopsis{ix} = removedoubleblanks(synopsistxt);
      end
      % deblank the function- or subfunction-name from name list
      names(iy,:) = space;
    end
    iy = iy + 1;
    while ((( iy <= numNames  ) && ...
	    ( currentStatement == statementNumbers(no(iy)) )))
      % parse variable names
      variablenames = char(variablenames,names(iy,:));
      iy = iy + 1;
    end
    if (ix>1) 
      i_H1lineSub = find(lineNumber(no(iy-1))+1==lineNumber);
      potentialH1line = ddeblank(string(i_H1lineSub));
      if any(potentialH1line) && (potentialH1line(1)=='%')
        subroutine.h1_lines{ix-1} = potentialH1line(2:end);
      else
        subroutine.h1_lines{ix-1} = '';
      end
    end  
  end 
  % Make cross reference info for subfunctions
  subroutine.hrefs = sparse(numFunctionStatements,numFunctionStatements);
  for k =2:numFunctionStatements
    % See if in the list
    res = strmatch(subroutine.names{k-1},names,'exact');
    for ix=1:length(res) 
      distance= (res(ix)-i_functions);
      distance(distance<0) = inf;
      if any(isfinite(distance))
        %main function: iy=1 
        % subfunction1: iy=2;
        %...
        % subfunctionN: iy=N+1
        [tmp,iy] = min(distance);
        subroutine.hrefs(iy,k) = 1; 
      end
    end;
  end
end

i_keyWords = find(ismember(lower(names), cellstr(matlabKeywords)));
nn =   setdiff(find(i_assg==0), statementNumbers(no(i_keyWords)));
if any(nn)
   % statements with no assignment and no matlab keywords:
   % i.e., statements like: 
   %    pause on
   % or
   %   interp(x,y,xi)
   %   error(nargchk(2,3,nargin))
   %
   % where first word on each statement is a function and 
   % the remaining words is a mix of variablenames and/or functions.
   nn = find(ismember(statementNumbers(no),nn));
 
   nn = nn((diff([0 statementNumbers(no(nn))])==0)); % index to variable names
   
   for ix=length(nn):-1:1
     thisName  = deblank(names(nn(ix),:));
     % 
     try
       status =  exist(thisName,'file') + exist(thisName,'builtin');
     catch
       status = exist(thisName);
     end
     if (2<=status) &&  (status<=6)
       % It is a function name in path: remove from list
       nn(ix) = [];
     end
   end
   if any(nn),
      variablenames = char(variablenames,...
                          setdiff(names(nn,:), subroutine.names));
   end
end
% Choose only those on the right-hand side
% which definitely is not a variable name

names = setdiff(names(~leftSide,:),cellstr(variablenames));
if any(i_exclamation)
  names = char(names,'system');
end
if ~isempty(names),
  % or a global variable
  names = setdiff(names,info1.global);
  if ~isempty(names),
    % or a persistent variable
    names = setdiff(names, info1.persistent);
    if ~isempty(names)
      % or subfunction name
      names = setdiff(names, subroutine.names);
    end
  end
end

% numSubroutines = length(subroutines);
% if numSubroutines>0
%   numStartBlocks = sum(ismember(names,startBlockKeywords,'rows'));
%   numEndBlocks   = sum(ismember(names,'end','rows'));
%   isFunctionNested = (numEndBlocks==numStartBlocks + numSubroutines+1);
%   
% end
 
i0       = ones(1,4);
i0(In_n) = 0;
if any(i0) && ~isempty(names),
  % Compile "trivial" or special names ............  
   trivialFunctionNames = cell(1,4);
  [trivialFunctionNames{:}] = fnames;
  %[trivia,math,controls,exim] = fnames;
  names2remove = char(trivialFunctionNames(find(i0)));
  % Remove standards: "trivia" such as "zeros", "ones"
  % and "controls" such as "for", "if", "end"
  names = setdiff(names, cellstr(names2remove));
end
if any(names(:)) && all(names(1,:)<=space)
  names(1,:) = [];
end
info1.functions_called = names;
%if length(subroutine.names)>0
  info1.subroutine = subroutine;
%end
return





function str = removedoubleblanks(str)
%REMOVEDOUBLEBLANKS 
%
% 


k = isspace(str);
if any(k)
isdoubleblank = [k(1:end-1) == k(2:end) & k(2:end)==1, k(end)];
if any(isdoubleblank)
  str(isdoubleblank) = [];
end
end
if any(str) && isspace(str(1))
  str(1) = [];
end
