function  [names,synopsis,subroutines,helpHeader,...
	   H1Line,todo,globalVariableNames,H1lineSubroutines,subHrefs] =...
    parsemfilestr(string,varargin)
%PARSEMFILESTR Returns function names used in m-file string.
%
% CALL: [names,synopsis,subroutines, HH, h1line,todo] = parsemfilestr(string,options);
%    
% names       = character array containing function names used in 
%                the m-file.
% synopsis    = cellarray with function synopsis line of the m-file and
%               its subfunctions.
% subroutines = cellarray with name of the subroutines in the m-file.
% HH          = help header
% h1line      = first comment line in help header of the m-file.
% todo        = structure containing information about potential todo
%               tags
% string      = m-file string to dissect.  
% options     = strings defining standard function names to be included in
%                output. Allowed options are: 
%                'trivia', 'math', 'controls' or 'exim' (see fnames)
%
%  PARSEMFILESTR parses STRING for all legal names
%  (excluding comments) and returns string matrix
%  Some simple functions and language commands
%  such as 'ones', 'rand', 'sin', 'for' , 'else',
%  etc. are normally excluded. To include them
%  the additional arguments such as 'trivia',
%  'math', 'controls' can be used.
%  See FNAMES function for details.
%
% Example: % Return functions (including all "math" functions) used in
%          % the dewhite function
% str =  freadtxt('dewhite'); 
% [names,synopsis,subroutines,HH,h1] = parsemfilestr(str,'math','trivia')
% assert(names, strvcat('cell', ¨'error','find','iscell','ischar','isempty',...
%       'isspace','max','nargchk','nargin','numel','size'))
%	 
% See also: fnames


% History:
% revised pab July 2006
% - added 
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
% BUGS % If Parts of code are never reached on execution, these parts may cause wrong names to be extracted as used in code.
% BUGS % functions called must be in path to be detected

% Handle input ..................................
%error(nargchk(1,6,nargin))
narginchk(1,6)
% initialise output
synopsis    = cell(1,0);
subroutines = cell(1,0);
H1lineSubroutines = cell(1,0);
names       = '';
helpHeader  = '';
H1Line      = '';
todo        = struct('line',[],'comment',{{}});
bugs        = todo;
globalVariableNames = cell(1,0);

LF          = char(10); % Linefeed character
space       = ' '; 

startBlockKeywords = strvcat('if','switch', 'try', 'while','for');
%- List of Matlab keywords (output from iskeyword)
matlabKeywords = strvcat(startBlockKeywords,'break', 'case', 'catch', ...
      'continue', 'elseif', 'else', 'end',  'function', 'global',  ...
      'otherwise', 'persistent', 'return');

In_str = strvcat('trivia','math','controls','exim');
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

string = [space string(:)' LF space];  % Add blank edges
string = dm2unix(string,space);     % convert to unix line endings.
lstr   = length(string);

 % Quotes and comments masks ...............
[mask_q,mask_c,i_LF,lineNumber] = findquot(string);

% Blank the quotes
string_q     = string;
string_q(mask_q) = space;

% Find help header and comments
nn = find(mask_c);
if any(nn),
  if nargout>3
    nn2 = [0,find(diff(nn)>1),length(nn)];
  
    helpHeader = string(nn(nn2(1)+1):nn(nn2(2)));
  
    indLF = find(helpHeader==LF);
    indP  = intersect(find(helpHeader=='%'),[1, indLF+1]);
    
    helpHeader(indP) = space; % blank the comment characters in HH
    if any(indLF),
      H1Line = helpHeader(1:indLF(1)-1);
    end
  
    if nargout>5
      no = findstr(string_q,'% TODO %');
      no = intersect(no,[1,i_LF+1]); % Dismiss all TODO comments not starting at beginning of line
      if any(no)
         todo.line    = lineNumber(no);
         nc           = i_LF(todo.line); % indices to end of line characters
         todoComments = extract(string,no+8,nc);
         Nc           = size(todoComments,1);
         todo.comment(1:Nc) = cellstr(todoComments);
      end
      no = findstr(string_q,'% BUGS %');
      no = intersect(no,[1,i_LF+1]); % Dismiss all BUGS comments not starting at beginning of line
      if any(no)
         bugs.line    = lineNumber(no);
         nc           = i_LF(bugs.line); % indices to end of line characters
         bugsComments = extract(string,no+8,nc);
         Nc           = size(bugsComments,1);
         bugs.comment(1:Nc) = cellstr(bugsComments);
      end
    end
  end % nargout > 3
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
  % If nothing left, return ......................
  return;
end

 % Determine bracketing .....
[brro,brrc,mask_brr] = brackets(string_q,'('); % parentheses
[brso,brsc,mask_brs] = brackets(string_q,'['); % Brackets
[brbo,brbc,mask_brb] = brackets(string_q,'{'); % Braces
 % Find end of statements, i.e., LF, semicolon or 
 % comma character which are not inside brackets.
i_end              = unique([i_LF,find(string_q==','|string_q==';'),lstr]);
notInsideBrackets  = (mask_brr(i_end)==0 & ...
                          mask_brs(i_end)==0 & ...
                          mask_brb(i_end)==0);
i_end              = i_end(notInsideBrackets);
i_start            = [1 i_end(1:end-1)+1];

statementList      = extract(string0,i_start,i_end+1);
nonEmptyStatements = find(any(statementList>space,2));
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
[names,no,nc] = getnames(string_q);

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


globalWords = strmatch('global', lower(names),'exact');
if any(globalWords)
  % find GLOBAL variable names on the right hand side of 
  % GLOBAL statements if any.
  
  globalStatements = statementNumbers(no(globalWords));
  
  numNames = size(names,1);
  
  for ix = 1:length(globalStatements)
    currentStatement = globalStatements(ix);
    % index to next word
    iy = globalWords(ix) + 1; 
    while ((( iy <= numNames  ) && ...
	    ( currentStatement == statementNumbers(no(iy)) )))
      % parse global variable names
      globalVariableNames{end+1} = deblank(names(iy,:));
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
      synopsistxt     = string(nc(i_functionWords(ix))+1:i_end(iz)-1);
      synopsistxt     = strrep(deblank(strrep(synopsistxt,'.','')),LF,'');
      synopsis{end+1} = rstrrep(synopsistxt,'  ',' ');
      i_functions(ix) = iy;
      if (ix>1)
        % Extract subroutine names.
        subroutines{end+1} = deblank(names(iy,:));
        %subroutines = strvcat(subroutines,names(iy,:));
      end
      % deblank the function- or subfunction-name from name list
      names(iy,:) = space;
    end
    iy = iy + 1;
    while ((( iy <= numNames  ) && ...
	    ( currentStatement == statementNumbers(no(iy)) )))
      % parse variable names
      variablenames = strvcat(variablenames,names(iy,:));
      iy = iy + 1;
    end
    if (ix>1) 
      i_H1lineSub = find(lineNumber(no(iy-1))+1==lineNumber);
      if any(i_H1lineSub) && (string(i_H1lineSub(1))=='%')
        H1lineSubroutines{ix-1} = deblank(string(i_H1lineSub(2:end)));
      else
        H1lineSubroutines{ix-1} = '';
      end
    end  
  end 
  % Make cross reference info for subfunctions
  subHrefs = sparse(numFunctionStatements,numFunctionStatements);
  for k =2:numFunctionStatements
    % See if in the list
    res = strmatch(subroutines{k-1},names,'exact');
    for ix=1:length(res) 
      distance= (res(ix)-i_functions);
      distance(distance<0) = inf;
      if any(isfinite(distance))
        %main function: iy=1 
        % subfunction1: iy=2;
        %...
        % subfunctionN: iy=N+1
        [tmp,iy] = min(distance);
        subHrefs(iy,k) = 1; 
      end
    end;
  end
  
  
  
end

[names, matlabKeywords] = space___padd(names, matlabKeywords);


i_keyWords = find(ismember(lower(names),matlabKeywords,'rows'));
nn =   setdiff(find(i_assg==0),statementNumbers(no(i_keyWords)));
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
 
   nn = nn(find(diff([0 statementNumbers(no(nn))])==0)); % index to variable names
   
   for ix=length(nn):-1:1
      status = exist(deblank(names(nn(ix),:)));
      if any(intersect(status,[2:6]))
         % It is a function name in path: remove from list            
         nn(ix) = []; 
      end
   end
   if any(nn),
      variablenames = strvcat(variablenames,...
			      setdiff(names(nn,:),char(subroutines),'rows'));
   end
end
% Choose only those on the right-hand side
% which definitely is not a variable name
nn    = find(~leftSide);
[names, variablenames] = space___padd(names, variablenames);
names = setdiff(names(nn,:),variablenames,'rows');
if isempty(names),
  return;
end
[names, globalVariableNames] = space___padd(names, globalVariableNames);
% or a global variable
names = setdiff(names,char(globalVariableNames),'rows');
if isempty(names),
  return;
end
[names, subroutines] = space___padd(names, char(subroutines));
% or subfunction name
names = setdiff(names,subroutines,'rows');

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
  names2remove = strvcat(trivialFunctionNames{find(i0)});
  % Remove standards: "trivia" such as "zeros", "ones"
  % and "controls" such as "for", "if", "end"
  [names, names2remove] = space___padd(names, names2remove);
  names = setdiff(names,names2remove,'rows');
end
return
end

function [a,b] = space___padd(a, b)
  n = size(a, 2)-size(b, 2);
  if n>0 
     b = cat(2, b, repmat(' ', size(b, 1), n));
  elseif n<0
     a = cat(2, a, repmat(' ', size(a, 1), -n));
  end
end
