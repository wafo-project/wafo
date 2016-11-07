function [options]= parseoptions(varargin)
%PARSEOPTIONS Create or alter a OPTIONS structure.
%
%  CALL:  options = parseoptions(legalFnames,options,funcname,opts1,opts2,...,
%                                par1,val1,par2,val2,...);
%         options = parseoptions(options,opts1,opts2,...,
%                                par1,val1,par2,val2,...);  
%
%   options    = OUT: options structure in which the named 
%                     parameters have the specified values.  
%                IN:  options structure giving the legal names and
%                     default options  
%   legalFnames = character array giving the names of functions for
%                 which default values for the options structure could be
%                 extracted.
%   funcname   = string giving the name of the function for which default
%                values for the options structure should be extracted.
%                (Must be equal to one of the names in LEGALFNAMES.)  
%   opts1,
%   opts2..    = options structures
%   par1,par2..= strings identifying the parameter to alter
%   val1,val2..= corresponding values the parameters are altered to.
%   
%   PARSEOPTIONS combines the default options for a function given by
%   FUNCNAME with new options structures (OPTS1,OPTS2,...) and/or with
%   the named parameters (PAR1,PAR2,...) with the corresponding values
%   (VAL1, VAL2,...).
%   The parameters are set in the same order as the input arguments.
%   Any parameters with non-empty values of the options struct overwrite
%   the corresponding old parameters. 
%   The input arguments can be given in any order with one exception:
%   PARx and VALx must be given in pairs in that order.
%   Any unspecified parameters for PARx are set to []. 
%   Parameters with value [] indicate to use the default value for that
%   parameter when OPTIONS is passed to the function. It is sufficient to
%   type only the first characters that uniquely identify the parameter
%   or function name.  Upper case letters for parameter names and values
%   that are strings are ignored. If an invalid string is provided, the
%   default is used.
%   
%   PARSEOPTIONS with no input and no output arguments displays this help 
%
% Examples:
%  defaultoptions = struct('test',[],'integrate',[] );
%  options = parseoptions(defaultoptions,'int','yes');
%  assert(options, struct('test',[], 'integrate', 'yes'))
%  opt = defaultoptions;
%  opt.test = 'yes';  
%  options2 = parseoptions(defaultoptions,'int','yes',opt);
%  assert(options2, struct('test','yes', 'integrate', 'yes'))
  
% History
% Revised pab april 2007
% -added possibility to have structarray, ie, numel(options)>1
% Revised pab Feb 2007
% - added possibility to parse custom classes different from a struct if this 
%   function resides in the custom class directory.
% by pab 14.05.2003
% based on MATLAB's optimset

% if (nargin == 0) 
%   % Display help
%   help('parseoptions')
%   return;
% end


defaultclass = 'struct';
% Initialization
nonValidClasses = {'double','single','logical','char','cell','function_handle','int8',...
      'uint8','int16','uint16','int32','uint32','int64','uint64'};    




fnames = '';

% Legal functions names
ix = 1;
if ((ix<=nargin) && (ischar(varargin{ix}))),
  fnames = varargin{ix};
  
  ix = 2;
end

try
  className = class(varargin{ix});
  if ismember(className,nonValidClasses) %(~isa(varargin{ix},defaultclass))),
      error('WAFO:PARSEOPTIONS','%s expected.',defaultclass)
  else
    defaultclass = className;
  end
  options = varargin{ix};
  ix = ix + 1;
catch
  error('WAFO:PARSEOPTIONS','%s expected.',defaultclass)
end

% Legal parameter names
namesc = fieldnames(options);
names = lower(namesc);
%names  = lower(strvcat(namesc{:}));

expectval = 0;         % start expecting a name or stucture, not a value
while ix <= nargin
  arg = varargin{ix};
  if expectval
    [options.(namesc{iy})]=deal(arg);
    expectval = 0;
  else
    switch lower(class(arg))
      case 'char',
        lowArg = strtok(lower(arg),'.');      % Downcase and strip .m extensions
        iy = findlname(lowArg,fnames,0);      % Find legal function name
        if length(iy)==1,
          opt1    = defaultopt(fnames(iy,:)); % Get default options for function
          options = combineopt(options,opt1); % Combine with old options
        else
          iy = findlname(lowArg,names,1);     % Find legal parameter name
          if isempty(iy),                     % Illegal parameter name
            ix = ix+1;                        % Skip next input
          else
            expectval = 1;                    % expect a value next
          end
        end
      case {'struct',defaultclass}
        options = combineopt(options,arg);
      otherwise
        error('WAFO:PARSEOPTIONS','Expected argument %d to be a string or a struct.', ix);
    end
  end
  ix = ix+1;
end

if expectval
   error('WAFO:PARSEOPTIONS','Expected value for parameter ''%s''.', arg);
end
return

function iy = findlname(arg,names,chkempty)
% FINDLNAME find index to legal parameter name  

iy = strmatch(arg,names);

Niy = length(iy);
if Niy<1 && chkempty                      % No matches
  msg = ['Unrecognized parameter name ' arg '.'];
  warning('WAFO:PARSEOPTIONS',msg);
elseif Niy > 1      % More than one match
  % Find exact matches in case any names are subsets of others
  k = strmatch(arg,names,'exact');
  if length(k) == 1
    iy = k;
  else
    if iscell(names)
      ambiguosNames = names(k);
    else
      ambiguosNames = cellstr(names(k,:));
    end
    ambTxt = sprintf('%s, ',ambiguosNames{:});
    ambTxt(end-1:end) = '';
    error('Ambiguous parameter name %s (%s)',arg,ambTxt);
  end
end
return

function options = combineopt(options,newopts)
%COMBINEOPT combines an existing options structure
%   OPTIONS with a new options structure NEWOPTS.  Any parameters in NEWOPTS
%   with non-empty values overwrite the corresponding old parameters in
%   OPTIONS. 

if ~isempty(newopts)                      % [] is a valid options argument
  newNames    = fieldnames(newopts);
  legalNames  = fieldnames(options);
  
  [commonNames, ia, ib] = intersect(lower(newNames),lower(legalNames));
  if isempty(ia),
    return;
  end
  cval = struct2cell(newopts);
  
  Numold = numel(options);
  Numnew = numel(newopts);
  if Numold>1
    % multidimensional struct
    if Numnew==1
      ind  = findNonEmptyCells(cval(ia));
      if any(ind)
        val  = struct2cell(options);
        val(ib(ind),:) = cval(ia(ind),ones(1,Numold));
        options = cell2struct(val,legalNames,1);
      end
    else
      if Numold~=Numnew
        error('WAFO:parseoptions','The new struct must match the number of elements in the old one!')
      end
      d = Numold;
      facta =  [0 length(newNames)*(1:d-1)];
      ia =  ia(:,ones(1,d))+facta(ones(length(ia),1),:);
      
      factb =  [0 length(legalNames)*(1:d-1)];
      ib =  ib(:,ones(1,d))+factb(ones(length(ib),1),:);
      ind  = findNonEmptyCells(cval(ia));
      if any(ind)
        val  = struct2cell(options);
        val(ib(ind)) = cval(ia(ind));
        options = cell2struct(val,legalNames,1);
      end
      
    end
  else
    if Numnew>1
      warning('WAFO:parseoptions','Multidimensional new options. Ignoring all elements but the first!')
    end
    ind  = findNonEmptyCells(cval(ia));
    if any(ind)
      val  = struct2cell(options);
      val(ib(ind)) = cval(ia(ind));
      options = cell2struct(val,legalNames,1);
    end
  end
end
return


function options = defaultopt(fname)
% DEFAULTOPT Get default OPTIONS structure from the function FNAME

fname = deblank(fname);

if ~exist(fname,'file')
  msg = ['No default options available: the function ' fname ' does not exist on the path.'];
  error(msg)
end
try 
  options = feval(fname,'defaults');
catch
  msg = ['No default options available for the function ' fname '.'];
  error(msg)
end
return

function ind = findNonEmptyCells(carray)
  try % matlab 5.3 or higher
    ind = find(~cellfun('isempty',carray)).';
  catch
    % Slow 
    n = length(carray);
    ind1 = zeros(1,n);
    for ix = 1:n
      ind1(ix) = isempty(carray{ix});
    end
    ind = find(~ind1);
  end
  return