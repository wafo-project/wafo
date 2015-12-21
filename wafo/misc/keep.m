function keep(varargin)
%KEEP Keeps workspace variables of your choice and clears the rest.
%	
%     CALL: keep var1 var2 var3 ...
% 
% KEEP keeps the workspace variables of your choice and clears the rest.
% It is like the CLEAR command except it only works for variables. 
% The wildcard character '*' can be used to keep variables that match a pattern.
% For instance, KEEP X* keeps all the variables in the current workspace
% that start with X.
%	
% Examples:
% [set1data, set2data, set3data, data,test,test2,test3] = deal(1);
% who                   % list your workspace variables
% keep *data
% who
% keep set1data set2data set3data
% who
%
% See also clear


% Revised pab 18 May 2005  
% added possibility to allow wildcards
%revised pab 12.11.2002
% fixed a bug.
%revised pab 21.06.2001
% -vectorized the code.
% Yoram Tal 5/7/98    yoramtal@internet-zahav.net
% MATLAB version 5.2
% Based on a program by Xiaoning (David) Yang


  
% Keep all
if isempty(varargin) || any(strcmp('*',varargin))
  return
end

% Find variables in workspace
workspace = 'caller';  % 'caller' or 'base'
vnames = evalin(workspace,'who');
vnames = vnames(:).';
% Remove variables in the "keep" list
del = getdeletelist(vnames,varargin{:});
if ~isempty(del)
  % Clear selected variables
  evalin(workspace,['clear ' del])
end
return
function deletelist = getdeletelist(vnames,varargin)
  
  index = zeros(1,length(vnames));
    
  % Create string of variables to keep.
  k = any(strvcat(varargin{:})=='*',2);

  k0 = find(~k(:).');
  if any(k0)
    [tmp,k10]  = intersect(vnames,varargin(k0));
    index(k10) = 1;
  end

  vnames1      = vnames;
  vnames1(2,:) = {'  '};
  deletelist   = [' ' vnames1{:,find(~index)}];

  k1 = find(k(:).');
  if any(k1)
    v   = version;
    ix  = find(v=='.');
    vnr = str2double(v(1:ix(min(2,length(ix)))-1));

    if (vnr>6.5)          
      for i = k1
        name    = ['\s' strrep(varargin{i},'*','\w*') '\s'];
        [lo,hi] = regexp(deletelist, name);
        idx     = idxlimits2idxvector(lo,hi);
        if any(idx)
          deletelist(idx) = ' ';
        end
      end
      
    else
      
      for i = 1:length(vnames)
        vnames{i} = [' ' vnames{i} ' '];
      end
      for i = k1
        name    = [' ' varargin{i} ' '];
        k = find(strwcmp(vnames,name));
        if any(k)
          index(k) = 1;
        end
      end
      deletelist = [' ' vnames1{:,find(~index)}];
    end
  end
  if all(isspace(deletelist))
    deletelist = '';
  end
  return
  
  
  function idx = idxlimits2idxvector(lo,hi)
    % IDXLIMITS2IDXVECTOR Convert from index limits to index vector
    %
    % Given the vectors LO and HI, this function creates the index vector
    %  idx = [lo(1):hi(1),lo(2):hi(2),.....]
    
    
    i  = find(lo > hi);
    if any(i)
      %Delete invalid entries.
      lo(i) = [];
      hi(i) = [];
    end
    
    m   = length(lo); %
    len = hi-lo+1;
    n   = sum(len);
    
    idx = ones(1,n);
    if n<1
      return
    end
    idx(1) = lo(1);
    len(1) = len(1) +1;
    idx(cumsum(len(1:m-1))) = lo(2:m)-hi(1:m-1);
    idx = cumsum(idx);
    
