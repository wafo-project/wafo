function [C,D,M,G2] = finddep(ldirs,edirs)
%FINDDEP Find dependencies and crossreferences to other functions 
% 
% CALL: [C,W,M,G] = finddep(ldirs,edirs);
% 
% C     = {C1 C2} cell array of cross reference matrices (size [n1 n2])
%         C1 : dependency matrix
%         C2 : cross references found in help header.
% W     = A structure array with the fields:
%         path : full path to directories given in ldirs
%         m    : cell array of m-file names in the ldirs directories.
% M     = character array of file names in directories ldirs and
%         edirs. (length n2)
% G     = character array of file names in directories ldirs (length n1).
% ldirs = character array of directories    (size [k1  m1])
%         which are fully checked for cross references.
%         (default current directory)
% edirs = character array of directories    (size [k2  m2]) 
%         which are only checked for files called 
%         by files in ldirs 
%
% Note: ldirs and edirs may contain partial pathnames.
% 
% FINDDEP finds the dependency  matrix between all the
% m-files in directories ldirs and edirs, i.e., which m-files 
% call one another and which function they refer to. 
% If only ldirs is given then G=M.
%
% Example
%    [C,W,M] = finddep(fullfile(waforoot, 'fileutil'));
%    assert(W.m(1:4), {'Contents.m','bindiff.m','cdtomfile.m','diary2m.m'});
% %   assert(size(C{1,1}), [52,52]);
%    assert(C{1,1}(1:4,1:4), sparse([2,3,4],[2,3,4], ones(1, 3)));
%
% See also: getnames, strmatch 

% History
% revised pab 21.06.2001
% -aaargh found a logical error
% revised pab 15.05.2001
% - simplified for loop.
% revised by pab 16.10.2000


error(nargchk(0,2,nargin))

ignoreU=1; % =1 if ignoring upper case letters
%strcmpi(computer,'pcwin'),

if nargin<1||isempty(ldirs)
  ldirs = pwd;
end
%ldirs
[k1,m]=size(ldirs);
G2=[];
for ix = 1:k1
  t = what(deblank(ldirs(ix,:)));
  D(ix) = t(1);
  G2 = strvcat(G2,D(ix).m{:});
end 

M = G2;
if nargin ==2&&~isempty(edirs),
  [k2,m]=size(edirs);
  for ix = 1:k2
    DD(ix) = what(deblank(edirs(ix,:)));
    M = strvcat(M,DD.m{:});
  end   
end;
if ignoreU
  G2=lower(G2);
  M=lower(M);
end
[n1,m1]=size(G2);
[n2,m2]=size(M);

%C = zeros(n1,n2);
C1 = sparse(n1,n2);
C2 = C1;
iy=0;
for ix=1:k1
  for j =1:length(D(ix).m),
    iy=iy+1;
    filename = fullfile(D(ix).path,D(ix).m{j});
    disp(['Checking file ',filename]);
   
    %HH = help(filename); % help header
    s = freadtxt(filename);
    if ignoreU,
      s  = lower(s);
     % HH = lower(HH);
    end
    if 1, % refer package
	  % Quotes and comments  masks...............
      [mask_q,mask_c] = findquot(s);
      space=' ';
      
      % Blank the openings of quotes (strings) ........
      %nn = find(diff(mask_q)==1)+1;
      %s(nn) = space(ones(size(nn)));
      
      % Blank the quotes
      nn    = find(mask_q);
      s(nn) = space(ones(size(nn)));
      
      %s1=s;
      % Find help header
      nn = find(mask_c);
      if any(nn)
        nn2 = find(diff(nn)>1);
        if isempty(nn2), nn2 = length(nn);end
        HH = s(nn(1):nn(nn2));
  
        % Blank the comments
        s(nn) = space(ones(size(nn)));
      else
        HH = '';
      end
      
      s = char(s);
    end
    % If nothing left, jump to next ......................
    % extract only the legal names (and remove duplicates)
    s = munique(getnames(s)); % More robust and faster
    HH = munique(getnames(HH)); % More robust and faster
    if any(s>32) || any(HH>32),
      s = cellstr(s);
      HH = cellstr(HH);
      for k =1:n2
        word = strtok(M(k,:),'.'); % remove .m-ending
        % See if in the list
        res = strmatch(word,s,'exact');
        if ~isempty(res),
          C1(iy,k) = 1;
        end;
        res = strmatch(word,HH,'exact');
        if ~isempty(res),
          C2(iy,k) = 1;
        end;
      end;
    end;
  end;
end;
C={C1 C2};

