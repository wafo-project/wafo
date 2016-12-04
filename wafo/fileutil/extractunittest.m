function str1 = extractunittest(mfile,method)
%extractunittest 
%
% CALL str1 = extractunittest(mfile,method)
%
% mfile = name of file
% method = 'example',
%          'unittest',
%          'all'
%
  str1 = cell(1,0);
  LF = char(10);

  if strncmpi(method,'all',1) || strncmpi(method,'example',1),
    str1 = extract_doctest(mfile);
  end
  if strncmpi(method,'all',1) || strncmpi(method,'unittest',1)  
    string = freadtxt([strtok(mfile,'.'),'.m']);
    LF  = char(10); % Linefeed character
    space = ' ';
    string = [space string(:)' LF space];  % Add blank edges
    string = dm2unix(string,space);     % convert to unix line endings.
    

    % Quotes and comments masks ...............
    [mask_q,mask_c,i_LF,lineNumber] = findquot(string);
    no = findstr(string,'% UNITTEST %');
    no = intersect(no,[1,i_LF+1]); % Dismiss all UNITTESTS not starting at beginning of line
    if any(no)
      no = i_LF(lineNumber(no))+1; % Index to the next line
      if mask_c(no)
        nc = i_LF(lineNumber(no));
        ic = find(mask_c(nc)& mask_c(nc+1));
        while any(ic)
          nc(ic) = i_LF(lineNumber(nc(ic))+1);
          ic = find(mask_c(nc) & mask_c(nc+1));
        end
        nco = findstr(string,'%');
        nco = intersect(nco,[1,i_LF+1]); % index to all comment charancters starting at beginning of line
        string(nco) = space; % replace with space
        for ix = 1:numel(no)
          str1{end+1} = string(no(ix),nc(ix));
        end
      end
    end
  end
end %% extract unittest

function str1 = extract_doctest(mfile)
  str1 = cell(1,0);
  LF = char(10);

  try  
     str = dm2unix(_help(mfile),' ');          % Extract help header
  
    if ~isempty(str)
      inl = findnl(str);             % Find newlines
      [names,p0 p1] = getnames(str); % Extract legal names
      names2 = lower(names);
      ind1 = max(strmatch('example',names2)); % Find token string
      if isempty(ind1)
        ind1 = max(strmatch('examples',names2));
      end
      if ~isempty(ind1) && any(str(p0(ind1)+(-4:-1) )==LF)
        % Find the first line after the token
        ind0 =  inl(find(inl-p0(ind1)>0,1,'first'))+1;
        if ~isempty(ind0)
          ind2 = strmatch('see',names2(ind1:end,:),'exact')+ind1-1;
          try
            ind2 = max(ind2(strmatch('also',names2(ind2+1,:),'exact')));
          catch
            ind2 = [];
          end
          if isempty(ind2),
             teststr = deblank(str(ind0:end));
          else
            ind1 = p0(ind2)-1;
            teststr= deblank(str(ind0:ind1));
          end
          if ~isempty(teststr)
            str1{end+1} = teststr;
          end
        end
      end
    end
  end
end

function s = _help(file)
  try
    s = help(file);
   catch
    [d, name, e] = fileparts(file);
    s = help(name);
  end
end
