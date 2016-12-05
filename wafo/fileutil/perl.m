function result = perl(varargin)
%PERL   Call Perl to execute input arguments.
%
%   PERL ARG1 ARG2 ... is like saying !PERL ARG1 ARG2 ... except that the
%   former takes special action on Windows systems to avoid using the
%   stripped-down Perl executable that is bundled with Matlab.
%
%   Note that certain arguments must be put inside single quotes to prevent
%   the string from being interpreted by matlab.
%
%   Examples:
%
%     perl -v                           % print version information
%
%     perl -le "print for @INC"         % print module search path
%     perl '-le' '"print for @INC"'     % equivalent
%     perl('-le', '"print for @INC"')   % ditto
%
%     perl -le '"print for split /;/, $ENV{PATH}"'      % print DOS path

%   Author:      Peter J. Acklam
%   Revised by   Per A. Brodtkorb 
%   Time-stamp:  2000-07-30 00:14:06
%   E-mail:      jacklam@math.uio.no
%   WWW URL:     http://www.math.uio.no/~jacklam

  
% build string with perl arguments
  
  perlargs = cell(2,nargin);
  perlargs(2,:) = {' '};
   
  for i=1:nargin
    thisArg = varargin{i};
    if (isempty(thisArg) || ~ischar(thisArg))
      error('All input arguments must be valid strings.');
    elseif (exist(thisArg, 'file')==2)
      % This is a valid file on the MATLAB path
      if isempty(dir(thisArg))
        % Not complete file specification
        % - file is not in current directory
        % - OR filename specified without extension
        % ==> get full file path
        thisArg = which(thisArg);
      end
    end
  
    % Wrap thisArg in double quotes if:
    % 1) it contains spaces
    % 2) and is not already put in double qoutes
    if (any(thisArg == ' ') && ((thisArg(1)~='"') || (thisArg(end)~='"')))
      thisArg = sprintf('"%s"', thisArg);
    end
    perlargs{1,i} = thisArg;
  end	
   
   perlArgString = [perlargs{1:end-1}];
      
   if isempty(perlArgString)
     error('No perl command(s) specified!')
   elseif strncmpi(computer,'PC',2), 
     % PC  
     perlCmd = fullfile(perlexepath,'perl.exe');
     cmdString = sprintf('"%s" %s',perlCmd,perlArgString);
     [status, result] = dos(cmdString);      
   else
     % unix
     [status, perlCmd] = unix('which perl');
     if (status == 0)
       cmdString = sprintf('perl %s', perlArgString);
       [status, result] = unix(cmdString);
     else
       error('Perl binary not found. Is Perl installed properly?');
     end
   end

   % Check for errors in shell command
   if (status~=0)
     errorTxt = sprintf('System error: %s \nCommand executed: %s\n',...
     result,  cmdString);
     error(errorTxt);
   end
   return