function  [mask_q,mask_c,inl,linenum] = findquot(string)
%FINDQUOT  Find quote and comment mask as well as new lines and line numbers.
%
% CALL: [maskq, maskc,inl,linenum] = findquot(string)
%
%  maskq, maskc = quote and comment mask, respectively. The masks
%                 includes  the quotes and the comment characters.
%  inl          = indices to new line characters of string.
%  linenum      = line number mask, (i.e., vector of line numbers).
%  string       = character vector to be searched.
%
% Note: MASKC includes the end-of-line characters.  
%  
% Example:
%  t = freadtxt('findquot.m');
%  [mq, mc] = findquot(t); 
%  assert(t(mq)(2:11), 'findquot.m')
%  assert(t(mc)(1:41), '%FINDQUOT  Find quote and comment mask as')
%
% See also: findnl 

% Tested on: Matlab 5.x
% History:
% revised 30.10.2003
%  added '{' as a legal character to precede quote character  
% revised pab 10.10.2003
% fixed a bug: added all white space characters to ch_str  
% revised pab 06.12.2000
%  - Fixed a bug: Now it does tackle quotes inside comment lines properly!
%  - made sure maskc and maskq is logical
%  - maskc now also contain the end of line character
% revised pab 10.11.2000
%  - Fixed a bug: forgot to remove the added ch_nl entries from the 
%      inl and linenum output
% revised pab 05.10.2000
%  - updated documentation
%  - fixed some bugs: 
%     1) the beginning and ending quote character is now included in maskq 
%     2) Comment lines on the last line produced an internal error in the
%        line marking the end of comments lines:
%           mask_c(inl(ind)) = -ones(size(ind));
%        This is fixed by adding ch_nl to string before the calculation
%        and remove the two last entries of maskq and maskc before
%        returning the result.
%  - reordered the output arguments so that inl and linenum come last
% by  Kirill K. Pankratov, kirill@plume.mit.edu  12/26/94

%error(nargchk(1,1,nargin))
narginchk(1,1)
% Determine the new line character .........................
ch_nl            = [13 10]; % New line characters
quoteCharacter   = '''';    % char(39);


strsize = size(string);


 % Make string a row vector and make sure the last character is a Newline
string = string';
string = [string(:)' ch_nl] ;

 % New lines and line numbering .........
[inl, linenum] = findnl(string);


 % Find explicit strings (inside single quotes) .....
i_sq   = find(string == quoteCharacter);   % Pointers to single quotes
num_sq = length(i_sq);


mask_c = zeros(size(string));
mask_q = mask_c ;

if num_sq>0, 
  
  whiteSpaceChars  = char(find(isspace(char(1:32)))); 
  % white space characters are spaces, newlines,
  % carriage returns, tabs, vertical tabs, and formfeeds

  % Legal char. before beginning of string
  %ch_str = [' '':;,+-*\/=><([',  ch_nl]'; 
  ch_str  = [':;,+-*\/=><([{',quoteCharacter, whiteSpaceChars].';
  
  ol    = ones(1,num_sq);
  o_str = ones(size(ch_str));
 
    
  % Sort indices to single quotes and newline characters.
  a0  = sort([i_sq inl+.5]);
  % Set quotes to 1 and newline characters to 0.
  ind = (floor(a0)==a0);
  
  
  % Find single quotes preceded by legal character
  a00 = any(string(o_str,i_sq-1)==ch_str(:,ol));
  
  a0  = zeros(size(a0));
  
  % Keep only single quotes preceeded by legal character
  a0(ind) = a00;
  a00     = a0;
  a1      = a0;
  
  % Find beginning of quotes
  ll = length(a0);
  stop = 0;
  while ~stop
    a1(2:ll) = a00(2:ll)&~a0(1:ll-1);
    stop     = all(a1==a0);
    a0       = a1;
  end

  a0  = a1(ind);   % Keep only single quotes
  ind = find(a0);  % Indices to beginnings of single quotes.


       % To be a quote, 2 single quotes on the same line is 
       % required. Check if this is the case. If not then it 
       % is a single quote inside a comment and must 
       % therefore be removed from ind.
  if ~isempty(ind) && any(num_sq == ind(end)),
    ind(end)=[];   % pab 06.12.2000
  end
  if ~isempty(ind)
    % Keep only quotes on the same line
    ind = ind(linenum(i_sq(ind))==linenum(i_sq(ind+1))); % pab 06.12.2000
  end
  
  % Calculate mask for quotes
  if ~isempty(ind)
    ol         = ones(size(ind));
    a0         = i_sq(ind);     % Indices to the opening quote
    mask_q(a0) = ol;              % Mark beginning of quote
    a1         = i_sq(ind+1);   % Indices to the closing quote
    mask_q(a1) = -ol;             % Mark end of quote
    mask_q     = cumsum(mask_q);  % Mask from start quote to end quote
    mask_q(a1) = ol;              % include end quote as well pab 06.12.2000
  end
end

 % Find comments ...................................
fnd_cm = find(string=='%');         % Find comments char.

if ~isempty(fnd_cm),
  fnd_cm = fnd_cm(~mask_q(fnd_cm)); % Exclude those inside strings
end
if ~isempty(fnd_cm),
  % Find lines where comment char. is found
  ind = linenum(fnd_cm);

  a0 = zeros(1,length(inl)+1);

  % Find first comment char. on each line + indices to the end of line char.
  a0(fliplr(ind)) = fliplr(fnd_cm);
  ind = find(a0);
  a1  = a0(ind);   % Indices to the first comment char.
  a0  = inl(ind);  % Indices to end of line char.
  
  % Mask for comments ..................
  ol         = ones(size(a1)); 
  mask_c(a1) = ol;             % Mark the first comment char on each line
  mask_c(a0) = -ol;            % Mark the end of line .....
  mask_c     = cumsum(mask_c); % Mask from '%' marker till the end of line
  mask_c(a0) = ol;             % Include end of line character in the mask
  
  % Remove strings inside comments .....
  %mask_q = (mask_q&~mask_c);
end

% Remove the mask to added newline characters
mask_q(end-1:end) = [];
mask_c(end-1:end) = [];

% Also remove added newline characters from inl and linenum
inl(end)           = [];
linenum(end-1:end) = [];


nsiz    = strsize([ 2 1 3:end]);
mask_q  = logical(reshape(mask_q,nsiz).');
mask_c  = logical(reshape(mask_c,nsiz).');


if prod(nsiz)~=max(nsiz)
  linenum = reshape(linenum,nsiz).';
  [I, J] = ind2sub(nsiz,inl);
  inl    = sub2ind(strsize,J,I);
end
return








