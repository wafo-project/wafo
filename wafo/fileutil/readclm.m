function  [outdata,head] = readclm(filename,nclm,skip,formt)
% READCLM Reads numerical data from a text file into a matrix.
%	Text file can begin with a header or comment block.
%	[DATA,HEAD] = READCLM(FILENAME,NCLM,SKIP,FORMAT)
%	Opens file FILENAME, skips first several lines specified
%	by SKIP number or beginning with comment '%'.
%	Then reads next several lines into a string matrix HEAD
%	until the first line with numerical data is encountered
%	(that is until first non-empty output of SSCANF).
%	Then reads the rest of the file into a numerical matrix
%	DATA in a format FORMAT with number of columns equal
%	to number of columns of the text file or specified by
%	number NCLM. If data does not match the size of the
%	matrix DATA, it is padded with NaN at the end.
%
%	READCLM(FILENAME) reads data from a text file FILENAME,
%	skipping only commented lines. It determines number of
%	columns by the length of the first data line and uses
%	the floating point format '%g';
%
%	READCLM uses FGETS to read the first lines and 	FSCANF
%	for reading data.

%  Kirill K. Pankratov,  kirill@plume.mit.edu
%  03/12/94, 01/10/95.
% revised pab 20.07.2002
%  -removed  ==[] statements with isempty


 % Defaults and  parameters ..............................
formt_dflt = '%g';  % Default format for fscanf
addn = NaN;         % Number to fill the end if necessary

 % Handle input ..........................................
if nargin<1, error('  File name is undefined'); end

if nargin<4||isempty(formt), formt = formt_dflt; end
if nargin<3||isempty(skip),  skip  = 0; end
if nargin<2||isempty(nclm),  nclm  = 0; end
 

 % Open file ............................
[fid,msg] = fopen(filename);
if fid<0, disp(msg), return, end

 % Find header and first  data line ......................
is_head = 1;
jl      = 0;
head    = ' ';

while (is_head), % Add lines to header.....
  s  = fgets(fid);           % Get next line
  jl = jl+1;
  %is_skip = (jl<=skip);
  is_skip = (jl<=skip) | (s(1)=='%');

  out1 = sscanf(s,formt);   % Try to read this line

   % If unreadable by SSCANF or skip, add to header
  is_head = (isempty(out1) | is_skip);

  if (is_head && ~is_skip)
    head = str2mat(head,s(1:length(s)-1));
  end
end
head = head(2:size(head,1),:);

 % Determine number of columns if not specified
out1 = out1(:)';
l1   = length(out1);

if (~nclm), nclm = l1; end

 % Read the rest of the file ..............................
if l1~=nclm  % First line format is different from ncolumns

  outdata = fscanf(fid,formt);
  lout    = length(outdata)+l1;
  ncu     = ceil(lout/nclm);
  lz      = nclm*ncu-lout;
  outdata = [out1'; outdata(:); ones(lz,1)*addn];
  outdata = reshape(outdata,nclm,ncu)';

else              % Regular case

  outdata = fscanf(fid,formt,[nclm inf]);
  outdata = [out1; outdata'];  % Add the first line

end

fclose (fid);     % Close file ..........

