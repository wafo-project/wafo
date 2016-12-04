function mlstrip(infile,outfile)
%MLSTRIP  Strip comments from a matlab m-file.
%
% CALL  mlstrip(original.m ,stripped.m)
%
%   See also LOOKFOR.


%   Installing
%   ----------
%   Put the Matlab .m file and the Perl .pl file in the same directory.
%   Make sure the directory is in Matlab's search path.
%
%   UNIX:     Make the Perl file executable, using something like
%             'chmod u+x *.pl'.  Make sure that the first line
%             (the "she-bang" line) of the Perl file contains the
%             correct path to perl.
%
%   Windows:  Put the Matlab .m file and the Perl .pl file in the same
%             directory.  Make sure the directory is in Matlab's
%             search path.  Edit the Matlab file (look for 'PC
%             section') and make sure that the path to the perl
%             executable is correct.
%
%   Interaction between Matlab and Perl
%   -----------------------------------
%   The user calls the Matlab file and gives the search expression as
%   input.  Matlab then creates a temporary file containing the Matlab
%   search path on the first line and the search expression on the
%   following lines.  Matlab then calls perl to execute the Perl file.
%   The only input argument to the Perl file is the name of the
%   temporary file.  Perl reads the information in this temporary file
%   and deletes the file immediately.  Then perl performs the actual
%   searching and displays the result in the Matlab command window.
%
%   How the matching is done
%   ------------------------
%   First the help text is extracted from the m-file.  All comment
%   markers ('%') are removed.  Then each sequence of whitespace
%   (tabs, spaces, newlines, carriage returns, etc.) is converted to a
%   single space character.  Finally the search expression is matched
%   against the help text.  For more information, see the Perl file.
%
%   Legal stuff
%   -----------
%   This program is free software; you can redistribute it and/or
%   modify it under the terms of the GNU General Public License as
%   published by the Free Software Foundation; either version 2 of
%   the License, or (at your option) any later version.
%
%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
%
%   The GNU General Public License can be found on the World Wide
%   Web at http://www.gnu.org/copyleft/gpl.html and it can be
%   obtained by writing to the Free Software Foundation, Inc.,
%   59 Temple Place - Suite 330,  Boston, MA 02111-1307, USA.

% History
% by pab
% based on MLOOKFOR.M
% Check arguments.
error(nargchk(2,2,nargin))

% Get the name of the perl-file and see if it exists.
fullname = which( mfilename );
[ dirpath, basename, suffix ] = fileparts(fullname);
perlfile = fullfile(dirpath, [ basename '.pl' ]);
if ~exist( perlfile, 'file' )
   error('Perl file "%s" not found.', perlfile);
end

% Write a temporary file that will be read by perl.
%tmpfile = tempname;
%fid = fopen( tmpfile, 'w' );
%if fid < 0
%   error(sprintf('Can''t open file "%s" for writing.', tmpfile));
%end
%fprintf(fid, '%s\n', matlabpath, varargin{:});
%fclose(fid);

% Let perl do the actual searching.
cmp = computer;
if isunix
   %
   % UNIX section.
   %
   unix([ perlfile ' ' infile ' > ' outfile ]);

elseif strcmp(cmp(1:2), 'PC')
   %
   % PC section.
   %
   dos([ perlexepath 'perl.exe ' perlfile ' ' infile ' > ' outfile ], '-echo');
elseif strcmp(cmp(1:2), 'MA')
   %
   % MAC section.
   %
   error('MLSTRIP is not supported on this platform.');
else
   error('Unknown platform.');
end

