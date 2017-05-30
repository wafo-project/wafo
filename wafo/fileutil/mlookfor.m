function mlookfor( varargin )
%MLOOKFOR Search the help text of all M-files.
%
%   MLOOKFOR SEARCH-EXPRESSION or MLOOKFOR( 'SEARCH-EXPRESSION' )
%   searches all m-files along the Matlab search path including
%   private subdirectories and class subdirectories.  If the help text
%   in an m-file matches SEARCH-EXPRESSION, the path of the m-file is
%   printed.  The whole help text is searched, not only the first help
%   line.  Contents.m and Readme.m files are not searched.
%
%   SEARCH-EXPRESSION might be any combination of globs, regular
%   expressions, logical operators and parentheses.
%
%     * Globs are given as 'plain' strings.  The metacharacters '*' and
%       '?' are supported, as are character classes.
%
%     * Regular expressions must be surrounded by slashes: /regex/.  The
%       regular expression syntax is that of Perl.  To include a slash
%       inside a regular expression, escape it with a backslash.
%
%     * Logical operators that are allowed are (ordered by increasing
%       precedence) 'or', 'and' and 'not'.
%
%      Note that
%
%         mlookfor inverse and pseudo or generalized
%
%      would not be the same, since it would also find m-files that
%      contain the word 'generalized' but possibly without 'inverse',
%      since 'or' has lower presedence than 'and'.  Also not that
%
%         mlookfor ( pseudo or generalized ) and inverse
%
%      would cause an error and would have to be rewritten as
%
%         mlookfor( '( pseudo or generalized ) and inverse' )
%
%   Some details
%   ------------
%   Default is to do a case-insensitive substring search.  This can
%   only be changed by using the appropriate regular expression:  To
%   do a case-sensitive search, prepend '(?-i)' to the regular
%   expression.  To match words, use word-boundary anchors '\b'.
%   Since a substring search is done, any leading or trailing glob or
%   regular expression that can match nothing is redundant.
%
%   Regular expression syntax
%   -------------------------
%   See the perlre(5) manual page for more details about the Perl
%   regular expression syntax.  Typing 'perldoc perlre' at a shell
%   command line should give you this page.

%
%   Examples
%   --------
% %   1) To find all m-files that contain 'pseudo' and 'inverse':
%
%         mlookfor pseudo and inverse
%
% %  2) To find all m-files that contain 'pseudo' or 'generalized' and
% %     'inverse':
%
%         mlookfor inverse and ( pseudo or generalized )
%
%
% %  3) To find all m-files that contain 'pseudoinverse' or
% %     'pseudo-inverse', use
%
%         mlookfor pseudoinverse or pseudo-inverse
%
% %     or a regular expression where the '-' is optional, like
%
%         mlookfor /pseudo-?inverse/
%
%
%   See also LOOKFOR.


%   Installing
%   ----------
%   Put the Matlab .m file and the Perl .pl file in the same directory.
%   Make sure the directory is in Matlab's search path.
%
%   UNIX:     Make the Perl file executable, using something like
%             'chmod u+x mlookfor.pl'.  Make sure that the first line
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

%   Author:      Peter J. Acklam
%   Time-stamp:  1999-11-15 16:12:28
%   E-mail:      jacklam@math.uio.no
%   WWW URL:     http://www.math.uio.no/~jacklam

% Check arguments.
if isempty(varargin)
   error('No seach expression specified.');
end

% Get the name of the perl-file and see if it exists.
fullname = which( mfilename );
[ dirpath, basename, suffix ] = fileparts(fullname);
perlfile = fullfile(dirpath, [ basename '.pl' ]);
if ~exist( perlfile, 'file' )
   error('Perl file "%s" not found.', perlfile);
end

% Write a temporary file that will be read by perl.
tmpfile = tempname;
fid = fopen( tmpfile, 'w' );
if fid < 0
   error('Can''t open file "%s" for writing.', tmpfile);
end
fprintf(fid, '%s\n', matlabpath, varargin{:});
fclose(fid);

% Let perl do the actual searching.
cmp = computer;
if isunix
   %
   % UNIX section.
   %
   unix([ perlfile ' ' tmpfile ]);
elseif strcmp(cmp(1:2), 'PC')
   %
   % PC section.
   %
   dos([fullfile(perlexepath ,'perl.exe') ' ' perlfile ' ' tmpfile ], '-echo');
elseif strcmp(cmp(1:2), 'MA')
   %
   % MAC section.
   %
   error('MLOOKFOR is not supported on this platform.');
else
   error('Unknown platform.');
end

