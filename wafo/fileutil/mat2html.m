function mat2html(varargin)
%MAT2HTML  Produce html files out of Matlab m-file.
%  Usage:
%
%  mat2html [-i] [-q] [-g] [-M matlab_dir] [-H html_dir]
%  mat2html [-i] [-q] [-g] [-p] [-H html_dir]
%
%mat2html reads a list of matlab .m files and/or directories from the
%standard input to produce a hypertext documentation suitable for
%browsing with a WWW browser such as Mosaic or Netscape.  An file
%html_dir/index.html in written. Subdirectories are written for matlab
%directory that is encountered containing
%.html files corresponding to every .m file.
%
%Help:
%
%  -h   Print this help message and exit.
%
%  -q   Quiet operation. Default is verbose.
%
%Output Options:
%
%  -H  Specify a directory to place the html files and subdirectories.
%      The default is the current directory. If necessary, an attempt
%      is made to create the directory. The file index.html is placed
%      in this directory.
%
%  -i  Include matlab source code in the html documentation
%
%  -g  Do global hypertext links, that is, hypertext links among
%      separate matlab directories. Default is to do hypertext
%      links only among functions in the same directory.
%
%Matlab Source Options:
%
%      The standard input is looked at first. If there is nothing there,
%      then we look in the current directory for .m files as if we did a
%      -M . option.
%
%  -M  Specify a root matlab directory to search. The standard input
%      is not read.
%
%  -p  Search the matlab path obtained by matlab -n. Options -M
%      and -p are incompatable.
%
%Typical usages are:
%
%  mat2html 
%      Produces a file index.html in the current directory, and a 
%      .html file corresponding to each .m file in the current directory.
%      An index.html file is produced in the current directory.
%
%  mat2html -M matlab -H html
%      Same as above, except it looks in subdirectory matlab of
%      the current directory for the .m files, and puts the results
%      in subdirectory html. Useful for separating the .m files
%      from the documentation.
%
%  mat2html -p -H html
%      Creates documentation for everything on the default matlab
%      path, and puts it into a directory html.
%
%  ls *.m | mat2html -H html
%      Index the .m files in the current directory
%
%  find . -name "*.m" -print | mat2html -H html -i
%      The find command recursively builds a list of all .m files
%      to be found the directory and its subdirectories. These
%      are then processed and the html files put in the directory
%      html. The matlab source code is included in the html files.
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
%% Check arguments.
%if length(varargin) == 0
%   error('No search expression specified.');
%end
% Get the name of the perl-file and see if it exists.
fullname = which( mfilename );
[ dirpath, basename, suffix ] = fileparts(fullname);
perlfile = fullfile(dirpath, [ basename '.pl' ]);
if ~exist( perlfile, 'file' )
   error(sprintf('Perl file "%s" not found.', perlfile));
end

args = varargin;
if ~isempty(args)
  [args{2,:}]=deal(' ');
end
% Write a temporary file that will be read by perl.
tmpfile = tempname;
fid = fopen( tmpfile, 'w' );
if fid < 0
   error('Can''t open file "%s" for writing.', tmpfile);
end
fprintf(fid, '%s ', varargin{:});
fclose(fid);
%tmpfile
% Let perl do the actual searching.
cmp = computer;
if isunix
   %
   % UNIX section.
   %
   unix([ perlfile ' ' tmpfile]);

elseif strcmp(cmp(1:2), 'PC')
   %
   % PC section.
   %
%   eval(['!' perlexepath 'perl.exe ' perlfile ' ' args{:}]);
   dos([ perlexepath 'perl.exe ' perlfile ' ' tmpfile ], '-echo');
elseif strcmp(cmp(1:2), 'MA')
   %
   % MAC section.
   %
   error('MAT2HTML is not supported on this platform.');
else
   error('Unknown platform.');
end

