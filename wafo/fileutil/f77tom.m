function f77tom( varargin )
%F77TOM Convert one or more F77 or F90 files to Matlab M-files.  
% 
% CALL f77tom(f1,f2,...)
%
% WARNING: This program is not up to version 1.0 as of July 1998.
%	There are a number of shortcomings and it is not fully tested.
%	RSVP to bob@math.umn.edu with all bug reports.
%	However, don't tell me about how it will not handle GO TOs in
% your Fortran code.  :^)
%	Note also the location of your perl binary may differ from what
% is on line number 1 here.  Adjust as necessary.
%	All development and testing done on a Pentium II with Linux 2.0.34.
% Your mileage may vary.
%
%   See also m2c.

%	If you have any questions, please contact Chris Cornuelle at
% bob@math.umn.edu, Minnesota Center for Industrial Mathematics
% at the University of Minnesota.  
%	Dr. C. Cornuelle
%	School of Mathematics/MCIM
%	206 Church Street SE
%	Minneapolis, MN 55455 USA



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

%history
% by pab 30.01.2001

% Check arguments.
if isempty(varargin)
   error('No search expression specified.');
end
mfilename = 'f77tom.pl';
% Get the name of the perl-file and see if it exists.
fullname = which( mfilename );
[ dirpath, basename, suffix ] = fileparts(fullname);
perlfile = fullfile(dirpath, [ basename '.pl' ]);
if ~exist( perlfile, 'file' )
   error('Perl file "%s" not found.', perlfile);
end

tmpfiles = varargin;
[tmpfiles{2,:}]=deal(' ');

% Let perl do the actual searching.
cmp = computer;
if isunix
   %
   % UNIX section.
   %
   unix([ perlfile ' ' tmpfiles{:} ]);
elseif strcmp(cmp(1:2), 'PC')
   %
   % PC section.
   %
   dos([ perlexepath 'perl.exe ' perlfile ' ' tmpfiles{:} ], '-echo');
elseif strcmp(cmp(1:2), 'MA')
   %
   % MAC section.
   %
   error('F77TOM is not supported on this platform.');
else
   error('Unknown platform.');
end

