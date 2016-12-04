#!/usr/local/bin/perl5 -w
#
# mlookfor - Search the help text of Matlab m-files.
# Copyright (c) 1999  Peter J. Acklam
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License as
# published by the Free Software Foundation; either version 2 of
# the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# The GNU General Public License can be found on the World Wide
# Web at http://www.gnu.org/copyleft/gpl.html and it can be
# obtained by writing to the Free Software Foundation, Inc.,
# 59 Temple Place - Suite 330,  Boston, MA 02111-1307, USA.
#
# Author:      Peter J. Acklam
# Time-stamp:  1999-11-12 21:29:31
# E-mail:      jacklam@math.uio.no
# URL:         http://www.math.uio.no/~jacklam

require 5.004;                  # Required version of Perl.

# Pragmatic modules.
use strict;                     # Restrict unsafe constructs.
#use diagnostics;                # Force verbose warning diagnostics.
use vars qw( $VERSION );        # Predeclare global variable names.

# Modules from the Standard Perl Library.
use File::Basename;     # Split a pathname into pieces.
use FileHandle;         # Object methods for filehandles.
use DirHandle;          # Object methods for directory handles.

# Global variables.
$VERSION = "0.001";

# Lexical (private) variables.
my $program = basename $0;                  # Name of this program.
my $pathsep = do {
                  if    ( $^O eq "MSWin32" ) { ";" }
                  elsif ( $^O eq "MacOS"   ) { ";" }
                  elsif ( $^O eq "VMS"     ) { "," }
                  else                       { ":" }
	      };
my $filesep = do {
                  if    ( $^O eq "MSWin32" ) { "\\" }
                  elsif ( $^O eq "MacOS"   ) { ":"  }
                  elsif ( $^O eq "VMS"     ) { "."  }
                  else                       { "/"  }
	      };

# Subroutine prototypes.
sub find_match ($);     # This subroutine is not def'ed at compile-time.

# Setting this to 1 prints some debugging information.
my $opt_debug = 0;

########################################################################
# glob_to_regex EXPR
#
# Given a filename glob EXPR, returns a regex.  Based on a subroutine
# found in Jeffrey E. F. Friedl's search.pl.
#
sub glob_to_regex ($) {
    my $glob = shift;
    my @parts = $glob
      =~ m{
                 \\.            # Any backslashed character
             |
                 [*?]           # Either * or ?
             |
                 \[]?[^]]*]     # A character class
             |
                 [^[\\*?]+      # Neither of the above
         }gx;

    foreach ( @parts ) {
        if ( $_ eq '*' or $_ eq '?' ) {
            $_ = ".$_";
        } elsif ( substr($_, 0, 1) ne '[' ) {
            s/^\\//;            # remove any leading backslash
            s/(\W)/\\$1/g;      # now quote anything dangerous
        }
    }
#    return join '', '^', @parts, '$';
    return join '', @parts;     # We do a substring match.
}

########################################################################
# is_valid_regex EXPR
#
# Return 1 if EXPR is a valid regex, and "" otherwise.  If EXPR is
# omitted, uses $_.  Suggested by Tom Christiansen
# <tchrist@mox.perl.com>.
#
sub is_valid_regex (;$) {
    my $pattern = @_ ? $_[0] : $_;
    my $string  = '';
    return eval { $string =~ /$pattern/, 1 };
}

########################################################################
# userexpr_to_perlexpr EXPR
#
# Convert search expression given by user to the Perl search expression.
#
sub userexpr_to_perlexpr ($) {
    my $userargs = shift;       # Arguments from user.
    $_ = lc $userargs;          # Convert to lowercase.

    my @procargs;               # List of processed arguments.
    my $parenlevel = 0;         # Parenthesis level.
    my $offset = 0;             # This variable is used when an error
                                # message is generated, indicating
                                # where the parse failed.

    while ( length ) {

        # whitespace
        if ( s{ ^               # anchor at beginning of string
                (               # start capturing into $1
                  \s+           # whitespace
                )               # end capturing into $1
             }{}x )
        {
            $offset += length $1;
        }

        # parenthesis
        if ( s{ ^               # anchor at beginning of string
                (               # start capturing into $1
                  [()]          # opening or closing parenthesis
                )               # end capturing into $1
             }{}x )
        {
            my $paren = $1;
            push @procargs, $paren;
            $paren eq "(" ? ++$parenlevel : --$parenlevel;
            # Give an error message if we found too many closing
            # parentheses.
            if ( $parenlevel < 0 ) {
                die "$program: invalid search expression\n",
                    "no opening parenthesis to match closing parenthesis\n",
                    "\n",
                    "$userargs\n",
                    " " x $offset, "^\n",
                    "\n";
            }
            $offset += length $paren;
        }

        # logical operator
        elsif ( s{ ^            # anchor at beginning of string
                   (            # start capturing into $1
                       and      # logical "and"
                     |          #   or
                       or       # logical "and"
                     |          #   or
                       not      # logical "not"
                   )            # end capturing into $1
                   # trailing context:
                   (?=          # start capturing into $2
                       [()]     # opening or closing parenthesis
                     |          #   or
                       \s       # whitespace
                     |          #   or
                       $        # end of string
                   )            # end capturing into $2
                }{}x )
        {
            push @procargs, $1;
            $offset += length $1;
        }

        # regular expression
        elsif ( s{ ^            # anchor at beginning of string
                   (            # start capturing into $1
                     /          # leading slash before regex
                     (?:        # non-capturing grouping parenthesis
                         \\.    # any escaped character
                       |        #   or
                         [^/]+  # anything not indicating end of regex
                     )*         # end grouping; repeat ad lib
                   )            # end capturing into $1
                   (            # start capturing into $2
                     /?         # trailing slash after regex
                   )            # end capturing into $2
                }{}x )
        {
            my $regex   = $1;
            my $closing = $2;
            $offset += length( $regex . $closing );
            $regex =~ s{^/}{};  # remove leading slash
            # Give an error message if we didn't found a slash
            # indicating the end of the regex.
            if ( $closing eq "" ) {
                die "$program: invalid search expression\n",
                    "missing / to end regular expression\n",
                    "\n",
                    "$userargs\n",
                    " " x $offset, "^\n",
                    "\n";
            }
            # Give an error message if the regex is invalid.
            unless ( is_valid_regex $regex ) {
                die "$program: invalid search expression\n",
                    "`$regex' is not a valid regular expression\n",
                    "\n",
                    "$userargs\n",
                    " " x $offset, "^\n",
                    "\n";
            }
            push @procargs, "/$regex/si";
        }

        # glob
        elsif ( s{ ^            # anchor at beginning of string
                   (            # start capturing into $1
                     \S+        # non-whitespace
                   )            # end capturing into $1
                }{}x )
        {
            my $glob = $1;
            my $regex = glob_to_regex $glob;    # convert to regex
            # Give an error message if the regex is invalid, assuming
            # that an invalid regex indicates an invalid glob.
            unless ( is_valid_regex $regex ) {
                die "$program: invalid search expression\n",
                    "`$glob' is not a valid glob\n",
                    "\n",
                    "$userargs\n",
                    " " x $offset, "^\n",
                    "\n";
            }
            push @procargs, "/$regex/si";
            $offset += length $glob;
        }

        # we only get here if the parsing failed
        else {
            die "$program: parse failed at indicated point\n",
                "\n",
                "$userargs\n",
                " " x $offset, "^\n",
                "\n";
        }
    }

    # Give an error message if we didn't find enough closing
    # parentheses.
    if ( $parenlevel > 0 ) {
        die "$program: invalid search expression\n",
            "missing closing parenthesis to match opening parenthesis\n",
            "\n",
            "$userargs\n",
            " " x $offset, "^\n",
            "\n";
    }

    my $perlexpr = join " ", @procargs;
    return $perlexpr
}

########################################################################
# Check how this program is called, and get the name of the file
# containing the information passed from Matlab.
#
die <<"EOF" unless @ARGV == 1;
Usage: $program FILE
FILE is a temporary file containing the information passed from Matlab.

This program should only be called by Matlab, not explicitly by a user.
EOF
my $tmpfile = shift;

########################################################################
# Read the file containing the information passed from Matlab.
#
my $tfh = new FileHandle $tmpfile, "r"
   or die "$program: can't open `$tmpfile' for reading: $!\n";
chomp( my $mpath    = $tfh->getline );
chomp( my $userexpr = do { local $/; $tfh->getline } );
$tfh->close;

# Delete the file.
unlink $tmpfile or warn "$program: can't delete `$tmpfile'\n";

# Trim the search expression given by the user.
foreach ( $userexpr ) {
    s/^\s+//;                   # Remove leading whitespace.
    s/\s+$//;                   # Remove trailing whitespace.
    s/\s+/ /g;                  # Convert whitespace to blanks.
}

########################################################################
# Convert the search expression from the user into a Perl expression.
#
my $searchexpr = userexpr_to_perlexpr $userexpr;

########################################################################
# Check the search expression.
#
$_ = '';
eval $searchexpr;
if ( $@ ) {
    die <<"EOF";
$program: Invalid search expression.

Search expression passed from user was:

    $userexpr

This was converted into the Perl expression:

    $searchexpr

EOF
}

########################################################################
# Print the search expression.
#
print <<"EOF" if $opt_debug;

Search expression passed from user was:

    $userexpr

This was converted into the Perl expression:

    $searchexpr

EOF

########################################################################
# Build the subroutine that does the matching.
#
eval <<"EOF";
sub find_match (\$) {
    local \$_ = shift;
    return ( $searchexpr );
}
EOF
if ( $@ ) {
    warn "$program: internal error: eval failed";
    die $@;
}

########################################################################
# Now scan through the directories and files.
#
my @dirs = split /\Q$pathsep\E/o, $mpath;
while ( defined( my $dir = shift @dirs ) ) {

    # Open the directory for reading.
    my $dh = new DirHandle $dir or do {
        warn "$program: can't open dir `$dir' for reading: $!\n";
        next;
    };

    # Scan through all files in the directory.
    foreach my $file ( $dh->read ) {

        # Look for private and class subdirectories.  This is a bit
        # too general, though.  Class directories can have private
        # subdirectories, but private directories can not have private
        # or class subdirectories.  We happily ignore that here.
        if ( ( $file eq "private" || $file =~ /^@[A-Za-z]\w*$/ )
             && -d ( my $subdir = $dir.$filesep.$file ) ) {
            unshift @dirs, $subdir;
            next;
        }

        # Skip everything unless it is a text file with a file name that
        # is valid for a Matlab m-file.
        next if $file eq "Contents.m" || $file eq "Readme.m";
        next unless $file =~ /^[A-Za-z]\w*\.m$/
          && -T ( my $fullfile = $dir.$filesep.$file );

        # Open file for reading, extract help text and close file.
        my $fh = new FileHandle $fullfile or do {
            warn "$program: can't open file `$fullfile' for reading: $!\n";
            next;
        };
        my @help;                       # List of help lines.
        my $fncnt = 0;                  # Function counter.
        while ( <$fh> ) {
            $fncnt += /^\s*function\b/ ? 1 : 0;
            last if $fncnt == 2;
            next unless s/^%//;         # Read next line unless it's help line.
            @help = $_;                 # Get first help line ("H1 line").
            while ( <$fh> ) {
                last unless s/^%//;     # Bail out unless it's a help line.
                push @help, $_;         # Get next help line.
            }
            last;                       # Bail out (end of help or file).
        }
        $fh->close;

        # Convert all sequences of whitespace to a single space.
        my $help = join " ", split /\s+/, join " ", @help;

        # See if we have a match.
        print "$fullfile\n" if find_match $help;

    }

}


