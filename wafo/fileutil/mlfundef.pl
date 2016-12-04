#!/local/bin/perl5 -w
#
# mlfundef - Print function definitions in a Matlab source file.
# Copyright (C) 1998  Peter J. Acklam
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or (at
# your option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
#
# The GNU General Public License can be found on the World Wide Web at
# http://www.gnu.org/copyleft/gpl.html and it can be obtained by writing
# to the Free Software Foundation, Inc., 59 Temple Place - Suite 330,
# Boston, MA 02111-1307, USA.
#
# Author:      Peter J. Acklam
# Time-stamp:  1998-10-23 21:33:38
# E-mail:      jacklam@math.uio.no (Internet)
# URL:         http://www.math.uio.no/~jacklam

########################################################################
# Miscellaneous "high-level" stuff.
########################################################################

require 5.004;                  # Required version of Perl.

#
# Pragmatic modules.
#
use strict;                     # Restrict unsafe constructs.
#use diagnostics;                # Force verbose warning diagnostics.
use vars qw( $VERSION );        # Predeclare global variable names.

#
# Modules from the Standard Perl Library.
#
use File::Basename;             # Split file into dir/file/suffix.
use Getopt::Long;               # Extended processing of options.
use FileHandle;                 # Object methods for filehandles.

#
# Global variables.
#
$VERSION = "1.005";

#
# Lexical (private) variables.
#
my $program = basename $0;      # Name of this program.
my $printusage;                 # See option part below.
my $printversion;               # See option part below.
#my $verbose;                    # See option part below.
my $printargnum;                # See option part below.
my $printarglist;               # See option part below.
my $printlinenum;               # See option part below.
my $printsynopsis;              # See option part below.

########################################################################
# Subroutines.
########################################################################

#
# Subroutine that prints program version number.
#
sub print_version () {
    print <<"EOF";

$program $VERSION

Copyright (c) 1998 Peter J. Acklam. All rights reserved.
This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License.

EOF

}

#
# Subroutines that prints program usage.
#
sub print_usage () {
    print <<"EOF";
Usage: $program [-hvV] [files...]
Print function definitions in a Matlab source file.

  -a, --numargs     print number of input and output arguments
  -n, --linenum     print line numbers
  -h, --help        print this help and exit
  -l, --listargs    list input and output arguments
  -s, --synopsis    print function synopsis
  -V, --version     print version information and exit

Report bugs to jacklam\@math.uio.no
EOF
}

#
# Subroutine that processes a line of Matlab code with a function
# declaration.
#
sub process_declaration ($\$\@\@) {

    my $declaration = shift;    # Function declaration string.

    my $name        = shift;    # Ref to scalar with function name.
    my $inputargs   = shift;    # Ref to list of input arguments.
    my $outputargs  = shift;    # Ref to list of output arguments.

    my ( $input, $output );
    ( $output, $$name, $input ) = $declaration =~
      m/
        ^                   # Anchor at beginning of line.
        \s*                 # Optional whitespace.
        function            # The function keyword.
        \b                  # Word boundary.
        \s*                 # Optional whitespace.
        (?:                 # Non-capturing grouping parenthesis.
          (                 # Capture output arguments.
              \[            # Either an opening left bracket
                [^\]]*      #   bracket contents (list of output args)
              \]            #   closing right bracket
            |               # or else
              \w+           #   output argument name without brackets.
          )                 # End capturing.
          \s*               # Optional whitespace.
          =                 # Assignment operator.
        )?                  # Optional. There need not be output args.
        \s*                 # Optional whitespace.
        (                   # Capture function name.
          \w+               # Function name.
        )                   # End capturing.
        \s*                 # Optional whitespace.
        (                   # Capture input arguments.
          \(                # Opening literal parenthesis
            [^)]*           #   bracket contents (list of input args).
          \)                # Closing literal parenthesis.
        )?                  # End capturing. Optional since there need
                            #   not be input args.
       /x;

    # Build a list of input arguments.
    if ( $input ) {
        $input =~ s/^\(?\s*//;  # Strip leading paren and whitespace.
        $input =~ s/\s*\)?$//;  # Strip trailing paren and whitespace.
        @$inputargs =  split /\s*,\s*/, $input;
    }

    # Build a list of output arguments.
    if ( $output ) {
        $output =~ s/^\[?\s*//; # Strip leading bracket and whitespace.
        $output =~ s/\s*\]?$//; # Strip trailing bracket and whitespace.
        @$outputargs =  split /\s*,\s*/, $output;
    }

}

#
# Subroutine that builds a "pretty-printed" function synopsis.
#
sub build_synopsis ($\@\@) {

    my $name        = shift;    # Scalar with function name.
    my $inputargs   = shift;    # Ref to list of input arguments.
    my $outputargs  = shift;    # Ref to list of output arguments.

    my $synopsis;               # String with function synopsis.

    $synopsis .=  join ", ", @$outputargs;
    $synopsis  =  "[ $synopsis ]" if @$outputargs > 1;
    $synopsis .=  " = " if @$outputargs;
    $synopsis .=  $name;
    $synopsis .=  "( " . join( ", ", @$inputargs ) . " )"
      if @$inputargs;

    return $synopsis;

}

#
# Subroutine that prints function declaration.
#
sub print_declaration ($$\@\@) {

    my $linenum     = shift;    # Line number where function is def'ed.
    my $name        = shift;    # Function name.
    my $inputargs   = shift;    # Ref to list of input arguments.
    my $outputargs  = shift;    # Ref to list of output arguments.

    my $synopsis = build_synopsis $name, @$inputargs, @$outputargs
      if $printsynopsis;

    if ( $printarglist ) {
        print "\n";
        print "Line:  $linenum\n" if $printlinenum;
        print "Name:  $name\n";
        print "In:    " . join( ", ", @$inputargs ) . "\n";
        print "Out:   " . join( ", ", @$outputargs ) . "\n";
        print "Syn.:  $synopsis\n" if $printsynopsis;
    } else {
        printf "%6d:  ", $linenum if $printlinenum;
        print $name;
        print " ", " " x ( 32 - length $name ) if $printargnum;
        print " (" . @$inputargs . " in, " . @$outputargs . " out)"
          if $printargnum;
        print "\n";
        print "         $synopsis\n" if $printsynopsis;
    }

}

#
# Subroutine that gets the user options.
#
sub get_options () {

    # Configure option handling.
    #
    # Getopt::Long::config is deprecated in favour of
    # Getopt::Long::Configure which was introduced in
    # $Getopt::Long::VERSION 2.17. We use the config routine for
    # backwards compatibility.
    #
    Getopt::Long::config( "no_ignore_case", "no_auto_abbrev", "bundling" );
#    Getopt::Long::Configure( "no_ignore_case", "no_auto_abbrev", "bundling" );

    # Read options and arguments and initialize variables.
    GetOptions(
        "a|numargs"  => \$printargnum,    # Print number of arguments.
        "h|help"     => \$printusage,     # Print help text.
        "l|listargs" => \$printarglist,   # List arguments.
        "n|linenum"  => \$printlinenum,   # Print line numbers.
        "s|synopsis" => \$printsynopsis,  # Print function synopsis.
#        "v|verbose"  => \$verbose,        # Verbose mode.
        "V|version"  => \$printversion,   # Print version number.
    );
}

#
# Subroutine that processes the user options.
#
sub process_options () {
    print_usage,   exit if $printusage;   # Print program usage.
    print_version, exit if $printversion; # Print version number.
}

########################################################################
# Main part.
########################################################################

get_options;
process_options;

while ( my $file = shift @ARGV ) {

    # Convert directory separators to UNIX style.
    $file =~ s!\\!/!g if $^O =~ /^MS(Win32|DOS)$/;

    # Perform glob expansion if necessary.
    if ( $file =~ m! (^|/)~ | [*?{}[\]^] !x ) {
        unshift @ARGV, glob $file;
        next;
    }

    # Add .m suffix if it is missing.
    $file .= ".m" unless $file =~ /\.[mM]$/;

    # Give a warning if the file does not exist or is not a file.
    warn( "$program: can't find $file\n" ),  next unless -e $file;
    warn( "$program: not a file: $file\n" ), next unless -f $file;

    # Process file.
    my $fh = new FileHandle $file
      or die "$program: can't open $file: $!\n";
    while ( <$fh> ) {

        # The function keyword can only be preceded by whitespace.
        next unless /^\s*function\b/;

        my $linenum     = $.;   # Line number where declation begins.
        my $declaration = "";   # Variable that holds declaration.
        my $continued   = 0;    # Indicates continued line.

        do {

            # Strip comment (easy here because strings are not allowed).
            s/%.*$//;

            # If line is continued (ends in three or more dots and maybe
            # some whitespace) replace with a blank.
            $continued = s/\.\.\.+.*$/ /;

            # Append this line to the previous if we are reading a
            # continued line.
            $declaration .= $_;

            # Get next line if line is continued.
            $_ = <$fh> if $continued;

        } while $continued;

        my ( $name, @inputargs, @outputargs );
        process_declaration $declaration, $name,
          @inputargs, @outputargs;
        print_declaration $linenum, $name, @inputargs, @outputargs;

    }

    $fh->close;

}
