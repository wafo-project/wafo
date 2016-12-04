#!/local/bin/perl5 -w
#
# mlcall - Call Matlab to execute commands and print result to stdout.
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
# Time-stamp:  1999-05-31 21:08:30
# E-mail:      jacklam@math.uio.no
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
use File::Basename;     # Split a pathname into pieces.
#use File::Path;         # Create or remove a series of directories.
#use File::Find;         # Traverse a file tree.
#use File::Copy;         # Copy files or filehandles.
#use Text::Tabs;         # Expand and unexpand tabs.
use Getopt::Long 2.10;  # Extended processing of command line options.
use FileHandle;         # Object methods for filehandles.
#use Cwd;                # Get pathname of current working directory.

#
# Global variables.
#
$VERSION = "2.000";

#
# Lexical (private) variables.
#
my $program = basename $0;              # Name of this program.
my $is_ms = $^O =~ /^(MS)?(Win32|DOS)$/i;

my $dryrun;             # See usage (and manual) part below.
my $printusage;         # See usage (and manual) part below.
my $printmanual;        # See usage (and manual) part below.
my $verbose;            # See usage (and manual) part below.
my $printversion;       # See usage (and manual) part below.

my $matlabexe;          # Location of matlab executable.

my $wrapperfun  = "tmpw____";          # Name of wrapper function.
my $usrcmdsfun  = "tmpu____";          # Name of user commands function.
my $wrapperfile = "$wrapperfun.m";     # Name of wrapper file.
my $usrcmdsfile = "$usrcmdsfun.m";     # Name of user commands file.
my $logfile     = "$wrapperfun.log";   # Name of Matlab output file.
my $cshscript   = "matcsh__";          # Name of csh script (UNIX only).

my $beginsep;           # User command output start separator.
my $endsep;             # User command output end separator.

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
# Subroutine that prints program usage.
#
sub print_usage () {
    print <<"EOF";
Usage: $program [-dhmvV] [commands...]
Call Matlab to run commands and print output in console window.

  -d, --dryrun      dryrun; print what Matlab would process and exit
  -h, --help        print usage and exit
  -m, --manual      print manual and exit
  -v, --verbose     verbose mode; tell what is going on
  -V, --version     print version information and exit

Report bugs to jacklam\@math.uio.no
EOF
}

#
# Subroutines that prints program manual.
#
sub print_manual () {

    my $pager = $ENV{PAGER} || "more";
    my $ofh;
    if ( open PIPE, "|$pager" ) {
        $ofh = \*PIPE;
    } else {
        $ofh = \*STDOUT;
    }
    print $ofh <<"EOF";

NAME
     $program - capture Matlab output to console window

SYNOPSIS
     $program [-dhmvV] [commands...]

DESCRIPTION
     $program takes any amount of Matlab code and passes it to Matlab.
     The output from Matlab is captured and printed on standard output.

     Only the output from processing the user commands is returned. Any
     output during Matlab's initialization (processing matlabrc.m,
     startup.m etc.) and exiting (printing flop count, processing
     finish.m etc.) is not included.

     The Matlab code may be given through standard input or on the
     command line. If commands are given both through standard input and
     on the command line, the code on standard input is ignored.

OPTIONS
     -d --dryrun
          Print what Matlab would process on standard output and exit
          successfully. This is useful to ensure that the commands given
          on the command line are fed to Matlab as intended.

     -h --help
          Print a usage message on standard output and exit
          successfully.

     -m --manual
          Print this manual page on standard output and exit
          successfully.

     -v --verbose
          Verbose mode. Print what is currently being done on standard
          output.

     -V --version
          Print version information on standard output then exit
          successfully.

EXAMPLES
     Execute the m-file foo.m

          $program foo

     Calculate pi

          $program '4*atan(1)'    (UNIX)
          $program "4*atan(1)"    (Windows)

     Invert a matrix (note the double quotes)

          $program 'format short g; inv( [ 2 3 ; 4 5 ] )'    (UNIX)
          $program "format short g; inv( [ 2 3 ; 4 5 ] )"    (Windows)

     The arguments may also be typed in like this. In this case the code
     is given through standard input so the code does not need to be
     surrounded by quotes.  The ending ^D is a Ctrl-D character (on
     Windows use Ctrl-Z in stead).

          $program
          format short g
          inv( [ 2 3 ; 4 5 ] )
          ^D

     Windows only:  If $program is a batch file (.bat file) that makes
     perl run $program.pl, commands may be given through standard input
     like this

          type foo.m | command.com /c $program

INPUT FILES
     The input must be given exactly as it would be given on the Matlab
     command line. Thus, any Matlab m-files must be given without the .m
     suffix. Explicit file names are invalid input.

OUTPUT FILES
     Temporary files are
          ./$wrapperfile      Name of wrapper file
          ./$usrcmdsfile      Name of user commands file
          ./$logfile          Name of Matlab output file
          ./$cshscript        Name of csh script (UNIX only)

KNOWN PROBLEMS
     None that are not mentioned above.

BUGS
     None known.

SEE ALSO
     The Matlab documentation.

NOTES
     Arguments passed on the command line will be pre-processed by the
     shell (command interpreter), which may lead to unexpected results.
     For instance, on Windows machines, command.com removes all commas
     and semicolons unless they are inside double quotes. If $program is
     fed [ 2 3 ; 4 5 ]^2, Matlab will receive [ 2 3 4 5 ]^2 which is not
     the same and will produce an error message because you can't square
     a vector. In this case the argument should be put inside double
     quotes. Quotes are not necessary if the arguments are given through
     standard input.

VERSION
     This man page documents $program $VERSION.

AUTHOR
     Peter J. Acklam <jacklam\@math.uio.no>

COPYRIGHT
     Copyright (c) 1998 Peter J. Acklam. All rights reserved.

LICENSE
     This program is free software; you can redistribute it and/or
     modify it under the terms of the GNU General Public License.
EOF

close $ofh;

}

#
# Subroutine that removes the temporary files that we have made.
#
sub cleanup () {
    unlink $usrcmdsfile
      or warn "$program: can't unlink $usrcmdsfile: $!\n"
        if -f $usrcmdsfile;
    unlink $wrapperfile
      or warn "$program: can't unlink $wrapperfile: $!\n"
        if -f $wrapperfile;
    unlink $logfile or warn "$program: can't unlink $logfile: $!\n"
      if -f $logfile;
    unlink $cshscript or warn "$program: can't unlink $cshscript: $!\n"
      if -f $cshscript;
}

#
# Subroutine that finds the location of the Matlab executable when this
# program is run under Windows.  This routine should perhaps be able to
# seach for the Matlab executable, but that is not implemented.
#
sub get_matlab_executable () {
    $matlabexe = 'x:\matlab\bin\matlab.exe';
    my $line = __LINE__ - 1;        # Line number of previous line.
    die <<"EOF"
$program: can't find $matlabexe

Please edit the file $0 line $line
and make the variable \$matlabexe point to the Matlab executable.
EOF
      unless -f $matlabexe;
}

#
# Subroutine that writes a Matlab m-file with the user commands.
#
sub print_user_commands_file () {

    print "Writing Matlab m-file with user commands.\n" if $verbose;

    $usrcmdsfile = "-" if $dryrun;      # Write on stdout if dryrun.
    my $ofh = new FileHandle ">$usrcmdsfile"
      or cleanup, die "$program: can't open $usrcmdsfile: $!\n";

    if ( @ARGV ) {                      # Matlab code on cmd line.
        foreach ( @ARGV ) { s/^"(.*)"$/$1/ if $is_ms }
        print $ofh join " ", @ARGV;
        print $ofh "\n";
    } else {                            # Matlab code through stdin.
        while ( <STDIN> ) { print $ofh $_ }
    }
    $ofh->close;

}

#
# Subroutine that writes a wrapper to the file with the user
# commands. This is necessary because we want Matlab to continue even if
# the file with the user commands is not syntactically correct.
#
sub print_wrapper_file () {

    my $ofh = new FileHandle ">$wrapperfile"
      or cleanup, die "$program: can't open $wrapperfile: $!\n";

    #
    # We don't want to include the output which is printed when Matlab
    # is initializing (processing matlabrc.m, startup.m, etc.) and
    # exiting (printing flop count etc.), so we include some separators.
    #
    $beginsep = "==> Output from user commands begins here <==";
    $endsep   = "==>  Output from user commands ends here  <==";

    print "Writing Matlab wrapper m-file.\n" if $verbose;

    print $ofh <<"EOF";
errortrap on;
fprintf( '\\n%s\\n', '$beginsep' );
$usrcmdsfun;
fprintf( '\\n%s\\n', '$endsep' );
exit;
EOF
    $ofh->close;

}

#
# Subroutine that calls Matlab to execute the wrapper file.
#
sub execute_wrapper_file () {
    print "Calling Matlab to execute the commands.\n" if $verbose;
    my $status;
    if ( $is_ms ) {
        $status = system <<"EOF";
$matlabexe -nosplash -minimize -r $wrapperfun -logfile $logfile
EOF
    } else {
        my $ofh = new FileHandle ">$cshscript"
          or cleanup, die "$program: can't open $cshscript: $!\n";
        print $ofh <<'EOF';
#!/bin/csh -f
# Clear the DISPLAY.
unsetenv DISPLAY
# Call MATLAB with the appropriate input and output,
# make it immune to hangups and quits using ''nohup''.
nohup matlab < $1 > $2
EOF
        $ofh->close;
        chmod 0700, $cshscript;
        $status = system <<"EOF";
$cshscript $wrapperfile $logfile
EOF
    }
    cleanup, die "$program: matlab exited with non-zero exit status\n"
      if $status;
}

#
# Subroutine that prints the Matlab output from the user commands.
#
sub print_results () {

    print "Processing Matlab output file.\n",
      "Here is the output from Matlab:\n" if $verbose;

    #
    # Print lines between the begin and end separators.  The end
    # separator might not be on a line by itself; it might be preceded
    # by Matlab output.
    #
    my $ifh = new FileHandle $logfile
      or cleanup, die "$program: can't open $logfile: $!\n";
    while ( <$ifh> ) {
        last if /^\Q$beginsep\E$/o;
    }
    while ( <$ifh> ) {
        last if s/^\Q$endsep\E$//o;
        print;
    }
    $ifh->close;

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
    GetOptions (
        "d|dryrun"      => \$dryrun,
        "h|help"        => \$printusage,
        "m|man"         => \$printmanual,
        "v|verbose"     => \$verbose,
        "V|version"     => \$printversion,
    );
}

#
# Subroutine that processes the user options.
#
sub process_options () {
    print_usage,   exit if $printusage;
    print_manual,  exit if $printmanual;
    print_version, exit if $printversion;
}

########################################################################
# Main part.
########################################################################

get_options;
process_options;
get_matlab_executable if $is_ms;
print_user_commands_file;
exit if $dryrun;
print_wrapper_file;
execute_wrapper_file;
print_results;
cleanup;
