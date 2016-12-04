#!/local/bin/perl5 -w
#
# mlfundefre - Print function definitions in a Matlab source file.
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
# Time-stamp:  1998-07-21 16:31:56
# E-mail:      jacklam@math.uio.no (Internet)
# URL:         http://www.math.uio.no/~jacklam

#
# Function definitions may span multiple lines, so read in paragraph
# mode.
#
$/ = "";

#
# Regex for matching a Matlab function definition.
#
$ml_fun_def_regex = q!^[ \t\f]*function\b(?:[ \t]*(?:\.\.\.+[^\n]*\n[ \t]*)*(?:\[[ \t]*(?:\.\.\.+[^\n]*\n[ \t]*)*(?:[a-zA-Z][a-zA-Z0-9_]*(?:[ \t]*(?:\.\.\.+[^\n]*\n[ \t]*)*,[ \t]*(?:\.\.\.+[^\n]*\n[ \t]*)*[a-zA-Z][a-zA-Z0-9_]*)*[ \t]*(?:\.\.\.+[^\n]*\n[ \t]*)*)?\]|[a-zA-Z][a-zA-Z0-9_]*)[ \t]*(?:\.\.\.+[^\n]*\n[ \t]*)*=)?[ \t]*(?:\.\.\.+[^\n]*\n[ \t]*)*[a-zA-Z][a-zA-Z0-9_]*(?:[ \t]*(?:\.\.\.+[^\n]*\n[ \t]*)*\([ \t]*(?:\.\.\.+[^\n]*\n[ \t]*)*(?:[a-zA-Z][a-zA-Z0-9_]*(?:[ \t]*(?:\.\.\.+[^\n]*\n[ \t]*)*,[ \t]*(?:\.\.\.+[^\n]*\n[ \t]*)*[a-zA-Z][a-zA-Z0-9_]*)*[ \t]*(?:\.\.\.+[^\n]*\n[ \t]*)*)?\))?(?=[ \t]*(?:\.\.\.+[^\n]*\n[ \t]*)*\r?[\n,;])!;

#
# Regex for matching a Matlab function definition.
#
$ml_fun_def_regex = <<'EOF';
^[ \t\f]*function\b(?:[ \t]*(?:\.\.\.+[^\n]*\n[ \t]*)*(?:\[[ \t]
*(?:\.\.\.+[^\n]*\n[ \t]*)*(?:[a-zA-Z][a-zA-Z0-9_]*(?:[ \t]*(?:\
.\.\.+[^\n]*\n[ \t]*)*,[ \t]*(?:\.\.\.+[^\n]*\n[ \t]*)*[a-zA-Z][
a-zA-Z0-9_]*)*[ \t]*(?:\.\.\.+[^\n]*\n[ \t]*)*)?\]|[a-zA-Z][a-zA
-Z0-9_]*)[ \t]*(?:\.\.\.+[^\n]*\n[ \t]*)*=)?[ \t]*(?:\.\.\.+[^\n
]*\n[ \t]*)*[a-zA-Z][a-zA-Z0-9_]*(?:[ \t]*(?:\.\.\.+[^\n]*\n[ \t
]*)*\([ \t]*(?:\.\.\.+[^\n]*\n[ \t]*)*(?:[a-zA-Z][a-zA-Z0-9_]*(?
:[ \t]*(?:\.\.\.+[^\n]*\n[ \t]*)*,[ \t]*(?:\.\.\.+[^\n]*\n[ \t]*
)*[a-zA-Z][a-zA-Z0-9_]*)*[ \t]*(?:\.\.\.+[^\n]*\n[ \t]*)*)?\))?(
?=[ \t]*(?:\.\.\.+[^\n]*\n[ \t]*)*\r?[\n,;])
EOF

$ml_fun_def_regex =~ s/^\s+//;          # Remove leading whitespace.
$ml_fun_def_regex =~ s/\s+$//;          # Remove trailing whitespace.
$ml_fun_def_regex =~ s/[\r\n]//g;       # Remove embedded newlines.

#
# Extract and print the function definitions.
#
while ( <> ) {
    foreach ( /$ml_fun_def_regex/gmo ) {

        print "$_\n";
        next;

        #
        # Comment out the two statements above to get the function
        # definitions pretty-printed.
        #
        s/\.\.\.+[^\n]*\n//g;       # Remove line continuations.
        s/([][()=,])/ $1 /g;        # Add space arount [ ] ( ) , =
        s/^\s+//;                   # Remove leading whitespace.
        s/\s+$//;                   # Remove trailing whitespace.
        s/\s+/ /g;                  # Compress whitespace.
        s/\s+([,(])/$1/g;           # Remove whitespace before , (
        s/\[\s+\]/[]/g;             # [ ] -> []
        print "$_\n";

    }
}
