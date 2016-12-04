#!/local/bin/perl5 -wp
#
# mlstrip - Strip comments from a matlab m-file.
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
# Time-stamp:  1999-07-01 23:20:29
# E-mail:      jacklam@math.uio.no
# URL:         http://www.math.uio.no/~jacklam
#
# Bugs
#
# Implicit comments (anything after a line continuation) are not
# removed.
#
# How to use
#
# UNIX:    mlstrip original.m > stripped.m
# MS-DOS:  perl mlstrip.pl original.m > stripped.m

#
# Remove a comment line.
#
s<^\s*%.*><>;

#
# Remove a comment on a line possibly containing Matlab code. It is done
# by finding the part of a line which is not a comment and then
# replacing the whole line (including the comment) with what was found.
# A comment is started by a percent sign, but a percent sign might be
# embedded in a string (in which case it does of course not start a
# comment).  Strings are delimited by single quotes, but finding the
# beginning of a string is non-trivial since single quotes are also used
# as transpose operator.
#
s<
  ^
  (                     # Grouping parenthesis (contents goes to $1).
    (?:                 # Non-backreferencing grouping parenthesis.

        # Match anything that is neither a comment nor a string.
        (?:             # Non-backreferencing grouping parenthesis.
            [])}\w.]    # Either a character followed by
            '+          #    one or more transpose operators
          |             # or else
            [^'%]       #   any character except single quote (which
                        #   starts a string) or a percent sign (which
                        #   starts a comment).
        )+              # Match one or more times.

      |

        # Match a string.
        '               # Opening single quote that starts the string.
          [^'\n]*       # Zero or more chars that are neither single
                        #   quotes (special) nor newlines (illegal).
          (?:           # Non-backreferencing grouping parenthesis.
            ''          # An embedded (literal) single quote character.
            [^'\n]*     # Again, zero or more chars that are neither
                        #   single quotes nor newlines.
          )*            # Match zero or more times.
        '               # Closing single quote that ends the string.

    )*          # Match zero or more times.
  )
  .*            # What remains must be a comment.
><$1>x;         # Replace everything with everything except the comment.

#
# Remove trailing whitespace from each line.
#
s<\s+$><\n>;
