#!/local/bin/perl5 -wn
#
# mlisfun - Returns 1 for Matlab functions and 0 for scripts.
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
# Time-stamp:  1998-07-21 16:38:36
# E-mail:      jacklam@math.uio.no (Internet)
# URL:         http://www.math.uio.no/~jacklam

s/(%|\.\.\.+).*//;
if ( /\S/ ) {
    print /^\s*function\b/ ? 1 : 0, "\n";
    exit;
}
