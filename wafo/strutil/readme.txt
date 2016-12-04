Below I present a set of string and text processing  functions  
and some other useful functions.

   MATLAB is certainly not a language for text processing,
Yet with growing popularity and diversity of MATLAB
applications at least some elementary text functions
are needed.
  Although MATLAB offers a convenient basis for making
such utilities, the current level of actual functions to
deal with text is inadequate.
  About the only routines that perform some text
manipulations (aside from very useful but trivial such as
STRCMP, STR2MAT) are FINDSTR and STRREP functions. Yet their
capabilities are severely limited even for the purposes which
they are designed for.

  Take for example STRREP function which replaces all
occurrences of one string withing another by a third string.
More often than not one needs a selective replacement,
such as only those occurrences which are names (words)
themselves and not parts of bigger words.
  For example one may want to replace all names "x1"
but not "x10", "x1_temp" or "dx1". Or, using example
from the FINDSTR finction:
 STRREP('How much wood would a woodchuck chuck?','chuck','pack')
should probably result in
'How much wood would a woodchuck pack?'     and not in
'How much wood would a woodpack pack?'.

  Similar considerations are applicable to FINDSTR.
Again one often needs to know whether found occurrence is a
name (word) or a part of a word.
  It is also quite obvious that such a function invites
a "wildcard" option - ability to match a specified number or
arbitrary number of any characters.
  Nothing like that is possible with FINDSTR.

  Over a time I wrote some routines which perform different
text-processing functions, including those described above.

  Firstly, there are substitutes for FINDSTR and STRREP -
FINDNAME and STRREPL functions. The function FINDNAME in
particular has some useful wilcard options:

>> s='How much wood would a woodchuck chuck?';
>> [names,p] = findname(s,'wo*')

names =

wood
would
woodchuck


p =

    10    15    23

>> [names,p] = findname(s,'wo??')

names =

wood


p =

    10


  Now let's go back to the original question.
  Here is the function FINDREF which can determine the functions
used in a file. Essentially it is similar to a MATLAB parser
written  in MATLAB itself.
  This is quite a non-trivial procedure.
One needs to analyze the file, check for comments, explicit
strings in quotes, determine beginnings and ends of statements,
extract all legal names, find which side of the assignment operator
they are on, find whether they are variables or function names.
In addition there are some non-trivial choices to be made: one
probably does not want to see all the trivia such as "zeros",
"ones", "find" etc. It is far from clear which function should be
excluded and which options should a user have in this respect.

  The function FINDREF attempts to accomplish this. Normally it
returns all function names used inside a given file excluding some
elementary language and math functions (written in accompanying
"library" file FNAMES.M). It also has the options to include
various classes of these elementary functions and also control
the exist status of returned function names (for example, return
only built-in functions).
FINDREF also returns the exist status of all found names. It also
has built-in potential to give much more information about the
processed file and all functions in it.
  It is fully vectorized and parses a file very quickly. About
80% of time is spent to check the exist status of potential
function names.
  For example, let us check which functions are used in the
FINDREF program itself:

>> findref findref
 Extracted all names. Checking EXIST status...

ans =

error    
nargchk  
str2mat  
int2str  
ischar   
lower    
strtok   
freadtxt 
fnames   
strvcat  
deblank  
findquot 
diff     
char     
brackets 
isempty  
getnames 
munique  
exist    

 There is also a routine READCLM which can read
an ASCII data file with header into a matrix.
It skips the commented or specified number of lines,
then extracts header (all lines until numerical
data, each line into a row of a string matrix)
then determines number of columns in a data matrix
and reads the rest of the file into it.

  Once you obtained a header (and you variable names)
you can use another routine GETNAMES which extracts
all legal MATLAB names from a string (each name into
a row of a string matrix).
  Then you can make necessary assignment with the EVAL
command.

  Here is an example.

  Suppose we have a data file TEXTCLM with header
(and comment lines):

>> type textclm

% This is a test for READCLM routine
A x_tmp, y zz
Now some data:

0.0727    0.1665    0.4940    0.4644
   0.6316   0.4865    0.2661    0.9410
    0.8847    0.8977    0.0907    0.0501
    0.2727    0.9092    0.9478    0.7615
    0.7665    0.9047    0.5007 0.8278



>> [outdata,head] = readclm('textclm')

outdata =

    0.0727    0.1665    0.4940    0.4644
    0.6316    0.4865    0.2661    0.9410
    0.8847    0.8977    0.0907    0.0501
    0.2727    0.9092    0.9478    0.7615
    0.7665    0.9047    0.5007    0.8278


head =

A x_tmp, y zz
Now some data:


  As you can see the program is capable to detect the
number of columns by itself.
  Now extract the names from the first header lines:

>> Names=getnames(head(1,:))

Names =

A
x_tmp
y
zz

  Now you can assign the data to the obtained variable
names with EVAL command:

  for jj=1:size(Names,1)
    eval([Names(jj,:) ' = outdata(:,jj);']);
  end

  Something similar can be done with the second row
of the header matrix.
 
  For more information see HELP for REFER and other accompanying
functions.

  Regards, Kirill

``````````````````````````````````````````````````````````
  Kirill Pankratov, Ph.D.

  Department of Earth, Atmospheric & Planetary Sciences,
  Massachusetts Institute of Technology
  Cambridge, MA, 02139

  kirill@plume.mit.edu
  Office (617)-253-5938
''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
Revised pab 13.12.2000
- updated readme to conform to the updated REFER Toolbox
- added readme.txt for readclm.m
