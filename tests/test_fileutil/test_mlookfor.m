function test_suite=test_mlookfor()
  initTestSuite;
end
function test_mlookfor_()
    -------- 
   1) To find all m-files that contain 'pseudo' and 'inverse': 
 
         mlookfor pseudo and inverse 
 
   2) To find all m-files that contain 'pseudo' or 'generalized' and 
      'inverse': 
 
         mlookfor inverse and ( pseudo or generalized ) 
 
      Note that 
 
         mlookfor inverse and pseudo or generalized 
 
      would not be the same, since it would also find m-files that 
      contain the word 'generalized' but possibly without 'inverse', 
      since 'or' has lower presedence than 'and'.  Also not that 
 
         mlookfor ( pseudo or generalized ) and inverse 
 
      would cause an error and would have to be rewritten as 
 
         mlookfor( '( pseudo or generalized ) and inverse' ) 
 
   3) To find all m-files that contain 'pseudoinverse' or 
      'pseudo-inverse', use 
 
         mlookfor pseudoinverse or pseudo-inverse 
 
      or a regular expression where the '-' is optional, like 
 
         mlookfor /pseudo-?inverse/ 
 
   Some details 
   ------------ 
   Default is to do a case-insensitive substring search.  This can 
   only be changed by using the appropriate regular expression:  To 
   do a case-sensitive search, prepend '(?-i)' to the regular 
   expression.  To match words, use word-boundary anchors '\b'. 
   Since a substring search is done, any leading or trailing glob or 
   regular expression that can match nothing is redundant. 
 
   Regular expression syntax 
   ------------------------- 
   See the perlre(5) manual page for more details about the Perl 
   regular expression syntax.  Typing 'perldoc perlre' at a shell 
   command line should give you this page.
end
