
 % ----------------------------------------------
 % HELP TEXT EXAMPLES FOR H:\matlab\strutil\BRACKETS.M : 
 
   s1 = '(How) much (wood (would a)) 
 % ----------------------------------------------
 % HELP TEXT EXAMPLES FOR H:\matlab\strutil\dm2unix.m : 
 
    LF = char(10); CR = char(13); CRLF = [ LF CR];
    str = cat(2,'test', CRLF, 'test2',CR,'test3',LF);
    str2 = dm2unix(str);
 
  
 % ----------------------------------------------
 % HELP TEXT EXAMPLES FOR H:\matlab\strutil\EXTRACT.M : 
 
   [names, p0, p1] = getnames('yy = foo(x1.^2,a_1*c,flag)');
   names2  = extract('yy = foo(x1.^2,a_1*c,flag)',p0,p1+1);   
 % ----------------------------------------------
 % HELP TEXT EXAMPLES FOR H:\matlab\strutil\FINDNAME.M : 
 
   s = 'How much wood would a Woodchuck chuck?';
   findname(s,'wo*'), findname(s,'wo*','ignore')  
 
  
 % ----------------------------------------------
 % HELP TEXT EXAMPLES FOR H:\matlab\strutil\FINDNL.M : 
 
   t = freadtxt('findnl.m');
   [inl,linenum] = findnl(t);
   t(1:inl(3))       
   t(linenum<=3)      
 % ----------------------------------------------
 % HELP TEXT EXAMPLES FOR H:\matlab\strutil\FINDQUOT.M : 
 
   t = freadtxt('findquot.m');
   [mq, mc] = findquot(t); 
   t(mq)
   t(mc)
 
  
 % ----------------------------------------------
 % HELP TEXT EXAMPLES FOR H:\matlab\strutil\GETNAMES.M : 
 
   names = getnames('yy = foo(x1.^2,a_1*c,flag)')
 
  
 % ----------------------------------------------
 % HELP TEXT EXAMPLES FOR H:\matlab\strutil\insert.m : 
 
  s1 = 'How much wood would a Woodchuck chuck?';
  no = [5 15]; s0 = strvcat('test','j');
  s2 = insert(s1,s0,no)       
 % ----------------------------------------------
 % HELP TEXT EXAMPLES FOR H:\matlab\strutil\munique.m : 
 
   t = freadtxt('munique.m');
   R = getnames(t); 
   [D,NR,C] = munique(R);
   [D(NR,:) R]  
 % ----------------------------------------------
 % HELP TEXT EXAMPLES FOR H:\matlab\strutil\rstrrep.m : 
 

 % ----------------------------------------------
 % HELP TEXT EXAMPLES FOR H:\matlab\strutil\sizestr.m : 
 
 
       sizestr(rand(3,1,4), '[', 'x')     returns '[3x1x4]'
       sizestr(rand(3,1,4), '', 'x')      returns '[3x1x4]'
       sizestr(rand(3,1,4), ' ', '-by-')  returns ' 3-by-1-by-4 '
       sizestr(rand(3,1,4), '( ', ' by ') returns '( 3 by 1 by 4 )'

 % ----------------------------------------------
 % HELP TEXT EXAMPLES FOR H:\matlab\strutil\STRIM.M : 
 
   s  = '    Testing testing 1-2-3   '
   s1 = strim(s)
 
  
 % ----------------------------------------------
 % HELP TEXT EXAMPLES FOR H:\matlab\strutil\strjoin.m : 
 
   strjoin('-by-', '2', '3', '4')  
 % ----------------------------------------------
 % HELP TEXT EXAMPLES FOR H:\matlab\strutil\strpat.m : 
 
 		strpat('aa_a_aaa_aa','aa','X') 
 		=>	X_a_aaa_X
 		strpat('aa_a_aaa_aa','a','XXX')
 		=>	aa_XXX_aaa_aa
 		strpat([[1:3] pi [5:7]],pi,nan)
 		=>	1 2 3 NaN 5 6 7
 		strpat(pi*[1:6],[4*pi 5*pi],nan)
 		=>      3.1416 6.2832 9.4248 NaN 18.85
 		strpat(pi*[1:6],[4*pi 5*pi],[1:3])
 		=>      3.1416 6.2832 9.4248 1 2 3 18.85

 % ----------------------------------------------
 % HELP TEXT EXAMPLES FOR H:\matlab\strutil\STRREPL.M : 
 
  s1 = 'How much wood would a Woodchuck chuck?';
  strrepl(s1,'chuck','pack') 
  strrepl(s1,'chuck','pack','all')  
  strrepl(s1,'Wood', 'food','all')
  strrepl(s1,'Wood', 'food','ignore','all')
 
  
 % ----------------------------------------------
 % HELP TEXT EXAMPLES FOR H:\matlab\strutil\STRWCMP.M : 
 
   s1 = 'How much wood would a Woodchuck chuck?';
   strwcmp(s1,'*chuck*') 
   strwcmp(s1,'how*')  
   strwcmp(s1,'how*','i')
 
  
 % ----------------------------------------------
 % HELP TEXT EXAMPLES FOR H:\matlab\strutil\substr.m : 
 
   s = 'How much wood would a Woodchuck chuck?';
   substr( s,  0, 1 )   
 % ----------------------------------------------
 % HELP TEXT EXAMPLES FOR H:\matlab\strutil\wordcase.m : 
 
  s1 = 'How much wood would a Woodchuck chuck?';
  wordcase(s1)
 
  