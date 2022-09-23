#! /bin/tcsh
#
#        Copyright (C) 2000-2020 the YAMBO team
#              http://www.yambo-code.org
#
# Authors (see AUTHORS file for details): HM
#
# This file is distributed under the terms of the GNU
# General Public License. You can redistribute it and/or
# modify it under the terms of the GNU General Public
# License as published by the Free Software Foundation;
# either version 2, or (at your option) any later version.
#
# This program is distributed in the hope that it will
# be useful, but WITHOUT ANY WARRANTY; without even the
# implied warranty of MERCHANTABILITY or FITNESS FOR A
# PARTICULAR PURPOSE.  See the GNU General Public License
# for more details.
#
# You should have received a copy of the GNU General Public
# License along with this program; if not, write to the Free
# Software Foundation, Inc., 59 Temple Place - Suite 330,Boston,
# MA 02111-1307, USA or visit http://www.gnu.org/copyleft/gpl.txt.
#
unalias mv rm cp

set OBJ="."
set action="clean"
set filter="changed"

#######################################################################
#
# clean_unused_variables.tcsh clean/list [dir/file] [all/changed]
#
if ( $#argv == 1 ) then
 if ( $argv[1] =~ "clean" ) then
  git ls-files --others | xargs rm -f
  exit 0
 endif
 set OBJ = $argv[1]
endif
if ( $#argv == 2 ) then
 set action = $argv[1]
 set OBJ = $argv[2]
endif
if ( $#argv == 3 ) then
 set action = $argv[1]
 set OBJ = $argv[2]
 set filter=$argv[3]
endif

echo "ACTION   :" $action
echo "FILTER   :" $filter
echo "OBJ/FILE :" $OBJ
#######################################################################

rm -f tmp  MODULE_* CLEAN LIST
if (-f $OBJ ) then
 rm -f "${OBJ}"_*
endif

######### AWK SECTION ####################

cat << EOF > AWK_split
{
 na = split (\$0,a)
 if (NR==1) {
  N=1
  start="no"
  contains="no"
 }
 MOD_NAME="MODULE_"
 if (N<10) MOD_NAME="MODULE_0"
 if (substr(a[1],1,1) == "!" ) {
  print \$0  >> MOD_NAME N
  next
 }
 if (a[1]=="module" || a[1]=="MODULE") {nextfile}
 if (a[1]=="contains") {
  contains="yes"
  print \$0 >>  MOD_NAME N
  next
 }
 if (a[1]=="function" || a[1]=="subroutine" || a[2] == "function" || a[2] == "subroutine") {
  if (contains=="yes") {
  print \$0 >>  MOD_NAME N
   next
  }
  if (start=="no") 
  {
   start="yes"
   #print N \$0  
  }
  else
  {
   start="no"
   print \$0 >>  MOD_NAME N
   N=N+1
   next
   #print N \$0  
  }
 }
 print \$0 >>  MOD_NAME N
}
EOF

cat << EOF > AWK_nl
{
 na = split (\$0,a)
 if (substr(a[1],1,1) == "!" ) {next}
 in_clean=\$0
 if (index(\$0,"!")>5) {in_clean=substr(\$0,1,index(\$0,"!")-1)}
 if (NR==1) {
  line=in_clean
  next
  na = split (in_clean,a)
 }
 if (a[1] == "&" ) {
  for (i = 2; i <= na; i++)  { 
    line = line " " a[i]
   }
 }else{
  na = split (line,a)
  if (na==0) {
   line=in_clean
   next
  }
  if( a[1]=="function" || a[1]=="subroutine" || a[2] == "function" || a[2] == "subroutine") {
    print line  >> "HEADER"
   }else{
    print line >> "PP"}
  line=in_clean
 }
}
EOF

cat << EOF > AWK_separate
{
 if (NR==1) found_call="no"
 line=\$0
 if (index(\$0,"!")>5) {line=substr(\$0,1,index(\$0,"!")-1)}
 #print line length(\$0) index(\$0,"!")
 na = split (line,a)
 var_line = "yes"
 if (index(line,"use")>0 && index (line,"only")==0 &&  index (line,"ONLY")==0 ) var_line = "no"
 if (index(line,"code_branch")>0) var_line = "no"
 if (index(line,"call")>0) found_call = "yes"
 if (index(line,"&")==0 && index(line,"real")==0  &&
     index(line,"integer")==0 && index(line,"complex")==0  &&
     index(line,"logical")==0 && index(line,"use")==0 && 
     index(line,"type")==0    && index(line,"character")==0 &&
     index(line,"REAL")==0  &&
     index(line,"INTEGER")==0 && index(line,"COMPLEX")==0  &&
     index(line,"LOGICAL")==0 && index(line,"CHARACTER")==0 ) {  var_line = "no" }
 if (var_line=="yes" && index(line,"data ")>0) { var_line = "no" }
 if (var_line=="yes" && index(line,"write ")>0) { var_line = "no" }
 if (var_line=="yes" && index(line,".not.")>0) { var_line = "no" }
 if (var_line=="yes" && index(line,".or.")>0) { var_line = "no" }
 if (var_line=="yes" && index(line,".and.")>0) { var_line = "no" }
 if (var_line=="yes" && index(line,"=")>0) { var_line = "no" }
 if (var_line=="yes" && index(line,"*")>0) { var_line = "no" }
 if (var_line=="yes" && found_call=="no") print line > "VARIABLES"
 if (var_line=="no" || found_call=="yes") print line > "BODY"
}
EOF

cat << EOF > AWK_analyze
{
 line=\$0
 if (NR ==1 ) { 
  NV=0
  V_not_found=""
 }
 if (index(line,"EOF")!=0) { 
  for (i = 1; i <= NV; i++) { 
    V_not_found=V_not_found" "V[i]
  }
  print V_not_found
 }
 if (index(line,"ONLY:")!=0) { line=substr(\$0,index(line,"ONLY:")+5) }
 if (index(line,"ONLY :")!=0) { line=substr(\$0,index(line,"ONLY :")+6) }
 if (index(line,"only:")!=0) { line=substr(\$0,index(line,"only:")+5) }
 if (index(line,"only :")!=0) { line=substr(\$0,index(line,"only :")+6) }
 if (index(line,"::")!=0) { 
   line=substr(\$0,index(line,"::")+2)
  }
 gsub(","," ",line)
 gsub("&","",line)
 na = split (line,a)
 for (i = 1; i <= na; i++)  { 
  if (index(a[i],"(")>0) { 
   a[i]=substr(a[i],1,index(a[i],"(")-1) 
   if (index(a[i+2],")")>0) {
    a[i+1]=" "
    a[i+2]=" "
   }
   if (index(a[i+1],")")>0) {a[i+1]=" "}
  }
  gsub("\\\("," ",a[i])
  gsub("\\\)"," ",a[i])
  gsub("\\\["," ",a[i])
  gsub("\\\]"," ",a[i])
  gsub("="," ",a[i])
  gsub("'"," ",a[i])
  gsub("/"," ",a[i])
  gsub("!"," ",a[i])
  gsub(":"," ",a[i])
  gsub("*"," ",a[i])
  gsub(" ","",a[i])
  tmp_var=a[i]
  for (j = 1; j <= 100; j++)  { 
   gsub("0","",tmp_var)
   gsub("1","",tmp_var)
   gsub("2","",tmp_var)
   gsub("3","",tmp_var)
   gsub("4","",tmp_var)
   gsub("5","",tmp_var)
   gsub("6","",tmp_var)
   gsub("7","",tmp_var)
   gsub("8","",tmp_var)
   gsub("9","",tmp_var)
  }
  is_var="yes"
  if (index(a[i],"%") > 0 ) {is_var="no"}
  if (length(a[i]) == 0 || length(tmp_var) ==0 ) {is_var="no"}
  if (is_var == "yes") 
  {
   NV++
   V[NV]=a[i]
   #print "|"a[i]"|"tmp_var"|"
  }
 }
}
EOF

######### AWK SECTION ####################
#

set FILES = (  )

if (-d $OBJ) then
 foreach kind ( "modified:" )  
  git status -uno $OBJ | grep $kind | grep -v '\.pl' | grep  -v '\.pm'  | grep -v '\.c' |grep -v '\.m4' | grep -v '\.git' |grep -v '\.\.\/' > LIST
  cat LIST | grep -v "mod_" | grep -v "Makefile" | grep -v "configure" | grep -v "\.h" | grep -v "\.object" | grep -v "\.tcsh" | grep -v "lib/yambo" > LIST
  sed -i -e 's/new file/new_file/g' LIST
  if ($kind =~ "renamed") then
   set FILES=($FILES `awk '{print $4}' "LIST"`)
  else
   set FILES=($FILES `awk '{print $2}' "LIST"`)
  endif
 end
endif
if ($filter =~ "all" && $OBJ =~ ".") then
 set FILES = `find src -name '*.F' `
 set FILES = ($FILES `find ypp -name '*.F' `)
 set FILES = ($FILES `find interfaces -name '*.F' `)
endif
if ($filter =~ "all" && $OBJ !~ ".") then
 set FILES = `find $OBJ -name '*.F'`
endif

if (-f $OBJ) then
 set FILES=$OBJ
endif

foreach file ($FILES)
 #echo "processing $file..."
 gawk -f AWK_split ${file}
 @ N_unused=0
 foreach source (MODULE_*)
  #echo "$file $source..."
  gawk -f AWK_nl $source
  if (! -f PP || ! -f HEADER) then
   echo "-----------------------------------------------------------------------------"
   echo "WARNING->"$file"/"$source": impossible to split"
   echo "-----------------------------------------------------------------------------"
   continue
  endif
  mv PP ${file}"_PP"
  mv HEADER ${file}"_HEADER"
  gawk -f AWK_separate ${file}"_PP"
  if (! -f VARIABLES) then
   echo "-----------------------------------------------------------------------------"
   echo "WARNING->"$file"/"$source": VARIABLES not found"
   echo "-----------------------------------------------------------------------------"
   continue
  endif
  mv VARIABLES ${file}"_VARIABLES"
  mv BODY ${file}"_BODY"
  echo "EOF" >>  ${file}"_VARIABLES"
  echo "EOF" >>  ${file}"_BODY"
  #gawk -f AWK_analyze "${file}"_VARIABLES
  set VARIABLES = `gawk -f AWK_analyze "${file}"_VARIABLES`
  set UNUSED = ( )
  foreach var ($VARIABLES)
   set var = `echo $var | sed 's/-//g'`
   set CHECK = `grep $var "${file}"_VARIABLES| wc -l`
   set FOUND = `grep $var "${file}"_BODY| wc -l`
   set HEAD  = `grep $var "${file}"_HEADER| wc -l`
   if ($FOUND == 0 && $CHECK == 1 && $HEAD == 0 && "$var" != "code_branch" ) then
    @ N_unused ++
    set UNUSED = ($UNUSED $var)
    #echo -n $var $source  $N_unused 
    #echo -n "..."$var"(removed) "

cat << EOF > AWK_replace
{
 na = split (\$0,a)
 if (substr(a[1],1,1) == "!" ) 
 {
  print \$0
  next
 }
 #
 gsub("\\\(","_PARL_")
 gsub("\\\)","_PARR_")
 # Isolate the real var
 PATTERN="${var}"
 line=\$0
 gsub(",${var}"," ${var}",line)
 na = split (line,a," ")
 for (i = 1; i <= na; i++)  { 
   if (index(a[i],"${var}_PARL_")>0) 
   {
    #print "S" i a[i]
    PATTERN=a[i]
    #print PATTERN,"PATT"
    #print \$0
    #gsub(PATTERN,"")
    #print \$0
   }
 }
 #
 if (index(\$0,PATTERN)==0)
 {
  gsub("_PARL_","(")
  gsub("_PARR_",")")
  print \$0
  next
 }
 gsub(PATTERN",&","&")
 gsub(PATTERN" ,&","&")
 gsub(PATTERN"  ,&","&")
 gsub(PATTERN",","")
 gsub(PATTERN" ,","")
 gsub(PATTERN"  ,","")
 gsub(","PATTERN,"")
 gsub(", "PATTERN,"")
 gsub(",  "PATTERN,"")
 gsub(":"PATTERN,"")
 gsub(": "PATTERN,"")
 gsub(":  "PATTERN,"")
 gsub(PATTERN",","")
 # Check if the USE call is empty
 print_the_var="yes"
 if (index(\$0,"use")>1 && index(\$0,":")==0) {print_the_var = "no" } 
 if (index(\$0,"logical")>1 && index(\$0,"::")==0) {print_the_var = "no" }
 if (index(\$0,"integer")>1 && index(\$0,"::")==0) {print_the_var = "no" }
 if (index(\$0,"character")>1 && index(\$0,"::")==0) {print_the_var = "no" }
 if (index(\$0,"real_PARL_")>1 && index(\$0,"::")==0) {print_the_var = "no" }
 if (index(\$0,"complex_PARL_")>1 && index(\$0,"::")==0) {print_the_var = "no" }
 gsub("_PARL_","(")
 gsub("_PARR_",")")
 if (print_the_var=="yes") {print \$0}
}
EOF
    gawk -f AWK_replace $source > CLEAN  
    mv CLEAN $source
   endif
  end
  #echo " " $#UNUSED "("$N_unused")"
  if ($#UNUSED>0) then
   echo "-----------------------------------------------------------------------------"
   echo $file"/"$source":"$UNUSED
   echo "-----------------------------------------------------------------------------"
  endif
 end
 if ($N_unused>0 && $action !~ "quiet") then
  cat MODULE_* > $file
 endif
 rm -f tmp "${file}"_* MODULE_* CLEAN LIST
end

rm -f AWK* 

exit 0
