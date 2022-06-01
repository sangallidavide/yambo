#!/bin/bash
#
#        Copyright (C) 2000-2022 the YAMBO team
#              http://www.yambo-code.org
#
# Authors (see AUTHORS file for details): AM
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
lock_files=`find $dir -name '*.lock'`
sorted_locks=$(echo "$lock_files"|tr " " "\n"|sort|uniq|tr "\n" " ")
#
# Locks -> string
#
lock_string=""
save_dir=""
for lock in $sorted_locks
do
 lock=`echo $lock | sed "s/.lock//"`
 lock=`basename $lock`
 lock_string="${lock} ${lock_string}"
 save_dir="${lock}_${save_dir}"
done
if [[ -z $save_dir ]] ; then
 return
fi
save_dir="${save_dir}.save"
#
# -D* -> string
#
flag_string=""
restore_dir=""
for flag in $sorted_precomps
do
 flag=`echo $flag | sed "s/\-D_//"`
 flag_string="${flag} ${flag_string}"
 restore_dir="${flag}_${restore_dir}"
done
restore_dir="${restore_dir}.save"
#
# Check for missing/new precomp flags
#
missing=`comm -23 <(tr ' ' $'\n' <<< $lock_string | sort) <(tr ' ' $'\n' <<< $flag_string | sort)`
new=`comm -23 <(tr ' ' $'\n' <<< $flag_string | sort) <(tr ' ' $'\n' <<< $lock_string | sort)`
#
unmatched="$missing $new"
#
# Now new nor missing -> exit
#
if [[ -z $new ]] && [[ -z $missing ]]; then
 return
fi
#
#echo "D" $dir 
#echo "L" $lock_string
#echo "F" $flag_string
#echo "M" $missing 
#echo "N" $new 
#echo "SAVE" $save_dir 
#echo "RESTORE" $restore_dir 
#
# SAVE (step #1) & RESTORE (step #2) PJ dependent objects (from .dep files)
#
step=1
while [ $step -le 2 ]
do
 for lock in $unmatched
 do
  #
  if test -f "$dir/${lock}_project.dep"; then
   deps=`cat $dir/${lock}_project.dep`
   for file in $deps
   do
     if [ $step == 1 ]; then source ./sbin/compilation/object_save_restore_remove.sh "save"; fi
     if [ $step == 2 ]; then source ./sbin/compilation/object_save_restore_remove.sh "restore"; fi
   done
  fi
  #
  # Remove the lock 
  #
  if [ "$VERB" == 1 ] && [ $step == 2 ] ; then
   echo "rm -f $dir/$lock.lock"
  else
   rm -f $dir/$lock.lock
  fi
  DIR_is_to_recompile=1
  #
 done
 ((step++))
done
#
#
