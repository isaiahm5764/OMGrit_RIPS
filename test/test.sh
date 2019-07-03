#!/bin/sh
#BHEADER**********************************************************************
#
# Copyright (c) 2013, Lawrence Livermore National Security, LLC. 
# Produced at the Lawrence Livermore National Laboratory. Written by 
# Jacob Schroder, Rob Falgout, Tzanio Kolev, Ulrike Yang, Veselin 
# Dobrev, et al. LLNL-CODE-660355. All rights reserved.
# 
# This file is part of XBraid. For support, post issues to the XBraid Github page.
# 
# This program is free software; you can redistribute it and/or modify it under
# the terms of the GNU General Public License (as published by the Free Software
# Foundation) version 2.1 dated February 1999.
# 
# This program is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the IMPLIED WARRANTY OF MERCHANTABILITY or FITNESS FOR A
# PARTICULAR PURPOSE. See the terms and conditions of the GNU General Public
# License for more details.
# 
# You should have received a copy of the GNU Lesser General Public License along
# with this program; if not, write to the Free Software Foundation, Inc., 59
# Temple Place, Suite 330, Boston, MA 02111-1307 USA
#
#EHEADER**********************************************************************


while [ "$*" ]
do
   case $1 in
      -h|-help)
         cat <<EOF

   $0 [options] {testname}.sh [{testname}.sh args]

    where: {testname} is the user-defined name for the test script

    with options:
       -h|-help       prints this usage information and exits
       -t|-trace      echo each command

   This script runs the Bourne shell test '{testname}.sh' and creates output
   files named '{testname}.err' and '{testname}.out' which capture the stderr
   and stdout output from the test.  The test script is run from the current
   directory, which should contain this script, '{test_name}.sh', and any other
   supporting files.

   A test is deemed to have passed when nothing is written to stderr.  A test
   may call other tests.  A test may take arguments, such as directories or
   files.  A test may also create output, which should be collected by the test
   in a directory named '{testname}.dir'.  A test may also require additional
   "filtering" in situations where information is erroneously written to stderr.
   Text identifying lines to be filtered are added to '{testname}.filters'.
   Usage documentation should appear at the top of each test.

   Example usage: $0 default.sh ..

EOF
         exit
         ;;
      -t|-trace)
         set -xv
         shift
         ;;
      *)
         break
         ;;
   esac
done

# Run the test and capture stdout, stderr
testname=`basename $1 .sh`
shift
echo "Running test [$testname]"
./$testname.sh $@ 1>"$testname.out" 2>"$testname.err"

