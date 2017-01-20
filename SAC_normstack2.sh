#!/bin/bash

normmac=$1

#Collect the files ready for stack and add them to the stack
#restack the chosen seismograms according to their assigned weights
#determine the SNR for the new stack (stack2.sac)


#make all the .stk2 files
sac << sacend
m mk_stk2.m
quit
sacend

#This block may take a long time to complete!!
echo "About to stack round 2"
sac << sacend
sss
m addstack.m
timewindow -10 10
sumstack normalization on
writestack stack2.sac
quitsub

m $normmac stack2.sac
quit
sacend
echo "completed stack round 2"

#Find the SNR of the stack and write to USER0 and USER5
sac << sacend
read stack2.sac
chnhdr b -10
wh
mtw -7 -2
rms noise off to USER0
mtw 0 5
rms noise off to USER5
write over
quit
sacend


#Compute correlations between all the files and the second stack
sac << sacend
read stack2.sac *stk2
correlate master stack2.sac normalized number 1 length 20 type rectangle
mtw -5 5
markptp length 10 to T8
write append .corr2
quit
sacend

#During the above process, we've correlated the stack with itself
cp stack2.sac.corr2 stack2_auto.sac
saclst depmax t9 f *.corr2 > correlation2.out  # filename, correlation co-efficient, XC location (correction)
