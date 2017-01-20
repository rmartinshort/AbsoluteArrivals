#!/bin/bash

# Prepare files before stacking
# Normalize the traces using macro normalize.m (Files *stk are overwritten)
# Calculate weighting factors used in noise assessment
# Save a phaseweighted stack
# Calculate the XC between the stack and each trace to check the alignment.

normmac=$1

sac << sacend
read *BHZ.p0
cuterr fillz 
cut t4 -10 10
read
synchronize
int
interpolate delta 0.025
taper width 0.3
rtrend
chnhdr b -10
write append .stk1
cut off
quit
sacend

#Do the normalization

for file in *BHZ.p0.stk1; do
echo "m $normmac "$file >> norm_trim.m
done
echo "quit" >> norm_trim.m

echo "m norm_trim.m" | sac


#calculate RMS value for 5s noise and signal windows
#do the stacking to form a stack and normalize it. Note that with 
#IRIS SAC we can only do a linear stack, but this probavbly won't make much of a
#difference to the final result - this should be checked

sac << sacend
read *BHZ.p0.stk1
mtw -7 -2
rms noise off to USER3
mtw 0 5
rms noise off to USER4
write append .rms

read *BHZ.p0.stk1
taper
sss
timewindow -10 10
sumstack normalization on
writestack stack.sac
quitsub
m $normmac stack.sac
quit
sacend

#Cross correlate all the trimmed files with the master stack
sac << sacend
read stack.sac *BHZ.p0.stk1
correlate master stack.sac normalized number 1 length 20 type rectangle
mtw -5 5
markptp length 10 to T8
write append .corr
quit
sacend

# Grab estimate for SNR for each trace
saclst user3 user4 f *.BHZ.p0.stk1.rms > rms_pre_arrival_signal.out # filename, rms of data in noise window, rms of data in signal window

# Grab correlation coefficients between stack and traces
saclst depmax t9 f *.BHZ.p0.stk1.corr | grep BHZ > correlation.out # filename, correlation co-efficient, XC location (correction)

