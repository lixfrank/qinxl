ROOT version: 6.04
gcc version: 4.9.3
we can use that version with the two different efficiency from KDE.



Maximum likelihood fit:
1. setup environment
   under the directory "plugins/Scripts/InitAnalysis.sh"  l2 & l11
2. change the directory of reading efficiency
   under the directory "plugins/ExtractYield.c" l779
3. choose type you want to fit
   under the directory "/python/ParameterFile.txt" 
4. compile the fitter program: make clean, make,  make ExtractYield
5. run the fit
   eg. ./ExtractYield fit-type sample-directory q^2bin
   (106: Reco MC, 206:GEN MC) 

Method of Moment:
1. change the directory of reading efficiency
   plugins/Moment.cc l371
2. change the directory of sample
   plugins/Moment.cc l204(gen) l384(reco)
3. compile the program: make Moment
4. run 
   if gen-level: ./Moment Paramater-type gen
   if reco-level: ./Moment Paramater-type reco event-tag(mis or good)



