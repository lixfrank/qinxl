if [ "$#" -ne 2 ] ; then
    export B0=/afs/cern.ch/user/l/llinwei/work2/qinxl/final/B0KstMuMu2
    source $B0/plugins/Scripts/InitAnalysis.sh $B0 0
    echo @@@ Error: Parameter missing @@@
    echo "Synopsis: InitAnalysis.sh analysis_path unset_display[0=false;1=true]"
else
    export ANALYPATH=$1
  
    export DATAYEAR=2017
    export DATADIR=/afs/cern.ch/user/l/llinwei/work2/qinxl/final/B0KstMuMu2/
    echo @@@ Analysis environment variable: $ANALYPATH @@@
    echo @@@ Directory with data: $DATADIR @@@
    echo @@@ Data year: $DATAYEAR @@@

    if [ "$2" -eq 1 ]; then
     unset DISPLAY
	echo @@@ Unset DISPLAY @@@
    fi
fi
