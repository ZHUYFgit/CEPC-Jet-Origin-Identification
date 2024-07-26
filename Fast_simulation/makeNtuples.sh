
path=`pwd`
dirs="bb cc uu dd ss"

for dir in $dirs
do
    OUTPUTDATADIR=$path/rootFiles/Herwig_test_onnx/${dir}bar
    #rm -rf $OUTPUTDATADIR
    mkdir -p $OUTPUTDATADIR

    INPUTDATADIR=$path/Herwig/${dir}

    Numjob=`ls ${INPUTDATADIR}/*.root |wc -l`
    numJob=$((10#${Numjob}))
    #numJob=280
    echo ${numJob}

    j=1
    while [ "$j" -le "${numJob}" ]
    #while [ "$j" -lt "2" ]
    do
        #job_sn=`printf "%04d" $j`
        job_sn=$j
        infile=$INPUTDATADIR/${job_sn}.root
        outfile=$OUTPUTDATADIR/${dir}bar_${j}.root
        filename=makeNtuplesZH_${dir}_${j}
        cp -f $path/makeNtuples.C $OUTPUTDATADIR/${filename}.C
        sed -i "s#INFILE#${infile}#g" $OUTPUTDATADIR/${filename}.C
        sed -i "s#OUTFILE#${outfile}#g" $OUTPUTDATADIR/${filename}.C
        sed -i "s#FILENAME#${filename}#g" $OUTPUTDATADIR/${filename}.C
        sed -i "s#RUNCLASS#${dir}#g" $OUTPUTDATADIR/${filename}.C

        echo \
            "#! /bin/bash
            source /cvmfs/sft.cern.ch/lcg/views/setupViews.sh LCG_100 x86_64-centos7-gcc8-opt
  
            # source /afs/ihep.ac.cn/users/z/zhuyf/cefs/workspace/GraphJetID/setup.sh
            # export PATH=/afs/ihep.ac.cn/users/z/zhuyf/cefs/workspace/GraphJetID/build/lib/libGraphJetID.so:\$PATH

            root -l $OUTPUTDATADIR/${filename}.C
            #> $OUTPUTDATADIR/makeNtuplesZH_${dir}_${j}.log

            " > $OUTPUTDATADIR/makeNtuplesZH_${dir}_${j}.sh
        #echo "$RecoWorkDir/$OUTPUTDATA/reco_${ojob}.sh"
        #export PATH=/afs/ihep.ac.cn/soft/common/sysgroup/hep_job/bin:$PATH
            chmod +x $OUTPUTDATADIR/makeNtuplesZH_${dir}_${j}.sh


        #hep_sub $OUTPUTDATADIR/makeNtuplesZH_${dir}_${j}.sh
        sh $OUTPUTDATADIR/makeNtuplesZH_${dir}_${j}.sh
        let "j+=1"
    done
done
