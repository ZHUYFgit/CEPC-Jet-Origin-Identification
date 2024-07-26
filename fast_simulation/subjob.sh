#!/bin/bash

path=`pwd`
cd $path

#htypes="X"
#ztypes="nn"
# ztypes="ll"
# ztypes="ww_l ww_sl ww_h zz_l zz_sl zz_h"
resos="cc"

export DelphesWorkDir=/afs/ihep.ac.cn/users/z/zhuyf/cefs/workspace/Delphes
for reso in ${resos}
do 
#    mkdir -p $DelphesWorkDir/Pythia6/${reso}
    OutputFileDir=$DelphesWorkDir/LLM/${reso}
    mkdir -p $OutputFileDir  
    tempcard=/afs/ihep.ac.cn/users/z/zhuyf/cefs/workspace/Delphes/yuexinCard/delphes_card_CEPC_4th.tcl

    ijob=1

    # INPUTDATADIR=/afs/ihep.ac.cn/users/z/zhuyf/cefs/workspace/madgraph/MG5_aMC_v2_9_18/${reso}/Events
    # Numjob=`ls -l ${INPUTDATADIR} | grep ^d | wc -l`
    # numJob=$((10#${Numjob}))

    # echo $numJob

    

    #while [ "$ijob" -le "$numJob" ]
    while [ "$ijob" -le "1000" ]
    do
        job_sn=`printf "%04d" $ijob`
        #InputFileDir=/cefs/data/stdhep/CEPC91/2fermions/wi_ISR_20220618_50M/2fermions/E91.2.Pbb.e0.p0.whizard195
        InputFileDir=/afs/ihep.ac.cn/users/z/zhuyf/cefs/workspace/whizard360/WhizardAis/data/higgs/E240.Pn2n2h_${reso}.e0.p0.whizard195
        #InputFileDir=${INPUTDATADIR}/run_${job_sn}
        # cd ${InputFileDir}
        # if [ -f "${InputFileDir}/tag_1_pythia8_events.hepmc.gz" ]; then
        #     gunzip ${InputFileDir}/tag_1_pythia8_events.hepmc.gz
        # fi
        # InputFiles_stdhep=/afs/ihep.ac.cn/users/z/zhuyf/cefs/workspace/madgraph/MG5_aMC_v2_9_18/${reso}/Events/run_${job_sn}/tag_1_pythia8_events.hepmc
        InputFiles_stdhep=${InputFileDir}/n2n2h_${reso}.e0.p0.${job_sn}.stdhep #bb.e0.p0.${job_sn}.stdhep
        OUTPUT_root=${job_sn}.root

        cp -f $DelphesWorkDir/run_delphes.sh $OutputFileDir/run_${reso}.e0.p0.${job_sn}.sh
        sed -i "s#OUTPUTFILE#${OutputFileDir}/${OUTPUT_root}#g" $OutputFileDir/run_${reso}.e0.p0.${job_sn}.sh
        sed -i "s#HEPFILE#${InputFiles_stdhep}#g" $OutputFileDir/run_${reso}.e0.p0.${job_sn}.sh
        sed -i "s#CARD#${tempcard}#g" $OutputFileDir/run_${reso}.e0.p0.${job_sn}.sh
        chmod +x $OutputFileDir/run_${reso}.e0.p0.${job_sn}.sh
#           source /afs/ihep.ac.cn/users/z/zhuyf/cefs/workspace/Delphes/setup.sh
        hep_sub  $OutputFileDir/run_${reso}.e0.p0.${job_sn}.sh
        #sh $OutputFileDir/run_${reso}.e0.p0.${job_sn}.sh


    let "ijob+=1"
    done
done
