source /afs/ihep.ac.cn/users/z/zhuyf/cefs/workspace/Delphes/setup.sh
#source /cvmfs/sft.cern.ch/lcg/views/LCG_99/x86_64-centos7-gcc10-opt/setup.sh

DelphesINSTALL=/afs/ihep.ac.cn/users/z/zhuyf/cefs/workspace/Delphes/Delphes_CEPC-1.3
#DelphesINSTALL=/cefs/higgs/gaoxu/delphes/delphes

#temppath=`pwd`

$DelphesINSTALL/DelphesSTDHEP CARD OUTPUTFILE HEPFILE
#$DelphesINSTALL/DelphesHepMC2 CARD OUTPUTFILE HEPFILE
#$DelphesINSTALL/DelphesHepMC2 /afs/ihep.ac.cn/users/z/zhuyf/cefs/workspace/Delphes/yuexinCard/delphes_card_CEPC_4th_4jets.tcl /afs/ihep.ac.cn/users/z/zhuyf/cefs/workspace/Delphes/eezh.root /afs/ihep.ac.cn/users/z/zhuyf/cefs/workspace/madgraph/MG5_aMC_v2_9_18/eezh/Events/run_01/eezh.hepmc


#$DelphesINSTALL/DelphesSTDHEP $DelphesINSTALL/cards/delphes_card_CEPC_4th.tcl OUTPUTFILE HEPFILE
#echo $DelphesINSTALL/cards/delphes_card_CEPC_4th.tcl

unset DelphesINSTALL
