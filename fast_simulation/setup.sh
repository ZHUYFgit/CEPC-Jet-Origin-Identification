source /cvmfs/sft.cern.ch/lcg/views/setupViews.sh LCG_100 x86_64-centos7-gcc8-opt
#source /cvmfs/sft.cern.ch/lcg/views/setupViews.sh LCG_103 x86_64-centos7-gcc11-opt

# fastjetpath=/cefs/higgs/wangshudong/Delphes_CEPC/external/fastjet
exrootpath=/afs/ihep.ac.cn/users/z/zhuyf/cefs/workspace/Delphes/Delphes_CEPC-1.3/external/ExRootAnalysis
# export LD_LIBRARY_PATH=$fastjetpath/lib:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=$exrootpath:$LD_LIBRARY_PATH
#export LD_LIBRARY_PATH=/cvmfs/sft.cern.ch/lcg/views/LCG_103/x86_64-centos7-gcc11-opt/lib:$LD_LIBRARY_PATH
delphespath=/afs/ihep.ac.cn/users/z/zhuyf/cefs/workspace/Delphes/Delphes_CEPC-1.3
export LD_LIBRARY_PATH=$delphespath:$LD_LIBRARY_PATH


#EOF
