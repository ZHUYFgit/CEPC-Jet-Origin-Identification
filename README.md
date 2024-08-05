# CEPC-Jet-Origin-Identification
from generator to Jet Origin Identification (JOI)

------

## Generator
 - Download madgraph from [http://madgraph.phys.ucl.ac.be](http://madgraph.phys.ucl.ac.be) and install it. And you need to install Pythia, ExRootAnalysis, and hepmc inside the madgraph.
 - To plot the JOI matrix shown in [https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.132.221802](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.132.221802), you need to generate the samples of e+e- -> ZH (Z->vv,H->bb/cc/uu/dd/ss/gg) with the center of mass energy of 240 GeV. Since the Standard Model can not generate H->ss/uu/dd, you need to import model HEL_UFO when generate samples with madgraph.
After launch madgraph, type import model HEL_UFO, then madgraph would automatically download this model and give back a command, you need to type this command.

## Simulation
#### Full Simulation
 - If you can access the computing resource from Institute of High Energy Physics, Chinese Academy of Sciences, you can do full simulation with CEPC Software.
 - Welcome to join CEPC.
 - The directory [[full_simulation]](full_simulation) provides the code used to extract the features from the reconstructed files (with postfix slcio).

#### Fast Simulation
 - Instead of full simulation, which need intensive computing resources, you can do fast simulation with Delphes.
 - Download a special version of Delphes from [https://github.com/oiunun/Delphes_CEPC/releases/tag/v1.0](https://github.com/oiunun/Delphes_CEPC/releases/tag/v1.0), unpack it and make. 
 - The Delphes card designed for 4th detector version of CEPC is delphes_card_CEPC_4th.tcl contained in directory [[fast_simulation]](fast_simulation).
 - The directory of [[fast_simulation]](fast_simulation) has the following files.
   * setup.sh set the path of Delphes installed
   * subjob.sh reads the file run_delphes.sh and set the input file (hepmc or stdhep) path, output file path, and card path.
   * makeNtuples.C illustrats how to get the data features used in model training.


## Install Miniconda3, weaver, and ParticleNet
 - Install Miniconda3 according to your OS, such as you can install it with the following commands. You need to change the path in env_conda.sh to your installed miniconda3 path.
 ```
$ wget https://repo.anaconda.com/miniconda/Miniconda3-latest-lnux-x86_64.sh
$ chmod +x Miniconda3-latest-Linux-x86_64.sh
$ ./Miniconda3-latest-Linux-x86_64.sh
$ source env_conda.sh
```
 - Create a virtual environment, activate the created environment, install pytorch (according to ou OS/CUDA version at [https://pytorch.org/get-started](https://pytorch.org/get-started)) and weaver with the following commands. 
```
$ conda create -n weaver python=3.10
$ conda activate weaver
$ conda install pytorch==1.13.1 torchvision==0.14.1 torchaudio==0.13.1 pytorch-cuda=11.6 -c pytorch -c nvidia
$ pip install weaver-core
```
 - Once you do something wrong with weaver env, you can delete the env with the following command and recreate the env with the above commands.
```
$ conda env remove --name weaver
```
#### Install ParticleNet
 - Download ParticleNet and Particle Transformer from github  https://github.com/jet-universe/particle_transformer. Once your analysis use the code from ParticleNet or Particle Transformer, you need to cite the papers listed in https://github.com/jet-universe/particle_transformer.
 - The director of ParticleNet (suppose that the directory name you downloaded is ParticleNet) has several files.
   * ParticleNet/env.sh: set the input directories of samples in this file (export DATADIR_JetClass=)
   * ParticleNet/data/JetClass/JetClass_full.yaml:
     * *new_variables* means you can construct new variables based on the variables stored in your generated root files
     * *Pt_points* has two variables used to calculate the distance between two particles in ParticleNet
     * *pf_features* are the features used in training the model
     * *pf_vectors* are four momentum of particles used to calculate the pair-wise features used in Particle Transformer
     * *labels* list the labels of your sample when you want to train a classfication model
     * *observers* list the variables do not used to train the model while keep them in the files after testing
     * *length* restrict the number of particle candidates within the jet. In proton-proton collision, the particles are sorted by the transver momentum, while in electron-positron collision, the particles are sorted by the energy. If the *length* is larger than the number of particles in the jet, the leading *length* particles are preserved. If the *length* is smaller than the number of particles in the jet, the program would add particles with all features equal to 0.
   * ParticleNet/train_JetClass.sh: set the detailed input paths, predicted output path, and other hyper parameters   


#### ParticleNet Usage method
 - If you follow the above steps producing the fast simulation files and extract the data featurew with [[fast_simulation/makeNtuples.C]](fast_simulation/makeNtuples.C), you can put the [[training/JetClass_M11.yaml]](training/JetClass_M11.yaml) into your ParticleNet/data/JetClass and [[training/train_JetClass_M11.sh]](training/train_JetClass_M11.sh) into your directory ParticleNet.
 - Set the path and parameters in [[training/train_JetClass_M11.sh]](training/train_JetClass_M11.sh), then sh [[training/train_JetClass_M11.sh]](training/train_JetClass_M11.sh) PN full, where PN means the ParticleNet model and full means nothing in this example.
 - You also can use sh [[training/train_JetClass_M11.sh]](training/train_JetClass_M11.sh) ParT full, where ParT means Particle Transformer Model.



## Acknowledgement

We extend our heartfelt thanks to Huilin Qu and Congqiao Li for their invaluable support in utilizing ParticleNet and Particle Transformer. Our gratitude also goes to Shudong Wang and Xu Gao for their guidance in Delphes, as well as to Sitian Qian and Ze Guan for their support with Herwig.
