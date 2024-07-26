# CEPC_PN-ParT_procedure
from generator to Jet Origin Identification (JOI)

------

**Generator**
 - Download madgraph from [http://madgraph.phys.ucl.ac.be](http://madgraph.phys.ucl.ac.be) and install it. And you need to install Pythia, ExRootAnalysis, and HEPTools inside the madgraph.
 - To plot the JOI matrix <p align="center">
    <img src="figures/ConfusionMatrix.pdf.pdf" width="800"\>
</p> shown in [https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.132.221802](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.132.221802), you need to generate the samples of e+e- -> ZH (Z->vv,H->bb/cc/uu/dd/ss/gg) with the center of mass energy of 240 GeV. Since the Standard Model can generate H->ss/uu/dd, you need to import model HEL_UFO when generate samples with madgraph.

**Full Simulation**
 - If you can access the computing resource from Institute of High Energy Physics, Chinese Academy of Sciences, you can do full simulation with CEPC Software.
 - Welcome to join CEPC.

**Fast Simulation**
 - Instead of full simulation, which need intensive computing resources, you can do fast simulation with Delphes.
 - Download a special version of Delphes from [https://github.com/oiunun/Delphes_CEPC](https://github.com/oiunun/Delphes_CEPC), unpack it and make.
 - The Delphes card designed for 4th detector version of CEPC is delphes_card_CEPC_4th.tcl contained in Fast_simulation.
 - The directory of Fast_simulation has the following files: subjob.sh reads the file run_delphes.sh and set the input file (hepmc or stdhep) path, output file path and card path.

**Prepare traing/validation/testing samples**
 - If you choose full simulation, you can get the reconstructed fills including the information of jets and reconstructed particles. Then you can get the data features, such as the four momentum, impact parameters, particle PIDs of particle candidates within the jet, and stored these features into thr root file. The directory of full simulation has two files illustrating how to get these features.
 - If you choose fast simulation, the file makeNtuples.C contained in Fast_simulation directory illustrating how to get the data features used in Machine Learning (ML) model training.

**ML**


[[Paper]](https://arxiv.org/abs/1711.11586) [[Code]](implementations/bicyclegan/bicyclegan.py)

<p align="center">
    <img src="assets/bicyclegan_architecture.jpg" width="800"\>
</p>

**[New] Keras/TensorFlow implemetation** 
 - [model](tf-keras/tf_keras_model.py)
 - Requires tensorflow>=2.0.0 or >=1.15rc2. 
 - A full training example is available in [tf-keras/keras_train.ipynb](tf-keras/keras_train.ipynb). 
    - The top tagging dataset can be obtained from [https://zenodo.org/record/2603256](https://zenodo.org/record/2603256) and converted with this [script](tf-keras/convert_dataset.ipynb). 

## How to use the model

#### MXNet model

The ParticleNet model can be obtained by calling the `get_particle_net` function in [particle_net.py](mxnet/particle_net.py), which can return either an MXNet `Symbol` or an MXNet Gluon `HybridBlock`. The model takes three input arrays:
 - `points`: the coordinates of the particles in the (eta, phi) space. It should be an array with a shape of (N, 2, P), where N is the batch size and P is the number of particles.
 - `features`: the features of the particles. It should be an array with a shape of (N, C, P), where N is the batch size, C is the number of features, and P is the number of particles.
 - `mask`: a mask array with a shape of (N, 1, P), taking a value of 0 for padded positions.

To have a simple implementation for batched training on GPUs, we use fixed-length input arrays for all the inputs, although in principle the  ParticleNet architecture can handle variable number of particles in each jet. Zero-padding is used for the `points` and `features` inputs such that they always have the same length, and a `mask` array is used to indicate if a position is occupied by a real particle or by a zero-padded value.

The implementation of a simplified model, ParticleNet-Lite, is also provided and can be accessed with the `get_particle_net_lite` function.

#### Keras/TensorFlow model

The use of the Keras/TensorFlow model is similar to the MXNet model. A full training example is available in [tf-keras/keras_train.ipynb](tf-keras/keras_train.ipynb).

## Citation
If you use ParticleNet in your research, please cite the paper:

	@article{Qu:2019gqs,
	      author         = "Qu, Huilin and Gouskos, Loukas",
	      title          = "{ParticleNet: Jet Tagging via Particle Clouds}",
	      year           = "2019",
	      eprint         = "1902.08570",
	      archivePrefix  = "arXiv",
	      primaryClass   = "hep-ph",
	      SLACcitation   = "%%CITATION = ARXIV:1902.08570;%%"
	}

## Acknowledgement
The ParticleNet model is developed based on the [Dynamic Graph CNN](https://arxiv.org/abs/1801.07829) model. The implementation of the EdgeConv operation in MXNet is adapted from the author's TensorFlow [implementation](https://github.com/WangYueFt/dgcnn), and also inspired by the MXNet [implementation](https://github.com/chinakook/PointCNN.MX) of PointCNN.
