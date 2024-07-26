####################################################################                                l
# CEPC 4th detector model
#
# Authors:
# email:
#####################################################################

## MOD2: set vtx mode timing to MC truth

set B 3.0

## Drift chamber coordinates
set DCHZMIN -2.98
set DCHZMAX  2.98
set DCHRMIN  0.6
set DCHRMAX  1.8

#######################################
# Order of execution of various modules
#######################################

set ExecutionPath {

  TruthVertexFinder
  ParticlePropagator

  ChargedHadronTrackingEfficiency
  ElectronTrackingEfficiency
  MuonTrackingEfficiency

  TrackMergerPre
  TrackSmearing
  ClusterCounting
  TimeSmearing
  TimeOfFlight

  TrackMerger
  IdentificationMap
  Calorimeter

  TimeSmearingNeutrals
  TimeOfFlightNeutralHadron
  EFlowMerger
  PhotonFilter
  PhotonEfficiency
  PhotonIsolation


  ElectronFilter
  ElectronEfficiency
  ElectronIsolation

  MuonFilter
  MuonEfficiency

  MissingET

  NeutrinoFilter
  GenJetFinder
  GenMissingET

  FastJetFinder

  JetEnergyScale

  JetFlavorAssociation

  BTagging
  TauTagging

  TreeWriter
}

#################################
# Truth Vertex Finder
#################################

module TruthVertexFinder TruthVertexFinder {

  ## below this distance two vertices are assumed to be merged
  set Resolution 1E-06
  set InputArray Delphes/stableParticles
  set VertexOutputArray vertices
}

#################################
# Propagate particles in cylinder
#################################

module ParticlePropagator ParticlePropagator {
  set InputArray Delphes/stableParticles
  set OutputArray stableParticles
  set ChargedHadronOutputArray chargedHadrons
  set ElectronOutputArray electrons
  set MuonOutputArray muons

  # inner radius of the solenoid, in m
  set Radius     2.25
  # half-length: z of the solenoid, in m
  set HalfLength 2.50
  # magnetic field, in T
  set Bz         $B
}

####################################
# Charged hadron tracking efficiency
####################################

module Efficiency ChargedHadronTrackingEfficiency {
    set InputArray ParticlePropagator/chargedHadrons
    set OutputArray chargedHadrons
    set UseMomentumVector true

    # We use only one efficiency, we set only 0 effincency out of eta bounds:

    set EfficiencyFormula {
        (abs(eta) > 3.0)                               * (0.000) +
        (pt >= 0.5) * (abs(eta) <= 3.0)                * (0.997) +
        (pt <  0.5 && pt >= 0.1) * (abs(eta) <= 3.0)   * (0.900) +
        (pt <  0.1) * (abs(eta) <= 3.0)                * (0.060)
    }

}

##############################
# Electron tracking efficiency
##############################

module Efficiency ElectronTrackingEfficiency {
    set InputArray ParticlePropagator/electrons
    set OutputArray electrons
    set UseMomentumVector true

    # Current full simulation with CLICdet provides for electrons:
    set EfficiencyFormula {
        (abs(eta) > 3.0)                               * (0.000) +
        (pt >= 0.5) * (abs(eta) <= 3.0)                * (0.997) +
        (pt <  0.5 && pt >= 0.1) * (abs(eta) <= 3.0)   * (0.900) +
        (pt <  0.1) * (abs(eta) <= 3.0)                * (0.060)
    }

}


##########################
# Muon tracking efficiency
##########################

module Efficiency MuonTrackingEfficiency {
    set InputArray ParticlePropagator/muons
    set OutputArray muons
    set UseMomentumVector true

    # Current full simulation with CLICdet provides for muons:
    set EfficiencyFormula {
        (abs(eta) > 3.0)                               * (0.000) +
        (pt >= 0.5) * (abs(eta) <= 3.0)                * (0.997) +
        (pt < 0.5 && pt >= 0.1) * (abs(eta) <= 3.0)    * (0.900) +
        (pt < 0.1) * (abs(eta) <= 3.0)                 * (0.060)
    }
}
##############
# Track merger
##############

module Merger TrackMergerPre {
# add InputArray InputArray
  add InputArray  ChargedHadronTrackingEfficiency/chargedHadrons
  add InputArray  ElectronTrackingEfficiency/electrons
  add InputArray  MuonTrackingEfficiency/muons
  set OutputArray tracks
}


########################################
# Smearing for charged tracks
########################################

module TrackCovariance TrackSmearing {

    set InputArray  TrackMergerPre/tracks
    set OutputArray tracks

    ## minimum number of hits to accept a track
    set NMinHits 6

    ## magnetic field
    set Bz $B

    ## uses https://raw.githubusercontent.com/selvaggi/FastTrackCovariance/master/GeoIDEA_BASE.txt
    set DetectorGeometry {

      # Layer type 1 = R (barrel) or 2 = z (forward/backward)
      # Layer label
      # Minimum dimension z for barrel or R for forward
      # Maximum dimension z for barrel or R for forward
      # R/z location of layer
      # Thickness (meters)
      # Radiation length (meters)
      # Number of measurements in layers (1D or 2D)
      # Stereo angle (rad) - 0(pi/2) = axial(z) layer - Upper side
      # Stereo angle (rad) - 0(pi/2) = axial(z) layer - Lower side
      # Resolution Upper side (meters) - 0 = no measurement
      # Resolution Lower side (meters) - 0 = no measurement
      # measurement flag = T, scattering only = F

      # barrel  name     zmin     zmax   r        w (m)      X0        n_meas  th_up (rad) th_down (rad)    reso_up (m)   reso_down (m)  flag
      1        PIPE     -3.0    3.0    0.014    0.00014   0.0937     0        0          0                0             0              0
      # Original
      # 1        VTX1A    -0.2      0.2      0.016   0.00014   0.0937     2       0          1.5708          3e-006        3e-006         1
      # 1        VTX1B    -0.2      0.2      0.018   0.00014   0.0937     2       0          1.5708          6e-006        6e-006         1
      # 1        VTX2A    -0.2      0.2      0.038   0.00014   0.0937     2       0          1.5708          4e-006        4e-006         1
      # 1        VTX2B    -0.2      0.2      0.040   0.00014   0.0937     2       0          1.5708          4e-006        4e-006         1
      # 1        VTX3A    -0.2      0.2      0.058   0.00014   0.0937     2       0          1.5708          4e-006        4e-006         1
      # 1        VTX3B    -0.2      0.2      0.060   0.00014   0.0937     2       0          1.5708          4e-006        4e-006         1

      # Baseline
      # 1        VTX1A    -0.03125  0.03125  0.016   0.00014   0.0937     2       0.2456     1.5708          2.8e-006      2.8e-006       1
      # 1        VTX1B    -0.03125  0.03125  0.018   0.00014   0.0937     2       0.2838     1.5708          6e-006        6e-006         1
      # 1        VTX2A    -0.0625   0.0625   0.037   0.00014   0.0937     2       0.2838     1.5708          4e-006        4e-006         1
      # 1        VTX2B    -0.0625   0.0625   0.039   0.00014   0.0937     2       0.3176     1.5708          4e-006        4e-006         1
      # 1        VTX3A    -0.0625   0.0625   0.058   0.00014   0.0937     2       0.4275     1.5708          4e-006        4e-006         1
      # 1        VTX3B    -0.0625   0.0625   0.060   0.00014   0.0937     2       0.4510     1.5708          4e-006        4e-006         1

      # Baseline but full acceptance
       1        VTX1A    -0.03125  0.03125  0.016   0.00014   0.0937     2       0     1.5708          2.8e-006      2.8e-006       1
       1        VTX1B    -0.03125  0.03125  0.018   0.00014   0.0937     2       0     1.5708          6e-006        6e-006         1
       1        VTX2A    -0.0625   0.0625   0.037   0.00014   0.0937     2       0     1.5708          4e-006        4e-006         1
       1        VTX2B    -0.0625   0.0625   0.039   0.00014   0.0937     2       0     1.5708          4e-006        4e-006         1
       1        VTX3A    -0.0625   0.0625   0.058   0.00014   0.0937     2       0     1.5708          4e-006        4e-006         1
       1        VTX3B    -0.0625   0.0625   0.060   0.00014   0.0937     2       0     1.5708          4e-006        4e-006         1

      # Test
      # 1        VTX1A    -0.3   0.3   0.016   0.00014   0.0937     2       0     1.5708          5e-006        5e-006         1
      # 1        VTX1B    -0.3   0.3   0.018   0.00014   0.0937     2       0     1.5708          6e-006        6e-006         1
      # 1        VTX2A    -0.3   0.3   0.037   0.00014   0.0937     2       0     1.5708          4e-006        4e-006         1
      # 1        VTX2B    -0.3   0.3   0.039   0.00014   0.0937     2       0     1.5708          4e-006        4e-006         1
      # 1        VTX3A    -0.3   0.3   0.058   0.00014   0.0937     2       0     1.5708          4e-006        4e-006         1
      # 1        VTX3B    -0.3   0.3   0.060   0.00014   0.0937     2       0     1.5708          4e-006        4e-006         1

      1        SHELL    -0.2   0.2   0.065    0.00014   0.0937     0        0          0                0             0              0
      1        SIT01    -0.241  0.241   0.12    0.00061   0.0937     2        0          1.5708           7e-006        86e-006        1
      1        SIT02    -0.455  0.455   0.27     0.00061   0.0937     2        0          1.5708           7e-006        86e-006        1
      1        SIT03    -0.721  0.721   0.42    0.00061   0.0937     2        0          1.5708           7e-006        86e-006        1
      1        SIT04    -0.988  0.988   0.57    0.00061   0.0937     2        0          1.5708           7e-006        86e-006        1

    # endcap  name     rmin   rmax       z       w (m)      X0       n_meas   th_up (rad)  th_down (rad)   reso_up (m)   reso_down (m) flag
      2        DSK1A     0.0295     0.12    0.241      0.00061   0.0937     2        0          1.5708           7e-006        86e-006        1
      2        DSK1B     0.0295     0.12    -0.241     0.00061   0.0937     2        0          1.5708           7e-006        86e-006        1
      2        DSK2A     0.0305     0.27    0.455    0.00061   0.0937     2        0          1.5708           7e-006        86e-006        1
      2        DSK2B     0.0305     0.27    -0.455   0.00061   0.0937     2        0          1.5708           7e-006        86e-006        1
      2        DSK3A     0.0325     0.42    0.721    0.00061   0.0937     2        0          1.5708           7e-006        86e-006        1
      2        DSK3B     0.0325     0.42    -0.721   0.00061   0.0937     2        0          1.5708           7e-006        86e-006        1
      2        DSK4A     0.0340     0.57     0.988     0.00061   0.0937     2        0          1.5708           7e-006        86e-006        1
      2        DSK4B     0.0340     0.57    -0.988    0.00061   0.0937     2        0          1.5708           7e-006        86e-006        1
      # Drift chamber
      1 DCHCANI  $DCHZMIN   $DCHZMAX    0.60478   0.00002   0.2797        0  0    0  0     0.0 0
      1 DCH       -2.98       2.98       0.605     0.015     $RLENGTH      1  0.05 0 0.0001 0 1
      1 DCH       -2.98       2.98       0.62     0.015     $RLENGTH      1  0.05 0 0.0001 0 1
      1 DCH       -2.98       2.98       0.635     0.015     $RLENGTH      1  0.05 0 0.0001 0 1
      1 DCH       -2.98       2.98       0.65     0.015     $RLENGTH      1  0.05 0 0.0001 0 1
      1 DCH       -2.98       2.98       0.665     0.015     $RLENGTH      1  0.05 0 0.0001 0 1
      1 DCH       -2.98       2.98       0.68     0.015     $RLENGTH      1  0.05 0 0.0001 0 1
      1 DCH       -2.98       2.98       0.695     0.015     $RLENGTH      1  0.05 0 0.0001 0 1
      1 DCH       -2.98       2.98       0.71     0.015     $RLENGTH      1  0.05 0 0.0001 0 1
      1 DCH       -2.98       2.98       0.725     0.015     $RLENGTH      1  0.05 0 0.0001 0 1
      1 DCH       -2.98       2.98       0.74     0.015     $RLENGTH      1  0.05 0 0.0001 0 1
      1 DCH       -2.98       2.98       0.755     0.015     $RLENGTH      1  0.05 0 0.0001 0 1
      1 DCH       -2.98       2.98       0.77     0.015     $RLENGTH      1  0.05 0 0.0001 0 1
      1 DCH       -2.98       2.98       0.785     0.015     $RLENGTH      1  0.05 0 0.0001 0 1
      1 DCH       -2.98       2.98       0.8     0.015     $RLENGTH      1  0.05 0 0.0001 0 1
      1 DCH       -2.98       2.98       0.815     0.015     $RLENGTH      1 -0.05 0 0.0001 0 1
      1 DCH       -2.98       2.98       0.83     0.015     $RLENGTH      1  0.05 0 0.0001 0 1
      1 DCH       -2.98       2.98       0.845     0.015     $RLENGTH      1 -0.05 0 0.0001 0 1
      1 DCH       -2.98       2.98       0.86     0.015     $RLENGTH      1  0.05 0 0.0001 0 1
      1 DCH       -2.98       2.98       0.875     0.015     $RLENGTH      1 -0.05 0 0.0001 0 1
      1 DCH       -2.98       2.98       0.89     0.015     $RLENGTH      1  0.05 0 0.0001 0 1
      1 DCH       -2.98       2.98       0.905     0.015     $RLENGTH      1 -0.05 0 0.0001 0 1
      1 DCH       -2.98       2.98       0.92     0.015     $RLENGTH      1  0.05 0 0.0001 0 1
      1 DCH       -2.98       2.98       0.935     0.015     $RLENGTH      1 -0.05 0 0.0001 0 1
      1 DCH       -2.98       2.98       0.95     0.015     $RLENGTH      1  0.05 0 0.0001 0 1
      1 DCH       -2.98       2.98       0.965     0.015     $RLENGTH      1 -0.05 0 0.0001 0 1
      1 DCH       -2.98       2.98       0.98     0.015     $RLENGTH      1  0.05 0 0.0001 0 1
      1 DCH       -2.98       2.98       0.995     0.015     $RLENGTH      1 -0.05 0 0.0001 0 1
      1 DCH       -2.98       2.98       1.01     0.015     $RLENGTH      1  0.05 0 0.0001 0 1
      1 DCH       -2.98       2.98       1.025     0.015     $RLENGTH      1 -0.05 0 0.0001 0 1
      1 DCH       -2.98       2.98       1.04     0.015     $RLENGTH      1  0.05 0 0.0001 0 1
      1 DCH       -2.98       2.98       1.055     0.015     $RLENGTH      1 -0.05 0 0.0001 0 1
      1 DCH       -2.98       2.98       1.07     0.015     $RLENGTH      1  0.05 0 0.0001 0 1
      1 DCH       -2.98       2.98       1.085     0.015     $RLENGTH      1 -0.05 0 0.0001 0 1
      1 DCH       -2.98       2.98       1.1     0.015     $RLENGTH      1  0.05 0 0.0001 0 1
      1 DCH       -2.98       2.98       1.115     0.015     $RLENGTH      1 -0.05 0 0.0001 0 1
      1 DCH       -2.98       2.98       1.13     0.015     $RLENGTH      1  0.05 0 0.0001 0 1
      1 DCH       -2.98       2.98       1.145     0.015     $RLENGTH      1 -0.05 0 0.0001 0 1
      1 DCH       -2.98       2.98       1.16     0.015     $RLENGTH      1  0.05 0 0.0001 0 1
      1 DCH       -2.98       2.98       1.175     0.015     $RLENGTH      1 -0.05 0 0.0001 0 1
      1 DCH       -2.98       2.98       1.19     0.015     $RLENGTH      1  0.05 0 0.0001 0 1
      1 DCH       -2.98       2.98       1.205     0.015     $RLENGTH      1 -0.05 0 0.0001 0 1
      1 DCH       -2.98       2.98       1.22     0.015     $RLENGTH      1  0.05 0 0.0001 0 1
      1 DCH       -2.98       2.98       1.235     0.015     $RLENGTH      1 -0.05 0 0.0001 0 1
      1 DCH       -2.98       2.98       1.25     0.015     $RLENGTH      1  0.05 0 0.0001 0 1
      1 DCH       -2.98       2.98       1.265     0.015     $RLENGTH      1 -0.05 0 0.0001 0 1
      1 DCH       -2.98       2.98       1.28     0.015     $RLENGTH      1  0.05 0 0.0001 0 1
      1 DCH       -2.98       2.98       1.295     0.015     $RLENGTH      1 -0.05 0 0.0001 0 1
      1 DCH       -2.98       2.98       1.31     0.015     $RLENGTH      1  0.05 0 0.0001 0 1
      1 DCH       -2.98       2.98       1.325     0.015     $RLENGTH      1 -0.05 0 0.0001 0 1
      1 DCH       -2.98       2.98       1.34     0.015     $RLENGTH      1  0.05 0 0.0001 0 1
      1 DCH       -2.98       2.98       1.355     0.015     $RLENGTH      1 -0.05 0 0.0001 0 1
      1 DCH       -2.98       2.98       1.37     0.015     $RLENGTH      1  0.05 0 0.0001 0 1
      1 DCH       -2.98       2.98       1.385     0.015     $RLENGTH      1 -0.05 0 0.0001 0 1
      1 DCH       -2.98       2.98       1.4     0.015     $RLENGTH      1  0.05 0 0.0001 0 1
      1 DCH       -2.98       2.98       1.415     0.015     $RLENGTH      1 -0.05 0 0.0001 0 1
      1 DCH       -2.98       2.98       1.43     0.015     $RLENGTH      1  0.05 0 0.0001 0 1
      1 DCH       -2.98       2.98       1.445     0.015     $RLENGTH      1 -0.05 0 0.0001 0 1
      1 DCH       -2.98       2.98       1.46     0.015     $RLENGTH      1  0.05 0 0.0001 0 1
      1 DCH       -2.98       2.98       1.475     0.015     $RLENGTH      1 -0.05 0 0.0001 0 1
      1 DCH       -2.98       2.98       1.49     0.015     $RLENGTH      1  0.05 0 0.0001 0 1
      1 DCH       -2.98       2.98       1.505     0.015     $RLENGTH      1 -0.05 0 0.0001 0 1
      1 DCH       -2.98       2.98       1.52     0.015     $RLENGTH      1  0.05 0 0.0001 0 1
      1 DCH       -2.98       2.98       1.535     0.015     $RLENGTH      1 -0.05 0 0.0001 0 1
      1 DCH       -2.98       2.98       1.55     0.015     $RLENGTH      1  0.05 0 0.0001 0 1
      1 DCH       -2.98       2.98       1.565     0.015     $RLENGTH      1 -0.05 0 0.0001 0 1
      1 DCH       -2.98       2.98       1.58     0.015     $RLENGTH      1  0.05 0 0.0001 0 1
      1 DCH       -2.98       2.98       1.595     0.015     $RLENGTH      1 -0.05 0 0.0001 0 1
      1 DCH       -2.98       2.98       1.61     0.015     $RLENGTH      1  0.05 0 0.0001 0 1
      1 DCH       -2.98       2.98       1.625     0.015     $RLENGTH      1 -0.05 0 0.0001 0 1
      1 DCH       -2.98       2.98       1.64     0.015     $RLENGTH      1  0.05 0 0.0001 0 1
      1 DCH       -2.98       2.98       1.655     0.015     $RLENGTH      1 -0.05 0 0.0001 0 1
      1 DCH       -2.98       2.98       1.67     0.015     $RLENGTH      1  0.05 0 0.0001 0 1
      1 DCH       -2.98       2.98       1.685     0.015     $RLENGTH      1 -0.05 0 0.0001 0 1
      1 DCH       -2.98       2.98       1.7     0.015     $RLENGTH      1  0.05 0 0.0001 0 1
      1 DCH       -2.98       2.98       1.715     0.015     $RLENGTH      1 -0.05 0 0.0001 0 1
      1 DCH       -2.98       2.98       1.73     0.015     $RLENGTH      1  0.05 0 0.0001 0 1
      1 DCH       -2.98       2.98       1.745     0.015     $RLENGTH      1 -0.05 0 0.0001 0 1
      1 DCH       -2.98       2.98       1.76     0.015     $RLENGTH      1  0.05 0 0.0001 0 1
      1 DCH       -2.98       2.98       1.775     0.015     $RLENGTH      1 -0.05 0 0.0001 0 1
      1 DCH       -2.98       2.98       1.79     0.01     $RLENGTH      1  0.05 0 0.0001 0 1

      1 DCHCANO $DCHZMIN  $DCHZMAX   1.80002      0.0028  0.2797         0  0    0 0      0.0 0
      # SET & ETD
      1 SET       -2.98      2.98        1.815     0.0061   0.0937         2  0  1.5708  7e-006 86e-006       1
      2 ETD        0.6       1.822       3.0      0.0061   0.0937         2  0  1.5708  7e-006 86e-006       1
      2 ETD        0.6       1.822      -3.0      0.0061   0.0937         2  0  1.5708  7e-006 86e-006       1

    }
}

###################
# Cluster Counting
###################

module ClusterCounting ClusterCounting {

  add InputArray TrackSmearing/tracks
  set OutputArray tracks

  set Bz $B

  ## check that these are consistent with DCHCANI/DCHNANO parameters in TrackCovariance module
  set Rmin $DCHRMIN
  set Rmax $DCHRMAX
  set Zmin $DCHZMIN
  set Zmax $DCHZMAX

  # gas mix option:
  # 0:  Helium 90% - Isobutane 10%
  # 1:  Helium 100%
  # 2:  Argon 50% - Ethane 50%
  # 3:  Argon 100%

  set GasOption 0

}


########################################
#   Time Smearing MIP
########################################

module TimeSmearing TimeSmearing {
  set InputArray ClusterCounting/tracks
  set OutputArray tracks

  # assume constant 30 ps resolution for now
  set TimeResolution {
                       (abs(eta) > 0.0 && abs(eta) <= 3.0)* 30E-12
                     }
}


########################################
#   Time Of Flight Measurement
########################################

module TimeOfFlight TimeOfFlight {
  set InputArray TimeSmearing/tracks
  set VertexInputArray TruthVertexFinder/vertices

  set OutputArray tracks

  # 0: assume vertex time tV from MC Truth (ideal case)
  # 1: assume vertex time tV = 0
  # 2: calculate vertex time as vertex TOF, assuming tPV=0

  set VertexTimeMode 0
}

##############
# Track merger
##############

module Merger TrackMerger {
# add InputArray InputArray
  add InputArray  TimeOfFlight/tracks
  set OutputArray tracks
}


#############
# Calorimeter
#############
module DualReadoutCalorimeter Calorimeter {
  set ParticleInputArray ParticlePropagator/stableParticles
  set TrackInputArray    IdentificationMap/tracks

  set TowerOutputArray   towers
  set PhotonOutputArray  photons

  set EFlowTrackOutputArray         eflowTracks
  set EFlowPhotonOutputArray        eflowPhotons
  set EFlowNeutralHadronOutputArray eflowNeutralHadrons

  set ECalEnergyMin 0.5
  set HCalEnergyMin 0.5
  set EnergyMin     0.5
  set ECalEnergySignificanceMin 1.0
  set HCalEnergySignificanceMin 1.0
  set EnergySignificanceMin     1.0

  set SmearTowerCenter true
  #set SmearTowerCenter false
  set pi [expr {acos(-1)}]

    # Lists of the edges of each tower in eta and phi;
    # each list starts with the lower edge of the first tower;
    # the list ends with the higher edged of the last tower.
    # Barrel:  deta=0.01 towers up to |eta| <= 0.88 ( up to 45°)
    # Endcaps: deta=0.01 towers up to |eta| <= 3.0 (8.6° = 100 mrad)
    # Cell size: about 2 cm x 3 cm

    #barrel:
    set PhiBins {}
    for {set i -120} {$i <= 120} {incr i} {
        add PhiBins [expr {$i * $pi/120}]
    }
    #deta=0.01 units for |eta| <= 0.88
    for {set i -44} {$i < 45} {incr i} {
        set eta [expr {$i * 0.02}]
        add EtaPhiBins $eta $PhiBins
    }

    #endcaps:
    set PhiBins {}
    for {set i -120} {$i <= 120} {incr i} {
        add PhiBins [expr {$i* $pi/120}]
    }
    #deta=0.01 units for 0.88 < |eta| <= 3.0
    #first, from -3.0 to -0.88
    for {set i 0} {$i <106} {incr i} {
        set eta [expr {-3.00 + $i * 0.02}]
        add EtaPhiBins $eta $PhiBins
    }
    #same for 0.88 to 3.0
    for  {set i 1} {$i <=106} {incr i} {
        set eta [expr {0.88 + $i * 0.02}]
        add EtaPhiBins $eta $PhiBins
    }

    # default energy fractions {abs(PDG code)} {Fecal Fhcal}
    add EnergyFraction {0} {0.0 1.0}
    # energy fractions for e, gamma and pi0
    add EnergyFraction {11} {1.0 0.0}
    add EnergyFraction {22} {1.0 0.0}
    add EnergyFraction {111} {1.0 0.0}
    # energy fractions for muon, neutrinos and neutralinos
    add EnergyFraction {12} {0.0 0.0}
    add EnergyFraction {13} {0.1 0.2}
    add EnergyFraction {14} {0.0 0.0}
    add EnergyFraction {16} {0.0 0.0}
    add EnergyFraction {1000022} {0.0 0.0}
    add EnergyFraction {1000023} {0.0 0.0}
    add EnergyFraction {1000025} {0.0 0.0}
    add EnergyFraction {1000035} {0.0 0.0}
    add EnergyFraction {1000045} {0.0 0.0}
    # energy fractions for K0short and Lambda
    add EnergyFraction {310} {0.3 0.7}
    add EnergyFraction {130} {0.3 0.7}
    add EnergyFraction {3122} {0.3 0.7}
    # energy fractions for pi
    add EnergyFraction {211} {0.3 0.7}
    add EnergyFraction {321} {0.3 0.7}
    add EnergyFraction {2212} {0.3 0.7}
    # set ECalResolutionFormula {resolution formula as a function of eta and energy}
    set ECalResolutionFormula {
    (abs(eta) <= 0.88 )                     * sqrt(energy^2*0.01^2 + energy*0.03^2)+
    (abs(eta) > 0.88 && abs(eta) <= 3.0)    * sqrt(energy^2*0.01^2 + energy*0.03^2)
    }

    # set HCalResolutionFormula {resolution formula as a function of eta and energy}
    set HCalResolutionFormula {
    (abs(eta) <= 0.88 )                     * sqrt(energy^2*0.02^2 + energy*0.40^2)+
    (abs(eta) > 0.88 && abs(eta) <= 3.0)    * sqrt(energy^2*0.02^2 + energy*0.40^2)
    }
}




####################################
# Charged hadron PID
####################################

module IdentificationMap IdentificationMap {
  set InputArray TrackMerger/tracks
  set OutputArray tracks

}
########################################
#   Time Smearing Neutrals
########################################

module TimeSmearing TimeSmearingNeutrals {
  set InputArray  Calorimeter/eflowNeutralHadrons
  set OutputArray eflowNeutralHadrons

  # assume constant 30 ps resolution for now
  set TimeResolution {
                       (abs(eta) > 0.0 && abs(eta) <= 3.0)* 30E-12
                     }
}

########################################
#   Time Of Flight Measurement
########################################

module TimeOfFlight TimeOfFlightNeutralHadron {
  set InputArray       TimeSmearingNeutrals/eflowNeutralHadrons
  set VertexInputArray TruthVertexFinder/vertices
  set OutputArray      eflowNeutralHadrons

  # 0: assume vertex time tV from MC Truth (ideal case)
  # 1: assume vertex time tV = 0
  # 2: calculate vertex time as vertex TOF, assuming tPV=0

  ## TBF (add option to take hard vertex time)
  set VertexTimeMode 1
}

####################
# Energy flow merger
####################

module Merger EFlowMerger {
  add InputArray  Calorimeter/eflowTracks
  add InputArray  Calorimeter/eflowPhotons
  add InputArray  TimeOfFlightNeutralHadron/eflowNeutralHadrons
  set OutputArray eflow
}


###################
# Photon efficiency
###################

module Efficiency PhotonEfficiency {
  set InputArray  Calorimeter/eflowPhotons
  set OutputArray photons

  # set EfficiencyFormula {efficiency formula as a function of eta and pt}
  # efficiency formula for photons
  set EfficiencyFormula {
        (energy <  2.0)                                       * (0.000) +
        (energy >= 2.0) * (abs(eta) <= 0.88)                  * (0.990) +
        (energy >= 2.0) * (abs(eta) >  0.88&&abs(eta) <= 3.0) * (0.990) +
        (abs(eta) > 3.0)                                      * (0.000)
  }
}

##################
# Photon isolation
##################

module Isolation PhotonIsolation {
  set CandidateInputArray PhotonEfficiency/photons
  set IsolationInputArray EFlowMerger/eflow

  set OutputArray photons
  set DeltaRMax 0.5
  set PTMin 0.5
  set PTRatioMax 999.
}

#################
# Photon filter
#################

module PdgCodeFilter PhotonFilter {
  set InputArray EFlowMerger/eflow
  set OutputArray1 WoPhotonpair
  set OutputArray2 Photonpair

  add EnMin   {35.0}
  add MassRes {125.25}
  add NP      {2}

  add PdgCode {22}
}

#################
# Electron filter
#################

module PdgCodeFilter ElectronFilter {
  set InputArray Calorimeter/eflowTracks
  set OutputArray1 electrons
  set Invert true
  add PdgCode { 11}
  add PdgCode {-11}
}



#####################
# Electron efficiency
#####################

module Efficiency ElectronEfficiency {
  set InputArray ElectronFilter/electrons
  set OutputArray electrons

  # set EfficiencyFormula {efficiency formula as a function of eta and pt}
  # efficiency formula for electrons
  set EfficiencyFormula {
        (energy <  2.0)                                       * (0.000)+
        (energy >= 2.0) * (abs(eta) <= 0.88)                  * (0.990)+
        (energy >= 2.0) * (abs(eta) >0.88 && abs(eta) <= 3.0) * (0.990)+
        (abs(eta)> 3.0)                                       * (0.000)
  }
}

####################
# Electron isolation
####################

module Isolation ElectronIsolation {
  set CandidateInputArray ElectronEfficiency/electrons
  set IsolationInputArray EFlowMerger/eflow
  set OutputArray         electrons

  set DeltaRMax  0.50
  set PTMin      0.50
  set PTRatioMax 0.12
}

#################
# Muon filter
#################

module PdgCodeFilter MuonFilter {
  set InputArray EFlowMerger/eflow
  set OutputArray1 Muon
  set Invert true
  set PTMin      0.5
#  set OutputArray2 MuonPair
#  add EnMin   {15.0}
#  add MassRes {91.18}
#  add NP      {2}

  add PdgCode {13}
  add PdgCode {-13}
}

#################
# Muon efficiency
#################

module Efficiency MuonEfficiency {
  set InputArray  MuonFilter/Muon
  set OutputArray muons

  # set EfficiencyFormula {efficiency as a function of eta and pt}
  # efficiency formula for muons
  set EfficiencyFormula {
        (energy < 2.0)                                         * (0.000)+
        (energy >= 2.0) * (abs(eta) <= 0.88)                   * (0.990)+
        (energy >= 2.0) * (abs(eta) >0.88 && abs(eta) <= 3.0)  * (0.990)+
        (abs(eta) > 3.0)                                       * (0.000)
  }
}

###################
# Missing ET merger
###################

module Merger MissingET {
  add InputArray EFlowMerger/eflow
  set MomentumOutputArray momentum
}

##################
# Scalar HT merger
##################

module Merger ScalarHT {
  add InputArray        UniqueObjectFinder/jets
  add InputArray        UniqueObjectFinder/electrons
  add InputArray        UniqueObjectFinder/photons
  add InputArray        UniqueObjectFinder/muons
  set EnergyOutputArray energy
}

#####################
# Neutrino Filter
#####################

module PdgCodeFilter NeutrinoFilter {

  set InputArray  Delphes/stableParticles
  set OutputArray1 filteredParticles
  set PTMin 0.0

  add PdgCode { 12}
  add PdgCode { 14}
  add PdgCode { 16}
  add PdgCode {-12}
  add PdgCode {-14}
  add PdgCode {-16}
}


#####################
# MC truth jet finder
#####################

module FastJetFinder GenJetFinder {
  set InputArray NeutrinoFilter/filteredParticles
  set OutputArray jets
  set ExclusiveClustering true

  set JetAlgorithm 10
  set ParameterR   2.5
  set ParameterP   1.0
  #set JetPTMin     1.0
  set Njets        2

}

#########################
# Gen Missing ET merger
########################

module Merger GenMissingET {
  add InputArray NeutrinoFilter/filteredParticles
  set MomentumOutputArray momentum
}

############
# Jet finder
############

module FastJetFinder FastJetFinder {
  set InputArray EFlowMerger/eflow
  set OutputArray jets
  set ExclusiveClustering true

  # algorithm: 1 CDFJetClu, 2 MidPoint, 3 SIScone, 4 kt, 5 Cambridge/Aachen, 6 antikt
  set JetAlgorithm 10
  set ParameterR   2.5
  set ParameterP   1.0
  set Njets        2
  #set JetPTMin     0.0
}

##################
# Jet Energy Scale
##################

module EnergyScale JetEnergyScale {
  set InputArray FastJetFinder/jets
  set OutputArray jets

  # scale formula for jets
  set ScaleFormula {1.025}
}

########################
# Jet Flavor Association
########################

module JetFlavorAssociation JetFlavorAssociation {

  set PartonInputArray Delphes/partons
  set ParticleInputArray Delphes/allParticles
  set ParticleLHEFInputArray Delphes/allParticlesLHEF
  set JetInputArray JetEnergyScale/jets

  set DeltaR       0.5
  set PartonPTMin  1.0
  set PartonEtaMax 3.0
}

###########
# b-tagging
###########

module BTagging BTagging {
  set JetInputArray JetEnergyScale/jets
  set BitNumber 0

  # add EfficiencyFormula {abs(PDG code)} {efficiency formula as a function of eta and pt}

  # default efficiency formula (misidentification rate)
  add EfficiencyFormula {0} {0.01}

  # efficiency formula for c-jets (misidentification rate)
  add EfficiencyFormula {4} {0.10}

  # efficiency formula for b-jets
  add EfficiencyFormula {5} {0.80}
}

#############
# tau-tagging
#############

module TauTagging TauTagging {
  set ParticleInputArray Delphes/allParticles
  set PartonInputArray Delphes/partons
  set JetInputArray JetEnergyScale/jets

  set DeltaR    0.5
  set TauPTMin  1.0
  set TauEtaMax 3.0

  # default efficiency formula (misidentification rate)
  add EfficiencyFormula { 0} {0.001}
  # efficiency formula for tau-jets
  add EfficiencyFormula {15} {0.600}
}


##################
# ROOT tree writer
##################

# Tracks, towers and eflow objects are not stored by default in the output.
# If needed (for jet constituent or other studies), uncomment the relevant
# "add Branch ..." lines.

module TreeWriter TreeWriter {
    # add Branch InputArray BranchName BranchClass
    add Branch Delphes/allParticles Particle GenParticle
    # add Branch TruthVertexFinder/vertices GenVertex Vertex

    # add Branch IdentificationMap/tracks Track Track
    #add Branch Calorimeter/towers Tower Tower

    #add Branch Calorimeter/eflowTracks EFlowTrack Track
    #add Branch Calorimeter/eflowPhotons EFlowPhoton Tower
    #add Branch TimeOfFlightNeutralHadron/eflowNeutralHadrons EFlowNeutralHadron Tower

    add Branch  EFlowMerger/eflow ParticleFlowCandidate ParticleFlowCandidate

    #add Branch Calorimeter/photons CaloPhoton Photon
    #add Branch PhotonEfficiency/photons PhotonEff Photon
    #add Branch PhotonIsolation/photons PhotonIso Photon
    #add Branch PhotonFilter/Photonpair Photonpair Photon

    #add Branch GenJetFinder/jets GenJet Jet
    #add Branch GenMissingET/momentum GenMissingET MissingET

    add Branch JetEnergyScale/jets Jet Jet
    #add Branch ElectronIsolation/electrons Electron Electron
    #add Branch PhotonIsolation/photons Photon Photon
    #add Branch MuonEfficiency/muons Muon Muon
    #add Branch MuonFilter/WoMuonPair WoMuonPair Muon

    #add Branch MissingET/momentum MissingET MissingET

    # add Info InfoName InfoValue
    add Info Bz $B
}
