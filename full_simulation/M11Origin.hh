#ifndef _M11Origin_hh_
#define _M11Origin_hh_

#include <string>
#include <iostream>
#include <fstream>
#include <marlin/Processor.h>
#include <TNtuple.h>
#include <TObject.h>
#include <TTree.h>
#include <TFile.h>
#include <TH1F.h>

#include <lcio.h>

#include <EVENT/ReconstructedParticle.h>
#include <IMPL/ReconstructedParticleImpl.h>
#include <EVENT/MCParticle.h>
#include <UTIL/LCRelationNavigator.h>

using namespace std ;
using namespace lcio ;
using namespace marlin ;

class TTree;

class M11Origin  : public marlin::Processor
{
	public:

		Processor*  newProcessor() { return new M11Origin ; }

		M11Origin();

		~M11Origin() {};

		void init();

		void processEvent( LCEvent * evtP );

		void end();

	protected:
		std::string _treeFileName;
		std::string _treeName;
		std::string _colName;
		std::string _colAdcVals;
//		TFile *tree_file;

		int _overwrite;
		TTree *_outputTree;

		unsigned int eventNr;
		int Num, count_jet;
    
   
               std::vector<double>  part_rdphi;
               std::vector<double>  part_rdtheta;
               std::vector<double>  part_rdlon;

               std::vector<double>  part_mom;

 
    
    std::vector<double > part_px;
    std::vector<double > part_py;
    std::vector<double > part_pz;
    std::vector<double> part_energy;
    std::vector<double> part_deta;
    std::vector<double> part_dphi;
    std::vector<double> part_pt;
    std::vector<double> part_d0val;
    std::vector<double> part_d0err;
    std::vector<double> part_dzval;
    std::vector<double> part_dzerr;
    std::vector<double> part_charge;

    std::vector<Int_t> part_isChargedHadron;
    std::vector<Int_t> part_isNeutralHadron;
    std::vector<Int_t> part_isPhoton;
    std::vector<Int_t> part_isElectron;
    std::vector<Int_t> part_isMuon;
    std::vector<Int_t> part_isPion;
    std::vector<Int_t> part_isChargedKaon;
    std::vector<Int_t> part_isProton;
    std::vector<Int_t> part_isKLong;
    std::vector<Int_t> part_isKShort;
    std::vector<Int_t>  part_isK0;
    std::vector<Int_t> part_isLambda;
    std::vector<Int_t> part_isSigma;
    std::vector<Int_t> part_isKxi;


    std::vector<Int_t> cpart_isChargedHadron;
    std::vector<Int_t> cpart_isNeutralHadron;
    std::vector<Int_t> cpart_isPhoton;
    std::vector<Int_t> cpart_isElectron;
    std::vector<Int_t> cpart_isMuon;
    std::vector<Int_t> cpart_isPion;
    std::vector<Int_t> cpart_isChargedKaon;
    std::vector<Int_t> cpart_isProton;

    std::vector<Int_t> BDaus;



    int KLong, KShort, K0, lambda, sigma, kxi;    
    float btag, ctag;
    bool label_bb, label_cc, label_gg;
    bool label_b, label_bbar, label_c, label_cbar, label_d, label_dbar, label_u, label_ubar, label_s, label_sbar, label_g;

    float quarks_Angle, quark1_En, quark2_En;    
    int jet_nparticles, charge_nparticles, BHadPDG;
    int jet_truthID;
    float jet_px, jet_py, jet_pz;
    float jet_costheta, jet_pt, jet_energy, jet_eta, leadMuonEn, leadElecEn, jet_charge1, jet_charge2, jet_charge3, jet_charge4, jet_charge5, jet_charge6;
    double thrust;
    
    float TM, TE, PtNeutrino, PtISR, miniThetaJet;
		std::string _fileName;
		std::ostream *_output;
		std::string _histFileName;
};

#endif


