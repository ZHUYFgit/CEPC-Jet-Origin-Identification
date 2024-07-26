#include <M11Origin.hh>
#include <EVENT/LCCollection.h>
#include <IMPL/LCCollectionVec.h>
#include <EVENT/LCFloatVec.h>
#include <EVENT/MCParticle.h>
#include <EVENT/Track.h>
#include <EVENT/ReconstructedParticle.h>
#include <IMPL/MCParticleImpl.h>
#include <EVENT/LCRelation.h>
#include <EVENT/Vertex.h>
#include <IMPL/LCRelationImpl.h>
#include <IMPL/LCFlagImpl.h>
#include <IMPL/ReconstructedParticleImpl.h>

#include <values.h>
#include <string>
#include <iostream>
#include <EVENT/LCFloatVec.h>
#include <EVENT/LCParameters.h>
#include <stdexcept>
#include <TFile.h>
#include <TTree.h>
#include <TVector3.h>
#include <TRandom.h>
#include <Rtypes.h>
#include <sstream>
#include <cmath>
#include <vector>
#include <TMath.h>
#include "TLorentzVector.h"
#include <UTIL/PIDHandler.h>

#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/Config.h"

using namespace std;

M11Origin a_M11Origin_instance;


M11Origin::M11Origin()
: Processor("M11Origin"),
_output(0)
{
    _description = "Print MC Truth" ;
    
    _treeFileName="MCTruth.root";
    registerProcessorParameter( "TreeOutputFile" ,
                               "The name of the file to which the ROOT tree will be written" ,
                               _treeFileName ,
                               _treeFileName);
    
    _treeName="tree";
    registerProcessorParameter( "TreeName" ,
                               "The name of the ROOT tree" ,
                               _treeName ,
                               _treeName);
    
//    registerOutputCollection( LCIO::RECONSTRUCTEDPARTICLE,
//                             "MCPSIMUFSP",
//                             "the final state particle after simulation",
//                             _outmcpsimufsp,
//                             std::string("ReconstructedParticle") );
//    
//    
//    registerOutputCollection( LCIO::LCRELATION,
//                             "mcpsimurelation",
//                             " relation between MCP and Reco after simulation",
//                             _outMCPSIMURelation,
//                             std::string("RelationusedSed") );
    
    _overwrite=0;
    registerProcessorParameter( "OverwriteFile" ,
                               "If zero an already existing file will not be overwritten." ,
                               _overwrite ,
                               _overwrite);
    
}

static int ISBbar(int mc);
static int ISB(int mc);
static int ISC(int mc);
static int ISCbar(int mc);
static int ISS(int mc);


static bool sortEn(ReconstructedParticle* a1, ReconstructedParticle* a2){
    return a1->getEnergy() >= a2->getEnergy();
}

static bool sortEnMC(MCParticle* a1, MCParticle* a2){
    return a1->getEnergy() >= a2->getEnergy();
}

typedef pair<ReconstructedParticle*, float> PAIR;
static float cmp(const PAIR& x, const PAIR& y){
    return x.second > y.second;
}

static double deltaPhi(double phi1, double phi2) { return TVector2::Phi_mpi_pi(phi1 - phi2); }
static void CalcuThrust(std::vector<TLorentzVector > UsedForThrust, std::vector<double> &result);
/*
float get_d0z0(float * pos, float * mom, float q, float B, float& d0return, float& z0return){
//    _referencePoint[0] = pos[0];
//    _referencePoint[1] = pos[1];
//    _referencePoint[2] = pos[2];
//    _momentum[0] = mom[0];
//    _momentum[1] = mom[1];
//    _momentum[2] = mom[2];
    float charge = q;
    float bField = B;
    float pxy = sqrt(mom[0]*mom[0]+mom[1]*mom[1]);
    float radius = _pxy / (_FCT*B);
    _omega = q/_radius;
    _tanLambda = mom[2]/_pxy;
    _phiMomRefPoint = atan2(mom[1],mom[0]);
    _xCentre = pos[0] + _radius*cos(_phiMomRefPoint-_const_pi2*q);
    _yCentre = pos[1] + _radius*sin(_phiMomRefPoint-_const_pi2*q);
    _phiRefPoint = atan2(pos[1]-_yCentre,pos[0]-_xCentre);
    _phiAtPCA = atan2(-_yCentre,-_xCentre);
    _phi0 = -_const_pi2*q + _phiAtPCA;
    while (_phi0<0) _phi0+=_const_2pi;
    while (_phi0>=_const_2pi) _phi0-=_const_2pi;
    _xAtPCA = _xCentre + _radius*cos(_phiAtPCA);
    _yAtPCA = _yCentre + _radius*sin(_phiAtPCA);
    //    _d0 = -_xAtPCA*sin(_phi0) + _yAtPCA*cos(_phi0);
    double pxy = double(_pxy);
    double radius = pxy/double(_FCT*B);
    double xCentre = double(pos[0]) + radius*double(cos(_phiMomRefPoint-_const_pi2*q));
    double yCentre = double(pos[1]) + radius*double(sin(_phiMomRefPoint-_const_pi2*q));

    double d0;

    if (q>0) {
      d0 = double(q)*radius - double(sqrt(xCentre*xCentre+yCentre*yCentre));
    }
    else {
      d0 = double(q)*radius + double(sqrt(xCentre*xCentre+yCentre*yCentre));
    }

    _d0 = FloatT(d0);

//     if (fabs(_d0)>0.001 ) {
//       std::cout << "New helix : " << std::endl;
//       std::cout << " Position : " << pos[0]
//         << " " << pos[1]
//         << " " << pos[2] << std::endl;
//       std::cout << " Radius = " << _radius << std::endl;
//       std::cout << " RC = " << sqrt(_xCentre*_xCentre+_yCentre*_yCentre) << std::endl;
//       std::cout << " D0 = " << _d0 << std::endl;
//     }

    _pxAtPCA = _pxy*cos(_phi0);
    _pyAtPCA = _pxy*sin(_phi0);
    FloatT deltaPhi = _phiRefPoint - _phiAtPCA;
    FloatT xCircles = -pos[2]*q/(_radius*_tanLambda) - deltaPhi;
    xCircles = xCircles/_const_2pi;
    int nCircles;
    int n1,n2;

    if (xCircles >= 0.) {
    n1 = int(xCircles);
    n2 = n1 + 1;
    }
    else {
    n1 = int(xCircles) - 1;
    n2 = n1 + 1;
    }

    if (fabs(n1-xCircles) < fabs(n2-xCircles)) {
    nCircles = n1;
    }
    else {
    nCircles = n2;
    }
    _z0 = pos[2] + _radius*_tanLambda*q*(deltaPhi + _const_2pi*nCircles);

    d0return = _d0;
    z0return = _z0;
    
}
*/

void M11Origin::init() {
    
    cout << "U are using this M11Origin !!!!!!!!!!" << endl;

    printParameters();
    TFile *tree_file=new TFile(_treeFileName.c_str(),(_overwrite ? "RECREATE" : "UPDATE"));
    if (!tree_file->IsOpen()) {
        delete tree_file;
        tree_file=new TFile(_treeFileName.c_str(),"NEW");
    }
    
    _outputTree = new TTree(_treeName.c_str(),_treeName.c_str());
    _outputTree->SetAutoSave(32*1024*1024);  // autosave every 32MB
    

    _outputTree->Branch("part_rdphi",     &part_rdphi);
    _outputTree->Branch("part_rdtheta",     &part_rdtheta);
    _outputTree->Branch("part_rdlon",     &part_rdlon);
    _outputTree->Branch("part_mom", &part_mom);
   
    
    _outputTree->Branch("part_px",     &part_px);
    _outputTree->Branch("part_py",     &part_py);
    _outputTree->Branch("part_pz",     &part_pz);
    _outputTree->Branch("part_energy", &part_energy);
    _outputTree->Branch("part_pt",     &part_pt);
    _outputTree->Branch("part_deta",   &part_deta);
    _outputTree->Branch("part_dphi",   &part_dphi);
    _outputTree->Branch("part_d0val",  &part_d0val);
    _outputTree->Branch("part_d0err",  &part_d0err);
    _outputTree->Branch("part_dzval",  &part_dzval);
    _outputTree->Branch("part_dzerr",  &part_dzerr);
    _outputTree->Branch("part_charge", &part_charge);
    
//part_isChargedHadron    
    _outputTree->Branch("part_isChargedHadron", &part_isChargedHadron);
    _outputTree->Branch("part_isNeutralHadron", &part_isNeutralHadron);
    _outputTree->Branch("part_isPhoton",        &part_isPhoton);
    _outputTree->Branch("part_isElectron",      &part_isElectron);
    _outputTree->Branch("part_isMuon",          &part_isMuon);
    _outputTree->Branch("part_isPion",          &part_isPion);
    _outputTree->Branch("part_isChargedKaon",   &part_isChargedKaon);
    _outputTree->Branch("part_isProton",        &part_isProton);   
    
    _outputTree->Branch("part_isKLong",         &part_isKLong);
    _outputTree->Branch("part_isKShort",        &part_isKShort);
 //   _outputTree->Branch("part_isK0",            &part_isK0);
//    _outputTree->Branch("part_isLambda",        &part_isLambda);
//    _outputTree->Branch("part_isSigma",         &part_isSigma);
//    _outputTree->Branch("part_isKxi",           &part_isKxi);


    _outputTree->Branch("cpart_isChargedHadron", &cpart_isChargedHadron);
    _outputTree->Branch("cpart_isNeutralHadron", &cpart_isNeutralHadron);
    _outputTree->Branch("cpart_isPhoton",        &cpart_isPhoton);
    _outputTree->Branch("cpart_isElectron",      &cpart_isElectron);
    _outputTree->Branch("cpart_isMuon",          &cpart_isMuon);
    _outputTree->Branch("cpart_isPion",          &cpart_isPion);
    _outputTree->Branch("cpart_isChargedKaon",   &cpart_isChargedKaon);
    _outputTree->Branch("cpart_isProton",        &cpart_isProton);


 
    _outputTree->Branch("EventNr",       &eventNr,        "EventNr/I");
    _outputTree->Branch("Num",           &Num,            "Num/I");
    
    _outputTree->Branch("BHadPDG",       &BHadPDG,        "BHadPDG/I");
    _outputTree->Branch("BDaus",         &BDaus);
    
    _outputTree->Branch("label_bb",     &label_bb,           "label_bb/O");
    _outputTree->Branch("label_cc",     &label_cc,           "label_cc/O");
    _outputTree->Branch("label_gg",     &label_gg,           "label_gg/O");

    _outputTree->Branch("label_b",     &label_b,           "label_b/O");
    _outputTree->Branch("label_bbar",  &label_bbar,        "label_bbar/O");
    _outputTree->Branch("label_c",     &label_c,           "label_c/O");
    _outputTree->Branch("label_cbar",  &label_cbar,        "label_cbar/O");
    _outputTree->Branch("label_u",     &label_u,           "label_u/O");
    _outputTree->Branch("label_ubar",  &label_ubar,        "label_ubar/O");
    _outputTree->Branch("label_d",     &label_d,           "label_d/O");
    _outputTree->Branch("label_dbar",  &label_dbar,        "label_dbar/O");
    _outputTree->Branch("label_s",     &label_s,           "label_s/O");
    _outputTree->Branch("label_sbar",  &label_sbar,        "label_sbar/O");
    _outputTree->Branch("label_g",     &label_g,           "label_g/O");


  
    _outputTree->Branch("jet_nparticles",&jet_nparticles,      "jet_nparticles/I");
    _outputTree->Branch("jet_costheta",  &jet_costheta,        "jet_costheta/F");
    _outputTree->Branch("jet_pt",        &jet_pt,              "jet_pt/F");
    _outputTree->Branch("jet_energy",    &jet_energy,          "jet_energy/F");
    _outputTree->Branch("jet_eta",       &jet_eta,             "jet_eta/F");
    _outputTree->Branch("jet_truthID",   &jet_truthID,         "jet_truthID/I");
    _outputTree->Branch("jet_px",        &jet_px,              "jet_px/F");
    _outputTree->Branch("jet_py",        &jet_py,              "jet_py/F");
    _outputTree->Branch("jet_pz",        &jet_pz,              "jet_pz/F");
    

    _outputTree->Branch("btag",          &btag,                "btag/F");
    _outputTree->Branch("ctag",          &ctag,                "ctag/F");
    _outputTree->Branch("thrust",        &thrust,              "thrust/D");
    
    _outputTree->Branch("charge_nparticles", &charge_nparticles, "charge_nparticles/I");
    _outputTree->Branch("leadMuonEn",        &leadMuonEn,        "leadMuonEn/F");
    _outputTree->Branch("leadElecEn",        &leadElecEn,        "leadElecEn/F");
    _outputTree->Branch("jet_charge1",        &jet_charge1,        "jet_charge1/F");
    _outputTree->Branch("TM",        &TM,        "TM/F");
    _outputTree->Branch("TE",        &TE,        "TE/F");
    _outputTree->Branch("PtISR",        &PtISR,        "PtISR/F");
    _outputTree->Branch("PtNeutrino",        &PtNeutrino,        "PtNeutrino/F");
    _outputTree->Branch("miniThetaJet",        &miniThetaJet,        "miniThetaJet/F");


    _outputTree->Branch("jet_charge2",        &jet_charge2,        "jet_charge2/F");
    _outputTree->Branch("jet_charge3",        &jet_charge3,        "jet_charge3/F");
    _outputTree->Branch("jet_charge4",        &jet_charge4,        "jet_charge4/F");
    _outputTree->Branch("jet_charge5",        &jet_charge5,        "jet_charge5/F");
    _outputTree->Branch("jet_charge6",        &jet_charge6,        "jet_charge6/F");
    _outputTree->Branch("quark1_En",        &quark1_En,        "quark1_En/F");
    _outputTree->Branch("quark2_En",        &quark2_En,        "quark2_En/F");
    _outputTree->Branch("quarks_Angle",        &quarks_Angle,        "quarks_Angle/F");	   

 
 
//    _outputTree->Branch("KLong",     &KLong,        "KLong/I");
//    _outputTree->Branch("KShort",    &KShort,       "KShort/I");
//    _outputTree->Branch("K0",        &K0,           "K0/I");
//    _outputTree->Branch("lambda",    &lambda,       "lambda/I");
//    _outputTree->Branch("sigma",     &sigma,        "sigma/I");
//    _outputTree->Branch("kxi",       &kxi,          "kxi/I");


//    count_jet = 0;
    Num = 0;
}

void M11Origin::processEvent( LCEvent * evtP )
{
    
    
    if (evtP)
    {
        try{
            
        //    cout<<"Next Event *******************************************************************************************************"<<endl;
            eventNr = evtP->getEventNumber();
        //    cout<<"eventNr : "<<eventNr<<" Num : "<<Num<<endl;
            
            label_bb = false, label_cc = false, label_gg = false;
            BHadPDG = 0;

	        quarks_Angle = 0, quark1_En = 0, quark2_En = 0;

	        miniThetaJet = 999, PtNeutrino = 999, PtISR = 999;  
	        TLorentzVector TLISR(0,0,0,0);
            TLorentzVector TLNeutrino(0,0,0,0);


            MCParticle* HDau1 = NULL;
            MCParticle* HDau2 = NULL;

            MCParticle* ZDau1 = NULL;
            MCParticle* ZDau2 = NULL;

            LCCollection* col_MCP = evtP->getCollection( "MCParticle" );
            int n_MCP = col_MCP->getNumberOfElements();
            
            std::vector<MCParticle*> BHadrons; BHadrons.clear();
            std::vector<MCParticle*> CHadrons; CHadrons.clear();
            std::vector<MCParticle*> SHadrons; SHadrons.clear();
            std::vector<MCParticle*> leadingHadrons; leadingHadrons.clear();
            int countQuark = 0;
            for(int i = 0; i<n_MCP; i++)
            {
                MCParticle* a_MCP = dynamic_cast<MCParticle*>(col_MCP->getElementAt(i));
                int NParents      = a_MCP->getParents().size();
                int NDaughters    = a_MCP->getDaughters().size();
                int PDG           = a_MCP->getPDG();
                TLorentzVector TLtemp(a_MCP->getMomentum(), a_MCP->getEnergy());

                if(NParents == 0 && abs(PDG) < 6 && abs(PDG) > 0){
                    if(countQuark == 0){ZDau1 = a_MCP; countQuark += 1;}
                    else{ZDau2 = a_MCP;}
                }
                if(PDG == 23 && NDaughters == 2){
                    MCParticle* ZDau = a_MCP->getDaughters()[0];
                    if( abs(ZDau->getPDG()) == 5){ label_bb = true; }
                    if( abs(ZDau->getPDG()) == 4){ label_cc = true; }
                    if( abs(ZDau->getPDG()) == 1 || abs(ZDau->getPDG()) == 2 || abs(ZDau->getPDG()) == 3 || abs(ZDau->getPDG()) == 21 ){ label_gg = true; }
                    if( abs(ZDau->getPDG()) > 0 && abs(ZDau->getPDG()) < 6 ){
                        ZDau1 = ZDau;
                        ZDau2 = a_MCP->getDaughters()[1];
                    }
                    
                }

                if(NParents == 0 && abs(PDG) == 22){TLISR += TLtemp;}
                if(a_MCP->getGeneratorStatus() == 1 && (abs(PDG) == 12 || abs(PDG) == 14 || abs(PDG) == 16) && NParents != 0){
		            MCParticle* parent = a_MCP;
  	                int PPDG = abs(parent->getPDG());
		            do{parent = parent->getParents()[0]; PPDG = abs(parent->getPDG());}
                    while(parent->getParents().size() != 0 && PPDG != 25);

		            if(PPDG == 25){

		                TLNeutrino += TLtemp;
		            }
		        }                

                if(NParents == 0 && abs(PDG) == 5){label_bb = true;}
                if(NParents == 0 && abs(PDG) == 4){label_cc = true;}
                if(NParents == 0 && (abs(PDG) == 1 || abs(PDG) == 2 || abs(PDG) == 3) ){label_gg = true;}
               
             

                if(PDG == 25 && NDaughters == 2){
                    MCParticle* HDau = a_MCP->getDaughters()[0];
                    if( abs(HDau->getPDG()) == 5){ label_bb = true; }
                    if( abs(HDau->getPDG()) == 4){ label_cc = true; }
                    if( abs(HDau->getPDG()) == 1 || abs(HDau->getPDG()) == 2 || abs(HDau->getPDG()) == 3 || abs(HDau->getPDG()) == 21 ){ label_gg = true; }
                

		            HDau1 = HDau;
                    HDau2 = a_MCP->getDaughters()[1];

		    // cout<<"HDau1->getPDG() "<<HDau1->getPDG()<<" HDau2->getPDG() "<<HDau2->getPDG()<<endl;

		            quark1_En = HDau1->getEnergy();
		            quark2_En = HDau2->getEnergy();
		            TLorentzVector TLq1(HDau1->getMomentum(), HDau1->getEnergy());
                    TVector3 TVq1 = TLq1.Vect();
		            TLorentzVector TLq2(HDau2->getMomentum(), HDau2->getEnergy());
                    TVector3 TVq2 = TLq2.Vect();
                    quarks_Angle = TVq1.Angle(TVq2);
			
		            miniThetaJet = abs(TLq1.CosTheta()) < abs(TLq2.CosTheta())? abs(TLq2.CosTheta()) : abs(TLq1.CosTheta());

                }

                if(PDG == 81){
                    for(int j = 0; j<NDaughters; j++){
                        MCParticle* Dau = a_MCP->getDaughters()[j];
                        if(ISBbar(Dau->getPDG()) || ISB(Dau->getPDG())){
                            BHadrons.push_back(Dau);
                        }
                        if(ISCbar(Dau->getPDG()) || ISC(Dau->getPDG())){
                            CHadrons.push_back(Dau);
                        }
                        if(ISS( abs(Dau->getPDG()) )){
                            SHadrons.push_back(Dau);
                        }
                    }
                }
                
                if(PDG == 92){
                    for(int j = 0; j<NDaughters; j++){
                        MCParticle* Dau = a_MCP->getDaughters()[j];
                        if( ISBbar(Dau->getPDG()) || ISB(Dau->getPDG()) ){
                            BHadrons.push_back(Dau);
                        }
                        if( ISCbar(Dau->getPDG()) || ISC(Dau->getPDG()) ){
                            CHadrons.push_back(Dau);
                        }
                        if( ISS( abs(Dau->getPDG()) ) ){
                            SHadrons.push_back(Dau);
                        }

                        // if( ISBbar(Dau->getPDG()) || ISB(Dau->getPDG()) || ISCbar(Dau->getPDG()) || ISC(Dau->getPDG()) ){
                        //     leadingHadrons.push_back( Dau );
                        // }

                    }
                    
                }

 
            }
            
            PtNeutrino = TLNeutrino.Perp();
            PtISR = TLISR.Perp();            
            
//            LCCollection* col_pfo = evtP->getCollection( "ArborPFOs" );
//            int n_pfo = col_pfo->getNumberOfElements();
//
//            for(int i = 0; i<n_pfo; i++)
//            {
//                ReconstructedParticle* a_pfo = dynamic_cast<ReconstructedParticle*>(col_pfo->getElementAt(i));
//                if(a_pfo->getCharge() != 0){
//                    Track* trk = a_pfo->getTracks()[0];
//                    cout<<"trk->getD0() : "<<trk->getD0()<<endl;
//                }
//
//            }


       
            LCCollection* col_Rela = evtP->getCollection( "RecoMCTruthLink" );
            int N_Link = col_Rela->getNumberOfElements();
 
            LCCollection* col_Jet = evtP->getCollection( "RefinedJets" );
            int num_jet = col_Jet->getNumberOfElements();

            // LCCollection* col_GJet = evtP->getCollection( "GraphJets" );
            // int num_Gjet = col_GJet->getNumberOfElements();

            // cout<<"num_Gjet : "<<num_Gjet<<endl;

           //	    if(num_jet < 2){ continue;}

       
        //    ReconstructedParticle* jeto = dynamic_cast<ReconstructedParticle*>(col_Jet->getElementAt(0));
        //    TLorentzVector TLjeto(jeto->getMomentum(), jeto->getEnergy());
        //    ReconstructedParticle* jett = dynamic_cast<ReconstructedParticle*>(col_Jet->getElementAt(1));
        //    TLorentzVector TLjett(jett->getMomentum(), jett->getEnergy());

//	        TM = 0;
  //          TM = (TLjeto + TLjett).M(); 
//	        TE = 0;
  //          TE = jeto->getEnergy() + jett->getEnergy();

            int count_jet = 0;
            for(int i = 0; i<num_jet; i++){
               
		        BDaus.clear();
                part_rdphi.clear();
                part_rdtheta.clear();
                part_rdlon.clear();
 
                part_mom.clear();
                part_px.clear();
                part_py.clear();
                part_pz.clear();
                part_energy.clear();
                part_pt.clear();
                part_deta.clear();
                part_dphi.clear();
                part_d0val.clear();
                part_d0err.clear();
                part_dzval.clear();
                part_dzerr.clear();
                part_charge.clear();
                
                part_isNeutralHadron.clear();
                part_isPhoton.clear();
                part_isElectron.clear();
                part_isMuon.clear();
                part_isPion.clear();
                part_isChargedKaon.clear();
                part_isProton.clear();  
                part_isChargedHadron.clear();
                part_isKLong.clear();
                part_isKShort.clear();


                cpart_isNeutralHadron.clear();
                cpart_isPhoton.clear();
                cpart_isElectron.clear();
                cpart_isMuon.clear();
                cpart_isPion.clear();
                cpart_isChargedKaon.clear();
                cpart_isProton.clear();
                cpart_isChargedHadron.clear();


 
                label_b = false, label_bbar = false, label_c = false, label_cbar = false, label_u = false, label_ubar = false;
                label_d = false, label_dbar = false, label_s = false, label_sbar = false, label_g = false;              
                
                btag = 0, ctag = 0;
                jet_nparticles = 0;
                jet_costheta = 0, jet_pt = 0, jet_energy = 0, jet_eta = 0;
                jet_truthID = 0, jet_px = 0, jet_py = 0, jet_pz = 0;
                
                ReconstructedParticle* jet1 = dynamic_cast<ReconstructedParticle*>(col_Jet->getElementAt(i));
                TLorentzVector TLjet1(jet1->getMomentum(), jet1->getEnergy());
                TVector3 TVjet1 = TLjet1.Vect();
                jet_px = jet1->getMomentum()[0]; jet_py = jet1->getMomentum()[1];
                jet_pz = jet1->getMomentum()[2];


                TVector3 ajet = TVjet1.Unit();
                TVector3 aphi = TVector3(0,0,1).Cross(ajet).Unit();
                TVector3 atheta = ajet.Cross(aphi).Unit();
 
                // if(BHadrons.size() > 1){
                //     sort(BHadrons.begin(), BHadrons.end(), sortEnMC);
                //     MCParticle* BHad1 = NULL;
                //     MCParticle* BHad2 = NULL;
                //     BHad1 = BHadrons.at(0);
                //     BHad2 = BHadrons.at(1);
                //     TLorentzVector TLBHad1(BHad1->getMomentum(), BHad1->getEnergy());
                //     TLorentzVector TLBHad2(BHad2->getMomentum(), BHad2->getEnergy());
                //     if( TVjet1.Angle( TLBHad1.Vect() ) < TVjet1.Angle( TLBHad2.Vect() ) ){
                //         BHadPDG = BHad1->getPDG();
                //     }
                //     else{
                //         BHadPDG = BHad2->getPDG();
                //     }
                    
                // }
                // else if(BHadrons.size() < 1 && CHadrons.size() > 1){
                //     sort(CHadrons.begin(), CHadrons.end(), sortEnMC);
                //     MCParticle* BHad1 = NULL;
                //     MCParticle* BHad2 = NULL;
                //     BHad1 = CHadrons.at(0);
                //     BHad2 = CHadrons.at(1);
                //     TLorentzVector TLBHad1(BHad1->getMomentum(), BHad1->getEnergy());
                //     TLorentzVector TLBHad2(BHad2->getMomentum(), BHad2->getEnergy());
                //     if( TVjet1.Angle( TLBHad1.Vect() ) < TVjet1.Angle( TLBHad2.Vect() ) ){
                //         BHadPDG = BHad1->getPDG();
                //     }
                //     else{
                //         BHadPDG = BHad2->getPDG();
                //     }
                // }
                // else if(BHadrons.size() < 1 && CHadrons.size() < 1 && SHadrons.size() > 1){
                //     sort(SHadrons.begin(), SHadrons.end(), sortEnMC);
                //     MCParticle* BHad1 = NULL;
                //     MCParticle* BHad2 = NULL;
                //     BHad1 = SHadrons.at(0);
                //     BHad2 = SHadrons.at(1);
                //     TLorentzVector TLBHad1(BHad1->getMomentum(), BHad1->getEnergy());
                //     TLorentzVector TLBHad2(BHad2->getMomentum(), BHad2->getEnergy());
                //     if( TVjet1.Angle( TLBHad1.Vect() ) < TVjet1.Angle( TLBHad2.Vect() ) ){
                //         BHadPDG = BHad1->getPDG();
                //     }
                //     else{
                //         BHadPDG = BHad2->getPDG();
                //     }
                // }

                if(HDau1 == NULL && HDau2 == NULL){
                    HDau1 = ZDau1;
                    HDau2 = ZDau2;
                }
                
                TLorentzVector TLHDau1(HDau1->getMomentum(), HDau1->getEnergy());
                TLorentzVector TLHDau2(HDau2->getMomentum(), HDau2->getEnergy());
                if( TVjet1.Angle( TLHDau1.Vect() ) < TVjet1.Angle( TLHDau2.Vect() ) ){
                    if(      HDau1->getPDG() == 5 ){ label_b = true; jet_truthID = 5;}
                    else if( HDau1->getPDG() == -5 ){ label_bbar = true; jet_truthID = -5;}
                    else if( HDau1->getPDG() == 4 ){ label_c = true; jet_truthID = 4;}
                    else if( HDau1->getPDG() == -4 ){ label_cbar = true; jet_truthID = -4;}
                    else if( HDau1->getPDG() == 1 ){ label_d = true; jet_truthID = 1;}
                    else if( HDau1->getPDG() == -1 ){ label_dbar = true;  jet_truthID = -1;}
                    else if( HDau1->getPDG() == 2 ){ label_u = true; jet_truthID = 2;}
                    else if( HDau1->getPDG() == -2 ){ label_ubar = true; jet_truthID = -2;}
                    else if( HDau1->getPDG() == 3 ){ label_s = true; jet_truthID = 3;}
                    else if( HDau1->getPDG() == -3 ){ label_sbar = true; jet_truthID = -3;}
                    else if( HDau1->getPDG() == 21 ){ label_g = true; jet_truthID = 21;}
                }
                else{


                    if(      HDau2->getPDG() == 5 ){ label_b = true; jet_truthID = 5;}
                    else if( HDau2->getPDG() == -5 ){ label_bbar = true; jet_truthID = -5;}
                    else if( HDau2->getPDG() == 4 ){ label_c = true; jet_truthID = 4;}
                    else if( HDau2->getPDG() == -4 ){ label_cbar = true; jet_truthID = -4;}
                    else if( HDau2->getPDG() == 1 ){ label_d = true; jet_truthID = 1;}
                    else if( HDau2->getPDG() == -1 ){ label_dbar = true; jet_truthID = -1;}
                    else if( HDau2->getPDG() == 2 ){ label_u = true; jet_truthID = 2;}
                    else if( HDau2->getPDG() == -2 ){ label_ubar = true; jet_truthID = -2;}
                    else if( HDau2->getPDG() == 3 ){ label_s = true; jet_truthID = 3;}
                    else if( HDau2->getPDG() == -3 ){ label_sbar = true; jet_truthID = -3;}
                    else if( HDau2->getPDG() == 21 ){ label_g = true; jet_truthID = 21;}


                }

                if(   abs(HDau1->getPDG()) == 5 ){leadingHadrons = BHadrons;}
                else if( abs(HDau1->getPDG()) == 4 ){leadingHadrons = CHadrons;}
                if(   abs(HDau1->getPDG()) == 3 ){leadingHadrons = SHadrons;}
                
                if(leadingHadrons.size() > 1){
                    sort(leadingHadrons.begin(), leadingHadrons.end(), sortEnMC);
                    MCParticle* BHad1 = NULL;
                    MCParticle* BHad2 = NULL;
                    BHad1 = leadingHadrons.at(0);
                    BHad2 = leadingHadrons.at(1);
                    TLorentzVector TLBHad1(BHad1->getMomentum(), BHad1->getEnergy());
                    TLorentzVector TLBHad2(BHad2->getMomentum(), BHad2->getEnergy());
                    if( TVjet1.Angle( TLBHad1.Vect() ) < TVjet1.Angle( TLBHad2.Vect() ) ){
                        BHadPDG = BHad1->getPDG();
                        
			            for(int k = 0; k<BHad1->getDaughters().size(); k++){
			                MCParticle* Bdau = BHad1->getDaughters()[k];
			                BDaus.push_back( Bdau->getPDG() );
			            }
                    }
                    else{
                        BHadPDG = BHad2->getPDG();
			            for(int k = 0; k<BHad2->getDaughters().size(); k++){
                            MCParticle* Bdau = BHad2->getDaughters()[k];
                            BDaus.push_back( Bdau->getPDG() );
                        }
			
                    }

			cout<<"test"<<endl;

	            //    cout<<"BHadPDG : "<<BHadPDG<<endl;			
		    //        for(int k = 0; k<BDaus.size(); k++){cout<<"BDaus.at(k) : "<<BDaus.at(k)<<endl;}

                }

           
                jet_energy = jet1->getEnergy(); jet_pt = TLjet1.Pt(); jet_costheta = TVjet1.CosTheta();
                jet_eta = TLjet1.Eta();
                
                PIDHandler pidh(col_Jet);
                int algo   = 0;
                int ibtag  =-1;
                int ictag  =-1;
                int ibctag =-1;
                int icat   =-1;
                algo      = pidh.getAlgorithmID("lcfiplus");
                ibtag     = pidh.getParameterIndex (algo,  "BTag");
                ictag     = pidh.getParameterIndex (algo,  "CTag");
                ibctag    = pidh.getParameterIndex (algo,  "BCTag");
                icat      = pidh.getParameterIndex (algo,  "Category");
                
                const ParticleID &pid1 = pidh.getParticleID(jet1, algo);
                btag  = pid1.getParameters()[ibtag];
                ctag  = pid1.getParameters()[ictag];
                                
                ReconstructedParticleVec comp = jet1->getParticles();
                jet_nparticles = comp.size();
 
                std::map<ReconstructedParticle*, float> recPt; recPt.clear();
                std::vector<TLorentzVector > recTL; recTL.clear();               
                recPt.clear();
                std::vector<ReconstructedParticle* > muons; muons.clear();
                std::vector<ReconstructedParticle* > electrons; electrons.clear();
                charge_nparticles = 0;
                jet_charge1 = 0;
                jet_charge2 = 0;
                jet_charge3 = 0;
                jet_charge4 = 0;
                jet_charge5 = 0;
                jet_charge6 = 0;



                for(int j = 0; j<jet_nparticles; j++){
                    ReconstructedParticle* rec = comp.at(j);
                    if(rec->getCharge() != 0){charge_nparticles += 1;}
                    TLorentzVector TLrec(rec->getMomentum(), rec->getEnergy());

                    TVector3 partV3 = TLrec.Vect();
                    TVector3 partV3_ = partV3.Unit();
                    float rdphi = TMath::ASin(partV3_.Dot(aphi));
                    float rdtheta = TMath::ASin(partV3_.Dot(atheta));
                    float rdlon = partV3_.Dot(ajet);
	    
 
                    part_rdphi.push_back(rdphi);
                    part_rdtheta.push_back(rdtheta);
                    part_rdlon.push_back(rdlon);

                    recTL.push_back(TLrec);
                    recPt[rec] = TLrec.E();
                    jet_charge1 += pow( (rec->getEnergy()/jet_energy), 0.1 ) * rec->getCharge();
                    jet_charge2 += pow( (rec->getEnergy()/jet_energy), 0.2 ) * rec->getCharge();
                    jet_charge3 += pow( (rec->getEnergy()/jet_energy), 0.3 ) * rec->getCharge();
                    jet_charge4 += pow( (rec->getEnergy()/jet_energy), 0.4 ) * rec->getCharge();
                    jet_charge5 += pow( (rec->getEnergy()/jet_energy), 0.5 ) * rec->getCharge();
                    jet_charge6 += pow( (rec->getEnergy()/jet_energy), 0.6 ) * rec->getCharge();
                    
                    if( abs(rec->getType()) == 11){electrons.push_back(rec);}
                    if( abs(rec->getType()) == 13){muons.push_back(rec);}
                }
                
                leadElecEn = 999;
                leadMuonEn = 999;
                if(electrons.size() != 0){
                    sort(electrons.begin(), electrons.end(), sortEn);
                    ReconstructedParticle* leadElec = electrons.at(0);
                    leadElecEn = leadElec->getEnergy();
                }
                if(muons.size() != 0){
                    sort(muons.begin(), muons.end(), sortEn);
                    ReconstructedParticle* leadMuon = muons.at(0);
                    leadMuonEn = leadMuon->getEnergy();
                }
                
 
                std::vector<double> thrustJet;
                thrust = 0;
                CalcuThrust(recTL, thrustJet);
                thrust = thrustJet.at(0);
                
                
                std::vector<PAIR> vec(recPt.begin(), recPt.end());
                sort(vec.begin(), vec.end(), cmp);
                
                std::map<MCParticle*, int> neutral_s_had; neutral_s_had.clear();

                for(int j = 0; j<vec.size(); j++){
                    ReconstructedParticle* rec = vec[j].first;
                    TLorentzVector TLrec(rec->getMomentum(), rec->getEnergy());
                    TVector3 TVrec = TLrec.Vect();
                  

                    if(rec->getCharge() != 0 && rec->getTracks().size()){
                        Track* trk = rec->getTracks()[0];
                        float d0 = trk->getD0();
                        float z0 = trk->getZ0();
                        if(abs(d0) < 1e-6){
                            if(d0 < 0){d0 = -1e-6;}
                            if(d0 > 0){d0 = 1e-6;}
                        }
                        if(abs(z0) < 1e-6){
                            if(z0 < 0){z0 = -1e-6;}
                            if(z0 > 0){z0 = 1e-6;}
                        }


                        part_d0val.push_back( 0.2 * ((d0 > 0) ? 1 : -1) * std::pow(std::abs(d0), 0.25)  );
                        float d0err = sqrt(trk->getCovMatrix()[0]);
                        float dzerr = sqrt(trk->getCovMatrix()[9]);
                        part_d0err.push_back( 0.2 * ((d0err > 0) ? 1 : -1) * std::pow(std::abs(d0err), 0.25) );
                        part_dzval.push_back( 0.2 * ((z0 > 0) ? 1 : -1) * std::pow(std::abs(z0), 0.25)  );
                        part_dzerr.push_back( 0.2 * ((dzerr > 0) ? 1 : -1) * std::pow(std::abs(dzerr), 0.25) );

                    }
                    else{
                        part_d0val.push_back(0);
                        part_d0err.push_back(0);
                        part_dzval.push_back(0);
                        part_dzerr.push_back(0);
                    
                    }
                    
                    part_px.push_back( rec->getMomentum()[0] );
                    part_py.push_back( rec->getMomentum()[1] );
                    part_pz.push_back( rec->getMomentum()[2] );
                    part_energy.push_back( rec->getEnergy() );
                    part_pt.push_back( TLrec.Pt() );
                    part_mom.push_back( pow( pow(rec->getMomentum()[0], 2) + pow(rec->getMomentum()[1], 2) + pow(rec->getMomentum()[2], 2)  ,0.5) );   

                    float deta = (TLjet1.Eta() > 0 ? 1 : -1) * (TLrec.Eta() - TLjet1.Eta());
                    float dphi = deltaPhi(TLrec.Phi(), TLjet1.Phi());
                    if(abs(dphi) < 1e-6){
                        if(dphi < 0){dphi = -1e-6;}
                        if(dphi > 0){dphi = 1e-6;}
                    }
                    if(abs(deta) < 1e-6){
                        if(deta < 0){deta = -1e-6;}
                        if(deta > 0){deta = 1e-6;}
                    }
                    part_deta.push_back( deta );
                    part_dphi.push_back( dphi );
                    part_charge.push_back(rec->getCharge());
                    

                   
                    cpart_isElectron.push_back(0);
                    cpart_isMuon.push_back(0);
                    cpart_isChargedHadron.push_back(0);
                    cpart_isPion.push_back(0);
                    cpart_isProton.push_back(0);
                    cpart_isPhoton.push_back(0);
                    cpart_isNeutralHadron.push_back(0);
                    cpart_isChargedKaon.push_back(0);
 
                    int type = abs(rec->getType());
                    if(rec->getCharge() != 0){
                        if(  type == 11 ){    cpart_isElectron.back() = 1; }
                        else if(type == 13){  cpart_isMuon.back() = 1; }
                        else if(type == 211){ cpart_isPion.back() = 1; cpart_isChargedHadron.back() = 1; }
                        else if(type == 321){ cpart_isChargedKaon.back() = 1; cpart_isChargedHadron.back() = 1; }
                        else if(type == 2212){ cpart_isPhoton.back() = 1; cpart_isChargedHadron.back() = 1; }
                        else{               cpart_isChargedHadron.back() = 1;}
                    }
                    else {
                        if(  type == 22 ){ cpart_isPhoton.back() = 1; }
                        else{             cpart_isNeutralHadron.back() = 1; }
                    }

 
                    part_isKLong.push_back(0);
                    part_isKShort.push_back(0);
//                    part_isK0.push_back(0);
//                    part_isLambda.push_back(0);
//                    part_isSigma.push_back(0);
//                    part_isKxi.push_back(0);
                    

                    int PDG = 0;
                    for(int k=0; k<N_Link; k++)
                    {
                        LCRelation* a_link = dynamic_cast<LCRelation*>(col_Rela->getElementAt(k));

                        if(a_link->getFrom() == rec){
                            MCParticle* a_mcp = dynamic_cast<MCParticle*>(a_link->getTo());
                            PDG = abs(a_mcp->getPDG());
                            MCParticle* parent = a_mcp;
                            int pPDG = abs(parent->getPDG());

                            if(pPDG == 130 || pPDG == 310 || pPDG == 311 || pPDG == 3122 || pPDG == 3212 || pPDG == 3322){
                                neutral_s_had[parent] += 1;
                                if(pPDG == 130){part_isKLong.back() = 1;}
                                else if(pPDG == 310){part_isKShort.back() = 1;}
                            }
                            else if(parent->getParents().size() != 0){
                                do{parent = parent->getParents()[0]; pPDG = abs(parent->getPDG()); }
                                while(parent->getParents().size() != 0 && pPDG != 130 && pPDG != 310 && pPDG != 311 && pPDG != 3122 && pPDG != 3212 && pPDG != 3322);
                                if(pPDG == 130 || pPDG == 310 || pPDG == 311 || pPDG == 3122 || pPDG == 3212 || pPDG == 3322){
                                    neutral_s_had[parent] += 1;
                                    if(pPDG == 130 || PDG == 130){part_isKLong.back() = 1;}
                                    else if(pPDG == 310){part_isKShort.back() = 1;}
           //                         else if(pPDG == 311){part_isK0.back() = 1;}
           //                         else if(pPDG == 3122){part_isLambda.back() = 1;}
           //                         else if(pPDG == 3212){part_isSigma.back() = 1;}
           //                         else if(pPDG == 3322){part_isKxi.back() = 1;}
                                }
                        //
                            }

                        }

                        
                        
                    }

                    if(rec->getCharge() != 0){
                        if(  PDG == 11 ){part_isElectron.push_back(1); part_isChargedHadron.push_back(0); part_isMuon.push_back(0); part_isPion.push_back(0); part_isChargedKaon.push_back(0); part_isProton.push_back(0);}
                        else if(PDG == 13){ part_isMuon.push_back(1); part_isChargedHadron.push_back(0); part_isElectron.push_back(0); part_isPion.push_back(0); part_isChargedKaon.push_back(0); part_isProton.push_back(0);}
                        else if(PDG == 211){ part_isPion.push_back(1); part_isChargedHadron.push_back(1); part_isElectron.push_back(0); part_isMuon.push_back(0); part_isChargedKaon.push_back(0); part_isProton.push_back(0);}
                        else if(PDG == 321){part_isMuon.push_back(0); part_isChargedHadron.push_back(1); part_isElectron.push_back(0); part_isPion.push_back(0); part_isChargedKaon.push_back(1); part_isProton.push_back(0);}
                        else if(PDG == 2212){part_isMuon.push_back(0); part_isChargedHadron.push_back(1); part_isElectron.push_back(0); part_isPion.push_back(0); part_isChargedKaon.push_back(0); part_isProton.push_back(1);}
                        else{part_isMuon.push_back(0); part_isChargedHadron.push_back(1); part_isElectron.push_back(0); part_isPion.push_back(0); part_isChargedKaon.push_back(0); part_isProton.push_back(0);}
                        part_isPhoton.push_back(0); part_isNeutralHadron.push_back(0);
                    }
                    else{
                        if(  PDG == 22 ){part_isPhoton.push_back(1); part_isNeutralHadron.push_back(0);}
                        else{ part_isPhoton.push_back(0); part_isNeutralHadron.push_back(1);}
                        part_isElectron.push_back(0); part_isMuon.push_back(0); part_isPion.push_back(0); part_isChargedKaon.push_back(0); part_isProton.push_back(0); part_isChargedHadron.push_back(0);
                    }
                }
               
                KLong = 0; KShort = 0; K0 = 0; lambda = 0; sigma = 0; kxi = 0;
                std::map<MCParticle*, int>::iterator it;
                for( it = neutral_s_had.begin(); it != neutral_s_had.end(); it++){
//                    cout<<"it->first : "<<it->first<<" it->second : "<<it->second<<endl;
                    MCParticle* temp = it->first;
                    int PDG = 0;
                    PDG = abs(temp->getPDG());
                    if(PDG == 130){KLong += 1;}
                    else if(PDG == 310){KShort += 1;}
                    else if(PDG == 311){K0 += 1;}
                    else if(PDG == 3122){lambda += 1;}
                    else if(PDG == 3212){sigma += 1;}
                    else if(PDG == 3322){kxi += 1;}
                }
                
                // for(int i = 0; i<part_px.size(); i++){
                //     cout<<"part_px.at(i) : "<<part_px.at(i)<<endl;
                // }

//                cout<<"label_bbar : "<<label_bbar<<endl;  // && abs(jet_costheta) < 0.8
                if( (label_bbar || label_cbar || label_dbar || label_ubar || label_sbar || label_g) && count_jet < 1 ){
                    count_jet += 1;
//                if(count_jet < 16000){
//                    cout<<"label_b : "<<label_b<<endl;
//                    cout<<"BHadPDG : "<<BHadPDG<<endl;
                        _outputTree->Fill();
//		           }
                }
                
            }
            

   

        }catch (lcio::DataNotAvailableException err) {  }

    }
    
    
    Num ++;
}



void M11Origin::end()
{
    
    if (_outputTree) {
        
        TFile *tree_file = _outputTree->GetCurrentFile(); //just in case we switched to a new file
        //tree_file->cd();
        tree_file->Write();
        delete tree_file;
        //tree_file->Close();
    }
    
}


static void CalcuThrust(std::vector<TLorentzVector > UsedForThrust, std::vector<double> &result){
    result.clear();
    double T = 0;
    double thetaMin = TMath::Pi(), phiMin = 2*TMath::Pi();
    double thetaMin2 = 0, phiMin2 = 0;
    double thetaL = 0, phiL = 0, thetaR = TMath::Pi(), phiR = 2*TMath::Pi();
    int iter = 0;
    double Told = 0;
    double Tnew = 0;
    double cut = 1;
    double thetaRange = 0, phiRange = 0;
    do{ 
        iter += 1; 
        if(iter == 1){ 
            thetaRange = thetaR - thetaL, phiRange = phiR - phiL;
        }
        else if(iter != 1){
            
            thetaRange = 0.1*(thetaR - thetaL);
            phiRange = 0.1*(phiR - phiL);
            
            thetaL =  thetaMin - thetaRange;
            thetaR = thetaMin + thetaRange;
            phiL = phiMin - phiRange;
            phiR = phiMin + phiRange;
            thetaRange = thetaR - thetaL, phiRange = phiR - phiL;
        }
        for(double theta = thetaL; theta <= thetaR; theta += 0.1*thetaRange){   //in this round, find the max T
            for(double phi = phiL; phi <= phiR; phi += 0.1*phiRange){
                
                double x = sin(theta)*cos(phi);
                double y = sin(theta)*sin(phi);
                double z = cos(theta);
                
                double denominator = 0;
                double numerator = 0;
                for(int i = 0; i<UsedForThrust.size(); i++){
                    TLorentzVector TLtemp = UsedForThrust.at(i);
                    TVector3 TVtemp = TLtemp.Vect();
                    denominator += TVtemp.Mag(); 
                    numerator += abs(x*TVtemp(0) + y*TVtemp(1) + z*TVtemp(2));
                }
                double Ttemp = numerator/denominator;
                if(Ttemp > T){
                    thetaMin = theta;   phiMin = phi; T = Ttemp;
                }
            }
        }
        if(iter == 1){Told = T; Tnew = T;}
        else if(T >= Tnew && iter != 1){
            Told = Tnew; Tnew = T; cut = (Tnew - Told)/Tnew;
        }
    }
    while(cut >= 0.2);
    result.push_back(T);
}


static int ISBbar(int mc){
    if(mc == 511 || mc == 521 || mc == 10511 || mc == 10521 || mc == 513 || mc == 523 || mc == 10513 || mc == 10523 || mc == 20513 || mc == 20523 || mc == 515 || mc == 525 || mc == 531 || mc == 10531 || mc == 533 || mc == 10533 || mc == 20533 || mc == 535 || mc == 541 || mc == 10541 || mc == 543 || mc == 10543 || mc == 20543 || mc == 545 || mc == -5122 || mc == -5112 || mc == -5212 || mc == -5222 || mc == -5114 || mc == -5214 || mc == -5224 || mc == -5132 || mc == -5232 || mc == -5312 || mc == -5322 || mc == -5314 || mc == -5324 || mc == -5332 || mc == -5334 || mc == -5142 || mc == -5242 || mc == -5412 || mc == -5422 || mc == -5414 || mc == -5424 || mc == -5342 || mc == -5432 || mc == -5434 || mc == -5442 || mc == -5444 || mc == -5512 || mc == -5522 || mc == -5514 || mc == -5524 || mc == -5532 || mc == -5534 || mc == -5542 || mc == -5544 || mc == -5544 ) {return 1;}
    else {return 0;}
}

static int ISB(int mc){
    if(mc == -511 || mc == -521 || mc == -10511 || mc == -10521 || mc == -513 || mc == -523 || mc == -10513 || mc == -10523 || mc == -20513 || mc == -20523 || mc == -515 || mc == -525 || mc == -531 || mc == -10531 || mc == -533 || mc == -10533 || mc == -20533 || mc == -535 || mc == -541 || mc == -10541 || mc == -543 || mc == -10543 || mc == -20543 || mc == -545 || mc == 5122 || mc == 5112 || mc == 5212 || mc == 5222 || mc == 5114 || mc == 5214 || mc == 5224 || mc == 5132 || mc == 5232 || mc == 5312 || mc == 5322 || mc == 5314 || mc == 5324 || mc == 5332 || mc == 5334 || mc == 5142 || mc == 5242 || mc == 5412 || mc == 5422 || mc == 5414 || mc == 5424 || mc == 5342 || mc == 5432 || mc == 5434 || mc == 5442 || mc == 5444 || mc == 5512 || mc == 5522 || mc == 5514 || mc == 5524 || mc == 5532 || mc == 5534 || mc == 5542 || mc == 5544 || mc == 5544 ) {return 1;}
    else {return 0;}

}

static int ISC(int mc){
    if(mc == 411 || mc == 421 || mc == 10411 || mc == 10421 || mc == 413 || mc == 423 || mc == 10413 || mc == 10423 || mc == 20413 || mc == 20423 || mc == 415 || mc == 425 || mc == 431 || mc == 10431 || mc == 433 || mc == 10433 || mc == 20433 || mc == 435      || mc == 4122 || mc == 4222 || mc == 4212 || mc == 4112 || mc == 4224 || mc == 4214 || mc == 4114 || mc == 4232 || mc == 4132 || mc == 4322 || mc == 4312 || mc == 4324 || mc == 4314 || mc == 4332 || mc == 4334 || mc == 4412 || mc == 4422 || mc == 4414 || mc == 4424 || mc == 4432 || mc == 4434 || mc == 4444) {return 1;}
    else {return 0;}
}

static int ISCbar(int mc){
    if(mc == -411 || mc == -421 || mc == -10411 || mc == -10421 || mc == -413 || mc == -423 || mc == -10413 || mc == -10423 || mc == -20413 || mc == -20423 || mc == -415 || mc == -425 || mc == -431 || mc == -10431 || mc == -433 || mc == -10433 || mc == -20433 || mc == -435      || mc == -4122 || mc == -4222 || mc == -4212 || mc == -4112 || mc == -4224 || mc == -4214 || mc == -4114 || mc == -4232 || mc == -4132 || mc == -4322 || mc == -4312 || mc == -4324 || mc == -4314 || mc == -4332 || mc == -4334 || mc == -4412 || mc == -4422 || mc == -4414 || mc == -4424 || mc == -4432 || mc == -4434 || mc == -4444) {return 1;}
    else {return 0;}
}

static int ISS(int mc){
    if(mc == 130 || mc == 310 || mc == 311 || mc == 321 || mc == 10311 || mc == 10321 || mc == 100311 || mc == 100321 || mc == 200311 || mc == 200321 || mc == 9000311 || mc == 9000321 || mc == 313 || mc == 323 || mc == 10313 || mc == 10323 || mc == 20313 || mc == 20323 || mc == 100313 || mc == 100323 || mc == 9000313 || mc == 9000323 || mc == 30313 || mc == 30323 || mc == 315 || mc == 325 || mc == 9000315 || mc == 9000325 || mc == 10315 || mc == 10325 || mc == 20315 || mc == 20325 || mc == 100315 || mc == 100325 || mc == 9010315 || mc == 9010325 || mc == 317 || mc == 327 || mc == 9010317 || mc == 9010327 || mc == 319 || mc == 329 || mc == 9000319 || mc == 9000329       || mc == 3122 || mc == 3222 || mc == 3212 || mc == 3112 || mc == 3224 || mc == 3214 || mc == 3114 || mc == 3322 || mc == 3312 || mc == 3324 || mc == 3314 || mc == 3334) {return 1;}
    else {return 0;}
}
