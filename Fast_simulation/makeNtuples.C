// #include <iostream>
// #include <unordered_set>
// #include <utility>
// #include "TClonesArray.h"
// #include "Math/LorentzVector.h"
// #include "classes/DelphesClasses.h"
// #include "external/ExRootAnalysis/ExRootTreeReader.h"

#ifdef __CLING__
R__LOAD_LIBRARY(/afs/ihep.ac.cn/users/z/zhuyf/cefs/workspace/Delphes/Delphes_CEPC-1.3/libDelphes.so)
#include "/afs/ihep.ac.cn/users/z/zhuyf/cefs/workspace/Delphes/Delphes_CEPC-1.3/classes/DelphesClasses.h"
#include "/afs/ihep.ac.cn/users/z/zhuyf/cefs/workspace/Delphes/Delphes_CEPC-1.3/external/ExRootAnalysis/ExRootTreeReader.h"
#else
class ExRootTreeReader;
#endif

#include <TVector3.h>
#include "TLorentzVector.h"
//#include "OnnxReader.h"

double deltaPhi(double phi1, double phi2) { return TVector2::Phi_mpi_pi(phi1 - phi2); }

double deltaR(double eta1, double phi1, double eta2, double phi2) {
  double deta = eta1 - eta2;
  double dphi = deltaPhi(phi1, phi2);
  return std::hypot(deta, dphi);
}

std::vector<float> Preprocessing(std::vector<float> input, size_t resize=50, float shift = 0, float scale = 1, float min = -1e9, float max = 1e9) {
    std::vector<float> res(resize);
    for (size_t i = 0; i < resize; i++ ) {
        if (i >= input.size()) { // reach the end of input vector, remaining containt of output should be zero.
            break;
        }
        res.at(i) = std::clamp((float(input.at(i)) - shift) * scale,
                               min, max);
    }

    return res;
};

template <class T1, class T2>
double deltaR(const T1 &a, const T2 &b) {
  return deltaR(a->Eta, a->Phi, b->Eta, b->Phi);
}

namespace ParticleID {
  enum PdgId {
    p_unknown,
    p_d,
    p_u,
    p_s,
    p_c,
    p_b,
    p_t,
    p_bprime,
    p_tprime,
    p_eminus = 11,
    p_nu_e,
    p_muminus,
    p_nu_mu,
    p_tauminus,
    p_nu_tau,
    p_tauprimeminus,
    p_nu_tauprime,
    p_g = 21,
    p_gamma,
    p_Z0,
    p_Wplus,
    p_h0,
    p_Zprime0 = 32,
    p_Zpprime0,
    p_Wprimeplus,
    p_H0,
    p_A0,
    p_Hplus,
    p_G = 39,
    p_R0 = 41,
    p_H30 = 45,
    p_A20 = 46,
    p_LQ,
    p_cluster = 91,
    p_string,
    p_pi0 = 111,
    p_rho0 = 113,
    p_klong = 130,
    p_piplus = 211,
    p_rhoplus = 213,
    p_eta = 221,
    p_omega = 223,
    p_kshort = 310,
    p_k0,
    p_kstar0 = 313,
    p_kplus = 321,
    p_kstarplus = 323,
    p_phi = 333,
    p_dplus = 411,
    p_d0 = 421,
    p_dsplus = 431,
    p_b0 = 511,
    p_bplus = 521,
    p_bs0 = 531,
    p_bcplus = 541,
    p_neutron = 2112,
    p_proton = 2212,
    p_sigmaminus = 3112,
    p_lambda0 = 3122,
    p_sigma0 = 3212,
    p_sigmaplus = 3222,
    p_ximinus = 3312,
    p_xi0 = 3322,
    p_omegaminus = 3334,
    p_sigmac0 = 4112,
    p_lambdacplus = 4122,
    p_xic0 = 4132,
    p_sigmacplus = 4212,
    p_sigmacpp = 4222,
    p_xicplus = 4232,
    p_omegac0 = 4332,
    p_sigmabminus = 5112,
    p_lambdab0 = 5122,
    p_xibminus = 5132,
    p_sigmab0 = 5212,
    p_sigmabplus = 5222,
    p_xib0 = 5232,
    p_omegabminus = 5332,
  };
}


std::vector<const GenParticle *> getLabel(const TClonesArray *branchParticle) {
  std::vector<const GenParticle *> genParticles;
  genParticles.clear();
  for (Int_t i = 0; i < branchParticle->GetEntriesFast(); ++i) {
    genParticles.push_back((GenParticle *)branchParticle->At(i));
  }
  return genParticles;
}

std::vector<const GenParticle *> getDaughterQuarks(std::vector<const GenParticle *> genParticles, const GenParticle *particle, bool allow_gluon = true) {
    std::vector<const GenParticle *> daughters;
    // for (int i = 0; i<particle->daughter().size(); ++i){
    //   const auto *dau = particle->daughter(i);
    //   if( (abs(dau->PID) > 0 && abs(dau->PID) < 6) || abs(dau->PID) == 22 ){
    //     daughters.push_back(dau);
    //   }
    // }
    for (int idau = particle->D1; idau <= particle->D2; ++idau) {
      const auto *dau = genParticles.at(idau);
      auto pdgid = std::abs(dau->PID);
      if (pdgid >= ParticleID::p_d && pdgid <= ParticleID::p_b) {
        daughters.push_back(dau);
      }
      if (allow_gluon && pdgid == ParticleID::p_g) {
        daughters.push_back(dau);
      }
    }
    return daughters;
  }


int ISBbar(int mc){
    if(mc == 511 || mc == 521 || mc == 10511 || mc == 10521 || mc == 513 || mc == 523 || mc == 10513 || mc == 10523 || mc == 20513 || mc == 20523 || mc == 515 || mc == 525 || mc == 531 || mc == 10531 || mc == 533 || mc == 10533 || mc == 20533 || mc == 535 || mc == 541 || mc == 10541 || mc == 543 || mc == 10543 || mc == 20543 || mc == 545 || mc == -5122 || mc == -5112 || mc == -5212 || mc == -5222 || mc == -5114 || mc == -5214 || mc == -5224 || mc == -5132 || mc == -5232 || mc == -5312 || mc == -5322 || mc == -5314 || mc == -5324 || mc == -5332 || mc == -5334 || mc == -5142 || mc == -5242 || mc == -5412 || mc == -5422 || mc == -5414 || mc == -5424 || mc == -5342 || mc == -5432 || mc == -5434 || mc == -5442 || mc == -5444 || mc == -5512 || mc == -5522 || mc == -5514 || mc == -5524 || mc == -5532 || mc == -5534 || mc == -5542 || mc == -5544 || mc == -5544 ) {return 1;}
    else {return 0;}
}

int ISB(int mc){
    if(mc == -511 || mc == -521 || mc == -10511 || mc == -10521 || mc == -513 || mc == -523 || mc == -10513 || mc == -10523 || mc == -20513 || mc == -20523 || mc == -515 || mc == -525 || mc == -531 || mc == -10531 || mc == -533 || mc == -10533 || mc == -20533 || mc == -535 || mc == -541 || mc == -10541 || mc == -543 || mc == -10543 || mc == -20543 || mc == -545 || mc == 5122 || mc == 5112 || mc == 5212 || mc == 5222 || mc == 5114 || mc == 5214 || mc == 5224 || mc == 5132 || mc == 5232 || mc == 5312 || mc == 5322 || mc == 5314 || mc == 5324 || mc == 5332 || mc == 5334 || mc == 5142 || mc == 5242 || mc == 5412 || mc == 5422 || mc == 5414 || mc == 5424 || mc == 5342 || mc == 5432 || mc == 5434 || mc == 5442 || mc == 5444 || mc == 5512 || mc == 5522 || mc == 5514 || mc == 5524 || mc == 5532 || mc == 5534 || mc == 5542 || mc == 5544 || mc == 5544 ) {return 1;}
    else {return 0;}

}


int ISC(int mc){
    if(mc == 411 || mc == 421 || mc == 10411 || mc == 10421 || mc == 413 || mc == 423 || mc == 10413 || mc == 10423 || mc == 20413 || mc == 20423 || mc == 415 || mc == 425 || mc == 431 || mc == 10431 || mc == 433 || mc == 10433 || mc == 20433 || mc == 435      || mc == 4122 || mc == 4222 || mc == 4212 || mc == 4112 || mc == 4224 || mc == 4214 || mc == 4114 || mc == 4232 || mc == 4132 || mc == 4322 || mc == 4312 || mc == 4324 || mc == 4314 || mc == 4332 || mc == 4334 || mc == 4412 || mc == 4422 || mc == 4414 || mc == 4424 || mc == 4432 || mc == 4434 || mc == 4444) {return 1;}
    else {return 0;}
}

int ISCbar(int mc){
    if(mc == -411 || mc == -421 || mc == -10411 || mc == -10421 || mc == -413 || mc == -423 || mc == -10413 || mc == -10423 || mc == -20413 || mc == -20423 || mc == -415 || mc == -425 || mc == -431 || mc == -10431 || mc == -433 || mc == -10433 || mc == -20433 || mc == -435      || mc == -4122 || mc == -4222 || mc == -4212 || mc == -4112 || mc == -4224 || mc == -4214 || mc == -4114 || mc == -4232 || mc == -4132 || mc == -4322 || mc == -4312 || mc == -4324 || mc == -4314 || mc == -4332 || mc == -4334 || mc == -4412 || mc == -4422 || mc == -4414 || mc == -4424 || mc == -4432 || mc == -4434 || mc == -4444) {return 1;}
    else {return 0;}
}

int ISS(int mc){
    if(mc == 130 || mc == 310 || mc == 311 || mc == 321 || mc == 10311 || mc == 10321 || mc == 100311 || mc == 100321 || mc == 200311 || mc == 200321 || mc == 9000311 || mc == 9000321 || mc == 313 || mc == 323 || mc == 10313 || mc == 10323 || mc == 20313 || mc == 20323 || mc == 100313 || mc == 100323 || mc == 9000313 || mc == 9000323 || mc == 30313 || mc == 30323 || mc == 315 || mc == 325 || mc == 9000315 || mc == 9000325 || mc == 10315 || mc == 10325 || mc == 20315 || mc == 20325 || mc == 100315 || mc == 100325 || mc == 9010315 || mc == 9010325 || mc == 317 || mc == 327 || mc == 9010317 || mc == 9010327 || mc == 319 || mc == 329 || mc == 9000319 || mc == 9000329       || mc == 3122 || mc == 3222 || mc == 3212 || mc == 3112 || mc == 3224 || mc == 3214 || mc == 3114 || mc == 3322 || mc == 3312 || mc == 3324 || mc == 3314 || mc == 3334) {return 1;}
    else {return 0;}
}






struct ParticleInfo {
    
    
  ParticleInfo(const GenParticle *particle) {
    /*
    pt = particle->PT;
    eta = particle->Eta;
    phi = particle->Phi;
    mass = particle->Mass;
    p4 = ROOT::Math::PtEtaPhiMVector(pt, eta, phi, mass);
    px = p4.px();
    py = p4.py();
    pz = p4.pz();
    energy = p4.energy();
    charge = particle->Charge;
    pid = particle->PID;
      
      if(charge != 0){
          if(abs(pid) == 11){isElectron = 1;}
          else if(abs(pid) == 13){isMuon = 1;}
          else{ isChargedHadron = 1; }
      }
      else{
          if(abs(pid) == 22){isPhoton = 1;}
          else{isNeutralHadron = 1;}
      }
      */
  }

    
  ParticleInfo(const ParticleFlowCandidate *particle) {
    pt = particle->PT;
    eta = particle->Eta;
    phi = particle->Phi;
    mass = particle->Mass;
    p4 = ROOT::Math::PtEtaPhiMVector(pt, eta, phi, mass);
    px = p4.px();
    py = p4.py();
    pz = p4.pz();
    energy = p4.energy();
    charge = particle->Charge;
    // pid = particle->PID;
    truth_pid = particle->truth_PID;
    d0 = particle->D0;
    d0err = particle->ErrorD0;
    dz = particle->DZ;
    dzerr = particle->ErrorDZ;
    
    // std::cerr<<"particle->PID : "<<particle->PID<<std::endl;
    MC = (GenParticle*) particle->Particles.At(0);
    pid = MC->PID;
    // std::cerr << "MC->PID: " << MC->PID << std::endl;

    //GenParticle *gPMother = gParticle->M1();
    
    //int MPID = gParticle->PID; 

    if(charge != 0){
      d0 = 0.2 * ((d0 > 0) ? 1 : -1) * std::pow(std::abs(d0), 0.25);
      d0err = 0.2 * ((d0err > 0) ? 1 : -1) * std::pow(std::abs(d0err), 0.25);
      dz = 0.2 * ((dz > 0) ? 1 : -1) * std::pow(std::abs(dz), 0.25);
      dzerr = 0.2 * ((dzerr > 0) ? 1 : -1) * std::pow(std::abs(dzerr), 0.25);
    }

      
      if(charge != 0){
          if(abs(truth_pid) == 11){isElectron = 1;}
          else if(abs(truth_pid) == 13){isMuon = 1;}
          else if(abs(truth_pid) == 211){isPion = 1; isChargedHadron = 1;}
          else if(abs(truth_pid) == 321){isChargedKaon = 1; isChargedHadron = 1;}
          else if(abs(truth_pid) == 2212){isProton = 1; isChargedHadron = 1;}
          else{ isChargedHadron = 1; }
      }
      else{
          if(abs(truth_pid) == 22){isPhoton = 1;}
          else{isNeutralHadron = 1;}
      }
      
      
      
  }

  GenParticle *MC = NULL;

  double pt;
  double eta;
  double phi;
  double mass;
  double px;
  double py;
  double pz;
  double energy;
  ROOT::Math::PtEtaPhiMVector p4;

  int charge;
  int pid;
  int truth_pid;

  float d0 = 0;
  float d0err = 0;
  float dz = 0;
  float dzerr = 0;
    
    bool isChargedHadron = false;
    bool isNeutralHadron = false;
    bool isPhoton = false;
    bool isElectron = false;
    bool isMuon = false;
    bool isPion = false;
    bool isChargedKaon = false;
    bool isProton = false;
    // bool isKLong = false;
    // bool isKShort = false;    

    
};

//------------------------------------------------------------------------------

void FILENAME() {
  gSystem->Load("libDelphes");

 //   OnnxReader* _reader = new OnnxReader(_onnx_file);

    std::vector<string> run_class; run_class.clear();
     run_class.push_back("RUNCLASS");
    //  run_class.push_back("cc"); 
    //  run_class.push_back("uu"); 
    // run_class.push_back("dd"); run_class.push_back("ss");   run_class.push_back("gg");
    for(int run = 0; run < run_class.size(); run++){
        string temp_run = run_class.at(run);
        for(int num_out = 0; num_out < 1; num_out++){
            
            //string outfile = "./rootFiles/higgs125/" + temp_run + "/" + temp_run + "_" + std::to_string(num_out) + ".root";
            string outfile = "OUTFILE";
            std::cerr << "outfile: " << outfile << std::endl;
            
            TFile *fout = new TFile(outfile.c_str(), "RECREATE");
            TTree *tree = new TTree("tree", "tree");
            
            
            
            std::map<TString, bool> labelVars;
            labelVars["label_bb"] = false;
            labelVars["label_cc"] = false;
            labelVars["label_gg"] = false;

            labelVars["label_b"] = false;
            labelVars["label_bbar"] = false;
            labelVars["label_c"] = false;
            labelVars["label_cbar"] = false;
            labelVars["label_u"] = false;
            labelVars["label_ubar"] = false;
            labelVars["label_d"] = false;
            labelVars["label_dbar"] = false;
            labelVars["label_s"] = false;
            labelVars["label_sbar"] = false;
            labelVars["label_g"] = false;
            
            
            
            // define branches
            std::map<TString, float> floatVars;
            floatVars["is_signal"] = 0;
            
            floatVars["gen_match"] = 0;
            floatVars["genpart_pt"] = 0;
            floatVars["genpart_eta"] = 0;
            floatVars["genpart_phi"] = 0;
            floatVars["genpart_pid"] = 0;
            
            floatVars["jet_pt"] = 0;
            floatVars["jet_eta"] = 0;
            floatVars["jet_phi"] = 0;
            floatVars["jet_energy"] = 0;
            floatVars["jet_nparticles"] = 0;
            floatVars["btag"] = 0;
            floatVars["ctag"] = 0;
            //  floatVars["jet_sdmass"] = 0;
            //  floatVars["jet_tau1"] = 0;
            //  floatVars["jet_tau2"] = 0;
            //  floatVars["jet_tau3"] = 0;
            //  floatVars["jet_tau4"] = 0;
            
            std::map<TString, std::vector<float>> arrayVars;
            arrayVars["part_px"];
            arrayVars["part_py"];
            arrayVars["part_pz"];
            arrayVars["part_energy"];
            arrayVars["part_pt"];
            arrayVars["part_deta"];
            arrayVars["part_dphi"];
            arrayVars["part_charge"];
            arrayVars["part_pid"];
            arrayVars["part_d0val"];
            arrayVars["part_d0err"];
            arrayVars["part_dzval"];
            arrayVars["part_dzerr"];
            arrayVars["part_deltaR"];
            arrayVars["part_logptrel"];
            arrayVars["part_logerel"];
            arrayVars["part_e_log"];
            arrayVars["part_pt_log"];
            arrayVars["part_d0"];
            arrayVars["part_dz"];
            arrayVars["part_vtxX"];
            arrayVars["part_vtxY"];
            arrayVars["part_vtxZ"];



         
            
            
            std::map<TString, std::vector<bool>> arrayBools;
            arrayBools["part_isChargedHadron"];
            arrayBools["part_isNeutralHadron"];
            arrayBools["part_isPhoton"];
            arrayBools["part_isElectron"];
            arrayBools["part_isMuon"];
            arrayBools["part_isPion"];
            arrayBools["part_isProton"];
            arrayBools["part_isChargedKaon"];
            arrayBools["part_isKLong"];
            arrayBools["part_isKShort"];   
            arrayBools["part_isPi0"];

         
            // book
            for (auto &v : floatVars) {
                tree->Branch(v.first.Data(), &v.second);
            }
            for (auto &v : arrayVars) {
                tree->Branch(v.first.Data(), &v.second, /*bufsize=*/1024000);
            }
            for (auto &v : arrayBools) {
                tree->Branch(v.first.Data(), &v.second);
            }
            for (auto &v : labelVars) {
                tree->Branch(v.first.Data(), &v.second);
            }
            
            
            
            std::vector<string> infiles; infiles.clear();
            //for(int i = num_out*25+1; i< (num_out+1)*25 +1; i++){
            for(int i = num_out; i< (num_out+1); i++){
                //infiles.push_back( "/afs/ihep.ac.cn/users/z/zhuyf/cefs/workspace/Delphes/Pythia6/" + temp_run + "/" + to_string(i) + ".root" );
                infiles.push_back( "INFILE" );
            }

//            infiles.push_back( "/cefs/higgs/PFAData/Software/Delphes/Data/Baseline/" + temp_run + "/" + to_string(num_out) + ".root" );
//            infiles.push_back( "/afs/ihep.ac.cn/users/z/zhuyf/cefs/workspace/Delphes/Pythia8/" + temp_run + "/" + to_string(num_out) + ".root" );

            // read input
            TChain *chain = new TChain("Delphes");
            //  chain->Add(inputFile);
            for(int i = 0; i<infiles.size(); i++){
                //        chain->Add( Form("./%s.root", names[i].Data()) );
                //        chain->Add( "./Zpole_bb.e0.p0.whizard195/1.root" );
                chain->Add( (infiles.at(i)).c_str() );
            }
            ExRootTreeReader *treeReader = new ExRootTreeReader(chain);
            Long64_t allEntries = treeReader->GetEntries();
            
            //  std::cerr << "** Input file: " << inputFile << std::endl;
            //  std::cerr << "** Jet branch: " << jetBranch << std::endl;
            std::cerr << "** Total events: " << allEntries << std::endl;
            
            // analyze
            TClonesArray *branchParticle = treeReader->UseBranch("Particle");
            TClonesArray *branchPFCand = treeReader->UseBranch("ParticleFlowCandidate");
            TClonesArray *branchJet = treeReader->UseBranch("Jet");
            
            //  FatJetMatching fjmatch(0.8);
            
            
            
            
            
            // Loop over all events
            int num_processed = 0;
            for (Long64_t entry = 0; entry < allEntries; ++entry) {
                if (entry % 1000 == 0) {
                  //  std::cerr << "processing " << entry << " of " << allEntries << " events." << std::endl;
                }
                
                // Load selected branches with data from specified event
                treeReader->ReadEntry(entry);

                std::vector<const GenParticle *> genParticles = getLabel(branchParticle);



                std::vector<const GenParticle *> BHadrons; BHadrons.clear();
                const GenParticle* hDau1 = NULL;
                const GenParticle* hDau2 = NULL;
                std::map<const GenParticle*, TVector3 > mcp_vtx; mcp_vtx.clear();
                for (Int_t i = 0; i < branchParticle->GetEntriesFast(); ++i) {
                    const GenParticle* particle = (GenParticle*)branchParticle->At(i);

                    if(particle->Status == 1 && particle->Charge != 0 && particle->M1 > 0){
                        int ip_mom = particle->M1;
                        GenParticle* p_mom = (GenParticle*)branchParticle->At(ip_mom);
                        int num_cdau = 0;
                        for(Int_t j = p_mom->D1; j < p_mom->D2; j++){
                            GenParticle* pdau = (GenParticle*)branchParticle->At(j);
                            if(pdau->Charge != 0 && pdau->Status == 1){num_cdau += 1;}
                        }
                        if( num_cdau > 1 ){
                        //if(num_cdau > 1 && (ISB(abs(p_mom->PID))==1 || ISC(abs(p_mom->PID))==1  ) ){
                            TVector3 TVp(particle->X, particle->Y, particle->Z);
                            mcp_vtx[p_mom] = TVp;
                        }
                    }


                    if (abs(particle->PID) == 25) {
                        auto hdaus = getDaughterQuarks(genParticles, particle);
                        if (hdaus.size() == 2){
                            hDau1 = hdaus.at(0);
                            hDau2 = hdaus.at(1);
                        }
                    }

		            if ( ISBbar(particle->PID) == 1 || ISB(particle->PID) == 1 ){
			            BHadrons.push_back(particle);
		            }


                }

                bool reverse = false;
                int jet1PID = 0;
                int jet2PID = 0;
                if( branchJet->GetEntriesFast() == 2 && hDau1 != NULL && hDau2 != NULL){
                //  std::cerr << "hDau1->PID " << hDau1->PID << " hDau2->PID " << hDau2->PID << std::endl;
                    const Jet *jet1 = (Jet *)branchJet->At(0);
                    const Jet *jet2 = (Jet *)branchJet->At(1);
                    double pt1 = jet1->PT; double eta1 = jet1->Eta; double phi1 = jet1->Phi; double mass1 = jet1->Mass;
                    ROOT::Math::PtEtaPhiMVector p41 = ROOT::Math::PtEtaPhiMVector(pt1, eta1, phi1, mass1);
                    double px1 = p41.px(); double py1 = p41.py(); double pz1 = p41.pz(); double energy1 = p41.energy();
                    double pt2 = jet2->PT; double eta2 = jet2->Eta; double phi2 = jet2->Phi; double mass2 = jet2->Mass;
                    ROOT::Math::PtEtaPhiMVector p42 = ROOT::Math::PtEtaPhiMVector(pt2, eta2, phi2, mass2);
                    double px2 = p42.px(); double py2 = p42.py(); double pz2 = p42.pz(); double energy2 = p42.energy();

                    TLorentzVector TLjet1(px1, py1, pz1, energy1);
                    TLorentzVector TLjet2(px2, py2, pz2, energy2);
                    TLorentzVector TLhDau1(hDau1->Px, hDau1->Py, hDau1->Pz, hDau1->E);
                    TLorentzVector TLhDau2(hDau2->Px, hDau2->Py, hDau2->Pz, hDau2->E);
                    TVector3 TVjet1 = TLjet1.Vect();
                    TVector3 TVjet2 = TLjet2.Vect();
                    TVector3 TVhDau1 = TLhDau1.Vect();
                    TVector3 TVhDau2 = TLhDau2.Vect();
                    if(TVjet1.Angle(TVhDau1) + TVjet2.Angle(TVhDau2) > TVjet1.Angle(TVhDau2) + TVjet2.Angle(TVhDau1)){
                        jet1PID = hDau2->PID;
                        jet2PID = hDau1->PID;
                    }
                    else{
                        jet1PID = hDau1->PID;
                        jet2PID = hDau2->PID;
                    }

                }


/*
	        if( BHadrons.size() == 2 && branchJet->GetEntriesFast() == 2 ){
		  const Jet *jet1 = (Jet *)branchJet->At(0);
                  const Jet *jet2 = (Jet *)branchJet->At(1);
                  const GenParticle* B1 = BHadrons.at(0);
                  const GenParticle* B2 = BHadrons.at(1);
                  TLorentzVector TLjet1(0,0,0,0);
                  TLorentzVector TLjet2(0,0,0,0);
                  TLorentzVector TLB1(0,0,0,0);
                  TLorentzVector TLB2(0,0,0,0);


		  TLjet1 = jet1->P4();
                  TLjet2 = jet2->P4();
                  TLB1   = B1->P4();
                  TLB2   = B2->P4();
	
		  TVector3 TVjet1 = TLjet1.Vect();
                  TVector3 TVjet2 = TLjet2.Vect();
                  TVector3 TVB1   = TLB1.Vect();
                  TVector3 TVB2   = TLB2.Vect();

		  if(TVjet1.Angle(TVB1) + TVjet2.Angle(TVB2) > TVjet1.Angle(TVB2) + TVjet2.Angle(TVB1)){
		     

		  }

		}
        */
                
                // Loop over all jets in event
                int count_jet = 0;
                for (Int_t i = 0; i < branchJet->GetEntriesFast(); ++i) {
                    const Jet *jet = (Jet *)branchJet->At(i);
                    
                    //      if (jet->PT < 500 || std::abs(jet->Eta) > 2) not useful at the CEPC
                    //        continue;
                    
                    for (auto &v : labelVars) {
                        v.second = false;
                    }
                    for (auto &v : floatVars) {
                        v.second = 0;
                    }
                    for (auto &v : arrayVars) {
                        v.second.clear();
                    }
                    for (auto &v : arrayBools) {
                        v.second.clear();
                    }
                    
                  
                    floatVars["btag"] = jet->BTag;
                    floatVars["ctag"] = 0.2;
                    floatVars["jet_pt"] = jet->PT;
                    floatVars["jet_eta"] = jet->Eta;
                    floatVars["jet_phi"] = jet->Phi;
                    floatVars["jet_energy"] = jet->P4().Energy();
                    
                    if(temp_run == "bb"){labelVars["label_bb"] = true; labelVars["label_cc"] = false; labelVars["label_gg"] = false;}
                    else if(temp_run == "cc"){labelVars["label_bb"] = false; labelVars["label_cc"] = true; labelVars["label_gg"] = false;}
                    else{labelVars["label_bb"] = false; labelVars["label_cc"] = false; labelVars["label_gg"] = true;}

                   

                    if(temp_run == "bb"){
                      if(i == 0 && jet1PID > 0){labelVars["label_b"] = true;}
                      else if(i == 0 && jet1PID < 0){labelVars["label_bbar"] = true;}
                      else if(i == 1 && jet2PID > 0){labelVars["label_b"] = true;}
                      else if(i == 1 && jet2PID < 0){labelVars["label_bbar"] = true;}
                     }
                    else if(temp_run == "cc"){
                      if(i == 0 && jet1PID > 0){labelVars["label_c"] = true;}
                      else if(i == 0 && jet1PID < 0){labelVars["label_cbar"] = true;}
                      else if(i == 1 && jet2PID > 0){labelVars["label_c"] = true;}
                      else if(i == 1 && jet2PID < 0){labelVars["label_cbar"] = true;}
                    }
                    else if(temp_run == "uu"){
                      if(i == 0 && jet1PID > 0){labelVars["label_u"] = true;}
                      else if(i == 0 && jet1PID < 0){labelVars["label_ubar"] = true;}
                      else if(i == 1 && jet2PID > 0){labelVars["label_u"] = true;}
                      else if(i == 1 && jet2PID < 0){labelVars["label_ubar"] = true;}
                    }
                    else if(temp_run == "dd"){
                      if(i == 0 && jet1PID > 0){labelVars["label_d"] = true;}
                      else if(i == 0 && jet1PID < 0){labelVars["label_dbar"] = true;}
                      else if(i == 1 && jet2PID > 0){labelVars["label_d"] = true;}
                      else if(i == 1 && jet2PID < 0){labelVars["label_dbar"] = true;}
                    }
                    else if(temp_run == "ss"){
                      if(i == 0 && jet1PID > 0){labelVars["label_s"] = true;}
                      else if(i == 0 && jet1PID < 0){labelVars["label_sbar"] = true;}
                      else if(i == 1 && jet2PID > 0){labelVars["label_s"] = true;}
                      else if(i == 1 && jet2PID < 0){labelVars["label_sbar"] = true;}
                    }
                    else if(temp_run == "gg"){
                      labelVars["label_g"] = true;
                    }

                    // Loop over all jet's constituents
                    std::vector<ParticleInfo> particles;
                    for (Int_t j = 0; j < jet->Constituents.GetEntriesFast(); ++j) {
                        const TObject *object = jet->Constituents.At(j);
                        
                        // Check if the constituent is accessible
                        if (!object)
                            continue;
                        
                        if (object->IsA() == GenParticle::Class()) {
                            particles.emplace_back((GenParticle *)object);
                        } else if (object->IsA() == ParticleFlowCandidate::Class()) {
                            particles.emplace_back((ParticleFlowCandidate *)object);
                        }

                    }
                    
                    // sort particles by pt
                    std::sort(particles.begin(), particles.end(), [](const auto &a, const auto &b) { return a.energy > b.energy; });
                    floatVars["jet_nparticles"] = particles.size();
                    for (const auto &p : particles) {

                        GenParticle* gen = p.MC;

                        //std::cerr<<"gen vertex : "<<gen->X<<" "<<gen->Y<<" "<<gen->Z<<endl;
                        //std::cerr<<"gen->PID : "<<gen->PID<<" gen->D1 : "<<gen->D1<<" gen->Status : "<<gen->Status<<std::endl;
                        // if(std::abs(gen->PID) == 130){std::cerr<<"KLong not decay"<<std::endl;}
                        // std::cerr<<"gen->PID : "<<gen->PID<<std::endl;
                        int MPID = 999;
                        //std::cerr<<"gen->M1 "<<gen->M1<<" "<<gen->M2<<std::endl;
                        int iMom = gen->M1;
                        GenParticle* Mom = (GenParticle*)branchParticle->At(iMom);
                        //const auto *Mom = genParticles.at(iMom);
                        MPID = std::abs(Mom->PID);
                        // std::cerr << "MPID: " << MPID << std::endl;
                        int MMPID = 999;
                        
                        if(Mom->M1 > 0){
                            //std::cerr<<"Mom->M1 : "<<Mom->M1<<"  "<<Mom->DecayPosition<<std::endl;
                            //std::cerr<<"Mom->Status : "<<Mom->Status<<std::endl;
                            int iMMom = Mom->M1;
                        //  std::cerr << "Mom->M1: " << Mom->M1 << std::endl;
                            GenParticle* MMom = (GenParticle*)branchParticle->At(iMMom);
                            MMPID = std::abs(MMom->PID);

                            // int num_charge_Daus = 0;
                            // for(Int_t k = MMom->D1; k <= MMom->D2; k++){
                            //     GenParticle* MDau = (GenParticle*)branchParticle->At(k);
                            //     if(MDau->Charge != 0){num_charge_Daus += 1;}
                            // }
                            auto it = mcp_vtx.find(Mom);
                            if(it != mcp_vtx.end()){ 
                                TVector3 v3 = mcp_vtx[Mom];
                                //std::cerr<<" "<<v3(0)<<" "<<v3(1)<<" "<<v3(2)<<endl;
                                arrayVars["part_vtxX"].push_back(v3(0)); arrayVars["part_vtxY"].push_back(v3(1)); arrayVars["part_vtxZ"].push_back(v3(2)); 
                            }
                            else{
                                arrayVars["part_vtxX"].push_back(0); arrayVars["part_vtxY"].push_back(0); arrayVars["part_vtxZ"].push_back(0);
                            }
                         // std::cerr<<"MMPID : "<<MMom->PID<<endl;
                        }
                        else{
                            arrayVars["part_vtxX"].push_back(0); arrayVars["part_vtxY"].push_back(0); arrayVars["part_vtxZ"].push_back(0);
                        }
                        

                        
                        

                        if(gen != NULL && std::abs(gen->PID) == 130){arrayBools["part_isKLong"].push_back(true); }
                        else{ arrayBools["part_isKLong"].push_back(false); }
                        if(gen != NULL && (MPID == 310 || MMPID == 310)){arrayBools["part_isKShort"].push_back(true); }
                        else{ arrayBools["part_isKShort"].push_back(false); }
                        if(gen != NULL && ( std::abs(p.pid) == 111 || std::abs(MPID) == 111 || MMPID == 111 ) ){
                          arrayBools["part_isPi0"].push_back(true);
                        }
                        else{ arrayBools["part_isPi0"].push_back(false); }
                        // else{ std::cerr<<"not KLong or KShort"<<std::endl; }



                        

                        arrayVars["part_px"].push_back(p.px);
                        arrayVars["part_py"].push_back(p.py);
                        arrayVars["part_pz"].push_back(p.pz);
                        arrayVars["part_energy"].push_back(p.energy);
                        arrayVars["part_pt"].push_back(p.pt);
                        arrayVars["part_deta"].push_back((jet->Eta > 0 ? 1 : -1) * (p.eta - jet->Eta));
                        arrayVars["part_dphi"].push_back(deltaPhi(p.phi, jet->Phi));
                        arrayVars["part_charge"].push_back(p.charge);
                        arrayVars["part_d0val"].push_back(p.d0);
                        arrayVars["part_d0err"].push_back(p.d0err);
                        arrayVars["part_dzval"].push_back(p.dz);
                        arrayVars["part_dzerr"].push_back(p.dzerr);
                        
                        arrayBools["part_isChargedHadron"].push_back(p.isChargedHadron);
                        arrayBools["part_isNeutralHadron"].push_back(p.isNeutralHadron);
                        arrayBools["part_isPhoton"].push_back(p.isPhoton);
                        arrayBools["part_isElectron"].push_back(p.isElectron);
                        arrayBools["part_isMuon"].push_back(p.isMuon);
                        arrayBools["part_isPion"].push_back(p.isPion);
                        arrayBools["part_isChargedKaon"].push_back(p.isChargedKaon);
                        arrayBools["part_isProton"].push_back(p.isProton);

                        arrayVars["part_deltaR"].push_back( std::hypot(  ((jet->Eta > 0 ? 1 : -1) * (p.eta - jet->Eta)),  deltaPhi(p.phi, jet->Phi)  )  );
                        arrayVars["part_logptrel"].push_back( std::log( p.pt ) );
                        arrayVars["part_logerel"].push_back( std::log( p.energy ) );
                        arrayVars["part_e_log"].push_back( std::log( p.energy/ jet->P4().Energy() ) );
                        arrayVars["part_pt_log"].push_back( std::log( p.pt/ jet->PT ) );
                        arrayVars["part_d0"].push_back( std::tanh(p.d0) );
                        arrayVars["part_dz"].push_back( std::tanh(p.dz) );
			                  
                        
                    }
                    
                    if( (labelVars["label_bbar"] || labelVars["label_cbar"] || labelVars["label_dbar"] || labelVars["label_ubar"] || labelVars["label_sbar"] || labelVars["label_g"]) && count_jet < 1 ){
                      count_jet += 1;
                      tree->Fill();
                    }


/*
                    std::map<std::string, std::vector<float>> map_input = {
                        {"part_pt_log",          Preprocessing(arrayVars["part_pt_log"], 50, -1.5, 1.0)},
                        {"part_e_log",           Preprocessing(arrayVars["part_e_log"], 50, -0.687, 1.0 )},
                        {"part_logptrel",        Preprocessing(arrayVars["part_logptrel"], 50, -4.7, 1.0)},
                        {"part_logerel",         Preprocessing(arrayVars["part_logerel"], 50, -4.473, 1.0)},
                        {"part_deltaR",          Preprocessing(arrayVars["part_deltaR"], 50, 2.1, 2.3)},
                        {"part_charge",          Preprocessing(arrayVars["part_charge"], 50)},
                        {"part_isChargedKaon",   Preprocessing(arrayVars["part_isChargedKaon"], 50)},
                        {"part_isPion",          Preprocessing(arrayVars["part_isPion"], 50)},
                        {"part_isProton",        Preprocessing(arrayVars["part_isProton"], 50)},
                        {"part_isElectron",      Preprocessing(arrayVars["part_isElectron"], 50)},
                        {"part_isMuon",          Preprocessing(arrayVars["part_isMuon"], 50)},
                        {"part_isNeutralHadron", Preprocessing(arrayVars["part_isNeutralHadron"], 50)},
                        {"part_isPhoton",        Preprocessing(arrayVars["part_isPhoton"], 50)},
                        {"part_d0",              Preprocessing(arrayVars["part_d0"], 50)},
                        {"part_dz",              Preprocessing(arrayVars["part_dz"], 50)},
                        {"part_deta",            Preprocessing(arrayVars["part_deta"], 50)},
                        {"part_dphi",            Preprocessing(arrayVars["part_dphi"], 50)},
                        {"part_mask",            std::vector<float>(50, 1.)}
                    };

                    std::vector<std::string> pf_points_name = {"part_dphi", "part_deta"};
                    std::vector<std::string> pf_features_name = {
                        "part_pt_log", "part_e_log", "part_logptrel", "part_logerel",
                        "part_deltaR", "part_charge", "part_isChargedKaon", "part_isPion",
                        "part_isProton", "part_isElectron", "part_isMuon", "part_isNeutralHadron",
                        "part_isPhoton", "part_d0", "part_dz", "part_deta",
                        "part_dphi"
                    };

                    std::vector<float> vec_pf_points;
                    for (const auto& k : pf_points_name) {
                        auto kv = map_input.find(k);
                        std::vector<float> a_feat;
                        if (kv == map_input.end()) {
                            std::cout << "no input: " << k << std::endl;
                            a_feat = std::vector<float>(50, 0.);
                        }
                        else {
                            a_feat = kv->second;
                        }

                        vec_pf_points.insert(vec_pf_points.end(), a_feat.begin(), a_feat.end());
                    }


                    std::vector<float> vec_pf_features;
                    for (const auto& k : pf_features_name) {
                        auto kv = map_input.find(k);
                        std::vector<float> a_feat;
                        if (kv == map_input.end()) {
                            std::cout << "no input: " << k << std::endl;
                            a_feat = std::vector<float>(50, 0.);
                        }
                        else {
                            a_feat = kv->second;
                        }

                        vec_pf_features.insert(vec_pf_features.end(), a_feat.begin(), a_feat.end());
                    }

                    std::vector<std::string> names = {
                        "pf_points", "pf_features", "pf_mask"
                    };
                    std::vector<std::vector<float>> input = {
                        vec_pf_points, vec_pf_features, std::vector<float>(50, 1.)
                    };
                    std::vector<std::vector<int64_t>> shapes = {
                        {1, 2, 50}, {1, 17, 50}, {1, 1, 50}
                    };

*/
                    // std::vector<float> output = reader.GetPrediction(names, input, shapes);

                    // for (const auto& val : output) {
                    //     std::cout << val << " ";
                    // }
                    // std::cout << std::endl;


                   // tree->Fill();
                    ++num_processed;
                }
            }
            
            tree->Write();
            //  std::cerr << TString::Format("** Written %d jets to output %s", num_processed, outputFile.Data()) << std::endl;
            
            delete treeReader;
            delete chain;
            delete fout;
            
            
            
        }
    }
    
    
}

//------------------------------------------------------------------------------
