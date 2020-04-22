#include <AsgTools/MessageCheck.h>
#include <LowPtAnalysis/TruthAnalysis.h>

#include <xAODEventInfo/EventInfo.h>

#include <xAODTruth/TruthParticleContainer.h>
#include <xAODTruth/xAODTruthHelpers.h>
#include <xAODTracking/TrackParticleContainer.h>

#include <cmath>

//ROOT includes
#include "TFile.h"

TruthAnalysis :: TruthAnalysis (const std::string& name,
                                  ISvcLocator *pSvcLocator)
    : EL::AnaAlgorithm (name, pSvcLocator)
{
  // Here you put any code for the base initialization of variables,
  // e.g. initialize all pointers to 0.  This is also where you
  // declare all properties for your algorithm.  Note that things like
  // resetting statistics variables or booking histograms should
  // rather go into the initialize() function.

  // Properties
  declareProperty("lep1_min_pt", lep1_min_pt = 27*GeV);
  declareProperty("lep2_min_pt", lep2_min_pt = 20*GeV);
  declareProperty("lep_max_eta", lep_max_eta = 2.5);
  declareProperty("dilep_min_pt", dilep_min_pt = 30.0*GeV);
  declareProperty("dilep_min_mass", dilep_min_mass = 20.0*GeV);
  declareProperty("tracks_min_pt", tracks_min_pt = 400*MeV,
		  "minimum reconstructible pT");
  declareProperty("tracks_max_eta", tracks_max_eta = 2.5,
		  "maximum reconstructible eta");
  declareProperty("tracks_max_n", tracks_max_n = 0,
		  "exclusivity selection");
  declareProperty("input_trk_eff_file", input_trk_eff_file = "",
		  "tracking efficiency input. if empty, use hard-coded numbers");
 declareProperty("input_pu_file", input_pu_file = "",
		 "PU info input. if empty, use hard-coded numbers");
 declareProperty("input_weight_file", input_weight_file = "/global/cfs/projectdirs/atlas/wmccorma/ExclWW/trk-exclusive-ww/analysis/weight_plots.root",
		 "weight info input. if empty, use hard-coded numbers");  
 declareProperty("filter_by_selections", filter_by_selections = false,
		  "if true, only events passing selections are stored in the output");
declareProperty("do_lpair_reweight", do_lpair_reweight = 0,
		  "if true, reweight events to match LPAIR ntrack distribution");
  declareProperty("random_seed", random_seed=29901,
		  "Random seed for tracking efficiency");

}



StatusCode TruthAnalysis :: initialize ()
{
  // Here you do everything that needs to be done at the very
  // beginning on each worker node, e.g. create histograms and output
  // trees.  This method gets called before any input files are
  // connected.

  ANA_MSG_INFO ("Initializing");

  // Random generator
  m_rnd = new TRandom3(random_seed);

  // Output tree/histograms
  ANA_CHECK (book (TTree ("analysis", "Truth analysis ntuple")));
  TTree* mytree = tree ("analysis");

  //Event-level
  mytree->Branch ("dilep_pt", &m_dilep_pt);
  mytree->Branch ("dilep_m", &m_dilep_m);
  mytree->Branch ("WW_pt", &m_WW_pt);
  mytree->Branch ("WW_m", &m_WW_m);
  mytree->Branch ("diphoton_m", &m_diphoton_m);
  m_pass_sel = new std::vector<char>(ncuts);
  mytree->Branch("pass_sel", m_pass_sel);
  mytree->Branch ("numPUtracks", &m_numPUtracks);
  mytree->Branch ("numHStracks", &m_numHStracks);
  mytree->Branch ("numHStracks2", &m_numHStracks2);
  mytree->Branch ("mcWeight", &m_mcWeight);
  mytree->Branch ("nTrackweight", &m_nTrackweight);
  //mytree->Branch ("mcWeight_0trk", &m_mcWeight_0trk);
  //mytree->Branch ("mcWeight_1trk", &m_mcWeight_1trk);
  //mytree->Branch ("mcWeight_2trk", &m_mcWeight_2trk);
  //mytree->Branch ("mcWeight_3trk", &m_mcWeight_3trk);
  //mytree->Branch ("rand", &m_rand);

  //Leptons
  m_lep_pt = new std::vector<float>();
  mytree->Branch ("lep_pt", &m_lep_pt);
  m_lep_eta = new std::vector<float>();
  mytree->Branch ("lep_eta", &m_lep_eta);
  m_lep_phi = new std::vector<float>();
  mytree->Branch ("lep_phi", &m_lep_phi);
  m_lep_charge = new std::vector<int>();
  mytree->Branch ("lep_charge", &m_lep_charge);
  m_lep_pdgid = new std::vector<int>();
  mytree->Branch ("lep_pdgid", &m_lep_pdgid);

  //photons
  m_photon_p = new std::vector<float>();
  mytree->Branch ("photon_p", &m_photon_p);

  m_W_pt = new std::vector<float>();
  mytree->Branch ("W_pt", &m_W_pt);
  m_W_m = new std::vector<float>();
  mytree->Branch ("W_m", &m_W_m);
  m_W_eta = new std::vector<float>();
  mytree->Branch ("W_eta", &m_W_eta);
  m_W_phi = new std::vector<float>();
  mytree->Branch ("W_phi", &m_W_phi);

  //Tracks
  m_trk_pt = new std::vector<float>();
  mytree->Branch ("trk_pt", &m_trk_pt);
  m_trk_eta = new std::vector<float>();
  mytree->Branch ("trk_eta", &m_trk_eta);
  m_trk_phi = new std::vector<float>();
  mytree->Branch ("trk_phi", &m_trk_phi);
  m_trk_charge = new std::vector<int>();
  mytree->Branch ("trk_charge", &m_trk_charge);
  m_trk_pdgid = new std::vector<int>();
  mytree->Branch ("trk_pdgid", &m_trk_pdgid);

  //////TRUTH STUFF no pseudo-reco
  mytree->Branch ("truedilep_pt", &m_truedilep_pt);
  mytree->Branch ("truedilep_m", &m_truedilep_m);
  mytree->Branch ("numTrueHSparticles", &m_numTrueHSparticles);
  m_truelep_pt = new std::vector<float>();
  mytree->Branch ("truelep_pt", &m_truelep_pt);
  m_truelep_eta = new std::vector<float>();
  mytree->Branch ("truelep_eta", &m_truelep_eta);
  m_truelep_phi = new std::vector<float>();
  mytree->Branch ("truelep_phi", &m_truelep_phi);
  m_truelep_charge = new std::vector<int>();
  mytree->Branch ("truelep_charge", &m_truelep_charge);
  m_truelep_pdgid = new std::vector<int>();
  mytree->Branch ("truelep_pdgid", &m_truelep_pdgid);
  m_truetrk_pt = new std::vector<float>();
  mytree->Branch ("truetrk_pt", &m_truetrk_pt);
  m_truetrk_eta = new std::vector<float>();
  mytree->Branch ("truetrk_eta", &m_truetrk_eta);
  m_truetrk_phi = new std::vector<float>();
  mytree->Branch ("truetrk_phi", &m_truetrk_phi);
  m_truetrk_charge = new std::vector<int>();
  mytree->Branch ("truetrk_charge", &m_truetrk_charge);
  m_truetrk_pdgid = new std::vector<int>();
  mytree->Branch ("truetrk_pdgid", &m_truetrk_pdgid);
  mytree->Branch ("true_twoleps", &m_true_twoleps);
  mytree->Branch ("true_twoOSleps", &m_true_twoOSleps);
  mytree->Branch ("true_absproduct_of_lepton_pdgidvals", &m_true_absproduct_of_lepton_pdgidvals);
  mytree->Branch ("true_leadlep_ptpass", &m_true_leadlep_ptpass);
  mytree->Branch ("true_subleadlep_ptpass", &m_true_subleadlep_ptpass);
  mytree->Branch ("true_dilep_masspass", &m_true_dilep_masspass);
  mytree->Branch ("true_dilep_ptpass", &m_true_dilep_ptpass);
  mytree->Branch ("true_numCh_and_PU", &m_true_numCh_and_PU);


  ANA_CHECK (book ( TH1F ("cutflow", "Cutflow", ncuts, -0.5, ncuts-0.5) ));
  ANA_CHECK (book ( TH1F ("num_fiducial_leptons", "Number of fiducial leptons", 10, 0, 10) ));
  ANA_CHECK (book ( TH1F ("dilep_m", "m(e#mu) (GeV)", 60, 0, 300.) ));
  ANA_CHECK (book ( TH1F ("dilep_pt", "p_{T} (e#mu) (GeV)", 60, 0, 300.) ));
  ANA_CHECK (book ( TH1F ("num_fiducial_tracks", "Number of fiducial tracks", 50, 0, 50) ));
  ANA_CHECK (book ( TH1F ("sr_dilep_pt", "p_{T} (e#mu) after all selections (GeV);p_{T}(e#mu) [GeV];Events/5 GeV", 60, 0, 300.) ));

  ANA_CHECK (book ( TH1D ("hCutFlow_Sum", "Cut Flow Sum", 3, 0, 3.) ));

  ANA_CHECK (book ( TH1F ("num_hs_tracks", "Number of HS Reco Tracks", 20, 0, 20) ));
  ANA_CHECK (book ( TH1F ("num_hs_truth_tracks", "Number of HS Truth Tracks", 20, 0, 20) ));
  ANA_CHECK (book ( TH1F ("num_hs_tracksRW", "Number of RW HS Reco Tracks", 20, 0, 20) ));

  std::vector<int> ptll_bins = {0,5,10,15,20,25,30,50,70,90,110,130};
  std::vector<int> mass_bins = {44, 60, 90, 110, 140, 180, 190, 200, 210, 220, 230, 250, 270, 290, 320, 350, 380, 410, 450, 500, 600};
  for(unsigned int i = 0; i<ptll_bins.size() - 1; i++){
    std::stringstream ss1;
    ss1<<ptll_bins.at(i);
    std::stringstream ss2;
    ss2<<ptll_bins.at(i+1);
    std::string s1;
    std::string s2;
    ss1>>s1;
    ss2>>s2;
    ANA_CHECK (book ( TH1F ( TString("num_hs_tracks_ptll_"+s1+"_"+s2), TString("Number of HS Reco Tracks for ptll "+s1+"_"+s2+"; N_{HStracks}^{PV}; Entries"), 20, 0, 20) ));
  }
    for(unsigned int i = 0; i < mass_bins.size() - 1; i++){
      std::stringstream ss1;
      ss1<<mass_bins.at(i);
      std::stringstream ss2;
      ss2<<mass_bins.at(i+1);
      std::string s1;
      std::string s2;
      ss1>>s1;
      ss2>>s2;
    ANA_CHECK (book ( TH1F ( TString("num_hs_tracks_mll_"+s1+"_"+s2), TString("Number of HS Reco Tracks for mll "+s1+"_"+s2+"; N_{HStracks}^{PV}; Entries"), 20, 0, 20) ));
    }


  //For Data/MC comparison similar to Filip's offline code
  MakeHistos();

  //set bin labels for cutflow
  for (int i=1; i<=ncuts;i++) {
    hist("cutflow")->GetXaxis()->SetBinLabel(i, cuts_labels[i-1].c_str());
  }

  //retrieve tracking efficiency, if needed
  h_trk_eff_pt_eta = nullptr;
  h_electron_eff= nullptr;
  h_muon_eff = nullptr;
  if (not input_trk_eff_file.empty()) {
    TFile *f_trk_eff = TFile::Open(input_trk_eff_file.c_str());
    h_trk_eff_pt_eta = static_cast<TProfile2D*>(f_trk_eff->Get("Tracking_Eff_2D"));
    h_electron_eff = static_cast<TProfile2D*>(f_trk_eff->Get("Electron_Eff_2D"));
    h_muon_eff = static_cast<TProfile2D*>(f_trk_eff->Get("Muon_Eff_2D"));
    if (h_trk_eff_pt_eta == nullptr) {
      ANA_MSG_ERROR("Error loading tracking efficiency from:" << input_trk_eff_file);
      return StatusCode::FAILURE;
    }
    if (h_electron_eff == nullptr) {
      ANA_MSG_ERROR("Error loading electron efficiency from:" << input_trk_eff_file);
      return StatusCode::FAILURE;
    }
    if (h_muon_eff == nullptr) {
      ANA_MSG_ERROR("Error loading muon efficiency from:" << input_trk_eff_file);
      return StatusCode::FAILURE;
    }
    ANA_MSG_INFO("Loaded tracking efficiency from " << input_trk_eff_file);
  }
  //retrieve tracking efficiency, if needed
  h_pu_info = nullptr;
  if (not input_pu_file.empty()) {
    TFile *f_pu = TFile::Open(input_pu_file.c_str());
    //h_pu_info = static_cast<TH1F*>(f_pu->Get("PUC_mumuOSEWW_MassZZmumu_PP8"));
    h_pu_info = static_cast<TH1F*>(f_pu->Get("PUC_mumuOSEWW_MassZdataFULL"));
    if (h_pu_info == nullptr) {
      ANA_MSG_ERROR("Error loading pu info from:" << input_pu_file);
      return StatusCode::FAILURE;
    }
    ANA_MSG_INFO("Loaded pu info from " << input_pu_file);
  }

  //h_weight_info = nullptr;
  if (not input_weight_file.empty()) {
    TFile *f_weight = TFile::Open(input_weight_file.c_str());
    std::vector<int> mass_bins = {44, 60, 90, 110, 140, 180, 190, 200, 210, 220, 230, 250, 270, 290, 320, 350, 380, 410, 450, 500, 600};
    for(unsigned int i = 0; i < mass_bins.size() - 1; i++){
      std::stringstream ss1;
      ss1<<mass_bins.at(i);
      std::stringstream ss2;
      ss2<<mass_bins.at(i+1);
      std::string s1;
      std::string s2;
      ss1>>s1;
      ss2>>s2;
      h_weight_info.push_back( static_cast<TH1F*>(f_weight->Get(TString("num_hs_tracks_mll_"+s1+"_"+s2+"_weight"))) );
    }
    if (h_weight_info.size() == 0) {
      ANA_MSG_ERROR("Error loading weight info from:" << input_weight_file);
      return StatusCode::FAILURE;
    }
    ANA_MSG_INFO("Loaded weight info from " << input_weight_file);
  }

  //print properties values
  ANA_MSG_INFO("Properties values:");
  ANA_MSG_INFO("lep1_min_pt = " << lep1_min_pt);
  ANA_MSG_INFO("lep2_min_pt = " << lep2_min_pt);
  ANA_MSG_INFO("lep_max_eta = " << lep_max_eta);
  ANA_MSG_INFO("dilep_min_pt = " << dilep_min_pt);
  ANA_MSG_INFO("tracks_min_pt = " << tracks_min_pt);
  ANA_MSG_INFO("tracks_max_eta = " << tracks_max_eta);
  ANA_MSG_INFO("filter_by_selections = " << filter_by_selections);
  ANA_MSG_INFO("input_trk_eff_file = " << input_trk_eff_file);
  ANA_MSG_INFO("input_pu_file = " << input_pu_file);

  return StatusCode::SUCCESS;
}

StatusCode TruthAnalysis :: execute ()
{
  //reset variables
  m_dilep_pt = -1.0;
  m_dilep_m = -1.0;
  m_WW_m = -1.0;
  m_WW_pt = -1.0;
  m_diphoton_m = -1.0;
  m_numPUtracks = 0;
  m_numHStracks = -1;
  m_numHStracks2 = -1;
  m_rand = -1.;
  m_mcWeight = -1;
  m_nTrackweight = -1;
  m_mcWeight_0trk = -1;
  m_mcWeight_1trk = -1;
  m_mcWeight_2trk = -1;
  m_mcWeight_3trk = -1;
  m_lep_pt->clear();
  m_lep_eta->clear();
  m_lep_phi->clear();
  m_lep_charge->clear();
  m_lep_pdgid->clear();
  m_photon_p->clear();
  m_W_pt->clear();
  m_W_m->clear();
  m_W_phi->clear();
  m_W_eta->clear();
  m_trk_pt->clear();
  m_trk_eta->clear();
  m_trk_phi->clear();
  m_trk_charge->clear();
  m_trk_pdgid->clear();

  m_truedilep_pt = -1.;
  m_truedilep_m = -1.;
  m_numTrueHSparticles = -1;
  m_truelep_pt->clear();
  m_truelep_eta->clear();
  m_truelep_phi->clear();
  m_truelep_charge->clear();
  m_truelep_pdgid->clear();
  m_truetrk_pt->clear();
  m_truetrk_eta->clear();
  m_truetrk_phi->clear();
  m_truetrk_charge->clear();
  m_truetrk_pdgid->clear();
  m_true_twoleps = false;
  m_true_twoOSleps = false;
  m_true_absproduct_of_lepton_pdgidvals = 0;
  m_true_leadlep_ptpass = false;
  m_true_subleadlep_ptpass = false;
  m_true_dilep_masspass = false;
  m_true_dilep_ptpass = false;
  m_true_numCh_and_PU = -1;

  for (int i=0; i<ncuts;i++)
    m_pass_sel->at(i)=false;


  //get event info
  const xAOD::EventInfo* eventInfo; 
  ANA_CHECK(evtStore()->retrieve(eventInfo,"EventInfo")); 
  double mcWeight = 1;
  const std::vector< float > weights = eventInfo->mcEventWeights();
  if( weights.size() > 0 ) mcWeight = weights[0];  

  hist("hCutFlow_Sum")->Fill(0.5, 1);
  hist("hCutFlow_Sum")->Fill(1.5, mcWeight);
  hist("hCutFlow_Sum")->Fill(2.5, mcWeight*mcWeight);

  m_mcWeight = mcWeight;
  
  //How many PU tracks in window?  Compare rand against the integral of the PDF of the track multiplicity
  if (h_pu_info != nullptr) {
    double rand = m_rnd->Rndm();
    m_rand = rand;
    double integral = 0;
    for(int b = 1; b < h_pu_info->GetNbinsX(); b++){
      integral += h_pu_info->GetBinContent(b);
      if(integral > rand){
	m_numPUtracks = b-1;
	//m_numHStracks = b-1;
	break;
      }
    }
  }
  
  
  // get truth particle container of interest
  const xAOD::TruthParticleContainer* truthParts = 0;
  ANA_CHECK (evtStore()->retrieve( truthParts, "TruthParticles"));
  ANA_MSG_DEBUG ("execute(): number of truth particles = " << truthParts->size());

  passCut(cut_nocut);

  // trigger, assumed 100% efficiency after other selections
  passCut(cut_trigger);


  int tmpsize = 0;
  int tmpsizetruth = 0;

  // loop over the particles in the container
  // and select fiducial leptons and tracks
  for (const xAOD::TruthParticle *part : *truthParts) {    
    //select fiducial truth particles
    ANA_MSG_VERBOSE ("Particle PDG " << part->auxdata<int>("pdgId") << ", pT=" << part->pt() << ", eta=" << part->eta() << ", charge=" << part->charge() << ", status=" << part->auxdata<int>("status") << ", barcode=" << part->auxdata<int>("barcode"));
    if( abs(part->auxdata<int>("pdgId")) == 24 || abs(part->auxdata<int>("pdgId")) == 11 || abs(part->auxdata<int>("pdgId")) == 13 || abs(part->auxdata<int>("pdgId")) == 15 ||  abs(part->auxdata<int>("pdgId")) == 23 || (part->auxdata<int>("status") < 25 && part->auxdata<int>("status") > 20) ){
      ANA_MSG_DEBUG ("Particle PDG " << part->auxdata<int>("pdgId") << ", pT=" << part->pt() << ", pZ=" << part->pz() << ", M=" << part->m() << ", eta=" << part->eta() << ", charge=" << part->charge() << ", status=" << part->auxdata<int>("status") << ", barcode=" << part->auxdata<int>("barcode"));
    }

    if(abs(part->auxdata<int>("pdgId")) == 22 && part->auxdata<int>("status") == 21 && ( part->auxdata<int>("barcode") == 3 || part->auxdata<int>("barcode") == 4 ))
    {
      	m_photon_p->push_back(part->pz());
    }
    if(abs(part->auxdata<int>("pdgId")) == 24 && part->auxdata<int>("status") == 22 && ( part->auxdata<int>("barcode") == 5 || part->auxdata<int>("barcode") == 6 ))
    {
      	m_W_pt->push_back(part->pt());
      	m_W_m->push_back(part->m());
      	m_W_eta->push_back(part->eta());
      	m_W_phi->push_back(part->phi());
    }
    if (part->pt() < std::min(tracks_min_pt,lep2_min_pt)) continue;
    if (part->abseta() > std::max(tracks_max_eta,lep_max_eta)) continue;
    if (part->charge() == 0) continue;
    if (part->auxdata<int>("status") != 1) continue;
    if (part->auxdata<int>("barcode") > 200000) continue;
    ANA_MSG_VERBOSE("Pass pre-fiducial cuts");

    int pdgid = part->auxdata<int>("pdgId");
    if ( (abs(pdgid) == 11) or (abs(pdgid)==13) ) {
      if ((part->pt() > lep2_min_pt) and
	  (part->abseta() < lep_max_eta) ) {      

	m_truelep_pt->push_back(part->pt());
	m_truelep_eta->push_back(part->eta());
	m_truelep_phi->push_back(part->phi());
	m_truelep_charge->push_back(static_cast<int>(part->charge()));
	m_truelep_pdgid->push_back(pdgid);

	if ( h_electron_eff != nullptr && (abs(pdgid) == 11) ) {
	  ANA_MSG_VERBOSE ("Candidate electron; is it reconstructed?");
	  int xbin = h_electron_eff->GetXaxis()->FindBin(part->pt());
	  int ybin = h_electron_eff->GetYaxis()->FindBin(part->eta());
	  float trk_eff = 0.0; //default if underflow
	  if ( xbin > 0 && ybin > 0 ) {
	    //use last bin if overflow
	    if (xbin > h_electron_eff->GetNbinsX()){
	      xbin = h_electron_eff->GetNbinsX();
	    }
	    trk_eff = h_electron_eff->GetBinContent(xbin, ybin);      
	    if ( std::abs(part->eta()) > 2.5 ){
	      trk_eff = 0.0;
	    }
	  }    
	  if (m_rnd->Rndm() > trk_eff){
	    continue;
	  }
	}

	if ( h_muon_eff != nullptr && (abs(pdgid) == 13) ) {
	  ANA_MSG_VERBOSE ("Candidate muon; is it reconstructed?");
	  int xbin = h_muon_eff->GetXaxis()->FindBin(part->pt());
	  int ybin = h_muon_eff->GetYaxis()->FindBin(part->eta());
	  float trk_eff = 0.0; //default if underflow
	  if ( xbin > 0 && ybin > 0 ) {
	    //use last bin if overflow
	    if (xbin > h_muon_eff->GetNbinsX()){
	      xbin = h_muon_eff->GetNbinsX();
	    }
	    trk_eff = h_muon_eff->GetBinContent(xbin, ybin);      
	    if ( std::abs(part->eta()) > 2.5 ){
	      trk_eff = 0.0;
	    }
	  }    
	  if (m_rnd->Rndm() > trk_eff)
	    continue;
	}


	ANA_MSG_VERBOSE ("It was recoed.  Store as lepton");
	m_lep_pt->push_back(part->pt());
	m_lep_eta->push_back(part->eta());
	m_lep_phi->push_back(part->phi());
	m_lep_charge->push_back(static_cast<int>(part->charge()));
	m_lep_pdgid->push_back(pdgid);
	continue; //do not consider it as track
      }
    }

    if (part->pt() < tracks_min_pt) continue;
    if (part->abseta() > tracks_max_eta) continue;
    ANA_MSG_VERBOSE("Candidate track");
    
    m_truetrk_pt->push_back(part->pt());
    m_truetrk_eta->push_back(part->eta());
    m_truetrk_phi->push_back(part->phi());
    m_truetrk_charge->push_back(static_cast<int>(part->charge()));
    m_truetrk_pdgid->push_back(pdgid);

    tmpsizetruth += 1;

    //apply parametrized tracking efficiency, if requested
    if (h_trk_eff_pt_eta != nullptr) {
      //int ibin = h_trk_eff_pt_eta->FindBin(part->pt(), part->eta());
      int xbin = h_trk_eff_pt_eta->GetXaxis()->FindBin(part->pt());
      int ybin = h_trk_eff_pt_eta->GetYaxis()->FindBin(part->eta());
      float trk_eff = 0.0; //default if underflow
      if ( xbin > 0 && ybin > 0 ) {
	//use last bin if overflow
	if (xbin > h_trk_eff_pt_eta->GetNbinsX()){
	  xbin = h_trk_eff_pt_eta->GetNbinsX();
	}
	trk_eff = h_trk_eff_pt_eta->GetBinContent(xbin, ybin);      
	if ( std::abs(part->eta()) > 2.5 ){
	  trk_eff = 0.0;
	}
      }    
      if (m_rnd->Rndm() > trk_eff)
	continue;
    }

    tmpsize += 1;

    //store basic kinematic and particle info
    ANA_MSG_VERBOSE("Store track");
    m_trk_pt->push_back(part->pt());
    m_trk_eta->push_back(part->eta());
    m_trk_phi->push_back(part->phi());
    m_trk_charge->push_back(static_cast<int>(part->charge()));
    m_trk_pdgid->push_back(pdgid);
    
  } // end for loop over truth particles



  ////////    TRUTH ANALYSIS STUFF with no pseudo-reco  /////////////////

  m_numTrueHSparticles = m_trk_pt->size();
  m_true_numCh_and_PU = m_numTrueHSparticles + m_numPUtracks;
  if (m_truelep_pt->size() == 2) {
    m_true_twoleps = true;
    if(m_truelep_charge->at(0)*m_truelep_charge->at(1) == -1){
      m_true_twoOSleps = true;
      m_true_absproduct_of_lepton_pdgidvals = abs(m_truelep_pdgid->at(0)*m_truelep_pdgid->at(1));

      if (std::min(m_truelep_pt->at(0), m_truelep_pt->at(1)) > lep2_min_pt) {
	m_true_subleadlep_ptpass = true;
      }
      if (std::max(m_truelep_pt->at(0), m_truelep_pt->at(1)) > lep1_min_pt) {
	m_true_leadlep_ptpass = true;
      }

      if( m_true_leadlep_ptpass && m_true_subleadlep_ptpass ){
	TLorentzVector truelep1, truelep2;
	int index_truelep1=0;
	int index_truelep2=1;
	if (m_truelep_pt->at(1) > m_truelep_pt->at(0)) {
	  index_truelep1=1;
	  index_truelep2=0;
	}
	truelep1.SetPtEtaPhiM(m_truelep_pt->at(index_truelep1), m_truelep_eta->at(index_truelep1), m_truelep_phi->at(index_truelep1), 0.0);
	truelep2.SetPtEtaPhiM(m_truelep_pt->at(index_truelep2), m_truelep_eta->at(index_truelep2), m_truelep_phi->at(index_truelep2), 0.0);
	m_truedilep_m = (truelep1+truelep2).M();
	m_truedilep_pt = (truelep1+truelep2).Pt();

	if(m_truedilep_pt < dilep_min_pt){
	  m_true_dilep_ptpass = true;
	}
	if (m_truedilep_m < dilep_min_mass) {
	  m_true_dilep_masspass = true;
	}

      }
    }
  }
  ///////////////////////////////////////////




  int num_track_in_window = m_trk_pt->size() + m_numPUtracks;
  if(num_track_in_window == 0) m_mcWeight_0trk = mcWeight;
  if(num_track_in_window == 1) m_mcWeight_1trk = mcWeight;
  if(num_track_in_window == 2) m_mcWeight_2trk = mcWeight;
  if(num_track_in_window == 3) m_mcWeight_3trk = mcWeight;

  if(m_photon_p->size() == 2){
    float p1 = abs(m_photon_p->at(0));
    float p2 = abs(m_photon_p->at(1));
    m_diphoton_m = sqrt((p1+p2)*(p1+p2) - (p1-p2)*(p1-p2));
    ANA_MSG_DEBUG("photons m"<<m_diphoton_m);
  }
  std::vector<int> mass_bins = {44, 60, 90, 110, 140, 180, 190, 200, 210, 220, 230, 250, 270, 290, 320, 350, 380, 410, 450, 500, 600};
  float additionalweight = 1;
  if(do_lpair_reweight){
    for(unsigned int i = 0; i < mass_bins.size() - 1; i++){
      if( ( (m_diphoton_m/1e3) < mass_bins.at(i+1) ) && ( (m_diphoton_m/1e3) > mass_bins.at(i) ) ){
	additionalweight = additionalweight*h_weight_info.at(i)->GetBinContent( num_track_in_window - m_numPUtracks + 1);
      }      
    }
  }
  m_nTrackweight = additionalweight;


  if(m_trk_pt->size() > 0){
    m_numHStracks = m_trk_pt->size();
  }
  else{
    m_numHStracks = 0;
  }

  hist("num_hs_tracks")->Fill(tmpsize);
  hist("num_hs_truth_tracks")->Fill(tmpsizetruth);

  // Now evaluate event-level selections  
  ANA_MSG_VERBOSE("Checking event selections");
  ANA_MSG_DEBUG("NLep = " << m_lep_pt->size() << ", NTracks = " << m_trk_pt->size());
  hist("num_fiducial_leptons")->Fill(m_lep_pt->size());

  //int tmpsize = m_trk_pt->size();
  //m_numHStracks = tmpsize;

  if (m_lep_pt->size() != 2) {saveTree(); return StatusCode::SUCCESS;}
  if (m_lep_charge->at(0)*m_lep_charge->at(1) != -1) {saveTree(); return StatusCode::SUCCESS;}
  if ( abs(m_lep_pdgid->at(0)*m_lep_pdgid->at(1)) != 11*13 && abs(m_lep_pdgid->at(0)*m_lep_pdgid->at(1)) != 13*13 && abs(m_lep_pdgid->at(0)*m_lep_pdgid->at(1)) != 11*11 ) {saveTree(); return StatusCode::SUCCESS;}
  passCut(cut_lep_ocof);
  ANA_MSG_VERBOSE("Pass cut_lep_ocof");

  if (std::min(m_lep_pt->at(0), m_lep_pt->at(1)) < lep2_min_pt) {saveTree(); return StatusCode::SUCCESS;}
  if (std::max(m_lep_pt->at(0), m_lep_pt->at(1)) < lep1_min_pt) {saveTree(); return StatusCode::SUCCESS;}
  passCut(cut_lep_minpt);
  ANA_MSG_VERBOSE("Pass cut_lep_minpt");

  TLorentzVector lep1, lep2;
  int index_lep1=0;
  int index_lep2=1;
  if (m_lep_pt->at(1) > m_lep_pt->at(0)) {
    index_lep1=1;
    index_lep2=0;
  }
  lep1.SetPtEtaPhiM(m_lep_pt->at(index_lep1), m_lep_eta->at(index_lep1), m_lep_phi->at(index_lep1), 0.0);
  lep2.SetPtEtaPhiM(m_lep_pt->at(index_lep2), m_lep_eta->at(index_lep2), m_lep_phi->at(index_lep2), 0.0);
  m_dilep_m = (lep1+lep2).M();
  m_dilep_pt = (lep1+lep2).Pt();
  ANA_MSG_DEBUG("ll m"<<m_dilep_m);
  ANA_MSG_DEBUG("ll pt"<<m_dilep_pt);
  /*
  if(m_photon_p->size() == 2){
    float p1 = abs(m_photon_p->at(0));
    float p2 = abs(m_photon_p->at(1));
    m_diphoton_m = sqrt((p1+p2)*(p1+p2) - (p1-p2)*(p1-p2));
    ANA_MSG_DEBUG("photons m"<<m_diphoton_m);
  }
  */
  if(m_W_pt->size() == 2){
    TLorentzVector W1, W2;
    W1.SetPtEtaPhiM(m_W_pt->at(0), m_W_eta->at(0), m_W_phi->at(0), m_W_m->at(0));
    W2.SetPtEtaPhiM(m_W_pt->at(1), m_W_eta->at(1), m_W_phi->at(1), m_W_m->at(1));
    m_WW_m = (W1+W2).M();
    m_WW_pt = (W1+W2).Pt();
    ANA_MSG_DEBUG("WW m"<<m_WW_m);
    ANA_MSG_DEBUG("WW pt"<<m_WW_pt);
  }

  std::vector<int> ptll_bins = {0,5,10,15,20,25,30,50,70,90,110,130};
  //std::vector<int> mass_bins = {44, 60, 90, 110, 140, 180, 190, 200, 210, 220, 230, 250, 270, 290, 320, 350, 380, 410, 450, 500, 600};
  for(unsigned int i = 0; i < ptll_bins.size() - 1; i++){
    if( ( (m_dilep_pt/1e3) < ptll_bins.at(i+1) ) && ( (m_dilep_pt/1e3) > ptll_bins.at(i) ) ){
      std::stringstream ss1;
      ss1<<ptll_bins.at(i);
      std::stringstream ss2;
      ss2<<ptll_bins.at(i+1);
      std::string s1;
      std::string s2;
      ss1>>s1;
      ss2>>s2;
      hist("num_hs_tracks_ptll_"+s1+"_"+s2)->Fill(tmpsize);
    }
  }
  for(unsigned int i = 0; i < mass_bins.size() - 1; i++){
    if( ( (m_dilep_m/1e3) < mass_bins.at(i+1) ) && ( (m_dilep_m/1e3) > mass_bins.at(i) ) ){
      std::stringstream ss1;
      ss1<<mass_bins.at(i);
      std::stringstream ss2;
      ss2<<mass_bins.at(i+1);
      std::string s1;
      std::string s2;
      ss1>>s1;
      ss2>>s2;
      hist("num_hs_tracks_mll_"+s1+"_"+s2)->Fill(tmpsize);
    }    
  }
  

  hist("dilep_m")->Fill(m_dilep_m/GeV);

  if(abs(m_lep_pdgid->at(0)*m_lep_pdgid->at(1)) == 11*13){
    FillHistos ("nominal/emuOS/EWW/Pre_Mll", mcWeight,  lep1, lep2, num_track_in_window, additionalweight);
  }
  else if(abs(m_lep_pdgid->at(0)*m_lep_pdgid->at(1)) == 13*13){
    FillHistos ("nominal/mumuOS/EWW/Pre_Mll", mcWeight,  lep1, lep2, num_track_in_window, additionalweight);
  }
  else if(abs(m_lep_pdgid->at(0)*m_lep_pdgid->at(1)) == 11*11){
    FillHistos ("nominal/eeOS/EWW/Pre_Mll", mcWeight,  lep1, lep2, num_track_in_window, additionalweight);
  }

  if (m_dilep_m < dilep_min_mass) {saveTree(); return StatusCode::SUCCESS;}
  passCut(cut_m_ll);
  ANA_MSG_VERBOSE("Pass cut_m_ll");

  
  if( ( abs(m_lep_pdgid->at(0)*m_lep_pdgid->at(1)) == 11*11 || abs(m_lep_pdgid->at(0)*m_lep_pdgid->at(1)) == 13*13 ) && m_dilep_m > 160000 && num_track_in_window == 0 ){
    FillHistos_kristin ("CR4_exclSF", mcWeight,  lep1, lep2, num_track_in_window, additionalweight);
  }
  if( ( abs(m_lep_pdgid->at(0)*m_lep_pdgid->at(1)) == 11*11 || abs(m_lep_pdgid->at(0)*m_lep_pdgid->at(1)) == 13*13 ) && m_dilep_m > 160000 && num_track_in_window == 1 ){
    FillHistos_kristin ("CR4_exclSF_1track", mcWeight,  lep1, lep2, num_track_in_window, additionalweight);
  }

  if( abs(m_lep_pdgid->at(0)*m_lep_pdgid->at(1)) == 11*13 && m_dilep_pt < dilep_min_pt && num_track_in_window == 0 ){
    FillHistos_kristin ("CR1", mcWeight,  lep1, lep2, num_track_in_window, additionalweight);
  }
  if( abs(m_lep_pdgid->at(0)*m_lep_pdgid->at(1)) == 11*13 && m_dilep_pt < dilep_min_pt && num_track_in_window == 1 ){
    FillHistos_kristin ("CR1_1track", mcWeight,  lep1, lep2, num_track_in_window, additionalweight);
    FillHistos_kristin ("CR3_1track", mcWeight,  lep1, lep2, num_track_in_window, additionalweight);
  }
  if( abs(m_lep_pdgid->at(0)*m_lep_pdgid->at(1)) == 11*13 && m_dilep_pt < dilep_min_pt && ( num_track_in_window > 0 && num_track_in_window < 5 ) ){
    FillHistos_kristin ("CR3", mcWeight,  lep1, lep2, num_track_in_window, additionalweight);
    if( num_track_in_window == 1 || num_track_in_window == 2 ){
      FillHistos_kristin ("CR3_1_to_2", mcWeight,  lep1, lep2, num_track_in_window, additionalweight);
    }
    if( num_track_in_window == 3 || num_track_in_window == 4 ){
      FillHistos_kristin ("CR3_3_to_4", mcWeight,  lep1, lep2, num_track_in_window, additionalweight);
    }
  }

  hist("dilep_pt")->Fill(m_dilep_pt/GeV);

  if(abs(m_lep_pdgid->at(0)*m_lep_pdgid->at(1)) == 11*13){
    FillHistos ("nominal/emuOS/EWW/Pre_Ptll", mcWeight,  lep1, lep2, num_track_in_window, additionalweight);
  }
  else if(abs(m_lep_pdgid->at(0)*m_lep_pdgid->at(1)) == 13*13){
    FillHistos ("nominal/mumuOS/EWW/Pre_Ptll", mcWeight,  lep1, lep2, num_track_in_window, additionalweight);
  }
  else if(abs(m_lep_pdgid->at(0)*m_lep_pdgid->at(1)) == 11*11){
    FillHistos ("nominal/eeOS/EWW/Pre_Ptll", mcWeight,  lep1, lep2, num_track_in_window, additionalweight);
  }

  if (m_dilep_pt < dilep_min_pt) {
    if( num_track_in_window > 0){
      if(abs(m_lep_pdgid->at(0)*m_lep_pdgid->at(1)) == 11*13){
	FillHistos ("nominal/emuOS/EWW/CR/Ztt_all_tracks", mcWeight,  lep1, lep2, num_track_in_window, additionalweight);
	if(num_track_in_window < 5){
	  FillHistos ("nominal/emuOS/EWW/CR/Ztt_1_to_4", mcWeight,  lep1, lep2, num_track_in_window, additionalweight);
	}
      }
    }
    saveTree(); return StatusCode::SUCCESS;
  }
  passCut(cut_pt_ll);
  ANA_MSG_VERBOSE("Pass cut_pt_ll");

  m_numHStracks2 = m_trk_pt->size();


  hist("num_hs_tracksRW")->Fill(tmpsize, additionalweight);

  if(abs(m_lep_pdgid->at(0)*m_lep_pdgid->at(1)) == 11*13){
    FillHistos ("nominal/emuOS/EWW/Ptll", mcWeight,  lep1, lep2, num_track_in_window, additionalweight);
  }
  else if(abs(m_lep_pdgid->at(0)*m_lep_pdgid->at(1)) == 13*13){
    FillHistos ("nominal/mumuOS/EWW/Ptll", mcWeight,  lep1, lep2, num_track_in_window, additionalweight);
  }
  else if(abs(m_lep_pdgid->at(0)*m_lep_pdgid->at(1)) == 11*11){
    FillHistos ("nominal/eeOS/EWW/Ptll", mcWeight,  lep1, lep2, num_track_in_window, additionalweight);
  }

  if( num_track_in_window > 0){
    if(abs(m_lep_pdgid->at(0)*m_lep_pdgid->at(1)) == 11*13){
      FillHistos ("nominal/emuOS/EWW/CR/InclWW_all_tracks", mcWeight,  lep1, lep2, num_track_in_window, additionalweight);
    }
    else if(abs(m_lep_pdgid->at(0)*m_lep_pdgid->at(1)) == 13*13){
      FillHistos ("nominal/mumuOS/EWW/CR/InclWW_all_tracks", mcWeight,  lep1, lep2, num_track_in_window, additionalweight);
    }
    else if(abs(m_lep_pdgid->at(0)*m_lep_pdgid->at(1)) == 11*11){
      FillHistos ("nominal/eeOS/EWW/CR/InclWW_all_tracks", mcWeight,  lep1, lep2, num_track_in_window, additionalweight);
    }

    if( std::min(m_lep_pt->at(0), m_lep_pt->at(1)) > 27*GeV && m_dilep_m > 55000 ){
      if(abs(m_lep_pdgid->at(0)*m_lep_pdgid->at(1)) == 11*13){
	FillHistos ("nominal/emuOS/EWW/CR/InclWW_alltrk_WWmeas", mcWeight,  lep1, lep2, num_track_in_window, additionalweight);
      }
      else if(abs(m_lep_pdgid->at(0)*m_lep_pdgid->at(1)) == 13*13){
	FillHistos ("nominal/mumuOS/EWW/CR/InclWW_alltrk_WWmeas", mcWeight,  lep1, lep2, num_track_in_window, additionalweight);
      }
      else if(abs(m_lep_pdgid->at(0)*m_lep_pdgid->at(1)) == 11*11){
	FillHistos ("nominal/eeOS/EWW/CR/InclWW_alltrk_WWmeas", mcWeight,  lep1, lep2, num_track_in_window, additionalweight);
      }
    }

  }

  if( num_track_in_window > 0 && num_track_in_window < 5 ){
    if(abs(m_lep_pdgid->at(0)*m_lep_pdgid->at(1)) == 11*13){
      FillHistos ("nominal/emuOS/EWW/CR/InclWW", mcWeight,  lep1, lep2, num_track_in_window, additionalweight);
      FillHistos_kristin ("CR2", mcWeight,  lep1, lep2, num_track_in_window, additionalweight);
      if( num_track_in_window == 1 ){
	FillHistos_kristin ("CR2_1track", mcWeight,  lep1, lep2, num_track_in_window, additionalweight);
	FillHistos_kristin ("SR_1track", mcWeight,  lep1, lep2, num_track_in_window, additionalweight);
      }
      if( num_track_in_window == 1 || num_track_in_window == 2 ){
	FillHistos_kristin ("CR2_1_to_2", mcWeight,  lep1, lep2, num_track_in_window, additionalweight);
      }
      if( num_track_in_window == 3 || num_track_in_window == 4 ){
	FillHistos_kristin ("CR2_3_to_4", mcWeight,  lep1, lep2, num_track_in_window, additionalweight);
      }
    }
    else if(abs(m_lep_pdgid->at(0)*m_lep_pdgid->at(1)) == 13*13){
      FillHistos ("nominal/mumuOS/EWW/CR/InclWW", mcWeight,  lep1, lep2, num_track_in_window, additionalweight);
    }
    else if(abs(m_lep_pdgid->at(0)*m_lep_pdgid->at(1)) == 11*11){
      FillHistos ("nominal/eeOS/EWW/CR/InclWW", mcWeight,  lep1, lep2, num_track_in_window, additionalweight);
    }
  }
  if( num_track_in_window == 0){
    if(abs(m_lep_pdgid->at(0)*m_lep_pdgid->at(1)) == 11*13){
      FillHistos ("nominal/emuOS/EWW/Excl", mcWeight,  lep1, lep2, num_track_in_window, additionalweight);
      FillHistos_kristin ("SR", mcWeight,  lep1, lep2, num_track_in_window, additionalweight);
    }
    else if(abs(m_lep_pdgid->at(0)*m_lep_pdgid->at(1)) == 13*13){
      FillHistos ("nominal/mumuOS/EWW/Excl", mcWeight,  lep1, lep2, num_track_in_window, additionalweight);
    }
    else if(abs(m_lep_pdgid->at(0)*m_lep_pdgid->at(1)) == 11*11){
      FillHistos ("nominal/eeOS/EWW/Excl", mcWeight,  lep1, lep2, num_track_in_window, additionalweight);
    }
  }

  hist("num_fiducial_tracks")->Fill(m_trk_pt->size());
  //if (m_trk_pt->size() > tracks_max_n) {
  //passCut(cut_exclusive);
  //}
  saveTree();
  //if (m_trk_pt->size() > tracks_max_n) {saveTree(); return StatusCode::SUCCESS;}
  if (m_trk_pt->size() > tracks_max_n) {return StatusCode::SUCCESS;}
  passCut(cut_exclusive);
  ANA_MSG_VERBOSE("Pass cut_exclusive");

  //Plots after all selections
  hist("sr_dilep_pt")->Fill(m_dilep_pt/GeV);

  // Fill the event into the tree
  saveTree();
  return StatusCode::SUCCESS;
}


StatusCode TruthAnalysis :: finalize ()
{
  // This method is the mirror image of initialize(), meaning it gets
  // called after the last event has been processed on the worker node
  // and allows you to finish up any objects you created in
  // initialize() before they are written to disk.  This is actually
  // fairly rare, since this happens separately for each worker node.
  // Most of the time you want to do your post-processing on the
  // submission node after all your histogram outputs have been
  // merged.
  return StatusCode::SUCCESS;
}


TruthAnalysis :: ~TruthAnalysis () {
  if (m_pass_sel) {
    delete m_lep_pt;
    delete m_lep_eta;
    delete m_lep_phi;
    delete m_lep_charge;
    delete m_lep_pdgid;
    
    delete m_trk_pt;
    delete m_trk_eta;
    delete m_trk_phi;
    delete m_trk_charge;
    delete m_trk_pdgid;
    
    delete m_pass_sel;
  }
}

void TruthAnalysis :: saveTree ()
{
  //save if requested or pass all selections
  if ((not filter_by_selections) or (m_pass_sel->at(ncuts-1)))
    tree ("analysis")->Fill ();
}

void TruthAnalysis :: passCut (int cut)
{
  hist("cutflow")->Fill(cut);
  m_pass_sel->at(cut)=true;
}

 StatusCode TruthAnalysis :: MakeHistos ()
{

  std::vector<TString> regions = {"nominal/emuOS/EWW/CR/InclWW/","nominal/emuOS/EWW/CR/InclWW_all_tracks/","nominal/emuOS/EWW/Excl/","nominal/emuOS/EWW/Ptll/","nominal/mumuOS/EWW/CR/InclWW/","nominal/mumuOS/EWW/CR/InclWW_all_tracks/","nominal/mumuOS/EWW/Excl/","nominal/mumuOS/EWW/Ptll/","nominal/eeOS/EWW/CR/InclWW/","nominal/eeOS/EWW/CR/InclWW_all_tracks/","nominal/eeOS/EWW/Excl/","nominal/eeOS/EWW/Ptll/","nominal/emuOS/EWW/CR/InclWW_alltrk_WWmeas/","nominal/mumuOS/EWW/CR/InclWW_alltrk_WWmeas/","nominal/eeOS/EWW/CR/InclWW_alltrk_WWmeas/", "nominal/emuOS/EWW/CR/Ztt_all_tracks/", "nominal/emuOS/EWW/CR/Ztt_1_to_4/", "nominal/emuOS/EWW/Pre_Ptll/", "nominal/emuOS/EWW/Pre_Mll/", "nominal/eeOS/EWW/Pre_Ptll/", "nominal/eeOS/EWW/Pre_Mll/", "nominal/mumuOS/EWW/Pre_Ptll/", "nominal/mumuOS/EWW/Pre_Mll/", "nominal/emuOS/EWW/CR/InclWW_all_tracks_near100/", "nominal/emuOS/EWW/CR/InclWW_near100/", "nominal/emuOS/EWW/CR/ttbar/", "nominal/emuOS/EWW/CR/ttbar_1j1b/", "nominal/emuOS/EWW/CR/ttbar_2j/", "nominal/emuOS/EWW/CR/ttbar_2Fj/", "nominal/emuOS/EWW/CR/ttbar_1j1b_pvcut/", "nominal/emuOS/EWW/CR/ttbar_2j_pvcut/", "nominal/emuOS/EWW/CR/ttbar_2Fj_pvcut/", "nominal/emuOS/EWW/CR/ttbar_1j1b_opp_pvcut/", "nominal/emuOS/EWW/CR/ttbar_2j_opp_pvcut/", "nominal/emuOS/EWW/CR/ttbar_2Fj_opp_pvcut/"};

  for(unsigned int i = 0; i < regions.size(); i++){
    //hNtrkPVW
    ANA_CHECK (book ( TH1D ( regions.at(i)+"NX/hNtrkPVW", "PV Track multiplicity weighted; N_{trk}^{PV,W}; Entries", 150, -0.5, 149.5) ));
    //hNtrkPVW_central
    ANA_CHECK (book ( TH1D ( regions.at(i)+"NX/hNtrkPVW_central", "Central PV Track multiplicity weighted; N_{centraltrk}^{PV,W}; Entries", 150, -0.5, 149.5) ));
    //hNtrkPVW_central
    ANA_CHECK (book ( TH1D ( regions.at(i)+"NX/hNtrkPVW_forward", "Forward PV Track multiplicity weighted; N_{forwardtrk}^{PV,W}; Entries", 150, -0.5, 149.5) ));
    //hNtrkPV
    ANA_CHECK (book ( TH1D ( regions.at(i)+"NX/hNtrkPV", "PV Track multiplicity; N_{trk}^{PV}; Entries", 150, -0.5, 149.5) ));
    //hDiLeptonMass
    ANA_CHECK (book ( TH1D ( regions.at(i)+"ll/hDiLeptonMass", "DiLepton mass;M_{ll}[GeV]; Entries", 2000, 0, 1000 ) ));
    //hDiLeptonPt
    ANA_CHECK (book ( TH1D ( regions.at(i)+"ll/hDiLeptonPt", "DiLepton p_{T}; p_{T,ll}; Entries", 200, 0, 200 ) ));
    //hDiLeptonAco
    ANA_CHECK (book ( TH1D ( regions.at(i)+"ll/hDiLeptonAco", "DiLepton Acoplanarity; Aco; Entries", 100, 0, 1 ) ));
    ANA_CHECK (book ( TH1D ( regions.at(i)+"ll/LogAcoplan", "Log DiLepton Acoplanarity; log Aco; Entries", 60, -5., 0. ) ));
    //hPtLead
    ANA_CHECK (book ( TH1D ( regions.at(i)+"Lepton/hPtLead", "p_{T} leading; p_{T}^{1}; Entries", 200, 0, 200 ) ));
    //hPtSubLead
    ANA_CHECK (book ( TH1D ( regions.at(i)+"Lepton/hPtSubLead", "p_{T} Subleading; p_{T}^{1}; Entries", 200, 0, 200 ) ));
    //hEtaLead
    ANA_CHECK (book ( TH1D ( regions.at(i)+"Lepton/hEtaLead", "#eta leading; #eta^{1}; Entries", 60, -3, 3 ) ));
    //hEtaSubLead
    ANA_CHECK (book ( TH1D ( regions.at(i)+"Lepton/hEtaSubLead", "#eta Subleading; #eta^{1}; Entries", 60, -3, 3 ) ));
    //hPhiLead
    ANA_CHECK (book ( TH1D ( regions.at(i)+"Lepton/hPhiLead", "#phi leading; #phi^{1}; Entries", 60, -TMath::Pi(), +TMath::Pi() ) ));
    //hPhiSubLead
    ANA_CHECK (book ( TH1D ( regions.at(i)+"Lepton/hPhiSubLead", "#phi Subleading; #phi^{1}; Entries", 60, -TMath::Pi(), +TMath::Pi() ) ));
    //hPtEtaLead
    ANA_CHECK (book ( TH2D ( regions.at(i)+"Lepton/hPtEtaLead", "p_{T}#eta leading; p_{T}^{1}; #eta^{1}; Entries", 200, 0, 200, 60, -3, 3 ) ));
    //hPtEtaSubLead
    ANA_CHECK (book ( TH2D ( regions.at(i)+"Lepton/hPtEtaSubLead", "p_{T}#eta Subleading; p_{T}^{1}; #eta^{1}; Entries", 200, 0, 200, 60, -3, 3 ) ));
    
    ANA_CHECK (book ( TH1D ( regions.at(i)+"Lepton/hLeadtopoetcone20_pT", "topoetcone20_pT^{1}; Entries", 100, -1, 1 ) ));
    ANA_CHECK (book ( TH1D ( regions.at(i)+"Lepton/hSubLeadtopoetcone20_pT", "topoetcone20_pT^{s1}; Entries", 100, -1, 1 ) ));
    ANA_CHECK (book ( TH1D ( regions.at(i)+"Lepton/hLeadMutopoetcone20_pT", "topoetcone20_pT^{1}; Entries", 100, -1, 1 ) ));
    ANA_CHECK (book ( TH1D ( regions.at(i)+"Lepton/hSubLeadMutopoetcone20_pT", "topoetcone20_pT^{s1}; Entries", 100, -1, 1 ) ));
    ANA_CHECK (book ( TH1D ( regions.at(i)+"Lepton/hLeadEltopoetcone20_pT", "topoetcone20_pT^{1}; Entries", 100, -1, 1 ) ));
    ANA_CHECK (book ( TH1D ( regions.at(i)+"Lepton/hSubLeadEltopoetcone20_pT", "topoetcone20_pT^{s1}; Entries", 100, -1, 1 ) ));
    ANA_CHECK (book ( TH1D ( regions.at(i)+"Lepton/hLeadptvarcone20_TightTTVA_pt1000_pT", "ptvarcone20_TightTTVA_pt1000_pT^{1}; Entries", 100, -1, 1 ) ));
    ANA_CHECK (book ( TH1D ( regions.at(i)+"Lepton/hSubLeadptvarcone20_TightTTVA_pt1000_pT", "ptvarcone20_TightTTVA_pt1000_pT^{s1}; Entries", 100, -1, 1 ) ));
    ANA_CHECK (book ( TH1D ( regions.at(i)+"Lepton/hLeadptvarcone30_TightTTVA_pt1000_pT", "ptvarcone30_TightTTVA_pt1000_pT^{1}; Entries", 100, -1, 1 ) ));
    ANA_CHECK (book ( TH1D ( regions.at(i)+"Lepton/hSubLeadptvarcone30_TightTTVA_pt1000_pT", "ptvarcone30_TightTTVA_pt1000_pT^{s1}; Entries", 100, -1, 1 ) ));
    ANA_CHECK (book ( TH1D ( regions.at(i)+"Lepton/hLeadptcone20_TightTTVA_pt1000_pT", "ptcone20_TightTTVA_pt1000_pT^{1}; Entries", 100, -1, 1 ) ));
    ANA_CHECK (book ( TH1D ( regions.at(i)+"Lepton/hSubLeadptcone20_TightTTVA_pt1000_pT", "ptcone20_TightTTVA_pt1000_pT^{s1}; Entries", 100, -1, 1 ) ));

    ANA_CHECK (book ( TH1D ( regions.at(i)+"ll/hDpvZ_lepZ", "| #Delta#z_{ll,PV}; | #Delta#z_{ll,PV} | [mm]; Entries", 100, 0, 5 ) ));
    ANA_CHECK (book ( TH1D ( regions.at(i)+"Jets/hNumJets", "Number of Jets; N_{Jets}; Entries", 20, 0, 20 ) ));
    ANA_CHECK (book ( TH1D ( regions.at(i)+"Jets/hNumBJets", "Number of B Jets; N_{BJets}; Entries", 10, 0, 10 ) ));
    ANA_CHECK (book ( TH1D ( regions.at(i)+"Jets/ hDphiJJ", "#Delta #Phi between leading jets; #Delta#phi_{JJ}; Entries", 10, 0, TMath::Pi() ) ));
    ANA_CHECK (book ( TH1D ( regions.at(i)+"Jets/ hDetaJJ", "| #Delta #Eta | between leading jets; #Delta#eta_{JJ}; Entries", 40, -10, 10 ) ));
    ANA_CHECK (book ( TH1D ( regions.at(i)+"Jets/ hDptJJ", "| #Delta Pt | between leading jets; #Delta# pt_{JJ} [GeV]; Entries", 10, 0, 100 ) ));
    ANA_CHECK (book ( TH1D ( regions.at(i)+"Jets/ hLeadJet_eta", "#Eta of leading jet; #eta_{lead jet}; Entries", 40, -10, 10 ) ));

  }

  std::vector<TString> regions_kristin = {"nominal/CR1/", "nominal/CR2/", "nominal/CR3/", "nominal/CR4_exclSF/", "nominal/SR/", "nominal/CR1_1track/", "nominal/CR2_1track/", "nominal/CR3_1track/", "nominal/CR4_exclSF_1track/", "nominal/SR_1track/", "nominal/CR2_1_to_2/", "nominal/CR3_1_to_2/", "nominal/CR2_3_to_4/", "nominal/CR3_3_to_4/"};

  for(unsigned int i = 0; i < regions_kristin.size(); i++){
    ANA_CHECK (book ( TH1D ( regions_kristin.at(i)+"nTrack_weighted", "PV Track multiplicity weighted; N_{trk}^{PV,W}; Entries", 150, -0.5, 149.5) ));
    ANA_CHECK (book ( TH1D ( regions_kristin.at(i)+"nTrack_noweight", "PV Track multiplicity; N_{trk}^{PV}; Entries", 150, -0.5, 149.5) ));
    ANA_CHECK (book ( TH1D ( regions_kristin.at(i)+"Mll", "DiLepton mass;M_{ll}[GeV]; Entries", 80, 0., 400. ) ));
    ANA_CHECK (book ( TH1D ( regions_kristin.at(i)+"Ptll", "DiLepton p_{T}; p_{T,ll}; Entries", 80, 0., 400. ) ));
    ANA_CHECK (book ( TH1D ( regions_kristin.at(i)+"LogAcoplan", "log10 of DiLepton Acoplanarity; log_{10}Aco; Entries", 60, -5., 0. ) ));
  }

  return StatusCode::SUCCESS;

}

void TruthAnalysis :: FillHistos (std::string region, double mcWeight,  TLorentzVector lep1, TLorentzVector lep2, int num_track_in_window, float additionalweight, int num_central, int num_forward )
{
  
  float dilep_pt;
  float dilep_m;
  dilep_m = (lep1+lep2).M();
  dilep_pt = (lep1+lep2).Pt();

  hist(region+"/NX/hNtrkPVW")->Fill(num_track_in_window, mcWeight*additionalweight);
  hist(region+"/NX/hNtrkPVW_central")->Fill(num_central, mcWeight*additionalweight);
  hist(region+"/NX/hNtrkPVW_forward")->Fill(num_forward, mcWeight*additionalweight);
  hist(region+"/NX/hNtrkPV")->Fill(num_track_in_window, mcWeight);
  hist(region+"/ll/hDiLeptonMass")->Fill(dilep_m/GeV, mcWeight*additionalweight);
  hist(region+"/ll/hDiLeptonPt")->Fill(dilep_pt/GeV, mcWeight*additionalweight);
  hist(region+"/ll/hDiLeptonAco")->Fill(1 - abs( lep1.Phi()-lep2.Phi() )/TMath::Pi(), mcWeight*additionalweight);
  hist(region+"/ll/LogAcoplan")->Fill( log10( 1 - abs( lep1.Phi()-lep2.Phi() )/TMath::Pi() ), mcWeight*additionalweight);
  hist(region+"/Lepton/hPtLead")->Fill(lep1.Pt()/GeV, mcWeight*additionalweight);
  hist(region+"/Lepton/hPtSubLead")->Fill(lep2.Pt()/GeV, mcWeight*additionalweight);
  hist(region+"/Lepton/hEtaLead")->Fill(lep1.Eta(), mcWeight*additionalweight);
  hist(region+"/Lepton/hEtaSubLead")->Fill(lep2.Eta(), mcWeight*additionalweight);
  hist(region+"/Lepton/hPhiLead")->Fill(lep1.Phi(), mcWeight*additionalweight);
  hist(region+"/Lepton/hPhiSubLead")->Fill(lep2.Phi(), mcWeight*additionalweight);
  ((TH2*) hist(region+"/Lepton/hPtEtaLead"))->Fill(lep1.Pt()/GeV, lep1.Eta(), mcWeight*additionalweight);
  ((TH2*) hist(region+"/Lepton/hPtEtaSubLead"))->Fill(lep2.Pt()/GeV, lep2.Eta(), mcWeight*additionalweight);
}

void TruthAnalysis :: FillHistos_kristin (std::string region, double mcWeight,  TLorentzVector lep1, TLorentzVector lep2, int num_track_in_window, float additionalweight)
{
  
  float dilep_pt;
  float dilep_m;
  dilep_m = (lep1+lep2).M();
  dilep_pt = (lep1+lep2).Pt();

  hist("nominal/"+region+"/nTrack_weighted")->Fill(num_track_in_window, mcWeight*additionalweight);
  hist("nominal/"+region+"/nTrack_noweight")->Fill(num_track_in_window, mcWeight);
  hist("nominal/"+region+"/Mll")->Fill(dilep_m/GeV, mcWeight*additionalweight);
  hist("nominal/"+region+"/Ptll")->Fill(dilep_pt/GeV, mcWeight*additionalweight);
  hist("nominal/"+region+"/LogAcoplan")->Fill( log10( 1 - abs( lep1.Phi()-lep2.Phi() )/TMath::Pi() ) , mcWeight*additionalweight);
}
