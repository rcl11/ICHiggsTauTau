#include "UserCode/ICHiggsTauTau/Analysis/HiggsTauTau/interface/HTTCategories.h"
#include "UserCode/ICHiggsTauTau/Analysis/HiggsTauTau/interface/HTTConfig.h"
#include "UserCode/ICHiggsTauTau/interface/PFJet.hh"
#include "UserCode/ICHiggsTauTau/Analysis/Utilities/interface/FnPredicates.h"
#include "UserCode/ICHiggsTauTau/Analysis/Utilities/interface/FnPairs.h"
#include "UserCode/ICHiggsTauTau/Analysis/HiggsTauTau/HHKinFit/include/HHKinFitMaster.h"
#include "UserCode/ICHiggsTauTau/Analysis/HiggsTauTau/HHKinFit/include/HHDiJetKinFitMaster.h"

#include "TMVA/Reader.h"
#include "TVector3.h"
#include "boost/format.hpp"
#include "TMath.h"
#include "TLorentzVector.h"

namespace ic {

  HTTCategories::HTTCategories(std::string const& name) : ModuleBase(name), 
      channel_(channel::et), 
      era_(era::data_2012_rereco),
      strategy_(strategy::paper2013) {
      ditau_label_ = "emtauCandidates";
      jets_label_ = "pfJetsPFlow";
      met_label_ = "pfMVAMetNoLeptons";
      mass_shift_ = 1.0;
      fs_ = NULL;
      write_tree_ = true;
      bjet_regression_ = false;
      make_sync_ntuple_ = false;
      sync_output_name_ = "SYNC.root";
      iso_study_=false;
      tau_id_study_=false;
      is_embedded_=false;
      is_data_=false;
      qcd_study_=false;
      kinfit_mode_ = 0; //0 = don't run, 1 = run simple 125,125 default fit, 2 = run extra masses default fit, 3 = run m_bb only fit
      systematic_shift_ = false;
      add_Hhh_variables_ = false; //set to include custom variables for the H->hh analysis
  }

  HTTCategories::~HTTCategories() {
    ;
  }

  int HTTCategories::PreAnalysis() {
    std::cout << "-------------------------------------" << std::endl;
    std::cout << "HTTCategories" << std::endl;
    std::cout << "-------------------------------------" << std::endl;    
      std::cout << boost::format(param_fmt()) % "channel"         % Channel2String(channel_);
      std::cout << boost::format(param_fmt()) % "strategy"        % Strategy2String(strategy_);
      std::cout << boost::format(param_fmt()) % "era"             % Era2String(era_);
      std::cout << boost::format(param_fmt()) % "dilepton_label"  % ditau_label_;
      std::cout << boost::format(param_fmt()) % "met_label"       % met_label_;
      std::cout << boost::format(param_fmt()) % "jets_label"      % jets_label_;
      std::cout << boost::format(param_fmt()) % "mass_shift"      % mass_shift_;
      std::cout << boost::format(param_fmt()) % "write_tree"      % write_tree_;
      std::cout << boost::format(param_fmt()) % "kinfit_mode"     % kinfit_mode_;
      std::cout << boost::format(param_fmt()) % "make_sync_ntuple" % make_sync_ntuple_;
      std::cout << boost::format(param_fmt()) % "bjet_regression" % bjet_regression_;

    if (fs_ && write_tree_) {
      outtree_ = fs_->make<TTree>("ntuple","ntuple");
      outtree_->Branch("event",             &event_);
      outtree_->Branch("wt",                &wt_);
      outtree_->Branch("mpt_1",              &mpt_1_);
      outtree_->Branch("tpt_1",              &tpt_1_);
      outtree_->Branch("meta_1",              &meta_1_);
      outtree_->Branch("teta_1",              &teta_1_);
      outtree_->Branch("mphi_1",              &mphi_1_);
      outtree_->Branch("tphi_1",              &tphi_1_);
      outtree_->Branch("mE_1",              &mE_1_);
      outtree_->Branch("tE_1",              &tE_1_);
      outtree_->Branch("mq_1",              &mq_1_);
      outtree_->Branch("tq_1",              &tq_1_);
      outtree_->Branch("met",              &mvamet_);
      outtree_->Branch("n_jets",            &n_jets_);
      outtree_->Branch("n_bjets",           &n_bjets_);
      outtree_->Branch("miso_1",             &iso_1_);
      outtree_->Branch("tiso_2",             &iso_2_);
      outtree_->Branch("mt_2",              &mt_2_);
      outtree_->Branch("mt_1",              &mt_1_);
      outtree_->Branch("pt_tt",             &pt_tt_);
      outtree_->Branch("m_vis",             &m_vis_);
      outtree_->Branch("mdxy_1",             &d0_1_);
      outtree_->Branch("tdxy_1",             &d0_2_);
      outtree_->Branch("mdz_1",              &dz_1_);
      outtree_->Branch("tdz_1",              &dz_2_);
      outtree_->Branch("jpt_1",             &jpt_1_);
      outtree_->Branch("jeta_1",            &jeta_1_);
      outtree_->Branch("jphi_1",            &jeta_1_);
      outtree_->Branch("bpt_1",             &bpt_1_);
      outtree_->Branch("beta_1",            &beta_1_);
      outtree_->Branch("bphi_1",            &beta_1_);
      outtree_->Branch("dphi",            &dphi_);
      outtree_->Branch("deta",            &deta_);
    }

    return 0;
  }

  int HTTCategories::Execute(TreeEvent *event) {

    // Get the objects we need from the event
    EventInfo const* eventInfo = event->GetPtr<EventInfo>("eventInfo");
    
    wt_ = {eventInfo->total_weight(), static_cast<float>(eventInfo->total_weight())};
    run_ = eventInfo->run();
    event_ = (unsigned long long) eventInfo->event();
    lumi_ = eventInfo->lumi_block();
 

    std::vector<CompositeCandidate *> const& ditau_vec = event->GetPtrVec<CompositeCandidate>(ditau_label_);
    CompositeCandidate const* ditau = ditau_vec.at(0);
    Candidate const* lep1 = ditau->GetCandidate("lepton1");
    Candidate const* lep2 = ditau->GetCandidate("lepton2");

    
    Met const* mets = NULL;
    mets = event->GetPtr<Met>(met_label_);

    std::vector<PFJet*> jets = event->GetPtrVec<PFJet>(jets_label_);
    std::vector<PFJet*> bjets = jets;
    ic::erase_if(bjets,!boost::bind(MinPtMaxEta, _1, 20.0, 2.4));
    double btag_wp;
    std::string btag_label;
    if(strategy_ == strategy::fall15) btag_label = "pfCombinedInclusiveSecondaryVertexV2BJetTags";
    if(strategy_ == strategy::fall15) btag_wp = 0.8;
    ic::erase_if(bjets, boost::bind(&PFJet::GetBDiscriminator, _1, btag_label) <btag_wp);
    
    Muon const* muon = dynamic_cast<Muon const*>(lep1);
    Tau const* tau = dynamic_cast<Tau const*>(lep2);

    mpt_1_ = muon->pt();
    meta_1_ = muon->eta();
    mphi_1_ = muon->phi();
    mE_1_ = muon->energy();
    mq_1_ = muon->charge();

    tpt_1_ = tau->pt();
    teta_1_ = tau->eta();
    tphi_1_ = tau->phi();
    tE_1_ = tau->energy();
    tq_1_ = tau->charge();

    mvamet_ = mets->pt();

    iso_1_ = PF04IsolationVal(muon, 0.5);
    iso_2_ = tau->GetTauID("byCombinedIsolationDeltaBetaCorrRaw3Hits");
    d0_1_ = muon->dxy_vertex();
    dz_1_ = muon->dz_vertex();
    d0_2_ = tau->lead_dxy_vertex();
    dz_2_ = tau->lead_dz_vertex();
  
    m_vis_ = ditau->M();
    pt_tt_ = (ditau->vector() + mets->vector()).pt();
    mt_1_ = MT(lep1, mets);
    mt_2_ = MT(lep2, mets);
    dphi_ = std::fabs(ROOT::Math::VectorUtil::DeltaPhi(lep1->vector(),lep2->vector()));
    deta_ = std::fabs(meta_1_ - teta_1_);

    n_jets_ = jets.size();
    n_bjets_ = bjets.size();

    if (n_jets_ >= 1) {
      jpt_1_ = jets[0]->pt();
      jeta_1_ = jets[0]->eta();
      jphi_1_ = jets[0]->phi();
    } else {
      jpt_1_ = -9999;
      jeta_1_ = -9999;
      jphi_1_ = -9999;
    }
    
if (n_bjets_ >= 1) {
      bpt_1_ = bjets[0]->pt();
      beta_1_ = bjets[0]->eta();
      bphi_1_ = bjets[0]->phi();
    } else {
      bpt_1_ = -9999;
      beta_1_ = -9999;
      bphi_1_ = -9999;
    }

    
    
    if (write_tree_ && fs_) outtree_->Fill();
    if (make_sync_ntuple_) synctree_->Fill();


    return 0;
  }

  int HTTCategories::PostAnalysis() {
    if(make_sync_ntuple_) {   
      lOFile->cd();
      synctree_->Write();
      lOFile->Close();
    }
    return 0;
  }

  void HTTCategories::PrintInfo() {
    ;
  }
}
