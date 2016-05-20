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
      outtree_->Branch("wt",                &wt_.var_double);
      outtree_->Branch("mpt_1",              &mpt_1_);
      outtree_->Branch("mpt_2",              &mpt_2_);
      outtree_->Branch("mpt_3",              &mpt_3_);
      outtree_->Branch("tpt_1",              &tpt_1_);
      outtree_->Branch("tpt_2",              &tpt_2_);
      outtree_->Branch("tpt_3",              &tpt_3_);
      outtree_->Branch("meta_1",              &meta_1_);
      outtree_->Branch("meta_2",              &meta_2_);
      outtree_->Branch("meta_3",              &meta_3_);
      outtree_->Branch("teta_1",              &teta_1_);
      outtree_->Branch("teta_2",              &teta_2_);
      outtree_->Branch("teta_3",              &teta_3_);
      outtree_->Branch("mphi_1",              &mphi_1_);
      outtree_->Branch("mphi_2",              &mphi_2_);
      outtree_->Branch("mphi_3",              &mphi_3_);
      outtree_->Branch("tphi_1",              &tphi_1_);
      outtree_->Branch("tphi_2",              &tphi_2_);
      outtree_->Branch("tphi_3",              &tphi_3_);
      outtree_->Branch("mE_1",              &mE_1_);
      outtree_->Branch("mE_2",              &mE_2_);
      outtree_->Branch("mE_3",              &mE_3_);
      outtree_->Branch("tE_1",              &tE_1_);
      outtree_->Branch("tE_2",              &tE_2_);
      outtree_->Branch("tE_3",              &tE_3_);
      outtree_->Branch("mq_1",              &mq_1_);
      outtree_->Branch("mq_2",              &mq_2_);
      outtree_->Branch("mq_3",              &mq_3_);
      outtree_->Branch("tq_1",              &tq_1_);
      outtree_->Branch("tq_2",              &tq_2_);
      outtree_->Branch("tq_3",              &tq_3_);
      outtree_->Branch("met",              &mvamet_);
        
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
  
    
    Met const* mets = NULL;
    mets = event->GetPtr<Met>(met_label_);
    
    std::vector<Muon *> const& muons = event->GetPtrVec<Muon>("sel_muons");
    Muon const* mu1 = muons.at(0);
    Muon const* mu2 = NULL;
    Muon const* mu3 = NULL;
    if(muons.size()>1) mu2 = muons.at(1);
    if(muons.size()>2) mu3 = muons.at(2);
   
    mpt_1_ = mu1->pt();
    meta_1_ = mu1->eta();
    mphi_1_ = mu1->phi();
    mE_1_ = mu1->energy();
    mq_1_ = mu1->charge();
    if(muons.size()>1) {
        mpt_2_ = mu2->pt();
        meta_2_ = mu2->eta();
        mphi_2_ = mu2->phi();
        mE_2_ = mu2->energy();
        mq_2_ = mu2->charge();
    } else {
        mpt_2_ = -9999;
        meta_2_ = -9999;
        mphi_2_ = -9999;
        mE_2_ = -9999;
        mq_2_ = -9999;
    }
    if(muons.size()>2) {
        mpt_3_ = mu3->pt();
        meta_3_ = mu3->eta();
        mphi_3_ = mu3->phi();
        mE_3_ = mu3->energy();
        mq_3_ = mu3->charge();
    } else {
        mpt_3_ = -9999;
        meta_3_ = -9999;
        mphi_3_ = -9999;
        mE_3_ = -9999;
        mq_3_ = -9999;
    }

    std::vector<Tau *> const& taus = event->GetPtrVec<Tau>("taus");
    Tau const* tau1 = taus.at(0);
    Tau const* tau2 = NULL;
    Tau const* tau3 = NULL;
    if(taus.size()>1) tau2 = taus.at(1);
    if(taus.size()>2) tau3 = taus.at(2);
   

    tpt_1_ = tau1->pt();
    teta_1_ = tau1->eta();
    tphi_1_ = tau1->phi();
    tE_1_ = tau1->energy();
    tq_1_ = tau1->charge();
    if(taus.size()>1) {
        tpt_2_ = tau2->pt();
        teta_2_ = tau2->eta();
        tphi_2_ = tau2->phi();
        tE_2_ = tau2->energy();
        tq_2_ = tau2->charge();
    } else {
        tpt_2_ = -9999;
        teta_2_ = -9999;
        tphi_2_ = -9999;
        tE_2_ = -9999;
        tq_2_ = -9999;
    }
    if(taus.size()>2) {
        tpt_3_ = tau3->pt();
        teta_3_ = tau3->eta();
        tphi_3_ = tau3->phi();
        tE_3_ = tau3->energy();
        tq_3_ = tau3->charge();
    } else {
        tpt_3_ = -9999;
        teta_3_ = -9999;
        tphi_3_ = -9999;
        tE_3_ = -9999;
        tq_3_ = -9999;
    }

    mvamet_ = mets->pt();
    
    
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
