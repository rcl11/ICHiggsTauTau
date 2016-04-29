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
      outtree_->Branch("wt_btag",           &wt_btag_);
      outtree_->Branch("os",                &os_);
      outtree_->Branch("m_sv",              &m_sv_.var_double);
      outtree_->Branch("mt_sv",              &mt_sv_.var_double);
      outtree_->Branch("m_vis",             &m_vis_.var_double);
      outtree_->Branch("pt_h",              &pt_h_.var_double);
      outtree_->Branch("pt_tt",             &pt_tt_.var_double);
      outtree_->Branch("mt_tot",            &mt_tot_.var_double);
      outtree_->Branch("mt_lep",            &mt_lep_.var_double);
      outtree_->Branch("mt_2",              &mt_2_.var_double);
      outtree_->Branch("mt_1",              &mt_1_.var_double);
      outtree_->Branch("m_2",               &m_2_.var_double);
      outtree_->Branch("pfmt_1",            &pfmt_1_.var_double);
      outtree_->Branch("puppimt_1",         &puppimt_1_.var_double);
      outtree_->Branch("pzeta",             &pzeta_.var_double);
      outtree_->Branch("pfpzeta",           &pfpzeta_.var_double);
      outtree_->Branch("puppipzeta",        &puppipzeta_.var_double);
      outtree_->Branch("iso_1",             &iso_1_.var_double);
      outtree_->Branch("iso_2",             &iso_2_.var_double);
      outtree_->Branch("iso_pho_sum_pt_2",  &lPhotonPtSum_2.var_double);
      outtree_->Branch("iso_pho_sum_pt_1",  &lPhotonPtSum_1.var_double);
      outtree_->Branch("antiele_1",         &antiele_1_);
      outtree_->Branch("antimu_1",          &antimu_1_);
      outtree_->Branch("antiele_2",         &antiele_2_);
      outtree_->Branch("antimu_2",          &antimu_2_);
      outtree_->Branch("leptonveto",        &lepton_veto_);
      outtree_->Branch("dilepton_veto",     &dilepton_veto_);
      outtree_->Branch("extraelec_veto",    &extraelec_veto_);
      outtree_->Branch("extramuon_veto",    &extramuon_veto_);
      outtree_->Branch("minimal_extraelec_veto",    &minimal_extraelec_veto_);
      outtree_->Branch("minimal_extramuon_veto",    &minimal_extramuon_veto_);
      outtree_->Branch("met",               &mvamet_.var_double);
      outtree_->Branch("pfmet",             &pfmet_.var_double);
      outtree_->Branch("n_jets",            &n_jets_);
      outtree_->Branch("n_bjets",           &n_bjets_);
      outtree_->Branch("n_loose_bjets",     &n_loose_bjets_);
      outtree_->Branch("mjj",               &mjj_.var_double);
      outtree_->Branch("n_jetsingap",       &n_jetsingap_);
      outtree_->Branch("jdeta",             &jdeta_.var_double);
      outtree_->Branch("n_lowpt_jets",      &n_lowpt_jets_);
      outtree_->Branch("n_jetsingap_lowpt", &n_jetsingap_lowpt_);
      outtree_->Branch("pt_2",              &pt_2_.var_double);
      outtree_->Branch("pt_1",              &pt_1_.var_double);
      outtree_->Branch("eta_1",             &eta_1_.var_double);
      outtree_->Branch("eta_2",             &eta_2_.var_double);
      outtree_->Branch("mjj_lowpt",         &mjj_lowpt_);
      outtree_->Branch("gen_match_1", &gen_match_1_);
      outtree_->Branch("gen_match_2", &gen_match_2_);
      outtree_->Branch("db_loose_1",&lbyLooseCombinedIsolation_1);
      outtree_->Branch("db_loose_2",&lbyLooseCombinedIsolation_2);
      outtree_->Branch("db_medium_1",&lbyMediumCombinedIsolation_1);
      outtree_->Branch("db_medium_2",&lbyMediumCombinedIsolation_2);
      outtree_->Branch("db_tight_1",&lbyTightCombinedIsolation_1);
      outtree_->Branch("db_tight_2",&lbyTightCombinedIsolation_2);
      outtree_->Branch("mva_olddm_vloose_1",&lbyVLooseIsolationMVArun2DBoldDMwLT_1);
      outtree_->Branch("mva_olddm_vloose_2",&lbyVLooseIsolationMVArun2DBoldDMwLT_2);
      outtree_->Branch("mva_olddm_loose_1",&lbyLooseIsolationMVArun2DBoldDMwLT_1);
      outtree_->Branch("mva_olddm_loose_2",&lbyLooseIsolationMVArun2DBoldDMwLT_2);
      outtree_->Branch("mva_olddm_medium_1",&lbyMediumIsolationMVArun2DBoldDMwLT_1);
      outtree_->Branch("mva_olddm_medium_2",&lbyMediumIsolationMVArun2DBoldDMwLT_2);
      outtree_->Branch("mva_olddm_tight_1",&lbyTightIsolationMVArun2DBoldDMwLT_1);
      outtree_->Branch("mva_olddm_tight_2",&lbyTightIsolationMVArun2DBoldDMwLT_2);
      outtree_->Branch("mva_olddm_vtight_1",&lbyVTightIsolationMVArun2DBoldDMwLT_1);
      outtree_->Branch("mva_olddm_vtight_2",&lbyVTightIsolationMVArun2DBoldDMwLT_2);
      outtree_->Branch("tau_decay_mode_2",    &tau_decay_mode_2_);
      outtree_->Branch("tau_decay_mode_1",    &tau_decay_mode_1_);

/*      outtree_->Branch("leading_lepton_match_pt", &leading_lepton_match_pt_);
      outtree_->Branch("subleading_lepton_match_pt",&subleading_lepton_match_pt_);
      outtree_->Branch("leading_lepton_match_DR", &leading_lepton_match_DR_);
      outtree_->Branch("subleading_lepton_match_DR",&subleading_lepton_match_DR_);*/

      outtree_->Branch("jdeta_lowpt",       &jdeta_lowpt_);
      if (channel_ == channel::em) {
        outtree_->Branch("em_gf_mva",         &em_gf_mva_);
        outtree_->Branch("wt_em_qcd",         &wt_em_qcd_);
        outtree_->Branch("wt_em_qcd_up",      &wt_em_qcd_up_);
        outtree_->Branch("wt_em_qcd_down",    &wt_em_qcd_down_);
        // outtree_->Branch("em_vbf_mva",        &em_vbf_mva_);
      }
      if(add_Hhh_variables_) { 
        outtree_->Branch("jet_csv_mjj",               &jet_csv_mjj_);
        outtree_->Branch("m_H_hh",     &m_H_hh_);
        outtree_->Branch("convergence_hh", &convergence_hh_);
        outtree_->Branch("mjj_tt",            &mjj_tt_);
        outtree_->Branch("n_jets_csv",        &n_jets_csv_);
        outtree_->Branch("n_bjets_csv",       &n_bjets_csv_);
        outtree_->Branch("jet_csvbcsv_1",     &jet_csvbcsv_1_);
        outtree_->Branch("jet_csvbcsv_2",     &jet_csvbcsv_2_);
      }
      if(iso_study_){
        //Add different isolation variables for if studying isolation
        outtree_->Branch("iso_1_db03", &iso_1_db03_);
        outtree_->Branch("iso_1_puw03", &iso_1_puw03_);
        outtree_->Branch("iso_1_puw04", &iso_1_puw04_);
        outtree_->Branch("iso_1_db03allch", &iso_1_db03allch_);
        outtree_->Branch("iso_1_db04allch", &iso_1_db04allch_);
        outtree_->Branch("iso_1_db04", &iso_1_db04_);
        outtree_->Branch("iso_1_ea03", &iso_1_ea03_);
        outtree_->Branch("iso_1_trk03", &iso_1_trk03_);
        outtree_->Branch("iso_2_db03", &iso_2_db03_);
        outtree_->Branch("iso_2_db03allch", &iso_2_db03allch_);
        outtree_->Branch("iso_2_db04allch", &iso_2_db04allch_);
        outtree_->Branch("iso_2_db04", &iso_2_db04_);
        outtree_->Branch("iso_2_ea03", &iso_2_ea03_);
        outtree_->Branch("iso_2_trk03", &iso_2_trk03_);
        outtree_->Branch("iso_2_puw03", &iso_2_puw03_);
        outtree_->Branch("iso_2_puw04", &iso_2_puw04_);
      }
 
      if(tau_id_study_){
       outtree_->Branch("mvadbnew_vloose_1",&lbyVLooseIsolationMVArun2DBnewDMwLT_1);
       outtree_->Branch("mvadbnew_vloose_2",&lbyVLooseIsolationMVArun2DBnewDMwLT_2);
       outtree_->Branch("mvadbnew_loose_1",&lbyLooseIsolationMVArun2DBnewDMwLT_1);
       outtree_->Branch("mvadbnew_loose_2",&lbyLooseIsolationMVArun2DBnewDMwLT_2);
       outtree_->Branch("mvadbnew_medium_1",&lbyMediumIsolationMVArun2DBnewDMwLT_1);
       outtree_->Branch("mvadbnew_medium_2",&lbyMediumIsolationMVArun2DBnewDMwLT_2);
       outtree_->Branch("mvadbnew_tight_1",&lbyTightIsolationMVArun2DBnewDMwLT_1);
       outtree_->Branch("mvadbnew_tight_2",&lbyTightIsolationMVArun2DBnewDMwLT_2);
       outtree_->Branch("mvadbnew_vtight_1",&lbyVTightIsolationMVArun2DBnewDMwLT_1);
       outtree_->Branch("mvadbnew_vtight_2",&lbyVTightIsolationMVArun2DBnewDMwLT_2);
       outtree_->Branch("mvadbnew_vvtight_1",&lbyVVTightIsolationMVArun2DBnewDMwLT_1);
       outtree_->Branch("mvadbnew_vvtight_2",&lbyVVTightIsolationMVArun2DBnewDMwLT_2);
       outtree_->Branch("mvadbold_vloose_1",&lbyVLooseIsolationMVArun2DBoldDMwLT_1);
       outtree_->Branch("mvadbold_vloose_2",&lbyVLooseIsolationMVArun2DBoldDMwLT_2);
       outtree_->Branch("mvadbold_loose_1",&lbyLooseIsolationMVArun2DBoldDMwLT_1);
       outtree_->Branch("mvadbold_loose_2",&lbyLooseIsolationMVArun2DBoldDMwLT_2);
       outtree_->Branch("mvadbold_medium_1",&lbyMediumIsolationMVArun2DBoldDMwLT_1);
       outtree_->Branch("mvadbold_medium_2",&lbyMediumIsolationMVArun2DBoldDMwLT_2);
       outtree_->Branch("mvadbold_tight_1",&lbyTightIsolationMVArun2DBoldDMwLT_1);
       outtree_->Branch("mvadbold_tight_2",&lbyTightIsolationMVArun2DBoldDMwLT_2);
       outtree_->Branch("mvadbold_vtight_1",&lbyVTightIsolationMVArun2DBoldDMwLT_1);
       outtree_->Branch("mvadbold_vtight_2",&lbyVTightIsolationMVArun2DBoldDMwLT_2);
       outtree_->Branch("mvadbold_vvtight_1",&lbyVVTightIsolationMVArun2DBoldDMwLT_1);
       outtree_->Branch("mvadbold_vvtight_2",&lbyVVTightIsolationMVArun2DBoldDMwLT_2);
       outtree_->Branch("mvapwnew_vloose_1",&lbyVLooseIsolationMVArun2PWnewDMwLT_1);
       outtree_->Branch("mvapwnew_vloose_2",&lbyVLooseIsolationMVArun2PWnewDMwLT_2);
       outtree_->Branch("mvapwnew_loose_1",&lbyLooseIsolationMVArun2PWnewDMwLT_1);
       outtree_->Branch("mvapwnew_loose_2",&lbyLooseIsolationMVArun2PWnewDMwLT_2);
       outtree_->Branch("mvapwnew_medium_1",&lbyMediumIsolationMVArun2PWnewDMwLT_1);
       outtree_->Branch("mvapwnew_medium_2",&lbyMediumIsolationMVArun2PWnewDMwLT_2);
       outtree_->Branch("mvapwnew_tight_1",&lbyTightIsolationMVArun2PWnewDMwLT_1);
       outtree_->Branch("mvapwnew_tight_2",&lbyTightIsolationMVArun2PWnewDMwLT_2);
       outtree_->Branch("mvapwnew_vtight_1",&lbyVTightIsolationMVArun2PWnewDMwLT_1);
       outtree_->Branch("mvapwnew_vtight_2",&lbyVTightIsolationMVArun2PWnewDMwLT_2);
       outtree_->Branch("mvapwnew_vvtight_1",&lbyVVTightIsolationMVArun2PWnewDMwLT_1);
       outtree_->Branch("mvapwnew_vvtight_2",&lbyVVTightIsolationMVArun2PWnewDMwLT_2);
       outtree_->Branch("mvapwold_vloose_1",&lbyVLooseIsolationMVArun2PWoldDMwLT_1);
       outtree_->Branch("mvapwold_vloose_2",&lbyVLooseIsolationMVArun2PWoldDMwLT_2);
       outtree_->Branch("mvapwold_loose_1",&lbyLooseIsolationMVArun2PWoldDMwLT_1);
       outtree_->Branch("mvapwold_loose_2",&lbyLooseIsolationMVArun2PWoldDMwLT_2);
       outtree_->Branch("mvapwold_medium_1",&lbyMediumIsolationMVArun2PWoldDMwLT_1);
       outtree_->Branch("mvapwold_medium_2",&lbyMediumIsolationMVArun2PWoldDMwLT_2);
       outtree_->Branch("mvapwold_tight_1",&lbyTightIsolationMVArun2PWoldDMwLT_1);
       outtree_->Branch("mvapwold_tight_2",&lbyTightIsolationMVArun2PWoldDMwLT_2);
       outtree_->Branch("mvapwold_vtight_1",&lbyVTightIsolationMVArun2PWoldDMwLT_1);
       outtree_->Branch("mvapwold_vtight_2",&lbyVTightIsolationMVArun2PWoldDMwLT_2);
       outtree_->Branch("mvapwold_vvtight_1",&lbyVVTightIsolationMVArun2PWoldDMwLT_1);
       outtree_->Branch("mvapwold_vvtight_2",&lbyVVTightIsolationMVArun2PWoldDMwLT_2);
       outtree_->Branch("puw_loose_1",&lbyLoosePileupWeightedIsolation_1);
       outtree_->Branch("puw_loose_2",&lbyLoosePileupWeightedIsolation_2);
       outtree_->Branch("puw_medium_1",&lbyMediumPileupWeightedIsolation_1);
       outtree_->Branch("puw_medium_2",&lbyMediumPileupWeightedIsolation_2);
       outtree_->Branch("puw_tight_1",&lbyTightPileupWeightedIsolation_1);
       outtree_->Branch("puw_tight_2",&lbyTightPileupWeightedIsolation_2);
       outtree_->Branch("antie_vloose_1",&lagainstElectronVLooseMVA_1);
       outtree_->Branch("antie_loose_1",&lagainstElectronLooseMVA_1);
       outtree_->Branch("antie_medium_1",&lagainstElectronMediumMVA_1); 
       outtree_->Branch("antie_tight_1",&lagainstElectronTightMVA_1);
       outtree_->Branch("antie_vtight_1",&lagainstElectronVTightMVA_1);
       outtree_->Branch("antimu_loose_1",&lagainstMuonLoose3_1);
       outtree_->Branch("antimu_tight_1",&lagainstMuonTight3_1);
       outtree_->Branch("antie_vloose_2",&lagainstElectronVLooseMVA_2);
       outtree_->Branch("antie_loose_2",&lagainstElectronLooseMVA_2);
       outtree_->Branch("antie_medium_2",&lagainstElectronMediumMVA_2); 
       outtree_->Branch("antie_tight_2",&lagainstElectronTightMVA_2);
       outtree_->Branch("antie_vtight_2",&lagainstElectronVTightMVA_2);
       outtree_->Branch("antimu_loose_2",&lagainstMuonLoose3_2);
       outtree_->Branch("antimu_tight_2",&lagainstMuonTight3_2);
       outtree_->Branch("isoPhoSumPt_2",&lPhotonPtSum_2.var_float);
       outtree_->Branch("isoPhoSumPt_1",&lPhotonPtSum_1.var_float);
       outtree_->Branch("iso_mvadb_new_1",&lbyIsolationMVArun2DBnewDMwLTraw_1.var_double);
       outtree_->Branch("iso_mvadb_old_1",&lbyIsolationMVArun2DBoldDMwLTraw_1.var_double);
       outtree_->Branch("iso_mvadb_new_2",&lbyIsolationMVArun2DBnewDMwLTraw_2.var_double);
       outtree_->Branch("iso_mvadb_old_2",&lbyIsolationMVArun2DBoldDMwLTraw_2.var_double);
       outtree_->Branch("iso_mvapw_new_1",&lbyIsolationMVArun2PWnewDMwLTraw_1.var_double);
       outtree_->Branch("iso_mvapw_old_1",&lbyIsolationMVArun2PWoldDMwLTraw_1.var_double);
       outtree_->Branch("iso_mvapw_new_2",&lbyIsolationMVArun2PWnewDMwLTraw_2.var_double);
       outtree_->Branch("iso_mvapw_old_2",&lbyIsolationMVArun2PWoldDMwLTraw_2.var_double);
       outtree_->Branch("olddm_1",&ldecayModeFindingOldDMs_1);
       outtree_->Branch("olddm_2",&ldecayModeFindingOldDMs_2);
      }
      if(qcd_study_){
        outtree_->Branch("jet_flav_1", &jet_flav_1_);
        outtree_->Branch("jet_flav_2", &jet_flav_2_);
      }

      if(channel_ == channel::tpzmm || channel_ == channel::tpzee){
        //Extra variables needed for tag and probe
        outtree_->Branch("id_1", &mva_1_.var_double);
        outtree_->Branch("id_2", &mva_2_.var_double);
        outtree_->Branch("q_1", &q_1_);
        outtree_->Branch("q_2", &q_2_);
        outtree_->Branch("dxy_1", &d0_1_.var_double);
        outtree_->Branch("dxy_2", &d0_2_.var_double);
        outtree_->Branch("dz_1", &dz_1_.var_double);
        outtree_->Branch("dz_2", &dz_2_.var_double);
        outtree_->Branch("trigger_match_1", &trigger_match_1_);
        outtree_->Branch("trigger_match_2", &trigger_match_2_);
      }
      //Variables needed for control plots need only be generated for central systematics
      if(!systematic_shift_) {
        //outtree_->Branch("wt_ggh_pt_up",      &wt_ggh_pt_up_);
        //outtree_->Branch("wt_ggh_pt_down",    &wt_ggh_pt_down_);
        //outtree_->Branch("wt_tau_fake_up",    &wt_tau_fake_up_);
        //outtree_->Branch("wt_tau_fake_down",  &wt_tau_fake_down_);
        outtree_->Branch("wt_tquark_up",      &wt_tquark_up_);
        outtree_->Branch("wt_tquark_down",    &wt_tquark_down_);
        outtree_->Branch("wt_zpt_up",         &wt_zpt_up_);
        outtree_->Branch("wt_zpt_down",       &wt_zpt_down_);
        outtree_->Branch("wt_tau_id_up",      &wt_tau_id_up_);
        outtree_->Branch("wt_tau_id_down",    &wt_tau_id_down_);
        outtree_->Branch("n_vtx",             &n_vtx_);
        outtree_->Branch("good_vtx",          &good_vtx_);
        outtree_->Branch("phi_1",             &phi_1_.var_double);
        outtree_->Branch("phi_2",             &phi_2_.var_double);
        if (channel_ != channel::em){
          outtree_->Branch("dphi",              &dphi_);
        }
        outtree_->Branch("E_1",               &E_1_);
        outtree_->Branch("E_2",               &E_2_);
        outtree_->Branch("z_2",               &z_2_);
        outtree_->Branch("met_phi",           &mvamet_phi_.var_double);
        outtree_->Branch("n_prebjets",        &n_prebjets_);
        outtree_->Branch("jpt_1",             &jpt_1_.var_double);
        outtree_->Branch("j1_dm",             &j1_dm_);
        outtree_->Branch("jpt_2",             &jpt_2_.var_double);
        outtree_->Branch("jeta_1",            &jeta_1_.var_double);
        outtree_->Branch("jeta_2",            &jeta_2_.var_double);
        outtree_->Branch("bpt_1",             &bpt_1_.var_double);
        outtree_->Branch("beta_1",            &beta_1_.var_double);
        outtree_->Branch("bcsv_1",            &bcsv_1_.var_double);
/*        outtree_->Branch("trigger_object_pt_1",&trigger_object_pt_1.var_double);
        outtree_->Branch("trigger_object_eta_1",&trigger_object_eta_1.var_double);
        outtree_->Branch("trigger_object_pt_2",&trigger_object_pt_2.var_double);
        outtree_->Branch("trigger_object_eta_2",&trigger_object_eta_2.var_double);
*/
        if (channel_ == channel::em) {
          outtree_->Branch("pzetavis",          &pzetavis_.var_double);
          outtree_->Branch("pzetamiss",         &pzetamiss_.var_double);
          outtree_->Branch("mt_ll",             &mt_ll_);
          outtree_->Branch("emu_dphi",          &dphi_);
          outtree_->Branch("emu_csv",           &emu_csv_);
          outtree_->Branch("emu_dxy_1",         &emu_dxy_1_);
          outtree_->Branch("emu_dxy_2",         &emu_dxy_2_);
        }
        if(add_Hhh_variables_) {
          outtree_->Branch("jet_csvpt_1",       &jet_csvpt_1_);
          outtree_->Branch("jet_csveta_1",      &jet_csveta_1_);
          outtree_->Branch("jet_csvpt_2",       &jet_csvpt_2_);
          outtree_->Branch("jet_csveta_2",      &jet_csveta_2_);
          outtree_->Branch("mjj_h",             &mjj_h_);
          outtree_->Branch("mbb_h",             &mbb_h_);
          if(kinfit_mode_ > 1) {
            outtree_->Branch("m_H_best",               &m_H_best_);
            outtree_->Branch("m_H_chi2_best",               &m_H_chi2_best_);
            outtree_->Branch("pull_balance_H_best", &pull_balance_H_best_);
            outtree_->Branch("convergence_H_best", &convergence_H_best_); 
            outtree_->Branch("m_H_hZ",          &m_H_hZ_);
            outtree_->Branch("m_H_hZ_chi2",     &m_H_hZ_chi2_);
            outtree_->Branch("pull_balance_hZ", &pull_balance_hZ_);
            outtree_->Branch("convergence_hZ", &convergence_hZ_);
            outtree_->Branch("m_H_Zh",          &m_H_Zh_);
            outtree_->Branch("m_H_Zh_chi2",     &m_H_Zh_chi2_);
            outtree_->Branch("pull_balance_Zh",  &pull_balance_Zh_);
            outtree_->Branch("convergence_Zh",  &convergence_Zh_);
            outtree_->Branch("m_H_hh_all",     &m_H_hh_all_);
            outtree_->Branch("m_H_hh_chi2",     &m_H_hh_chi2_);
            outtree_->Branch("pull_balance_hh", &pull_balance_hh_);
            outtree_->Branch("m_bb",     &m_bb_);
            outtree_->Branch("m_bb_chi2",     &m_bb_chi2_);
            outtree_->Branch("pull_balance_bb", &pull_balance_bb_);
            outtree_->Branch("convergence_bb", &convergence_bb_);
          }
        }
      }
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
    std::vector<PileupInfo *> puInfo;
    float true_int = -1;

    if (event->Exists("pileupInfo") || strategy_ == strategy::phys14 || ((strategy_==strategy::spring15||strategy_==strategy::fall15) && !is_data_) ) {
     puInfo = event->GetPtrVec<PileupInfo>("pileupInfo");
      for (unsigned i = 0; i < puInfo.size(); ++i) {
        if (puInfo[i]->bunch_crossing() == 0)
          true_int = puInfo[i]->true_num_interactions();
      }
    }
    n_pu_ = true_int;
    rho_ = eventInfo->jet_rho();
    if(event->Exists("gen_match_1")) gen_match_1_ = MCOrigin2UInt(event->Get<ic::mcorigin>("gen_match_1"));
    if(event->Exists("gen_match_2")) gen_match_2_ = MCOrigin2UInt(event->Get<ic::mcorigin>("gen_match_2"));
  
    
 //   Met const* mets = NULL;
 //   mets = event->GetPtr<Met>(met_label_);
    
    std::vector<Muon *> const& muons = event->GetPtrVec<Muon>("sel_muons");
    Muon const* lep1 = muons.at(0);
    Muon const* lep2 = NULL;
    if(muons.size()>1) lep2 = muons.at(1);
   


    std::vector<PFJet*> jets = event->GetPtrVec<PFJet>(jets_label_);
    std::vector<PFJet*> uncleaned_jets = event->GetPtrVec<PFJet>(jets_label_+"UnFiltered");
    std::vector<PFJet*> corrected_jets;
    if(bjet_regression_) corrected_jets = event->GetPtrVec<PFJet>(jets_label_+"Corrected");
    std::sort(jets.begin(), jets.end(), bind(&Candidate::pt, _1) > bind(&Candidate::pt, _2));
    std::vector<PFJet*> lowpt_jets = jets;
    ic::erase_if(jets,!boost::bind(MinPtMaxEta, _1, 30.0, 4.7));
    ic::erase_if(lowpt_jets,!boost::bind(MinPtMaxEta, _1, 20.0, 4.7));
    std::vector<PFJet*> prebjets = lowpt_jets;
    ic::erase_if(prebjets,!boost::bind(MinPtMaxEta, _1, 20.0, 2.4));
    std::vector<PFJet*> bjets = prebjets;
    std::vector<PFJet*> loose_bjets = prebjets;
    std::string btag_label="combinedSecondaryVertexBJetTags";
    double btag_wp =  0.679;
    double loose_btag_wp = 0.244;
    if(strategy_ == strategy::phys14) btag_label = "combinedInclusiveSecondaryVertexV2BJetTags";
    if(strategy_ == strategy::phys14) btag_wp = 0.814 ;
    if(strategy_ == strategy::spring15) btag_label = "pfCombinedInclusiveSecondaryVertexV2BJetTags";
    if(strategy_ == strategy::spring15) btag_wp = 0.89 ;
    if(strategy_ == strategy::fall15) btag_label = "pfCombinedInclusiveSecondaryVertexV2BJetTags";
    if(strategy_ == strategy::fall15) btag_wp = 0.8;
    if(strategy_ == strategy::fall15) loose_btag_wp = 0.46;


    //Sort out the loose (em,mt,et) or medium (tt) b-jets
    if(channel_!= channel::tt){
      ic::erase_if(loose_bjets, boost::bind(&PFJet::GetBDiscriminator, _1, btag_label) < loose_btag_wp);
    } else {
      ic::erase_if(bjets, boost::bind(&PFJet::GetBDiscriminator, _1, btag_label) <btag_wp);
    }


    // Instead of changing b-tag value in the promote/demote method we look for a map of bools
    // that say whether a jet should pass the WP or not
    if (event->Exists("retag_result")) {
      auto const& retag_result = event->Get<std::map<std::size_t,bool>>("retag_result"); 
       if(channel_ != channel::tt){
          ic::erase_if(bjets, !boost::bind(IsReBTagged, _1, retag_result));
       } else {
          ic::erase_if(loose_bjets, !boost::bind(IsReBTagged, _1, retag_result));
      }
    } else{ 
      if(channel_ != channel::tt){
        ic::erase_if(bjets, boost::bind(&PFJet::GetBDiscriminator, _1, btag_label) < btag_wp);
      } else {
        ic::erase_if(loose_bjets, boost::bind(&PFJet::GetBDiscriminator, _1, btag_label) < loose_btag_wp);
      }
    } 
    
    // Define event properties
    // IMPORTANT: Make sure each property is re-set
    // for each new event
/*    if (PairOppSign(ditau)) {
      os_ = true;
    } else {
      os_ = false;
    }*/


    pt_1_ = lep1->pt();
    if(muons.size()>1) pt_2_ = lep2->pt();
    else pt_2_ = -9999;
/*    eta_1_ = lep1->eta();
    eta_2_ = lep2->eta();
    phi_1_ = lep1->phi();
    phi_2_ = lep2->phi();
    dphi_ = std::fabs(ROOT::Math::VectorUtil::DeltaPhi(lep1->vector(),lep2->vector()));
    E_1_ = lep1->energy();
    E_2_ = lep2->energy();
    m_1_ = lep1->M();
    m_2_ = lep2->M();
    q_1_ = lep1->charge();
    q_2_ = lep2->charge();
    if(make_sync_ntuple_){
      event->Exists("genpX") ? gen_px_ = event->Get<double>("genpX") : 0.;
      event->Exists("genpY") ? gen_py_ = event->Get<double>("genpY") : 0.;
      event->Exists("vispX") ? vis_px_ = event->Get<double>("vispX") : 0.;
      event->Exists("vispY") ? vis_py_ = event->Get<double>("vispY") : 0.;
    }
    mvamet_ = mets->pt();
    mvamet_phi_ = mets->phi();
    
    emu_dxy_1_ = 0.0;
    emu_dxy_2_ = 0.0;
    
    antiele_1_ = true;
    antimu_1_ = true;
    antiele_2_ = true;
    antimu_2_ = true;
    
    
    if (channel_ == channel::mt || channel_ == channel::mtmet) {
      Muon const* muon1 = dynamic_cast<Muon const*>(lep1);
      Muon const* muon2 = dynamic_cast<Muon const*>(lep2);
      d0_1_ = muon1->dxy_vertex();
      dz_1_ = muon1->dz_vertex();
      d0_2_ = muon2->dxy_vertex();
      dz_2_ = muon2->dz_vertex();
    }
*/

    n_jets_ = jets.size();
    n_lowpt_jets_ = lowpt_jets.size();
    n_bjets_ = bjets.size();
    n_prebjets_ = prebjets.size();
    n_loose_bjets_ = loose_bjets.size();


    if (n_lowpt_jets_ >= 1) {
      jpt_1_ = lowpt_jets[0]->pt();
      jeta_1_ = lowpt_jets[0]->eta();
      jphi_1_ = lowpt_jets[0]->phi();
      jrawf_1_ = lowpt_jets[0]->uncorrected_energy()/jets[0]->energy();//* (jets[0]->pt() / jets[0]->energy());
      jptunc_1_ = 0.0;
      jmva_1_ = lowpt_jets[0]->pu_id_mva_value();
      jlrm_1_ = lowpt_jets[0]->linear_radial_moment();
      jctm_1_ = lowpt_jets[0]->charged_multiplicity_nopu();
      std::vector<ic::Tau *> taus = event->GetPtrVec<Tau>("taus");
      std::vector<ic::Jet *> leadjet = { jets[0] };
      std::vector<std::pair<ic::Jet *, ic::Tau *>> matches = MatchByDR(leadjet, taus, 0.5, true, true);
      if (matches.size() == 1) {
        j1_dm_ = matches[0].second->decay_mode();
      } else {
        j1_dm_ = -1;
      }
    } else {
      jpt_1_ = -9999;
      jeta_1_ = -9999;
      jphi_1_ = -9999;
      jrawf_1_ = -9999;
      jptunc_1_ = -9999;
      jmva_1_ = -9999;
      jlrm_1_ = -9999;
      jctm_1_ = -9999;
    }

    if (n_lowpt_jets_ >= 2) {
      jpt_2_ = lowpt_jets[1]->pt();
      jeta_2_ = lowpt_jets[1]->eta();
      jphi_2_ = lowpt_jets[1]->phi();
      jrawf_2_ = lowpt_jets[1]->uncorrected_energy()/jets[1]->energy();// * (jets[1]->pt() / jets[1]->energy());
      jptunc_2_ = 0.0;
      jmva_2_ = lowpt_jets[1]->pu_id_mva_value();
      jlrm_2_ = lowpt_jets[1]->linear_radial_moment();
      jctm_2_ = lowpt_jets[1]->charged_multiplicity_nopu();
      mjj_ = (lowpt_jets[0]->vector() + lowpt_jets[1]->vector()).M();
      jdeta_ = fabs(lowpt_jets[0]->eta() - lowpt_jets[1]->eta());
      jdphi_ =  std::fabs(ROOT::Math::VectorUtil::DeltaPhi(lowpt_jets[0]->vector(), lowpt_jets[1]->vector()));
      double eta_high = (lowpt_jets[0]->eta() > lowpt_jets[1]->eta()) ? lowpt_jets[0]->eta() : lowpt_jets[1]->eta();
      double eta_low = (lowpt_jets[0]->eta() > lowpt_jets[1]->eta()) ? lowpt_jets[1]->eta() : lowpt_jets[0]->eta();
      n_jetsingap_ = 0;
      n_jetsingap20_ = 0;
      if (n_lowpt_jets_ > 2) {
        for (unsigned i = 2; i < lowpt_jets.size(); ++i) {
         if (lowpt_jets[i]->pt() > 30.0 &&  lowpt_jets[i]->eta() > eta_low && lowpt_jets[i]->eta() < eta_high) ++n_jetsingap_;
         if (lowpt_jets[i]->pt() > 20.0 &&  lowpt_jets[i]->eta() > eta_low && lowpt_jets[i]->eta() < eta_high) ++n_jetsingap20_;
        }
      }
    } else {
      jpt_2_ = -9999;
      jeta_2_ = -9999;
      jphi_2_ = -9999;
      mjj_ = -9999;
      jdeta_ = -9999;
      jdphi_ = -9999;
      jrawf_2_ = -9999;
      jptunc_2_ = -9999;
      jmva_2_ = -9999;
      jlrm_2_ = -9999;
      jctm_2_ = -9999;
      n_jetsingap_ = 9999;
      n_jetsingap20_ = 9999;
    }

    if (n_lowpt_jets_ >= 2) {
      mjj_lowpt_ = (lowpt_jets[0]->vector() + lowpt_jets[1]->vector()).M();
      jdeta_lowpt_ = fabs(lowpt_jets[0]->eta() - lowpt_jets[1]->eta());
      double eta_high = (lowpt_jets[0]->eta() > lowpt_jets[1]->eta()) ? lowpt_jets[0]->eta() : lowpt_jets[1]->eta();
      double eta_low = (lowpt_jets[0]->eta() > lowpt_jets[1]->eta()) ? lowpt_jets[1]->eta() : lowpt_jets[0]->eta();
      n_jetsingap_lowpt_ = 0;
      if (n_lowpt_jets_ > 2) {
        for (unsigned i = 2; i < lowpt_jets.size(); ++i) {
         if (lowpt_jets[i]->pt() > 30.0 &&  lowpt_jets[i]->eta() > eta_low && lowpt_jets[i]->eta() < eta_high) ++n_jetsingap_lowpt_;
        }
      }
    } else {
      mjj_lowpt_ = -9999;
      jdeta_lowpt_ = -9999;
      n_jetsingap_lowpt_ = 9999;
    }

    if (channel_ != channel::tt){
      if (n_bjets_ >= 1) {
        bpt_1_ = bjets[0]->pt();
        brawf_1_ = bjets[0]->uncorrected_energy()/bjets[0]->energy();//* (jets[0]->pt() / jets[0]->energy());
        beta_1_ = bjets[0]->eta();
        bphi_1_ = bjets[0]->phi();
        bmva_1_ = bjets[0]->pu_id_mva_value();
      
      } else {
        bpt_1_ = -9999;
        brawf_1_ = -9999;
        beta_1_ = -9999;
        bphi_1_ = -9999;
        bmva_1_ = -9999;
      }

      if (n_bjets_ >= 2) {
        bpt_2_ = bjets[1]->pt();
        brawf_2_ = bjets[1]->uncorrected_energy()/bjets[1]->energy();//* (jets[0]->pt() / jets[0]->energy());
        beta_2_ = bjets[1]->eta();
        bphi_2_ = bjets[1]->phi();
        bmva_2_ = bjets[1]->pu_id_mva_value();
      
      } else {
        bpt_2_ = -9999;
        brawf_2_ = -9999;
        beta_2_ = -9999;
        bphi_2_ = -9999;
        bmva_2_ = -9999;
      }
    } else {//We use the loose CSV wp for fully hadronic, so adjust definitions accordingly
      if (n_loose_bjets_ >= 1) {
        bpt_1_ = loose_bjets[0]->pt();
        brawf_1_ = loose_bjets[0]->uncorrected_energy()/loose_bjets[0]->energy();//* (jets[0]->pt() / jets[0]->energy());
        beta_1_ = loose_bjets[0]->eta();
        bphi_1_ = loose_bjets[0]->phi();
        bmva_1_ = loose_bjets[0]->pu_id_mva_value();
      
      } else {
        bpt_1_ = -9999;
        brawf_1_ = -9999;
        beta_1_ = -9999;
        bphi_1_ = -9999;
        bmva_1_ = -9999;
      }

      if (n_loose_bjets_ >= 2) {
        bpt_2_ = loose_bjets[1]->pt();
        brawf_2_ = loose_bjets[1]->uncorrected_energy()/loose_bjets[1]->energy();//* (jets[0]->pt() / jets[0]->energy());
        beta_2_ = loose_bjets[1]->eta();
        bphi_2_ = loose_bjets[1]->phi();
        bmva_2_ = loose_bjets[1]->pu_id_mva_value();
      
      } else {
        bpt_2_ = -9999;
        brawf_2_ = -9999;
        beta_2_ = -9999;
        bphi_2_ = -9999;
        bmva_2_ = -9999;
      }

    }



    if (n_prebjets_ >= 1) {
      bcsv_1_ = prebjets[0]->GetBDiscriminator(btag_label);
    } else {
      bcsv_1_ = -9999;
    }
    if (n_prebjets_ >= 2) {
      bcsv_2_ = prebjets[1]->GetBDiscriminator(btag_label);
    } else {
      bcsv_2_ = -9999;
    }

    emu_csv_ = (bcsv_1_.var_double > 0.244) ? bcsv_1_.var_double : -1.0;


    
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
