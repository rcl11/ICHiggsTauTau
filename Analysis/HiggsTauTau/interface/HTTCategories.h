#ifndef ICHiggsTauTau_HiggsTauTau_HTTCategories_h
#define ICHiggsTauTau_HiggsTauTau_HTTCategories_h
#include "PhysicsTools/FWLite/interface/TFileService.h"
#include "UserCode/ICHiggsTauTau/Analysis/Core/interface/TreeEvent.h"
#include "UserCode/ICHiggsTauTau/Analysis/Core/interface/ModuleBase.h"
#include "UserCode/ICHiggsTauTau/Analysis/Utilities/interface/HTTPlots.h"
#include "UserCode/ICHiggsTauTau/Analysis/HiggsTauTau/interface/HTTConfig.h"
#include "UserCode/ICHiggsTauTau/Analysis/Utilities/interface/HistoSet.h"

#include <string>

namespace ic {

class HTTCategories : public ModuleBase {

 private:
  CLASS_MEMBER(HTTCategories, std::string, ditau_label)
  CLASS_MEMBER(HTTCategories, std::string, met_label)
  CLASS_MEMBER(HTTCategories, std::string, jets_label)
  CLASS_MEMBER(HTTCategories, double, mass_shift)
  CLASS_MEMBER(HTTCategories, ic::channel, channel)
  CLASS_MEMBER(HTTCategories, ic::era, era)
  CLASS_MEMBER(HTTCategories, ic::strategy, strategy)
  CLASS_MEMBER(HTTCategories, bool, write_tree)
  CLASS_MEMBER(HTTCategories, bool, bjet_regression)
  CLASS_MEMBER(HTTCategories, bool, make_sync_ntuple)
  CLASS_MEMBER(HTTCategories, bool, is_embedded)
  CLASS_MEMBER(HTTCategories, bool, is_data)
  CLASS_MEMBER(HTTCategories, bool, systematic_shift)
  CLASS_MEMBER(HTTCategories, bool, add_Hhh_variables)
  CLASS_MEMBER(HTTCategories, std::string, sync_output_name)
  CLASS_MEMBER(HTTCategories, bool, iso_study)
  CLASS_MEMBER(HTTCategories, bool, tau_id_study)
  CLASS_MEMBER(HTTCategories, bool, qcd_study)
  CLASS_MEMBER(HTTCategories, int, kinfit_mode )
  CLASS_MEMBER(HTTCategories, fwlite::TFileService*, fs)
 
  TTree *outtree_;
  TTree *synctree_;
  TFile *lOFile;

  struct branch_var {
      double var_double;  
      float var_float;
      branch_var& operator=(double doub){
          var_double = doub;
          var_float = static_cast<float>(doub);
          return *this;
      }
      branch_var& operator*(double sf){
          var_double*=sf;
          var_float*=sf;
          return *this;
      }
  };

  // Event Properties
  branch_var wt_;
  int run_;
  unsigned long long event_;
  int lumi_;
  float rho_;
  double mpt_1_;
  double tpt_1_;
  double meta_1_;
  double teta_1_;
  double mphi_1_;
  double tphi_1_;
  double mE_1_;
  double tE_1_;
  double mq_1_;
  double tq_1_;
  double mvamet_;
  int n_jets_, n_bjets_;
  double iso_1_, iso_2_;
  double mt_1_, mt_2_;
  double pt_tt_;
  double m_vis_;
  double dxy_1_, dxy_2_;
  double dz_1_, dz_2_;
  double jpt_1_, jpt_2_;
  double bpt_1_, bpt_2_;
  double jeta_1_, jeta_2_;
  double beta_1_, beta_2_;
  double jphi_1_, jphi_2_;
  double bphi_1_, bphi_2_;
  double deta_, dphi_;
  double d0_1_, d0_2_;

 public:
  HTTCategories(std::string const& name);
  virtual ~HTTCategories();

  virtual int PreAnalysis();
  virtual int Execute(TreeEvent *event);
  virtual int PostAnalysis();
  virtual void PrintInfo();






};

}


#endif
