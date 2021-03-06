#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include "boost/lexical_cast.hpp"
#include "boost/algorithm/string.hpp"
#include "boost/format.hpp"
#include "boost/program_options.hpp"
#include "UserCode/ICHiggsTauTau/Analysis/Core/interface/Plot.h"
#include "UserCode/ICHiggsTauTau/Analysis/Core/interface/TextElement.h"
#include "UserCode/ICHiggsTauTau/Analysis/Utilities/interface/SimpleParamParser.h"
#include "TH1.h"
#include "TFile.h"
#include "UserCode/ICHiggsTauTau/Analysis/HiggsNuNu/interface/HinvConfig.h"
#include "UserCode/ICHiggsTauTau/Analysis/HiggsNuNu/interface/Utilities.h"

namespace po = boost::program_options;

using namespace std;
using namespace ic;

void SumHistograms(std::string lName, ic::TH1PlotElement & lEle,std::string lStr, ic::TH1PlotElement & lHist){
  if (lName.find(lStr) != lName.npos) {
    if (!lHist.hist_ptr()) lHist.SetElement((TH1F*)(lEle.hist_ptr()->Clone()));
    else (lHist.hist_ptr())->Add((TH1F*)(lEle.hist_ptr()->Clone()));
  }
}

string Token(string const& file, string const& selection) {
	return (file+"_"+selection);
}

double Integral(TH1F const* hist) {
  if (hist) {
    double ltmp =hist->Integral(0, hist->GetNbinsX() + 1); 
    if (ltmp<0 || ltmp != ltmp) {
      std::cout << " -- Warning: integral is " << ltmp << ". Setting to 0." << std::endl;
      ltmp=0;
    }
    return ltmp;
  }
  else return 0;
}

double Error(TH1F const* hist) {
  double err = 0.0;
  if (hist) {
    hist->IntegralAndError(0, hist->GetNbinsX() + 1, err);
    if (err<0 || err != err) {
      std::cout << " -- Warning: error on integral is " << err << ". Setting to 0." << std::endl;
      err=0;
    }
  }
  return err;
}

double sumSqErr(double n1, double n2, 
		double n3 = 0, double n4 = 0, 
		double n5 = 0, double n6 = 0,
		double n7 = 0, double n8 = 0){
  return sqrt(pow(n1,2)+pow(n2,2)+
	      pow(n3,2)+pow(n4,2)+
	      pow(n5,2)+pow(n6,2)+
	      pow(n7,2)+pow(n8,2)
	      );
}

void SetBkgStyle(ic::TH1PlotElement & ele, unsigned color) {
  ele.set_marker_color(color);
  ele.set_line_color(color);
  ele.set_fill_color(color);
  ele.set_fill_style(1001);
  ele.set_draw_fill(true);
  ele.set_draw_marker(false);
  ele.set_draw_line(false);
  ele.set_line_width(2);
  ele.set_draw_stat_error_y(false);
  ele.set_in_stack(true);
  return;
}


void SetSignalStyle(ic::TH1PlotElement & ele, unsigned color) {
  ele.set_marker_color(color);
  ele.set_marker_style(1);
  ele.set_line_color(color);
  ele.set_line_style(2);
  ele.set_fill_color(color);
  ele.set_fill_style(0);
  ele.set_draw_fill(true);
  ele.set_draw_marker(false);
  ele.set_draw_line(true);
  ele.set_smooth_curve(false);
  ele.set_line_width(3);
  ele.set_draw_stat_error_y(false);
  ele.set_in_stack(false);
  return;
}


void SetDataStyle(ic::TH1PlotElement & ele) {
  ele.set_marker_color(1);
  ele.set_line_color(1);
  ele.set_fill_color(1);
  ele.set_fill_style(0);
  ele.set_draw_fill(false);
  ele.set_line_width(3);
  ele.set_draw_marker(true);
  ele.set_draw_line(true);
  ele.set_marker_style(20);
  ele.set_draw_stat_error_y(true);
  ele.set_marker_size(1.3);

  //ele.set_draw_normalised(true);

  //ele.set_in_stack(true);
  return;
}

void SetStyle(ic::RatioPlotElement & ele, unsigned color) {
  ele.set_marker_color(color);
  ele.set_line_color(color);
  ele.set_draw_stat_error_y(true);
  ele.set_draw_line(false);
  ele.set_draw_marker(true);
  ele.set_marker_style(20);

  return;
}


int main(int argc, char* argv[]){

	// Configurable parameters
  string cfg;			     // The configuration file
  string paramfile;		     // The parameters files									

  string plot_dir;                   // Name of the folder comtaining the output plots
  string plot_name;		     // Name of the plot to draw
  string folder;		     // Folder containing input files
  bool is_2012;			     // false = 7 TeV, true = 8 TeV
  bool no_plot;			     // Don't actually generate image files if true
  bool verbose;			     // Verbose output, useful for diagnostic purposes

  bool saveroot;                     // Save plots as a root file

  // Plotting options
  string x_axis_label;		     // Label for the X-axis
  string x_axis_bin_labels;          // of the form "a:b:c"
  unsigned rebin;		     // Rebin factor
  bool custom_x_axis_range;	     // Can optionally specify an x-axis range
  double x_axis_min;		     // If custom_x_axis_range is true, use this as min
  double x_axis_max;		     // If custom_x_axis_range is true, use this as max
  bool custom_y_axis_range;	     // Can optionally specify a y-axis range
  double y_axis_min;		     // If custom_y_axis_range is true, use this as min
  double y_axis_max;		     // If custom_y_axis_range is true, use this as max
  double y_ratio_min;		     // Use this as min for ratio plot
  double y_ratio_max;		     // Use this as max for ratio plot
  double extra_pad;		     // Expand the y-axis by an additional factor
  bool blind;			     // Blind some region of the data
  double x_blind_min;		     // If blind is true, use this as min x for blinding.
                                     // If min>max blind everything
  double x_blind_max;		     // If blind is true, use this as max x for blinding
  bool log_y;			     // Draw plot in log scale
  bool norm_bins;		     // Normalise using bin width
  bool signal_no_stack;		     // Don't stack the signal contributions on the backgrounds
  bool draw_ratio;                   // Draw a ratio box
  bool plot_qcd;                     // Include QCD in plots
  bool plot_data_qcd;                // Include QCD from data in plots
  bool plot_wjets_comp;              // Separate Wjets components in plots
  bool use_embedded_data;            // use embedded data instead of data
  bool use_rereco_data;            // use rereco data instead of prompt data

  bool addOverflows;                 //add two bins to histo for underflows and overflows.
  bool dolumixsweight;               //dolumixsweights in controlplots

  // Options to manually shift backgrounds and draw uncertainty bands
  bool shift_backgrounds = false;
  double qcd_shift = 1.0;
  double top_shift = 1.0;
  double ztt_shift = 1.0;
  double zll_shift = 1.0;
  double w_shift = 1.0;
  double vv_shift = 1.0;
  bool draw_band_on_stack = false;
  double band_size_fractional = 0.0;

  po::options_description preconfig("Pre-Configuration");
  preconfig.add_options()("cfg", po::value<std::string>(&cfg)->required());
  po::variables_map vm;
  po::store(po::command_line_parser(argc, argv).
            options(preconfig).allow_unregistered().run(), vm);
  po::notify(vm);
  po::options_description config("Configuration");
  config.add_options()
    ("paramfile",           po::value<string>(&paramfile)->required())
    ("plot_dir",            po::value<string>(&plot_dir)->default_value("PLOTS"))
    ("plot_name",           po::value<string>(&plot_name)->required())
    ("folder",              po::value<string>(&folder)->required())
    ("is_2012",             po::value<bool>(&is_2012)->required())
    ("draw_ratio",          po::value<bool>(&draw_ratio)->default_value(false))
    ("no_plot",             po::value<bool>(&no_plot)->default_value(false))
    ("verbose",             po::value<bool>(&verbose)->default_value(false))
    ("saveroot",            po::value<bool>(&saveroot)->default_value(false))
    ("x_axis_label",        po::value<string>(&x_axis_label)->required())
    ("x_axis_bin_labels",   po::value<string>(&x_axis_bin_labels)->default_value(""))
    ("rebin",               po::value<unsigned>(&rebin)->default_value(1))
    ("custom_x_axis_range", po::value<bool>(&custom_x_axis_range)->default_value(false))
    ("x_axis_min",          po::value<double>(&x_axis_min)->default_value(0))
    ("x_axis_max",          po::value<double>(&x_axis_max)->default_value(0))
    ("custom_y_axis_range", po::value<bool>(&custom_y_axis_range)->default_value(false))
    ("y_axis_min",          po::value<double>(&y_axis_min)->default_value(0))
    ("y_axis_max",          po::value<double>(&y_axis_max)->default_value(0))
    ("y_ratio_min",         po::value<double>(&y_ratio_min)->default_value(0))
    ("y_ratio_max",         po::value<double>(&y_ratio_max)->default_value(2.))
    ("extra_pad",           po::value<double>(&extra_pad)->default_value(1.0))
    ("blind",               po::value<bool>(&blind)->default_value(false))
    ("x_blind_min",         po::value<double>(&x_blind_min)->default_value(0))
    ("x_blind_max",         po::value<double>(&x_blind_max)->default_value(0))
    ("log_y",               po::value<bool>(&log_y)->default_value(false))
    ("norm_bins",           po::value<bool>(&norm_bins)->default_value(false))
    ("signal_no_stack",     po::value<bool>(&signal_no_stack)->default_value(false))
    ("plot_qcd",            po::value<bool>(&plot_qcd)->default_value(true))
    ("plot_data_qcd",       po::value<bool>(&plot_data_qcd)->default_value(false))
    ("plot_wjets_comp",     po::value<bool>(&plot_wjets_comp)->default_value(true))
    ("use_embedded_data",   po::value<bool>(&use_embedded_data)->default_value(false))
    ("use_rereco_data",   po::value<bool>(&use_rereco_data)->default_value(false))
    ("addOverflows",        po::value<bool>(&addOverflows)->default_value(false))
    ("dolumixsweight",      po::value<bool>(&dolumixsweight)->default_value(false))
    ("shift_backgrounds",   po::value<bool>(&shift_backgrounds)->default_value(false))
    ("draw_band_on_stack",  po::value<bool>(&draw_band_on_stack)->default_value(false))
    ("qcd_shift",           po::value<double>(&qcd_shift)->default_value(1.0))
    ("top_shift",           po::value<double>(&top_shift)->default_value(1.0))
    ("ztt_shift",           po::value<double>(&ztt_shift)->default_value(1.0))
    ("zll_shift",           po::value<double>(&zll_shift)->default_value(1.0))
    ("w_shift",             po::value<double>(&w_shift)->default_value(1.0))
    ("vv_shift",            po::value<double>(&vv_shift)->default_value(1.0))
    ("band_size_fractional", po::value<double>(&band_size_fractional)->default_value(0.0));
  po::store(po::command_line_parser(argc, argv).options(config).allow_unregistered().run(), vm);
  po::store(po::parse_config_file<char>(cfg.c_str(), config), vm);
  po::notify(vm);

  std::cout << "**** ControlPlots *****" << std::endl;

  string param_fmt = "%-25s %-40s\n";
  std::cout << boost::format(param_fmt) % "paramfile" 	% paramfile;
  std::cout << boost::format(param_fmt) % "plot_dir" 	% plot_dir;
  std::cout << boost::format(param_fmt) % "plot_name" 	% plot_name;
  std::cout << boost::format(param_fmt) % "folder" 		% folder;
  std::cout << boost::format(param_fmt) % "is_2012" 		% is_2012;
  std::cout << boost::format(param_fmt) % "no_plot" 		% no_plot;
  std::cout << boost::format(param_fmt) % "plot_qcd" 		% plot_qcd;
  std::cout << boost::format(param_fmt) % "plot_data_qcd" 	% plot_data_qcd;
  std::cout << boost::format(param_fmt) % "plot_wjets_comp" 	% plot_wjets_comp;
  std::cout << boost::format(param_fmt) % "use_embedded_data" 	% use_embedded_data ;


  // Parse the parameter file
  SimpleParamParser parser;
  std::cout << "** Parsing parameter file... **" << paramfile << std::endl;
  parser.ParseFile(paramfile);


  bool draw_signal = parser.GetParam<bool>("DRAW_SIGNAL");
  unsigned draw_signal_factor = parser.GetParam<unsigned>("DRAW_SIGNAL_FACTOR");
  string year_label = parser.GetParam<string>("YEAR_LABEL");
  double data_lumi = parser.GetParam<double>("LUMI_DATA");
  double data_embedded_lumi = 9964.0;//parser.GetParam<double>("LUMI_DATA_EMBEDDED");
  string lumi_data_label = parser.GetParam<string>("LUMI_DATA_LABEL");


  //List of input files
  vector<string> files;
  //files.push_back("Data_MET-2012A-13Jul2012-v1_0_795305fb");
  if (use_embedded_data) {
    files.push_back("DataEmbedded_METembedded-all");
  }
  else if(use_rereco_data){
    files.push_back("MET_MET-2012A-22Jan2013-v1");
    files.push_back("MET_MET-2012B-22Jan2013-v1");
    files.push_back("MET_MET-2012C-22Jan2013-v1");
    files.push_back("MET_METParked-2012D-22Jan2013-v1");
  }
  else{
    files.push_back("Data_MET-2012A-13Jul2012-v1");
    files.push_back("Data_MET-2012A-06Aug2012-v1");
    files.push_back("Data_MET-2012B-13Jul2012-v1");
    files.push_back("Data_MET-2012C-24Aug2012-v1");
    files.push_back("Data_MET-2012C-11Dec2012-v1"); 
    files.push_back("Data_MET-2012C-PromptReco-v2");
    files.push_back("Data_MET-2012D-PromptReco-v1");
    //files.push_back("Data_MET-2012D-PromptReco-v1_7_269fb");
    //files.push_back("Data_MET-2012D-PromptReco-v1_7_193fb");
    if (plot_data_qcd) files.push_back("Data_MET0to120");
  }
  files.push_back("MC_QCD-Pt-30to50-pythia6");
  files.push_back("MC_QCD-Pt-50to80-pythia6");
  files.push_back("MC_QCD-Pt-80to120-pythia6");
  files.push_back("MC_QCD-Pt-120to170-pythia6");
  files.push_back("MC_QCD-Pt-170to300-pythia6");
  files.push_back("MC_QCD-Pt-300to470-pythia6");
  files.push_back("MC_QCD-Pt-470to600-pythia6");
  files.push_back("MC_QCD-Pt-600to800-pythia6");
  files.push_back("MC_QCD-Pt-800to1000-pythia6");
  files.push_back("MC_QCD-Pt-1000to1400-pythia6");
  files.push_back("MC_QCD-Pt-1400to1800-pythia6");
  files.push_back("MC_QCD-Pt-1800-pythia6");
  files.push_back("MC_TTJets");
  //powheg samples
  files.push_back("MC_TT-v1");
  files.push_back("MC_TT-v2");
  //
  files.push_back("MC_T-tW");
  files.push_back("MC_Tbar-tW");
  files.push_back("MC_SingleT-s-powheg-tauola");
  files.push_back("MC_SingleTBar-s-powheg-tauola");
  files.push_back("MC_SingleT-t-powheg-tauola");
  files.push_back("MC_SingleTBar-t-powheg-tauola");
  files.push_back("MC_WW-pythia6-tauola");
  files.push_back("MC_WZ-pythia6-tauola");
  files.push_back("MC_ZZ-pythia6-tauola");
  files.push_back("MC_W1JetsToLNu_enu");
  files.push_back("MC_W2JetsToLNu_enu");
  files.push_back("MC_W3JetsToLNu_enu");
  files.push_back("MC_W4JetsToLNu_enu");
  files.push_back("MC_WJetsToLNu-v1_enu");
  files.push_back("MC_WJetsToLNu-v2_enu");
  files.push_back("MC_W1JetsToLNu_munu");
  files.push_back("MC_W2JetsToLNu_munu");
  files.push_back("MC_W3JetsToLNu_munu");
  files.push_back("MC_W4JetsToLNu_munu");
  files.push_back("MC_WJetsToLNu-v1_munu");
  files.push_back("MC_WJetsToLNu-v2_munu");
  files.push_back("MC_W1JetsToLNu_taunu");
  files.push_back("MC_W2JetsToLNu_taunu");
  files.push_back("MC_W3JetsToLNu_taunu");
  files.push_back("MC_W4JetsToLNu_taunu");
  files.push_back("MC_WJetsToLNu-v1_taunu");
  files.push_back("MC_WJetsToLNu-v2_taunu");
  files.push_back("MC_DYJetsToLL");
  files.push_back("MC_DY1JetsToLL");
  files.push_back("MC_DY2JetsToLL");
  files.push_back("MC_DY3JetsToLL");
  files.push_back("MC_DY4JetsToLL");
  files.push_back("MC_DYJetsToLL_PtZ-100-madgraph");//High pt sample
  files.push_back("MC_ZJetsToNuNu_100_HT_200");
  files.push_back("MC_ZJetsToNuNu_200_HT_400");
  files.push_back("MC_ZJetsToNuNu_400_HT_inf");
  files.push_back("MC_ZJetsToNuNu_50_HT_100");
  files.push_back("MC_GJets-HT-200To400-madgraph");
  files.push_back("MC_GJets-HT-400ToInf-madgraph");
  files.push_back("MC_VBF_HToZZTo4Nu_M-120");
  files.push_back("MC_EWK-Z2j");
  //files.push_back("MC_EWK-Z2jiglep");
  files.push_back("MC_EWK-W2jminus_enu");
  files.push_back("MC_EWK-W2jplus_enu");
  files.push_back("MC_EWK-W2jminus_munu");
  files.push_back("MC_EWK-W2jplus_munu");
  files.push_back("MC_EWK-W2jminus_taunu");
  files.push_back("MC_EWK-W2jplus_taunu");
  files.push_back("MC_WGamma");
  


  //build a list of selections
  vector<string> selections;
  selections.push_back("HLTMetClean");
  selections.push_back("LeptonVeto");
  selections.push_back("WSelection");
  selections.push_back("ZSelection");
  selections.push_back("IgnoreLeptons");
  selections.push_back("JetPair");
  selections.push_back("AN");
  selections.push_back("DEta");
  selections.push_back("MET");
  selections.push_back("LooseMjj");
  selections.push_back("TightMjj");
  selections.push_back("DPhiSIGNAL_noCJV");
  selections.push_back("DPhiQCD_noCJV");
  selections.push_back("CJVfail");
  selections.push_back("DPhiSIGNAL_CJVfail");
  selections.push_back("DPhiQCD_CJVfail");
  selections.push_back("CJVpass");
  selections.push_back("DPhiSIGNAL_CJVpass");
  selections.push_back("DPhiQCD_CJVpass");
  selections.push_back("MtZepp");
  selections.push_back("VBF");
  selections.push_back("GenNoFilters");
  selections.push_back("GenZMassFiltered");
  
  vector<string> latex;
  latex.push_back("HLTMetClean");
  latex.push_back("LeptonVeto");
  latex.push_back("WSelection");
  latex.push_back("ZSelection");
  latex.push_back("IgnoreLeptons");
  latex.push_back("JetPair");
  latex.push_back("AN");
  latex.push_back("DEta");
  latex.push_back("MET");
  latex.push_back("LooseMjj");
  latex.push_back("TightMjj");
  latex.push_back("SIG noCJV");
  latex.push_back("QCD noCJV");
  latex.push_back("CJVfail");
  latex.push_back("SIG CJVfail");
  latex.push_back("QCD CJVfail");
  latex.push_back("CJVpass");
  latex.push_back("SIG CJVpass");
  latex.push_back("QCD CJVpass");
  latex.push_back("MtZepp");
  latex.push_back("VBF");
  latex.push_back("GenNoFilters");
  latex.push_back("GenZMassFiltered");
  

  vector<string> selectionsdir = selections;

//   if(dopusyst){//If do PU uncertainty get n_jets hist from systematics directory
//     if((plot_name.find("n_jets") != plot_name.npos && plot_name.find("ingap") == plot_name.npos) || (plot_name.find("n_vtx") != plot_name.npos) || (plot_name.find("dphijj") != plot_name.npos)){
//       for(int m=0;unsigned(m)<selections.size();m++){
// 	selectionsdir[m]+="/systematics";
//       }
//       if(puupordown){//Do PUUP error
// 	plot_name+="_puUp";
//       }
//       else if(!puupordown){//Do PUDOWN error
// 	plot_name+="_puDown";
//       }
//     }
//   }

  if (plot_name.find("Up") != plot_name.npos  || plot_name.find("Down") != plot_name.npos) {
    for(int m=0;unsigned(m)<selections.size();m++){
      selectionsdir[m]+="/systematics";
    }
  }
  else if (plot_name.find("_pu") != plot_name.npos || plot_name.find("_noW") != plot_name.npos) {
    for(int m=0;unsigned(m)<selections.size();m++){
      selectionsdir[m]+="/weights";
    }
  }

  std::string lSuffix = "";
  if (plot_name.find("n_jets")!= plot_name.npos) lSuffix = plot_name.substr(plot_name.find("n_jets")+6,plot_name.npos);

  std::map<std::string, TFile *> tfiles;
  for (unsigned i = 0; i < files.size(); ++i) {
    std::string filename = (files[i]+".root");
    TFile * tmp = new TFile((folder+"/"+filename).c_str());
    if (!tmp) {
      std::cerr << "Warning, file " << filename << " could not be opened." << std::endl;
    } else {
      tfiles[files[i]] = tmp;      
    }
  }

  // Get Plots and Scale
  map<string, ic::TH1PlotElement> plots;
  bool skip[selectionsdir.size()];
  for (unsigned k = 0; k < selectionsdir.size(); ++k) {
    skip[k] = false;
  }
  for (unsigned i = 0; i < files.size(); ++i) {
    for (unsigned k = 0; k < selectionsdir.size(); ++k) {
      
      string f = files[i];
      string s = selectionsdir[k];

      std::string nm = Token(f, s);
      if (plots.count(nm)) continue;
      if (verbose) std::cout << "Adding plot: " << nm << std::endl;

      if (tfiles[f]->GetDirectory(("/"+s).c_str())) {
	plots[nm] = ic::TH1PlotElement(nm, tfiles[f], "/"+s, plot_name);
	if (!plots[nm].hist_ptr()) {
	  if (!skip[k]) std::cerr << " Histogram " << plot_name << " not found for " << nm << ". Skipping..." << std::endl;
	  skip[k] = true;
	  continue;
	}
	if (plots[nm].hist_ptr()->GetSumw2N() == 0) {
	  cout << "Warning: Plot " << nm << " does not have a weights structure" << endl;
	  plots[nm].hist_ptr()->Sumw2();
	}

	if(dolumixsweight){
	  string lookup = f;
	  double sample_events = parser.GetParam<double>("EVT_"+lookup);
	  double sample_xs = parser.GetParam<double>("XS_"+lookup);
	  double sample_lumi = sample_events / sample_xs;
	  double sample_scale = 1;
	  if (sample_xs > 0) sample_scale = data_lumi / sample_lumi;
	  if (use_embedded_data) sample_scale = data_lumi / data_embedded_lumi;
	  //std::cout<<"Doing lumixsweight: "<<sample_scale<<std::endl;
	  plots[nm].hist_ptr()->Scale(sample_scale);
	}

	plots[nm].hist_ptr()->Rebin(rebin);
	
	if (addOverflows) plots[nm].AddOverflows();

	//if (k==0) std::cout << f << " " << sample_scale << std::endl;

      }
      else {
	if (!skip[k]) std::cerr << " Directory " << s << " not found. Skipping..." << std::endl;
	skip[k] = true;
      }

    }//loop on selections
  }//loop on files

  
  std::ofstream lTexOutput;
  std::ofstream lDatOutput[selections.size()];

  bool lFillSummaryTable = plot_name.find("n_jets") != plot_name.npos && plot_name.find("ingap") == plot_name.npos;

  //output a table with number of events selected
  if (lFillSummaryTable){
     //lTexOutput << "\\begin{tabular}{|l|c|c|c|c|c||c|c||c|}" << std::endl
     lTexOutput.open(plot_dir+"/SummaryTable"+lSuffix+".txt",std::ios_base::out);
     lTexOutput << "\\begin{tabular}{|l|p{0.07\\textwidth}|p{0.07\\textwidth}|p{0.07\\textwidth}|p{0.07\\textwidth}|p{0.07\\textwidth}|p{0.07\\textwidth}||p{0.07\\textwidth}|c||p{0.07\\textwidth}|}" << std::endl
		<<"\\hline" << std::endl
		<< "Step & QCD & $\\gamma$+jets & Top & W+jets & Z+jets & VV & SumMC & Data & Signal 120 \\\\"
		<< std::endl
		<< "\\hline" << std::endl;
   }


  for (unsigned k = 0; k < selections.size(); ++k) {
    if (skip[k]) continue;
    if (lFillSummaryTable) lDatOutput[k].open(plot_dir+"/SummaryTable_"+selections[k]+lSuffix+".dat",std::ios_base::out);
    ic::Plot plot;
    plot.output_filename = plot_dir+"/"+plot_name + "_" + selections[k] + ".pdf";
    if (log_y) plot.output_filename = plot_dir+"/"+plot_name + "_" + selections[k]  + "_log.pdf";
  
    plot.x_bin_labels_ = x_axis_bin_labels;

    ic::TH1PlotElement data_hist = ic::TH1PlotElement("Data");//+selections[k]);
    ic::TH1PlotElement data_qcd_hist = ic::TH1PlotElement("DataQCD");
    ic::TH1PlotElement signal_hist = ic::TH1PlotElement("Signal");
    ic::TH1PlotElement qcd_hist = ic::TH1PlotElement("QCD");
    ic::TH1PlotElement top_hist = ic::TH1PlotElement("Top");
    ic::TH1PlotElement ttbar_hist = ic::TH1PlotElement("TTBar");
    ic::TH1PlotElement ttpowheg_hist = ic::TH1PlotElement("TTBarPowheg");
    ic::TH1PlotElement singletop_hist = ic::TH1PlotElement("SingleTop");
    ic::TH1PlotElement tW_hist = ic::TH1PlotElement("tW");
    ic::TH1PlotElement WJets_hist = ic::TH1PlotElement("WJets");
    ic::TH1PlotElement WJets_enu_hist = ic::TH1PlotElement("WJets_enu");
    ic::TH1PlotElement WJets_munu_hist = ic::TH1PlotElement("WJets_munu");
    ic::TH1PlotElement WJets_taunu_hist = ic::TH1PlotElement("WJets_taunu");
    ic::TH1PlotElement ZJetsToLL_hist = ic::TH1PlotElement("ZJetsToLL");
    ic::TH1PlotElement ZJetsToNuNu_hist = ic::TH1PlotElement("ZJetsToNuNu");
    ic::TH1PlotElement VBFZ_hist = ic::TH1PlotElement("VBFZ");
    ic::TH1PlotElement VBFZnunu_hist = ic::TH1PlotElement("VBFZnunu");
    ic::TH1PlotElement VBFW_hist = ic::TH1PlotElement("VBFW");
    ic::TH1PlotElement VBFW_enu_hist = ic::TH1PlotElement("VBFW_enu");
    ic::TH1PlotElement VBFW_munu_hist = ic::TH1PlotElement("VBFW_munu");
    ic::TH1PlotElement VBFW_taunu_hist = ic::TH1PlotElement("VBFW_taunu");
    ic::TH1PlotElement GJets_hist = ic::TH1PlotElement("GJets");
    ic::TH1PlotElement VV_hist = ic::TH1PlotElement("Dibosons");

 
    for (unsigned i = 0; i < files.size(); ++i) {
      string f = files[i];
      string s = selectionsdir[k];

      std::string nm = Token(f, s);

      if (f.find("0to120") == f.npos){
	if(use_rereco_data){
	  SumHistograms(f,plots[nm],"MET",data_hist);
	}
	else{
	  SumHistograms(f,plots[nm],"Data",data_hist);
	}
      }
      SumHistograms(f,plots[nm],"Data_MET0to120",data_qcd_hist);
      SumHistograms(f,plots[nm],"VBF_H",signal_hist);
      //qcd
      SumHistograms(f,plots[nm],"MC_QCD",qcd_hist);
      SumHistograms(f,plots[nm],"GJets",GJets_hist);
      //top
      SumHistograms(f,plots[nm],"MC_TTJets",top_hist);
      SumHistograms(f,plots[nm],"MC_T-tW",top_hist);
      SumHistograms(f,plots[nm],"MC_Tbar-tW",top_hist);
      SumHistograms(f,plots[nm],"MC_SingleT",top_hist);
      //top-separated
      SumHistograms(f,plots[nm],"MC_TTJets",ttbar_hist);
      SumHistograms(f,plots[nm],"MC_TT-v",ttpowheg_hist);
      SumHistograms(f,plots[nm],"MC_SingleT",singletop_hist);
      SumHistograms(f,plots[nm],"tW",tW_hist);
      //ewk V+2j
      SumHistograms(f,plots[nm],"EWK-Z2j",VBFZ_hist);
      SumHistograms(f,plots[nm],"EWK-Z2jiglep",VBFZnunu_hist);
      SumHistograms(f,plots[nm],"EWK-W2j",VBFW_hist);
      SumHistograms(f,plots[nm],"EWK-W2jminus_enu",VBFW_enu_hist);
      SumHistograms(f,plots[nm],"EWK-W2jplus_enu",VBFW_enu_hist);
      SumHistograms(f,plots[nm],"EWK-W2jminus_munu",VBFW_munu_hist);
      SumHistograms(f,plots[nm],"EWK-W2jplus_munu",VBFW_munu_hist);
      SumHistograms(f,plots[nm],"EWK-W2jminus_taunu",VBFW_taunu_hist);
      SumHistograms(f,plots[nm],"EWK-W2jplus_taunu",VBFW_taunu_hist);
      //Z and W
      SumHistograms(f,plots[nm],"JetsToLL",ZJetsToLL_hist);
      SumHistograms(f,plots[nm],"JetsToLNu",WJets_hist);
      //SumHistograms(f,plots[nm],"EWK-W2j",WJets_hist);
      if (f.find("JetsToLNu") != f.npos){
	SumHistograms(f,plots[nm],"enu",WJets_enu_hist);
	SumHistograms(f,plots[nm],"munu",WJets_munu_hist);
	SumHistograms(f,plots[nm],"taunu",WJets_taunu_hist);
      }
      SumHistograms(f,plots[nm],"ZJetsToNuNu",ZJetsToNuNu_hist);
      //VV
      SumHistograms(f,plots[nm],"MC_WW",VV_hist);
      SumHistograms(f,plots[nm],"MC_WZ",VV_hist);
      SumHistograms(f,plots[nm],"MC_ZZ",VV_hist);
      SumHistograms(f,plots[nm],"MC_WGamma",VV_hist);
      
      
 
    }//loop on files

    (signal_hist.hist_ptr())->Scale(draw_signal_factor);

    // std::cout << " DEBUG: " << WJets_hist.hist_ptr()->Integral() << " " << WJets_enu_hist.hist_ptr()->Integral()+WJets_munu_hist.hist_ptr()->Integral()+WJets_taunu_hist .hist_ptr()->Integral()<< std::endl;


//     SetSignalStyle(signal_hist,2);
//     SetDataStyle(data_hist);
//     SetBkgStyle(qcd_hist,8);
//     SetBkgStyle(top_hist,2);
//     SetBkgStyle(WJets_hist,4);
//     SetBkgStyle(ZJetsToLL_hist,7);
//     SetBkgStyle(ZJetsToNuNu_hist,7);
//     SetBkgStyle(VBFZ_hist,7);
//     SetBkgStyle(GJets_hist,8);
//     SetBkgStyle(VV_hist,5);

 
    SetSignalStyle(signal_hist,1);
    SetDataStyle(data_hist);
    SetBkgStyle(qcd_hist,7);
    SetBkgStyle(data_qcd_hist,7);
    SetBkgStyle(top_hist,5);
    SetBkgStyle(WJets_hist,6);
    SetBkgStyle(VBFW_hist,kViolet-6);
    SetBkgStyle(WJets_enu_hist,2);
    SetBkgStyle(WJets_munu_hist,kOrange);
    SetBkgStyle(WJets_taunu_hist,6);
    SetBkgStyle(VBFW_enu_hist,kRed);
    SetBkgStyle(VBFW_munu_hist,kOrange+8);
    SetBkgStyle(VBFW_taunu_hist,kViolet-6);
    SetBkgStyle(ZJetsToLL_hist,3);
    SetBkgStyle(ZJetsToNuNu_hist,3);
    SetBkgStyle(VBFZ_hist,3);
    SetBkgStyle(GJets_hist,7);
    SetBkgStyle(VV_hist,4);

    data_hist.set_legend_text("Data");
    if (draw_signal) {
      if (draw_signal_factor>1) 
	signal_hist.set_legend_text("VBF m_{H}=120 GeV #times"+boost::lexical_cast<std::string>(draw_signal_factor));
      else signal_hist.set_legend_text("VBF m_{H}=120 GeV");
    }
    if (plot_qcd || plot_data_qcd) {
      qcd_hist.set_legend_text("QCD,#gamma+jets");
      data_qcd_hist.set_legend_text("QCD");
      //GJets_hist.set_legend_text("#gamma + jets");
    }
    ZJetsToNuNu_hist.set_legend_text("Z+jets,EWK Z");
    if (!plot_wjets_comp) {
      WJets_hist.set_legend_text("QCD W+jets");
      VBFW_hist.set_legend_text("EWK W+jets");
    }
    else {
      WJets_taunu_hist.set_legend_text("W#rightarrow#tau#nu+jets");
      WJets_munu_hist.set_legend_text("W#rightarrow#mu#nu+jets");
      WJets_enu_hist.set_legend_text("W#rightarrow e#nu+jets");
    }
    top_hist.set_legend_text("t#bar{t},t,tW");
    VV_hist.set_legend_text("Dibosons");
    //VBFZ_hist.set_legend_text("VBF Z+2j");
    //ZJetsToLL_hist.set_legend_text("Z#rightarrow ll + jets");
    //ZJetsToNuNu_hist.set_legend_text("Z#rightarrow #nu#nu + jets");


    bool lBlind = blind && (
			    selections[k].find("DPhiSIGNAL") != selections[k].npos ||
			    (selections[k].find("CJV") != selections[k].npos && 
			     plot_name.find("dphijj") != plot_name.npos )
			    );

    if (lBlind) {
      for (int j = 0; j < data_hist.hist_ptr()->GetNbinsX(); ++j) {
	if (x_blind_min < x_blind_max){
	  double low_edge = data_hist.hist_ptr()->GetBinLowEdge(j+1);
	  double high_edge = data_hist.hist_ptr()->GetBinWidth(j+1)+data_hist.hist_ptr()->GetBinLowEdge(j+1);
	  if ((low_edge > x_blind_min && low_edge < x_blind_max) || (high_edge > x_blind_min && high_edge < x_blind_max)) {
	    data_hist.hist_ptr()->SetBinContent(j+1,0);
	    data_hist.hist_ptr()->SetBinError(j+1,0);
	  }
	}
	else {
	  data_hist.hist_ptr()->SetBinContent(j+1,0);
	  data_hist.hist_ptr()->SetBinError(j+1,0);
	}
      }
    }

    plot.AddTH1PlotElement(VV_hist);
    plot.AddTH1PlotElement(top_hist);
    if (!plot_wjets_comp) {
      plot.AddTH1PlotElement(WJets_hist);
      plot.AddTH1PlotElement(VBFW_hist);
    }
    else {
      plot.AddTH1PlotElement(WJets_enu_hist);
      plot.AddTH1PlotElement(VBFW_enu_hist);
      plot.AddTH1PlotElement(WJets_munu_hist);
      plot.AddTH1PlotElement(VBFW_munu_hist);
      plot.AddTH1PlotElement(WJets_taunu_hist);
      plot.AddTH1PlotElement(VBFW_taunu_hist);
    }
    plot.AddTH1PlotElement(ZJetsToNuNu_hist);
    plot.AddTH1PlotElement(ZJetsToLL_hist);
    
    //moving signal before data screws up everything !!
    plot.AddTH1PlotElement(VBFZ_hist);
    if (plot_qcd){
      plot.AddTH1PlotElement(GJets_hist);
      plot.AddTH1PlotElement(qcd_hist);
    }
    else if (plot_data_qcd){
      plot.AddTH1PlotElement(data_qcd_hist);
    }
   
    plot.AddTH1PlotElement(data_hist);
 
    if (draw_signal) {
      plot.AddTH1PlotElement(signal_hist);
    }

    plot.extra_pad = extra_pad;
    plot.custom_x_axis_range = custom_x_axis_range;
    if (custom_x_axis_range){
      plot.x_axis_min = x_axis_min;
      plot.x_axis_max = x_axis_max;
    }
    if (log_y) plot.y_axis_log = true;
    if (y_axis_min > 0) plot.y_axis_min = y_axis_min;
    if (custom_y_axis_range){
      plot.y_axis_min = y_axis_min;
      plot.y_axis_max = y_axis_max;
    }
    plot.legend_height = 0.045;
    plot.x_axis_title = x_axis_label;
    plot.y_axis_title = "Events";
    plot.title_left = "CMS Preliminary " + year_label +", #sqrt{s} = " + (is_2012 ? "8":"7") +" TeV, "+ lumi_data_label;
    plot.title_right = selections[k];
    

    //ratio plot

    plot.draw_ratio_hist = draw_ratio;
    plot.draw_signif = false;
    
    string background_list = "VV+Top+WJets_enu+WJets_munu+WJets_taunu+VBFW_enu+VBFW_munu+VBFW_taunu+ZJetsToLL+VBFZ+ZJetsToNuNu";
    if (!plot_wjets_comp) background_list = "VV+Top+WJets+VBFW+ZJetsToLL+VBFZ+ZJetsToNuNu";
    if (plot_qcd) background_list += "+GJets+QCD";
    else if (plot_data_qcd) background_list += "+DataQCD";
    ic::RatioPlotElement ratio("DataOverMC","Data",background_list);
    
    plot.band_size_fractional_ = band_size_fractional;
    plot.draw_band_on_stack_ = draw_band_on_stack;
    plot.samples_for_band_ = background_list;
    
    SetStyle(ratio,1);
    ratio.set_multi_mode(true);
    //plot.ratio_y_axis_title = "#frac{(N_{obs}-N_{exp})}{#sqrt{N_{obs}}}";
    plot.ratio_y_axis_title = "Data/MC";
    plot.AddRatioPlotElement(ratio);
    plot.custom_ratio_y_axis_range = true;
    plot.ratio_y_axis_min = y_ratio_min;
    plot.ratio_y_axis_max = y_ratio_max;
    
    
    if (!no_plot){
      std::vector<string> types={"pdf","png"};
      if(saveroot){
	types.push_back("root");
	types.push_back("C");
      }
      plot.GeneratePlot(types);
    }

    if (lFillSummaryTable) {

      Utilities n_qcd = Utilities(Integral(qcd_hist.hist_ptr()),Error(qcd_hist.hist_ptr()));
      Utilities n_VV = Utilities(Integral(VV_hist.hist_ptr()),Error(VV_hist.hist_ptr()));
      Utilities n_gjets = Utilities(Integral(GJets_hist.hist_ptr()),Error(GJets_hist.hist_ptr()));
      //top
      Utilities n_top = Utilities(Integral(top_hist.hist_ptr()),Error(top_hist.hist_ptr()));
      Utilities n_ttbar = Utilities(Integral(ttbar_hist.hist_ptr()),Error(ttbar_hist.hist_ptr()));
      Utilities n_ttpowheg = Utilities(Integral(ttpowheg_hist.hist_ptr()),Error(ttpowheg_hist.hist_ptr()));
      Utilities n_singletop = Utilities(Integral(singletop_hist.hist_ptr()),Error(singletop_hist.hist_ptr()));
      Utilities n_tW = Utilities(Integral(tW_hist.hist_ptr()),Error(tW_hist.hist_ptr()));
      //W+jets
      Utilities n_WJets = Utilities(Integral(WJets_hist.hist_ptr()),Error(WJets_hist.hist_ptr()));
      Utilities n_EWKW = Utilities(Integral(VBFW_hist.hist_ptr()),Error(VBFW_hist.hist_ptr()));
      Utilities n_WJets_enu = Utilities(Integral(WJets_enu_hist.hist_ptr()),Error(WJets_enu_hist.hist_ptr()));
      Utilities n_WJets_munu = Utilities(Integral(WJets_munu_hist.hist_ptr()),Error(WJets_munu_hist.hist_ptr()));
      Utilities n_WJets_taunu = Utilities(Integral(WJets_taunu_hist.hist_ptr()),Error(WJets_taunu_hist.hist_ptr()));
      Utilities n_EWKW_enu = Utilities(Integral(VBFW_enu_hist.hist_ptr()),Error(VBFW_enu_hist.hist_ptr()));
      Utilities n_EWKW_munu = Utilities(Integral(VBFW_munu_hist.hist_ptr()),Error(VBFW_munu_hist.hist_ptr()));
      Utilities n_EWKW_taunu = Utilities(Integral(VBFW_taunu_hist.hist_ptr()),Error(VBFW_taunu_hist.hist_ptr()));
      Utilities n_WJets_tot = Utilities(n_WJets.rawNumber()+n_EWKW.rawNumber(),
					sumSqErr(n_WJets.rawError(),n_EWKW.rawError()));
      //Z+jets
      Utilities n_ZJets_nunu = Utilities(Integral(ZJetsToNuNu_hist.hist_ptr()),Error(ZJetsToNuNu_hist.hist_ptr()));
      Utilities n_ZJets_ll   = Utilities(Integral(ZJetsToLL_hist.hist_ptr()),Error(ZJetsToLL_hist.hist_ptr()));
      Utilities n_ZJets_vbf  = Utilities(Integral(VBFZ_hist.hist_ptr()),Error(VBFZ_hist.hist_ptr()));
      Utilities n_ZJets_vbf_nunu  = Utilities(Integral(VBFZnunu_hist.hist_ptr()),Error(VBFZnunu_hist.hist_ptr()));
      Utilities n_ZJets = Utilities(n_ZJets_ll.rawNumber()+n_ZJets_nunu.rawNumber()+n_ZJets_vbf.rawNumber(),
				    sumSqErr(n_ZJets_ll.rawError(),n_ZJets_nunu.rawError(),n_ZJets_vbf.rawError()));
      //others  
      Utilities n_others = Utilities(n_gjets.rawNumber()+n_VV.rawNumber(),
				     sumSqErr(n_gjets.rawError(),n_VV.rawError()));
      Utilities n_data = Utilities(Integral(data_hist.hist_ptr()),Error(data_hist.hist_ptr()));
      Utilities n_data_qcd = Utilities(Integral(data_qcd_hist.hist_ptr()),Error(data_qcd_hist.hist_ptr()));
      Utilities n_signal = Utilities(Integral(signal_hist.hist_ptr()),Error(signal_hist.hist_ptr()));
      Utilities n_Tot = Utilities(n_qcd.rawNumber()+n_top.rawNumber()+
				  n_WJets.rawNumber()+n_EWKW.rawNumber()+
				  n_ZJets.rawNumber()+n_others.rawNumber(),
       				  sumSqErr(n_qcd.rawError(),n_top.rawError(),
					   n_WJets.rawError(),n_EWKW.rawError(),
					   n_ZJets.rawError(),n_others.rawError()));
       
      //std::cout.precision(2);
      //std::cout << std::scientific;

      lTexOutput << latex[k] << " & " 
		 << n_qcd.roundedResult() << " & " 
		 << n_gjets.roundedResult() << " & "
		 << n_top.roundedResult() << " & "
		 << n_WJets_tot.roundedResult() << " & "
	//<< n_WJets_enu.roundedResult() << " & "
	//<< n_WJets_munu.roundedResult() << " & "
	//<< n_WJets_taunu.roundedResult() << " & "
		 << n_ZJets.roundedResult() << " & "
		 << n_VV.roundedResult() << " & "
		 << n_Tot.roundedResult() << " & ";

      if (!lBlind)
	lTexOutput << n_data.roundedResult(false) << " & ";
      else lTexOutput << "XXX" <<  " & ";

      lTexOutput << n_signal.roundedResult()
		 << " \\\\ "
		 << std::endl;

      lDatOutput[k] << "QCD " << n_qcd.rawNumber() << " " << n_qcd.rawError() << std::endl
		    << "GJets " << n_gjets.rawNumber() << " " << n_gjets.rawError() << std::endl
		    << "Top " << n_top.rawNumber() << " " << n_top.rawError() <<  std::endl
		    << "TTbar " << n_ttbar.rawNumber() << " " << n_ttbar.rawError() <<  std::endl
		    << "TTbarPowheg " << n_ttpowheg.rawNumber() << " " << n_ttpowheg.rawError() <<  std::endl
		    << "SingleTop " << n_singletop.rawNumber() << " " << n_singletop.rawError() <<  std::endl
		    << "TW " << n_tW.rawNumber() << " " << n_tW.rawError() <<  std::endl
		    << "WJets " << n_WJets_tot.rawNumber() << " " << n_WJets_tot.rawError() <<  std::endl
		    << "QCD_WJets " << n_WJets.rawNumber() << " " << n_WJets.rawError() <<  std::endl
		    << "EWK_WJets " << n_EWKW.rawNumber() << " " << n_EWKW.rawError() <<  std::endl
		    << "QCD_WJets_enu " << n_WJets_enu.rawNumber() << " " << n_WJets_enu.rawError() <<  std::endl
		    << "EWK_WJets_enu " << n_EWKW_enu.rawNumber() << " " << n_EWKW_enu.rawError() <<  std::endl
		    << "QCD_WJets_munu " << n_WJets_munu.rawNumber() << " " << n_WJets_munu.rawError() <<  std::endl
		    << "EWK_WJets_munu " << n_EWKW_munu.rawNumber() << " " << n_EWKW_munu.rawError() <<  std::endl
		    << "QCD_WJets_taunu " << n_WJets_taunu.rawNumber() << " " << n_WJets_taunu.rawError() <<  std::endl
		    << "EWK_WJets_taunu " << n_EWKW_taunu.rawNumber() << " " << n_EWKW_taunu.rawError() <<  std::endl
		    << "ZJets " << n_ZJets.rawNumber() << " " << n_ZJets.rawError() <<  std::endl
		    << "ZJets_ll " << n_ZJets_ll.rawNumber() << " " << n_ZJets_ll.rawError() <<  std::endl
		    << "ZJets_nunu " << n_ZJets_nunu.rawNumber() << " " << n_ZJets_nunu.rawError() <<  std::endl
		    << "ZJets_vbf " << n_ZJets_vbf.rawNumber() << " " << n_ZJets_vbf.rawError() <<  std::endl
		    << "ZJets_vbf_nunu " << n_ZJets_vbf_nunu.rawNumber() << " " << n_ZJets_vbf_nunu.rawError() <<  std::endl
		    << "VV " << n_VV.rawNumber() << " " << n_VV.rawError() <<  std::endl;
      if (plot_data_qcd)
	lDatOutput[k] << "Data_QCD " << n_data_qcd.rawNumber() << " " << n_data_qcd.rawError() <<  std::endl;
      lDatOutput[k] << "Data " << n_data.rawNumber() << " " << n_data.rawError() <<  std::endl
		    << "Signal " << n_signal.rawNumber() << " " << n_signal.rawError() 
		    << std::endl;

      lDatOutput[k].close();
 
    }



  }//loop on selection

  if (lFillSummaryTable){
     lTexOutput << "\\hline" << std::endl
		<< "\\end{tabular}" << std::endl;
		
     lTexOutput.close();
   }


   return 0;
   
}//main
