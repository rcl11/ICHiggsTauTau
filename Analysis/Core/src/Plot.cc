#include "Core/interface/Plot.h"

#include "Core/interface/TH1PlotElement.h"
#include "Core/interface/RatioPlotElement.h"
#include "Core/interface/TextElement.h"
#include <vector>
#include <iostream>
#include <cmath>
#include "TAxis.h"
#include "TH1.h"
#include "TROOT.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TLegend.h"
#include "THStack.h"
#include "TPad.h"
#include "TLine.h"
#include "TGraphErrors.h"
#include "boost/algorithm/string.hpp"

namespace ic {

  bool Plot::AddTH1PlotElement(ic::TH1PlotElement element) {
    if (element_map_.count(element.name()) == 0) {
      if (!element.hist_ptr()) {
        std::cerr << "Error: Cannot add TH1PlotElement with name \"" <<
          element.name() << "\", TH1F pointer is null." << std::endl;
        return false;
      } else {
        elements_.push_back(element);
        element_map_[element.name()] = element;
        if (element.in_stack()) ++n_in_stack;
      }
      return true;
    } else {
      std::cerr << "Error: Cannot add TH1PlotElement with name \"" <<
        element.name() << "\", name already in use." << std::endl;
      return false;
    }
  }

  Plot::~Plot() {

  }

  bool Plot::AddRatioPlotElement(ic::RatioPlotElement element) {
    ratios_.push_back(element);
    return true;
  }

  void Plot::AddTextElement(ic::TextElement element) {
    texts_.push_back(element);
  }
  
  void Plot::SetRatioYaxisRange(double min, double max) {
    custom_ratio_y_axis_range = true;
    ratio_y_axis_min = min;
    ratio_y_axis_max = max;
  }

  void Plot::SetXaxisRange(double min, double max) { 
    custom_x_axis_range = true;
    x_axis_min = min;
    x_axis_max = max;
  }

  void Plot::SetYaxisRange(double min, double max) { 
    custom_y_axis_range = true;
    y_axis_min = min;
    y_axis_max = max;
  }

  int Plot::GeneratePlot(std::vector<std::string> types) {

    if (use_htt_style) {
      SetHTTStyle();
    } else {
      SetTdrStyle();      
    }
    if (draw_ratio_hist) {
      //gStyle->SetCanvasDefH(900); //Height of canvas
      if (use_htt_style) {
        gStyle->SetLabelSize(0.045, "XYZ");
        gStyle->SetTitleSize(0.055, "XYZ");
      } else {
        gStyle->SetLabelSize(0.035, "XYZ");
        gStyle->SetTitleSize(0.045, "XYZ");
      }
    }

    unsigned n_elements = elements_.size();
    unsigned n_legend = 0;
    double y_pos_max = use_htt_style ? 0.88 : 0.92;
    //double y_pos_current = 0.90;
    THStack thstack("stack","stack");
    unsigned n_added_to_stack = 0;

    //1st Loop through elements - information gathering
    for (unsigned i = 0; i < n_elements; ++i) {
      TH1PlotElement const& ele = elements_[i];
      //Count number of elements with a legend entry
      if (ele.legend_text() != "") ++n_legend;
      //Overide whatever default style has been applied to the histograms
      ele.hist_ptr()->UseCurrentStyle();
    }

    //Create a TLegend object, if necessary
    TLegend *legend = 0;
    if (n_legend > 0) {
      //Allow 5% of the Canvas height for each 'row'
      double legend_y1 = y_pos_max - (n_legend * legend_height);
      //y_pos_current = legend_y1 - 0.03;
      double x1 = (legend_pos == 0) ? 0.60 : 0.16;
      if (use_htt_style) x1 = (legend_pos == 0) ? 0.55 : 0.16;
      if (legend_left > 0) x1 = legend_left;
      double x2 = (legend_pos == 0) ? 0.92 : 0.48;
      legend = new TLegend(x1,legend_y1,x2,y_pos_max,"","brNDC");
      if (use_htt_style) legend->SetTextFont(82);
    }
    //Set up canvas (and optionally pads)
    TCanvas* canv = new TCanvas("canv","canv");
    if (use_htt_style) {
      canv->SetFillColor      (0);
      canv->SetBorderMode     (0);
      canv->SetBorderSize     (10);
      // Set margins to reasonable defaults
      canv->SetLeftMargin     (0.18);
      canv->SetLeftMargin     (0.17);
      canv->SetRightMargin    (0.05);
      canv->SetTopMargin      (0.12);
      if (!draw_ratio_hist) canv->SetTopMargin      (0.08);
      canv->SetBottomMargin   (0.14);
      // Setup a frame which makes sense
      canv->SetFrameFillStyle (0);
      canv->SetFrameLineStyle (0);
      canv->SetFrameBorderMode(0);
      canv->SetFrameBorderSize(10);
      canv->SetFrameFillStyle (0);
      canv->SetFrameLineStyle (0);
      canv->SetFrameBorderMode(0);
      canv->SetFrameBorderSize(10);
    }
    canv->cd();
    TPad* upper = nullptr;
    TPad* lower = nullptr;
    if (draw_ratio_hist) {
      upper = new TPad("upper","pad",0, 0.26 ,1 ,1);
      lower = new TPad("lower","pad",0, 0   ,1 ,0.26);
      upper->SetBottomMargin(0.02);
      if (use_htt_style) upper->SetTopMargin(0.09);
      upper->Draw();
      upper->cd();
    }

    double max_bin_content = 0.0;
    bool first_drawn = false;
    for (unsigned i = 0; i < n_elements; ++i) {
      TH1PlotElement & ele = elements_[i];
      //This is a fix for histograms with negative bin content
      //for which Sumw2 has not already been called
      // for (int j = 0; j < ele.hist_ptr()->GetNbinsX(); ++j) {
      //   if (ele.hist_ptr()->GetBinContent(j) < 0) 
      //     ele.hist_ptr()->SetBinContent(j,0.0);
      // }
      if (ele.hist_ptr()->GetSumw2N() == 0) ele.hist_ptr()->Sumw2();
      //If draw_normalised is set, this takes priority over any scale factor
      if (!ele.draw_normalised()) {
        ele.hist_ptr()->Scale(ele.scale_factor());
      } else {
        //I THINK THIS IGNORES UNDER/OVERFLOW BINS - BUT SHOULD
        //REALLY BE BETTER DEFINED
        ele.hist_ptr()->Scale(1./ele.hist_ptr()->Integral());
      }

      ele.hist_ptr()->Rebin(ele.rebin_factor());
      ele.hist_ptr()->SetLineColor(ele.line_color());
      ele.hist_ptr()->SetLineWidth(ele.line_width());
      ele.hist_ptr()->SetLineStyle(ele.line_style());
      ele.hist_ptr()->SetMarkerColor(ele.marker_color());
      ele.hist_ptr()->SetMarkerSize(ele.marker_size());
      ele.hist_ptr()->SetMarkerStyle(ele.marker_style());
      ele.hist_ptr()->SetFillColor(ele.fill_color());
      ele.hist_ptr()->SetFillStyle(ele.fill_style());
      //UPDATE: MAKE THIS AN OPTION
      ele.hist_ptr()->SetStats(false);


      // if (ele.legend_text() != ""){
      //   std::string leg_options = "";
      //   if (ele.draw_fill() && ele.draw_fill_in_legend()) leg_options += "f";
      //   if (ele.draw_line()) leg_options += "l";
      //   if (ele.draw_marker() || ele.draw_stat_error_y()) leg_options += "p";
      //   legend->AddEntry(ele.hist_ptr(),
      //       ele.legend_text().c_str(),leg_options.c_str());
      // }
      std::string draw_options = "";
      if (ele.draw_fill()) draw_options += "HIST";
      if (ele.smooth_curve()) draw_options += "C";
      if (ele.draw_marker()) draw_options += "P";
      if (ele.draw_stat_error_y()) draw_options += "E1";
      if (ele.draw_options() != "") draw_options = ele.draw_options();
      //ele.hist_ptr()->Smooth(10);

      //Loop through bins of each histo (respecting any custom
      //x-axis range that is applied) to find the largest bin
      for (int j = 0; j < ele.hist_ptr()->GetNbinsX(); ++j) {
        bool bin_in_range = true;
        if (custom_x_axis_range) {
          if ((ele.hist_ptr()->GetBinLowEdge(j) < x_axis_min)
              || (ele.hist_ptr()->GetXaxis()->GetBinUpEdge(j) > x_axis_max))
            bin_in_range = false;
        }
        if (bin_in_range) {
          if (ele.hist_ptr()->GetBinContent(j) > max_bin_content) {
            max_bin_content = ele.hist_ptr()->GetBinContent(j);
          }
        }
      }

      //Special actions for first histogram drawn
      if (i == 0) {
        ele.hist_ptr()->GetXaxis()->SetTitle(x_axis_title.c_str());
        if (use_htt_style) {
          ele.hist_ptr()->SetTitleSize  (0.055,"Y");
          ele.hist_ptr()->SetTitleOffset(1.400,"Y");
          if (!draw_ratio_hist) ele.hist_ptr()->SetTitleOffset(1.500,"Y");
          ele.hist_ptr()->SetLabelOffset(0.010,"Y");
          ele.hist_ptr()->SetLabelSize  (0.040,"Y");
          ele.hist_ptr()->SetLabelFont  (42   ,"Y");
          ele.hist_ptr()->SetTitleSize  (0.055,"X");
          ele.hist_ptr()->SetTitleOffset(1.100,"X");
          ele.hist_ptr()->SetLabelOffset(0.010,"X");
          ele.hist_ptr()->SetLabelSize  (0.040,"X");
          ele.hist_ptr()->SetLabelFont  (42   ,"X");
          // ele.hist_ptr()->SetMarkerStyle(20);
          // ele.hist_ptr()->SetMarkerColor(color);
          // ele.hist_ptr()->SetMarkerSize (0.6);
          ele.hist_ptr()->GetYaxis()->SetTitleFont(62);
          ele.hist_ptr()->GetXaxis()->SetTitleFont(62);
        }
        if (x_axis_min > -9999 && x_axis_max > 0){
          ele.hist_ptr()->GetXaxis()->SetRangeUser(x_axis_min,x_axis_max);
        }
        if (y_axis_min > 0 && y_axis_max > 0){
          ele.hist_ptr()->GetYaxis()->SetRangeUser(y_axis_min,y_axis_max);
        }

        if (x_bin_labels_ != "") {
          std::vector<std::string> bin_labels;
          boost::split(bin_labels, x_bin_labels_, boost::is_any_of(":"));

          if (long(bin_labels.size()) <= long(ele.hist_ptr()->GetNbinsX())) {
            for (unsigned binl = 0; binl < bin_labels.size(); ++binl) {
              ele.hist_ptr()->GetXaxis()->SetBinLabel(binl+1, bin_labels[binl].c_str());
            }
          }

        }
  
        if (draw_ratio_hist) {
          ele.hist_ptr()->GetXaxis()->SetLabelSize(0.0);
          ele.hist_ptr()->GetXaxis()->SetTitleSize(0.0);
        }

        ele.hist_ptr()->GetYaxis()->SetTitle(y_axis_title.c_str());
        //ele.hist_ptr()->SetName("Madgraph");
        if (!ele.in_stack()) {
          ele.hist_ptr()->Draw(draw_options.c_str());
          first_drawn = true;
        } else {
          thstack.Add(ele.hist_ptr(), draw_options.c_str());
          if (n_added_to_stack == 0) thstack.SetHistogram((TH1F*)ele.hist_ptr()->Clone());
          ++n_added_to_stack;
          if (n_added_to_stack == n_in_stack) {
            thstack.Draw();
            first_drawn = true;
          }
        }

        canv->Update();
        /*
           TPaveStats *st = (TPaveStats*)ele.hist_ptr()->FindObject("stats");
           st->SetX1NDC(0.65); 
           st->SetX2NDC(0.89); 
           st->SetY1NDC(0.6); 
           st->SetY2NDC(0.75); 
           */
      } else {
        //ele.hist_ptr()->SetName("aMC@NLO");
        if (!ele.in_stack()) {
          std::string same_opt = "";
          if (first_drawn) same_opt = "SAME";
          first_drawn = true;
          ele.hist_ptr()->Draw((same_opt+draw_options).c_str());
        } else {
          thstack.Add(ele.hist_ptr());
          ++n_added_to_stack;
          if (n_added_to_stack == n_in_stack) {
            std::string same_opt = "";
            if (first_drawn) same_opt = "SAME";
            first_drawn = true;
            thstack.Draw((same_opt+"HIST").c_str());
          }
        }
        canv->Update();
        /*
           TPaveStats *st = (TPaveStats*)ele.hist_ptr()->FindObject("stats");
           st->SetX1NDC(0.65); 
           st->SetX2NDC(0.89); 
           st->SetY1NDC(0.45); 
           st->SetY2NDC(0.6); 
           */
      }
    }
    unsigned j = n_elements -1;
    for (unsigned i = 0; i < n_elements; ++i, --j) {
      TH1PlotElement & ele = elements_[j];
      if (ele.legend_text() != ""){
        std::string leg_options = "";
        if (ele.draw_fill() && ele.draw_fill_in_legend()) leg_options += "f";
        if (ele.draw_line()) leg_options += "l";
        if (ele.draw_marker() || ele.draw_stat_error_y()) leg_options += "p";
        legend->AddEntry(ele.hist_ptr(),
            ele.legend_text().c_str(),leg_options.c_str());
      }
    }

    if (draw_band_on_stack_) {
      std::vector<std::string> samples;
      boost::split(samples, samples_for_band_, boost::is_any_of("+"));
      //std::cout << "Found " << samples.size() << " samples to sum for band" << std::endl;
      TH1F* error_band = (TH1F*)element_map_[samples[0]].hist_ptr()->Clone();
      for (unsigned i = 1; i < samples.size(); ++i) {
        //std::cout << "Adding sample: " << samples[i] << std::endl;
        error_band->Add(element_map_[samples[i]].hist_ptr());
      }
      for (unsigned i = 1; i <= unsigned(error_band->GetNbinsX()); ++i) {
        error_band->SetBinError(i, error_band->GetBinContent(i)*band_size_fractional_);
      }
      error_band->SetMarkerSize(0);
      error_band->SetFillColor(1);
      error_band->SetFillStyle(3013);
      error_band->SetLineWidth(1);
      error_band->Draw("E2SAME");
      legend->AddEntry(error_band, "bkg. uncertainty" , "F" );
    }


    TLatex *title_latex = new TLatex();
    title_latex->SetNDC();
    if (use_htt_style) {
      title_latex->SetTextSize(0.04);
      if (draw_ratio_hist) title_latex->SetTextSize(0.05);
      title_latex->SetTextFont(62);
      title_latex->SetTextAlign(31);
      double height = draw_ratio_hist ? 0.93 : 0.94;
      title_latex->DrawLatex(0.95,height+0.01,title_right.c_str());
      title_latex->SetTextAlign(11);
      title_latex->DrawLatex(0.17,height,title_left.c_str());
    } else {
      title_latex->SetTextSize(0.03);
      title_latex->SetTextAlign(31);
      title_latex->DrawLatex(0.95,0.965,title_right.c_str());
      title_latex->SetTextAlign(11);
      title_latex->DrawLatex(0.16,0.965,title_left.c_str());
    }

    for (unsigned te = 0; te < texts_.size(); ++te) {
      title_latex->SetTextSize(texts_[te].size());
      title_latex->SetTextAlign(11);
      title_latex->DrawLatex(texts_[te].x_pos(), texts_[te].y_pos(), texts_[te].text().c_str());
    }


    //Adjust the y-axis maximum to account for largest bin
    if (y_axis_log) {
      double tgt_max = max_bin_content*1.5*extra_pad;
      elements_[0].hist_ptr()->SetMaximum(tgt_max);
      if (y_axis_min > 0) {
          elements_[0].hist_ptr()->SetMinimum(y_axis_min);
      }
      if (thstack.GetHistogram()) {
        thstack.SetMaximum(tgt_max);
        // This is absolutely crazy. Turns out ROOT ignores
        // the min and max values to a THStack when drawing
        // with a log scale, but rather adjusts them to some
        // wider range for you.  This code ensures we get the
        // actual range we want. See here for a better long-
        // term fix: http://root.cern.ch/phpBB3/viewtopic.php?f=3&t=12702
        // (basically figure out axis ranges first and draw with the pad)
        if (y_axis_min > 0) {
            thstack.SetMaximum(tgt_max/(1+0.2*std::log10(tgt_max/y_axis_min)));
            thstack.SetMinimum(y_axis_min*(1+0.5*std::log10(tgt_max/y_axis_min)));
        }
      }
    } else {
      if (!custom_y_axis_range) {
        elements_[0].hist_ptr()->SetMaximum(max_bin_content*1.1*extra_pad);
        if (thstack.GetHistogram()) thstack.SetMaximum(thstack.GetMaximum()*1.1*extra_pad);
      } else {
        elements_[0].hist_ptr()->SetMaximum(y_axis_max);
        if (thstack.GetHistogram()) thstack.SetMaximum(y_axis_max);
        elements_[0].hist_ptr()->SetMinimum(y_axis_min);
        if (thstack.GetHistogram()) thstack.SetMinimum(y_axis_min);
      }
    }

    /* 
    // TLine *line = new TLine(30., 0., 30., elements_[0].hist_ptr()->GetMaximum()); // mT
    TLine *line = new TLine(-20., 0., -20., 12000);   // Dzeta
    // TLine *line = new TLine(-0.5, 0., -0.5, 5000000);   // BDT
    TLine *line2 = new TLine(70., 0., 70., elements_[0].hist_ptr()->GetMaximum());
    TLatex *ex_latex = new TLatex();
    ex_latex->SetTextSize(0.035);
    ex_latex->SetTextFont(62);
    ex_latex->SetTextAlign(31);
    //ex_latex->DrawLatex(27,15000,"#splitline{Baseline}{selection}"); // mT
    ex_latex->DrawLatex(51,9500,"#splitline{MSSM Baseline}{selection}"); // Dzeta
    //ex_latex->DrawLatex(0.1,6000,"#splitline{SM Baseline}{selection}"); // BDT
    ex_latex->SetTextAlign(11);
    //ex_latex->DrawLatex(100,6000,"#splitline{High-m_{T}}{control region}");
    line->SetLineWidth(3);
    line->SetLineStyle(3);
    line->SetLineColor(1);
    line2->SetLineWidth(3);
    line2->SetLineStyle(7);
    line2->SetLineColor(1);
    canv->Update();
    line->Draw();
    //line2->Draw();
    */
    
    canv->RedrawAxis();
    canv->Update();


    //Apply legend style options
    if (n_legend > 0) {
      legend->SetBorderSize(1);
      legend->SetTextFont(42);
      if (use_htt_style) legend->SetTextFont(62);
      legend->SetLineColor(0);
      legend->SetLineStyle(1);
      legend->SetLineWidth(1);
      legend->SetFillColor(0);
      legend->SetFillStyle(1001);
      legend->Draw();
    }


    if (draw_ratio_hist) {
      canv->cd();
      lower->SetTopMargin(0.026);
      lower->SetBottomMargin(0.24);
      if (use_htt_style) lower->SetBottomMargin(0.35);
      lower->Draw();
      lower->cd();
      std::vector<TH1F*> ratio_ele;
      ratio_ele.resize(ratios_.size());
      for (unsigned k = 0; k < ratios_.size(); ++k) {
        TH1F* num = 0;
        TH1F* den = 0;
        if (ratios_[k].multi_mode()) {
          std::vector<std::string> num_strings;
          std::vector<std::string> den_strings;
          std::string num_str = ratios_[k].hist_numerator();
          std::string den_str = ratios_[k].hist_denominator();
          boost::split(num_strings, num_str, boost::is_any_of("+"));
          boost::split(den_strings, den_str, boost::is_any_of("+"));
          for (unsigned n = 0; n < num_strings.size(); ++n) {
            for (unsigned m = 0; m < elements_.size(); ++m) {
              if (elements_[m].name() == num_strings[n]) {
                if (!num) {
                  num = (TH1F*)elements_[m].hist_ptr()->Clone();
                } else {
                  num->Add(elements_[m].hist_ptr());
                }
              }
            }
          }
          for (unsigned n = 0; n < den_strings.size(); ++n) {
            for (unsigned m = 0; m < elements_.size(); ++m) {
              if (elements_[m].name() == den_strings[n]) {
                if (!den) {
                  den = (TH1F*)elements_[m].hist_ptr()->Clone();
                } else {
                  den->Add(elements_[m].hist_ptr());
                }
              }
            }
          }
        } else {
          for (unsigned m = 0; m < elements_.size(); ++m) {
            if (elements_[m].name() == ratios_[k].hist_numerator())
              num = elements_[m].hist_ptr();
            if (elements_[m].name() == ratios_[k].hist_denominator())
              den = elements_[m].hist_ptr();
          }
        }
        if (num && den) {

          ratio_ele[k] = (TH1F*)num->Clone();
          if (x_axis_min > -9999 && x_axis_max > 0){
            ratio_ele[k]->GetXaxis()->SetRangeUser(x_axis_min,x_axis_max);
          }

          TGraphErrors *grErr = nullptr;

          if (!draw_signif) {
            ratio_ele[k]->Divide(den);
          } else {
            TH1F *h1 = (TH1F*)ratio_ele[k]->Clone("data1");
            TH1F *h2 = (TH1F*)ratio_ele[k]->Clone("data2");
            TH1F *h3 = (TH1F*)ratio_ele[k]->Clone("data3");
            TH1F *h4 = (TH1F*)ratio_ele[k]->Clone("data4");
            h3->Add(den,-1);
            int nbins = num->GetNbinsX();
            for(int i = 0; i < nbins; i++) {
              double err=0;
              if(h1->GetBinContent(i+1)==0 ) err=0;
              else err=(0.5*h1->GetBinError(i+1)/sqrt(h1->GetBinContent(i+1)));
              h2->SetBinContent(i+1,sqrt(h1->GetBinContent(i+1)));
              h2->SetBinError(i+1,err);
            }
            h3->Divide(h2);
            for(int i = 0; i < nbins; i++) {
              h4->SetBinContent(i+1,0);
              //h4->SetBinError(i+1,h3->GetBinError(i+1));
              h4->SetBinError(i+1,num->GetBinError(i+1) / den->GetBinError(i+1));
              h3->SetBinError(i+1,0);
              if(h3->GetBinContent(i+1)==0) h3->SetBinContent(i+1,-10);      // hack, because I can't figure out how to suppress the points for bins with zeroes
            }
            grErr = new TGraphErrors(h4);
            grErr->SetFillColor(kPink-5);
            grErr->SetFillStyle(3002);
            ratio_ele[k] = h3;
          }

          ratio_ele[k]->SetLineColor(ratios_[k].line_color());
          ratio_ele[k]->SetLineWidth(ratios_[k].line_width());
          ratio_ele[k]->SetLineStyle(ratios_[k].line_style());
          ratio_ele[k]->SetMarkerColor(ratios_[k].marker_color());
          ratio_ele[k]->SetMarkerSize(ratios_[k].marker_size());
          ratio_ele[k]->SetMarkerStyle(ratios_[k].marker_style());
          ratio_ele[k]->SetFillColor(ratios_[k].fill_color());
          ratio_ele[k]->SetFillStyle(ratios_[k].fill_style());
           std::string draw_options = "";
          if (ratios_[k].draw_fill()) draw_options += "HIST";
          if (ratios_[k].draw_marker()) draw_options += "P";
          if (ratios_[k].draw_stat_error_y()) {
            draw_options += use_htt_style ? "E0" : "E1";
          }
          if (ratios_[k].draw_options() != "") draw_options = ratios_[k].draw_options();

          if (k == 0)
          {
            if (use_htt_style) {
              ratio_ele[k]->GetXaxis()->SetLabelSize(0.11);
              ratio_ele[k]->GetXaxis()->SetTitleFont(62);
              ratio_ele[k]->GetYaxis()->SetTitleFont(62);
              ratio_ele[k]->GetXaxis()->SetTitleSize(0.16);
              ratio_ele[k]->GetYaxis()->SetTitleSize(0.16);
              ratio_ele[k]->GetYaxis()->SetTitleOffset(0.460);
              ratio_ele[k]->GetXaxis()->SetTitleOffset(0.900);
              if (x_bin_labels_ != "") ratio_ele[k]->GetXaxis()->SetLabelSize(0.20);
            } else {
              ratio_ele[k]->GetXaxis()->SetLabelSize(0.10);
              ratio_ele[k]->GetXaxis()->SetTitleSize(0.12);
              ratio_ele[k]->GetYaxis()->SetTitleSize(0.12);
              ratio_ele[k]->GetYaxis()->SetTitleOffset(0.55);
            }
            if (draw_signif) ratio_ele[k]->GetYaxis()->SetTitleOffset(0.45);
            if (x_bin_labels_ != "") {
              std::vector<std::string> bin_labels;
              boost::split(bin_labels, x_bin_labels_, boost::is_any_of(":"));
              if (long(bin_labels.size()) <= long(ratio_ele[k]->GetNbinsX())) {
                for (unsigned binl = 0; binl < bin_labels.size(); ++binl) {
                  ratio_ele[k]->GetXaxis()->SetBinLabel(binl+1, bin_labels[binl].c_str());
                }
              }

            }
            ratio_ele[k]->GetYaxis()->SetLabelSize(0.10);
            if (use_htt_style) ratio_ele[k]->GetYaxis()->SetLabelSize(0.11);
            ratio_ele[k]->GetXaxis()->SetTitle(x_axis_title.c_str());
            ratio_ele[k]->GetYaxis()->SetNdivisions(505);
            if (custom_ratio_y_axis_range) {
              ratio_ele[k]->SetMaximum(ratio_y_axis_max);
              ratio_ele[k]->SetMinimum(ratio_y_axis_min);
            }
            ratio_ele[k]->GetYaxis()->SetTitle(ratio_y_axis_title.c_str());
            ratio_ele[k]->SetStats(false);
            if (draw_y_gridlines) lower->SetGrid(0,1);
            ratio_ele[k]->Draw(draw_options.c_str());
          } else {
            ratio_ele[k]->Draw(("SAME"+draw_options).c_str());
          }
          //if (draw_signif) grErr->Draw("2");
        }
      }
    }

    canv->cd();

    if (y_axis_log) {
      if (draw_ratio_hist) {
        upper->SetLogy(true);
      } else {
        canv->SetLogy(true);
      }
    }
    canv->Update();
    std::size_t pos = output_filename.find(".pdf");
    if (pos != output_filename.npos) {
      output_filename = output_filename.substr(0, pos);
    }

    for (auto const& type : types) {
      if (type == "pdf") {
        if (append > 0) {
          if (append == 1) canv->Print((output_filename+".pdf"+"(").c_str(),"pdf");
          if (append == 2) canv->Print((output_filename+".pdf").c_str(),"pdf");
          if (append == 3) canv->Print((output_filename+".pdf"+")").c_str(),"pdf");
        } else {
          canv->SaveAs((output_filename+".pdf").c_str());
        }
      } else {
        canv->SaveAs((output_filename+"."+type).c_str());
      }
    }
    delete legend;
    delete title_latex;
    delete upper;
    delete lower;
    delete canv;

    
    return 0;
  }

  void Plot::SetTdrStyle() {
    gStyle->SetCanvasBorderMode(0);
    gStyle->SetCanvasColor(kWhite);
    gStyle->SetCanvasDefH(700); //Height of canvas
    gStyle->SetCanvasDefW(700); //Width of canvas

    gStyle->SetCanvasDefX(0);   //POsition on screen
    gStyle->SetCanvasDefY(0);

    // For the Pad:
    gStyle->SetPadBorderMode(0);
    // gStyle->SetPadBorderSize(Width_t size = 1);
    gStyle->SetPadColor(kWhite);
    gStyle->SetPadGridX(false);
    gStyle->SetPadGridY(false);
    gStyle->SetGridColor(0);
    gStyle->SetGridStyle(3);
    gStyle->SetGridWidth(1);

    // For the frame:
    gStyle->SetFrameBorderMode(0);
    gStyle->SetFrameBorderSize(10);
    gStyle->SetFrameFillColor(0);
    gStyle->SetFrameFillStyle(0);
    gStyle->SetFrameLineColor(1);
    gStyle->SetFrameLineStyle(0);
    gStyle->SetFrameLineWidth(2);

    // For the histo:
    // gStyle->SetHistFillColor(1);
    // gStyle->SetHistFillStyle(0);
    gStyle->SetHistLineColor(1);
    gStyle->SetHistLineStyle(0);
    gStyle->SetHistLineWidth(short(0.1));
    // gStyle->SetLegoInnerR(Float_t rad = 0.5);
    // gStyle->SetNumberContours(Int_t number = 20);

    gStyle->SetEndErrorSize(0);
    //gStyle->SetErrorMarker(20);
    gStyle->SetErrorX(0.5);

    gStyle->SetMarkerStyle(20);

    //For the fit/function:
    gStyle->SetOptFit(1);
    gStyle->SetFitFormat("5.4g");
    gStyle->SetFuncColor(2);
    gStyle->SetFuncStyle(1);
    gStyle->SetFuncWidth(1);

    //For the date:
    gStyle->SetOptDate(0);
    // gStyle->SetDateX(Float_t x = 0.01);
    // gStyle->SetDateY(Float_t y = 0.01);

    // For the statistics box:
    //  gStyle->SetOptFile(0);
    gStyle->SetOptStat(1);
    gStyle->SetOptStat("MRn");
    gStyle->SetStatColor(kWhite);
    gStyle->SetStatFont(42);
    gStyle->SetStatFontSize(0.04);///---> gStyle->SetStatFontSize(0.025);
    gStyle->SetStatTextColor(1);
    gStyle->SetStatFormat("6.4g");
    gStyle->SetStatBorderSize(1);
    gStyle->SetStatH(0.1);
    gStyle->SetStatW(0.2);///---> gStyle->SetStatW(0.15);

    // gStyle->SetStatStyle(Style_t style = 1001);
    // gStyle->SetStatX(Float_t x = 0);
    // gStyle->SetStatY(Float_t y = 0);

    // Margins:
    gStyle->SetPadTopMargin(0.05);
    gStyle->SetPadBottomMargin(0.10);
    gStyle->SetPadLeftMargin(0.16);
    gStyle->SetPadRightMargin(0.05);

    // For the Global title:

    gStyle->SetOptTitle(0);
    gStyle->SetTitleFont(42);
    gStyle->SetTitleColor(1);
    gStyle->SetTitleTextColor(1);
    gStyle->SetTitleFillColor(10);
    gStyle->SetTitleFontSize(0.05);
    // gStyle->SetTitleH(0); // Set the height of the title box
    // gStyle->SetTitleW(0); // Set the width of the title box
    // gStyle->SetTitleX(0); // Set the position of the title box
    // gStyle->SetTitleY(0.985); // Set the position of the title box
    // gStyle->SetTitleStyle(Style_t style = 1001);
    // gStyle->SetTitleBorderSize(2);

    // For the axis titles:

    gStyle->SetTitleColor(1, "XYZ");
    gStyle->SetTitleFont(42, "XYZ");
    gStyle->SetTitleSize(0.05, "XYZ");
    // gStyle->SetTitleXSize(Float_t size = 0.02); // Another way to set the size?
    // gStyle->SetTitleYSize(Float_t size = 0.02);
    gStyle->SetTitleXOffset(0.9);
    gStyle->SetTitleYOffset(1.55);
    // gStyle->SetTitleOffset(1.1, "Y"); // Another way to set the Offset

    // For the axis labels:

    gStyle->SetLabelColor(1, "XYZ");
    gStyle->SetLabelFont(42, "XYZ");
    gStyle->SetLabelOffset(0.007, "XYZ");
    gStyle->SetLabelSize(0.04, "XYZ");

    // For the axis:

    gStyle->SetAxisColor(1, "XYZ");
    gStyle->SetStripDecimals(kTRUE);
    gStyle->SetTickLength(0.02, "XYZ");
    gStyle->SetNdivisions(510, "XYZ");
    gStyle->SetPadTickX(1);  // To get tick marks on the opposite side of the frame
    gStyle->SetPadTickY(1);

    // Change for log plots:
    gStyle->SetOptLogx(0);
    gStyle->SetOptLogy(0);
    gStyle->SetOptLogz(0);
    gStyle->SetLineStyleString(11,"20 10");

    // Postscript options:
    //  gStyle->SetPaperSize(20.,20.);
    // gStyle->SetLineScalePS(Float_t scale = 3);
    // gStyle->SetLineStyleString(Int_t i, const char* text);
    // gStyle->SetHeaderPS(const char* header);
    // gStyle->SetTitlePS(const char* pstitle);

    // gStyle->SetBarOffset(Float_t baroff = 0.5);
    // gStyle->SetBarWidth(Float_t barwidth = 0.5);
    // gStyle->SetPaintTextFormat(const char* format = "g");
    // gStyle->SetPalette(Int_t ncolors = 0, Int_t* colors = 0);
    // gStyle->SetTimeOffset(Double_t toffset);
    // gStyle->SetHistMinimumZero(kTRUE);

    //   gStyle->cd();

    gROOT->ForceStyle();
  }

  void Plot::SetHTTStyle() {
    TStyle *HttStyle = new TStyle("Htt-Style","The Perfect Style for Plots ;-)");
    gStyle = HttStyle;

    // Canvas
    HttStyle->SetCanvasColor     (0);
    HttStyle->SetCanvasBorderSize(10);
    HttStyle->SetCanvasBorderMode(0);
    HttStyle->SetCanvasDefH      (600);
    HttStyle->SetCanvasDefW      (600);
    HttStyle->SetCanvasDefX      (100);
    HttStyle->SetCanvasDefY      (100);

    // color palette for 2D temperature plots
    HttStyle->SetPalette(1,0);
    
    // Pads
    HttStyle->SetPadColor       (0);
    HttStyle->SetPadBorderSize  (10);
    HttStyle->SetPadBorderMode  (0);
    HttStyle->SetPadBottomMargin(0.13);
    HttStyle->SetPadTopMargin   (0.08);
    HttStyle->SetPadLeftMargin  (0.17);
    HttStyle->SetPadRightMargin (0.05);
    HttStyle->SetPadGridX       (0);
    HttStyle->SetPadGridY       (0);
    HttStyle->SetPadTickX       (1);
    HttStyle->SetPadTickY       (1);

    // Frames
    HttStyle->SetLineWidth(3);
    HttStyle->SetFrameFillStyle ( 0);
    HttStyle->SetFrameFillColor ( 0);
    HttStyle->SetFrameLineColor ( 1);
    HttStyle->SetFrameLineStyle ( 0);
    HttStyle->SetFrameLineWidth ( 2);
    HttStyle->SetFrameBorderSize(10);
    HttStyle->SetFrameBorderMode( 0);

    // Histograms
    HttStyle->SetHistFillColor(2);
    HttStyle->SetHistFillStyle(0);
    HttStyle->SetHistLineColor(1);
    HttStyle->SetHistLineStyle(0);
    HttStyle->SetHistLineWidth(3);
    HttStyle->SetNdivisions(510);

    // Functions
    HttStyle->SetFuncColor(1);
    HttStyle->SetFuncStyle(0);
    HttStyle->SetFuncWidth(2);

    // Various
    HttStyle->SetMarkerStyle(20);
    HttStyle->SetMarkerColor(kBlack);
    HttStyle->SetMarkerSize (1.1);

    HttStyle->SetTitleBorderSize(0);
    HttStyle->SetTitleFillColor (0);
    HttStyle->SetTitleX         (0.2);

    HttStyle->SetTitleSize  (0.055,"X");
    HttStyle->SetTitleOffset(1.200,"X");
    HttStyle->SetLabelOffset(0.005,"X");
    HttStyle->SetLabelSize  (0.040,"X");
    HttStyle->SetLabelFont  (42   ,"X");

    HttStyle->SetStripDecimals(kFALSE);
    HttStyle->SetLineStyleString(11,"20 10");

    HttStyle->SetTitleSize  (0.055,"Y");
    HttStyle->SetTitleOffset(1.600,"Y");
    HttStyle->SetLabelOffset(0.010,"Y");
    HttStyle->SetLabelSize  (0.040,"Y");
    HttStyle->SetLabelFont  (42   ,"Y");

    HttStyle->SetOptTitle(0);


    HttStyle->SetTextSize   (0.055);
    HttStyle->SetTextFont   (42);

    HttStyle->SetStatFont   (42);
    HttStyle->SetTitleFont  (42);
    HttStyle->SetTitleFont  (42,"X");
    HttStyle->SetTitleFont  (42,"Y");
    HttStyle->SetEndErrorSize(0);


    HttStyle->SetOptStat    (0);

    gROOT->ForceStyle();
    return;

  }




}
