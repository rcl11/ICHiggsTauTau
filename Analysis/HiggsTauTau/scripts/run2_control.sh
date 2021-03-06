
if (( "$#" != "4" ))
then
  echo "Usage: new_plot_control.sh [PROD] [YEAR] [ANALYSIS] [POSTFIT]"
  exit
fi


PROD=$1
YEAR=$2
ANA=$3
POSTFIT=$4
PARAMS=./scripts/"$PROD"_params_$YEAR.dat
SVSMBINS="0,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180,190,200"
SMBINS="0,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150"
#SMBINS="0,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180,190,200,225,250,275,300,325,350"
#MSSMBINS="0,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180,190,200,225,250,275,300,325,350,400,450,500,550,600,650,700,750,800,850,900,950,1000,1100,1200,1300,1400,1500"
MSSMBINS="0,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180,190,200,225,250,275,300,325,350,400,500,700,900"

if [ "$POSTFIT" -gt "0" ]
  then
    echo "Applying SF's ..."
    ET_INC_SHIFT=" --draw_error_band=true --auto_error_band=0.01036 --shift_backgrounds=QCD:0.97707,TT:0.98411,VV:0.98437,W:1.47733,ZJ:0.92786,ZL:1.37040,ZTT:0.81560 "
    MT_INC_SHIFT=" --draw_error_band=true --auto_error_band=0.00711 --shift_backgrounds=QCD:0.90291,TT:0.99461,VV:0.98942,W:1.25993,ZJ:0.93176,ZL:1.87241,ZTT:0.88709 "
    EM_INC_SHIFT=" --draw_error_band=true --auto_error_band=0.00684 --shift_backgrounds=QCD:1.02079,TT:0.97237,VV:0.98234,W:0.99838,ZLL:0.95286,ZTT:0.91173 "
    TT_INC_SHIFT=" --draw_error_band=true --auto_error_band=0.02456 --shift_backgrounds=QCD:1.00606,TT:0.93180,VV:0.97842,W:1.01697,ZJ:0.98546,ZL:0.96595,ZTT:0.92505 "
    ET_BAND_ONLY=" --draw_error_band=true --auto_error_band=0.01036 "
    MT_BAND_ONLY=" --draw_error_band=true --auto_error_band=0.00711 "
    EM_BAND_ONLY=" --draw_error_band=true --auto_error_band=0.00684 "
    TT_BAND_ONLY=" --draw_error_band=true --auto_error_band=0.02456 "
  else
   echo "Setting prefit uncertainty..."
    ET_INC_SHIFT=" --draw_error_band=true --auto_error_band=0.07559 " 
    MT_INC_SHIFT=" --draw_error_band=true --auto_error_band=0.08397 "
    EM_INC_SHIFT=" --draw_error_band=true --auto_error_band=0.07061 "
    TT_INC_SHIFT=" --draw_error_band=true --auto_error_band=0.10421 "
    ET_BAND_ONLY=$ET_INC_SHIFT
    MT_BAND_ONLY=$MT_INC_SHIFT
    EM_BAND_ONLY=$EM_INC_SHIFT
    TT_BAND_ONLY=$TT_INC_SHIFT
fi


################################################################################
##### Inclusive selection plots
################################################################################
#### SVFit Mass
'''
#./bin/HiggsTauTauPlot5 --cfg=scripts/new_plot_"$ANA"_"$YEAR".cfg --channel=et --set_alias="sel:mt_1<30." \
#  --method=8 --var="m_sv"["$SMBINS"] --cat="inclusive" \
#  --x_axis_label="M_{#tau#tau} [GeV]" --y_axis_label="dN/dm_{#tau#tau} [1/GeV]"\
#  --norm_bins=true --datacard="inclusive" $ET_INC_SHIFT \
#  --background_scheme="et_default"

#./bin/HiggsTauTauPlot5 --cfg=scripts/new_plot_"$ANA"_"$YEAR".cfg --channel=mt --set_alias="sel:mt_1<30." \
#  --method=8 --var="m_sv"["$SMBINS"] --cat="inclusive" \
#  --x_axis_label="M_{#tau#tau} [GeV]" --y_axis_label="dN/dm_{#tau#tau} [1/GeV]" $MT_INC_SHIFT \
#  --norm_bins=true --datacard="inclusive" \
#  --background_scheme="mt_with_zmm"

#./bin/HiggsTauTauPlot5 --cfg=scripts/new_plot_"$ANA"_"$YEAR".cfg --channel=em --set_alias="sel:pzeta>-20." \
#  --method=8 --var="m_sv"["$SMBINS"] --cat="inclusive" \
#  --x_axis_label="M_{#tau#tau} [GeV]" --y_axis_label="dN/dm_{#tau#tau} [1/GeV]" $EM_INC_SHIFT \
#  --norm_bins=true --datacard="inclusive" \
#  --background_scheme="em_default"

#### SVFit Mass MSSM

#./bin/HiggsTauTauPlot5 --cfg=scripts/new_plot_"$ANA"_"$YEAR".cfg --channel=et --set_alias="sel:mt_1<30." \
#  --method=8 --var="m_sv"["$MSSMBINS"] --cat="inclusive" \
#  --x_axis_label="M_{#tau#tau} [GeV]" $ET_INC_SHIFT \
#  --norm_bins=true --datacard="inclusive" --log_y=true --draw_ratio=true \
#  --custom_y_axis_min=true --y_axis_min=0.0099 \
#  --blind=false --x_blind_min=100 --x_blind_max=2000 \
#  --background_scheme="et_default"

#./bin/HiggsTauTauPlot5 --cfg=scripts/new_plot_"$ANA"_"$YEAR".cfg --channel=mt --set_alias="sel:mt_1<30." \
#  --method=8 --var="m_sv"["$MSSMBINS"] --cat="inclusive" \
#  --x_axis_label="M_{#tau#tau} [GeV]" $MT_INC_SHIFT \
#  --norm_bins=true --datacard="inclusive" --log_y=true --draw_ratio=true \
#  --custom_y_axis_min=true --y_axis_min=0.0099 \
#  --blind=false --x_blind_min=100 --x_blind_max=2000 \
#  --background_scheme="mt_with_zmm"
#
#./bin/HiggsTauTauPlot5 --cfg=scripts/new_plot_"$ANA"_"$YEAR".cfg --channel=em --set_alias="sel:pzeta>-20." \
#  --method=8 --var="m_sv"["$MSSMBINS"] --cat="inclusive" \
#  --x_axis_label="M_{#tau#tau} [GeV]" $EM_INC_SHIFT \
#  --norm_bins=true --datacard="inclusive" --log_y=true --draw_ratio=true \
#  --custom_y_axis_min=true --y_axis_min=0.0099 \
#  --blind=false --x_blind_min=100 --x_blind_max=2000 \
#  --background_scheme="em_default"
##
#### Visible Mass

./bin/HiggsTauTauPlot5 --cfg=scripts/new_plot_"$ANA"_"$YEAR".cfg --channel=et --set_alias="sel:mt_1<40." \
  --method=12 --var="m_vis"["$SMBINS"] --cat="inclusive" \
  --x_axis_label="M_{#tau#tau}^{vis} [GeV]" $ET_INC_SHIFT \
  --norm_bins=true --datacard="inclusive"\
  --background_scheme="et_default" --extra_pad=0.2 --draw_error_band=true --auto_error_band=0.00001 --qcd_os_ss_factor=1.0

./bin/HiggsTauTauPlot5 --cfg=scripts/new_plot_"$ANA"_"$YEAR".cfg --channel=mt --set_alias="sel:mt_1<30." \
  --method=12 --var="m_vis"["$SMBINS"] --cat="inclusive" \
  --x_axis_label="M_{#tau#tau}^{vis} [GeV]" $MT_INC_SHIFT \
  --norm_bins=true --datacard="inclusive"\
  --background_scheme="mt_with_zmm" --extra_pad=0.2 --draw_error_band=true --auto_error_band=0.00001 --qcd_os_ss_factor=1.17

./bin/HiggsTauTauPlot5 --cfg=scripts/new_plot_"$ANA"_"$YEAR".cfg --channel=em --set_alias="sel:pzeta>-20." \
  --method=8 --var="m_vis"["$SMBINS"] --cat="inclusive" \
 --x_axis_label="M_{#tau#tau}^{vis} [GeV]" $EM_INC_SHIFT \
  --norm_bins=true --datacard="inclusive"\
  --background_scheme="em_default" --extra_pad=0.2 --draw_error_band=true --auto_error_band=0.00001 --qcd_os_ss_factor=2

./bin/HiggsTauTauPlot5 --cfg=scripts/new_plot_"$ANA"_"$YEAR".cfg --channel=tt --set_alias="sel:1" \
  --method=8 --var="m_vis"["$SMBINS"] --cat="inclusive" \
 --x_axis_label="M_{#tau#tau}^{vis} [GeV]" $EM_INC_SHIFT \
  --norm_bins=true --datacard="inclusive"\
  --background_scheme="tt_default" --extra_pad=0.2 --draw_error_band=true --auto_error_band=0.00001

#### Total Transverse Mass

./bin/HiggsTauTauPlot5 --cfg=scripts/new_plot_"$ANA"_"$YEAR".cfg --channel=et --set_alias="sel:mt_1<40." \
  --method=12 --var="mt_tot"["$SVSMBINS"] --cat="inclusive" \
  --x_axis_label="M_{T}^{tot} [GeV]" $ET_INC_SHIFT \
  --norm_bins=true --datacard="inclusive"\
  --background_scheme="et_default" --extra_pad=0.2 --draw_error_band=true --auto_error_band=0.00001 --qcd_os_ss_factor=1.0

./bin/HiggsTauTauPlot5 --cfg=scripts/new_plot_"$ANA"_"$YEAR".cfg --channel=mt --set_alias="sel:mt_1<30." \
  --method=12 --var="mt_tot"["$SVSMBINS"] --cat="inclusive" \
  --x_axis_label="M_{T}^{tot} [GeV]" $MT_INC_SHIFT \
  --norm_bins=true --datacard="inclusive"\
  --background_scheme="mt_with_zmm" --extra_pad=0.2 --draw_error_band=true --auto_error_band=0.00001 --qcd_os_ss_factor=1.17

./bin/HiggsTauTauPlot5 --cfg=scripts/new_plot_"$ANA"_"$YEAR".cfg --channel=em --set_alias="sel:pzeta>-20." \
  --method=8 --var="mt_tot"["$SVSMBINS"] --cat="inclusive" \
 --x_axis_label="M_{T}^{tot} [GeV]" $EM_INC_SHIFT \
  --norm_bins=true --datacard="inclusive"\
  --background_scheme="em_default" --extra_pad=0.2 --draw_error_band=true --auto_error_band=0.00001 --qcd_os_ss_factor=2

./bin/HiggsTauTauPlot5 --cfg=scripts/new_plot_"$ANA"_"$YEAR".cfg --channel=tt --set_alias="sel:1" \
  --method=8 --var="mt_tot"["$SVSMBINS"] --cat="inclusive" \
 --x_axis_label="M_{T}^{tot} [GeV]" $EM_INC_SHIFT \
  --norm_bins=true --datacard="inclusive"\
  --background_scheme="tt_default" --extra_pad=0.2 --draw_error_band=true --auto_error_band=0.00001

#### SVFit Mass

./bin/HiggsTauTauPlot5 --cfg=scripts/new_plot_"$ANA"_"$YEAR".cfg --channel=et --set_alias="sel:mt_1<40." \
  --method=12 --var="m_sv"["$SVSMBINS"] --cat="inclusive" \
  --x_axis_label="M_{#tau#tau} [GeV]" $ET_INC_SHIFT \
  --norm_bins=true --datacard="inclusive"\
  --background_scheme="et_default" --extra_pad=0.2 --draw_error_band=true --auto_error_band=0.00001 --qcd_os_ss_factor=1.0

./bin/HiggsTauTauPlot5 --cfg=scripts/new_plot_"$ANA"_"$YEAR".cfg --channel=mt --set_alias="sel:mt_1<30." \
  --method=12  --var="m_sv"["$SVSMBINS"] --cat="inclusive" \
  --x_axis_label="M_{#tau#tau} [GeV]" $MT_INC_SHIFT \
  --norm_bins=true --datacard="inclusive"\
  --background_scheme="mt_with_zmm" --extra_pad=0.2 --draw_error_band=true --auto_error_band=0.00001 --qcd_os_ss_factor=1.17

./bin/HiggsTauTauPlot5 --cfg=scripts/new_plot_"$ANA"_"$YEAR".cfg --channel=em --set_alias="sel:pzeta>-20." \
  --method=8 --var="m_sv"["$SVSMBINS"] --cat="inclusive" \
 --x_axis_label="M_{#tau#tau} [GeV]" $EM_INC_SHIFT \
  --norm_bins=true --datacard="inclusive"\
  --background_scheme="em_default" --extra_pad=0.2 --draw_error_band=true --auto_error_band=0.00001 --qcd_os_ss_factor=2

./bin/HiggsTauTauPlot5 --cfg=scripts/new_plot_"$ANA"_"$YEAR".cfg --channel=tt --set_alias="sel:1" \
  --method=8 --var="m_sv"["$SVSMBINS"] --cat="inclusive" \
 --x_axis_label="M_{#tau#tau} [GeV]" $EM_INC_SHIFT \
  --norm_bins=true --datacard="inclusive"\
  --background_scheme="tt_default" --extra_pad=0.2 --draw_error_band=true --auto_error_band=0.00001
'''

#### SVFit Mass
'''
./bin/HiggsTauTauPlot5 --cfg=scripts/new_plot_"$ANA"_"$YEAR".cfg --channel=et --set_alias="sel:mt_1<40." \
  --method=12 --var="mt_sv"["$SVSMBINS"] --cat="inclusive" \
  --x_axis_label="M_{T,#tau#tau} [GeV]" $ET_INC_SHIFT \
  --norm_bins=true --datacard="inclusive"\
  --background_scheme="et_default" --extra_pad=0.2 #--draw_error_band=true --auto_error_band=0.00001 --qcd_os_ss_factor=1.0

./bin/HiggsTauTauPlot5 --cfg=scripts/new_plot_"$ANA"_"$YEAR".cfg --channel=mt --set_alias="sel:mt_1<30." \
  --method=12  --var="mt_sv"["$SVSMBINS"] --cat="inclusive" \
  --x_axis_label="M_{T,#tau#tau} [GeV]" $MT_INC_SHIFT \
  --norm_bins=true --datacard="inclusive"\
  --background_scheme="mt_with_zmm" --extra_pad=0.2 #--draw_error_band=true --auto_error_band=0.00001 --qcd_os_ss_factor=1.17

./bin/HiggsTauTauPlot5 --cfg=scripts/new_plot_"$ANA"_"$YEAR".cfg --channel=em --set_alias="sel:pzeta>-20." \
  --method=8 --var="mt_sv"["$SVSMBINS"] --cat="inclusive" \
 --x_axis_label="M_{T,#tau#tau} [GeV]" $EM_INC_SHIFT \
  --norm_bins=true --datacard="inclusive"\
  --background_scheme="em_default" --extra_pad=0.2 #--draw_error_band=true --auto_error_band=0.00001 --qcd_os_ss_factor=2

./bin/HiggsTauTauPlot5 --cfg=scripts/new_plot_"$ANA"_"$YEAR".cfg --channel=tt --set_alias="sel:1" \
  --method=8 --var="mt_sv"["$SVSMBINS"] --cat="inclusive" \
 --x_axis_label="M_{T,#tau#tau} [GeV]" $TT_INC_SHIFT \
  --norm_bins=true --datacard="inclusive"\
  --background_scheme="tt_default" --extra_pad=0.2 #--draw_error_band=true --auto_error_band=0.00001
'''
#### Visible Mass MSSM

#./bin/HiggsTauTauPlot5 --cfg=scripts/new_plot_"$ANA"_"$YEAR".cfg --channel=et --set_alias="sel:mt_1<30." \
#  --method=8 --var="m_vis"["$MSSMBINS"] --cat="inclusive" \
#  --x_axis_label="M_{#tau#tau}^{vis} [GeV]" $ET_INC_SHIFT \
#  --norm_bins=true --datacard="inclusive" --log_y=true --draw_ratio=true \
#  --custom_y_axis_min=true --y_axis_min=0.0099 \
#  --blind=false --x_blind_min=150 --x_blind_max=2000 \
#  --background_scheme="et_default" --extra_pad=0.3 --draw_error_band=true --auto_error_band=0.00001

#./bin/HiggsTauTauPlot5 --cfg=scripts/new_plot_"$ANA"_"$YEAR".cfg --channel=mt --set_alias="sel:mt_1<30." \
#  --method=8 --var="m_vis"["$MSSMBINS"] --cat="inclusive" \
#  --x_axis_label="M_{#tau#tau}^{vis} [GeV]" $MT_INC_SHIFT \
#  --norm_bins=true --datacard="inclusive" --log_y=true --draw_ratio=true \
#  --custom_y_axis_min=true --y_axis_min=0.0099 \
#  --blind=false --x_blind_min=150 --x_blind_max=2000 \
#  --background_scheme="mt_with_zmm" --extra_pad=0.3 --draw_error_band=true --auto_error_band=0.00001

#./bin/HiggsTauTauPlot5 --cfg=scripts/new_plot_"$ANA"_"$YEAR".cfg --channel=em --set_alias="sel:pzeta>-20." \
#  --method=8 --var="m_vis"["$MSSMBINS"] --cat="inclusive" \
#  --x_axis_label="M_{#tau#tau}^{vis} [GeV]" $EM_INC_SHIFT \
#  --norm_bins=true --datacard="inclusive" --log_y=true --draw_ratio=true \
#  --custom_y_axis_min=true --y_axis_min=0.0099 \
#  --blind=false --x_blind_min=150 --x_blind_max=2000 \
#  --background_scheme="em_default" --extra_pad=0.3 --draw_error_band=true --auto_error_band=0.00001

#./bin/HiggsTauTauPlot5 --cfg=scripts/new_plot_"$ANA"_"$YEAR".cfg --channel=tt --set_alias="sel:1" \
#  --method=8 --var="m_vis"["$MSSMBINS"] --cat="inclusive" \
#  --x_axis_label="M_{#tau#tau}^{vis} [GeV]" $TT_INC_SHIFT \
#  --norm_bins=true --datacard="inclusive" --log_y=true --draw_ratio=true \
#  --custom_y_axis_min=true --y_axis_min=0.0099 \
#  --blind=false --x_blind_min=150 --x_blind_max=2000 \
#  --background_scheme="tt_default" --extra_pad=0.3 --draw_error_band=true --auto_error_band=0.00001
#
#### MET

./bin/HiggsTauTauPlot5 --cfg=scripts/new_plot_"$ANA"_"$YEAR".cfg --channel=et --set_alias="sel:mt_1<40." \
  --method=12 --cat="inclusive" --var="met(20,0,100)" --x_axis_label="MVA E_{T}^{miss} [GeV]" \
  --norm_bins=true --datacard="inclusive" --extra_pad=0.2 \
  --background_scheme="et_default" $ET_INC_SHIFT  #--draw_error_band=true --auto_error_band=0.00001 --qcd_os_ss_factor=1.0

./bin/HiggsTauTauPlot5 --cfg=scripts/new_plot_"$ANA"_"$YEAR".cfg --channel=mt --set_alias="sel:mt_1<30." \
  --method=12 --cat="inclusive" --var="met(20,0,100)" --x_axis_label="MVA E_{T}^{miss} [GeV]" \
  --norm_bins=true --datacard="inclusive" --extra_pad=0.2 \
  --background_scheme="mt_with_zmm" $MT_INC_SHIFT #--draw_error_band=true --auto_error_band=0.00001 --qcd_os_ss_factor=1.17

./bin/HiggsTauTauPlot5 --cfg=scripts/new_plot_"$ANA"_"$YEAR".cfg --channel=em --set_alias="sel:pzeta>-20." \
  --method=15 --cat="inclusive" --var="met(20,0,100)" --x_axis_label="MVA E_{T}^{miss} [GeV]" \
  --custom_x_axis_range=true --x_axis_min=0 --x_axis_max=100 \
  --norm_bins=true --datacard="inclusive" --extra_pad=0.2 \
  --background_scheme="em_default" $EM_INC_SHIFT #--draw_error_band=true --auto_error_band=0.00001 --qcd_os_ss_factor=2

./bin/HiggsTauTauPlot5 --cfg=scripts/new_plot_"$ANA"_"$YEAR".cfg --channel=tt --set_alias="sel:1" \
  --method=8 --cat="inclusive" --var="met(20,0,100)" --x_axis_label="MVA E_{T}^{miss} [GeV]" \
  --custom_x_axis_range=true --x_axis_min=0 --x_axis_max=100 \
  --norm_bins=true --datacard="inclusive" --extra_pad=0.2 \
  --background_scheme="tt_default" $TT_INC_SHIFT #--draw_error_band=true --auto_error_band=0.00001

#./bin/HiggsTauTauPlot5 --cfg=scripts/new_plot_"$ANA"_"$YEAR".cfg --channel=et --set_alias="sel:mt_1<40." \
#  --method=8 --cat="inclusive" --var="pfmet(20,0,100)" --x_axis_label="pf E_{T}^{miss} [GeV]" \
#  --norm_bins=true --datacard="inclusive" --extra_pad=0.2 \
#  --background_scheme="et_default" $ET_INC_SHIFT  --draw_error_band=true --auto_error_band=0.00001 --qcd_os_ss_factor=1.0
#
#./bin/HiggsTauTauPlot5 --cfg=scripts/new_plot_"$ANA"_"$YEAR".cfg --channel=mt --set_alias="sel:mt_1<30." \
#  --method=8 --cat="inclusive" --var="pfmet(20,0,100)" --x_axis_label="pf E_{T}^{miss} [GeV]" \
#  --norm_bins=true --datacard="inclusive" --extra_pad=0.2 \
#  --background_scheme="mt_with_zmm" $MT_INC_SHIFT --draw_error_band=true --auto_error_band=0.00001 --qcd_os_ss_factor=1.17
#
#./bin/HiggsTauTauPlot5 --cfg=scripts/new_plot_"$ANA"_"$YEAR".cfg --channel=em --set_alias="sel:pzeta>-20." \
#  --method=8 --cat="inclusive" --var="pfmet(20,0,100)" --x_axis_label="pf E_{T}^{miss} [GeV]" \
#  --custom_x_axis_range=true --x_axis_min=0 --x_axis_max=100 \
#  --norm_bins=true --datacard="inclusive" --extra_pad=0.2 \
#  --background_scheme="em_default" $EM_INC_SHIFT --draw_error_band=true --auto_error_band=0.00001 --qcd_os_ss_factor=2
#
#./bin/HiggsTauTauPlot5 --cfg=scripts/new_plot_"$ANA"_"$YEAR".cfg --channel=tt --set_alias="sel:1" \
#  --method=8 --cat="inclusive" --var="pfmet(20,0,100)" --x_axis_label="pf E_{T}^{miss} [GeV]" \
#  --custom_x_axis_range=true --x_axis_min=0 --x_axis_max=100 \
#  --norm_bins=true --datacard="inclusive" --extra_pad=0.2 \
#  --background_scheme="tt_default" $EM_INC_SHIFT --draw_error_band=true --auto_error_band=0.00001

#### pt_1

./bin/HiggsTauTauPlot5 --cfg=scripts/new_plot_"$ANA"_"$YEAR".cfg --channel=et --set_alias="sel:mt_1<40." \
  --method=12 --cat="inclusive" --var="pt_1(25,0,100)" \
  --x_axis_label="Electron p_{T} [GeV]" --datacard="inclusive" \
  --background_scheme="et_default" $ET_INC_SHIFT #--draw_error_band=true --auto_error_band=0.00001 --qcd_os_ss_factor=1.0
  
./bin/HiggsTauTauPlot5 --cfg=scripts/new_plot_"$ANA"_"$YEAR".cfg --channel=mt --set_alias="sel:mt_1<30." \
  --method=12 --cat="inclusive" --var="pt_1(25,0,100)" \
  --x_axis_label="Muon p_{T} [GeV]" --datacard="inclusive" \
  --background_scheme="mt_with_zmm" $MT_INC_SHIFT #--draw_error_band=true --auto_error_band=0.00001 --qcd_os_ss_factor=1.17

./bin/HiggsTauTauPlot5 --cfg=scripts/new_plot_"$ANA"_"$YEAR".cfg --channel=em --set_alias="sel:pzeta>-20." \
  --method=15 --cat="inclusive" --var="pt_1(25,0,100)" \
  --x_axis_label="Electron p_{T} [GeV]" --datacard="inclusive" \
  --background_scheme="em_default" $EM_INC_SHIFT #--draw_error_band=true --auto_error_band=0.00001 --qcd_os_ss_factor=2

./bin/HiggsTauTauPlot5 --cfg=scripts/new_plot_"$ANA"_"$YEAR".cfg --channel=tt --set_alias="sel:1" \
  --method=8 --cat="inclusive" --var="pt_1(25,0,100)" \
  --x_axis_label="Leading #tau_{h} p_{T} [GeV]" --datacard="inclusive" \
  --background_scheme="tt_default" $TT_INC_SHIFT #--draw_error_band=true --auto_error_band=0.00001


#### pt_2

./bin/HiggsTauTauPlot5 --cfg=scripts/new_plot_"$ANA"_"$YEAR".cfg --channel=et --set_alias="sel:mt_1<40." \
  --method=12 --cat="inclusive" --var="pt_2(25,0,100)" \
  --x_axis_label="p_{T}(#tau_{h}) [GeV]" --datacard="inclusive" \
  --background_scheme="et_default" $ET_INC_SHIFT #--draw_error_band=true --auto_error_band=0.00001 --qcd_os_ss_factor=1.0

./bin/HiggsTauTauPlot5 --cfg=scripts/new_plot_"$ANA"_"$YEAR".cfg --channel=mt --set_alias="sel:mt_1<30." \
  --method=12 --cat="inclusive" --var="pt_2(25,0,100)" \
  --x_axis_label="p_{T}(#tau_{h}) [GeV]" --datacard="inclusive" \
  --background_scheme="mt_with_zmm" $MT_INC_SHIFT #--draw_error_band=true --auto_error_band=0.00001 --qcd_os_ss_factor=1.17

./bin/HiggsTauTauPlot5 --cfg=scripts/new_plot_"$ANA"_"$YEAR".cfg --channel=em --set_alias="sel:pzeta>-20." \
  --method=15 --cat="inclusive" --var="pt_2(25,0,100)" \
  --x_axis_label="Muon p_{T} [GeV]" --datacard="inclusive" \
  --background_scheme="em_default" $EM_INC_SHIFT #--draw_error_band=true --auto_error_band=0.00001 --qcd_os_ss_factor=2

./bin/HiggsTauTauPlot5 --cfg=scripts/new_plot_"$ANA"_"$YEAR".cfg --channel=tt --set_alias="sel:1" \
  --method=8 --cat="inclusive" --var="pt_2(25,0,100)" \
  --x_axis_label="Trailing #tau_{h} p_{T} [GeV]" --datacard="inclusive" \
  --background_scheme="tt_default" $TT_INC_SHIFT #--draw_error_band=true --auto_error_band=0.00001


#### pt_tt
./bin/HiggsTauTauPlot5 --cfg=scripts/new_plot_"$ANA"_"$YEAR".cfg --channel=et --set_alias="sel:mt_1<40." \
  --method=12 --cat="inclusive" --var="pt_tt(30,0,300)" --extra_pad=0.2 \
  --x_axis_label="p_{T}^{#tau#tau} [GeV]" --datacard="inclusive" \
  --background_scheme="et_default" $ET_INC_SHIFT #--draw_error_band=true --auto_error_band=0.00001 --qcd_os_ss_factor=1.0

./bin/HiggsTauTauPlot5 --cfg=scripts/new_plot_"$ANA"_"$YEAR".cfg --channel=mt --set_alias="sel:mt_1<30." \
  --method=12 --cat="inclusive" --var="pt_tt(30,0,300)" --extra_pad=0.2  \
  --x_axis_label="p_{T}^{#tau#tau} [GeV]" --datacard="inclusive" \
  --background_scheme="mt_with_zmm" $MT_INC_SHIFT #--draw_error_band=true --auto_error_band=0.00001 --qcd_os_ss_factor=1.17

./bin/HiggsTauTauPlot5 --cfg=scripts/new_plot_"$ANA"_"$YEAR".cfg --channel=em --set_alias="sel:pzeta>-20." \
  --method=15 --cat="inclusive" --var="pt_tt(30,0,300)" --extra_pad=0.3  \
  --x_axis_label="p_{T}^{#tau#tau} [GeV]" --datacard="inclusive" \
  --background_scheme="em_default" $EM_INC_SHIFT #--draw_error_band=true --auto_error_band=0.00001 --qcd_os_ss_factor=2

./bin/HiggsTauTauPlot5 --cfg=scripts/new_plot_"$ANA"_"$YEAR".cfg --channel=tt --set_alias="sel:1" \
  --method=8 --cat="inclusive" --var="pt_tt(30,0,300)" --extra_pad=0.3  \
  --x_axis_label="p_{T}^{#tau#tau} [GeV]" --datacard="inclusive" \
  --background_scheme="tt_default" $TT_INC_SHIFT #--draw_error_band=true --auto_error_band=0.00001

#### pt_tt log
./bin/HiggsTauTauPlot5 --cfg=scripts/new_plot_"$ANA"_"$YEAR".cfg --channel=et --set_alias="sel:mt_1<40." \
  --method=12 --cat="inclusive" --var="pt_tt(30,0,300)" --extra_pad=0.2  \
  --x_axis_label="p_{T}^{#tau#tau} [GeV]" --datacard="inclusive" --log_y=true --draw_ratio=true \
  --background_scheme="et_default" $ET_INC_SHIFT #--draw_error_band=true --auto_error_band=0.00001 --qcd_os_ss_factor=1.0

./bin/HiggsTauTauPlot5 --cfg=scripts/new_plot_"$ANA"_"$YEAR".cfg --channel=mt --set_alias="sel:mt_1<30." \
  --method=12 --cat="inclusive" --var="pt_tt(30,0,300)"  --extra_pad=0.2 \
  --x_axis_label="p_{T}^{#tau#tau} [GeV]" --datacard="inclusive" --log_y=true --draw_ratio=true \
  --background_scheme="mt_with_zmm" $MT_INC_SHIFT #--draw_error_band=true --auto_error_band=0.00001 --qcd_os_ss_factor=1.17

./bin/HiggsTauTauPlot5 --cfg=scripts/new_plot_"$ANA"_"$YEAR".cfg --channel=em --set_alias="sel:pzeta>-20." \
  --method=12 --cat="inclusive" --var="pt_tt(30,0,300)"  --extra_pad=0.3 \
  --x_axis_label="p_{T}^{#tau#tau} [GeV]" --datacard="inclusive" --log_y=true --draw_ratio=true \
  --background_scheme="em_default" $EM_INC_SHIFT #--draw_error_band=true --auto_error_band=0.00001 --qcd_os_ss_factor=2

./bin/HiggsTauTauPlot5 --cfg=scripts/new_plot_"$ANA"_"$YEAR".cfg --channel=tt --set_alias="sel:1" \
  --method=8 --cat="inclusive" --var="pt_tt(30,0,300)"  --extra_pad=0.3 \
  --x_axis_label="p_{T}^{#tau#tau} [GeV]" --datacard="inclusive" --log_y=true --draw_ratio=true \
  --background_scheme="tt_default" $TT_INC_SHIFT #--draw_error_band=true --auto_error_band=0.00001

#### eta_1

./bin/HiggsTauTauPlot5 --cfg=scripts/new_plot_"$ANA"_"$YEAR".cfg --channel=et --set_alias="sel:mt_1<40." \
  --method=12 --cat="inclusive" --var="eta_1(30,-3,3)" \
  --x_axis_label="Electron #eta" --datacard="inclusive" \
  --background_scheme="et_default" $ET_INC_SHIFT #--draw_error_band=true --auto_error_band=0.00001 --qcd_os_ss_factor=1.0

./bin/HiggsTauTauPlot5 --cfg=scripts/new_plot_"$ANA"_"$YEAR".cfg --channel=mt --set_alias="sel:mt_1<30." \
  --method=12 --cat="inclusive" --var="eta_1(30,-3,3)" \
  --x_axis_label="Muon #eta" --datacard="inclusive" \
  --background_scheme="mt_with_zmm" $MT_INC_SHIFT #--draw_error_band=true --auto_error_band=0.00001 --qcd_os_ss_factor=1.17

./bin/HiggsTauTauPlot5 --cfg=scripts/new_plot_"$ANA"_"$YEAR".cfg --channel=em --set_alias="sel:pzeta>-20."\
  --method=15 --cat="inclusive" --var="eta_1(30,-3,3)" \
  --x_axis_label="Electron #eta" --datacard="inclusive" \
  --background_scheme="em_default" $EM_INC_SHIFT #--draw_error_band=true --auto_error_band=0.00001 --qcd_os_ss_factor=2

./bin/HiggsTauTauPlot5 --cfg=scripts/new_plot_"$ANA"_"$YEAR".cfg --channel=tt --set_alias="sel:1"\
  --method=8 --cat="inclusive" --var="eta_1(30,-3,3)" \
  --x_axis_label="Leading #tau #eta" --datacard="inclusive" \
  --background_scheme="tt_default" $TT_INC_SHIFT #--draw_error_band=true --auto_error_band=0.00001


#### eta_2

./bin/HiggsTauTauPlot5 --cfg=scripts/new_plot_"$ANA"_"$YEAR".cfg --channel=et --set_alias="sel:mt_1<40." \
  --method=12 --cat="inclusive" --var="eta_2(30,-3,3)" \
  --x_axis_label="Tau #eta" --datacard="inclusive" \
  --background_scheme="et_default" $ET_INC_SHIFT #--draw_error_band=true --auto_error_band=0.00001 --qcd_os_ss_factor=1.0

./bin/HiggsTauTauPlot5 --cfg=scripts/new_plot_"$ANA"_"$YEAR".cfg --channel=mt --set_alias="sel:mt_1<30." \
  --method=12 --cat="inclusive" --var="eta_2(30,-3,3)" \
  --x_axis_label="Tau #eta" --datacard="inclusive" \
  --background_scheme="mt_with_zmm" $MT_INC_SHIFT #--draw_error_band=true --auto_error_band=0.00001 --qcd_os_ss_factor=1.17

./bin/HiggsTauTauPlot5 --cfg=scripts/new_plot_"$ANA"_"$YEAR".cfg --channel=em --set_alias="sel:pzeta>-20."\
  --method=15 --cat="inclusive" --var="eta_2(30,-3,3)" \
  --x_axis_label="Muon #eta" --datacard="inclusive" \
  --background_scheme="em_default" $EM_INC_SHIFT #--draw_error_band=true --auto_error_band=0.00001 --qcd_os_ss_factor=2

./bin/HiggsTauTauPlot5 --cfg=scripts/new_plot_"$ANA"_"$YEAR".cfg --channel=tt --set_alias="sel:1"\
  --method=8 --cat="inclusive" --var="eta_2(30,-3,3)" \
  --x_axis_label="Subleading #tau #eta" --datacard="inclusive" \
  --background_scheme="tt_default" $TT_INC_SHIFT #--draw_error_band=true --auto_error_band=0.00001

#
#### m_2
./bin/HiggsTauTauPlot5 --cfg=scripts/new_plot_"$ANA"_"$YEAR".cfg --channel=et --set_alias="sel:mt_1<40." \
  --method=12 --cat="inclusive" --var="m_2(20,0,2)" --extra_pad=0.2 \
  --x_axis_label="Tau Mass [GeV]" --datacard="inclusive" \
  --background_scheme="et_default" $ET_INC_SHIFT #--draw_error_band=true --auto_error_band=0.00001 --qcd_os_ss_factor=1.0

./bin/HiggsTauTauPlot5 --cfg=scripts/new_plot_"$ANA"_"$YEAR".cfg --channel=tt --set_alias="sel:1." \
  --method=8 --cat="inclusive" --var="m_2(20,0,2)" --extra_pad=0.2 \
  --x_axis_label="Tau Mass [GeV]" --datacard="inclusive" \
  --background_scheme="tt_default" $TT_INC_SHIFT #--draw_error_band=true --auto_error_band=0.00001

./bin/HiggsTauTauPlot5 --cfg=scripts/new_plot_"$ANA"_"$YEAR".cfg --channel=mt --set_alias="sel:mt_1<30." \
  --method=12 --cat="inclusive" --var="m_2(20,0,2)" --extra_pad=0.2\
  --x_axis_label="Tau Mass [GeV]" --datacard="inclusive" \
  --background_scheme="mt_with_zmm" $MT_INC_SHIFT #--draw_error_band=true --auto_error_band=0.00001 --qcd_os_ss_factor=1.17

./bin/HiggsTauTauPlot5 --cfg=scripts/new_plot_"$ANA"_"$YEAR".cfg --channel=et --set_alias="sel:mt_1<40." \
  --method=12 --cat="pt_2>30. && tau_decay_mode_2==1" --var="m_2(20,0,2)" --extra_pad=0.5 \
  --x_axis_label="Tau Mass [GeV]" --datacard="1prong" \
  --background_scheme="et_default" $ET_INC_SHIFT #--draw_error_band=true --auto_error_band=0.00001 --qcd_os_ss_factor=1.0

#./bin/HiggsTauTauPlot5 --cfg=scripts/new_plot_"$ANA"_"$YEAR".cfg --channel=tt --set_alias="sel:1." \
#  --method=8 --cat="pt_2>30. && tau_decay_mode_1==1" --var="m_2(20,0,2)" --extra_pad=0.5 \
#  --x_axis_label="Tau Mass [GeV]" --datacard="1prong" \
#  --background_scheme="tt_default" $ET_INC_SHIFT --draw_error_band=true --auto_error_band=0.00001

./bin/HiggsTauTauPlot5 --cfg=scripts/new_plot_"$ANA"_"$YEAR".cfg --channel=mt --set_alias="sel:mt_1<30." \
  --method=12 --cat="pt_2>30. && tau_decay_mode_2==1" --var="m_2(20,0,2)" --extra_pad=0.5\
  --x_axis_label="Tau Mass [GeV]" --datacard="1prong" \
  --background_scheme="mt_with_zmm" $MT_INC_SHIFT #--draw_error_band=true --auto_error_band=0.00001 --qcd_os_ss_factor=1.17

./bin/HiggsTauTauPlot5 --cfg=scripts/new_plot_"$ANA"_"$YEAR".cfg --channel=et --set_alias="sel:mt_1<40." \
  --method=12 --cat="pt_2>30. && tau_decay_mode_2==10" --var="m_2(20,0,2)" --extra_pad=0.5 \
  --x_axis_label="Tau Mass [GeV]" --datacard="3prong" \
  --background_scheme="et_default" $ET_INC_SHIFT #�--draw_error_band=true --auto_error_band=0.00001 --qcd_os_ss_factor=1.0

#./bin/HiggsTauTauPlot5 --cfg=scripts/new_plot_"$ANA"_"$YEAR".cfg --channel=tt --set_alias="sel:1." \
#  --method=8 --cat="pt_2>30. && tau_decay_mode==10" --var="m_2(20,0,2)" --extra_pad=0.5 \
#  --x_axis_label="Tau Mass [GeV]" --datacard="3prong" \
#  --background_scheme="tt_default" $ET_INC_SHIFT

./bin/HiggsTauTauPlot5 --cfg=scripts/new_plot_"$ANA"_"$YEAR".cfg --channel=mt --set_alias="sel:mt_1<30." \
  --method=12 --cat="pt_2>30. && tau_decay_mode_2==10" --var="m_2(20,0,2)" --extra_pad=0.5 \
  --x_axis_label="Tau Mass [GeV]" --datacard="3prong" \
  --background_scheme="mt_with_zmm" $MT_INC_SHIFT #--draw_error_band=true --auto_error_band=0.00001 --qcd_os_ss_factor=1.17

#### n_jets (log)

./bin/HiggsTauTauPlot5 --cfg=scripts/new_plot_"$ANA"_"$YEAR".cfg --channel=et --set_alias="sel:mt_1<40." \
  --method=12 --cat="inclusive" --var="n_jets(10,-0.5,9.5)"  --x_axis_label="Number of Jets" \
  --draw_ratio=true --log_y=true --extra_pad=0.3 --datacard="inclusive" \
  --background_scheme="et_default" $ET_INC_SHIFT \
  --custom_y_axis_min=true --y_axis_min=0.99 #--draw_error_band=true --auto_error_band=0.00001 --qcd_os_ss_factor=1.0

./bin/HiggsTauTauPlot5 --cfg=scripts/new_plot_"$ANA"_"$YEAR".cfg --channel=tt --set_alias="sel:1." \
  --method=8 --cat="inclusive" --var="n_jets(10,-0.5,9.5)"  --x_axis_label="Number of Jets" \
  --draw_ratio=true --log_y=true --extra_pad=0.3 --datacard="inclusive" \
  --background_scheme="tt_default" $ET_INC_SHIFT \
  --custom_y_axis_min=true --y_axis_min=0.99 #--draw_error_band=true --auto_error_band=0.00001


./bin/HiggsTauTauPlot5 --cfg=scripts/new_plot_"$ANA"_"$YEAR".cfg --channel=mt --set_alias="sel:mt_1<30." \
  --method=12 --cat="inclusive" --var="n_jets(10,-0.5,9.5)"  --x_axis_label="Number of Jets" \
  --draw_ratio=true --log_y=true --extra_pad=0.3 --datacard="inclusive" \
  --background_scheme="mt_with_zmm" $MT_INC_SHIFT \
  --custom_y_axis_min=true --y_axis_min=0.99 #--draw_error_band=true --auto_error_band=0.00001 --qcd_os_ss_factor=1.17

./bin/HiggsTauTauPlot5 --cfg=scripts/new_plot_"$ANA"_"$YEAR".cfg --channel=em --set_alias="sel:pzeta>-20." \
  --method=15 --cat="inclusive" --var="n_jets(10,-0.5,9.5)"  --x_axis_label="Number of Jets" \
  --draw_ratio=true --log_y=true --extra_pad=0.3 --datacard="inclusive" \
  --background_scheme="em_default" $EM_INC_SHIFT --qcd_os_ss_factor=2 \
  --custom_y_axis_min=true --y_axis_min=0.99 #--draw_error_band=true --auto_error_band=0.00001
#
#### n_bjets (log)

./bin/HiggsTauTauPlot5 --cfg=scripts/new_plot_"$ANA"_"$YEAR".cfg --channel=et --set_alias="sel:mt_1<40." \
  --method=12 --cat="inclusive" --var="n_bjets(5,-0.5,4.5)"  --x_axis_label="Number of b-tagged Jets" \
  --draw_ratio=true --log_y=true --extra_pad=0.3 --datacard="inclusive" \
  --background_scheme="et_default" $ET_INC_SHIFT \
  --custom_y_axis_min=true --y_axis_min=0.99 #--draw_error_band=true --auto_error_band=0.00001 --qcd_os_ss_factor=1.0

./bin/HiggsTauTauPlot5 --cfg=scripts/new_plot_"$ANA"_"$YEAR".cfg --channel=tt --set_alias="sel:1." \
  --method=8 --cat="inclusive" --var="n_bjets(5,-0.5,4.5)"  --x_axis_label="Number of b-tagged Jets (medium) " \
  --draw_ratio=true --log_y=true --extra_pad=0.3 --datacard="inclusive" \
  --background_scheme="tt_default" $TT_INC_SHIFT \
  --custom_y_axis_min=true --y_axis_min=0.99 #--draw_error_band=true --auto_error_band=0.00001


./bin/HiggsTauTauPlot5 --cfg=scripts/new_plot_"$ANA"_"$YEAR".cfg --channel=tt --set_alias="sel:1." \
  --method=8 --cat="inclusive" --var="n_loose_bjets(5,-0.5,4.5)"  --x_axis_label="Number of b-tagged Jets" \
  --draw_ratio=true --log_y=true --extra_pad=0.3 --datacard="inclusive" \
  --background_scheme="tt_default" $TT_INC_SHIFT \
  --custom_y_axis_min=true --y_axis_min=0.99 #--draw_error_band=true --auto_error_band=0.00001


./bin/HiggsTauTauPlot5 --cfg=scripts/new_plot_"$ANA"_"$YEAR".cfg --channel=mt --set_alias="sel:mt_1<30." \
  --method=12 --cat="inclusive" --var="n_bjets(5,-0.5,4.5)"  --x_axis_label="Number of b-tagged Jets" \
  --draw_ratio=true --log_y=true --extra_pad=0.3 --datacard="inclusive" \
  --background_scheme="mt_with_zmm" $MT_INC_SHIFT \
  --custom_y_axis_min=true --y_axis_min=0.99 #--draw_error_band=true --auto_error_band=0.00001 --qcd_os_ss_factor=1.17

./bin/HiggsTauTauPlot5 --cfg=scripts/new_plot_"$ANA"_"$YEAR".cfg --channel=em --set_alias="sel:pzeta>-20." \
  --method=15 --cat="inclusive" --var="n_bjets(5,-0.5,4.5)"  --x_axis_label="Number of b-tagged Jets" \
  --draw_ratio=true --log_y=true --extra_pad=0.3 --datacard="inclusive" \
  --background_scheme="em_default" $EM_INC_SHIFT --qcd_os_ss_factor=2 \
  --custom_y_axis_min=true --y_axis_min=0.99 #--draw_error_band=true --auto_error_band=0.00001
#
#### mt_1 / pzeta

./bin/HiggsTauTauPlot5 --cfg=scripts/new_plot_"$ANA"_"$YEAR".cfg --channel=et  \
  --method=12 --cat="inclusive" --set_alias="sel:1" --var="mt_1(20,0,160)" \
  --x_axis_label="m_{T} [GeV]" --datacard="inclusive" --extra_pad=0.2 $ET_BAND_ONLY \
  --background_scheme="et_default" #--draw_error_band=true --auto_error_band=0.00001 --qcd_os_ss_factor=1.0

./bin/HiggsTauTauPlot5 --cfg=scripts/new_plot_"$ANA"_"$YEAR".cfg --channel=mt \
  --method=12 --cat="inclusive" --set_alias="sel:1" --var="mt_1(20,0,160)" \
  --x_axis_label="m_{T} [GeV]" --datacard="inclusive" --extra_pad=0.2 $MT_BAND_ONLY \
  --background_scheme="mt_with_zmm" #--draw_error_band=true --auto_error_band=0.00001 --qcd_os_ss_factor=1.17
'''
./bin/HiggsTauTauPlot5 --cfg=scripts/new_plot_"$ANA"_"$YEAR".cfg --channel=et  \
  --method=12 --cat="inclusive" --set_alias="sel:1" --var="pfmt_1(20,0,160)" \
  --x_axis_label="m_{T} (pfmet) [GeV]" --datacard="inclusive" --extra_pad=0.2 $ET_BAND_ONLY \
  --background_scheme="et_default" #--draw_error_band=true --auto_error_band=0.00001

./bin/HiggsTauTauPlot5 --cfg=scripts/new_plot_"$ANA"_"$YEAR".cfg --channel=mt \
  --method=12 --cat="inclusive" --set_alias="sel:1" --var="pfmt_1(20,0,160)" \
  --x_axis_label="m_{T} (pfmet) [GeV]" --datacard="inclusive" --extra_pad=0.2 $MT_BAND_ONLY \
  --background_scheme="mt_with_zmm" #--draw_error_band=true --auto_error_band=0.00001
'''
./bin/HiggsTauTauPlot5 --cfg=scripts/new_plot_"$ANA"_"$YEAR".cfg --channel=em --set_alias="sel:1" \
  --method=15 --cat="inclusive" --var="pzeta(48,-100,100)" --extra_pad=0.3 $EM_BAND_ONLY \
  --x_axis_label="D_{#zeta} [GeV]" --datacard="inclusive" \
  --background_scheme="em_default" #--draw_error_band=true --auto_error_band=0.00001 --qcd_os_ss_factor=2
'''
./bin/HiggsTauTauPlot5 --cfg=scripts/new_plot_"$ANA"_"$YEAR".cfg --channel=em --set_alias="sel:1" \
  --method=15 --cat="inclusive" --var="pfpzeta(48,-100,100)" --extra_pad=0.3 $EM_BAND_ONLY \
  --x_axis_label="D_{#zeta} (pfmet) [GeV]" --datacard="inclusive" \
  --background_scheme="em_default" #--draw_error_band=true --auto_error_band=0.00001 --qcd_os_ss_factor=2
'''
./bin/HiggsTauTauPlot5 --cfg=scripts/new_plot_"$ANA"_"$YEAR".cfg --channel=em --set_alias="sel:1" \
  --method=15 --cat="inclusive" --var="pzetavis(25,0,50)"  $EM_BAND_ONLY \
  --x_axis_label="P_{#zeta}^{vis} [GeV]" --datacard="inclusive" \
  --background_scheme="em_default" --extra_pad=0.2 #--draw_error_band=true --auto_error_band=0.00001 --qcd_os_ss_factor=2

 ./bin/HiggsTauTauPlot5 --cfg=scripts/new_plot_"$ANA"_"$YEAR".cfg --channel=em --set_alias="sel:1" \
  --method=15 --cat="inclusive" --var="pzetamiss(25,-100,100)"  $EM_BAND_ONLY \
  --x_axis_label="#slash{P}_{#zeta} [GeV]" --datacard="inclusive" \
  --background_scheme="em_default" --extra_pad=0.3 #--draw_error_band=true --auto_error_band=0.00001 --qcd_os_ss_factor=2
#
# ./bin/HiggsTauTauPlot5 --cfg=scripts/new_plot_"$ANA"_"$YEAR".cfg --channel=em --set_alias="sel:pzeta>-20." \
#  --method=8 --cat="inclusive" --var="em_gf_mva(40,-1,1)" --extra_pad=0.2 $EM_BAND_ONLY \
#  --x_axis_label="BDT output" --datacard="inclusive" \
#  --background_scheme="em_default"
# 
# ./bin/HiggsTauTauPlot5 --cfg=scripts/new_plot_"$ANA"_"$YEAR".cfg --channel=em --set_alias="sel:pzeta>-20." \
#  --method=8 --cat="inclusive" --var="em_gf_mva(40,-1,1)" --extra_pad=0.3 $EM_BAND_ONLY \
#  --x_axis_label="BDT output" --datacard="inclusive" --log_y=true --draw_ratio=true \
#  --background_scheme="em_default"
# 
#./bin/HiggsTauTauPlot5 --cfg=scripts/new_plot_"$ANA"_"$YEAR".cfg --channel=em --set_alias="sel:pzeta>-20." \
#  --method=8 --cat="inclusive" --var="mt_ll(40,0,160)" $EM_BAND_ONLY \
#  --x_axis_label="m_{T}(ll) [GeV]" --datacard="inclusive" \
#  --background_scheme="em_default" --extra_pad=0.3
#
#./bin/HiggsTauTauPlot5 --cfg=scripts/new_plot_"$ANA"_"$YEAR".cfg --channel=em --set_alias="sel:pzeta>-20." \
#  --method=8 --cat="inclusive" --var="emu_dphi(40,0,3.15)" --extra_pad=0.3 $EM_BAND_ONLY \
#  --x_axis_label="#Delta#phi_{e,#mu}" --datacard="inclusive" \
#  --background_scheme="em_default"
#
#./bin/HiggsTauTauPlot5 --cfg=scripts/new_plot_"$ANA"_"$YEAR".cfg --channel=em --set_alias="sel:pzeta>-20." \
#  --method=8 --cat="inclusive" --var="emu_dxy_1(40,-0.02,0.02)" --extra_pad=0.3 --log_y=true --draw_ratio=true \
#  --x_axis_label="Electron d_{0}^{vtx} [cm]" --datacard="inclusive" $EM_BAND_ONLY \
#  --background_scheme="em_default"
#
#### n_vtx

./bin/HiggsTauTauPlot5 --cfg=scripts/new_plot_"$ANA"_"$YEAR".cfg --channel=et --set_alias="sel:mt_1<40." \
  --method=12 --cat="inclusive" --var="n_vtx(30,0,30)" $ET_INC_SHIFT \
  --x_axis_label="Number of Vertices" --datacard="inclusive" --extra_pad=0.4 \
  --background_scheme="et_default" #--draw_error_band=true --auto_error_band=0.00001 --qcd_os_ss_factor=1.0

./bin/HiggsTauTauPlot5 --cfg=scripts/new_plot_"$ANA"_"$YEAR".cfg --channel=tt --set_alias="sel:1." \
  --method=8 --cat="inclusive" --var="n_vtx(30,0,30)" $ET_INC_SHIFT \
  --x_axis_label="Number of Vertices" --datacard="inclusive" --extra_pad=0.4 \
  --background_scheme="tt_default" #--draw_error_band=true --auto_error_band=0.00001

./bin/HiggsTauTauPlot5 --cfg=scripts/new_plot_"$ANA"_"$YEAR".cfg --channel=mt --set_alias="sel:mt_1<30." \
  --method=12 --cat="inclusive" --var="n_vtx(30,0,30)" $MT_INC_SHIFT \
  --x_axis_label="Number of Vertices" --datacard="inclusive" --extra_pad=0.4 \
  --background_scheme="mt_with_zmm" #--draw_error_band=true --auto_error_band=0.00001 --qcd_os_ss_factor=1.17

./bin/HiggsTauTauPlot5 --cfg=scripts/new_plot_"$ANA"_"$YEAR".cfg --channel=em --set_alias="sel:pzeta>-20." \
  --method=15 --cat="inclusive" --var="n_vtx(30,0,30)" $EM_INC_SHIFT \
  --x_axis_label="Number of Vertices" --datacard="inclusive" --extra_pad=0.4 \
  --background_scheme="em_default" #--draw_error_band=true --auto_error_band=0.00001 --qcd_os_ss_factor=2

#### tau_decay_mode

#./bin/HiggsTauTauPlot5 --cfg=scripts/new_plot_"$ANA"_"$YEAR".cfg --channel=et --set_alias="sel:mt_1<30." \
#  --method=8 --cat="inclusive" --var="tau_decay_mode(4,-3,13)" $ET_INC_SHIFT \
#  --x_axis_label="Tau Decay Mode" --datacard="inclusive" --extra_pad=0.5 \
#  --background_scheme="et_default"\
#  --x_axis_bin_labels="1 Prong 0 #pi^{0}:1 Prong 1 #pi^{0}:2 Prong:3 Prong" 

#./bin/HiggsTauTauPlot5 --cfg=scripts/new_plot_"$ANA"_"$YEAR".cfg --channel=mt --set_alias="sel:mt_1<30." \
#  --method=8 --cat="inclusive" --var="tau_decay_mode(4,-3,13)" $MT_INC_SHIFT \
#  --x_axis_label="Tau Decay Mode" --datacard="inclusive" --extra_pad=0.5 \
#  --background_scheme="mt_with_zmm"\
#  --x_axis_bin_labels="1 Prong 0 #pi^{0}:1 Prong 1 #pi^{0}:2 Prong:3 Prong" 

./bin/HiggsTauTauPlot5 --cfg=scripts/new_plot_"$ANA"_"$YEAR".cfg --channel=et --set_alias="sel:mt_1<40." \
  --method=12 --cat="inclusive" --var="tau_decay_mode_2(14,0,14)" $ET_INC_SHIFT \
  --x_axis_label="Tau Decay Mode" --datacard="inclusive" --extra_pad=0.5 \
  --background_scheme="et_default" #--draw_error_band=true --auto_error_band=0.00001 --qcd_os_ss_factor=1.0

./bin/HiggsTauTauPlot5 --cfg=scripts/new_plot_"$ANA"_"$YEAR".cfg --channel=tt --set_alias="sel:1." \
  --method=8 --cat="inclusive" --var="tau_decay_mode_1(14,0,14)" $TT_INC_SHIFT \
  --x_axis_label="Leading Tau Decay Mode" --datacard="inclusive" --extra_pad=0.5 \
  --background_scheme="tt_default" #--draw_error_band=true --auto_error_band=0.00001 \


./bin/HiggsTauTauPlot5 --cfg=scripts/new_plot_"$ANA"_"$YEAR".cfg --channel=tt --set_alias="sel:1." \
  --method=8 --cat="inclusive" --var="tau_decay_mode_2(14,0,14)" $TT_INC_SHIFT \
  --x_axis_label="Subleading Tau Decay Mode" --datacard="inclusive" --extra_pad=0.5 \
  --background_scheme="tt_default" #--draw_error_band=true --auto_error_band=0.00001 \


./bin/HiggsTauTauPlot5 --cfg=scripts/new_plot_"$ANA"_"$YEAR".cfg --channel=mt --set_alias="sel:mt_1<30." \
  --method=12 --cat="inclusive" --var="tau_decay_mode_2(14,0,14)" $MT_INC_SHIFT \
  --x_axis_label="Tau Decay Mode" --datacard="inclusive" --extra_pad=0.5 \
  --background_scheme="mt_with_zmm" #--draw_error_band=true --auto_error_band=0.00001 --qcd_os_ss_factor=1.17 \

'''
#### bpt_1

./bin/HiggsTauTauPlot5 --cfg=scripts/new_plot_"$ANA"_"$YEAR".cfg --channel=et --set_alias="sel:mt_1<40." \
  --method=12 --cat="btag" --var="bpt_1(20,0,200)" $ET_BTAG_SHIFT \
  --x_axis_label="Leading b-tagged jet p_{T} [GeV]" --datacard="btag" \
  --background_scheme="et_default" --extra_pad=0.3 --draw_error_band=true --auto_error_band=0.00001 --qcd_os_ss_factor=1.0

./bin/HiggsTauTauPlot5 --cfg=scripts/new_plot_"$ANA"_"$YEAR".cfg --channel=mt --set_alias="sel:mt_1<30." \
  --method=12 --cat="btag" --var="bpt_1(20,0,200)" \
  --x_axis_label="Leading b-tagged jet p_{T} [GeV]" --datacard="btag" $MT_BTAG_SHIFT \
  --background_scheme="mt_with_zmm" --extra_pad=0.3 --draw_error_band=true --auto_error_band=0.00001 --qcd_os_ss_factor=1.17

./bin/HiggsTauTauPlot5 --cfg=scripts/new_plot_"$ANA"_"$YEAR".cfg --channel=tt --set_alias="sel:1" \
  --method=8 --cat="btag" --var="bpt_1(20,0,200)" \
  --x_axis_label="Leading b-tagged jet p_{T} [GeV]" --datacard="btag" $MT_BTAG_SHIFT \
  --background_scheme="tt_default" --extra_pad=0.3 --draw_error_band=true --auto_error_band=0.00001

./bin/HiggsTauTauPlot5 --cfg=scripts/new_plot_"$ANA"_"$YEAR".cfg --channel=em --set_alias="sel:pzeta>-20." \
  --method=8 --cat="btag" --var="bpt_1(20,0,200)" $EM_BTAG_SHIFT \
  --x_axis_label="Leading b-tagged jet p_{T} [GeV]" --datacard="btag" \
  --background_scheme="em_default" --extra_pad=0.3 --draw_error_band=true --auto_error_band=0.00001 --qcd_os_ss_factor=2
#
#### beta_1
BETA_BINS="beta_1(10,-3,3)"
#if [ "$YEAR" = "2012" ]
#then
#  BETA_BINS="beta_1(15,-3,3)"
#fi
./bin/HiggsTauTauPlot5 --cfg=scripts/new_plot_"$ANA"_"$YEAR".cfg --channel=et --set_alias="sel:mt_1<40." \
  --method=12 --cat="btag" --var=$BETA_BINS --extra_pad=0.3 $ET_BTAG_SHIFT \
  --x_axis_label="Leading b-tagged jet #eta" --datacard="btag" \
  --background_scheme="et_default" --draw_error_band=true --auto_error_band=0.00001 --qcd_os_ss_factor=1.0

./bin/HiggsTauTauPlot5 --cfg=scripts/new_plot_"$ANA"_"$YEAR".cfg --channel=tt --set_alias="sel:1" \
  --method=8 --cat="btag" --var=$BETA_BINS --extra_pad=0.3 $ET_BTAG_SHIFT \
  --x_axis_label="Leading b-tagged jet #eta" --datacard="btag" \
  --background_scheme="tt_default" --draw_error_band=true --auto_error_band=0.00001

./bin/HiggsTauTauPlot5 --cfg=scripts/new_plot_"$ANA"_"$YEAR".cfg --channel=mt --set_alias="sel:mt_1<30." \
  --method=12 --cat="btag" --var=$BETA_BINS --extra_pad=0.3 \
  --x_axis_label="Leading b-tagged jet #eta" --datacard="btag" $MT_BTAG_SHIFT \
  --background_scheme="mt_with_zmm" --draw_error_band=true --auto_error_band=0.00001 --qcd_os_ss_factor=1.17

./bin/HiggsTauTauPlot5 --cfg=scripts/new_plot_"$ANA"_"$YEAR".cfg --channel=em --set_alias="sel:pzeta>-20." \
  --method=8 --cat="btag" --var=$BETA_BINS --extra_pad=0.3 \
  --x_axis_label="Leading b-tagged jet #eta" --datacard="btag" $EM_BTAG_SHIFT \
  --background_scheme="em_default" --draw_error_band=true --auto_error_band=0.00001 --qcd_os_ss_factor=2
#
#### bcsv

./bin/HiggsTauTauPlot5 --cfg=scripts/new_plot_"$ANA"_"$YEAR".cfg --channel=et --set_alias="sel:mt_1<40. && tau_decay_mode_2!=5 && tau_decay_mode_2!=6" \
  --method=8 --cat="prebtag" --var="bcsv_1(25,0,1)" \
  --x_axis_label="CSV discriminator" --datacard="prebtag"  $ET_BTAG_SHIFT \
  --background_scheme="et_default" --extra_pad=0.2 --draw_error_band=true --auto_error_band=0.00001 --qcd_os_ss_factor=1.0

./bin/HiggsTauTauPlot5 --cfg=scripts/new_plot_"$ANA"_"$YEAR".cfg --channel=tt --set_alias="sel:tau_decay_mode_2!=5 && tau_decay_mode_2!=6" \
  --method=8 --cat="prebtag" --var="bcsv_1(25,0,1)" \
  --x_axis_label="CSV discriminator" --datacard="prebtag"  $ET_BTAG_SHIFT \
  --background_scheme="tt_default" --extra_pad=0.2 --draw_error_band=true --auto_error_band=0.00001

./bin/HiggsTauTauPlot5 --cfg=scripts/new_plot_"$ANA"_"$YEAR".cfg --channel=mt --set_alias="sel:mt_1<30. && tau_decay_mode_1!=5 && tau_decay_mode_1!=6 && tau_decay_mode_2!=5 && tau_decay_mode_2!=6" \
  --method=8 --cat="prebtag" --var="bcsv_1(25,0,1)" \
  --x_axis_label="CSV discriminator" --datacard="prebtag"  $MT_BTAG_SHIFT \
  --background_scheme="mt_with_zmm" --extra_pad=0.2 --draw_error_band=true --auto_error_band=0.00001 --qcd_os_ss_factor=1.17

./bin/HiggsTauTauPlot5 --cfg=scripts/new_plot_"$ANA"_"$YEAR".cfg --channel=em --set_alias="sel:pzeta>-20." \
  --method=8 --cat="prebtag" --var="bcsv_1(25,0,1)" \
  --x_axis_label="CSV discriminator" --datacard="prebtag"  $EM_BTAG_SHIFT \
  --background_scheme="em_default" --extra_pad=0.2 --draw_error_band=true --auto_error_band=0.00001 --qcd_os_ss_factor=2
'''
#### mjj

#./bin/HiggsTauTauPlot5 --cfg=scripts/new_plot_"$ANA"_"$YEAR".cfg --channel=et --set_alias="sel:mt_1<30." \
#  --method=8 --cat="twojet" --var="mjj(20,0,1000)" $ET_2JET_SHIFT \
#  --x_axis_label="M_{jj} [GeV]" --datacard="twojet" \
#  --background_scheme="et_default" --extra_pad=0.2

#./bin/HiggsTauTauPlot5 --cfg=scripts/new_plot_"$ANA"_"$YEAR".cfg --channel=mt --set_alias="sel:mt_1<30." \
#  --method=8 --cat="twojet" --var="mjj(20,0,1000)" $MT_2JET_SHIFT \
#  --x_axis_label="M_{jj} [GeV]" --datacard="twojet" \
#  --background_scheme="mt_with_zmm" --extra_pad=0.2

#./bin/HiggsTauTauPlot5 --cfg=scripts/new_plot_"$ANA"_"$YEAR".cfg --channel=em --set_alias="sel:pzeta>-20." \
#  --method=8 --cat="twojet" --var="mjj(20,0,1000)" $EM_2JET_SHIFT \
#  --x_axis_label="M_{jj} [GeV]" --datacard="twojet" \
#  --background_scheme="em_default" --extra_pad=0.2
#
#### jdeta

#./bin/HiggsTauTauPlot5 --cfg=scripts/new_plot_"$ANA"_"$YEAR".cfg --channel=et --set_alias="sel:mt_1<30." \
#  --method=8 --cat="twojet" --var="jdeta(20,0,10)" $ET_2JET_SHIFT \
#  --x_axis_label="#Delta#eta_{jj}" \
#  --datacard="twojet" \
#  --background_scheme="et_default" --extra_pad=0.2

#./bin/HiggsTauTauPlot5 --cfg=scripts/new_plot_"$ANA"_"$YEAR".cfg --channel=mt --set_alias="sel:mt_1<30." \
#  --method=8 --cat="twojet" --var="jdeta(20,0,10)" $MT_2JET_SHIFT \
#  --x_axis_label="#Delta#eta_{jj}" \
#  --datacard="twojet" \
#  --background_scheme="mt_with_zmm" --extra_pad=0.2

#./bin/HiggsTauTauPlot5 --cfg=scripts/new_plot_"$ANA"_"$YEAR".cfg --channel=em --set_alias="sel:pzeta>-20." \
#  --method=8 --cat="twojet" --var="jdeta(20,0,10)" $EM_2JET_SHIFT \
#  --x_axis_label="#Delta#eta_{jj}" \
#  --datacard="twojet" \
#  --background_scheme="em_default" --extra_pad=0.2

#### Z->ee / Z->mumu and W->munu plots

#./bin/ZllPlot --cfg=scripts/new_plot_zll_2015.cfg --channel=zmm --set_alias="sel:1"   --method=8 \
#   --cat="inclusive" --var="m_vis(40,60,140)"  --x_axis_label="m_{#mu#mu} [GeV]" --datacard="inclusive" \
#   --background_scheme="zmm_default" --extra_pad=0.2 --verbosity=0 

#./bin/ZllPlot --cfg=scripts/new_plot_zll_2015.cfg --channel=zmm --set_alias="sel:1"   --method=8 \
#   --cat="inclusive" --var="m_vis(40,60,140)"  --x_axis_label="m_{#mu#mu} [GeV]" --datacard="inclusive" \
#   --background_scheme="zmm_default" --extra_pad=0.2 --verbosity=0 --log_y=true

#./bin/ZllPlot --cfg=scripts/new_plot_zll_2015.cfg --channel=zee --set_alias="sel:1"   --method=8 \
#   --cat="inclusive" --var="m_vis(40,60,140)"  --x_axis_label="m_{#mu#mu} [GeV]" --datacard="inclusive" \
#   --background_scheme="zee_default" --extra_pad=0.2 --verbosity=0 

#./bin/ZllPlot --cfg=scripts/new_plot_zll_2015.cfg --channel=zee --set_alias="sel:1"   --method=8 \
#   --cat="inclusive" --var="m_vis(40,60,140)"  --x_axis_label="m_{#mu#mu} [GeV]" --datacard="inclusive" \
#   --background_scheme="zee_default" --extra_pad=0.2 --verbosity=0 --log_y=true

#./bin/ZllPlot --cfg=scripts/new_plot_zll_2015.cfg --channel=wmnu --set_alias="sel:met>40" \
#--method=8 --cat="inclusive" --var="mt_1(40,0,160)"  --x_axis_label="m_{T} [GeV]" --datacard="inclusive" \
#--background_scheme="zmm_default" --extra_pad=0.2 --verbosity=0 --folder=output/July27_NewGolden_WMuNu/
