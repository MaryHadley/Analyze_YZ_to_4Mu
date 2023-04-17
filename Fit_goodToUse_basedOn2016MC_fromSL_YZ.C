// how to run
// >> root -l Fit.C++

#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "RooAddPdf.h"
#include "RooArgusBG.h"
#include "RooBinning.h"
#include "RooCBShape.h"
#include "RooCategory.h"
#include "RooChebychev.h"
#include "RooConstVar.h"
#include "RooDLLSignificanceMCSModule.h"
#include "RooDataHist.h"
#include "RooDataSet.h"
#include "RooDecay.h"
#include "RooErrorVar.h"
#include "RooExponential.h"
#include "RooExtendPdf.h"
#include "RooFFTConvPdf.h"
#include "RooFitResult.h"
#include "RooFormulaVar.h"
#include "RooGaussModel.h"
#include "RooGaussian.h"
#include "RooGenericPdf.h"
#include "RooKeysPdf.h"
#include "RooLandau.h"
#include "RooMCStudy.h"
#include "RooNDKeysPdf.h"
#include "RooNumConvPdf.h"
#include "RooPlot.h"
#include "RooPolynomial.h"
#include "RooProdPdf.h"
#include "RooProfileLL.h"
#include "RooRealVar.h"
#include "RooSimultaneous.h"
#include "RooStats/HypoTestResult.h"
#include "RooStats/LikelihoodIntervalPlot.h"
#include "RooStats/MaxLikelihoodEstimateTestStat.h"
#include "RooStats/ModelConfig.h"
#include "RooStats/NumberCountingPdfFactory.h"
#include "RooStats/ProfileLikelihoodCalculator.h"
#include "RooStats/RatioOfProfiledLikelihoodsTestStat.h"
#include "RooStats/SPlot.h"
#include "RooVoigtian.h"
#include "RooWorkspace.h"
#include "TArrow.h"
#include "TAxis.h"
#include "TCanvas.h"
#include "TColor.h"
#include "TCut.h"
#include "TF1.h"
#include "TFile.h"
#include "TGraph.h"
#include "TH1D.h"
#include "TH3.h"
#include "TLegend.h"
#include "TMath.h"
#include "TPaveLabel.h"
#include "TRandom3.h"
#include "TStyle.h"
#include "TTree.h"

  #define DrawPulls
//#define CalculatePulls

using namespace RooFit ;
using namespace RooStats ;

void fit();
void Fit_goodToUse_basedOn2016MC_fromSL_YZ() { fit(); }
void fit() {
  
  //Constants
  double PI = 3.14159265358;
  
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  
  //Variables we fit to get yields to determine the sWeights, i.e. the discriminating variables in sPlot terminology.
  // We need to calculate these yields from the extended maximum likelihood fit in order to get sWeights
  RooRealVar upsi_mass  ("upsi_mass", "m[#mu#mu] [GeV]",  8.5, 11.); //8, 12
  RooRealVar Z_mass  ("Z_mass", "m[#mu#mu] [GeV]",  66., 116.);

 //Upsi variables
  RooRealVar upsi_pT ("upsi_pT", "upsi_pT", 0, 250); // this was originally a dummy variable, hence this comment, which I leave because it amuses me
  //although it is no longer accurate, the upsi_pT is now no longer a dummy variable and doesn't need to be updated
  //this needs to be updated of course ;) //Stefanos' comment, which made me laugh because of the winking face
  RooRealVar upsi_RAPIDITY("upsi_RAPIDITY", "upsi_RAPIDITY", -2.5, 2.5);
  RooRealVar upsi_phi("upsi_phi", "upsi_phi", -PI, PI);
  RooRealVar upsi_eta("upsi_eta", "upsi_eta", -10, 10);
  
 //Z variables
  RooRealVar Z_pT("Z_pT","Z_pT",0, 250); //Needed to raise the allowed range for the Z_pT to 1100 in order to not have any MC events excluded
  RooRealVar Z_RAPIDITY("Z_RAPIDITY", "Z_RAPIDITY", -2.5, 2.5);
  RooRealVar Z_phi("Z_phi", "Z_phi",-PI,PI);
  RooRealVar Z_eta("Z_eta", "Z_eta", -12, 12);
 
 //Lead pT mu from Z Variables
 RooRealVar lead_pT_mu_from_Z_pT("lead_pT_mu_from_Z_pT", "lead_pT_mu_from_Z_pT", 20, 200); 
 RooRealVar lead_pT_mu_from_Z_RAPIDITY("lead_pT_mu_from_Z_RAPIDITY", "lead_pT_mu_from_Z_RAPIDITY", -2.5, 2.5 );
 RooRealVar lead_pT_mu_from_Z_eta("lead_pT_mu_from_Z_eta", "lead_pT_mu_from_Z_eta",-2.5,2.5); //BE CAREFUL WITH YOUR COMMA, IF IT IS NOT OUTSIDE THE QUOTES, YOU GET NONSENSE
 RooRealVar lead_pT_mu_from_Z_phi("lead_pT_mu_from_Z_phi","lead_pT_mu_from_Z_phi", -PI,PI);
 
 //Sublead pT mu from Z Variables
 RooRealVar sublead_pT_mu_from_Z_pT("sublead_pT_mu_from_Z_pT", "sublead_pT_mu_from_Z_pT",10,150);
 RooRealVar sublead_pT_mu_from_Z_RAPIDITY("sublead_pT_mu_from_Z_RAPIDITY", "sublead_pT_mu_from_Z_RAPIDITY", -2.5, 2.5);
 RooRealVar sublead_pT_mu_from_Z_eta("sublead_pT_mu_from_Z_eta", "sublead_pT_mu_from_Z_eta", -2.5, 2.5);
 RooRealVar sublead_pT_mu_from_Z_phi("sublead_pT_mu_from_Z_phi", "sublead_pT_mu_from_Z_phi", -PI, PI);
 
 //Lead pT mu from upsi variables
 RooRealVar lead_pT_mu_from_upsi_pT ("lead_pT_mu_from_upsi_pT", "lead_pT_mu_from_upsi_pT", 0, 120);
 RooRealVar lead_pT_mu_from_upsi_RAPIDITY("lead_pT_mu_from_upsi_RAPIDITY", "lead_pT_mu_from_upsi_RAPIDITY", -2.5, 2.5);
 RooRealVar lead_pT_mu_from_upsi_eta("lead_pT_mu_from_upsi_eta", "lead_pT_mu_from_upsi_eta", -2.5, 2.5);
 RooRealVar lead_pT_mu_from_upsi_phi("lead_pT_mu_from_upsi_phi", "lead_pT_mu_from_upsi_phi", -PI, PI);
 
 //Sublead pT mu from upsi variables
 RooRealVar sublead_pT_mu_from_upsi_pT("sublead_pT_mu_from_upsi_pT", "sublead_pT_mu_from_upsi_pT", 0, 100);
 RooRealVar sublead_pT_mu_from_upsi_RAPIDITY("sublead_pT_mu_from_upsi_RAPIDITY", "sublead_pT_mu_from_upsi_RAPIDITY", -2.5, 2.5);
 RooRealVar sublead_pT_mu_from_upsi_eta("sublead_pT_mu_from_upsi_eta", "sublead_pT_mu_from_upsi_eta", -2.5, 2.5);
 RooRealVar sublead_pT_mu_from_upsi_phi("sublead_pT_mu_from_upsi_phi", "sublead_pT_mu_from_upsi_phi", -PI, PI);
 
 //Add upsi_type so that we can cut on it
  RooRealVar upsi_type("upsi_type", "upsi_type", -5, 5); //The way Stefanos suggested of trying to trick the code
 //by calling upsi_type a RooRealVar and then doing upsi_type > 0.5 and upsi_type < 1.5 if we were trying to pick
 //out the upsi_type ==1 case for example, trying the approach again now that upsi_type is a double
 
 //  RooCategory upsi_type("upsi_type", "upsi_type");
//   upsi_type.defineType("Upsi1", 1);
//   upsi_type.defineType("Upsi2", 2);
//   upsi_type.defineType("Upsi3", 3);
//   
  //Adding this as a test
//  RooRealVar big4MuVtxProb("big4MuVtxProb", "big4MuVtxProb", 0, 1);
 
 /////////////////////////////////////////////////////////////////////////////////
//  TFile *ntuple_data  = new TFile("ntuple_skimmed_maryTest_12July2021.root");
 // TFile *ntuple_data   = new TFile( "big4MuVtxCutRemoved_Runs2016-2017-2018_Total.root");
 // TFile *ntuple_data   = new TFile("ntuple_skimmed_big4muVtxCutRemoved_25Jan2022_inputFileIs_maryTest_27Jan2022_inputFileIs_MC_DPS_2016_YZ_00623223-2B20-AB42-A456-670F9B3875D5.root");
 // TFile *ntuple_data  = new TFile("ntuple_skimmed_big4muVtxCutRemoved_DimuonVtxCutsRemoved_10Feb2022_inputFileIs_maryTest_27Jan2022_inputFileIs_MC_DPS_2016_YZ_00623223-2B20-AB42-A456-670F9B3875D5.root");
//  TFile *ntuple_data = new TFile("test_Runs2016-2017-2018_Total.root");
//  TFile *ntuple_data = new TFile("ntuple_skimmed_2016-2017-2018_Total_9May2022.root");
 // TFile *ntuple_data = new TFile("testVersion_2016-2017-2018_Total.root");
  //TFile *ntuple_data = new TFile("ntuple_skimmed_inputFileIs_ZYto4Mu_Zto4Mu_pTCut3_Bjorn_MC_inputFileIs_MC_DPS_2018_YZ_04A4F969-2F02-F24D-9BA7-2FAB6D708CB6.root");
//  TFile *ntuple_data   = new TFile("bjornStyle_ntuple_Runs2016-2017-2018_Total_noRecoToTrigMuMatching.root");
//  TFile *ntuple_data   = new TFile("12July2022_skim_total_ZplusY_trigsApplied_noRecoToTrigMuMatching_JanCutsRemoved.root");
//  TFile *ntuple_data     = new TFile("12July2022_total_skim_JanCutsIncl.root");
//  TTree* tree_data    = (TTree*) ntuple_data->Get("tree");
 //   TFile *ntuple_data = new TFile("WIP_Total_ZplusY_2016_2017_2018_big4MuVtxProbCutIncluded_YMuIsoRemoved.root");
 //   TFile *ntuple_data = new TFile("2Oct2022_Runs_2016_2017_2018_Total_isoAppliedToUpsiMu_4MuVtxProbCutApplied.root");
//    TFile *ntuple_data = new TFile("2Oct2022_Runs_2016_2017_2018_Total_iso_0p85_AppliedToUpsiMu_4MuVtxProbCutApplied.root");
//    TFile *ntuple_data = new TFile("Run2Total_pfIso0p35_for_Zmu_pfIso2p0_for_upsiMu.root");
//    TFile *ntuple_data = new TFile("13Oct2022_Data_allYears_SL-like.root");
//   TFile *ntuple_data  = new TFile("2016_2017_2018_pfIso_0p35_forZmu_0p7_forUpsiMu.root");
//    TFile *ntuple_data = new TFile("ntuple_2016_2017_2018_pfIso0p7forUpsiMu_0p2forZMu.root");
//    TFile *ntuple_data = new TFile("30March2023_ntuple_2016_2017_2018_Run2_Data_Total.root"); //should be the same as 2016_2017_2018_pfIso_0p35_forZmu_0p7_forUpsiMu.root
//    TFile *ntuple_data   = new TFile("ntuple_2016_upsi_type_double.root"); //this is for testing only!
    TFile *ntuple_data = new TFile("12April2023_ntuple_v3_pfIso0p35forZmu_0p7forUpsiMu_2016_2017_2018_Total_Data.root"); // this version of the root file as upsi_type saved as a double
    //which is critical for getting it added to our RooArgSet //v3 indicates that upsi_type is now a double
//    TFile *ntuple_data = new TFile("MC_DPS_Weighted_Run2_Total_YZ_v3.root");
    TTree* tree_data = (TTree*) ntuple_data->Get("tree");
    
 // TFile *ntuple_mc    = new TFile("ntuple_skimmed_maryTest_12July2021.root");
  //TFile *ntuple_mc    = new TFile(" big4MuVtxCutRemoved_Runs2016-2017-2018_Total.root");
 // TFile *ntuple_mc    = new TFile("ntuple_skimmed_big4muVtxCutRemoved_25Jan2022_inputFileIs_maryTest_27Jan2022_inputFileIs_MC_DPS_2016_YZ_00623223-2B20-AB42-A456-670F9B3875D5.root");
 // TFile *ntuple_mc    = new TFile("ntuple_skimmed_big4muVtxCutRemoved_DimuonVtxCutsRemoved_10Feb2022_inputFileIs_maryTest_27Jan2022_inputFileIs_MC_DPS_2016_YZ_00623223-2B20-AB42-A456-670F9B3875D5.root");
 // TFile *ntuple_mc    = new TFile("test_Runs2016-2017-2018_Total.root");
 // TFile *ntuple_mc    = new TFile("ntuple_skimmed_2016-2017-2018_Total_9May2022.root");
//  TFile *ntuple_mc    = new TFile("testVersion_2016-2017-2018_Total.root");
//  TFile *ntuple_mc    = new TFile("ntuple_skimmed_inputFileIs_ZYto4Mu_Zto4Mu_pTCut3_Bjorn_MC_inputFileIs_MC_DPS_2018_YZ_04A4F969-2F02-F24D-9BA7-2FAB6D708CB6.root");
//  TFile *ntuple_mc      = new TFile("bjornStyle_ntuple_Runs2016-2017-2018_Total_noRecoToTrigMuMatching.root");
//  TFile *ntuple_mc     = new TFile("12July2022_skim_total_ZplusY_trigsApplied_noRecoToTrigMuMatching_JanCutsRemoved.root");
 // TFile *ntuple_mc    =  new TFile("12July2022_total_skim_JanCutsIncl.root");
//  TFile *ntuple_mc      = new TFile("WIP_Total_ZplusY_2016_2017_2018_big4MuVtxProbCutIncluded_YMuIsoRemoved.root");
//  TFile *ntuple_mc     = new TFile("2Oct2022_Runs_2016_2017_2018_Total_isoAppliedToUpsiMu_4MuVtxProbCutApplied.root");
//  TFile *ntuple_mc    = new TFile("2Oct2022_Runs_2016_2017_2018_Total_iso_0p85_AppliedToUpsiMu_4MuVtxProbCutApplied.root");
//  TFile *ntuple_mc = new TFile("Run2Total_pfIso0p35_for_Zmu_pfIso2p0_for_upsiMu.root");
//  TFile *ntuple_mc = new TFile("13Oct2022_Data_allYears_SL-like.root");
//  TFile *ntuple_mc  = new TFile("2016_2017_2018_pfIso_0p35_forZmu_0p7_forUpsiMu.root");
//  TFile *ntuple_mc = new TFile("ntuple_2016_2017_2018_pfIso0p7forUpsiMu_0p2forZMu.root");
//  TFile *ntuple_mc = new TFile("30March2023_ntuple_2016_2017_2018_Run2_Data_Total.root"); //should be the same as 2016_2017_2018_pfIso_0p35_forZmu_0p7_forUpsiMu.root
 // TFile *ntuple_mc    = new TFile("MC_Weighted_Run2_Total_YZ.root");
 //   TFile *ntuple_mc    = new TFile("ntuple_2016_upsi_type_double.root"); //this is for testing only!
//    TFile *ntuple_mc = new TFile("12April2023_ntuple_v3_pfIso0p35forZmu_0p7forUpsiMu_2016_2017_2018_Total_Data.root");
  TFile *ntuple_mc = new TFile("MC_DPS_Weighted_Run2_Total_YZ_v3.root"); //this version of the root file, the v3, has upsi_type as a double, which turns out to be 
  //critical for being able to add it our RooArgSet
  TTree* tree_mc      = (TTree*) ntuple_mc->Get("tree");

  RooArgSet Variables(upsi_mass, Z_mass, upsi_pT); //If you try to put too many in here, things will break, so better to do Variables.add(blah blah blah) as I do below when you want to add something
 //  RooArgSet Variables(upsi_mass, Z_mass); //upsi_pT
 /// /adding Z variables
     Variables.add(Z_pT);
     Variables.add(Z_RAPIDITY);
     Variables.add(Z_phi);
     Variables.add(Z_eta);
// // //   
// //   //adding upsi variables
     Variables.add(upsi_RAPIDITY);
     Variables.add(upsi_phi);
     Variables.add(upsi_eta);
// //   
// //   //adding lead pT mu from Z variables
     Variables.add(lead_pT_mu_from_Z_pT);
    Variables.add(lead_pT_mu_from_Z_RAPIDITY);
    Variables.add(lead_pT_mu_from_Z_eta);
    Variables.add(lead_pT_mu_from_Z_phi);
// //   
// //   //adding sublead pT mu from Z variables
     Variables.add(sublead_pT_mu_from_Z_pT);
    Variables.add(sublead_pT_mu_from_Z_RAPIDITY);
    Variables.add(sublead_pT_mu_from_Z_eta);
    Variables.add(sublead_pT_mu_from_Z_phi);
// //   
// //   //adding lead pT mu from upsi variables
     Variables.add(lead_pT_mu_from_upsi_pT);
    Variables.add(lead_pT_mu_from_upsi_RAPIDITY);
    Variables.add(lead_pT_mu_from_upsi_eta);
    Variables.add(lead_pT_mu_from_upsi_phi);
// //   
// //   //adding sublead pT mu from upsi variables
     Variables.add(sublead_pT_mu_from_upsi_pT);
    Variables.add(sublead_pT_mu_from_upsi_RAPIDITY);
    Variables.add(sublead_pT_mu_from_upsi_eta);
    Variables.add(sublead_pT_mu_from_upsi_phi);
    
    //Add upsi_type as a variable so that we can cut on it
      Variables.add(upsi_type);
// //   
//Adding this as a test
//   Variables.add(big4MuVtxProb);
  //////////////////////////////////////////////////////////////
  
  RooDataSet *data    = new RooDataSet("data", "data", tree_data, Variables);
  RooDataSet *mc      = new RooDataSet("mc",   "mc",   tree_mc,   Variables);
  
   std::cout << "Checkpoint 0" << std::endl;
   std::cout << mc->sumEntries() << std::endl;
// that's how you apply cuts if needed
//  TCut SelectionCut = "1.";
//  RooDataSet *cut_data      = (RooDataSet*)data     ->reduce(SelectionCut);

    TCut SelectionCut1 = "upsi_type > 0.5 && upsi_type < 1.5";
    TCut SelectionCut2 = "upsi_type > 1.5 && upsi_type < 2.5";
    TCut SelectionCut3 = "upsi_type > 2.5 && upsi_type < 3.5 ";
   
    RooDataSet *cut_mc1 = (RooDataSet*)mc->reduce(SelectionCut1);
    RooDataSet *cut_mc2 = (RooDataSet*)mc->reduce(SelectionCut2);
    RooDataSet *cut_mc3 = (RooDataSet*)mc->reduce(SelectionCut3);
    
    std::cout << "Checkpoint 1" << std::endl;
    std::cout << cut_mc1->sumEntries() << std::endl;
    std::cout << cut_mc2->sumEntries() << std::endl;
    std::cout << cut_mc3->sumEntries() << std::endl;

  RooRealVar mean_mUpsilon1  ("mean_mUpsilon1","mean of gaussian", 9.46);//, 9.26, 9.66); //fix mean at PDG value
  RooRealVar sigma_mUpsilon1 ("sigma_mUpsilon1","Scale Factor 1",  0.101, 0.098114, 0.103886);  //0.1004, .0853, 0.1154
  RooGaussian mSigUpsilon1   ("mSigUpsilon1","signal p.d.f.", upsi_mass, mean_mUpsilon1, sigma_mUpsilon1);

  RooRealVar mean_mUpsilon2  ("mean_mUpsilon2","mean of gaussian", 10.02);//, 9.82, 10.22); //fix mean at PDG value
  RooRealVar sigma_mUpsilon2 ("sigma_mUpsilon2","Scale Factor 1",  0.113, 0.09806, 0.12794); 
  RooGaussian mSigUpsilon2   ("mSigUpsilon2","signal p.d.f.", upsi_mass, mean_mUpsilon2, sigma_mUpsilon2);

  RooRealVar mean_mUpsilon3  ("mean_mUpsilon3","mean of gaussian", 10.36);//, 10.16, 10.56); //fix mean at PDG value
  RooRealVar sigma_mUpsilon3 ("sigma_mUpsilon3","Scale Factor 1",  0.0883, 0.08215, 0.09445);
  RooGaussian mSigUpsilon3   ("mSigUpsilon3","signal p.d.f.", upsi_mass, mean_mUpsilon3, sigma_mUpsilon3);  

  RooRealVar mean_mZ  ("mean_mZ","mean of gaussian", 91.2);
  RooRealVar width_mZ ("width_mZ","Scale Factor Z",  2.4952); //fix Z width to its PDG value
  RooRealVar sigma_mZ ("sigma_mZ","Scale Factor Z",  1.28, 1.1786, 1.3814); 
//  RooGaussian mSigZ   ("mSigZ","signal p.d.f.", Z_mass, mean_mZ, sigma_mZ);
  RooVoigtian mSigZ   ("mSigZ","signal p.d.f.", Z_mass, mean_mZ, width_mZ, sigma_mZ); //RooVoigtian is a convolution of BW and Gaussian //see: https://root.cern/doc/master/classRooVoigtian.html //width needs to come first, followed by sigma 

  RooRealVar        cUpsilon("cUpsilon", "cUpsilon", -0.267, -.5412, 0.0072); // -1.04, -2.54, 0// -.5, -20, 0
  RooExponential mBkgUpsilon("mBkgUpsilon", "exponential", upsi_mass, cUpsilon);

  RooRealVar        cZ("cZ", "cZ", -0.0476, -0.07592, -0.01928); // -.043, -.0902, 0 //-.5, -20, 0
  RooExponential mBkgZ("mBkgZ", "exponential", Z_mass, cZ);

  RooProdPdf Upsi1_S_Z_B ("Upsi1_S_Z_B", "Upsi1_S_Z_B", RooArgSet(mSigUpsilon1, mBkgZ));
  RooProdPdf Upsi2_S_Z_B ("Upsi2_S_Z_B", "Upsi2_S_Z_B", RooArgSet(mSigUpsilon2, mBkgZ));
  RooProdPdf Upsi3_S_Z_B ("Upsi3_S_Z_B", "Upsi3_S_Z_B", RooArgSet(mSigUpsilon3, mBkgZ));
  RooProdPdf UpsiX_B_Z_B ("UpsiX_B_Z_B", "UpsiX_B_Z_B", RooArgSet(mBkgUpsilon,  mBkgZ));
  RooProdPdf UpsiX_B_Z_S ("UpsiX_B_Z_S", "UpsiX_B_Z_S", RooArgSet(mBkgUpsilon,  mSigZ));
  RooProdPdf Upsi1_S_Z_S ("Upsi1_S_Z_S", "Upsi1_S_Z_S", RooArgSet(mSigUpsilon1, mSigZ));
  RooProdPdf Upsi2_S_Z_S ("Upsi2_S_Z_S", "Upsi2_S_Z_S", RooArgSet(mSigUpsilon2, mSigZ));
  RooProdPdf Upsi3_S_Z_S ("Upsi3_S_Z_S", "Upsi3_S_Z_S", RooArgSet(mSigUpsilon3, mSigZ));
 
 
 //  RooRealVar N_Upsi1_S_Z_B ("N_Upsi1_S_Z_B", "N_Upsi1_S_Z_B", 0); //10., 0., 10000.);
//   RooRealVar N_Upsi2_S_Z_B ("N_Upsi2_S_Z_B", "N_Upsi2_S_Z_B", 0); //10., 0., 10000.);
//   RooRealVar N_Upsi3_S_Z_B ("N_Upsi3_S_Z_B", "N_Upsi3_S_Z_B", 0); //10., 0., 10000.);
//   RooRealVar N_UpsiX_B_Z_B ("N_UpsiX_B_Z_B", "N_UpsiX_B_Z_B", 0); //10., 0., 10000.);
//   RooRealVar N_UpsiX_B_Z_S ("N_UpsiX_B_Z_S", "N_UpsiX_B_Z_S", 0); //10., 0., 10000.);
//   
  //Let these background variables above float freely instead of having them fixed to 0 
  RooRealVar N_Upsi1_S_Z_B ("N_Upsi1_S_Z_B", "N_Upsi1_S_Z_B", 10., 0., 10000.);
  RooRealVar N_Upsi2_S_Z_B ("N_Upsi2_S_Z_B", "N_Upsi2_S_Z_B", 10., 0., 10000.);
  RooRealVar N_Upsi3_S_Z_B ("N_Upsi3_S_Z_B", "N_Upsi3_S_Z_B", 10., 0., 10000.);
  RooRealVar N_UpsiX_B_Z_B ("N_UpsiX_B_Z_B", "N_UpsiX_B_Z_B", 10., 0., 10000.);
  RooRealVar N_UpsiX_B_Z_S ("N_UpsiX_B_Z_S", "N_UpsiX_B_Z_S", 10., 0., 10000.);
  RooRealVar N_Upsi1_S_Z_S ("N_Upsi1_S_Z_S", "N_Upsi1_S_Z_S", 100., 0., 1000.);
  RooRealVar N_Upsi2_S_Z_S ("N_Upsi2_S_Z_S", "N_Upsi2_S_Z_S", 50.,  0., 1000.);
  RooRealVar N_Upsi3_S_Z_S ("N_Upsi3_S_Z_S", "N_Upsi3_S_Z_S", 20.,  0., 1000.);

  RooExtendPdf e_Upsi1_S_Z_B ("e_Upsi1_S_Z_B", "e_Upsi1_S_Z_B", Upsi1_S_Z_B, N_Upsi1_S_Z_B);
  RooExtendPdf e_Upsi2_S_Z_B ("e_Upsi2_S_Z_B", "e_Upsi2_S_Z_B", Upsi2_S_Z_B, N_Upsi2_S_Z_B);
  RooExtendPdf e_Upsi3_S_Z_B ("e_Upsi3_S_Z_B", "e_Upsi3_S_Z_B", Upsi3_S_Z_B, N_Upsi3_S_Z_B);
  RooExtendPdf e_UpsiX_B_Z_B ("e_UpsiX_B_Z_B", "e_UpsiX_B_Z_B", UpsiX_B_Z_B, N_UpsiX_B_Z_B);
  RooExtendPdf e_UpsiX_B_Z_S ("e_UpsiX_B_Z_S", "e_UpsiX_B_Z_S", UpsiX_B_Z_S, N_UpsiX_B_Z_S);
  RooExtendPdf e_Upsi1_S_Z_S ("e_Upsi1_S_Z_S", "e_Upsi1_S_Z_S", Upsi1_S_Z_S, N_Upsi1_S_Z_S);
  RooExtendPdf e_Upsi2_S_Z_S ("e_Upsi2_S_Z_S", "e_Upsi2_S_Z_S", Upsi2_S_Z_S, N_Upsi2_S_Z_S);
  RooExtendPdf e_Upsi3_S_Z_S ("e_Upsi3_S_Z_S", "e_Upsi3_S_Z_S", Upsi3_S_Z_S, N_Upsi3_S_Z_S);

  RooAddPdf eSum ("eSum", "eSum", RooArgList(e_Upsi1_S_Z_B, e_Upsi2_S_Z_B, e_Upsi3_S_Z_B, e_UpsiX_B_Z_B, e_UpsiX_B_Z_S, e_Upsi1_S_Z_S, e_Upsi2_S_Z_S, e_Upsi3_S_Z_S));

  // fit
  RooFitResult *fr = eSum.fitTo(*data, NumCPU(4, kTRUE), Save(), Extended());
  
  //Upsi (1S) + Z case
  double nll_0sig, nll_obs, p0_nosyst;
  RooAbsReal* nll = eSum.createNLL(*data,NumCPU(4)) ; //Only need to do this step once!

  double yields_NOniaSigZSig =N_Upsi1_S_Z_S.getVal();
  N_Upsi1_S_Z_S.setVal(0);                   nll_0sig = nll->getVal();
  N_Upsi1_S_Z_S.setVal(yields_NOniaSigZSig); nll_obs  = nll->getVal();
  p0_nosyst = TMath::Prob( 2.*(nll_0sig-nll_obs), 1) / 2.;
  cout << "Significance for Y(1S) + Z case: " << p0_nosyst << " " << RooStats::PValueToSignificance(p0_nosyst) << endl;

//Upsi(2S) + Z case
double nll_0sig_2S, nll_obs_2S, p0_nosyst_2S;
double yields_NOniaSigZSig_2S = N_Upsi2_S_Z_S.getVal();
N_Upsi2_S_Z_S.setVal(0); nll_0sig_2S = nll->getVal();
N_Upsi2_S_Z_S.setVal(yields_NOniaSigZSig_2S); nll_obs_2S = nll->getVal();
p0_nosyst_2S = TMath::Prob( 2.*(nll_0sig_2S-nll_obs_2S), 1) / 2.;
std::cout << "Significance for Y(2S) + Z case: " << p0_nosyst_2S << " " << RooStats::PValueToSignificance(p0_nosyst_2S) << std::endl;

//Upsi(3S) + Z case
double nll_0sig_3S, nll_obs_3S, p0_nosyst_3S;
double yields_NOniaSigZSig_3S = N_Upsi3_S_Z_S.getVal();
N_Upsi3_S_Z_S.setVal(0); nll_0sig_3S = nll->getVal();
N_Upsi3_S_Z_S.setVal(yields_NOniaSigZSig_3S); nll_obs_3S = nll->getVal();
p0_nosyst_3S = TMath::Prob( 2.*(nll_0sig_3S-nll_obs_3S), 1) / 2.;
std::cout << "Significance for Y(3S) + Z case: " << p0_nosyst_3S << " " << RooStats::PValueToSignificance(p0_nosyst_3S) << std::endl;



  
//   #ifdef CalculatePulls
//   int Ntoys = 500;
//   RooMCStudy *SimToy = new RooMCStudy(eSum, RooArgSet(upsi_mass,Z_mass), Binned(kFALSE), Silence(), Extended(), FitOptions(Save(kTRUE)), PrintEvalErrors(0), Minos(kTRUE));
// 
//   SimToy->generateAndFit(Ntoys);
// 
//   RooPlot *N_Upsi1_S_Z_S_pull_frame = SimToy->plotPull(N_Upsi1_S_Z_S, -4, 4, 8, kFALSE);
//   RooPlot *N_Upsi2_S_Z_S_pull_frame = SimToy->plotPull(N_Upsi2_S_Z_S, -4, 4, 8, kFALSE);
//   RooPlot *N_Upsi3_S_Z_S_pull_frame = SimToy->plotPull(N_Upsi3_S_Z_S, -4, 4, 8, kFALSE);
// 
//   TCanvas* cPulls = new TCanvas("cPulls","toymc", 1200, 400); cPulls->Divide(3,1);
//   cPulls->cd(1); N_Upsi1_S_Z_S_pull_frame->Draw(); N_Upsi1_S_Z_S_pull_frame->SetXTitle("Z+Y(1S) signal yield");
//   cPulls->cd(2); N_Upsi2_S_Z_S_pull_frame->Draw(); N_Upsi2_S_Z_S_pull_frame->SetXTitle("Z+Y(2S) signal yield");
//   cPulls->cd(3); N_Upsi3_S_Z_S_pull_frame->Draw(); N_Upsi3_S_Z_S_pull_frame->SetXTitle("Z+Y(3S) signal yield");
//   cPulls->SaveAs("c_pulls.pdf");
// #endif

  // plot part
  int Upsi_bins = 25; //20
  int Z_bins = 25; //20

  RooPlot *frame_main_fit1 = upsi_mass.frame(Title("mass 1 fit"), Bins(Upsi_bins));
  data->plotOn(frame_main_fit1, XErrorSize(0), Name("plotdata"));
  eSum.plotOn(frame_main_fit1, LineColor(kGray+2), Name("totalpdf"));
  eSum.plotOn(frame_main_fit1, Components(e_Upsi1_S_Z_S), LineStyle(kDashed), LineColor(kBlue), Name("e_Upsi1_S_Z_S"));
  eSum.plotOn(frame_main_fit1, Components(e_Upsi2_S_Z_S), LineStyle(kDashed), LineColor(kRed), Name("e_Upsi2_S_Z_S"));
  eSum.plotOn(frame_main_fit1, Components(e_Upsi3_S_Z_S), LineStyle(kDashed), LineColor(kOrange), Name("e_Upsi3_S_Z_S"));
  
  eSum.plotOn(frame_main_fit1, Components(e_Upsi1_S_Z_B), LineStyle(kDashed), LineColor(kMagenta), Name ("e_Upsi1_S_Z_B"));
  eSum.plotOn(frame_main_fit1, Components(e_Upsi2_S_Z_B), LineStyle(kDashed), LineColor(kViolet), Name("e_Upsi2_S_Z_B"));
  eSum.plotOn(frame_main_fit1, Components(e_Upsi3_S_Z_B), LineStyle(kDashed), LineColor(kCyan), Name("e_Upsi3_S_Z_B"));
  
  eSum.plotOn(frame_main_fit1, Components(e_UpsiX_B_Z_B), LineStyle(kDashed), LineColor(kYellow), Name("e_UpsiX_B_Z_B"));
  eSum.plotOn(frame_main_fit1, Components(e_UpsiX_B_Z_S), LineStyle(kDashed), LineColor(kGreen), Name("e_UpsiX_B_Z_S"));
  
  data->plotOn(frame_main_fit1, XErrorSize(0));//redraw the data in case it got covered by the other curves
  TCanvas *c_mass_1 = new TCanvas("c_mass_1", "c_mass_1", 900, 900); c_mass_1->cd();
  
  #ifdef DrawPulls
  RooPlot* dummy_frame_jpsi = upsi_mass.frame(Title("dummy frame to extract pulls"), Bins(Upsi_bins));
  data->plotOn(dummy_frame_jpsi, XErrorSize(0)); eSum.plotOn(dummy_frame_jpsi);

  RooHist* h_pulls_mass_jpsi = dummy_frame_jpsi->pullHist();
  RooPlot* frame_pulls_mass_jpsi = upsi_mass.frame(Title("Pull Distribution #mu^{+}#mu^{#font[122]{\55}} mass"));
  frame_pulls_mass_jpsi->addPlotable(h_pulls_mass_jpsi, "P");

  TPad *pad1_jpsi = new TPad("pad1_jpsi", "The pad 80% of the height",0.0,0.2,1.0,1.0,21);
  TPad *pad2_jpsi = new TPad("pad2_jpsi", "The pad 20% of the height",0.0,0.0,1.0,0.2,22);
  pad1_jpsi->SetFillColor(0); pad2_jpsi->SetFillColor(0);
  pad1_jpsi->Draw();
  pad2_jpsi->Draw();

  pad2_jpsi->cd();
  frame_pulls_mass_jpsi->Draw();
  TLine *l1 = new TLine(8.5, 0, 11, 0); l1->Draw("");
  pad1_jpsi->cd();
#endif
  
  frame_main_fit1->Draw();
  TLegend *leg1 = new TLegend (0.5, 0.6, 0.8,0.85);
  leg1->AddEntry("plotdata", "Data", "P");
  leg1->AddEntry("totalpdf", "Fit", "LP");
  leg1->AddEntry("e_Upsi1_S_Z_S", "Y (1S) + Z Signal", "LP");
  leg1->AddEntry("e_Upsi2_S_Z_S", "Y (2S) + Z Signal", "LP");
  leg1->AddEntry("e_Upsi3_S_Z_S", "Y (3S) + Z Signal", "LP");
  leg1->AddEntry("e_Upsi1_S_Z_B", "Y(1S) Sig, Z Bkg", "LP");
  leg1->AddEntry("e_Upsi2_S_Z_B", "Y(2S) Sig, Z Bkg", "LP");
  leg1->AddEntry("e_Upsi3_S_Z_B", "Y(3S) Sig, Z Bkg", "LP");
  leg1->AddEntry("e_UpsiX_B_Z_B", "Y(nS) Bkg, Z Bkg", "LP");
  leg1->AddEntry("e_UpsiX_B_Z_S", "Y(nS) Bkg, Z Sig", "LP");
  leg1->Draw();
  c_mass_1->SaveAs("c_mass_1.pdf");

  RooPlot *frame_main_fit2 = Z_mass.frame(Title("mass 2 fit"), Bins(Z_bins));
  data->plotOn(frame_main_fit2, XErrorSize(0), Name("plotdata"));
  eSum.plotOn(frame_main_fit2, LineColor(kGray+2), Name("Fit_Z"));
  
  eSum.plotOn(frame_main_fit2, Components(e_Upsi1_S_Z_S), LineStyle(kDashed), LineColor(kBlue), Name("e_Upsi1_S_Z_S"));
  eSum.plotOn(frame_main_fit2, Components(e_Upsi2_S_Z_S), LineStyle(kDashed), LineColor(kRed), Name("e_Upsi2_S_Z_S"));
  eSum.plotOn(frame_main_fit2, Components(e_Upsi3_S_Z_S), LineStyle(kDashed), LineColor(kOrange), Name("e_Upsi3_S_Z_S"));
  eSum.plotOn(frame_main_fit2, Components(e_UpsiX_B_Z_S), LineStyle(kDashed), LineColor(kGreen), Name("e_UpsiX_B_Z_S"));
  
  eSum.plotOn(frame_main_fit2, Components(e_Upsi1_S_Z_B), LineStyle(kDashed), LineColor(kMagenta), Name ("e_Upsi1_S_Z_B"));
  eSum.plotOn(frame_main_fit2, Components(e_Upsi2_S_Z_B), LineStyle(kDashed), LineColor(kViolet), Name("e_Upsi2_S_Z_B"));
  eSum.plotOn(frame_main_fit2, Components(e_Upsi3_S_Z_B), LineStyle(kDashed), LineColor(kCyan), Name("e_Upsi3_S_Z_B"));
  
  eSum.plotOn(frame_main_fit2, Components(e_UpsiX_B_Z_B), LineStyle(kDashed), LineColor(kYellow), Name("e_UpsiX_B_Z_B"));
  
  data->plotOn(frame_main_fit2, XErrorSize(0)); //redraw the data in case it got covered by the other curves
  TCanvas *c_mass_2 = new TCanvas("c_mass_2", "c_mass_2", 900, 900); c_mass_2->cd();
  
  #ifdef DrawPulls
   RooPlot* dummy_frame_z = Z_mass.frame(Title("dummy frame to extract pulls"), Bins(Upsi_bins));
  data->plotOn(dummy_frame_z, XErrorSize(0)); eSum.plotOn(dummy_frame_z);

  RooHist* h_pulls_mass_z = dummy_frame_z->pullHist();
  RooPlot* frame_pulls_mass_z = Z_mass.frame(Title("Pull Distribution #mu^{+}#mu^{#font[122]{\55}} mass"));
  frame_pulls_mass_z->addPlotable(h_pulls_mass_z, "P");

  TPad *pad1_z = new TPad("pad1_z", "The pad 80% of the height",0.0,0.2,1.0,1.0,21);
  TPad *pad2_z = new TPad("pad2_z", "The pad 20% of the height",0.0,0.0,1.0,0.2,22);
  pad1_z->SetFillColor(0); pad2_z->SetFillColor(0);
  pad1_z->Draw();
  pad2_z->Draw();

  pad2_z->cd();
  frame_pulls_mass_z->Draw();
  TLine *l2 = new TLine(66, 0., 116, 0.); l2->Draw("");
  pad1_z->cd();
#endif

  frame_main_fit2->Draw();
  TLegend *leg2 = new TLegend(0.6, 0.6, 0.8,0.85);
  leg2->AddEntry("plotdata", "Data", "P");
  leg2->AddEntry("Fit_Z", "Fit", "LP");
  leg2->AddEntry("e_Upsi1_S_Z_S", "Y (1S) + Z Signal", "LP");
  leg2->AddEntry("e_Upsi2_S_Z_S", "Y (2S) + Z Signal", "LP");
  leg2->AddEntry("e_Upsi3_S_Z_S", "Y (3S) + Z Signal", "LP");
  leg2->AddEntry("e_Upsi1_S_Z_B", "Y(1S) Sig, Z Bkg", "LP");
  leg2->AddEntry("e_Upsi2_S_Z_B", "Y(2S) Sig, Z Bkg", "LP");
  leg2->AddEntry("e_Upsi3_S_Z_B", "Y(3S) Sig, Z Bkg", "LP");
  leg2->AddEntry("e_UpsiX_B_Z_B", "Y(nS) Bkg, Z Bkg", "LP");
  leg2->AddEntry("e_UpsiX_B_Z_S", "Y (nS) Bkg, Z Sig", "LP");
  leg2->Draw();
  
  c_mass_2->SaveAs("c_mass_2.pdf");



  // splot part

  RooWorkspace* ws = new RooWorkspace("myWS");
  ws->import(eSum);
  ws->import(*data, Rename("sPlotdata"));

  RooAbsPdf* sPloteSum = ws->pdf("eSum");

  RooRealVar* sPlot_N_Upsi1_S_Z_B = ws->var("N_Upsi1_S_Z_B");
  RooRealVar* sPlot_N_Upsi2_S_Z_B = ws->var("N_Upsi2_S_Z_B");
  RooRealVar* sPlot_N_Upsi3_S_Z_B = ws->var("N_Upsi3_S_Z_B");
  RooRealVar* sPlot_N_UpsiX_B_Z_B = ws->var("N_UpsiX_B_Z_B");
  RooRealVar* sPlot_N_UpsiX_B_Z_S = ws->var("N_UpsiX_B_Z_S");
  RooRealVar* sPlot_N_Upsi1_S_Z_S = ws->var("N_Upsi1_S_Z_S");
  RooRealVar* sPlot_N_Upsi2_S_Z_S = ws->var("N_Upsi2_S_Z_S");
  RooRealVar* sPlot_N_Upsi3_S_Z_S = ws->var("N_Upsi3_S_Z_S");

  sPlot_N_Upsi1_S_Z_B->setConstant(); 
  sPlot_N_Upsi2_S_Z_B->setConstant();
  sPlot_N_Upsi3_S_Z_B->setConstant();
  sPlot_N_UpsiX_B_Z_B->setConstant();
  sPlot_N_UpsiX_B_Z_S->setConstant();
  sPlot_N_Upsi1_S_Z_S->setConstant();
  sPlot_N_Upsi2_S_Z_S->setConstant();
  sPlot_N_Upsi3_S_Z_S->setConstant();

  RooStats::SPlot* sData = new RooStats::SPlot("sData","An SPlot", *data, sPloteSum, RooArgList(*sPlot_N_Upsi1_S_Z_B, *sPlot_N_Upsi2_S_Z_B, *sPlot_N_Upsi3_S_Z_B, *sPlot_N_UpsiX_B_Z_B,
                                                                                                *sPlot_N_UpsiX_B_Z_S, *sPlot_N_Upsi1_S_Z_S, *sPlot_N_Upsi2_S_Z_S, *sPlot_N_Upsi3_S_Z_S));
//Make the new weighted datasets for all variables 
  RooDataSet * data_weighted_1 = new RooDataSet(data->GetName(), data->GetTitle(), data, *data->get(), 0, "N_Upsi1_S_Z_S_sw"); // this selected the Z+Y(1S). Pick the right variable and add _sw
  RooDataSet * data_weighted_2 = new RooDataSet(data->GetName(), data->GetTitle(), data, *data->get(), 0, "N_Upsi2_S_Z_S_sw"); 
  RooDataSet * data_weighted_3 = new RooDataSet(data->GetName(), data->GetTitle(), data, *data->get(), 0, "N_Upsi3_S_Z_S_sw"); 
  
  //Attempt to read out the sWeights 
  
   for(Int_t i=0; i < data->numEntries(); i++) {
    const RooArgSet*  obs = data->get(i);
    double weight = sData->GetSWeight(i,"N_Upsi1_S_Z_S_sw");
    std::cout << "sWeight for event i = " << i << "," << " related to how likely this event is to be Upsi1_S_Z_S:  " << weight << std::endl;
    
//    std::cout << "///////////////////////////////////////////" << std::endl;
    
    double weight2 = sData->GetSWeight(i, "N_Upsi2_S_Z_S_sw");
    std::cout << "sWeight for event i = " << i << "," << " related to how likely this event is to be Upsi2_S_Z_S:  " << weight2 << std::endl;
    
//    std::cout << "///////////////////////////////////////////" << std::endl; 
    
    double weight3 = sData->GetSWeight(i, "N_Upsi3_S_Z_S_sw");
    std::cout << "sWeight for event i = " << i << "," << " related to how likely this event is to be Upsi3_S_Z_S:  " << weight3 << std::endl;
    
    std::cout << "///////////////////////////////////////////" << std::endl; 
}
 
 int upsi_pT_bins = 5;
 int upsi_RAP_bins = 5;
 int upsi_phi_bins = 5;
 int upsi_eta_bins = 8;
 
 int Z_pT_bins = 5;
 int Z_RAP_bins = 5;
 int Z_phi_bins = 5;
 int Z_eta_bins = 8;
 
 int lead_pT_mu_from_Z_pT_bins = 5;
 int lead_pT_mu_from_Z_RAP_bins = 5;
 int lead_pT_mu_from_Z_eta_bins = 5;
 int lead_pT_mu_from_Z_phi_bins = 5;
 
 int sublead_pT_mu_from_Z_pT_bins = 5;
 int sublead_pT_mu_from_Z_RAP_bins = 5;
 int sublead_pT_mu_from_Z_eta_bins = 5;
 int sublead_pT_mu_from_Z_phi_bins = 5;
 
 int lead_pT_mu_from_upsi_pT_bins = 5;
 int lead_pT_mu_from_upsi_RAP_bins = 5;
 int lead_pT_mu_from_upsi_eta_bins = 5;
 int lead_pT_mu_from_upsi_phi_bins = 5;
 
 int sublead_pT_mu_from_upsi_pT_bins = 5;
 int sublead_pT_mu_from_upsi_RAP_bins = 5;
 int sublead_pT_mu_from_upsi_eta_bins = 5;
 int sublead_pT_mu_from_upsi_phi_bins = 5;

 //Make plots of the sWeighted distributions for each variable 
  RooPlot *frame_pt_upsilon1 = upsi_pT.frame(Title("pT of upsi(1S) produced in association with Z"), Bins(upsi_pT_bins));
  data_weighted_1->plotOn(frame_pt_upsilon1, XErrorSize(0));

  RooPlot *frame_pt_upsilon2 = upsi_pT.frame(Title("pT of upsi(2S) produced in association with Z"), Bins(upsi_pT_bins));
  data_weighted_2->plotOn(frame_pt_upsilon2, XErrorSize(0));

  RooPlot *frame_pt_upsilon3 = upsi_pT.frame(Title("pT of upsi(3S) produced in association with Z"), Bins(upsi_pT_bins));
  data_weighted_3->plotOn(frame_pt_upsilon3, XErrorSize(0));
  
  RooPlot *frame_RAP_upsi1 = upsi_RAPIDITY.frame(Title("Rapidity of upsi(1S) produced in association with Z"), Bins(upsi_RAP_bins));
  data_weighted_1->plotOn(frame_RAP_upsi1, XErrorSize(0));
  
  RooPlot *frame_RAP_upsi2 = upsi_RAPIDITY.frame(Title("Rapidity of upsi(2S) produced in association with Z"), Bins(upsi_RAP_bins));
  data_weighted_2->plotOn(frame_RAP_upsi2, XErrorSize(0));
  
  RooPlot *frame_RAP_upsi3 = upsi_RAPIDITY.frame(Title("Rapidity of upsi(3S) produced in association with Z"), Bins(upsi_RAP_bins));
  data_weighted_3->plotOn(frame_RAP_upsi3, XErrorSize(0));
  
  RooPlot *frame_phi_upsi1 = upsi_phi.frame(Title("Phi of upsi(1S) produced in association with Z"),Bins(upsi_phi_bins) );
  data_weighted_1->plotOn(frame_phi_upsi1, XErrorSize(0)); 
  
  RooPlot *frame_phi_upsi2 = upsi_phi.frame(Title("Phi of upsi(2S) produced in association with Z"),Bins(upsi_phi_bins) );
  data_weighted_2->plotOn(frame_phi_upsi2, XErrorSize(0)); 
  
  RooPlot *frame_phi_upsi3 = upsi_phi.frame(Title("Phi of upsi(3S) produced in association with Z"),Bins(upsi_phi_bins) );
  data_weighted_3->plotOn(frame_phi_upsi3, XErrorSize(0)); 
  
  RooPlot *frame_eta_upsi1 = upsi_eta.frame(Title("Eta of upsi(1S) produced in association with Z"), Bins(upsi_eta_bins));
  data_weighted_1->plotOn(frame_eta_upsi1, XErrorSize(0));
  
  RooPlot *frame_eta_upsi2 = upsi_eta.frame(Title("Eta of upsi(2S) produced in association with Z"), Bins(upsi_eta_bins));
  data_weighted_2->plotOn(frame_eta_upsi2, XErrorSize(0));
  
  RooPlot *frame_eta_upsi3 = upsi_eta.frame(Title("Eta of upsi(3S) produced in association with Z"), Bins(upsi_eta_bins));
  data_weighted_3->plotOn(frame_eta_upsi3, XErrorSize(0));
  
  RooPlot *frame_pt_Z_upsi1 = Z_pT.frame(Title("pT of Z produced in association with upsi(1S)"), Bins(Z_pT_bins));
  data_weighted_1->plotOn(frame_pt_Z_upsi1, XErrorSize(0));
  
  RooPlot *frame_pt_Z_upsi2 = Z_pT.frame(Title("pT of Z produced in association with upsi(2S)"), Bins(Z_pT_bins));
  data_weighted_2->plotOn(frame_pt_Z_upsi2, XErrorSize(0));
  
  RooPlot *frame_pt_Z_upsi3 = Z_pT.frame(Title("pT of Z produced in association with upsi(3S)"), Bins(Z_pT_bins));
  data_weighted_3->plotOn(frame_pt_Z_upsi3, XErrorSize(0));
  
  RooPlot *frame_Z_RAP_upsi1 = Z_RAPIDITY.frame(Title("Rapidity of Z produced in association with upsi(1S)"), Bins(Z_RAP_bins));
  data_weighted_1->plotOn(frame_Z_RAP_upsi1, XErrorSize(0));
  
  RooPlot *frame_Z_RAP_upsi2 = Z_RAPIDITY.frame(Title("Rapidity of Z produced in association with upsi(2S)"), Bins(Z_RAP_bins));
  data_weighted_2->plotOn(frame_Z_RAP_upsi2, XErrorSize(0));
  
  RooPlot *frame_Z_RAP_upsi3 = Z_RAPIDITY.frame(Title("Rapidity of Z produced in association with upsi(3S)"), Bins(Z_RAP_bins));
  data_weighted_3->plotOn(frame_Z_RAP_upsi3, XErrorSize(0));
  
  RooPlot *frame_Z_phi_upsi1 = Z_phi.frame(Title("Phi of Z produced in association with upsi(1S)"), Bins(Z_phi_bins));
  data_weighted_1->plotOn(frame_Z_phi_upsi1, XErrorSize(0));
  
  RooPlot *frame_Z_phi_upsi2 = Z_phi.frame(Title("Phi of Z produced in association with upsi(2S)"), Bins(Z_phi_bins));
  data_weighted_2->plotOn(frame_Z_phi_upsi2, XErrorSize(0));
  
  RooPlot *frame_Z_phi_upsi3 = Z_phi.frame(Title("Phi of Z produced in association with upsi(3S)"), Bins(Z_phi_bins));
  data_weighted_3->plotOn(frame_Z_phi_upsi3, XErrorSize(0));
  
  RooPlot *frame_Z_eta_upsi1 = Z_eta.frame(Title("Eta of Z produced in association with upsi(1S)"), Bins(Z_eta_bins));
  data_weighted_1->plotOn(frame_Z_eta_upsi1, XErrorSize(0));
  
  RooPlot *frame_Z_eta_upsi2 = Z_eta.frame(Title("Eta of Z produced in association with upsi(2S)"), Bins(Z_eta_bins));
  data_weighted_2->plotOn(frame_Z_eta_upsi2, XErrorSize(0));
  
  RooPlot *frame_Z_eta_upsi3 = Z_eta.frame(Title("Eta of Z produced in association with upsi(3S)"), Bins(Z_eta_bins));
  data_weighted_3->plotOn(frame_Z_eta_upsi3, XErrorSize(0));
  
  RooPlot *frame_lead_pT_mu_from_Z_pT_upsi1 = lead_pT_mu_from_Z_pT.frame(Title("pT of lead pT mu from Z produced in association with upsi(1S)"), Bins(lead_pT_mu_from_Z_pT_bins));
  data_weighted_1->plotOn(frame_lead_pT_mu_from_Z_pT_upsi1, XErrorSize(0));
  
  RooPlot *frame_lead_pT_mu_from_Z_pT_upsi2 = lead_pT_mu_from_Z_pT.frame(Title("pT of lead pT mu from Z produced in association with upsi(2S)"), Bins(lead_pT_mu_from_Z_pT_bins));
  data_weighted_2->plotOn(frame_lead_pT_mu_from_Z_pT_upsi2, XErrorSize(0));
  
  RooPlot *frame_lead_pT_mu_from_Z_pT_upsi3 = lead_pT_mu_from_Z_pT.frame(Title("pT of lead pT mu from Z produced in association with upsi(3S)"), Bins(lead_pT_mu_from_Z_pT_bins));
  data_weighted_3->plotOn(frame_lead_pT_mu_from_Z_pT_upsi3, XErrorSize(0));
  
  RooPlot *frame_lead_pT_mu_from_Z_RAP_upsi1 = lead_pT_mu_from_Z_RAPIDITY.frame(Title("Rapidity of lead pT mu from Z produced in association with upsi(1S)"), Bins(lead_pT_mu_from_Z_RAP_bins));
  data_weighted_1->plotOn(frame_lead_pT_mu_from_Z_RAP_upsi1, XErrorSize(0));
  
  RooPlot *frame_lead_pT_mu_from_Z_RAP_upsi2 = lead_pT_mu_from_Z_RAPIDITY.frame(Title("Rapidity of lead pT mu from Z produced in association with upsi(2S)"), Bins(lead_pT_mu_from_Z_RAP_bins));
  data_weighted_2->plotOn(frame_lead_pT_mu_from_Z_RAP_upsi2, XErrorSize(0));
  
  RooPlot *frame_lead_pT_mu_from_Z_RAP_upsi3 = lead_pT_mu_from_Z_RAPIDITY.frame(Title("Rapidity of lead pT mu from Z produced in association with upsi(3S)"), Bins(lead_pT_mu_from_Z_RAP_bins));
  data_weighted_3->plotOn(frame_lead_pT_mu_from_Z_RAP_upsi3, XErrorSize(0));
  
  RooPlot *frame_lead_pT_mu_from_Z_eta_upsi1 = lead_pT_mu_from_Z_eta.frame(Title("Eta of lead pT mu from Z produced in association with upsi(1S)"), Bins(lead_pT_mu_from_Z_eta_bins));
  data_weighted_1->plotOn(frame_lead_pT_mu_from_Z_eta_upsi1, XErrorSize(0));
   
  RooPlot *frame_lead_pT_mu_from_Z_eta_upsi2 = lead_pT_mu_from_Z_eta.frame(Title("Eta of lead pT mu from Z produced in association with upsi(2S)"), Bins(lead_pT_mu_from_Z_eta_bins));
  data_weighted_2->plotOn(frame_lead_pT_mu_from_Z_eta_upsi2, XErrorSize(0));
  
  RooPlot *frame_lead_pT_mu_from_Z_eta_upsi3 = lead_pT_mu_from_Z_eta.frame(Title("Eta of lead pT mu from Z produced in association with upsi(3S)"), Bins(lead_pT_mu_from_Z_eta_bins));
  data_weighted_3->plotOn(frame_lead_pT_mu_from_Z_eta_upsi3, XErrorSize(0));
  
  RooPlot *frame_lead_pT_mu_from_Z_phi_upsi1 = lead_pT_mu_from_Z_phi.frame(Title("Phi of lead pT mu from Z produced in association with upsi(1S)"),Bins(lead_pT_mu_from_Z_phi_bins));
  data_weighted_1->plotOn(frame_lead_pT_mu_from_Z_phi_upsi1, XErrorSize(0));   
  
  RooPlot *frame_lead_pT_mu_from_Z_phi_upsi2 = lead_pT_mu_from_Z_phi.frame(Title("Phi of lead pT mu from Z produced in association with upsi(2S)"),Bins(lead_pT_mu_from_Z_phi_bins));
  data_weighted_2->plotOn(frame_lead_pT_mu_from_Z_phi_upsi2, XErrorSize(0));
  
  RooPlot *frame_lead_pT_mu_from_Z_phi_upsi3 = lead_pT_mu_from_Z_phi.frame(Title("Phi of lead pT mu from Z produced in association with upsi(3S)"),Bins(lead_pT_mu_from_Z_phi_bins));
  data_weighted_3->plotOn(frame_lead_pT_mu_from_Z_phi_upsi3, XErrorSize(0));
  
  RooPlot *frame_sublead_pT_mu_from_Z_pT_upsi1 = sublead_pT_mu_from_Z_pT.frame(Title("pT of sublead pT mu from Z produced in association with upsi(1S)"), Bins(sublead_pT_mu_from_Z_pT_bins));
  data_weighted_1->plotOn(frame_sublead_pT_mu_from_Z_pT_upsi1, XErrorSize(0));
  
  RooPlot *frame_sublead_pT_mu_from_Z_pT_upsi2 = sublead_pT_mu_from_Z_pT.frame(Title("pT of sublead pT mu from Z produced in association with upsi(2S)"), Bins(sublead_pT_mu_from_Z_pT_bins));
  data_weighted_2->plotOn(frame_sublead_pT_mu_from_Z_pT_upsi2, XErrorSize(0));
  
  RooPlot *frame_sublead_pT_mu_from_Z_pT_upsi3 = sublead_pT_mu_from_Z_pT.frame(Title("pT of sublead pT mu from Z produced in association with upsi(3S)"), Bins(sublead_pT_mu_from_Z_pT_bins));
  data_weighted_3->plotOn(frame_sublead_pT_mu_from_Z_pT_upsi3, XErrorSize(0));
  
  RooPlot *frame_sublead_pT_mu_from_Z_RAP_upsi1 = sublead_pT_mu_from_Z_RAPIDITY.frame(Title("Rapidity of sublead pT mu from Z produced in association with uspi(1S)"), Bins(sublead_pT_mu_from_Z_RAP_bins));
  data_weighted_1->plotOn(frame_sublead_pT_mu_from_Z_RAP_upsi1, XErrorSize(0));
  
  RooPlot *frame_sublead_pT_mu_from_Z_RAP_upsi2 = sublead_pT_mu_from_Z_RAPIDITY.frame(Title("Rapidity of sublead pT mu from Z produced in association with uspi(2S)"), Bins(sublead_pT_mu_from_Z_RAP_bins));
  data_weighted_2->plotOn(frame_sublead_pT_mu_from_Z_RAP_upsi2, XErrorSize(0));
  
  RooPlot *frame_sublead_pT_mu_from_Z_RAP_upsi3 = sublead_pT_mu_from_Z_RAPIDITY.frame(Title("Rapidity of sublead pT mu from Z produced in association with uspi(3S)"), Bins(sublead_pT_mu_from_Z_RAP_bins));
  data_weighted_3->plotOn(frame_sublead_pT_mu_from_Z_RAP_upsi3, XErrorSize(0));
  
  RooPlot *frame_sublead_pT_mu_from_Z_eta_upsi1 = sublead_pT_mu_from_Z_eta.frame(Title("Eta of sublead pT mu from Z produced in association with upsi(1S)"), Bins(sublead_pT_mu_from_Z_eta_bins));
  data_weighted_1->plotOn(frame_sublead_pT_mu_from_Z_eta_upsi1, XErrorSize(0));
  
  RooPlot *frame_sublead_pT_mu_from_Z_eta_upsi2 = sublead_pT_mu_from_Z_eta.frame(Title("Eta of sublead pT mu from Z produced in association with upsi(2S)"), Bins(sublead_pT_mu_from_Z_eta_bins));
  data_weighted_2->plotOn(frame_sublead_pT_mu_from_Z_eta_upsi2, XErrorSize(0));
  
  RooPlot *frame_sublead_pT_mu_from_Z_eta_upsi3 = sublead_pT_mu_from_Z_eta.frame(Title("Eta of sublead pT mu from Z produced in association with upsi(3S)"), Bins(sublead_pT_mu_from_Z_eta_bins));
  data_weighted_3->plotOn(frame_sublead_pT_mu_from_Z_eta_upsi3, XErrorSize(0));
  
  RooPlot *frame_sublead_pT_mu_from_Z_phi_upsi1 = sublead_pT_mu_from_Z_phi.frame(Title("Phi of sublead pT mu from Z produced in association with upsi(1S)"), Bins(sublead_pT_mu_from_Z_phi_bins));
  data_weighted_1->plotOn(frame_sublead_pT_mu_from_Z_phi_upsi1, XErrorSize(0));
  
  RooPlot *frame_sublead_pT_mu_from_Z_phi_upsi2 = sublead_pT_mu_from_Z_phi.frame(Title("Phi of sublead pT mu from Z produced in association with upsi(2S)"), Bins(sublead_pT_mu_from_Z_phi_bins));
  data_weighted_2->plotOn(frame_sublead_pT_mu_from_Z_phi_upsi2, XErrorSize(0));
  
  RooPlot *frame_sublead_pT_mu_from_Z_phi_upsi3 = sublead_pT_mu_from_Z_phi.frame(Title("Phi of sublead pT mu from Z produced in association with upsi(3S)"), Bins(sublead_pT_mu_from_Z_phi_bins));
  data_weighted_3->plotOn(frame_sublead_pT_mu_from_Z_phi_upsi3, XErrorSize(0));
  
  RooPlot *frame_lead_pT_mu_from_upsi_pT_upsi1 = lead_pT_mu_from_upsi_pT.frame(Title("pT of lead pT mu from upsi(1S) produced in association with Z"), Bins(lead_pT_mu_from_upsi_pT_bins));
  data_weighted_1->plotOn(frame_lead_pT_mu_from_upsi_pT_upsi1, XErrorSize(0));
  
  RooPlot *frame_lead_pT_mu_from_upsi_pT_upsi2 = lead_pT_mu_from_upsi_pT.frame(Title("pT of lead pT mu from upsi(2S) produced in association with Z"), Bins(lead_pT_mu_from_upsi_pT_bins));
  data_weighted_2->plotOn(frame_lead_pT_mu_from_upsi_pT_upsi2, XErrorSize(0));
  
  RooPlot *frame_lead_pT_mu_from_upsi_pT_upsi3 = lead_pT_mu_from_upsi_pT.frame(Title("pT of lead pT mu from upsi(3S) produced in association with Z"), Bins(lead_pT_mu_from_upsi_pT_bins));
  data_weighted_3->plotOn(frame_lead_pT_mu_from_upsi_pT_upsi3, XErrorSize(0));
  
  RooPlot *frame_lead_pT_mu_from_upsi_RAP_upsi1 = lead_pT_mu_from_upsi_RAPIDITY.frame(Title("Rapidity of lead pT mu from upsi(1S) produced in association with Z"), Bins(lead_pT_mu_from_upsi_RAP_bins));
  data_weighted_1->plotOn(frame_lead_pT_mu_from_upsi_RAP_upsi1, XErrorSize(0));
  
  RooPlot *frame_lead_pT_mu_from_upsi_RAP_upsi2 = lead_pT_mu_from_upsi_RAPIDITY.frame(Title("Rapidity of lead pT mu from upsi(2S) produced in association with Z"), Bins(lead_pT_mu_from_upsi_RAP_bins));
  data_weighted_2->plotOn(frame_lead_pT_mu_from_upsi_RAP_upsi2, XErrorSize(0));
  
  RooPlot *frame_lead_pT_mu_from_upsi_RAP_upsi3 = lead_pT_mu_from_upsi_RAPIDITY.frame(Title("Rapidity of lead pT mu from upsi(3S) produced in association with Z"), Bins(lead_pT_mu_from_upsi_RAP_bins));
  data_weighted_3->plotOn(frame_lead_pT_mu_from_upsi_RAP_upsi3, XErrorSize(0));
  
  RooPlot *frame_lead_pT_mu_from_upsi_eta_upsi1 = lead_pT_mu_from_upsi_eta.frame(Title("Eta of lead pT mu from upsi(1S) produced in association with Z"), Bins(lead_pT_mu_from_upsi_eta_bins));
  data_weighted_1->plotOn(frame_lead_pT_mu_from_upsi_eta_upsi1, XErrorSize(0));
  
  RooPlot *frame_lead_pT_mu_from_upsi_eta_upsi2 = lead_pT_mu_from_upsi_eta.frame(Title("Eta of lead pT mu from upsi(2S) produced in association with Z"), Bins(lead_pT_mu_from_upsi_eta_bins));
  data_weighted_2->plotOn(frame_lead_pT_mu_from_upsi_eta_upsi2, XErrorSize(0));
  
  RooPlot *frame_lead_pT_mu_from_upsi_eta_upsi3 = lead_pT_mu_from_upsi_eta.frame(Title("Eta of lead pT mu from upsi(3S) produced in association with Z"), Bins(lead_pT_mu_from_upsi_eta_bins));
  data_weighted_3->plotOn(frame_lead_pT_mu_from_upsi_eta_upsi3, XErrorSize(0));
  
  RooPlot *frame_lead_pT_mu_from_upsi_phi_upsi1 = lead_pT_mu_from_upsi_phi.frame(Title("Phi of lead pT mu from upsi(1S) produced in association wtih Z"), Bins(lead_pT_mu_from_upsi_phi_bins));
  data_weighted_1->plotOn(frame_lead_pT_mu_from_upsi_phi_upsi1, XErrorSize(0));
  
  RooPlot *frame_lead_pT_mu_from_upsi_phi_upsi2 = lead_pT_mu_from_upsi_phi.frame(Title("Phi of lead pT mu from upsi(2S) produced in association wtih Z"), Bins(lead_pT_mu_from_upsi_phi_bins));
  data_weighted_2->plotOn(frame_lead_pT_mu_from_upsi_phi_upsi2, XErrorSize(0));
  
  RooPlot *frame_lead_pT_mu_from_upsi_phi_upsi3 = lead_pT_mu_from_upsi_phi.frame(Title("Phi of lead pT mu from upsi(3S) produced in association wtih Z"), Bins(lead_pT_mu_from_upsi_phi_bins));
  data_weighted_3->plotOn(frame_lead_pT_mu_from_upsi_phi_upsi3, XErrorSize(0));
  
  RooPlot *frame_sublead_pT_mu_from_upsi_pT_upsi1 = sublead_pT_mu_from_upsi_pT.frame(Title("pT of sublead pT mu from upsi(1S) produced in association with Z"), Bins(sublead_pT_mu_from_upsi_pT_bins));
  data_weighted_1->plotOn(frame_sublead_pT_mu_from_upsi_pT_upsi1, XErrorSize(0));
  
  RooPlot *frame_sublead_pT_mu_from_upsi_pT_upsi2 = sublead_pT_mu_from_upsi_pT.frame(Title("pT of sublead pT mu from upsi(2S) produced in association with Z"), Bins(sublead_pT_mu_from_upsi_pT_bins));
  data_weighted_2->plotOn(frame_sublead_pT_mu_from_upsi_pT_upsi2, XErrorSize(0));
  
  RooPlot *frame_sublead_pT_mu_from_upsi_pT_upsi3 = sublead_pT_mu_from_upsi_pT.frame(Title("pT of sublead pT mu from upsi(3S) produced in association with Z"), Bins(sublead_pT_mu_from_upsi_pT_bins));
  data_weighted_3->plotOn(frame_sublead_pT_mu_from_upsi_pT_upsi3, XErrorSize(0));
  
  RooPlot *frame_sublead_pT_mu_from_upsi_RAP_upsi1 = sublead_pT_mu_from_upsi_RAPIDITY.frame(Title("Rapidity of sublead pT mu from upsi(1S) produced in association with Z"), Bins(sublead_pT_mu_from_upsi_RAP_bins));
  data_weighted_1->plotOn(frame_sublead_pT_mu_from_upsi_RAP_upsi1, XErrorSize(0));
  
  RooPlot *frame_sublead_pT_mu_from_upsi_RAP_upsi2 = sublead_pT_mu_from_upsi_RAPIDITY.frame(Title("Rapidity of sublead pT mu from upsi(2S) produced in association with Z"), Bins(sublead_pT_mu_from_upsi_RAP_bins));
  data_weighted_2->plotOn(frame_sublead_pT_mu_from_upsi_RAP_upsi2, XErrorSize(0));
  
  RooPlot *frame_sublead_pT_mu_from_upsi_RAP_upsi3 = sublead_pT_mu_from_upsi_RAPIDITY.frame(Title("Rapidity of sublead pT mu from upsi(3S) produced in association with Z"), Bins(sublead_pT_mu_from_upsi_RAP_bins));
  data_weighted_3->plotOn(frame_sublead_pT_mu_from_upsi_RAP_upsi3, XErrorSize(0));
  
  RooPlot *frame_sublead_pT_mu_from_upsi_eta_upsi1 = sublead_pT_mu_from_upsi_eta.frame(Title("Eta of sublead pT mu from upsi(1S) produced in association with Z"), Bins(sublead_pT_mu_from_upsi_eta_bins));
  data_weighted_1->plotOn(frame_sublead_pT_mu_from_upsi_eta_upsi1, XErrorSize(0));
  
  RooPlot *frame_sublead_pT_mu_from_upsi_eta_upsi2 = sublead_pT_mu_from_upsi_eta.frame(Title("Eta of sublead pT mu from upsi(2S) produced in association with Z"), Bins(sublead_pT_mu_from_upsi_eta_bins));
  data_weighted_2->plotOn(frame_sublead_pT_mu_from_upsi_eta_upsi2, XErrorSize(0));
  
  RooPlot *frame_sublead_pT_mu_from_upsi_eta_upsi3 = sublead_pT_mu_from_upsi_eta.frame(Title("Eta of sublead pT mu from upsi(3S) produced in association with Z"), Bins(sublead_pT_mu_from_upsi_eta_bins));
  data_weighted_3->plotOn(frame_sublead_pT_mu_from_upsi_eta_upsi3, XErrorSize(0));
  
  RooPlot *frame_sublead_pT_mu_from_upsi_phi_upsi1 = sublead_pT_mu_from_upsi_phi.frame(Title("Phi of sublead pT mu from uspi(1S) produced in association with Z"), Bins(sublead_pT_mu_from_upsi_phi_bins));
  data_weighted_1->plotOn(frame_sublead_pT_mu_from_upsi_phi_upsi1, XErrorSize(0));
  
  RooPlot *frame_sublead_pT_mu_from_upsi_phi_upsi2 = sublead_pT_mu_from_upsi_phi.frame(Title("Phi of sublead pT mu from uspi(2S) produced in association with Z"), Bins(sublead_pT_mu_from_upsi_phi_bins));
  data_weighted_2->plotOn(frame_sublead_pT_mu_from_upsi_phi_upsi2, XErrorSize(0));
  
  RooPlot *frame_sublead_pT_mu_from_upsi_phi_upsi3 = sublead_pT_mu_from_upsi_phi.frame(Title("Phi of sublead pT mu from uspi(3S) produced in association with Z"), Bins(sublead_pT_mu_from_upsi_phi_bins));
  data_weighted_3->plotOn(frame_sublead_pT_mu_from_upsi_phi_upsi3, XErrorSize(0));
  
  //Won't cause a bug, nothing will be drawn because I added the lines in the wrong place, they need to go as shown in the example of c_weighted 
 
//TH1* tmp6 = mc->createHistogram("Z_mass",Z_bins); tmp6->Scale(N_Upsi1_S_Z_S.getVal() / mc->sumEntries()); tmp6->Draw("hesame"); tmp6->SetLineColor(kRed); tmp6->SetMarkerSize(0);
  
 //Canvases, draw on them and save them
 
 //upper and lower limits for Y axis of ratio plots defined here
 double lowerLimitYForRatioPlot = 0.;
 double upperLimitYForRatioPlot = 5.;
 
 //Canvas 1
  TCanvas *c_weighted = new TCanvas("c_weighted", "c_weighted", 1200, 400); c_weighted->Divide(3,1);
  
  c_weighted->cd(1); 
  TPad *pad1 = new TPad("pad1", "The pad 80% of the height",0.0,0.2,1.0,1.0,21);
  TPad *pad2 = new TPad("pad2", "The pad 20% of the height",0.0,0.0,1.0,0.2,22);
  pad1->SetFillColor(0); pad2->SetFillColor(0);
  pad1->Draw();
  pad2->Draw();
  
  pad1->cd();
  std::cout << "CHECKPOINT A" << std::endl;
  frame_pt_upsilon1->Draw(); //Draw the sWeighted data
  TH1* tmp1 = cut_mc1->createHistogram("upsi_pT", upsi_pT_bins); //Draw and scale the upsi pT of upsilons that have been tagged as upsi 1
  tmp1->Scale(N_Upsi1_S_Z_S.getVal()/cut_mc1->sumEntries());
  tmp1->Draw("hesame");
  tmp1->SetLineColor(kRed);
  tmp1->SetMarkerSize(0);
  c_weighted->cd(1);
  
  pad2->cd();
  std::cout << "CHECKPOINT A1" << std::endl;
 // frame_pt_upsilon1->Sumw2(); //Doesn't work, there is no member Sumw2 in RooPlot
  std::cout << "CHECKPOINT B" << std::endl;
 // TH1* h1 = (TH1*)frame_pt_upsilon1->Clone("h1"); //this doesn't cause a compilation error, but creates a silent monster that then causes seg faults when you try
 //to do things with it, such as call Sumw2 (you get a seg fault) or Divide (another seg fault). Looks like you cannot convert a RooPlot object to a TH1 in this way
 //Unfortunately, ROOT doesn't help you out by giving a useful error when you try to do this
  TH1* h1 = data_weighted_1->createHistogram("upsi_pT", upsi_pT_bins);
  std::cout << "CHECKPOINT C" << std::endl;
  std::cout << h1->GetBinContent(1) << std::endl;
  std::cout << h1->GetBinContent(2) << std::endl;
  std::cout << tmp1->GetBinContent(1) << std::endl;
  std::cout << tmp1->GetBinContent(2) << std::endl;
  h1->Sumw2();
  std::cout << "CHECKPOINT D"  << std::endl;
  h1->Divide(tmp1);
  std::cout << "CHECKPOINT E" << std::endl;
  h1->Draw();
  h1->GetYaxis()->SetRangeUser(lowerLimitYForRatioPlot,upperLimitYForRatioPlot)  ;
  c_weighted->cd(1);
  
  c_weighted->cd(2); 
  TPad *pad3 = new TPad("pad3", "The pad 80% of the height",0.0,0.2,1.0,1.0,21);
  TPad *pad4 = new TPad("pad4", "The pad 20% of the height",0.0,0.0,1.0,0.2,22);
  pad3->SetFillColor(0); pad4->SetFillColor(0);
  pad3->Draw();
  pad4->Draw();
  
  pad3->cd();
  frame_pt_upsilon2->Draw();
  TH1* tmp2 = cut_mc2->createHistogram("upsi_pT", upsi_pT_bins);
  tmp2->Scale(N_Upsi2_S_Z_S.getVal()/cut_mc2->sumEntries());
  tmp2->Draw("hesame");
  tmp2->SetLineColor(kRed);
  tmp2->SetMarkerSize(0);
  c_weighted->cd(2);
  
  pad4->cd();
  TH1* h2 = data_weighted_2->createHistogram("upsi_pT", upsi_pT_bins);
  h2->Sumw2();
  h2->Divide(tmp2);
  h2->Draw();
  h2->GetYaxis()->SetRangeUser(lowerLimitYForRatioPlot,upperLimitYForRatioPlot)  ;
  
  c_weighted->cd(3); 
  TPad *pad5 = new TPad("pad5", "The pad 80% of the height",0.0,0.2,1.0,1.0,21);
  TPad *pad6 = new TPad("pad6", "The pad 20% of the height",0.0,0.0,1.0,0.2,22);
  pad5->SetFillColor(0); pad6->SetFillColor(0);
  pad5->Draw();
  pad6->Draw();
  
  pad5->cd();
  frame_pt_upsilon3->Draw();
  TH1* tmp3= cut_mc3->createHistogram("upsi_pT", upsi_pT_bins);
  tmp3->Scale(N_Upsi3_S_Z_S.getVal()/cut_mc3->sumEntries());
  tmp3->Draw("hesame");
  tmp3->SetLineColor(kRed);
  tmp3->SetMarkerSize(0);
  c_weighted->cd(3);
  
  pad6->cd();
  TH1* h3 = data_weighted_3->createHistogram("upsi_pT", upsi_pT_bins);
  h3->Sumw2();
  std::cout << "CHECKPOINT F" << std::endl;
  std::cout << h3->GetBinContent(2) << std::endl;
  std::cout << tmp3->GetBinContent(2) << std::endl;
  h3->Divide(tmp3);
  h3->Draw();
  h3->GetYaxis()->SetRangeUser(lowerLimitYForRatioPlot,upperLimitYForRatioPlot)  ;
  
  c_weighted->SaveAs("c_weighted_upsi_pT.pdf");

  
//Canvas 2
  TCanvas *c_weighted2 = new TCanvas("c_weighted2", "c_weighted2", 1200,400); c_weighted2->Divide(3,1);
  
  c_weighted2->cd(1); 
  TPad *pad7 = new TPad("pad7", "The pad 80% of the height",0.0,0.2,1.0,1.0,21);
  TPad *pad8 = new TPad("pad8", "The pad 20% of the height",0.0,0.0,1.0,0.2,22);
  pad7->SetFillColor(0); pad8->SetFillColor(0);
  pad7->Draw();
  pad8->Draw();
  
  pad7->cd();
  frame_pt_Z_upsi1->Draw();
  TH1* tmp4 = cut_mc1->createHistogram("Z_pT", Z_pT_bins);
  tmp4->Scale(N_Upsi1_S_Z_S.getVal()/cut_mc1->sumEntries());
  tmp4->Draw("hesame");
  tmp4->SetLineColor(kRed);
  tmp4->SetMarkerSize(0);
  c_weighted2->cd(1);
  
  pad8->cd();
  TH1* h4 = data_weighted_1->createHistogram("Z_pT", Z_pT_bins);
  h4->Sumw2();
  std::cout << "CHECKPOINT G" << std::endl;
  std::cout << h4->GetBinContent(2) << std::endl;
  std::cout << tmp4->GetBinContent(2) << std::endl;
  std::cout << h4->GetBinContent(3) << std::endl;
  std::cout << tmp4->GetBinContent(3) << std::endl;
  h4->Divide(tmp4);
  h4->Draw();
  h4->GetYaxis()->SetRangeUser(lowerLimitYForRatioPlot,upperLimitYForRatioPlot)  ;
  
  c_weighted2->cd(2); 
  TPad *pad9 = new TPad("pad9", "The pad 80% of the height",0.0,0.2,1.0,1.0,21);
  TPad *pad10 = new TPad("pad10", "The pad 20% of the height",0.0,0.0,1.0,0.2,22);
  pad9->SetFillColor(0); pad10->SetFillColor(0);
  pad9->Draw();
  pad10->Draw();
  
  pad9->cd();
  frame_pt_Z_upsi2->Draw();
  TH1* tmp5 = cut_mc2->createHistogram("Z_pT", Z_pT_bins);
  tmp5->Scale(N_Upsi2_S_Z_S.getVal()/cut_mc2->sumEntries());
  tmp5->Draw("hesame");
  tmp5->SetLineColor(kRed);
  tmp5->SetMarkerSize(0);
  c_weighted2->cd(2);
  
  pad10->cd();
  TH1* h5 = data_weighted_2->createHistogram("Z_pT", Z_pT_bins);
  h5->Sumw2();
  h5->Divide(tmp5);
  h5->Draw();
  h5->GetYaxis()->SetRangeUser(lowerLimitYForRatioPlot,upperLimitYForRatioPlot)  ;
  
  c_weighted2->cd(3); 
  TPad *pad11 = new TPad("pad11", "The pad 80% of the height",0.0,0.2,1.0,1.0,21);
  TPad *pad12 = new TPad("pad12", "The pad 20% of the height",0.0,0.0,1.0,0.2,22);
  pad11->SetFillColor(0); pad12->SetFillColor(0);
  pad11->Draw();
  pad12->Draw();
  
  pad11->cd();
  frame_pt_Z_upsi3->Draw();
  TH1* tmp6 = cut_mc3->createHistogram("Z_pT", Z_pT_bins);
  tmp6->Scale(N_Upsi3_S_Z_S.getVal()/cut_mc3->sumEntries());
  tmp6->Draw("hesame");
  tmp6->SetLineColor(kRed);
  tmp6->SetMarkerSize(0);
  c_weighted2->cd(3);
  
  pad12->cd();
  TH1* h6 = data_weighted_3->createHistogram("Z_pT", Z_pT_bins);
  h6->Sumw2();
  h6->Divide(tmp6);
  h6->Draw();
  h6->GetYaxis()->SetRangeUser(lowerLimitYForRatioPlot,upperLimitYForRatioPlot)  ;
  
  c_weighted2->SaveAs("c_weighted_Z_pT.pdf");
  
  //Canvas 3
  TCanvas *c_weighted3 = new TCanvas("c_weighted3", "c_weighted3", 1200, 400); c_weighted3->Divide(3,1);
  
  c_weighted3->cd(1); 
  TPad *pad13 = new TPad("pad13", "The pad 80% of the height",0.0,0.2,1.0,1.0,21);
  TPad *pad14 = new TPad("pad14", "The pad 20% of the height",0.0,0.0,1.0,0.2,22);
  pad13->SetFillColor(0); pad14->SetFillColor(0);
  pad13->Draw();
  pad14->Draw();
  
  pad13->cd();
  frame_Z_RAP_upsi1->Draw();
  TH1* tmp7 = cut_mc1->createHistogram("Z_RAPIDITY", Z_RAP_bins);
  tmp7->Scale(N_Upsi1_S_Z_S.getVal()/cut_mc1->sumEntries());
  tmp7->Draw("hesame");
  tmp7->SetLineColor(kRed);
  tmp7->SetMarkerSize(0);
  c_weighted3->cd(1);
  
  pad14->cd();
  TH1* h7 = data_weighted_1->createHistogram("Z_RAPIDITY", Z_RAP_bins);
  h7->Sumw2();
  h7->Divide(tmp7);
  h7->Draw();
  h7->GetYaxis()->SetRangeUser(lowerLimitYForRatioPlot,upperLimitYForRatioPlot)  ;
  
  c_weighted3->cd(2); 
  TPad *pad15 = new TPad("pad15", "The pad 80% of the height",0.0,0.2,1.0,1.0,21);
  TPad *pad16 = new TPad("pad16", "The pad 20% of the height",0.0,0.0,1.0,0.2,22);
  pad15->SetFillColor(0); pad16->SetFillColor(0);
  pad15->Draw();
  pad16->Draw();
  
  pad15->cd();
  frame_Z_RAP_upsi2->Draw();
  TH1* tmp8 = cut_mc2->createHistogram("Z_RAPIDITY", Z_RAP_bins);
  tmp8->Scale(N_Upsi2_S_Z_S.getVal()/cut_mc2->sumEntries());
  tmp8->Draw("hesame");
  tmp8->SetLineColor(kRed);
  tmp8->SetMarkerSize(0);
  c_weighted3->cd(2);
  
  pad16->cd();
  TH1* h8 = data_weighted_2->createHistogram("Z_RAPIDITY", Z_RAP_bins);
  h8->Sumw2();
  h8->Divide(tmp8);
  h8->Draw();
  h8->GetYaxis()->SetRangeUser(lowerLimitYForRatioPlot,upperLimitYForRatioPlot)  ;
  
  c_weighted3->cd(3); 
  TPad *pad17 = new TPad("pad17", "The pad 80% of the height",0.0,0.2,1.0,1.0,21);
  TPad *pad18 = new TPad("pad18", "The pad 20% of the height",0.0,0.0,1.0,0.2,22);
  pad17->SetFillColor(0); pad18->SetFillColor(0);
  pad17->Draw();
  pad18->Draw();
  
  pad17->cd();
  frame_Z_RAP_upsi3->Draw();
  TH1* tmp9 = cut_mc3->createHistogram("Z_RAPIDITY", Z_RAP_bins);
  tmp9->Scale(N_Upsi3_S_Z_S.getVal()/cut_mc3->sumEntries());
  tmp9->Draw("hesame");
  tmp9->SetLineColor(kRed);
  tmp9->SetMarkerSize(0);
  c_weighted3->cd(3);
  
  pad18->cd();
  TH1* h9 = data_weighted_3->createHistogram("Z_RAPIDITY", Z_RAP_bins);
  h9->Sumw2();
  h9->Divide(tmp9);
  h9->Draw();
  h9->GetYaxis()->SetRangeUser(lowerLimitYForRatioPlot,upperLimitYForRatioPlot)  ;
  
  c_weighted3->SaveAs("c_weighted_Z_RAP.pdf");
  
  //Canvas 4
  TCanvas *c_weighted4 = new TCanvas("c_weighted4", "c_weighted4", 1200, 400); c_weighted4->Divide(3,1);
  
  c_weighted4->cd(1);
  
  TPad *pad19 = new TPad("pad19", "The pad 80% of the height",0.0,0.2,1.0,1.0,21);
  TPad *pad20 = new TPad("pad20", "The pad 20% of the height",0.0,0.0,1.0,0.2,22);
  pad19->SetFillColor(0); pad20->SetFillColor(0);
  pad19->Draw();
  pad20->Draw();
  
  pad19->cd();
  frame_Z_phi_upsi1->Draw();
  TH1* tmp10 = cut_mc1->createHistogram("Z_phi", Z_phi_bins);
  tmp10->Scale(N_Upsi1_S_Z_S.getVal()/cut_mc1->sumEntries());
  tmp10->Draw("hesame");
  tmp10->SetLineColor(kRed);
  tmp10->SetMarkerSize(0);
  c_weighted4->cd(1);
  
  pad20->cd();
  TH1* h10 = data_weighted_1->createHistogram("Z_phi", Z_phi_bins);
  h10->Sumw2();
  h10->Divide(tmp10);
  h10->Draw();
  h10->GetYaxis()->SetRangeUser(lowerLimitYForRatioPlot,upperLimitYForRatioPlot)  ;
  
  c_weighted4->cd(2); 
  
  TPad *pad21 = new TPad("pad21", "The pad 80% of the height",0.0,0.2,1.0,1.0,21);
  TPad *pad22 = new TPad("pad22", "The pad 20% of the height",0.0,0.0,1.0,0.2,22);
  pad21->SetFillColor(0); pad22->SetFillColor(0);
  pad21->Draw();
  pad22->Draw();
  
  pad21->cd();
  frame_Z_phi_upsi2->Draw();
  TH1* tmp11 = cut_mc2->createHistogram("Z_phi", Z_phi_bins);
  tmp11->Scale(N_Upsi2_S_Z_S.getVal()/cut_mc2->sumEntries());
  tmp11->Draw("hesame");
  tmp11->SetLineColor(kRed);
  tmp11->SetMarkerSize(0);
  c_weighted4->cd(2);
  
  pad22->cd();
  TH1* h11 = data_weighted_2->createHistogram("Z_phi", Z_phi_bins);
  h11->Sumw2();
  h11->Divide(tmp11);
  h11->Draw();
  h11->GetYaxis()->SetRangeUser(lowerLimitYForRatioPlot,upperLimitYForRatioPlot)  ;
  
  c_weighted4->cd(3); 
   
  TPad *pad23 = new TPad("pad23", "The pad 80% of the height",0.0,0.2,1.0,1.0,23); //these last numbers, the 23 here and the 24 below, appear not to make a difference, and S.L. told me he wasn't sure what they do...
  TPad *pad24 = new TPad("pad24", "The pad 20% of the height",0.0,0.0,1.0,0.2,24);
  pad23->SetFillColor(0); pad24->SetFillColor(0);
  pad23->Draw();
  pad24->Draw();
  
  pad23->cd();
  frame_Z_phi_upsi3->Draw();
  TH1* tmp12 = cut_mc3->createHistogram("Z_phi", Z_phi_bins);
  tmp12->Scale(N_Upsi3_S_Z_S.getVal()/cut_mc3->sumEntries());
  tmp12->Draw("hesame");
  tmp12->SetLineColor(kRed);
  tmp12->SetMarkerSize(0);
  c_weighted4->cd(3);
  
  pad24->cd();
  TH1* h12 = data_weighted_3->createHistogram("Z_phi", Z_phi_bins);
  h12->Sumw2();
  h12->Divide(tmp12);
  h12->Draw();
  h12->GetYaxis()->SetRangeUser(lowerLimitYForRatioPlot,upperLimitYForRatioPlot)  ;
  
  c_weighted4->SaveAs("c_weighted_Z_phi.pdf");
  
  //Canvas 5
  TCanvas *c_weighted5 = new TCanvas("c_weighted5", "c_weighted5", 1200, 400); c_weighted5->Divide(3,1);
  
  c_weighted5->cd(1); 
  
  TPad *pad25 = new TPad("pad25", "The pad 80% of the height",0.0,0.2,1.0,1.0,23);
  TPad *pad26 = new TPad("pad26", "The pad 20% of the height",0.0,0.0,1.0,0.2,24);
  pad25->SetFillColor(0); pad26->SetFillColor(0);
  pad25->Draw();
  pad26->Draw();
  
  pad25->cd();
  frame_Z_eta_upsi1->Draw();
  TH1* tmp13 = cut_mc1->createHistogram("Z_eta", Z_eta_bins);
  std::cout << "CHECKPOINT 1BIS" << std::endl;
  std::cout << N_Upsi1_S_Z_S.getVal() << std::endl;
  std::cout << cut_mc1->sumEntries() << std::endl;
  tmp13->Scale(N_Upsi1_S_Z_S.getVal()/cut_mc1->sumEntries());
  tmp13->Draw("hesame");
  tmp13->SetLineColor(kRed);
  tmp13->SetMarkerSize(0);
  c_weighted5->cd(1);
  
  std::cout << "CHECKPOINT H" << std::endl;
  std::cout << pad26 << std::endl;
  pad26->cd();
  std::cout << "CHECKPOINT H1" << std::endl;
  TH1* h13 = data_weighted_1->createHistogram("Z_eta", Z_eta_bins);
  std::cout << "CHECKPOINT I" << std::endl;
  h13->Sumw2();
  h13->Divide(tmp13);
  h13->Draw();
  h13->GetYaxis()->SetRangeUser(lowerLimitYForRatioPlot,upperLimitYForRatioPlot)  ;
  
  //Start here 14 April 2023
  c_weighted5->cd(2); 
 
  TPad *pad27 = new TPad("pad27", "The pad 80% of the height",0.0,0.2,1.0,1.0,23);
  TPad *pad28 = new TPad("pad28", "The pad 20% of the height",0.0,0.0,1.0,0.2,24);
  pad27->SetFillColor(0); pad28->SetFillColor(0);
  pad27->Draw();
  pad28->Draw();
  
  pad27->cd();
  frame_Z_eta_upsi2->Draw();
  TH1* tmp14 = cut_mc2->createHistogram("Z_eta", Z_eta_bins);
  tmp14->Scale(N_Upsi2_S_Z_S.getVal()/cut_mc2->sumEntries());
  tmp14->Draw("hesame");
  tmp14->SetLineColor(kRed);
  tmp14->SetMarkerSize(0);
  c_weighted5->cd(2);
  
  pad28->cd();
  TH1* h14 = data_weighted_2->createHistogram("Z_eta", Z_eta_bins);
  h14->Sumw2();
  h14->Divide(tmp14);
  h14->Draw();
  h14->GetYaxis()->SetRangeUser(lowerLimitYForRatioPlot,upperLimitYForRatioPlot)  ;
  
  c_weighted5->cd(3); 
  
  TPad *pad29 = new TPad("pad29", "The pad 80% of the height",0.0,0.2,1.0,1.0,23);
  TPad *pad30 = new TPad("pad30", "The pad 20% of the height",0.0,0.0,1.0,0.2,24);
  pad29->SetFillColor(0); pad30->SetFillColor(0);
  pad29->Draw();
  pad30->Draw();
  
  pad29->cd();
  frame_Z_eta_upsi3->Draw();
  TH1* tmp15 = cut_mc3->createHistogram("Z_eta", Z_eta_bins);
  tmp15->Scale(N_Upsi3_S_Z_S.getVal()/cut_mc3->sumEntries());
  tmp15->Draw("hesame");
  tmp15->SetLineColor(kRed);
  tmp15->SetMarkerSize(0);
  c_weighted5->cd(3);
  
  pad30->cd();
  TH1* h15 = data_weighted_3->createHistogram("Z_eta", Z_eta_bins);
  h15->Sumw2();
  h15->Divide(tmp15);
  h15->Draw();
  h15->GetYaxis()->SetRangeUser(lowerLimitYForRatioPlot,upperLimitYForRatioPlot)  ;
  
  c_weighted5->SaveAs("c_weighted_Z_eta.pdf");
  
  //Canvas 6 
  TCanvas *c_weighted6 = new TCanvas("c_weighted6", "c_weighted6", 1200, 400); c_weighted6->Divide(3,1);
  
  c_weighted6->cd(1); 
  TPad *pad31 = new TPad("pad31", "The pad 80% of the height",0.0,0.2,1.0,1.0,23);
  TPad *pad32 = new TPad("pad32", "The pad 20% of the height",0.0,0.0,1.0,0.2,24);
  pad31->SetFillColor(0); pad32->SetFillColor(0);
  pad31->Draw();
  pad32->Draw();
  
  pad31->cd();
  frame_RAP_upsi1->Draw();
  TH1* tmp16 = cut_mc1->createHistogram("upsi_RAPIDITY", upsi_RAP_bins);
  std::cout << "CHECKPOINT 2" << std::endl;
  std::cout << N_Upsi1_S_Z_S.getVal() << std::endl;
  std::cout << cut_mc1->sumEntries() << std::endl; 
  std::cout << N_Upsi1_S_Z_S.getVal()/cut_mc1->sumEntries() << std::endl;
  tmp16->Scale(N_Upsi1_S_Z_S.getVal()/cut_mc1->sumEntries());
  tmp16->Draw("hesame");
  tmp16->SetLineColor(kRed);
  tmp16->SetMarkerSize(0);
  c_weighted6->cd(1);
  
  pad32->cd();
  TH1* h16 = data_weighted_1->createHistogram("upsi_RAPIDITY", upsi_RAP_bins);
  h16->Sumw2();
  h16->Divide(tmp16);
  h16->Draw();
  h16->GetYaxis()->SetRangeUser(lowerLimitYForRatioPlot,upperLimitYForRatioPlot)  ;
  
  c_weighted6->cd(2); 
  TPad *pad33 = new TPad("pad33", "The pad 80% of the height",0.0,0.2,1.0,1.0,23);
  TPad *pad34 = new TPad("pad34", "The pad 20% of the height",0.0,0.0,1.0,0.2,24);
  pad33->SetFillColor(0); pad34->SetFillColor(0);
  pad33->Draw();
  pad34->Draw();
  
  pad33->cd();
  frame_RAP_upsi2->Draw();
  TH1* tmp17 = cut_mc2->createHistogram("upsi_RAPIDITY", upsi_RAP_bins);
  tmp17->Scale(N_Upsi2_S_Z_S.getVal()/cut_mc2->sumEntries());
  tmp17->Draw("hesame");
  tmp17->SetLineColor(kRed);
  tmp17->SetMarkerSize(0);
  c_weighted6->cd(2);
  
  pad34->cd();
  TH1* h17 = data_weighted_2->createHistogram("upsi_RAPIDITY", upsi_RAP_bins);
  h17->Sumw2();
  h17->Divide(tmp17);
  h17->Draw();
  h17->GetYaxis()->SetRangeUser(lowerLimitYForRatioPlot,upperLimitYForRatioPlot)  ;
  
  c_weighted6->cd(3); 
  TPad *pad35 = new TPad("pad35", "The pad 80% of the height",0.0,0.2,1.0,1.0,23);
  TPad *pad36 = new TPad("pad36", "The pad 20% of the height",0.0,0.0,1.0,0.2,24);
  pad35->SetFillColor(0); pad36->SetFillColor(0);
  pad35->Draw();
  pad36->Draw();
  
  pad35->cd();
  frame_RAP_upsi3->Draw();
  TH1* tmp18 = cut_mc3->createHistogram("upsi_RAPIDITY", upsi_RAP_bins);
  tmp18->Scale(N_Upsi3_S_Z_S.getVal()/cut_mc3->sumEntries());
  tmp18->Draw("hesame");
  tmp18->SetLineColor(kRed);
  tmp18->SetMarkerSize(0);
  c_weighted6->cd(3);
  
  pad36->cd();
  TH1* h18 = data_weighted_3->createHistogram("upsi_RAPIDITY", upsi_RAP_bins);
  h18->Sumw2();
  h18->Divide(tmp18);
  h18->Draw();
  h18->GetYaxis()->SetRangeUser(lowerLimitYForRatioPlot,upperLimitYForRatioPlot)  ;
  
  c_weighted6->SaveAs("c_weighted_upsi_RAP.pdf");
  
  //Canvas 7
  TCanvas *c_weighted7 = new TCanvas("c_weighted7", "c_weighted7", 1200, 400); c_weighted7->Divide(3,1);
  
  c_weighted7->cd(1); 
  TPad *pad37 = new TPad("pad37", "The pad 80% of the height",0.0,0.2,1.0,1.0,23);
  TPad *pad38 = new TPad("pad38", "The pad 20% of the height",0.0,0.0,1.0,0.2,24);
  pad37->SetFillColor(0); pad38->SetFillColor(0);
  pad37->Draw();
  pad38->Draw();
  
  pad37->cd();
  frame_phi_upsi1->Draw();
  TH1* tmp19 = cut_mc1->createHistogram("upsi_phi", upsi_phi_bins);
  tmp19->Scale(N_Upsi1_S_Z_S.getVal()/cut_mc1->sumEntries());
  tmp19->Draw("hesame");
  tmp19->SetLineColor(kRed);
  tmp19->SetMarkerSize(0);
  c_weighted7->cd(1);
  
  pad38->cd();
  TH1* h19 = data_weighted_1->createHistogram("upsi_phi", upsi_phi_bins);
  h19->Sumw2();
  h19->Divide(tmp19);
  h19->Draw();
  h19->GetYaxis()->SetRangeUser(lowerLimitYForRatioPlot,upperLimitYForRatioPlot)  ;
  
  c_weighted7->cd(2); 
  TPad *pad39 = new TPad("pad39", "The pad 80% of the height",0.0,0.2,1.0,1.0,23);
  TPad *pad40 = new TPad("pad40", "The pad 20% of the height",0.0,0.0,1.0,0.2,24);
  pad39->SetFillColor(0); pad40->SetFillColor(0);
  pad39->Draw();
  pad40->Draw();
  
  pad39->cd();
  frame_phi_upsi2->Draw();
  TH1* tmp20 = cut_mc2->createHistogram("upsi_phi", upsi_phi_bins);
  tmp20->Scale(N_Upsi2_S_Z_S.getVal()/cut_mc2->sumEntries());
  tmp20->Draw("hesame");
  tmp20->SetLineColor(kRed);
  tmp20->SetMarkerSize(0);
  c_weighted7->cd(2);
  
  pad40->cd();
  TH1* h20 = data_weighted_2->createHistogram("upsi_phi", upsi_phi_bins);
  h20->Sumw2();
  h20->Divide(tmp20);
  h20->Draw();
  h20->GetYaxis()->SetRangeUser(lowerLimitYForRatioPlot,upperLimitYForRatioPlot)  ;
  
  c_weighted7->cd(3); 
  TPad *pad41 = new TPad("pad41", "The pad 80% of the height",0.0,0.2,1.0,1.0,23);
  TPad *pad42 = new TPad("pad42", "The pad 20% of the height",0.0,0.0,1.0,0.2,24);
  pad41->SetFillColor(0); pad42->SetFillColor(0);
  pad41->Draw();
  pad42->Draw();
  
  pad41->cd();
  frame_phi_upsi3->Draw();
  TH1* tmp21 = cut_mc3->createHistogram("upsi_phi", upsi_phi_bins);
  tmp21->Scale(N_Upsi3_S_Z_S.getVal()/cut_mc3->sumEntries());
  tmp21->Draw("hesame");
  tmp21->SetLineColor(kRed);
  tmp21->SetMarkerSize(0);
  c_weighted7->cd(3);
  
  pad42->cd();
  TH1* h21 = data_weighted_3->createHistogram("upsi_phi", upsi_phi_bins);
  h21->Sumw2();
  h21->Divide(tmp21);
  h21->Draw();
  h21->GetYaxis()->SetRangeUser(lowerLimitYForRatioPlot,upperLimitYForRatioPlot)  ;
  
  c_weighted7->SaveAs("c_weighted_upsi_phi.pdf");
  
  //Canvas 8
  TCanvas *c_weighted8 = new TCanvas("c_weighted8", "c_weighted8", 1200, 400); c_weighted8->Divide(3,1);
  
  c_weighted8->cd(1); 
  TPad *pad43 = new TPad("pad43", "The pad 80% of the height",0.0,0.2,1.0,1.0,23);
  TPad *pad44 = new TPad("pad44", "The pad 20% of the height",0.0,0.0,1.0,0.2,24);
  pad43->SetFillColor(0); pad44->SetFillColor(0);
  pad43->Draw();
  pad44->Draw();
  
  pad43->cd();
  frame_eta_upsi1->Draw();
  TH1* tmp22 = cut_mc1->createHistogram("upsi_eta", upsi_eta_bins);
  tmp22->Scale(N_Upsi1_S_Z_S.getVal()/cut_mc1->sumEntries());
  tmp22->Draw("hesame");
  tmp22->SetLineColor(kRed);
  tmp22->SetMarkerSize(0);
  c_weighted8->cd(1);
  
  pad44->cd();
  TH1* h22 = data_weighted_1->createHistogram("upsi_eta",upsi_eta_bins);
  h22->Sumw2();
  h22->Divide(tmp22);
  h22->Draw();
  h22->GetYaxis()->SetRangeUser(lowerLimitYForRatioPlot,upperLimitYForRatioPlot)  ;
  
  c_weighted8->cd(2); 
  TPad *pad45 = new TPad("pad45", "The pad 80% of the height",0.0,0.2,1.0,1.0,23);
  TPad *pad46 = new TPad("pad46", "The pad 20% of the height",0.0,0.0,1.0,0.2,24);
  pad45->SetFillColor(0); pad46->SetFillColor(0);
  pad45->Draw();
  pad46->Draw();
  
  pad45->cd();
  frame_eta_upsi2->Draw();
  TH1* tmp23 = cut_mc2->createHistogram("upsi_eta", upsi_eta_bins);
  tmp23->Scale(N_Upsi2_S_Z_S.getVal()/cut_mc2->sumEntries());
  tmp23->Draw("hesame");
  tmp23->SetLineColor(kRed);
  tmp23->SetMarkerSize(0);
  c_weighted8->cd(2);
  
  pad46->cd();
  TH1* h23 = data_weighted_2->createHistogram("upsi_eta", upsi_eta_bins);
  h23->Sumw2();
  h23->Divide(tmp23);
  h23->Draw();
  h23->GetYaxis()->SetRangeUser(lowerLimitYForRatioPlot,upperLimitYForRatioPlot)  ;
  
  c_weighted8->cd(3); 
  TPad *pad47 = new TPad("pad47", "The pad 80% of the height",0.0,0.2,1.0,1.0,23);
  TPad *pad48 = new TPad("pad48", "The pad 20% of the height",0.0,0.0,1.0,0.2,24);
  pad47->SetFillColor(0); pad48->SetFillColor(0);
  pad47->Draw();
  pad48->Draw();
  
  pad47->cd();
  frame_eta_upsi3->Draw();
  TH1* tmp24 = cut_mc3->createHistogram("upsi_eta", upsi_eta_bins);
  tmp24->Scale(N_Upsi3_S_Z_S.getVal()/cut_mc3->sumEntries());
  std::cout << "CHECKPOINT I1" << std::endl;
  std::cout << tmp24->GetBinContent(5) << std::endl;
  std::cout << tmp24->GetBinContent(4) << std::endl;
  tmp24->Draw("hesame");
  tmp24->SetLineColor(kRed);
  tmp24->SetMarkerSize(0);
  c_weighted8->cd(3);
  
  pad48->cd();
  TH1* h24 = data_weighted_3->createHistogram("upsi_eta", upsi_eta_bins);
  std::cout << "CHECKPOINT J" << std::endl;
  std::cout << h24->GetBinContent(5) << std::endl;
  std::cout << h24->GetBinContent(4) << std::endl;
  h24->Sumw2();
  h24->Divide(tmp24);
  h24->Draw();
  h24->GetYaxis()->SetRangeUser(lowerLimitYForRatioPlot,upperLimitYForRatioPlot)  ;
  
  
  c_weighted8->SaveAs("c_weighted_upsi_eta.pdf");
  
  //Canvas 9
  TCanvas *c_weighted9 = new TCanvas("c_weighted9", "c_weighted9", 1200, 400); c_weighted9->Divide(3,1);
  
  c_weighted9->cd(1); 
  TPad *pad49 = new TPad("pad49", "The pad 80% of the height",0.0,0.2,1.0,1.0,23);
  TPad *pad50 = new TPad("pad50", "The pad 20% of the height",0.0,0.0,1.0,0.2,24);
  pad49->SetFillColor(0); pad50->SetFillColor(0);
  pad49->Draw();
  pad50->Draw();
  
  pad49->cd();
  frame_lead_pT_mu_from_Z_pT_upsi1->Draw();
  TH1* tmp25 = cut_mc1->createHistogram("lead_pT_mu_from_Z_pT", lead_pT_mu_from_Z_pT_bins);
  tmp25->Scale(N_Upsi1_S_Z_S.getVal()/cut_mc1->sumEntries());
  tmp25->Draw("hesame");
  tmp25->SetLineColor(kRed);
  tmp25->SetMarkerSize(0);
  std::cout << "CHECKPOINT K" << std::endl;
  std::cout << tmp25->GetBinContent(4) << std::endl;
  c_weighted9->cd(1);
  
  pad50->cd();
  TH1* h25 = data_weighted_1->createHistogram("lead_pT_mu_from_Z_pT", lead_pT_mu_from_Z_pT_bins);
  h25->Sumw2();
  std::cout << h25->GetBinContent(4) << std::endl;
  h25->Divide(tmp25);
  h25->Draw();
  h25->GetYaxis()->SetRangeUser(lowerLimitYForRatioPlot,upperLimitYForRatioPlot)  ;
  
  c_weighted9->cd(2); 
  TPad *pad51 = new TPad("pad51", "The pad 80% of the height",0.0,0.2,1.0,1.0,23);
  TPad *pad52 = new TPad("pad52", "The pad 20% of the height",0.0,0.0,1.0,0.2,24);
  pad51->SetFillColor(0); pad52->SetFillColor(0);
  pad51->Draw();
  pad52->Draw();
  
  pad51->cd();
  frame_lead_pT_mu_from_Z_pT_upsi2->Draw();
  TH1* tmp26 = cut_mc2->createHistogram("lead_pT_mu_from_Z_pT", lead_pT_mu_from_Z_pT_bins);
  tmp26->Scale(N_Upsi2_S_Z_S.getVal()/cut_mc2->sumEntries());
  tmp26->Draw("hesame");
  tmp26->SetLineColor(kRed);
  tmp26->SetMarkerSize(0);
  c_weighted9->cd(2);
  
  pad52->cd();
  TH1* h26 = data_weighted_2->createHistogram("lead_pT_mu_from_Z_pT", lead_pT_mu_from_Z_pT_bins);
  h26->Sumw2();
  h26->Divide(tmp26);
  h26->Draw();
  h26->GetYaxis()->SetRangeUser(lowerLimitYForRatioPlot,upperLimitYForRatioPlot)  ;
  
  c_weighted9->cd(3); 
  TPad *pad53 = new TPad("pad53", "The pad 80% of the height",0.0,0.2,1.0,1.0,23);
  TPad *pad54 = new TPad("pad54", "The pad 20% of the height",0.0,0.0,1.0,0.2,24);
  pad53->SetFillColor(0); pad54->SetFillColor(0);
  pad53->Draw();
  pad54->Draw();
  
  pad53->cd();
  frame_lead_pT_mu_from_Z_pT_upsi3->Draw();
  TH1* tmp27 = cut_mc3->createHistogram("lead_pT_mu_from_Z_pT", lead_pT_mu_from_Z_pT_bins);
  tmp27->Scale(N_Upsi3_S_Z_S.getVal()/cut_mc3->sumEntries());
  tmp27->Draw("hesame");
  tmp27->SetLineColor(kRed);
  tmp27->SetMarkerSize(0);
  c_weighted9->cd(3);
  
  pad54->cd();
  TH1* h27 = data_weighted_3->createHistogram("lead_pT_mu_from_Z_pT", lead_pT_mu_from_Z_pT_bins);
  h27->Sumw2();
  h27->Divide(tmp27);
  h27->Draw();
  h27->GetYaxis()->SetRangeUser(lowerLimitYForRatioPlot,upperLimitYForRatioPlot)  ;
  
  c_weighted9->SaveAs("c_weighted_lead_pT_mu_from_Z_pT.pdf");
  
  //Canvas 10
  TCanvas *c_weighted10 = new TCanvas("c_weighted10", "c_weighted10", 1200, 400); c_weighted10->Divide(3,1);
  
  c_weighted10->cd(1); 
  TPad *pad55 = new TPad("pad55", "The pad 80% of the height",0.0,0.2,1.0,1.0,23);
  TPad *pad56 = new TPad("pad56", "The pad 20% of the height",0.0,0.0,1.0,0.2,24);
  pad55->SetFillColor(0); pad56->SetFillColor(0);
  pad55->Draw();
  pad56->Draw();
  
  pad55->cd();
  frame_lead_pT_mu_from_Z_RAP_upsi1->Draw();
  TH1* tmp28 = cut_mc1->createHistogram("lead_pT_mu_from_Z_RAPIDITY", lead_pT_mu_from_Z_RAP_bins);
  tmp28->Scale(N_Upsi1_S_Z_S.getVal()/cut_mc1->sumEntries());
  tmp28->Draw("hesame");
  tmp28->SetLineColor(kRed);
  tmp28->SetMarkerSize(0);
  c_weighted10->cd(1);
  
  pad56->cd();
  TH1* h28 = data_weighted_1->createHistogram("lead_pT_mu_from_Z_RAPIDITY", lead_pT_mu_from_Z_RAP_bins);
  h28->Sumw2();
  h28->Divide(tmp28);
  h28->Draw();
  h28->GetYaxis()->SetRangeUser(lowerLimitYForRatioPlot,upperLimitYForRatioPlot)  ;
 
  c_weighted10->cd(2); 
  TPad *pad57 = new TPad("pad57", "The pad 80% of the height",0.0,0.2,1.0,1.0,23);
  TPad *pad58 = new TPad("pad58", "The pad 20% of the height",0.0,0.0,1.0,0.2,24);
  pad57->SetFillColor(0); pad58->SetFillColor(0);
  pad57->Draw();
  pad58->Draw();
  
  pad57->cd();
  frame_lead_pT_mu_from_Z_RAP_upsi2->Draw();
  TH1* tmp29 = cut_mc2->createHistogram("lead_pT_mu_from_Z_RAPIDITY", lead_pT_mu_from_Z_RAP_bins);
  tmp29->Scale(N_Upsi2_S_Z_S.getVal()/cut_mc2->sumEntries());
  tmp29->Draw("hesame");
  tmp29->SetLineColor(kRed);
  tmp29->SetMarkerSize(0);
  c_weighted10->cd(2);
  
  pad58->cd();
  TH1* h29 = data_weighted_2->createHistogram("lead_pT_mu_from_Z_RAPIDITY", lead_pT_mu_from_Z_RAP_bins);
  h29->Sumw2();
  h29->Divide(tmp29);
  h29->Draw();
  h29->GetYaxis()->SetRangeUser(lowerLimitYForRatioPlot,upperLimitYForRatioPlot)  ;
 
  c_weighted10->cd(3); 
  TPad *pad59 = new TPad("pad59", "The pad 80% of the height",0.0,0.2,1.0,1.0,23);
  TPad *pad60 = new TPad("pad60", "The pad 20% of the height",0.0,0.0,1.0,0.2,24);
  pad59->SetFillColor(0); pad60->SetFillColor(0);
  pad59->Draw();
  pad60->Draw();
 
  pad59->cd();
  frame_lead_pT_mu_from_Z_RAP_upsi3->Draw();
  TH1* tmp30 = cut_mc3->createHistogram("lead_pT_mu_from_Z_RAPIDITY", lead_pT_mu_from_Z_RAP_bins);
  tmp30->Scale(N_Upsi3_S_Z_S.getVal()/cut_mc3->sumEntries());
  tmp30->Draw("hesame");
  tmp30->SetLineColor(kRed);
  tmp30->SetMarkerSize(0);
  c_weighted10->cd(3);
  
  pad60->cd();
  TH1* h30 = data_weighted_3->createHistogram("lead_pT_mu_from_Z_RAPIDITY", lead_pT_mu_from_Z_RAP_bins);
  h30->Sumw2();
  h30->Divide(tmp30);
  h30->Draw();
  h30->GetYaxis()->SetRangeUser(lowerLimitYForRatioPlot,upperLimitYForRatioPlot)  ;
  
  c_weighted10->SaveAs("c_weighted_lead_pT_mu_from_Z_RAP.pdf");
  
  //Canvas 11
  TCanvas *c_weighted11 = new TCanvas("c_weighted11", "c_weighted11", 1200, 400); c_weighted11->Divide(3,1);
  
  c_weighted11->cd(1); 
  TPad *pad61 = new TPad("pad61", "The pad 80% of the height",0.0,0.2,1.0,1.0,23);
  TPad *pad62 = new TPad("pad62", "The pad 20% of the height",0.0,0.0,1.0,0.2,24);
  pad61->SetFillColor(0); pad62->SetFillColor(0);
  pad61->Draw();
  pad62->Draw();
 
  pad61->cd();
  frame_lead_pT_mu_from_Z_eta_upsi1->Draw();
  TH1* tmp31 = cut_mc1->createHistogram("lead_pT_mu_from_Z_eta", lead_pT_mu_from_Z_eta_bins);
  tmp31->Scale(N_Upsi1_S_Z_S.getVal()/cut_mc1->sumEntries());
  tmp31->Draw("hesame");
  tmp31->SetLineColor(kRed);
  tmp31->SetMarkerSize(0);
  c_weighted11->cd(1);
  
  pad62->cd();
  TH1* h31 = data_weighted_1->createHistogram("lead_pT_mu_from_Z_eta", lead_pT_mu_from_Z_eta_bins);
  h31->Sumw2();
  h31->Divide(tmp31);
  h31->Draw();
  h31->GetYaxis()->SetRangeUser(lowerLimitYForRatioPlot,upperLimitYForRatioPlot)  ;
  
  c_weighted11->cd(2); 
  TPad *pad63 = new TPad("pad63", "The pad 80% of the height",0.0,0.2,1.0,1.0,23);
  TPad *pad64 = new TPad("pad64", "The pad 20% of the height",0.0,0.0,1.0,0.2,24);
  pad63->SetFillColor(0); pad64->SetFillColor(0);
  pad63->Draw();
  pad64->Draw();
  
  pad63->cd();
  frame_lead_pT_mu_from_Z_eta_upsi2->Draw();
  TH1* tmp32 = cut_mc2->createHistogram("lead_pT_mu_from_Z_eta", lead_pT_mu_from_Z_eta_bins);
  tmp32->Scale(N_Upsi2_S_Z_S.getVal()/cut_mc2->sumEntries());
  tmp32->Draw("hesame");
  tmp32->SetLineColor(kRed);
  tmp32->SetMarkerSize(0);
  c_weighted11->cd(2);
  
  pad64->cd();
  TH1* h32 = data_weighted_2->createHistogram("lead_pT_mu_from_Z_eta", lead_pT_mu_from_Z_eta_bins);
  h32->Sumw2();
  h32->Divide(tmp32);
  h32->Draw();
  h32->GetYaxis()->SetRangeUser(lowerLimitYForRatioPlot,upperLimitYForRatioPlot)  ;
  
  c_weighted11->cd(3); 
  TPad *pad65 = new TPad("pad65", "The pad 80% of the height",0.0,0.2,1.0,1.0,23);
  TPad *pad66 = new TPad("pad66", "The pad 20% of the height",0.0,0.0,1.0,0.2,24);
  pad65->SetFillColor(0); pad66->SetFillColor(0);
  pad65->Draw();
  pad66->Draw();
  
  pad65->cd();
  frame_lead_pT_mu_from_Z_eta_upsi3->Draw();
  TH1* tmp33 = cut_mc3->createHistogram("lead_pT_mu_from_Z_eta", lead_pT_mu_from_Z_eta_bins);
  tmp33->Scale(N_Upsi3_S_Z_S.getVal()/cut_mc3->sumEntries());
  tmp33->Draw("hesame");
  tmp33->SetLineColor(kRed);
  tmp33->SetMarkerSize(0);
  c_weighted11->cd(3);
  
  pad66->cd();
  TH1* h33 = data_weighted_3->createHistogram("lead_pT_mu_from_Z_eta", lead_pT_mu_from_Z_eta_bins);
  h33->Sumw2();
  h33->Divide(tmp33);
  h33->Draw();
  h33->GetYaxis()->SetRangeUser(lowerLimitYForRatioPlot,upperLimitYForRatioPlot)  ;
  
  c_weighted11->SaveAs("c_weighted_lead_pT_mu_from_Z_eta.pdf");
  
  
  //Canvas 12
  TCanvas *c_weighted12 = new TCanvas("c_weighted12", "c_weighted12", 1200, 400); c_weighted12->Divide(3,1);
  
  c_weighted12->cd(1); 
  TPad *pad67 = new TPad("pad67", "The pad 80% of the height",0.0,0.2,1.0,1.0,23);
  TPad *pad68 = new TPad("pad68", "The pad 20% of the height",0.0,0.0,1.0,0.2,24);
  pad67->SetFillColor(0); pad68->SetFillColor(0);
  pad67->Draw();
  pad68->Draw();
  
  pad67->cd();
  frame_lead_pT_mu_from_Z_phi_upsi1->Draw();
  TH1* tmp34 = cut_mc1->createHistogram("lead_pT_mu_from_Z_phi", lead_pT_mu_from_Z_phi_bins);
  tmp34->Scale(N_Upsi1_S_Z_S.getVal()/cut_mc1->sumEntries());
  tmp34->Draw("hesame");
  tmp34->SetLineColor(kRed);
  tmp34->SetMarkerSize(0);
  c_weighted12->cd(1);
  
  pad68->cd();
  TH1* h34 = data_weighted_1->createHistogram("lead_pT_mu_from_Z_phi", lead_pT_mu_from_Z_phi_bins);
  h34->Sumw2();
  h34->Divide(tmp34);
  h34->Draw();
  h34->GetYaxis()->SetRangeUser(lowerLimitYForRatioPlot,upperLimitYForRatioPlot)  ;
  
  c_weighted12->cd(2); 
  TPad *pad69 = new TPad("pad69", "The pad 80% of the height",0.0,0.2,1.0,1.0,23);
  TPad *pad70 = new TPad("pad70", "The pad 20% of the height",0.0,0.0,1.0,0.2,24);
  pad69->SetFillColor(0); pad70->SetFillColor(0);
  pad69->Draw();
  pad70->Draw();
  
  pad69->cd();
  frame_lead_pT_mu_from_Z_phi_upsi2->Draw();
  TH1* tmp35 = cut_mc2->createHistogram("lead_pT_mu_from_Z_phi", lead_pT_mu_from_Z_phi_bins);
  tmp35->Scale(N_Upsi2_S_Z_S.getVal()/cut_mc2->sumEntries());
  tmp35->Draw("hesame");
  tmp35->SetLineColor(kRed);
  tmp35->SetMarkerSize(0);
  c_weighted12->cd(2);
  
  pad70->cd();
  TH1* h35 = data_weighted_2->createHistogram("lead_pT_mu_from_Z_phi", lead_pT_mu_from_Z_phi_bins);
  h35->Sumw2();
  h35->Divide(tmp35);
  h35->Draw();
  h35->GetYaxis()->SetRangeUser(lowerLimitYForRatioPlot,upperLimitYForRatioPlot)  ;
  
  c_weighted12->cd(3); 
  TPad *pad71 = new TPad("pad71", "The pad 80% of the height",0.0,0.2,1.0,1.0,23);
  TPad *pad72 = new TPad("pad72", "The pad 20% of the height",0.0,0.0,1.0,0.2,24);
  pad71->SetFillColor(0); pad72->SetFillColor(0);
  pad71->Draw();
  pad72->Draw();
  
  pad71->cd();
  frame_lead_pT_mu_from_Z_phi_upsi3->Draw();
  TH1* tmp36 = cut_mc3->createHistogram("lead_pT_mu_from_Z_phi", lead_pT_mu_from_Z_phi_bins);
  tmp36->Scale(N_Upsi3_S_Z_S.getVal()/cut_mc3->sumEntries());
  tmp36->Draw("hesame");
  tmp36->SetLineColor(kRed);
  tmp36->SetMarkerSize(0);
  c_weighted12->cd(3);
  
  pad72->cd();
  TH1* h36 = data_weighted_3->createHistogram("lead_pT_mu_from_Z_phi", lead_pT_mu_from_Z_phi_bins);
  h36->Sumw2();
  h36->Divide(tmp36);
  h36->Draw();
  h36->GetYaxis()->SetRangeUser(lowerLimitYForRatioPlot,upperLimitYForRatioPlot)  ;
  
  c_weighted12->SaveAs("c_weighted_lead_pT_mu_from_Z_phi.pdf");
  
  //Canvas 13
  TCanvas *c_weighted13 = new TCanvas("c_weighted13", "c_weighted13", 1200, 400); c_weighted13->Divide(3,1);
  
  c_weighted13->cd(1); 
  TPad *pad73 = new TPad("pad73", "The pad 80% of the height",0.0,0.2,1.0,1.0,23);
  TPad *pad74 = new TPad("pad74", "The pad 20% of the height",0.0,0.0,1.0,0.2,24);
  pad73->SetFillColor(0); pad74->SetFillColor(0);
  pad73->Draw();
  pad74->Draw();
  
  pad73->cd();
  frame_sublead_pT_mu_from_Z_pT_upsi1->Draw();
  TH1* tmp37 = cut_mc1->createHistogram("sublead_pT_mu_from_Z_pT", sublead_pT_mu_from_Z_pT_bins);
  tmp37->Scale(N_Upsi1_S_Z_S.getVal()/cut_mc1->sumEntries());
  tmp37->Draw("hesame");
  tmp37->SetLineColor(kRed);
  tmp37->SetMarkerSize(0);
  c_weighted13->cd(1);
  
  pad74->cd();
  TH1* h37 = data_weighted_1->createHistogram("sublead_pT_mu_from_Z_pT", sublead_pT_mu_from_Z_pT_bins);
  h37->Sumw2();
  h37->Divide(tmp37);
  h37->Draw();
  h37->GetYaxis()->SetRangeUser(lowerLimitYForRatioPlot,upperLimitYForRatioPlot)  ;
  
  c_weighted13->cd(2); 
  TPad *pad75 = new TPad("pad75", "The pad 80% of the height",0.0,0.2,1.0,1.0,23);
  TPad *pad76 = new TPad("pad76", "The pad 20% of the height",0.0,0.0,1.0,0.2,24);
  pad75->SetFillColor(0); pad76->SetFillColor(0);
  pad75->Draw();
  pad76->Draw();
  
  pad75->cd();
  frame_sublead_pT_mu_from_Z_pT_upsi2->Draw();
  TH1* tmp38 = cut_mc2->createHistogram("sublead_pT_mu_from_Z_pT", sublead_pT_mu_from_Z_pT_bins);
  tmp38->Scale(N_Upsi2_S_Z_S.getVal()/cut_mc2->sumEntries());
  tmp38->Draw("hesame");
  tmp38->SetLineColor(kRed);
  tmp38->SetMarkerSize(0);
  c_weighted13->cd(2);
  
  pad76->cd();
  TH1* h38 = data_weighted_2->createHistogram("sublead_pT_mu_from_Z_pT", sublead_pT_mu_from_Z_pT_bins);
  h38->Sumw2();
  h38->Divide(tmp38);
  h38->Draw();
  h38->GetYaxis()->SetRangeUser(lowerLimitYForRatioPlot,upperLimitYForRatioPlot)  ;
  
  c_weighted13->cd(3); 
  TPad *pad77 = new TPad("pad77", "The pad 80% of the height",0.0,0.2,1.0,1.0,23);
  TPad *pad78 = new TPad("pad78", "The pad 20% of the height",0.0,0.0,1.0,0.2,24);
  pad77->SetFillColor(0); pad78->SetFillColor(0);
  pad77->Draw();
  pad78->Draw();
  
  pad77->cd();
  frame_sublead_pT_mu_from_Z_pT_upsi3->Draw();
  TH1* tmp39 = cut_mc3->createHistogram("sublead_pT_mu_from_Z_pT", sublead_pT_mu_from_Z_pT_bins);
  tmp39->Scale(N_Upsi3_S_Z_S.getVal()/cut_mc3->sumEntries());
  tmp39->Draw("hesame");
  tmp39->SetLineColor(kRed);
  tmp39->SetMarkerSize(0);
  c_weighted13->cd(3);
  
  pad78->cd();
  TH1* h39 = data_weighted_3->createHistogram("sublead_pT_mu_from_Z_pT", sublead_pT_mu_from_Z_pT_bins);
  h39->Sumw2();
  h39->Divide(tmp39);
  h39->Draw();
  h39->GetYaxis()->SetRangeUser(lowerLimitYForRatioPlot,upperLimitYForRatioPlot)  ;
  
  c_weighted13->SaveAs("c_weighted_sublead_pT_mu_from_Z_pT.pdf");
  
  //Canvas 14
  TCanvas *c_weighted14 = new TCanvas("c_weighted14", "c_weighted14", 1200,400); c_weighted14->Divide(3,1);
  
  c_weighted14->cd(1); 
  TPad *pad79 = new TPad("pad79", "The pad 80% of the height",0.0,0.2,1.0,1.0,23);
  TPad *pad80 = new TPad("pad80", "The pad 20% of the height",0.0,0.0,1.0,0.2,24);
  pad79->SetFillColor(0); pad80->SetFillColor(0);
  pad79->Draw();
  pad80->Draw();
  
  pad79->cd();
  frame_sublead_pT_mu_from_Z_RAP_upsi1->Draw();
  TH1* tmp40 = cut_mc1->createHistogram("sublead_pT_mu_from_Z_RAPIDITY", sublead_pT_mu_from_Z_RAP_bins);
  tmp40->Scale(N_Upsi1_S_Z_S.getVal()/cut_mc1->sumEntries());
  tmp40->Draw("hesame");
  tmp40->SetLineColor(kRed);
  tmp40->SetMarkerSize(0);
  c_weighted14->cd(1);
  
  pad80->cd();
  TH1* h40 = data_weighted_1->createHistogram("sublead_pT_mu_from_Z_RAPIDITY", sublead_pT_mu_from_Z_RAP_bins);
  h40->Sumw2();
  h40->Divide(tmp40);
  h40->Draw();
  h40->GetYaxis()->SetRangeUser(lowerLimitYForRatioPlot,upperLimitYForRatioPlot)  ;
  
  c_weighted14->cd(2); 
  TPad *pad81 = new TPad("pad81", "The pad 80% of the height",0.0,0.2,1.0,1.0,23);
  TPad *pad82 = new TPad("pad82", "The pad 20% of the height",0.0,0.0,1.0,0.2,24);
  pad81->SetFillColor(0); pad82->SetFillColor(0);
  pad81->Draw();
  pad82->Draw();
  
  pad81->cd();
  frame_sublead_pT_mu_from_Z_RAP_upsi2->Draw();
  TH1* tmp41 = cut_mc2->createHistogram("sublead_pT_mu_from_Z_RAPIDITY", sublead_pT_mu_from_Z_RAP_bins);
  tmp41->Scale(N_Upsi2_S_Z_S.getVal()/cut_mc2->sumEntries());
  tmp41->Draw("hesame");
  tmp41->SetLineColor(kRed);
  tmp41->SetMarkerSize(0);
  c_weighted14->cd(2);
  
  pad82->cd();
  TH1* h41 = data_weighted_2->createHistogram("sublead_pT_mu_from_Z_RAPIDITY", sublead_pT_mu_from_Z_RAP_bins);
  h41->Sumw2();
  h41->Divide(tmp41);
  h41->Draw();
  h41->GetYaxis()->SetRangeUser(lowerLimitYForRatioPlot,upperLimitYForRatioPlot)  ;
  
  c_weighted14->cd(3); 
  TPad *pad83 = new TPad("pad83", "The pad 80% of the height",0.0,0.2,1.0,1.0,23);
  TPad *pad84 = new TPad("pad84", "The pad 20% of the height",0.0,0.0,1.0,0.2,24);
  pad83->SetFillColor(0); pad84->SetFillColor(0);
  pad83->Draw();
  pad84->Draw();
  
  pad83->cd();
  frame_sublead_pT_mu_from_Z_RAP_upsi3->Draw();
  TH1* tmp42 = cut_mc3->createHistogram("sublead_pT_mu_from_Z_RAPIDITY", sublead_pT_mu_from_Z_RAP_bins);
  tmp42->Scale(N_Upsi3_S_Z_S.getVal()/cut_mc3->sumEntries());
  tmp42->Draw("hesame");
  tmp42->SetLineColor(kRed);
  tmp42->SetMarkerSize(0);
  c_weighted14->cd(3);
  
  pad84->cd();
  TH1* h42 = data_weighted_3->createHistogram("sublead_pT_mu_from_Z_RAPIDITY", sublead_pT_mu_from_Z_RAP_bins);
  h42->Sumw2();
  h42->Divide(tmp42);
  h42->Draw();
  h42->GetYaxis()->SetRangeUser(lowerLimitYForRatioPlot,upperLimitYForRatioPlot)  ;
  
 
  c_weighted14->SaveAs("c_weighted_sublead_pT_mu_from_Z_RAP.pdf");
  
  //Canvas 15
  TCanvas *c_weighted15 = new TCanvas("c_weighted15", "c_weighted15", 1200, 400); c_weighted15->Divide(3,1);
  
  c_weighted15->cd(1);
  TPad *pad85 = new TPad("pad85", "The pad 80% of the height",0.0,0.2,1.0,1.0,23);
  TPad *pad86 = new TPad("pad86", "The pad 20% of the height",0.0,0.0,1.0,0.2,24);
  pad85->SetFillColor(0); pad86->SetFillColor(0);
  pad85->Draw();
  pad86->Draw();
  
  pad85->cd();
  frame_sublead_pT_mu_from_Z_eta_upsi1->Draw();
  TH1* tmp43 = cut_mc1->createHistogram("sublead_pT_mu_from_Z_eta", sublead_pT_mu_from_Z_eta_bins);
  tmp43->Scale(N_Upsi1_S_Z_S.getVal()/cut_mc1->sumEntries());
  tmp43->Draw("hesame");
  tmp43->SetLineColor(kRed);
  tmp43->SetMarkerSize(0);
  c_weighted15->cd(1);
  
  pad86->cd();
  TH1* h43 = data_weighted_1->createHistogram("sublead_pT_mu_from_Z_eta", sublead_pT_mu_from_Z_eta_bins);
  h43->Sumw2();
  h43->Divide(tmp43);
  h43->Draw();
  h43->GetYaxis()->SetRangeUser(lowerLimitYForRatioPlot,upperLimitYForRatioPlot)  ;
  
  
  c_weighted15->cd(2); 
  TPad *pad87 = new TPad("pad87", "The pad 80% of the height",0.0,0.2,1.0,1.0,23);
  TPad *pad88 = new TPad("pad88", "The pad 20% of the height",0.0,0.0,1.0,0.2,24);
  pad87->SetFillColor(0); pad88->SetFillColor(0);
  pad87->Draw();
  pad88->Draw();
  
  pad87->cd();
  frame_sublead_pT_mu_from_Z_eta_upsi2->Draw();
  TH1* tmp44 = cut_mc2->createHistogram("sublead_pT_mu_from_Z_eta", sublead_pT_mu_from_Z_eta_bins);
  tmp44->Scale(N_Upsi2_S_Z_S.getVal()/cut_mc2->sumEntries());
  tmp44->Draw("hesame");
  tmp44->SetLineColor(kRed);
  tmp44->SetMarkerSize(0);
  c_weighted15->cd(2);
  
  pad88->cd();
  TH1* h44 = data_weighted_2->createHistogram("sublead_pT_mu_from_Z_eta", sublead_pT_mu_from_Z_eta_bins);
  h44->Sumw2();
  h44->Divide(tmp44);
  h44->Draw();
  h44->GetYaxis()->SetRangeUser(lowerLimitYForRatioPlot,upperLimitYForRatioPlot)  ;
  
  
  c_weighted15->cd(3); 
  TPad *pad89 = new TPad("pad89", "The pad 80% of the height",0.0,0.2,1.0,1.0,23);
  TPad *pad90 = new TPad("pad90", "The pad 20% of the height",0.0,0.0,1.0,0.2,24);
  pad89->SetFillColor(0); pad90->SetFillColor(0);
  pad89->Draw();
  pad90->Draw();
  
  pad89->cd();
  frame_sublead_pT_mu_from_Z_eta_upsi3->Draw();
  TH1* tmp45 = cut_mc3->createHistogram("sublead_pT_mu_from_Z_eta", sublead_pT_mu_from_Z_eta_bins);
  tmp45->Scale(N_Upsi3_S_Z_S.getVal()/cut_mc3->sumEntries());
  tmp45->Draw("hesame");
  tmp45->SetLineColor(kRed);
  tmp45->SetMarkerSize(0);
  c_weighted15->cd(3);
  
  pad90->cd();
  TH1* h45 = data_weighted_3->createHistogram("sublead_pT_mu_from_Z_eta", sublead_pT_mu_from_Z_eta_bins);
  h45->Sumw2();
  h45->Divide(tmp45);
  h45->Draw();
  h45->GetYaxis()->SetRangeUser(lowerLimitYForRatioPlot,upperLimitYForRatioPlot)  ;
  
   c_weighted15->SaveAs("c_weighted_sublead_pT_mu_from_Z_eta.pdf");
  
  //Canvas 16
  TCanvas *c_weighted16 = new TCanvas("c_weighted16", "c_weighted16", 1200,400); c_weighted16->Divide(3,1);
  
  c_weighted16->cd(1); 
  TPad *pad91 = new TPad("pad91", "The pad 80% of the height",0.0,0.2,1.0,1.0,23);
  TPad *pad92 = new TPad("pad92", "The pad 20% of the height",0.0,0.0,1.0,0.2,24);
  pad91->SetFillColor(0); pad92->SetFillColor(0);
  pad91->Draw();
  pad92->Draw();
  
  pad91->cd();
  frame_sublead_pT_mu_from_Z_phi_upsi1->Draw();
  TH1* tmp46 = cut_mc1->createHistogram("sublead_pT_mu_from_Z_phi", sublead_pT_mu_from_Z_phi_bins);
  tmp46->Scale(N_Upsi1_S_Z_S.getVal()/cut_mc1->sumEntries());
  tmp46->Draw("hesame");
  tmp46->SetLineColor(kRed);
  tmp46->SetMarkerSize(0);
  c_weighted16->cd(1);
  
  pad92->cd();
  TH1* h46 = data_weighted_1->createHistogram("sublead_pT_mu_from_Z_phi", sublead_pT_mu_from_Z_phi_bins);
  h46->Sumw2();
  h46->Divide(tmp46);
  h46->Draw();
  h46->GetYaxis()->SetRangeUser(lowerLimitYForRatioPlot,upperLimitYForRatioPlot)  ;
  
  c_weighted16->cd(2); 
  TPad *pad93 = new TPad("pad93", "The pad 80% of the height",0.0,0.2,1.0,1.0,23);
  TPad *pad94 = new TPad("pad94", "The pad 20% of the height",0.0,0.0,1.0,0.2,24);
  pad93->SetFillColor(0); pad94->SetFillColor(0);
  pad93->Draw();
  pad94->Draw();
  
  pad93->cd();
  frame_sublead_pT_mu_from_Z_phi_upsi2->Draw();
  TH1* tmp47 = cut_mc2->createHistogram("sublead_pT_mu_from_Z_phi", sublead_pT_mu_from_Z_phi_bins);
  tmp47->Scale(N_Upsi2_S_Z_S.getVal()/cut_mc2->sumEntries());
  tmp47->Draw("hesame");
  tmp47->SetLineColor(kRed);
  tmp47->SetMarkerSize(0);
  c_weighted16->cd(2);
  
  pad94->cd();
  TH1* h47 = data_weighted_2->createHistogram("sublead_pT_mu_from_Z_phi", sublead_pT_mu_from_Z_phi_bins);
  h47->Sumw2();
  h47->Divide(tmp47);
  h47->Draw();
  h47->GetYaxis()->SetRangeUser(lowerLimitYForRatioPlot,upperLimitYForRatioPlot)  ;
  
  c_weighted16->cd(3); 
  TPad *pad95 = new TPad("pad95", "The pad 80% of the height",0.0,0.2,1.0,1.0,23);
  TPad *pad96 = new TPad("pad96", "The pad 20% of the height",0.0,0.0,1.0,0.2,24);
  pad95->SetFillColor(0); pad96->SetFillColor(0);
  pad95->Draw();
  pad96->Draw();
  
  pad95->cd();
  frame_sublead_pT_mu_from_Z_phi_upsi3->Draw();
  TH1* tmp48 = cut_mc3->createHistogram("sublead_pT_mu_from_Z_phi", sublead_pT_mu_from_Z_phi_bins);
  tmp48->Scale(N_Upsi3_S_Z_S.getVal()/cut_mc3->sumEntries());
  tmp48->Draw("hesame");
  tmp48->SetLineColor(kRed);
  tmp48->SetMarkerSize(0);
  c_weighted16->cd(3);
  
  pad96->cd();
  TH1* h48 = data_weighted_3->createHistogram("sublead_pT_mu_from_Z_phi", sublead_pT_mu_from_Z_phi_bins);
  h48->Sumw2();
  h48->Divide(tmp48);
  h48->Draw();
  h48->GetYaxis()->SetRangeUser(lowerLimitYForRatioPlot,upperLimitYForRatioPlot)  ;
  
  c_weighted16->SaveAs("c_weighted_sublead_pT_mu_from_Z_phi.pdf");
  
  //Canvas 17
  TCanvas *c_weighted17 = new TCanvas("c_weighted17", "c_weighted17", 1200, 400); c_weighted17->Divide(3,1);
  
  c_weighted17->cd(1); 
  TPad *pad97 = new TPad("pad97", "The pad 80% of the height",0.0,0.2,1.0,1.0,23);
  TPad *pad98 = new TPad("pad98", "The pad 20% of the height",0.0,0.0,1.0,0.2,24);
  pad97->SetFillColor(0); pad98->SetFillColor(0);
  pad97->Draw();
  pad98->Draw();
  
  pad97->cd();
  frame_lead_pT_mu_from_upsi_pT_upsi1->Draw();
  TH1* tmp49 = cut_mc1->createHistogram("lead_pT_mu_from_upsi_pT", lead_pT_mu_from_upsi_pT_bins);
  tmp49->Scale(N_Upsi1_S_Z_S.getVal()/cut_mc1->sumEntries());
  tmp49->Draw("hesame");
  tmp49->SetLineColor(kRed);
  tmp49->SetMarkerSize(0);
  c_weighted17->cd();
  
  pad98->cd();
  TH1* h49 = data_weighted_1->createHistogram("lead_pT_mu_from_upsi_pT", lead_pT_mu_from_upsi_pT_bins);
  h49->Sumw2();
  h49->Divide(tmp49);
  h49->Draw();
  h49->GetYaxis()->SetRangeUser(lowerLimitYForRatioPlot,upperLimitYForRatioPlot)  ;
  
  c_weighted17->cd(2); 
  TPad *pad99 = new TPad("pad99", "The pad 80% of the height",0.0,0.2,1.0,1.0,23);
  TPad *pad100 = new TPad("pad100", "The pad 20% of the height",0.0,0.0,1.0,0.2,24);
  pad99->SetFillColor(0); pad100->SetFillColor(0);
  pad99->Draw();
  pad100->Draw();
  
  pad99->cd();
  frame_lead_pT_mu_from_upsi_pT_upsi2->Draw();
  TH1* tmp50 = cut_mc2->createHistogram("lead_pT_mu_from_upsi_pT", lead_pT_mu_from_upsi_pT_bins);
  tmp50->Scale(N_Upsi2_S_Z_S.getVal()/cut_mc2->sumEntries());
  tmp50->Draw("hesame");
  tmp50->SetLineColor(kRed);
  tmp50->SetMarkerSize(0);
  c_weighted17->cd(2);
  
  pad100->cd();
  TH1* h50 = data_weighted_2->createHistogram("lead_pT_mu_from_upsi_pT", lead_pT_mu_from_upsi_pT_bins);
  h50->Sumw2();
  h50->Divide(tmp50);
  h50->Draw();
  h50->GetYaxis()->SetRangeUser(lowerLimitYForRatioPlot,upperLimitYForRatioPlot)  ;
  
  
  c_weighted17->cd(3); 
  TPad *pad101 = new TPad("pad101", "The pad 80% of the height",0.0,0.2,1.0,1.0,23);
  TPad *pad102 = new TPad("pad102", "The pad 20% of the height",0.0,0.0,1.0,0.2,24);
  pad101->SetFillColor(0); pad102->SetFillColor(0);
  pad101->Draw();
  pad102->Draw();
  
  pad101->cd();
  frame_lead_pT_mu_from_upsi_pT_upsi3->Draw();
  TH1* tmp51 = cut_mc3->createHistogram("lead_pT_mu_from_upsi_pT", lead_pT_mu_from_upsi_pT_bins);
  tmp51->Scale(N_Upsi3_S_Z_S.getVal()/cut_mc3->sumEntries());
  tmp51->Draw("hesame");
  tmp51->SetLineColor(kRed);
  tmp51->SetMarkerSize(0);
  c_weighted17->cd(3);
  
  pad102->cd();
  TH1* h51 = data_weighted_3->createHistogram("lead_pT_mu_from_upsi_pT", lead_pT_mu_from_upsi_pT_bins);
  h51->Sumw2();
  h51->Divide(tmp51);
  h51->Draw();
  h51->GetYaxis()->SetRangeUser(lowerLimitYForRatioPlot,upperLimitYForRatioPlot)  ;
  
  c_weighted17->SaveAs("c_weighted_lead_pT_mu_from_upsi_pT.pdf");
  
  //Canvas 18
  TCanvas *c_weighted18 = new TCanvas("c_weighted18", "c_weighted18", 1200, 400); c_weighted18->Divide(3,1);
  
  c_weighted18->cd(1); 
  TPad *pad103 = new TPad("pad103", "The pad 80% of the height",0.0,0.2,1.0,1.0,23);
  TPad *pad104 = new TPad("pad104", "The pad 20% of the height",0.0,0.0,1.0,0.2,24);
  pad103->SetFillColor(0); pad104->SetFillColor(0);
  pad103->Draw();
  pad104->Draw();
  
  pad103->cd();
  frame_lead_pT_mu_from_upsi_RAP_upsi1->Draw();
  TH1* tmp52 = cut_mc1->createHistogram("lead_pT_mu_from_upsi_RAPIDITY", lead_pT_mu_from_upsi_RAP_bins);
  tmp52->Scale(N_Upsi1_S_Z_S.getVal()/cut_mc1->sumEntries());
  tmp52->Draw("hesame");
  tmp52->SetLineColor(kRed);
  tmp52->SetMarkerSize(0);
  c_weighted18->cd(1);
  
  pad104->cd();
  TH1* h52 = data_weighted_1->createHistogram("lead_pT_mu_from_upsi_RAPIDITY", lead_pT_mu_from_upsi_RAP_bins);
  h52->Sumw2();
  h52->Divide(tmp52);
  h52->Draw();
  h52->GetYaxis()->SetRangeUser(lowerLimitYForRatioPlot,upperLimitYForRatioPlot)  ;
  
  c_weighted18->cd(2); 
  
  TPad *pad105 = new TPad("pad105", "The pad 80% of the height",0.0,0.2,1.0,1.0,23);
  TPad *pad106 = new TPad("pad106", "The pad 20% of the height",0.0,0.0,1.0,0.2,24);
  pad105->SetFillColor(0); pad106->SetFillColor(0);
  pad105->Draw();
  pad106->Draw();
  
  pad105->cd();
  frame_lead_pT_mu_from_upsi_RAP_upsi2->Draw();
  TH1* tmp53 = cut_mc2->createHistogram("lead_pT_mu_from_upsi_RAPIDITY", lead_pT_mu_from_upsi_RAP_bins);
  tmp53->Scale(N_Upsi2_S_Z_S.getVal()/cut_mc2->sumEntries());
  tmp53->Draw("hesame");
  tmp53->SetLineColor(kRed);
  tmp53->SetMarkerSize(0);
  c_weighted18->cd(2);
  
  pad106->cd();
  TH1* h53 = data_weighted_2->createHistogram("lead_pT_mu_from_upsi_RAPIDITY", lead_pT_mu_from_upsi_RAP_bins);
  h53->Sumw2();
  h53->Divide(tmp53);
  h53->Draw();
  h53->GetYaxis()->SetRangeUser(lowerLimitYForRatioPlot,upperLimitYForRatioPlot)  ;
  
  c_weighted18->cd(3); 
  
  TPad *pad107 = new TPad("pad107", "The pad 80% of the height",0.0,0.2,1.0,1.0,23);
  TPad *pad108 = new TPad("pad108", "The pad 20% of the height",0.0,0.0,1.0,0.2,24);
  pad107->SetFillColor(0); pad108->SetFillColor(0);
  pad107->Draw();
  pad108->Draw();
  
  pad107->cd();
  frame_lead_pT_mu_from_upsi_RAP_upsi3->Draw();
  TH1* tmp54 = cut_mc3->createHistogram("lead_pT_mu_from_upsi_RAPIDITY", lead_pT_mu_from_upsi_RAP_bins);
  tmp54->Scale(N_Upsi3_S_Z_S.getVal()/cut_mc3->sumEntries());
  tmp54->Draw("hesame");
  tmp54->SetLineColor(kRed);
  tmp54->SetMarkerSize(0);
  c_weighted18->cd(3);
  
  pad108->cd();
  TH1* h54 = data_weighted_3->createHistogram("lead_pT_mu_from_upsi_RAPIDITY", lead_pT_mu_from_upsi_RAP_bins);
  h54->Sumw2();
  h54->Divide(tmp54);
  h54->Draw();
  h54->GetYaxis()->SetRangeUser(lowerLimitYForRatioPlot,upperLimitYForRatioPlot)  ;
  
  c_weighted18->SaveAs("c_weighted_lead_pT_mu_from_upsi_RAP.pdf");
  
  //Canvas 19
  TCanvas *c_weighted19 = new TCanvas("c_weighted19", "c_weighted19", 1200, 400); c_weighted19->Divide(3,1);
  
  c_weighted19->cd(1); 
  TPad *pad109 = new TPad("pad109", "The pad 80% of the height",0.0,0.2,1.0,1.0,23);
  TPad *pad110 = new TPad("pad110", "The pad 20% of the height",0.0,0.0,1.0,0.2,24);
  pad109->SetFillColor(0); pad110->SetFillColor(0);
  pad109->Draw();
  pad110->Draw();
  
  pad109->cd();
  frame_lead_pT_mu_from_upsi_eta_upsi1->Draw();
  TH1* tmp55 = cut_mc1->createHistogram("lead_pT_mu_from_upsi_eta", lead_pT_mu_from_upsi_eta_bins);
  tmp55->Scale(N_Upsi1_S_Z_S.getVal()/cut_mc1->sumEntries());
  tmp55->Draw("hesame");
  tmp55->SetLineColor(kRed);
  tmp55->SetMarkerSize(0);
  c_weighted19->cd(1);
  
  pad110->cd();
  TH1* h55 = data_weighted_1->createHistogram("lead_pT_mu_from_upsi_eta", lead_pT_mu_from_upsi_eta_bins);
  h55->Sumw2();
  h55->Divide(tmp55);
  h55->Draw();
  h55->GetYaxis()->SetRangeUser(lowerLimitYForRatioPlot,upperLimitYForRatioPlot)  ;
  
  c_weighted19->cd(2); 
  TPad *pad111 = new TPad("pad111", "The pad 80% of the height",0.0,0.2,1.0,1.0,23);
  TPad *pad112 = new TPad("pad112", "The pad 20% of the height",0.0,0.0,1.0,0.2,24);
  pad111->SetFillColor(0); pad112->SetFillColor(0);
  pad111->Draw();
  pad112->Draw();
  
  pad111->cd();
  frame_lead_pT_mu_from_upsi_eta_upsi2->Draw();
  TH1* tmp56 = cut_mc2->createHistogram("lead_pT_mu_from_upsi_eta", lead_pT_mu_from_upsi_eta_bins);
  tmp56->Scale(N_Upsi2_S_Z_S.getVal()/cut_mc2->sumEntries());
  tmp56->Draw("hesame");
  tmp56->SetLineColor(kRed);
  tmp56->SetMarkerSize(0);
  c_weighted19->cd(2);
  
  pad112->cd();
  TH1* h56 = data_weighted_2->createHistogram("lead_pT_mu_from_upsi_eta", lead_pT_mu_from_upsi_eta_bins);
  h56->Sumw2();
  h56->Divide(tmp56);
  h56->Draw();
  h56->GetYaxis()->SetRangeUser(lowerLimitYForRatioPlot,upperLimitYForRatioPlot)  ;
  
  c_weighted19->cd(3); 
  TPad *pad113 = new TPad("pad113", "The pad 80% of the height",0.0,0.2,1.0,1.0,23);
  TPad *pad114 = new TPad("pad114", "The pad 20% of the height",0.0,0.0,1.0,0.2,24);
  pad113->SetFillColor(0); pad114->SetFillColor(0);
  pad113->Draw();
  pad114->Draw();
  
  pad113->cd();
  frame_lead_pT_mu_from_upsi_eta_upsi3->Draw();
  TH1* tmp57 = cut_mc3->createHistogram("lead_pT_mu_from_upsi_eta", lead_pT_mu_from_upsi_eta_bins);
  tmp57->Scale(N_Upsi3_S_Z_S.getVal()/cut_mc3->sumEntries());
  tmp57->Draw("hesame");
  tmp57->SetLineColor(kRed);
  tmp57->SetMarkerSize(0);
  c_weighted19->cd(3);
  
  pad114->cd();
  TH1* h57 = data_weighted_3->createHistogram("lead_pT_mu_from_upsi_eta", lead_pT_mu_from_upsi_eta_bins);
  h57->Sumw2();
  h57->Divide(tmp57);
  h57->Draw();
  h57->GetYaxis()->SetRangeUser(lowerLimitYForRatioPlot,upperLimitYForRatioPlot)  ;
  
  c_weighted19->SaveAs("c_weighted_lead_pT_mu_from_upsi_eta.pdf");
  
  //Canvas 20
  TCanvas *c_weighted20 = new TCanvas("c_weighted20", "c_weighted20", 1200,400); c_weighted20->Divide(3,1);
  
  c_weighted20->cd(1); 
  TPad *pad115 = new TPad("pad115", "The pad 80% of the height",0.0,0.2,1.0,1.0,23);
  TPad *pad116 = new TPad("pad116", "The pad 20% of the height",0.0,0.0,1.0,0.2,24);
  pad115->SetFillColor(0); pad116->SetFillColor(0);
  pad115->Draw();
  pad116->Draw();
  
  pad115->cd();
  frame_lead_pT_mu_from_upsi_phi_upsi1->Draw();
  TH1* tmp58 = cut_mc1->createHistogram("lead_pT_mu_from_upsi_phi", lead_pT_mu_from_upsi_phi_bins);
  tmp58->Scale(N_Upsi1_S_Z_S.getVal()/cut_mc1->sumEntries());
  tmp58->Draw("hesame");
  tmp58->SetLineColor(kRed);
  tmp58->SetMarkerSize(0);
  c_weighted20->cd(1);
  
  pad116->cd();
  TH1* h58 = data_weighted_1->createHistogram("lead_pT_mu_from_upsi_phi", lead_pT_mu_from_upsi_phi_bins);
  h58->Sumw2();
  h58->Divide(tmp58);
  h58->Draw();
  h58->GetYaxis()->SetRangeUser(lowerLimitYForRatioPlot,upperLimitYForRatioPlot)  ;
  
  c_weighted20->cd(2); 
  
  TPad *pad117 = new TPad("pad117", "The pad 80% of the height",0.0,0.2,1.0,1.0,23);
  TPad *pad118 = new TPad("pad118", "The pad 20% of the height",0.0,0.0,1.0,0.2,24);
  pad117->SetFillColor(0); pad118->SetFillColor(0);
  pad117->Draw();
  pad118->Draw();
  
  pad117->cd();
  frame_lead_pT_mu_from_upsi_phi_upsi2->Draw();
  TH1* tmp59 = cut_mc2->createHistogram("lead_pT_mu_from_upsi_phi", lead_pT_mu_from_upsi_phi_bins);
  tmp59->Scale(N_Upsi2_S_Z_S.getVal()/cut_mc2->sumEntries());
  tmp59->Draw("hesame");
  tmp59->SetLineColor(kRed);
  tmp59->SetMarkerSize(0);
  c_weighted20->cd(2);
  
  pad118->cd();
  TH1* h59 = data_weighted_2->createHistogram("lead_pT_mu_from_upsi_phi", lead_pT_mu_from_upsi_phi_bins);
  h59->Sumw2();
  h59->Divide(tmp59);
  h59->Draw();
  h59->GetYaxis()->SetRangeUser(lowerLimitYForRatioPlot,upperLimitYForRatioPlot)  ;
  
  c_weighted20->cd(3); 
  TPad *pad119 = new TPad("pad119", "The pad 80% of the height",0.0,0.2,1.0,1.0,23);
  TPad *pad120 = new TPad("pad120", "The pad 20% of the height",0.0,0.0,1.0,0.2,24);
  pad119->SetFillColor(0); pad120->SetFillColor(0);
  pad119->Draw();
  pad120->Draw();
  
  pad119->cd();
  frame_lead_pT_mu_from_upsi_phi_upsi3->Draw();
  TH1* tmp60 = cut_mc3->createHistogram("lead_pT_mu_from_upsi_phi", lead_pT_mu_from_upsi_phi_bins);
  tmp60->Scale(N_Upsi3_S_Z_S.getVal()/cut_mc3->sumEntries());
  tmp60->Draw("hesame");
  tmp60->SetLineColor(kRed);
  tmp60->SetMarkerSize(0);
  c_weighted20->cd(3);
  
  pad120->cd();
  TH1* h60 = data_weighted_3->createHistogram("lead_pT_mu_from_upsi_phi", lead_pT_mu_from_upsi_phi_bins);
  h60->Sumw2();
  h60->Divide(tmp60);
  h60->Draw();
  h60->GetYaxis()->SetRangeUser(lowerLimitYForRatioPlot,upperLimitYForRatioPlot)  ;
  
  c_weighted20->SaveAs("c_weighted_lead_pT_mu_from_upsi_phi.pdf");
  
  //Canvas 21
  TCanvas *c_weighted21 = new TCanvas("c_weighted21", "c_weighted21", 1200,400); c_weighted21->Divide(3,1);
  
  c_weighted21->cd(1); 
  TPad *pad121 = new TPad("pad121", "The pad 80% of the height",0.0,0.2,1.0,1.0,23);
  TPad *pad122 = new TPad("pad122", "The pad 20% of the height",0.0,0.0,1.0,0.2,24);
  pad121->SetFillColor(0); pad122->SetFillColor(0);
  pad121->Draw();
  pad122->Draw();
  
  pad121->cd();
  frame_sublead_pT_mu_from_upsi_pT_upsi1->Draw();
  TH1* tmp61 = cut_mc1->createHistogram("sublead_pT_mu_from_upsi_pT", sublead_pT_mu_from_upsi_pT_bins);
  tmp61->Scale(N_Upsi1_S_Z_S.getVal()/cut_mc1->sumEntries());
  tmp61->Draw("hesame");
  tmp61->SetLineColor(kRed);
  tmp61->SetMarkerSize(0);
  c_weighted21->cd(1);
  
  pad122->cd();
  TH1* h61 = data_weighted_1->createHistogram("sublead_pT_mu_from_upsi_pT", sublead_pT_mu_from_upsi_pT_bins);
  h61->Sumw2();
  h61->Divide(tmp61);
  h61->Draw();
  h61->GetYaxis()->SetRangeUser(lowerLimitYForRatioPlot,upperLimitYForRatioPlot)  ;
  
  c_weighted21->cd(2);
  TPad *pad123 = new TPad("pad123", "The pad 80% of the height",0.0,0.2,1.0,1.0,23);
  TPad *pad124 = new TPad("pad124", "The pad 20% of the height",0.0,0.0,1.0,0.2,24);
  pad123->SetFillColor(0); pad124->SetFillColor(0);
  pad123->Draw();
  pad124->Draw();
  
  pad123->cd();
  frame_sublead_pT_mu_from_upsi_pT_upsi2->Draw();
  TH1* tmp62 = cut_mc2->createHistogram("sublead_pT_mu_from_upsi_pT", sublead_pT_mu_from_upsi_pT_bins);
  tmp62->Scale(N_Upsi2_S_Z_S.getVal()/cut_mc2->sumEntries());
  tmp62->Draw("hesame");
  tmp62->SetLineColor(kRed);
  tmp62->SetMarkerSize(0);
  c_weighted21->cd(2);
  
  pad124->cd();
  TH1* h62 = data_weighted_2->createHistogram("sublead_pT_mu_from_upsi_pT", sublead_pT_mu_from_upsi_pT_bins);
  h62->Sumw2();
  h62->Divide(tmp62);
  h62->Draw();
  h62->GetYaxis()->SetRangeUser(lowerLimitYForRatioPlot,upperLimitYForRatioPlot)  ;
  
  c_weighted21->cd(3); 
  TPad *pad125 = new TPad("pad125", "The pad 80% of the height",0.0,0.2,1.0,1.0,23);
  TPad *pad126 = new TPad("pad126", "The pad 20% of the height",0.0,0.0,1.0,0.2,24);
  pad125->SetFillColor(0); pad126->SetFillColor(0);
  pad125->Draw();
  pad126->Draw();
   
  pad125->cd(); 
  frame_sublead_pT_mu_from_upsi_pT_upsi3->Draw();
  TH1* tmp63 = cut_mc3->createHistogram("sublead_pT_mu_from_upsi_pT", sublead_pT_mu_from_upsi_pT_bins);
  tmp63->Scale(N_Upsi3_S_Z_S.getVal()/cut_mc3->sumEntries());
  tmp63->Draw("hesame");
  tmp63->SetLineColor(kRed);
  tmp63->SetMarkerSize(0);
  c_weighted21->cd(3);
  
  pad126->cd();
  TH1* h63 = data_weighted_3->createHistogram("sublead_pT_mu_from_upsi_pT", sublead_pT_mu_from_upsi_pT_bins);
  h63->Sumw2();
  h63->Divide(tmp63);
  h63->Draw();
  h63->GetYaxis()->SetRangeUser(lowerLimitYForRatioPlot,upperLimitYForRatioPlot)  ;
  
  c_weighted21->SaveAs("c_weighted_sublead_pT_mu_from_upsi_pT.pdf");
  
  //Canvas 22
  TCanvas *c_weighted22 = new TCanvas("c_weighted22", "c_weighted22", 1200,400); c_weighted22->Divide(3,1);
  
  c_weighted22->cd(1); 
  TPad *pad127 = new TPad("pad127", "The pad 80% of the height",0.0,0.2,1.0,1.0,23);
  TPad *pad128 = new TPad("pad128", "The pad 20% of the height",0.0,0.0,1.0,0.2,24);
  pad127->SetFillColor(0); pad128->SetFillColor(0);
  pad127->Draw();
  pad128->Draw();
   
  pad127->cd(); 
  frame_sublead_pT_mu_from_upsi_RAP_upsi1->Draw();
  TH1* tmp64 = cut_mc1->createHistogram("sublead_pT_mu_from_upsi_RAPIDITY", sublead_pT_mu_from_upsi_RAP_bins);
  tmp64->Scale(N_Upsi1_S_Z_S.getVal()/cut_mc1->sumEntries());
  tmp64->Draw("hesame");
  tmp64->SetLineColor(kRed);
  tmp64->SetMarkerSize(0);
  c_weighted22->cd(1);
  
  pad128->cd();
  TH1* h64 = data_weighted_1->createHistogram("sublead_pT_mu_from_upsi_RAPIDITY", sublead_pT_mu_from_upsi_RAP_bins);
  h64->Sumw2();
  h64->Divide(tmp64);
  h64->Draw();
  h64->GetYaxis()->SetRangeUser(lowerLimitYForRatioPlot,upperLimitYForRatioPlot)  ;
  
  c_weighted22->cd(2); 
  TPad *pad129 = new TPad("pad129", "The pad 80% of the height",0.0,0.2,1.0,1.0,23);
  TPad *pad130 = new TPad("pad130", "The pad 20% of the height",0.0,0.0,1.0,0.2,24);
  pad129->SetFillColor(0); pad130->SetFillColor(0);
  pad129->Draw();
  pad130->Draw();
   
  pad129->cd(); 
  frame_sublead_pT_mu_from_upsi_RAP_upsi2->Draw();
  TH1* tmp65 = cut_mc2->createHistogram("sublead_pT_mu_from_upsi_RAPIDITY", sublead_pT_mu_from_upsi_RAP_bins);
  tmp65->Scale(N_Upsi2_S_Z_S.getVal()/cut_mc2->sumEntries());
  tmp65->Draw("hesame");
  tmp65->SetLineColor(kRed);
  tmp65->SetMarkerSize(0);
  c_weighted22->cd(2);
  
  pad130->cd();
  TH1* h65 = data_weighted_2->createHistogram("sublead_pT_mu_from_upsi_RAPIDITY", sublead_pT_mu_from_upsi_RAP_bins);
  h65->Sumw2();
  h65->Divide(tmp65);
  h65->Draw();
  h65->GetYaxis()->SetRangeUser(lowerLimitYForRatioPlot,upperLimitYForRatioPlot)  ;
  
  c_weighted22->cd(3);
  TPad *pad131 = new TPad("pad131", "The pad 80% of the height",0.0,0.2,1.0,1.0,23);
  TPad *pad132 = new TPad("pad132", "The pad 20% of the height",0.0,0.0,1.0,0.2,24);
  pad131->SetFillColor(0); pad132->SetFillColor(0);
  pad131->Draw();
  pad132->Draw();
   
  pad131->cd();  
  frame_sublead_pT_mu_from_upsi_RAP_upsi3->Draw();
  TH1* tmp66 = cut_mc3->createHistogram("sublead_pT_mu_from_upsi_RAPIDITY", sublead_pT_mu_from_upsi_RAP_bins);
  tmp66->Scale(N_Upsi3_S_Z_S.getVal()/cut_mc3->sumEntries());
  tmp66->Draw("hesame");
  tmp66->SetLineColor(kRed);
  tmp66->SetMarkerSize(0);
  c_weighted22->cd(3);
  
  pad132->cd();
  TH1* h66 = data_weighted_3->createHistogram("sublead_pT_mu_from_upsi_RAPIDITY", sublead_pT_mu_from_upsi_RAP_bins);
  h66->Sumw2();
  h66->Divide(tmp66);
  h66->Draw();
  h66->GetYaxis()->SetRangeUser(lowerLimitYForRatioPlot,upperLimitYForRatioPlot)  ;
  
  c_weighted22->SaveAs("c_weighted_sublead_pT_mu_from_upsi_RAP.pdf");
  
  //Canvas 23
  TCanvas *c_weighted23 = new TCanvas("c_weighted23", "c_weighted23",1200, 400); c_weighted23->Divide(3,1);
  
  c_weighted23->cd(1); 
  TPad *pad133 = new TPad("pad133", "The pad 80% of the height",0.0,0.2,1.0,1.0,23);
  TPad *pad134 = new TPad("pad134", "The pad 20% of the height",0.0,0.0,1.0,0.2,24);
  pad133->SetFillColor(0); pad134->SetFillColor(0);
  pad133->Draw();
  pad134->Draw();
  
  pad133->cd();
  frame_sublead_pT_mu_from_upsi_eta_upsi1->Draw();
  TH1* tmp67 = cut_mc1->createHistogram("sublead_pT_mu_from_upsi_eta", sublead_pT_mu_from_upsi_eta_bins);
  tmp67->Scale(N_Upsi1_S_Z_S.getVal()/cut_mc1->sumEntries());
  tmp67->Draw("hesame");
  tmp67->SetLineColor(kRed);
  tmp67->SetMarkerSize(0);
  c_weighted23->cd(1);
  
  pad134->cd();
  TH1* h67 = data_weighted_1->createHistogram("sublead_pT_mu_from_upsi_eta", sublead_pT_mu_from_upsi_eta_bins);
  h67->Sumw2();
  h67->Divide(tmp67);
  h67->Draw();
  h67->GetYaxis()->SetRangeUser(lowerLimitYForRatioPlot,upperLimitYForRatioPlot)  ;
  
  c_weighted23->cd(2); 
  TPad *pad135 = new TPad("pad135", "The pad 80% of the height",0.0,0.2,1.0,1.0,23);
  TPad *pad136 = new TPad("pad136", "The pad 20% of the height",0.0,0.0,1.0,0.2,24);
  pad135->SetFillColor(0); pad136->SetFillColor(0);
  pad135->Draw();
  pad136->Draw();
  
  pad135->cd();
  frame_sublead_pT_mu_from_upsi_eta_upsi2->Draw();
  TH1* tmp68 = cut_mc2->createHistogram("sublead_pT_mu_from_upsi_eta", sublead_pT_mu_from_upsi_eta_bins);
  tmp68->Scale(N_Upsi2_S_Z_S.getVal()/cut_mc2->sumEntries());
  tmp68->Draw("hesame");
  tmp68->SetLineColor(kRed);
  tmp68->SetMarkerSize(0);
  c_weighted23->cd(2);
  
  pad136->cd();
  TH1* h68 = data_weighted_2->createHistogram("sublead_pT_mu_from_upsi_eta", sublead_pT_mu_from_upsi_eta_bins);
  h68->Sumw2();
  h68->Divide(tmp68);
  h68->Draw();
  h68->GetYaxis()->SetRangeUser(lowerLimitYForRatioPlot,upperLimitYForRatioPlot)  ;
  
  c_weighted23->cd(3); 
  TPad *pad137 = new TPad("pad137", "The pad 80% of the height",0.0,0.2,1.0,1.0,23);
  TPad *pad138 = new TPad("pad138", "The pad 20% of the height",0.0,0.0,1.0,0.2,24);
  pad137->SetFillColor(0); pad138->SetFillColor(0);
  pad137->Draw();
  pad138->Draw();
  
  pad137->cd();
  frame_sublead_pT_mu_from_upsi_eta_upsi3->Draw();
  TH1* tmp69 = cut_mc3->createHistogram("sublead_pT_mu_from_upsi_eta", sublead_pT_mu_from_upsi_eta_bins);
  tmp69->Scale(N_Upsi3_S_Z_S.getVal()/cut_mc3->sumEntries());
  tmp69->Draw("hesame");
  tmp69->SetLineColor(kRed);
  tmp69->SetMarkerSize(0);
  c_weighted23->cd(3);
  
  pad138->cd();
  TH1* h69 = data_weighted_3->createHistogram("sublead_pT_mu_from_upsi_eta", sublead_pT_mu_from_upsi_eta_bins);
  h69->Sumw2();
  h69->Divide(tmp69);
  h69->Draw();
  h69->GetYaxis()->SetRangeUser(lowerLimitYForRatioPlot,upperLimitYForRatioPlot)  ;
  
  c_weighted23->SaveAs("c_weighted_sublead_pT_mu_from_upsi_eta.pdf");
  
  //Canvas 24
  TCanvas *c_weighted24 = new TCanvas("c_weighted24", "c_weighted24", 1200, 400); c_weighted24->Divide(3,1);
  
  c_weighted24->cd(1); 
  TPad *pad139 = new TPad("pad139", "The pad 80% of the height",0.0,0.2,1.0,1.0,23);
  TPad *pad140 = new TPad("pad140", "The pad 20% of the height",0.0,0.0,1.0,0.2,24);
  pad139->SetFillColor(0); pad140->SetFillColor(0);
  pad139->Draw();
  pad140->Draw();
  
  pad139->cd();
  frame_sublead_pT_mu_from_upsi_phi_upsi1->Draw();
  TH1* tmp70 = cut_mc1->createHistogram("sublead_pT_mu_from_upsi_phi", sublead_pT_mu_from_upsi_phi_bins);
  tmp70->Scale(N_Upsi1_S_Z_S.getVal()/cut_mc1->sumEntries());
  tmp70->Draw("hesame");
  tmp70->SetLineColor(kRed);
  tmp70->SetMarkerSize(0);
  c_weighted24->cd(1);
  
  pad140->cd();
  TH1* h70 = data_weighted_1->createHistogram("sublead_pT_mu_from_upsi_phi", sublead_pT_mu_from_upsi_phi_bins);
  h70->Sumw2();
  h70->Divide(tmp70);
  h70->Draw();
  h70->GetYaxis()->SetRangeUser(lowerLimitYForRatioPlot,upperLimitYForRatioPlot)  ;
  
  c_weighted24->cd(2); 
  TPad *pad141 = new TPad("pad141", "The pad 80% of the height",0.0,0.2,1.0,1.0,23);
  TPad *pad142 = new TPad("pad142", "The pad 20% of the height",0.0,0.0,1.0,0.2,24);
  pad141->SetFillColor(0); pad142->SetFillColor(0);
  pad141->Draw();
  pad142->Draw();
  
  pad141->cd();
  frame_sublead_pT_mu_from_upsi_phi_upsi2->Draw();
  TH1* tmp71 = cut_mc2->createHistogram("sublead_pT_mu_from_upsi_phi", sublead_pT_mu_from_upsi_phi_bins);
  tmp71->Scale(N_Upsi2_S_Z_S.getVal()/cut_mc2->sumEntries());
  tmp71->Draw("hesame");
  tmp71->SetLineColor(kRed);
  tmp71->SetMarkerSize(0);
  c_weighted24->cd(2);
  
  pad142->cd();
  TH1* h71 = data_weighted_2->createHistogram("sublead_pT_mu_from_upsi_phi", sublead_pT_mu_from_upsi_phi_bins);
  h71->Sumw2();
  h71->Divide(tmp71);
  h71->Draw();
  h71->GetYaxis()->SetRangeUser(lowerLimitYForRatioPlot,upperLimitYForRatioPlot)  ;
  
  c_weighted24->cd(3); 
  TPad *pad143 = new TPad("pad143", "The pad 80% of the height",0.0,0.2,1.0,1.0,23);
  TPad *pad144 = new TPad("pad144", "The pad 20% of the height",0.0,0.0,1.0,0.2,24);
  pad143->SetFillColor(0); pad144->SetFillColor(0);
  pad143->Draw();
  pad144->Draw();
  
  pad143->cd();
  frame_sublead_pT_mu_from_upsi_phi_upsi3->Draw();
  TH1* tmp72 = cut_mc3->createHistogram("sublead_pT_mu_from_upsi_phi", sublead_pT_mu_from_upsi_phi_bins);
  tmp72->Scale(N_Upsi3_S_Z_S.getVal()/cut_mc3->sumEntries());
  tmp72->Draw("hesame");
  tmp72->SetLineColor(kRed);
  tmp72->SetMarkerSize(0);
  c_weighted24->cd(3);
  
  pad144->cd();
  TH1* h72 = data_weighted_3->createHistogram("sublead_pT_mu_from_upsi_phi", sublead_pT_mu_from_upsi_phi_bins);
  h72->Sumw2();
  h72->Divide(tmp72);
  h72->Draw();
  h72->GetYaxis()->SetRangeUser(lowerLimitYForRatioPlot,upperLimitYForRatioPlot)  ;
  
  c_weighted24->SaveAs("c_weighted_sublead_pT_mu_from_upsi_phi.pdf");
  
   #ifdef CalculatePulls
  int Ntoys = 500;
  RooMCStudy *SimToy = new RooMCStudy(eSum, RooArgSet(upsi_mass,Z_mass), Binned(kFALSE), Silence(), Extended(), FitOptions(Save(kTRUE)), PrintEvalErrors(0), Minos(kTRUE));

  SimToy->generateAndFit(Ntoys);

  RooPlot *N_Upsi1_S_Z_S_pull_frame = SimToy->plotPull(N_Upsi1_S_Z_S, -4, 4, 8, kFALSE);
  RooPlot *N_Upsi2_S_Z_S_pull_frame = SimToy->plotPull(N_Upsi2_S_Z_S, -4, 4, 8, kFALSE);
  RooPlot *N_Upsi3_S_Z_S_pull_frame = SimToy->plotPull(N_Upsi3_S_Z_S, -4, 4, 8, kFALSE);

  TCanvas* cPulls = new TCanvas("cPulls","toymc", 1200, 400); cPulls->Divide(3,1);
  cPulls->cd(1); N_Upsi1_S_Z_S_pull_frame->Draw(); N_Upsi1_S_Z_S_pull_frame->SetXTitle("Z+Y(1S) signal yield");
  cPulls->cd(2); N_Upsi2_S_Z_S_pull_frame->Draw(); N_Upsi2_S_Z_S_pull_frame->SetXTitle("Z+Y(2S) signal yield");
  cPulls->cd(3); N_Upsi3_S_Z_S_pull_frame->Draw(); N_Upsi3_S_Z_S_pull_frame->SetXTitle("Z+Y(3S) signal yield");
  cPulls->SaveAs("c_pulls.pdf");
#endif


}

