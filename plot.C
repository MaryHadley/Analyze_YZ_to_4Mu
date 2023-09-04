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
#include "TLatex.h"
#include "TMath.h"
#include "TPaveLabel.h"
#include "TRandom3.h"
#include "TStyle.h"
#include "TTree.h"

  #define DrawPulls
//#define CalculatePulls
//#define datadriven
#define DPSSPSfit
//#define getkernels
//#define getSPSkernels
#define makeplots

//#define ConstraintY2SY3S

using namespace RooFit ;
using namespace RooStats ;

void fit();
void plot() { fit(); }
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
  RooRealVar Z_pT("Z_pT","Z_pT",0, 500); //Needed to raise the allowed range for the Z_pT to 1100 in order to not have any MC events excluded
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
 
 //Add variables to use for SPS vs. DPS determination
 RooRealVar deltaPhi_upsiZ("deltaPhi_upsiZ", "#Delta#phi(Z, Y(nS)", 0., PI);
 RooRealVar deltaRAP_upsiZ("deltaRAP_upsiZ", "#Deltay(Z, Y(nS)", 0., 6.);
 
 
 TFile *ntuple_data = new TFile("ntuple_data_2016_2017_2018_useGlobalMuonsFalse_upsiMuPtCut3_newVariablesIncluded.root");
 TTree* tree_data = (TTree*) ntuple_data->Get("tree");
    
// TFile *ntuple_mc = new TFile("MC_YZ_10KEvents.root");
 TFile *ntuple_mc = new TFile("MC_DPS_lumiWeighted_Run2_Total_YZ.root");
 TTree* tree_mc      = (TTree*) ntuple_mc->Get("tree");

//  TFile *ntuple_mc_sps  = new TFile("sps.root");
 TFile *ntuple_mc_sps  = new TFile("SPS_2018_Y1SZ_test_with_weights.root");
 TTree* tree_mc_sps    = (TTree*) ntuple_mc_sps->Get("tree");

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
      
    //Add variables that are useful for SPS vs. DPS separation
    Variables.add(deltaPhi_upsiZ);
    Variables.add(deltaRAP_upsiZ);
// //   
//Adding this as a test
//   Variables.add(big4MuVtxProb);
  //////////////////////////////////////////////////////////////

    RooArgSet VariablesSPS(upsi_mass, Z_mass, upsi_pT); //If you try to put too many in here, things will break, so better to do Variables.add(blah blah blah) as I do below when you want to add something
    VariablesSPS.add(Z_pT);
    VariablesSPS.add(Z_RAPIDITY);
    VariablesSPS.add(Z_phi);
    VariablesSPS.add(Z_eta);
    VariablesSPS.add(upsi_RAPIDITY);
    VariablesSPS.add(upsi_phi);
    VariablesSPS.add(upsi_eta);
    VariablesSPS.add(lead_pT_mu_from_Z_pT);
    VariablesSPS.add(lead_pT_mu_from_Z_RAPIDITY);
    VariablesSPS.add(lead_pT_mu_from_Z_eta);
    VariablesSPS.add(lead_pT_mu_from_Z_phi);
    VariablesSPS.add(sublead_pT_mu_from_Z_pT);
    VariablesSPS.add(sublead_pT_mu_from_Z_RAPIDITY);
    VariablesSPS.add(sublead_pT_mu_from_Z_eta);
    VariablesSPS.add(sublead_pT_mu_from_Z_phi);
    VariablesSPS.add(lead_pT_mu_from_upsi_pT);
    VariablesSPS.add(lead_pT_mu_from_upsi_RAPIDITY);
    VariablesSPS.add(lead_pT_mu_from_upsi_eta);
    VariablesSPS.add(lead_pT_mu_from_upsi_phi);
    VariablesSPS.add(sublead_pT_mu_from_upsi_pT);
    VariablesSPS.add(sublead_pT_mu_from_upsi_RAPIDITY);
    VariablesSPS.add(sublead_pT_mu_from_upsi_eta);
    VariablesSPS.add(sublead_pT_mu_from_upsi_phi);
    VariablesSPS.add(upsi_type);
    VariablesSPS.add(deltaPhi_upsiZ);
    VariablesSPS.add(deltaRAP_upsiZ);
 
  RooDataSet *data    = new RooDataSet("data", "data", tree_data, Variables);
  RooDataSet *mc      = new RooDataSet("mc",   "mc",   tree_mc,   Variables);
//  RooDataSet *mc_sps  = new RooDataSet("mc_sps",   "mc_sps",   tree_mc_sps,   Variables);

  RooRealVar weightToApply("weightToApply", "weightToApply", -1000, 10000);
    VariablesSPS.add(weightToApply);
  RooFormulaVar n1Func("n1Func","1.*@0",RooArgList(weightToApply));
  RooDataSet *mc0  = new RooDataSet("mc0",   "mc0",   tree_mc_sps,   VariablesSPS);
  RooRealVar* w = (RooRealVar*) mc0->addColumn(n1Func) ;
  RooDataSet mc1 (mc0->GetName(),mc0->GetTitle(),mc0,*mc0->get(),0,w->GetName()) ;
  RooDataSet *mc_sps = (RooDataSet*)mc1.reduce("1");
  
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
  
  RooDataSet *cut_mc1_sps = (RooDataSet*)mc_sps->reduce(SelectionCut1);

  std::cout << "Checkpoint 1" << std::endl;
  std::cout << cut_mc1->sumEntries() << std::endl;
  std::cout << cut_mc2->sumEntries() << std::endl;
  std::cout << cut_mc3->sumEntries() << std::endl;

  // plot part
  int Upsi_bins = 25; //20
  int Z_bins = 25; //20

 int bin_factor = 5;

 int upsi_pT_bins  = bin_factor*5;
 int upsi_RAP_bins = bin_factor*5;
 int upsi_phi_bins = bin_factor*5;
 int upsi_eta_bins = bin_factor*8;
 
 int Z_pT_bins  = bin_factor*5;
 int Z_RAP_bins = bin_factor*5;
 int Z_phi_bins = bin_factor*5;
 int Z_eta_bins = bin_factor*8;
 
 int lead_pT_mu_from_Z_pT_bins  = bin_factor*5;
 int lead_pT_mu_from_Z_RAP_bins = bin_factor*5;
 int lead_pT_mu_from_Z_eta_bins = bin_factor*5;
 int lead_pT_mu_from_Z_phi_bins = bin_factor*5;
 
 int sublead_pT_mu_from_Z_pT_bins  = bin_factor*5;
 int sublead_pT_mu_from_Z_RAP_bins = bin_factor*5;
 int sublead_pT_mu_from_Z_eta_bins = bin_factor*5;
 int sublead_pT_mu_from_Z_phi_bins = bin_factor*5;
 
 int lead_pT_mu_from_upsi_pT_bins  = bin_factor*5;
 int lead_pT_mu_from_upsi_RAP_bins = bin_factor*5;
 int lead_pT_mu_from_upsi_eta_bins = bin_factor*5;
 int lead_pT_mu_from_upsi_phi_bins = bin_factor*5;
 
 int sublead_pT_mu_from_upsi_pT_bins  = bin_factor*5;
 int sublead_pT_mu_from_upsi_RAP_bins = bin_factor*5;
 int sublead_pT_mu_from_upsi_eta_bins = bin_factor*5;
 int sublead_pT_mu_from_upsi_phi_bins = bin_factor*5;
 
 int deltaPhi_upsiZ_bins = bin_factor*6; //Changed from 3
 int deltaRAP_upsiZ_bins = bin_factor*6;

 //Make plots of the sWeighted distributions for each variable 
  
  //Won't cause a bug, nothing will be drawn because I added the lines in the wrong place, they need to go as shown in the example of c_weighted 
 
 //Canvases, draw on them and save them
 
 //upper and lower limits for Y axis of ratio plots defined here
 double lowerLimitYForRatioPlot = 0.;
 double upperLimitYForRatioPlot = 5.;

 double DPS_scale_factor = 1. / cut_mc1->sumEntries();
 double SPS_scale_factor = 1. / cut_mc1_sps->sumEntries();
  
//Canvas 1
  TCanvas *c_weighted = new TCanvas("c_weighted", "c_weighted", 1200, 400); c_weighted->Divide(3,1);
  
  c_weighted->cd(1); 
  TH1* tmp1 = cut_mc1->createHistogram("upsi_pT", upsi_pT_bins); //Draw and scale the upsi pT of upsilons that have been tagged as upsi 1
  TH1* tmp1SPS = cut_mc1_sps->createHistogram("upsi_pT", upsi_pT_bins); //Draw and scale the upsi pT of upsilons that have been tagged as upsi 1
  tmp1->SetLineColor(kRed);
  tmp1SPS->SetLineColor(kBlue);
  tmp1->SetMarkerSize(0);
  tmp1SPS->SetMarkerSize(0);
  tmp1->DrawNormalized("he");
  tmp1SPS->DrawNormalized("hesame");
 
  c_weighted->cd(2); 
  TH1* tmp2 = cut_mc2->createHistogram("upsi_pT", upsi_pT_bins);
  tmp2->SetLineColor(kRed);
  tmp2->SetMarkerSize(0);
  tmp2->DrawNormalized("he");
  
  c_weighted->cd(3); 
  TH1* tmp3= cut_mc3->createHistogram("upsi_pT", upsi_pT_bins);
  tmp3->SetLineColor(kRed);
  tmp3->SetMarkerSize(0);
  tmp3->DrawNormalized("he");
  
  c_weighted->SaveAs("c_weighted_upsi_pT_norm.pdf");

#ifdef makeplots 
//Canvas 2
  TCanvas *c_weighted2 = new TCanvas("c_weighted2", "c_weighted2", 1200,400); c_weighted2->Divide(3,1);
  
  c_weighted2->cd(1); 

  TH1* tmp4 = cut_mc1->createHistogram("Z_pT", Z_pT_bins);
  TH1* tmp4SPS = cut_mc1_sps->createHistogram("Z_pT", Z_pT_bins);
  tmp4->SetLineColor(kRed);
  tmp4SPS->SetLineColor(kBlue);
  tmp4->SetMarkerSize(0);
  tmp4SPS->SetMarkerSize(0);
  tmp4->SetMarkerSize(0);
  tmp4->DrawNormalized("he");
  tmp4SPS->DrawNormalized("hesame");
  
  c_weighted2->cd(2); 
  
  TH1* tmp5 = cut_mc2->createHistogram("Z_pT", Z_pT_bins);
  tmp5->SetLineColor(kRed);
  tmp5->SetMarkerSize(0);
  tmp5->DrawNormalized("he");
  
  c_weighted2->cd(3); 
  TH1* tmp6 = cut_mc3->createHistogram("Z_pT", Z_pT_bins);
  tmp6->SetLineColor(kRed);
  tmp6->SetMarkerSize(0);
  tmp6->DrawNormalized("he");
  
  c_weighted2->SaveAs("c_weighted_Z_pT_norm.pdf");
  
  //Canvas 3
  TCanvas *c_weighted3 = new TCanvas("c_weighted3", "c_weighted3", 1200, 400); c_weighted3->Divide(3,1);
  
  c_weighted3->cd(1); 

  TH1* tmp7 = cut_mc1->createHistogram("Z_RAPIDITY", Z_RAP_bins);
  TH1* tmp7SPS = cut_mc1_sps->createHistogram("Z_RAPIDITY", Z_RAP_bins);
  tmp7->SetLineColor(kRed);
  tmp7SPS->SetLineColor(kBlue);
  tmp7->SetMarkerSize(0);
  tmp7SPS->SetMarkerSize(0);
  tmp7SPS->DrawNormalized("he");
  tmp7->DrawNormalized("hesame");

  c_weighted3->cd(2); 
  TH1* tmp8 = cut_mc2->createHistogram("Z_RAPIDITY", Z_RAP_bins);
  tmp8->SetLineColor(kRed);
  tmp8->SetMarkerSize(0);
  tmp8->DrawNormalized("he");
  
  c_weighted3->cd(3); 
  TH1* tmp9 = cut_mc3->createHistogram("Z_RAPIDITY", Z_RAP_bins);
  tmp9->SetLineColor(kRed);
  tmp9->SetMarkerSize(0);
  tmp9->DrawNormalized("he");
  
  c_weighted3->SaveAs("c_weighted_Z_RAP_norm.pdf");
  
  //Canvas 4
  TCanvas *c_weighted4 = new TCanvas("c_weighted4", "c_weighted4", 1200, 400); c_weighted4->Divide(3,1);
  
  c_weighted4->cd(1);
  TH1* tmp10 = cut_mc1->createHistogram("Z_phi", Z_phi_bins);
  TH1* tmp10SPS = cut_mc1_sps->createHistogram("Z_phi", Z_phi_bins);
  tmp10->SetLineColor(kRed);
  tmp10SPS->SetLineColor(kBlue);
  tmp10->SetMarkerSize(0);
  tmp10SPS->SetMarkerSize(0);
  tmp10SPS->DrawNormalized("he");
  tmp10->DrawNormalized("hesame");
  
  c_weighted4->cd(2); 
  TH1* tmp11 = cut_mc2->createHistogram("Z_phi", Z_phi_bins);
  tmp11->SetLineColor(kRed);
  tmp11->SetMarkerSize(0);
  tmp11->DrawNormalized("he");
  
  c_weighted4->cd(3); 
  TH1* tmp12 = cut_mc3->createHistogram("Z_phi", Z_phi_bins);
  tmp12->SetLineColor(kRed);
  tmp12->SetMarkerSize(0);
  tmp12->DrawNormalized("he");
  
  c_weighted4->SaveAs("c_weighted_Z_phi_norm.pdf");
  
  //Canvas 5
  TCanvas *c_weighted5 = new TCanvas("c_weighted5", "c_weighted5", 1200, 400); c_weighted5->Divide(3,1);
  
  c_weighted5->cd(1); 
  
  TH1* tmp13 = cut_mc1->createHistogram("Z_eta", Z_eta_bins);
  TH1* tmp13SPS = cut_mc1_sps->createHistogram("Z_eta", Z_eta_bins);
  tmp13->SetLineColor(kRed);
  tmp13SPS->SetLineColor(kBlue);
  tmp13->SetMarkerSize(0);
  tmp13SPS->SetMarkerSize(0);
  tmp13SPS->DrawNormalized("he");
  tmp13->DrawNormalized("hesame");
 
  //Start here 14 April 2023
  c_weighted5->cd(2); 
 
  
  TH1* tmp14 = cut_mc2->createHistogram("Z_eta", Z_eta_bins);
  tmp14->SetLineColor(kRed);
  tmp14->SetMarkerSize(0);
  tmp14->DrawNormalized("he");
  
  c_weighted5->cd(3); 
  
  
  TH1* tmp15 = cut_mc3->createHistogram("Z_eta", Z_eta_bins);
  tmp15->SetLineColor(kRed);
  tmp15->SetMarkerSize(0);
  tmp15->DrawNormalized("he");
  c_weighted5->cd(3);
  
  c_weighted5->SaveAs("c_weighted_Z_eta_norm.pdf");
  
  //Canvas 6 
  TCanvas *c_weighted6 = new TCanvas("c_weighted6", "c_weighted6", 1200, 400); c_weighted6->Divide(3,1);
  
  c_weighted6->cd(1); 
  
  TH1* tmp16 = cut_mc1->createHistogram("upsi_RAPIDITY", upsi_RAP_bins);
  TH1* tmp16SPS = cut_mc1_sps->createHistogram("upsi_RAPIDITY", upsi_RAP_bins);
  tmp16->SetLineColor(kRed);
  tmp16SPS->SetLineColor(kBlue);
  tmp16->SetMarkerSize(0);
  tmp16SPS->SetMarkerSize(0);
  tmp16SPS->DrawNormalized("he");
  tmp16->DrawNormalized("hesame");
  
  c_weighted6->cd(2); 
  
  TH1* tmp17 = cut_mc2->createHistogram("upsi_RAPIDITY", upsi_RAP_bins);
  tmp17->SetLineColor(kRed);
  tmp17->SetMarkerSize(0);
  tmp17->DrawNormalized("he");
  c_weighted6->cd(2);
  
  
  c_weighted6->cd(3); 
  
  TH1* tmp18 = cut_mc3->createHistogram("upsi_RAPIDITY", upsi_RAP_bins);
  tmp18->SetLineColor(kRed);
  tmp18->SetMarkerSize(0);
  tmp18->DrawNormalized("he");
  c_weighted6->cd(3);
  
  c_weighted6->SaveAs("c_weighted_upsi_RAP_norm.pdf");
  
  //Canvas 7
  TCanvas *c_weighted7 = new TCanvas("c_weighted7", "c_weighted7", 1200, 400); c_weighted7->Divide(3,1);
  
  c_weighted7->cd(1); 
  
  TH1* tmp19 = cut_mc1->createHistogram("upsi_phi", upsi_phi_bins);
  TH1* tmp19SPS = cut_mc1_sps->createHistogram("upsi_phi", upsi_phi_bins);
  tmp19->SetLineColor(kRed);
  tmp19SPS->SetLineColor(kBlue);
  tmp19->SetMarkerSize(0);
  tmp19SPS->SetMarkerSize(0);
  tmp19SPS->DrawNormalized("he");
  tmp19->DrawNormalized("hesame");

  c_weighted7->cd(1);
  
  c_weighted7->cd(2); 
  
  TH1* tmp20 = cut_mc2->createHistogram("upsi_phi", upsi_phi_bins);
  tmp20->SetLineColor(kRed);
  tmp20->SetMarkerSize(0);
  tmp20->DrawNormalized("he");
  c_weighted7->cd(2);
  
  c_weighted7->cd(3); 
  
  TH1* tmp21 = cut_mc3->createHistogram("upsi_phi", upsi_phi_bins);
  tmp21->SetLineColor(kRed);
  tmp21->SetMarkerSize(0);
  tmp21->DrawNormalized("he");
  c_weighted7->cd(3);
  
  c_weighted7->SaveAs("c_weighted_upsi_phi_norm.pdf");
  
  //Canvas 8
  TCanvas *c_weighted8 = new TCanvas("c_weighted8", "c_weighted8", 1200, 400); c_weighted8->Divide(3,1);
  
  c_weighted8->cd(1); 
  
  TH1* tmp22 = cut_mc1->createHistogram("upsi_eta", upsi_eta_bins);
  TH1* tmp22SPS = cut_mc1_sps->createHistogram("upsi_eta", upsi_eta_bins);
  tmp22->SetLineColor(kRed);
  tmp22SPS->SetLineColor(kBlue);
  tmp22->SetMarkerSize(0);
  tmp22SPS->SetMarkerSize(0);
  tmp22SPS->DrawNormalized("he");
  tmp22->DrawNormalized("hesame");
 
 
  c_weighted8->cd(2); 
  
  TH1* tmp23 = cut_mc2->createHistogram("upsi_eta", upsi_eta_bins);
  tmp23->SetLineColor(kRed);
  tmp23->SetMarkerSize(0);
  tmp23->DrawNormalized("he");
  
  c_weighted8->cd(3); 
  TH1* tmp24 = cut_mc3->createHistogram("upsi_eta", upsi_eta_bins);
  tmp24->SetLineColor(kRed);
  tmp24->SetMarkerSize(0);
  tmp24->DrawNormalized("he");
  c_weighted8->cd(3);
  
  c_weighted8->SaveAs("c_weighted_upsi_eta_norm.pdf");
  
  //Canvas 9
  TCanvas *c_weighted9 = new TCanvas("c_weighted9", "c_weighted9", 1200, 400); c_weighted9->Divide(3,1);
  
  c_weighted9->cd(1); 
  
  TH1* tmp25 = cut_mc1->createHistogram("lead_pT_mu_from_Z_pT", lead_pT_mu_from_Z_pT_bins);
  TH1* tmp25SPS = cut_mc1_sps->createHistogram("lead_pT_mu_from_Z_pT", lead_pT_mu_from_Z_pT_bins);
  tmp25->SetLineColor(kRed);
  tmp25SPS->SetLineColor(kBlue);
  tmp25->SetMarkerSize(0);
  tmp25SPS->SetMarkerSize(0);
  std::cout << "CHECKPOINT K" << std::endl;
   tmp25->DrawNormalized("he");
  tmp25SPS->DrawNormalized("hesame");
 std::cout << tmp25->GetBinContent(4) << std::endl;
  
  c_weighted9->cd(2); 
  
  TH1* tmp26 = cut_mc2->createHistogram("lead_pT_mu_from_Z_pT", lead_pT_mu_from_Z_pT_bins);
  tmp26->SetLineColor(kRed);
  tmp26->SetMarkerSize(0);
  tmp26->DrawNormalized("he");
  c_weighted9->cd(2);
  
  c_weighted9->cd(3); 
  
  TH1* tmp27 = cut_mc3->createHistogram("lead_pT_mu_from_Z_pT", lead_pT_mu_from_Z_pT_bins);
  tmp27->SetLineColor(kRed);
  tmp27->SetMarkerSize(0);
  tmp27->DrawNormalized("he");
  c_weighted9->cd(3);
  
  c_weighted9->SaveAs("c_weighted_lead_pT_mu_from_Z_pT_norm.pdf");
  
  //Canvas 10
  TCanvas *c_weighted10 = new TCanvas("c_weighted10", "c_weighted10", 1200, 400); c_weighted10->Divide(3,1);
  
  c_weighted10->cd(1); 
  
  TH1* tmp28 = cut_mc1->createHistogram("lead_pT_mu_from_Z_RAPIDITY", lead_pT_mu_from_Z_RAP_bins);
  TH1* tmp28SPS = cut_mc1_sps->createHistogram("lead_pT_mu_from_Z_RAPIDITY", lead_pT_mu_from_Z_RAP_bins);
  tmp28->SetLineColor(kRed);
  tmp28SPS->SetLineColor(kBlue);
  tmp28->SetMarkerSize(0);
  tmp28SPS->SetMarkerSize(0);
  tmp28SPS->DrawNormalized("he");
  tmp28->DrawNormalized("hesame");

  c_weighted10->cd(2); 
  
  TH1* tmp29 = cut_mc2->createHistogram("lead_pT_mu_from_Z_RAPIDITY", lead_pT_mu_from_Z_RAP_bins);
  tmp29->SetLineColor(kRed);
  tmp29->SetMarkerSize(0);
  tmp29->DrawNormalized("he");
  c_weighted10->cd(2);
  
  c_weighted10->cd(3); 
 
  TH1* tmp30 = cut_mc3->createHistogram("lead_pT_mu_from_Z_RAPIDITY", lead_pT_mu_from_Z_RAP_bins);
  tmp30->SetLineColor(kRed);
  tmp30->SetMarkerSize(0);
  tmp30->DrawNormalized("he");
  c_weighted10->cd(3);
  
  c_weighted10->SaveAs("c_weighted_lead_pT_mu_from_Z_RAP_norm.pdf");
  
  //Canvas 11
  TCanvas *c_weighted11 = new TCanvas("c_weighted11", "c_weighted11", 1200, 400); c_weighted11->Divide(3,1);
  
  c_weighted11->cd(1); 
 
  TH1* tmp31 = cut_mc1->createHistogram("lead_pT_mu_from_Z_eta", lead_pT_mu_from_Z_eta_bins);
  TH1* tmp31SPS = cut_mc1_sps->createHistogram("lead_pT_mu_from_Z_eta", lead_pT_mu_from_Z_eta_bins);
  tmp31->SetLineColor(kRed);
  tmp31SPS->SetLineColor(kBlue);
  tmp31->SetMarkerSize(0);
  tmp31SPS->SetMarkerSize(0);
  tmp31SPS->DrawNormalized("he");
  tmp31->DrawNormalized("hesame");
 
  c_weighted11->cd(2); 
  
  TH1* tmp32 = cut_mc2->createHistogram("lead_pT_mu_from_Z_eta", lead_pT_mu_from_Z_eta_bins);
  tmp32->SetLineColor(kRed);
  tmp32->SetMarkerSize(0);
  tmp32->DrawNormalized("he");
  c_weighted11->cd(2);
  
  c_weighted11->cd(3); 
  
  TH1* tmp33 = cut_mc3->createHistogram("lead_pT_mu_from_Z_eta", lead_pT_mu_from_Z_eta_bins);
  tmp33->SetLineColor(kRed);
  tmp33->SetMarkerSize(0);
  tmp33->DrawNormalized("he");
  c_weighted11->cd(3);
  
  c_weighted11->SaveAs("c_weighted_lead_pT_mu_from_Z_eta_norm.pdf");
  
  
  //Canvas 12
  TCanvas *c_weighted12 = new TCanvas("c_weighted12", "c_weighted12", 1200, 400); c_weighted12->Divide(3,1);
  
  c_weighted12->cd(1); 
  
  TH1* tmp34 = cut_mc1->createHistogram("lead_pT_mu_from_Z_phi", lead_pT_mu_from_Z_phi_bins);
  TH1* tmp34SPS = cut_mc1_sps->createHistogram("lead_pT_mu_from_Z_phi", lead_pT_mu_from_Z_phi_bins);
  tmp34->SetLineColor(kRed);
  tmp34SPS->SetLineColor(kBlue);
  tmp34->SetMarkerSize(0);
  tmp34SPS->SetMarkerSize(0);
  tmp34SPS->DrawNormalized("he");
  tmp34->DrawNormalized("hesame");
 
  c_weighted12->cd(2); 
  
  TH1* tmp35 = cut_mc2->createHistogram("lead_pT_mu_from_Z_phi", lead_pT_mu_from_Z_phi_bins);
  tmp35->SetLineColor(kRed);
  tmp35->SetMarkerSize(0);
  tmp35->DrawNormalized("he");
  c_weighted12->cd(2);
  
  c_weighted12->cd(3); 
  
  TH1* tmp36 = cut_mc3->createHistogram("lead_pT_mu_from_Z_phi", lead_pT_mu_from_Z_phi_bins);
  tmp36->SetLineColor(kRed);
  tmp36->SetMarkerSize(0);
  tmp36->DrawNormalized("he");
  c_weighted12->cd(3);
  
  c_weighted12->SaveAs("c_weighted_lead_pT_mu_from_Z_phi_norm.pdf");
  
  //Canvas 13
  TCanvas *c_weighted13 = new TCanvas("c_weighted13", "c_weighted13", 1200, 400); c_weighted13->Divide(3,1);
  
  c_weighted13->cd(1); 
  
  TH1* tmp37 = cut_mc1->createHistogram("sublead_pT_mu_from_Z_pT", sublead_pT_mu_from_Z_pT_bins);
  TH1* tmp37SPS = cut_mc1_sps->createHistogram("sublead_pT_mu_from_Z_pT", sublead_pT_mu_from_Z_pT_bins);
  tmp37->SetLineColor(kRed);
  tmp37SPS->SetLineColor(kBlue);
  tmp37->SetMarkerSize(0);
  tmp37SPS->SetMarkerSize(0);
  tmp37->DrawNormalized("he");
  tmp37SPS->DrawNormalized("hesame");
 
  c_weighted13->cd(2); 
  
  TH1* tmp38 = cut_mc2->createHistogram("sublead_pT_mu_from_Z_pT", sublead_pT_mu_from_Z_pT_bins);
  tmp38->SetLineColor(kRed);
  tmp38->SetMarkerSize(0);
  tmp38->DrawNormalized("he");
  c_weighted13->cd(2);
  
  c_weighted13->cd(3); 
  
  TH1* tmp39 = cut_mc3->createHistogram("sublead_pT_mu_from_Z_pT", sublead_pT_mu_from_Z_pT_bins);
  tmp39->SetLineColor(kRed);
  tmp39->SetMarkerSize(0);
  tmp39->DrawNormalized("he");
  c_weighted13->cd(3);
  
  c_weighted13->SaveAs("c_weighted_sublead_pT_mu_from_Z_pT_norm.pdf");
  
  //Canvas 14
  TCanvas *c_weighted14 = new TCanvas("c_weighted14", "c_weighted14", 1200,400); c_weighted14->Divide(3,1);
  
  c_weighted14->cd(1); 
  
  TH1* tmp40 = cut_mc1->createHistogram("sublead_pT_mu_from_Z_RAPIDITY", sublead_pT_mu_from_Z_RAP_bins);
  TH1* tmp40SPS = cut_mc1_sps->createHistogram("sublead_pT_mu_from_Z_RAPIDITY", sublead_pT_mu_from_Z_RAP_bins);
  tmp40->SetLineColor(kRed);
  tmp40SPS->SetLineColor(kBlue);
  tmp40->SetMarkerSize(0);
  tmp40SPS->SetMarkerSize(0);
  tmp40SPS->DrawNormalized("he");
  tmp40->DrawNormalized("hesame");

  c_weighted14->cd(2); 
  
  TH1* tmp41 = cut_mc2->createHistogram("sublead_pT_mu_from_Z_RAPIDITY", sublead_pT_mu_from_Z_RAP_bins);
  tmp41->SetLineColor(kRed);
  tmp41->SetMarkerSize(0);
  tmp41->DrawNormalized("he");
  c_weighted14->cd(2);
  
  c_weighted14->cd(3); 
  
  TH1* tmp42 = cut_mc3->createHistogram("sublead_pT_mu_from_Z_RAPIDITY", sublead_pT_mu_from_Z_RAP_bins);
  tmp42->SetLineColor(kRed);
  tmp42->SetMarkerSize(0);
  tmp42->DrawNormalized("he");
  c_weighted14->cd(3);
  
  c_weighted14->SaveAs("c_weighted_sublead_pT_mu_from_Z_RAP_norm.pdf");
  
  //Canvas 15
  TCanvas *c_weighted15 = new TCanvas("c_weighted15", "c_weighted15", 1200, 400); c_weighted15->Divide(3,1);
  
  c_weighted15->cd(1);
  
  TH1* tmp43 = cut_mc1->createHistogram("sublead_pT_mu_from_Z_eta", sublead_pT_mu_from_Z_eta_bins);
  TH1* tmp43SPS = cut_mc1_sps->createHistogram("sublead_pT_mu_from_Z_eta", sublead_pT_mu_from_Z_eta_bins);
  tmp43->SetLineColor(kRed);
  tmp43SPS->SetLineColor(kBlue);
  tmp43->SetMarkerSize(0);
  tmp43SPS->SetMarkerSize(0);
  tmp43SPS->DrawNormalized("he");
  tmp43->DrawNormalized("hesame");
 
  c_weighted15->cd(2); 
  
  TH1* tmp44 = cut_mc2->createHistogram("sublead_pT_mu_from_Z_eta", sublead_pT_mu_from_Z_eta_bins);
  tmp44->SetLineColor(kRed);
  tmp44->SetMarkerSize(0);
  tmp44->DrawNormalized("he");
  c_weighted15->cd(2);
  
  c_weighted15->cd(3); 
  
  TH1* tmp45 = cut_mc3->createHistogram("sublead_pT_mu_from_Z_eta", sublead_pT_mu_from_Z_eta_bins);
  tmp45->SetLineColor(kRed);
  tmp45->SetMarkerSize(0);
  tmp45->DrawNormalized("he");
  c_weighted15->cd(3);
  
  c_weighted15->SaveAs("c_weighted_sublead_pT_mu_from_Z_eta_norm.pdf");
  
  //Canvas 16
  TCanvas *c_weighted16 = new TCanvas("c_weighted16", "c_weighted16", 1200,400); c_weighted16->Divide(3,1);
  
  c_weighted16->cd(1); 
  
  TH1* tmp46 = cut_mc1->createHistogram("sublead_pT_mu_from_Z_phi", sublead_pT_mu_from_Z_phi_bins);
  TH1* tmp46SPS = cut_mc1_sps->createHistogram("sublead_pT_mu_from_Z_phi", sublead_pT_mu_from_Z_phi_bins);
  tmp46->SetLineColor(kRed);
  tmp46SPS->SetLineColor(kBlue);
  tmp46->SetMarkerSize(0);
  tmp46SPS->SetMarkerSize(0);
  tmp46SPS->DrawNormalized("he");
  tmp46->DrawNormalized("hesame");
 
  c_weighted16->cd(2); 
  
  TH1* tmp47 = cut_mc2->createHistogram("sublead_pT_mu_from_Z_phi", sublead_pT_mu_from_Z_phi_bins);
  tmp47->SetLineColor(kRed);
  tmp47->SetMarkerSize(0);
  tmp47->DrawNormalized("he");
  c_weighted16->cd(2);
  
  c_weighted16->cd(3); 
  
  TH1* tmp48 = cut_mc3->createHistogram("sublead_pT_mu_from_Z_phi", sublead_pT_mu_from_Z_phi_bins);
  tmp48->SetLineColor(kRed);
  tmp48->SetMarkerSize(0);
  tmp48->DrawNormalized("he");
  c_weighted16->cd(3);
  
  c_weighted16->SaveAs("c_weighted_sublead_pT_mu_from_Z_phi_norm.pdf");
  
  //Canvas 17
  TCanvas *c_weighted17 = new TCanvas("c_weighted17", "c_weighted17", 1200, 400); c_weighted17->Divide(3,1);
  
  c_weighted17->cd(1); 
  
  TH1* tmp49 = cut_mc1->createHistogram("lead_pT_mu_from_upsi_pT", lead_pT_mu_from_upsi_pT_bins);
  TH1* tmp49SPS = cut_mc1_sps->createHistogram("lead_pT_mu_from_upsi_pT", lead_pT_mu_from_upsi_pT_bins);
  tmp49->SetLineColor(kRed);
  tmp49SPS->SetLineColor(kBlue);
  tmp49->SetMarkerSize(0);
  tmp49SPS->SetMarkerSize(0);
  tmp49->DrawNormalized("he");
  tmp49SPS->DrawNormalized("hesame");
  
  c_weighted17->cd(2); 
  
  TH1* tmp50 = cut_mc2->createHistogram("lead_pT_mu_from_upsi_pT", lead_pT_mu_from_upsi_pT_bins);
  tmp50->SetLineColor(kRed);
  tmp50->SetMarkerSize(0);
  tmp50->DrawNormalized("he");
  c_weighted17->cd(2);
  
  c_weighted17->cd(3); 
  
  TH1* tmp51 = cut_mc3->createHistogram("lead_pT_mu_from_upsi_pT", lead_pT_mu_from_upsi_pT_bins);
  tmp51->SetLineColor(kRed);
  tmp51->SetMarkerSize(0);
  tmp51->DrawNormalized("he");
  c_weighted17->cd(3);
  
  c_weighted17->SaveAs("c_weighted_lead_pT_mu_from_upsi_pT_norm.pdf");
  
  //Canvas 18
  TCanvas *c_weighted18 = new TCanvas("c_weighted18", "c_weighted18", 1200, 400); c_weighted18->Divide(3,1);
  
  c_weighted18->cd(1); 
  
  TH1* tmp52 = cut_mc1->createHistogram("lead_pT_mu_from_upsi_RAPIDITY", lead_pT_mu_from_upsi_RAP_bins);
  TH1* tmp52SPS = cut_mc1_sps->createHistogram("lead_pT_mu_from_upsi_RAPIDITY", lead_pT_mu_from_upsi_RAP_bins);
  tmp52->SetLineColor(kRed);
  tmp52SPS->SetLineColor(kBlue);
  tmp52->SetMarkerSize(0);
  tmp52SPS->SetMarkerSize(0);
  tmp52SPS->DrawNormalized("he");
  tmp52->DrawNormalized("hesame");

  c_weighted18->cd(2); 
  
  
  TH1* tmp53 = cut_mc2->createHistogram("lead_pT_mu_from_upsi_RAPIDITY", lead_pT_mu_from_upsi_RAP_bins);
  tmp53->SetLineColor(kRed);
  tmp53->SetMarkerSize(0);
  tmp53->DrawNormalized("he");
  c_weighted18->cd(2);
  
  c_weighted18->cd(3); 
  
  
  TH1* tmp54 = cut_mc3->createHistogram("lead_pT_mu_from_upsi_RAPIDITY", lead_pT_mu_from_upsi_RAP_bins);
  tmp54->SetLineColor(kRed);
  tmp54->SetMarkerSize(0);
  tmp54->DrawNormalized("he");
  c_weighted18->cd(3);
  
  c_weighted18->SaveAs("c_weighted_lead_pT_mu_from_upsi_RAP_norm.pdf");
  
  //Canvas 19
  TCanvas *c_weighted19 = new TCanvas("c_weighted19", "c_weighted19", 1200, 400); c_weighted19->Divide(3,1);
  
  c_weighted19->cd(1); 
  
  TH1* tmp55 = cut_mc1->createHistogram("lead_pT_mu_from_upsi_eta", lead_pT_mu_from_upsi_eta_bins);
  TH1* tmp55SPS = cut_mc1_sps->createHistogram("lead_pT_mu_from_upsi_eta", lead_pT_mu_from_upsi_eta_bins);
  tmp55->SetLineColor(kRed);
  tmp55SPS->SetLineColor(kBlue);
  tmp55->SetMarkerSize(0);
  tmp55SPS->SetMarkerSize(0);
  tmp55SPS->DrawNormalized("he");
  tmp55->DrawNormalized("hesame");

  c_weighted19->cd(2); 
  
  TH1* tmp56 = cut_mc2->createHistogram("lead_pT_mu_from_upsi_eta", lead_pT_mu_from_upsi_eta_bins);
  tmp56->SetLineColor(kRed);
  tmp56->SetMarkerSize(0);
  tmp56->DrawNormalized("he");
  c_weighted19->cd(2);
  
  c_weighted19->cd(3); 
  
  TH1* tmp57 = cut_mc3->createHistogram("lead_pT_mu_from_upsi_eta", lead_pT_mu_from_upsi_eta_bins);
  tmp57->SetLineColor(kRed);
  tmp57->SetMarkerSize(0);
  tmp57->DrawNormalized("he");
  c_weighted19->cd(3);
  
  c_weighted19->SaveAs("c_weighted_lead_pT_mu_from_upsi_eta_norm.pdf");
  
  //Canvas 20
  TCanvas *c_weighted20 = new TCanvas("c_weighted20", "c_weighted20", 1200,400); c_weighted20->Divide(3,1);
  
  c_weighted20->cd(1); 
  TH1* tmp58 = cut_mc1->createHistogram("lead_pT_mu_from_upsi_phi", lead_pT_mu_from_upsi_phi_bins);
  TH1* tmp58SPS = cut_mc1_sps->createHistogram("lead_pT_mu_from_upsi_phi", lead_pT_mu_from_upsi_phi_bins);
  tmp58->SetLineColor(kRed);
  tmp58SPS->SetLineColor(kBlue);
  tmp58->SetMarkerSize(0);
  tmp58SPS->SetMarkerSize(0);
  tmp58SPS->DrawNormalized("he");
  tmp58->DrawNormalized("hesame");

  c_weighted20->cd(2); 
  TH1* tmp59 = cut_mc2->createHistogram("lead_pT_mu_from_upsi_phi", lead_pT_mu_from_upsi_phi_bins);
  tmp59->SetLineColor(kRed);
  tmp59->SetMarkerSize(0);
  tmp59->DrawNormalized("he");
  
  c_weighted20->cd(3); 
  
  TH1* tmp60 = cut_mc3->createHistogram("lead_pT_mu_from_upsi_phi", lead_pT_mu_from_upsi_phi_bins);
  tmp60->SetLineColor(kRed);
  tmp60->SetMarkerSize(0);
  tmp60->DrawNormalized("he");
  c_weighted20->cd(3);
  
  c_weighted20->SaveAs("c_weighted_lead_pT_mu_from_upsi_phi_norm.pdf");
  
  //Canvas 21
  TCanvas *c_weighted21 = new TCanvas("c_weighted21", "c_weighted21", 1200,400); c_weighted21->Divide(3,1);
  
  c_weighted21->cd(1); 
  
  TH1* tmp61 = cut_mc1->createHistogram("sublead_pT_mu_from_upsi_pT", sublead_pT_mu_from_upsi_pT_bins);
  TH1* tmp61SPS = cut_mc1_sps->createHistogram("sublead_pT_mu_from_upsi_pT", sublead_pT_mu_from_upsi_pT_bins);
  tmp61->SetLineColor(kRed);
  tmp61SPS->SetLineColor(kBlue);
  tmp61->SetMarkerSize(0);
  tmp61SPS->SetMarkerSize(0);
   tmp61->DrawNormalized("he");
  tmp61SPS->DrawNormalized("hesame");
 
  c_weighted21->cd(2);
  
  TH1* tmp62 = cut_mc2->createHistogram("sublead_pT_mu_from_upsi_pT", sublead_pT_mu_from_upsi_pT_bins);
  tmp62->SetLineColor(kRed);
  tmp62->SetMarkerSize(0);
  tmp62->DrawNormalized("he");
  c_weighted21->cd(2);
  
  c_weighted21->cd(3); 
   
  TH1* tmp63 = cut_mc3->createHistogram("sublead_pT_mu_from_upsi_pT", sublead_pT_mu_from_upsi_pT_bins);
  tmp63->SetLineColor(kRed);
  tmp63->SetMarkerSize(0);
  tmp63->DrawNormalized("he");
  c_weighted21->cd(3);
  
  c_weighted21->SaveAs("c_weighted_sublead_pT_mu_from_upsi_pT_norm.pdf");
  
  //Canvas 22
  TCanvas *c_weighted22 = new TCanvas("c_weighted22", "c_weighted22", 1200,400); c_weighted22->Divide(3,1);
  
  c_weighted22->cd(1); 
   
  TH1* tmp64 = cut_mc1->createHistogram("sublead_pT_mu_from_upsi_RAPIDITY", sublead_pT_mu_from_upsi_RAP_bins);
  TH1* tmp64SPS = cut_mc1_sps->createHistogram("sublead_pT_mu_from_upsi_RAPIDITY", sublead_pT_mu_from_upsi_RAP_bins);
  tmp64->SetLineColor(kRed);
  tmp64SPS->SetLineColor(kBlue);
  tmp64->SetMarkerSize(0);
  tmp64SPS->SetMarkerSize(0);
  tmp64SPS->DrawNormalized("he");
  tmp64->DrawNormalized("hesame");

  c_weighted22->cd(2); 
   
  TH1* tmp65 = cut_mc2->createHistogram("sublead_pT_mu_from_upsi_RAPIDITY", sublead_pT_mu_from_upsi_RAP_bins);
  tmp65->SetLineColor(kRed);
  tmp65->SetMarkerSize(0);
  tmp65->DrawNormalized("he");
  c_weighted22->cd(2);
  
  c_weighted22->cd(3);
   
  TH1* tmp66 = cut_mc3->createHistogram("sublead_pT_mu_from_upsi_RAPIDITY", sublead_pT_mu_from_upsi_RAP_bins);
  tmp66->SetLineColor(kRed);
  tmp66->SetMarkerSize(0);
  tmp66->DrawNormalized("he");
  c_weighted22->cd(3);
  
  c_weighted22->SaveAs("c_weighted_sublead_pT_mu_from_upsi_RAP_norm.pdf");
  
  //Canvas 23
  TCanvas *c_weighted23 = new TCanvas("c_weighted23", "c_weighted23",1200, 400); c_weighted23->Divide(3,1);
  
  c_weighted23->cd(1); 
  
  TH1* tmp67 = cut_mc1->createHistogram("sublead_pT_mu_from_upsi_eta", sublead_pT_mu_from_upsi_eta_bins);
  TH1* tmp67SPS = cut_mc1_sps->createHistogram("sublead_pT_mu_from_upsi_eta", sublead_pT_mu_from_upsi_eta_bins);
  tmp67->SetLineColor(kRed);
  tmp67SPS->SetLineColor(kBlue);
  tmp67->SetMarkerSize(0);
  tmp67SPS->SetMarkerSize(0);
  tmp67SPS->DrawNormalized("he");
  tmp67->DrawNormalized("hesame");

  c_weighted23->cd(2); 
  
  TH1* tmp68 = cut_mc2->createHistogram("sublead_pT_mu_from_upsi_eta", sublead_pT_mu_from_upsi_eta_bins);
  tmp68->SetLineColor(kRed);
  tmp68->SetMarkerSize(0);
  tmp68->DrawNormalized("he");
  c_weighted23->cd(2);
  
  c_weighted23->cd(3); 
  
  TH1* tmp69 = cut_mc3->createHistogram("sublead_pT_mu_from_upsi_eta", sublead_pT_mu_from_upsi_eta_bins);
  tmp69->SetLineColor(kRed);
  tmp69->SetMarkerSize(0);
  tmp69->DrawNormalized("he");
  c_weighted23->cd(3);
  
  c_weighted23->SaveAs("c_weighted_sublead_pT_mu_from_upsi_eta_norm.pdf");
  
  //Canvas 24
  /* 
TCanvas *c_weighted24 = new TCanvas("c_weighted24", "c_weighted24", 1200, 400); c_weighted24->Divide(3,1);
  
  c_weighted24->cd(1); 
  
  TH1* tmp70 = cut_mc1->createHistogram("sublead_pT_mu_from_upsi_phi", sublead_pT_mu_from_upsi_phi_bins);
  tmp70->DrawNormalized("he");
  tmp70->SetLineColor(kRed);
  tmp70->SetMarkerSize(0);
  
  TH1* h70 = data_weighted_1->createHistogram("sublead_pT_mu_from_upsi_phi", sublead_pT_mu_from_upsi_phi_bins);
  h70->Sumw2();
  h70->Divide(tmp70);
  h70->DrawNormalized();
  h70->GetYaxis()->SetRangeUser(lowerLimitYForRatioPlot,upperLimitYForRatioPlot);
  h70->GetYaxis()->SetTitle("Data/MC");
  
  c_weighted24->cd(2); 
  
  TH1* tmp71 = cut_mc2->createHistogram("sublead_pT_mu_from_upsi_phi", sublead_pT_mu_from_upsi_phi_bins);
  tmp71->DrawNormalized("hesame");
  tmp71->SetLineColor(kRed);
  tmp71->SetMarkerSize(0);
  c_weighted24->cd(2);
  
  TH1* h71 = data_weighted_2->createHistogram("sublead_pT_mu_from_upsi_phi", sublead_pT_mu_from_upsi_phi_bins);
  h71->Sumw2();
  h71->Divide(tmp71);
  h71->DrawNormalized();
  h71->GetYaxis()->SetRangeUser(lowerLimitYForRatioPlot,upperLimitYForRatioPlot);
  h71->GetYaxis()->SetTitle("Data/MC");
  
  c_weighted24->cd(3); 
  
  TH1* tmp72 = cut_mc3->createHistogram("sublead_pT_mu_from_upsi_phi", sublead_pT_mu_from_upsi_phi_bins);
  tmp72->DrawNormalized("he");
  tmp72->SetLineColor(kRed);
  tmp72->SetMarkerSize(0);
  c_weighted24->cd(3);
  
  TH1* h72 = data_weighted_3->createHistogram("sublead_pT_mu_from_upsi_phi", sublead_pT_mu_from_upsi_phi_bins);
  h72->Sumw2();
  h72->Divide(tmp72);
  h72->DrawNormalized();
  h72->GetYaxis()->SetRangeUser(lowerLimitYForRatioPlot,upperLimitYForRatioPlot);
  h72->GetYaxis()->SetTitle("Data/MC");
  
  c_weighted24->SaveAs("c_weighted_sublead_pT_mu_from_upsi_phi_norm.pdf");
  
 */

#endif

  //Canvas 25
  TCanvas *c_weighted25 = new TCanvas("c_weighted25", "c_weighted25", 1200, 400); c_weighted25->Divide(3,1);
  
  c_weighted25->cd(1);
  TH1 * tmp73 = cut_mc1->createHistogram("deltaPhi_upsiZ", deltaPhi_upsiZ_bins);
  TH1 * tmp73SPS = cut_mc1_sps->createHistogram("deltaPhi_upsiZ", deltaPhi_upsiZ_bins);
  tmp73->SetLineColor(kRed);
  tmp73SPS->SetLineColor(kBlue);
  tmp73->SetMarkerSize(0);
  tmp73SPS->SetMarkerSize(0);
  tmp73SPS->DrawNormalized("he");
  tmp73->DrawNormalized("hesame");
  
  c_weighted25->cd(2);
  
  TH1* tmp74 = cut_mc2->createHistogram("deltaPhi_upsiZ", deltaPhi_upsiZ_bins);
  tmp74->SetLineColor(kRed);
  tmp74->SetMarkerSize(0);
  tmp74->DrawNormalized("he");
  
  c_weighted25->cd(3);
  TH1* tmp75 = cut_mc3->createHistogram("deltaPhi_upsiZ", deltaPhi_upsiZ_bins);
  tmp75->SetLineColor(kRed);
  tmp75->SetMarkerSize(0);
  tmp75->DrawNormalized("he");
  
 
  c_weighted25->SaveAs("c_weighted_deltaPhi_upsiZ_norm.pdf");
  
  //Canvas 26
  TCanvas *c_weighted26 = new TCanvas("c_weighted26", "c_weighted26", 1200, 400); c_weighted26->Divide(3,1);
  
  c_weighted26->cd(1);
  TH1* tmp76 = cut_mc1->createHistogram("deltaRAP_upsiZ", deltaRAP_upsiZ_bins);
  TH1* tmp76SPS = cut_mc1_sps->createHistogram("deltaRAP_upsiZ", deltaRAP_upsiZ_bins);
  tmp76->SetLineColor(kRed);
  tmp76SPS->SetLineColor(kBlue);
  tmp76->SetMarkerSize(0);
  tmp76SPS->SetMarkerSize(0);
  tmp76->DrawNormalized("he");
  tmp76SPS->DrawNormalized("hesame");
  
  c_weighted26->cd(2);
  TH1* tmp77 = cut_mc2->createHistogram("deltaRAP_upsiZ", deltaRAP_upsiZ_bins);
  tmp77->SetLineColor(kRed);
  tmp77->SetMarkerSize(0);
  tmp77->DrawNormalized("he");
  
  c_weighted26->cd(3);
  TH1* tmp78 = cut_mc3->createHistogram("deltaRAP_upsiZ", deltaRAP_upsiZ_bins);
  tmp78->SetLineColor(kRed);
  tmp78->SetMarkerSize(0);
  tmp78->DrawNormalized("he");
  
  c_weighted26->SaveAs("c_weighted_deltaRAP_upsiZ_norm.pdf");
  
}

