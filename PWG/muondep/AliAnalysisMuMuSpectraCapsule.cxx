/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

#include "AliAnalysisMuMuSpectraCapsule.h"

ClassImp(AliAnalysisMuMuSpectraCapsule)


#include "AliLog.h"
#include "TObject.h"
#include <TString.h>
#include <iostream>
#include <string>
#include "AliAnalysisMuMuResult.h"
#include "AliAnalysisMuMuJpsiResult.h"
#include "AliAnalysisMuMuBinning.h"
#include "TObjArray.h"
#include "TStyle.h"
#include "TH1F.h"
#include "TLine.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "TString.h"
#include "TPaveText.h"
#include "TPaveStats.h"
#include "TLatex.h"
#include "TList.h"
#include <fstream>
#include <string>

using std::cout;
using std::endl;
using std::ifstream;


//_____________________________________________________________________________
AliAnalysisMuMuSpectraCapsule::AliAnalysisMuMuSpectraCapsule() : TObject()
{
  /// Default ctor
}

//_____________________________________________________________________________
AliAnalysisMuMuSpectraCapsule::~AliAnalysisMuMuSpectraCapsule()
{
  // dtor
}

//_____________________________________________________________________________
Bool_t AliAnalysisMuMuSpectraCapsule::SetConstantFromExternFile(const char* file, Double_t* constantArray, const TString* spectraName)
{
  /// Set member constants depending on centrality bin from an ewternfile.
  /// If values are empty and can be obtained from a graph provided by the AliAnalysisMuMu Framework, this value is set by default
  /// For the PP capsule ine could be :
  ///
  /// #centrality Low  High  lumi.    lumi (stat)  lumi (syst. %)  Trigg  Trigg (local board)   Traj.err.(%)  MC Input (%)  Matching(%)  AccEff  dAccEff  NofJpsi Stat.Jpsi SystJpsi
  /// PP          PP   PP    0.         0.0         0.0            0.0    0.0                    04           03            01           0.0     0.0      0.0     0.0        0.0
  ///
  /// Note that for the PP case, the centrality limits are irrelevant
  ///
  /// For PbPb capsule :
  /// #centrality Low  High  <npart>    d<npart>    TAA           dTAA      sys.AP(%)           Traj.err.(%)   Trigg.err.(%) Matching(%)   AccEff   dAccEff  NofJpsi Stat.Jpsi SystJpsi
  /// V0M         00   10    359        31.2        23.4          0.351     2.00                04             03             01            0.1297   0.00040  105159  1693      488

    // Reset on fConstant
    for (int i = 0; i < 13; ++i) constantArray[i]=0.;
    Bool_t ok= kFALSE;
    AliDebug(1,Form("Reading from file %s",file));

    //________Open file
    ifstream infile(file,std::ios::in);
    TString line;
    TObjArray* lineArray;

    if (infile){
        AliDebug(1, " ==== opening file ==== ");
        // Loop until end of file is reached
        while(infile.eof()!=kTRUE){

            //read the line
            line.ReadLine(infile,kFALSE);
            if (line.BeginsWith("#"))continue;
            AliDebug(1,Form(" Read line : %s",line.Data()));

            // Put the line in a TObjArray
            lineArray = line.Tokenize(" ");


            // Select the good interval. Since interval is written in <binAsString>, just need them to match
            TString centrality   =  static_cast<TObjString*>(lineArray->At(0))->String().Data();
            TString intervalLow  =  TString::Format("%.2f",static_cast<TObjString*>(lineArray->At(1))->String().Atof());
            TString intervalHigh =  TString::Format("%.2f",static_cast<TObjString*>(lineArray->At(2))->String().Atof());
            AliDebug(1,Form(" --__--__-- interval low = %s",intervalLow.Data()));
            AliDebug(1,Form(" --__--__-- interval high = %s",intervalHigh.Data()));
            if (intervalLow.EqualTo("0.00")) intervalLow ="00.00";

            // Select the good interval for PbPb case. Since interval is written in <binAsString>, just need them to match
            if(spectraName->Contains(Form("%s",centrality.Data()))&& spectraName->Contains(Form("%s_%s",intervalLow.Data(),intervalHigh.Data())) && spectraName->Contains(Form("%s_%s",centrality.Data(),intervalLow.Data()))){
                AliDebug(1,Form(" spectraName = %s",spectraName->Data()));
                AliDebug(1,Form(" -- line selected -- "));
                ok = kTRUE;
                break;
            }
            // PP case
            else if(centrality.Contains("PP")){
                AliDebug(1,Form(" spectraName = %s",spectraName->Data()));
                AliDebug(1,Form(" -- line selected -- "));
                ok = kTRUE;
                break;
            }
            else continue;
        }
        infile.close();
        AliDebug(1, " ==== closing file ==== ");

        // Store the value
        for (int i =0 ; i<13 ; i++) {
            constantArray[i]= static_cast<TObjString*>(lineArray->At(i+3))->String().Atof();
        }
        return ok;
    }
    else return ok;
}


//_____________________________________________________________________________
void AliAnalysisMuMuSpectraCapsule::PrintNofWhat(const char* what) const
{
  /// Print whar number for each results on terminal.


  //Check point
  if(!GetSpectra() || strcmp(what,"")==1 )
    {
      AliError("No Spectra or no arguments given !");
      return ;
    }

  // Pointers to handle results and subresults and binning
  AliAnalysisMuMuResult    * result;
  AliAnalysisMuMuJpsiResult* subresult;
  AliAnalysisMuMuResult    * sr;
  AliAnalysisMuMuBinning   ::Range* r;

  TString swhat(what);
  // Array to store bins for the while loop
  TObjArray * bins=GetSpectra()->Binning()->CreateBinObjArray();// (intrinseque 'new')
  if (!bins)
  {
    AliError(Form("Cannot find bins"));
    return;
  }
  TCanvas *se[5];
  TCanvas *schi2[5];
  for(Int_t can=0; can<5; can++){
    se[can] = new TCanvas;
    schi2[can] = new TCanvas;
  }
  // se->Divide(1,bins->GetEntries());
  //Counters and Iterator for bin
  Int_t nofResult = 0;
  TIter nextBin(bins);
  nextBin.Reset();

  //Draw syst
  // TGraphErrors *gr[bins->GetEntries()];
  // TCanvas *cmult[bins->GetEntries()];
  // TCanvas *sysDis =new TCanvas("Systematics", Form("Systematics for %s",what),1200,800);
  // Int_t nx = 1;
  // Int_t ny = 1;
  // // if(bins->GetEntries()>1) {
  // //   ny = (bins->GetEntries()-1)/3+1;
  // //   nx = bins->GetEntries()/ny+1;
  // // }
  // AliInfo(Form("Dividing canvas into %d x %d (%d entries)",nx,ny,bins->GetEntries()));
  // sysDis->Divide(bins->GetEntries(),1);

  // Loop on bins
  //==============================================================================
  while ((r = static_cast<AliAnalysisMuMuBinning::Range*>(nextBin())))
  {
    // Make bin a MuMuResult
    result = GetSpectra()->GetResultForBin(*r);
    if (!result)
    {
      AliError(Form("Cannot find result "));
      return;
    }
    AliDebug(1, Form("result(%s) = %p ",result->GetName(),result));

    //plotting
    // gr[nofResult] = new TGraphErrors(result->SubResults()->GetEntries());


    Int_t nofSubResult = 0; // Counter for subresult
    TIter nextSubResult(result->SubResults());// Iterator for subresults
    nextSubResult.Reset();

    //Some variables
    TString  binAsString(r->AsString());// Usefull for the coming loop
     // To store subresults values
    Double_t subNofWhat[result->SubResults()->GetEntries()];
    Double_t subNofWhatStatError[result->SubResults()->GetEntries()];
    Double_t subChi2[result->SubResults()->GetEntries()];
    const char * srName[result->SubResults()->GetEntries()];

    TString srToExclude("");
    Int_t nofExcludedSr=0;

    cout << Form(" -_-_-_-_- %s_%s -_-_-_-_- ",binAsString.Data(),GetSpectraName().Data()) << endl;
    // Loop on subresults
    //==============================================================================
    while ((sr = static_cast<AliAnalysisMuMuResult*>(nextSubResult())))
    {
      // Get our final result
      subresult = static_cast<AliAnalysisMuMuJpsiResult*>(result->SubResult(Form("%s",sr->GetName())));
      if (!subresult)
      {
        AliError(Form("Cannot find subresult "));
        return;
      }
      AliDebug(1,Form("subresult(%s) = %p",sr->GetName(),subresult));

      //Get quantities
      Double_t NofWhat      = subresult->GetValue(what);
      Double_t NofWhatErrorStat = subresult->GetErrorStat(what);
      Double_t chi2 = subresult->GetValue("FitChi2PerNDF");
      subNofWhat[nofSubResult]          = NofWhat;
      subNofWhatStatError[nofSubResult] = NofWhatErrorStat;
      subChi2[nofSubResult] = chi2;
      srName[nofSubResult] =sr->GetName();
      // if(!TString(srName[nofSubResult]).Contains("_2.2_4.5")) srName[nofSubResult] = "";

               // else  srName[nofSubResult] = "";      if(TString(srName[nofSubResult]).Contains("CB2VWG2_BKGMV2POL2_2.2_4.5(Sig:2.2_4.5_SP1.1)")) srName[nofSubResult] = "CB2VWG2_pol2_2.2";
      // else if(TString(srName[nofSubResult]).Contains("CB2VWG2_BKGMV2POL2_2.3_4.6(Sig:2.2_4.5_SP1.1)")) srName[nofSubResult] = "CB2VWG2_pol2_2.3";
      // else if(TString(srName[nofSubResult]).Contains("CB2VWG2_BKGMV2POL2_2.4_4.7(Sig:2.2_4.5_SP1.1)")) srName[nofSubResult] = "CB2VWG2_pol2_2.4";
      // else if(TString(srName[nofSubResult]).Contains("CB2VWG2_BKGMV2POLEXP_2.2_4.5(Sig:2.2_4.5_SP1.1)")) srName[nofSubResult] = "CB2VWG2_polEx_2.2";
      // else if(TString(srName[nofSubResult]).Contains("CB2VWG2_BKGMV2POLEXP_2.3_4.6(Sig:2.2_4.5_SP1.1)")) srName[nofSubResult] = "CB2VWG2_polEx_2.3";
      // else if(TString(srName[nofSubResult]).Contains("CB2VWG2_BKGMV2POLEXP_2.4_4.7(Sig:2.2_4.5_SP1.1)")) srName[nofSubResult] = "CB2VWG2_polEx_2.4";
      // else if(TString(srName[nofSubResult]).Contains("CB2VWG2_BKGMV2POL4Cheb_2.2_4.5(Sig:2.2_4.5_SP1.1)")) srName[nofSubResult] = "CB2VWG2_pol4C_2.2";
      // else if(TString(srName[nofSubResult]).Contains("CB2VWG2_BKGMV2POL4Cheb_2.3_4.6(Sig:2.2_4.5_SP1.1)")) srName[nofSubResult] = "CB2VWG2_pol4C_2.3";
      // else if(TString(srName[nofSubResult]).Contains("CB2VWG2_BKGMV2POL4Cheb_2.4_4.7(Sig:2.2_4.5_SP1.1)")) srName[nofSubResult] = "CB2VWG2_pol4C_2.4";
      // else if(TString(srName[nofSubResult]).Contains("CB2POL2POL3_BKGMV2POL2_2.2_4.5(Sig:2.2_4.5_SP1.1)")) srName[nofSubResult] = "CB2POL_pol2_2.2";
      // else if(TString(srName[nofSubResult]).Contains("CB2POL2POL3_BKGMV2POL2_2.3_4.6(Sig:2.2_4.5_SP1.1)")) srName[nofSubResult] = "CB2POL_pol2_2.3";
      // else if(TString(srName[nofSubResult]).Contains("CB2POL2POL3_BKGMV2POL2_2.4_4.7(Sig:2.2_4.5_SP1.1)")) srName[nofSubResult] = "CB2POL_pol2_2.4";
      // else if(TString(srName[nofSubResult]).Contains("CB2POL2POL3_BKGMV2POLEXP_2.2_4.5(Sig:2.2_4.5_SP1.1)")) srName[nofSubResult] = "CB2POL_polEx_2.2";
      // else if(TString(srName[nofSubResult]).Contains("CB2POL2POL3_BKGMV2POLEXP_2.3_4.6(Sig:2.2_4.5_SP1.1)")) srName[nofSubResult] = "CB2POL_polEx_2.3";
      // else if(TString(srName[nofSubResult]).Contains("CB2POL2POL3_BKGMV2POLEXP_2.4_4.7(Sig:2.2_4.5_SP1.1)")) srName[nofSubResult] = "CB2POL_polEx_2.4";
      // else if(TString(srName[nofSubResult]).Contains("CB2POL2POL3_BKGMV2POL4Cheb_2.2_4.5(Sig:2.2_4.5_SP1.1)")) srName[nofSubResult] = "CB2POL_pol4C_2.2";
      // else if(TString(srName[nofSubResult]).Contains("CB2POL2POL3_BKGMV2POL4Cheb_2.3_4.6(Sig:2.2_4.5_SP1.1)")) srName[nofSubResult] = "CB2POL_pol4C_2.3";
      // else if(TString(srName[nofSubResult]).Contains("CB2POL2POL3_BKGMV2POL4Cheb_2.4_4.7(Sig:2.2_4.5_SP1.1)")) srName[nofSubResult] = "CB2POL_pol4C_2.4";
      // else if(TString(srName[nofSubResult]).Contains("NA60NEWVWG2_BKGMV2POL2_2.2_4.5(Sig:2.2_4.5_SP1.1)")) srName[nofSubResult] = "NA60VWG2_pol2_2.2";
      // else if(TString(srName[nofSubResult]).Contains("NA60NEWVWG2_BKGMV2POL2_2.3_4.6(Sig:2.2_4.5_SP1.1)")) srName[nofSubResult] = "NA60VWG2_pol2_2.3";
      // else if(TString(srName[nofSubResult]).Contains("NA60NEWVWG2_BKGMV2POL2_2.4_4.7(Sig:2.2_4.5_SP1.1)")) srName[nofSubResult] = "NA60VWG2_pol2_2.4";
      // else if(TString(srName[nofSubResult]).Contains("NA60NEWVWG2_BKGMV2POLEXP_2.2_4.5(Sig:2.2_4.5_SP1.1)")) srName[nofSubResult] = "NA60VWG2_polEx_2.2";
      // else if(TString(srName[nofSubResult]).Contains("NA60NEWVWG2_BKGMV2POLEXP_2.3_4.6(Sig:2.2_4.5_SP1.1)")) srName[nofSubResult] = "NA60VWG2_polEx_2.3";
      // else if(TString(srName[nofSubResult]).Contains("NA60NEWVWG2_BKGMV2POLEXP_2.4_4.7(Sig:2.2_4.5_SP1.1)")) srName[nofSubResult] = "NA60VWG2_polEx_2.4";
      // else if(TString(srName[nofSubResult]).Contains("NA60NEWVWG2_BKGMV2POL4Cheb_2.2_4.5(Sig:2.2_4.5_SP1.1)")) srName[nofSubResult] = "NA60VWG2_pol4C_2.2";
      // else if(TString(srName[nofSubResult]).Contains("NA60NEWVWG2_BKGMV2POL4Cheb_2.3_4.6(Sig:2.2_4.5_SP1.1)")) srName[nofSubResult] = "NA60VWG2_pol4C_2.3";
      // else if(TString(srName[nofSubResult]).Contains("NA60NEWVWG2_BKGMV2POL4Cheb_2.4_4.7(Sig:2.2_4.5_SP1.1)")) srName[nofSubResult] = "NA60VWG2_pol4C_2.4";
      // else if(TString(srName[nofSubResult]).Contains("NA60NEWPOL2POL3_BKGMV2POL2_2.2_4.5(Sig:2.2_4.5_SP1.1)")) srName[nofSubResult] = "NA60POL_pol2_2.2";
      // else if(TString(srName[nofSubResult]).Contains("NA60NEWPOL2POL3_BKGMV2POL2_2.3_4.6(Sig:2.2_4.5_SP1.1)")) srName[nofSubResult] = "NA60POL_pol2_2.3";
      // else if(TString(srName[nofSubResult]).Contains("NA60NEWPOL2POL3_BKGMV2POL2_2.4_4.7(Sig:2.2_4.5_SP1.1)")) srName[nofSubResult] = "NA60POL_pol2_2.4";
      // else if(TString(srName[nofSubResult]).Contains("NA60NEWPOL2POL3_BKGMV2POLEXP_2.2_4.5(Sig:2.2_4.5_SP1.1)")) srName[nofSubResult] = "NA60POL_polEx_2.2";
      // else if(TString(srName[nofSubResult]).Contains("NA60NEWPOL2POL3_BKGMV2POLEXP_2.3_4.6(Sig:2.2_4.5_SP1.1)")) srName[nofSubResult] = "NA60POL_polEx_2.3";
      // else if(TString(srName[nofSubResult]).Contains("NA60NEWPOL2POL3_BKGMV2POLEXP_2.4_4.7(Sig:2.2_4.5_SP1.1)")) srName[nofSubResult] = "NA60POL_polEx_2.4";
      // else if(TString(srName[nofSubResult]).Contains("NA60NEWPOL2POL3_BKGMV2POL4Cheb_2.2_4.5(Sig:2.2_4.5_SP1.1)")) srName[nofSubResult] = "NA60POL_pol4C_2.2";
      // else if(TString(srName[nofSubResult]).Contains("NA60NEWPOL2POL3_BKGMV2POL4Cheb_2.3_4.6(Sig:2.2_4.5_SP1.1)")) srName[nofSubResult] = "NA60POL_pol4C_2.3";
      // else if(TString(srName[nofSubResult]).Contains("NA60NEWPOL2POL3_BKGMV2POL4Cheb_2.4_4.7(Sig:2.2_4.5_SP1.1)")) srName[nofSubResult] = "NA60POL_pol4C_2.4";
      // else  srName[nofSubResult] = "";

      if(sr->GetValue("FitStatus")!=0 || subresult->GetValue("FitChi2PerNDF")>2.){
      // || (swhat.Contains("v2") && (sr->GetValue("<v2>JPsi")>0.1|| sr->GetValue("<v2>JPsi") < -0.1))){
        srToExclude += Form("%s,",sr->GetName());
        nofExcludedSr++;
        nofSubResult++;
        continue;
      }

      //Output messages
      AliDebug(0,Form(" -------- "));

      if(swhat.Contains("v2"))AliInfo(Form(" -- subresult %s :  %.4f +/- %.4f, FitStatus :%.0f, FitChi2PerNDF :%.1f",sr->GetName(),NofWhat,NofWhatErrorStat,sr->GetValue("FitStatus"),subresult->GetValue("FitChi2PerNDF")));
      else AliInfo(Form(" -- subresult %s :  %.0f +/- %.0f, FitStatus :%.0f, FitChi2PerNDF :%.1f ",sr->GetName(),NofWhat,NofWhatErrorStat,sr->GetValue("FitStatus"),subresult->GetValue("FitChi2PerNDF")));
      nofSubResult++;
    }
    AliInfo(Form("%d excluded fits : %s",nofExcludedSr,srToExclude.Data()));

    result->Exclude(srToExclude);
    AliInfo(Form(" -------- "));
    if(swhat.Contains("v2")) AliInfo(Form(" ------ Mean :  %.4f +/- %.4f (%.1f %%) +/- %.4f (%.1f %%) ------ \n",
      result->GetValue(what),result->GetErrorStat(what),100*result->GetErrorStat(what)/result->GetValue(what),result->GetRMS(what),100*result->GetRMS(what)/result->GetValue(what)));
    else AliInfo(Form(" ------ Mean :  %.1f +/- %.1f (%.1f %%) +/- %.1f (%.1f %%) ------ \n",
      result->GetValue(what),result->GetErrorStat(what),100*result->GetErrorStat(what)/result->GetValue(what),result->GetRMS(what),100*result->GetRMS(what)/result->GetValue(what)));
    // AliDebug(0,"");

    // Plot the histograms
    TH1F * h_test = new TH1F(Form("%s_%s",what,r->AsString().Data()),Form("%s_%s",what,r->AsString().Data()),result->SubResults()->GetEntries(),0,result->SubResults()->GetEntries());
    TH1F * h_chi2 = new TH1F(Form("Chi2_%s_%s",what,r->AsString().Data()),Form("Chi2_%s_%s",what,r->AsString().Data()),result->SubResults()->GetEntries(),0,result->SubResults()->GetEntries());
    if(swhat.Contains("v2"))
    {
      h_test->SetMarkerColor(kAzure-2);
      h_test->SetLineColor(kAzure-2);
    }
    else {
      h_test->SetMarkerColor(kRed+1);
      h_test->SetLineColor(kRed+1);
    }
    for (int i = 0; i < result->SubResults()->GetEntries(); ++i)
    {
    //     if(!TString(srName[nofSubResult]).Contains("CB2VWG2_BKGMV2POL2_2.2_4.5(Sig:2.2_4.5_SP1.1)")&&!TString(srName[nofSubResult]).Contains("CB2VWG2_BKGMV2POL2_2.3_4.6(Sig:2.3_4.6_SP1.1)")&&!TString(srName[nofSubResult]).Contains("CB2VWG2_BKGMV2POL2_2.4_4.7(Sig:2.4_4.7_SP1.1)")
    //     &&!TString(srName[nofSubResult]).Contains("CB2VWG2_BKGMV2POLEXP_2.2_4.5(Sig:2.2_4.5_SP1.1)")&&!TString(srName[nofSubResult]).Contains("CB2VWG2_BKGMV2POLEXP_2.3_4.6(Sig:2.3_4.6_SP1.1)")&&!TString(srName[nofSubResult]).Contains("CB2VWG2_BKGMV2POLEXP_2.4_4.7(Sig:2.4_4.7_SP1.1)")
    //     &&!TString(srName[nofSubResult]).Contains("CB2VWG2_BKGMV2POL4Cheb_2.2_4.5(Sig:2.2_4.5_SP1.1)")&&!TString(srName[nofSubResult]).Contains("CB2VWG2_BKGMV2POL4Cheb_2.3_4.6(Sig:2.3_4.6_SP1.1)")&&!TString(srName[nofSubResult]).Contains("CB2VWG2_BKGMV2POL4Cheb_2.4_4.7(Sig:2.4_4.7_SP1.1)")
    //     &&!TString(srName[nofSubResult]).Contains("CB2POL2POL3_BKGMV2POL2_2.2_4.5(Sig:2.2_4.5_SP1.1)"&&!TString(srName[nofSubResult]).Contains("CB2POL2POL3_BKGMV2POL2_2.3_4.6(Sig:2.3_4.6_SP1.1)"&&!TString(srName[nofSubResult]).Contains("CB2POL2POL3_BKGMV2POL2_2.4_4.7(Sig:2.4_4.7_SP1.1)"
    //     &&!TString(srName[nofSubResult]).Contains("CB2POL2POL3_BKGMV2POLEXP_2.2_4.5(Sig:2.2_4.5_SP1.1)")&&!TString(srName[nofSubResult]).Contains("CB2POL2POL3_BKGMV2POLEXP_2.3_4.6(Sig:2.3_4.6_SP1.1)")&&!TString(srName[nofSubResult]).Contains("CB2POL2POL3_BKGMV2POLEXP_2.4_4.7(Sig:2.4_4.7_SP1.1)")
    //     &&!TString(srName[nofSubResult]).Contains("CB2POL2POL3_BKGMV2POL4Cheb_2.2_4.5(Sig:2.2_4.5_SP1.1)")&&!TString(srName[nofSubResult]).Contains("CB2POL2POL3_BKGMV2POL4Cheb_2.3_4.6(Sig:2.3_4.6_SP1.1)")&&!TString(srName[nofSubResult]).Contains("CB2POL2POL3_BKGMV2POL4Cheb_2.4_4.7(Sig:2.4_4.7_SP1.1)")
    //     &&!TString(srName[nofSubResult]).Contains("NA60NEWVWG2_BKGMV2POL2_2.2_4.5(Sig:2.2_4.5_SP1.1)")&&!TString(srName[nofSubResult]).Contains("NA60NEWVWG2_BKGMV2POL2_2.3_4.6(Sig:2.3_4.6_SP1.1)")&&!TString(srName[nofSubResult]).Contains("NA60NEWVWG2_BKGMV2POL2_2.4_4.7(Sig:2.4_4.7_SP1.1)")
    //     &&!TString(srName[nofSubResult]).Contains("NA60NEWVWG2_BKGMV2POLEXP_2.2_4.5(Sig:2.2_4.5_SP1.1)")) &&!TString(srName[nofSubResult]).Contains("NA60NEWVWG2_BKGMV2POLEXP_2.3_4.6(Sig:2.3_4.6_SP1.1)")) &&!TString(srName[nofSubResult]).Contains("NA60NEWVWG2_BKGMV2POLEXP_2.4_4.7(Sig:2.4_4.7_SP1.1)")
    //     &&!TString(srName[nofSubResult]).Contains("NA60NEWVWG2_BKGMV2POL4Cheb_2.2_4.5(Sig:2.2_4.5_SP1.1)")) &&!TString(srName[nofSubResult]).Contains("NA60NEWVWG2_BKGMV2POL4Cheb_2.3_4.6(Sig:2.3_4.6_SP1.1)") &&!TString(srName[nofSubResult]).Contains("NA60NEWVWG2_BKGMV2POL4Cheb_2.4_4.7(Sig:2.4_4.7_SP1.1)")
    //     &&!TString(srName[nofSubResult]).Contains("NA60NEWPOL2POL3_BKGMV2POL2_2.2_4.5(Sig:2.2_4.5_SP1.1)")&&!TString(srName[nofSubResult]).Contains("NA60NEWPOL2POL3_BKGMV2POL2_2.3_4.6(Sig:2.3_4.6_SP1.1)")&&!TString(srName[nofSubResult]).Contains("NA60NEWPOL2POL3_BKGMV2POL2_2.4_4.7(Sig:2.4_4.7_SP1.1)")
    //     &&!TString(srName[nofSubResult]).Contains("NA60NEWPOL2POL3_BKGMV2POLEXP_2.2_4.5(Sig:2.2_4.5_SP1.1)")&&!TString(srName[nofSubResult]).Contains("NA60NEWPOL2POL3_BKGMV2POLEXP_2.3_4.6(Sig:2.3_4.6_SP1.1)")&&!TString(srName[nofSubResult]).Contains("NA60NEWPOL2POL3_BKGMV2POLEXP_2.4_4.7(Sig:2.4_4.7_SP1.1)")
    //     &&!TString(srName[nofSubResult]).Contains("NA60NEWPOL2POL3_BKGMV2POL4Cheb_2.2_4.5(Sig:2.2_4.5_SP1.1)")&&!TString(srName[nofSubResult]).Contains("NA60NEWPOL2POL3_BKGMV2POL4Cheb_2.3_4.6(Sig:2.3_4.6_SP1.1)")&&!TString(srName[nofSubResult]).Contains("NA60NEWPOL2POL3_BKGMV2POL4Cheb_2.4_4.7(Sig:2.4_4.7_SP1.1)"))
    //      {continue;}
        // the histo we plot
        h_test->SetBinContent(i+1,subNofWhat[i]);
        h_test->SetBinError(i+1,subNofWhatStatError[i]);

        // Here we change the label names
        // if(nofResult==(bins->GetEntries()-1))
        h_test->GetXaxis()->SetBinLabel(i+1,Form("%s",srName[i]));
        h_chi2->GetXaxis()->SetBinLabel(i+1,Form("%s",srName[i]));

        h_chi2->SetBinContent(i+1,subChi2[i]);
    }

    // --- Here we draw ---
    se[nofResult]->cd();
    // se->cd(nofResult+1);
    se[nofResult]->SetBottomMargin(4.);
    h_test->GetYaxis()->SetTitle(what);
    h_test->SetStats(0);
    h_test->DrawCopy();

    TLine *line1 = new TLine(0,result->GetValue(what),result->SubResults()->GetEntries(),result->GetValue(what));
    line1->SetLineColor(kGray+3);
    line1->SetLineWidth(3);

    TLine *line2 = new TLine(0,result->GetValue(what)-result->GetRMS(what),result->SubResults()->GetEntries(),result->GetValue(what)-result->GetRMS(what));
    line2->SetLineColor(kGray+3);
    line2->SetLineWidth(3);
    line2->SetLineStyle(8);

    TLine *line3 = new TLine(0,result->GetValue(what)+result->GetRMS(what),result->SubResults()->GetEntries(),result->GetValue(what)+result->GetRMS(what));
    line3->SetLineColor(kGray+3);
    line3->SetLineWidth(3);
    line3->SetLineStyle(8);
    line1->Draw("same");
    line2->Draw("same");
    line3->Draw("same");

    TPaveText* pt = new TPaveText(0.16,0.8,0.47,0.85,"nbNDC");
    pt->AddText(Form(" [%0.f < #it{p}_{T} < %0.f GeV/#it{c}] Excluded fits : %d / %d",r->Xmin(), r->Xmax(), nofExcludedSr,nofSubResult));
    pt->Draw("same");
    // m0.SetLineStyle(1);
    // s1.SetLineStyle(5);
    // s2.SetLineStyle(8);
    // TPaveText ex(0.16,0.8,0.47,0.85,"nbNDC");
    // ex.SetTextAlign(11);
    // ex.SetBorderSize(0);
    // ex.SetFillStyle(0);
    // ex.SetTextColor(kGray+3);
    // ex.SetTextFont(42);
    // // ex.SetTextSize(gStyle->GetTextSize()*0.9);

    schi2[nofResult]->cd();
    h_chi2->DrawCopy("");
    nofResult++;
  }

  // delete result;
  // delete subresult;
  // delete sr;
  // delete r;
  // delete bins;
  // delete gr;
  // delete ex;
}
