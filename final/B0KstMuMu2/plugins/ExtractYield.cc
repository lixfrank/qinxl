// ################################################################################
// # Program to perform the full angular analysis of the decay B0 --> K*0 mu+ mu- #
// ################################################################################

#include <TROOT.h>
#include <TApplication.h>
#include <TSystem.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <Riostream.h>
#include <TH3.h>
#include <TEnv.h>
#include <TSystem.h>
#include <TF2.h>
#include <TStyle.h>
#include <TFile.h>
#include <TTree.h>
#include <TGaxis.h>
#include <TLegend.h>
#include <TPaveText.h>
#include <TLatex.h>
#include <TCutG.h>
#include <Math/Functor.h>

#include <TMinuit.h>
#include <RooMinimizer.h>
#include <RooNLLVar.h>
#include <RooRealVar.h>
#include <RooGaussian.h>
#include <RooAbsPdf.h>
#include <RooHistPdf.h>
#include <RooEffProd.h>
#include <RooAddPdf.h>
#include <RooDataSet.h>
#include <RooGenericPdf.h>
#include <RooPlot.h>
#include <RooHistFunc.h>
#include <RooArgSet.h>
#include <RooFitResult.h>
#include <RooPolynomial.h>
#include <RooMCStudy.h>
#include <RooMinuit.h>
#include <RooNLLVar.h>
#include <RooWorkspace.h>
#include <RooConstVar.h>
#include <RooRandom.h>
#include <RooDataHist.h>
#include <RooFunctorBinding.h>
#include <RooStats/RooStatsUtils.h>

#include <TBenchmark.h>
#include <time.h>
#include <ctime>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <utility>
#include <sstream>

#include "DecayRate.h"
#include "PdfSigAng.h"
#include "PdfRT.h"
#include "PdfWT.h"
#include "ReadParameters.h"
#include "Utils.h"
#include "B0KstMuMuTreeContent.h"
using std::cout;
using std::endl;
using std::string;
using std::stringstream;
using std::vector;
using std::ios_base;
using std::pair;
using std::fstream;
using std::ifstream;
using std::make_pair;
using namespace RooFit;


// ####################
// # Global constants #
// ####################
#define NBINS        20
#define MULTYIELD     1. // Multiplication factor to the number of entry in toy-MC
#define NCOEFFPOLYBKG 5  // Maximum number of coefficients (= degree) of the polynomial describing the combinatorial bkg

#define _USE_MATH_DEFINES

// ##########################################
// # Internal flags to control the workflow #
// ##########################################
#define PLOT          true  //[true= plot the results]
#define SAVEPLOT      true   //2015-12-17
#define MINIMIZER     "Minuit" // Minimizer type for 3D MODEL actual fit ["Minuit"; "Minuit2"]

#define doTransf       false    // limits
#define ROOMINIMIZER   false     // if true: RooMinimizer, if false: Roofit
// ##################
// # External files #
// ##################
#define PARAMETERFILEIN  "/python/ParameterFile.txt"
#define PARAMETERFILEOUT "ParameterFileOut.txt"


// ############################################
// # Global variables from configuration file #
// ############################################
double LUMI;

string CTRLfitWRKflow;
string ParameterFILE;

vector<double> q2Bins;
vector<double> cosThetaKBins;
vector<double> cosThetaLBins;
vector<double> phiBins;

vector<vector<string>*>             fitParam;    // Vector containing the pointers to the vectors containing the starting values for the fit
vector<vector<unsigned int>*>       configParam; // Vector containing the pointers to the vectors containing the configuration parameters for the fit

// ####################
// # Global variables #
// ####################
TTree* theTree;
Utils* Utility;
B0KstMuMuTreeContent* NTuple;

double* q2BinsHisto;

// ####################################
// # Useful variables from the NTuple #
// ####################################
RooDataSet* toy;
RooDataSet* toy2;
RooDataSet* SingleCandNTuple_JPsi;
RooDataSet* SingleCandNTuple_PsiP;
RooDataSet* SingleCandNTuple_RejectPsi;
RooDataSet* SingleCandNTuple;
RooRealVar* B0mass;
RooRealVar* mumuMass;
RooRealVar* mumuMassE;
RooRealVar* ctK;                //CosThetaKArb
RooRealVar* ctL;                //CosThetaMuArb
RooRealVar* phi;
RooRealVar* tagB0;
RooRealVar* genSignal;


// #################################
// # Variables and pdf for the fit #
// #################################

// ##################
// # Signal B0 mass #
// ##################

// #################
// # Signal angles #
// #################
RooRealVar* Fl;

RooRealVar* P4p;
RooRealVar* P5p;
RooRealVar* P6p;
RooRealVar* P8p;
RooRealVar* P1;
RooRealVar* P2;
RooRealVar* P3;              

RooRealVar* S3S;
RooRealVar* S4S;
RooRealVar* S5S;
RooRealVar* AFBS;
RooRealVar* S7S;
RooRealVar* S8S;
RooRealVar* S9S;

RooAbsPdf*  AngleGood;
RooEffProd*  AngleS;

// ####################
// # Total signal pdf #
// ####################
RooAbsPdf* Signal;
RooAbsPdf* SignalT;

// ####################################
// # Combinatorial background B0 mass #
// ####################################

// #########################
// # Mistag signal B0 mass #
// #########################
RooRealVar* fracMisTag;
RooAbsPdf*  AngleMisTag;
// ##############################
// # Peaking background B0 mass #
// ##############################

// #####################
// # Background angles #
// #####################

// ########################
// # Total background pdf #
// ########################
RooRealVar* nMisTagFrac;

// ##################################
// # Total pdf for B0 --> K*0 mu mu #
// ##################################
RooAbsPdf* TotalPDFRejectPsi;

// ##############################################
// # Vector containing the Gaussian constraints #
// ##############################################
RooArgSet vecConstr;


// #######################
// # Function Definition #
// #######################

struct MyProdPdf
{
public:
  MyProdPdf (RooAbsPdf& pdf1, RooAbsPdf& pdf2) : _pdf1(pdf1), _pdf2(pdf2)
  {
    const RooArgSet* allvar1 = pdf1.getVariables();
    const RooArgSet* allvar2 = pdf2.getVariables();

    _vars.add(*allvar1);
//    _vars.add(*allvar2,kFALSE);

    delete allvar1;
    delete allvar2;
  }

  int ndim ()
  {
    return _vars.getSize();
  }

  const RooArgList& vars() const
  {
    return _vars;
  }
  double operator() (const double* v)
  {
    for (int i = 0; i < ndim(); ++i) ((RooRealVar&)_vars[i]).setVal(v[i]);

   return _pdf1.getVal() * _pdf2.getVal();
  }


private:
  RooAbsPdf& _pdf1;
  RooAbsPdf& _pdf2;
  RooArgList _vars;

};



bool CheckGoodFit               (RooFitResult* fitResult, TPaveText* paveText = NULL);
RooRealVar* GetVar              (RooAbsPdf* pdf, string varName);
bool GetValueAndErrors          (RooAbsPdf* pdf, string varName, stringstream* myString, double& val, double& errLo, double& errHi);
void SetValueAndErrors          (RooAbsPdf* pdf, string varName, double multi, stringstream* myString, double* val, double* errLo, double* errHi);
void PrintVariables             (RooArgSet* setVar, string type);
void ClearVars                  (RooArgSet* vecVars);

void SetStyle                   ();
string MakeName                 (const RooDataSet* data, int ID);
void DrawString                 (double Lumi, RooPlot* myFrame = NULL);
void AddGaussConstraint         (RooArgSet* vecConstr, RooAbsPdf* pdf, string varName);

double PrintFitResults    (RooAbsPdf* pdf, RooFitResult* fitResult);

void MakeDatasets               (B0KstMuMuTreeContent* NTuple, unsigned int FitType, int SampleType);

//###############
// 3-D models   #
//###############

void Instantiate3AnglesFit     (RooAbsPdf** TotalPDF,
				RooRealVar* y, RooRealVar* z,RooRealVar* p,
                                unsigned int FitType,
				unsigned int q2BinIndx);


RooFitResult* Make3AnglesFit (RooDataSet* dataSet, RooAbsPdf** TotalPDF, RooRealVar* y, RooRealVar* z,RooRealVar* p,  unsigned int FitType, RooArgSet* vecConstr, TCanvas* Canv, int specBin);

void Iterative3AnglesFitq2Bins (RooDataSet* dataSet,
                                    RooRealVar* y, RooRealVar* z,RooRealVar* p,
                                    int specBin,
                                    unsigned int FitType,
                                    vector<double>* q2Bins,
                                    RooArgSet* vecConstr);


// ###########################
// # Function Implementation #
// ###########################
bool CheckGoodFit (RooFitResult* fitResult, TPaveText* paveText)
// ####################################################
// # Covariance matrix quality:                       #
// # -1 : "Unknown, matrix was externally provided"   #
// #  0 : "Not calculated at all"                     #
// #  1 : "Approximation only, not accurate"          #
// #  2 : "Full matrix, but forced positive-definite" #
// #  3 : "Full, accurate covariance matrix"          #
// ####################################################
{
  if (fitResult != NULL)
    {
      if ((fitResult->covQual() == 3) && (fitResult->status() == 0))
	{
          cout << "status from fitresult = " <<  fitResult->status() << endl;
	  if (paveText != NULL) paveText->AddText("Fit status: GOOD");
	  return true;
	}
      else
	{
	  if (paveText != NULL) paveText->AddText("Fit status: BAD");
	  return false;
	}
    }


  return false;
}


RooRealVar* GetVar (RooAbsPdf* pdf, string varName)
{
  return (RooRealVar*)(pdf->getVariables()->find(varName.c_str()));
}


bool GetValueAndErrors (RooAbsPdf* pdf, string varName, stringstream* myString, double& val, double& errLo, double& errHi)
{
  if (GetVar(pdf,varName.c_str()) != NULL)
    {
      val   = GetVar(pdf,varName.c_str())->getVal();
      errLo = GetVar(pdf,varName.c_str())->getErrorLo();
      errHi = GetVar(pdf,varName.c_str())->getErrorHi();

      (*myString) << val << "   " << errHi << "   " << errLo << "   ";
      return true;
    }
  else (*myString) << "0.0   0.0   0.0";


  return false;
}


void SetValueAndErrors (RooAbsPdf* pdf, string varName, double multi, stringstream* myString, double* val, double* errLo, double* errHi)
// #############################################################################
// # If the error is an empty string --> setAsymError = -1/+1 and setError = 1 #
// # If both errLo and errHi are 0.0 --> setAsymError = -1/+1 and setError = 1 #
// #############################################################################
{
  string tmpStr;


  if (myString->str().empty() == false)
    {
      tmpStr.clear();
      (*myString) >> tmpStr;
      *val = atof(tmpStr.c_str()) * multi;

      tmpStr.clear();
      (*myString) >> tmpStr;
      if (tmpStr.empty() == true) *errLo = -1.0;
      else *errLo = atof(tmpStr.c_str()) * multi;

      tmpStr.clear();
      (*myString) >> tmpStr;
      if (tmpStr.empty() == true) *errHi = 1.0;
      else *errHi = atof(tmpStr.c_str()) * multi;
    }

  if ((pdf != NULL) && (GetVar(pdf,varName) != NULL))
    {
      pdf->getVariables()->setRealValue(varName.c_str(),*val);
      if ((*errLo == 0.0) && (*errHi == 0.0))
      	{
	  GetVar(pdf,varName)->setAsymError(0.0,0.0);
	  GetVar(pdf,varName)->setError(0.0);
      	}
      else
      	{
      	  GetVar(pdf,varName)->setAsymError(*errLo,*errHi);
	  GetVar(pdf,varName)->setError((*errHi - *errLo) / 2.);
      	}
    }
}



void PrintVariables (RooArgSet* setVar, string type)
{
  RooRealVar* tmpVar;
  int nEleSet = setVar->getSize();


  if (type == "vars")
    {
      cout << "\n@@@@@@@@@@@@@@@@@@@@@@@@@@" << endl;
      cout <<   "@@@ Printing variables @@@" << endl;
      cout <<   "@@@@@@@@@@@@@@@@@@@@@@@@@@" << endl;
      
      if (setVar != NULL)
	{
	  TIterator* it = setVar->createIterator();
	  for (int i = 0; i < nEleSet; i++)
	    {
	      tmpVar = (RooRealVar*)it->Next();
	      cout << "Variable: " << i;
	      cout << "\tname: "   << tmpVar->GetName();
	      cout << "\tvalue: "  << tmpVar->getVal();
	      cout << "\terr: "    << tmpVar->getError();
	      cout << "\tErrLo: "  << tmpVar->getErrorLo();
	      cout << "\tErrHi: "  << tmpVar->getErrorHi() << endl;
	    }
	}
    }
  else if (type == "cons")
    {
      cout << "\n@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@" << endl;
      cout <<   "@@@@@@@@@ Printing constraints @@@@@@@@@" << endl;
      cout <<   "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@" << endl;

      if (setVar != NULL)
	{
	  TIterator* it = setVar->createIterator();
	  for (int i = 0; i < nEleSet; i++) PrintVariables(((RooAbsPdf*)it->Next())->getVariables(),"vars");
	}
    }
  else
    {
      cout << "[ExtractYield::PrintVariables]\tWrong parameter: " << type << endl;
      exit (EXIT_FAILURE);
    }
}


void ClearVars (RooArgSet* vecVars)
{
  if (vecVars != NULL)
    {
      int nEle = vecVars->getSize();
      
      TIterator* it = vecVars->createIterator();
      for (int i = 0; i < nEle; i++) delete it->Next();
      
      vecVars->removeAll();
    }
}


void SetStyle ()
{
  gROOT->SetStyle("Plain");
  gROOT->ForceStyle();
  gStyle->SetTextFont(42);

  gStyle->SetOptFit(1112);
  gStyle->SetOptStat(1110);
  gStyle->SetOptTitle(0);

  gStyle->SetPadRightMargin(0.02);
  gStyle->SetPadTopMargin(0.11);
  gStyle->SetPadBottomMargin(0.12);

  gStyle->SetTitleFont(42,"x");
  gStyle->SetTitleFont(42,"y");
  gStyle->SetTitleOffset(1.05,"x");
  gStyle->SetTitleOffset(0.95,"y");
  gStyle->SetTitleSize(0.05,"x");
  gStyle->SetTitleSize(0.05,"y");

  gStyle->SetLabelFont(42,"x");
  gStyle->SetLabelFont(42,"y");
  gStyle->SetLabelSize(0.05,"x");
  gStyle->SetLabelSize(0.05,"y");

  TGaxis::SetMaxDigits(3);
  gStyle->SetStatY(0.9);
}


string MakeName (const RooDataSet* data,  int ID)
{
  stringstream myString;


  myString.clear(); myString.str("");
  myString << data->GetName() << "_" << ID;


  return myString.str();
}


void DrawString (double Lumi, RooPlot* myFrame)
{
  stringstream myString;
  double scaleRespect2CMS = 0.75;


  myString.clear(); myString.str("");
  myString << "CMS";
  TLatex* LumiTex1 = new TLatex(0.18,0.9,myString.str().c_str());
  LumiTex1->SetTextFont(61);
  LumiTex1->SetTextSize(0.05);
  LumiTex1->SetTextColor(kBlack);
  LumiTex1->SetNDC(true);
  if (myFrame == NULL) LumiTex1->DrawLatex(0.18,0.9,myString.str().c_str());
  else
    {
      LumiTex1->Paint();
      myFrame->addObject(LumiTex1);
    }


  myString.clear(); myString.str("");
  myString << "#it{Preliminary}";
  TLatex* LumiTex2 = new TLatex(0.265,0.9,myString.str().c_str());
  LumiTex2->SetTextFont(42);
  LumiTex2->SetTextSize(0.05 * scaleRespect2CMS);
  LumiTex2->SetTextColor(kBlack);
  LumiTex2->SetNDC(true);
  if (myFrame == NULL) LumiTex2->DrawLatex(0.265,0.9,myString.str().c_str());
  else
    {
      LumiTex2->Paint();
      myFrame->addObject(LumiTex2);
    }


  myString.clear(); myString.str("");
  myString << Lumi <<  " fb#lower[0.4]{^{#font[122]{\55}1}} (8 TeV)";
  TLatex* LumiTex3 = new TLatex(0.8,0.9,myString.str().c_str());
  LumiTex3->SetTextFont(42);
  LumiTex3->SetTextSize(0.05 * scaleRespect2CMS);
  LumiTex3->SetTextColor(kBlack);
  LumiTex3->SetNDC(true);
  if (myFrame == NULL) LumiTex3->DrawLatex(0.8,0.9,myString.str().c_str());
  else
    {
      LumiTex3->Paint();
      myFrame->addObject(LumiTex3);
    }
}

void AddGaussConstraint (RooArgSet* vecConstr, RooAbsPdf* pdf, string varName)
{
  stringstream myString;


  myString.clear(); myString.str("");
  myString << varName;

  if (GetVar(pdf,myString.str().c_str())->isConstant() == false)
    {
      RooRealVar* varConstr = GetVar(pdf,myString.str().c_str());
      double mean  = GetVar(pdf,myString.str().c_str())->getVal();
      double sigma = (varConstr->getErrorHi() - varConstr->getErrorLo()) / 2.;
      
      myString << "_constr";
      
      RooGaussian* newConstr = new RooGaussian(myString.str().c_str(), myString.str().c_str(), *varConstr, RooConst(mean), RooConst(sigma));
      vecConstr->add(*newConstr);
    }
}

double PrintFitResults (RooAbsPdf* pdf, RooFitResult* fitResult)
{
  double varVal, varValELo, varValEHi;

  cout << "[ExtractYield::PrintFitResults]\t Print Fitting Results" << endl;
  if (fitResult != NULL)
     {
      if (GetVar(pdf,"FlS") != NULL)
        {
          cout << " # FlS +/- err" << endl;
          cout << varVal << " +/- " << (varValEHi - varValELo) / 2. << " (" << varValEHi << "/" << varValELo << ")" << endl; 
        }      
      if (GetVar(pdf,"P1") != NULL)
        {
          cout << " # P1 +/- err" << endl;
          cout << GetVar(pdf,"P1")->getVal() << " +/- " << GetVar(pdf,"P1")->getError() <<  " (" << GetVar(pdf,"P1")->getErrorHi() << "/" << GetVar(pdf,"P1")->getErrorLo() << ")" << endl;
        }
      if (GetVar(pdf,"P2") != NULL)
        {
          cout << " # P2 +/- err" << endl;
          cout << GetVar(pdf,"P2")->getVal() << " +/- " << GetVar(pdf,"P2")->getError() <<  " (" << GetVar(pdf,"P2")->getErrorHi() << "/" << GetVar(pdf,"P2")->getErrorLo() << ")" << endl;
        }
      if (GetVar(pdf,"P3") != NULL)
        {
          cout << " # P3 +/- err" << endl;
          cout << GetVar(pdf,"P3")->getVal() << " +/- " << GetVar(pdf,"P3")->getError() <<  " (" << GetVar(pdf,"P3")->getErrorHi() << "/" << GetVar(pdf,"P3")->getErrorLo() << ")" << endl;
        }
      if (GetVar(pdf,"P4p") != NULL)
        {
          cout << " # P4p +/- err" << endl;
          cout << GetVar(pdf,"P4p")->getVal() << " +/- " << GetVar(pdf,"P4p")->getError() <<  " (" << GetVar(pdf,"P4p")->getErrorHi() << "/" << GetVar(pdf,"P4p")->getErrorLo() << ")" << endl;
        }
      if (GetVar(pdf,"P5p") != NULL)
        {
          cout << " # P5p +/- err" << endl;
          cout << GetVar(pdf,"P5p")->getVal() << " +/- " << GetVar(pdf,"P5p")->getError() <<  " (" << GetVar(pdf,"P5p")->getErrorHi() << "/" << GetVar(pdf,"P5p")->getErrorLo() << ")" << endl;
        }
      if (GetVar(pdf,"P6p") != NULL)
        {
          cout << " # P6p +/- err" << endl;
          cout << GetVar(pdf,"P6p")->getVal() << " +/- " << GetVar(pdf,"P6p")->getError() <<  " (" << GetVar(pdf,"P6p")->getErrorHi() << "/" << GetVar(pdf,"P6p")->getErrorLo() << ")" << endl;
        }
      if (GetVar(pdf,"P8p") != NULL)
        {
          cout << " # P8p +/- err" << endl;
          cout << GetVar(pdf,"P8p")->getVal() << " +/- " << GetVar(pdf,"P8p")->getError() <<  " (" << GetVar(pdf,"P8p")->getErrorHi() << "/" << GetVar(pdf,"P8p")->getErrorLo() << ")" << endl;
        } 
      if (GetVar(pdf,"fracMisTag") != NULL)
        {
          cout << " # fracMisTag +/- err" << endl;
          cout << GetVar(pdf,"fracMisTag")->getVal() << " +/- " << GetVar(pdf,"fracMisTag")->getError() <<  " (" << GetVar(pdf,"fracMisTag")->getErrorHi() << "/" << GetVar(pdf,"fracMisTag")->getErrorLo() << ")" << endl;
        }
     }
  return 1;
}


void MakeDatasets (B0KstMuMuTreeContent* NTuple, unsigned int FitType, int SampleType)
{
  stringstream myString;

  // ###########################
  // # Define useful variables #
  // ###########################
   B0mass  = new RooRealVar("B0mass","#font[12]{m}(#font[122]{K}#kern[0.1]{#lower[0.4]{^{#font[122]{+}}}}#kern[-0.3]{#pi}#kern[-0.3]{#lower[0.6]{^{#font[122]{\55}}}}#mu#kern[-0.9]{#lower[0.6]{^{#font[122]{+}}}}#kern[-0.1]{#mu}#kern[-1.3]{#lower[0.6]{^{#font[122]{\55}}}})",Utility->B0Mass - atof(Utility->GetGenericParam("B0MassIntervalLeft").c_str()),Utility->B0Mass + atof(Utility->GetGenericParam("B0MassIntervalRight").c_str()),"GeV");
  mumuMass           = new RooRealVar("mumuMass","#mu#kern[-0.9]{#lower[0.6]{^{#font[122]{+}}}}#kern[-0.1]{#mu}#kern[-1.3]{#lower[0.6]{^{#font[122]{\55}}}} inv. mass",0.0,6.0,"GeV");
  mumuMassE          = new RooRealVar("mumuMassE","#mu#kern[-0.9]{#lower[0.6]{^{#font[122]{+}}}}#kern[-0.1]{#mu}#kern[-1.3]{#lower[0.6]{^{#font[122]{\55}}}} inv. mass error",0.0,0.5,"GeV");
  ctK     = new RooRealVar("ctK","cos#theta_{K}",-1.0,1.0);
  ctL     = new RooRealVar("ctL","cos#theta_{L}",-1.0,1.0);
  phi     = new RooRealVar("phi","#phi", - Utility->PI ,Utility->PI);
  tagB0   = new RooRealVar("tagB0","Tagged Events",-0.1,1.1);
  genSignal     = new RooRealVar("genSignal","gen Signal",0.0,3.0);

  if ( (FitType == 205) || (FitType == 206) ||
       (FitType == 105) || (FitType == 106)
     )
    {
      RooArgSet Vars;
      Vars.add(*B0mass);
      Vars.add(*mumuMass);
      Vars.add(*mumuMassE);
      Vars.add(*ctK);
      Vars.add(*ctL);
      Vars.add(*phi);
      Vars.add(*tagB0);
      Vars.add(*genSignal);

      SingleCandNTuple           = new RooDataSet("SingleCandNTuple"          ,"SingleCandNTuple"          ,Vars);
      SingleCandNTuple_RejectPsi = new RooDataSet("SingleCandNTuple_RejectPsi","SingleCandNTuple_RejectPsi",Vars);

      // #############################
      // # Load values from the tree #
      // #############################
      NTuple->ClearNTuple( SampleType );
      NTuple->SetBranchAddresses(theTree, SampleType);
      int nEntries = theTree->GetEntries();
      for (int entry = 0;  entry <nEntries; entry++)
	{
	  theTree->GetEntry(entry);
	 
	  if (
              (FitType == 105) || (FitType == 106) || (FitType == 206) || (FitType == 205)
	     )
	    { 
              if ((FitType == 206))
                {
                   Vars.setRealValue("mumuMass",          NTuple->genQ);
                   Vars.setRealValue("ctK",               NTuple->gen_cos_theta_k);
                   Vars.setRealValue("ctL",               NTuple->gen_cos_theta_l);
                   Vars.setRealValue("phi",               NTuple->gen_phi_kst_mumu);
                }
              else if (((FitType == 105) || (FitType == 106)) 
                      )
                {
                   Vars.setRealValue("B0mass",            NTuple->tagged_mass);
        	   Vars.setRealValue("mumuMass",          NTuple->mumuMass);
                   Vars.setRealValue("mumuMassE",         NTuple->mumuMassE);
	           Vars.setRealValue("ctK",               NTuple->cos_theta_k);
                   Vars.setRealValue("ctL",               NTuple->cos_theta_l);
                   Vars.setRealValue("phi",               NTuple->phi_kst_mumu);
                   Vars.setRealValue("tagB0",             NTuple->tagB0);
                   Vars.setRealValue("genSignal",         NTuple->genSignal);                         
                }

              // ########################
              // # NTuple with all data #
              // ########################
	      SingleCandNTuple->add(Vars);
          if (
               (((FitType == 105) || (FitType == 106)) &&
                (((strcmp(CTRLfitWRKflow.c_str(),"trueGoodtag") == 0) && ((NTuple->tagB0 ==1 && NTuple->genSignal ==1) || (NTuple->tagB0 == 0 && NTuple->genSignal ==2))) ||
                 ((strcmp(CTRLfitWRKflow.c_str(),"trueMistag")  == 0) && (!((NTuple->tagB0 ==1 && NTuple->genSignal ==1) || (NTuple->tagB0 ==0 && NTuple->genSignal ==2))))||
                  (strcmp(CTRLfitWRKflow.c_str(),"trueAll")     == 0) ||
                  (strcmp(CTRLfitWRKflow.c_str(),"allEvts")     == 0)) &&
                (Utility->PsiRejection(NTuple->tagged_mass,NTuple->mumuMass,NTuple->mumuMassE,"rejectPsi",true) == true))
             )
           SingleCandNTuple_RejectPsi->add(Vars);             
         }	
       }
      cout << "\n[ExtractYield::MakeDatasets]\t@@@ NTuple with all data @@@" << endl;
      SingleCandNTuple->Print("v");

      cout << "[ExtractYield::MakeDatasets]\t@@@ NTuple without J/psi and psi(2S) regions @@@" << endl;
      SingleCandNTuple_RejectPsi->Print("v");
      }
  // ####################################################
  // # Setting initial values for independent variables #
  // ####################################################
   B0mass->setVal(Utility->B0Mass);
   ctK ->setVal(0.0);
   ctL->setVal(0.0);
   phi->setVal(0.0);
 }



     //#############################
     //# GEN level fitting  PDF #
     //#############################
void Instantiate3AnglesFit    (RooAbsPdf** TotalPDF,
				RooRealVar* y, RooRealVar* z,RooRealVar* p,
				unsigned int FitType,
			        unsigned int q2BinIndx)
// #########################
// # y: angle cos(theta_l) #
// # z: angle cos(theta_K) #
// # p: angle phi          #
// #########################
{

  // ###################
  // # Local variables #
  // ###################
    stringstream myString;
    stringstream effString;

  RooRealVar* Fl = new RooRealVar("Fl","F_{L}",0.5,0,1);
  RooRealVar* P1 = new RooRealVar("P1","P_{1}",0,-1,1);
  RooRealVar* P2 = new RooRealVar("P2","P_{2}",0,-1,1);
  RooRealVar* P3 = new RooRealVar("P3","P_{3}",0,-1,1);
  RooRealVar* P4p = new RooRealVar("P4p","P'_{4}",0,-1*sqrt(2),sqrt(2));
  RooRealVar* P5p = new RooRealVar("P5p","P'_{5}",0,-1*sqrt(2),sqrt(2));
  RooRealVar* P6p = new RooRealVar("P6p","P'_{6}",0,-1*sqrt(2),sqrt(2));
  RooRealVar* P8p = new RooRealVar("P8p","P'_{8}",0,-1*sqrt(2),sqrt(2));
  RooRealVar*  fracMisTag = new RooRealVar("fracMisTag","Fraction of mistag",1);
  fracMisTag->setConstant();

  if (FitType == 206)
      *TotalPDF = new DecayRate("TotalPDF","TotalPDF",*z,*y,*p,*Fl,*P1,*P2,*P3,*P4p,*P5p,*P6p,*P8p);
  else if (FitType == 106)
    {
        // ################################
        //  # Read configuration variables #
        //  ################################
       // unsigned int useSignal = configParam->operator[](Utility->GetConfigParamIndx("SigType"))->operator[](q2BinIndx);
        int parity = 1;
        RooArgList vars (*z, *y, *p);
        //import KDE efficiency histograms
        myString.clear(); myString.str("");
        myString <<"/afs/cern.ch/user/l/llinwei/work2/B0KstMuMufull/B0KstMuMu2/efficiency/KDE/latestVersion-integrals/KDEeff_b" <<  q2BinIndx << "_od.root";
        std::cout << myString.str().c_str() << std::endl;
        TFile* efffile= new TFile(myString.str().c_str(),"READ");

        // two tag
        TH3D* effCHist = (TH3D*)efffile->Get(Form("effCHist_b%ip%i",q2BinIndx,parity));
        TH3D* effWHist = (TH3D*)efffile->Get(Form("effWHist_b%ip%i",q2BinIndx,parity));
        RooDataHist* effCData = new RooDataHist("effCData","effCData",vars,effCHist);
        RooDataHist* effWData = new RooDataHist("effWData","effWData",vars,effWHist);
        RooAbsReal* effC = new RooHistFunc("effC","effC",vars,*effCData,1);
        RooAbsReal* effW = new RooHistFunc("effW","effW",vars,*effWData,1);

        // import precomputed integrals
        vector<double> intCVec (0);
        vector<double> intWVec (0);
        TH1D* intCHist = (TH1D*)efffile->Get(Form("MCint_b%ip%it1",q2BinIndx,parity));
        TH1D* intWHist = (TH1D*)efffile->Get(Form("MCint_b%ip%it0",q2BinIndx,parity));

        if ( !intCHist || intCHist->IsZombie() || !intWHist || intWHist->IsZombie() ) {
            cout<<"Integral histograms not found in file: "<<endl<<"Using rooFit integration"<<endl;
            intCVec.push_back(0);
            intWVec.push_back(0);
        } else if ( strcmp( intCHist->GetTitle(), effCHist->GetTitle() ) || strcmp( intWHist->GetTitle(), effWHist->GetTitle() ) ) {
            cout<<"Integral histograms are incoherent with efficiency in file: "<<endl;
            cout<<"Efficiency (CT) conf: "<<effCHist->GetTitle()<<endl;
            cout<<"Integral (CT) conf: "<<intCHist->GetTitle()<<endl;
            cout<<"Efficiency (WT) conf: "<<effWHist->GetTitle()<<endl;
            cout<<"Integral (WT) conf: "<<intWHist->GetTitle()<<endl;
            cout<<"Using rooFit integration"<<endl;
            intCVec.push_back(0);
            intWVec.push_back(0);
        } else {
            for (int i=1; i<=intCHist->GetNbinsX(); ++i) intCVec.push_back(intCHist->GetBinContent(i));
            for (int i=1; i<=intWHist->GetNbinsX(); ++i) intWVec.push_back(intWHist->GetBinContent(i));
        }


        if (strcmp(CTRLfitWRKflow.c_str(),"trueGoodtag") == 0)    *TotalPDF = new PdfRT("TotalPDF","TotalPDF",*z,*y,*p,*Fl,*P1,*P2,*P3,*P4p,*P5p,*P6p,*P8p,*effC,intCVec);
        else if (strcmp(CTRLfitWRKflow.c_str(),"trueMistag")  == 0)     *TotalPDF = new PdfWT("TotalPDF","TotalPDF",*z,*y,*p,*Fl,*P1,*P2,*P3,*P4p,*P5p,*P6p,*P8p,*effW,intWVec);
        else if ((strcmp(CTRLfitWRKflow.c_str(),"trueAll") == 0) || (strcmp(CTRLfitWRKflow.c_str(),"allEvts") == 0) )   *TotalPDF = new PdfSigAng("TotalPDF","TotalPDF",*z,*y,*p,*Fl,*P1,*P2,*P3,*P4p,*P5p,*P6p,*P8p,*fracMisTag,*effC,*effW,intCVec,intWVec);
    }
   else
    {
      cout << "[ExtractYield::Instantiate3AnglesFit]\tIncorrect configuration sequence : useSignal = " << "useSignal" <<  endl;
      exit (EXIT_FAILURE);
    }

}

RooFitResult* Make3AnglesFit (RooDataSet* dataSet, RooAbsPdf** TotalPDF, RooRealVar* y, RooRealVar* z,RooRealVar* p,  unsigned int FitType, RooArgSet* vecConstr, TCanvas* Canv, int specBin)
{
  // ###################
  // # Local variables #
  // ###################
  RooFitResult* fitResult = NULL;
  stringstream myString;

  unsigned int nElements   = 0;
  unsigned int it          = 0;
  TString legNames[6]; 
  TLegend*   legY             = NULL;
  TLegend*   legZ             = NULL;
  TLegend*   legP             = NULL;

  if ((FitType == 206) || (FitType == 205) || (FitType ==105) || (FitType == 106))
    {
      clock_t fitstart = clock();
      TBenchmark* Fittime = new TBenchmark();
      Fittime->Start("test");
      // ###################
      // # Make actual fit #
      // ###################
      if (ROOMINIMIZER == true )
       {
         clock_t nlls = clock();
         TBenchmark* nll = new TBenchmark();
         nll->Start("nll");
         RooAbsReal* MC_nll = (*TotalPDF)->createNLL(*dataSet,RooFit::NumCPU(10));    
         RooMinuit Minuit(*MC_nll) ;
         Minuit.setPrintLevel(-1);
         clock_t nlle = clock();
         double nlltime = (double)(nlle - nlls)/ CLOCKS_PER_SEC;
         cout <<  "[ExtractYield::Make3AnglesFit]\t @@@ Time for nll clock_t \t" << nlltime << "s\t @@@" << endl;
         nll->Stop("nll");
         nll->Show("nll");

         clock_t migrads = clock();
         TBenchmark* migrad = new TBenchmark();
         migrad->Start("migrad");
         Minuit.migrad() ;
         clock_t migrade = clock();
         double migradtime = (double)(migrade - migrads)/ CLOCKS_PER_SEC;
         cout <<  "[ExtractYield::Make3AnglesFit]\t @@@ Time for Migrad clock_t \t" << migradtime << "s\t @@@" << endl;
         migrad->Stop("migrad");
         migrad->Show("migrad");
        
         if (Minuit.migrad() == 0) {
         clock_t hesses = clock();
         TBenchmark* hesse = new TBenchmark();
         hesse->Start("hesse");
         Minuit.hesse();
         clock_t hessee = clock();
         double hessetime = (double)(hessee - hesses)/ CLOCKS_PER_SEC;
         cout <<  "[ExtractYield::Make3AnglesFit]\t @@@ Time for Hesse clock_t \t" << hessetime << "s\t @@@" << endl;
         hesse->Stop("hesse");
         hesse->Show("hesse");
         fitResult = Minuit.save(); 
         }
        else exit(EXIT_FAILURE);
       }
      else 
	fitResult = (*TotalPDF)->fitTo(*dataSet,Save(true),Minos(0),Timer(true),Minimizer(MINIMIZER),NumCPU(10));
      cout <<  "[ExtractYield::Make3AnglesFit]\t @@@ status = " << fitResult->status() << endl;
      clock_t fitend = clock();
      double fittime = (double) (fitend - fitstart)/ CLOCKS_PER_SEC;
       cout <<  "[ExtractYield::Make3AnglesFit]\t @@@ Time for fitting \t" << fittime << "s\t @@@" << endl; 
       
       Fittime->Stop("test");
       Fittime->Show("test");
      // ###################################################
      // # Set p.d.f. independent variables to known point #
      // ###################################################
      if (GetVar(*TotalPDF,y->getPlotLabel()) != NULL) (*TotalPDF)->getVariables()->setRealValue(y->getPlotLabel(),0.0);
      if (GetVar(*TotalPDF,z->getPlotLabel()) != NULL) (*TotalPDF)->getVariables()->setRealValue(z->getPlotLabel(),0.0);
      if (GetVar(*TotalPDF,p->getPlotLabel()) != NULL) (*TotalPDF)->getVariables()->setRealValue(p->getPlotLabel(),0.0);
      if (fitResult != NULL) fitResult->Print("v");
      
      // ###########################
      // # Costheta-l plot results #
      // ###########################
       clock_t plotstart = clock();
      if ((PLOT==true ) &&  CheckGoodFit(fitResult) == true)
      {
      Canv->cd(1);
      RooPlot* myFrameY = y->frame(NBINS);

      dataSet->plotOn(myFrameY, Name(MakeName(dataSet,specBin).c_str()));
      if ((FitType == 206) || (FitType == 205))  legNames[nElements++] = "GEN-MC";
      else if ((FitType == 106) || (FitType == 105))  legNames[nElements++] = "Reco-MC";
      (*TotalPDF)->plotOn(myFrameY, Name((*TotalPDF)->getPlotLabel()), LineColor(kBlack), Project(RooArgSet(*z,*p)));
      if ((FitType == 206) || (FitType == 205))  legNames[nElements++] = "Gen p.d.f.";
      else if ((FitType == 106) || (FitType == 105))  legNames[nElements++] = "Reco p.d.f.";
      TPaveText* paveTextY = new TPaveText(0.12,0.78,0.4,0.86,"NDC");
      paveTextY->AddText(Form("%s%.2f","#chi#lower[0.4]{^{2}}/DoF = ",myFrameY->chiSquare((*TotalPDF)->getPlotLabel(),MakeName(dataSet,specBin).c_str())));
      CheckGoodFit(fitResult,paveTextY);
      paveTextY->SetTextAlign(11);
      paveTextY->SetBorderSize(0.0);
      paveTextY->SetFillStyle(0);
      paveTextY->SetTextSize(0.04);
      paveTextY->Paint();
      myFrameY->addObject(paveTextY);

      legY = new TLegend(0.78, 0.65, 0.97, 0.88, "");
      it = 0;
      for (unsigned int i = 0; i < 2*nElements; i++)
	{
	  TString objName = myFrameY->nameOf(i);
	  if ((objName == "") || (objName.Contains("paramBox") == true) || (objName.Contains("TPave") == true) || ((i > 0) && (objName == myFrameY->nameOf(i-1)))) continue;
	  TObject* obj = myFrameY->findObject(objName.Data());
	  legY->AddEntry(obj,legNames[it++],"PL");
	  legY->SetTextFont(42);
	}     
      legY->SetFillStyle(0);
      legY->SetFillColor(0);
      legY->SetTextSize(0.04);
      legY->SetBorderSize(0);
     
      myFrameY->SetMaximum(1.4*myFrameY->GetMaximum());
      myFrameY->Draw();
      legY->Draw("same");
     
      // ###########################
      // # Costheta-k plot results #
      // ###########################
      Canv->cd(2);
      RooPlot* myFrameZ = z->frame(NBINS);

      dataSet->plotOn(myFrameZ, Name(MakeName(dataSet,specBin).c_str()));
      if ((FitType == 206) || (FitType == 205))  legNames[nElements++] = "GEN-MC";
      else if ((FitType == 106) || (FitType == 105))  legNames[nElements++] = "Reco-MC";
      (*TotalPDF)->plotOn(myFrameZ, Name((*TotalPDF)->getPlotLabel()), LineColor(kBlack), Project(RooArgSet(*y,*p)));
      if ((FitType == 206) || (FitType == 205))  legNames[nElements++] = "Gen p.d.f.";
      else if ((FitType == 106) || (FitType == 105))  legNames[nElements++] = "Reco p.d.f";
 
      TPaveText* paveTextZ = new TPaveText(0.12,0.82,0.4,0.86,"NDC");
      paveTextZ->AddText(Form("%s%.2f","#chi#lower[0.4]{^{2}}/DoF = ",myFrameZ->chiSquare((*TotalPDF)->getPlotLabel(),MakeName(dataSet,specBin).c_str())));
      paveTextZ->SetTextAlign(11);
      paveTextZ->SetBorderSize(0.0);
      paveTextZ->SetFillStyle(0);
      paveTextZ->SetTextSize(0.04);
      paveTextZ->Paint();
      myFrameZ->addObject(paveTextZ);

      legZ = new TLegend(0.78, 0.65, 0.97, 0.88, "");
      it = 0;
      for (unsigned int i = 0; i < 2*nElements; i++)
	{
	  TString objName = myFrameZ->nameOf(i);
	  if ((objName == "") || (objName.Contains("paramBox") == true) || (objName.Contains("TPave") == true) || ((i > 0) && (objName == myFrameZ->nameOf(i-1)))) continue;
	  TObject* obj = myFrameZ->findObject(objName.Data());
	  legZ->AddEntry(obj,legNames[it++],"PL");
	  legZ->SetTextFont(42);
	}     
      legZ->SetFillStyle(0);
      legZ->SetFillColor(0);
      legZ->SetTextSize(0.04);
      legZ->SetBorderSize(0);
     
      myFrameZ->SetMaximum(1.4*myFrameZ->GetMaximum());
      myFrameZ->Draw();
      legZ->Draw("same");


      // ###########################
      // # Costheta-l plot results #
      // ###########################
      Canv->cd(3);
      RooPlot* myFrameP = p->frame(NBINS);

      dataSet->plotOn(myFrameP, Name(MakeName(dataSet,specBin).c_str()));
      if ((FitType == 206) || (FitType == 205))  legNames[nElements++] = "GEN-MC";
      else if ((FitType == 106) || (FitType == 105))  legNames[nElements++] = "Reco-MC";
      (*TotalPDF)->plotOn(myFrameP, Name((*TotalPDF)->getPlotLabel()), LineColor(kBlack), Project(RooArgSet(*y,*z)));
      if ((FitType == 206) || (FitType == 205))  legNames[nElements++] = "Gen p.d.f.";
      else if ((FitType == 106) || (FitType == 105))  legNames[nElements++] = "Reco p.d.f";      

      TPaveText* paveTextP = new TPaveText(0.12,0.82,0.4,0.86,"NDC");
      paveTextP->AddText(Form("%s%.2f","#chi#lower[0.4]{^{2}}/DoF = ",myFrameP->chiSquare((*TotalPDF)->getPlotLabel(),MakeName(dataSet,specBin).c_str())));
      paveTextP->SetTextAlign(11);
      paveTextP->SetBorderSize(0.0);
      paveTextP->SetFillStyle(0);
      paveTextP->SetTextSize(0.04);
      paveTextP->Paint();
      myFrameP->addObject(paveTextP);

      legP = new TLegend(0.78, 0.65, 0.97, 0.88, "");
      it = 0;
      for (unsigned int i = 0; i < 2*nElements; i++)
	{
	  TString objName = myFrameP->nameOf(i);
	  if ((objName == "") || (objName.Contains("paramBox") == true) || (objName.Contains("TPave") == true) || ((i > 0) && (objName == myFrameY->nameOf(i-1)))) continue;
	  TObject* obj = myFrameP->findObject(objName.Data());
	  legP->AddEntry(obj,legNames[it++],"PL");
	  legP->SetTextFont(42);
	}     
      legP->SetFillStyle(0);
      legP->SetFillColor(0);
      legP->SetTextSize(0.04);
      legP->SetBorderSize(0);
     
      myFrameP->SetMaximum(1.4*myFrameP->GetMaximum());
      myFrameP->Draw();
      legP->Draw("same");
    }

    Canv->Modified();
    Canv->Update();

    if (SAVEPLOT == true)
    {
     myString.clear(); myString.str("");
     myString << (*TotalPDF)->getPlotLabel() << "_Canv" << specBin << ".pdf";
     Canv->Print(myString.str().c_str());

     myString.clear(); myString.str("");
     myString << (*TotalPDF)->getPlotLabel() << "_Canv" << specBin  << ".root";
     Canv->Print(myString.str().c_str());
    }
   clock_t plotend = clock(); 
   double plottime = (double) (plotend - plotstart)/ CLOCKS_PER_SEC;
   cout <<  "[ExtractYield::Make3AnglesFit]\t @@@ Time for plotting \t" << plottime << "s\t @@@" << endl;
  } 
    return fitResult;
}


void Iterative3AnglesFitq2Bins (RooDataSet* dataSet,
				    RooRealVar* y, RooRealVar* z,RooRealVar* p,
				    int specBin,
				    unsigned int FitType,
				    vector<double>* q2Bins,
				    RooArgSet* vecConstr)
{
  // ###################
  // # Local variables #
  // ###################
  stringstream myString;   

  if ( (FitType == 105) || (FitType == 106) || (FitType == 205) || (FitType == 206) )
    {
      RooAbsPdf*  TotalPDFq2Bins[q2Bins->size()-1];
      RooFitResult* fitResult = NULL;

      RooDataSet* dataSet_q2Bins[q2Bins->size()-1];
      for (unsigned int i = (specBin == -1 ? 0 : specBin); i < (specBin == -1 ? q2Bins->size()-1 : specBin+1); i++)
       {
         myString.clear(); myString.str("");
         myString << "(mumuMass*mumuMass) > " << q2Bins->operator[](i) << " && (mumuMass*mumuMass) <= " << q2Bins->operator[](i+1);
         cout << "\n[ExtractYield::IterativeAnglesFitq2Bins]\tCut string: " << myString.str() << endl;
         dataSet_q2Bins[i] = (RooDataSet*)dataSet->reduce(myString.str().c_str());
         cout << "[ExtractYield::IterativeAnglesFitq2Bins]\tNumber of events : " << dataSet_q2Bins[i]->sumEntries() << endl;

      if ( FitType == 206)
        {
          TCanvas*    Gen[q2Bins->size()-1];
          Gen[i] = new TCanvas("Gen","Gen",10,10,700,500);
          Gen[i]->Divide(2,2);

          myString.clear(); myString.str("");
          myString << "GenTotalPDFq2Bin_" << i;
          Instantiate3AnglesFit(&TotalPDFq2Bins[i],y,z,p,FitType,i);

          // ########################################
          // # Initialize p.d.f & Apply constraints #
          // ########################################
          clock_t copystart = clock();
          ClearVars(vecConstr);
 
          fitResult = Make3AnglesFit(dataSet_q2Bins[i],&TotalPDFq2Bins[i],y,z,p,FitType,vecConstr,Gen[i],specBin);
          if (CheckGoodFit(fitResult) == true) cout << "\n[ExtractYield::Iterative3AnglesFitq2Bins]\t@@@ Fit converged ! @@@" << endl;
          else                                 cout << "\n[ExtractYield::Iterative3AnglesFitq2Bins]\t@@@ Fit didn't converge ! @@@" << endl;
                
        }
      else if (  FitType == 105 || FitType == 106)
        {
          TCanvas*    Reco[q2Bins->size()-1]; 
          Reco[i] = new TCanvas("Reco","Reco",10,10,700,500);
          Reco[i]->Divide(2,2);

          unsigned int countMisTag  = 0;
          unsigned int countGoodTag = 0;
          for (int j = 0; j < static_cast<int>(dataSet_q2Bins[i]->sumEntries()); j++)
           {
             if ( ((dataSet_q2Bins[i]->get(j)->getRealValue("tagB0") == 1) && (dataSet_q2Bins[i]->get(j)->getRealValue("genSignal") == 1)) || ((dataSet_q2Bins[i]->get(j)->getRealValue("tagB0") == 0) && (dataSet_q2Bins[i]->get(j)->getRealValue("genSignal") == 2)))   
             countGoodTag++;
             else
              countMisTag++;
           } 
          cout << "[ExtractYield::IterativeMassFitq2Bins]\tDynamic mis-tag fraction : " << static_cast<double>(countMisTag) / static_cast<double>(countMisTag + countGoodTag) << " = (" << countMisTag << "/(" << countMisTag << "+" << countGoodTag << "))" << endl;
          cout << countGoodTag << endl;
          cout << countMisTag << endl; 
          myString.clear(); myString.str("");
          myString << "RecoTotalPDFq2Bin_" << i;
          Instantiate3AnglesFit(&TotalPDFq2Bins[i],y,z,p,FitType,i);
          ClearVars(vecConstr);

      // ###################
      // # Perform the fit #
      // ##################
      fitResult = Make3AnglesFit(dataSet_q2Bins[i],&TotalPDFq2Bins[i],y,z,p,FitType,vecConstr,Reco[i],specBin);
      if (CheckGoodFit(fitResult) == true) cout << "\n[ExtractYield::Iterative3AnglesFitq2Bins]\t@@@ Fit converged ! @@@" << endl;
      else                                 cout << "\n[ExtractYield::Iterative3AnglesFitq2Bins]\t@@@ Fit didn't converge ! @@@" << endl;
        }
      PrintFitResults(*&TotalPDFq2Bins[i], fitResult);
    }
  }
}



int main(int argc, char** argv)
{
  
  if (argc >= 4)
    {
      // ##################
      // # Main variables #
      // ##################
      stringstream myString;
      
      string fileName           = "";
      string correct4Efficiency = "";
      string tmpFileName        = "";

      int specBin               = -1;
      unsigned int FitType      = atoi(argv[1]);
      int SampleType            = -1;
      TFile* NtplFile           = NULL;

      if (((FitType == 106) ||
           (FitType == 206))
          && (argc >= 3))
	  
	{
	  ParameterFILE = Utility->MakeAnalysisPATH(PARAMETERFILEIN).c_str();


 	  // ###################
	  // # Read parameters #
 	  // ###################
	  Utility = new Utils(false);
	  Utility->ReadAllBins(ParameterFILE,&q2Bins);//,&cosThetaKBins,&cosThetaLBins,&phiBins);

	  // #################################
	  // # Check that FitType is correct #
	  // #################################
	  if ((FitType == 106) ||
	      (FitType == 206))
	    {
	      fileName           = argv[2];
	    }
          if ( (FitType == 106))
          SampleType         = 1;
          else if ((FitType == 206) )
          SampleType         = 0;

  
          // #################################
          // # Read the q^2 bin and the rest #
          // #################################
          if (argc >= 4) specBin = atoi(argv[3]);


          cout << "\n[ExtractYield::main]\t@@@ Input variables from command line @@@" << endl;
          cout << "- input/outputFile.root = " << fileName.c_str() << endl;
          cout << "- specBin = "               << specBin << endl;
          cout << "- FitType = "               << FitType << endl;
          cout << "- SampleType = "            << SampleType << endl;
          cout << "- ParameterFILE = "         << ParameterFILE.c_str() << endl;

          cout << "\n[ExtractYield::main]\t@@@ Internal settings @@@" << endl;
          cout << "NBINS = "         << NBINS << endl;
          cout << "MULTYIELD = "     << MULTYIELD << endl;
          cout << "NCOEFFPOLYBKG = " << NCOEFFPOLYBKG << endl;

          cout << "\nPARAMETERFILEIN = " << PARAMETERFILEIN << endl;
          cout << "PARAMETERFILEOUT = "  << PARAMETERFILEOUT << endl;
 
          // ##########################
          // # Set histo layout style #
          // ##########################
          SetStyle();

          // ###################
          // # Read parameters #
          // ###################       
          Utility->ReadGenericParam(ParameterFILE);
          Utility->ReadFitStartingValues(ParameterFILE,&fitParam,&configParam,Utility->ParFileBlockN("fitValBins"));

          CTRLfitWRKflow = Utility->GetGenericParam("CtrlFitWrkFlow");

          cout << "Use MINOS: "                          << Utility->GetGenericParam("UseMINOS").c_str() << " (0 = false; 1 = true)" << endl;
          cout << "Control fit workflow: "               << CTRLfitWRKflow.c_str() << endl;
          // ###############################################################################################
          // # Read other parameters : this also allow to understand if the parameter file is well written #
          // ###############################################################################################
          LUMI = Utility->ReadLumi(ParameterFILE);
          if (Utility->WhatIsThis(ParameterFILE) == 0) cout << "\n[ExtractYield::main]\t@@@ I recognize that this is a DATA file @@@" << endl;
          else                                         cout << "\n[ExtractYield::main]\t@@@ I recognize that this is a Monte Carlo file @@@" << endl;


         // ###################
          // # Select fit type #
          // ###################
          if (
              (FitType == 106) ||
              (FitType == 206)
             )
            {
              NtplFile = new TFile(fileName.c_str(),"READ");
              if ( (FitType == 106))      theTree  = (TTree*) NtplFile->Get("ntuple");
              else if ((FitType == 206))      theTree  = (TTree*) NtplFile->Get("ntuple");
              NTuple   = new B0KstMuMuTreeContent( SampleType);
              NTuple->Init(SampleType);

	      // #################
	      // # Make datasets #
	      // #################
	      cout << "\n[ExtractYield::main]\t@@@ Making datasets @@@" << endl;
              clock_t start = clock();
	      MakeDatasets(NTuple,FitType,SampleType);
              clock_t end = clock();
              double time = (double) (end - start)/ CLOCKS_PER_SEC;
              cout << "\n[ExtractYield::main]\t@@@ Time Reading tree is\n" << time <<"s @@@" << endl;
         
	      // ##############################
	      // # Select the proper fit type #
	      // ##############################
	     
              // #############################
              // # 3D-fit P5P-Fl per q^2 bin #
	      // #############################
		  cout << "\n[ExtractYield::main]\t@@@ Now fit invariant mass, cos(theta_K) and cos(theta_l) per mumu q^2 bins @@@" << endl;

               if ((FitType == 206) )                      Iterative3AnglesFitq2Bins(SingleCandNTuple,
		  							                        	  ctL,
		  								                          ctK,
		  								                          phi,
		  								                          specBin,
		  								                          FitType,
		  								                          &q2Bins,
		  								                          &vecConstr);
               
	        else if ((FitType == 106)  )                Iterative3AnglesFitq2Bins(SingleCandNTuple_RejectPsi,
                                                                                                          ctL,
                                                                                                          ctK,
                                                                                                          phi,
                                                                                                          specBin,
                                                                                                          FitType,
                                                                                                          &q2Bins,
                                                                                                          &vecConstr);
	    }

	}
      else
	{
	  cout << "Wrong parameter: " << endl;
	  cout << "./ExtractYield [FitType] [input/output[if toy-MC]File.root] [noEffCorr yesEffCorr]" << endl;
	  cout << "               [q^2 bin to fit (0 - ...)]" << endl;

	 
	  return EXIT_FAILURE;
	}
    }
  else
    {
          cout << "Parameter missing: " << endl;
	  cout << "./ExtractYield [FitType] [input/output[if toy-MC]File.root] [noEffCorr yesEffCorr]" << endl;
	  cout << "               [q^2 bin to fit (0 - ...)]" << endl;

      cout << "\n@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@  Signa  @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@" << endl;
      cout << "FitType = 206: 3D P5p-Fl (cos(theta_K), cos(theta_l)) per q^2 bin on Gen-level" << endl;
      cout << "FitType = 205: 3D Si-AFB (cos(theta_K), cos(theta_l)) per q^2 bin on Gen-level" << endl;
      cout << "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@" << endl;
      
      return EXIT_FAILURE;
    }

  return 0;
}
