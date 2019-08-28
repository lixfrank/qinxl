#include <RooGenericPdf.h>
#include "RooRealVar.h"
#include "RooEffProd.h"
#include "RooDataSet.h"
#include "RooProdPdf.h"
#include "RooWorkspace.h"
#include "RooPolynomial.h"
#include "RooExponential.h"
#include "RooRandom.h"

#include <fstream>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <utility>
#include <sstream>

using std::cout;
using std::endl;
using std::string;
using std::stringstream;
using std::vector;
using namespace RooFit ;

#include <TROOT.h>
#include "TFile.h"
#include "TMath.h"

#include "Utils.h"

#define PARAMETERFILEIN "/python/ParameterFile.txt"
Utils* Utility;


int Enumber = 200;
float fl[9] = {0.70106, 0.80106, 0.72858, 0.60421, 0, 0.44971, 0, 0.36234, 0.33799};
float p1[9] = { -0.011948, -0.11812, -0.052155, -0.18852, 0, -0.22149, 0, -0.38442, -0.61077};
float p2[9] = {-0.34202, -0.19378, 0.2103, 0.41643, 0, 0.46365, 0, 0.43995, 0.37644};
float p3[9] = {-0.0066693, 0.020142, 0.0078555, -0.003102, 0,  -0.013429, 0, 0.010475, 0.014601};
float p4p[9] = {-0.099372, 0.057716, -0.064411, -0.001402, 0, 0.036453, 0, -0.033719, 0.016814};
float p5p[9] = {0.048159, -0.0098887, -0.030408, -0.0033749, 0,  0.026821, 0, -0.014805, 0.027373};
float p6p[9] = {-0.0026263, -0.018793,-0.044764, 0.0026991, 0, -0.0031762, 0, -0.0076682, -0.0078121 };
float p8p[9] = {0.019295, -0.024787, -0.034604, -0.0010124, 0, -0.034945, 0, -0.018051, 0.063683};

RooRealVar ctK("ctK","cos#theta_{K}",-1,1);
RooRealVar ctL("ctL","cos#theta_{L}",-1,1);
RooRealVar phi("phi","#phi",-TMath::Pi(),TMath::Pi());

RooArgSet VarsAng;

int sigYield[9] = {84,145,119,225,0,361,0,222,239};
int bkgYield[9] = {91,289,216,343,0,567,0,178,82};

vector<double> q2Bins;

void loadBinning() {
  Utility->Readq2Bins(Utility->MakeAnalysisPATH(PARAMETERFILEIN).c_str(),&q2Bins);
}

void GenerateSig ( int toyIndx, unsigned int q2BinIndx )
{
    RooAbsPdf* AnglesPDF = NULL;
    RooRealVar FlS("FlS","F_{L}", fl[q2BinIndx-1]);
    RooRealVar P1S("P1S","P_{1}", p1[q2BinIndx-1]);
    RooRealVar P2S("P2S","P_{2}", p2[q2BinIndx-1]);
    RooRealVar P3S("P3S","P_{3}", p3[q2BinIndx-1]);
    RooRealVar P4pS("P4pS","P_{4p}", p4p[q2BinIndx-1]);
    RooRealVar P5pS("P5pS","P_{5p}", p5p[q2BinIndx-1]);
    RooRealVar P6pS("P6pS","P_{6p}", p6p[q2BinIndx-1]);
    RooRealVar P8pS("P8pS","P_{8p}", p8p[q2BinIndx-1]);

    cout << "Fl:" << fl[q2BinIndx-1] << endl;
    cout << "P1:" << p1[q2BinIndx-1] << endl;
    cout << "P5p:" << p5p[q2BinIndx-1] << endl;


    VarsAng.add(FlS);
    VarsAng.add(P1S);
    VarsAng.add(P2S);
    VarsAng.add(P3S);
    VarsAng.add(P4pS);
    VarsAng.add(P5pS);
    VarsAng.add(P6pS);
    VarsAng.add(P8pS);
    VarsAng.add(ctL);
    VarsAng.add(ctK);
    VarsAng.add(phi);

    // #####################
    // # P-wave decay rate #
    // #####################
    stringstream myString;
    myString.clear(); myString.str("");
    myString << "(9/(32 *" << Utility->PI << ") * ( 1/2 * (1-" << "FlS" << ") * (1-" << "ctK" << "*" << "ctK" << ") * ( 1+ " << "ctL" << "*" << "ctL" << ") + 2 * " << "FlS" << "* " << "ctK" << "*" << "ctK" << " *( 1- " << "ctL" << "*" << "ctL" << ") +";
    myString << "0.5 * P1S" << " * (1-" << "FlS" << ") * (1- " << "ctK" << "*" << "ctK" << ") * ( 1- " << "ctL" << "*" << "ctL" << ") * cos (2 *"<< "phi"<< ") + ";
    myString << "P4pS" << "* 2 * cos(" << "phi" << ") * " << "ctK" << "*" << "ctL" << " * sqrt(" << "FlS" << "* (1-" << "FlS" << ") * ( 1-" << "ctL" << "*" << "ctL" << ") * (1-" << "ctK" << "*" << "ctK" << ")) + " ;
    myString << "P5pS" << "* 2 * cos(" << "phi" << ") * " << "ctK" << "* sqrt (" << "FlS" << "* (1-" << "FlS" << ") * ( 1-" << "ctL" << "*" << "ctL" << ") * (1-" << "ctK" << "*" << "ctK" << ")) - " ;
    myString << "P6pS" << "* 2 * sin(" << "phi" << ") * " << "ctK" << "* sqrt (" << "FlS" << "* (1-" << "FlS" << ") * ( 1-" << "ctL" << "*" << "ctL" << ") * (1-" << "ctK" << "*" << "ctK" << ")) + " ;
    myString << "P8pS" << "* 2 * sin(" << "phi" << ") * " << "ctK" << "*" << "ctL" << " * sqrt(" << "FlS" << "* (1-" << "FlS" << ") * ( 1-" << "ctL" << "*" << "ctL" << ") * (1-" << "ctK" << "*" << "ctK" << ")) + " ;
    myString << "P2S" << "* 2 * (1-" << "FlS" << ")*" << "ctL" << "* (1-" <<  "ctK" << "*" << "ctK" << ") -" ;
    myString << "P3S" << " * sin( 2 *" << "phi" << ")  * ( 1-" << "FlS" << ") * (1-" << "ctL" << "*" << "ctL" << ") * (1-" << "ctK" << "*" << "ctK" << ")))" ;
    RooGenericPdf* _AnglesPDF = new RooGenericPdf("_AnglesPDF",myString.str().c_str(),RooArgSet(VarsAng));

    // ############################
    // # Read efficiency function #
    // ############################
    RooAbsPdf* EffPdf_R = Utility->ReadRTEffPDF(q2BinIndx -1, 4);
    cout << "\n[ExtractYield::MakeAngWithEffPDF]\t@@@ The efficiency function from Alessio  @@@" << endl;

    AnglesPDF  = new RooEffProd("AngleS","Signal * Efficiency", *_AnglesPDF, *EffPdf_R);

    stringstream fileNameOutput;
    fileNameOutput.clear(); fileNameOutput.str("");
    fileNameOutput <<  "toyDatasets" << q2BinIndx -1 << "_" << Enumber << ".root";
    TFile* fout = new TFile(fileNameOutput.str().c_str(), "new");
    RooWorkspace* ws = (RooWorkspace*)fout->Get("ws");
    RooDataSet* dataSet[toyIndx+1];
    // #####################
    // # Generate Datasets #
    // #####################
    for(int toyid = 0; toyid < toyIndx; toyid++)
      {
        RooRandom::randomGenerator()->SetSeed(toyid+1);
        dataSet[toyid+1] = AnglesPDF->generate(RooArgSet(ctK,ctL,phi), Enumber);
        dataSet[toyid+1]->SetName(Form("toy%i%i",toyid+1,q2BinIndx));
        dataSet[toyid+1]->SetTitle(Form("Dataset for toy %i (bin %i)",toyid+1,q2BinIndx));
        dataSet[toyid+1]->Print("v");

        if ( ws==0 ) ws = new RooWorkspace("ws");
        ws->import(*dataSet[toyid+1]);
        ws->Write(0);
      }
    fout->Close();
    delete fout;
//    delete dataSet;


}


int main(int argc, char** argv)
{
  if (argc > 0)
  {
    int toyIndx=0;
    int q2BinIndx=0;

    if ( argc > 1 ) toyIndx = atoi(argv[1]);
    if ( argc > 2 ) q2BinIndx = atoi(argv[2]);


    Utility = new Utils(false);
    loadBinning();

    if (toyIndx<0) {
      cout<<"toy index must be greater than 0. FAILURE!"<<endl;
      return EXIT_FAILURE;
    }

    // do all q2Bins at once
    if (q2BinIndx==0)
      for (q2BinIndx=1; q2BinIndx<10; ++q2BinIndx) {
	    if (q2BinIndx==5 || q2BinIndx==7) continue;
	    GenerateSig(toyIndx, q2BinIndx);
      }


    else if (q2BinIndx>0 && q2BinIndx<10 && q2BinIndx!=5 && q2BinIndx!=7) {
	  GenerateSig(toyIndx, q2BinIndx);

    }
    else {
      cout<<"q2Bin must be greater than 0 and smaller than 10 and no 5 nor 7. FAILURE!"<<endl;
      return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
  }
  else
  {
    cout << "Parameter missing: " << endl;
    cout << "./" << argv[0] << " <toy index> " << " <q2Bin, 0:all (default 0)> " << endl;

    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
