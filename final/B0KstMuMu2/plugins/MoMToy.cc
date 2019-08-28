#include <TCanvas.h>
#include <TTree.h>
#include <TFile.h>
#include <TGaxis.h>
#include <TH1.h>
#include <Math/Functor.h>
#include <TMath.h>
#include <RooAbsPdf.h>
#include <RooWorkspace.h>
#include <RooDataSet.h>
#include <RooArgSet.h>
#include <RooRealVar.h>
#include <TROOT.h>
#include <TApplication.h>
#include <TSystem.h>

#include <sstream>
#include <time.h>
#include <ctime>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <stdio.h>

using std::cout;
using std::endl;
using std::stringstream;
using std::string;
using namespace RooFit;

// #############
// # Variables #
// #############
double q2Min = 0.;
double q2Max = 0.;
int q2Bin = -1;
double mumuMass2 = -1.;

double wight;
double M_6s;
double M_6c;
double f_1s;
double f_3;
double f_4;
double f_5;
double f_6s;
double f_6c;
double f_7;
double f_8;
double f_9;

double Mi;
// ##############
// # Parameters #
// ##############
double value;
double error;
double Fl = 0.;
double AFB = 0.;
double S3 = 0.;
double S4 = 0.;
double S5 = 0.;
double S6 = 0.;
double S7 = 0.;
double S8 = 0.;
double S9 = 0.;

double P1 = 0.;
double P2 = 0.;
double P3 = 0.;
double P4p = 0.;
double P5p = 0.;
double P6p = 0.;
double P8p = 0.;

// ####################
// # angular varables #
// ####################
double cos_theta_k = 0;
double cos_theta_l = 0;
double phi_kst_mumu = 0;

int idn = 501;

void CalValue(int q2Bin, string Paratype) {

    RooDataSet* toy[idn];
    stringstream myString;
    wight = 0.0;
    M_6s = 0.;
    M_6c = 0.;
    f_1s = 0.;
    f_3 = 0.;
    f_4 = 0.;
    f_5 = 0.;
    f_6s = 0.;
    f_6c = 0.;
    f_7 = 0.;
    f_8 = 0.;
    f_9 = 0.;

    TH1D *f1s = new TH1D("f1s", "f1s", 100, -100., 100.0);
    TH1D *m6s = new TH1D("m6s", "m6s", 200, -100.0, 100.0);
    TH1D *m6c = new TH1D("m6c", "m6c", 200, -100.0, 100.0);
    TH1D *f3 = new TH1D("f3", "f3", 200, -100.0, 100.0);
    TH1D *f4 = new TH1D("f4", "f4", 200, -100.0, 100.0);
    TH1D *f5 = new TH1D("f5", "f5", 200, -100.0, 100.0);
    TH1D *f7 = new TH1D("f7", "f7", 200, -100.0, 100.0);
    TH1D *f8 = new TH1D("f8", "f8", 200, -100.0, 100.0);
    TH1D *f9 = new TH1D("f9", "f9", 200, -100.0, 100.0);

    TH1D *fv = new TH1D("fv", "fv", 100, 0.65, 0.85);

    // #############
    // # open file #
    // #############
    myString.clear();
    myString.str("");
    myString << "q2bin" << q2Bin << Paratype << ".txt";
    FILE* fp = fopen(myString.str().c_str(),"w");

    // ############################
    // # Read Effcicency Function #
    // ############################
    int test = 4;

    myString.clear();
    myString.str("");
    myString << "../efficiency/projection/effProjection_sh" << test << "o_b" << q2Bin << "ct_25_25_25.root";
    cout << "\n[Moment::ReCalValue]\tTry to open" << myString.str().c_str() << endl;

    TFile *Eff = new TFile(myString.str().c_str(), "READ");
    RooWorkspace *w = (RooWorkspace *) Eff->Get("ws");
    RooAbsPdf *EffPdf = (RooAbsPdf *) w->function("projectedFunc");
    RooArgSet *param = EffPdf->getVariables();

    // ##################
    // # Read Data & MC #
    // ##################
    cout << "\n[Moment::CalValue]\t @@@ Making datasets @@@ " << endl;
    TFile *f = new TFile("/afs/cern.ch/user/l/llinwei/work2/B0KstMuMufull/B0KstMuMu2/plugins/toyDatasets2_200.root");
    RooWorkspace *te = (RooWorkspace *) f->Get("ws");
    for(int toyid = 1; toyid < idn; toyid++) {

        // # selecte the toyid #
        myString.clear();
        myString.str("");
        myString << "toy" << toyid << q2Bin +1 ;
        toy[toyid] = (RooDataSet*) te->data(myString.str().c_str());
       // toy[toyid]->Print("v");

        Int_t entries = toy[toyid]->sumEntries();
        cout << "\n[Moment::CalValue]\tTotal number of events in the sample: " << entries << " @@@" << endl;
        for (Int_t i = 0; i < entries; i++) {

            // ##################
            // # complete wight #
            // ##################
            RooRealVar *k = (RooRealVar *) toy[toyid]->get(i)->find("ctK");
            RooRealVar *l = (RooRealVar *) toy[toyid]->get(i)->find("ctL");
            RooRealVar *p = (RooRealVar *) toy[toyid]->get(i)->find("phi");
            cos_theta_k = k->getValV();
            cos_theta_l = l->getValV();
            phi_kst_mumu = p->getValV();
            param->setRealValue("ctK", cos_theta_k);
            param->setRealValue("ctL", cos_theta_l);
            param->setRealValue("phi", phi_kst_mumu);

            // ###############################
            // # define orthogonal functions #
            // ###############################
            f_1s = 1 - cos_theta_k * cos_theta_k;
            f_6s = (1 - cos_theta_k * cos_theta_k) * cos_theta_l;
            f_6c = cos_theta_k * cos_theta_k * cos_theta_l;
            f_3 = (1 - cos_theta_k * cos_theta_k) * (1 - cos_theta_l * cos_theta_l) * TMath::Cos(2 * phi_kst_mumu);
            f_4 = 4 * cos_theta_k * cos_theta_l * TMath::Cos(phi_kst_mumu) *
                  TMath::Sqrt((1 - cos_theta_k * cos_theta_k) * (1 - cos_theta_l * cos_theta_l));
            f_5 = 2 * cos_theta_k * TMath::Cos(phi_kst_mumu) *
                  TMath::Sqrt((1 - cos_theta_k * cos_theta_k) * (1 - cos_theta_l * cos_theta_l));
            f_7 = 2 * cos_theta_k * TMath::Sin(phi_kst_mumu) *
                  TMath::Sqrt((1 - cos_theta_k * cos_theta_k) * (1 - cos_theta_l * cos_theta_l));
            f_8 = 4 * cos_theta_k * cos_theta_l * TMath::Sin(phi_kst_mumu) *
                  TMath::Sqrt((1 - cos_theta_k * cos_theta_k) * (1 - cos_theta_l * cos_theta_l));
            f_9 = (1 - cos_theta_k * cos_theta_k) * (1 - cos_theta_l * cos_theta_l) * TMath::Sin(2 * phi_kst_mumu);
            M_6s = 3.0 * f_6s - 2.0 * f_6c;
            M_6c = 2.0 * (4 * f_6c - f_6s);
            if (EffPdf->getValV() < 0.0)
                cout << "Eff = " << EffPdf->getValV() << " " << cos_theta_k << "" << cos_theta_l << "" << phi_kst_mumu
                     << endl;
            wight = std::max(1 / (EffPdf->getVal()), 0.);
            f1s->Fill(f_1s, wight);
            m6s->Fill(M_6s, wight);
            m6c->Fill(M_6c, wight);
            f3->Fill(f_3, wight);
            f4->Fill(f_4, wight);
            f5->Fill(f_5, wight);
            f7->Fill(f_7, wight);
            f8->Fill(f_8, wight);
            f9->Fill(f_9, wight);
        }

        if (Paratype == "FlS") {
            value = 2.0 - 2.5 * f1s->GetMean();
            error = 2.5 * f1s->GetRMS() / TMath::Sqrt(f1s->GetEntries() - 1);
        } else if (Paratype == "AFBS") {
            value = 3.0 / 4.0 * (3.0 * m6s->GetMean() - 2.0 * m6c->GetMean());
            error = 3.0 / 4.0 * TMath::Sqrt(4.0 * TMath::Power(m6c->GetRMS(), 2) / (m6c->GetEntries() - 1) +
                                            9.0 * TMath::Power(m6s->GetRMS(), 2) / (m6s->GetEntries() - 1));
        } else if (Paratype == "S3S") {
            value = 25.0 / 8.0 * f3->GetMean();
            error = 25.0 / 8.0 * f3->GetRMS() / TMath::Sqrt(f3->GetEntries() - 1);
        } else if (Paratype == "S4S") {
            value = 25.0 / 8.0 * f4->GetMean();
            error = 25.0 / 8.0 * f4->GetRMS() / TMath::Sqrt(f4->GetEntries() - 1);
        } else if (Paratype == "S5S") {
            value = 2.5 * f5->GetMean();
            error = 2.5 * f5->GetRMS() / TMath::Sqrt(f5->GetEntries() - 1);
        } else if (Paratype == "S7S") {
            value = 2.5 * f7->GetMean();
            error = 2.5 * f7->GetRMS() / TMath::Sqrt(f7->GetEntries() - 1);
        } else if (Paratype == "S8S") {
            value = 25.0 / 8.0 * f8->GetMean();
            error = 25.0 / 8.0 * f8->GetRMS() / TMath::Sqrt(f8->GetEntries() - 1);
        } else if (Paratype == "S9S") {
            value = 25.0 / 8.0 * f9->GetMean();
            error = 25.0 / 8.0 * f9->GetRMS() / TMath::Sqrt(f9->GetEntries() - 1);
        } else if (Paratype == "P1S") {
            value = 25. / 4. * f3->GetMean() / (1. - (2.0 - 2.5 * f1s->GetMean()));
            error = TMath::Sqrt(TMath::Power(
                    2 / (1.0 - (2.0 - 2.5 * f1s->GetMean())) *
                    (25 / 8 * f3->GetRMS() / TMath::Sqrt(f3->GetEntries() - 1)),
                    2) + TMath::Power((25 / 4 * f3->GetMean() / TMath::Power((1.0 - (2.0 - 2.5 * f1s->GetMean())), 2)) *
                                      (2.5 * f1s->GetRMS() / TMath::Sqrt(f1s->GetEntries() - 1)), 2));
        } else if (Paratype == "P2S") {
            value = 0.5 * m6s->GetMean() / (1. - (2.0 - 2.5 * f1s->GetMean()));
            error = TMath::Sqrt(TMath::Power(1 / (2 * (1 - (2.0 - 2.5 * f1s->GetMean()))), 2) *
                                (4 * TMath::Power(m6c->GetRMS(), 2) / (m6c->GetEntries() - 1) +
                                 9 * TMath::Power(m6s->GetRMS(), 2) / (m6s->GetEntries() - 1)) + TMath::Power(
                    (3.0 * m6s->GetMean() - 2.0 * m6c->GetMean()) /
                    (2 * TMath::Power((1 - (2.0 - 2.5 * f1s->GetMean())), 2)), 2) * (TMath::Power(2.5 * f1s->GetRMS(),
                                                                                                  2) /
                                                                                     (f1s->GetEntries() - 1)));
        } else if (Paratype == "P3S") {
            value = -25.0 / 8.0 * f9->GetMean() / (1.0 - (2.0 - 2.5 * f1s->GetMean()));
            double a = TMath::Power(2.5 * f1s->GetRMS(), 2) / (f1s->GetEntries() - 1);
            double b = TMath::Power(25 / 8 * f9->GetRMS(), 2) / (f9->GetEntries() - 1);
            error = TMath::Sqrt(TMath::Power(1.0 / (1.0 - (2.0 - 2.5 * f1s->GetMean())), 2) * b +
                                TMath::Power(
                                        25 / 8 * f9->GetMean() / TMath::Power((1.0 - (2.0 - 2.5 * f1s->GetMean())), 2),
                                        2) * a);
        } else if (Paratype == "P4pS") {
            value = 4 * 25 / 8 * f4->GetMean() /
                    TMath::Sqrt((2.0 - 2.5 * f1s->GetMean()) * (1.0 - (2.0 - 2.5 * f1s->GetMean())));
            double a = TMath::Power(2.5 * f1s->GetRMS(), 2) / (f1s->GetEntries() - 1);
            double b = TMath::Power(25 / 8 * f4->GetRMS(), 2) / (f4->GetEntries() - 1);
            error = 4 * TMath::Sqrt(1 / ((2.0 - 2.5 * f1s->GetMean()) * (1 - (2.0 - 2.5 * f1s->GetMean()))) * b +
                                    TMath::Power((1 - 2 * (2.0 - 2.5 * f1s->GetMean())) * 25 / 8 * f4->GetMean(), 2) *
                                    a /
                                    (4 * TMath::Power((1 - (2.0 - 2.5 * f1s->GetMean())), 3)));
        } else if (Paratype == "P5pS") {
            double a =  2.0 - 2.5 * f1s->GetMean();
            value = 2.5 * f5->GetMean() / TMath::Sqrt(a * (1 - a));
            double b = TMath::Power(2.5 * f1s->GetRMS(), 2) / (f1s->GetEntries() - 1);
            double c = TMath::Power(2.5 * f5->GetRMS(), 2) / (f5->GetEntries() - 1);
            error = TMath::Sqrt(1 / (a * (1 - a)) * c +
                                TMath::Power((1 - 2 * a) * 2.5 * f5->GetMean(), 2) * b /
                                (4 * TMath::Power((1 - a), 3)));
        } else if (Paratype == "P6pS") {
            double a =  2.0 - 2.5 * f1s->GetMean();
            value = - 2.5 * f7->GetMean() / TMath::Sqrt(a * (1 - a));
            double b = TMath::Power(2.5 * f1s->GetRMS(), 2) / (f1s->GetEntries() - 1);
            double c = TMath::Power(2.5 * f7->GetRMS(), 2) / (f7->GetEntries() - 1);
            error = TMath::Sqrt(1 / (a * (1 - a)) * c +
                                TMath::Power((1 - 2 * a) * 2.5 * f7->GetMean(), 2) * b /
                                (4 * TMath::Power((1 - a), 3)));
        } else if (Paratype == "P8pS") {
            double a = 2.0 - 2.5 * f1s->GetMean();
            value = 2* 25 / 8 * f8->GetMean() / TMath::Sqrt(a * (1 - a));
            double b = TMath::Power(2.5 * f1s->GetRMS(), 2) / (f1s->GetEntries() - 1);
            double c = TMath::Power(25 / 8 * f8->GetRMS(), 2) / (f8->GetEntries() - 1);
            error = 2 * TMath::Sqrt(1 / (a * (1 - a)) * c +
                                TMath::Power((1 - 2 * a) * 25 / 8 * f8->GetMean(), 2) * b /
                                (4 * TMath::Power((1 - a), 3)));
        }

        // #########################
        // # write value and error #
        // #########################
        fprintf(fp,"%f %f\n",value,error);
        fv->Fill(value);
    }
    fclose(fp);

    // #########
    // # plots #
    // #########
    myString.clear();
    myString.str("");
    myString << Paratype;

    TCanvas* c = new TCanvas("c","c");
    c->cd();
    fv->GetXaxis()->SetTitle(myString.str().c_str());
    fv->Draw();
    myString << ".pdf";
    c->Update();
    c->Print(myString.str().c_str());
}

int main(int argc, char **argv) {

    string ParaType = "";
    int q2BinIn = -1;

    // ###################################
    // # Check that Parameter is correct #
    // ###################################
    ParaType = argv[1];
    if ((ParaType != "FlS") && (ParaType != "AFBS") && (ParaType != "P1S") && (ParaType != "P2S") &&
        (ParaType != "P3S") && (ParaType != "P4pS") && (ParaType != "P5pS") && (ParaType != "P6pS") &&
        (ParaType != "P8pS") && (ParaType != "S3S") && (ParaType != "S4S") && (ParaType != "S5S") &&
        (ParaType != "S6S") && (ParaType != "S7S") && (ParaType != "S8S") && (ParaType != "S9S")) {
        cout << "\n[Moment::main]\tIncorrect Parameters\t" << ParaType << endl;
        return EXIT_FAILURE;
    }
    else
    cout << "\n[Moment::main]\tParameter Type = \t" << ParaType << endl;

    // ##########################
    // # Check Q2bin is correct #
    // ##########################
    q2BinIn = atoi(argv[2]);
    if ((q2BinIn<0) || (q2BinIn>8) || (q2BinIn ==4) || (q2BinIn == 6))
    {
        cout << "\n[Moment::main]\tIncorrect q2bin\t" << q2BinIn << endl;
        return EXIT_FAILURE;
    }
    else
        cout << "\n[Moment::main]\tq2bin = \t" << q2BinIn << endl;
    CalValue(q2BinIn, ParaType);


}
