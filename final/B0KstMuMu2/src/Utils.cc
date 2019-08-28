
#include "../interface/Utils.h"
#include "../interface/ReadParameters.h"

#include <TAxis.h>
#include <TMath.h>
#include <TFile.h>

#include <RooWorkspace.h>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <fstream>

// ####################
// # Global constants #
// ####################
#define YvalueOutsideLimits 20.0 // Value given to bins with zero error in order not to show them
#define ANALYPATH "ANALYPATH"


Utils::Utils (bool rightFlavorTag)
{
  muonMass     = 0.10565837;
  pionMass     = 0.13957018;
  kaonMass     = 0.493677;
  kstMass      = 0.896;
  B0Mass       = 5.27958;
  JPsiMass     = 3.096916;
  PsiPMass     = 3.686109;

  JPsiBF       =  7.87e-5; // B0 --> J/psi(mu+mu-) K*0          (1.32+/-0.06 * 5.96+/-0.03)
  JPsiKpiBF    =  5.25e-5; // B0 --> J/psi(mu+mu-) K*0(K+pi-)   (1.32+/-0.06 * 5.96+/-0.03 * 2/3)

  KstMuMuBF    =  1.05e-6; // B0 --> K*0 mu+mu-
  KstKpiMuMuBF =  7e-7;    // B0 --> K*0(K+pi-) mu+mu-          (1.05+/-0.1 * 2/3)

  PsiPBF       = 47.72e-7; // B0 --> psi(2S)(mu+mu-) K*0        (6.04+/-0.4 * 7.9+/-0.9)
  PsiPKpiBF    = 31.81e-7; // B0 --> psi(2S)(mu+mu-) K*0(K+pi-) (6.04+/-0.4 * 7.9+/-0.9 * 2/3)

  muonMassErr  = 3.5e-9;
  pionMassErr  = 3.5e-7;
  kaonMassErr  = 1.6e-5;
  B0MassErr    = 1.7e-4;
  kstSigma     = 0.05;

  nFitParam    = 16;
  nConfigParam = 4;
  nFitObserv   = 8; // FL P5p P1 P2 P3 P4p P6p P8p S3 S4 S5 AFB S7 S8 S9 

  NcoeffThetaL = 6;
  NcoeffThetaK = 4;
  NcoeffPhi    = 4;

  PI = 3.141592653589793;

  ProbThreshold = 0.0; // Confidence Level for accepting the null hypothesis: "the two mass hypothesis are statistically indistinguishable"
  // if (F-test < ProbThreshold) --> accept the null hypothesis
  // if (F-test > ProbThreshold) --> reject the null hypothesis
  scrambleFraction = 0.0; // Fraction of events with random CP-tagging
  KstMassShape = new TF1("KstMassShape",
			 "2*sqrt(2)*[0]*[1]* sqrt([0]*[0] * ([0]*[0] + [1]*[1])) / (TMath::Pi()*sqrt([0]*[0] + sqrt([0]*[0] * ([0]*[0] + [1]*[1])))) / ((x*x - [0]*[0]) * (x*x - [0]*[0]) + [0]*[0]*[1]*[1])",
			 0.0,kstMass*2.);
  // Breit-Wigner distribution:
  // [0]: mass of the resonance
  // [1]: width of the resonance
  KstMassShape->SetParameter(0,kstMass);
  KstMassShape->SetParameter(1,kstSigma);
  KstMassShape->SetParNames("Mean","Width");

  // Define whether to compute the efficiency with good-tagged or mis-tagged events
  RIGHTflavorTAG = rightFlavorTag;

  // ###############################################
  // # ===> Define codes to identify MC type <===  #
  // #                                             #
  // # 1 = B0    --> K*0(K+pi-) mu+mu-             #
  // # 2 = B0bar --> K*0bar(K-pi+) mu+mu-          #
  // #                                             #
  // # 3 = B0    --> K*0(K+pi-) J/psi(mu+mu-)      #
  // # 4 = B0bar --> K*0bar(K-pi+) J/psi(mu+mu-)   #
  // #                                             #
  // # 5 = B0    --> K*0(K-pi+) psi(2S)(mu+mu-)    #
  // # 6 = B0bar --> K*0bar(K-pi+) psi(2S)(mu+mu-) #
  // ###############################################
  B0ToKstMuMu  = 1;
  B0ToJPsiKst  = 3;
  B0ToPsi2SKst = 5;

  // ################################
  // # Print out internal variables #
  // ################################
  std::cout << "\n@@@@@@ Utils class settings : private @@@@@@" << std::endl;
  std::cout << "Analysis environment variable: " << ANALYPATH << std::endl;
  std::cout << "nFitObserv: "                << nFitObserv << std::endl;
  std::cout << "ProbThreshold: "             << ProbThreshold << std::endl;
  std::cout << "scrambleFraction: "          << scrambleFraction << std::endl;

  std::cout << "\n@@@@@@ Utils class settings : public  @@@@@@" << std::endl;
  std::cout << "NcoeffThetaL: "              << NcoeffThetaL << std::endl;
  std::cout << "NcoeffThetaK: "              << NcoeffThetaK << std::endl;
  std::cout << "NcoeffPhi: "                 << NcoeffPhi << std::endl;
  std::cout << "RIGHTflavorTAG: "            << RIGHTflavorTAG << std::endl;
  std::cout << "B0ToKstMuMu: "               << B0ToKstMuMu << std::endl;
  std::cout << "B0ToJPsiKst: "               << B0ToJPsiKst << std::endl;
  std::cout << "B0ToPsi2SKst: "              << B0ToPsi2SKst << std::endl;
  std::cout << "nFitParam: "                 << nFitParam << std::endl;
  std::cout << "nConfigParam: "              << nConfigParam << std::endl;
}

double Utils::computeInvMass (double Px1,
                              double Py1,
                              double Pz1,
                              double mass1,
                              double Px2,
                              double Py2,
                              double Pz2,
                              double mass2,
                              double Px3,
                              double Py3,
                              double Pz3,
                              double mass3)
{
  double Energy1 = sqrt(Px1*Px1 + Py1*Py1 + Pz1*Pz1 + mass1*mass1);
  double Energy2 = sqrt(Px2*Px2 + Py2*Py2 + Pz2*Pz2 + mass2*mass2);
  double Energy3 = sqrt(Px3*Px3 + Py3*Py3 + Pz3*Pz3 + mass3*mass3);
  return sqrt((Energy1+Energy2+Energy3) * (Energy1+Energy2+Energy3) - ((Px1+Px2+Px3) * (Px1+Px2+Px3) + (Py1+Py2+Py3) * (Py1+Py2+Py3) + (Pz1+Pz2+Pz3) * (Pz1+Pz2+Pz3)));
}

void Utils::ReadAllBins (std::string fileName, std::vector<double>* q2Bins)
// ##########################
// # signalType = "goodtag" #
// # signalType = "mistag"  #
// ##########################
{
  std::vector<std::string> ParVector;
  ReadParameters* ParameterFile = new ReadParameters(fileName.c_str());


  // #################
  // # Read q^2 bins #
  // #################
  std::cout << "\n[Utils::ReadAllBins]\tAll bins from file : " << fileName.c_str() << std::endl;
  ParVector.clear();
  ParameterFile->ReadFromFile(ParFileBlockN("q2"),&ParVector);
  for (unsigned int i = 0; i < ParVector.size(); i++)
    {
      q2Bins->push_back(atof(ParVector[i].c_str()));
      std::cout << "Bin " << i << "\t" << q2Bins->back() << std::endl;
    }

  ParVector.clear();
  delete ParameterFile;
}

void Utils::Readq2Bins (std::string fileName, std::vector<double>* q2Bins)
{
  std::vector<std::string>* ParVector = NULL;
  ReadParVsq2Bins(fileName,"q2",&ParVector);


  std::cout << "\n[Utils::Readq2Bins]\tq^2 bins from file : " << fileName.c_str() << std::endl;
  for (unsigned int i = 0; i < ParVector->size(); i++)
    {
      q2Bins->push_back(atof(ParVector->operator[](i).c_str()));
      std::cout << "Bin " << i << "\t" << q2Bins->back() << std::endl;
    }
  
    
  ParVector->clear();
  delete ParVector;
}

double Utils::ReadLumi (std::string fileName)
{
  double val = 0.0;
  std::vector<std::string> ParVector;
  ReadParameters* ParameterFile = new ReadParameters(fileName.c_str());


  // ##############################
  // # Read integrated luminosity #
  // ##############################
  ParameterFile->ReadFromFile(ParFileBlockN("lumi"),&ParVector);
  for (unsigned int i = 0; i < ParVector.size(); i++) val = val + atof(ParVector[i].c_str());
  std::cout << "\n[Utils::ReadLumi]\t@@@ Recorded luminosity: " << val << " fb-1 @@@" << std::endl;


  ParVector.clear();
  delete ParameterFile;
  return val;
}

void Utils::ReadFitStartingValues (std::string fileName, std::vector<std::vector<std::string>*>* vecParam, std::vector<std::vector<unsigned int>*>* configParam, const unsigned int dataBlockN)
{
  std::vector<std::string> ParVector;
  ReadParameters* ParameterFile = new ReadParameters(fileName.c_str());


  // ###################
  // # Clear vector(s) #
  // ###################
  vecParam->clear();
  configParam->clear();


  // ############################
  // # Read fit-starting values #
  // ############################
  ParameterFile->ReadFromFile(dataBlockN,&ParVector);

  for (unsigned int j = 0; j < nFitParam; j++) vecParam->push_back(new std::vector<std::string>);
  for (unsigned int j = 0; j < nConfigParam; j++) configParam->push_back(new std::vector<unsigned int>);
  
  for (unsigned int i = 0; i < ParVector.size(); i=i+nFitParam+nConfigParam)
    {

      std::cout << "\nRead set-" << static_cast<int>(static_cast<double>(i)/static_cast<double>(nFitParam+nConfigParam)) << " of fit-parameter starting values" << std::endl;

      for (unsigned int j = 0; j < nFitParam; j++)
	{
	  std::stringstream rawString(ParVector[i+j]);
	  vecParam->operator[](j)->push_back(rawString.str());
	  std::cout << "Fit parameter-" << j << " value: " << vecParam->operator[](j)->back() << std::endl;
	}
      for (unsigned int j = 0; j < nConfigParam; j++)
	{
	  configParam->operator[](j)->push_back(atoi(ParVector[i+nFitParam+j].c_str()));
	  std::cout << "Config. parameter-" << j << " value: " << configParam->operator[](j)->back() << std::endl;
	}
    }

  
  ParVector.clear();
  delete ParameterFile;
}

void Utils::ReadParVsq2Bins (std::string fileName, std::string praName, std::vector<std::string>** vecParam)
{
  ReadParameters* ParameterFile = new ReadParameters(fileName.c_str());

  if (*vecParam == NULL) *vecParam = new std::vector<std::string>;
  else                   (*vecParam)->clear();

  ParameterFile->ReadFromFile(ParFileBlockN(praName.c_str()),*vecParam);

  std::cout << "\n[Utils::ReadParVsq2Bins]\tReading parameters vs q^2 from file : " << fileName << std::endl;
  for (unsigned int i = 0; i < (*vecParam)->size(); i++)
    std::cout << "Parameter value and errors for q2 bin " << i << ": " << (*vecParam)->operator[](i) << std::endl;
  
  delete ParameterFile;
}

int Utils::WhatIsThis (std::string fileName)
{
  int val = 0;
  std::vector<std::string> ParVector;
  ReadParameters* ParameterFile = new ReadParameters(fileName.c_str());


  // ###############################
  // # Read data type (Data or MC) #
  // ###############################
  ParameterFile->ReadFromFile(ParFileBlockN("dtype"),&ParVector);
  for (unsigned int i = 0; i < ParVector.size(); i++) val = atoi(ParVector[i].c_str());


  ParVector.clear();
  delete ParameterFile;
  return val;
}

void Utils::SaveFitValues (std::string fileName, std::vector<std::string>* vecParStr, int indx, std::string howOpen, std::string str)
// #################################################
// # If indx == -1 then use str within default str #
// # If indx == -2 then use only str               #
// #################################################
{
  std::stringstream myString;
  
  myString.clear(); myString.str("");
  if      (indx == -1) myString << "#@@@@@@@@@@@@@@@@@@@@@@@@@@@" << str.c_str()   <<  "@@@@@@@@@@@@@@@@@@@@@@@@@@@#";
  else if (indx == -2) myString << str.c_str();
  else                 myString << "#@@@@@@@@@@@@@@@@@@@@@@@@@@@ q^2 bin " << indx << " @@@@@@@@@@@@@@@@@@@@@@@@@@@#";
  vecParStr->insert(vecParStr->begin(),myString.str());

  vecParStr->insert(vecParStr->end(),"");
  vecParStr->insert(vecParStr->end(),"");

  ReadParameters* ParameterFile = new ReadParameters(fileName.c_str(),howOpen.c_str());
  ParameterFile->WriteToFile(vecParStr);
  delete ParameterFile;
}


std::string Utils::MakeAnalysisPATH (std::string relativePath)
{
  std::stringstream myString;
  char* absolutePath = getenv(ANALYPATH);

  if (absolutePath == NULL)
    {
      std::cout << "[Utils::MakeAnalysisPATH]\tAnalysis environment variable " << ANALYPATH << " not defined" << std::endl;
      exit (EXIT_FAILURE);
    }
  myString.clear(); myString.str("");
  myString << absolutePath << relativePath;

  return myString.str();
}


unsigned int Utils::ParFileBlockN (std::string blockName)
{
  // #########################
  // # Parameter file blocks #
  // #########################
  if      (blockName == "q2")             return 1;

  else if (blockName == "fitValBins")     return 2;

  else if (blockName == "genericpar")     return 3;
  else if (blockName == "lumi")           return 4;
  else if (blockName == "dtype")          return 5;

  std::cout << "[Utils::ParFileBlockN]\tError wrong index name : " << blockName << std::endl;
  exit (EXIT_FAILURE);
}

unsigned int Utils::GetFitParamIndx (std::string varName)
{
  if      (varName == "Fl")            return 0;
  else if (varName == "P1")            return 1;
  else if (varName == "P2")            return 2;
  else if (varName == "P3")            return 3;
  else if (varName == "P4p")           return 4;
  else if (varName == "P5p")           return 5;
  else if (varName == "P6p")           return 6;
  else if (varName == "P8p")           return 7;
 
  else if (varName == "AFBS")           return 8;
  else if (varName == "S3S")            return 9;
  else if (varName == "S4S")            return 10;
  else if (varName == "S5S")            return 11;
  else if (varName == "S7S")            return 12;
  else if (varName == "S8S")            return 13;
  else if (varName == "S9S")            return 14;

  else if (varName == "fracMisTag")    return 15;
  std::cout << "[Utils::GetFitParamIndx]\tError wrong index name : " << varName << std::endl;
  exit (EXIT_FAILURE);
}

unsigned int Utils::GetConfigParamIndx (std::string varName)
{
  if      (varName == "SigType")      return 0;
  else if (varName == "PeakBkgType")  return 1;
  else if (varName == "CombBkgType")  return 2;
  else if (varName == "MistagType")   return 3;

  std::cout << "[Utils::GetConfigParamIndx]\tError wrong index name : " << varName << std::endl;
  exit (EXIT_FAILURE);
}

void Utils::ReadGenericParam (std::string fileName)
// #########################################
// # UseMINOS            = GenericPars[0]  #
// # B0MassIntervalLeft  = GenericPars[1] #
// # B0MassIntervalRight = GenericPars[2] #
// # S-P Wave            = GenericPars[3] #
// # CtrlFitWrkFlow      = GenericPars[4]  #
// #########################################
{
  std::vector<std::string> ParVector;
  ReadParameters* ParameterFile = new ReadParameters(fileName.c_str());


  // ###########################
  // # Read generic parameters #
  // ###########################
  std::cout << "\n[Utils::ReadGenericParam]\tGeneric parameters from file : " << fileName.c_str() << std::endl;
  ParVector.clear();
  ParameterFile->ReadFromFile(ParFileBlockN("genericpar"),&ParVector);
  for (unsigned int i = 0; i < ParVector.size(); i++)
    {
      GenericPars.push_back(ParVector[i]);
      std::cout << "Generic parameter #" << i << " from config file : " << GenericPars[i].c_str() << std::endl;
    }


  ParVector.clear();
  delete ParameterFile;
}

bool Utils::SetGenericParam (std::string parName, std::string val)
{
  if      (parName == "UseMINOS")            GenericPars[0]  = val;
  else if (parName == "B0MassIntervalLeft")  GenericPars[1] = val;
  else if (parName == "B0MassIntervalRight") GenericPars[2] = val;
  else if (parName == "UseSPwave")           GenericPars[3] = val;
  else if (parName == "CtrlFitWrkFlow")      GenericPars[4]  = val;
  else if (parName == "NSigmaB0")            GenericPars[5] = val;
  else if (parName == "NSigmaB0S")           GenericPars[6] = val;
  else if (parName == "NSigmaB0B")           GenericPars[7] = val;
  else if (parName == "NSigmaPsi")           GenericPars[8] = val;
  else if (parName == "B&psiMassJpsiLo")     GenericPars[9] = val;
  else if (parName == "B&psiMassJpsiHi")     GenericPars[10] = val;
  else if (parName == "B&psiMassPsiPLo")     GenericPars[11] = val;
  else if (parName == "B&psiMassPsiPHi")     GenericPars[12] = val;
  else if (parName == "CtrlMisTagWrkFlow")   GenericPars[13]  = val;
  else return false;

  return true;
}

std::string Utils::GetGenericParam (std::string parName)
{
  if      (parName == "UseMINOS")            return GenericPars[0];
  else if (parName == "B0MassIntervalLeft")  return GenericPars[1];
  else if (parName == "B0MassIntervalRight") return GenericPars[2];
  else if (parName == "UseSPwave")           return GenericPars[3];
  else if (parName == "CtrlFitWrkFlow")      return GenericPars[4];
  else if (parName == "NSigmaB0")            return GenericPars[5];
  else if (parName == "NSigmaB0S")           return GenericPars[6];
  else if (parName == "NSigmaB0B")           return GenericPars[7];
  else if (parName == "NSigmaPsi")           return GenericPars[8];
  else if (parName == "B&psiMassJpsiLo")     return GenericPars[9];
  else if (parName == "B&psiMassJpsiHi")     return GenericPars[10];
  else if (parName == "B&psiMassPsiPLo")     return GenericPars[11];
  else if (parName == "B&psiMassPsiPHi")     return GenericPars[12];
  else if (parName == "CtrlMisTagWrkFlow")   return GenericPars[13];
  else
    {
      std::cout << "[Utils::GetGenericParam]\tGeneric parameter not valid : " << parName << std::endl;
      exit (EXIT_FAILURE);
    }
}

bool Utils::PsiRejection (double myB0Mass, double myMuMuMass, double myMuMuMassE, std::string seleType, bool B0andPsiCut)
  // ###########################
  // # seleType == "keepJpsi"  #
  // # seleType == "keepPsiP"  #
  // # seleType == "rejectPsi" #
  // # seleType == "keepPsi"   #
  // ###########################
{
  // ####################################################################
  // # This method is used together with the method: "ReadGenericParam" #
  // ####################################################################
  if (seleType == "keepJpsi")
     {
       if ((fabs(myMuMuMass - JPsiMass) < atof(GetGenericParam("NSigmaPsi").c_str()) * myMuMuMassE) && (B0andPsiCut == false))
       return true;
     }
  else if (seleType == "keepPsiP")
     {
       if ((fabs(myMuMuMass - PsiPMass) < atof(GetGenericParam("NSigmaPsi").c_str()) * myMuMuMassE) && (B0andPsiCut == false))
       return true;
     }
  else if (seleType == "rejectPsi")
     {
       if ((fabs(myMuMuMass - JPsiMass) > atof(GetGenericParam("NSigmaPsi").c_str()) * myMuMuMassE) &&
           (fabs(myMuMuMass - PsiPMass) > atof(GetGenericParam("NSigmaPsi").c_str()) * myMuMuMassE) &&

           ((B0andPsiCut == false) ||
            ((B0andPsiCut == true) &&
             (((myMuMuMass < JPsiMass) &&
               (!((fabs((myB0Mass - B0Mass) - (myMuMuMass - JPsiMass)) < atof(GetGenericParam("B&psiMassJpsiLo").c_str()))    ||
                  (fabs((myB0Mass - B0Mass) - (myMuMuMass - PsiPMass)) < atof(GetGenericParam("B&psiMassPsiPLo").c_str()))))) ||
               
              ((myMuMuMass > PsiPMass) &&
               (!((fabs((myB0Mass - B0Mass) - (myMuMuMass - JPsiMass)) < atof(GetGenericParam("B&psiMassJpsiHi").c_str()))    ||
                  (fabs((myB0Mass - B0Mass) - (myMuMuMass - PsiPMass)) < atof(GetGenericParam("B&psiMassPsiPHi").c_str()))))) ||

              ((myMuMuMass > JPsiMass) && (myMuMuMass < PsiPMass) &&
               (!((fabs((myB0Mass - B0Mass) - (myMuMuMass - JPsiMass)) < atof(GetGenericParam("B&psiMassJpsiHi").c_str()))     ||
                  (fabs((myB0Mass - B0Mass) - (myMuMuMass - PsiPMass)) < atof(GetGenericParam("B&psiMassPsiPLo").c_str())))))))))

         return true;
     }
   else if (seleType == "keepPsi")
    {
       if (((fabs(myMuMuMass - JPsiMass) < atof(GetGenericParam("NSigmaPsi").c_str()) * myMuMuMassE)  ||
            (fabs(myMuMuMass - PsiPMass) < atof(GetGenericParam("NSigmaPsi").c_str()) * myMuMuMassE)) && (B0andPsiCut == false))
       return true;
    }
   else
    {
       std::cout << "[Utils::PsiRejection]\tSelection type not valid : " << seleType << std::endl;
       exit (EXIT_FAILURE);
    }

   return false;

}

double* Utils::MakeBinning (std::vector<double>* STLvec)
{
  double* vec = new double[STLvec->size()];

  for (unsigned int i = 0; i < STLvec->size(); i++) vec[i] = STLvec->operator[](i);

  return vec;
}

#if ROOFIT
using namespace RooFit;
/*
RooBernsteinEffi* Utils::GetFitReff   (unsigned int q2Indx)
{  
  std::stringstream myString;
  RooBernsteinEffi* Effi =NULL; 
 
  myString.clear(); myString.str("");
  myString << "/afs/cern.ch/user/l/llinwei/work2/B0KstMuMufull/B0KstMuMu2/efficiency/global/saveEffi-2016-Q2Bin-" << q2Indx  <<"-w.root";
  std::cout << myString.str().c_str() << std::endl;
  TFile* file= new TFile(myString.str().c_str(),"READ");
  RooWorkspace* ws = (RooWorkspace*)file->Get("w"); 
  Effi = (RooBernsteinEffi*) ws->function("Effi");
  return Effi;
}
*/
#endif
