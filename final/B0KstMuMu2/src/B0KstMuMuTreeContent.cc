#include "../interface/B0KstMuMuTreeContent.h"
#include <iostream>

// ###################################
// # Gen-level MC, SampleType == 0   #
// # Reco MC & DATA, SampleType == 1 #
// ###################################


B0KstMuMuTreeContent::B0KstMuMuTreeContent ( int SampleType )
{
  if (SampleType == 0)
   {
    // ########################
    // # Gen-Level parameters #
    // ########################

    genQ  = 0;
    gen_cos_theta_l   =0;
    gen_cos_theta_k   =0;
    gen_phi_kst_mumu  =0;
    }
  else if (SampleType == 1)
    {
     // #############################
     // # reco MC & DATA parameters #
     // #############################
     tagB0 = 0;
     genSignal = 0;
     mumuMass  = 0;
     mumuMassE  =0;
     cos_theta_l = 0;
     cos_theta_k = 0;
     phi_kst_mumu = 0;
     tagged_mass = 0;
   }
}

void B0KstMuMuTreeContent::Init ( int SampleType )
{
 if (SampleType == 0)
   {
     // ########################
     // # Gen-level parameters #
     // ########################
     genQ  = 0;

     gen_cos_theta_l  = 0;
     gen_cos_theta_k  = 0;
     gen_phi_kst_mumu     = 0;
    }
  else if (SampleType == 1 )
    {
      // #############################
      // # Reco MC & DATA parameters #
      // #############################
      // ### tag information ###
      tagB0 = 0;
      genSignal = 0;

       // ### Three angulars ###
      cos_theta_l = 0;
      cos_theta_k = 0;
      phi_kst_mumu = 0;

      mumuMass = 0;
      mumuMassE = 0;
      tagged_mass = 0; 
  }

}
void B0KstMuMuTreeContent::ClearNTuple  (int SampleType )
{
  if (SampleType == 0)
    {
       // ########################
       // # Gen-Level Parameters #
       // ########################
       // ### B0 Parameter ###
       genQ     =0;
       gen_cos_theta_l  = 0;
       gen_cos_theta_k  = 0;
       gen_phi_kst_mumu = 0;
    }
  
  // #############################
  // # reco MC & DATA Parameters #
  // #############################
  else if (SampleType ==1)
    {
       cos_theta_l = 0;
       cos_theta_k = 0;
       phi_kst_mumu = 0;
       tagged_mass = 0;
       tagB0 =0;
       genSignal =0;       
       mumuMass = 0;
       mumuMassE = 0;
   }
}

void B0KstMuMuTreeContent::MakeTreeBranches (TTree* theTree, int SampleType)
{
 // ##########
 // # Gen MC #
 // ########## 
 if (SampleType == 0)
 { // ### dumion mass ###
  theTree->Branch("genQ",    &genQ,    "genQ");

  // ### angular parameter ###
  theTree->Branch("gen_cos_theta_l",    &gen_cos_theta_l,   "gen_cos_theta_l");
  theTree->Branch("gen_cos_theta_k",    &gen_cos_theta_k,   "gen_cos_theta_k");
  theTree->Branch("gen_phi_kst_mumu",   &gen_phi_kst_mumu,  "gen_phi_kst_mumu");
  }
 // ##################
 // # Reco MC & DATA #
 // ##################
 else if (SampleType == 1)
  {
   theTree->Branch("cos_theta_l", &cos_theta_l, "cos_theta_l");
   theTree->Branch("cos_theta_k", &cos_theta_k, "cos_theta_k");
   theTree->Branch("phi_kst_mumu", &phi_kst_mumu, "phi_kst_mumu");
   theTree->Branch("tagged_mass", &tagged_mass, "tagged_mass");
   theTree->Branch("mumuMass", &mumuMass, "mumuMass");
   theTree->Branch("mumuMassE", &mumuMassE, "mumuMassE");
   theTree->Branch("tagB0", &tagB0, "tagB0");
   theTree->Branch("genSiganl", &genSignal, "genSiganl");
  }
}

void B0KstMuMuTreeContent::SetBranchAddresses (TTree* theTree, int SampleType)
{
 // #############
 // # Gen-Level #
 // ############
 if (SampleType == 0)
  {
  theTree->SetBranchAddress("genQ", &genQ);

  theTree->SetBranchAddress("gen_cos_theta_l", &gen_cos_theta_l);
  theTree->SetBranchAddress("gen_cos_theta_k", &gen_cos_theta_k);
  theTree->SetBranchAddress("gen_phi_kst_mumu", &gen_phi_kst_mumu);
  }

  // ##################
  // # reco MC & DATA #
  // ##################
  else if (SampleType == 1)
   {
  theTree->SetBranchAddress("cos_theta_l", &cos_theta_l);
  theTree->SetBranchAddress("cos_theta_k", &cos_theta_k);
  theTree->SetBranchAddress("phi_kst_mumu", &phi_kst_mumu);
  theTree->SetBranchAddress("tagB0", &tagB0);
  theTree->SetBranchAddress("genSignal", &genSignal);
  theTree->SetBranchAddress("mumuMass", &mumuMass);
  theTree->SetBranchAddress("mumuMassE", &mumuMassE);
  theTree->SetBranchAddress("tagged_mass", &tagged_mass);
 } 
}
void B0KstMuMuTreeContent::CopyCandidate (B0KstMuMuTreeContent* NTupleIn, int SampleType)
{ 
  CopyScalars(NTupleIn, SampleType);
}

void B0KstMuMuTreeContent::CopyScalars (B0KstMuMuTreeContent* NTupleIn, int SampleType )
{
 // ################
 // # Gen-Level MC #
 // ################
 if (SampleType == 0)
  {
  genQ           = NTupleIn->genQ;

  gen_cos_theta_l     = NTupleIn->gen_cos_theta_l;
  gen_cos_theta_k     = NTupleIn->gen_cos_theta_k;
  gen_phi_kst_mumu    = NTupleIn->gen_phi_kst_mumu;
  }

  // ##################
  // # reco MC & DATA #
  // ##################
  else if (SampleType == 1)
   {
  cos_theta_l          = NTupleIn->cos_theta_l;
  cos_theta_k          = NTupleIn->cos_theta_k;
  phi_kst_mumu         = NTupleIn->phi_kst_mumu;
  tagB0               = NTupleIn->tagB0;
  genSignal           = NTupleIn->genSignal;tagged_mass          = NTupleIn->tagged_mass;
  mumuMass            = NTupleIn->mumuMass;
  mumuMassE           = NTupleIn->mumuMassE;
  }
}






