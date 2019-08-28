#ifndef B0KSTMUMUTREECONTENT_H
#define B0KSTMUMUTREECONTENT_H

#include "TTree.h"


class B0KstMuMuTreeContent
{
 public:
  
  B0KstMuMuTreeContent (int SampleType);

  void Init ( int SampleType );
  void ClearNTuple ( int SampleType );
  void MakeTreeBranches (TTree* theTree, int SampleType );
  void SetBranchAddresses (TTree* theTree, int SampleType );
  void CopyCandidate (B0KstMuMuTreeContent* NTupleIn, int SampleType);

  double  genQ;

  double tagged_mass;
  double tagB0, genSignal;
  double cos_theta_l, cos_theta_k, phi_kst_mumu;
  double mumuMass, mumuMassE;
 
 double gen_cos_theta_l, gen_cos_theta_k, gen_phi_kst_mumu; 
 private:

  void CopyScalars (B0KstMuMuTreeContent* NTupleIn, int SampleType);
};

#endif
