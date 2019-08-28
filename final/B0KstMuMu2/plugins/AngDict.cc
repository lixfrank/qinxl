// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME AngDict

/*******************************************************************/
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#define G__DICTIONARY
#include "RConfig.h"
#include "TClass.h"
#include "TDictAttributeMap.h"
#include "TInterpreter.h"
#include "TROOT.h"
#include "TBuffer.h"
#include "TMemberInspector.h"
#include "TInterpreter.h"
#include "TVirtualMutex.h"
#include "TError.h"

#ifndef G__ROOT
#define G__ROOT
#endif

#include "RtypesImp.h"
#include "TIsAProxy.h"
#include "TFileMergeInfo.h"
#include <algorithm>
#include "TCollectionProxyInfo.h"
/*******************************************************************/

#include "TDataMember.h"

// Since CINT ignores the std namespace, we need to do so in this file.
namespace std {} using namespace std;

// Header files passed as explicit arguments
#include "../interface/PdfRT.h"
#include "../interface/PdfWT.h"
#include "../interface/DecayRate.h"
#include "../interface/PdfSigAng.h"

// Header files passed via #pragma extra_include

namespace ROOT {
   static void *new_PdfRT(void *p = 0);
   static void *newArray_PdfRT(Long_t size, void *p);
   static void delete_PdfRT(void *p);
   static void deleteArray_PdfRT(void *p);
   static void destruct_PdfRT(void *p);
   static void streamer_PdfRT(TBuffer &buf, void *obj);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::PdfRT*)
   {
      ::PdfRT *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::PdfRT >(0);
      static ::ROOT::TGenericClassInfo 
         instance("PdfRT", ::PdfRT::Class_Version(), "interface/PdfRT.h", 23,
                  typeid(::PdfRT), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::PdfRT::Dictionary, isa_proxy, 16,
                  sizeof(::PdfRT) );
      instance.SetNew(&new_PdfRT);
      instance.SetNewArray(&newArray_PdfRT);
      instance.SetDelete(&delete_PdfRT);
      instance.SetDeleteArray(&deleteArray_PdfRT);
      instance.SetDestructor(&destruct_PdfRT);
      instance.SetStreamerFunc(&streamer_PdfRT);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::PdfRT*)
   {
      return GenerateInitInstanceLocal((::PdfRT*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::PdfRT*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_PdfWT(void *p = 0);
   static void *newArray_PdfWT(Long_t size, void *p);
   static void delete_PdfWT(void *p);
   static void deleteArray_PdfWT(void *p);
   static void destruct_PdfWT(void *p);
   static void streamer_PdfWT(TBuffer &buf, void *obj);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::PdfWT*)
   {
      ::PdfWT *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::PdfWT >(0);
      static ::ROOT::TGenericClassInfo 
         instance("PdfWT", ::PdfWT::Class_Version(), "interface/PdfWT.h", 23,
                  typeid(::PdfWT), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::PdfWT::Dictionary, isa_proxy, 16,
                  sizeof(::PdfWT) );
      instance.SetNew(&new_PdfWT);
      instance.SetNewArray(&newArray_PdfWT);
      instance.SetDelete(&delete_PdfWT);
      instance.SetDeleteArray(&deleteArray_PdfWT);
      instance.SetDestructor(&destruct_PdfWT);
      instance.SetStreamerFunc(&streamer_PdfWT);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::PdfWT*)
   {
      return GenerateInitInstanceLocal((::PdfWT*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::PdfWT*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_DecayRate(void *p = 0);
   static void *newArray_DecayRate(Long_t size, void *p);
   static void delete_DecayRate(void *p);
   static void deleteArray_DecayRate(void *p);
   static void destruct_DecayRate(void *p);
   static void streamer_DecayRate(TBuffer &buf, void *obj);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::DecayRate*)
   {
      ::DecayRate *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::DecayRate >(0);
      static ::ROOT::TGenericClassInfo 
         instance("DecayRate", ::DecayRate::Class_Version(), "interface/DecayRate.h", 22,
                  typeid(::DecayRate), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::DecayRate::Dictionary, isa_proxy, 16,
                  sizeof(::DecayRate) );
      instance.SetNew(&new_DecayRate);
      instance.SetNewArray(&newArray_DecayRate);
      instance.SetDelete(&delete_DecayRate);
      instance.SetDeleteArray(&deleteArray_DecayRate);
      instance.SetDestructor(&destruct_DecayRate);
      instance.SetStreamerFunc(&streamer_DecayRate);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::DecayRate*)
   {
      return GenerateInitInstanceLocal((::DecayRate*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::DecayRate*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_PdfSigAng(void *p = 0);
   static void *newArray_PdfSigAng(Long_t size, void *p);
   static void delete_PdfSigAng(void *p);
   static void deleteArray_PdfSigAng(void *p);
   static void destruct_PdfSigAng(void *p);
   static void streamer_PdfSigAng(TBuffer &buf, void *obj);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::PdfSigAng*)
   {
      ::PdfSigAng *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::PdfSigAng >(0);
      static ::ROOT::TGenericClassInfo 
         instance("PdfSigAng", ::PdfSigAng::Class_Version(), "interface/PdfSigAng.h", 23,
                  typeid(::PdfSigAng), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::PdfSigAng::Dictionary, isa_proxy, 16,
                  sizeof(::PdfSigAng) );
      instance.SetNew(&new_PdfSigAng);
      instance.SetNewArray(&newArray_PdfSigAng);
      instance.SetDelete(&delete_PdfSigAng);
      instance.SetDeleteArray(&deleteArray_PdfSigAng);
      instance.SetDestructor(&destruct_PdfSigAng);
      instance.SetStreamerFunc(&streamer_PdfSigAng);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::PdfSigAng*)
   {
      return GenerateInitInstanceLocal((::PdfSigAng*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::PdfSigAng*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

//______________________________________________________________________________
atomic_TClass_ptr PdfRT::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *PdfRT::Class_Name()
{
   return "PdfRT";
}

//______________________________________________________________________________
const char *PdfRT::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::PdfRT*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int PdfRT::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::PdfRT*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *PdfRT::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::PdfRT*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *PdfRT::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::PdfRT*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr PdfWT::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *PdfWT::Class_Name()
{
   return "PdfWT";
}

//______________________________________________________________________________
const char *PdfWT::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::PdfWT*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int PdfWT::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::PdfWT*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *PdfWT::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::PdfWT*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *PdfWT::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::PdfWT*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr DecayRate::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *DecayRate::Class_Name()
{
   return "DecayRate";
}

//______________________________________________________________________________
const char *DecayRate::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::DecayRate*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int DecayRate::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::DecayRate*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *DecayRate::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::DecayRate*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *DecayRate::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::DecayRate*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr PdfSigAng::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *PdfSigAng::Class_Name()
{
   return "PdfSigAng";
}

//______________________________________________________________________________
const char *PdfSigAng::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::PdfSigAng*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int PdfSigAng::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::PdfSigAng*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *PdfSigAng::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::PdfSigAng*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *PdfSigAng::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::PdfSigAng*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
void PdfRT::Streamer(TBuffer &R__b)
{
   // Stream an object of class PdfRT.

   UInt_t R__s, R__c;
   if (R__b.IsReading()) {
      Version_t R__v = R__b.ReadVersion(&R__s, &R__c); if (R__v) { }
      RooAbsPdf::Streamer(R__b);
      ctK.Streamer(R__b);
      ctL.Streamer(R__b);
      phi.Streamer(R__b);
      Fl.Streamer(R__b);
      P1.Streamer(R__b);
      P2.Streamer(R__b);
      P3.Streamer(R__b);
      P4p.Streamer(R__b);
      P5p.Streamer(R__b);
      P6p.Streamer(R__b);
      P8p.Streamer(R__b);
      Eff.Streamer(R__b);
      {
         vector<double> &R__stl =  intPart;
         R__stl.clear();
         int R__i, R__n;
         R__b >> R__n;
         R__stl.reserve(R__n);
         for (R__i = 0; R__i < R__n; R__i++) {
            double R__t;
            R__b >> R__t;
            R__stl.push_back(R__t);
         }
      }
      R__b.CheckByteCount(R__s, R__c, PdfRT::IsA());
   } else {
      R__c = R__b.WriteVersion(PdfRT::IsA(), kTRUE);
      RooAbsPdf::Streamer(R__b);
      ctK.Streamer(R__b);
      ctL.Streamer(R__b);
      phi.Streamer(R__b);
      Fl.Streamer(R__b);
      P1.Streamer(R__b);
      P2.Streamer(R__b);
      P3.Streamer(R__b);
      P4p.Streamer(R__b);
      P5p.Streamer(R__b);
      P6p.Streamer(R__b);
      P8p.Streamer(R__b);
      Eff.Streamer(R__b);
      {
         vector<double> &R__stl =  intPart;
         int R__n=int(R__stl.size());
         R__b << R__n;
         if(R__n) {
            vector<double>::iterator R__k;
            for (R__k = R__stl.begin(); R__k != R__stl.end(); ++R__k) {
            R__b << (*R__k);
            }
         }
      }
      R__b.SetByteCount(R__c, kTRUE);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_PdfRT(void *p) {
      return  p ? new(p) ::PdfRT : new ::PdfRT;
   }
   static void *newArray_PdfRT(Long_t nElements, void *p) {
      return p ? new(p) ::PdfRT[nElements] : new ::PdfRT[nElements];
   }
   // Wrapper around operator delete
   static void delete_PdfRT(void *p) {
      delete ((::PdfRT*)p);
   }
   static void deleteArray_PdfRT(void *p) {
      delete [] ((::PdfRT*)p);
   }
   static void destruct_PdfRT(void *p) {
      typedef ::PdfRT current_t;
      ((current_t*)p)->~current_t();
   }
   // Wrapper around a custom streamer member function.
   static void streamer_PdfRT(TBuffer &buf, void *obj) {
      ((::PdfRT*)obj)->::PdfRT::Streamer(buf);
   }
} // end of namespace ROOT for class ::PdfRT

//______________________________________________________________________________
void PdfWT::Streamer(TBuffer &R__b)
{
   // Stream an object of class PdfWT.

   UInt_t R__s, R__c;
   if (R__b.IsReading()) {
      Version_t R__v = R__b.ReadVersion(&R__s, &R__c); if (R__v) { }
      RooAbsPdf::Streamer(R__b);
      ctK.Streamer(R__b);
      ctL.Streamer(R__b);
      phi.Streamer(R__b);
      Fl.Streamer(R__b);
      P1.Streamer(R__b);
      P2.Streamer(R__b);
      P3.Streamer(R__b);
      P4p.Streamer(R__b);
      P5p.Streamer(R__b);
      P6p.Streamer(R__b);
      P8p.Streamer(R__b);
      Eff.Streamer(R__b);
      {
         vector<double> &R__stl =  intPart;
         R__stl.clear();
         int R__i, R__n;
         R__b >> R__n;
         R__stl.reserve(R__n);
         for (R__i = 0; R__i < R__n; R__i++) {
            double R__t;
            R__b >> R__t;
            R__stl.push_back(R__t);
         }
      }
      R__b.CheckByteCount(R__s, R__c, PdfWT::IsA());
   } else {
      R__c = R__b.WriteVersion(PdfWT::IsA(), kTRUE);
      RooAbsPdf::Streamer(R__b);
      ctK.Streamer(R__b);
      ctL.Streamer(R__b);
      phi.Streamer(R__b);
      Fl.Streamer(R__b);
      P1.Streamer(R__b);
      P2.Streamer(R__b);
      P3.Streamer(R__b);
      P4p.Streamer(R__b);
      P5p.Streamer(R__b);
      P6p.Streamer(R__b);
      P8p.Streamer(R__b);
      Eff.Streamer(R__b);
      {
         vector<double> &R__stl =  intPart;
         int R__n=int(R__stl.size());
         R__b << R__n;
         if(R__n) {
            vector<double>::iterator R__k;
            for (R__k = R__stl.begin(); R__k != R__stl.end(); ++R__k) {
            R__b << (*R__k);
            }
         }
      }
      R__b.SetByteCount(R__c, kTRUE);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_PdfWT(void *p) {
      return  p ? new(p) ::PdfWT : new ::PdfWT;
   }
   static void *newArray_PdfWT(Long_t nElements, void *p) {
      return p ? new(p) ::PdfWT[nElements] : new ::PdfWT[nElements];
   }
   // Wrapper around operator delete
   static void delete_PdfWT(void *p) {
      delete ((::PdfWT*)p);
   }
   static void deleteArray_PdfWT(void *p) {
      delete [] ((::PdfWT*)p);
   }
   static void destruct_PdfWT(void *p) {
      typedef ::PdfWT current_t;
      ((current_t*)p)->~current_t();
   }
   // Wrapper around a custom streamer member function.
   static void streamer_PdfWT(TBuffer &buf, void *obj) {
      ((::PdfWT*)obj)->::PdfWT::Streamer(buf);
   }
} // end of namespace ROOT for class ::PdfWT

//______________________________________________________________________________
void DecayRate::Streamer(TBuffer &R__b)
{
   // Stream an object of class DecayRate.

   UInt_t R__s, R__c;
   if (R__b.IsReading()) {
      Version_t R__v = R__b.ReadVersion(&R__s, &R__c); if (R__v) { }
      RooAbsPdf::Streamer(R__b);
      ctK.Streamer(R__b);
      ctL.Streamer(R__b);
      phi.Streamer(R__b);
      Fl.Streamer(R__b);
      P1.Streamer(R__b);
      P2.Streamer(R__b);
      P3.Streamer(R__b);
      P4p.Streamer(R__b);
      P5p.Streamer(R__b);
      P6p.Streamer(R__b);
      P8p.Streamer(R__b);
      R__b.CheckByteCount(R__s, R__c, DecayRate::IsA());
   } else {
      R__c = R__b.WriteVersion(DecayRate::IsA(), kTRUE);
      RooAbsPdf::Streamer(R__b);
      ctK.Streamer(R__b);
      ctL.Streamer(R__b);
      phi.Streamer(R__b);
      Fl.Streamer(R__b);
      P1.Streamer(R__b);
      P2.Streamer(R__b);
      P3.Streamer(R__b);
      P4p.Streamer(R__b);
      P5p.Streamer(R__b);
      P6p.Streamer(R__b);
      P8p.Streamer(R__b);
      R__b.SetByteCount(R__c, kTRUE);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_DecayRate(void *p) {
      return  p ? new(p) ::DecayRate : new ::DecayRate;
   }
   static void *newArray_DecayRate(Long_t nElements, void *p) {
      return p ? new(p) ::DecayRate[nElements] : new ::DecayRate[nElements];
   }
   // Wrapper around operator delete
   static void delete_DecayRate(void *p) {
      delete ((::DecayRate*)p);
   }
   static void deleteArray_DecayRate(void *p) {
      delete [] ((::DecayRate*)p);
   }
   static void destruct_DecayRate(void *p) {
      typedef ::DecayRate current_t;
      ((current_t*)p)->~current_t();
   }
   // Wrapper around a custom streamer member function.
   static void streamer_DecayRate(TBuffer &buf, void *obj) {
      ((::DecayRate*)obj)->::DecayRate::Streamer(buf);
   }
} // end of namespace ROOT for class ::DecayRate

//______________________________________________________________________________
void PdfSigAng::Streamer(TBuffer &R__b)
{
   // Stream an object of class PdfSigAng.

   UInt_t R__s, R__c;
   if (R__b.IsReading()) {
      Version_t R__v = R__b.ReadVersion(&R__s, &R__c); if (R__v) { }
      RooAbsPdf::Streamer(R__b);
      ctK.Streamer(R__b);
      ctL.Streamer(R__b);
      phi.Streamer(R__b);
      Fl.Streamer(R__b);
      P1.Streamer(R__b);
      P2.Streamer(R__b);
      P3.Streamer(R__b);
      P4p.Streamer(R__b);
      P5p.Streamer(R__b);
      P6p.Streamer(R__b);
      P8p.Streamer(R__b);
      mFrac.Streamer(R__b);
      EffC.Streamer(R__b);
      EffW.Streamer(R__b);
      {
         vector<double> &R__stl =  intCPart;
         R__stl.clear();
         int R__i, R__n;
         R__b >> R__n;
         R__stl.reserve(R__n);
         for (R__i = 0; R__i < R__n; R__i++) {
            double R__t;
            R__b >> R__t;
            R__stl.push_back(R__t);
         }
      }
      {
         vector<double> &R__stl =  intWPart;
         R__stl.clear();
         int R__i, R__n;
         R__b >> R__n;
         R__stl.reserve(R__n);
         for (R__i = 0; R__i < R__n; R__i++) {
            double R__t;
            R__b >> R__t;
            R__stl.push_back(R__t);
         }
      }
      R__b.CheckByteCount(R__s, R__c, PdfSigAng::IsA());
   } else {
      R__c = R__b.WriteVersion(PdfSigAng::IsA(), kTRUE);
      RooAbsPdf::Streamer(R__b);
      ctK.Streamer(R__b);
      ctL.Streamer(R__b);
      phi.Streamer(R__b);
      Fl.Streamer(R__b);
      P1.Streamer(R__b);
      P2.Streamer(R__b);
      P3.Streamer(R__b);
      P4p.Streamer(R__b);
      P5p.Streamer(R__b);
      P6p.Streamer(R__b);
      P8p.Streamer(R__b);
      mFrac.Streamer(R__b);
      EffC.Streamer(R__b);
      EffW.Streamer(R__b);
      {
         vector<double> &R__stl =  intCPart;
         int R__n=int(R__stl.size());
         R__b << R__n;
         if(R__n) {
            vector<double>::iterator R__k;
            for (R__k = R__stl.begin(); R__k != R__stl.end(); ++R__k) {
            R__b << (*R__k);
            }
         }
      }
      {
         vector<double> &R__stl =  intWPart;
         int R__n=int(R__stl.size());
         R__b << R__n;
         if(R__n) {
            vector<double>::iterator R__k;
            for (R__k = R__stl.begin(); R__k != R__stl.end(); ++R__k) {
            R__b << (*R__k);
            }
         }
      }
      R__b.SetByteCount(R__c, kTRUE);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_PdfSigAng(void *p) {
      return  p ? new(p) ::PdfSigAng : new ::PdfSigAng;
   }
   static void *newArray_PdfSigAng(Long_t nElements, void *p) {
      return p ? new(p) ::PdfSigAng[nElements] : new ::PdfSigAng[nElements];
   }
   // Wrapper around operator delete
   static void delete_PdfSigAng(void *p) {
      delete ((::PdfSigAng*)p);
   }
   static void deleteArray_PdfSigAng(void *p) {
      delete [] ((::PdfSigAng*)p);
   }
   static void destruct_PdfSigAng(void *p) {
      typedef ::PdfSigAng current_t;
      ((current_t*)p)->~current_t();
   }
   // Wrapper around a custom streamer member function.
   static void streamer_PdfSigAng(TBuffer &buf, void *obj) {
      ((::PdfSigAng*)obj)->::PdfSigAng::Streamer(buf);
   }
} // end of namespace ROOT for class ::PdfSigAng

namespace ROOT {
   static TClass *vectorlEdoublegR_Dictionary();
   static void vectorlEdoublegR_TClassManip(TClass*);
   static void *new_vectorlEdoublegR(void *p = 0);
   static void *newArray_vectorlEdoublegR(Long_t size, void *p);
   static void delete_vectorlEdoublegR(void *p);
   static void deleteArray_vectorlEdoublegR(void *p);
   static void destruct_vectorlEdoublegR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<double>*)
   {
      vector<double> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<double>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<double>", -2, "vector", 214,
                  typeid(vector<double>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEdoublegR_Dictionary, isa_proxy, 0,
                  sizeof(vector<double>) );
      instance.SetNew(&new_vectorlEdoublegR);
      instance.SetNewArray(&newArray_vectorlEdoublegR);
      instance.SetDelete(&delete_vectorlEdoublegR);
      instance.SetDeleteArray(&deleteArray_vectorlEdoublegR);
      instance.SetDestructor(&destruct_vectorlEdoublegR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<double> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<double>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEdoublegR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<double>*)0x0)->GetClass();
      vectorlEdoublegR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEdoublegR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEdoublegR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<double> : new vector<double>;
   }
   static void *newArray_vectorlEdoublegR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<double>[nElements] : new vector<double>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEdoublegR(void *p) {
      delete ((vector<double>*)p);
   }
   static void deleteArray_vectorlEdoublegR(void *p) {
      delete [] ((vector<double>*)p);
   }
   static void destruct_vectorlEdoublegR(void *p) {
      typedef vector<double> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<double>

namespace {
  void TriggerDictionaryInitialization_AngDict_Impl() {
    static const char* headers[] = {
"interface/PdfRT.h",
"interface/PdfWT.h",
"interface/DecayRate.h",
"interface/PdfSigAng.h",
0
    };
    static const char* includePaths[] = {
"/usr/include/root",
"/afs/cern.ch/user/l/llinwei/private/eff-KDE-working-fitWorkflow/",
0
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "AngDict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_Autoloading_Map;
class __attribute__((annotate(R"ATTRDUMP(PDF for (angular decay rate x efficiency) description)ATTRDUMP"))) __attribute__((annotate("$clingAutoload$interface/PdfRT.h")))  PdfRT;
class __attribute__((annotate(R"ATTRDUMP(PDF for (angular decay rate x efficiency) description)ATTRDUMP"))) __attribute__((annotate("$clingAutoload$interface/PdfWT.h")))  PdfWT;
class __attribute__((annotate(R"ATTRDUMP(PDF for angular decay rate description)ATTRDUMP"))) __attribute__((annotate("$clingAutoload$interface/DecayRate.h")))  DecayRate;
class __attribute__((annotate(R"ATTRDUMP(PDF for (angular decay rate x efficiency) of both correctly-tagged and wrongly-tagged events)ATTRDUMP"))) __attribute__((annotate("$clingAutoload$interface/PdfSigAng.h")))  PdfSigAng;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "AngDict dictionary payload"

#ifndef G__VECTOR_HAS_CLASS_ITERATOR
  #define G__VECTOR_HAS_CLASS_ITERATOR 1
#endif

#define _BACKWARD_BACKWARD_WARNING_H
#include "interface/PdfRT.h"
#include "interface/PdfWT.h"
#include "interface/DecayRate.h"
#include "interface/PdfSigAng.h"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[]={
"DecayRate", payloadCode, "@",
"PdfRT", payloadCode, "@",
"PdfSigAng", payloadCode, "@",
"PdfWT", payloadCode, "@",
nullptr};

    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("AngDict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_AngDict_Impl, {}, classesHeaders, /*has no C++ module*/false);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_AngDict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_AngDict() {
  TriggerDictionaryInitialization_AngDict_Impl();
}
