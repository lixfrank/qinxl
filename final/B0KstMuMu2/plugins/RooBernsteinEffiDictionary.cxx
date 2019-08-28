// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME RooBernsteinEffiDictionary

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
#include "RooBernsteinEffi.h"

// Header files passed via #pragma extra_include

namespace ROOT {
   static void *new_RooBernsteinEffi(void *p = 0);
   static void *newArray_RooBernsteinEffi(Long_t size, void *p);
   static void delete_RooBernsteinEffi(void *p);
   static void deleteArray_RooBernsteinEffi(void *p);
   static void destruct_RooBernsteinEffi(void *p);
   static void streamer_RooBernsteinEffi(TBuffer &buf, void *obj);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::RooBernsteinEffi*)
   {
      ::RooBernsteinEffi *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::RooBernsteinEffi >(0);
      static ::ROOT::TGenericClassInfo 
         instance("RooBernsteinEffi", ::RooBernsteinEffi::Class_Version(), "RooBernsteinEffi.h", 20,
                  typeid(::RooBernsteinEffi), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::RooBernsteinEffi::Dictionary, isa_proxy, 16,
                  sizeof(::RooBernsteinEffi) );
      instance.SetNew(&new_RooBernsteinEffi);
      instance.SetNewArray(&newArray_RooBernsteinEffi);
      instance.SetDelete(&delete_RooBernsteinEffi);
      instance.SetDeleteArray(&deleteArray_RooBernsteinEffi);
      instance.SetDestructor(&destruct_RooBernsteinEffi);
      instance.SetStreamerFunc(&streamer_RooBernsteinEffi);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::RooBernsteinEffi*)
   {
      return GenerateInitInstanceLocal((::RooBernsteinEffi*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::RooBernsteinEffi*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

//______________________________________________________________________________
atomic_TClass_ptr RooBernsteinEffi::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *RooBernsteinEffi::Class_Name()
{
   return "RooBernsteinEffi";
}

//______________________________________________________________________________
const char *RooBernsteinEffi::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::RooBernsteinEffi*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int RooBernsteinEffi::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::RooBernsteinEffi*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *RooBernsteinEffi::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::RooBernsteinEffi*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *RooBernsteinEffi::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::RooBernsteinEffi*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
void RooBernsteinEffi::Streamer(TBuffer &R__b)
{
   // Stream an object of class RooBernsteinEffi.

   UInt_t R__s, R__c;
   if (R__b.IsReading()) {
      Version_t R__v = R__b.ReadVersion(&R__s, &R__c); if (R__v) { }
      RooAbsReal::Streamer(R__b);
      _x.Streamer(R__b);
      _y.Streamer(R__b);
      _z.Streamer(R__b);
      _coefList.Streamer(R__b);
      R__b >> _maxDegree1;
      R__b >> _maxDegree2;
      R__b >> _maxDegree3;
      R__b >> _maxDegreeV;
      R__b.CheckByteCount(R__s, R__c, RooBernsteinEffi::IsA());
   } else {
      R__c = R__b.WriteVersion(RooBernsteinEffi::IsA(), kTRUE);
      RooAbsReal::Streamer(R__b);
      _x.Streamer(R__b);
      _y.Streamer(R__b);
      _z.Streamer(R__b);
      _coefList.Streamer(R__b);
      R__b << _maxDegree1;
      R__b << _maxDegree2;
      R__b << _maxDegree3;
      R__b << _maxDegreeV;
      R__b.SetByteCount(R__c, kTRUE);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_RooBernsteinEffi(void *p) {
      return  p ? new(p) ::RooBernsteinEffi : new ::RooBernsteinEffi;
   }
   static void *newArray_RooBernsteinEffi(Long_t nElements, void *p) {
      return p ? new(p) ::RooBernsteinEffi[nElements] : new ::RooBernsteinEffi[nElements];
   }
   // Wrapper around operator delete
   static void delete_RooBernsteinEffi(void *p) {
      delete ((::RooBernsteinEffi*)p);
   }
   static void deleteArray_RooBernsteinEffi(void *p) {
      delete [] ((::RooBernsteinEffi*)p);
   }
   static void destruct_RooBernsteinEffi(void *p) {
      typedef ::RooBernsteinEffi current_t;
      ((current_t*)p)->~current_t();
   }
   // Wrapper around a custom streamer member function.
   static void streamer_RooBernsteinEffi(TBuffer &buf, void *obj) {
      ((::RooBernsteinEffi*)obj)->::RooBernsteinEffi::Streamer(buf);
   }
} // end of namespace ROOT for class ::RooBernsteinEffi

namespace {
  void TriggerDictionaryInitialization_RooBernsteinEffiDictionary_Impl() {
    static const char* headers[] = {
"RooBernsteinEffi.h",
0
    };
    static const char* includePaths[] = {
"/usr/include/root",
"/afs/cern.ch/user/l/llinwei/private/B0KstMuMufull/B0KstMuMu2/plugins/",
0
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "RooBernsteinEffiDictionary dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_Autoloading_Map;
class __attribute__((annotate(R"ATTRDUMP(Your description goes here...)ATTRDUMP"))) __attribute__((annotate("$clingAutoload$RooBernsteinEffi.h")))  RooBernsteinEffi;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "RooBernsteinEffiDictionary dictionary payload"

#ifndef G__VECTOR_HAS_CLASS_ITERATOR
  #define G__VECTOR_HAS_CLASS_ITERATOR 1
#endif

#define _BACKWARD_BACKWARD_WARNING_H
#include "RooBernsteinEffi.h"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[]={
"RooBernsteinEffi", payloadCode, "@",
nullptr};

    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("RooBernsteinEffiDictionary",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_RooBernsteinEffiDictionary_Impl, {}, classesHeaders, /*has no C++ module*/false);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_RooBernsteinEffiDictionary_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_RooBernsteinEffiDictionary() {
  TriggerDictionaryInitialization_RooBernsteinEffiDictionary_Impl();
}
