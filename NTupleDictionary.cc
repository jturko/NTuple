// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME NTupleDictionary

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
#include "Griffin.hh"

// Header files passed via #pragma extra_include

namespace ROOT {
   static void *new_Detector(void *p = 0);
   static void *newArray_Detector(Long_t size, void *p);
   static void delete_Detector(void *p);
   static void deleteArray_Detector(void *p);
   static void destruct_Detector(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::Detector*)
   {
      ::Detector *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::Detector >(0);
      static ::ROOT::TGenericClassInfo 
         instance("Detector", ::Detector::Class_Version(), "Griffin.hh", 6,
                  typeid(::Detector), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::Detector::Dictionary, isa_proxy, 4,
                  sizeof(::Detector) );
      instance.SetNew(&new_Detector);
      instance.SetNewArray(&newArray_Detector);
      instance.SetDelete(&delete_Detector);
      instance.SetDeleteArray(&deleteArray_Detector);
      instance.SetDestructor(&destruct_Detector);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::Detector*)
   {
      return GenerateInitInstanceLocal((::Detector*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::Detector*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

//______________________________________________________________________________
atomic_TClass_ptr Detector::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *Detector::Class_Name()
{
   return "Detector";
}

//______________________________________________________________________________
const char *Detector::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::Detector*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int Detector::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::Detector*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *Detector::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::Detector*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *Detector::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::Detector*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
void Detector::Streamer(TBuffer &R__b)
{
   // Stream an object of class Detector.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(Detector::Class(),this);
   } else {
      R__b.WriteClassBuffer(Detector::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_Detector(void *p) {
      return  p ? new(p) ::Detector : new ::Detector;
   }
   static void *newArray_Detector(Long_t nElements, void *p) {
      return p ? new(p) ::Detector[nElements] : new ::Detector[nElements];
   }
   // Wrapper around operator delete
   static void delete_Detector(void *p) {
      delete ((::Detector*)p);
   }
   static void deleteArray_Detector(void *p) {
      delete [] ((::Detector*)p);
   }
   static void destruct_Detector(void *p) {
      typedef ::Detector current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::Detector

namespace ROOT {
   static TClass *vectorlEDetectorgR_Dictionary();
   static void vectorlEDetectorgR_TClassManip(TClass*);
   static void *new_vectorlEDetectorgR(void *p = 0);
   static void *newArray_vectorlEDetectorgR(Long_t size, void *p);
   static void delete_vectorlEDetectorgR(void *p);
   static void deleteArray_vectorlEDetectorgR(void *p);
   static void destruct_vectorlEDetectorgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<Detector>*)
   {
      vector<Detector> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<Detector>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<Detector>", -2, "vector", 447,
                  typeid(vector<Detector>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEDetectorgR_Dictionary, isa_proxy, 4,
                  sizeof(vector<Detector>) );
      instance.SetNew(&new_vectorlEDetectorgR);
      instance.SetNewArray(&newArray_vectorlEDetectorgR);
      instance.SetDelete(&delete_vectorlEDetectorgR);
      instance.SetDeleteArray(&deleteArray_vectorlEDetectorgR);
      instance.SetDestructor(&destruct_vectorlEDetectorgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<Detector> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<Detector>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEDetectorgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<Detector>*)0x0)->GetClass();
      vectorlEDetectorgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEDetectorgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEDetectorgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<Detector> : new vector<Detector>;
   }
   static void *newArray_vectorlEDetectorgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<Detector>[nElements] : new vector<Detector>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEDetectorgR(void *p) {
      delete ((vector<Detector>*)p);
   }
   static void deleteArray_vectorlEDetectorgR(void *p) {
      delete [] ((vector<Detector>*)p);
   }
   static void destruct_vectorlEDetectorgR(void *p) {
      typedef vector<Detector> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<Detector>

namespace {
  void TriggerDictionaryInitialization_NTupleDictionary_Impl() {
    static const char* headers[] = {
"Griffin.hh",
0
    };
    static const char* includePaths[] = {
"/Users/JTurko/programs/root_v6-14-00/build_me/include",
"/Users/JTurko/programs/CommandLineInterface",
"/Users/JTurko/programs/root_v6-14-00/build_me/include",
"/Users/JTurko/programs/NTuple_TISTAR/",
0
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "NTupleDictionary dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_Autoloading_Map;
class __attribute__((annotate("$clingAutoload$Griffin.hh")))  Detector;
namespace std{inline namespace __1{template <class _Tp> class __attribute__((annotate("$clingAutoload$iosfwd")))  __attribute__((annotate("$clingAutoload$string")))  allocator;
}}
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "NTupleDictionary dictionary payload"

#ifndef G__VECTOR_HAS_CLASS_ITERATOR
  #define G__VECTOR_HAS_CLASS_ITERATOR 1
#endif

#define _BACKWARD_BACKWARD_WARNING_H
#include "Griffin.hh"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[]={
"Detector", payloadCode, "@",
nullptr};

    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("NTupleDictionary",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_NTupleDictionary_Impl, {}, classesHeaders, /*has no C++ module*/false);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_NTupleDictionary_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_NTupleDictionary() {
  TriggerDictionaryInitialization_NTupleDictionary_Impl();
}
