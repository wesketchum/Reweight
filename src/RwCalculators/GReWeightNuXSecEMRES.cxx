//____________________________________________________________________________
/*
 Copyright (c) 2003-2018, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org

 Authors: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          University of Liverpool & STFC Rutherford Appleton Lab

          Jim Dobson <J.Dobson07 \at imperial.ac.uk>
          Imperial College London
*/
//____________________________________________________________________________

#include <cassert>

#include <TMath.h>

// GENIE/Generator includes
#include "Framework/Algorithm/AlgFactory.h"
#include "Framework/Algorithm/AlgConfigPool.h"
#include "Framework/EventGen/XSecAlgorithmI.h"
#include "Framework/Conventions/Units.h"
#include "Framework/Conventions/Controls.h"
#include "Framework/EventGen/EventRecord.h"
#include "Framework/GHEP/GHepParticle.h"
#include "Framework/Interaction/Interaction.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/Registry/Registry.h"

// GENIE/Reweight includes
#include "RwCalculators/GReWeightNuXSecEMRES.h"
#include "RwFramework/GSystSet.h"
#include "RwFramework/GSystUncertainty.h"

using namespace genie;
using namespace genie::rew;

const int GReWeightNuXSecEMRES::kModeMaMv;
const int GReWeightNuXSecEMRES::kModeNormAndMaMvShape;

//_______________________________________________________________________________________
GReWeightNuXSecEMRES::GReWeightNuXSecEMRES() :
GReWeightModel("EMRES"),
fManualModelName(),
fManualModelType()
{
  this->Init();
}
//_______________________________________________________________________________________
GReWeightNuXSecEMRES::GReWeightNuXSecEMRES(std::string model, std::string type) :
GReWeightModel("EMRES"),
fManualModelName(model),
fManualModelType(type)
{
  this->Init();
}
//_______________________________________________________________________________________
GReWeightNuXSecEMRES::~GReWeightNuXSecEMRES()
{

}
//_______________________________________________________________________________________
bool GReWeightNuXSecEMRES::IsHandled(GSyst_t syst) const
{
   bool handle;

   switch(syst) {

     case ( kXSecTwkDial_NormEMRES    ) :
     case ( kXSecTwkDial_MvEMRESshape ) :
       if(fMode==kModeNormAndMaMvShape) {
          handle = true;
       } else {
          handle = false;
       }
       break;

     case ( kXSecTwkDial_MvEMRES ) :
       if(fMode==kModeMaMv) {
          handle = true;
       } else {
          handle = false;
       }
       break;

     default:
          handle = false;
          break;
   }

   return handle;
}
//_______________________________________________________________________________________
bool GReWeightNuXSecEMRES::AppliesTo(ScatteringType_t type, bool is_cc) const
{
  //if (type==kScResonant && !is_cc) {
  //  return true;
  //}
  if (type==kScResonant) {
    return true;
  }
  return false;
}
//_______________________________________________________________________________________
void GReWeightNuXSecEMRES::SetSystematic(GSyst_t syst, double twk_dial)
{
  if(!this->IsHandled(syst)) return;

  switch(syst) {
    case ( kXSecTwkDial_NormEMRES ) :
      fNormTwkDial = twk_dial;
      break;
    case ( kXSecTwkDial_MvEMRESshape ) :
    case ( kXSecTwkDial_MvEMRES ) :
      fMvTwkDial = twk_dial;
      break;
    default:
      break;
  }
}
//_______________________________________________________________________________________
void GReWeightNuXSecEMRES::Reset(void)
{
  fNormTwkDial = 0.;
  fNormCurr    = fNormDef;
  fMaTwkDial   = 0.;
  fMaCurr      = fMaDef;
  fMvTwkDial   = 0.;
  fMvCurr      = fMvDef;

  this->Reconfigure();
}
//_______________________________________________________________________________________
void GReWeightNuXSecEMRES::Reconfigure(void)
{
  GSystUncertainty * fracerr = GSystUncertainty::Instance();

  if(fMode==kModeMaMv) {
     double fracerr_mv = fracerr->OneSigmaErr(kXSecTwkDial_MvEMRES);
     fMvCurr = fMvDef * (1. + fMvTwkDial * fracerr_mv);
///     std::cout<< "fMvCurr " << fMvCurr<<" fMvDef " << fMvDef << " fracerr_mv " <<fracerr_mv<< " fMvTwkDial "<< fMvTwkDial <<std::endl;
  }
  else
  if(fMode==kModeNormAndMaMvShape) {
     double fracerr_norm = fracerr->OneSigmaErr(kXSecTwkDial_NormEMRES);
     double fracerr_mvsh = fracerr->OneSigmaErr(kXSecTwkDial_MvEMRESshape);
     fNormCurr = fNormDef * (1. + fNormTwkDial * fracerr_norm);
     fMvCurr   = fMvDef   * (1. + fMvTwkDial   * fracerr_mvsh);
  }

  fNormCurr = TMath::Max(0., fNormCurr);
  fMvCurr   = TMath::Max(0., fMvCurr  );

  Registry & r = *fXSecModelConfig;

  r.Set(fMvPath, fMvCurr);
  fXSecModel->Configure(r);

//LOG("ReW, pDEBUG) << *fXSecModel;
}
//_______________________________________________________________________________________
double GReWeightNuXSecEMRES::CalcWeight(const genie::EventRecord & event)
{
  bool is_res = event.Summary()->ProcInfo().IsResonant();
  bool is_nc  = event.Summary()->ProcInfo().IsWeakNC();
  bool is_em  = event.Summary()->ProcInfo().IsEM();
  //if(!is_res || !is_nc) return 1.;
  if(!is_res || !is_em) return 1.;

  int nupdg = event.Probe()->Pdg();
  // Let's just try reweighting all probes
  //if(nupdg==kPdgNuMu     && !fRewNumu   ) return 1.;
  //if(nupdg==kPdgAntiNuMu && !fRewNumubar) return 1.;
  //if(nupdg==kPdgNuE      && !fRewNue    ) return 1.;
  //if(nupdg==kPdgAntiNuE  && !fRewNuebar ) return 1.;

  if(fMode==kModeMaMv) {
     double wght = this->CalcWeightMaMv(event);
     return wght;
  }
  else
  if(fMode==kModeNormAndMaMvShape) {
     double wght =
         this->CalcWeightNorm      (event) *
         this->CalcWeightMaMvShape (event);
     return wght;
  }

  return 1.;
}
//_______________________________________________________________________________________
void GReWeightNuXSecEMRES::Init(void)
{
  AlgConfigPool * conf_pool = AlgConfigPool::Instance();
  Registry * gpl = conf_pool->GlobalParameterList();
  RgAlg xsec_alg = gpl->GetAlg("XSecModel@genie::EventGenerator/RES-EM");

  AlgId id(xsec_alg);

  AlgId twk_id(id);
  if (fManualModelName.size()) {
    twk_id = AlgId(fManualModelName,fManualModelType);
  }

  AlgFactory * algf = AlgFactory::Instance();

  Algorithm * algdef = algf->AdoptAlgorithm(id);
  fXSecModelDef = dynamic_cast<XSecAlgorithmI*>(algdef);
  fXSecModelDef->AdoptSubstructure();

  Algorithm * alg = algf->AdoptAlgorithm(twk_id);
  fXSecModel = dynamic_cast<XSecAlgorithmI*>(alg);
  fXSecModel->AdoptSubstructure();

  fXSecModelConfig = new Registry(fXSecModel->GetConfig());
//LOG("ReW", pNOTICE) << *fXSecModelConfig;

  this->SetMode(kModeNormAndMaMvShape);

  this->RewNue    (true);
  this->RewNuebar (true);
  this->RewNumu   (true);
  this->RewNumubar(true);

  this->SetMaPath("RES-Ma");
  this->SetMvPath("RES-Mv");

  fNormTwkDial = 0.;
  fNormDef     = 1.;
  fNormCurr    = fNormDef;
  fMaTwkDial   = 0.;
  fMaDef       = fXSecModelConfig->GetDouble(fMaPath);
  fMaCurr      = fMaDef;
  fMvTwkDial   = 0.;
  fMvDef       = fXSecModelConfig->GetDouble(fMvPath);
  fMvCurr      = fMvDef;
}
//_______________________________________________________________________________________
double GReWeightNuXSecEMRES::CalcWeightNorm(const genie::EventRecord & /*event*/)
{
  bool tweaked = (TMath::Abs(fNormTwkDial) > controls::kASmallNum);
  if(!tweaked) return 1.0;

  double wght = fNormCurr;
  return wght;
}
//_______________________________________________________________________________________
double GReWeightNuXSecEMRES::CalcWeightMaMv(const genie::EventRecord & event)
{
  bool tweaked =
     (TMath::Abs(fMaTwkDial) > controls::kASmallNum) ||
     (TMath::Abs(fMvTwkDial) > controls::kASmallNum);
  if(!tweaked) return 1.0;

  Interaction * interaction = event.Summary();

  interaction->KinePtr()->UseSelectedKinematics();

  const KinePhaseSpace_t phase_space = kPSWQ2fE;

  double old_xsec   = event.DiffXSec();
  double calc_old_xsec = fXSecModelDef->XSec(interaction, phase_space);
  if (std::abs(calc_old_xsec - old_xsec)/old_xsec > controls::kASmallNum) {
        LOG("ReW",pWARN) << "Warning - default dxsec does not match dxsec saved in tree. Does the config match?";
      }
///  std::cout<< "old_xsec: " << old_xsec << " calc_old_xsec: " << calc_old_xsec <<std::endl;
  
//  double old_xsec   = event.DiffXSec();
//  if (!fUseOldWeightFromFile || fNWeightChecksDone < fNWeightChecksToDo) {
//    double calc_old_xsec = fXSecModelDef->XSec(interaction, phase_space);
//    if (fNWeightChecksDone < fNWeightChecksToDo) {
//      if (std::abs(calc_old_xsec - old_xsec)/old_xsec > controls::kASmallNum) {
//        LOG("ReW",pWARN) << "Warning - default dxsec does not match dxsec saved in tree. Does the config match?";
//      }
//      fNWeightChecksDone++;
//    }
//    if(!fUseOldWeightFromFile) {
//      old_xsec = calc_old_xsec;

  double old_weight = event.Weight();
  double new_xsec   = fXSecModel->XSec(interaction, phase_space);
  double new_weight = old_weight * (new_xsec/old_xsec);
//testing weights
///  std::cout << "old_weight" << "old_xsec" << "new_weight" << "new_xsec" <<std::endl;
///  std::cout << old_weight <<" "<<  old_xsec <<" "<< new_weight <<" "<< new_xsec <<std::endl;

  interaction->KinePtr()->ClearRunningValues();

  return new_weight;
}
//_______________________________________________________________________________________
double GReWeightNuXSecEMRES::CalcWeightMaMvShape(const genie::EventRecord & event)
{
  bool tweaked =
     (TMath::Abs(fMaTwkDial) > controls::kASmallNum) ||
     (TMath::Abs(fMvTwkDial) > controls::kASmallNum);
  if(!tweaked) return 1.0;

  Interaction * interaction = event.Summary();

  interaction->KinePtr()->UseSelectedKinematics();

  const KinePhaseSpace_t phase_space = kPSWQ2fE;

  double old_xsec   = event.DiffXSec();
  if (!fUseOldWeightFromFile || fNWeightChecksDone < fNWeightChecksToDo) {
    double calc_old_xsec = fXSecModelDef->XSec(interaction, phase_space);
    if (fNWeightChecksDone < fNWeightChecksToDo) {
      if (std::abs(calc_old_xsec - old_xsec)/old_xsec > controls::kASmallNum) {
        LOG("ReW",pWARN) << "Warning - default dxsec does not match dxsec saved in tree. Does the config match?";
      }
      fNWeightChecksDone++;
    }
    if(!fUseOldWeightFromFile) {
      old_xsec = calc_old_xsec;
    }
  }

  double old_weight = event.Weight();
  double new_xsec   = fXSecModel->XSec(interaction, phase_space);
  double new_weight = old_weight * (new_xsec/old_xsec);

//LOG("ReW", pDEBUG) << "differential cross section (old) = " << old_xsec;
//LOG("ReW", pDEBUG) << "differential cross section (new) = " << new_xsec;
//LOG("ReW", pDEBUG) << "event generation weight = " << old_weight;
//LOG("ReW", pDEBUG) << "new weight = " << new_weight;

//double old_integrated_xsec = event.XSec();
  double old_integrated_xsec = fXSecModelDef -> Integral(interaction);
  double twk_integrated_xsec = fXSecModel    -> Integral(interaction);
  assert(twk_integrated_xsec > 0);
  new_weight *= (old_integrated_xsec/twk_integrated_xsec);

//LOG("ReW", pDEBUG) << "integrated cross section (old) = " << old_integrated_xsec;
//LOG("ReW", pDEBUG) << "integrated cross section (twk) = " << twk_integrated_xsec;
//LOG("ReW", pDEBUG) << "new weight (normalized to const integral) = " << new_weight;

  interaction->KinePtr()->ClearRunningValues();

  return new_weight;
}
//_______________________________________________________________________________________
