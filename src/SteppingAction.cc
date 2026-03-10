//
// ********************************************************************
// * License and Disclaimer                                           *
// ********************************************************************
//
/// \file SteppingAction.cc
/// \brief Implementation of the B1::SteppingAction class

#include "SteppingAction.hh"

#include "DetectorConstruction.hh"
#include "EventAction.hh"

#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4LogicalVolume.hh"
#include "G4RunManager.hh"
#include "G4Step.hh"
#include "G4StepPoint.hh"
#include "G4SystemOfUnits.hh"
#include "G4Track.hh"
#include "G4VPhysicalVolume.hh"
#include "G4VProcess.hh"
#include "G4ThreeVector.hh"
#include "G4ios.hh"

namespace B1
{

SteppingAction::SteppingAction(EventAction* eventAction)
  : fEventAction(eventAction)
{}

void SteppingAction::UserSteppingAction(const G4Step* step)
{
  if (!fScoringVolume) {
    const auto detConstruction = static_cast<const DetectorConstruction*>(
      G4RunManager::GetRunManager()->GetUserDetectorConstruction());
    fScoringVolume = detConstruction->GetScoringVolume();
  }

  const G4StepPoint* preStepPoint  = step->GetPreStepPoint();
  const G4StepPoint* postStepPoint = step->GetPostStepPoint();
  const G4Track* track             = step->GetTrack();

  if (!preStepPoint || !postStepPoint || !track) return;

  G4VPhysicalVolume* prePV  = preStepPoint->GetPhysicalVolume();
  G4VPhysicalVolume* postPV = postStepPoint->GetPhysicalVolume();

  if (!prePV) return;

  const G4String preVolName  = prePV->GetName();
  const G4String postVolName = (postPV) ? postPV->GetName() : "OutOfWorld";

  const G4String particleName = track->GetDefinition()->GetParticleName();
  const G4double kineticEnergyPre = preStepPoint->GetKineticEnergy();
  const G4double edepStep = step->GetTotalEnergyDeposit();

  const G4VProcess* process = postStepPoint->GetProcessDefinedStep();
  const G4String processName = (process) ? process->GetProcessName() : "NoProcess";

  const G4int eventID =
    G4EventManager::GetEventManager()->GetConstCurrentEvent()->GetEventID();

  // Original phantom scoring
  G4LogicalVolume* volume =
    preStepPoint->GetTouchableHandle()->GetVolume()->GetLogicalVolume();

  if (volume == fScoringVolume) {
    fEventAction->AddEdep(edepStep);
  }

  const G4bool isNewVolumeEntry = (preVolName != postVolName);

  // Detector / analysis volumes
  const G4bool isAnalysisVolume =
      (postVolName == "Cherenkov"   ||
       postVolName == "BremsTarget" ||
       postVolName == "Phantom"     ||
       postVolName == "Collimator"  ||
       postVolName == "Scint1"      ||
       postVolName == "DWC1"        ||
       postVolName == "DWC2"        ||
       postVolName == "Scint2");

  if (isNewVolumeEntry && isAnalysisVolume) {
    VolumeEntryRecord rec;
    rec.eventID    = eventID;
    rec.volume     = postVolName;
    rec.particle   = particleName;
    rec.energyMeV  = kineticEnergyPre / MeV;
    rec.xmm        = postStepPoint->GetPosition().x() / mm;
    rec.ymm        = postStepPoint->GetPosition().y() / mm;
    rec.zmm        = postStepPoint->GetPosition().z() / mm;
    rec.process    = processName;
    rec.fromVolume = preVolName;
    fEventAction->AddVolumeEntry(rec);

    DetectorSpectrumRecord srec;
    srec.eventID   = eventID;
    srec.volume    = postVolName;
    srec.particle  = particleName;
    srec.energyMeV = kineticEnergyPre / MeV;
    srec.xmm       = postStepPoint->GetPosition().x() / mm;
    srec.ymm       = postStepPoint->GetPosition().y() / mm;
    srec.zmm       = postStepPoint->GetPosition().z() / mm;
    fEventAction->AddDetectorSpectrumRecord(srec);

    if (particleName == "gamma") {
      const G4ThreeVector dir = postStepPoint->GetMomentumDirection();

      AngularRecord arec;
      arec.eventID   = eventID;
      arec.volume    = postVolName;
      arec.particle  = particleName;
      arec.energyMeV = kineticEnergyPre / MeV;
      arec.thetaDeg  = dir.theta() / deg;
      arec.phiDeg    = dir.phi() / deg;
      fEventAction->AddAngularRecord(arec);
    }
  }

  if (edepStep > 0.) {
    if (preVolName == "Phantom" ||
        preVolName == "Scint1"  ||
        preVolName == "Scint2"  ||
        preVolName == "BremsTarget") {

      EdepRecord rec;
      rec.eventID  = eventID;
      rec.volume   = preVolName;
      rec.particle = particleName;
      rec.edepMeV  = edepStep / MeV;
      fEventAction->AddEdepRecord(rec);
    }

    if (preVolName == "Phantom") {
      PhantomDoseRecord drec;
      drec.eventID  = eventID;
      drec.particle = particleName;
      drec.xmm      = preStepPoint->GetPosition().x() / mm;
      drec.ymm      = preStepPoint->GetPosition().y() / mm;
      drec.zmm      = preStepPoint->GetPosition().z() / mm;
      drec.edepMeV  = edepStep / MeV;
      fEventAction->AddPhantomDoseRecord(drec);
    }
  }

  // Optional terminal debug
  if (preVolName == "BremsTarget" && processName == "eBrem") {
    G4cout << "[Brems] event=" << eventID
           << " particle=" << particleName
           << " process=" << processName
           << " E=" << kineticEnergyPre / MeV << " MeV"
           << " volume=" << preVolName
           << G4endl;
  }

  if (preVolName == "Phantom" && processName == "compt") {
    G4cout << "[Compton] event=" << eventID
           << " particle=" << particleName
           << " process=" << processName
           << " E=" << kineticEnergyPre / MeV << " MeV"
           << " volume=" << preVolName
           << G4endl;
  }
}

}  // namespace B1