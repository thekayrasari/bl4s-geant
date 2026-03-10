//
// ********************************************************************
// * License and Disclaimer                                           *
// ********************************************************************
//
/// \file RunAction.cc
/// \brief Implementation of the B1::RunAction class

#include "RunAction.hh"

#include "DetectorConstruction.hh"
#include "PrimaryGeneratorAction.hh"

#include "G4AccumulableManager.hh"
#include "G4LogicalVolume.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleGun.hh"
#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"

#include <cmath>

namespace B1
{

RunAction::RunAction()
{
  const G4double milligray = 1.e-3 * gray;
  const G4double microgray = 1.e-6 * gray;
  const G4double nanogray  = 1.e-9 * gray;
  const G4double picogray  = 1.e-12 * gray;

  new G4UnitDefinition("milligray", "milliGy", "Dose", milligray);
  new G4UnitDefinition("microgray", "microGy", "Dose", microgray);
  new G4UnitDefinition("nanogray", "nanoGy", "Dose", nanogray);
  new G4UnitDefinition("picogray", "picoGy", "Dose", picogray);

  G4AccumulableManager* accumulableManager = G4AccumulableManager::Instance();
  accumulableManager->Register(fEdep);
  accumulableManager->Register(fEdep2);
}

RunAction::~RunAction()
{
  if (fVolumeEntriesFile.is_open()) fVolumeEntriesFile.close();
  if (fEdepFile.is_open()) fEdepFile.close();
  if (fDetectorSpectrumFile.is_open()) fDetectorSpectrumFile.close();
  if (fPhotonSpectrumFile.is_open()) fPhotonSpectrumFile.close();
  if (fPhantomDoseFile.is_open()) fPhantomDoseFile.close();
  if (fAngularFile.is_open()) fAngularFile.close();
}

void RunAction::BeginOfRunAction(const G4Run*)
{
  G4RunManager::GetRunManager()->SetRandomNumberStore(false);

  G4AccumulableManager* accumulableManager = G4AccumulableManager::Instance();
  accumulableManager->Reset();

  if (IsMaster()) return;

  fVolumeEntriesFile.open("volume_entries.csv", std::ios::out);
  fEdepFile.open("volume_edep.csv", std::ios::out);
  fDetectorSpectrumFile.open("detector_spectrum.csv", std::ios::out);
  fPhotonSpectrumFile.open("photon_spectrum.csv", std::ios::out);
  fPhantomDoseFile.open("phantom_dose_distribution.csv", std::ios::out);
  fAngularFile.open("angular_distribution.csv", std::ios::out);

  if (fVolumeEntriesFile.is_open()) {
    fVolumeEntriesFile
      << "eventID,volume,particle,energy_MeV,x_mm,y_mm,z_mm,process,fromVolume\n";
  }

  if (fEdepFile.is_open()) {
    fEdepFile
      << "eventID,volume,particle,edep_MeV\n";
  }

  if (fDetectorSpectrumFile.is_open()) {
    fDetectorSpectrumFile
      << "eventID,volume,particle,energy_MeV,x_mm,y_mm,z_mm\n";
  }

  if (fPhotonSpectrumFile.is_open()) {
    fPhotonSpectrumFile
      << "eventID,volume,energy_MeV,x_mm,y_mm,z_mm\n";
  }

  if (fPhantomDoseFile.is_open()) {
    fPhantomDoseFile
      << "eventID,particle,x_mm,y_mm,z_mm,edep_MeV\n";
  }

  if (fAngularFile.is_open()) {
    fAngularFile
      << "eventID,volume,particle,energy_MeV,theta_deg,phi_deg\n";
  }
}

void RunAction::EndOfRunAction(const G4Run* run)
{
  G4int nofEvents = run->GetNumberOfEvent();
  if (nofEvents == 0) return;

  G4AccumulableManager* accumulableManager = G4AccumulableManager::Instance();
  accumulableManager->Merge();

  G4double edep  = fEdep.GetValue();
  G4double edep2 = fEdep2.GetValue();

  G4double rms = edep2 - edep * edep / nofEvents;
  if (rms > 0.) rms = std::sqrt(rms);
  else rms = 0.;

  const auto detConstruction = static_cast<const DetectorConstruction*>(
    G4RunManager::GetRunManager()->GetUserDetectorConstruction());
  G4double mass = detConstruction->GetScoringVolume()->GetMass();
  G4double dose = edep / mass;
  G4double rmsDose = rms / mass;

  const auto generatorAction = static_cast<const PrimaryGeneratorAction*>(
    G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction());
  G4String runCondition;
  if (generatorAction) {
    const G4ParticleGun* particleGun = generatorAction->GetParticleGun();
    runCondition += particleGun->GetParticleDefinition()->GetParticleName();
    runCondition += " of ";
    G4double particleEnergy = particleGun->GetParticleEnergy();
    runCondition += G4BestUnit(particleEnergy, "Energy");
  }

  if (IsMaster()) {
    G4cout << G4endl << "--------------------End of Global Run-----------------------";
  } else {
    G4cout << G4endl << "--------------------End of Local Run------------------------";
  }

  G4cout << G4endl << " The run is " << nofEvents << " " << runCondition << G4endl << G4endl;
  G4cout << "  --> cumulated edep per run in scoring volume = "
         << G4BestUnit(edep, "Energy")
         << " = " << edep / joule << " joule" << G4endl;
  G4cout << "  --> mass of scoring volume = " << G4BestUnit(mass, "Mass") << G4endl
         << G4endl;
  G4cout << " Absorbed dose per run in scoring volume = edep/mass = "
         << G4BestUnit(dose, "Dose")
         << "; rms = " << G4BestUnit(rmsDose, "Dose") << G4endl
         << "------------------------------------------------------------" << G4endl << G4endl;

  if (!IsMaster()) {
    if (fVolumeEntriesFile.is_open()) fVolumeEntriesFile.close();
    if (fEdepFile.is_open()) fEdepFile.close();
    if (fDetectorSpectrumFile.is_open()) fDetectorSpectrumFile.close();
    if (fPhotonSpectrumFile.is_open()) fPhotonSpectrumFile.close();
    if (fPhantomDoseFile.is_open()) fPhantomDoseFile.close();
    if (fAngularFile.is_open()) fAngularFile.close();
  }
}

void RunAction::AddEdep(G4double edep)
{
  fEdep += edep;
  fEdep2 += edep * edep;
}

void RunAction::WriteVolumeEntry(const VolumeEntryRecord& rec)
{
  if (IsMaster()) return;
  if (!fVolumeEntriesFile.is_open()) return;

  fVolumeEntriesFile
    << rec.eventID << ","
    << rec.volume << ","
    << rec.particle << ","
    << rec.energyMeV << ","
    << rec.xmm << ","
    << rec.ymm << ","
    << rec.zmm << ","
    << rec.process << ","
    << rec.fromVolume << "\n";
}

void RunAction::WriteEdepRecord(const EdepRecord& rec)
{
  if (IsMaster()) return;
  if (!fEdepFile.is_open()) return;

  fEdepFile
    << rec.eventID << ","
    << rec.volume << ","
    << rec.particle << ","
    << rec.edepMeV << "\n";
}

void RunAction::WriteDetectorSpectrumRecord(const DetectorSpectrumRecord& rec)
{
  if (IsMaster()) return;
  if (!fDetectorSpectrumFile.is_open()) return;

  fDetectorSpectrumFile
    << rec.eventID << ","
    << rec.volume << ","
    << rec.particle << ","
    << rec.energyMeV << ","
    << rec.xmm << ","
    << rec.ymm << ","
    << rec.zmm << "\n";

  if (rec.particle == "gamma" && fPhotonSpectrumFile.is_open()) {
    fPhotonSpectrumFile
      << rec.eventID << ","
      << rec.volume << ","
      << rec.energyMeV << ","
      << rec.xmm << ","
      << rec.ymm << ","
      << rec.zmm << "\n";
  }
}

void RunAction::WritePhantomDoseRecord(const PhantomDoseRecord& rec)
{
  if (IsMaster()) return;
  if (!fPhantomDoseFile.is_open()) return;

  fPhantomDoseFile
    << rec.eventID << ","
    << rec.particle << ","
    << rec.xmm << ","
    << rec.ymm << ","
    << rec.zmm << ","
    << rec.edepMeV << "\n";
}

void RunAction::WriteAngularRecord(const AngularRecord& rec)
{
  if (IsMaster()) return;
  if (!fAngularFile.is_open()) return;

  fAngularFile
    << rec.eventID << ","
    << rec.volume << ","
    << rec.particle << ","
    << rec.energyMeV << ","
    << rec.thetaDeg << ","
    << rec.phiDeg << "\n";
}

}  // namespace B1