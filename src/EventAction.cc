//
// ********************************************************************
// * License and Disclaimer                                           *
// ********************************************************************
//
/// \file EventAction.cc
/// \brief Implementation of the B1::EventAction class

#include "EventAction.hh"

#include "RunAction.hh"

#include "G4Event.hh"

namespace B1
{

EventAction::EventAction(RunAction* runAction)
  : fRunAction(runAction)
{}

void EventAction::BeginOfEventAction(const G4Event*)
{
  fEdep = 0.;
  fVolumeEntries.clear();
  fEdepRecords.clear();
  fDetectorSpectrumRecords.clear();
  fPhantomDoseRecords.clear();
  fAngularRecords.clear();
}

void EventAction::EndOfEventAction(const G4Event*)
{
  fRunAction->AddEdep(fEdep);

  for (const auto& rec : fVolumeEntries) {
    fRunAction->WriteVolumeEntry(rec);
  }

  for (const auto& rec : fEdepRecords) {
    fRunAction->WriteEdepRecord(rec);
  }

  for (const auto& rec : fDetectorSpectrumRecords) {
    fRunAction->WriteDetectorSpectrumRecord(rec);
  }

  for (const auto& rec : fPhantomDoseRecords) {
    fRunAction->WritePhantomDoseRecord(rec);
  }

  for (const auto& rec : fAngularRecords) {
    fRunAction->WriteAngularRecord(rec);
  }
}

void EventAction::AddVolumeEntry(const VolumeEntryRecord& rec)
{
  fVolumeEntries.push_back(rec);
}

void EventAction::AddEdepRecord(const EdepRecord& rec)
{
  fEdepRecords.push_back(rec);
}

void EventAction::AddDetectorSpectrumRecord(const DetectorSpectrumRecord& rec)
{
  fDetectorSpectrumRecords.push_back(rec);
}

void EventAction::AddPhantomDoseRecord(const PhantomDoseRecord& rec)
{
  fPhantomDoseRecords.push_back(rec);
}

void EventAction::AddAngularRecord(const AngularRecord& rec)
{
  fAngularRecords.push_back(rec);
}

}  // namespace B1