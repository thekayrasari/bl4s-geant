//
// ********************************************************************
// * License and Disclaimer                                           *
// ********************************************************************
//
/// \file RunAction.hh
/// \brief Definition of the B1::RunAction class

#ifndef B1RunAction_h
#define B1RunAction_h 1

#include "G4UserRunAction.hh"
#include "G4Accumulable.hh"
#include "globals.hh"

#include "EventAction.hh"

#include <fstream>

class G4Run;

namespace B1
{

class RunAction : public G4UserRunAction
{
  public:
    RunAction();
    ~RunAction() override;

    void BeginOfRunAction(const G4Run*) override;
    void EndOfRunAction(const G4Run*) override;

    void AddEdep(G4double edep);

    void WriteVolumeEntry(const VolumeEntryRecord& rec);
    void WriteEdepRecord(const EdepRecord& rec);
    void WriteDetectorSpectrumRecord(const DetectorSpectrumRecord& rec);
    void WritePhantomDoseRecord(const PhantomDoseRecord& rec);
    void WriteAngularRecord(const AngularRecord& rec);

  private:
    G4Accumulable<G4double> fEdep = 0.;
    G4Accumulable<G4double> fEdep2 = 0.;

    std::ofstream fVolumeEntriesFile;
    std::ofstream fEdepFile;
    std::ofstream fDetectorSpectrumFile;
    std::ofstream fPhotonSpectrumFile;
    std::ofstream fPhantomDoseFile;
    std::ofstream fAngularFile;
};

}  // namespace B1

#endif