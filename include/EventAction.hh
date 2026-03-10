//
// ********************************************************************
// * License and Disclaimer                                           *
// ********************************************************************
//
/// \file EventAction.hh
/// \brief Definition of the B1::EventAction class

#ifndef B1EventAction_h
#define B1EventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"

#include <vector>

class G4Event;

namespace B1
{

class RunAction;

struct VolumeEntryRecord
{
  G4int eventID = -1;
  G4String volume;
  G4String particle;
  G4double energyMeV = 0.0;
  G4double xmm = 0.0;
  G4double ymm = 0.0;
  G4double zmm = 0.0;
  G4String process;
  G4String fromVolume;
};

struct EdepRecord
{
  G4int eventID = -1;
  G4String volume;
  G4String particle;
  G4double edepMeV = 0.0;
};

struct DetectorSpectrumRecord
{
  G4int eventID = -1;
  G4String volume;
  G4String particle;
  G4double energyMeV = 0.0;
  G4double xmm = 0.0;
  G4double ymm = 0.0;
  G4double zmm = 0.0;
};

struct PhantomDoseRecord
{
  G4int eventID = -1;
  G4String particle;
  G4double xmm = 0.0;
  G4double ymm = 0.0;
  G4double zmm = 0.0;
  G4double edepMeV = 0.0;
};

struct AngularRecord
{
  G4int eventID = -1;
  G4String volume;
  G4String particle;
  G4double energyMeV = 0.0;
  G4double thetaDeg = 0.0;
  G4double phiDeg = 0.0;
};

class EventAction : public G4UserEventAction
{
  public:
    EventAction(RunAction* runAction);
    ~EventAction() override = default;

    void BeginOfEventAction(const G4Event* event) override;
    void EndOfEventAction(const G4Event* event) override;

    void AddEdep(G4double edep) { fEdep += edep; }

    void AddVolumeEntry(const VolumeEntryRecord& rec);
    void AddEdepRecord(const EdepRecord& rec);
    void AddDetectorSpectrumRecord(const DetectorSpectrumRecord& rec);
    void AddPhantomDoseRecord(const PhantomDoseRecord& rec);
    void AddAngularRecord(const AngularRecord& rec);

  private:
    RunAction* fRunAction = nullptr;
    G4double fEdep = 0.;

    std::vector<VolumeEntryRecord> fVolumeEntries;
    std::vector<EdepRecord> fEdepRecords;
    std::vector<DetectorSpectrumRecord> fDetectorSpectrumRecords;
    std::vector<PhantomDoseRecord> fPhantomDoseRecords;
    std::vector<AngularRecord> fAngularRecords;
};

}  // namespace B1

#endif