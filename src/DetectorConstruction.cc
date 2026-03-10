//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
/// \file DetectorConstruction.cc
/// \brief Implementation of the B1::DetectorConstruction class

#include "DetectorConstruction.hh"

#include "G4Box.hh"
#include "G4Element.hh"
#include "G4LogicalVolume.hh"
#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4PVPlacement.hh"
#include "G4SubtractionSolid.hh"
#include "G4SystemOfUnits.hh"
#include "G4Tubs.hh"

namespace B1
{

G4VPhysicalVolume* DetectorConstruction::Construct()
{
  G4NistManager* nist = G4NistManager::Instance();
  G4bool checkOverlaps = true;

  //
  // Materials
  //
  auto worldMat  = nist->FindOrBuildMaterial("G4_AIR");
  auto tungsten  = nist->FindOrBuildMaterial("G4_W");
  auto tissue    = nist->FindOrBuildMaterial("G4_TISSUE_SOFT_ICRP");
  auto lead      = nist->FindOrBuildMaterial("G4_Pb");
  auto argonGas  = nist->FindOrBuildMaterial("G4_Ar");

  // Plastic scintillator
  auto scintMat = nist->FindOrBuildMaterial("G4_POLYSTYRENE");

  //
  // Cherenkov gas radiator (simple CO2 gas)
  //
  auto elC = nist->FindOrBuildElement("C");
  auto elO = nist->FindOrBuildElement("O");

  G4double temperature = 293.15 * kelvin;
  G4double pressure    = 1.0 * atmosphere;
  G4double density     = 1.842e-3 * g / cm3;

  auto cherenkovGas =
    new G4Material("CherenkovGasCO2", density, 2, kStateGas, temperature, pressure);
  cherenkovGas->AddElement(elC, 1);
  cherenkovGas->AddElement(elO, 2);

  //
  // World
  //
  G4double worldXY = 120.0 * cm;
  G4double worldZ  = 300.0 * cm;

  auto solidWorld = new G4Box("World", worldXY / 2., worldXY / 2., worldZ / 2.);

  auto logicWorld = new G4LogicalVolume(solidWorld, worldMat, "World");

  auto physWorld = new G4PVPlacement(nullptr,
                                     G4ThreeVector(),
                                     logicWorld,
                                     "World",
                                     nullptr,
                                     false,
                                     0,
                                     checkOverlaps);

  //
  // Common transverse dimensions
  //
  G4double detXY = 10.0 * cm;

  //
  // Longitudinal thicknesses
  //
  G4double cherZ    = 10.0 * cm;

  // Tungsten radiation length X0 ≈ 3.5 mm
  // chosen target thickness = 2.5 X0 ≈ 8.75 mm
  G4double bremsTargetThickness = 8.75 * mm;

  G4double phantomZ = 3.0 * cm;
  G4double collZ    = 5.0 * cm;
  G4double scintZ   = 5.0 * mm;
  G4double dwcXY    = 20.0 * cm;
  G4double dwcZ     = 20.0 * cm;

  //
  // 1) Threshold Cherenkov detector (gas radiator)
  //
  auto solidCher = new G4Box("Cherenkov", detXY / 2., detXY / 2., cherZ / 2.);

  auto logicCher = new G4LogicalVolume(solidCher, cherenkovGas, "Cherenkov");

  new G4PVPlacement(nullptr,
                    G4ThreeVector(0, 0, -90.0 * cm),
                    logicCher,
                    "Cherenkov",
                    logicWorld,
                    false,
                    0,
                    checkOverlaps);

  //
  // 2) Tungsten brems target
  //
  auto solidBrems =
    new G4Box("BremsTarget", detXY / 2., detXY / 2., bremsTargetThickness / 2.);

  auto logicBrems = new G4LogicalVolume(solidBrems, tungsten, "BremsTarget");

  new G4PVPlacement(nullptr,
                    G4ThreeVector(0, 0, -75.0 * cm),
                    logicBrems,
                    "BremsTarget",
                    logicWorld,
                    false,
                    0,
                    checkOverlaps);

  //
  // 3) Tissue phantom
  //
  auto solidPhantom = new G4Box("Phantom", detXY / 2., detXY / 2., phantomZ / 2.);

  auto logicPhantom = new G4LogicalVolume(solidPhantom, tissue, "Phantom");

  new G4PVPlacement(nullptr,
                    G4ThreeVector(0, 0, -55.0 * cm),
                    logicPhantom,
                    "Phantom",
                    logicWorld,
                    false,
                    0,
                    checkOverlaps);

  //
  // 4) Collimator (lead block with cylindrical hole)
  //
  auto solidCollBlock = new G4Box("CollBlock", 6.0 * cm, 6.0 * cm, collZ / 2.);
  auto solidCollHole =
    new G4Tubs("CollHole", 0., 1.0 * cm, collZ / 2. + 1.0 * mm, 0., 360. * deg);

  auto solidCollimator =
    new G4SubtractionSolid("Collimator", solidCollBlock, solidCollHole);

  auto logicCollimator = new G4LogicalVolume(solidCollimator, lead, "Collimator");

  new G4PVPlacement(nullptr,
                    G4ThreeVector(0, 0, -35.0 * cm),
                    logicCollimator,
                    "Collimator",
                    logicWorld,
                    false,
                    0,
                    checkOverlaps);

  //
  // 5) Scintillator 1
  //
  auto solidScint1 = new G4Box("Scint1", detXY / 2., detXY / 2., scintZ / 2.);

  auto logicScint1 = new G4LogicalVolume(solidScint1, scintMat, "Scint1");

  new G4PVPlacement(nullptr,
                    G4ThreeVector(0, 0, -20.0 * cm),
                    logicScint1,
                    "Scint1",
                    logicWorld,
                    false,
                    0,
                    checkOverlaps);

  //
  // 6) DWC1
  //
  auto solidDWC1 = new G4Box("DWC1", dwcXY / 2., dwcXY / 2., dwcZ / 2.);

  auto logicDWC1 = new G4LogicalVolume(solidDWC1, argonGas, "DWC1");

  new G4PVPlacement(nullptr,
                    G4ThreeVector(0, 0, 5.0 * cm),
                    logicDWC1,
                    "DWC1",
                    logicWorld,
                    false,
                    0,
                    checkOverlaps);

  //
  // 7) DWC2
  //
  auto solidDWC2 = new G4Box("DWC2", dwcXY / 2., dwcXY / 2., dwcZ / 2.);

  auto logicDWC2 = new G4LogicalVolume(solidDWC2, argonGas, "DWC2");

  new G4PVPlacement(nullptr,
                    G4ThreeVector(0, 0, 30.0 * cm),
                    logicDWC2,
                    "DWC2",
                    logicWorld,
                    false,
                    0,
                    checkOverlaps);

  //
  // 8) Scintillator 2
  //
  auto solidScint2 = new G4Box("Scint2", detXY / 2., detXY / 2., scintZ / 2.);

  auto logicScint2 = new G4LogicalVolume(solidScint2, scintMat, "Scint2");

  new G4PVPlacement(nullptr,
                    G4ThreeVector(0, 0, 60.0 * cm),
                    logicScint2,
                    "Scint2",
                    logicWorld,
                    false,
                    0,
                    checkOverlaps);

  //
  // Scoring volume: phantom
  //
  fScoringVolume = logicPhantom;

  return physWorld;
}

}  // namespace B1