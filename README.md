# Compton Scattering TOF Resolution — Geant4 Simulation

> **BL4S (Beam Line for Schools) 2026 — TED Ege College, Aydın, Turkey**  
> *Characterisation of Compton Scattering Geometry and Its Contribution to Timing Resolution in Time-of-Flight Detection Systems*

---

## Overview

This repository contains a Geant4 simulation of the experimental beamline proposed for CERN's **Beam Line for Schools 2026** competition. The simulation models a sequential detection chain designed to quantify how the geometric variability of Compton-scattered electron trajectories contributes to the integral timing resolution of a time-of-flight (TOF) system — a question directly relevant to the design of next-generation TOF-PET scanners used in oncological imaging.

The central research question is:

> *To what extent do scattering geometry and path-length variation contribute to the integral timing resolution of the system, independently of intrinsic detector jitter?*

---

## Beamline Layout

The simulated geometry follows the proposed experimental setup in order of particle traversal:

| Position (z) | Component | Material | Purpose |
|---|---|---|---|
| −120 cm | Electron gun | — | 0.5 GeV e⁻ source |
| −90 cm | Cherenkov detector | CO₂ gas | Beam tagging & energy verification |
| −75 cm | Bremsstrahlung target | Tungsten (2.5 X₀ ≈ 8.75 mm) | Gamma photon production |
| −55 cm | Tissue phantom | ICRP soft tissue | Clinical reference medium |
| −35 cm | Collimator | Lead (5 cm, 1 cm aperture) | Beam definition |
| −20 cm | Scintillator 1 (t₁) | Polystyrene | First TOF timestamp |
| +5 cm | DWC1 | Argon gas | Track reconstruction |
| +30 cm | DWC2 | Argon gas | Track reconstruction |
| +60 cm | Scintillator 2 (t₂) | Polystyrene | Second TOF timestamp |

The tungsten target thickness of 2.5 X₀ was chosen to maximise Bremsstrahlung photon yield while minimising photon self-absorption via the photoelectric effect and pair production.

---

## Physics

The simulation uses the **QBBC** physics list, which includes:

- Full electromagnetic physics (Livermore model for low energies)
- Hadronic interaction models (FTFP + Bertini cascade)
- Compton scattering via the Klein–Nishina model
- Bremsstrahlung with LPM effect above 1 GeV

The Compton kinematic relation governing the experiment is:

```
Δλ = λ' − λ = (h / mₑc)(1 − cosθ)
```

---

## Output Data

Each simulation run produces six CSV files:

| File | Contents |
|---|---|
| `volume_entries.csv` | Every particle entering each detector volume (eventID, particle, energy, position, process, origin volume) |
| `volume_edep.csv` | Energy depositions per volume per step |
| `detector_spectrum.csv` | Particle spectra at each detector |
| `photon_spectrum.csv` | Gamma-only subset of detector spectrum |
| `phantom_dose_distribution.csv` | 3D dose map within the tissue phantom |
| `angular_distribution.csv` | θ and φ of gamma photons entering each volume |

These outputs are designed to allow the geometric timing contribution (σ²_geo) and the intrinsic detector jitter contribution (σ²_det) to be separated via:

```
σ²_total = σ²_geo + σ²_det
```

---

## Prerequisites

- [Geant4](https://geant4.web.cern.ch/) ≥ 11.0 (built with Qt and OpenGL for visualisation)
- CMake ≥ 3.16
- A C++17-compatible compiler

Ensure the following Geant4 datasets are installed: `G4LEDATA`, `G4LEVELGAMMADATA`, `G4NEUTRONXSDATA`, `G4SAIDXSDATA`, `G4ENSDFSTATEDATA`.

---

## Build Instructions

```bash
# Clone the repository
git clone https://github.com/<your-org>/<repo-name>.git
cd <repo-name>

# Create a build directory
mkdir build && cd build

# Configure
cmake ..

# Build
make -j$(nproc)
```

---

## Running the Simulation

### Interactive mode (with visualisation)
```bash
./exampleB1
```
The visualiser opens automatically. Trajectories are colour-coded: **red** = e⁻, **green** = γ, **blue** = e⁺.

### Batch mode (recommended for data collection)
```bash
# Short test run (5 gamma + 1 proton events, verbose tracking)
./exampleB1 run1.mac

# Full statistics run (10,000 events per particle type)
./exampleB1 run2.mac

# Predefined test with output redirect
./exampleB1 exampleB1.in > exampleB1.out
```

---

## Project Structure

```
.
├── include/
│   ├── ActionInitialization.hh
│   ├── DetectorConstruction.hh
│   ├── EventAction.hh
│   ├── PrimaryGeneratorAction.hh
│   ├── RunAction.hh
│   └── SteppingAction.hh
├── src/
│   ├── ActionInitialization.cc
│   ├── DetectorConstruction.cc
│   ├── EventAction.cc
│   ├── PrimaryGeneratorAction.cc
│   ├── RunAction.cc
│   └── SteppingAction.cc
├── exampleB1.cc        # Main entry point
├── CMakeLists.txt
├── run1.mac            # Interactive test macro
├── run2.mac            # High-statistics batch macro
├── exampleB1.in        # Predefined batch input
├── vis.mac             # Visualisation settings
└── README.md
```

---

## Scientific Context

Compton scattering accounts for 70–80% of photon interactions in the energy range used in radiotherapy. In TOF-PET systems, the precision with which an annihilation vertex can be localised is bounded by timing resolution, expressed as:

```
Δx = c · Δt / 2
```

Improved timing resolution directly improves the signal-to-noise ratio and enables more targeted dose delivery to tumour tissue. This simulation provides a Monte Carlo baseline for separating the geometric path-length contribution from detector jitter — data that can inform the optimisation of sub-100 ps TOF-PET detectors.

---

## Team

**TED Ege College, Aydın, Turkey — BL4S 2026**

Arhan Hasan Ünsal · Atakan Korkmaz · Beren Duygu Yılmaz · Doruk Turan · Doruk Utku Tarim · Leman Ece Genclesen · Neva Yildizli · Kayra Sari · Arda Genc · Duru Sefa

Supervisor: Eda Erdogan

---

## References

1. Surti S, Karp JS. Advances in time-of-flight PET. *Phys Med Biol.* 2016.
2. Joseph S, Matthijs V, et al. Bremsstrahlung, Synchrotron Radiation, and Compton Scattering. *Rev. Mod. Phys.* 1970; 42:237.
3. Surti S, et al. Impact of TOF on PET/CT imaging. *J Nucl Med.* 2010; 51(2):237–245.
4. Conti M, et al. Impact of TOF on PET tumor detection. *J Nucl Med.* 2009; 50(8):1315–1323.
5. Tajima M, et al. Carbon target in 90° Compton spectroscopy. *J Nucl Sci Technol.* 2008; 45(8):760–765.
6. Llosá G, et al. Noise evaluation of Compton camera imaging. *Phys Med Biol.* 2015; 60(5):1845.
7. Shimazoe K, et al. TOF-PET using Compton scattering by plastic scintillators. *VCI Proceedings.* 2016.
