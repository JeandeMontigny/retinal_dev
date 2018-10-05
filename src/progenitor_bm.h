#ifndef PROGENITOR_BM_
#define PROGENITOR_BM_

#include "biodynamo.h"
#include "rgc_soma_bm.h"

namespace bdm {
using namespace std;

struct Progenitor_behaviour_BM : public BaseBiologyModule {
  Progenitor_behaviour_BM() : BaseBiologyModule(gNullEventId) {}

  /// Default event constructor
  template <typename TEvent, typename TBm>
  Progenitor_behaviour_BM(const TEvent& event, TBm* other, uint64_t new_oid = 0) {
  }

  /// Default event handler (exising biology module won't be modified on
  /// any event)
  template <typename TEvent, typename... TBms>
  void EventHandler(const TEvent&, TBms*...) {
  }

  template <typename T, typename TSimulation = Simulation<>>
  void Run(T* cell) {
    auto* sim = TSimulation::GetActive();
    auto* rm = sim->GetResourceManager();
    auto* random = sim->GetRandom();

    if (!init_) {
      dg_substanceGuide_ = rm->GetDiffusionGrid("progenitors_guide");
      init_ = true;
    }

    array<double, 3> gradient;
    dg_substanceGuide_->GetGradient(cell->GetPosition(), &gradient);
    double concentration = dg_substanceGuide_->GetConcentration(cell->GetPosition());

    if (cell->GetDiameter() < 7) {
      cell->ChangeVolume(100);
    }

    if (!rgc_generated_) {
      if (concentration < 0.004) {
        cell->UpdatePosition(gradient);
      }
      if (concentration > 0.00039) {
        auto&& daughterRGC = cell->Divide();
        daughterRGC->SetCellType(-1);
        daughterRGC->SetInternalClock(-200);
        daughterRGC->AddBiologyModule(RGC_axial_migration_BM());
        rgc_generated_ = true;
      }
    }
    else if (!migrated_back_) {
      cell->UpdatePosition(Math::ScalarMult(-1, gradient));
      if (concentration < 0.0001 + random->Uniform(-0.00002, 0.00002)) {
        migrated_back_ = true;
      }
    }


  }  // end Run()

 private:
  bool init_ = false;
  bool rgc_generated_ = false;
  bool migrated_back_ = false;
  DiffusionGrid* dg_substanceGuide_ = nullptr;
  ClassDefNV(Progenitor_behaviour_BM, 1);
};  // end biologyModule RGC_mosaic_BM

} // end namespace bdm

#endif
