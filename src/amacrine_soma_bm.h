#ifndef AMACRINE_SOMA_BM_
#define AMACRINE_SOMA_BM_

#include "amacrine_dendrite_bm.h"
#include "extended_objects.h"

namespace bdm {
using namespace std;

// Define Amacrine cell behavior for neurite creation
struct Amacrine_Neurite_creation_BM: public BaseBiologyModule {
  Amacrine_Neurite_creation_BM() : BaseBiologyModule(gNullEventId) {}

  /// Default event constructor
  template <typename TEvent, typename TBm>
  Amacrine_Neurite_creation_BM(const TEvent& event, TBm* other, uint64_t new_oid = 0) {}

  template <typename TEvent, typename... TBms>
  void EventHandler(const TEvent&, TBms*...) {}

  template <typename T, typename TSimulation = Simulation<>>
  void Run(T* soma) {
    if (soma->GetDaughters().size()==0) {
      auto* sim = TSimulation::GetActive();
      auto* random = sim->GetRandom();
      //TODO: add correct number of dendrites: 6.4, std 1.6
      for (int i = 0; i <= (int)random->Uniform(2, 7); i++) {
        auto&& ne = soma->ExtendNewNeurite({0, 0, -1});
        ne->AddBiologyModule(Amacrine_dendrite_growth_BM());
        ne->SetMySoma(soma->GetSoPtr());
      }
      //TODO: remove this BM
    }
  } // end run

private:
  ClassDefNV(Amacrine_Neurite_creation_BM, 1);
}; // end Amacrine_Neurite_creation_BM


// Define cell behavior for RGC migration to GCL
struct Amacrine_axial_migration_BM: public BaseBiologyModule {
  Amacrine_axial_migration_BM() : BaseBiologyModule(gNullEventId) {}

  /// Default event constructor
  template <typename TEvent, typename TBm>
  Amacrine_axial_migration_BM(const TEvent& event, TBm* other, uint64_t new_oid = 0) {}

  template <typename TEvent, typename... TBms>
  void EventHandler(const TEvent&, TBms*...) {}

  template <typename T, typename TSimulation = Simulation<>>
  void Run(T* cell) {
    if (!migrated_) {
      auto* sim = TSimulation::GetActive();
      auto* random = sim->GetRandom();

      if (!init_) {
        auto* rm = sim->GetResourceManager();
        dg_0_ = rm->GetDiffusionGrid("on_diffusion");
        dg_1_ = rm->GetDiffusionGrid("off_diffusion");
        dg_2_ = rm->GetDiffusionGrid("on_off_diffusion");
        init_ = true;
      }

      auto& position = cell->GetPosition();
      double concentration = dg_0_->GetConcentration(position)
        + dg_1_->GetConcentration(position) + dg_2_->GetConcentration(position);

      array<double, 3> gradient_0_;
      array<double, 3> gradient_1_;
      array<double, 3> gradient_2_;
      dg_0_->GetGradient(position, &gradient_0_);
      dg_1_->GetGradient(position, &gradient_1_);
      dg_2_->GetGradient(position, &gradient_2_);
      double zDirection = gradient_0_[2] + gradient_1_[2] + gradient_2_[2];

      if (zDirection < 0) {
        int multiple = 0;
        double temps = zDirection;
        while (temps > -0.1) {
          temps = temps * 10;
          multiple ++;
        }
        zDirection = zDirection * pow(10, multiple);
      }

      if (cell->GetDiameter() < 8 && random->Uniform(0, 1) < 0.01) {
        cell->ChangeVolume(2000);
      }

      // migrate to INL ~ 350
      if (concentration < 8e-3) {
        cell->UpdatePosition({random->Uniform(-0.1, 0.1), random->Uniform(-0.1, 0.1), zDirection});
      }
      else {
        //TODO: remove that BM
        migrated_ = true;
      }

    } // end if not migrated_
    else {
      //TODO: seg fault
      // cell->AddBiologyModule(Amacrine_Neurite_creation_BM());
    }
  } // end run

private:
  bool init_ = false;
  bool migrated_ = false;
  DiffusionGrid* dg_0_ = nullptr;
  DiffusionGrid* dg_1_ = nullptr;
  DiffusionGrid* dg_2_ = nullptr;
  ClassDefNV(Amacrine_axial_migration_BM, 1);
}; // Amacrine_axial_migration_BM


} // end namespace bdm

#endif
