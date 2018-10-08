#ifndef AMACRINE_SOMA_BM_
#define AMACRINE_SOMA_BM_

namespace bdm {
using namespace std;

// Define cell behavior for neurite creation
struct Amacrine_Neurite_creation_BM: public BaseBiologyModule {
  Amacrine_Neurite_creation_BM() : BaseBiologyModule(gNullEventId) {}

  /// Default event constructor
  template <typename TEvent, typename TBm>
  Amacrine_Neurite_creation_BM(const TEvent& event, TBm* other, uint64_t new_oid = 0) {}

  template <typename TEvent, typename... TBms>
  void EventHandler(const TEvent&, TBms*...) {}

  template <typename T, typename TSimulation = Simulation<>>
  void Run(T* soma) {

    bool createDendrites = false;

    if (createDendrites) {
      //TODO: add correct number of dendrites
      soma->ExtendNewNeurite({0,0,-1});
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
      while (temps > -1) {
        temps = temps * 10;
        multiple ++;
      }
      zDirection = zDirection * pow(10, multiple-1);
    }

    // migrate to INL
    if (concentration < 1e-10 + random->Uniform(-2e-11, 2e-11)) {
      cell->UpdatePosition({0, 0, zDirection});
    }
    else {
      //TODO: seg fault
      // cell->AddBiologyModule(Amacrine_Neurite_creation_BM());
      //TODO: remove that BM
    }
  } // end run

private:
  bool init_ = false;
  DiffusionGrid* dg_0_ = nullptr;
  DiffusionGrid* dg_1_ = nullptr;
  DiffusionGrid* dg_2_ = nullptr;
  ClassDefNV(Amacrine_axial_migration_BM, 1);
}; // Amacrine_axial_migration_BM


} // end namespace bdm

#endif
