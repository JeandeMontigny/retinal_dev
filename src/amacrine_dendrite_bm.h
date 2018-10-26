#ifndef AMACRINE_DENDRITE_BM_
#define AMACRINE_DENDRITE_BM_

#include "biodynamo.h"

namespace bdm {
using namespace std;

// Define dendrites behavior for Amacrine dendritic growth
struct Amacrine_dendrite_growth_BM : public BaseBiologyModule {
  Amacrine_dendrite_growth_BM() : BaseBiologyModule({
    experimental::neuroscience::NeuriteBifurcationEvent::kEventId
  }) {}

  /// Default event constructor
  template <typename TEvent, typename TBm>
  Amacrine_dendrite_growth_BM(const TEvent& event, TBm* other,
    uint64_t new_oid = 0) : BaseBiologyModule(event, other, new_oid) {
  }

  template <typename T, typename TSimulation = Simulation<>>
  void Run(T* ne) {

    if (ne->IsTerminal() && ne->GetDiameter() >= 0.5) {
      auto* sim = TSimulation::GetActive();
      auto* random = sim->GetRandom();

      if (!init_) {
        auto* rm = sim->GetResourceManager();
        dg_0_ = rm->GetDiffusionGrid("on_diffusion");
        dg_1_ = rm->GetDiffusionGrid("off_diffusion");
        dg_2_ = rm->GetDiffusionGrid("on_off_diffusion");
        init_ = true;
      }

      array<double, 3> position = ne->GetPosition();

      array<double, 3> gradient_0_;
      array<double, 3> gradient_1_;
      array<double, 3> gradient_2_;
      dg_0_->GetGradient(position, &gradient_0_);
      dg_1_->GetGradient(position, &gradient_1_);
      dg_2_->GetGradient(position, &gradient_2_);

      array<double, 3> gradient_combined = Math::Add(
          Math::Add(gradient_0_, gradient_1_), gradient_2_);

      double gradientWeight = 0.2;
      double randomnessWeight = 0.2;
      double oldDirectionWeight = 1.6;
      array<double, 3> random_axis = {random->Uniform(-1, 1),
                                      random->Uniform(-1, 1),
                                      random->Uniform(-1, 1)};
      auto oldDirection =
          Math::ScalarMult(oldDirectionWeight, ne->GetSpringAxis());
      auto gradDirection = Math::ScalarMult(
          gradientWeight, Math::Normalize(gradient_combined));
      auto randomDirection =
          Math::ScalarMult(randomnessWeight, random_axis);
      array<double, 3> newStepDirection = Math::Add(
          Math::Add(oldDirection, randomDirection), gradDirection);

      ne->ElongateTerminalEnd(25, newStepDirection);
      ne->SetDiameter(ne->GetDiameter()-0.00075);
    }

  }  // end run

 private:
   bool init_ = false;
   DiffusionGrid* dg_0_ = nullptr;
   DiffusionGrid* dg_1_ = nullptr;
   DiffusionGrid* dg_2_ = nullptr;
   ClassDefNV(Amacrine_dendrite_growth_BM, 1);
}; // end Amacrine_dendrite_growth_BM

} // end namespace bdm

#endif
