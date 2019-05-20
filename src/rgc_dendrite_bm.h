#ifndef RGC_DENDRITE_BM_
#define RGC_DENDRITE_BM_

#include "biodynamo.h"

namespace bdm {
using namespace std;

  // Define dendrites behavior for RGC dendritic growth
  struct RGC_dendrite_growth_BM : public BaseBiologyModule {
    RGC_dendrite_growth_BM() : BaseBiologyModule({
     experimental::neuroscience::NeuriteBifurcationEvent::kEventId
    }) {}

    /// Default event constructor
    template <typename TEvent, typename TBm>
    RGC_dendrite_growth_BM(const TEvent& event, TBm* other,
                           uint64_t new_oid = 0) : BaseBiologyModule(event, other, new_oid) {
    }

    template <typename T, typename TSimulation = Simulation<>>
    void Run(T* ne) {
      auto* sim = TSimulation::GetActive();
      auto* random = sim->GetRandom();
      auto* rm = sim->GetResourceManager();

      if (!init_) {
        dg_on_RGCguide_ =
          rm->GetDiffusionGrid("on_substance_RGC_guide");
        dg_off_RGCguide_ =
          rm->GetDiffusionGrid("off_substance_RGC_guide");
        init_ = true;
      }

      if (ne->IsTerminal() && ne->GetDiameter() >= 0.5) {

        int cellClock = ne->GetMySoma()->GetInternalClock();

        // -------------------------------- //

        array<double, 3> gradient_RGCguide;
        double concentration;

        double conc_on =
          dg_on_RGCguide_->GetConcentration(ne->GetPosition());
        double conc_off =
          dg_off_RGCguide_->GetConcentration(ne->GetPosition());
        if (conc_on*0.2 > conc_off) {
          concentration = conc_on;
          dg_on_RGCguide_->GetGradient(ne->GetPosition(),
            &gradient_RGCguide);
          gradient_RGCguide = Math::ScalarMult(1, gradient_RGCguide);
        } else {
          concentration = conc_off;
          dg_off_RGCguide_->GetGradient(ne->GetPosition(),
            &gradient_RGCguide);
        }


        if (cellClock < 500) {

          // ---------------- //

          // array<double, 3> gradient_on;
          // array<double, 3> gradient_off;
          // array<double, 3> gradient_sum;
          //
          // double conc_on =
          //   dg_on_RGCguide_->GetConcentration(ne->GetPosition());
          // double conc_off =
          //   dg_off_RGCguide_->GetConcentration(ne->GetPosition());
          // double conc_sum = conc_on + conc_off;
          //
          // dg_on_RGCguide_->GetGradient(ne->GetPosition(),
          //   &gradient_on);
          // dg_off_RGCguide_->GetGradient(ne->GetPosition(),
          //   &gradient_off);
          //   for (int i = 0; i < 3; i++) {
          //     gradient_sum[i] = gradient_on[i] + gradient_off[i];
          //     // gradient_sum[i] = max(gradient_on[i], gradient_off[i]);
          //   }
          //
          //   double bifurcProba = 0.015*ne->GetDiameter();
          //
          //   double gradientWeight = 0.5; // 0.2
          //   double randomnessWeight = 0.2; // 0.5
          //   double oldDirectionWeight = 1; // 4.5
          //   array<double, 3> random_axis = {random->Uniform(-1, 1),
          //                                   random->Uniform(-1, 1),
          //                                   random->Uniform(-1, 1)};
          //
          //   auto oldDirection =
          //     Math::ScalarMult(oldDirectionWeight, ne->GetSpringAxis());
          //   auto gradDirection = Math::ScalarMult(
          //     gradientWeight, Math::Normalize(gradient_sum));
          //   auto randomDirection =
          //     Math::ScalarMult(randomnessWeight, random_axis);
          //   array<double, 3> newStepDirection = Math::Add(
          //     Math::Add(oldDirection, randomDirection), gradDirection);
          //
          //   ne->ElongateTerminalEnd(25, newStepDirection);
          //   ne->SetDiameter(ne->GetDiameter()-0.0007);
          //
          //   if (conc_sum > 0.06 && random->Uniform() < bifurcProba) {
          //     ne->SetDiameter(ne->GetDiameter()-0.005);
          //     ne->Bifurcate();
          //   }

            // ---------------- //

          double bifurcProba = 0.012*ne->GetDiameter();

          double gradientWeight = 0.15; // 0.2
          double randomnessWeight = 0.6; // 0.5
          double oldDirectionWeight = 4.5; // 4.5
          array<double, 3> random_axis = {random->Uniform(-1, 1),
                                          random->Uniform(-1, 1),
                                          random->Uniform(-1, 1)};
          auto oldDirection =
            Math::ScalarMult(oldDirectionWeight, ne->GetSpringAxis());
          auto gradDirection = Math::ScalarMult(
            gradientWeight, Math::Normalize(gradient_RGCguide));
          auto randomDirection =
            Math::ScalarMult(randomnessWeight, random_axis);
          array<double, 3> newStepDirection = Math::Add(
            Math::Add(oldDirection, randomDirection), gradDirection);

          ne->ElongateTerminalEnd(25, newStepDirection);
          ne->SetDiameter(ne->GetDiameter()-0.0007);

          if (concentration > 0.04 && random->Uniform() < bifurcProba) {
            ne->SetDiameter(ne->GetDiameter()-0.005);
            ne->Bifurcate();
          }

        } // end if cellClock

        else if (cellClock > 600) {

          if (conc_on < 0.065 && conc_off < 0.02) {
            ne->RetractTerminalEnd(40);
          }

          // if (concentration > 0.04 && random->Uniform() < 0.001
          //     && ne->GetDaughterRight() == nullptr) {
          //   ne->Branch(ne->GetDiameter()-0.005);
          // }

        }


        // -------------------------------- //

        // array<double, 3> gradient_RGCguide;
        // double concentration = 0;
        // // initialise the correct substance as guide depending on cell type
        // if (ne->GetSubtype()/100 == 0) {
        //   dg_on_RGCguide_->GetGradient(ne->GetPosition(),
        //     &gradient_RGCguide);
        //   concentration =
        //     dg_on_RGCguide_->GetConcentration(ne->GetPosition());
        // }
        // if (ne->GetSubtype()/100 == 1) {
        //   dg_off_RGCguide_->GetGradient(ne->GetPosition(),
        //     &gradient_RGCguide);
        //   concentration =
        //     dg_off_RGCguide_->GetConcentration(ne->GetPosition());
        // }
        // if (ne->GetSubtype()/100 == 2) {
        //   double conc_on =
        //     dg_on_RGCguide_->GetConcentration(ne->GetPosition());
        //   double conc_off =
        //     dg_off_RGCguide_->GetConcentration(ne->GetPosition());
        //   if (conc_on > conc_off) {
        //     concentration = conc_on;
        //     dg_on_RGCguide_->GetGradient(ne->GetPosition(),
        //       &gradient_RGCguide);
        //   } else {
        //     concentration = conc_off;
        //     dg_off_RGCguide_->GetGradient(ne->GetPosition(),
        //       &gradient_RGCguide);
        //   }
        // }
        //
        // // if neurite doesn't have to retract
        // if (!ne->GetHasToRetract()) {
        //   double bifurcProba = 0.01*ne->GetDiameter();
        //
        //   double gradientWeight = 0.2;
        //   double randomnessWeight = 0.5;
        //   double oldDirectionWeight = 4.5;
        //   array<double, 3> random_axis = {random->Uniform(-1, 1),
        //                                   random->Uniform(-1, 1),
        //                                   random->Uniform(-1, 1)};
        //   auto oldDirection =
        //     Math::ScalarMult(oldDirectionWeight, ne->GetSpringAxis());
        //   auto gradDirection = Math::ScalarMult(
        //     gradientWeight, Math::Normalize(gradient_RGCguide));
        //   auto randomDirection =
        //     Math::ScalarMult(randomnessWeight, random_axis);
        //   array<double, 3> newStepDirection = Math::Add(
        //     Math::Add(oldDirection, randomDirection), gradDirection);
        //
        //   ne->ElongateTerminalEnd(25, newStepDirection);
        //   ne->SetDiameter(ne->GetDiameter()-0.0007);
        //
        //   if (concentration > 0.04 && random->Uniform() < bifurcProba) {
        //     ne->SetDiameter(ne->GetDiameter()-0.005);
        //     ne->Bifurcate();
        //   }
        //
        //   // homo-type interaction
        //   double squared_radius = 1.1; // 1.2
        //   int sameType = 0;
        //   int otherType = 0;
        //   // lambda updating counters for neurites neighbours
        //   auto countNeighbours = [&](const auto* neighbor) {
        //     // if neighbor is a NeuriteElement
        //     if (neighbor->template IsSoType<MyNeurite>()) {
        //       auto n_soptr = neighbor->template
        //         ReinterpretCast<MyNeurite>()->GetSoPtr();
        //       // if neurites have not the same soma
        //       if (!(n_soptr->GetMySoma() == ne->GetMySoma())) {
        //         // if neurites got the same type
        //         if (n_soptr->GetSubtype() == ne->GetSubtype()) {
        //           sameType++;
        //         }
        //         else {
        //           otherType++;
        //         }
        //       }
        //       else {
        //         sameType--;
        //       }
        //     }
        //   }; // end lambda
        //
        //   auto* ctxt = TSimulation::GetActive()->GetExecutionContext();
        //   ctxt->ForEachNeighborWithinRadius(countNeighbours, *ne, squared_radius);
        //   // if is surrounded by homotype dendrites
        //   if (sameType > otherType) {
        //     ne->SetHasToRetract(true);
        //     ne->SetDiamBeforeRetraction(ne->GetDiameter());
        //   }
        //
        //   // if neurite is going too far away from guide
        //   if (concentration < 0.01 && ne->GetDiameter() < 0.9) {
  			// 		ne->SetHasToRetract(true);
  			// 		ne->SetBeyondThreshold(true);
  			// 	}
        //
        // } // if ! has to retract
        //
        // // if neurite has to retract
        // else {
        //   ne->RetractTerminalEnd(40);
        //   // if neurite has retracted enough because of interactions
        // 	if (!ne->GetBeyondThreshold()
        //     && ne->GetDiameter() > ne->GetDiamBeforeRetraction()+0.02) {
        // 	  ne->SetSleepMode(true);
        //      // ne->RemoveBiologyModule(
        //      //   ne->GetBiologyModules<RGC_dendrite_growth_BM>()[0]);
        // 	}
        //   // if neurite is back to higher concentration
        // 	if (ne->GetBeyondThreshold() && concentration>0.02) {
        // 		ne->SetBeyondThreshold(false);
        // 		ne->SetHasToRetract(false);
        // 	}
        // } // end has to retract

      } // if is terminal

      // else {
      //   if (ne->GetMySoma()->GetInternalClock() > 600) {
      //     array<double, 3> gradient_RGCguide;
      //     double conc_on =
      //       dg_on_RGCguide_->GetConcentration(ne->GetPosition());
      //     double conc_off =
      //       dg_off_RGCguide_->GetConcentration(ne->GetPosition());
      //
      //     if (conc_on*0.4 > conc_off) {
      //       dg_on_RGCguide_->GetGradient(ne->GetPosition(),
      //         &gradient_RGCguide);
      //       gradient_RGCguide = Math::ScalarMult(1, gradient_RGCguide);
      //     } else {
      //       dg_off_RGCguide_->GetGradient(ne->GetPosition(),
      //         &gradient_RGCguide);
      //     }
      //
      //     if (conc_on < 0.075 && conc_off < 0.03) {
      //      ne->MovePointMass(10, gradient_RGCguide);
      //     }
      //   }
      // } // end if ne is not terminal

    }  // end run

   private:
    bool init_ = false;
    DiffusionGrid* dg_on_RGCguide_ = nullptr;
    DiffusionGrid* dg_off_RGCguide_ = nullptr;
    ClassDefNV(RGC_dendrite_growth_BM, 1);
  }; // end RGC_dendrite_growth_BM

} // end namespace bdm

#endif
