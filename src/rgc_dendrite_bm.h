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

      if (ne->IsTerminal() && ne->GetDiameter() > ne->GetDiamLimit()) {

        int cellClock = ne->GetMySoma()->GetInternalClock();
        array<double, 3> gradient_RGCguide;
        double concentration = 0;
        double conc_on =
          dg_on_RGCguide_->GetConcentration(ne->GetPosition());
        double conc_off =
          dg_off_RGCguide_->GetConcentration(ne->GetPosition());
        if (conc_on*0.8 > conc_off) { // 0.4
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
          // bi-stratified proportion: 2/3 at the begining -> 1/3 at the end
          // 50% of on and off cells will be mono-stratified
          if (ne->GetSubtype()/100 == 0 && ne->GetSubtype()%2 == 0) {
            dg_on_RGCguide_->GetGradient(ne->GetPosition(),
              &gradient_RGCguide);
            concentration =
              dg_on_RGCguide_->GetConcentration(ne->GetPosition());
          } // end if on cell
          if (ne->GetSubtype()/100 == 1 && ne->GetSubtype()%2 == 0) {
            dg_off_RGCguide_->GetGradient(ne->GetPosition(),
              &gradient_RGCguide);
            concentration =
              dg_off_RGCguide_->GetConcentration(ne->GetPosition());
          } // end if off cell

          double bifurcProba = 0.0135*ne->GetDiameter(); // 0.0133

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

          if (concentration > 0.04 && ne->GetDiameter() > 0.55
              && random->Uniform() < bifurcProba) {
            ne->Bifurcate();
          }

        } // end if cellClock

        else if (cellClock < 1000) {

            array<double, 3> gradient_RGCguide;
            double concentration = 0;
            // initialise the correct substance as guide depending on cell type
            if (ne->GetSubtype()/100 == 0) {
              dg_on_RGCguide_->GetGradient(ne->GetPosition(),
                &gradient_RGCguide);
              concentration =
                dg_on_RGCguide_->GetConcentration(ne->GetPosition());
              // shrinkage
              if (conc_on < 0.055) { // 0.06
                ne->SetHasToRetract(true);
                ne->SetDiamBeforeRetraction(ne->GetDiameter());
              }
            } // end if on cell

            if (ne->GetSubtype()/100 == 1) {
              dg_off_RGCguide_->GetGradient(ne->GetPosition(),
                &gradient_RGCguide);
              concentration =
                dg_off_RGCguide_->GetConcentration(ne->GetPosition());
              // shrinkage
              if (conc_off < 0.04) {
                ne->SetHasToRetract(true);
                ne->SetDiamBeforeRetraction(ne->GetDiameter());
              }
            } // end if off cell

            if (ne->GetSubtype()/100 == 2) {
              double conc_on =
                dg_on_RGCguide_->GetConcentration(ne->GetPosition());
              double conc_off =
                dg_off_RGCguide_->GetConcentration(ne->GetPosition());
              if (conc_on*0.8 > conc_off) {
                concentration = conc_on;
                dg_on_RGCguide_->GetGradient(ne->GetPosition(),
                  &gradient_RGCguide);
              } else {
                concentration = conc_off;
                dg_off_RGCguide_->GetGradient(ne->GetPosition(),
                  &gradient_RGCguide);
              }
              // shrinkage
              if (conc_on < 0.055 && conc_off < 0.04) { // 0.06 ; 0.04
                ne->SetHasToRetract(true);
                ne->SetDiamBeforeRetraction(ne->GetDiameter());
              }
            } // end if on-off cell


            // branching (not bifurcation)
            // if (random->Uniform() < 0.001 && ne->GetDaughterRight() == nullptr) {
            //   ne->Branch(ne->GetDiameter()-0.005);
            // }

            // normal grwoth/retract routine

            // if neurite doesn't have to retract
            if (!ne->GetHasToRetract()) {
              double bifurcProba = 0.011*ne->GetDiameter();

              double gradientWeight = 0.3; // 0.2
              double randomnessWeight = 0.5;
              double oldDirectionWeight = 4.5;
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

              if ((conc_on > 0.04 || conc_off > 0.016) && ne->GetDiameter() > 0.55
                  && random->Uniform() < bifurcProba) {
                ne->Bifurcate();
              }

              // homo-type interaction
              double squared_radius = 1; // 1.2
              int sameType = 0;
              int otherType = 0;
              // lambda updating counters for neurites neighbours
              auto countNeighbours = [&](const auto* neighbor) {
                // if neighbor is a NeuriteElement
                if (neighbor->template IsSoType<MyNeurite>()) {
                  auto n_soptr = neighbor->template
                    ReinterpretCast<MyNeurite>()->GetSoPtr();
                  // if neurites have not the same soma
                  if (!(n_soptr->GetMySoma() == ne->GetMySoma())) {
                    // if neurites got the same type
                    if (n_soptr->GetSubtype() == ne->GetSubtype()) {
                      sameType++;
                    }
                    else {
                      otherType++;
                    }
                  }
                  else {
                    sameType--;
                  }
                }
              }; // end lambda

              auto* ctxt = TSimulation::GetActive()->GetExecutionContext();
              ctxt->ForEachNeighborWithinRadius(countNeighbours, *ne, squared_radius);
              // if is surrounded by homotype dendrites
              if (sameType > otherType) {
                ne->SetHasToRetract(true);
                ne->SetDiamBeforeRetraction(ne->GetDiameter());
              }

              // if neurite is going too far away from guide
              // if (concentration < 0.01 && ne->GetDiameter() < 0.9) {
      				// 	ne->SetHasToRetract(true);
      				// 	ne->SetBeyondThreshold(true);
      				// }

            } // if ! has to retract

            // if neurite has to retract
            else {
              ne->RetractTerminalEnd(40);
              // if neurite has retracted enough because of interactions
            	if (!ne->GetBeyondThreshold()
                && ne->GetDiameter() > ne->GetDiamBeforeRetraction()+0.02) {
            	  ne->SetSleepMode(true);
                 // ne->RemoveBiologyModule(
                 //   ne->GetBiologyModules<RGC_dendrite_growth_BM>()[0]);
            	}
              // if neurite is back to higher concentration
            	if (ne->GetBeyondThreshold() && concentration > 0.02) {
            		ne->SetBeyondThreshold(false);
            		ne->SetHasToRetract(false);
            	}
            } // end has to retract


          } // end if cellClock > 600


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
