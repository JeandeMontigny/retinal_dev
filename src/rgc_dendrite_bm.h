#ifndef RGC_DENDRITE_BM_
#define RGC_DENDRITE_BM_

#include "biodynamo.h"

namespace bdm {
using namespace std;

  // Define dendrites behavior for RGC dendritic growth
  struct RGC_dendrite_growth_BM : public BaseBiologyModule {
    BDM_STATELESS_BM_HEADER(RGC_dendrite_growth_BM, BaseBiologyModule, 1);

  public:
    RGC_dendrite_growth_BM() : BaseBiologyModule({
     experimental::neuroscience::NeuriteBifurcationEvent::kEventId
    }) {}

    void Run(SimObject* so) override {
      if (auto* ne = so->As<MyNeurite>()) {
        auto* sim = Simulation::GetActive();
        auto* random = sim->GetRandom();
        auto* rm = sim->GetResourceManager();

        if (ne->IsTerminal() && ne->GetDiameter() >= 0.5) {
          DiffusionGrid* dg_on_RGCguide_ =
            rm->GetDiffusionGrid("on_substance_RGC_guide");;
          DiffusionGrid* dg_off_RGCguide_ =
            rm->GetDiffusionGrid("off_substance_RGC_guide");;

          array<double, 3> gradient_RGCguide;
          double concentration = 0;
          array<double, 3> dendritePosition = ne->GetPosition();
          // initialise the correct substance as guide depending on cell type
          if (ne->GetSubtype()/100 == 0) {
            dg_on_RGCguide_->GetGradient(dendritePosition,
              &gradient_RGCguide);
            concentration =
              dg_on_RGCguide_->GetConcentration(dendritePosition);
          }
          if (ne->GetSubtype()/100 == 1) {
            dg_off_RGCguide_->GetGradient(dendritePosition,
              &gradient_RGCguide);
            concentration =
              dg_off_RGCguide_->GetConcentration(dendritePosition);
          }
          if (ne->GetSubtype()/100 == 2) {
            double conc_on =
              dg_on_RGCguide_->GetConcentration(dendritePosition);
            double conc_off =
              dg_off_RGCguide_->GetConcentration(dendritePosition);
            if (conc_on > conc_off) {
              concentration = conc_on;
              dg_on_RGCguide_->GetGradient(dendritePosition,
                &gradient_RGCguide);
            } else {
              concentration = conc_off;
              dg_off_RGCguide_->GetGradient(dendritePosition,
                &gradient_RGCguide);
            }
          }

          // if neurite doesn't have to retract
          if (!ne->GetHasToRetract()) {
            double bifurcProba = 0.0108*ne->GetDiameter();

            double gradientWeight = 0.2;
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
            ne->SetDiameter(ne->GetDiameter()-0.00075);

            if (concentration > 0.04 && random->Uniform() < bifurcProba) {
              ne->SetDiameter(ne->GetDiameter()-0.005);
              ne->Bifurcate();
            }

            // homo-type interaction
            double squared_radius = 1;
            int sameType = 0;
            int otherType = 0;
            // lambda updating counters for neurites neighbours
            auto countNeighbours = [&](const auto* neighbor) {
              // if neighbor is a NeuriteElement
              if (auto* neighbor_asMN = neighbor->template As<MyNeurite>()) {
                auto n_soptr = neighbor_asMN->GetSoPtr();
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

            auto* ctxt = Simulation::GetActive()->GetExecutionContext();
            ctxt->ForEachNeighborWithinRadius(countNeighbours, *ne, squared_radius);
            // if is surrounded by homotype dendrites
            if (sameType > otherType) {
              ne->SetHasToRetract(true);
              ne->SetDiamBeforeRetraction(ne->GetDiameter());
            }

            // if neurite is going too far away from guide
            if (concentration < 0.01 && ne->GetDiameter() < 0.9) {
    					ne->SetHasToRetract(true);
    					ne->SetBeyondThreshold(true);
    				}

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
          	if (ne->GetBeyondThreshold() && concentration>0.02) {
          		ne->SetBeyondThreshold(false);
          		ne->SetHasToRetract(false);
          	}
          } // end has to retract

        } // if is terminal
      }
    }  // end run

  }; // end RGC_dendrite_growth_BM

} // end namespace bdm

#endif
