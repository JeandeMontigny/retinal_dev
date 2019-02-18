#ifndef RGC_SOMA_BM_
#define RGC_SOMA_BM_

#include "rgc_dendrite_bm.h"
#include "extended_objects.h"

namespace bdm {
using namespace std;

  // Define cell behavior for mosaic formation
  struct RGC_mosaic_BM : public BaseBiologyModule {
    BDM_STATELESS_BM_HEADER(RGC_mosaic_BM, BaseBiologyModule, 1)

  public:
    RGC_mosaic_BM() : BaseBiologyModule(gNullEventId) {}

    void Run(SimObject* so) override {
      if (auto* cell = so->As<MyCell>()) {
        auto* sim = Simulation::GetActive();
        auto* random = sim->GetRandom();
        auto* rm = sim->GetResourceManager();
        DiffusionGrid* dg = nullptr;

        // run tangential migration
        bool withMovement = true;
        double movementThreshold = 1.21;
        // run cell death
        bool withDeath = true;
        double deathThreshold = 1.216;

        // initialise substance diffusions
        if (cell->GetCellType() == 0) {
          dg = rm->GetDiffusionGrid("on_diffusion");
        }
        else if (cell->GetCellType() == 1) {
          dg = rm->GetDiffusionGrid("off_diffusion");
        }
        else if (cell->GetCellType() == 2) {
          dg = rm->GetDiffusionGrid("on_off_diffusion");
        }

        auto& position = cell->GetPosition();
        array<double, 3> gradient_;
        array<double, 3> diff_gradient;
        array<double, 3> gradient_z;
        double concentration = 0;
        int cellClock = cell->GetInternalClock();

        // get concentration and gradient of correct substance
        if (dg != nullptr) {
          dg->GetGradient(position, &gradient_);
          gradient_z = Math::ScalarMult(0.2, gradient_);
          gradient_z[0] = 0;
          gradient_z[1] = 0;
          diff_gradient = Math::ScalarMult(-0.1, gradient_);
          diff_gradient[2] = 0;
          concentration = dg->GetConcentration(position);
        }

        if (cellClock < 700) {
          // // add small random movements
          // cell->UpdatePosition(
          //     {random->Uniform(-0.01, 0.01), random->Uniform(-0.01, 0.01), 0});
          // cell growth
          if (cell->GetDiameter() < 14 && random->Uniform(0, 1) < 0.02) {
            cell->ChangeVolume(5500);
          }
        }

        /* -- cell movement -- */
        if (withMovement && cellClock >= 100 && cellClock < 900
            && concentration >= movementThreshold) {
          // cell movement based on homotype substance gradient
            cell->UpdatePosition(diff_gradient);
            // update distance travelled by this cell
            array<double, 3> previousPosition = cell->GetPreviousPosition();
            array<double, 3> currentPosition = cell->GetPosition();
            cell->SetDistanceTravelled(cell->GetDistanceTravelled() +
              (sqrt(pow(currentPosition[0] - previousPosition[0], 2) +
                   pow(currentPosition[1] - previousPosition[1], 2))));
            cell->SetPreviousPosition(cell->GetPosition());
        }  // end tangential migration

        /* -- cell death -- */
        if (cell->GetCellType() != -1 && withDeath
            && cellClock >= 100 && cellClock < 700) {
          // add vertical migration as the multi layer colapse in just on layer
          cell->UpdatePosition(gradient_z);
          // cell death depending on homotype substance concentration
          if (concentration > deathThreshold && random->Uniform(0, 1) < 0.05) {
            cell->RemoveFromSimulation();
          }
        } // end cell death

        /* -- cell fate -- */
        // cell type attribution depending on concentrations
        // if cell in undifferentiated; random so don't change simultaneously
        if (cell->GetCellType() == -1 && random->Uniform(0, 1) < 0.05) {
          DiffusionGrid* dg_0_ = rm->GetDiffusionGrid("on_diffusion");;
          DiffusionGrid* dg_1_ = rm->GetDiffusionGrid("off_diffusion");;
          DiffusionGrid* dg_2_ = rm->GetDiffusionGrid("on_off_diffusion");;
          array<double, 3> gradient_0_;
          array<double, 3> gradient_1_;
          array<double, 3> gradient_2_;

          dg_0_->GetGradient(position, &gradient_0_);
          double concentration_0 = dg_0_->GetConcentration(position);
          dg_1_->GetGradient(position, &gradient_1_);
          double concentration_1 = dg_1_->GetConcentration(position);
          dg_2_->GetGradient(position, &gradient_2_);
          double concentration_2 = dg_2_->GetConcentration(position);

          map<int, double> concentrationMap;
          concentrationMap[0] = concentration_0;
          concentrationMap[1] = concentration_1;
          concentrationMap[2] = concentration_2;

          vector<int> possibleCellType;
          int nbOfZero = 0;
          double smallestValue = 1e10;
          int smallestConcentrationType = -1;

          for (auto it = concentrationMap.begin(); it != concentrationMap.end();
               ++it) {
            if (it->second == 0) {
              possibleCellType.push_back(it->first);
              smallestConcentrationType = it->first;
              nbOfZero++;
            }
            if (it->second < smallestValue) {
              smallestValue = it->second;
              smallestConcentrationType = it->first;
            }
          }

          if (nbOfZero < 2) {
            cell->SetCellType(smallestConcentrationType);
          } else {
            cell->SetCellType(
                possibleCellType[random->Uniform(0, possibleCellType.size())]);
          }
        }  // end cell type = -1

        /* -- internal clock -- */
        // probability to increase internal clock
        if (random->Uniform(0, 1) < 0.96) {
          // update cell internal clock
          cell->SetInternalClock(cell->GetInternalClock() + 1);
        } // end update cell internal clock
      }
    }  // end Run()

  }; // end biologyModule RGC_mosaic_BM


  // Define cell behavior for substance secretion
  struct Substance_secretion_BM : public BaseBiologyModule {
    BDM_STATELESS_BM_HEADER(Substance_secretion_BM, BaseBiologyModule, 1);

    public:
      Substance_secretion_BM() : BaseBiologyModule(gNullEventId) {}

    void Run(SimObject* so) override {
      if (auto* cell = so->As<MyCell>()) {
        auto* sim = Simulation::GetActive();
        auto* rm = sim->GetResourceManager();
        DiffusionGrid* dg = nullptr;

        if (cell->GetInternalClock()%2==0) {
          auto& secretion_position = cell->GetPosition();
          if (cell->GetCellType() == 0) {
            dg = rm->GetDiffusionGrid("on_diffusion");
          }
          else if (cell->GetCellType() == 1) {
            dg = rm->GetDiffusionGrid("off_diffusion");
          }
          else if (cell->GetCellType() == 2) {
            dg = rm->GetDiffusionGrid("on_off_diffusion");
          }
          // if dg has been correctly initialised
          if (dg != nullptr) {
            dg->IncreaseConcentrationBy(secretion_position, 1);
          }
        }
      }
    } // end Run()

  }; // end biologyModule Substance_secretion_BM


  // Define cell behavior for neurite creation
  struct Neurite_creation_BM: public BaseBiologyModule {
    BDM_STATELESS_BM_HEADER(Neurite_creation_BM, BaseBiologyModule, 1);

  public:
    Neurite_creation_BM() : BaseBiologyModule(gNullEventId) {}

    void Run(SimObject* so) override {
      if (auto* soma = so->As<MyCell>()) {
        auto* sim = Simulation::GetActive();
        // if we want to create dendrites after mosaics formation
        bool createDendrites = true;

        // if mosaic formation is over and cell got a type
        if (createDendrites && soma->GetInternalClock() == 1000
            && soma->GetCellType() != -1) {
          auto* random = sim->GetRandom();
          // dendrite per cell: average=4.5; std=1.2
          int thisSubType = soma->GetCellType()*100 + (int)random->Uniform(0, 20);
          for (int i = 0; i <= (int)random->Uniform(2, 7); i++) {
            auto&& ne = soma->ExtendNewNeurite({0, 0, 1});
            ne->AddBiologyModule(new RGC_dendrite_growth_BM());
            ne->SetHasToRetract(false);
            ne->SetSleepMode(false);
            ne->SetBeyondThreshold(false);
            ne->SetSubtype(thisSubType);
            ne->SetMySoma(soma->GetSoPtr());
          }
          // remove BM that are not needed anymore
          // soma->RemoveBiologyModule(soma->template GetBiologyModules<RGC_mosaic_BM>()[0]);
  //        soma->RemoveBiologyModule(soma->template GetBiologyModules<Neurite_creation_BM>()[0]);
        }
      }
    } // end run

  }; // endNeurite_creation_BM

} // end namespace bdm

#endif
