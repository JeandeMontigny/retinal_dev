#ifndef RGC_SOMA_BM_
#define RGC_SOMA_BM_

#include "rgc_dendrite_bm.h"
#include "extended_objects.h"

namespace bdm {
using namespace std;

  // Define cell behavior for mosaic formation
  struct RGC_mosaic_BM : public BaseBiologyModule {
    RGC_mosaic_BM() : BaseBiologyModule(gNullEventId) {}

    // Default event constructor
    template <typename TEvent, typename TBm>
    RGC_mosaic_BM(const TEvent& event, TBm* other, uint64_t new_oid = 0) {
    }

    // Default event handler
    template <typename TEvent, typename... TBms>
    void EventHandler(const TEvent&, TBms*...) {
    }

    template <typename T, typename TSimulation = Simulation<>>
    void Run(T* cell) {
      auto* sim = TSimulation::GetActive();
      auto* random = sim->GetRandom();

      // run tangential migration
      bool withMovement = true;
      double movementThreshold = 1.21;
      bool withDeath = true;
      double deathThreshold = 1.216;

      // if not initialised, initialise substance diffusions
      if (!init_) {
        auto* rm = sim->GetResourceManager();
        dg_0_ = rm->GetDiffusionGrid("on_diffusion");
        dg_1_ = rm->GetDiffusionGrid("off_diffusion");
        dg_2_ = rm->GetDiffusionGrid("on_off_diffusion");
        init_ = true;
      }

      auto& position = cell->GetPosition();
      array<double, 3> gradient_0_;
      array<double, 3> gradient_1_;
      array<double, 3> gradient_2_;
      array<double, 3> diff_gradient;
      array<double, 3> gradient_z;
      double concentration = 0;
      int cellClock = cell->GetInternalClock();

      // if cell is type 0, concentration and gradient are substance 0
      if (cell->GetCellType() == 0) {
        dg_0_->GetGradient(position, &gradient_0_);
        gradient_z = Math::ScalarMult(0.2, gradient_0_);
        gradient_z[0] = 0;
        gradient_z[1] = 0;
        diff_gradient = Math::ScalarMult(-0.1, gradient_0_);
        diff_gradient[2] = 0;
        concentration = dg_0_->GetConcentration(position);
      }
      // if cell is type 1, concentration and gradient are substance 1
      if (cell->GetCellType() == 1) {
        dg_1_->GetGradient(position, &gradient_1_);
        gradient_z = Math::ScalarMult(0.2, gradient_1_);
        gradient_z[0] = 0;
        gradient_z[1] = 0;
        diff_gradient = Math::ScalarMult(-0.1, gradient_1_);
        diff_gradient[2] = 0;
        concentration = dg_1_->GetConcentration(position);
      }
      // if cell is type 2, concentration and gradient are substance 2
      if (cell->GetCellType() == 2) {
        dg_2_->GetGradient(position, &gradient_2_);
        gradient_z = Math::ScalarMult(0.2, gradient_2_);
        gradient_z[0] = 0;
        gradient_z[1] = 0;
        diff_gradient = Math::ScalarMult(-0.1, gradient_2_);
        diff_gradient[2] = 0;
        concentration = dg_2_->GetConcentration(position);
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
      if (withDeath && cellClock >= 100 && cellClock < 700) {
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

    }  // end Run()

   private:
    bool init_ = false;
    DiffusionGrid* dg_0_ = nullptr;
    DiffusionGrid* dg_1_ = nullptr;
    DiffusionGrid* dg_2_ = nullptr;
    ClassDefNV(RGC_mosaic_BM, 1);
  }; // end biologyModule RGC_mosaic_BM


  // Define cell behavior for substance secretion
  struct Substance_secretion_BM : public BaseBiologyModule {
    // Daughter cells inherit this biology module
    Substance_secretion_BM() : BaseBiologyModule(gNullEventId) {}

    /// Default event constructor
    template <typename TEvent, typename TBm>
    Substance_secretion_BM(const TEvent& event, TBm* other, uint64_t new_oid = 0) {
    }

    /// Default event handler (exising biology module won't be modified on
    /// any event)
    template <typename TEvent, typename... TBms>
    void EventHandler(const TEvent&, TBms*...) {
    }

    template <typename T, typename TSimulation = Simulation<>>
    void Run(T* cell) {
      auto* sim = TSimulation::GetActive();

      if (!init_) {
        auto* rm = sim->GetResourceManager();
        dg_0_ = rm->GetDiffusionGrid("on_diffusion");
        dg_1_ = rm->GetDiffusionGrid("off_diffusion");
        dg_2_ = rm->GetDiffusionGrid("on_off_diffusion");
        init_ = true;
      }

      if (cell->GetInternalClock()%2==0) {
        auto& secretion_position = cell->GetPosition();
        // if on cell, secrete on cells substance
        if (cell->GetCellType() == 0) {
          dg_0_->IncreaseConcentrationBy(secretion_position, 1);
        }
        // is off cell, secrete off cells substance
        else if (cell->GetCellType() == 1) {
          dg_1_->IncreaseConcentrationBy(secretion_position, 1);
        }
        // is on-off cell, secrete on-off cells substance
        else if (cell->GetCellType() == 2) {
          dg_2_->IncreaseConcentrationBy(secretion_position, 1);
        }
      }
    }

   private:
    bool init_ = false;
    DiffusionGrid* dg_0_ = nullptr;
    DiffusionGrid* dg_1_ = nullptr;
    DiffusionGrid* dg_2_ = nullptr;
    ClassDefNV(Substance_secretion_BM, 1);
  }; // end biologyModule Substance_secretion_BM


  // Define cell behavior for neurite creation
  struct Neurite_creation_BM: public BaseBiologyModule {
    Neurite_creation_BM() : BaseBiologyModule(gNullEventId) {}

    /// Default event constructor
    template <typename TEvent, typename TBm>
    Neurite_creation_BM(const TEvent& event, TBm* other, uint64_t new_oid = 0) {}

    template <typename TEvent, typename... TBms>
    void EventHandler(const TEvent&, TBms*...) {}

    template <typename T, typename TSimulation = Simulation<>>
    void Run(T* soma) {
      auto* sim = TSimulation::GetActive();
      auto* random = sim->GetRandom();

      /* -- internal clock -- */
      // probability to increase internal clock
      if (random->Uniform(0, 1) < 0.96) {
        // update cell internal clock
        soma->SetInternalClock(soma->GetInternalClock() + 1);
      } // end update cell internal clock

      if (createDendrites_ && soma->GetInternalClock() == 1 && soma->GetCellType() != -1) {

        // dendrite per cell: average=4.5; std=1.2
        int thisSubType = soma->GetCellType()*100 + (int)random->Uniform(0, 20);

        // vector< array<double, 3>> root_list  = {{0, 0, 1}, {0.2, 0, 1}, {-0.2, 0, 1},
        //   {0, 0.2, 1}, {0, -0.2, 1}, {0.2, 0.2, 1}, {-0.2, 0.2, 1}, {0.2, -0.5, 1}};

        for (int i = 0; i <= (int)random->Uniform(3, 6); i++) {
          // if all dendrites are created at the same location, setting mechanical
          // interaction to true crashes simu
          //TODO: cleaner fix
          // array<double, 3> dendrite_root = root_list[i];
          array<double, 3> dendrite_root = {0,0,1};

          auto&& ne = soma->ExtendNewNeurite(dendrite_root);
          ne->AddBiologyModule(RGC_dendrite_growth_BM());
          ne->SetHasToRetract(false);
          ne->SetSleepMode(false);
          ne->SetBeyondThreshold(false);
          ne->SetSubtype(thisSubType);
          ne->SetDiamLimit(double(random->Uniform(0.4, 0.6))); // 0.5
          ne->SetMySoma(soma->GetSoPtr());
        }
        // remove BM that are not needed anymore
        // soma->RemoveBiologyModule(soma->template GetBiologyModules<RGC_mosaic_BM>()[0]);
        // soma->RemoveBiologyModule(soma->template GetBiologyModules<Neurite_creation_BM>()[0]);
      createDendrites_ = false;
      }
      if (soma->GetInternalClock() > 100 && soma->GetDaughters().size()==0) {
	soma->RemoveFromSimulation();
      }
    } // end run

  private:
    ClassDefNV(Neurite_creation_BM, 1);
    bool createDendrites_ = true;
  }; // endNeurite_creation_BM

} // end namespace bdm

#endif
