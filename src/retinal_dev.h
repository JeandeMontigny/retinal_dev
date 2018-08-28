#ifndef RETINAL_DEV_
#define RETINAL_DEV_

#include "biodynamo.h"
#include "random.h"
#include "substance_initializers.h"
#include "model_initializer.h"

#include "neuroscience/compile_time_param.h"
#include "neuroscience/neuron_soma.h"
#include "neuroscience/neurite_element.h"

namespace bdm {
// cell type: 0=on; 1=off; -1=not differentiated yet
// std::cout << setprecision (15) <<

// TODO: add external biology modules

/* TODO
- other comparison needed. This is not just one mosaic, but two. have to take
that into account. RI doesn't necessary reflect the mosaic (ie two mosaic
without overlaping would have high R) but just the regularity. check for cells
pairs (pair of diff types). cell fate + movement should have a better result
than movement or cell fate alone.
*/

using namespace std;

  enum Substances {
    on_diffusion, off_diffusion,
    on_substance_RGC_guide, off_substance_RGC_guide
  };

// Define my custom cell MyCell extending NeuronSoma
BDM_SIM_OBJECT(MyCell, experimental::neuroscience::NeuronSoma) {
  BDM_SIM_OBJECT_HEADER(MyCellExt, 1, cell_type_, internal_clock_);

 public:
  MyCellExt() {}

  MyCellExt(const array<double, 3>& position) : Base(position) {}

  void SetCellType(int t) { cell_type_[kIdx] = t; }
  int GetCellType() const { return cell_type_[kIdx]; }
  // This function is used by ParaView for coloring the cells by their type
  int* GetCellTypePtr() { return cell_type_.data(); }

  void SetInternalClock(int t) { internal_clock_[kIdx] = t; }
  int GetInternalClock() const { return internal_clock_[kIdx]; }

 private:
  vec<int> cell_type_;
  vec<int> internal_clock_;
};

// Define my custom neurite MyNeurite, which extends NeuriteElement
BDM_SIM_OBJECT(MyNeurite, experimental::neuroscience::NeuriteElement) {
  BDM_SIM_OBJECT_HEADER(MyNeuriteExt, 1, has_to_retract_, beyond_threshold_, sleep_mode_, diam_before_retract_);

public:
  MyNeuriteExt() {}
  MyNeuriteExt(const array<double, 3>& position) : Base(position) {}

  void SetHasToRetract(int r) { has_to_retract_[kIdx] = r; }
  bool GetHasToRetract() const { return has_to_retract_[kIdx]; }

  void SetBeyondThreshold(int r) { beyond_threshold_[kIdx] = r; }
  bool GetBeyondThreshold() const { return beyond_threshold_[kIdx]; }

  void SetSleepMode(int r) { sleep_mode_[kIdx] = r; }
  bool GetSleepMode() const { return sleep_mode_[kIdx]; }

  void SetDiamBeforeRetraction(double d) { diam_before_retract_[kIdx] = d; }
  double GetDiamBeforeRetraction() const { return diam_before_retract_[kIdx]; }

 private:
  vec<bool> has_to_retract_;
  vec<bool> beyond_threshold_;
  vec<bool> sleep_mode_;
  vec<int> diam_before_retract_;
};


template <typename TSimulation = Simulation<>>
struct RGC_dendrite_growth_test : public BaseBiologyModule {
  RGC_dendrite_growth_test() : BaseBiologyModule(gAllBmEvents) {}

  template <typename T>
  void Run(T* sim_object) {
    auto* sim = TSimulation::GetActive();
    auto* random = sim->GetRandom();
//    auto* param = sim->GetParam();
    auto* rm = sim->GetResourceManager();

    if (sim_object->template IsSoType<MyNeurite>()) {
      auto&& sim_objectNe = sim_object->template ReinterpretCast<MyNeurite>();
      auto ne = sim_objectNe.GetSoPtr();

      if (ne->IsTerminal() && ne->GetDiameter() >= 0.6) {
        if (!init_) {
          dg_on_RGCguide_ = rm->GetDiffusionGrid("on_substance_RGC_guide");
          dg_off_RGCguide_ = rm->GetDiffusionGrid("off_substance_RGC_guide");
          init_ = true;
        }
        array<double, 3> gradient_RGCguide;
        double concentration = 0;
        // initialise the correct substance as guide depending on cell type
        if (ne->GetNeuronSomaOfNeurite()->GetCellType() == 0) {
          dg_on_RGCguide_->GetGradient(ne->GetPosition(), &gradient_RGCguide);
          concentration = dg_on_RGCguide_->GetConcentration(ne->GetPosition());
        }
        if (ne->GetNeuronSomaOfNeurite()->GetCellType() == 1) {
          dg_off_RGCguide_->GetGradient(ne->GetPosition(), &gradient_RGCguide);
          concentration = dg_off_RGCguide_->GetConcentration(ne->GetPosition());
        }

        double gradientWeight = 0.5;
        double randomnessWeight = 0.6;
        double oldDirectionWeight = 2.0;
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
        ne->SetDiameter(ne->GetDiameter()-0.0005);

        if (concentration > 0.04 && random->Uniform() < 0.01) {
          ne->SetDiameter(ne->GetDiameter()-0.01);
          ne->Bifurcate();
        }

      } // if is terminal
    } // if MyNeurite object
  } // end run

  private:
    bool init_ = false;
    DiffusionGrid* dg_on_RGCguide_ = nullptr;
    DiffusionGrid* dg_off_RGCguide_ = nullptr;
    ClassDefNV(RGC_dendrite_growth_test, 1);
  };


// Define dendrites behavior for RGC dendritic growth
template <typename TSimulation = Simulation<>>
struct RGC_dendrite_growth : public BaseBiologyModule {
  RGC_dendrite_growth() : BaseBiologyModule(gAllBmEvents) {}

  template <typename T>
  void Run(T* sim_object) {
    auto* sim = TSimulation::GetActive();
    auto* random = sim->GetRandom();
    auto* rm = sim->GetResourceManager();

    if (sim_object->template IsSoType<MyNeurite>()) {
      auto&& sim_objectNe = sim_object->template ReinterpretCast<MyNeurite>();

      if (!init_) {
        dg_on_RGCguide_ = rm->GetDiffusionGrid("on_substance_RGC_guide");
        dg_off_RGCguide_ = rm->GetDiffusionGrid("off_substance_RGC_guide");
        init_ = true;
      }

      double diamReducSpeed = 0.00052;
      double branchingFactor = 0.02; // 0.002
      // inline double getBranchingFactor(NeuriteElement ne) {
      //   return branchingFactor(4.5/ne->GetDiameter()*ne->GetDiameter());
      // }

      double threshold_1 = 0.8;
      double threshold_2 = 0.6;
      double growthSpeed = 10;
      array<double, 3> gradient_RGCguide;
      double concentration = 0;

      auto ne = sim_objectNe.GetSoPtr();
      // if neurite is not in sleep & is terminal & mosaic is over
      if (!ne->GetSleepMode() && ne->IsTerminal() &&
          ne->GetNeuronSomaOfNeurite()->GetCellType() != -1 &&
          ne->GetNeuronSomaOfNeurite()->GetInternalClock() > 400) {

        // if neurite does not have to retract
        if (!ne->GetHasToRetract()) {
          auto& positionNeurite = ne->GetPosition();

          double oldDirectionWeight = 0.6; // cylinder axis
          double gradientWeight = 0.8;
          double randomnessWeight = 0.6;

          // initialise the correct substance as guide depending on cell type
          if (ne->GetNeuronSomaOfNeurite()->GetCellType() == 0) {
            dg_on_RGCguide_->GetGradient(positionNeurite, &gradient_RGCguide);
            concentration = dg_on_RGCguide_->GetConcentration(positionNeurite);
          }
          if (ne->GetNeuronSomaOfNeurite()->GetCellType() == 1) {
            dg_off_RGCguide_->GetGradient(positionNeurite, &gradient_RGCguide);
            concentration = dg_off_RGCguide_->GetConcentration(positionNeurite);
          }

          // branching factor should depend on diameter
          if (ne->GetDiameter() > 0.6 &&
              random->Uniform(0, 1) < branchingFactor) {
            ne->SetDiameter(ne->GetDiameter() - 0.0016);
            ne->Bifurcate();
            // auto ne2 = ne->Bifurcate()[1];
            // ne2->AddBiologyModule(RGC_dendrite_growth<>());
          }

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

          // ne->SetDiameter(ne->GetDiameter()-diamReducSpeed);
          // TODO: if elongation create new cylinder, transfer biologyModule
          // module to that new cylinder. make sure to deleted old one.
          ne->ElongateTerminalEnd(growthSpeed, newStepDirection);
          ne->SetDiameter(ne->GetDiameter()-diamReducSpeed);

          // homo-type interaction
          // lambda updating external counters for neighbor neurites
          int ownType = 0;
          int otherType = 0;
          // auto countNeighbours = [&](auto&& neighbor, SoHandle neighbor_handle) {
          //   // if neighbor is a NeuriteElement
          //   if (neighbor->GetSoPtr()->IsSoType(ne)) {
          //     auto&& neighbor_rc = neighbor->template
          //       ReinterpretCast<MyNeurite>();
          //     auto n_soptr = neighbor_rc->GetSoPtr();
          //     // if it is a direct relative
          //     if (n_soptr->GetNeuronSomaOfNeurite() !=
          //         ne->GetNeuronSomaOfNeurite() &&
          //         n_soptr->GetNeuronSomaOfNeurite()->GetCellType() ==
          //         ne->GetNeuronSomaOfNeurite()->GetCellType()) {
          //       ownType++;
          //     }
          //     else {
          //       otherType++;
          //     }
          //   }
          // };
          // auto* grid = sim->GetGrid();
          // grid->ForEachNeighborWithinRadius(
          //   countNeighbours, ne, ne->GetSoHandle(), 4);

          if (ownType > otherType) {
            ne->SetHasToRetract(true);
            ne->SetDiamBeforeRetraction(ne->GetDiameter());
          }

          if (concentration < threshold_2 && ne->GetDiameter() < 0.8) {
            ne->SetHasToRetract(true);
            ne->SetBeyondThreshold(true);
            ne->SetDiamBeforeRetraction(ne->GetDiameter());
          }
        }  // end not has to retract

        // if has to retract
        else {
//          ne->RetractTerminalEnd(40);

          // what to do in case of retraction due to interaction
          if (!ne->GetBeyondThreshold() &&
              ne->GetDiamBeforeRetraction() > ne->GetDiameter() + 0.15) {
            ne->SetHasToRetract(false);
            ne->SetSleepMode(true);
            // TODO: remove biologyModule for this neurite
            // ne->RemoveLocalBiologyModule(RGC_dendrite_growth<>());
          }

          // if retraction due to threshold
          if (ne->GetBeyondThreshold()) {
            // initialise the correct substance as guide depending on cell type
            auto& positionNeurite = ne->GetPosition();
            if (ne->GetNeuronSomaOfNeurite()->GetCellType() == 0) {
              dg_on_RGCguide_->GetGradient(positionNeurite, &gradient_RGCguide);
              concentration =
                  dg_on_RGCguide_->GetConcentration(positionNeurite);
            }
            if (ne->GetNeuronSomaOfNeurite()->GetCellType() == 1) {
              dg_off_RGCguide_->GetGradient(positionNeurite,
                                            &gradient_RGCguide);
              concentration =
                  dg_off_RGCguide_->GetConcentration(positionNeurite);
            }
            if (concentration > threshold_1) {
              ne->SetBeyondThreshold(false);
              ne->SetHasToRetract(false);
              ne->ElongateTerminalEnd(1, {random->Uniform(-1, 1),
                random->Uniform(-1, 1), random->Uniform(-1, 1)});
            }
          } // end if beyound threshold
        } // end has to retract

      } // end if neurite not in sleep & cell is not neutral type
    } // end if is MyNeurite element
  } // end run

 private:
  bool init_ = false;
  DiffusionGrid* dg_on_RGCguide_ = nullptr;
  DiffusionGrid* dg_off_RGCguide_ = nullptr;
  ClassDefNV(RGC_dendrite_growth, 1);
}; // end biologyModule RGC_dendrite_growth


// Define cell behavior for mosaic formation
template <typename TSimulation = Simulation<>>
struct Chemotaxis : public BaseBiologyModule {
  Chemotaxis() : BaseBiologyModule(gAllBmEvents) {}

  template <typename T>
  void Run(T* sim_object) {
    auto* sim = TSimulation::GetActive();
    auto* random = sim->GetRandom();

    if (sim_object->template IsSoType<MyCell>()) {
      auto&& cell = sim_object->template ReinterpretCast<MyCell>();

      bool withCellDeath = false; // run cell death mechanism
      bool withMovement = false; // run tangential migration

      // if not initialised, initialise substance diffusions
      if (!init_) {
        auto* rm = sim->GetResourceManager();
        dg_0_ = rm->GetDiffusionGrid("on_diffusion");
        dg_1_ = rm->GetDiffusionGrid("off_diffusion");
        init_ = true;
      }

      auto& position = cell->GetPosition();
      array<double, 3> gradient_0_;
      array<double, 3> gradient_1_;
      array<double, 3> diff_gradient;
      array<double, 3> gradient_z;
      double concentration = 0;
      int cellClock = cell->GetInternalClock();

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

      // add small random movements
      // cell->UpdatePosition({random->Uniform(-0.1, 0.1),
      //                       random->Uniform(-0.1, 0.1), 0});
      // cell growth
      if (cell->GetDiameter() < 15 && random->Uniform(0, 1) < 0.01) {
        cell->ChangeVolume(0.1);
      }

      /* -- cell movement -- */
      if (withMovement && cellClock >= 400) {
        // cell movement based on homotype substance gradient
        // 0. with cell death - 0. without
        // 10100755032 too much - 1010075504 not enough
        if (concentration >= 0.10100755034) {
          cell->UpdatePosition(diff_gradient);
        }
      } // end tangential migration

      /* -- cell death -- */
      if (withCellDeath && cellClock >= 400 && cellClock < 2000) {
        // add vertical migration as the multi layer colapse in just on layer
        cell->UpdatePosition(gradient_z);
        // cell death depending on homotype substance concentration
        // with cell fate: 0.11
        // with cell movement:
        // with cell fate and cell movement:
        if (concentration >= 0.11 && random->Uniform(0, 1) < 0.1) {
          cell->RemoveFromSimulation();
        }

        auto killNeighbor = [&](auto&& neighbor, SoHandle neighbor_handle) {
          if (random->Uniform(0, 1) < 0.2) {
            neighbor->RemoveFromSimulation();
          }
        };
        auto* grid = sim->GetGrid();
        grid->ForEachNeighborWithinRadius(
          killNeighbor, cell, cell->GetSoHandle(), cell->GetDiameter()*1.5);
      } // end cell death

      /* -- cell fate -- */
      // cell type attribution depending on concentrations
      if (cell->GetCellType() == -1) {  // if cell type is not on or off
        dg_1_->GetGradient(position, &gradient_1_);
        double concentration_1 = dg_1_->GetConcentration(position);
        dg_0_->GetGradient(position, &gradient_0_);
        double concentration_0 = dg_0_->GetConcentration(position);
        // if no substances
        // random so all cell types doesn't create all randomly at step 1
        if (concentration_1 == 0 &&
          concentration_0 == 0 && random->Uniform(0, 1) < 0.0001) {
          // random attribution of a cell type
          cell->SetCellType((int)random->Uniform(0, 2));
        }
        // if [Off substance] > [On substance] - 1e-8 for diffusion = 0.02
        // random so all cell doesn't choose their type in the same time
        if (concentration_1 > 1e-8 && concentration_1 > concentration_0
          && random->Uniform(0, 1) < 0.1) {
            cell->SetCellType(0); // cell become on
        }
        // if [On substance] > [Off substance]
        // random so cells don't choose their type in the same time
        if (concentration_0 > 1e-8 && concentration_0 > concentration_1
          && random->Uniform(0, 1) < 0.1) {
            cell->SetCellType(1); // cell become off
        }
      }  // end cell type = -1

      // probability to increase internal clock
      if (random->Uniform(0, 1) < 0.96) {
        // update cell internal clock
        cell->SetInternalClock(cell->GetInternalClock() + 1);
      } // end update cell internal clock

    }  // end of if neuron soma
  }    // end Run()

 private:
  bool init_ = false;
  DiffusionGrid* dg_0_ = nullptr;
  DiffusionGrid* dg_1_ = nullptr;
  ClassDefNV(Chemotaxis, 1);
};  // end biologyModule Chemotaxis


// 1b. Define secretion behavior:
template <typename TSimulation = Simulation<>>
struct SubstanceSecretion : public BaseBiologyModule {
  // Daughter cells inherit this biology module
  SubstanceSecretion() : BaseBiologyModule(gAllBmEvents) {}

  template <typename T>
  void Run(T* sim_object) {
    auto* sim = TSimulation::GetActive();

    if (sim_object->template IsSoType<MyCell>()) {
      auto&& cell = sim_object->template ReinterpretCast<MyCell>();

      if (!init_) {
        auto* rm = sim->GetResourceManager();
        dg_0_ = rm->GetDiffusionGrid("on_diffusion");
        dg_1_ = rm->GetDiffusionGrid("off_diffusion");
        init_ = true;
      }
      auto& secretion_position = cell->GetPosition();
      // if off cell, secrete off cells substance
      if (cell->GetCellType() == 1) {
        dg_1_->IncreaseConcentrationBy(secretion_position, 0.1);
      }
      // is on cell, secrete on cells substance
      else if (cell->GetCellType() == 0) {
        dg_0_->IncreaseConcentrationBy(secretion_position, 0.1);
      }
    }
  }

 private:
  bool init_ = false;
  DiffusionGrid* dg_0_ = nullptr;
  DiffusionGrid* dg_1_ = nullptr;
  ClassDefNV(SubstanceSecretion, 1);
};  // end biologyModule SubstanceSecretion


// define compile time parameter
template <typename Backend>
struct CompileTimeParam : public DefaultCompileTimeParam<Backend> {
  using BiologyModules =
      Variant<Chemotaxis<>, SubstanceSecretion<>, RGC_dendrite_growth<>, RGC_dendrite_growth_test<>>;
  using AtomicTypes = VariadicTypedef<MyCell, MyNeurite>;
  using NeuronSoma = MyCell;
  using NeuriteElement = MyNeurite;
};

// define my cell creator
template <typename TSimulation = Simulation<>>
static void CellCreator(double min, double max, int num_cells, int cellType) {
  auto* sim = TSimulation::GetActive();
  auto* rm = sim->GetResourceManager();
  auto* random = sim->GetRandom();

  auto* container = rm->template Get<MyCell>();
  container->reserve(num_cells);

  for (int i = 0; i < num_cells; i++) {
    double x = random->Uniform(min + 20, max - 20);
    double y = random->Uniform(min + 20, max - 20);
    double z = random->Uniform(min + 20, 40);  // 24
    std::array<double, 3> position = {x, y, z};

    auto&& cell = rm->template New<MyCell>(position);
    cell.SetDiameter(random->Uniform(7, 8)); // random diameter
    cell.SetCellType(cellType);
    cell.SetInternalClock(0);
    cell.AddBiologyModule(SubstanceSecretion<>());
    cell.AddBiologyModule(Chemotaxis<>());

    auto&& ne = cell.ExtendNewNeurite({0, 0, 1});
    ne->GetSoPtr()->AddBiologyModule(RGC_dendrite_growth<>());
    ne->GetSoPtr()->SetHasToRetract(false);
    ne->GetSoPtr()->SetSleepMode(false);
    ne->GetSoPtr()->SetBeyondThreshold(false);

    // container->push_back(new_simulation_object);
  }
  container->Commit();
}  // end CellCreator

// position exporteur
template <typename TSimulation = Simulation<>>
inline void position_exporteur(int i) {
  int seed = TSimulation::GetActive()->GetRandom()->GetSeed();
  ofstream positionFileOn;
  ofstream positionFileOff;
  ofstream positionFileAll;
  stringstream sstmOn;
  stringstream sstmOff;
  stringstream sstmAll;
  sstmOn << "on_t" << i << "_seed" << seed << ".txt";
  sstmOff << "off_t" << i << "_seed" << seed << ".txt";
  sstmAll << "all_t" << i << "_seed" << seed << ".txt";

  string fileNameOn = sstmOn.str();
  string fileNameOff = sstmOff.str();
  string fileNameAll = sstmAll.str();
  positionFileOn.open(fileNameOn);
  positionFileOff.open(fileNameOff);
  positionFileAll.open(fileNameAll);

  auto* sim = TSimulation::GetActive();
  auto* rm = sim->GetResourceManager();
  auto my_cells = rm->template Get<MyCell>();
  int numberOfCells = my_cells->size();

  for (int cellNum = 0; cellNum < numberOfCells; cellNum++) {
    auto thisCellType = (*my_cells)[cellNum].GetCellType();
    auto position = (*my_cells)[cellNum].GetPosition();
    if (thisCellType == 0) {
      positionFileOn << position[0] << " " << position[1] << "\n";
      positionFileAll << position[0] << " " << position[1] << " "
                      << position[2] << " on\n";
    }
    else if (thisCellType == 1) {
      positionFileOff << position[0] << " " << position[1] << "\n";
      positionFileAll << position[0] << " " << position[1] << " "
                      << position[2] << " off\n";
    }
    else {
    positionFileAll << position[0] << " " << position[1] << " "
                    << position[2] << " nd\n";
    }
  }
  positionFileOn.close();
  positionFileOff.close();
  positionFileAll.close();
}  // end position_exporteur

// RI computation
inline double computeRI(vector<array<double, 3>> coordList) {
  vector<double> shortestDistList;
  for (unsigned int i = 0; i < coordList.size();
       i++) {  // for each cell of same type in the simulation
    array<double, 3> cellPosition = coordList[i];

    double shortestDist = 9999;
    // for each other cell of same type in the simulation
    for (unsigned int j = 0; j < coordList.size();  j++) {
      array<double, 3> otherCellPosition = coordList[j];

      // get the distance between those 2 cells (x-y plan only)
      double tempsDistance =
          sqrt(pow(cellPosition[0] - otherCellPosition[0], 2) +
               pow(cellPosition[1] - otherCellPosition[1], 2));
      // if cell is closer and is not itself
      if (tempsDistance < shortestDist && tempsDistance != 0) {
        shortestDist = tempsDistance; // updade closest neighbour distance
      }
    }
    shortestDistList.push_back(
        shortestDist);  // save the shortest distance to neighbour
  }
  // compute mean
  double temps_sum = 0;
  for (unsigned int i = 0; i < shortestDistList.size(); i++) {
    temps_sum += shortestDistList[i];
  }
  double aveShortestDist = temps_sum / (double)shortestDistList.size();
  // compute std
  double temp = 0;
  for (unsigned int i = 0; i < shortestDistList.size(); i++) {
    double val = shortestDistList[i];
    double squrDiffToMean = pow(val - aveShortestDist, 2);
    temp += squrDiffToMean;
  }
  double meanOfDiffs = temp / (double)(shortestDistList.size());
  double std = sqrt(meanOfDiffs);

  return aveShortestDist / std;  // return RI
}  // end computeRI

// get cells coordinate of same cell_type_ to call computeRI
template <typename TSimulation = Simulation<>>
inline double getRI(int desiredCellType) {
  auto* sim = TSimulation::GetActive();
  auto* rm = sim->GetResourceManager();
  auto my_cells = rm->template Get<MyCell>(); // get cell list
  vector<array<double, 3>> coordList; // list of coordinate
  int numberOfCells = my_cells->size();
  // for each cell in simulation
  for (int cellNum = 0; cellNum < numberOfCells; cellNum++) {
    auto thisCellType = (*my_cells)[cellNum].GetCellType();
    if (thisCellType == desiredCellType) { // if cell is of the desired type
      auto position = (*my_cells)[cellNum].GetPosition(); // get its position
      coordList.push_back(position); // put cell coord in the list
    }
  }
  return computeRI(coordList);  // return RI for desired cell type
}  // end getRI

/* -------- simulate -------- */
template <typename TSimulation = Simulation<>>
inline int Simulate(int argc, const char** argv) {
  Simulation<> simulation(argc, argv);
  auto* rm = simulation.GetResourceManager();
  auto* random = simulation.GetRandom();
  auto* scheduler = simulation.GetScheduler();
  auto* param = simulation.GetParam();

  // number of simulation steps
  int maxStep = 2500;
  // Create an artificial bounds for the simulation space
  int cubeDim = 500;
  int num_cells = 4400;
  double cellDensity = (double)num_cells * 1e6 / (cubeDim * cubeDim);
  cout << "cell density: " << cellDensity << " cells per cm2" << endl;

  // cell are created with +20 to min and -20 to max
  param->bound_space_ = true;
  param->min_bound_ = 0;
  param->max_bound_ = cubeDim + 40;
  // set min and max length for neurite segments
  param->kNeuriteMinLength = 1.0;
  param->kNeuriteMaxLength = 2.0;
//  param->kOutputDir = "myOutput";

  int mySeed = rand() % 10000;
  mySeed = 9784;  // 9784
  random->SetSeed(mySeed);
  cout << "modelling with seed " << mySeed << endl;

  // min position, max position, number of cells , cell type
  CellCreator(param->min_bound_, param->max_bound_, num_cells/2, 0);
  cout << "on cells created" << endl;
  CellCreator(param->min_bound_, param->max_bound_, num_cells/2, 1);
  cout << "off cells created" << endl;
  CellCreator(param->min_bound_, param->max_bound_, 0, -1); // num_cells
  cout << "undifferentiated cells created" << endl;

  // 3. Define substances
  // Order: substance_name, diffusion_coefficient, decay_constant, resolution
  // if diffusion_coefficient is low, diffusion distance is short
  // if decay_constant is high, diffusion distance is short
  // resolution is number of point in one domaine dimension
  ModelInitializer::DefineSubstance(0, "on_diffusion",
                                    0.05, 0.99, param->max_bound_/10);
  ModelInitializer::DefineSubstance(1, "off_diffusion",
                                    0.05, 0.99, param->max_bound_/10);

  // define substances for RGC dendrite attraction
  ModelInitializer::DefineSubstance(2, "on_substance_RGC_guide",
                                    0, 0, param->max_bound_/5);
  ModelInitializer::DefineSubstance(3, "off_substance_RGC_guide",
                                    0, 0, param->max_bound_/5);
  // create substance gaussian distribution for RGC dendrite attraction
  // average peak distance for ON cells: 15.959 with std of 5.297;
  ModelInitializer::InitializeSubstance(2, "on_substance_RGC_guide",
                                        GaussianBand(46, 6, Axis::kZAxis));
  // average peak distance for OFF cells: 40.405 with std of 8.39;
  ModelInitializer::InitializeSubstance(3, "off_substance_RGC_guide",
                                        GaussianBand(70, 8, Axis::kZAxis));
  cout << "substances initialised" << endl;

  // set write output param
  // if you want to write file for RI and cell position
  bool writeRI = false;
  bool writePositionExport = false;
  // create cell position files every outputFrequence steps
  int outputFrequence = 100;
  ofstream outputFile;

  if (writeRI) {
    outputFile.open("RI_" + to_string(mySeed) + ".txt");
  }

  // 4. Run simulation for maxStep timesteps
  auto my_cells = rm->template Get<MyCell>();
  int numberOfCells = my_cells->size();
  int numberOfCells0 = 0;
  int numberOfCells1 = 0;

  for (int i = 0; i < maxStep; i++) {
    scheduler->Simulate(1);

    if (i % 10 == 0) {  // write RI in file
      double RIon = getRI(0);
      double RIoff = getRI(1);
      //        cout << "RI on: " << RIon << " ; RI off: " << RIoff << endl;
      if (writeRI) {
        outputFile << RIon << " " << RIoff << "\n";
      }

      if (i % 100 == 0) {  // print
        // get cell list size
        rm = simulation.GetResourceManager();
        my_cells = rm->template Get<MyCell>();
        numberOfCells = my_cells->size();
        numberOfCells0 = 0;
        numberOfCells1 = 0;

        for (int cellNum = 0; cellNum < numberOfCells;
             cellNum++) {  // for each cell in simulation
          auto thisCellType = (*my_cells)[cellNum].GetCellType();
          if (thisCellType == 0) {
            numberOfCells0++;
          } else if (thisCellType == 1) {
            numberOfCells1++;
          }
        }
        cout << "step " << i << " out of " << maxStep << "\n"
             << numberOfCells << " cells in simulation: "
             << (1 - ((double)numberOfCells / num_cells)) * 100
             << "% of cell death\n"
             << numberOfCells0 << " cells are type 0 (on) ; "
             << numberOfCells1 << " cells are type 1 (off)\n"
             << "RI on: " << RIon << " ; RI off: " << RIoff
             << " ; mean: " << (RIon + RIoff) / 2 << endl;
      } // end every 100 simu steps
    } // end every 10 simu steps

    if (writePositionExport && i % outputFrequence == 0) {
      position_exporteur(i); // export cell position
    }

  } // end for Simulate

  outputFile.close();

  return 0;
}

}  // namespace bdm

#endif // RETINAL_DEV_
