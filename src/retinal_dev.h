#ifndef RETINAL_DEV_
#define RETINAL_DEV_

#include "biodynamo.h"
#include "random.h"
#include "substance_initializers.h"
#include "model_initializer.h"

namespace bdm {
// cell type: 0=on; 1=off; -1=not differentiated yet
// std::cout << setprecision (15) <<

/* TODO
- other comparison needed. This is not just one mosaic, but two. have to take
that into account. RI doesn't necessary reflect the mosaic (ie two mosaic
without overlaping would have high R) but just the regularity. check for cells
pairs (pair of diff types). cell fate + movement should have a better result
than movement or cell fate alone.
*/

using namespace std;

// TODO: add external biology modules: bioM_*
  enum Substances { on_substance_RGC_guide, off_substance_RGC_guide };

// Define my custom cell MyCell extending Cell
BDM_SIM_OBJECT(MyCell, Cell) {
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

      bool withCellDeath = true; // run cell death mechanism
      bool withMovement = true; // run tangential migration

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
      cell->UpdatePosition({random->Uniform(-0.1, 0.1), random->Uniform(-0.1, 0.1), 0});
      if (cell->GetDiameter() < 15 && random->Uniform(0, 1) < 0.01) {
        cell->ChangeVolume(0.1);
      }

      /* -- cell movement -- */
      if (withMovement && cellClock >= 400) {
        // cell movement based on homotype substance gradient
        // 0. with cell death - 0. without
        // 1010075 too much - 1010076 not enough
        if (concentration >= 0.1010076) {
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
        grid->ForEachNeighborWithinRadius(killNeighbor, cell, cell->GetSoHandle(), cell->GetDiameter()*1.5);
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
        if (concentration_1 == 0 && concentration_0 == 0 && random->Uniform(0, 1) < 0.0001) {
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
      Variant<Chemotaxis<>, SubstanceSecretion<>>;
  using AtomicTypes = VariadicTypedef<MyCell>;
  // using Cell = MyCell;
};

// define my cell creator
template <typename Function, typename TSimulation = Simulation<>>
static void CellCreator(double min, double max, int num_cells,
                        Function cell_builder) {
  auto* sim = TSimulation::GetActive();
  auto* rm = sim->GetResourceManager();
  auto* random = sim->GetRandom();

  // Determine simulation object type which is returned by the cell_builder
  using FunctionReturnType = decltype(cell_builder({0, 0, 0}));

  auto container = rm->template Get<FunctionReturnType>();
  container->reserve(num_cells);

  for (int i = 0; i < num_cells; i++) {
    double x = random->Uniform(min + 20, max - 20);
    double y = random->Uniform(min + 20, max - 20);
    double z = random->Uniform(min + 20, 40);  // 24
    auto new_simulation_object = cell_builder({x, y, z});
    container->push_back(new_simulation_object);
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
  int maxStep = 2000;
  // Create an artificial bounds for the simulation space
  int cubeDim = 500;
  int num_cells = 4400;
  double cellDensity = (double)num_cells * 1e6 / (cubeDim * cubeDim);
  cout << "cell density: " << cellDensity << " cells per cm2" << endl;

  param->bound_space_ = true;
  param->min_bound_ = 0;
  // cell are created with +20 to min and -20 to max. so physical cube has to be cubeDim+40
  param->max_bound_ = cubeDim + 40;
  param->run_mechanical_interactions_ = true;

  int mySeed = rand() % 10000;
  mySeed = 9784;  // 9784
  random->SetSeed(mySeed);
  cout << "modelling with seed " << mySeed << endl;

  // Construct cells of on cells (type 0)
  auto construct_on = [](const array<double, 3>& position) {
    auto* simulation = TSimulation::GetActive();
    auto* random = simulation->GetRandom();
    MyCell cell(position);
    cell.SetDiameter(random->Uniform(8, 9));  // random diameter between 8 and 9
    cell.SetCellType(0);
    cell.SetInternalClock(0);
    cell.AddBiologyModule(SubstanceSecretion<>());
    cell.AddBiologyModule(Chemotaxis<>());
    return cell;
  };
  CellCreator(param->min_bound_, param->max_bound_, 0, construct_on); // num_cells/2
  cout << "on cells created" << endl;

  // Construct cells of off cells (type 1)
  auto construct_off = [](const array<double, 3>& position) {
    auto* simulation = TSimulation::GetActive();
    auto* random = simulation->GetRandom();
    MyCell cell(position);
    cell.SetDiameter(random->Uniform(8, 9)); // random diameter
    cell.SetCellType(1);
    cell.SetInternalClock(0);
    cell.AddBiologyModule(SubstanceSecretion<>());
    cell.AddBiologyModule(Chemotaxis<>());
    return cell;
  };
  CellCreator(param->min_bound_, param->max_bound_, 0, construct_off); // num_cells/2
  cout << "off cells created" << endl;

  // construct neutral cells (type -1)
  auto construct_nonType = [](const array<double, 3>& position) {
    auto* simulation = TSimulation::GetActive();
    auto* random = simulation->GetRandom();
    MyCell cell(position);
    cell.SetDiameter(random->Uniform(8, 9));
    cell.SetCellType(-1);
    cell.SetInternalClock(0);
    cell.AddBiologyModule(SubstanceSecretion<>());
    cell.AddBiologyModule(Chemotaxis<>());
    return cell;
  };
  CellCreator(param->min_bound_, param->max_bound_, num_cells, construct_nonType); // num_cells
  cout << "undifferentiated cells created" << endl;

  // 3. Define the substances that cells may secrete
  // Order: substance_name, diffusion_coefficient, decay_constant, resolution
  // if diffusion_coefficient is low, diffusion distance is short
  // if decay_constant is high, diffusion distance is short
  // resolution is number of point in one domaine dimension
  ModelInitializer::DefineSubstance(0, "on_diffusion", 0.05, 0.99, param->max_bound_/10);
  ModelInitializer::DefineSubstance(1, "off_diffusion", 0.05, 0.99, param->max_bound_/10);

  // set write output param
  // if you want to write file for RI and cell position
  bool writeRI = true;
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
