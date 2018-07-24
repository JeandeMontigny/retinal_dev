#ifndef RETINAL_DEV_
#define RETINAL_DEV_

#include "biodynamo.h"
#include "random.h"
#include "substance_initializers.h"
#include "model_initializer.h"

namespace bdm {
// cell type: 0=on; 1=off; -1=not attributed yet
// setprecision (15) <<

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

// Define my custom cell MyCell, which extends Cell
BDM_SIM_OBJECT(MyCell, Cell) {
  //  BDM_SIM_OBJECT(MyCell, Cell) {
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
  int* GetInternalClockPtr() { return internal_clock_.data(); }

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

      bool withCellDeath = false;
      bool withMovement = false;

      // if not initialised, initialise substance diffusions
      if (!init_) {
        auto* rm = sim->GetResourceManager();
        dg_0_ = rm->GetDiffusionGrid("on_diffusion");
        dg_1_ = rm->GetDiffusionGrid("off_diffusion");
        init_ = true;
      }

      auto& position = cell->GetPosition();
      array<double, 3> diff_gradient;
      array<double, 3> gradient_z;
      double concentration = 0;
      int cellClock = cell->GetInternalClock();

      // if cell is type 1, concentration and gradient are the one of substance
      // 1
      if (cell->GetCellType() == 1) {
        dg_1_->GetGradient(position, &gradient_1_);
        gradient_z = Math::ScalarMult(0.08, gradient_1_);
        gradient_z[0] = 0;
        gradient_z[1] = 0;
        diff_gradient = Math::ScalarMult(-0.1, gradient_1_);
        diff_gradient[2] = 0;
        concentration = dg_1_->GetConcentration(position);
      }
      // else if cell is type 2, concentration and gradient are the one of
      // substance 2
      if (cell->GetCellType() == 0) {
        dg_0_->GetGradient(position, &gradient_0_);
        gradient_z = Math::ScalarMult(0.08, gradient_0_);
        gradient_z[0] = 0;
        gradient_z[1] = 0;
        diff_gradient = Math::ScalarMult(-0.1, gradient_0_);
        diff_gradient[2] = 0;
        concentration = dg_0_->GetConcentration(position);
      }

      /* -- cell movement -- */
      if (withMovement && cellClock >= 200) {
        // cell movement based on homotype substance gradient
        // 0. for high density - 0. for normal density // 0. for cell death with layer collapse
        if (concentration >= 0.10475) {
          cell->UpdatePosition(diff_gradient);
          cell->UpdatePosition(gradient_z);
          cell->SetPosition(cell->GetPosition());
        }
        // random movement
        //        cell->UpdatePosition({random->Uniform(-1, 1),
        //        random->Uniform(-1, 1), 0});
      }

      /* -- cell death -- */
      if (withCellDeath && cellClock >= 200) {
        // add small random movements if tangential dispersion is off
        if (!withMovement && cellClock < 400 && concentration >= 0.104725) {
          cell->UpdatePosition({random->Uniform(-0.2, 0.2), random->Uniform(-0.2, 0.2), 0});
          cell->SetPosition(cell->GetPosition());
        }
        // add vertical migration as the multi layer colapse in just on layer
        if (concentration >= 0.1045) {  // 104
          cell->UpdatePosition(gradient_z);
          //          cell->UpdatePosition(diff_gradient);
          cell->SetPosition(cell->GetPosition());
        }
        // cell death depending on homotype substance concentration

        // die if concentration is too high; proba so all cells don't die simultaneously
        // 0.1047 for already existing mosaic - 0.1048 for random position - 0.1047x for inbetween
        if (concentration >= 0.10485 && random->Uniform(0, 1) < 0.1) {
          cell->RemoveFromSimulation();
        }
        // randomly kill ~60% cells (over 250 steps)
        //      if (random->Uniform(0, 1) < 0.004) {
        //        cell->RemoveFromSimulation();
        //      }
      }

      /* -- cell fate -- */
      // cell type attribution depending on concentrations
      if (cell->GetCellType() == -1) {  // if cell type is not on or off
        dg_1_->GetGradient(position, &gradient_1_);
        double concentration_1 = dg_1_->GetConcentration(position);
        dg_0_->GetGradient(position, &gradient_0_);
        double concentration_0 = dg_0_->GetConcentration(position);

        // if no substances
        // random so all cell types doesn't create all randomly at step 1
        // if (concentration_1 == 0 && concentration_0 == 0 && random->Uniform(0, 1) < 0.0001) {
        //   // random attribution of a cell type
        //   cell->SetCellType((int)random->Uniform(0, 2));
        // }

        // if [Off substance] > [On substance]
        // random so all cell doesn't choose their type in the same time
        if (concentration_1 > concentration_0 && concentration_1 > 0.01 // 0.000001
          && random->Uniform(0, 1) < 0.1) {
            cell->SetCellType(0);  // cell become on
        }

        // if [On substance] > [Off substance]
        // random so cells don't choose their type in the same time
        if (concentration_0 > concentration_1 && concentration_0 > 0.01
          && random->Uniform(0, 1) < 0.1) {
          cell->SetCellType(1);  // cell become off
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
  array<double, 3> gradient_0_;
  array<double, 3> gradient_1_;
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
        dg_1_->IncreaseConcentrationBy(secretion_position, 1);  // 0.1
      }
      // is on cell, secrete on cells substance
      else if (cell->GetCellType() == 0) {
        dg_0_->IncreaseConcentrationBy(secretion_position, 1); // 0.1
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
  ofstream positionFileOn;
  ofstream positionFileOff;
  ofstream positionFileAll;
  stringstream sstmOn;
  stringstream sstmOff;
  stringstream sstmAll;
  sstmOn << "on_t" << i << ".txt";
  sstmOff << "off_t" << i << ".txt";
  sstmAll << "all_t" << i << ".txt";

  string fileNameOn = sstmOn.str();
  string fileNameOff = sstmOff.str();
  string fileNameAll = sstmAll.str();
  positionFileOn.open(fileNameOn);  // TODO: include seed number to file name
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
      positionFileAll << position[0] << " " << position[1] << " " << position[2]
                      << " on\n";
    } else {
      positionFileOff << position[0] << " " << position[1] << "\n";
      positionFileAll << position[0] << " " << position[1] << " " << position[2]
                      << " off\n";
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
    for (unsigned int j = 0; j < coordList.size();
         j++) {  // for each other cell of same type in the simulation
      array<double, 3> otherCellPosition = coordList[j];

      double tempsDistance =
          sqrt(pow(cellPosition[0] - otherCellPosition[0], 2) +
               pow(cellPosition[1] - otherCellPosition[1],
                   2));  // get the distance between those 2 cells
      if (tempsDistance < shortestDist &&
          tempsDistance != 0) {        // if cell is not itself
        shortestDist = tempsDistance;  // updade closest neighbour distance
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

template <typename TSimulation = Simulation<>>
inline double getRI(int desiredCellType) {
  auto* sim = TSimulation::GetActive();
  auto* rm = sim->GetResourceManager();
  auto my_cells = rm->template Get<MyCell>();  // get cell list
  vector<array<double, 3>> coordList;          // list of coordinate
  int numberOfCells = my_cells->size();

  for (int cellNum = 0; cellNum < numberOfCells;
       cellNum++) {  // for each cell in simulation
    auto thisCellType = (*my_cells)[cellNum].GetCellType();
    if (thisCellType == desiredCellType) {  // if cell is of the desired type
      auto position = (*my_cells)[cellNum].GetPosition();  // get its position
      coordList.push_back(position);  // put cell coord in the list
    }
  }
  //    cout << coordList.size() << " cells of type " << desiredCellType <<
  //    endl;
  return computeRI(coordList);  // return RI for desired cell type
}  //; end getRI

/* -------- simulate -------- */
// template <typename TResourceManager = ResourceManager<>>
template <typename TSimulation = Simulation<>>
inline int Simulate(int argc, const char** argv) {

  double diffusionCoef = 0.0;
  double decayConst = 0.0;

  for (int eachDiffus = 0; eachDiffus < 10; eachDiffus++) {

    decayConst = 0;
    diffusionCoef += 0.01;

    for (int eachDecay = 0; eachDecay < 10; eachDecay++) {

      decayConst += 0.01;

      auto simulation = new Simulation<>(argc, argv);
      simulation->Activate();

      auto* rm = simulation->GetResourceManager();
      auto* random = simulation->GetRandom();
      auto* scheduler = simulation->GetScheduler();
      auto* param = simulation->GetParam();

      // Create an artificial bounds for the simulation space
      int cubeDim = 40; // 500
      param->bound_space_ = true;
      param->min_bound_ = 0;
      // cells are created with +20 to min and -20 to max. so physical cube has to be cubeDim+40
      param->max_bound_ = cubeDim;

      int mySeed = rand() % 10000;
      mySeed = 8274;  // 9784
      random->SetSeed(mySeed);

      auto* cells = rm->template Get<MyCell>();

      MyCell cell0({10,25,25});
      cell0.SetDiameter(8);
      cell0.SetCellType(0);
      cell0.SetInternalClock(0);
      cell0.AddBiologyModule(SubstanceSecretion<>());
      cell0.AddBiologyModule(Chemotaxis<>());
      cells->push_back(cell0);

      MyCell cell1({20,25,25});
      cell1.SetDiameter(8);
      cell1.SetCellType(-1);
      cell1.SetInternalClock(0);
      cell1.AddBiologyModule(Chemotaxis<>());
      cells->push_back(cell1);

      MyCell cell2({30,50,50});
      cell2.SetDiameter(8);
      cell2.SetCellType(-1);
      cell2.SetInternalClock(0);
      cell2.AddBiologyModule(Chemotaxis<>());
      cells->push_back(cell2);

      cells->Commit();

      // 3. Define the substances that cells may secrete
      // Order: substance_name, diffusion_coefficient, decay_constant, resolution
      // if diffusion_coefficient is low, diffusion distance is short
      // if decay_constant is high, diffusion distance is short
      // resolution is number of point in one domaine dimension
      ModelInitializer::DefineSubstance(0, "on_diffusion", diffusionCoef, decayConst, param->max_bound_);  // 0.5, 0.1
      ModelInitializer::DefineSubstance(1, "off_diffusion", 1, 0.5, param->max_bound_);  // 0.5, 0.1

      // set some param
      // number of simulation steps // 1201
      int maxStep = 1000;
      // if you want to write file for RI and cell position
      // create cell position files every outputFrequence steps
      ofstream outputFile;

      outputFile.open("diffus=" + to_string(diffusionCoef) + "_decay=" + to_string(decayConst) + ".txt");

      auto dg_0_ = rm->GetDiffusionGrid("on_diffusion");

      for (int i = 0; i < maxStep; i++) {
        scheduler->Simulate(1);

        if (i % 10 == 0) {  // write RI in file
          outputFile << dg_0_->GetConcentration(cell1.GetPosition()) << " " << dg_0_->GetConcentration(cell2.GetPosition()) << "\n";
        }

      } // end for i to maxStep
      outputFile.close();

      std::cout << "simulation " << diffusionCoef << " " << decayConst << " done" << std::endl;

      delete simulation;

    } // end for eachDecay
  } // end eachDiffus

  return 0;

} // end Simulate

}  // namespace bdm

#endif // RETINAL_DEV_
