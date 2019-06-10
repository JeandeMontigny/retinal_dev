#ifndef RETINAL_DEV_
#define RETINAL_DEV_
#include "biodynamo.h"
#include "neuroscience/neuroscience.h"
#include "substance_initializers.h"

#include "extended_objects.h"
#include "rgc_soma_bm.h"
#include "util_methods.h"

namespace bdm {
using namespace std;

enum Substances {
  on_diffusion,
  off_diffusion,
  on_off_diffusion,
};

// define compile time parameter
BDM_CTPARAM(experimental::neuroscience) {
  BDM_CTPARAM_HEADER(experimental::neuroscience);

  using SimObjectTypes = CTList<MyCell>;
  using NeuronSoma = MyCell;

  BDM_CTPARAM_FOR(bdm, MyCell) {
    using BiologyModules = CTList<Vertical_migration_BM, RGC_mosaic_BM,
                            Substance_secretion_BM>;
  };

};


/* -------- simulate -------- */
template <typename TSimulation = Simulation<>>
inline int Simulate(int argc, const char** argv) {

  // check if output file doesn't exists
  if (access("param_RI_study_density.txt", F_OK ) == -1) {
    // create output file and return true if failed
    if (system("touch ./param_RI_study_density.txt")) {
      std::cout << "could not create output file" << std::endl;
    }
  }

  int start_num_cells = 50;
  int max_num_cells = 4450; // 4450
  int step_num_cells = 100;

  int numberOfIteration = 5; // 5
  bool writePositionExport = true;
  bool writeIndividualRI = true;

    for (int num_cells = start_num_cells; num_cells <= max_num_cells;
         num_cells += step_num_cells) {

      cout << "modelling with parameters "
           << num_cells << " out of " << max_num_cells <<" (step of "
           << step_num_cells << ")" << endl;

      for (int iteration = 0; iteration < numberOfIteration; iteration++) {
        int maxStep = 100;
        int cubeDim = 500;

        auto set_param = [&](auto* param) {
          param->bound_space_ = true;
          param->min_bound_ = 0;
          param->max_bound_ = cubeDim;
      	  // set output directory name
      	  param->output_dir_ = "output/" + to_string(num_cells);
        };

        Simulation<> simulation(argc, argv, set_param);
        auto* random = simulation.GetRandom();
        auto* scheduler = simulation.GetScheduler();
        auto* rm = simulation.GetResourceManager();
        auto* param = simulation.GetParam();

        int mySeed = rand() % 10000;
        //  mySeed = 9784;
        random->SetSeed(mySeed);

        // min position, max position, number of cells , cell type
        // CellCreator<MyCell>(param->min_bound_, param->max_bound_, num_cells, 0);
        // CellCreator<MyCell>(param->min_bound_, param->max_bound_, num_cells/3, 1);
        // CellCreator<MyCell>(param->min_bound_, param->max_bound_, num_cells/3, 2);


        vector<array<double, 3>> cell_pos_list;

        rm->template Reserve<MyCell>(num_cells);

        for (int i = 0; i < num_cells; i++) {
          double x = random->Uniform(param->min_bound_, param->max_bound_);
          double y = random->Uniform(param->min_bound_, param->max_bound_);
          // RGCL thickness before cell death ~24
          double z = random->Uniform(20, 35); // 20, 35
          std::array<double, 3> position = {x, y, z};

          if (!conflict(position, 12)) {
            MyCell cell(position);
            cell.SetDiameter(12);
            cell.SetCellType(0);
            cell.SetInternalClock(0);
            cell.SetDensity(0.001);
            cell.SetPreviousPosition(position);
            cell.AddBiologyModule(Vertical_migration_BM());
            // cell.AddBiologyModule(Substance_secretion_BM());
            // cell.AddBiologyModule(RGC_mosaic_BM());
            rm->push_back(cell);
            cell_pos_list.push_back(position);
          }
        }
        //
        // for (int i = 0; i < num_cells/3; i++) {
        //   double x = random->Uniform(param->min_bound_, param->max_bound_);
        //   double y = random->Uniform(param->min_bound_, param->max_bound_);
        //   // RGCL thickness before cell death ~24
        //   double z = random->Uniform(20, 20); // 20, 35
        //   std::array<double, 3> position = {x, y, z};
        //
        //   if (!conflict(position, 12)) {
        //     MyCell cell(position);
        //     cell.SetDiameter(12);
        //     cell.SetCellType(1);
        //     cell.SetInternalClock(0);
        //     cell.SetDensity(0.001);
        //     cell.SetPreviousPosition(position);
        //     // cell.AddBiologyModule(Substance_secretion_BM());
        //     // cell.AddBiologyModule(RGC_mosaic_BM());
        //     rm->push_back(cell);
        //     cell_pos_list.push_back(position);
        //   }
        // }
        //
        // for (int i = 0; i < num_cells/3; i++) {
        //   double x = random->Uniform(param->min_bound_, param->max_bound_);
        //   double y = random->Uniform(param->min_bound_, param->max_bound_);
        //   // RGCL thickness before cell death ~24
        //   double z = random->Uniform(20, 20); // 20, 35
        //   std::array<double, 3> position = {x, y, z};
        //
        //   if (!conflict(position, 12)) {
        //     MyCell cell(position);
        //     cell.SetDiameter(12);
        //     cell.SetCellType(2);
        //     cell.SetInternalClock(0);
        //     cell.SetDensity(0.001);
        //     cell.SetPreviousPosition(position);
        //     // cell.AddBiologyModule(Substance_secretion_BM());
        //     // cell.AddBiologyModule(RGC_mosaic_BM());
        //     rm->push_back(cell);
        //     cell_pos_list.push_back(position);
        //   }
        // }

        // 4. Run simulation for maxStep timesteps
        double cellDensity;
        if (writeIndividualRI) {
          cellDensity = (double)num_cells * 1e6 / (cubeDim * cubeDim);
          ofstream ri_outputFile;
          ri_outputFile.open(Concat("individual/individual_RI_density_", cellDensity, ".txt").c_str(), std::ios::app);
          for (int step=0; step<maxStep; step++) {
            scheduler->Simulate(1);
            ri_outputFile << cellDensity << " " << getRI(0) << "\n";
          }
        }
        else {
          scheduler->Simulate(maxStep);
        }

      	// export position of every cells
        if (writePositionExport) {
      	  position_exporteur(iteration);
      	}

        cellDensity = (double)num_cells * 1e6 / (cubeDim * cubeDim);

        double ri0 = getRI(0);
        // double ri1 = getRI(1);
        // double ri2 = getRI(2);

        ofstream param_outputFile;
        param_outputFile.open("param_RI_study_density.txt", std::ios::app);

        // param_outputFile << cellDensity << " " << (ri0+ri1+ri2)/3 << "\n";
        param_outputFile << cellDensity << " " << ri0 << "\n";

        param_outputFile.close();
    } // end for iteration
} // end for num_cells

  return 0;

} // end Simulate

}  // end namespace bdm

#endif  // RETINAL_DEV_
