#ifndef RETINAL_DEV_
#define RETINAL_DEV_

#include "biodynamo.h"
#include "neuroscience/neuroscience.h"
#include "substance_initializers.h"

#include "extended_objects.h"
#include "rgc_soma_bm.h"
#include "rgc_dendrite_bm.h"
#include "util_methods.h"

namespace bdm {
using namespace std;

enum Substances {
  on_diffusion,
  off_diffusion,
  on_substance_RGC_guide,
  off_substance_RGC_guide
};


// define compile time parameter
BDM_CTPARAM(experimental::neuroscience) {
  BDM_CTPARAM_HEADER(experimental::neuroscience);

  using SimObjectTypes = CTList<MyCell, MyNeurite>;
  using NeuronSoma = MyCell;
  using NeuriteElement = MyNeurite;

  BDM_CTPARAM_FOR(bdm, MyCell) {
    using BiologyModules = CTList<RGC_mosaic_BM, Substance_secretion_BM, Neurite_creation_BM>;
  };

  BDM_CTPARAM_FOR(bdm, MyNeurite) {
    using BiologyModules =
        CTList<RGC_dendrite_growth_BM>;
  };
};


/* -------- simulate -------- */
template <typename TSimulation = Simulation<>>
inline int Simulate(int argc, const char** argv) {

  ofstream param_outputFile;
  param_outputFile.open("param_RI_study.txt");

  for (int parameter = 0; parameter < 0.1; parameter = parameter + 0.1) {

    // auto simulation = new Simulation<>(argc, argv);
    // simulation->Activate();

    // number of simulation steps
    int maxStep = 1800;
    // Create an artificial bounds for the simulation space
    int cubeDim = 500;
    int num_cells = 4400;
    double cellDensity = (double)num_cells * 1e6 / (cubeDim * cubeDim);
    cout << "cell density: " << cellDensity << " cells per cm2" << endl;

    // set write output param
    // if you want to write file for RI and cell position
    bool writeRI = true;
    bool writePositionExport = false;
    bool writeSWC = false;
    bool writeMigrationDistance = true;
    // create cell position files every outputFrequence steps
    int outputFrequence = 100;

    auto set_param = [&](auto* param) {
      // cell are created with +100 to min and -100 to max
      param->bound_space_ = true;
      param->min_bound_ = 0;
      param->max_bound_ = cubeDim + 100;
    };

    Simulation<> simulation(argc, argv, set_param);
    auto* rm = simulation.GetResourceManager();
    auto* random = simulation.GetRandom();
    auto* scheduler = simulation.GetScheduler();
    auto* param = simulation.GetParam();

    int mySeed = rand() % 10000;
  //  mySeed = 9784;  // 9784
    random->SetSeed(mySeed);
    cout << "modelling with seed " << mySeed << endl;

    // min position, max position, number of cells , cell type
    CellCreator(param->min_bound_, param->max_bound_, 0, 0); //num_cells/3
    cout << "on cells created" << endl;
    CellCreator(param->min_bound_, param->max_bound_, 0, 1);
    cout << "off cells created" << endl;
    CellCreator(param->min_bound_, param->max_bound_, 0, 2);
    cout << "on-off cells created" << endl;
    CellCreator(param->min_bound_, param->max_bound_, num_cells, -1);
    cout << "undifferentiated cells created" << endl;

    // 3. Define substances
    // Order: substance_name, diffusion_coefficient, decay_constant, resolution
    // if diffusion_coefficient is low, diffusion distance is short
    // if decay_constant is high, diffusion distance is short
    // resolution is number of point in one domaine dimension
    ModelInitializer::DefineSubstance(0, "on_diffusion", 0.65, 0,
                                      param->max_bound_ / 2);
    ModelInitializer::DefineSubstance(1, "off_diffusion", 0.65, 0,
                                      param->max_bound_ / 2);
    ModelInitializer::DefineSubstance(2, "on_off_diffusion", 0.65, 0,
                                      param->max_bound_ / 2);

    cout << "substances initialised" << endl;

    // prepare export
    ofstream outputFile;
    // delete previous simulation export
    if (writePositionExport && system(
      Concat("mkdir -p ", param->output_dir_, "/cells_position").c_str())) {
      cout << "error in " << param->output_dir_
           << "/cells_position folder creation" << endl;
    }
    // create RI export file
    if (writeRI && !writePositionExport &&
      system(Concat("mkdir -p ", param->output_dir_, "/cells_position").c_str())){
      cout << "error in " << param->output_dir_
           << "/cells_position folder creation" << endl;
    } else {
      outputFile.open(Concat(param->output_dir_, "/cells_position/RI_" +
                      to_string(mySeed) + ".txt"));
    }

    // 4. Run simulation for maxStep timesteps
    for (int i = 0; i <= maxStep; i++) {
      scheduler->Simulate(1);

      if (i % 10 == 0) {  // write RI in file
        double RIon = getRI(0);
        double RIoff = getRI(1);
        double RIonOff = getRI(2);
        if (writeRI) {
          outputFile << RIon << " " << RIoff << " " << RIonOff << "\n";
        }

        // print
        if (i % 100 == 0) {
          // get cell list size
          rm = simulation.GetResourceManager();
          auto my_cells = rm->template Get<MyCell>();
          int numberOfCells = my_cells->size();
          // TODO: vector for unknow number of cell type
          int numberOfCells0 = 0;
          int numberOfCells1 = 0;
          int numberOfCells2 = 0;
          int numberOfDendrites = 0;

          for (int cellNum = 0; cellNum < numberOfCells; cellNum++) {
            numberOfDendrites += (*my_cells)[cellNum].GetDaughters().size();
            auto thisCellType = (*my_cells)[cellNum].GetCellType();
            if (thisCellType == 0) { numberOfCells0++; }
            else if (thisCellType == 1) { numberOfCells1++; }
            else if (thisCellType == 2) { numberOfCells2++; }
          }

          cout << "-- step " << i << " out of " << maxStep << " --\n"
               << numberOfCells << " cells in simulation: "
               << (1 - ((double)numberOfCells / num_cells)) * 100
               << "% of cell death\n"
               << numberOfCells0 << " cells are type 0 (on) ; " << numberOfCells1
               << " cells are type 1 (off) ; " << numberOfCells2
               << " cells are type 2 (on-off) ; "
               << (double)(numberOfCells0 + numberOfCells1 + numberOfCells2) /
                      numberOfCells * 100
               << "% got type\n"
               << numberOfDendrites << " apical dendrites in simulation: "
               << (double)numberOfDendrites / numberOfCells << " dendrites per cell\n"
               << "RI on: " << RIon << " ; RI off: " << RIoff
               << " ; RI on-off: " << RIonOff
               << " ; mean: " << (RIon + RIoff + RIonOff) / 3 << endl;
        }  // end every 100 simu steps
      }    // end every 10 simu steps

      if (writePositionExport && i % outputFrequence == 0) {
        // export cell position
        position_exporteur(i);
      }

    }  // end for Simulate

    outputFile.close();

    // remove previous swc files and export neurons morphology
    if (writeSWC &&
        !system(Concat("mkdir -p ", param->output_dir_, "/swc_files").c_str())) {
      if (system(Concat("rm -r ", param->output_dir_, "/swc_files/*").c_str())) {
        cout << "could not delete previous swc files" << endl;
      };
      morpho_exporteur();
    }

    if (writeMigrationDistance) {
      exportMigrationDist();
    }


    // delete simulation;

} // end for parameter

  param_outputFile.close();


  return 0;

} // end Simulate

}  // end namespace bdm

#endif  // RETINAL_DEV_
