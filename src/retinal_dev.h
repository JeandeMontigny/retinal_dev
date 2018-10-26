#ifndef RETINAL_DEV_
#define RETINAL_DEV_

#include "biodynamo.h"
#include "neuroscience/neuroscience.h"
#include "substance_initializers.h"

#include "extended_objects.h"

#include "progenitor_bm.h"
#include "rgc_soma_bm.h"
#include "rgc_dendrite_bm.h"
#include "amacrine_soma_bm.h"
#include "amacrine_dendrite_bm.h"

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
  // using NeuronSoma = Progenitor;
  using NeuronSoma = MyCell;
  using NeuriteElement = MyNeurite;

  // BDM_CTPARAM_FOR(bdm, Progenitor) {
  //   using BiologyModules = CTList<Progenitor_behaviour_BM>;
  // };

  BDM_CTPARAM_FOR(bdm, MyCell) {
    using BiologyModules =
      CTList<Progenitor_behaviour_BM, RGC_axial_migration_BM,
      Amacrine_axial_migration_BM, Amacrine_Neurite_creation_BM,
      RGC_mosaic_BM, Substance_secretion_BM, Rgc_Neurite_creation_BM>;
  };

  BDM_CTPARAM_FOR(bdm, MyNeurite) {
    using BiologyModules =
      CTList<Amacrine_dendrite_growth_BM, RGC_dendrite_growth_BM>;
  };
};


/* -------- simulate -------- */
template <typename TSimulation = Simulation<>>
inline int Simulate(int argc, const char** argv) {
  // number of simulation steps
  int maxStep = 2200;
  // Create an artificial bounds for the simulation space
  int cubeDim = 500;
  int num_cells = 4400;
  double cellDensity = (double)num_cells * 1e6 / (cubeDim * cubeDim);
  cout << "cell density: " << cellDensity << " cells per cm2" << endl;

  // set write output param
  // if you want to write file for RI and cell position
  bool writeRI = false;
  bool writePositionExport = false;
  bool writeSWC = false;
  bool writeMigrationDistance = false;
  // terminal print and create cell position files every outputFrequence steps
  int outputFrequence = 100;

  auto set_param = [&](auto* param) {
    // cell are created with +100 to min and -100 to max
    param->bound_space_ = true;
    param->min_bound_ = 0;
    param->max_bound_ = cubeDim + 200;
    // set min and max length for neurite segments
    param->neurite_min_length_ = 1.0;
    param->neurite_max_length_ = 2.0;
  };

  Simulation<> simulation(argc, argv, set_param);
  auto* random = simulation.GetRandom();
  auto* scheduler = simulation.GetScheduler();
  auto* param = simulation.GetParam();

  int mySeed = rand() % 10000;
  mySeed = 9784;  // 9784
  random->SetSeed(mySeed);
  cout << "modelling with seed " << mySeed << endl;

  ProgenitorsCreator(param->min_bound_, param->max_bound_, num_cells);

  // min position, max position, number of cells , cell type
  MyCellCreator(param->min_bound_, param->max_bound_, 0, 0);
  // cout << "on cells created" << endl;
  MyCellCreator(param->min_bound_, param->max_bound_, 0, 1);
  // cout << "off cells created" << endl;
  MyCellCreator(param->min_bound_, param->max_bound_, 0, 2);
  // cout << "on-off cells created" << endl;
  MyCellCreator(param->min_bound_, param->max_bound_, 0, -1);
  // cout << "undifferentiated cells created" << endl;

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

  // gaussian substances
  // create substance for Progenitor guide
  ModelInitializer::DefineSubstance(5, "progenitors_guide", 0, 0,
                                    param->max_bound_ / 2);
  ModelInitializer::InitializeSubstance(5, "progenitors_guide",
                                        GaussianBand(param->min_bound_,
                                          param->max_bound_ / 5, Axis::kZAxis));

  // define substances for RGC dendrite attraction
  ModelInitializer::DefineSubstance(3, "on_substance_RGC_guide", 0, 0,
                                    param->max_bound_ / 2);
  ModelInitializer::DefineSubstance(4, "off_substance_RGC_guide", 0, 0,
                                    param->max_bound_ / 2);
  // create substance gaussian distribution for RGC dendrite attraction
  // average peak distance for ON cells: 15.959 with std of 5.297;
  ModelInitializer::InitializeSubstance(3, "on_substance_RGC_guide",
                                        GaussianBand(45, 6, Axis::kZAxis));
  // average peak distance for OFF cells: 40.405 with std of 8.39;
  ModelInitializer::InitializeSubstance(4, "off_substance_RGC_guide",
                                        GaussianBand(69, 8, Axis::kZAxis));
  cout << "substances created" << endl;

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
  int stepNb = 0;
  if (writeRI) {
    maxStep = (int)maxStep / 10;
    stepNb = 10;
  }
  else {
    maxStep = (int)maxStep / outputFrequence;
    stepNb = outputFrequence;
  }

  for (int i = 0; i <= maxStep; i++) {
    // simulate for stepNb steps
    scheduler->Simulate(stepNb);

    if (writeRI) {
      double RIon = getRI(0);
      double RIoff = getRI(1);
      double RIonOff = getRI(2);
      outputFile << RIon << " " << RIoff << " " << RIonOff << "\n";

      if (i % 10 == 0) {
        PrintTerminalOutput(i*stepNb, maxStep*10, num_cells);
      }
    } // end if writeRI
    else {
      PrintTerminalOutput(i*stepNb, maxStep*outputFrequence, num_cells);
    }

    if (writePositionExport && i*stepNb % outputFrequence == 0) {
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

  return 0;
} // end Simulate

}  // end namespace bdm

#endif  // RETINAL_DEV_
