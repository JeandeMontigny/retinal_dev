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
    using BiologyModules = CTList<Neurite_creation_BM>;
  };

  BDM_CTPARAM_FOR(bdm, MyNeurite) {
    using BiologyModules =
        CTList<RGC_dendrite_growth_BM>;
  };
};


/* -------- simulate -------- */
template <typename TSimulation = Simulation<>>
inline int Simulate(int argc, const char** argv) {
  // number of simulation steps
  int maxStep = 200;
  // Create an artificial bounds for the simulation space
  int cubeDim = 1000; // 500
  int num_cells = 17600; // 4400
  double cellDensity = (double)num_cells * 1e6 / (cubeDim * cubeDim);
  cout << "cell density: " << cellDensity << " cells per cm2" << endl;

  // set write output param
  bool writeSWC = true;

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
//  mySeed = 9784;  // 9784
  random->SetSeed(mySeed);
  cout << "modelling with seed " << mySeed << endl;

  // min position, max position, number of cells , cell type
  CellCreator(param->min_bound_, param->max_bound_, num_cells/3, 0);
  cout << "on cells created" << endl;
  CellCreator(param->min_bound_, param->max_bound_, num_cells/3, 1);
  cout << "off cells created" << endl;
  CellCreator(param->min_bound_, param->max_bound_, num_cells/3, 2);
  cout << "on-off cells created" << endl;
  CellCreator(param->min_bound_, param->max_bound_, 0, -1);
  cout << "undifferentiated cells created" << endl;

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
  cout << "substances initialised" << endl;

  // 4. Run simulation for maxStep timesteps
  scheduler->Simulate(maxStep);

  // remove previous swc files and export neurons morphology
  if (writeSWC &&
      !system(Concat("mkdir -p ", param->output_dir_, "/swc_files").c_str())) {
    if (system(Concat("rm -r ", param->output_dir_, "/swc_files/*").c_str())) {
      cout << "could not delete previous swc files" << endl;
    };
    morpho_exporteur();
  }

  return 0;
} // end Simulate

}  // end namespace bdm

#endif  // RETINAL_DEV_
