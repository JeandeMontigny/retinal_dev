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
    using BiologyModules = CTList<RGC_mosaic_BM,
                            Substance_secretion_BM>;
  };

};


/* -------- simulate -------- */
template <typename TSimulation = Simulation<>>
inline int Simulate(int argc, const char** argv) {

  double startParameterStudyDeath = 1.21;
  double maxParameterStudyDeath = 1.23;
  double ParamStepDeath = 0.01;

  double startParameterStudyMovement = 1.20;
  double maxParameterStudyMovement = 1.22;
  double ParamStepMovement = 0.01;

  int numberOfIteration = 5;

  // check if output file doesn't exists
  if (access("param_RI_study.txt", F_OK ) == -1) {
    // create output file and return true if failed
    if (system("touch ./param_RI_study.txt")) {
      std::cout << "could not create output file" << std::endl;
    }
  }

  for (double deathThreshold = startParameterStudyDeath;
       deathThreshold < maxParameterStudyDeath;
       deathThreshold += ParamStepMovement) {

    for (double movementThreshold = startParameterStudyMovement;
         movementThreshold < maxParameterStudyMovement;
         movementThreshold += ParamStepMovement) {

           cout << "\nmodelling with parameters "
           << movementThreshold << " / " << deathThreshold
           << " out of " << maxParameterStudyMovement << " / " << maxParameterStudyDeath << " (step of " << ParamStepMovement << " / " << ParamStepDeath << ")" << endl;

      for (int iteration = 0; iteration < numberOfIteration; iteration++) {

        // number of simulation steps
        int maxStep = 1800;
        // Create an artificial bounds for the simulation space
        int cubeDim = 500;
        int num_cells = 4400;
        double diffusion_coef = 0.65;
        double decay_const = 0.1;

        double cellDensity = (double)num_cells * 1e6 / (cubeDim * cubeDim);

        auto set_param = [&](auto* param) {
          // cell are created with +100 to min and -100 to max
          param->bound_space_ = true;
          param->min_bound_ = 0;
          param->max_bound_ = cubeDim + 200;
          // neuroscience/param.h
          param->my_parameter_1_ = movementThreshold;
          param->my_parameter_2_ = deathThreshold;
        };

        Simulation<> simulation(argc, argv, set_param);
        auto* rm = simulation.GetResourceManager();
        auto* random = simulation.GetRandom();
        auto* scheduler = simulation.GetScheduler();
        auto* param = simulation.GetParam();

        int mySeed = rand() % 10000;
      //  mySeed = 9784;
        random->SetSeed(mySeed);

        cout << "iteration " << iteration+1 << " out of "
             << numberOfIteration << endl;

        // min position, max position, number of cells , cell type
        CellCreator(param->min_bound_, param->max_bound_, 0, 0);
        CellCreator(param->min_bound_, param->max_bound_, 0, 1);
        CellCreator(param->min_bound_, param->max_bound_, 0, 2);
        CellCreator(param->min_bound_, param->max_bound_, num_cells, -1);

        // 3. Define substances
        ModelInitializer::DefineSubstance(0, "on_diffusion", diffusion_coef, decay_const,
                                          param->max_bound_ / 2);
        ModelInitializer::DefineSubstance(1, "off_diffusion", diffusion_coef, decay_const,
                                          param->max_bound_ / 2);
        ModelInitializer::DefineSubstance(2, "on_off_diffusion", diffusion_coef, decay_const,
                                          param->max_bound_ / 2);

        // 4. Run simulation for maxStep timesteps
        scheduler->Simulate(maxStep);

        auto my_cells = rm->template Get<MyCell>();
        int numberOfCells = my_cells->size();

        int tempsMigrationDist = 0;
        vector<array<double, 3>> coordList;

        for (int cellNum = 0; cellNum < numberOfCells; cellNum++) {
          array<double, 3> positionAtCreation = (*my_cells)[cellNum].GetPreviousPosition();
          array<double, 3> currentPosition = (*my_cells)[cellNum].GetPosition();
          double distance = sqrt(pow(currentPosition[0] - positionAtCreation[0], 2) +
                                 pow(currentPosition[1] - positionAtCreation[1], 2));
          tempsMigrationDist += distance;
        }


        double ri0 = getRI(0);
        double ri1 = getRI(1);
        double ri2 = getRI(2);
        // mean ri
        double meanRi = (ri0+ri1+ri2)/3;
        // ri std
        double stdRi = sqrt((pow(ri0 - meanRi, 2) + pow(ri1 - meanRi, 2)
                          + pow(ri2 - meanRi, 2)) / 3);
        double cellDeath = (1 - ((double)numberOfCells / num_cells)) * 100;

        cout << "  average ri: " << meanRi << " with std: " << stdRi << " ; "
             << "cell death: " << cellDeath << endl;

        cout << movementThreshold << " " << movementThreshold << " "
             << ri0 << " " << ri1 << " " << ri2 << " "
             << cellDensity << " " << cellDeath << " "
              << (double)tempsMigrationDist/numberOfCells << endl;

        ofstream param_outputFile;
        param_outputFile.open("param_RI_study.txt", std::ios::app);

        param_outputFile << movementThreshold << " " << movementThreshold << " "
                         << ri0 << " " << ri1 << " " << ri2 << " "
                         << cellDensity << " " << cellDeath << " "
                         << (double)tempsMigrationDist/numberOfCells
                         << "\n";

        param_outputFile.close();
    } // end for iteration
  } // end for parameter movement
} // end fot parameter death

  return 0;

} // end Simulate

}  // end namespace bdm

#endif  // RETINAL_DEV_
