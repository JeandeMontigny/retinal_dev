#ifndef RETINAL_DEV_
#define RETINAL_DEV_

#include "biodynamo.h"
#include "neuroscience/neuroscience.h"
#include "substance_initializers.h"

namespace bdm {
// cell type: 0=on; 1=off; -1=not differentiated yet
// std::cout << setprecision (15) <<


using namespace std;

enum Substances {
  on_diffusion,
};

// Define secretion behavior:
struct Substance_secretion_BM : public BaseBiologyModule {
  // Daughter cells inherit this biology module
  Substance_secretion_BM() : BaseBiologyModule(gAllEventIds) {}

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
      init_ = true;
    }
    auto& secretion_position = cell->GetPosition();
    dg_0_->IncreaseConcentrationBy(secretion_position, 1);
  }

 private:
  bool init_ = false;
  DiffusionGrid* dg_0_ = nullptr;
  ClassDefNV(Substance_secretion_BM, 1);
};  // end biologyModule Substance_secretion_BM

// define compile time parameter
BDM_CTPARAM(experimental::neuroscience) {
  BDM_CTPARAM_HEADER(experimental::neuroscience);

  using SimObjectTypes = CTList<Cell>;

  BDM_CTPARAM_FOR(bdm, Cell) {
    using BiologyModules = CTList<Substance_secretion_BM>;
  };

};

/* -------- simulate -------- */
template <typename TSimulation = Simulation<>>
inline int Simulate(int argc, const char** argv) {

  auto set_param = [&](auto* param) {
    // cell are created with +100 to min and -100 to max
    param->bound_space_ = true;
    param->min_bound_ = 0;
    param->max_bound_ = 100;
  };

  Simulation<> simulation(argc, argv, set_param);
  auto* rm = simulation.GetResourceManager();
  auto* scheduler = simulation.GetScheduler();
  auto* param = simulation.GetParam();

  std::array<double, 3> pos = {35,50,50};
  // auto&& soma = rm->New<Cell>(pos);
  // soma.SetDiameter(12);
  // soma.AddBiologyModule(Substance_secretion_BM());

  pos = {50,50,50};
  auto&& soma2 = rm->New<Cell>(pos);
  soma2.SetDiameter(12);
  soma2.AddBiologyModule(Substance_secretion_BM());

  // pos = {65,50,50};
  // auto&& soma3 = rm->New<Cell>(pos);
  // soma3.SetDiameter(12);
  // soma3.AddBiologyModule(Substance_secretion_BM());

  // 3. Define substances
  ModelInitializer::DefineSubstance(0, "on_diffusion", 0.65, 0, param->max_bound_/2);

  // 4. Run simulation for maxStep timesteps
  scheduler->Simulate(100);

  // print diffusion_study.txt
  DiffusionGrid* dg_0_ = rm->GetDiffusionGrid("on_diffusion");
  ofstream diffu_outputFile;
  diffu_outputFile.open("diffusion_study.txt");

  for (double dist=0; dist < 100; dist=dist+0.1) {
    double concentration = dg_0_->GetConcentration({dist, 50, 50});
//    if (concentration > 0.1) { concentration=0.1; }
    diffu_outputFile << dist << " " << concentration << "\n";
  }
  diffu_outputFile.close();
  cout << "done" << endl;

  return 0;
} // end Simulate

}  // end namespace bdm

#endif  // RETINAL_DEV_
