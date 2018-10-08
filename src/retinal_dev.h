#ifndef RETINAL_DEV_
#define RETINAL_DEV_

#include "biodynamo.h"
#include "neuroscience/neuroscience.h"
#include "substance_initializers.h"

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
  on_diffusion,
};

// Define my custom cell MyCell extending NeuronSoma
BDM_SIM_OBJECT(MyCell, experimental::neuroscience::NeuronSoma) {
  BDM_SIM_OBJECT_HEADER(MyCellExt, 1, cell_type_, internal_clock_, labelSWC_);

 public:
  MyCellExt() {}

  MyCellExt(const array<double, 3>& position) : Base(position) {}

  /// Default event constructor
  template <typename TEvent, typename TOther>
  MyCellExt(const TEvent& event, TOther* other, uint64_t new_oid = 0) {
    // TODO(jean) implement
  }

  /// Default event handler (exising biology module won't be modified on
  /// any event)
  template <typename TEvent, typename... TOthers>
  void EventHandler(const TEvent&, TOthers*...) {
    // TODO(jean) implement
  }

  void SetCellType(int t) { cell_type_[kIdx] = t; }
  int GetCellType() const { return cell_type_[kIdx]; }
  // This function is used by ParaView for coloring the cells by their type
  int* GetCellTypePtr() { return cell_type_.data(); }

  void SetInternalClock(int t) { internal_clock_[kIdx] = t; }
  int GetInternalClock() const { return internal_clock_[kIdx]; }

  inline void SetLabel(int label) { labelSWC_[kIdx] = label; }
  inline int GetLabel() { return labelSWC_[kIdx]; }
  inline void IncreaseLabel() { labelSWC_[kIdx] = labelSWC_[kIdx] + 1; }

 private:
  vec<int> cell_type_;
  vec<int> internal_clock_;
  vec<int> labelSWC_;
};

// Define my custom neurite MyNeurite, which extends NeuriteElement
BDM_SIM_OBJECT(MyNeurite, experimental::neuroscience::NeuriteElement) {
  BDM_SIM_OBJECT_HEADER(MyNeuriteExt, 1, has_to_retract_, beyond_threshold_,
                        sleep_mode_, diam_before_retract_, subtype_, its_soma_);

 public:
  MyNeuriteExt() {}
  MyNeuriteExt(const array<double, 3>& position) : Base(position) {}

  using NeuronSoma = typename TCompileTimeParam::NeuronSoma;
  using NeuronSomaSoPtr = ToSoPtr<NeuronSoma>;

  /// Default event constructor
  template <typename TEvent, typename TOther>
  MyNeuriteExt(const TEvent& event, TOther* other, uint64_t new_oid = 0) {
    subtype_[kIdx] = other->subtype_[other->kIdx];
    its_soma_[kIdx] = other->its_soma_[kIdx];
  }

  template <typename TOther>
  MyNeuriteExt(const experimental::neuroscience::NewNeuriteExtensionEvent& event, TOther* other, uint64_t new_oid = 0) {
    its_soma_[kIdx] = other->GetSoPtr();
  }

  /// Default event handler
  template <typename TEvent, typename... TOthers>
  void EventHandler(const TEvent&, TOthers*...) {}

  void SetHasToRetract(int r) { has_to_retract_[kIdx] = r; }
  bool GetHasToRetract() const { return has_to_retract_[kIdx]; }

  void SetBeyondThreshold(int r) { beyond_threshold_[kIdx] = r; }
  bool GetBeyondThreshold() const { return beyond_threshold_[kIdx]; }

  void SetSleepMode(int r) { sleep_mode_[kIdx] = r; }
  bool GetSleepMode() const { return sleep_mode_[kIdx]; }

  void SetDiamBeforeRetraction(double d) { diam_before_retract_[kIdx] = d; }
  double GetDiamBeforeRetraction() const { return diam_before_retract_[kIdx]; }

  void SetSubtype(int st) { subtype_[kIdx] = st; }
  int GetSubtype() { return subtype_[kIdx]; }
  // ParaView
  NeuronSomaSoPtr* GetSubtypePtr() { return subtype_.data(); }

  void SetMySoma(NeuronSomaSoPtr soma) { its_soma_[kIdx] = soma; }
  NeuronSomaSoPtr GetMySoma() { return its_soma_[kIdx]; }

 private:
  vec<bool> has_to_retract_;
  vec<bool> beyond_threshold_;
  vec<bool> sleep_mode_;
  vec<int> diam_before_retract_;
  vec<int> subtype_;
  vec<NeuronSomaSoPtr> its_soma_;
};

// Define dendrites behavior for RGC dendritic growth
struct RGC_dendrite_growth_BM : public BaseBiologyModule {
  RGC_dendrite_growth_BM() : BaseBiologyModule(gAllEventIds) {}

  /// Default event constructor
  template <typename TEvent, typename TBm>
  RGC_dendrite_growth_BM(const TEvent& event, TBm* other,
                           uint64_t new_oid = 0) {
  }

  /// Default event handler (exising biology module won't be modified on
  /// any event)
  template <typename TEvent, typename... TBms>
  void EventHandler(const TEvent&, TBms*...) {
  }

  template <typename T, typename TSimulation = Simulation<>>
  void Run(T* ne) {

  }  // end run

  ClassDefNV(RGC_dendrite_growth_BM, 1);
}; // end RGC_dendrite_growth_BM

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
    // if on cell, secrete on cells substance
    // dg_0_->IncreaseConcentrationBy(secretion_position, 1e5/pow(dg_0_->GetBoxLength(), 3));
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

  using SimObjectTypes = CTList<MyCell, MyNeurite>;
  using NeuronSoma = MyCell;
  using NeuriteElement = MyNeurite;

  BDM_CTPARAM_FOR(bdm, MyCell) {
    using BiologyModules = CTList<Substance_secretion_BM>;
  };

  BDM_CTPARAM_FOR(bdm, MyNeurite) {
    using BiologyModules =
        CTList<RGC_dendrite_growth_BM>;
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
  // auto&& soma = rm->New<MyCell>(pos);
  // soma.SetDiameter(12);
  // soma.SetCellType(0);
  // soma.AddBiologyModule(Substance_secretion_BM());

  pos = {50,50,50};
  auto&& soma2 = rm->New<MyCell>(pos);
  soma2.SetDiameter(12);
  soma2.SetCellType(0);
  soma2.AddBiologyModule(Substance_secretion_BM());

  // pos = {65,50,50};
  // auto&& soma3 = rm->New<MyCell>(pos);
  // soma3.SetDiameter(12);
  // soma3.SetCellType(0);
  // soma3.AddBiologyModule(Substance_secretion_BM());

  // 3. Define substances
  ModelInitializer::DefineSubstance(0, "on_diffusion", 1, 0.5,
                                    param->max_bound_);

  // 4. Run simulation for maxStep timesteps
  scheduler->Simulate(100);

  // print diffusion_study.txt
  DiffusionGrid* dg_0_ = rm->GetDiffusionGrid("on_diffusion");
  ofstream diffu_outputFile;
  diffu_outputFile.open("diffusion_study.txt");

  for (double dist=0; dist < 100; dist=dist+0.1) {
    double concentration = dg_0_->GetConcentration({dist, 50, 50});
    diffu_outputFile << dist << " " << concentration << "\n";
  }
  diffu_outputFile.close();
  cout << "done" << endl;

  return 0;
} // end Simulate

}  // end namespace bdm

#endif  // RETINAL_DEV_
