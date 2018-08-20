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

  template <typename TSimulation = Simulation<>>
  struct CellMovement : public BaseBiologyModule {
    CellMovement() : BaseBiologyModule(gAllBmEvents) {}

    template <typename T>
    void Run(T* sim_object) {
      auto&& cell = sim_object->template ReinterpretCast<experimental::neuroscience::NeuronSoma>();

      cell->UpdatePosition({1, 0, 0});
    }

    ClassDefNV(CellMovement, 1);
  };


  // Define compile time parameter
  template <typename Backend>
  struct CompileTimeParam : public DefaultCompileTimeParam<Backend> {
    using BiologyModules = Variant<CellMovement<>>;
    using AtomicTypes = VariadicTypedef<experimental::neuroscience::NeuronSoma, experimental::neuroscience::NeuriteElement>;
    using NeuronSoma = experimental::neuroscience::NeuronSoma;
    using NeuriteElement = experimental::neuroscience::NeuriteElement;
  };

  template <typename TSimulation = Simulation<>>
  inline int Simulate(int argc, const char** argv) {
    Simulation<> simulation(argc, argv);
    auto* rm = simulation.GetResourceManager();
    auto* scheduler = simulation.GetScheduler();

    experimental::neuroscience::NeuronSoma cell({20, 50, 50});
    cell.SetDiameter(10);
    cell.AddBiologyModule(CellMovement<>());
    auto ne = cell.ExtendNewNeurite({0, 0, 1});
    ne->GetSoPtr()->SetDiameter(1);
    rm->push_back(cell);

    scheduler->Simulate(100);

    std::cout << ne->GetSoPtr()->GetPosition()[1] << std::endl;

    return 0;
    }
}  // namespace bdm

#endif // RETINAL_DEV_
