#ifndef RGC_SOMA_BM_
#define RGC_SOMA_BM_

#include "rgc_dendrite_bm.h"
#include "extended_objects.h"

namespace bdm {
using namespace std;

  // Define cell behavior for neurite creation
  struct Neurite_creation_BM: public BaseBiologyModule {
    Neurite_creation_BM() : BaseBiologyModule(gNullEventId) {}

    /// Default event constructor
    template <typename TEvent, typename TBm>
    Neurite_creation_BM(const TEvent& event, TBm* other, uint64_t new_oid = 0) {}

    template <typename TEvent, typename... TBms>
    void EventHandler(const TEvent&, TBms*...) {}

    template <typename T, typename TSimulation = Simulation<>>
    void Run(T* soma) {

      if (!init_) {
        auto* sim = TSimulation::GetActive();
        auto* random = sim->GetRandom();
        // dendrite per cell: average=4.5; std=1.2
        int thisSubType = soma->GetCellType()*100 + (int)random->Uniform(0, 20);
        for (int i = 0; i <= (int)random->Uniform(2, 7); i++) {
          auto&& ne = soma->ExtendNewNeurite({0, 0, 1});
          ne->AddBiologyModule(RGC_dendrite_growth_BM());
          ne->SetHasToRetract(false);
          ne->SetSleepMode(false);
          ne->SetBeyondThreshold(false);
          ne->SetSubtype(thisSubType);
          ne->SetMySoma(soma->GetSoPtr());

          init_ = true;
        }
      }
    } // end run

  private:
    bool init_ = false;
    ClassDefNV(Neurite_creation_BM, 1);
  }; // endNeurite_creation_BM

} // end namespace bdm

#endif
