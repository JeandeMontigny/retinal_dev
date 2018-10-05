#ifndef AMACRINE_SOMA_BM_
#define AMACRINE_SOMA_BM_

namespace bdm {
using namespace std;

// Define cell behavior for neurite creation
struct Amacrine_Neurite_creation_BM: public BaseBiologyModule {
  Amacrine_Neurite_creation_BM() : BaseBiologyModule(gNullEventId) {}

  /// Default event constructor
  template <typename TEvent, typename TBm>
  Amacrine_Neurite_creation_BM(const TEvent& event, TBm* other, uint64_t new_oid = 0) {}

  template <typename TEvent, typename... TBms>
  void EventHandler(const TEvent&, TBms*...) {}

  template <typename T, typename TSimulation = Simulation<>>
  void Run(T* soma) {

    bool createDendrites = false;

    if (createDendrites) {
    }
  } // end run

private:
  ClassDefNV(Amacrine_Neurite_creation_BM, 1);
}; // end Amacrine_Neurite_creation_BM


// Define cell behavior for RGC migration to GCL
struct Amacrine_axial_migration_BM: public BaseBiologyModule {
  Amacrine_axial_migration_BM() : BaseBiologyModule(gNullEventId) {}

  /// Default event constructor
  template <typename TEvent, typename TBm>
  Amacrine_axial_migration_BM(const TEvent& event, TBm* other, uint64_t new_oid = 0) {}

  template <typename TEvent, typename... TBms>
  void EventHandler(const TEvent&, TBms*...) {}

  template <typename T, typename TSimulation = Simulation<>>
  void Run(T* cell) {
    // TODO: migrate to GCL -> attracted by RGC
  } // end run

private:
  ClassDefNV(Amacrine_axial_migration_BM, 1);
}; // Amacrine_axial_migration_BM


} // end namespace bdm

#endif
