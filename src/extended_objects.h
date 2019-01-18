#ifndef EXTENDED_OBJECTS_
#define EXTENDED_OBJECTS_

#include "biodynamo.h"
#include "neuroscience/neuroscience.h"

namespace bdm {
using namespace std;

  // Define my custom cell MyCell extending NeuronSoma
  BDM_SIM_OBJECT(MyCell, experimental::neuroscience::NeuronSoma) {
    BDM_SIM_OBJECT_HEADER(MyCell, experimental::neuroscience::NeuronSoma, 1, cell_type_, internal_clock_,
      labelSWC_, previous_position_, distance_travelled_);

   public:
    MyCellExt() {}

    MyCellExt(const array<double, 3>& position) : Base(position) {}

    /// Default event constructor
    template <typename TEvent, typename TOther>
    MyCellExt(const TEvent& event, TOther* other, uint64_t new_oid = 0) {
    }

    /// Default event handler (exising biology module won't be modified on
    /// any event)
    template <typename TEvent, typename... TOthers>
    void EventHandler(const TEvent& event, TOthers*... others) {
      Base::EventHandler(event, others...);
    }

    void SetCellType(int t) { cell_type_[kIdx] = t; }
    int GetCellType() const { return cell_type_[kIdx]; }
    // This function is used by ParaView for coloring the cells by their type
    int* GetCellTypePtr() { return cell_type_.data(); }

    void SetInternalClock(int t) { internal_clock_[kIdx] = t; }
    int GetInternalClock() const { return internal_clock_[kIdx]; }

    inline void SetLabel(int label) { labelSWC_[kIdx] = label; }
    inline int GetLabel() const { return labelSWC_[kIdx]; }
    inline void IncreaseLabel() { labelSWC_[kIdx] = labelSWC_[kIdx] + 1; }

    void SetPreviousPosition(array<double, 3> position) { previous_position_[kIdx] = position; }
    array<double, 3> GetPreviousPosition() const { return previous_position_[kIdx]; }

    void SetDistanceTravelled(double distance) { distance_travelled_[kIdx] = distance; }
    double GetDistanceTravelled() const {return distance_travelled_[kIdx]; }

   private:
    vec<int> cell_type_;
    vec<int> internal_clock_;
    vec<int> labelSWC_;
    vec<array<double, 3>> previous_position_;
    vec<double> distance_travelled_;
  };


} // end namespace bdm

#endif
