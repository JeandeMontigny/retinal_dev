#ifndef EXTENDED_OBJECTS_
#define EXTENDED_OBJECTS_

#include "biodynamo.h"
#include "neuroscience/neuroscience.h"

namespace bdm {
using namespace std;

  // Define my custom cell Progenitor extending NeuronSoma
  // BDM_SIM_OBJECT(Progenitor, experimental::neuroscience::NeuronSoma) {
  //   BDM_SIM_OBJECT_HEADER(ProgenitorExt, 1, internal_clock_);
  //
  //  public:
  //   ProgenitorExt() {}
  //
  //   ProgenitorExt(const array<double, 3>& position) : Base(position) {}
  //
  //   /// Default event constructor
  //   template <typename TEvent, typename TOther>
  //   ProgenitorExt(const TEvent& event, TOther* other, uint64_t new_oid = 0) {
  //   }
  //
  //   /// Default event handler (exising biology module won't be modified on
  //   /// any event)
  //   template <typename TEvent, typename... TOthers>
  //   void EventHandler(const TEvent& event, TOthers*... others) {
  //     Base::EventHandler(event, others...);
  //   }
  //
  //   void SetInternalClock(int t) { internal_clock_[kIdx] = t; }
  //   int GetInternalClock() const { return internal_clock_[kIdx]; }
  //
  //  private:
  //   vec<int> internal_clock_;
  // };


  // Define my custom cell MyCell extending NeuronSoma
  BDM_SIM_OBJECT(MyCell, experimental::neuroscience::NeuronSoma) {
    BDM_SIM_OBJECT_HEADER(MyCellExt, 1, cell_type_, internal_clock_,
      labelSWC_, previous_position_, distance_travelled_);

   public:
    MyCellExt() {}

    MyCellExt(const array<double, 3>& position) : Base(position) {}

    /// Default event constructor
    template <typename TEvent, typename TOther>
    MyCellExt(const TEvent& event, TOther* other, uint64_t new_oid = 0) : Base(event, other, new_oid){
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
    inline int GetLabel() { return labelSWC_[kIdx]; }
    inline void IncreaseLabel() { labelSWC_[kIdx] = labelSWC_[kIdx] + 1; }

    void SetPreviousPosition(array<double, 3> position) { previous_position_[kIdx] = position; }
    array<double, 3> GetPreviousPosition() { return previous_position_[kIdx]; }

    void SetDistanceTravelled(double distance) { distance_travelled_[kIdx] = distance; }
    double GetDistanceTravelled() {return distance_travelled_[kIdx]; }

   private:
    vec<int> cell_type_;
    vec<int> internal_clock_;
    vec<int> labelSWC_;
    vec<array<double, 3>> previous_position_;
    vec<double> distance_travelled_;
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
    MyNeuriteExt(const TEvent& event, TOther* other, uint64_t new_oid = 0) : Base(event, other, new_oid) {
      subtype_[kIdx] = other->subtype_[other->kIdx];
      its_soma_[kIdx] = other->its_soma_[kIdx];
    }

    template <typename TOther>
    MyNeuriteExt(const experimental::neuroscience::NewNeuriteExtensionEvent& event,
      TOther* other, uint64_t new_oid = 0) : Base(event, other, new_oid) {
      its_soma_[kIdx] = other->GetSoPtr();
    }

    // MyNeuriteExt(const experimental::neuroscience::SplitNeuriteElementEvent& event,
    //   TNeurite* proximal) {
    //     //TODO: delete BM of the mother
    //     // bdm::experimental::neuroscience::SplitNeuriteElementEvent::kEventId
    //   }

    /// Default event handler
    template <typename TEvent, typename... TOthers>
    void EventHandler(const TEvent& event, TOthers*... others) {
      Base::EventHandler(event, others...);
    }

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

} // end namespace bdm

#endif
