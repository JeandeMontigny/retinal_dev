#ifndef EXTENDED_OBJECTS_
#define EXTENDED_OBJECTS_

#include "biodynamo.h"
#include "neuroscience/neuroscience.h"

namespace bdm {
using namespace std;

  // Define my custom cell MyCell extending NeuronSoma
  class MyCell : public experimental::neuroscience::NeuronSoma {
    BDM_SIM_OBJECT_HEADER(MyCell, experimental::neuroscience::NeuronSoma, 1, cell_type_,
      internal_clock_, labelSWC_, previous_position_, distance_travelled_);

   public:
    MyCell() {}
    explicit MyCell(const array<double, 3>& position) : Base(position) {}

    /// Default event constructor
    MyCell(const Event& event, SimObject* other, uint64_t new_oid = 0)
      : Base(event, other, new_oid) {}

    // Default event handler
    void EventHandler(const Event& event, SimObject* other1,
                      SimObject* other2 = nullptr) override {
      Base::EventHandler(event, other1, other2);
    }

    void SetCellType(int t) { cell_type_ = t; }
    int GetCellType() const { return cell_type_; }
    // This function is used by ParaView for coloring the cells by their type
    // int* GetCellTypePtr() { return cell_type_.data(); }

    void SetInternalClock(int t) { internal_clock_ = t; }
    int GetInternalClock() const { return internal_clock_; }

    inline void SetLabel(int label) { labelSWC_ = label; }
    inline int GetLabel() const { return labelSWC_; }
    inline void IncreaseLabel() { labelSWC_ = labelSWC_ + 1; }

    void SetPreviousPosition(array<double, 3> position) { previous_position_ = position; }
    const array<double, 3>& GetPreviousPosition() const { return previous_position_; }

    void SetDistanceTravelled(double distance) { distance_travelled_ = distance; }
    double GetDistanceTravelled() const {return distance_travelled_; }

   private:
    int cell_type_;
    int internal_clock_;
    int labelSWC_;
    array<double, 3> previous_position_;
    double distance_travelled_;
  };


  // Define my custom neurite MyNeurite, which extends NeuriteElement
  class MyNeurite : public experimental::neuroscience::NeuriteElement {
    BDM_SIM_OBJECT_HEADER(MyNeurite, experimental::neuroscience::NeuriteElement, 1,
      has_to_retract_, beyond_threshold_, sleep_mode_,
      diam_before_retract_, subtype_, its_soma_);

   public:
    MyNeurite() {}
    explicit MyNeurite(const array<double, 3>& position) : Base(position) {}

    // Default event constructor
    MyNeurite(const Event& event, SimObject* other,
              uint64_t new_oid = 0) : Base(event, other, new_oid) {
      this->subtype_ = other->subtype_;
      this->its_soma_ = other->its_soma_;
    }

    MyNeurite(const experimental::neuroscience::NewNeuriteExtensionEvent& event,
      SimObject* other, uint64_t new_oid = 0) : Base(event, other, new_oid) {
      this->its_soma_ = other->GetSoPtr();
    }

    // Default event handler
    void EventHandler(const Event& event, SimObject* other1,
                      SimObject* other2 = nullptr) override {
      Base::EventHandler(event, other1, other2);
    }

    void SetHasToRetract(int r) { has_to_retract_ = r; }
    bool GetHasToRetract() const { return has_to_retract_; }

    void SetBeyondThreshold(int r) { beyond_threshold_ = r; }
    bool GetBeyondThreshold() const { return beyond_threshold_; }

    void SetSleepMode(int r) { sleep_mode_ = r; }
    bool GetSleepMode() const { return sleep_mode_; }

    void SetDiamBeforeRetraction(double d) { diam_before_retract_ = d; }
    double GetDiamBeforeRetraction() const { return diam_before_retract_; }

    void SetSubtype(int st) { subtype_ = st; }
    int GetSubtype() { return subtype_; }
    // ParaView
    // TNeuronSomaSoPtr* GetSubtypePtr() { return subtype_.data(); }

    void SetMySoma(NeuronSomaSoPtr soma) { its_soma_ = soma; }
    TNeuronSomaSoPtr GetMySoma() { return its_soma_; }

   private:
    bool has_to_retract_;
    bool beyond_threshold_;
    bool sleep_mode_;
    int diam_before_retract_;
    int subtype_;
    TNeuronSomaSoPtr its_soma_;
  };

} // end namespace bdm

#endif
