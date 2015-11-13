#include "../interface/Tau.hh"
#include "../interface/city.h"
#include "boost/format.hpp"

namespace ic {

Tau::Tau()
    : ic::Candidate(),
      decay_mode_(0),
      lead_ecal_energy_(0.),
      lead_hcal_energy_(0.),
      lead_p_(0.),
      lead_dxy_vertex_(0.),
      lead_dz_vertex_(0.) {}

Tau::~Tau() {}

void Tau::Print() const { 
  std::cout << boost::format("%-17s | ")   % "p4";
  Candidate::Print();
  std::cout << boost::format("%-17s | %10i\n")   % "decay_mode" % decay_mode_;
}

void Tau::SetTauID(std::string const& name, float const& value) {
  tau_ids_[CityHash64(name)] = value;
}

float Tau::GetTauID(std::string const& name) const {
  UFmap::const_iterator iter = tau_ids_.find(CityHash64(name));
  if (iter != tau_ids_.end()) {
    return iter->second;
  } else {
    std::cerr << "Warning in <Tau::GetTauID>: Algorithm \"" << name
              << "\" not found" << std::endl;
    return 0.0;
  }
}

bool Tau::HasTauID(std::string const& name) const {
  return tau_ids_.count(CityHash64(name)) > 0;
}
}
