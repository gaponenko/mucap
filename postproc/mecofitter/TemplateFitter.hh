// Andrei Gaponenko, 2012

#ifndef postproc_mecofitter_TemplateFitter_hh
#define postproc_mecofitter_TemplateFitter_hh

#include "fhiclcpp/ParameterSet.h"

#include "FitRegion.hh"
#include "FitData.hh"
#include "FitMC.hh"

class TH2;
class TH3;

namespace mucap {
  namespace fitter {
    class TemplateFitter {
      fhicl::ParameterSet pset_;

      FitRegion region_;

      FitRegion computeFitRegion(const fhicl::ParameterSet& cuts,
                                 const TH2 *data,
                                 const TH2 *mcproj) const;

      TH2 *visualizeFitRegion(const TH2 *binningTemplate, const FitRegion &r) const;

      FitData getFitData(const TH2 *data, const FitRegion& r) const;
      FitMC getFitMC(const TH3 *hmc, const FitRegion& r) const;

    public:
      explicit TemplateFitter(const fhicl::ParameterSet& pset);

    };
  } // fitter
} // mucap

#endif/*postproc_mecofitter_TemplateFitter_hh*/
