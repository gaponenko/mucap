// Andrei Gaponenko, 2012

#ifndef postproc_mecofitter_TemplateFitter_hh
#define postproc_mecofitter_TemplateFitter_hh

#include <memory>

#include "fhiclcpp/ParameterSet.h"

#include "Math/Minimizer.h"

#include "FitRegion.hh"
#include "FitData.hh"
#include "FitMC.hh"
#include "TargetFunctionMECO.hh"

class TH2;
class TH3;

namespace mucap {
  namespace fitter {
    class TemplateFitter {
      fhicl::ParameterSet pset_;

      std::unique_ptr<const TH2> origData_;
      std::unique_ptr<const TH3> origMC_;
      FitRegion region_;
      TargetFunctionMECO func_;
      std::unique_ptr<ROOT::Math::Minimizer> min_;

      TH2 *visualizeFitRegion(const TH2 *binningTemplate, const FitRegion &r) const;
      TH2 *visualizeFitData(const std::string& hname,
                            const FitData& data,
                            const TH2 *binningTemplate,
                            const FitRegion &r) const;


      static TH2 *readDataHisto(const fhicl::ParameterSet& pset);
      static TH3 *readMCHisto(const fhicl::ParameterSet& pset);
      static FitRegion computeFitRegion(const fhicl::ParameterSet& cuts,
                                        const TH2 *data,
                                        const TH3 *mc);

      static FitData extractFitData(const TH2 *data, const FitRegion& r);
      static FitMC extractFitMC(const TH3 *hmc, const FitRegion& r);



      static ROOT::Math::Minimizer *createMinimizer(const fhicl::ParameterSet& pset);

    public:

      explicit TemplateFitter(const fhicl::ParameterSet& pset);

      ROOT::Math::Minimizer& getMinimizer() { return *min_; }

      const TargetFunctionMECO& getFunction() const { return func_; }

      TH2 *visualizeFitRegion() const {
        return visualizeFitRegion(origData_.get(), region_);
      }

      TH2 *visualizeFitData(const std::string& name, const FitData& data) const {
        return visualizeFitData(name, data, origData_.get(), region_);
      }
    };
  } // fitter
} // mucap

#endif/*postproc_mecofitter_TemplateFitter_hh*/
