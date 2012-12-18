// Andrei Gaponenko, 2012

#include "MuCapG4/inc/MuCapMaterials.hh"

#include <iostream>

#include "cetlib/exception.h"

#include "CLHEP/Units/SystemOfUnits.h"

#include "G4Element.hh"
#include "G4Material.hh"
#include "G4NistManager.hh"

#include "Mu2eG4/inc/findMaterialOrThrow.hh"

namespace mu2e {

  void MuCapMaterials::construct(){
    using CLHEP::mg;
    using CLHEP::cm3;

    int ncomp=0, natoms=0;

    const double temperature0 = 273.15;
    const double absDetectorTemperature =  temperature0 + pset_.get<double>("TCelsius");
    const double gasDensityFactorFrom15C = (15 + temperature0)/absDetectorTemperature;

    //----------------------------------------------------------------
    const double cf4_density_at_15C = 3.72 * mg/cm3;
    G4Material *cf4 =
      new G4Material("MUCAP_CF4", cf4_density_at_15C * gasDensityFactorFrom15C, ncomp=2);

    cf4->AddElement(getElementOrThrow("C"), natoms=1);
    cf4->AddElement(getElementOrThrow("F"), natoms=4);

    //----------------------------------------------------------------
    const double dme_density_at_15C = 1.97 * mg/cm3;
    G4Material *dme =
      new G4Material("MUCAP_DME", dme_density_at_15C * gasDensityFactorFrom15C, ncomp=3);

    dme->AddElement(getElementOrThrow("O"), natoms=1);
    dme->AddElement(getElementOrThrow("C"), natoms=2);
    dme->AddElement(getElementOrThrow("H"), natoms=6);

    //----------------------------------------------------------------
    const double isobutane_density_at_15C = 2.51 * mg/cm3;
    G4Material *isobutane =
      new G4Material("MUCAP_ISOBUTANE", isobutane_density_at_15C * gasDensityFactorFrom15C, ncomp=2);

    isobutane->AddElement(getElementOrThrow("C"), natoms=4);
    isobutane->AddElement(getElementOrThrow("H"), natoms=10);

    //----------------------------------------------------------------
    const double cf4MassFraction = pset_.get<double>("pcGasCF4MassFraction");

    const double pcGasDensity = 1./
      (cf4MassFraction/cf4->GetDensity() + (1-cf4MassFraction)/isobutane->GetDensity());

    G4Material *pcgas = new G4Material("MUCAP_PC_GAS", pcGasDensity, ncomp=2);

    pcgas->AddMaterial(cf4, cf4MassFraction);
    pcgas->AddMaterial(isobutane, (1-cf4MassFraction));

    //----------------------------------------------------------------
    G4Material *G4He = findMaterialOrThrow("G4_He");
    G4Material *G4N2 = findMaterialOrThrow("G4_N");

    const double n2VolumeFraction = pset_.get<double>("cradleGasN2VolumeFraction");
    const double n2MassFraction = n2VolumeFraction * G4N2->GetDensity()
      /(n2VolumeFraction * G4N2->GetDensity() + (1-n2VolumeFraction)*G4He->GetDensity());

    const double cradleGasDensity =
      n2VolumeFraction * G4N2->GetDensity() + (1-n2VolumeFraction)*G4He->GetDensity();

    std::cout<<"MUCAP_CRADLE_GAS: N2 volume fraction = "<<n2VolumeFraction
             <<", mass fraction = "<<n2MassFraction
             <<", density = "<<(cradleGasDensity)/(mg/cm3)
             <<" mg/cm3"<<std::endl;

    G4Material *cradlegas = new G4Material("MUCAP_CRADLE_GAS", cradleGasDensity, ncomp=2);

    cradlegas->AddMaterial(G4N2, n2MassFraction);
    cradlegas->AddMaterial(G4He, 1-n2MassFraction);

    //----------------------------------------------------------------
    // Print element table, if requested.
    if(pset_.get<bool>("printElements")) {
      std::cout<<"MuCapMaterials printout of elements"<<std::endl;
      std::cout<< *G4Element::GetElementTable()<<std::endl;
    }

    // Print material table, if requested.
    if(pset_.get<bool>("printMaterials")) {
      std::cout<<"MuCapMaterials printout of materials"<<std::endl;
      std::cout<<*G4Material::GetMaterialTable()<<std::endl;
    }
  }

  //================================================================
  G4Element* MuCapMaterials::getElementOrThrow(const G4String& name) {
    G4NistManager* nistMan = G4NistManager::Instance();
    G4Element* answer = nistMan->FindOrBuildElement(name,true);
    if(!answer){
      throw cet::exception("GEOM")
        << "mu2e::MuCapMaterials: "
        << "Could not load predefined G4 element named: "
        << name
        << "\n";
    }
    return answer;
  }

} // end namespace mu2e
