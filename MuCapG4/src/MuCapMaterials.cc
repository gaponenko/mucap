// Andrei Gaponenko, 2012

#include "MuCapG4/inc/MuCapMaterials.hh"

#include <iostream>

#include "cetlib/exception.h"

#include "CLHEP/Units/SystemOfUnits.h"

#include "G4Element.hh"
#include "G4Material.hh"
#include "G4NistManager.hh"

#include "Mu2eG4/inc/findMaterialOrThrow.hh"

namespace mucap {

  void MuCapMaterials::construct(){
    using CLHEP::g;
    using CLHEP::mg;
    using CLHEP::cm3;

    int ncomp=0, natoms=0;

    const double P_STP = pset_.get<double>("P_STP");
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
    G4Material *tec_dme =
      new G4Material("MUCAP_TEC_GAS",
                     (pset_.get<double>("tecGasPressure")/P_STP)*
                     dme_density_at_15C * gasDensityFactorFrom15C,
                     ncomp=3);

    tec_dme->AddElement(getElementOrThrow("O"), natoms=1);
    tec_dme->AddElement(getElementOrThrow("C"), natoms=2);
    tec_dme->AddElement(getElementOrThrow("H"), natoms=6);

    //----------------------------------------------------------------
    const double isobutane_density_at_15C = 2.51 * mg/cm3;
    G4Material *isobutane =
      new G4Material("MUCAP_ISOBUTANE", isobutane_density_at_15C * gasDensityFactorFrom15C, ncomp=2);

    isobutane->AddElement(getElementOrThrow("C"), natoms=4);
    isobutane->AddElement(getElementOrThrow("H"), natoms=10);

    //----------------------------------------------------------------
    const double cf4VolumeFraction = pset_.get<double>("pcGasCF4VolumeFraction");
    const double cf4MassFraction = cf4VolumeFraction * cf4->GetDensity()
      /(cf4VolumeFraction * cf4->GetDensity() + (1-cf4VolumeFraction)*isobutane->GetDensity());

    const double pcGasDensity =
      cf4VolumeFraction*cf4->GetDensity() + (1-cf4VolumeFraction)*isobutane->GetDensity();

    std::cout<<"MUCAP_PC_GAS: CF4 volume fraction = "<<cf4VolumeFraction
             <<", mass fraction = "<<cf4MassFraction
             <<", density = "<<(pcGasDensity)/(mg/cm3)<<" mg/cm3"<<std::endl;

    G4Material *pcgas = new G4Material("MUCAP_PC_GAS", pcGasDensity, ncomp=2);

    pcgas->AddMaterial(cf4, cf4MassFraction);
    pcgas->AddMaterial(isobutane, (1-cf4MassFraction));

    //----------------------------------------------------------------
    G4Material *G4He = mu2e::findMaterialOrThrow("G4_He");
    G4Material *G4N2 = mu2e::findMaterialOrThrow("G4_N");

    const double densityHe = G4He->GetDensity() * G4He->GetTemperature()/absDetectorTemperature;
    const double densityN2 = G4N2->GetDensity() * G4N2->GetTemperature()/absDetectorTemperature;

    const double n2VolumeFraction = pset_.get<double>("cradleGasN2VolumeFraction");
    const double n2MassFraction = n2VolumeFraction * densityN2
      /(n2VolumeFraction * densityN2 + (1-n2VolumeFraction)*densityHe);

    const double cradleGasDensity = n2VolumeFraction * densityN2 + (1-n2VolumeFraction)*densityHe;

    std::cout<<"MUCAP_CRADLE_GAS: N2 volume fraction = "<<n2VolumeFraction
             <<", mass fraction = "<<n2MassFraction
             <<", density = "<<(cradleGasDensity)/(mg/cm3)
             <<" mg/cm3"<<std::endl;

    G4Material *cradlegas = new G4Material("MUCAP_CRADLE_GAS", cradleGasDensity, ncomp=2);

    cradlegas->AddMaterial(G4N2, n2MassFraction);
    cradlegas->AddMaterial(G4He, 1-n2MassFraction);

    //----------------------------------------------------------------
    G4Material *G4CO2 = mu2e::findMaterialOrThrow("G4_CARBON_DIOXIDE");

    G4double densityCO2  = G4CO2->GetDensity() * G4CO2->GetTemperature()/absDetectorTemperature;
    G4double co2VolumeFraction  = pset_.get<double>("gabsGasCO2VolumeFraction");

    const double co2MassFraction = co2VolumeFraction * densityCO2
      /(co2VolumeFraction * densityCO2 + (1-co2VolumeFraction)*densityHe);

    double gabs_density = co2VolumeFraction*densityCO2 + (1.0-co2VolumeFraction)*densityHe;

    std::cout<<"MUCAP_GABS_GAS: CO2 volume fraction = "<<co2VolumeFraction
             <<", mass fraction = "<<co2MassFraction
             <<", density = "<<(gabs_density)/(mg/cm3)
             <<" mg/cm3"<<std::endl;

    G4Material *gabs_gas = new G4Material("MUCAP_GABS_GAS", gabs_density, ncomp=2);

    gabs_gas->AddMaterial(G4CO2, co2MassFraction);
    gabs_gas->AddMaterial(G4He, 1-co2MassFraction);

    //----------------------------------------------------------------
    // FR4 was sometimes referred to as "G10" in TWIST, but is a
    // different newer material.  This is a fiberglass-epoxy
    // composite.  For the lack of better information I'll use
    // pre-defined G4 plastic and glass materials to build it.

    G4Material *fr4 = new G4Material("MUCAP_FR4", 1.85*g/cm3, ncomp=2);
    fr4->AddMaterial(mu2e::findMaterialOrThrow("G4_BAKELITE"), 0.44);
    fr4->AddMaterial(mu2e::findMaterialOrThrow("G4_Pyrex_Glass"), 0.56);

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

} // end namespace mucap
