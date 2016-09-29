/* 
   +----------------------------------------------------------------------+
   |                                 HEP                                  |
   |                             --- HEP ---                              |
   |                             Module File                              |
   +----------------------------------------------------------------------+
   HEP coherent system of Units AND HEP coherent Physical Constants
  
   This file has been provided to CLHEP by Geant4 (simulation toolkit for HEP).
  
   The basic units are :
   millimeter              (millimeter)
   nanosecond              (nanosecond)
   Mega electron Volt      (MeV)
   positron charge         (eplus)
   degree Kelvin           (kelvin)
   the amount of substance (mole)
   luminous intensity      (candela)
   radian                  (radian)
   steradian               (steradian)
  
   Below is a non exhaustive list of derived and pratical units
   (i.e. mostly the SI units).
   You can add your own units.
  
   The SI numerical value of the positron charge is defined here,
   as it is needed for conversion factor : positron charge = e_SI (coulomb)
  
  
   Below is a non exhaustive list of Physical CONSTANTS,
   computed in the Internal HEP System Of Units.
  
   Most of them are extracted from the Particle Data Book :
          Phys. Rev. D  volume 50 3-1 (1994) page 1233
   
          ...with a meaningful (?) name ...
  
   You can add your own constants.

   +----------------------------------------------------------------------+
   | JavaScript                                                           |
   +----------------------------------------------------------------------+
   F. Quinonez - Created 2016-01-30                  
               - Geant4 File
   +----------------------------------------------------------------------+
   | C++                                                                  |
   +----------------------------------------------------------------------+
   PhysicalUnits.h:
   Authors: M.Maire, S.Giani
  
   History:
  
   06.02.96   Created.
   28.03.96   Added miscellaneous constants.
   05.12.97   E.Tcherniaev: Redefined pascal (to avoid warnings on WinNT)
   20.05.98   names: meter, second, gram, radian, degree
              (from Brian.Lasiuk@yale.edu (STAR)). Added luminous units.
   05.08.98   angstrom, picobarn, microsecond, picosecond, petaelectronvolt
   01.03.01   parsec    
   31.01.06   kilogray, milligray, microgray    
   29.04.08   use PDG 2006 value of e_SI
   03.11.08   use PDG 2008 value of e_SI
   19.08.15   added liter and its sub units (mma)

   PhysicalConstants.h: 
   Author: M.Maire
  
   History:
  
   23.02.96 Created
   26.03.96 Added constants for standard conditions of temperature
            and pressure; also added Gas threshold.
   29.04.08   use PDG 2006 values
   03.11.08   use PDG 2008 values
*/

  "use strict";

  //constructor: HEP,
  function HEP(){
  }


  HEP.pi  = 3.14159265358979323846;
  HEP.twopi  = 2 * HEP.pi;
  HEP.halfpi  = HEP.pi / 2;
  HEP.pi2 = HEP.pi * HEP.pi;

  // 
  // Length [L]
  //
  HEP.millimeter  = 1.;                        
  HEP.millimeter2 = HEP.millimeter * HEP.millimeter;
  HEP.millimeter3 = HEP.millimeter2 * HEP.millimeter;

  HEP.centimeter  = 10. * HEP.millimeter;   
  HEP.centimeter2 = HEP.centimeter * HEP.centimeter;
  HEP.centimeter3 = HEP.centimeter2 * HEP.centimeter;

  HEP.meter  = 1000. * HEP.millimeter;                  
  HEP.meter2 = HEP.meter * HEP.meter;
  HEP.meter3 = HEP.meter2 * HEP.meter;

  HEP.kilometer = 1000. * HEP.meter;                   
  HEP.kilometer2 = HEP.kilometer * HEP.kilometer;
  HEP.kilometer3 = HEP.kilometer2 * HEP.kilometer;

  HEP.parsec = 3.0856775807e+16 * HEP.meter;

  HEP.micrometer = 1.e-6 * HEP.meter;             
  HEP.nanometer = 1.e-9 * HEP.meter;
  HEP.angstrom  = 1.e-10 * HEP.meter;
  HEP.fermi     = 1.e-15 * HEP.meter;

  HEP.barn = 1.e-28 * HEP.meter2;
  HEP.millibarn = 1.e-3 * HEP.barn;
  HEP.microbarn = 1.e-6 * HEP.barn;
  HEP.nanobarn = 1.e-9 * HEP.barn;
  HEP.picobarn = 1.e-12 * HEP.barn;

  // symbols
  HEP.nm  = HEP.nanometer;                        
  HEP.um  = HEP.micrometer;                        

  HEP.mm  = HEP.millimeter;                        
  HEP.mm2 = HEP.millimeter2;
  HEP.mm3 = HEP.millimeter3;

  HEP.cm  = HEP.centimeter;   
  HEP.cm2 = HEP.centimeter2;
  HEP.cm3 = HEP.centimeter3;

  HEP.liter = 1.e+3 * HEP.cm3;
  HEP.L = HEP.liter;
  HEP.dL = 1.e-1 * HEP.liter;
  HEP.cL = 1.e-2 * HEP.liter;
  HEP.mL = 1.e-3 * HEP.liter;       

  HEP.m  = HEP.meter;                  
  HEP.m2 = HEP.meter2;
  HEP.m3 = HEP.meter3;

  HEP.km  = HEP.kilometer;                   
  HEP.km2 = HEP.kilometer2;
  HEP.km3 = HEP.kilometer3;

  HEP.pc = HEP.parsec;

  //
  // Angle
  //
  HEP.radian      = 1.;                  
  HEP.milliradian = 1.e-3 * HEP.radian;
  HEP.degree = (HEP.pi/180.0) * HEP.radian;

  HEP.steradian = 1.;
  
  // symbols
  HEP.rad  = HEP.radian;
  HEP.mrad = HEP.milliradian;
  HEP.sr   = HEP.steradian;
  HEP.deg  = HEP.degree;

  //
  // Time [T]
  //
  HEP.nanosecond  = 1.;
  HEP.second      = 1.e+9  * HEP.nanosecond;
  HEP.millisecond = 1.e-3  * HEP.second;
  HEP.microsecond = 1.e-6  * HEP.second;
  HEP.picosecond = 1.e-12 * HEP.second;

  HEP.hertz = 1. / HEP.second;
  HEP.kilohertz = 1.e+3 * HEP.hertz;
  HEP.megahertz = 1.e+6 * HEP.hertz;

  // symbols
  HEP.ns = HEP.nanosecond;
  HEP.s = HEP.second;
  HEP.ms = HEP.millisecond;

  //
  // Electric charge [Q]
  //
  HEP.eplus = 1. ;// positron charge
  HEP.e_SI  = 1.602176487e-19;// positron charge in coulomb
  HEP.coulomb = HEP.eplus / HEP.e_SI;// coulomb = 6.24150 e+18  *  eplus

  //
  // Energy [E]
  //
  HEP.megaelectronvolt = 1. ;
  HEP.electronvolt = 1.e-6 * HEP.megaelectronvolt;
  HEP.kiloelectronvolt = 1.e-3 * HEP.megaelectronvolt;
  HEP.gigaelectronvolt = 1.e+3 * HEP.megaelectronvolt;
  HEP.teraelectronvolt = 1.e+6 * HEP.megaelectronvolt;
  HEP.petaelectronvolt = 1.e+9 * HEP.megaelectronvolt;

  HEP.joule = HEP.electronvolt / HEP.e_SI;// joule = 6.24150 e+12  *  MeV

  // symbols
  HEP.MeV = HEP.megaelectronvolt;
  HEP.eV = HEP.electronvolt;
  HEP.keV = HEP.kiloelectronvolt;
  HEP.GeV = HEP.gigaelectronvolt;
  HEP.TeV = HEP.teraelectronvolt;
  HEP.PeV = HEP.petaelectronvolt;

  //
  // Mass [E][T^2][L^-2]
  //
  HEP.kilogram = HEP.joule * HEP.second * HEP.second/(HEP.meter * HEP.meter);   
  HEP.gram = 1.e-3 * HEP.kilogram;
  HEP.milligram = 1.e-3 * HEP.gram;

  // symbols
  HEP.kg = HEP.kilogram;
  HEP.g = HEP.gram;
  HEP.mg = HEP.milligram;

  //
  // Power [E][T^-1]
  //
  HEP.watt = HEP.joule / HEP.second;// watt = 6.24150 e+3  *  MeV/ns

  //
  // Force [E][L^-1]
  //
  HEP.newton = HEP.joule / HEP.meter; 
  // newton = 6.24150 e+9  *  MeV/mm

  //
  // Pressure [E][L^-3]
  //
  HEP.pascal = HEP.newton/HEP.m2;   // pascal = 6.24150 e+3  *  MeV/mm3
  HEP.bar        = 100000 * HEP.pascal; // bar    = 6.24150 e+8  *  MeV/mm3
  HEP.atmosphere = 101325 * HEP.pascal; // atm    = 6.32420 e+8  *  MeV/mm3

  //
  // Electric current [Q][T^-1]
  //
  HEP.ampere = HEP.coulomb/HEP.second; // ampere = 6.24150 e+9  *  eplus/ns
  HEP.milliampere = 1.e-3 * HEP.ampere;
  HEP.microampere = 1.e-6 * HEP.ampere;
  HEP.nanoampere = 1.e-9 * HEP.ampere;

  //
  // Electric potential [E][Q^-1]
  //
  HEP.megavolt = HEP.megaelectronvolt/HEP.eplus;
  HEP.kilovolt = 1.e-3 * HEP.megavolt;
  HEP.volt = 1.e-6 * HEP.megavolt;

  //
  // Electric resistance [E][T][Q^-2]
  //
  HEP.ohm = HEP.volt / HEP.ampere;// ohm = 1.60217e-16 * (MeV/eplus)/(eplus/ns)

  //
  // Electric capacitance [Q^2][E^-1]
  //
  HEP.farad = HEP.coulomb / HEP.volt;// farad = 6.24150e+24  *  eplus/Megavolt
  HEP.millifarad = 1.e-3 * HEP.farad;
  HEP.microfarad = 1.e-6 * HEP.farad;
  HEP.nanofarad = 1.e-9 * HEP.farad;
  HEP.picofarad = 1.e-12 * HEP.farad;

  //
  // Magnetic Flux [T][E][Q^-1]
  //
  HEP.weber = HEP.volt * HEP.second;// weber = 1000 * megavolt * ns

  //
  // Magnetic Field [T][E][Q^-1][L^-2]
  //
  HEP.tesla     = HEP.volt * HEP.second / HEP.meter2;// tesla =0.001 * megavolt * ns/mm2

  HEP.gauss     = 1.e-4 * HEP.tesla;
  HEP.kilogauss = 1.e-1 * HEP.tesla;

  //
  // Inductance [T^2][E][Q^-2]
  //
  HEP.henry = HEP.weber / HEP.ampere;// henry = 1.60217e-7 * MeV * (ns/eplus) *  * 2

  //
  // Temperature
  //
  HEP.kelvin = 1.;

  //
  // Amount of substance
  //
  HEP.mole = 1.;

  //
  // Activity [T^-1]
  //
  HEP.becquerel = 1./HEP.second ;
  HEP.curie = 3.7e+10  *  HEP.becquerel;
  HEP.kilobecquerel = 1.e+3 * HEP.becquerel;
  HEP.megabecquerel = 1.e+6 * HEP.becquerel;
  HEP.gigabecquerel = 1.e+9 * HEP.becquerel;
  HEP.millicurie = 1.e-3 * HEP.curie;
  HEP.microcurie = 1.e-6 * HEP.curie;
  HEP.Bq = HEP.becquerel;
  HEP.kBq = HEP.kilobecquerel;
  HEP.MBq = HEP.megabecquerel;
  HEP.GBq = HEP.gigabecquerel;
  HEP.Ci = HEP.curie;
  HEP.mCi = HEP.millicurie;
  HEP.uCi = HEP.microcurie;

  //
  // Absorbed dose [L^2][T^-2]
  //
  HEP.gray = HEP.joule / HEP.kilogram ;
  HEP.kilogray = 1.e+3 * HEP.gray;
  HEP.milligray = 1.e-3 * HEP.gray;
  HEP.microgray = 1.e-6 * HEP.gray;

  //
  // Luminous intensity [I]
  //
  HEP.candela = 1.;

  //
  // Luminous flux [I]
  //
  HEP.lumen = HEP.candela * HEP.steradian;

  //
  // Illuminance [I][L^-2]
  //
  HEP.lux = HEP.lumen / HEP.meter2;

  //
  // Miscellaneous
  //
  HEP.perCent     = 0.01 ;
  HEP.perThousand = 0.001;
  HEP.perMillion  = 0.000001;


//
// 
//
HEP.Avogadro = 6.02214179e+23 / HEP.mole;

//
// c   = 299.792458 mm/ns
// c^2 = 898.7404 (mm/ns)^2 
//
HEP.c_light   = 2.99792458e+8  *  HEP.m / HEP.s;
HEP.c_squared = HEP.c_light  *  HEP.c_light;

//
// h     = 4.13566e-12 MeV * ns
// hbar  = 6.58212e-13 MeV * ns
// hbarc = 197.32705e-12 MeV * mm
//
HEP.h_Planck      = 6.62606896e-34  *  HEP.joule * HEP.s;
HEP.hbar_Planck   = HEP.h_Planck / HEP.twopi;
HEP.hbarc         = HEP.hbar_Planck  *  HEP.c_light;
HEP.hbarc_squared = HEP.hbarc  *  HEP.hbarc;

//
//
//
HEP.electron_charge = - HEP.eplus; // see SystemOfUnits.h
HEP.e_squared = HEP.eplus  *  HEP.eplus;

//
// amu_c2 - atomic equivalent mass unit
//        - AKA, unified atomic mass unit (u)
// amu    - atomic mass unit
//
HEP.electron_mass_c2 = 0.510998910  *  HEP.MeV;
HEP.proton_mass_c2 = 938.272013  *  HEP.MeV;
HEP.neutron_mass_c2 = 939.56536  *  HEP.MeV;
HEP.amu_c2 = 931.494028  *  HEP.MeV;
HEP.amu = HEP.amu_c2 / HEP.c_squared;

//
// permeability of free space mu0    = 2.01334e-16 Mev * (ns * eplus)^2/mm
// permittivity of free space epsil0 = 5.52636e+10 eplus^2/(MeV * mm)
//
HEP.mu0      = 4 * HEP.pi * 1.e-7  *  HEP.henry/HEP.m;
HEP.epsilon0 = 1./(HEP.c_squared * HEP.mu0);

//
// electromagnetic coupling = 1.43996e-12 MeV * mm/(eplus^2)
//
HEP.elm_coupling           = HEP.e_squared / (4 * HEP.pi * HEP.epsilon0);
HEP.fine_structure_const   = HEP.elm_coupling / HEP.hbarc;
HEP.classic_electr_radius  = HEP.elm_coupling / HEP.electron_mass_c2;
HEP.electron_Compton_length = HEP.hbarc / HEP.electron_mass_c2;
HEP.Bohr_radius = HEP.electron_Compton_length / HEP.fine_structure_const;

HEP.alpha_rcl2 = HEP.fine_structure_const * HEP.classic_electr_radius * HEP.classic_electr_radius;

HEP.twopi_mc2_rcl2 = HEP.twopi * HEP.electron_mass_c2 * HEP.classic_electr_radius * HEP.classic_electr_radius;
//
//
//
HEP.k_Boltzmann = 8.617343e-11  *  HEP.MeV / HEP.kelvin;

//
//
//
HEP.STP_Temperature = 273.15 * HEP.kelvin;
HEP.STP_Pressure    = 1. * HEP.atmosphere;
HEP.kGasThreshold   = 10. * HEP.mg / HEP.cm3;

//
//
//
HEP.universe_mean_density = 1.e-25 * HEP.g / HEP.cm3;



module.exports = HEP;



