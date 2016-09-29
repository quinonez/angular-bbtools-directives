# HEPhys Module
This module written in JavaScript has the units and constants used everyday in High Energy Physics
calculations.

HEPhys is based on CLHEP - A Class Library for High Energy Physics, which is written in C++.
This HEPhys specifically presents the same units and constants that are defined in files PhysicalUnits.h and PhysicalConstants.h.
 

One way to use it is:

```javascript
    var HEP = require('hephys');
    var radiuscircle = 3.4 * HEP.cm; 
    var perimeter = radiuscircle * HEP.twopi;
    var area      = Math.pow( radiuscircle, 2) * HEP.pi;
    console.log( "area = %d cm2", area / HEP.cm2 );
    console.log( "area = %d m2", area / HEP.m2 );
    console.log( "perimeter = %d cm", perimeter / HEP.cm );
    console.log( "perimeter = %d m", perimeter / HEP.m );
```

## Public Functions


## Static Public Variables

#### HEP.pi
#### HEP.pi2
#### HEP.halfpi
#### HEP.twopi
### Length Units
#### HEP.millimeter = 1
#### HEP.mm = 1
#### HEP.millimeter2
#### HEP.millimeter3
#### HEP.centimeter
#### HEP.cm
#### HEP.centimeter2
#### HEP.cm2
#### HEP.centimeter3
#### HEP.cm3
#### HEP.meter
#### HEP.m
#### HEP.meter2
#### HEP.m2
#### HEP.meter3
#### HEP.m3
#### HEP.liter 
#### HEP.L
#### HEP.dL 
#### HEP.cL
#### HEP.mL 
#### HEP.kilometer
#### HEP.km
#### HEP.kilometer2
#### HEP.km2
#### HEP.kilometer3
#### HEP.km3
#### HEP.parsec
#### HEP.pc
#### HEP.micrometer
#### HEP.um
#### HEP.nanometer
#### HEP.nm
#### HEP.angstrom
#### HEP.fermi
#### HEP.barn
#### HEP.millibarn
#### HEP.microbarn
#### HEP.picobarn
### Angle Units
#### HEP.radian = 1
#### HEP.rad = 1
#### HEP.milliradian
#### HEP.mrad
#### HEP.degree
#### HEP.deg
#### HEP.steradian = 1
#### HEP.sr = 1
### Time Units
#### HEP.nanosecond = 1
#### HEP.ns = 1
#### HEP.second
#### HEP.s
#### HEP.millisecond
#### HEP.ms
#### HEP.microsecond
#### HEP.us
#### HEP.picosecond
### Electric Charge Units
#### HEP.eplus = 1
#### HEP.e\_SI
#### HEP.coulomb
### Energy Units
#### HEP.megaelectronvolt = 1
#### HEP.MeV
#### HEP.electronvolt
#### HEP.eV
#### HEP.gigaelectronvolt
#### HEP.GeV
#### HEP.teraelectronvolt
#### HEP.TeV
#### HEP.petaelectronvolt
#### HEP.PeV
#### HEP.joule
### Mass Units
#### HEP.kilogram
#### HEP.kg
#### HEP.gram
#### HEP.g
#### HEP.milligram
#### HEP.mg
### Power Units
#### HEP.watt
### Force Units
#### HEP.newton
### Pressure Units
#### HEP.hep\_pascal
#### HEP.bar
#### HEP.atmosphere
### Electric Current Units
#### HEP.ampere
#### HEP.milliampere
#### HEP.microampere
#### HEP.nanoampere
### Electric Voltage Units
#### HEP.megavolt
#### HEP.kilovolt
#### HEP.volt
### Electric Resistance Units
#### HEP.ohm
### Capacitance Units
#### HEP.farad
#### HEP.millifarad
#### HEP.microfarad
#### HEP.nanofarad
#### HEP.picofarad
### Magnetic Field Units
#### HEP.tesla
#### HEP.gauss
#### HEP.kilogauss
### Magnetic Flux Units
#### HEP.weber
### Inductance Units
#### HEP.henry
### Temperature Units
#### HEP.kelvin = 1
### Amount of Substance Units
#### HEP.mole = 1
### Activity Units
#### HEP.hertz
#### HEP.kilohertz
#### HEP.megahertz
#### HEP.becquerel
#### HEP.Bq
#### HEP.kilobecquerel
#### HEP.kBq
#### HEP.megabecquerel
#### HEP.MBq
#### HEP.gigabecquerel
#### HEP.GBq
#### HEP.curie
#### HEP.Ci
#### HEP.millicurie
#### HEP.mCi
#### HEP.microcurie
#### HEP.uCi
### Absorbed Dose
#### HEP.gray
#### HEP.kilogray
#### HEP.milligray
#### HEP.microgray
### Luminous Intensity
#### HEP.candela = 1
### Luminous Flux
#### HEP.lumen = 1
### Illuminance
#### HEP.lux
### Miscelanea
#### HEP.perCent = 0.01
#### HEP.perThousand = 0.001
#### HEP.perMillion = 0.000001
### Constants
#### HEP.Avogadro
#### HEP.c\_light
#### HEP.c\_squared
#### HEP.h\_Planck
#### HEP.hbar\_Planck
#### HEP.hbarc
#### HEP.hbarc\_squared
#### HEP.electron\_charge
#### HEP.e\_squared
#### HEP.electron\_mass\_c2
#### HEP.proton\_mass\_c2
#### HEP.neutron\_mass\_c2
#### HEP.amu\_c2
#### HEP.amu
#### HEP.mu0
#### HEP.epsilon0
#### HEP.elm\_coupling
#### HEP.fine\_structure\_const
#### HEP.classic\_electr\_radius
#### HEP.electron\_Compton\_length
#### HEP.Bohr\_radius
#### HEP.alpha\_rcl2
#### HEP.twopi\_mc2\_rcl2
#### HEP.k\_Boltzmann
#### HEP.STP\_Temperature
#### HEP.STP\_Pressure
#### HEP.kGasThreshold
#### HEP.universe\_mean\_density
