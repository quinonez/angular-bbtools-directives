/* 
   +----------------------------------------------------------------------+
   |                            HEP Random                                |
   |                       ---  JamesRandom ---                           |
   |                            Module File                               |
   +----------------------------------------------------------------------+
  
   This algorithm implements the original universal random number generator
   as proposed by Marsaglia & Zaman in report FSU-SCRI-87-50 and coded
   in FORTRAN77 by Fred James as the RANMAR generator, part of the MATHLIB
   HEP library. 

   +----------------------------------------------------------------------+
   | JavaScript                                                           |
   +----------------------------------------------------------------------+
   F. Quinonez - Created 2015-10-19                  
               - Geant4 File
   +----------------------------------------------------------------------+
   | C++                                                                  |
   +----------------------------------------------------------------------+
   Gabriele Cosmo - Created: 5th September 1995
                  - Fixed a bug in setSeed(): 26th February 1996
                  - Minor corrections: 31st October 1996
                  - Added methods for engine status: 19th November 1996
                  - Fixed bug in setSeeds(): 15th September 1997
   J.Marraffino   - Added stream operators and related constructor.
                    Added automatic seed selection from seed table and
                    engine counter: 16th Feb 1998
   Ken Smith      - Added conversion operators:  6th Aug 1998
   J. Marraffino  - Remove dependence on hepString class  13 May 1999
   V. Innocente   - changed pointers to indices     3 may 2000
   M. Fischler    - In restore, checkFile for file not found  03 Dec 2004
   M. Fischler    - Methods for distrib. instacne save/restore  12/8/04    
   M. Fischler    - split get() into tag validation and 
                    getState() for anonymous restores 12/27/04    
   M. Fischler    - Enforcement that seeds be non-negative
  		    (lest the sequence be non-random) 2/14/05    
   M. Fischler    - put/get for vectors of ulongs	3/14/05
   M. Fischler    - State-saving using only ints, for portability 4/12/05
*/

  "use strict";

  // Number of instances with automatic seed selection
  var numberOfEngines = 0;
  // Maximum index into the seed table
  var maxIndex = 215;

  //constructor: JamesRandom,
  function JamesRandom( args ){
    this.u = [];
    this.c;
    this.cd;
    this.cm;
    this.i97;
    this.j97; 

    this.Name = function(){ return "HepJamesRandom"; };
    
    // Returns a pseudo random number between 0 and 1 (excluding the end points)
    this.Flat = function(){
      var uni;

      do {
         uni = this.u[this.i97] - this.u[this.j97];
         if ( uni < 0.0 ) uni++;
         this.u[this.i97] = uni;
         
         if (this.i97 == 0) this.i97 = 96;
         else this.i97--;
         
         if (this.j97 == 0) this.j97 = 96;
         else this.j97--;
         
         this.c -= this.cd;
         if (this.c < 0.0) this.c += this.cm;
         
         uni -= this.c;
         if (uni < 0.0) uni += 1.0;
      } while ( uni <= 0.0 || uni >= 1.0 );
      
      return uni;
         
    };

    // Fills the array "vect" of specified size with flat random values.
    this.FlatArray = function( size, vect ){
      var i;

      for (i=0; i<size; ++i) {
        vect[i] = this.Flat();
      }   

    };

    // Sets the state of the algorithm according to seed.
    this.SetSeed = function( seed, dum ){
      // The input value for "seed" should be within the range [0,900000000]
      //
      // Negative seeds result in serious flaws in the randomness;
      // seeds above 900000000 are OK because of the %177 in the expression for i,
      // but may have the same effect as other seeds below 900000000.

      var m, n, mm; // integers
      var s, t;

      if( seed < 0 ){
        var msg = "Seed for HepJamesRandom must be non-negative. Seed value supplied was" + seed + ". Using its absolute value instead.";
        console.log( msg );
        seed = -seed;
      }

      var ij = seed / 30082;
      var kl = seed - ( 30082 * ij );
      var i = ( ( ij / 177 ) % 177 ) + 2;
      var j = ( ij % 177 ) + 2;
      var k = ( ( kl / 169 ) % 178 ) + 1;
      var l = kl % 169;
      
      this.theSeed = seed;

      for( n = 1 ; n < 98; n++ ) {
        s = 0.0;
        t = 0.5;
        for( m = 1 ; m < 25; m++ ) {
          mm = ( ( (i*j) % 179 ) * k ) % 179;
          i = j;
          j = k;
          k = mm;
          l = ( 53 * l + 1 ) % 169;
          if ( (l*mm % 64 ) >= 32 )
            s += t;
          t *= 0.5;
        }
        this.u.push( s );
      }
      this.c = 362436.0 / 16777216.0;
      this.cd = 7654321.0 / 16777216.0;
      this.cm = 16777213.0 / 16777216.0;
    
      this.i97 = 96;
      this.j97 = 32;
    };

    // Sets the state of the zero terminated array of seeds. Only the first seed is used.
    this.SetSeeds = function( seeds, dum ){
      this.SetSeed( seeds? seeds: 19780503, 0 );
      this.theSeeds = seeds;
    };

    this.theSeed = args.seed || new Date().getDate(); 

    this.SetSeed( this.theSeed, 0 )
    this.SetSeeds( this.theSeed, 0 );
  }

  JamesRandom.MarkerLen = 64; // Enough room to hold a begin or end marker.

  JamesRandom.EngineName = function(){ return "HepJamesRandom"; };
  JamesRandom.BeginTag = function(){ return "JamesRandom-Begin"; };




module.exports = JamesRandom;

