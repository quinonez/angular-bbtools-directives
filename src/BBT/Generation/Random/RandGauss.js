/* 
   +----------------------------------------------------------------------+
   |                            HEP Random                                |
   |                          --- RandGauss ---                           |
   |                            Module File                               |
   +----------------------------------------------------------------------+
   Module defining methods for shooting gaussian distributed random values,
   given a mean (default=0) or specifying also a deviation (default=1).
   Gaussian random numbers are generated two at the time, so every
   other time shoot is called the number returned is the one generated the
   time before.
   +----------------------------------------------------------------------+
   | JavaScript                                                           |
   +----------------------------------------------------------------------+
   F. Quinonez - Created 2015-10-30                  
   +----------------------------------------------------------------------+
   | C++                                                                  |
   +----------------------------------------------------------------------+
   Gabriele Cosmo - Created: 5th September 1995
                  - Added methods to shoot arrays: 28th July 1997
   J.Marraffino   - Added default arguments as attributes and
                    operator() with arguments. Introduced method normal()
                    for computation in fire(): 16th Feb 1998
   Gabriele Cosmo - Relocated static data from HepRandom: 5th Jan 1999
   M Fischler     - Copy constructor should supply right engine to HepRandom:
  		    1/26/00.
   M Fischler     - Workaround for problem of non-reproducing saveEngineStatus
  		    by saving cached gaussian.  March 2000.
   M Fischler     - Avoiding hang when file not found in restoreEngineStatus 
                    12/3/04
   M Fischler     - put and get to/from streams 12/8/04
   M Fischler     - save and restore dist to streams 12/20/04
   M Fischler	  - put/get to/from streams uses pairs of ulongs when
  		    storing doubles avoid problems with precision.
  		    Similarly for saveEngineStatus and RestoreEngineStatus
  		    and for save/restore distState
  		    Care was taken that old-form output can still be read back.
  			4/14/05
*/

  "use strict";
  var JamesRandom = require('jamesrandom');

  function RandGauss( args ){
    this.fmean = args.mean || 0.0;
    this.fstdDev = args.stdDev || 1.0;
    this.fengine = args.engine || new JamesRandom({});
    this.fset = false;
    this.fnextGauss;

    this.Fire = function(){
      return this.Normal * this.fstdDev + this.fmean;
    };

    this.FireArray = function( /* size of vect */ size, /* Array */ vect ){
      for( var i = 0; i < size; ++i ){
        vect.push( this.Fire() );
      }
    };

    this.Normal = function(){
      if( this.fset ){
        this.fset = false;
        return this.fnextGauss;
      }

      var r,v1,v2,fac,val; // double

      do {
        v1 = 2.0 * this.fengine.Flat() - 1.0;
        v2 = 2.0 * this.fengine.Flat() - 1.0;
        r = v1 * v1 + v2 * v2;
      } while( r > 1.0 );

      fac = Math.sqrt( -2.0 * Math.log( r ) / r );
      val = v1 * fac;
      this.fnextGauss = val;
      this.fset = true;
      return v2 * fac;
    };
  } 

  RandGauss.snextGauss = 0.0;
  RandGauss.sset = false;

  RandGauss.GetFlag = function(){
    return RandGauss.sset;
  };

  RandGauss.GetVal = function(){ 
    return RandGauss.snextGauss; 
  };

  RandGauss.SetFlag = function( value ){
    RandGauss.sset = value; 
  };

  RandGauss.SetVal = function( value ){ 
    RandGauss.snextGauss = value; 
  };

  RandGauss.ShootAux = function( sengine ){
    // Gaussian random numbers are generated two at the time, so every other
    // time this is called we just return a number generated the time before.
    if ( RandGauss.GetFlag() ) {
      RandGauss.SetFlag( false );
      var x = RandGauss.GetVal();
      return x; 
    } 

    var r,v1,v2,fac,val; // double

    do {
      v1 = 2.0 * sengine.Flat() - 1.0;
      v2 = 2.0 * sengine.Flat() - 1.0;
      r = v1 * v1 + v2 * v2;
    } while( r > 1.0 );

    fac = Math.sqrt( -2.0 * Math.log( r ) / r );
    val = v1 * fac;
    RandGauss.SetVal( val );
    RandGauss.SetFlag( true );
    return ( v2 * fac );
  };

  RandGauss.Shoot = function( args ){
    var smean = args.mean || 0.0;
    var sstdDev = args.stdDev || 1.0;
    var sengine = args.engine || new JamesRandom({});
    var stdValue = RandGauss.ShootAux( sengine );
    console.log(stdValue*sstdDev+smean);
    return  stdValue * sstdDev + smean;
  };

  RandGauss.ShootArray = function( args ){
    var ssize = args.size || 1;
    var smean = args.mean || 0.0;
    console.log( "mean: %f", smean );
    var sstdDev = args.stdDev || 1.0;
    console.log( "stdDev: %f", sstdDev );
    var sengine = args.engine || new JamesRandom({});
    // var svect = args.vect;

    var argsShoot = { mean: smean, stdDev: sstdDev, engine: sengine };

    for( var i = 0; i < ssize; ++i ){
      args.vect.push( RandGauss.Shoot( argsShoot ) );
    }

  };

  module.exports = RandGauss;
