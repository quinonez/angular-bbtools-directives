/* 
   +----------------------------------------------------------------------+
   |                            HEP Random                                |
   |                        --- RandPoisson ---                           |
   |                            Module File                               |
   +----------------------------------------------------------------------+
   Module defining methods for shooting numbers according to the Poisson
   distribution, given a mean (Algorithm taken from "W.H.Press et al.,
   Numerical Recipes in C, Second Edition".
   Default mean value is set to 1.
   +----------------------------------------------------------------------+
   | JavaScript                                                           |
   +----------------------------------------------------------------------+
   F. Quinonez - Created 2015-10-30                  
   +----------------------------------------------------------------------+
   | C++                                                                  |
   +----------------------------------------------------------------------+
   Gabriele Cosmo - Created: 5th September 1995
                  - Added not static Shoot() method: 17th May 1996
                  - Algorithm now operates on doubles: 31st Oct 1996
                  - Added methods to shoot arrays: 28th July 1997
                  - Added check in case xm=-1: 4th February 1998
   J.Marraffino   - Added default mean as attribute and
                    operator() with mean: 16th Feb 1998
   Gabriele Cosmo - Relocated static data from HepRandom: 5th Jan 1999
   M Fischler     - put and get to/from streams 12/15/04
   M Fischler	      - put/get to/from streams uses pairs of ulongs when
  			+ storing doubles avoid problems with precision 
  			4/14/05
   Mark Fischler  - Repaired BUG - when mean > 2 billion, was returning
                    mean instead of the proper value.  01/13/06
*/

  "use strict";
  var JamesRandom = require('jamesrandom');

  function gammln( xx ){
    // Returns the value ln(Gamma(xx) for xx > 0.  Full accuracy is obtained for 
    // xx > 1. For 0 < xx < 1. the reflection formula (6.1.4) can be used first.
    // (Adapted from Numerical Recipes in C)
    var cof = [ 76.18009172947146,
                -86.50532032941677,
                24.01409824083091, 
                -1.231739572450155,
                0.1208650973866179e-2, 
                -0.5395239384953e-5 ];
    var j;
    var x = xx - 1.0;
    var tmp = x + 5.5;
    tmp -= ( x + 0.5 ) * Math.log( tmp );
    var ser = 1.000000000190015;

    for( j = 0; j < cof.length; ++j ){
      x += 1.0;
      ser += cof[ j ] / x;
    }

    return -tmp + Math.log( 2.5066282746310005 * ser );
  };

  function normal( sengine ){
    var r,v1,v2,fac; // double

    do {
      v1 = 2.0 * sengine.Flat() - 1.0;
      v2 = 2.0 * sengine.Flat() - 1.0;
      r = v1 * v1 + v2 * v2;
    } while( r > 1.0 );

    fac = Math.sqrt( -2.0 * Math.log( r ) / r );
    return v2 * fac;
  };




  function RandPoisson( args ){
    this.fmean = args.mean || 1.0;
    this.fstatus = [ 0.0, 0.0, 0.0 ];
    this.foldm = undefined;
    this.fmeanMax;
    this.fengine = args.engine || new JamesRandom({});

    this.Fire = function(){
      // Returns as a floating-point number an integer value that is a random
      // deviation drawn from a Poisson distribution of mean xm, using flat()
      // as a source of uniform random numbers.
      // (Adapted from Numerical Recipes in C)

      var em, t, y;
      var sq, alxm, g1;

      sq = this.fstatus[0];
      alxm = this.fstatus[1];
      g1 = this.fstatus[2];

      if( this.fmean === -1 ) return 0;
      if( this.fmean < 12.0 ){
        if( this.fmean != this.foldm ){
          this.foldm = this.fmean;
          g1 = Math.exp( -this.fmean );
        }
        em = -1;
        t = 1.0;
        do {
          em += 1.0;
          t *= this.fengine.Flat();
        } while( t > g1 );
      } else if ( this.fmean < this.fmeanMax ) {
        if ( this.fmean != this.foldm ) {
          oldm = this.fmean;
          sq = Math.sqrt( 2.0 * this.fmean );
          alxm = Math.log( this.fmean );
          g1 = this.fmean*alxm - gammln( this.fmean + 1.0 );
        }
        do {
          do {
      	    y = Math.tan( Math.PI * this.fengine.Flat() );
            em = sq * y + this.fmean;
          } while( em < 0.0 );
          em = Math.floor( em );
          t = 0.9 * ( 1.0 + y * y ) * Math.exp( em * alxm - gammln( em + 1.0 ) - g1 );
        } while( this.fengine.Flat() > t );
      } else {
        em = this.fmean + Math.sqrt( this.fmean ) * normal( sengine );
        if ( ( em | 0 ) < 0 ) 
          em = ( this.fmean | 0 ) >= 0? this.fmean: GetMaxMean();
      }    
      this.fstatus[0] = sq; 
      this.fstatus[1] = alxm; 
      this.fstatus[2] = g1;

      return ( em | 0 );    
    };

    this.FireArray = function( /* size of vect */ size, /* Array */ vect ){
      for( var i = 0; i < size; ++i ){
        vect.push( this.Fire( this.fmean ) );
      }
    };
  }



  RandPoisson.sstatus = [ undefined, undefined, undefined ];

  RandPoisson.soldMean = undefined;

  RandPoisson.smeanMax = undefined;

  RandPoisson.GetMaxMean = function(){
    return RandPoisson.smeanMax;    
  };

  RandPoisson.GetOldMean = function(){
    return RandPoisson.soldMean;    
  };

  RandPoisson.GetSStatus = function(){
    return RandPoisson.sstatus;    
  };

  RandPoisson.SetOldMean = function( value ){
    RandPoisson.soldMean = value;
  };


  RandPoisson.Shoot = function( args ){
    var smean = args.mean || 1.0;
    var poisson = new RandPoisson( args ); 
    var sengine = args.engine || new JamesRandom({});

    // Returns as a floating-point number an integer value that is a random
    // deviation drawn from a Poisson distribution of mean xm, using flat()
    // as a source of uniform random numbers.
    // (Adapted from Numerical Recipes in C)

    var em, t, y;
    var sq, alxm, g1;
    var om = RandPoisson.GetOldMean();

    sq = poisson.fstatus[0];
    alxm = poisson.fstatus[1];
    g1 = poisson.fstatus[2];

    if( smean === -1 ) return 0;
    if( smean < 12.0 ){
      if( smean != om ){
        RandPoisson.SetOldMean( smean );
        g1 = Math.exp( -smean );
      }
      em = -1;
      t = 1.0;
      do {
        em += 1.0;
        t *= sengine.Flat();
      } while( t > g1 );
    } else if ( smean < RandPoisson.GetMaxMean() ) {
      if ( smean != om ) {
        RandPoisson.SetOldMean( smean );
        sq = Math.sqrt( 2.0 * smean );
        alxm = Math.log( smean );
        g1 = smean * alxm - gammln( smean + 1.0 );
      }
      do {
        do {
    	  y = Math.tan( Math.PI * sengine.Flat() );
          em = sq * y + smean;
        } while( em < 0.0 );
        em = Math.floor( em );
        t = 0.9 * ( 1.0 + y * y ) * Math.exp( em * alxm - gammln( em + 1.0 ) - g1 );
      } while( sengine.Flat() > t );
    } else {
      em = smean + Math.sqrt( smean ) * normal( sengine );
      if ( ( em | 0 ) < 0 ) 
        em = ( smean | 0 ) >= 0? smean: RandPoisson.GetMaxMean();
    }    

    RandPoisson.SetSStatus( sq, alxm, g1 );

    return ( em | 0 );    
  };

  RandPoisson.ShootArray = function( args ){
    var ssize = args.size || 1;
    var smean = args.mean || 1.0;
    var sengine = args.engine || new JamesRandom({});
    // var svect = args.vect;

    var argsShoot = { mean: smean, engine: sengine };

    for( var i = 0; i < ssize; ++i ){
      args.vect.push( RandPoisson.Shoot( argsShoot ) );
    }
  };


  module.exports = RandPoisson;

