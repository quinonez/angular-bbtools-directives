/* 
   +----------------------------------------------------------------------+
   |                            HEP Random                                |
   |                       --- RandStudentT ---                           |
   |                            Module File                               |
   +----------------------------------------------------------------------+
   Module defining methods for shooting Student's t- distributed random 
   values, given a number of degrees of freedom a (default=1.0).
  
   Valid input values are a > 0.  When invalid values are presented, the
   code silently returns Number.MAX_VALUE.
   +----------------------------------------------------------------------+
   | JavaScript                                                           |
   +----------------------------------------------------------------------+
   F. Quinonez - Created 2015-10-30                  
   +----------------------------------------------------------------------+
   | C++                                                                  |
   +----------------------------------------------------------------------+
   John Marraffino - Created: 12th May 1998
   G.Cosmo         - Fixed minor bug on inline definition for shoot()
                     methods : 20th Aug 1998
   M Fischler      - put and get to/from streams 12/13/04
   M Fischler	      - put/get to/from streams uses pairs of ulongs when
  			+ storing doubles avoid problems with precision 
  			4/14/05
*/

  "use strict";
  var JamesRandom = require('jamesrandom');



  function RandStudentT( args ){
    this.fa = args.a || 1.0;
    this.fengine = args.engine || new JamesRandom({});

    this.Fire = function(){

      var u,v,w;

      do {
        u = 2.0 * this.fengine.Flat() - 1.0;
        v = 2.0 * this.fengine.Flat() - 1.0;
      } while( ( w = u * u + v * v) > 1.0 );

      return ( u * Math.sqrt( this.fa * ( Math.exp( -2.0 / this.fa * Math.log( w ) ) - 1.0 ) / w ) );
    };

    this.FireArray = function( /* size of vect */ size, /* Array */ vect ){
      for( var i = 0; i < size; i++ ){
        vect.push( this.Fire() );  
      }
    };

  } 



  RandStudentT.Shoot = function( args ){
    /******************************************************************
     *                                                                *
     *           Student-t Distribution - Polar Method                *
     *                                                                *
     ******************************************************************
     *                                                                *
     * The polar method of Box/Muller for generating Normal variates  *
     * is adapted to the Student-t distribution. The two generated    *
     * variates are not independent and the expected no. of uniforms  *
     * per variate is 2.5464.                                         *
     *                                                                *
     ******************************************************************
     *                                                                *
     * FUNCTION :   - tpol  samples a random number from the          *
     *                Student-t distribution with parameter a > 0.    *
     * REFERENCE :  - R.W. Bailey (1994): Polar generation of random  *
     *                variates with the t-distribution, Mathematics   *
     *                of Computation 62, 779-781.                     *
     * SUBPROGRAM : -  ... (0,1)-Uniform generator                    *
     *                                                                *
     *                                                                *
     * Implemented by F. Niederl, 1994                                *
     ******************************************************************/
 
    var sa = args.a;
    var sengine = args.engine || new JamesRandom({});

    var u,v,w;

    // check for valid input value
    if( sa < 0.0 ) return Number.MAX_VALUE;

    do {
      u = 2.0 * sengine.Flat() - 1.0;
      v = 2.0 * sengine.Flat() - 1.0;
    } while( ( w = u * u + v * v ) > 1.0 );

      return ( u * Math.sqrt( sa * ( Math.exp( -2.0 / sa * Math.log( w ) ) - 1.0 ) / w ) );

  };


  RandStudentT.ShootArray = function( args ){
    var ssize = args.size || 1;
    var sa = args.a || 1.0;
    var sengine = args.engine || new JamesRandom({});
    // var svect = args.vect;

    var argsShoot = { a: sa, engine: sengine };

    for( var i = 0; i < ssize; ++i ){
      args.vect.push( RandStudentT.Shoot( argsShoot ) );
    }
  };


  module.exports = RandStudentT;
