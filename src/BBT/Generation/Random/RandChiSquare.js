/* 
   +----------------------------------------------------------------------+
   |                            HEP Random                                |
   |                      --- RandChiSquare ---                           |
   |                            Module File                               |
   +----------------------------------------------------------------------+
   Module defining methods for shooting Chi^2 distributed random values,
   given a number of degrees of freedom a (default=1.0).

   Valid values of a satisfy a > 1. When invalid values are presented,
   the code silently returns -1.0.
   +----------------------------------------------------------------------+
   | JavaScript                                                           |
   +----------------------------------------------------------------------+
   F. Quinonez - Created 2015-10-30                  
   +----------------------------------------------------------------------+
   | C++                                                                  |
   +----------------------------------------------------------------------+
   John Marraffino - Created: 12th May 1998
   M Fischler     - put and get to/from streams 12/10/04
   M Fischler	      - put/get to/from streams uses pairs of ulongs when
  			+ storing doubles avoid problems with precision 
  			4/14/05
*/

  "use strict";
  var JamesRandom = require('jamesrandom');

  function RandChiSquare( args ){
    this.fa = args.a || 1.0;
    this.fengine = args.engine || new JamesRandom({});

    this.Fire = function(){
      return RandChiSquare.GenChiSquare( this.fengine, this.fa );
    };

    this.FireArray = function( /* size of vect */ size, /* Array */ vect ){
      for( var i = 0; i < size; ++i ){
        vect.push( this.Fire() );  
      }
    };
  } 



  RandChiSquare.Shoot = function( args ){
    var sa = args.a || 1.0;
    var sengine = args.engine || new JamesRandom({});
    return RandChiSquare.GenChiSquare( sengine, sa );
  };



  RandChiSquare.ShootArray = function( args ){
    var ssize = args.size || 1;
    // var svect = args.vect;
    var sa = args.a || 1.0;
    var sengine = args.engine || new JamesRandom({});
    var argsShoot = { a: sa, engine: sengine };
    for( var i = 0; i < ssize; ++i ){
      args.vect.push( RandChiSquare.Shoot( argsShoot ) );
    }
  };



  RandChiSquare.GenChiSquare = function( sengine, sa ){
    /******************************************************************
     *                                                                *
     *        Chi Distribution - Ratio of Uniforms  with shift        *
     *                                                                *
     ******************************************************************
     *                                                                *
     * FUNCTION :   - chru samples a random number from the Chi       *
     *                distribution with parameter  a > 1.             *
     * REFERENCE :  - J.F. Monahan (1987): An algorithm for           *
     *                generating chi random variables, ACM Trans.     *
     *                Math. Software 13, 168-172.                     *
     * SUBPROGRAM : - anEngine  ... pointer to a (0,1)-Uniform        *
     *                engine                                          *
     *                                                                *
     * Implemented by R. Kremer, 1990                                 *
     ******************************************************************/
    var a_in = -1.0;
    var b,vm,vp,vd; // CLHEP_THREAD_LOCAL double's
    var u,v,z,zz,r; // double's

    //console.log("a = ",sa);
    // chech for invalid input value
    if( sa < 1 ) return -1.0;

    //var pasada = 0;

    if( sa == 1 ){
      for(;;){
        u = sengine.Flat();
        v = sengine.Flat() * 0.857763884960707;
        z = v / u;
        if( z < 0 ) continue;
        zz = z * z;
        r = 2.5 - zz;
        if( z < 0.0 ) r += zz * z / ( 3.0 * z );
        if( u < r * 0.3894003915 ) return( z * z );
        if( zz > ( 1.036961043 / u + 1.4 ) ) continue;
        if( 2 * Math.log( u ) < ( - zz * 0.5 ) ) return( z * z );
      }
    } else {
      if( sa != a_in ){
        b = Math.sqrt( sa - 1.0 );
        // console.log("b = %f",b);
        vm = - 0.6065306597 * ( 1.0 - 0.25 / ( b * b + 1.0 ) );
        vm = ( -b > vm )? -b: vm;
        vp = 0.6065306597 * (0.7071067812 + b) / (0.5 + b);
        vd = vp - vm;
        a_in = sa;
      }
      for(;;){
        //console.log(pasada);
        u = sengine.Flat();
        v = sengine.Flat() * vd + vm;
        z = v / u;
        // console.log("u: %f",u);
        // console.log("vm: %f",vm);
        // console.log("vp: %f",vp);
        // console.log("v: %f",v);
        // console.log("z: %f",z);
        if( z < -b ) continue;
        zz = z * z;
        r = 2.5 - zz;
        // console.log("zz: %f",zz);
        if( z < 0.0 ) r += zz * z / ( 3.0 * ( z + b ) );
        // console.log("r: %f",r);


        if( u < r * 0.3894003915 ){
          // console.log("A");
          return( ( z + b ) * ( z + b ) );
        }
        if( zz > ( 1.036961043 / u + 1.4 ) ){
          // console.log("B");
          continue;
        }
        if( 2 * Math.log( u ) < ( Math.log( 1.0 + z / b ) * b * b - zz * 0.5 - z * b ) ){
          // console.log("C");
          return ( ( z + b ) * ( z + b ) );
        }
        //pasada++;
        //if( pasada === 10 ) return;
      }       
    }
  };



  module.exports = RandChiSquare;
