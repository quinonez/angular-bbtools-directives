/* 
   +----------------------------------------------------------------------+
   |                            BBT Random                                |
   |                           --- RandBit ---                            |
   |                            Module File                               |
   +----------------------------------------------------------------------+
   Module defining methods for shooting Flat or Bit random numbers, 
   double or integers.
   It provides methods to fill with double flat values arrays of
   specified size, as well as methods for shooting sequences of 0,1 (bits).
   +----------------------------------------------------------------------+
   | JavaScript                                                           |
   +----------------------------------------------------------------------+
   F. Quinonez - Created 2015-10-30                  
   +----------------------------------------------------------------------+
   | C++                                                                  |
   +----------------------------------------------------------------------+
   This is derived from RandFlat and is a drop-in replacement.  However
   the shootBit() and fireBit() methods are stateless (which makes them
   an order of magnitude slower, but allows save/restore engine status
   to work correctly).

   M Fischler     - Created from RandFlat.cc, deleting almost all the 
                    content since inheritance takes care of it.  2/15/00
   M Fischler     - put and get to/from streams 12/10/04
*/

  "use strict";
  var JamesRandom = require('jamesrandom');

  function RandBit( args ){
    this.fengine = args.engine || new JamesRandom({});

    this.Fire = function(){
      var x = this.fengine.Flat();
      var bit = ( x > 0.5 )? 1: 0;
      return bit;
    };

    this.FireArray = function( /* size of vect */ size, /* Array */ vect ){
      for( var i = 0; i < size; ++i ){
        vect.push( this.Fire() );
      }
    };
  } 

  RandBit.Shoot = function( args ){
    var sengine = args.engine || new JamesRandom({});
    var x = sengine.Flat();
    var bit = ( x > 0.5 )? 1: 0;
    return bit;
  };

  RandBit.ShootArray = function( args ){
    var ssize = args.size || 1;
    var sengine = args.engine || new JamesRandom({});
    // var svect = args.vect;

    var argsShoot = { engine: sengine };

    for( var i = 0; i < ssize; ++i ){
      args.vect.push( RandBit.Shoot( argsShoot ) );
    }
  };



  module.exports = RandBit;


