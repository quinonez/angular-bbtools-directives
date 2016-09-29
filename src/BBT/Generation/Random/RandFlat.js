/* 
   +----------------------------------------------------------------------+
   |                            BBT Random                                |
   |                          --- RandFlat ---                            |
   |                            Module File                               |
   +----------------------------------------------------------------------+
   Module defining methods for shooting flat random numbers, double or
   integers.
   It provides methods to fill with double flat values arrays of
   specified size, as well as methods for shooting sequences of 0,1 (bits).
   +----------------------------------------------------------------------+
   | JavaScript                                                           |
   +----------------------------------------------------------------------+
   F. Quinonez - Created 2015-10-30                  
   +----------------------------------------------------------------------+
   | C++                                                                  |
   +----------------------------------------------------------------------+
   Gabriele Cosmo - Created: 5th September 1995
   Peter Urban    - ShootBit() and related stuff added: 5th Sep 1996
   Gabriele Cosmo - Added operator() and additional methods to fill
                    arrays specifying boundaries: 24th Jul 1997 
   J.Marraffino   - Added default arguments as attributes and
                    operator() with arguments: 16th Feb 1998
   M. Fischler    - Moved copy constructor to protected so that
  		    Rderived RandBit can get at it.
   M Fischler      - put and get to/from streams 12/10/04
*/

  "use strict";
  var JamesRandom = require('jamesrandom');

  function RandFlat( args ){
    this.fa = args.a || 0;
    this.fb = args.b || 1;
    this.fwidth = args.width || ( args.b - args.a );
    this.fengine = args.engine || new JamesRandom({});

    this.Fire = function(){
      return ( this.fb - this.fa ) * fengine.Flat() + this.fa;
    };

    this.FireArray = function( /* size of vect */ size, /* Array */ vect ){
      for( var i = 0; i < size; ++i ){
        vect.push( this.Fire() );
      }
    };
  } 


  RandFlat.Shoot = function( args ){
    var sa = args.a || 0;
    var sb = args.b || 1;
    var swidth = args.width || ( args.b - args.a );  
    var sengine = args.engine || new JamesRandom({});
    return ( swidth * sengine.Flat() + sa ); 
  };

  RandFlat.ShootArray = function( args ){
    var ssize = args.size || 1;
    var sa = args.a || 0;
    var sb = args.b || 1;
    var swidth = args.width || ( args.b - args.a );  
    var sengine = args.engine || new JamesRandom({});
    // var svect = args.vect;

    var argsShoot = { a: sa, b: sb, width: swidth, engine: sengine };

    for( var i = 0; i < ssize; ++i ){
      args.vect.push( RandFlat.Shoot( argsShoot ) );
    }

  };


  module.exports = RandFlat;

