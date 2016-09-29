/* 
   +----------------------------------------------------------------------+
   |                            HEP Random                                |
   |                       --- RandBreitWigner ---                        |
   |                            Module File                               |
   +----------------------------------------------------------------------+
   Module defining methods for shooting numbers according to the
   Breit-Wigner distribution algorithms (plain or mean^2).
   Default values are set: mean=1, gamma=.2, cut=undefined.
   Plain algorithm is used for shootArray() and fireArray().
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
                    operator() with arguments: 16th Feb 1998
   M Fischler      - put and get to/from streams 12/10/04
*/
  "use strict";
  var JamesRandom = require('jamesrandom');


  function RandBreitWigner( args ){
    this.fmean = args.mean || 1.0;
    this.fgamma = args.gamma || 0.2;
    this.fcut = args.cut || undefined;
    this.fengine = args.engine || new JamesRandom({});

    this.Fire = function(){
      if( this.fgamma == 0 ) return this.fmean;
      var rval, displ;
      rval  = 2.0 * this.fengine.Flat() - 1.0;

      if( this.fcut === undefined ){
        // Without cut
        displ = 0.5 * this.fgamma * Math.tan( rval * Math.PI * 0.5 ); 
      } else {
        // With cut
        var val;
        val = Math.atan( 2.0 * this.fcut / this.fgamma );        
        displ = 0.5 * gamma * Math.tan( rval * val ); 
      }
      return ( this.fmean + displ );
    };

    this.FireM2 = function(){
      var val, rval, displ;
      if( this.fgamma == 0.0 ) return this.fmean;

      if( this.fcut === undefined ){
        // Without cut
        val = Math.atan( -this.fmean / this.fgamma );
        rval = RandBreitWigner.ShootFlat( {  a: val, b: Math.PI / 2, engine: this.fengine  } );
        displ = this.fgamma * Math.tan( rval );
        return Math.sqrt( Math.pow( this.fmean, 2 ) + this.fmean * displ );          
      } else {
        // With cut
        var lower, upper, tmp;

        tmp = Math.max( 0.0, this.fmean - this.fcut );
        lower = Math.atan( ( tmp * tmp - Math.pow( this.fmean, 2) ) / ( this.fmean * this.fgamma ) );

        upper = Math.atan( ( Math.pow( this.fmean + this.fcut, 2 ) - Math.pow( this.fmean, 2 ) ) / ( this.fmean * this.fgamma ) );

        rval = RandBreitWigner.ShootFlat( { a: lower, b: upper, engine: this.fengine } );

        displ = this.fgamma * Math.tan( rval );

        return Math.sqrt( Math.max( 0.0, Math.pow( this.fmean, 2 ) + this.fmean * displ ) );
      }
    };

    this.FireArray = function( /* size of vect */ size, /* Array */ vect ){
      for( var i = 0; i < size; ++i ){
        vect.push( this.Fire() );
      }
    };
  } 

  RandBreitWigner.Shoot = function( args ){
    var smean = args.mean || 1.0;
    var sgamma = args.gamma || 0.2;
    var scut = args.cut || undefined;
    var sengine = args.engine || new JamesRandom({});

    if( sgamma == 0 ) return smean;
    
    var rval, displ;
    rval  = 2.0 * sengine.Flat() - 1.0;
    if( scut === undefined ){
      // Without cut
      displ = 0.5 * sgamma * Math.tan( rval * Math.PI * 0.5 ); 
    } else {
      // With cut
      var val;
      val = Math.atan( 2.0 * scut / sgamma );        
      displ = 0.5 * gamma * Math.tan( rval * val ); 
    }
    return ( smean + displ );
  };

  RandBreitWigner.ShootArray = function( args ){
    var ssize = args.size || 1;
    var smean = args.mean || 1.0;
    var sgamma = args.gamma || 0.2;
    var scut = args.cut || undefined;
    var sengine = args.engine || new JamesRandom({});
    // var svect = args.vect;

    var argsShoot = { mean: smean, gamma: sgamma, cut: scut, engine: sengine };

    for( var i = 0; i < ssize; ++i ){
      args.vect.push( RandBreitWigner.Shoot( argsShoot ) );
    }
    
  };

  // Tweak CLHEP C++ that used the RandFlat::shoot static method.
  // Now this distribution has its own Shoot static function.  
  RandBreitWigner.ShootFlat = function( args ){
    var sa = args.a || 0;
    var sb = args.b || 1;
    var swidth = args.width || ( args.b - args.a );  
    var sengine = args.engine || new JamesRandom({});
    return ( swidth * sengine.Flat() + sa ); 
  };

  RandBreitWigner.ShootM2 = function( args ){
    var sengine = args.engine || new JamesRandom({});
    var smean = args.mean || 1.0;
    var sgamma = args.gamma || 0.2;
    var scut = args.cut || undefined;

    var val, rval, displ;
    if( sgamma == 0.0 ) return smean;

    if( scut === undefined ){
      // Without cut
      val = Math.atan( -smean / sgamma );
      rval = RandBreitWigner.ShootFlat( {  a: val, b: Math.PI / 2, engine: sengine  } );
      displ = sgamma * Math.tan( rval );
      return Math.sqrt( Math.pow( smean, 2 ) + smean * displ );          
    } else {
      // With cut
      var lower, upper, tmp;

      tmp = Math.max( 0.0, smean - scut );
      lower = Math.atan( ( tmp * tmp - Math.pow( smean, 2) ) / ( smean * sgamma ) );

      upper = Math.atan( ( Math.pow( smean + scut, 2 ) - Math.pow( smean, 2 ) ) / ( smean * sgamma ) );

      rval = RandBreitWigner.ShootFlat( { a: lower, b: upper, engine: sengine } );

      displ = sgamma * Math.tan( rval );

      return Math.sqrt( Math.max( 0.0, Math.pow( smean, 2 ) + smean * displ ) );
    }
  };


  module.exports = RandBreitWigner;
