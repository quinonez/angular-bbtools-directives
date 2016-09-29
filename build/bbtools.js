  "use strict";

  function BBT(){

    this.JamesRandom = require('jamesrandom');
    this.HEP         = require('hephys');

    this.RandBinomial    = require('./BBT/Generation/Random/RandBinomial.js'); 
    this.RandBit         = require('./BBT/Generation/Random/RandBit.js'); 
    this.RandBreitWigner = require('./BBT/Generation/Random/RandBreitWigner.js'); 
    this.RandChiSquare   = require('./BBT/Generation/Random/RandChiSquare.js'); 
    this.RandExponential = require('./BBT/Generation/Random/RandExponential.js'); 
    this.RandFlat        = require('./BBT/Generation/Random/RandFlat.js'); 
    this.RandGamma       = require('./BBT/Generation/Random/RandGamma.js'); 
    this.RandGauss       = require('./BBT/Generation/Random/RandGauss.js'); 
    this.RandLandau      = require('./BBT/Generation/Random/RandLandau.js'); 
    this.RandPoisson     = require('./BBT/Generation/Random/RandPoisson.js'); 
    this.RandStudentT    = require('./BBT/Generation/Random/RandStudentT.js'); 


    this.H1 = require('./BBT/Visualization/H1.js'); 



  }




  module.exports = BBT;
;/* 
   +----------------------------------------------------------------------+
   |                            BBT Random                                |
   |                       --- RandBinomial ---                           |
   |                            Module File                               |
   +----------------------------------------------------------------------+
  
   Module defining methods for shooting binomial distributed random values,
   given a sample size n (default=1) and a probability p (default=0.5).
  
   Valid input values satisfy the relation n*min(p,1-p) > 0. When invalid
   values are presented, the code silently returns -1.0.
  
   +----------------------------------------------------------------------+
   | JavaScript                                                           |
   +----------------------------------------------------------------------+
   F. Quinonez - Created 2015-10-30                  
   +----------------------------------------------------------------------+
   | C++                                                                  |
   +----------------------------------------------------------------------+
   John Marraffino - Created: 12th May 1998  Based on the C-Rand package
                     by Ernst Stadlober and Franz Niederl of the Technical
                     University of Graz, Austria.
   Gabriele Cosmo  - Removed useless methods and data: 5th Jan 1999
   M Fischler      - put and get to/from streams 12/10/04
*/

  "use strict";
  var JamesRandom = require('jamesrandom');

  function StirlingCorrection( k ){
    /*
    +----------------------------------------------------------------------+
    | StirlingCorrection                                                   |
    +----------------------------------------------------------------------+

    Correction term of the Stirling approximation for std::log(k!)          
    (series in 1/k, or table values for small k)                         
    with long int parameter k                                            
                                                                         
                                                                         
    log k! = (k + 1/2)log(k + 1) - (k + 1) + (1/2)log(2Pi) +              
             StirlingCorrection(k + 1)                                    
                                                                         
    log k! = (k + 1/2)log(k)     -  k      + (1/2)log(2Pi) +              
             StirlingCorrection(k)                                        
                                                                         
    */
    var C1 = 8.33333333333333333e-02;     //  +1/12 
    var C3 = -2.77777777777777778e-03;     //  -1/360
    var C5 = 7.93650793650793651e-04;     //  +1/1260
    var C7 = -5.95238095238095238e-04;     //  -1/1680
 
    var  c = new Array(31);
    c = [   0.0,
 			     8.106146679532726e-02, 4.134069595540929e-02,
 			     2.767792568499834e-02, 2.079067210376509e-02,
 			     1.664469118982119e-02, 1.387612882307075e-02,
 			     1.189670994589177e-02, 1.041126526197209e-02,
 			     9.255462182712733e-03, 8.330563433362871e-03,
 			     7.573675487951841e-03, 6.942840107209530e-03,
 			     6.408994188004207e-03, 5.951370112758848e-03,
 			     5.554733551962801e-03, 5.207655919609640e-03,
 			     4.901395948434738e-03, 4.629153749334029e-03,
 			     4.385560249232324e-03, 4.166319691996922e-03,
 			     3.967954218640860e-03, 3.787618068444430e-03,
 			     3.622960224683090e-03, 3.472021382978770e-03,
 			     3.333155636728090e-03, 3.204970228055040e-03,
 			     3.086278682608780e-03, 2.976063983550410e-03,
 			     2.873449362352470e-03, 2.777674929752690e-03,
    ];
    var r, rr;
 
    if( k > 30 ){
      r = 1.0 / k;
      rr = r * r;
      return( r * ( C1 + rr * ( C3 + rr * ( C5 + rr * C7 ) ) ) );
    } else return( c[ k ] );
  };



  function RandBinomial( args ){
    this.fn = args.n || 1;
    this.fp = args.p || 0.5;
    this.fengine = args.engine || new JamesRandom({});

    this.Fire = function( ){
      return RandBinomial.GenBinomial( this.fengine, this.fn, this.fp );    
    };

    this.FireArray = function( /* size of vect */ size, /* Array */ vect ){
      for( var i = 0; i < size; i++ ){
        vect.push( this.Fire() );  
      }
    };
  }


  RandBinomial.GenBinomial = function( sengine, n, p ){
    /*
    +----------------------------------------------------------------------+
    |     Binomial-Distribution - Acceptance Rejection/Inversion           |
    +----------------------------------------------------------------------+
                                                                     
    Acceptance Rejection method combined with Inversion for generating 
    Binomial random numbers with parameters n (number of trials) and p 
    (probability of success).           
    
    For  min(n*p,n*(1-p)) < 10  the Inversion method is applied: The random 
    numbers are generated via sequential search, starting at the lowest 
    index k=0. 
    The cumulative probabilities are avoided by using the technique of 
    chop-down.               

    For  min(n*p,n*(1-p)) >= 10  Acceptance Rejection is used: The algorithm
    is based on a hat-function which is uniform in the centre region and 
    exponential in the tails.                
    A triangular immediate acceptance region in the centre speeds up the 
    generation of binomial variates.                       
    If candidate k is near the mode, f(k) is computed recursively starting 
    at the mode m.                                        
    The acceptance test by Stirling's formula is modified according to 

    W. Hoermann (1992): The generation of binomial    
    random variates, to appear in J. Statist. Comput. Simul.       

    If  p < .5  the algorithm is applied to parameters n, p.       
    Otherwise p is replaced by 1-p, and k is replaced by n - k.    
                                                                   
                                                                   
    FUNCTION:    - btpec samples a random number from the binomial 
                   distribution with parameters n and p  and is    
                   valid for  n*min(p,1-p)  >  0.                  
    REFERENCE:   - V. Kachitvichyanukul, B.W. Schmeiser (1988):    
                   Binomial random variate generation,             
                   Communications of the ACM 31, 216-222.          
    SUBPROGRAMS: - StirlingCorrection()                            
                               ... Correction term of the Stirling 
                                   approximation for std::log(k!)  
                                   (series in 1/k or table values  
                                   for small k) with long int k    
                 - sengine    ... Pointer to a (0,1)-Uniform       
                                   engine                          
                                                                   
    Implemented by H. Zechner and P. Busswald, September 1992      

    */
    
    var C1_3 = 0.33333333333333333;
    var C5_8 = 0.62500000000000000;
    var C1_6 = 0.16666666666666667;
    var DMAX_KM = 20;

    
    var n_last = -1, n_prev = -1; // static CLHEP_THREAD_LOCAL long int
    var p_last = -1.0, p_prev = -1.0; // static CLHEP_THREAD_LOCAL double

    var par, np, p0, q, pq, rc, ss, xm, xl, xr, ll, lr, c, p1, p2, p3, p4, ch; // static CLHEP_THREAD_LOCAL double
    var f, rm, U, V, X, T, E; // double
    var b, m, nm; // static CLHEP_THREAD_LOCAL long 

    var bh, i, K, Km, nK; // long

    if( n != n_last || p != p_last ){
      // set-up 
      n_last = n;
      p_last = p;
      par = Math.min( p, 1.0 - p );
      q = 1.0 - par;
      np = n * par;
    
      // Check for invalid input values
    
      if( np <= 0.0 ) return ( -1.0 );
    
      rm = np + par;
      m  = rm;                      // mode, integer 
      if( np < 10 ){
    	p0 = Math.exp( n * Math.log( q ) );       // Chop-down
    	bh =  np + 10.0 * Math.sqrt( np * q );
    	b = Math.min( n, bh );
      } else {
    	rc = (n + 1.0) * ( pq = par / q );          // recurr. relat.
    	ss = np * q;                              // variance  
    	i  = 2.195 * Math.sqrt( ss ) - 4.6 * q ; // i = p1 - 0.5
    	xm = m + 0.5;
    	xl = m - i;                    // limit left 
    	xr = m + i + 1;               // limit right
    	f  = ( rm - xl ) / ( rm - xl * par );  
        ll = f * ( 1.0 + 0.5 * f );
    	f  = ( xr - rm ) / ( xr * q );     
        lr = f * ( 1.0 + 0.5 * f );
    	c  = 0.134 + 20.5 / ( 15.3 + m );    // parallelogram
    						  // height
    	p1 = i + 0.5;
    	p2 = p1 * ( 1.0 + c + c );                  // probabilities
    	p3 = p2 + c / ll;                           // of regions 1-4
    	p4 = p3 + c / lr;
      }
    }
    if( np <= 0.0 ) return ( -1.0 );
    if( np < 10 ){                                      
      //Inversion Chop-down
      var pk;
    
      K = 0;
      pk = p0;
      U = sengine.Flat();

      while( U > pk ){
    	++K;
    	if( K > b ){
    	  U = sengine.Flat();
    	  K = 0;
    	  pk = p0;
    	} else {
    	  U -= pk;
    	  pk = ( ( n - K + 1 ) * par * pk ) / ( K * q );
    	}
      }

      return ( ( p > 0.5 )? ( n - K ): K );
    } 
    
    for( ; ; ){
      V = sengine.Flat();
      if( ( U = sengine.Flat() * p4 ) <= p1 ){
        // triangular region
    	K = xm - U + p1 * V;
    	return ( ( p > 0.5)? ( n - K ): K );  // immediate accept
      }
      if (U <= p2){
        // parallelogram
        X = xl + ( U - p1 ) / c;
    	if( ( V = V * c + 1.0 - Math.abs( xm - X ) / p1 ) >= 1.0 ) continue;
        K = X;
      } else if (U <= p3){
        // left tail
        if( ( X = xl + Math.log( V ) / ll ) < 0.0 ) continue;
          K = X;
    	  V *= ( U - p2 ) * ll;
      } else {
        // right tail
        if( ( K = xr - Math.log( V ) / lr ) > n )  continue;
        V *= (U - p3) * lr;
      }
    
      // acceptance test :  two cases, depending on |K - m|
      // console.log(abs( K - m ));
      if( ( ( Km = Math.abs( K - m ) ) <= DMAX_KM )  || ( (Km + Km + 2) >= ss ) ){
        // computation of p(K) via recurrence relationship from the mode
    	f = 1.0; // f(m)
    	if( m < K ){
    	  for( i = m; i < K; ){
            if( ( f *= ( rc / ++i - pq ) ) < V )  break; // multiply  f
    	  }
    	} else {
    	  for( i = K; i < m; ){
    	    if( ( V *= ( rc / ++i - pq ) ) > f )  break; // multiply  V
          }
    	}
    	if (V <= f)  break; // acceptance test
      } else {
        // lower and upper squeeze tests, based on lower bounds for log p(K)
  	V = Math.log( V );
    	T = - Km * Km / ( ss + ss );
    	E = ( Km / ss ) * ( ( Km * ( Km * C1_3 + C5_8 ) + C1_6 ) / ss + 0.5 );
        if( V <= T - E ) break;
    	if( V <= T + E ){
    	  if( n != n_prev || par != p_prev ){
    	    n_prev = n;
    	    p_prev = par;  
    	    nm = n - m + 1;
    	    ch = xm * Math.log( ( m + 1.0 ) / ( pq * nm ) ) + StirlingCorrection( m + 1 ) + StirlingCorrection( nm );
    	  } 
    	  nK = n - K + 1;
    
          // computation of log f(K) via Stirling's formula
          // final acceptance-rejection test
    	  if( V <= ( ch  
                       + ( n + 1.0 ) * Math.log( nm / nK )
                       + ( K + 0.5 ) * Math.log( nK * pq / ( K + 1.0 ) )
                       - StirlingCorrection( K + 1 ) 
                       - StirlingCorrection( nK ) 
                     ) )  break;
    	} 
      }
    } // end for 
    return ( ( p > 0.5 )? ( n - K ): K );
  };


  RandBinomial.Shoot = function( args ){
    var sn = args.n || 1;
    var sp = args.p || 0.5;
    var sengine = args.engine || new JamesRandom({});
    return RandBinomial.GenBinomial( sengine, sn, sp );    
  };

  RandBinomial.ShootArray = function( args ){
    var ssize = args.size || 1;
    var sn = args.n || 1;
    var sp = args.p || 0.5;
    var sengine = args.engine || new JamesRandom({});
    // var svect = args.vect;

    var argsShoot = { n: sn, p: sp, engine: sengine };

    for( var i = 0; i < ssize; ++i ){
      args.vect.push( RandBinomial.Shoot( argsShoot ) );
    }
  };


  module.exports = RandBinomial;;

;/* 
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


;/* 
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
;/* 
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
;/* 
   +----------------------------------------------------------------------+
   |                            HEP Random                                |
   |                      --- RandExponential ---                         |
   |                            Module File                               |
   +----------------------------------------------------------------------+
   Class defining methods for shooting exponential distributed random
   values, given a mean (default mean = 1).
   +----------------------------------------------------------------------+
   | JavaScript                                                           |
   +----------------------------------------------------------------------+
   F. Quinonez - Created 2015-10-30                  
   +----------------------------------------------------------------------+
   | C++                                                                  |
   +----------------------------------------------------------------------+
   Gabriele Cosmo - Created: 17th May 1996
                  - Added methods to shoot arrays: 28th July 1997
   J.Marraffino   - Added default mean as attribute and
                    operator() with mean: 16th Feb 1998
   M Fischler      - put and get to/from streams 12/15/04
   M Fischler	      - put/get to/from streams uses pairs of ulongs when
  			+ storing doubles avoid problems with precision 
  			4/14/05
*/

  "use strict";
  var JamesRandom = require('jamesrandom');

  function RandExponential( args ){
    this.fmean = args.mean || 1.0;
    this.fengine = args.engine || new JamesRandom({});
    //this.fnextGauss;
    //this.set;

    this.Fire = function( ){
      return -Math.log( this.fengine.Flat() ) * this.fmean;
    };

    this.FireArray = function( /* size of vect */ size, /* Array */ vect ){
      for( var i = 0; i < size; ++i ){
        vect.push( this.Fire() );
      }
    };
  } 




  RandExponential.Shoot = function( args ){
    var smean = args.mean || 1.0;
    var sengine = args.engine || new JamesRandom({});
    return -Math.log( sengine.Flat() ) * smean; 
  };




  RandExponential.ShootArray = function( args ){
    var ssize = args.size || 1;
    var smean = args.mean || 1.0;
    var sengine = args.engine || new JamesRandom({});
    // var svect = args.vect;

    var argsShoot = { mean: smean, engine: sengine };

    for( var i = 0; i < ssize; ++i ){
      args.vect.push( RandExponential.Shoot( argsShoot ) );
    }

  };



  module.exports = RandExponential;

;/* 
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

;/* 
   +----------------------------------------------------------------------+
   |                            HEP Random                                |
   |                         --- RandGamma ---                            |
   |                            Module File                               |
   +----------------------------------------------------------------------+
   Module defining methods for shooting gamma distributed random values,
   given a k (default=1) and specifying also a lambda (default=1).
   Default values are used for operator()().
  
   Valid input values are k > 0 and lambda > 0.  When invalid values are
   presented, the code silently returns -1.0.
   +----------------------------------------------------------------------+
   | JavaScript                                                           |
   +----------------------------------------------------------------------+
   F. Quinonez - Created 2015-10-30                  
   +----------------------------------------------------------------------+
   | C++                                                                  |
   +----------------------------------------------------------------------+
   John Marraffino - Created: 12th May 1998
   M Fischler      - put and get to/from streams 12/13/04
   M Fischler	      - put/get to/from streams uses pairs of ulongs when
  			+ storing doubles avoid problems with precision 
  			4/14/05
*/

  "use strict";
  var JamesRandom = require('jamesrandom');

  function RandGamma( args ){
    this.fk = args.k || 1.0;
    this.flambda = args.lambda || 1.0;
    this.fengine = args.engine || new JamesRandom({});

    this.Fire = function(){
      return RandGamma.GenGamma( this.fengine, this.fk, this.flambda );
    };

    this.FireArray = function( /* size of vect */ size, /* Array */ vect ){
      for( var i = 0; i < size; ++i ){
        vect.push( this.Fire() );
      }
    };
  } 

  RandGamma.GenGamma = function( sengine, sk, slambda ){
    /*************************************************************************
    *         Gamma Distribution - Rejection algorithm gs combined with     *
    *                              Acceptance complement method gd          *
    *************************************************************************/

    var aa = -1.0, aaa = -1.0, b, c, d, e, r, s, si, ss, q0;
    var q1 = 0.0416666664, q2 =  0.0208333723, q3 = 0.0079849875,
        q4 = 0.0015746717, q5 = -0.0003349403, q6 = 0.0003340332,
        q7 = 0.0006053049, q8 = -0.0004701849, q9 = 0.0001710320,
        a1 = 0.333333333,  a2 = -0.249999949,  a3 = 0.199999867,
        a4 =-0.166677482,  a5 =  0.142873973,  a6 =-0.124385581,
        a7 = 0.110368310,  a8 = -0.112750886,  a9 = 0.104089866,
        e1 = 1.000000000,  e2 =  0.499999994,  e3 = 0.166666848,
        e4 = 0.041664508,  e5 =  0.008345522,  e6 = 0.001353826,
        e7 = 0.000247453;

    var gds, p, q, t, sign_u, u, v, w, x;
    var v1, v2, v12;

    // Check for invalid input values

    if( sk <= 0.0 ) return (-1.0);
    if( slambda <= 0.0 ) return (-1.0);

    if ( sk < 1.0 ){
      // CASE A: Acceptance rejection algorithm gs
      b = 1.0 + 0.36788794412 * sk;       // Step 1
      for(;;){
        p = b * sengine.Flat();
        if( p <= 1.0 ){                            // Step 2. Case gds <= 1
	  gds = Math.exp( Math.log( p ) / sk );
	  if ( Math.log( sengine.Flat() ) <= -gds) return( gds / slambda );
	} else {                            // Step 3. Case gds > 1
	  gds = - Math.log ( ( b - p ) / sk );
	  if( Math.log( sengine.Flat() ) <= ( ( sk - 1.0 ) * Math.log( gds ) ) ) return( gds / slambda );
	}
      }
    } else {          // CASE B: Acceptance complement algorithm gd
      if( sk != aa ){                               // Step 1. Preparations
        aa = sk;
        ss = sk - 0.5;
        s = Math.sqrt( ss );
        d = 5.656854249 - 12.0 * s;
      }
                                              // Step 2. Normal deviate
      do {
        v1 = 2.0 * sengine.Flat() - 1.0;
        v2 = 2.0 * sengine.Flat() - 1.0;
        v12 = v1 * v1 + v2 * v2;
      } while( v12 > 1.0 );
      t = v1 * Math.sqrt( -2.0 * Math.log( v12 ) / v12 );
      x = s + 0.5 * t;
      gds = x * x;
      if(t >= 0.0 ) return( gds / slambda);         // Immediate acceptance

      u = sengine.Flat();            // Step 3. Uniform random number
      if( d * u <= t * t * t ) return( gds / slambda ); // Squeeze acceptance

      if( sk != aaa ){                               // Step 4. Set-up for hat case
        aaa = sk;
	r = 1.0 / sk;
	q0 = (((((((( q9 * r + q8 ) * r + q7 ) * r + q6 ) * r + q5 ) * r + q4 ) * r + q3 ) * r + q2 ) * r + q1 ) * r;
	if( sk > 3.686 ){
	  if( sk > 13.022 ){
	    b = 1.77;
	    si = 0.75;
	    c = 0.1515 / s;
	  } else {
            b = 1.654 + 0.0076 * ss;
	    si = 1.68 / s + 0.275;
	    c = 0.062 / s + 0.024;
	  }
	} else {
	  b = 0.463 + s - 0.178 * ss;
	  si = 1.235;
	  c = 0.195 / s - 0.079 + 0.016 * s;
	}
      }
      if( x > 0.0 ){                // Step 5. Calculation of q
	v = t / ( s + s );               // Step 6.
	if( Math.abs( v ) > 0.25 ){
	  q = q0 - s * t + 0.25 * t * t + ( ss + ss ) * Math.log( 1.0 + v );
	} else {
	  q = q0 + 0.5 * t * t * (((((((( a9 * v + a8 ) * v + a7 ) * v + a6) * v + a5) * v + a4 ) * v + a3 ) * v + a2 ) * v + a1 ) * v;
	}                // Step 7. Quotient acceptance
	if( Math.log( 1.0 - u ) <= q ) return( gds / slambda );
      }

      for(;;){                    // Step 8. Double exponential deviate t
        do {
	  e = -Math.log( sengine.Flat() );
	  u = sengine.Flat();
	  u = u + u - 1.0;
	  sign_u = ( u > 0 )? 1.0: -1.0;
	  t = b + (e * si) * sign_u;
	} while( t <= -0.71874483771719 );   // Step 9. Rejection of t
	v = t / ( s + s );                  // Step 10. New q(t)
	if( Math.abs( v ) > 0.25 ){
	    q = q0 - s * t + 0.25 * t * t + ( ss + ss ) * Math.log( 1.0 + v );
	} else {
	  q = q0 + 0.5 * t * t * (((((((( a9 * v + a8 ) * v + a7 ) * v + a6 ) * v + a5 ) * v + a4 ) * v + a3 ) * v + a2 ) * v + a1 ) * v;
	}
	if ( q <= 0.0 ) continue;           // Step 11.
	if( q > 0.5 ){
	  w = Math.exp( q ) - 1.0;
	} else {
	  w = (((((( e7 * q + e6 ) * q + e5 ) * q + e4 ) * q + e3 ) * q + e2 ) * q + e1 ) * q;
	}                    // Step 12. Hat acceptance
	if( c * u * sign_u <= w * Math.exp( e - 0.5 * t * t ) ){
	  x = s + 0.5 * t;
	  return( x * x / slambda );
	}
      }
    }
  }; 

  RandGamma.Shoot = function( args ){
    var sk = args.k || 1.0;
    var slambda = args.lambda || 1.0;
    var sengine = args.engine || new JamesRandom({});
    return RandGamma.GenGamma( sengine, sk, slambda );
  };

  RandGamma.ShootArray = function( args ){
    var ssize = args.size || 1;
    var sk = args.k || 1.0;
    var slambda = args.lambda || 1.0;
    var sengine = args.engine || new JamesRandom({});
    var argsShoot = { k: sk, lambda: slambda, engine: sengine };

    for( var i = 0; i < ssize; ++i ){
      args.vect.push( RandGamma.Shoot( argsShoot ) );
    }

  };


  module.exports = RandGamma;
;/* 
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
;/* 
   +----------------------------------------------------------------------+
   |                            HEP Random                                |
   |                         --- RandLandau ---                           |
   |                            Module File                               |
   +----------------------------------------------------------------------+
   Module defining methods for shooting or firing Landau distributed 
   random values.   
   
   The Landau distribution is parameterless and describes the fluctuations
   in energy loss of a particle, making certain assumptions.  For 
   definitions and algorithms, the following papers could be read:
  
   Landau, Jour Phys VIII, No. 4, p. 201 (1944)
   Borsh-Supan, Jour Res. of NBS 65B NO. 4 p. 245 (1961)
   Kolbig & Schorr Comp Phys Comm 31 p. 97 (1984)
  
   The algorithm implemented comes form RANLAN in CERNLIB.
   +----------------------------------------------------------------------+
   | JavaScript                                                           |
   +----------------------------------------------------------------------+
   F. Quinonez - Created 2015-10-30                  
   +----------------------------------------------------------------------+
   | C++                                                                  |
   +----------------------------------------------------------------------+
   M Fischler	  - Created 1/6/2000.
  
  		    The key Transform() method uses the algorithm in CERNLIB.
  		    This is because I trust that RANLAN routine more than 
  		    I trust the Bukin-Grozina inverseLandau, which is not
  		    claimed to be better than 1% accurate.  
  
   M Fischler      - put and get to/from streams 12/13/04
*/

  "use strict";
  var JamesRandom = require('jamesrandom');

  var TABLE_INTERVAL = 0.001;
  var TABLE_END = 982;
  var TABLE_MULTIPLIER = 1.0 / TABLE_INTERVAL;
  var inverseLandau = [
    0.0,							     // .000
    0.0, 	    0.0, 	0.0, 	    0.0, 	0.0, 	     // .001 - .005
    -2.244733, -2.204365, -2.168163, -2.135219, -2.104898,  // .006 - .010
    -2.076740, -2.050397, -2.025605, -2.002150, -1.979866,
    -1.958612, -1.938275, -1.918760, -1.899984, -1.881879,  // .020
    -1.864385, -1.847451, -1.831030, -1.815083, -1.799574,
    -1.784473, -1.769751, -1.755383, -1.741346, -1.727620,  // .030
    -1.714187, -1.701029, -1.688130, -1.675477, -1.663057, 
    -1.650858, -1.638868, -1.627078, -1.615477, -1.604058,  // .040
    -1.592811, -1.581729, -1.570806, -1.560034, -1.549407,
    -1.538919, -1.528565, -1.518339, -1.508237, -1.498254,  // .050
    -1.488386, -1.478628, -1.468976, -1.459428, -1.449979,
    -1.440626, -1.431365, -1.422195, -1.413111, -1.404112,  // .060
    -1.395194, -1.386356, -1.377594, -1.368906, -1.360291,
    -1.351746, -1.343269, -1.334859, -1.326512, -1.318229,  // .070
    -1.310006, -1.301843, -1.293737, -1.285688, -1.277693,
    -1.269752, -1.261863, -1.254024, -1.246235, -1.238494,  // .080
    -1.230800, -1.223153, -1.215550, -1.207990, -1.200474,
    -1.192999, -1.185566, -1.178172, -1.170817, -1.163500,  // .090
    -1.156220, -1.148977, -1.141770, -1.134598, -1.127459,
    -1.120354, -1.113282, -1.106242, -1.099233, -1.092255,  // .100
    
    -1.085306, -1.078388, -1.071498, -1.064636, -1.057802,
    -1.050996, -1.044215, -1.037461, -1.030733, -1.024029,
    -1.017350, -1.010695, -1.004064, -.997456,  -.990871, 
    -.984308, -.977767, -.971247, -.964749, -.958271, 
    -.951813, -.945375, -.938957, -.932558, -.926178, 
    -.919816, -.913472, -.907146, -.900838, -.894547,
    -.888272, -.882014, -.875773, -.869547, -.863337, 
    -.857142, -.850963, -.844798, -.838648, -.832512, 
    -.826390, -.820282, -.814187, -.808106, -.802038, 
    -.795982, -.789940, -.783909, -.777891, -.771884, 	// .150
    -.765889, -.759906, -.753934, -.747973, -.742023, 
    -.736084, -.730155, -.724237, -.718328, -.712429, 
    -.706541, -.700661, -.694791, -.688931, -.683079, 
    -.677236, -.671402, -.665576, -.659759, -.653950, 
    -.648149, -.642356, -.636570, -.630793, -.625022, 
    -.619259, -.613503, -.607754, -.602012, -.596276, 
    -.590548, -.584825, -.579109, -.573399, -.567695, 
    -.561997, -.556305, -.550618, -.544937, -.539262,
    -.533592, -.527926, -.522266, -.516611, -.510961, 
    -.505315, -.499674, -.494037, -.488405, -.482777,	// .200
    
    -.477153, -.471533, -.465917, -.460305, -.454697, 
    -.449092, -.443491, -.437893, -.432299, -.426707, 
    -.421119, -.415534, -.409951, -.404372, -.398795, 
    -.393221, -.387649, -.382080, -.376513, -.370949, 
    -.365387, -.359826, -.354268, -.348712, -.343157, 
    -.337604, -.332053, -.326503, -.320955, -.315408,
    -.309863, -.304318, -.298775, -.293233, -.287692,
    -.282152, -.276613, -.271074, -.265536, -.259999, 
    -.254462, -.248926, -.243389, -.237854, -.232318, 
    -.226783, -.221247, -.215712, -.210176, -.204641, 	// .250
    -.199105, -.193568, -.188032, -.182495, -.176957, 
    -.171419, -.165880, -.160341, -.154800, -.149259,
    -.143717, -.138173, -.132629, -.127083, -.121537, 
    -.115989, -.110439, -.104889, -.099336, -.093782, 
    -.088227, -.082670, -.077111, -.071550, -.065987, 
    -.060423, -.054856, -.049288, -.043717, -.038144, 
    -.032569, -.026991, -.021411, -.015828, -.010243, 
    -.004656,  .000934,  .006527,  .012123,  .017722,
    .023323, .028928,  .034535,  .040146,  .045759, 
    .051376, .056997,  .062620,  .068247,  .073877,	// .300
    
    .079511,  .085149,  .090790,  .096435,  .102083,  
    .107736,  .113392,  .119052,  .124716,  .130385,  
    .136057,  .141734,  .147414,  .153100,  .158789,  
    .164483,  .170181,  .175884,  .181592,  .187304,  
    .193021,  .198743,  .204469,  .210201,  .215937,  
    .221678,  .227425,  .233177,  .238933,  .244696,  
    .250463,  .256236,  .262014,  .267798,  .273587,  
    .279382,  .285183,  .290989,  .296801,  .302619,  
    .308443,  .314273,  .320109,  .325951,  .331799,  
    .337654,  .343515,  .349382,  .355255,  .361135,	// .350  
    .367022,  .372915,  .378815,  .384721,  .390634,  
    .396554,  .402481,  .408415,  .414356,  .420304,
    .426260,  .432222,  .438192,  .444169,  .450153,  
    .456145,  .462144,  .468151,  .474166,  .480188,  
    .486218,  .492256,  .498302,  .504356,  .510418,  
    .516488,  .522566,  .528653,  .534747,  .540850,  
    .546962,  .553082,  .559210,  .565347,  .571493,  
    .577648,  .583811,  .589983,  .596164,  .602355,
    .608554,  .614762,  .620980,  .627207,  .633444,  
    .639689,  .645945,  .652210,  .658484,  .664768,	// .400
    
    .671062,  .677366,  .683680,  .690004,  .696338,  
    .702682,  .709036,  .715400,  .721775,  .728160,  
    .734556,  .740963,  .747379,  .753807,  .760246,  
    .766695,  .773155,  .779627,  .786109,  .792603,  
    .799107,  .805624,  .812151,  .818690,  .825241,  
    .831803,  .838377,  .844962,  .851560,  .858170,
    .864791,  .871425,  .878071,  .884729,  .891399,  
    .898082,  .904778,  .911486,  .918206,  .924940,  
    .931686,  .938446,  .945218,  .952003,  .958802,  
    .965614,  .972439,  .979278,  .986130,  .992996, 	// .450 
    .999875,  1.006769, 1.013676, 1.020597, 1.027533, 
    1.034482, 1.041446, 1.048424, 1.055417, 1.062424,
    1.069446, 1.076482, 1.083534, 1.090600, 1.097681, 
    1.104778, 1.111889, 1.119016, 1.126159, 1.133316, 
    1.140490, 1.147679, 1.154884, 1.162105, 1.169342, 
    1.176595, 1.183864, 1.191149, 1.198451, 1.205770, 
    1.213105, 1.220457, 1.227826, 1.235211, 1.242614, 
    1.250034, 1.257471, 1.264926, 1.272398, 1.279888,
    1.287395, 1.294921, 1.302464, 1.310026, 1.317605, 
    1.325203, 1.332819, 1.340454, 1.348108, 1.355780,	// .500
    
    1.363472, 1.371182, 1.378912, 1.386660, 1.394429, 
    1.402216, 1.410024, 1.417851, 1.425698, 1.433565, 
    1.441453, 1.449360, 1.457288, 1.465237, 1.473206, 
    1.481196, 1.489208, 1.497240, 1.505293, 1.513368, 
    1.521465, 1.529583, 1.537723, 1.545885, 1.554068, 
    1.562275, 1.570503, 1.578754, 1.587028, 1.595325,
    1.603644, 1.611987, 1.620353, 1.628743, 1.637156, 
    1.645593, 1.654053, 1.662538, 1.671047, 1.679581, 
    1.688139, 1.696721, 1.705329, 1.713961, 1.722619, 
    1.731303, 1.740011, 1.748746, 1.757506, 1.766293, 	// .550
    1.775106, 1.783945, 1.792810, 1.801703, 1.810623, 
    1.819569, 1.828543, 1.837545, 1.846574, 1.855631,
    1.864717, 1.873830, 1.882972, 1.892143, 1.901343, 
    1.910572, 1.919830, 1.929117, 1.938434, 1.947781, 
    1.957158, 1.966566, 1.976004, 1.985473, 1.994972, 
    2.004503, 2.014065, 2.023659, 2.033285, 2.042943, 
    2.052633, 2.062355, 2.072110, 2.081899, 2.091720, 
    2.101575, 2.111464, 2.121386, 2.131343, 2.141334,
    2.151360, 2.161421, 2.171517, 2.181648, 2.191815, 
    2.202018, 2.212257, 2.222533, 2.232845, 2.243195,	// .600
    
    2.253582, 2.264006, 2.274468, 2.284968, 2.295507, 
    2.306084, 2.316701, 2.327356, 2.338051, 2.348786, 
    2.359562, 2.370377, 2.381234, 2.392131, 2.403070, 
    2.414051, 2.425073, 2.436138, 2.447246, 2.458397, 
    2.469591, 2.480828, 2.492110, 2.503436, 2.514807, 
    2.526222, 2.537684, 2.549190, 2.560743, 2.572343,
    2.583989, 2.595682, 2.607423, 2.619212, 2.631050, 
    2.642936, 2.654871, 2.666855, 2.678890, 2.690975, 
    2.703110, 2.715297, 2.727535, 2.739825, 2.752168, 
    2.764563, 2.777012, 2.789514, 2.802070, 2.814681,	// .650 
    2.827347, 2.840069, 2.852846, 2.865680, 2.878570, 
    2.891518, 2.904524, 2.917588, 2.930712, 2.943894,
    2.957136, 2.970439, 2.983802, 2.997227, 3.010714,
    3.024263, 3.037875, 3.051551, 3.065290, 3.079095, 
    3.092965, 3.106900, 3.120902, 3.134971, 3.149107, 
    3.163312, 3.177585, 3.191928, 3.206340, 3.220824, 
    3.235378, 3.250005, 3.264704, 3.279477, 3.294323, 
    3.309244, 3.324240, 3.339312, 3.354461, 3.369687,
    3.384992, 3.400375, 3.415838, 3.431381, 3.447005, 
    3.462711, 3.478500, 3.494372, 3.510328, 3.526370,	// .700
    
    3.542497, 3.558711, 3.575012, 3.591402, 3.607881, 
    3.624450, 3.641111, 3.657863, 3.674708, 3.691646, 
    3.708680, 3.725809, 3.743034, 3.760357, 3.777779, 
    3.795300, 3.812921, 3.830645, 3.848470, 3.866400, 
    3.884434, 3.902574, 3.920821, 3.939176, 3.957640, 
    3.976215, 3.994901, 4.013699, 4.032612, 4.051639,
    4.070783, 4.090045, 4.109425, 4.128925, 4.148547,  
    4.168292, 4.188160, 4.208154, 4.228275, 4.248524, 
    4.268903, 4.289413, 4.310056, 4.330832, 4.351745, 
    4.372794, 4.393982, 4.415310, 4.436781, 4.458395, 
    4.480154, 4.502060, 4.524114, 4.546319, 4.568676,	// .750 
    4.591187, 4.613854, 4.636678, 4.659662, 4.682807,
    4.706116, 4.729590, 4.753231, 4.777041, 4.801024, 
    4.825179, 4.849511, 4.874020, 4.898710, 4.923582, 
    4.948639, 4.973883, 4.999316, 5.024942, 5.050761, 
    5.076778, 5.102993, 5.129411, 5.156034, 5.182864, 
    5.209903, 5.237156, 5.264625, 5.292312, 5.320220, 
    5.348354, 5.376714, 5.405306, 5.434131, 5.463193,
    5.492496, 5.522042, 5.551836, 5.581880, 5.612178, 
    5.642734, 5.673552, 5.704634, 5.735986, 5.767610,	// .800
    
    5.799512, 5.831694, 5.864161, 5.896918, 5.929968, 
    5.963316, 5.996967, 6.030925, 6.065194, 6.099780, 
    6.134687, 6.169921, 6.205486, 6.241387, 6.277630, 
    6.314220, 6.351163, 6.388465, 6.426130, 6.464166, 
    6.502578, 6.541371, 6.580553, 6.620130, 6.660109, 
    6.700495, 6.741297, 6.782520, 6.824173, 6.866262,
    6.908795, 6.951780, 6.995225, 7.039137, 7.083525, 
    7.128398, 7.173764, 7.219632, 7.266011, 7.312910, 
    7.360339, 7.408308, 7.456827, 7.505905, 7.555554, 
    7.605785, 7.656608, 7.708035, 7.760077, 7.812747, 	// .850
    7.866057, 7.920019, 7.974647, 8.029953, 8.085952, 
    8.142657, 8.200083, 8.258245, 8.317158, 8.376837,
    8.437300, 8.498562, 8.560641, 8.623554, 8.687319, 
    8.751955, 8.817481, 8.883916, 8.951282, 9.019600, 
    9.088889, 9.159174, 9.230477, 9.302822, 9.376233, 
    9.450735, 9.526355, 9.603118, 9.681054, 9.760191, 
     9.840558,  9.922186, 10.005107, 10.089353, 10.174959,
    10.261958, 10.350389, 10.440287, 10.531693, 10.624646,
    10.719188, 10.815362, 10.913214, 11.012789, 11.114137,
    11.217307, 11.322352, 11.429325, 11.538283, 11.649285,	// .900
    
    11.762390, 11.877664, 11.995170, 12.114979, 12.237161, 
    12.361791, 12.488946, 12.618708, 12.751161, 12.886394, 
    13.024498, 13.165570, 13.309711, 13.457026, 13.607625, 
    13.761625, 13.919145, 14.080314, 14.245263, 14.414134, 
    14.587072, 14.764233, 14.945778, 15.131877, 15.322712, 
    15.518470, 15.719353, 15.925570, 16.137345, 16.354912, 
    16.578520, 16.808433, 17.044929, 17.288305, 17.538873, 
    17.796967, 18.062943, 18.337176, 18.620068, 18.912049, 
    19.213574, 19.525133, 19.847249, 20.180480, 20.525429, 
    20.882738, 21.253102, 21.637266, 22.036036, 22.450278, 	// .950
    22.880933, 23.329017, 23.795634, 24.281981, 24.789364, 
    25.319207, 25.873062, 26.452634, 27.059789, 27.696581,     // .960
    28.365274, 29.068370, 29.808638, 30.589157, 31.413354, 
    32.285060, 33.208568, 34.188705, 35.230920, 36.341388,     // .970
    37.527131, 38.796172, 40.157721, 41.622399, 43.202525, 
    44.912465, 46.769077, 48.792279, 51.005773, 53.437996,     // .980
    56.123356, 59.103894, 					// .982
  ];
    


  function RandLandau( args ){
    var sengine = args.engine || new JamesRandom({});

    this.Fire = function(){
      return RandLandau.Transform( this.fengine.Flat() );
    };

    this.FireArray = function( /* size of vect */ size, /* Array */ vect ){
      for( var i = 0; i < size; ++i ){
        vect.push( this.Fire() );
      }
    };

  } 

    
  RandLandau.Shoot = function( args ){
    var sengine = args.engine || new JamesRandom({});
    return RandLandau.Transform( sengine.Flat() );

  };

  RandLandau.ShootArray = function( args ){
    var ssize = args.size || 1;
    var sengine = args.engine || new JamesRandom({});
    // var svect = args.vect;

    var argsShoot = { engine: sengine };

    for( var i = 0; i < ssize; ++i ){
      args.vect.push( RandLandau.Shoot( argsShoot ) );
    }
  };

  RandLandau.Transform = function( r ){
    var  u = r * TABLE_MULTIPLIER; 
    var index = ( u | 0 );
    var du = u - index;
  
    // du is scaled such that the we dont have to multiply by TABLE_INTERVAL
    // when interpolating.
  
    // Five cases:
    // A) Between .070 and .800 the function is so smooth, straight
    //	linear interpolation is adequate.
    // B) Between .007 and .070, and between .800 and .980, quadratic
    //    interpolation is used.  This requires the same 4 points as
    //	a cubic spline (thus we need .006 and .981 and .982) but
    //	the quadratic interpolation is accurate enough and quicker.
    // C) Below .007 an asymptotic expansion for low negative lambda 
    //    (involving two logs) is used; there is a pade-style correction 
    //	factor.
    // D) Above .980, a simple pade approximation is made (asymptotic to
    //    1/(1-r)), but...
    // E) the coefficients in that pade are different above r=.999.
  
    if( index >= 70 && index <= 800 ){		// (A)
  
      var f0 = inverseLandau[ index ];
      var f1 = inverseLandau[ index+1 ];
      return f0 + du * ( f1 - f0 );
  
    } else if (  index >= 7 && index <= 980  ) {	// ( B )
  
      var f_1 = inverseLandau[ index-1 ];
      var f0  = inverseLandau[ index ];
      var f1  = inverseLandau[ index+1 ];
      var f2  = inverseLandau[ index+2 ];
  
      return f0 + du * ( f1 - f0 - 0.25 * ( 1 - du ) * ( f2 -f1 - f0 + f_1 )  );
  
    } else if (  index < 7  ) {			// ( C )
  
      var n0 =  0.99858950;
      var n1 = 34.5213058;	var d1 = 34.1760202;
      var n2 = 17.0854528;	var d2 =  4.01244582;
  
      var logr = Math.log( r );
      var x    = 1 / logr;
      var x2   = x * x;
  
      var pade = ( n0 + n1 * x + n2 * x2 ) / ( 1.0 + d1 * x + d2 * x2 );
  
      return ( - Math.log ( -0.91893853 - logr ) -1  ) * pade;
  
    } else if (  index <= 999  ) {			// ( D )
  
      var n0 =  1.00060006;
      var n1 =  263.991156;	var d1 =  257.368075;
      var n2 = 4373.20068;	var d2 = 3414.48018;
  
      var x = 1 - r;
      var x2 = x * x;
  
      return ( n0 + n1 * x + n2 * x2 ) / ( x * ( 1.0 + d1 * x + d2 * x2 ) );
  
    } else { 					// ( E )
  
      var n0 =   1.00001538;
      var n1 =   6075.14119;	var d1 = 6065.11919;
      var n2 = 734266.409;	var d2 = 694021.044;
  
      var x = 1 - r;
      var x2   = x * x;
  
      return ( n0 + n1 * x + n2 * x2 ) / ( x * ( 1.0 + d1 * x + d2 * x2 ) );
  
    }        
  };

  module.exports = RandLandau;

;/* 
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

;/* 
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
;/* 
   +----------------------------------------------------------------------+
   |                          BBT Visualization                           |
   |                            --- Axis ---                              |
   |                            Module File                               |
   +----------------------------------------------------------------------+
  
   Graphical Module Axis.
   This module has all functionalities for graphical display of axes for 
   user data Histograms.
  
   +----------------------------------------------------------------------+
   | JavaScript                                                           |
   +----------------------------------------------------------------------+
   F. Quinonez - Created 2015-10-30 
   +----------------------------------------------------------------------+
   | C++                                                                  |
   +----------------------------------------------------------------------+
   Inspired by Class TAxis https://root.cern.ch/doc/master/classTAxis.html 
   of ROOT https://root.cern.ch/ 

   */

  "use strict";
  var d3 = require('../../bower_components/d3/d3.js')


  // Object Constructor BBH1 Building Block Histogram 1D. 
 
  function Axis( name, title ){
    this.name = name;
    this.title = title;

    this.fNdivisions;
    this.fAxisColor;
    this.fLabelColor;
    this.fLabelFont;
    this.fLabelOffset;
    this.fLabelOffset;
    this.fLabelSize;
    this.fTitleOffset;
    this.d3Axis = d3.svg.axis();


    // ****************************************************
    // * Functions coming from ROOT mother class TAxisAtt *
    // ****************************************************

    // 
    /*SaveAttributes: function( out, name, subname ){
      
      return this;
    },

    // Set color of the line axis and tick marks.
    SetAxisColor: function( color ){
      this.d3Axis.style( "stroke", color );     

      return this;
    },

    // Set color of labels.
    SetLabelColor: function( color ){

      return this;
    },

    // Set labels' font.
    SetLabelFont: function( font ){

      return this;
    },

    // Set distance between the axis and the labels.
    // In ROOT: The distance is expressed in per cent of the pad width.
    SetLabelOffset: function( offset  ){

      return this;
    },

    // Set size of axis labels.
    // In ROOT: The size is expressed in per cent of the pad width.
    SetLabelSize: function( size  ){

      return this;
    },

    /* Set the number of divisions for this axis.

    if optim = kTRUE (default), the number of divisions will be
                   optimized around the specified value.
    if optim = kFALSE, or n < 0, the axis will be forced to use
                   exactly n divisions.

    n = n1 + 100*n2 + 10000*n3

    Where n1 is the number of primary divisions,
    n2 is the number of second order divisions and
    n3 is the number of third order divisions.

    e.g. 512 means 12 primary and 5 secondary divisions.

    If the number of divisions is "optimized" (see above) n1, n2, n3 are
    maximum values. 
    SetNdivisions: function( n, optim  ){
      
      return this;
    },
    SetNdivisions: function( n1, n2, n3, optim){
      
      return this;
    },

    // Set tick mark length
    // The length is expressed in per cent of the pad width
    SetTickLength: function( lenght ){
      
      return this;
    },

    /* Set distance between the axis and the axis title
    Offset is a correction factor with respect to the "standard" value.
    offset = 1  uses the default position that is computed in function
    of the label offset and size.
    offset = 1.2 will add 20 per cent more to the default offset. 
    SetTitleOffset: function( offset  ){
    
      return this;
    },
    
    // Set size of axis title.
    // The size is expressed in per cent of the pad width.
    SetTitleSize: function( size  ){
    
      return this;
    },

    // Set color of axis title
    SetTitleColor: function( color ){

      return this;
    },

    // Set the title font.
    SetTitleFont: function( font ){
      
      return this;
    },

    GetNdivisions: function(){
      return this.fNdivisions;
    },

    GetAxisColor: function(){
      return this.fAxisColor;
    },

    GetLabelColor: function(){
      return this.fLabelColor;
    },

    GetLabelFont: function(){
      return this.fLabelFont;
    },

    GetLabelOffset: function(){
      return this.fLabelOffset;
    },

    GetLabelSize: function(){
      return this.fLabelSize;
    },

    GetTitleOffset: function(){
      return this.fTitleOffset;
    },

    GetTitleSize: function(){
      return this.fTitleSize;
    },

    GetTitleColor: function(){
      return this.fTitleColor;
    },

    GetTitleFont: function(){
      return this.fTitleFont;
    },
*/
  }

  module.exports = Axis;
;/* 
   +----------------------------------------------------------------------+
   |                          BBT Visualization                           |
   |                             --- H1 ---                               |
   |                            Module File                               |
   +----------------------------------------------------------------------+
  
   Graphical Module H1.
   This module has all functionalities for graphical display of user data 
   by using 1D Histograms.
  
   +----------------------------------------------------------------------+
   | JavaScript                                                           |
   +----------------------------------------------------------------------+
   F. Quinonez - Created 2015-10-30 
               - Developed functions: H1, Fill, GetXaxis, GetYaxis, 
                 FillRandom

   A. J. Hernandez Goez - Developed functions: Draw.
   +----------------------------------------------------------------------+
   | C++                                                                  |
   +----------------------------------------------------------------------+
   Inspired by Class TH1 https://root.cern.ch/doc/master/classTH1.html 
   of ROOT https://root.cern.ch/ 

   */

  'use strict';

  var JamesRandom     = require('jamesrandom');

  var d3              = require('../../bower_components/d3/d3.js');

  var Axis            = require('./Axis.js');

  var RandBinomial    = require('../Generation/Random/RandBinomial.js'); 
  var RandBit         = require('../Generation/Random/RandBit.js'); 
  var RandBreitWigner = require('../Generation/Random/RandBreitWigner.js'); 
  var RandChiSquare   = require('../Generation/Random/RandChiSquare.js'); 
  var RandExponential = require('../Generation/Random/RandExponential.js'); 
  var RandFlat        = require('../Generation/Random/RandFlat.js');
  var RandGamma       = require('../Generation/Random/RandGamma.js'); 
  var RandGauss       = require('../Generation/Random/RandGauss.js'); 
  var RandLandau      = require('../Generation/Random/RandLandau.js'); 
  var RandPoisson     = require('../Generation/Random/RandPoisson.js'); 
  var RandStudentT    = require('../Generation/Random/RandStudentT.js');
 

  // BBT Histogram 1D Constructor. 
  function H1( name, title, nbinsx, xmin, xmax ){
    this.fDimension = 1; // Dimension of the plot data.
    this.name = name;
    this.title = title;
    this.nbinsx = nbinsx;
    this.xmin = xmin;
    this.xmax = xmax;
    this.rawData = []; // data array filled by Fill function.
    this.freqData = new Array( this.nbinsx + 2 ); // data array of frequencies for each bin. Two bins extra have been added, one for underflow and the another one for overflow.

    this.fXaxis = new Axis("xaxishisto", "Eje x en [u]");
    this.fYaxis = new Axis("xaxishisto", "Eje x en [u]");
	//this.fPainter = Object.create( Painter() );

    this.Fill = function( value ) {
      return this.rawData.push( value );
    };

    this.FillRandom = function( args ){
      var siz3 = args.size || 1;
      var n = args.n || siz3;
      siz3 = n;
      var v3ct = args.vect || [];
      var pdf = args.pdf || "Gauss";
      var s = args.seed || new Date().getDate();

      var pdf_mean = args.pdf_mean || undefined; // Used in: BreitWigner, Gauss, Poisson, Exponential
      console.log( "mean: %f", pdf_mean );
      var pdf_stdDev = args.pdf_stdDev || undefined; // Used in: Gauss
      console.log( "stdDev: %f", pdf_stdDev );
      var pdf_a = args.pdf_a || undefined; // Used in: ChiSquare, StudentT, Flat
      var pdf_b = args.pdf_b || undefined; // Used in: Flat
      var pdf_width = pdf_b - pdf_a; // Used in: Flat
      var pdf_cut = args.pdf_cut || undefined; // Used in: BreitWigner
      var pdf_n = args.pdf_n || undefined; // Used in: Binomial
      var pdf_p = args.pdf_p || undefined; // Used in: Binomial
      var pdf_k = args.pdf_k || undefined; // Used in: Gamma
      var pdf_lambda = args.pdf_lambda || undefined; // Used in: Gamma
      var pdf_gamma = args.pdf_gamma || undefined; // Used in: BreitWigner

      switch( pdf ){
        case "Bit":
          var engin3 = new JamesRandom( { seed: s } );
          var argum = { size: siz3, vect: v3ct, engine: engin3 };
          RandBit.ShootArray( argum );
          this.rawData = v3ct;
          break;

        case "Binomial":
          var engin3 = new JamesRandom( { seed: s } );
          var argum = { size: siz3, vect: v3ct, engine: engin3, n: pdf_n, p: pdf_p };
          RandBinomial.ShootArray( argum );
          this.rawData = v3ct;
          break;

        case "BreitWigner":
          engin3 = new JamesRandom( { seed: s } );
          var argum = { size: siz3, vect: v3ct, engine: engin3, mean: pdf_mean, gamma: pdf_gamma, cut: pdf_cut };
          RandBreitWigner.ShootArray( argum );
          this.rawData = v3ct;
          break;

        case "ChiSquare":
          var engin3 = new JamesRandom( { seed: s } );
          var argum = { size: siz3, vect: v3ct, engine: engin3, a: pdf_a };
          RandChiSquare.ShootArray( argum );
          this.rawData = v3ct;
          break;

        case "Exponential":
          var engin3 = new JamesRandom( { seed: s } );
          var argum = { size: siz3, vect: v3ct, engine: engin3, mean: pdf_mean };
          RandExponential.ShootArray( argum );
          this.rawData = v3ct;
          break;

        case "Flat":
          var engin3 = new JamesRandom( { seed: s } );
          var argum = { size: siz3, vect: v3ct, engine: engin3, a: pdf_a, b: pdf_b, width: pdf_width };
          RandFlat.ShootArray( argum );
          this.rawData = v3ct;
          break;

        case "Gamma":
          var engin3 = new JamesRandom( { seed: s } );
          var argum = { size: siz3, vect: v3ct, engine: engin3, k: pdf_k, lamda: pdf_lambda };
          RandGamma.ShootArray( argum );
          this.rawData = v3ct;
          break;

        case "Gauss":
          var engin3 = new JamesRandom( { seed: s } );
          var argum = { size: siz3, vect: v3ct, engine: engin3, mean: pdf_mean, stdDev: pdf_stdDev };
          RandGauss.ShootArray( argum );
          this.rawData = v3ct;
          break;

        case "Landau":
          var engin3 = new JamesRandom( { seed: s } );
          var argum = { size: siz3, vect: v3ct, engine: engin3 };
          RandLandau.ShootArray( argum );
          this.rawData = v3ct;
          break;

        case "Poisson":
          var engin3 = new JamesRandom( { seed: s } );
          var argum = { size: siz3, vect: v3ct, engine: engin3, mean: pdf_mean };
          RandPoisson.ShootArray( argum );
          this.rawData = v3ct;
          break;

        case "StudentT":
          var engin3 = new JamesRandom( { seed: s } );
          var argum = { size: siz3, vect: v3ct, engine: engin3, a: pdf_a };
          RandStudentT.ShootArray( argum );
          this.rawData = v3ct;
          break;
      }
    };
 
    // Deattach data preparation
    this.Prepare = function(){
      var binxlow = this.xmin;
      var binwidth = ( this.xmax - this.xmin ) / this.nbinsx;
      var binxup = binxlow + binwidth;
  
      // Initialize frequencies to zero.
      for( var i = 0; i < this.freqData.length; i++ ){
        this.freqData[ i ] = 0;
      }
  
      // First fill array of frequencies for each bin: freqData.
      for( var i = 0; i < this.rawData.length; i++ ){
        binxlow = this.xmin;
        binxup = binxlow + binwidth;
        if( this.rawData[ i ] < this.xmin ) this.freqData[ 0 ] += 1; // Underflow.
        if( this.rawData[ i ] >= this.xmax ) this.freqData[ this.nbinsx + 1 ] += 1; // Overflow.
        // From index 1 til index nbinsx-2 of freqData Array.
        for( var j = 1; j < this.freqData.length - 1; j++ ){
          if( this.rawData[ i ] >= binxlow && this.rawData[ i ] < binxup ){
            this.freqData[ j ] += 1;
            break;
          }
          // Go to the next bin
          binxlow = binxup;
          binxup += binwidth;   
        }
      }
      console.log(this.freqData);
       
    };
    
    this.Draw = function(){
      this.Prepare();
      var binwidth = ( this.xmax - this.xmin ) / this.nbinsx;
      // **************************************************************
      // Now start visualization
      // **************************************************************
      var margin = { top:20, right:40, bottom: 30, left:60 };
      var width = 640 -  margin.left - margin.right ;
      var height = 400 -  margin.top - margin.bottom ;
      var svg = d3.select( "#c1" ).append( "svg" )
        .attr( "width", width + margin.left + margin.right )
        .attr( "height", height + margin.top + margin.bottom );
      // Drawing a chart inside the svg
      var chart = svg.append( "g" )
        .attr( "transform", "translate(" + margin.left + "," + margin.top + ")" );
      // Scaling the data
      var xScale = d3.scale.linear().domain([ this.xmin, this.xmax ]).range( [ 0, width ] );
      var yScale = d3.scale.linear().domain([ 0, d3.max( this.freqData ) ]).range( [ height, 0 ] );
      // Defining the plot's domain
      //xScale..nice();
      //yScale.nice();

      var xAxis = d3.svg.axis()
        .scale( xScale )
        .orient( "bottom" );

      var yAxis = d3.svg.axis()
        .scale( yScale )
        .orient( "left" );

      var xAxisGroup = chart.append( "g" )
        .attr( "transform", "translate(0," + height + ")" );

      var yAxisGroup = chart.append( "g" )
        .attr( "transform", "translate(0,0)" );

      xAxis( xAxisGroup );
      yAxis( yAxisGroup );
      // Definining text labels for axis.
      // Don't do it because that need to be accessed by the user

      // Displaying data    
      this.freqData.forEach(function(d,i,a){
		        if( i > 0 && i < a.length - 1 ) chart.append( "rect" )
          .attr( "x", function(){return (i-1)*xScale(binwidth);} )
          .attr( "width", xScale(binwidth)-2)
          .attr( "y", yScale(d) )
          .attr( "height", function(){return height-yScale(d);}  )
		        .attr("fill", "steelblue");
      });

    }; // Ends function Draw

    this.GetXaxis = function(){
      return this.fXaxis;
    };

    this.GetYaxis = function(){
      return this.fYaxis;
    };

  }

  module.exports = H1;
