/* 
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

