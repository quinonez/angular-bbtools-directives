/* 
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
