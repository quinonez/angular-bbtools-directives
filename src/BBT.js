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
