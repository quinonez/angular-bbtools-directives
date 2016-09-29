'use strict';

/* App Module */

var BBTApp = angular.module('BBTApp', [
  'ngRoute',
  'BBTAnimations',
  'BBTControllers',
  'BBTFilters',
  'BBTServices',
  'BBTDirectives'
]);

BBTApp.config(['$routeProvider',
  function($routeProvider) {
    $routeProvider
     .when('/', {
        templateUrl: 'index.html',
        controller: 'StartCtrl'
      })
     .when('/random', {
        templateUrl: 'demo/random.html',
        controller: 'RandomCtrl'
      })
     .when('/pdf', {
        templateUrl: 'demo/pdf.html',
        controller: 'PDFCtrl'
      })
     .when('/h1', {
        templateUrl: 'demo/h1.html',
        controller: 'H1Ctrl'
      })

     ;
  }]);
