'use strict';

/* Services */

var BBTServices = angular.module('BBTServices', ['ngResource']);

BBTServices.factory('UserService', [function() {
  var sdo = {
    isLogged: false,
    user: '',
    profile:'',
    id:''
  };
  return sdo;
}]);
BBTServices.factory('Simulacro', [function() {
  var sdo = {
    materia:''
  };
  return sdo;
}]);

BBTServices.factory('QuestionService', [function() {
  var sdo = [];
  return sdo;
}]);
BBTServices.factory('LC11Srvc', ['$resource',
  function($resource){
    return $resource('phones/:phoneId.json', {}, {
      phpfilltable_lc11q: {method:'PUT', params:{phoneId:'phones'}, isArray:true}
    });
  }]);


BBTServices.factory('GlobalesSesion', ['$resource',
  function($resource){
    return $resource('users/:userId.json', {userId:'users'}, {get: { method: "GET", isArray: true }});
  }]);

