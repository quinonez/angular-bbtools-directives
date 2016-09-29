# JamesRandom Module
This pseudo random number generator have been written in JavaScript in order to bring 
a better mathematical behavior that it is present in Math.Random built in function.

JamesRandom is based on CLHEP - A Class Library for High Energy Physics, which is written in C++.
This JamesRandom written in JS does not have all the functions that have JamesRandom CLHEP, but it does what it must to do, i.e. 
generates random numbers (by default between 0 and 1).
 

One way to use it is:

```javascript
    var JamesRandom = require('jamesrandom');

    var gen = new JamesRandom({seed: 12345});
    var A = [];
    var sizeA = 40;
    gen.FlatArray( sizeA, A );
    console.log(A);
```



## Public Functions

### Name()
Returns String "HepJamesRandom".
### Flat()
Returns one pseudo random number between 0 and 1.
### FlatArray( size, vect ) 
Where ```size``` is an integer and is intended to be the size of the array ```vect``` of pseudorandom numbers. 
### SetSeed( seed, dum )
Where ```seed``` is an integer and ```dum``` is another integer generally 0.
### SetSeeds( seed, dum )
Where ```seed``` is an integer and ```dum``` is another integer generally 0.


## Static Public Functions

### JamesRandom.EngineName()
Returns String "HepJamesRandom".

### JamesRandom.BeginTag()
Returns String "JamesRandom-Begin".

## Public Variables

Variables that keep the present status of the generator

### u
This is an array of size 97.
### c
Double.
### cd
Double.
### cm
Double.
### i97
Integer.
### j97
Integer.

## Static Public Variables

### JamesRandom.SetMarkerLen = 64




