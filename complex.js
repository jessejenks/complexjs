/**
 * A library of complex variables and functions, emphasizing Riemann Zeta
 * @author Jesse E Jenks
 * @author Sam Burdick
 */

// constants
var EPSILON = 0.0000001;
var PI = new Complex(Math.PI,0);

// a complex variable

function Complex(real, imag) {
  this.real = real;
  this.imag = imag;
  this.toString = function() { // TODO: decide how to print if re or im == 0
    var op;
    var im = this.imag;
    if (this.imag >= 0){
      op = '+';
    } else {
      op = '-';
    }
    if(this.imag === 1 || this.imag === -1) {
      im = '';
    }
    return this.real + ' ' + op + ' ' + im + 'i';
  }
};

// Elementary operations

/*
* Examples:
* var w = new Complex(2,5); // 2+5i
* var z = new Complex(1,1); // 1+2i
* var w_x_z = w.multiply(z); // (2*1 - 5*2) + (2*2 + 5*1)i
*/

Complex.prototype.multiply = function(z) {
  return new Complex(this.real*z.real - this.imag*z.imag, this.real*z.imag + this.imag*z.real);
};
Complex.prototype.scalarMult = function(r) {
  return new Complex(this.real*r, this.imag*r);
};
Complex.prototype.power = function(n) {
  var zNew = new Complex(this.real, this.imag);
  for (var i = 0; i<n-1; i++) {
    zNew = zNew.multiply(this);
  }
  return zNew;
};
Complex.prototype.add = function(z) {
  return new Complex(this.real+z.real, this.imag+z.imag);
};
Complex.prototype.add_real = function(r){
 return new Complex(this.real+ r, this.imag);
};
Complex.prototype.subtract = function(z) {
  return new Complex(this.real-z.real, this.imag-z.imag);
};
/*
* w/z = (a+bi)/(c+di) = (a+bi)(c+di)/(c+di)^2 = (ac - bd)/(c^2 - d^2) + (ad + bc)/(c^2 - d^2) i
*/
Complex.prototype.divide = function(z) {
  return new Complex( (this.real*z.real + this.imag*z.imag)/ (z.real*z.real + z.imag*z.imag), ( z.real*this.imag - this.real*z.imag )/(z.real*z.real + z.imag*z.imag));
};
Complex.prototype.arg = function () {
  return Math.atan2(this.imag, this.real);
};
Complex.prototype.raise_to = function (c) {
  var theta = this.imag*Math.log(c);
  var coeff = Math.pow(c,this.real);
  return new Complex(coeff*Math.cos(theta),coeff*Math.sin(theta));
};
Complex.prototype.e_to_the = function() {
  return new Complex(Math.exp(this.real)*Math.cos(this.imag),Math.exp(this.real)*Math.sin(this.imag));
};
Complex.prototype.distance_sqr = function(z){
  return (this.real-z.real)*(this.real-z.real)+(this.imag-z.imag)*(this.imag-z.imag);
};
Complex.prototype.magnitude = function() {
  return sqrt(this.real*this.real + this.imag*this.imag);
};
Complex.prototype.mag_sqr = function(){
  return this.real*this.real + this.imag*this.imag;
};
Complex.prototype.sub_real = function(x) {
  return new Complex(this.real - x, this.imag);
};
Complex.prototype.add_real = function(x) {
  return new Complex(this.real + x, this.imag);
};
Complex.prototype.raise_to_z = function(z) {
  var a = this.real;
  var b = this.imag;
  var c = z.real;
  var d = z.imag;
  var arg = this.arg();
  var asbs = (a*a + b*b);
  var param = c*arg + (1/2) * d * Math.log(asbs);
  var to_return = new Complex ( Math.cos(param), Math.sin(param));
  var sc = Math.pow(asbs , c / 2) * Math.pow(Math.E , -d*arg);
  return to_return.scalarMult(sc);
};

// Taylor evaluations of trig fns

Complex.prototype.sin_taylor = function() {
// Complex.prototype.sin_taylor = function(k) {
  // if (k === 1) return this;
  // TODO: return a degree-k aproximation
  // I don't know if that should mean at most the k-th power
  // or 2*k + 1-st power
  return this.subtract(this.power(3).scalarMult(1/6)).add(this.power(5).scalarMult(1/120)).subtract(this.power(7).scalarMult(1/5040)).add(this.power(9).scalarMult(1/362880));
};
Complex.prototype.cos_taylor = function() {
  // z.power(2).multiply().add_real(1);
  return this.power(2).scalarMult(-1/2).add(this.power(4).scalarMult(1/24)).subtract(this.power(6).scalarMult(1/720)).add(this.power(8).scalarMult(1/40320)).add_real(1);
};

// Gamma function

// polynomial coeffs
var P = [ 676.5203681218851,   -1259.1392167224028,
          771.32342877765313, -176.61502916214059,
          12.507343278686905, -0.13857109526572012,
          9.9843695780195716e-6, 1.5056327351493116e-7];

/**
 * Lanczos approximation of Γ(z) implementation based on https://en.wikipedia.org/wiki/Lanczos_approximation
 * @param Complex z
 */
function gamma(z) {
  var result;
  var h = 0.5;
  if (z.real < h) {
    result = PI.divide( z.scalarMult(Math.PI)
                         .sin_taylor()
                         .multiply( gamma( z.scalarMult(-1).add_real(1) ) ) );
  } else {
    z.real -= 1;
    var x = new Complex(0.99999999999980993 , 0);
    var d = new Complex(1,1);
    for(var i = 0; i < P.length; i++) {
      x = x.add( new Complex(P[i] , 0)
            .divide( z.add_real(i+1) ) );
    }
    var t = z.add_real(P.length)
             .sub_real(h);
    result = x.scalarMult( Math.sqrt(2 * Math.PI) )
              .multiply( t.raise_to_z( z.add_real(h) ) )
              .multiply( t.scalarMult(-1).e_to_the() );
  }
  if (within_epsilon(result.imag)) {
    return result.real;
  }
  return result;
}

/**
 * @param real x
 */
function within_epsilon(x) {
  return Math.abs(x) <= EPSILON;
}

// Riemann Zeta function

Complex.prototype.zeta_n = function () {
  // d/ds zeta_N = - sum_{n=1}^{N} log(n)/n^s
};
/*
* based on work from http://people.math.sfu.ca/~pmenz/thesis.pdf
* http://numbers.computation.free.fr/Constants/Miscellaneous/zetaevaluations.pdf
*/
Complex.prototype.zeta = function(epsilon) {
  var p, q;
  if (this.real === 0 && this.imag === 0) return new Complex(-1/2, 0);
  if (this.real === 1 && this.imag === 0) return new Complex(Infinity, 0);
  if (this.real < 0) {
    var news = this.scalarMult(-1).add_real(1);
    return news.scalarMult(-1).raise_to(2*Math.PI)
    .multiply(this.scalarMult(Math.PI/2).sin_taylor())
    .multiply(news.gamma())
    .multiply(news.zeta(epsilon))
    .scalarMult(2);
  }
  // else {
  //   var alpha;
  //   if (this.imag === 0) {
  //
  //   } else {
  //
  //   }
  // }
};
