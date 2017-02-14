(function e(t,n,r){function s(o,u){if(!n[o]){if(!t[o]){var a=typeof require=="function"&&require;if(!u&&a)return a(o,!0);if(i)return i(o,!0);var f=new Error("Cannot find module '"+o+"'");throw f.code="MODULE_NOT_FOUND",f}var l=n[o]={exports:{}};t[o][0].call(l.exports,function(e){var n=t[o][1][e];return s(n?n:e)},l,l.exports,e,t,n,r)}return n[o].exports}var i=typeof require=="function"&&require;for(var o=0;o<r.length;o++)s(r[o]);return s})({1:[function(require,module,exports){
// ==ClosureCompiler==
// @output_file_name fit-curve.min.js
// @compilation_level SIMPLE_OPTIMIZATIONS
// ==/ClosureCompiler==

/**
 *  @preserve  JavaScript implementation of
 *  Algorithm for Automatically Fitting Digitized Curves
 *  by Philip J. Schneider
 *  "Graphics Gems", Academic Press, 1990
 *
 *  The MIT License (MIT)
 *
 *  https://github.com/soswow/fit-curves
 */

/**
 * Fit one or more Bezier curves to a set of points.
 *
 * @param {Array<Array<Number>>} points - Array of digitized points, e.g. [[5,5],[5,50],[110,140],[210,160],[320,110]]
 * @param {Number} maxError - Tolerance, squared error between points and fitted curve
 * @returns {Array<Array<Array<Number>>>} Array of Bezier curves, where each element is [first-point, control-point-1, control-point-2, second-point] and points are [x, y]
 */
function fitCurve(points, maxError, progressCallback) {
    if (!Array.isArray(points)) {
        throw new TypeError("First argument should be an array");
    }
    points.forEach((point) => {
        if(!Array.isArray(point) || point.length !== 2
        || typeof point[0] !== 'number' || typeof point[1] !== 'number'){
            throw Error("Each point should be an array of two numbers")
        }
    });
    // Remove duplicate points
    points = points.filter((point, i) =>
        i === 0 || !(point[0] === points[i-1][0] && point[1] === points[i-1][1])
    );

    if (points.length < 2) {
        return [];
    }

    const len = points.length;
    const leftTangent = createTangent(points[1], points[0]);
    const rightTangent = createTangent(points[len - 2], points[len - 1]);

    return fitCubic(points, leftTangent, rightTangent, maxError, progressCallback);
}

/**
 * Fit a Bezier curve to a (sub)set of digitized points.
 * Your code should not call this function directly. Use {@link fitCurve} instead.
 *
 * @param {Array<Array<Number>>} points - Array of digitized points, e.g. [[5,5],[5,50],[110,140],[210,160],[320,110]]
 * @param {Array<Number>} leftTangent - Unit tangent vector at start point
 * @param {Array<Number>} rightTangent - Unit tangent vector at end point
 * @param {Number} error - Tolerance, squared error between points and fitted curve
 * @returns {Array<Array<Array<Number>>>} Array of Bezier curves, where each element is [first-point, control-point-1, control-point-2, second-point] and points are [x, y]
 */
function fitCubic(points, leftTangent, rightTangent, error, progressCallback) {
    const MaxIterations = 20;   //Max times to try iterating (to find an acceptable curve)

    var bezCurve,               //Control points of fitted Bezier curve
        u,                      //Parameter values for point
        uPrime,                 //Improved parameter values
        maxError, prevErr,      //Maximum fitting error
        splitPoint, prevSplit,  //Point to split point set at if we need more than one curve
        centerVector, toCenterTangent, fromCenterTangent,  //Unit tangent vector(s) at splitPoint
        beziers,                //Array of fitted Bezier curves if we need more than one curve
        dist, i;

    //console.log('fitCubic, ', points.length);

    //Use heuristic if region only has two points in it
    if (points.length === 2) {
        dist = maths.vectorLen(maths.subtract(points[0], points[1])) / 3.0;
        bezCurve = [
            points[0],
            maths.addArrays(points[0], maths.mulItems(leftTangent,  dist)),
            maths.addArrays(points[1], maths.mulItems(rightTangent, dist)),
            points[1]
        ];
        return [bezCurve];
    }

    //Parameterize points, and attempt to fit curve
    u = chordLengthParameterize(points);
    [bezCurve, maxError, splitPoint] = generateAndReport(points, u, u, leftTangent, rightTangent, progressCallback)

    if (maxError < error) {
        return [bezCurve];
    }
    //If error not too large, try some reparameterization and iteration
    if (maxError < (error*error)) {

        uPrime = u;
        prevErr = maxError;
        prevSplit = splitPoint;

        for (i = 0; i < MaxIterations; i++) {

            uPrime = reparameterize(bezCurve, points, uPrime);
            [bezCurve, maxError, splitPoint] = generateAndReport(points, u, uPrime, leftTangent, rightTangent, progressCallback);

            if (maxError < error) {
                return [bezCurve];
            }
            //If the development of the fitted curve grinds to a halt,
            //we abort this attempt (and try a shorter curve):
            else if(splitPoint === prevSplit) {
                let errChange = maxError/prevErr;
                if((errChange > .9999) && (errChange < 1.0001)) {
                    break;
                }
            }

            prevErr = maxError;
            prevSplit = splitPoint;
        }
    }

    //Fitting failed -- split at max error point and fit recursively
    beziers = [];

    //To create a smooth transition from one curve segment to the next,
    //we calculate the tangent of the points directly before and after the center,
    //and use that same tangent both to and from the center point.
    centerVector = maths.subtract(points[splitPoint - 1], points[splitPoint + 1]);
    //However, should those two points be equal, the normal tangent calculation will fail.
    //Instead, we calculate the tangent from that "double-point" to the center point, and rotate 90deg.
    if((centerVector[0] === 0) && (centerVector[1] === 0)) {
        //toCenterTangent = createTangent(points[splitPoint - 1], points[splitPoint]);
        //fromCenterTangent = createTangent(points[splitPoint + 1], points[splitPoint]);

        //[x,y] -> [-y,x]: http://stackoverflow.com/a/4780141/1869660
        centerVector = maths.subtract(points[splitPoint - 1], points[splitPoint])
                            .reverse();
        centerVector[0] = -centerVector[0];
    }
    toCenterTangent = maths.normalize(centerVector);
    //To and from need to point in opposite directions:
    fromCenterTangent = maths.mulItems(toCenterTangent, -1);

    /*
    Note: An alternative to this "divide and conquer" recursion could be to always
          let new curve segments start by trying to go all the way to the end,
          instead of only to the end of the current subdivided polyline.
          That might let many segments fit a few points more, reducing the number of total segments.

          However, a few tests have shown that the segment reduction is insignificant
          (240 pts, 100 err: 25 curves vs 27 curves. 140 pts, 100 err: 17 curves on both),
          and the results take twice as many steps and milliseconds to finish,
          without looking any better than what we already have.
    */
    beziers = beziers.concat(fitCubic(points.slice(0, splitPoint + 1), leftTangent, toCenterTangent,    error, progressCallback));
    beziers = beziers.concat(fitCubic(points.slice(splitPoint),        fromCenterTangent, rightTangent, error, progressCallback));
    return beziers;
};

function generateAndReport(points, paramsOrig, paramsPrime, leftTangent, rightTangent, progressCallback) {
    var bezCurve, maxError, splitPoint;

    bezCurve = generateBezier(points, paramsPrime, leftTangent, rightTangent, progressCallback);
    //Find max deviation of points to fitted curve.
    //Here we always use the original parameters (from chordLengthParameterize()),
    //because we need to compare the current curve to the actual source polyline,
    //and not the currently iterated parameters which reparameterize() & generateBezier() use,
    //as those have probably drifted far away and may no longer be in ascending order.
    [maxError, splitPoint] = computeMaxError(points, bezCurve, paramsOrig);

    if(progressCallback) {
        progressCallback({
            bez: bezCurve,
            points: points,
            params: paramsOrig,
            maxErr: maxError,
            maxPoint: splitPoint,
        });
    }

    return [bezCurve, maxError, splitPoint];
}

/**
 * Use least-squares method to find Bezier control points for region.
 *
 * @param {Array<Array<Number>>} points - Array of digitized points
 * @param {Array<Number>} parameters - Parameter values for region
 * @param {Array<Number>} leftTangent - Unit tangent vector at start point
 * @param {Array<Number>} rightTangent - Unit tangent vector at end point
 * @returns {Array<Array<Number>>} Approximated Bezier curve: [first-point, control-point-1, control-point-2, second-point] where points are [x, y]
 */
function generateBezier(points, parameters, leftTangent, rightTangent) {
    var bezCurve,                       //Bezier curve ctl pts
        A, a,                           //Precomputed rhs for eqn
        C, X,                           //Matrices C & X
        det_C0_C1, det_C0_X, det_X_C1,  //Determinants of matrices
        alpha_l, alpha_r,               //Alpha values, left and right

        epsilon, segLength,
        i, len, tmp, u, ux,
        firstPoint = points[0],
        lastPoint = points[points.length-1];

    bezCurve = [firstPoint, null, null, lastPoint];
    //console.log('gb', parameters.length);

    //Compute the A's
    A = maths.zeros_Xx2x2(parameters.length);
    for (i = 0, len = parameters.length; i < len; i++) {
        u = parameters[i];
        ux = 1 - u;
        a = A[i];

        a[0] = maths.mulItems(leftTangent,  3 * u  * (ux*ux));
        a[1] = maths.mulItems(rightTangent, 3 * ux * (u*u));
    }

    //Create the C and X matrices
    C = [[0,0], [0,0]];
    X = [0,0];
    for (i = 0, len = points.length; i < len; i++) {
        u = parameters[i];
        a = A[i];

        C[0][0] += maths.dot(a[0], a[0]);
        C[0][1] += maths.dot(a[0], a[1]);
        C[1][0] += maths.dot(a[0], a[1]);
        C[1][1] += maths.dot(a[1], a[1]);

        tmp = maths.subtract(points[i], bezier.q([firstPoint, firstPoint, lastPoint, lastPoint], u));

        X[0] += maths.dot(a[0], tmp);
        X[1] += maths.dot(a[1], tmp);
    }

    //Compute the determinants of C and X
    det_C0_C1 = (C[0][0] * C[1][1]) - (C[1][0] * C[0][1]);
    det_C0_X  = (C[0][0] * X[1]   ) - (C[1][0] * X[0]   );
    det_X_C1  = (X[0]    * C[1][1]) - (X[1]    * C[0][1]);

    //Finally, derive alpha values
    alpha_l = det_C0_C1 === 0 ? 0 : det_X_C1 / det_C0_C1;
    alpha_r = det_C0_C1 === 0 ? 0 : det_C0_X / det_C0_C1;

    //If alpha negative, use the Wu/Barsky heuristic (see text).
    //If alpha is 0, you get coincident control points that lead to
    //divide by zero in any subsequent NewtonRaphsonRootFind() call.
    segLength = maths.vectorLen(maths.subtract(firstPoint, lastPoint));
    epsilon = 1.0e-6 * segLength;
    if (alpha_l < epsilon || alpha_r < epsilon) {
        //Fall back on standard (probably inaccurate) formula, and subdivide further if needed.
        bezCurve[1] = maths.addArrays(firstPoint, maths.mulItems(leftTangent,  segLength / 3.0));
        bezCurve[2] = maths.addArrays(lastPoint,  maths.mulItems(rightTangent, segLength / 3.0));
    } else {
        //First and last control points of the Bezier curve are
        //positioned exactly at the first and last data points
        //Control points 1 and 2 are positioned an alpha distance out
        //on the tangent vectors, left and right, respectively
        bezCurve[1] = maths.addArrays(firstPoint, maths.mulItems(leftTangent,  alpha_l));
        bezCurve[2] = maths.addArrays(lastPoint,  maths.mulItems(rightTangent, alpha_r));
    }

    return bezCurve;
};

/**
 * Given set of points and their parameterization, try to find a better parameterization.
 *
 * @param {Array<Array<Number>>} bezier - Current fitted curve
 * @param {Array<Array<Number>>} points - Array of digitized points
 * @param {Array<Number>} parameters - Current parameter values
 * @returns {Array<Number>} New parameter values
 */
function reparameterize(bezier, points, parameters) {
    /*
    var j, len, point, results, u;
    results = [];
    for (j = 0, len = points.length; j < len; j++) {
        point = points[j], u = parameters[j];

        results.push(newtonRaphsonRootFind(bezier, point, u));
    }
    return results;
    //*/
    return parameters.map((p, i) => newtonRaphsonRootFind(bezier, points[i], p));
};

/**
 * Use Newton-Raphson iteration to find better root.
 *
 * @param {Array<Array<Number>>} bez - Current fitted curve
 * @param {Array<Number>} point - Digitized point
 * @param {Number} u - Parameter value for "P"
 * @returns {Number} New u
 */
function newtonRaphsonRootFind(bez, point, u) {
    /*
        Newton's root finding algorithm calculates f(x)=0 by reiterating
        x_n+1 = x_n - f(x_n)/f'(x_n)
        We are trying to find curve parameter u for some point p that minimizes
        the distance from that point to the curve. Distance point to curve is d=q(u)-p.
        At minimum distance the point is perpendicular to the curve.
        We are solving
        f = q(u)-p * q'(u) = 0
        with
        f' = q'(u) * q'(u) + q(u)-p * q''(u)
        gives
        u_n+1 = u_n - |q(u_n)-p * q'(u_n)| / |q'(u_n)**2 + q(u_n)-p * q''(u_n)|
    */

    var d = maths.subtract(bezier.q(bez, u), point),
        qprime = bezier.qprime(bez, u),
        numerator = /*sum(*/maths.mulMatrix(d, qprime)/*)*/,
        denominator = maths.sum(maths.addItems( maths.squareItems(qprime), maths.mulMatrix(d, bezier.qprimeprime(bez, u)) ));

    if (denominator === 0) {
        return u;
    } else {
        return u - (numerator/denominator);
    }
};

/**
 * Assign parameter values to digitized points using relative distances between points.
 *
 * @param {Array<Array<Number>>} points - Array of digitized points
 * @returns {Array<Number>} Parameter values
 */
function chordLengthParameterize(points) {
    var u = [], currU, prevU, prevP;

    points.forEach((p, i) => {
        currU = i ? prevU + maths.vectorLen(maths.subtract(p, prevP))
                  : 0;
        u.push(currU);

        prevU = currU;
        prevP = p;
    })
    u = u.map(x => x/prevU);

    return u;
};

/**
 * Find the maximum squared distance of digitized points to fitted curve.
 *
 * @param {Array<Array<Number>>} points - Array of digitized points
 * @param {Array<Array<Number>>} bez - Fitted curve
 * @param {Array<Number>} parameters - Parameterization of points
 * @returns {Array<Number>} Maximum error (squared) and point of max error
 */
function computeMaxError(points, bez, parameters) {
    var dist,       //Current error
        maxDist,    //Maximum error
        splitPoint, //Point of maximum error
        v,          //Vector from point to curve
        i, count, point, t;

    maxDist = 0;
    splitPoint = points.length / 2;

    const t_distMap = mapTtoRelativeDistances(bez, 10);

    for (i = 0, count = points.length; i < count; i++) {
        point = points[i];
        //Find 't' for a point on the bez curve that's as close to 'point' as possible:
        t = find_t(bez, parameters[i], t_distMap, 10);

        v = maths.subtract(bezier.q(bez, t), point);
        dist = v[0]*v[0] + v[1]*v[1];

        if (dist > maxDist) {
            maxDist = dist;
            splitPoint = i;
        }
    }

    return [maxDist, splitPoint];
};

//Sample 't's and map them to relative distances along the curve:
var mapTtoRelativeDistances = function (bez, B_parts) {
    var B_t_curr;
    var B_t_dist = [0];
    var B_t_prev = bez[0];
    var sumLen = 0;

    for (var i=1; i<=B_parts; i++) {
      B_t_curr = bezier.q(bez, i/B_parts);

      sumLen += maths.vectorLen(maths.subtract(B_t_curr, B_t_prev));

      B_t_dist.push(sumLen);
      B_t_prev = B_t_curr;
    }

    //Normalize B_length to the same interval as the parameter distances; 0 to 1:
    B_t_dist = B_t_dist.map(x => x/sumLen);
    return B_t_dist;
};

function find_t(bez, param, t_distMap, B_parts) {
    if(param < 0) { return 0; }
    if(param > 1) { return 1; }

    /*
        'param' is a value between 0 and 1 telling us the relative position
        of a point on the source polyline (linearly from the start (0) to the end (1)).
        To see if a given curve - 'bez' - is a close approximation of the polyline,
        we compare such a poly-point to the point on the curve that's the same
        relative distance along the curve's length.

        But finding that curve-point takes a little work:
        There is a function "B(t)" to find points along a curve from the parametric parameter 't'
        (also relative from 0 to 1: http://stackoverflow.com/a/32841764/1869660
                                    http://pomax.github.io/bezierinfo/#explanation),
        but 't' isn't linear by length (http://gamedev.stackexchange.com/questions/105230).

        So, we sample some points along the curve using a handful of values for 't'.
        Then, we calculate the length between those samples via plain euclidean distance;
        B(t) concentrates the points around sharp turns, so this should give us a good-enough outline of the curve.
        Thus, for a given relative distance ('param'), we can now find an upper and lower value
        for the corresponding 't' by searching through those sampled distances.
        Finally, we just use linear interpolation to find a better value for the exact 't'.

        More info:
            http://gamedev.stackexchange.com/questions/105230/points-evenly-spaced-along-a-bezier-curve
            http://stackoverflow.com/questions/29438398/cheap-way-of-calculating-cubic-bezier-length
            http://steve.hollasch.net/cgindex/curves/cbezarclen.html
            https://github.com/retuxx/tinyspline
    */
    var lenMax, lenMin, tMax, tMin, t;

    //Find the two t-s that the current param distance lies between,
    //and then interpolate a somewhat accurate value for the exact t:
    for(var i = 1; i <= B_parts; i++) {

        if(param <= t_distMap[i]) {
            tMin   = (i-1) / B_parts;
            tMax   = i / B_parts;
            lenMin = t_distMap[i-1];
            lenMax = t_distMap[i];

            t = (param-lenMin)/(lenMax-lenMin) * (tMax-tMin) + tMin;
            break;
        }
    }
    return t;
}

/**
 * Creates a vector of length 1 which shows the direction from B to A
 */
function createTangent(pointA, pointB) {
    return maths.normalize(maths.subtract(pointA, pointB));
}

/*
    Simplified versions of what we need from math.js
    Optimized for our input, which is only numbers and 1x2 arrays (i.e. [x, y] coordinates).
*/
class maths {
    //zeros = logAndRun(math.zeros);
    static zeros_Xx2x2(x) {
        var zs = [];
        while(x--) { zs.push([0,0]); }
        return zs
    }

    //multiply = logAndRun(math.multiply);
    static mulItems(items, multiplier) {
        //return items.map(x => x*multiplier);
        return [items[0]*multiplier, items[1]*multiplier];
    }
    static mulMatrix(m1, m2) {
        //https://en.wikipedia.org/wiki/Matrix_multiplication#Matrix_product_.28two_matrices.29
        //Simplified to only handle 1-dimensional matrices (i.e. arrays) of equal length:
        //  return m1.reduce((sum,x1,i) => sum + (x1*m2[i]),
        //                   0);
        return (m1[0]*m2[0]) + (m1[1]*m2[1]);
    }

    //Only used to subract to points (or at least arrays):
    //  subtract = logAndRun(math.subtract);
    static subtract(arr1, arr2) {
        //return arr1.map((x1, i) => x1 - arr2[i]);
        return [arr1[0]-arr2[0], arr1[1]-arr2[1]];
    }

    //add = logAndRun(math.add);
    static addArrays(arr1, arr2) {
        //return arr1.map((x1, i) => x1 + arr2[i]);
        return [arr1[0]+arr2[0], arr1[1]+arr2[1]];
    }
    static addItems(items, addition) {
        //return items.map(x => x+addition);
        return [items[0]+addition, items[1]+addition];
    }

    //var sum = logAndRun(math.sum);
    static sum(items) {
        return items.reduce((sum,x) => sum + x);
    }

    //chain = math.chain;

    //Only used on two arrays. The dot product is equal to the matrix product in this case:
    //  dot = logAndRun(math.dot);
    static dot(m1, m2) {
        return maths.mulMatrix(m1, m2);
    }

    //https://en.wikipedia.org/wiki/Norm_(mathematics)#Euclidean_norm
    //  var norm = logAndRun(math.norm);
    static vectorLen(v) {
        var a = v[0], b = v[1];
        return Math.sqrt(a*a + b*b);
    }

    //math.divide = logAndRun(math.divide);
    static divItems(items, divisor) {
        //return items.map(x => x/divisor);
        return [items[0]/divisor, items[1]/divisor];
    }

    //var dotPow = logAndRun(math.dotPow);
    static squareItems(items) {
        //return items.map(x => x*x);
        var a = items[0], b = items[1];
        return [a*a, b*b];
    }

    static normalize(v) {
        return this.divItems(v, this.vectorLen(v));
    }

    //Math.pow = logAndRun(Math.pow);
}


class bezier {
    //Evaluates cubic bezier at t, return point
    static q(ctrlPoly, t) {
        var tx = 1.0 - t;
        var pA = maths.mulItems( ctrlPoly[0],      tx * tx * tx ),
            pB = maths.mulItems( ctrlPoly[1],  3 * tx * tx *  t ),
            pC = maths.mulItems( ctrlPoly[2],  3 * tx *  t *  t ),
            pD = maths.mulItems( ctrlPoly[3],       t *  t *  t );
        return maths.addArrays(maths.addArrays(pA, pB), maths.addArrays(pC, pD));
    }

    //Evaluates cubic bezier first derivative at t, return point
    static qprime(ctrlPoly, t) {
        var tx = 1.0 - t;
        var pA = maths.mulItems( maths.subtract(ctrlPoly[1], ctrlPoly[0]),  3 * tx * tx ),
            pB = maths.mulItems( maths.subtract(ctrlPoly[2], ctrlPoly[1]),  6 * tx *  t ),
            pC = maths.mulItems( maths.subtract(ctrlPoly[3], ctrlPoly[2]),  3 *  t *  t );
        return maths.addArrays(maths.addArrays(pA, pB), pC);
    }

    //Evaluates cubic bezier second derivative at t, return point
    static qprimeprime(ctrlPoly, t) {
        return maths.addArrays(maths.mulItems( maths.addArrays(maths.subtract(ctrlPoly[2], maths.mulItems(ctrlPoly[1], 2)), ctrlPoly[0]),  6 * (1.0 - t) ),
                               maths.mulItems( maths.addArrays(maths.subtract(ctrlPoly[3], maths.mulItems(ctrlPoly[2], 2)), ctrlPoly[1]),  6 *        t  ));
    }
}

module.exports = fitCurve;

},{}],2:[function(require,module,exports){
/*!
* svg.js - A lightweight library for manipulating and animating SVG.
* @version 2.3.7
* https://svgdotjs.github.io/
*
* @copyright Wout Fierens <wout@mick-wout.com>
* @license MIT
*
* BUILT: Sat Jan 14 2017 07:23:18 GMT+0100 (CET)
*/;
(function(root, factory) {
  if (typeof define === 'function' && define.amd) {
    define(function(){
      return factory(root, root.document)
    })
  } else if (typeof exports === 'object') {
    module.exports = root.document ? factory(root, root.document) : function(w){ return factory(w, w.document) }
  } else {
    root.SVG = factory(root, root.document)
  }
}(typeof window !== "undefined" ? window : this, function(window, document) {

// The main wrapping element
var SVG = this.SVG = function(element) {
  if (SVG.supported) {
    element = new SVG.Doc(element)

    if(!SVG.parser.draw)
      SVG.prepare()

    return element
  }
}

// Default namespaces
SVG.ns    = 'http://www.w3.org/2000/svg'
SVG.xmlns = 'http://www.w3.org/2000/xmlns/'
SVG.xlink = 'http://www.w3.org/1999/xlink'
SVG.svgjs = 'http://svgjs.com/svgjs'

// Svg support test
SVG.supported = (function() {
  return !! document.createElementNS &&
         !! document.createElementNS(SVG.ns,'svg').createSVGRect
})()

// Don't bother to continue if SVG is not supported
if (!SVG.supported) return false

// Element id sequence
SVG.did  = 1000

// Get next named element id
SVG.eid = function(name) {
  return 'Svgjs' + capitalize(name) + (SVG.did++)
}

// Method for element creation
SVG.create = function(name) {
  // create element
  var element = document.createElementNS(this.ns, name)

  // apply unique id
  element.setAttribute('id', this.eid(name))

  return element
}

// Method for extending objects
SVG.extend = function() {
  var modules, methods, key, i

  // Get list of modules
  modules = [].slice.call(arguments)

  // Get object with extensions
  methods = modules.pop()

  for (i = modules.length - 1; i >= 0; i--)
    if (modules[i])
      for (key in methods)
        modules[i].prototype[key] = methods[key]

  // Make sure SVG.Set inherits any newly added methods
  if (SVG.Set && SVG.Set.inherit)
    SVG.Set.inherit()
}

// Invent new element
SVG.invent = function(config) {
  // Create element initializer
  var initializer = typeof config.create == 'function' ?
    config.create :
    function() {
      this.constructor.call(this, SVG.create(config.create))
    }

  // Inherit prototype
  if (config.inherit)
    initializer.prototype = new config.inherit

  // Extend with methods
  if (config.extend)
    SVG.extend(initializer, config.extend)

  // Attach construct method to parent
  if (config.construct)
    SVG.extend(config.parent || SVG.Container, config.construct)

  return initializer
}

// Adopt existing svg elements
SVG.adopt = function(node) {
  // check for presence of node
  if (!node) return null

  // make sure a node isn't already adopted
  if (node.instance) return node.instance

  // initialize variables
  var element

  // adopt with element-specific settings
  if (node.nodeName == 'svg')
    element = node.parentNode instanceof SVGElement ? new SVG.Nested : new SVG.Doc
  else if (node.nodeName == 'linearGradient')
    element = new SVG.Gradient('linear')
  else if (node.nodeName == 'radialGradient')
    element = new SVG.Gradient('radial')
  else if (SVG[capitalize(node.nodeName)])
    element = new SVG[capitalize(node.nodeName)]
  else
    element = new SVG.Element(node)

  // ensure references
  element.type  = node.nodeName
  element.node  = node
  node.instance = element

  // SVG.Class specific preparations
  if (element instanceof SVG.Doc)
    element.namespace().defs()

  // pull svgjs data from the dom (getAttributeNS doesn't work in html5)
  element.setData(JSON.parse(node.getAttribute('svgjs:data')) || {})

  return element
}

// Initialize parsing element
SVG.prepare = function() {
  // Select document body and create invisible svg element
  var body = document.getElementsByTagName('body')[0]
    , draw = (body ? new SVG.Doc(body) :  new SVG.Doc(document.documentElement).nested()).size(2, 0)

  // Create parser object
  SVG.parser = {
    body: body || document.documentElement
  , draw: draw.style('opacity:0;position:fixed;left:100%;top:100%;overflow:hidden')
  , poly: draw.polyline().node
  , path: draw.path().node
  , native: SVG.create('svg')
  }
}

SVG.parser = {
  native: SVG.create('svg')
}

document.addEventListener('DOMContentLoaded', function() {
  if(!SVG.parser.draw)
    SVG.prepare()
}, false)

// Storage for regular expressions
SVG.regex = {
  // Parse unit value
  numberAndUnit:    /^([+-]?(\d+(\.\d*)?|\.\d+)(e[+-]?\d+)?)([a-z%]*)$/i

  // Parse hex value
, hex:              /^#?([a-f\d]{2})([a-f\d]{2})([a-f\d]{2})$/i

  // Parse rgb value
, rgb:              /rgb\((\d+),(\d+),(\d+)\)/

  // Parse reference id
, reference:        /#([a-z0-9\-_]+)/i

  // Parse matrix wrapper
, matrix:           /matrix\(|\)/g

  // Elements of a matrix
, matrixElements:   /,*\s+|,/

  // Whitespace
, whitespace:       /\s/g

  // Test hex value
, isHex:            /^#[a-f0-9]{3,6}$/i

  // Test rgb value
, isRgb:            /^rgb\(/

  // Test css declaration
, isCss:            /[^:]+:[^;]+;?/

  // Test for blank string
, isBlank:          /^(\s+)?$/

  // Test for numeric string
, isNumber:         /^[+-]?(\d+(\.\d*)?|\.\d+)(e[+-]?\d+)?$/i

  // Test for percent value
, isPercent:        /^-?[\d\.]+%$/

  // Test for image url
, isImage:          /\.(jpg|jpeg|png|gif|svg)(\?[^=]+.*)?/i

  // The following regex are used to parse the d attribute of a path

  // Replaces all negative exponents
, negExp:           /e\-/gi

  // Replaces all comma
, comma:            /,/g

  // Replaces all hyphens
, hyphen:           /\-/g

  // Replaces and tests for all path letters
, pathLetters:      /[MLHVCSQTAZ]/gi

  // yes we need this one, too
, isPathLetter:     /[MLHVCSQTAZ]/i

  // split at whitespaces
, whitespaces:      /\s+/

  // matches X
, X:                /X/g
}

SVG.utils = {
  // Map function
  map: function(array, block) {
    var i
      , il = array.length
      , result = []

    for (i = 0; i < il; i++)
      result.push(block(array[i]))

    return result
  }

  // Filter function
, filter: function(array, block) {
    var i
      , il = array.length
      , result = []

    for (i = 0; i < il; i++)
      if (block(array[i]))
        result.push(array[i])

    return result
  }

  // Degrees to radians
, radians: function(d) {
    return d % 360 * Math.PI / 180
  }

  // Radians to degrees
, degrees: function(r) {
    return r * 180 / Math.PI % 360
  }

, filterSVGElements: function(nodes) {
    return this.filter( nodes, function(el) { return el instanceof SVGElement })
  }

}

SVG.defaults = {
  // Default attribute values
  attrs: {
    // fill and stroke
    'fill-opacity':     1
  , 'stroke-opacity':   1
  , 'stroke-width':     0
  , 'stroke-linejoin':  'miter'
  , 'stroke-linecap':   'butt'
  , fill:               '#000000'
  , stroke:             '#000000'
  , opacity:            1
    // position
  , x:                  0
  , y:                  0
  , cx:                 0
  , cy:                 0
    // size
  , width:              0
  , height:             0
    // radius
  , r:                  0
  , rx:                 0
  , ry:                 0
    // gradient
  , offset:             0
  , 'stop-opacity':     1
  , 'stop-color':       '#000000'
    // text
  , 'font-size':        16
  , 'font-family':      'Helvetica, Arial, sans-serif'
  , 'text-anchor':      'start'
  }

}
// Module for color convertions
SVG.Color = function(color) {
  var match

  // initialize defaults
  this.r = 0
  this.g = 0
  this.b = 0

  if(!color) return

  // parse color
  if (typeof color === 'string') {
    if (SVG.regex.isRgb.test(color)) {
      // get rgb values
      match = SVG.regex.rgb.exec(color.replace(/\s/g,''))

      // parse numeric values
      this.r = parseInt(match[1])
      this.g = parseInt(match[2])
      this.b = parseInt(match[3])

    } else if (SVG.regex.isHex.test(color)) {
      // get hex values
      match = SVG.regex.hex.exec(fullHex(color))

      // parse numeric values
      this.r = parseInt(match[1], 16)
      this.g = parseInt(match[2], 16)
      this.b = parseInt(match[3], 16)

    }

  } else if (typeof color === 'object') {
    this.r = color.r
    this.g = color.g
    this.b = color.b

  }

}

SVG.extend(SVG.Color, {
  // Default to hex conversion
  toString: function() {
    return this.toHex()
  }
  // Build hex value
, toHex: function() {
    return '#'
      + compToHex(this.r)
      + compToHex(this.g)
      + compToHex(this.b)
  }
  // Build rgb value
, toRgb: function() {
    return 'rgb(' + [this.r, this.g, this.b].join() + ')'
  }
  // Calculate true brightness
, brightness: function() {
    return (this.r / 255 * 0.30)
         + (this.g / 255 * 0.59)
         + (this.b / 255 * 0.11)
  }
  // Make color morphable
, morph: function(color) {
    this.destination = new SVG.Color(color)

    return this
  }
  // Get morphed color at given position
, at: function(pos) {
    // make sure a destination is defined
    if (!this.destination) return this

    // normalise pos
    pos = pos < 0 ? 0 : pos > 1 ? 1 : pos

    // generate morphed color
    return new SVG.Color({
      r: ~~(this.r + (this.destination.r - this.r) * pos)
    , g: ~~(this.g + (this.destination.g - this.g) * pos)
    , b: ~~(this.b + (this.destination.b - this.b) * pos)
    })
  }

})

// Testers

// Test if given value is a color string
SVG.Color.test = function(color) {
  color += ''
  return SVG.regex.isHex.test(color)
      || SVG.regex.isRgb.test(color)
}

// Test if given value is a rgb object
SVG.Color.isRgb = function(color) {
  return color && typeof color.r == 'number'
               && typeof color.g == 'number'
               && typeof color.b == 'number'
}

// Test if given value is a color
SVG.Color.isColor = function(color) {
  return SVG.Color.isRgb(color) || SVG.Color.test(color)
}
// Module for array conversion
SVG.Array = function(array, fallback) {
  array = (array || []).valueOf()

  // if array is empty and fallback is provided, use fallback
  if (array.length == 0 && fallback)
    array = fallback.valueOf()

  // parse array
  this.value = this.parse(array)
}

SVG.extend(SVG.Array, {
  // Make array morphable
  morph: function(array) {
    this.destination = this.parse(array)

    // normalize length of arrays
    if (this.value.length != this.destination.length) {
      var lastValue       = this.value[this.value.length - 1]
        , lastDestination = this.destination[this.destination.length - 1]

      while(this.value.length > this.destination.length)
        this.destination.push(lastDestination)
      while(this.value.length < this.destination.length)
        this.value.push(lastValue)
    }

    return this
  }
  // Clean up any duplicate points
, settle: function() {
    // find all unique values
    for (var i = 0, il = this.value.length, seen = []; i < il; i++)
      if (seen.indexOf(this.value[i]) == -1)
        seen.push(this.value[i])

    // set new value
    return this.value = seen
  }
  // Get morphed array at given position
, at: function(pos) {
    // make sure a destination is defined
    if (!this.destination) return this

    // generate morphed array
    for (var i = 0, il = this.value.length, array = []; i < il; i++)
      array.push(this.value[i] + (this.destination[i] - this.value[i]) * pos)

    return new SVG.Array(array)
  }
  // Convert array to string
, toString: function() {
    return this.value.join(' ')
  }
  // Real value
, valueOf: function() {
    return this.value
  }
  // Parse whitespace separated string
, parse: function(array) {
    array = array.valueOf()

    // if already is an array, no need to parse it
    if (Array.isArray(array)) return array

    return this.split(array)
  }
  // Strip unnecessary whitespace
, split: function(string) {
    return string.trim().split(/\s+/)
  }
  // Reverse array
, reverse: function() {
    this.value.reverse()

    return this
  }

})
// Poly points array
SVG.PointArray = function(array, fallback) {
  this.constructor.call(this, array, fallback || [[0,0]])
}

// Inherit from SVG.Array
SVG.PointArray.prototype = new SVG.Array

SVG.extend(SVG.PointArray, {
  // Convert array to string
  toString: function() {
    // convert to a poly point string
    for (var i = 0, il = this.value.length, array = []; i < il; i++)
      array.push(this.value[i].join(','))

    return array.join(' ')
  }
  // Convert array to line object
, toLine: function() {
    return {
      x1: this.value[0][0]
    , y1: this.value[0][1]
    , x2: this.value[1][0]
    , y2: this.value[1][1]
    }
  }
  // Get morphed array at given position
, at: function(pos) {
    // make sure a destination is defined
    if (!this.destination) return this

    // generate morphed point string
    for (var i = 0, il = this.value.length, array = []; i < il; i++)
      array.push([
        this.value[i][0] + (this.destination[i][0] - this.value[i][0]) * pos
      , this.value[i][1] + (this.destination[i][1] - this.value[i][1]) * pos
      ])

    return new SVG.PointArray(array)
  }
  // Parse point string
, parse: function(array) {
    var points = []

    array = array.valueOf()

    // if already is an array, no need to parse it
    if (Array.isArray(array)) return array

    // parse points
    array = array.trim().split(/\s+|,/)

    // validate points - https://svgwg.org/svg2-draft/shapes.html#DataTypePoints
    // Odd number of coordinates is an error. In such cases, drop the last odd coordinate.
    if (array.length % 2 !== 0) array.pop()

    // wrap points in two-tuples and parse points as floats
    for(var i = 0, len = array.length; i < len; i = i + 2)
      points.push([ parseFloat(array[i]), parseFloat(array[i+1]) ])

    return points
  }
  // Move point string
, move: function(x, y) {
    var box = this.bbox()

    // get relative offset
    x -= box.x
    y -= box.y

    // move every point
    if (!isNaN(x) && !isNaN(y))
      for (var i = this.value.length - 1; i >= 0; i--)
        this.value[i] = [this.value[i][0] + x, this.value[i][1] + y]

    return this
  }
  // Resize poly string
, size: function(width, height) {
    var i, box = this.bbox()

    // recalculate position of all points according to new size
    for (i = this.value.length - 1; i >= 0; i--) {
      this.value[i][0] = ((this.value[i][0] - box.x) * width)  / box.width  + box.x
      this.value[i][1] = ((this.value[i][1] - box.y) * height) / box.height + box.y
    }

    return this
  }
  // Get bounding box of points
, bbox: function() {
    SVG.parser.poly.setAttribute('points', this.toString())

    return SVG.parser.poly.getBBox()
  }

})
// Path points array
SVG.PathArray = function(array, fallback) {
  this.constructor.call(this, array, fallback || [['M', 0, 0]])
}

// Inherit from SVG.Array
SVG.PathArray.prototype = new SVG.Array

SVG.extend(SVG.PathArray, {
  // Convert array to string
  toString: function() {
    return arrayToString(this.value)
  }
  // Move path string
, move: function(x, y) {
    // get bounding box of current situation
    var box = this.bbox()

    // get relative offset
    x -= box.x
    y -= box.y

    if (!isNaN(x) && !isNaN(y)) {
      // move every point
      for (var l, i = this.value.length - 1; i >= 0; i--) {
        l = this.value[i][0]

        if (l == 'M' || l == 'L' || l == 'T')  {
          this.value[i][1] += x
          this.value[i][2] += y

        } else if (l == 'H')  {
          this.value[i][1] += x

        } else if (l == 'V')  {
          this.value[i][1] += y

        } else if (l == 'C' || l == 'S' || l == 'Q')  {
          this.value[i][1] += x
          this.value[i][2] += y
          this.value[i][3] += x
          this.value[i][4] += y

          if (l == 'C')  {
            this.value[i][5] += x
            this.value[i][6] += y
          }

        } else if (l == 'A')  {
          this.value[i][6] += x
          this.value[i][7] += y
        }

      }
    }

    return this
  }
  // Resize path string
, size: function(width, height) {
    // get bounding box of current situation
    var i, l, box = this.bbox()

    // recalculate position of all points according to new size
    for (i = this.value.length - 1; i >= 0; i--) {
      l = this.value[i][0]

      if (l == 'M' || l == 'L' || l == 'T')  {
        this.value[i][1] = ((this.value[i][1] - box.x) * width)  / box.width  + box.x
        this.value[i][2] = ((this.value[i][2] - box.y) * height) / box.height + box.y

      } else if (l == 'H')  {
        this.value[i][1] = ((this.value[i][1] - box.x) * width)  / box.width  + box.x

      } else if (l == 'V')  {
        this.value[i][1] = ((this.value[i][1] - box.y) * height) / box.height + box.y

      } else if (l == 'C' || l == 'S' || l == 'Q')  {
        this.value[i][1] = ((this.value[i][1] - box.x) * width)  / box.width  + box.x
        this.value[i][2] = ((this.value[i][2] - box.y) * height) / box.height + box.y
        this.value[i][3] = ((this.value[i][3] - box.x) * width)  / box.width  + box.x
        this.value[i][4] = ((this.value[i][4] - box.y) * height) / box.height + box.y

        if (l == 'C')  {
          this.value[i][5] = ((this.value[i][5] - box.x) * width)  / box.width  + box.x
          this.value[i][6] = ((this.value[i][6] - box.y) * height) / box.height + box.y
        }

      } else if (l == 'A')  {
        // resize radii
        this.value[i][1] = (this.value[i][1] * width)  / box.width
        this.value[i][2] = (this.value[i][2] * height) / box.height

        // move position values
        this.value[i][6] = ((this.value[i][6] - box.x) * width)  / box.width  + box.x
        this.value[i][7] = ((this.value[i][7] - box.y) * height) / box.height + box.y
      }

    }

    return this
  }
  // Test if the passed path array use the same path data commands as this path array
, equalCommands: function(pathArray) {
    var i, il, equalCommands

    pathArray = new SVG.PathArray(pathArray)

    equalCommands = this.value.length === pathArray.value.length
    for(i = 0, il = this.value.length; equalCommands && i < il; i++) {
      equalCommands = this.value[i][0] === pathArray.value[i][0]
    }

    return equalCommands
  }
  // Make path array morphable
, morph: function(pathArray) {
    pathArray = new SVG.PathArray(pathArray)

    if(this.equalCommands(pathArray)) {
      this.destination = pathArray
    } else {
      this.destination = null
    }

    return this
  }
  // Get morphed path array at given position
, at: function(pos) {
    // make sure a destination is defined
    if (!this.destination) return this

    var sourceArray = this.value
      , destinationArray = this.destination.value
      , array = [], pathArray = new SVG.PathArray()
      , i, il, j, jl

    // Animate has specified in the SVG spec
    // See: https://www.w3.org/TR/SVG11/paths.html#PathElement
    for (i = 0, il = sourceArray.length; i < il; i++) {
      array[i] = [sourceArray[i][0]]
      for(j = 1, jl = sourceArray[i].length; j < jl; j++) {
        array[i][j] = sourceArray[i][j] + (destinationArray[i][j] - sourceArray[i][j]) * pos
      }
      // For the two flags of the elliptical arc command, the SVG spec say:
      // Flags and booleans are interpolated as fractions between zero and one, with any non-zero value considered to be a value of one/true
      // Elliptical arc command as an array followed by corresponding indexes:
      // ['A', rx, ry, x-axis-rotation, large-arc-flag, sweep-flag, x, y]
      //   0    1   2        3                 4             5      6  7
      if(array[i][0] === 'A') {
        array[i][4] = +(array[i][4] != 0)
        array[i][5] = +(array[i][5] != 0)
      }
    }

    // Directly modify the value of a path array, this is done this way for performance
    pathArray.value = array
    return pathArray
  }
  // Absolutize and parse path to array
, parse: function(array) {
    // if it's already a patharray, no need to parse it
    if (array instanceof SVG.PathArray) return array.valueOf()

    // prepare for parsing
    var i, x0, y0, s, seg, arr
      , x = 0
      , y = 0
      , paramCnt = { 'M':2, 'L':2, 'H':1, 'V':1, 'C':6, 'S':4, 'Q':4, 'T':2, 'A':7 }

    if(typeof array == 'string'){

      array = array
        .replace(SVG.regex.negExp, 'X')         // replace all negative exponents with certain char
        .replace(SVG.regex.pathLetters, ' $& ') // put some room between letters and numbers
        .replace(SVG.regex.hyphen, ' -')        // add space before hyphen
        .replace(SVG.regex.comma, ' ')          // unify all spaces
        .replace(SVG.regex.X, 'e-')             // add back the expoent
        .trim()                                 // trim
        .split(SVG.regex.whitespaces)           // split into array

      // at this place there could be parts like ['3.124.854.32'] because we could not determine the point as seperator till now
      // we fix this elements in the next loop
      for(i = array.length; --i;){
        if(array[i].indexOf('.') != array[i].lastIndexOf('.')){
          var split = array[i].split('.') // split at the point
          var first = [split.shift(), split.shift()].join('.') // join the first number together
          array.splice.apply(array, [i, 1].concat(first, split.map(function(el){ return '.'+el }))) // add first and all other entries back to array
        }
      }

    }else{
      array = array.reduce(function(prev, curr){
        return [].concat.apply(prev, curr)
      }, [])
    }

    // array now is an array containing all parts of a path e.g. ['M', '0', '0', 'L', '30', '30' ...]

    var arr = []

    do{

      // Test if we have a path letter
      if(SVG.regex.isPathLetter.test(array[0])){
        s = array[0]
        array.shift()
      // If last letter was a move command and we got no new, it defaults to [L]ine
      }else if(s == 'M'){
        s = 'L'
      }else if(s == 'm'){
        s = 'l'
      }

      // add path letter as first element
      seg = [s.toUpperCase()]

      // push all necessary parameters to segment
      for(i = 0; i < paramCnt[seg[0]]; ++i){
        seg.push(parseFloat(array.shift()))
      }

      // upper case
      if(s == seg[0]){

        if(s == 'M' || s == 'L' || s == 'C' || s == 'Q' || s == 'S' || s == 'T'){
          x = seg[paramCnt[seg[0]]-1]
          y = seg[paramCnt[seg[0]]]
        }else if(s == 'V'){
          y = seg[1]
        }else if(s == 'H'){
          x = seg[1]
        }else if(s == 'A'){
          x = seg[6]
          y = seg[7]
        }

      // lower case
      }else{

        // convert relative to absolute values
        if(s == 'm' || s == 'l' || s == 'c' || s == 's' || s == 'q' || s == 't'){

          seg[1] += x
          seg[2] += y

          if(seg[3] != null){
            seg[3] += x
            seg[4] += y
          }

          if(seg[5] != null){
            seg[5] += x
            seg[6] += y
          }

          // move pointer
          x = seg[paramCnt[seg[0]]-1]
          y = seg[paramCnt[seg[0]]]

        }else if(s == 'v'){
          seg[1] += y
          y = seg[1]
        }else if(s == 'h'){
          seg[1] += x
          x = seg[1]
        }else if(s == 'a'){
          seg[6] += x
          seg[7] += y
          x = seg[6]
          y = seg[7]
        }

      }

      if(seg[0] == 'M'){
        x0 = x
        y0 = y
      }

      if(seg[0] == 'Z'){
        x = x0
        y = y0
      }

      arr.push(seg)

    }while(array.length)

    return arr

  }
  // Get bounding box of path
, bbox: function() {
    SVG.parser.path.setAttribute('d', this.toString())

    return SVG.parser.path.getBBox()
  }

})

// Module for unit convertions
SVG.Number = SVG.invent({
  // Initialize
  create: function(value, unit) {
    // initialize defaults
    this.value = 0
    this.unit  = unit || ''

    // parse value
    if (typeof value === 'number') {
      // ensure a valid numeric value
      this.value = isNaN(value) ? 0 : !isFinite(value) ? (value < 0 ? -3.4e+38 : +3.4e+38) : value

    } else if (typeof value === 'string') {
      unit = value.match(SVG.regex.numberAndUnit)

      if (unit) {
        // make value numeric
        this.value = parseFloat(unit[1])

        // normalize
        if (unit[5] == '%')
          this.value /= 100
        else if (unit[5] == 's')
          this.value *= 1000

        // store unit
        this.unit = unit[5]
      }

    } else {
      if (value instanceof SVG.Number) {
        this.value = value.valueOf()
        this.unit  = value.unit
      }
    }

  }
  // Add methods
, extend: {
    // Stringalize
    toString: function() {
      return (
        this.unit == '%' ?
          ~~(this.value * 1e8) / 1e6:
        this.unit == 's' ?
          this.value / 1e3 :
          this.value
      ) + this.unit
    }
  , toJSON: function() {
      return this.toString()
    }
  , // Convert to primitive
    valueOf: function() {
      return this.value
    }
    // Add number
  , plus: function(number) {
      return new SVG.Number(this + new SVG.Number(number), this.unit)
    }
    // Subtract number
  , minus: function(number) {
      return this.plus(-new SVG.Number(number))
    }
    // Multiply number
  , times: function(number) {
      return new SVG.Number(this * new SVG.Number(number), this.unit)
    }
    // Divide number
  , divide: function(number) {
      return new SVG.Number(this / new SVG.Number(number), this.unit)
    }
    // Convert to different unit
  , to: function(unit) {
      var number = new SVG.Number(this)

      if (typeof unit === 'string')
        number.unit = unit

      return number
    }
    // Make number morphable
  , morph: function(number) {
      this.destination = new SVG.Number(number)

      return this
    }
    // Get morphed number at given position
  , at: function(pos) {
      // Make sure a destination is defined
      if (!this.destination) return this

      // Generate new morphed number
      return new SVG.Number(this.destination)
          .minus(this)
          .times(pos)
          .plus(this)
    }

  }
})

SVG.Element = SVG.invent({
  // Initialize node
  create: function(node) {
    // make stroke value accessible dynamically
    this._stroke = SVG.defaults.attrs.stroke

    // initialize data object
    this.dom = {}

    // create circular reference
    if (this.node = node) {
      this.type = node.nodeName
      this.node.instance = this

      // store current attribute value
      this._stroke = node.getAttribute('stroke') || this._stroke
    }
  }

  // Add class methods
, extend: {
    // Move over x-axis
    x: function(x) {
      return this.attr('x', x)
    }
    // Move over y-axis
  , y: function(y) {
      return this.attr('y', y)
    }
    // Move by center over x-axis
  , cx: function(x) {
      return x == null ? this.x() + this.width() / 2 : this.x(x - this.width() / 2)
    }
    // Move by center over y-axis
  , cy: function(y) {
      return y == null ? this.y() + this.height() / 2 : this.y(y - this.height() / 2)
    }
    // Move element to given x and y values
  , move: function(x, y) {
      return this.x(x).y(y)
    }
    // Move element by its center
  , center: function(x, y) {
      return this.cx(x).cy(y)
    }
    // Set width of element
  , width: function(width) {
      return this.attr('width', width)
    }
    // Set height of element
  , height: function(height) {
      return this.attr('height', height)
    }
    // Set element size to given width and height
  , size: function(width, height) {
      var p = proportionalSize(this, width, height)

      return this
        .width(new SVG.Number(p.width))
        .height(new SVG.Number(p.height))
    }
    // Clone element
  , clone: function(parent) {
      // clone element and assign new id
      var clone = assignNewId(this.node.cloneNode(true))

      // insert the clone in the given parent or after myself
      if(parent) parent.add(clone)
      else this.after(clone)

      return clone
    }
    // Remove element
  , remove: function() {
      if (this.parent())
        this.parent().removeElement(this)

      return this
    }
    // Replace element
  , replace: function(element) {
      this.after(element).remove()

      return element
    }
    // Add element to given container and return self
  , addTo: function(parent) {
      return parent.put(this)
    }
    // Add element to given container and return container
  , putIn: function(parent) {
      return parent.add(this)
    }
    // Get / set id
  , id: function(id) {
      return this.attr('id', id)
    }
    // Checks whether the given point inside the bounding box of the element
  , inside: function(x, y) {
      var box = this.bbox()

      return x > box.x
          && y > box.y
          && x < box.x + box.width
          && y < box.y + box.height
    }
    // Show element
  , show: function() {
      return this.style('display', '')
    }
    // Hide element
  , hide: function() {
      return this.style('display', 'none')
    }
    // Is element visible?
  , visible: function() {
      return this.style('display') != 'none'
    }
    // Return id on string conversion
  , toString: function() {
      return this.attr('id')
    }
    // Return array of classes on the node
  , classes: function() {
      var attr = this.attr('class')

      return attr == null ? [] : attr.trim().split(/\s+/)
    }
    // Return true if class exists on the node, false otherwise
  , hasClass: function(name) {
      return this.classes().indexOf(name) != -1
    }
    // Add class to the node
  , addClass: function(name) {
      if (!this.hasClass(name)) {
        var array = this.classes()
        array.push(name)
        this.attr('class', array.join(' '))
      }

      return this
    }
    // Remove class from the node
  , removeClass: function(name) {
      if (this.hasClass(name)) {
        this.attr('class', this.classes().filter(function(c) {
          return c != name
        }).join(' '))
      }

      return this
    }
    // Toggle the presence of a class on the node
  , toggleClass: function(name) {
      return this.hasClass(name) ? this.removeClass(name) : this.addClass(name)
    }
    // Get referenced element form attribute value
  , reference: function(attr) {
      return SVG.get(this.attr(attr))
    }
    // Returns the parent element instance
  , parent: function(type) {
      var parent = this

      // check for parent
      if(!parent.node.parentNode) return null

      // get parent element
      parent = SVG.adopt(parent.node.parentNode)

      if(!type) return parent

      // loop trough ancestors if type is given
      while(parent && parent.node instanceof SVGElement){
        if(typeof type === 'string' ? parent.matches(type) : parent instanceof type) return parent
        parent = SVG.adopt(parent.node.parentNode)
      }
    }
    // Get parent document
  , doc: function() {
      return this instanceof SVG.Doc ? this : this.parent(SVG.Doc)
    }
    // return array of all ancestors of given type up to the root svg
  , parents: function(type) {
      var parents = [], parent = this

      do{
        parent = parent.parent(type)
        if(!parent || !parent.node) break

        parents.push(parent)
      } while(parent.parent)

      return parents
    }
    // matches the element vs a css selector
  , matches: function(selector){
      return matches(this.node, selector)
    }
    // Returns the svg node to call native svg methods on it
  , native: function() {
      return this.node
    }
    // Import raw svg
  , svg: function(svg) {
      // create temporary holder
      var well = document.createElement('svg')

      // act as a setter if svg is given
      if (svg && this instanceof SVG.Parent) {
        // dump raw svg
        well.innerHTML = '<svg>' + svg.replace(/\n/, '').replace(/<(\w+)([^<]+?)\/>/g, '<$1$2></$1>') + '</svg>'

        // transplant nodes
        for (var i = 0, il = well.firstChild.childNodes.length; i < il; i++)
          this.node.appendChild(well.firstChild.firstChild)

      // otherwise act as a getter
      } else {
        // create a wrapping svg element in case of partial content
        well.appendChild(svg = document.createElement('svg'))

        // write svgjs data to the dom
        this.writeDataToDom()

        // insert a copy of this node
        svg.appendChild(this.node.cloneNode(true))

        // return target element
        return well.innerHTML.replace(/^<svg>/, '').replace(/<\/svg>$/, '')
      }

      return this
    }
  // write svgjs data to the dom
  , writeDataToDom: function() {

      // dump variables recursively
      if(this.each || this.lines){
        var fn = this.each ? this : this.lines();
        fn.each(function(){
          this.writeDataToDom()
        })
      }

      // remove previously set data
      this.node.removeAttribute('svgjs:data')

      if(Object.keys(this.dom).length)
        this.node.setAttribute('svgjs:data', JSON.stringify(this.dom)) // see #428

      return this
    }
  // set given data to the elements data property
  , setData: function(o){
      this.dom = o
      return this
    }
  , is: function(obj){
      return is(this, obj)
    }
  }
})

SVG.easing = {
  '-': function(pos){return pos}
, '<>':function(pos){return -Math.cos(pos * Math.PI) / 2 + 0.5}
, '>': function(pos){return  Math.sin(pos * Math.PI / 2)}
, '<': function(pos){return -Math.cos(pos * Math.PI / 2) + 1}
}

SVG.morph = function(pos){
  return function(from, to) {
    return new SVG.MorphObj(from, to).at(pos)
  }
}

SVG.Situation = SVG.invent({

  create: function(o){
    this.init = false
    this.reversed = false
    this.reversing = false

    this.duration = new SVG.Number(o.duration).valueOf()
    this.delay = new SVG.Number(o.delay).valueOf()

    this.start = +new Date() + this.delay
    this.finish = this.start + this.duration
    this.ease = o.ease

    // this.loop is incremented from 0 to this.loops
    // it is also incremented when in an infinite loop (when this.loops is true)
    this.loop = 0
    this.loops = false

    this.animations = {
      // functionToCall: [list of morphable objects]
      // e.g. move: [SVG.Number, SVG.Number]
    }

    this.attrs = {
      // holds all attributes which are not represented from a function svg.js provides
      // e.g. someAttr: SVG.Number
    }

    this.styles = {
      // holds all styles which should be animated
      // e.g. fill-color: SVG.Color
    }

    this.transforms = [
      // holds all transformations as transformation objects
      // e.g. [SVG.Rotate, SVG.Translate, SVG.Matrix]
    ]

    this.once = {
      // functions to fire at a specific position
      // e.g. "0.5": function foo(){}
    }

  }

})


SVG.FX = SVG.invent({

  create: function(element) {
    this._target = element
    this.situations = []
    this.active = false
    this.situation = null
    this.paused = false
    this.lastPos = 0
    this.pos = 0
    // The absolute position of an animation is its position in the context of its complete duration (including delay and loops)
    // When performing a delay, absPos is below 0 and when performing a loop, its value is above 1
    this.absPos = 0
    this._speed = 1
  }

, extend: {

    /**
     * sets or returns the target of this animation
     * @param o object || number In case of Object it holds all parameters. In case of number its the duration of the animation
     * @param ease function || string Function which should be used for easing or easing keyword
     * @param delay Number indicating the delay before the animation starts
     * @return target || this
     */
    animate: function(o, ease, delay){

      if(typeof o == 'object'){
        ease = o.ease
        delay = o.delay
        o = o.duration
      }

      var situation = new SVG.Situation({
        duration: o || 1000,
        delay: delay || 0,
        ease: SVG.easing[ease || '-'] || ease
      })

      this.queue(situation)

      return this
    }

    /**
     * sets a delay before the next element of the queue is called
     * @param delay Duration of delay in milliseconds
     * @return this.target()
     */
  , delay: function(delay){
      // The delay is performed by an empty situation with its duration
      // attribute set to the duration of the delay
      var situation = new SVG.Situation({
        duration: delay,
        delay: 0,
        ease: SVG.easing['-']
      })

      return this.queue(situation)
    }

    /**
     * sets or returns the target of this animation
     * @param null || target SVG.Element which should be set as new target
     * @return target || this
     */
  , target: function(target){
      if(target && target instanceof SVG.Element){
        this._target = target
        return this
      }

      return this._target
    }

    // returns the absolute position at a given time
  , timeToAbsPos: function(timestamp){
      return (timestamp - this.situation.start) / (this.situation.duration/this._speed)
    }

    // returns the timestamp from a given absolute positon
  , absPosToTime: function(absPos){
      return this.situation.duration/this._speed * absPos + this.situation.start
    }

    // starts the animationloop
  , startAnimFrame: function(){
      this.stopAnimFrame()
      this.animationFrame = requestAnimationFrame(function(){ this.step() }.bind(this))
    }

    // cancels the animationframe
  , stopAnimFrame: function(){
      cancelAnimationFrame(this.animationFrame)
    }

    // kicks off the animation - only does something when the queue is currently not active and at least one situation is set
  , start: function(){
      // dont start if already started
      if(!this.active && this.situation){
        this.active = true
        this.startCurrent()
      }

      return this
    }

    // start the current situation
  , startCurrent: function(){
      this.situation.start = +new Date + this.situation.delay/this._speed
      this.situation.finish = this.situation.start + this.situation.duration/this._speed
      return this.initAnimations().step()
    }

    /**
     * adds a function / Situation to the animation queue
     * @param fn function / situation to add
     * @return this
     */
  , queue: function(fn){
      if(typeof fn == 'function' || fn instanceof SVG.Situation)
        this.situations.push(fn)

      if(!this.situation) this.situation = this.situations.shift()

      return this
    }

    /**
     * pulls next element from the queue and execute it
     * @return this
     */
  , dequeue: function(){
      // stop current animation
      this.situation && this.situation.stop && this.situation.stop()

      // get next animation from queue
      this.situation = this.situations.shift()

      if(this.situation){
        if(this.situation instanceof SVG.Situation) {
          this.startCurrent()
        } else {
          // If it is not a SVG.Situation, then it is a function, we execute it
          this.situation.call(this)
        }
      }

      return this
    }

    // updates all animations to the current state of the element
    // this is important when one property could be changed from another property
  , initAnimations: function() {
      var i
      var s = this.situation

      if(s.init) return this

      for(i in s.animations){

        if(i == 'viewbox'){
          s.animations[i] = this.target().viewbox().morph(s.animations[i])
        }else{

          // TODO: this is not a clean clone of the array. We may have some unchecked references
          s.animations[i].value = (i == 'plot' ? this.target().array().value : this.target()[i]())

          // sometimes we get back an object and not the real value, fix this
          if(s.animations[i].value.value){
            s.animations[i].value = s.animations[i].value.value
          }

          if(s.animations[i].relative)
            s.animations[i].destination.value = s.animations[i].destination.value + s.animations[i].value

        }

      }

      for(i in s.attrs){
        if(s.attrs[i] instanceof SVG.Color){
          var color = new SVG.Color(this.target().attr(i))
          s.attrs[i].r = color.r
          s.attrs[i].g = color.g
          s.attrs[i].b = color.b
        }else{
          s.attrs[i].value = this.target().attr(i)// + s.attrs[i].value
        }
      }

      for(i in s.styles){
        s.styles[i].value = this.target().style(i)
      }

      s.initialTransformation = this.target().matrixify()

      s.init = true
      return this
    }
  , clearQueue: function(){
      this.situations = []
      return this
    }
  , clearCurrent: function(){
      this.situation = null
      return this
    }
    /** stops the animation immediately
     * @param jumpToEnd A Boolean indicating whether to complete the current animation immediately.
     * @param clearQueue A Boolean indicating whether to remove queued animation as well.
     * @return this
     */
  , stop: function(jumpToEnd, clearQueue){
      if(!this.active) this.start()

      if(clearQueue){
        this.clearQueue()
      }

      this.active = false

      if(jumpToEnd && this.situation){
        this.atEnd()
      }

      this.stopAnimFrame()

      return this.clearCurrent()
    }

    /** resets the element to the state where the current element has started
     * @return this
     */
  , reset: function(){
      if(this.situation){
        var temp = this.situation
        this.stop()
        this.situation = temp
        this.atStart()
      }
      return this
    }

    // Stop the currently-running animation, remove all queued animations, and complete all animations for the element.
  , finish: function(){

      this.stop(true, false)

      while(this.dequeue().situation && this.stop(true, false));

      this.clearQueue().clearCurrent()

      return this
    }

    // set the internal animation pointer at the start position, before any loops, and updates the visualisation
  , atStart: function() {
    return this.at(0, true)
  }

    // set the internal animation pointer at the end position, after all the loops, and updates the visualisation
  , atEnd: function() {
    if (this.situation.loops === true) {
      // If in a infinite loop, we end the current iteration
      return this.at(this.situation.loop+1, true)
    } else if(typeof this.situation.loops == 'number') {
      // If performing a finite number of loops, we go after all the loops
      return this.at(this.situation.loops, true)
    } else {
      // If no loops, we just go at the end
      return this.at(1, true)
    }
  }

    // set the internal animation pointer to the specified position and updates the visualisation
    // if isAbsPos is true, pos is treated as an absolute position
  , at: function(pos, isAbsPos){
      var durDivSpd = this.situation.duration/this._speed

      this.absPos = pos
      // If pos is not an absolute position, we convert it into one
      if (!isAbsPos) {
        if (this.situation.reversed) this.absPos = 1 - this.absPos
        this.absPos += this.situation.loop
      }

      this.situation.start = +new Date - this.absPos * durDivSpd
      this.situation.finish = this.situation.start + durDivSpd

      return this.step(true)
    }

    /**
     * sets or returns the speed of the animations
     * @param speed null || Number The new speed of the animations
     * @return Number || this
     */
  , speed: function(speed){
      if (speed === 0) return this.pause()

      if (speed) {
        this._speed = speed
        // We use an absolute position here so that speed can affect the delay before the animation
        return this.at(this.absPos, true)
      } else return this._speed
    }

    // Make loopable
  , loop: function(times, reverse) {
      var c = this.last()

      // store total loops
      c.loops = (times != null) ? times : true
      c.loop = 0

      if(reverse) c.reversing = true
      return this
    }

    // pauses the animation
  , pause: function(){
      this.paused = true
      this.stopAnimFrame()

      return this
    }

    // unpause the animation
  , play: function(){
      if(!this.paused) return this
      this.paused = false
      // We use an absolute position here so that the delay before the animation can be paused
      return this.at(this.absPos, true)
    }

    /**
     * toggle or set the direction of the animation
     * true sets direction to backwards while false sets it to forwards
     * @param reversed Boolean indicating whether to reverse the animation or not (default: toggle the reverse status)
     * @return this
     */
  , reverse: function(reversed){
      var c = this.last()

      if(typeof reversed == 'undefined') c.reversed = !c.reversed
      else c.reversed = reversed

      return this
    }


    /**
     * returns a float from 0-1 indicating the progress of the current animation
     * @param eased Boolean indicating whether the returned position should be eased or not
     * @return number
     */
  , progress: function(easeIt){
      return easeIt ? this.situation.ease(this.pos) : this.pos
    }

    /**
     * adds a callback function which is called when the current animation is finished
     * @param fn Function which should be executed as callback
     * @return number
     */
  , after: function(fn){
      var c = this.last()
        , wrapper = function wrapper(e){
            if(e.detail.situation == c){
              fn.call(this, c)
              this.off('finished.fx', wrapper) // prevent memory leak
            }
          }

      this.target().on('finished.fx', wrapper)
      return this
    }

    // adds a callback which is called whenever one animation step is performed
  , during: function(fn){
      var c = this.last()
        , wrapper = function(e){
            if(e.detail.situation == c){
              fn.call(this, e.detail.pos, SVG.morph(e.detail.pos), e.detail.eased, c)
            }
          }

      // see above
      this.target().off('during.fx', wrapper).on('during.fx', wrapper)

      return this.after(function(){
        this.off('during.fx', wrapper)
      })
    }

    // calls after ALL animations in the queue are finished
  , afterAll: function(fn){
      var wrapper = function wrapper(e){
            fn.call(this)
            this.off('allfinished.fx', wrapper)
          }

      // see above
      this.target().off('allfinished.fx', wrapper).on('allfinished.fx', wrapper)
      return this
    }

    // calls on every animation step for all animations
  , duringAll: function(fn){
      var wrapper = function(e){
            fn.call(this, e.detail.pos, SVG.morph(e.detail.pos), e.detail.eased, e.detail.situation)
          }

      this.target().off('during.fx', wrapper).on('during.fx', wrapper)

      return this.afterAll(function(){
        this.off('during.fx', wrapper)
      })
    }

  , last: function(){
      return this.situations.length ? this.situations[this.situations.length-1] : this.situation
    }

    // adds one property to the animations
  , add: function(method, args, type){
      this.last()[type || 'animations'][method] = args
      setTimeout(function(){this.start()}.bind(this), 0)
      return this
    }

    /** perform one step of the animation
     *  @param ignoreTime Boolean indicating whether to ignore time and use position directly or recalculate position based on time
     *  @return this
     */
  , step: function(ignoreTime){

      // convert current time to an absolute position
      if(!ignoreTime) this.absPos = this.timeToAbsPos(+new Date)

      // This part convert an absolute position to a position
      if(this.situation.loops !== false) {
        var absPos, absPosInt, lastLoop

        // If the absolute position is below 0, we just treat it as if it was 0
        absPos = Math.max(this.absPos, 0)
        absPosInt = Math.floor(absPos)

        if(this.situation.loops === true || absPosInt < this.situation.loops) {
          this.pos = absPos - absPosInt
          lastLoop = this.situation.loop
          this.situation.loop = absPosInt
        } else {
          this.absPos = this.situation.loops
          this.pos = 1
          // The -1 here is because we don't want to toggle reversed when all the loops have been completed
          lastLoop = this.situation.loop - 1
          this.situation.loop = this.situation.loops
        }

        if(this.situation.reversing) {
          // Toggle reversed if an odd number of loops as occured since the last call of step
          this.situation.reversed = this.situation.reversed != Boolean((this.situation.loop - lastLoop) % 2)
        }

      } else {
        // If there are no loop, the absolute position must not be above 1
        this.absPos = Math.min(this.absPos, 1)
        this.pos = this.absPos
      }

      // while the absolute position can be below 0, the position must not be below 0
      if(this.pos < 0) this.pos = 0

      if(this.situation.reversed) this.pos = 1 - this.pos


      // apply easing
      var eased = this.situation.ease(this.pos)

      // call once-callbacks
      for(var i in this.situation.once){
        if(i > this.lastPos && i <= eased){
          this.situation.once[i].call(this.target(), this.pos, eased)
          delete this.situation.once[i]
        }
      }

      // fire during callback with position, eased position and current situation as parameter
      if(this.active) this.target().fire('during', {pos: this.pos, eased: eased, fx: this, situation: this.situation})

      // the user may call stop or finish in the during callback
      // so make sure that we still have a valid situation
      if(!this.situation){
        return this
      }

      // apply the actual animation to every property
      this.eachAt()

      // do final code when situation is finished
      if((this.pos == 1 && !this.situation.reversed) || (this.situation.reversed && this.pos == 0)){

        // stop animation callback
        this.stopAnimFrame()

        // fire finished callback with current situation as parameter
        this.target().fire('finished', {fx:this, situation: this.situation})

        if(!this.situations.length){
          this.target().fire('allfinished')
          this.target().off('.fx') // there shouldnt be any binding left, but to make sure...
          this.active = false
        }

        // start next animation
        if(this.active) this.dequeue()
        else this.clearCurrent()

      }else if(!this.paused && this.active){
        // we continue animating when we are not at the end
        this.startAnimFrame()
      }

      // save last eased position for once callback triggering
      this.lastPos = eased
      return this

    }

    // calculates the step for every property and calls block with it
  , eachAt: function(){
      var i, at, self = this, target = this.target(), s = this.situation

      // apply animations which can be called trough a method
      for(i in s.animations){

        at = [].concat(s.animations[i]).map(function(el){
          return typeof el !== 'string' && el.at ? el.at(s.ease(self.pos), self.pos) : el
        })

        target[i].apply(target, at)

      }

      // apply animation which has to be applied with attr()
      for(i in s.attrs){

        at = [i].concat(s.attrs[i]).map(function(el){
          return typeof el !== 'string' && el.at ? el.at(s.ease(self.pos), self.pos) : el
        })

        target.attr.apply(target, at)

      }

      // apply animation which has to be applied with style()
      for(i in s.styles){

        at = [i].concat(s.styles[i]).map(function(el){
          return typeof el !== 'string' && el.at ? el.at(s.ease(self.pos), self.pos) : el
        })

        target.style.apply(target, at)

      }

      // animate initialTransformation which has to be chained
      if(s.transforms.length){

        // get initial initialTransformation
        at = s.initialTransformation
        for(i = 0, len = s.transforms.length; i < len; i++){

          // get next transformation in chain
          var a = s.transforms[i]

          // multiply matrix directly
          if(a instanceof SVG.Matrix){

            if(a.relative){
              at = at.multiply(new SVG.Matrix().morph(a).at(s.ease(this.pos)))
            }else{
              at = at.morph(a).at(s.ease(this.pos))
            }
            continue
          }

          // when transformation is absolute we have to reset the needed transformation first
          if(!a.relative)
            a.undo(at.extract())

          // and reapply it after
          at = at.multiply(a.at(s.ease(this.pos)))

        }

        // set new matrix on element
        target.matrix(at)
      }

      return this

    }


    // adds an once-callback which is called at a specific position and never again
  , once: function(pos, fn, isEased){

      if(!isEased)pos = this.situation.ease(pos)

      this.situation.once[pos] = fn

      return this
    }

  }

, parent: SVG.Element

  // Add method to parent elements
, construct: {
    // Get fx module or create a new one, then animate with given duration and ease
    animate: function(o, ease, delay) {
      return (this.fx || (this.fx = new SVG.FX(this))).animate(o, ease, delay)
    }
  , delay: function(delay){
      return (this.fx || (this.fx = new SVG.FX(this))).delay(delay)
    }
  , stop: function(jumpToEnd, clearQueue) {
      if (this.fx)
        this.fx.stop(jumpToEnd, clearQueue)

      return this
    }
  , finish: function() {
      if (this.fx)
        this.fx.finish()

      return this
    }
    // Pause current animation
  , pause: function() {
      if (this.fx)
        this.fx.pause()

      return this
    }
    // Play paused current animation
  , play: function() {
      if (this.fx)
        this.fx.play()

      return this
    }
    // Set/Get the speed of the animations
  , speed: function(speed) {
      if (this.fx)
        if (speed == null)
          return this.fx.speed()
        else
          this.fx.speed(speed)

      return this
    }
  }

})

// MorphObj is used whenever no morphable object is given
SVG.MorphObj = SVG.invent({

  create: function(from, to){
    // prepare color for morphing
    if(SVG.Color.isColor(to)) return new SVG.Color(from).morph(to)
    // prepare number for morphing
    if(SVG.regex.numberAndUnit.test(to)) return new SVG.Number(from).morph(to)

    // prepare for plain morphing
    this.value = 0
    this.destination = to
  }

, extend: {
    at: function(pos, real){
      return real < 1 ? this.value : this.destination
    },

    valueOf: function(){
      return this.value
    }
  }

})

SVG.extend(SVG.FX, {
  // Add animatable attributes
  attr: function(a, v, relative) {
    // apply attributes individually
    if (typeof a == 'object') {
      for (var key in a)
        this.attr(key, a[key])

    } else {
      // the MorphObj takes care about the right function used
      this.add(a, new SVG.MorphObj(null, v), 'attrs')
    }

    return this
  }
  // Add animatable styles
, style: function(s, v) {
    if (typeof s == 'object')
      for (var key in s)
        this.style(key, s[key])

    else
      this.add(s, new SVG.MorphObj(null, v), 'styles')

    return this
  }
  // Animatable x-axis
, x: function(x, relative) {
    if(this.target() instanceof SVG.G){
      this.transform({x:x}, relative)
      return this
    }

    var num = new SVG.Number().morph(x)
    num.relative = relative
    return this.add('x', num)
  }
  // Animatable y-axis
, y: function(y, relative) {
    if(this.target() instanceof SVG.G){
      this.transform({y:y}, relative)
      return this
    }

    var num = new SVG.Number().morph(y)
    num.relative = relative
    return this.add('y', num)
  }
  // Animatable center x-axis
, cx: function(x) {
    return this.add('cx', new SVG.Number().morph(x))
  }
  // Animatable center y-axis
, cy: function(y) {
    return this.add('cy', new SVG.Number().morph(y))
  }
  // Add animatable move
, move: function(x, y) {
    return this.x(x).y(y)
  }
  // Add animatable center
, center: function(x, y) {
    return this.cx(x).cy(y)
  }
  // Add animatable size
, size: function(width, height) {
    if (this.target() instanceof SVG.Text) {
      // animate font size for Text elements
      this.attr('font-size', width)

    } else {
      // animate bbox based size for all other elements
      var box

      if(!width || !height){
        box = this.target().bbox()
      }

      if(!width){
        width = box.width / box.height  * height
      }

      if(!height){
        height = box.height / box.width  * width
      }

      this.add('width' , new SVG.Number().morph(width))
          .add('height', new SVG.Number().morph(height))

    }

    return this
  }
  // Add animatable plot
, plot: function(p) {
    return this.add('plot', this.target().array().morph(p))
  }
  // Add leading method
, leading: function(value) {
    return this.target().leading ?
      this.add('leading', new SVG.Number().morph(value)) :
      this
  }
  // Add animatable viewbox
, viewbox: function(x, y, width, height) {
    if (this.target() instanceof SVG.Container) {
      this.add('viewbox', new SVG.ViewBox(x, y, width, height))
    }

    return this
  }
, update: function(o) {
    if (this.target() instanceof SVG.Stop) {
      if (typeof o == 'number' || o instanceof SVG.Number) {
        return this.update({
          offset:  arguments[0]
        , color:   arguments[1]
        , opacity: arguments[2]
        })
      }

      if (o.opacity != null) this.attr('stop-opacity', o.opacity)
      if (o.color   != null) this.attr('stop-color', o.color)
      if (o.offset  != null) this.attr('offset', o.offset)
    }

    return this
  }
})

SVG.BBox = SVG.invent({
  // Initialize
  create: function(element) {
    // get values if element is given
    if (element) {
      var box

      // yes this is ugly, but Firefox can be a bitch when it comes to elements that are not yet rendered
      try {

        // the element is NOT in the dom, throw error
        if(!document.documentElement.contains(element.node)) throw new Exception('Element not in the dom')

        // find native bbox
        box = element.node.getBBox()
      } catch(e) {
        if(element instanceof SVG.Shape){
          var clone = element.clone(SVG.parser.draw).show()
          box = clone.bbox()
          clone.remove()
        }else{
          box = {
            x:      element.node.clientLeft
          , y:      element.node.clientTop
          , width:  element.node.clientWidth
          , height: element.node.clientHeight
          }
        }
      }

      // plain x and y
      this.x = box.x
      this.y = box.y

      // plain width and height
      this.width  = box.width
      this.height = box.height
    }

    // add center, right and bottom
    fullBox(this)
  }

  // Define Parent
, parent: SVG.Element

  // Constructor
, construct: {
    // Get bounding box
    bbox: function() {
      return new SVG.BBox(this)
    }
  }

})

SVG.TBox = SVG.invent({
  // Initialize
  create: function(element) {
    // get values if element is given
    if (element) {
      var t   = element.ctm().extract()
        , box = element.bbox()

      // width and height including transformations
      this.width  = box.width  * t.scaleX
      this.height = box.height * t.scaleY

      // x and y including transformations
      this.x = box.x + t.x
      this.y = box.y + t.y
    }

    // add center, right and bottom
    fullBox(this)
  }

  // Define Parent
, parent: SVG.Element

  // Constructor
, construct: {
    // Get transformed bounding box
    tbox: function() {
      return new SVG.TBox(this)
    }
  }

})


SVG.RBox = SVG.invent({
  // Initialize
  create: function(element) {
    if (element) {
      var e    = element.doc().parent()
        , box  = element.node.getBoundingClientRect()
        , zoom = 1

      // get screen offset
      this.x = box.left
      this.y = box.top

      // subtract parent offset
      this.x -= e.offsetLeft
      this.y -= e.offsetTop

      while (e = e.offsetParent) {
        this.x -= e.offsetLeft
        this.y -= e.offsetTop
      }

      // calculate cumulative zoom from svg documents
      e = element
      while (e.parent && (e = e.parent())) {
        if (e.viewbox) {
          zoom *= e.viewbox().zoom
          this.x -= e.x() || 0
          this.y -= e.y() || 0
        }
      }

      // recalculate viewbox distortion
      this.width  = box.width  /= zoom
      this.height = box.height /= zoom
    }

    // add center, right and bottom
    fullBox(this)

    // offset by window scroll position, because getBoundingClientRect changes when window is scrolled
    this.x += window.pageXOffset
    this.y += window.pageYOffset
  }

  // define Parent
, parent: SVG.Element

  // Constructor
, construct: {
    // Get rect box
    rbox: function() {
      return new SVG.RBox(this)
    }
  }

})

// Add universal merge method
;[SVG.BBox, SVG.TBox, SVG.RBox].forEach(function(c) {

  SVG.extend(c, {
    // Merge rect box with another, return a new instance
    merge: function(box) {
      var b = new c()

      // merge boxes
      b.x      = Math.min(this.x, box.x)
      b.y      = Math.min(this.y, box.y)
      b.width  = Math.max(this.x + this.width,  box.x + box.width)  - b.x
      b.height = Math.max(this.y + this.height, box.y + box.height) - b.y

      return fullBox(b)
    }

  })

})

SVG.Matrix = SVG.invent({
  // Initialize
  create: function(source) {
    var i, base = arrayToMatrix([1, 0, 0, 1, 0, 0])

    // ensure source as object
    source = source instanceof SVG.Element ?
      source.matrixify() :
    typeof source === 'string' ?
      stringToMatrix(source) :
    arguments.length == 6 ?
      arrayToMatrix([].slice.call(arguments)) :
    typeof source === 'object' ?
      source : base

    // merge source
    for (i = abcdef.length - 1; i >= 0; --i)
      this[abcdef[i]] = source && typeof source[abcdef[i]] === 'number' ?
        source[abcdef[i]] : base[abcdef[i]]
  }

  // Add methods
, extend: {
    // Extract individual transformations
    extract: function() {
      // find delta transform points
      var px    = deltaTransformPoint(this, 0, 1)
        , py    = deltaTransformPoint(this, 1, 0)
        , skewX = 180 / Math.PI * Math.atan2(px.y, px.x) - 90

      return {
        // translation
        x:        this.e
      , y:        this.f
      , transformedX:(this.e * Math.cos(skewX * Math.PI / 180) + this.f * Math.sin(skewX * Math.PI / 180)) / Math.sqrt(this.a * this.a + this.b * this.b)
      , transformedY:(this.f * Math.cos(skewX * Math.PI / 180) + this.e * Math.sin(-skewX * Math.PI / 180)) / Math.sqrt(this.c * this.c + this.d * this.d)
        // skew
      , skewX:    -skewX
      , skewY:    180 / Math.PI * Math.atan2(py.y, py.x)
        // scale
      , scaleX:   Math.sqrt(this.a * this.a + this.b * this.b)
      , scaleY:   Math.sqrt(this.c * this.c + this.d * this.d)
        // rotation
      , rotation: skewX
      , a: this.a
      , b: this.b
      , c: this.c
      , d: this.d
      , e: this.e
      , f: this.f
      , matrix: new SVG.Matrix(this)
      }
    }
    // Clone matrix
  , clone: function() {
      return new SVG.Matrix(this)
    }
    // Morph one matrix into another
  , morph: function(matrix) {
      // store new destination
      this.destination = new SVG.Matrix(matrix)

      return this
    }
    // Get morphed matrix at a given position
  , at: function(pos) {
      // make sure a destination is defined
      if (!this.destination) return this

      // calculate morphed matrix at a given position
      var matrix = new SVG.Matrix({
        a: this.a + (this.destination.a - this.a) * pos
      , b: this.b + (this.destination.b - this.b) * pos
      , c: this.c + (this.destination.c - this.c) * pos
      , d: this.d + (this.destination.d - this.d) * pos
      , e: this.e + (this.destination.e - this.e) * pos
      , f: this.f + (this.destination.f - this.f) * pos
      })

      // process parametric rotation if present
      if (this.param && this.param.to) {
        // calculate current parametric position
        var param = {
          rotation: this.param.from.rotation + (this.param.to.rotation - this.param.from.rotation) * pos
        , cx:       this.param.from.cx
        , cy:       this.param.from.cy
        }

        // rotate matrix
        matrix = matrix.rotate(
          (this.param.to.rotation - this.param.from.rotation * 2) * pos
        , param.cx
        , param.cy
        )

        // store current parametric values
        matrix.param = param
      }

      return matrix
    }
    // Multiplies by given matrix
  , multiply: function(matrix) {
      return new SVG.Matrix(this.native().multiply(parseMatrix(matrix).native()))
    }
    // Inverses matrix
  , inverse: function() {
      return new SVG.Matrix(this.native().inverse())
    }
    // Translate matrix
  , translate: function(x, y) {
      return new SVG.Matrix(this.native().translate(x || 0, y || 0))
    }
    // Scale matrix
  , scale: function(x, y, cx, cy) {
      // support uniformal scale
      if (arguments.length == 1) {
        y = x
      } else if (arguments.length == 3) {
        cy = cx
        cx = y
        y = x
      }

      return this.around(cx, cy, new SVG.Matrix(x, 0, 0, y, 0, 0))
    }
    // Rotate matrix
  , rotate: function(r, cx, cy) {
      // convert degrees to radians
      r = SVG.utils.radians(r)

      return this.around(cx, cy, new SVG.Matrix(Math.cos(r), Math.sin(r), -Math.sin(r), Math.cos(r), 0, 0))
    }
    // Flip matrix on x or y, at a given offset
  , flip: function(a, o) {
      return a == 'x' ? this.scale(-1, 1, o, 0) : this.scale(1, -1, 0, o)
    }
    // Skew
  , skew: function(x, y, cx, cy) {
      // support uniformal skew
      if (arguments.length == 1) {
        y = x
      } else if (arguments.length == 3) {
        cy = cx
        cx = y
        y = x
      }

      // convert degrees to radians
      x = SVG.utils.radians(x)
      y = SVG.utils.radians(y)

      return this.around(cx, cy, new SVG.Matrix(1, Math.tan(y), Math.tan(x), 1, 0, 0))
    }
    // SkewX
  , skewX: function(x, cx, cy) {
      return this.skew(x, 0, cx, cy)
    }
    // SkewY
  , skewY: function(y, cx, cy) {
      return this.skew(0, y, cx, cy)
    }
    // Transform around a center point
  , around: function(cx, cy, matrix) {
      return this
        .multiply(new SVG.Matrix(1, 0, 0, 1, cx || 0, cy || 0))
        .multiply(matrix)
        .multiply(new SVG.Matrix(1, 0, 0, 1, -cx || 0, -cy || 0))
    }
    // Convert to native SVGMatrix
  , native: function() {
      // create new matrix
      var matrix = SVG.parser.native.createSVGMatrix()

      // update with current values
      for (var i = abcdef.length - 1; i >= 0; i--)
        matrix[abcdef[i]] = this[abcdef[i]]

      return matrix
    }
    // Convert matrix to string
  , toString: function() {
      return 'matrix(' + this.a + ',' + this.b + ',' + this.c + ',' + this.d + ',' + this.e + ',' + this.f + ')'
    }
  }

  // Define parent
, parent: SVG.Element

  // Add parent method
, construct: {
    // Get current matrix
    ctm: function() {
      return new SVG.Matrix(this.node.getCTM())
    },
    // Get current screen matrix
    screenCTM: function() {
      return new SVG.Matrix(this.node.getScreenCTM())
    }

  }

})

SVG.Point = SVG.invent({
  // Initialize
  create: function(x,y) {
    var i, source, base = {x:0, y:0}

    // ensure source as object
    source = Array.isArray(x) ?
      {x:x[0], y:x[1]} :
    typeof x === 'object' ?
      {x:x.x, y:x.y} :
    x != null ?
      {x:x, y:(y != null ? y : x)} : base // If y has no value, then x is used has its value

    // merge source
    this.x = source.x
    this.y = source.y
  }

  // Add methods
, extend: {
    // Clone point
    clone: function() {
      return new SVG.Point(this)
    }
    // Morph one point into another
  , morph: function(x, y) {
      // store new destination
      this.destination = new SVG.Point(x, y)

      return this
    }
    // Get morphed point at a given position
  , at: function(pos) {
      // make sure a destination is defined
      if (!this.destination) return this

      // calculate morphed matrix at a given position
      var point = new SVG.Point({
        x: this.x + (this.destination.x - this.x) * pos
      , y: this.y + (this.destination.y - this.y) * pos
      })

      return point
    }
    // Convert to native SVGPoint
  , native: function() {
      // create new point
      var point = SVG.parser.native.createSVGPoint()

      // update with current values
      point.x = this.x
      point.y = this.y

      return point
    }
    // transform point with matrix
  , transform: function(matrix) {
      return new SVG.Point(this.native().matrixTransform(matrix.native()))
    }

  }

})

SVG.extend(SVG.Element, {

  // Get point
  point: function(x, y) {
    return new SVG.Point(x,y).transform(this.screenCTM().inverse());
  }

})

SVG.extend(SVG.Element, {
  // Set svg element attribute
  attr: function(a, v, n) {
    // act as full getter
    if (a == null) {
      // get an object of attributes
      a = {}
      v = this.node.attributes
      for (n = v.length - 1; n >= 0; n--)
        a[v[n].nodeName] = SVG.regex.isNumber.test(v[n].nodeValue) ? parseFloat(v[n].nodeValue) : v[n].nodeValue

      return a

    } else if (typeof a == 'object') {
      // apply every attribute individually if an object is passed
      for (v in a) this.attr(v, a[v])

    } else if (v === null) {
        // remove value
        this.node.removeAttribute(a)

    } else if (v == null) {
      // act as a getter if the first and only argument is not an object
      v = this.node.getAttribute(a)
      return v == null ?
        SVG.defaults.attrs[a] :
      SVG.regex.isNumber.test(v) ?
        parseFloat(v) : v

    } else {
      // BUG FIX: some browsers will render a stroke if a color is given even though stroke width is 0
      if (a == 'stroke-width')
        this.attr('stroke', parseFloat(v) > 0 ? this._stroke : null)
      else if (a == 'stroke')
        this._stroke = v

      // convert image fill and stroke to patterns
      if (a == 'fill' || a == 'stroke') {
        if (SVG.regex.isImage.test(v))
          v = this.doc().defs().image(v, 0, 0)

        if (v instanceof SVG.Image)
          v = this.doc().defs().pattern(0, 0, function() {
            this.add(v)
          })
      }

      // ensure correct numeric values (also accepts NaN and Infinity)
      if (typeof v === 'number')
        v = new SVG.Number(v)

      // ensure full hex color
      else if (SVG.Color.isColor(v))
        v = new SVG.Color(v)

      // parse array values
      else if (Array.isArray(v))
        v = new SVG.Array(v)

      // store parametric transformation values locally
      else if (v instanceof SVG.Matrix && v.param)
        this.param = v.param

      // if the passed attribute is leading...
      if (a == 'leading') {
        // ... call the leading method instead
        if (this.leading)
          this.leading(v)
      } else {
        // set given attribute on node
        typeof n === 'string' ?
          this.node.setAttributeNS(n, a, v.toString()) :
          this.node.setAttribute(a, v.toString())
      }

      // rebuild if required
      if (this.rebuild && (a == 'font-size' || a == 'x'))
        this.rebuild(a, v)
    }

    return this
  }
})
SVG.extend(SVG.Element, {
  // Add transformations
  transform: function(o, relative) {
    // get target in case of the fx module, otherwise reference this
    var target = this
      , matrix

    // act as a getter
    if (typeof o !== 'object') {
      // get current matrix
      matrix = new SVG.Matrix(target).extract()

      return typeof o === 'string' ? matrix[o] : matrix
    }

    // get current matrix
    matrix = new SVG.Matrix(target)

    // ensure relative flag
    relative = !!relative || !!o.relative

    // act on matrix
    if (o.a != null) {
      matrix = relative ?
        // relative
        matrix.multiply(new SVG.Matrix(o)) :
        // absolute
        new SVG.Matrix(o)

    // act on rotation
    } else if (o.rotation != null) {
      // ensure centre point
      ensureCentre(o, target)

      // apply transformation
      matrix = relative ?
        // relative
        matrix.rotate(o.rotation, o.cx, o.cy) :
        // absolute
        matrix.rotate(o.rotation - matrix.extract().rotation, o.cx, o.cy)

    // act on scale
    } else if (o.scale != null || o.scaleX != null || o.scaleY != null) {
      // ensure centre point
      ensureCentre(o, target)

      // ensure scale values on both axes
      o.scaleX = o.scale != null ? o.scale : o.scaleX != null ? o.scaleX : 1
      o.scaleY = o.scale != null ? o.scale : o.scaleY != null ? o.scaleY : 1

      if (!relative) {
        // absolute; multiply inversed values
        var e = matrix.extract()
        o.scaleX = o.scaleX * 1 / e.scaleX
        o.scaleY = o.scaleY * 1 / e.scaleY
      }

      matrix = matrix.scale(o.scaleX, o.scaleY, o.cx, o.cy)

    // act on skew
    } else if (o.skew != null || o.skewX != null || o.skewY != null) {
      // ensure centre point
      ensureCentre(o, target)

      // ensure skew values on both axes
      o.skewX = o.skew != null ? o.skew : o.skewX != null ? o.skewX : 0
      o.skewY = o.skew != null ? o.skew : o.skewY != null ? o.skewY : 0

      if (!relative) {
        // absolute; reset skew values
        var e = matrix.extract()
        matrix = matrix.multiply(new SVG.Matrix().skew(e.skewX, e.skewY, o.cx, o.cy).inverse())
      }

      matrix = matrix.skew(o.skewX, o.skewY, o.cx, o.cy)

    // act on flip
    } else if (o.flip) {
      matrix = matrix.flip(
        o.flip
      , o.offset == null ? target.bbox()['c' + o.flip] : o.offset
      )

    // act on translate
    } else if (o.x != null || o.y != null) {
      if (relative) {
        // relative
        matrix = matrix.translate(o.x, o.y)
      } else {
        // absolute
        if (o.x != null) matrix.e = o.x
        if (o.y != null) matrix.f = o.y
      }
    }

    return this.attr('transform', matrix)
  }
})

SVG.extend(SVG.FX, {
  transform: function(o, relative) {
    // get target in case of the fx module, otherwise reference this
    var target = this.target()
      , matrix

    // act as a getter
    if (typeof o !== 'object') {
      // get current matrix
      matrix = new SVG.Matrix(target).extract()

      return typeof o === 'string' ? matrix[o] : matrix
    }

    // ensure relative flag
    relative = !!relative || !!o.relative

    // act on matrix
    if (o.a != null) {
      matrix = new SVG.Matrix(o)

    // act on rotation
    } else if (o.rotation != null) {
      // ensure centre point
      ensureCentre(o, target)

      // apply transformation
      matrix = new SVG.Rotate(o.rotation, o.cx, o.cy)

    // act on scale
    } else if (o.scale != null || o.scaleX != null || o.scaleY != null) {
      // ensure centre point
      ensureCentre(o, target)

      // ensure scale values on both axes
      o.scaleX = o.scale != null ? o.scale : o.scaleX != null ? o.scaleX : 1
      o.scaleY = o.scale != null ? o.scale : o.scaleY != null ? o.scaleY : 1

      matrix = new SVG.Scale(o.scaleX, o.scaleY, o.cx, o.cy)

    // act on skew
    } else if (o.skewX != null || o.skewY != null) {
      // ensure centre point
      ensureCentre(o, target)

      // ensure skew values on both axes
      o.skewX = o.skewX != null ? o.skewX : 0
      o.skewY = o.skewY != null ? o.skewY : 0

      matrix = new SVG.Skew(o.skewX, o.skewY, o.cx, o.cy)

    // act on flip
    } else if (o.flip) {
      matrix = new SVG.Matrix().morph(new SVG.Matrix().flip(
        o.flip
      , o.offset == null ? target.bbox()['c' + o.flip] : o.offset
      ))

    // act on translate
    } else if (o.x != null || o.y != null) {
      matrix = new SVG.Translate(o.x, o.y)
    }

    if(!matrix) return this

    matrix.relative = relative

    this.last().transforms.push(matrix)

    setTimeout(function(){this.start()}.bind(this), 0)

    return this
  }
})

SVG.extend(SVG.Element, {
  // Reset all transformations
  untransform: function() {
    return this.attr('transform', null)
  },
  // merge the whole transformation chain into one matrix and returns it
  matrixify: function() {

    var matrix = (this.attr('transform') || '')
      // split transformations
      .split(/\)\s*,?\s*/).slice(0,-1).map(function(str){
        // generate key => value pairs
        var kv = str.trim().split('(')
        return [kv[0], kv[1].split(SVG.regex.matrixElements).map(function(str){ return parseFloat(str) })]
      })
      // calculate every transformation into one matrix
      .reduce(function(matrix, transform){

        if(transform[0] == 'matrix') return matrix.multiply(arrayToMatrix(transform[1]))
        return matrix[transform[0]].apply(matrix, transform[1])

      }, new SVG.Matrix())

    return matrix
  },
  // add an element to another parent without changing the visual representation on the screen
  toParent: function(parent) {
    if(this == parent) return this
    var ctm = this.screenCTM()
    var temp = parent.rect(1,1)
    var pCtm = temp.screenCTM().inverse()
    temp.remove()

    this.addTo(parent).untransform().transform(pCtm.multiply(ctm))

    return this
  },
  // same as above with parent equals root-svg
  toDoc: function() {
    return this.toParent(this.doc())
  }

})

SVG.Transformation = SVG.invent({

  create: function(source, inversed){

    if(arguments.length > 1 && typeof inversed != 'boolean'){
      return this.create([].slice.call(arguments))
    }

    if(typeof source == 'object'){
      for(var i = 0, len = this.arguments.length; i < len; ++i){
        this[this.arguments[i]] = source[this.arguments[i]]
      }
    }

    if(Array.isArray(source)){
      for(var i = 0, len = this.arguments.length; i < len; ++i){
        this[this.arguments[i]] = source[i]
      }
    }

    this.inversed = false

    if(inversed === true){
      this.inversed = true
    }

  }

, extend: {

    at: function(pos){

      var params = []

      for(var i = 0, len = this.arguments.length; i < len; ++i){
        params.push(this[this.arguments[i]])
      }

      var m = this._undo || new SVG.Matrix()

      m = new SVG.Matrix().morph(SVG.Matrix.prototype[this.method].apply(m, params)).at(pos)

      return this.inversed ? m.inverse() : m

    }

  , undo: function(o){
      for(var i = 0, len = this.arguments.length; i < len; ++i){
        o[this.arguments[i]] = typeof this[this.arguments[i]] == 'undefined' ? 0 : o[this.arguments[i]]
      }

      // The method SVG.Matrix.extract which was used before calling this
      // method to obtain a value for the parameter o doesn't return a cx and
      // a cy so we use the ones that were provided to this object at its creation
      o.cx = this.cx
      o.cy = this.cy

      this._undo = new SVG[capitalize(this.method)](o, true).at(1)

      return this
    }

  }

})

SVG.Translate = SVG.invent({

  parent: SVG.Matrix
, inherit: SVG.Transformation

, create: function(source, inversed){
    if(typeof source == 'object') this.constructor.call(this, source, inversed)
    else this.constructor.call(this, [].slice.call(arguments))
  }

, extend: {
    arguments: ['transformedX', 'transformedY']
  , method: 'translate'
  }

})

SVG.Rotate = SVG.invent({

  parent: SVG.Matrix
, inherit: SVG.Transformation

, create: function(source, inversed){
    if(typeof source == 'object') this.constructor.call(this, source, inversed)
    else this.constructor.call(this, [].slice.call(arguments))
  }

, extend: {
    arguments: ['rotation', 'cx', 'cy']
  , method: 'rotate'
  , at: function(pos){
      var m = new SVG.Matrix().rotate(new SVG.Number().morph(this.rotation - (this._undo ? this._undo.rotation : 0)).at(pos), this.cx, this.cy)
      return this.inversed ? m.inverse() : m
    }
  , undo: function(o){
      this._undo = o
    }
  }

})

SVG.Scale = SVG.invent({

  parent: SVG.Matrix
, inherit: SVG.Transformation

, create: function(source, inversed){
    if(typeof source == 'object') this.constructor.call(this, source, inversed)
    else this.constructor.call(this, [].slice.call(arguments))
  }

, extend: {
    arguments: ['scaleX', 'scaleY', 'cx', 'cy']
  , method: 'scale'
  }

})

SVG.Skew = SVG.invent({

  parent: SVG.Matrix
, inherit: SVG.Transformation

, create: function(source, inversed){
    if(typeof source == 'object') this.constructor.call(this, source, inversed)
    else this.constructor.call(this, [].slice.call(arguments))
  }

, extend: {
    arguments: ['skewX', 'skewY', 'cx', 'cy']
  , method: 'skew'
  }

})

SVG.extend(SVG.Element, {
  // Dynamic style generator
  style: function(s, v) {
    if (arguments.length == 0) {
      // get full style
      return this.node.style.cssText || ''

    } else if (arguments.length < 2) {
      // apply every style individually if an object is passed
      if (typeof s == 'object') {
        for (v in s) this.style(v, s[v])

      } else if (SVG.regex.isCss.test(s)) {
        // parse css string
        s = s.split(';')

        // apply every definition individually
        for (var i = 0; i < s.length; i++) {
          v = s[i].split(':')
          this.style(v[0].replace(/\s+/g, ''), v[1])
        }
      } else {
        // act as a getter if the first and only argument is not an object
        return this.node.style[camelCase(s)]
      }

    } else {
      this.node.style[camelCase(s)] = v === null || SVG.regex.isBlank.test(v) ? '' : v
    }

    return this
  }
})
SVG.Parent = SVG.invent({
  // Initialize node
  create: function(element) {
    this.constructor.call(this, element)
  }

  // Inherit from
, inherit: SVG.Element

  // Add class methods
, extend: {
    // Returns all child elements
    children: function() {
      return SVG.utils.map(SVG.utils.filterSVGElements(this.node.childNodes), function(node) {
        return SVG.adopt(node)
      })
    }
    // Add given element at a position
  , add: function(element, i) {
      if (i == null)
        this.node.appendChild(element.node)
      else if (element.node != this.node.childNodes[i])
        this.node.insertBefore(element.node, this.node.childNodes[i])

      return this
    }
    // Basically does the same as `add()` but returns the added element instead
  , put: function(element, i) {
      this.add(element, i)
      return element
    }
    // Checks if the given element is a child
  , has: function(element) {
      return this.index(element) >= 0
    }
    // Gets index of given element
  , index: function(element) {
      return [].slice.call(this.node.childNodes).indexOf(element.node)
    }
    // Get a element at the given index
  , get: function(i) {
      return SVG.adopt(this.node.childNodes[i])
    }
    // Get first child
  , first: function() {
      return this.get(0)
    }
    // Get the last child
  , last: function() {
      return this.get(this.node.childNodes.length - 1)
    }
    // Iterates over all children and invokes a given block
  , each: function(block, deep) {
      var i, il
        , children = this.children()

      for (i = 0, il = children.length; i < il; i++) {
        if (children[i] instanceof SVG.Element)
          block.apply(children[i], [i, children])

        if (deep && (children[i] instanceof SVG.Container))
          children[i].each(block, deep)
      }

      return this
    }
    // Remove a given child
  , removeElement: function(element) {
      this.node.removeChild(element.node)

      return this
    }
    // Remove all elements in this container
  , clear: function() {
      // remove children
      while(this.node.hasChildNodes())
        this.node.removeChild(this.node.lastChild)

      // remove defs reference
      delete this._defs

      return this
    }
  , // Get defs
    defs: function() {
      return this.doc().defs()
    }
  }

})

SVG.extend(SVG.Parent, {

  ungroup: function(parent, depth) {
    if(depth === 0 || this instanceof SVG.Defs) return this

    parent = parent || (this instanceof SVG.Doc ? this : this.parent(SVG.Parent))
    depth = depth || Infinity

    this.each(function(){
      if(this instanceof SVG.Defs) return this
      if(this instanceof SVG.Parent) return this.ungroup(parent, depth-1)
      return this.toParent(parent)
    })

    this.node.firstChild || this.remove()

    return this
  },

  flatten: function(parent, depth) {
    return this.ungroup(parent, depth)
  }

})
SVG.Container = SVG.invent({
  // Initialize node
  create: function(element) {
    this.constructor.call(this, element)
  }

  // Inherit from
, inherit: SVG.Parent

})

SVG.ViewBox = SVG.invent({

  create: function(source) {
    var i, base = [0, 0, 0, 0]

    var x, y, width, height, box, view, we, he
      , wm   = 1 // width multiplier
      , hm   = 1 // height multiplier
      , reg  = /[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:e[+-]?\d+)?/gi

    if(source instanceof SVG.Element){

      we = source
      he = source
      view = (source.attr('viewBox') || '').match(reg)
      box = source.bbox

      // get dimensions of current node
      width  = new SVG.Number(source.width())
      height = new SVG.Number(source.height())

      // find nearest non-percentual dimensions
      while (width.unit == '%') {
        wm *= width.value
        width = new SVG.Number(we instanceof SVG.Doc ? we.parent().offsetWidth : we.parent().width())
        we = we.parent()
      }
      while (height.unit == '%') {
        hm *= height.value
        height = new SVG.Number(he instanceof SVG.Doc ? he.parent().offsetHeight : he.parent().height())
        he = he.parent()
      }

      // ensure defaults
      this.x      = 0
      this.y      = 0
      this.width  = width  * wm
      this.height = height * hm
      this.zoom   = 1

      if (view) {
        // get width and height from viewbox
        x      = parseFloat(view[0])
        y      = parseFloat(view[1])
        width  = parseFloat(view[2])
        height = parseFloat(view[3])

        // calculate zoom accoring to viewbox
        this.zoom = ((this.width / this.height) > (width / height)) ?
          this.height / height :
          this.width  / width

        // calculate real pixel dimensions on parent SVG.Doc element
        this.x      = x
        this.y      = y
        this.width  = width
        this.height = height

      }

    }else{

      // ensure source as object
      source = typeof source === 'string' ?
        source.match(reg).map(function(el){ return parseFloat(el) }) :
      Array.isArray(source) ?
        source :
      typeof source == 'object' ?
        [source.x, source.y, source.width, source.height] :
      arguments.length == 4 ?
        [].slice.call(arguments) :
        base

      this.x = source[0]
      this.y = source[1]
      this.width = source[2]
      this.height = source[3]
    }


  }

, extend: {

    toString: function() {
      return this.x + ' ' + this.y + ' ' + this.width + ' ' + this.height
    }
  , morph: function(v){

      var v = arguments.length == 1 ?
        [v.x, v.y, v.width, v.height] :
        [].slice.call(arguments)

      this.destination = new SVG.ViewBox(v)

      return this

    }

  , at: function(pos) {

    if(!this.destination) return this

    return new SVG.ViewBox([
        this.x + (this.destination.x - this.x) * pos
      , this.y + (this.destination.y - this.y) * pos
      , this.width + (this.destination.width - this.width) * pos
      , this.height + (this.destination.height - this.height) * pos
    ])

    }

  }

  // Define parent
, parent: SVG.Container

  // Add parent method
, construct: {

    // get/set viewbox
    viewbox: function(v) {
      if (arguments.length == 0)
        // act as a getter if there are no arguments
        return new SVG.ViewBox(this)

      // otherwise act as a setter
      v = arguments.length == 1 ?
        [v.x, v.y, v.width, v.height] :
        [].slice.call(arguments)

      return this.attr('viewBox', v)
    }

  }

})
// Add events to elements
;[  'click'
  , 'dblclick'
  , 'mousedown'
  , 'mouseup'
  , 'mouseover'
  , 'mouseout'
  , 'mousemove'
  // , 'mouseenter' -> not supported by IE
  // , 'mouseleave' -> not supported by IE
  , 'touchstart'
  , 'touchmove'
  , 'touchleave'
  , 'touchend'
  , 'touchcancel' ].forEach(function(event) {

  // add event to SVG.Element
  SVG.Element.prototype[event] = function(f) {
    var self = this

    // bind event to element rather than element node
    this.node['on' + event] = typeof f == 'function' ?
      function() { return f.apply(self, arguments) } : null

    return this
  }

})

// Initialize listeners stack
SVG.listeners = []
SVG.handlerMap = []
SVG.listenerId = 0

// Add event binder in the SVG namespace
SVG.on = function(node, event, listener, binding) {
  // create listener, get object-index
  var l     = listener.bind(binding || node.instance || node)
    , index = (SVG.handlerMap.indexOf(node) + 1 || SVG.handlerMap.push(node)) - 1
    , ev    = event.split('.')[0]
    , ns    = event.split('.')[1] || '*'


  // ensure valid object
  SVG.listeners[index]         = SVG.listeners[index]         || {}
  SVG.listeners[index][ev]     = SVG.listeners[index][ev]     || {}
  SVG.listeners[index][ev][ns] = SVG.listeners[index][ev][ns] || {}

  if(!listener._svgjsListenerId)
    listener._svgjsListenerId = ++SVG.listenerId

  // reference listener
  SVG.listeners[index][ev][ns][listener._svgjsListenerId] = l

  // add listener
  node.addEventListener(ev, l, false)
}

// Add event unbinder in the SVG namespace
SVG.off = function(node, event, listener) {
  var index = SVG.handlerMap.indexOf(node)
    , ev    = event && event.split('.')[0]
    , ns    = event && event.split('.')[1]

  if(index == -1) return

  if (listener) {
    if(typeof listener == 'function') listener = listener._svgjsListenerId
    if(!listener) return

    // remove listener reference
    if (SVG.listeners[index][ev] && SVG.listeners[index][ev][ns || '*']) {
      // remove listener
      node.removeEventListener(ev, SVG.listeners[index][ev][ns || '*'][listener], false)

      delete SVG.listeners[index][ev][ns || '*'][listener]
    }

  } else if (ns && ev) {
    // remove all listeners for a namespaced event
    if (SVG.listeners[index][ev] && SVG.listeners[index][ev][ns]) {
      for (listener in SVG.listeners[index][ev][ns])
        SVG.off(node, [ev, ns].join('.'), listener)

      delete SVG.listeners[index][ev][ns]
    }

  } else if (ns){
    // remove all listeners for a specific namespace
    for(event in SVG.listeners[index]){
        for(namespace in SVG.listeners[index][event]){
            if(ns === namespace){
                SVG.off(node, [event, ns].join('.'))
            }
        }
    }

  } else if (ev) {
    // remove all listeners for the event
    if (SVG.listeners[index][ev]) {
      for (namespace in SVG.listeners[index][ev])
        SVG.off(node, [ev, namespace].join('.'))

      delete SVG.listeners[index][ev]
    }

  } else {
    // remove all listeners on a given node
    for (event in SVG.listeners[index])
      SVG.off(node, event)

    delete SVG.listeners[index]

  }
}

//
SVG.extend(SVG.Element, {
  // Bind given event to listener
  on: function(event, listener, binding) {
    SVG.on(this.node, event, listener, binding)

    return this
  }
  // Unbind event from listener
, off: function(event, listener) {
    SVG.off(this.node, event, listener)

    return this
  }
  // Fire given event
, fire: function(event, data) {

    // Dispatch event
    if(event instanceof Event){
        this.node.dispatchEvent(event)
    }else{
        this.node.dispatchEvent(new CustomEvent(event, {detail:data}))
    }

    return this
  }
})

SVG.Defs = SVG.invent({
  // Initialize node
  create: 'defs'

  // Inherit from
, inherit: SVG.Container

})
SVG.G = SVG.invent({
  // Initialize node
  create: 'g'

  // Inherit from
, inherit: SVG.Container

  // Add class methods
, extend: {
    // Move over x-axis
    x: function(x) {
      return x == null ? this.transform('x') : this.transform({ x: x - this.x() }, true)
    }
    // Move over y-axis
  , y: function(y) {
      return y == null ? this.transform('y') : this.transform({ y: y - this.y() }, true)
    }
    // Move by center over x-axis
  , cx: function(x) {
      return x == null ? this.gbox().cx : this.x(x - this.gbox().width / 2)
    }
    // Move by center over y-axis
  , cy: function(y) {
      return y == null ? this.gbox().cy : this.y(y - this.gbox().height / 2)
    }
  , gbox: function() {

      var bbox  = this.bbox()
        , trans = this.transform()

      bbox.x  += trans.x
      bbox.x2 += trans.x
      bbox.cx += trans.x

      bbox.y  += trans.y
      bbox.y2 += trans.y
      bbox.cy += trans.y

      return bbox
    }
  }

  // Add parent method
, construct: {
    // Create a group element
    group: function() {
      return this.put(new SVG.G)
    }
  }
})

// ### This module adds backward / forward functionality to elements.

//
SVG.extend(SVG.Element, {
  // Get all siblings, including myself
  siblings: function() {
    return this.parent().children()
  }
  // Get the curent position siblings
, position: function() {
    return this.parent().index(this)
  }
  // Get the next element (will return null if there is none)
, next: function() {
    return this.siblings()[this.position() + 1]
  }
  // Get the next element (will return null if there is none)
, previous: function() {
    return this.siblings()[this.position() - 1]
  }
  // Send given element one step forward
, forward: function() {
    var i = this.position() + 1
      , p = this.parent()

    // move node one step forward
    p.removeElement(this).add(this, i)

    // make sure defs node is always at the top
    if (p instanceof SVG.Doc)
      p.node.appendChild(p.defs().node)

    return this
  }
  // Send given element one step backward
, backward: function() {
    var i = this.position()

    if (i > 0)
      this.parent().removeElement(this).add(this, i - 1)

    return this
  }
  // Send given element all the way to the front
, front: function() {
    var p = this.parent()

    // Move node forward
    p.node.appendChild(this.node)

    // Make sure defs node is always at the top
    if (p instanceof SVG.Doc)
      p.node.appendChild(p.defs().node)

    return this
  }
  // Send given element all the way to the back
, back: function() {
    if (this.position() > 0)
      this.parent().removeElement(this).add(this, 0)

    return this
  }
  // Inserts a given element before the targeted element
, before: function(element) {
    element.remove()

    var i = this.position()

    this.parent().add(element, i)

    return this
  }
  // Insters a given element after the targeted element
, after: function(element) {
    element.remove()

    var i = this.position()

    this.parent().add(element, i + 1)

    return this
  }

})
SVG.Mask = SVG.invent({
  // Initialize node
  create: function() {
    this.constructor.call(this, SVG.create('mask'))

    // keep references to masked elements
    this.targets = []
  }

  // Inherit from
, inherit: SVG.Container

  // Add class methods
, extend: {
    // Unmask all masked elements and remove itself
    remove: function() {
      // unmask all targets
      for (var i = this.targets.length - 1; i >= 0; i--)
        if (this.targets[i])
          this.targets[i].unmask()
      this.targets = []

      // remove mask from parent
      this.parent().removeElement(this)

      return this
    }
  }

  // Add parent method
, construct: {
    // Create masking element
    mask: function() {
      return this.defs().put(new SVG.Mask)
    }
  }
})


SVG.extend(SVG.Element, {
  // Distribute mask to svg element
  maskWith: function(element) {
    // use given mask or create a new one
    this.masker = element instanceof SVG.Mask ? element : this.parent().mask().add(element)

    // store reverence on self in mask
    this.masker.targets.push(this)

    // apply mask
    return this.attr('mask', 'url("#' + this.masker.attr('id') + '")')
  }
  // Unmask element
, unmask: function() {
    delete this.masker
    return this.attr('mask', null)
  }

})

SVG.ClipPath = SVG.invent({
  // Initialize node
  create: function() {
    this.constructor.call(this, SVG.create('clipPath'))

    // keep references to clipped elements
    this.targets = []
  }

  // Inherit from
, inherit: SVG.Container

  // Add class methods
, extend: {
    // Unclip all clipped elements and remove itself
    remove: function() {
      // unclip all targets
      for (var i = this.targets.length - 1; i >= 0; i--)
        if (this.targets[i])
          this.targets[i].unclip()
      this.targets = []

      // remove clipPath from parent
      this.parent().removeElement(this)

      return this
    }
  }

  // Add parent method
, construct: {
    // Create clipping element
    clip: function() {
      return this.defs().put(new SVG.ClipPath)
    }
  }
})

//
SVG.extend(SVG.Element, {
  // Distribute clipPath to svg element
  clipWith: function(element) {
    // use given clip or create a new one
    this.clipper = element instanceof SVG.ClipPath ? element : this.parent().clip().add(element)

    // store reverence on self in mask
    this.clipper.targets.push(this)

    // apply mask
    return this.attr('clip-path', 'url("#' + this.clipper.attr('id') + '")')
  }
  // Unclip element
, unclip: function() {
    delete this.clipper
    return this.attr('clip-path', null)
  }

})
SVG.Gradient = SVG.invent({
  // Initialize node
  create: function(type) {
    this.constructor.call(this, SVG.create(type + 'Gradient'))

    // store type
    this.type = type
  }

  // Inherit from
, inherit: SVG.Container

  // Add class methods
, extend: {
    // Add a color stop
    at: function(offset, color, opacity) {
      return this.put(new SVG.Stop).update(offset, color, opacity)
    }
    // Update gradient
  , update: function(block) {
      // remove all stops
      this.clear()

      // invoke passed block
      if (typeof block == 'function')
        block.call(this, this)

      return this
    }
    // Return the fill id
  , fill: function() {
      return 'url(#' + this.id() + ')'
    }
    // Alias string convertion to fill
  , toString: function() {
      return this.fill()
    }
    // custom attr to handle transform
  , attr: function(a, b, c) {
      if(a == 'transform') a = 'gradientTransform'
      return SVG.Container.prototype.attr.call(this, a, b, c)
    }
  }

  // Add parent method
, construct: {
    // Create gradient element in defs
    gradient: function(type, block) {
      return this.defs().gradient(type, block)
    }
  }
})

// Add animatable methods to both gradient and fx module
SVG.extend(SVG.Gradient, SVG.FX, {
  // From position
  from: function(x, y) {
    return (this._target || this).type == 'radial' ?
      this.attr({ fx: new SVG.Number(x), fy: new SVG.Number(y) }) :
      this.attr({ x1: new SVG.Number(x), y1: new SVG.Number(y) })
  }
  // To position
, to: function(x, y) {
    return (this._target || this).type == 'radial' ?
      this.attr({ cx: new SVG.Number(x), cy: new SVG.Number(y) }) :
      this.attr({ x2: new SVG.Number(x), y2: new SVG.Number(y) })
  }
})

// Base gradient generation
SVG.extend(SVG.Defs, {
  // define gradient
  gradient: function(type, block) {
    return this.put(new SVG.Gradient(type)).update(block)
  }

})

SVG.Stop = SVG.invent({
  // Initialize node
  create: 'stop'

  // Inherit from
, inherit: SVG.Element

  // Add class methods
, extend: {
    // add color stops
    update: function(o) {
      if (typeof o == 'number' || o instanceof SVG.Number) {
        o = {
          offset:  arguments[0]
        , color:   arguments[1]
        , opacity: arguments[2]
        }
      }

      // set attributes
      if (o.opacity != null) this.attr('stop-opacity', o.opacity)
      if (o.color   != null) this.attr('stop-color', o.color)
      if (o.offset  != null) this.attr('offset', new SVG.Number(o.offset))

      return this
    }
  }

})

SVG.Pattern = SVG.invent({
  // Initialize node
  create: 'pattern'

  // Inherit from
, inherit: SVG.Container

  // Add class methods
, extend: {
    // Return the fill id
    fill: function() {
      return 'url(#' + this.id() + ')'
    }
    // Update pattern by rebuilding
  , update: function(block) {
      // remove content
      this.clear()

      // invoke passed block
      if (typeof block == 'function')
        block.call(this, this)

      return this
    }
    // Alias string convertion to fill
  , toString: function() {
      return this.fill()
    }
    // custom attr to handle transform
  , attr: function(a, b, c) {
      if(a == 'transform') a = 'patternTransform'
      return SVG.Container.prototype.attr.call(this, a, b, c)
    }

  }

  // Add parent method
, construct: {
    // Create pattern element in defs
    pattern: function(width, height, block) {
      return this.defs().pattern(width, height, block)
    }
  }
})

SVG.extend(SVG.Defs, {
  // Define gradient
  pattern: function(width, height, block) {
    return this.put(new SVG.Pattern).update(block).attr({
      x:            0
    , y:            0
    , width:        width
    , height:       height
    , patternUnits: 'userSpaceOnUse'
    })
  }

})
SVG.Doc = SVG.invent({
  // Initialize node
  create: function(element) {
    if (element) {
      // ensure the presence of a dom element
      element = typeof element == 'string' ?
        document.getElementById(element) :
        element

      // If the target is an svg element, use that element as the main wrapper.
      // This allows svg.js to work with svg documents as well.
      if (element.nodeName == 'svg') {
        this.constructor.call(this, element)
      } else {
        this.constructor.call(this, SVG.create('svg'))
        element.appendChild(this.node)
        this.size('100%', '100%')
      }

      // set svg element attributes and ensure defs node
      this.namespace().defs()
    }
  }

  // Inherit from
, inherit: SVG.Container

  // Add class methods
, extend: {
    // Add namespaces
    namespace: function() {
      return this
        .attr({ xmlns: SVG.ns, version: '1.1' })
        .attr('xmlns:xlink', SVG.xlink, SVG.xmlns)
        .attr('xmlns:svgjs', SVG.svgjs, SVG.xmlns)
    }
    // Creates and returns defs element
  , defs: function() {
      if (!this._defs) {
        var defs

        // Find or create a defs element in this instance
        if (defs = this.node.getElementsByTagName('defs')[0])
          this._defs = SVG.adopt(defs)
        else
          this._defs = new SVG.Defs

        // Make sure the defs node is at the end of the stack
        this.node.appendChild(this._defs.node)
      }

      return this._defs
    }
    // custom parent method
  , parent: function() {
      return this.node.parentNode.nodeName == '#document' ? null : this.node.parentNode
    }
    // Fix for possible sub-pixel offset. See:
    // https://bugzilla.mozilla.org/show_bug.cgi?id=608812
  , spof: function(spof) {
      var pos = this.node.getScreenCTM()

      if (pos)
        this
          .style('left', (-pos.e % 1) + 'px')
          .style('top',  (-pos.f % 1) + 'px')

      return this
    }

      // Removes the doc from the DOM
  , remove: function() {
      if(this.parent()) {
        this.parent().removeChild(this.node);
      }

      return this;
    }
  }

})

SVG.Shape = SVG.invent({
  // Initialize node
  create: function(element) {
    this.constructor.call(this, element)
  }

  // Inherit from
, inherit: SVG.Element

})

SVG.Bare = SVG.invent({
  // Initialize
  create: function(element, inherit) {
    // construct element
    this.constructor.call(this, SVG.create(element))

    // inherit custom methods
    if (inherit)
      for (var method in inherit.prototype)
        if (typeof inherit.prototype[method] === 'function')
          this[method] = inherit.prototype[method]
  }

  // Inherit from
, inherit: SVG.Element

  // Add methods
, extend: {
    // Insert some plain text
    words: function(text) {
      // remove contents
      while (this.node.hasChildNodes())
        this.node.removeChild(this.node.lastChild)

      // create text node
      this.node.appendChild(document.createTextNode(text))

      return this
    }
  }
})


SVG.extend(SVG.Parent, {
  // Create an element that is not described by SVG.js
  element: function(element, inherit) {
    return this.put(new SVG.Bare(element, inherit))
  }
  // Add symbol element
, symbol: function() {
    return this.defs().element('symbol', SVG.Container)
  }

})
SVG.Use = SVG.invent({
  // Initialize node
  create: 'use'

  // Inherit from
, inherit: SVG.Shape

  // Add class methods
, extend: {
    // Use element as a reference
    element: function(element, file) {
      // Set lined element
      return this.attr('href', (file || '') + '#' + element, SVG.xlink)
    }
  }

  // Add parent method
, construct: {
    // Create a use element
    use: function(element, file) {
      return this.put(new SVG.Use).element(element, file)
    }
  }
})
SVG.Rect = SVG.invent({
  // Initialize node
  create: 'rect'

  // Inherit from
, inherit: SVG.Shape

  // Add parent method
, construct: {
    // Create a rect element
    rect: function(width, height) {
      return this.put(new SVG.Rect()).size(width, height)
    }
  }
})
SVG.Circle = SVG.invent({
  // Initialize node
  create: 'circle'

  // Inherit from
, inherit: SVG.Shape

  // Add parent method
, construct: {
    // Create circle element, based on ellipse
    circle: function(size) {
      return this.put(new SVG.Circle).rx(new SVG.Number(size).divide(2)).move(0, 0)
    }
  }
})

SVG.extend(SVG.Circle, SVG.FX, {
  // Radius x value
  rx: function(rx) {
    return this.attr('r', rx)
  }
  // Alias radius x value
, ry: function(ry) {
    return this.rx(ry)
  }
})

SVG.Ellipse = SVG.invent({
  // Initialize node
  create: 'ellipse'

  // Inherit from
, inherit: SVG.Shape

  // Add parent method
, construct: {
    // Create an ellipse
    ellipse: function(width, height) {
      return this.put(new SVG.Ellipse).size(width, height).move(0, 0)
    }
  }
})

SVG.extend(SVG.Ellipse, SVG.Rect, SVG.FX, {
  // Radius x value
  rx: function(rx) {
    return this.attr('rx', rx)
  }
  // Radius y value
, ry: function(ry) {
    return this.attr('ry', ry)
  }
})

// Add common method
SVG.extend(SVG.Circle, SVG.Ellipse, {
    // Move over x-axis
    x: function(x) {
      return x == null ? this.cx() - this.rx() : this.cx(x + this.rx())
    }
    // Move over y-axis
  , y: function(y) {
      return y == null ? this.cy() - this.ry() : this.cy(y + this.ry())
    }
    // Move by center over x-axis
  , cx: function(x) {
      return x == null ? this.attr('cx') : this.attr('cx', x)
    }
    // Move by center over y-axis
  , cy: function(y) {
      return y == null ? this.attr('cy') : this.attr('cy', y)
    }
    // Set width of element
  , width: function(width) {
      return width == null ? this.rx() * 2 : this.rx(new SVG.Number(width).divide(2))
    }
    // Set height of element
  , height: function(height) {
      return height == null ? this.ry() * 2 : this.ry(new SVG.Number(height).divide(2))
    }
    // Custom size function
  , size: function(width, height) {
      var p = proportionalSize(this, width, height)

      return this
        .rx(new SVG.Number(p.width).divide(2))
        .ry(new SVG.Number(p.height).divide(2))
    }
})
SVG.Line = SVG.invent({
  // Initialize node
  create: 'line'

  // Inherit from
, inherit: SVG.Shape

  // Add class methods
, extend: {
    // Get array
    array: function() {
      return new SVG.PointArray([
        [ this.attr('x1'), this.attr('y1') ]
      , [ this.attr('x2'), this.attr('y2') ]
      ])
    }
    // Overwrite native plot() method
  , plot: function(x1, y1, x2, y2) {
      if (typeof y1 !== 'undefined')
        x1 = { x1: x1, y1: y1, x2: x2, y2: y2 }
      else
        x1 = new SVG.PointArray(x1).toLine()

      return this.attr(x1)
    }
    // Move by left top corner
  , move: function(x, y) {
      return this.attr(this.array().move(x, y).toLine())
    }
    // Set element size to given width and height
  , size: function(width, height) {
      var p = proportionalSize(this, width, height)

      return this.attr(this.array().size(p.width, p.height).toLine())
    }
  }

  // Add parent method
, construct: {
    // Create a line element
    line: function(x1, y1, x2, y2) {
      return this.put(new SVG.Line).plot(x1, y1, x2, y2)
    }
  }
})

SVG.Polyline = SVG.invent({
  // Initialize node
  create: 'polyline'

  // Inherit from
, inherit: SVG.Shape

  // Add parent method
, construct: {
    // Create a wrapped polyline element
    polyline: function(p) {
      return this.put(new SVG.Polyline).plot(p)
    }
  }
})

SVG.Polygon = SVG.invent({
  // Initialize node
  create: 'polygon'

  // Inherit from
, inherit: SVG.Shape

  // Add parent method
, construct: {
    // Create a wrapped polygon element
    polygon: function(p) {
      return this.put(new SVG.Polygon).plot(p)
    }
  }
})

// Add polygon-specific functions
SVG.extend(SVG.Polyline, SVG.Polygon, {
  // Get array
  array: function() {
    return this._array || (this._array = new SVG.PointArray(this.attr('points')))
  }
  // Plot new path
, plot: function(p) {
    return this.attr('points', (this._array = new SVG.PointArray(p)))
  }
  // Move by left top corner
, move: function(x, y) {
    return this.attr('points', this.array().move(x, y))
  }
  // Set element size to given width and height
, size: function(width, height) {
    var p = proportionalSize(this, width, height)

    return this.attr('points', this.array().size(p.width, p.height))
  }

})
// unify all point to point elements
SVG.extend(SVG.Line, SVG.Polyline, SVG.Polygon, {
  // Define morphable array
  morphArray:  SVG.PointArray
  // Move by left top corner over x-axis
, x: function(x) {
    return x == null ? this.bbox().x : this.move(x, this.bbox().y)
  }
  // Move by left top corner over y-axis
, y: function(y) {
    return y == null ? this.bbox().y : this.move(this.bbox().x, y)
  }
  // Set width of element
, width: function(width) {
    var b = this.bbox()

    return width == null ? b.width : this.size(width, b.height)
  }
  // Set height of element
, height: function(height) {
    var b = this.bbox()

    return height == null ? b.height : this.size(b.width, height)
  }
})
SVG.Path = SVG.invent({
  // Initialize node
  create: 'path'

  // Inherit from
, inherit: SVG.Shape

  // Add class methods
, extend: {
    // Define morphable array
    morphArray:  SVG.PathArray
    // Get array
  , array: function() {
      return this._array || (this._array = new SVG.PathArray(this.attr('d')))
    }
    // Plot new poly points
  , plot: function(p) {
      return this.attr('d', (this._array = new SVG.PathArray(p)))
    }
    // Move by left top corner
  , move: function(x, y) {
      return this.attr('d', this.array().move(x, y))
    }
    // Move by left top corner over x-axis
  , x: function(x) {
      return x == null ? this.bbox().x : this.move(x, this.bbox().y)
    }
    // Move by left top corner over y-axis
  , y: function(y) {
      return y == null ? this.bbox().y : this.move(this.bbox().x, y)
    }
    // Set element size to given width and height
  , size: function(width, height) {
      var p = proportionalSize(this, width, height)

      return this.attr('d', this.array().size(p.width, p.height))
    }
    // Set width of element
  , width: function(width) {
      return width == null ? this.bbox().width : this.size(width, this.bbox().height)
    }
    // Set height of element
  , height: function(height) {
      return height == null ? this.bbox().height : this.size(this.bbox().width, height)
    }

  }

  // Add parent method
, construct: {
    // Create a wrapped path element
    path: function(d) {
      return this.put(new SVG.Path).plot(d)
    }
  }
})
SVG.Image = SVG.invent({
  // Initialize node
  create: 'image'

  // Inherit from
, inherit: SVG.Shape

  // Add class methods
, extend: {
    // (re)load image
    load: function(url) {
      if (!url) return this

      var self = this
        , img  = document.createElement('img')

      // preload image
      img.onload = function() {
        var p = self.parent(SVG.Pattern)

        if(p === null) return

        // ensure image size
        if (self.width() == 0 && self.height() == 0)
          self.size(img.width, img.height)

        // ensure pattern size if not set
        if (p && p.width() == 0 && p.height() == 0)
          p.size(self.width(), self.height())

        // callback
        if (typeof self._loaded === 'function')
          self._loaded.call(self, {
            width:  img.width
          , height: img.height
          , ratio:  img.width / img.height
          , url:    url
          })
      }

      img.onerror = function(e){
        if (typeof self._error === 'function'){
            self._error.call(self, e)
        }
      }

      return this.attr('href', (img.src = this.src = url), SVG.xlink)
    }
    // Add loaded callback
  , loaded: function(loaded) {
      this._loaded = loaded
      return this
    }

  , error: function(error) {
      this._error = error
      return this
    }
  }

  // Add parent method
, construct: {
    // create image element, load image and set its size
    image: function(source, width, height) {
      return this.put(new SVG.Image).load(source).size(width || 0, height || width || 0)
    }
  }

})
SVG.Text = SVG.invent({
  // Initialize node
  create: function() {
    this.constructor.call(this, SVG.create('text'))

    this.dom.leading = new SVG.Number(1.3)    // store leading value for rebuilding
    this._rebuild = true                      // enable automatic updating of dy values
    this._build   = false                     // disable build mode for adding multiple lines

    // set default font
    this.attr('font-family', SVG.defaults.attrs['font-family'])
  }

  // Inherit from
, inherit: SVG.Shape

  // Add class methods
, extend: {
    // Move over x-axis
    x: function(x) {
      // act as getter
      if (x == null)
        return this.attr('x')

      // move lines as well if no textPath is present
      if (!this.textPath)
        this.lines().each(function() { if (this.dom.newLined) this.x(x) })

      return this.attr('x', x)
    }
    // Move over y-axis
  , y: function(y) {
      var oy = this.attr('y')
        , o  = typeof oy === 'number' ? oy - this.bbox().y : 0

      // act as getter
      if (y == null)
        return typeof oy === 'number' ? oy - o : oy

      return this.attr('y', typeof y === 'number' ? y + o : y)
    }
    // Move center over x-axis
  , cx: function(x) {
      return x == null ? this.bbox().cx : this.x(x - this.bbox().width / 2)
    }
    // Move center over y-axis
  , cy: function(y) {
      return y == null ? this.bbox().cy : this.y(y - this.bbox().height / 2)
    }
    // Set the text content
  , text: function(text) {
      // act as getter
      if (typeof text === 'undefined'){
        var text = ''
        var children = this.node.childNodes
        for(var i = 0, len = children.length; i < len; ++i){

          // add newline if its not the first child and newLined is set to true
          if(i != 0 && children[i].nodeType != 3 && SVG.adopt(children[i]).dom.newLined == true){
            text += '\n'
          }

          // add content of this node
          text += children[i].textContent
        }

        return text
      }

      // remove existing content
      this.clear().build(true)

      if (typeof text === 'function') {
        // call block
        text.call(this, this)

      } else {
        // store text and make sure text is not blank
        text = text.split('\n')

        // build new lines
        for (var i = 0, il = text.length; i < il; i++)
          this.tspan(text[i]).newLine()
      }

      // disable build mode and rebuild lines
      return this.build(false).rebuild()
    }
    // Set font size
  , size: function(size) {
      return this.attr('font-size', size).rebuild()
    }
    // Set / get leading
  , leading: function(value) {
      // act as getter
      if (value == null)
        return this.dom.leading

      // act as setter
      this.dom.leading = new SVG.Number(value)

      return this.rebuild()
    }
    // Get all the first level lines
  , lines: function() {
      var node = (this.textPath && this.textPath() || this).node

      // filter tspans and map them to SVG.js instances
      var lines = SVG.utils.map(SVG.utils.filterSVGElements(node.childNodes), function(el){
        return SVG.adopt(el)
      })

      // return an instance of SVG.set
      return new SVG.Set(lines)
    }
    // Rebuild appearance type
  , rebuild: function(rebuild) {
      // store new rebuild flag if given
      if (typeof rebuild == 'boolean')
        this._rebuild = rebuild

      // define position of all lines
      if (this._rebuild) {
        var self = this
          , blankLineOffset = 0
          , dy = this.dom.leading * new SVG.Number(this.attr('font-size'))

        this.lines().each(function() {
          if (this.dom.newLined) {
            if (!this.textPath)
              this.attr('x', self.attr('x'))

            if(this.text() == '\n') {
              blankLineOffset += dy
            }else{
              this.attr('dy', dy + blankLineOffset)
              blankLineOffset = 0
            }
          }
        })

        this.fire('rebuild')
      }

      return this
    }
    // Enable / disable build mode
  , build: function(build) {
      this._build = !!build
      return this
    }
    // overwrite method from parent to set data properly
  , setData: function(o){
      this.dom = o
      this.dom.leading = new SVG.Number(o.leading || 1.3)
      return this
    }
  }

  // Add parent method
, construct: {
    // Create text element
    text: function(text) {
      return this.put(new SVG.Text).text(text)
    }
    // Create plain text element
  , plain: function(text) {
      return this.put(new SVG.Text).plain(text)
    }
  }

})

SVG.Tspan = SVG.invent({
  // Initialize node
  create: 'tspan'

  // Inherit from
, inherit: SVG.Shape

  // Add class methods
, extend: {
    // Set text content
    text: function(text) {
      if(text == null) return this.node.textContent + (this.dom.newLined ? '\n' : '')

      typeof text === 'function' ? text.call(this, this) : this.plain(text)

      return this
    }
    // Shortcut dx
  , dx: function(dx) {
      return this.attr('dx', dx)
    }
    // Shortcut dy
  , dy: function(dy) {
      return this.attr('dy', dy)
    }
    // Create new line
  , newLine: function() {
      // fetch text parent
      var t = this.parent(SVG.Text)

      // mark new line
      this.dom.newLined = true

      // apply new hy¡n
      return this.dy(t.dom.leading * t.attr('font-size')).attr('x', t.x())
    }
  }

})

SVG.extend(SVG.Text, SVG.Tspan, {
  // Create plain text node
  plain: function(text) {
    // clear if build mode is disabled
    if (this._build === false)
      this.clear()

    // create text node
    this.node.appendChild(document.createTextNode(text))

    return this
  }
  // Create a tspan
, tspan: function(text) {
    var node  = (this.textPath && this.textPath() || this).node
      , tspan = new SVG.Tspan

    // clear if build mode is disabled
    if (this._build === false)
      this.clear()

    // add new tspan
    node.appendChild(tspan.node)

    return tspan.text(text)
  }
  // Clear all lines
, clear: function() {
    var node = (this.textPath && this.textPath() || this).node

    // remove existing child nodes
    while (node.hasChildNodes())
      node.removeChild(node.lastChild)

    return this
  }
  // Get length of text element
, length: function() {
    return this.node.getComputedTextLength()
  }
})

SVG.TextPath = SVG.invent({
  // Initialize node
  create: 'textPath'

  // Inherit from
, inherit: SVG.Parent

  // Define parent class
, parent: SVG.Text

  // Add parent method
, construct: {
    // Create path for text to run on
    path: function(d) {
      // create textPath element
      var path  = new SVG.TextPath
        , track = this.doc().defs().path(d)

      // move lines to textpath
      while (this.node.hasChildNodes())
        path.node.appendChild(this.node.firstChild)

      // add textPath element as child node
      this.node.appendChild(path.node)

      // link textPath to path and add content
      path.attr('href', '#' + track, SVG.xlink)

      return this
    }
    // Plot path if any
  , plot: function(d) {
      var track = this.track()

      if (track)
        track.plot(d)

      return this
    }
    // Get the path track element
  , track: function() {
      var path = this.textPath()

      if (path)
        return path.reference('href')
    }
    // Get the textPath child
  , textPath: function() {
      if (this.node.firstChild && this.node.firstChild.nodeName == 'textPath')
        return SVG.adopt(this.node.firstChild)
    }
  }
})
SVG.Nested = SVG.invent({
  // Initialize node
  create: function() {
    this.constructor.call(this, SVG.create('svg'))

    this.style('overflow', 'visible')
  }

  // Inherit from
, inherit: SVG.Container

  // Add parent method
, construct: {
    // Create nested svg document
    nested: function() {
      return this.put(new SVG.Nested)
    }
  }
})
SVG.A = SVG.invent({
  // Initialize node
  create: 'a'

  // Inherit from
, inherit: SVG.Container

  // Add class methods
, extend: {
    // Link url
    to: function(url) {
      return this.attr('href', url, SVG.xlink)
    }
    // Link show attribute
  , show: function(target) {
      return this.attr('show', target, SVG.xlink)
    }
    // Link target attribute
  , target: function(target) {
      return this.attr('target', target)
    }
  }

  // Add parent method
, construct: {
    // Create a hyperlink element
    link: function(url) {
      return this.put(new SVG.A).to(url)
    }
  }
})

SVG.extend(SVG.Element, {
  // Create a hyperlink element
  linkTo: function(url) {
    var link = new SVG.A

    if (typeof url == 'function')
      url.call(link, link)
    else
      link.to(url)

    return this.parent().put(link).put(this)
  }

})
SVG.Marker = SVG.invent({
  // Initialize node
  create: 'marker'

  // Inherit from
, inherit: SVG.Container

  // Add class methods
, extend: {
    // Set width of element
    width: function(width) {
      return this.attr('markerWidth', width)
    }
    // Set height of element
  , height: function(height) {
      return this.attr('markerHeight', height)
    }
    // Set marker refX and refY
  , ref: function(x, y) {
      return this.attr('refX', x).attr('refY', y)
    }
    // Update marker
  , update: function(block) {
      // remove all content
      this.clear()

      // invoke passed block
      if (typeof block == 'function')
        block.call(this, this)

      return this
    }
    // Return the fill id
  , toString: function() {
      return 'url(#' + this.id() + ')'
    }
  }

  // Add parent method
, construct: {
    marker: function(width, height, block) {
      // Create marker element in defs
      return this.defs().marker(width, height, block)
    }
  }

})

SVG.extend(SVG.Defs, {
  // Create marker
  marker: function(width, height, block) {
    // Set default viewbox to match the width and height, set ref to cx and cy and set orient to auto
    return this.put(new SVG.Marker)
      .size(width, height)
      .ref(width / 2, height / 2)
      .viewbox(0, 0, width, height)
      .attr('orient', 'auto')
      .update(block)
  }

})

SVG.extend(SVG.Line, SVG.Polyline, SVG.Polygon, SVG.Path, {
  // Create and attach markers
  marker: function(marker, width, height, block) {
    var attr = ['marker']

    // Build attribute name
    if (marker != 'all') attr.push(marker)
    attr = attr.join('-')

    // Set marker attribute
    marker = arguments[1] instanceof SVG.Marker ?
      arguments[1] :
      this.doc().marker(width, height, block)

    return this.attr(attr, marker)
  }

})
// Define list of available attributes for stroke and fill
var sugar = {
  stroke: ['color', 'width', 'opacity', 'linecap', 'linejoin', 'miterlimit', 'dasharray', 'dashoffset']
, fill:   ['color', 'opacity', 'rule']
, prefix: function(t, a) {
    return a == 'color' ? t : t + '-' + a
  }
}

// Add sugar for fill and stroke
;['fill', 'stroke'].forEach(function(m) {
  var i, extension = {}

  extension[m] = function(o) {
    if (typeof o == 'undefined')
      return this
    if (typeof o == 'string' || SVG.Color.isRgb(o) || (o && typeof o.fill === 'function'))
      this.attr(m, o)

    else
      // set all attributes from sugar.fill and sugar.stroke list
      for (i = sugar[m].length - 1; i >= 0; i--)
        if (o[sugar[m][i]] != null)
          this.attr(sugar.prefix(m, sugar[m][i]), o[sugar[m][i]])

    return this
  }

  SVG.extend(SVG.Element, SVG.FX, extension)

})

SVG.extend(SVG.Element, SVG.FX, {
  // Map rotation to transform
  rotate: function(d, cx, cy) {
    return this.transform({ rotation: d, cx: cx, cy: cy })
  }
  // Map skew to transform
, skew: function(x, y, cx, cy) {
    return arguments.length == 1  || arguments.length == 3 ?
      this.transform({ skew: x, cx: y, cy: cx }) :
      this.transform({ skewX: x, skewY: y, cx: cx, cy: cy })
  }
  // Map scale to transform
, scale: function(x, y, cx, cy) {
    return arguments.length == 1  || arguments.length == 3 ?
      this.transform({ scale: x, cx: y, cy: cx }) :
      this.transform({ scaleX: x, scaleY: y, cx: cx, cy: cy })
  }
  // Map translate to transform
, translate: function(x, y) {
    return this.transform({ x: x, y: y })
  }
  // Map flip to transform
, flip: function(a, o) {
    return this.transform({ flip: a, offset: o })
  }
  // Map matrix to transform
, matrix: function(m) {
    return this.attr('transform', new SVG.Matrix(m))
  }
  // Opacity
, opacity: function(value) {
    return this.attr('opacity', value)
  }
  // Relative move over x axis
, dx: function(x) {
    return this.x((this instanceof SVG.FX ? 0 : this.x()) + x, true)
  }
  // Relative move over y axis
, dy: function(y) {
    return this.y((this instanceof SVG.FX ? 0 : this.y()) + y, true)
  }
  // Relative move over x and y axes
, dmove: function(x, y) {
    return this.dx(x).dy(y)
  }
})

SVG.extend(SVG.Rect, SVG.Ellipse, SVG.Circle, SVG.Gradient, SVG.FX, {
  // Add x and y radius
  radius: function(x, y) {
    var type = (this._target || this).type;
    return type == 'radial' || type == 'circle' ?
      this.attr('r', new SVG.Number(x)) :
      this.rx(x).ry(y == null ? x : y)
  }
})

SVG.extend(SVG.Path, {
  // Get path length
  length: function() {
    return this.node.getTotalLength()
  }
  // Get point at length
, pointAt: function(length) {
    return this.node.getPointAtLength(length)
  }
})

SVG.extend(SVG.Parent, SVG.Text, SVG.FX, {
  // Set font
  font: function(o) {
    for (var k in o)
      k == 'leading' ?
        this.leading(o[k]) :
      k == 'anchor' ?
        this.attr('text-anchor', o[k]) :
      k == 'size' || k == 'family' || k == 'weight' || k == 'stretch' || k == 'variant' || k == 'style' ?
        this.attr('font-'+ k, o[k]) :
        this.attr(k, o[k])

    return this
  }
})

SVG.Set = SVG.invent({
  // Initialize
  create: function(members) {
    // Set initial state
    Array.isArray(members) ? this.members = members : this.clear()
  }

  // Add class methods
, extend: {
    // Add element to set
    add: function() {
      var i, il, elements = [].slice.call(arguments)

      for (i = 0, il = elements.length; i < il; i++)
        this.members.push(elements[i])

      return this
    }
    // Remove element from set
  , remove: function(element) {
      var i = this.index(element)

      // remove given child
      if (i > -1)
        this.members.splice(i, 1)

      return this
    }
    // Iterate over all members
  , each: function(block) {
      for (var i = 0, il = this.members.length; i < il; i++)
        block.apply(this.members[i], [i, this.members])

      return this
    }
    // Restore to defaults
  , clear: function() {
      // initialize store
      this.members = []

      return this
    }
    // Get the length of a set
  , length: function() {
      return this.members.length
    }
    // Checks if a given element is present in set
  , has: function(element) {
      return this.index(element) >= 0
    }
    // retuns index of given element in set
  , index: function(element) {
      return this.members.indexOf(element)
    }
    // Get member at given index
  , get: function(i) {
      return this.members[i]
    }
    // Get first member
  , first: function() {
      return this.get(0)
    }
    // Get last member
  , last: function() {
      return this.get(this.members.length - 1)
    }
    // Default value
  , valueOf: function() {
      return this.members
    }
    // Get the bounding box of all members included or empty box if set has no items
  , bbox: function(){
      var box = new SVG.BBox()

      // return an empty box of there are no members
      if (this.members.length == 0)
        return box

      // get the first rbox and update the target bbox
      var rbox = this.members[0].rbox()
      box.x      = rbox.x
      box.y      = rbox.y
      box.width  = rbox.width
      box.height = rbox.height

      this.each(function() {
        // user rbox for correct position and visual representation
        box = box.merge(this.rbox())
      })

      return box
    }
  }

  // Add parent method
, construct: {
    // Create a new set
    set: function(members) {
      return new SVG.Set(members)
    }
  }
})

SVG.FX.Set = SVG.invent({
  // Initialize node
  create: function(set) {
    // store reference to set
    this.set = set
  }

})

// Alias methods
SVG.Set.inherit = function() {
  var m
    , methods = []

  // gather shape methods
  for(var m in SVG.Shape.prototype)
    if (typeof SVG.Shape.prototype[m] == 'function' && typeof SVG.Set.prototype[m] != 'function')
      methods.push(m)

  // apply shape aliasses
  methods.forEach(function(method) {
    SVG.Set.prototype[method] = function() {
      for (var i = 0, il = this.members.length; i < il; i++)
        if (this.members[i] && typeof this.members[i][method] == 'function')
          this.members[i][method].apply(this.members[i], arguments)

      return method == 'animate' ? (this.fx || (this.fx = new SVG.FX.Set(this))) : this
    }
  })

  // clear methods for the next round
  methods = []

  // gather fx methods
  for(var m in SVG.FX.prototype)
    if (typeof SVG.FX.prototype[m] == 'function' && typeof SVG.FX.Set.prototype[m] != 'function')
      methods.push(m)

  // apply fx aliasses
  methods.forEach(function(method) {
    SVG.FX.Set.prototype[method] = function() {
      for (var i = 0, il = this.set.members.length; i < il; i++)
        this.set.members[i].fx[method].apply(this.set.members[i].fx, arguments)

      return this
    }
  })
}




SVG.extend(SVG.Element, {
  // Store data values on svg nodes
  data: function(a, v, r) {
    if (typeof a == 'object') {
      for (v in a)
        this.data(v, a[v])

    } else if (arguments.length < 2) {
      try {
        return JSON.parse(this.attr('data-' + a))
      } catch(e) {
        return this.attr('data-' + a)
      }

    } else {
      this.attr(
        'data-' + a
      , v === null ?
          null :
        r === true || typeof v === 'string' || typeof v === 'number' ?
          v :
          JSON.stringify(v)
      )
    }

    return this
  }
})
SVG.extend(SVG.Element, {
  // Remember arbitrary data
  remember: function(k, v) {
    // remember every item in an object individually
    if (typeof arguments[0] == 'object')
      for (var v in k)
        this.remember(v, k[v])

    // retrieve memory
    else if (arguments.length == 1)
      return this.memory()[k]

    // store memory
    else
      this.memory()[k] = v

    return this
  }

  // Erase a given memory
, forget: function() {
    if (arguments.length == 0)
      this._memory = {}
    else
      for (var i = arguments.length - 1; i >= 0; i--)
        delete this.memory()[arguments[i]]

    return this
  }

  // Initialize or return local memory object
, memory: function() {
    return this._memory || (this._memory = {})
  }

})
// Method for getting an element by id
SVG.get = function(id) {
  var node = document.getElementById(idFromReference(id) || id)
  return SVG.adopt(node)
}

// Select elements by query string
SVG.select = function(query, parent) {
  return new SVG.Set(
    SVG.utils.map((parent || document).querySelectorAll(query), function(node) {
      return SVG.adopt(node)
    })
  )
}

SVG.extend(SVG.Parent, {
  // Scoped select method
  select: function(query) {
    return SVG.select(query, this.node)
  }

})
function is(el, obj){
  return el instanceof obj
}

// tests if a given selector matches an element
function matches(el, selector) {
  return (el.matches || el.matchesSelector || el.msMatchesSelector || el.mozMatchesSelector || el.webkitMatchesSelector || el.oMatchesSelector).call(el, selector);
}

// Convert dash-separated-string to camelCase
function camelCase(s) {
  return s.toLowerCase().replace(/-(.)/g, function(m, g) {
    return g.toUpperCase()
  })
}

// Capitalize first letter of a string
function capitalize(s) {
  return s.charAt(0).toUpperCase() + s.slice(1)
}

// Ensure to six-based hex
function fullHex(hex) {
  return hex.length == 4 ?
    [ '#',
      hex.substring(1, 2), hex.substring(1, 2)
    , hex.substring(2, 3), hex.substring(2, 3)
    , hex.substring(3, 4), hex.substring(3, 4)
    ].join('') : hex
}

// Component to hex value
function compToHex(comp) {
  var hex = comp.toString(16)
  return hex.length == 1 ? '0' + hex : hex
}

// Calculate proportional width and height values when necessary
function proportionalSize(element, width, height) {
  if (width == null || height == null) {
    var box = element.bbox()

    if (width == null)
      width = box.width / box.height * height
    else if (height == null)
      height = box.height / box.width * width
  }

  return {
    width:  width
  , height: height
  }
}

// Delta transform point
function deltaTransformPoint(matrix, x, y) {
  return {
    x: x * matrix.a + y * matrix.c + 0
  , y: x * matrix.b + y * matrix.d + 0
  }
}

// Map matrix array to object
function arrayToMatrix(a) {
  return { a: a[0], b: a[1], c: a[2], d: a[3], e: a[4], f: a[5] }
}

// Parse matrix if required
function parseMatrix(matrix) {
  if (!(matrix instanceof SVG.Matrix))
    matrix = new SVG.Matrix(matrix)

  return matrix
}

// Add centre point to transform object
function ensureCentre(o, target) {
  o.cx = o.cx == null ? target.bbox().cx : o.cx
  o.cy = o.cy == null ? target.bbox().cy : o.cy
}

// Convert string to matrix
function stringToMatrix(source) {
  // remove matrix wrapper and split to individual numbers
  source = source
    .replace(SVG.regex.whitespace, '')
    .replace(SVG.regex.matrix, '')
    .split(SVG.regex.matrixElements)

  // convert string values to floats and convert to a matrix-formatted object
  return arrayToMatrix(
    SVG.utils.map(source, function(n) {
      return parseFloat(n)
    })
  )
}

// Calculate position according to from and to
function at(o, pos) {
  // number recalculation (don't bother converting to SVG.Number for performance reasons)
  return typeof o.from == 'number' ?
    o.from + (o.to - o.from) * pos :

  // instance recalculation
  o instanceof SVG.Color || o instanceof SVG.Number || o instanceof SVG.Matrix ? o.at(pos) :

  // for all other values wait until pos has reached 1 to return the final value
  pos < 1 ? o.from : o.to
}

// PathArray Helpers
function arrayToString(a) {
  for (var i = 0, il = a.length, s = ''; i < il; i++) {
    s += a[i][0]

    if (a[i][1] != null) {
      s += a[i][1]

      if (a[i][2] != null) {
        s += ' '
        s += a[i][2]

        if (a[i][3] != null) {
          s += ' '
          s += a[i][3]
          s += ' '
          s += a[i][4]

          if (a[i][5] != null) {
            s += ' '
            s += a[i][5]
            s += ' '
            s += a[i][6]

            if (a[i][7] != null) {
              s += ' '
              s += a[i][7]
            }
          }
        }
      }
    }
  }

  return s + ' '
}

// Deep new id assignment
function assignNewId(node) {
  // do the same for SVG child nodes as well
  for (var i = node.childNodes.length - 1; i >= 0; i--)
    if (node.childNodes[i] instanceof SVGElement)
      assignNewId(node.childNodes[i])

  return SVG.adopt(node).id(SVG.eid(node.nodeName))
}

// Add more bounding box properties
function fullBox(b) {
  if (b.x == null) {
    b.x      = 0
    b.y      = 0
    b.width  = 0
    b.height = 0
  }

  b.w  = b.width
  b.h  = b.height
  b.x2 = b.x + b.width
  b.y2 = b.y + b.height
  b.cx = b.x + b.width / 2
  b.cy = b.y + b.height / 2

  return b
}

// Get id from reference string
function idFromReference(url) {
  var m = url.toString().match(SVG.regex.reference)

  if (m) return m[1]
}

// Create matrix array for looping
var abcdef = 'abcdef'.split('')
// Add CustomEvent to IE9 and IE10
if (typeof CustomEvent !== 'function') {
  // Code from: https://developer.mozilla.org/en-US/docs/Web/API/CustomEvent
  var CustomEvent = function(event, options) {
    options = options || { bubbles: false, cancelable: false, detail: undefined }
    var e = document.createEvent('CustomEvent')
    e.initCustomEvent(event, options.bubbles, options.cancelable, options.detail)
    return e
  }

  CustomEvent.prototype = window.Event.prototype

  window.CustomEvent = CustomEvent
}

// requestAnimationFrame / cancelAnimationFrame Polyfill with fallback based on Paul Irish
(function(w) {
  var lastTime = 0
  var vendors = ['moz', 'webkit']

  for(var x = 0; x < vendors.length && !window.requestAnimationFrame; ++x) {
    w.requestAnimationFrame = w[vendors[x] + 'RequestAnimationFrame']
    w.cancelAnimationFrame  = w[vendors[x] + 'CancelAnimationFrame'] ||
                              w[vendors[x] + 'CancelRequestAnimationFrame']
  }

  w.requestAnimationFrame = w.requestAnimationFrame ||
    function(callback) {
      var currTime = new Date().getTime()
      var timeToCall = Math.max(0, 16 - (currTime - lastTime))

      var id = w.setTimeout(function() {
        callback(currTime + timeToCall)
      }, timeToCall)

      lastTime = currTime + timeToCall
      return id
    }

  w.cancelAnimationFrame = w.cancelAnimationFrame || w.clearTimeout;

}(window))

return SVG

}));
},{}],3:[function(require,module,exports){
'use strict';

Object.defineProperty(exports, "__esModule", {
	value: true
});

var _fitCurve = require('fit-curve');

var _fitCurve2 = _interopRequireDefault(_fitCurve);

function _interopRequireDefault(obj) { return obj && obj.__esModule ? obj : { default: obj }; }

var error = 100;

function PaintControl(pannel) {
	var rawPointData = [];
	var paintingPolyLine = undefined;

	this.start = function (point) {
		rawPointData.push(point);
		paintingPolyLine = pannel.polyline().fill('none').stroke({ width: 1 });
	};
	this.update = function (point) {
		rawPointData.push(point);
		updateLines(paintingPolyLine, rawPointData);
	};

	this.end = function () {
		var smoothBizer = (0, _fitCurve2.default)(rawPointData, error);
		var pathString = fittedCurveToPathString(smoothBizer);

		// draw magnetic curve
		drawOnPannel(pannel, pathString);
		clearRawData();
	};

	function updateLines(paintingPolyLine, rawPointData) {
		paintingPolyLine.plot(rawPointData);
	}
	function fittedCurveToPathString(fittedLineData) {
		var str = '';
		//bezier : [ [c0], [c1], [c2], [c3] ]
		fittedLineData.map(function (bezier, i) {
			if (i == 0) {
				str += 'M ' + bezier[0][0] + ' ' + bezier[0][1];
			}

			str += 'C ' + bezier[1][0] + ' ' + bezier[1][1] + ', ' + bezier[2][0] + ' ' + bezier[2][1] + ', ' + bezier[3][0] + ' ' + bezier[3][1] + ' ';
		});

		return str;
	}
	function drawOnPannel(pannel, pathString) {
		pannel.path(pathString).fill('none').stroke({ width: 3 }).stroke('#f06');
	}
	function clearRawData() {
		rawPointData.length = 0;
		paintingPolyLine.remove();
	}
}

exports.default = PaintControl;

},{"fit-curve":1}],4:[function(require,module,exports){
'use strict';

var _PaintControl = require('./Controls/PaintControl');

var _PaintControl2 = _interopRequireDefault(_PaintControl);

var _svg = require('svg.js');

var _svg2 = _interopRequireDefault(_svg);

function _interopRequireDefault(obj) { return obj && obj.__esModule ? obj : { default: obj }; }

var draw = (0, _svg2.default)('drawing').size(300, 300);

setControl(draw);

function setControl(_container) {
	var isMouseDown = false;
	var currnetControl = new _PaintControl2.default(draw);
	var top = draw.node.getBoundingClientRect().top;
	var left = draw.node.getBoundingClientRect().left;

	_container.on('mousedown', function (e) {
		var point = [e.clientX - top, e.clientY - left];
		isMouseDown = true;
		currnetControl.start(point);
	});
	_container.on('mouseup', function () {
		isMouseDown = false;
		currnetControl.end();
	});
	_container.on('mousemove', function (e) {
		var x = e.offsetX;
		var y = e.offsetY;
		if (isMouseDown) {
			currnetControl.update([x, y]);
		}
	});
}

},{"./Controls/PaintControl":3,"svg.js":2}]},{},[4])
//# sourceMappingURL=data:application/json;charset=utf-8;base64,eyJ2ZXJzaW9uIjozLCJzb3VyY2VzIjpbIm5vZGVfbW9kdWxlcy9icm93c2VyLXBhY2svX3ByZWx1ZGUuanMiLCJub2RlX21vZHVsZXMvZml0LWN1cnZlL3NyYy9maXQtY3VydmUuanMiLCJub2RlX21vZHVsZXMvc3ZnLmpzL2Rpc3Qvc3ZnLmpzIiwic3JjL0NvbnRyb2xzL1BhaW50Q29udHJvbC5qcyIsInNyYy9tYWluLmpzIl0sIm5hbWVzIjpbXSwibWFwcGluZ3MiOiJBQUFBO0FDQUE7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7O0FDMWpCQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7Ozs7Ozs7O0FDLzFLQTs7Ozs7O0FBRUEsSUFBTSxRQUFRLEdBQWQ7O0FBRUEsU0FBUyxZQUFULENBQXNCLE1BQXRCLEVBQThCO0FBQzdCLEtBQUksZUFBZSxFQUFuQjtBQUNBLEtBQUksbUJBQW1CLFNBQXZCOztBQUVBLE1BQUssS0FBTCxHQUFhLFVBQVUsS0FBVixFQUFrQjtBQUM5QixlQUFhLElBQWIsQ0FBbUIsS0FBbkI7QUFDQSxxQkFBbUIsT0FBTyxRQUFQLEdBQWtCLElBQWxCLENBQXVCLE1BQXZCLEVBQStCLE1BQS9CLENBQXNDLEVBQUUsT0FBTyxDQUFULEVBQXRDLENBQW5CO0FBRUEsRUFKRDtBQUtBLE1BQUssTUFBTCxHQUFjLFVBQVUsS0FBVixFQUFrQjtBQUMvQixlQUFhLElBQWIsQ0FBbUIsS0FBbkI7QUFDQSxjQUFhLGdCQUFiLEVBQStCLFlBQS9CO0FBQ0EsRUFIRDs7QUFLQSxNQUFLLEdBQUwsR0FBVyxZQUFXO0FBQ3JCLE1BQUksY0FBYyx3QkFBVSxZQUFWLEVBQXdCLEtBQXhCLENBQWxCO0FBQ0EsTUFBSSxhQUFhLHdCQUF3QixXQUF4QixDQUFqQjs7QUFFQTtBQUNBLGVBQWEsTUFBYixFQUFxQixVQUFyQjtBQUNBO0FBQ0EsRUFQRDs7QUFTQSxVQUFTLFdBQVQsQ0FBcUIsZ0JBQXJCLEVBQXVDLFlBQXZDLEVBQXFEO0FBQ3BELG1CQUFpQixJQUFqQixDQUF1QixZQUF2QjtBQUNBO0FBQ0QsVUFBUyx1QkFBVCxDQUFpQyxjQUFqQyxFQUFpRDtBQUNoRCxNQUFJLE1BQU0sRUFBVjtBQUNBO0FBQ0EsaUJBQWUsR0FBZixDQUFtQixVQUFVLE1BQVYsRUFBa0IsQ0FBbEIsRUFBcUI7QUFDdkMsT0FBSSxLQUFLLENBQVQsRUFBWTtBQUNYLFdBQU8sT0FBTyxPQUFPLENBQVAsRUFBVSxDQUFWLENBQVAsR0FBc0IsR0FBdEIsR0FBNEIsT0FBTyxDQUFQLEVBQVUsQ0FBVixDQUFuQztBQUNBOztBQUVELFVBQU8sT0FBTyxPQUFPLENBQVAsRUFBVSxDQUFWLENBQVAsR0FBc0IsR0FBdEIsR0FBNEIsT0FBTyxDQUFQLEVBQVUsQ0FBVixDQUE1QixHQUEyQyxJQUEzQyxHQUNQLE9BQU8sQ0FBUCxFQUFVLENBQVYsQ0FETyxHQUNRLEdBRFIsR0FDYyxPQUFPLENBQVAsRUFBVSxDQUFWLENBRGQsR0FDNkIsSUFEN0IsR0FFUCxPQUFPLENBQVAsRUFBVSxDQUFWLENBRk8sR0FFUSxHQUZSLEdBRWMsT0FBTyxDQUFQLEVBQVUsQ0FBVixDQUZkLEdBRTZCLEdBRnBDO0FBSUEsR0FURDs7QUFXQSxTQUFPLEdBQVA7QUFDQTtBQUNELFVBQVMsWUFBVCxDQUFzQixNQUF0QixFQUE4QixVQUE5QixFQUF5QztBQUN4QyxTQUFPLElBQVAsQ0FBYSxVQUFiLEVBQTBCLElBQTFCLENBQStCLE1BQS9CLEVBQXVDLE1BQXZDLENBQThDLEVBQUUsT0FBTyxDQUFULEVBQTlDLEVBQTRELE1BQTVELENBQW1FLE1BQW5FO0FBQ0E7QUFDRCxVQUFTLFlBQVQsR0FBdUI7QUFDdEIsZUFBYSxNQUFiLEdBQXNCLENBQXRCO0FBQ0EsbUJBQWlCLE1BQWpCO0FBQ0E7QUFDRDs7a0JBRWMsWTs7Ozs7QUN2RGY7Ozs7QUFDQTs7Ozs7O0FBRUEsSUFBSSxPQUFPLG1CQUFJLFNBQUosRUFBZSxJQUFmLENBQW9CLEdBQXBCLEVBQXlCLEdBQXpCLENBQVg7O0FBR0EsV0FBVyxJQUFYOztBQUVBLFNBQVMsVUFBVCxDQUFvQixVQUFwQixFQUFnQztBQUMvQixLQUFJLGNBQWMsS0FBbEI7QUFDQSxLQUFJLGlCQUFpQiwyQkFBaUIsSUFBakIsQ0FBckI7QUFDQSxLQUFNLE1BQU0sS0FBSyxJQUFMLENBQVUscUJBQVYsR0FBa0MsR0FBOUM7QUFDQSxLQUFNLE9BQU8sS0FBSyxJQUFMLENBQVUscUJBQVYsR0FBa0MsSUFBL0M7O0FBRUEsWUFBVyxFQUFYLENBQWMsV0FBZCxFQUEyQixVQUFVLENBQVYsRUFBYTtBQUN2QyxNQUFNLFFBQVEsQ0FDYixFQUFFLE9BQUYsR0FBWSxHQURDLEVBRWIsRUFBRSxPQUFGLEdBQVksSUFGQyxDQUFkO0FBSUEsZ0JBQWMsSUFBZDtBQUNBLGlCQUFlLEtBQWYsQ0FBcUIsS0FBckI7QUFFQSxFQVJEO0FBU0EsWUFBVyxFQUFYLENBQWMsU0FBZCxFQUF5QixZQUFZO0FBQ3BDLGdCQUFjLEtBQWQ7QUFDQSxpQkFBZSxHQUFmO0FBQ0EsRUFIRDtBQUlBLFlBQVcsRUFBWCxDQUFjLFdBQWQsRUFBMkIsVUFBVSxDQUFWLEVBQWE7QUFDdkMsTUFBSSxJQUFJLEVBQUUsT0FBVjtBQUNBLE1BQUksSUFBSSxFQUFFLE9BQVY7QUFDQSxNQUFJLFdBQUosRUFBaUI7QUFDaEIsa0JBQWUsTUFBZixDQUFzQixDQUFDLENBQUQsRUFBSSxDQUFKLENBQXRCO0FBQ0E7QUFDRCxFQU5EO0FBT0EiLCJmaWxlIjoiZ2VuZXJhdGVkLmpzIiwic291cmNlUm9vdCI6IiIsInNvdXJjZXNDb250ZW50IjpbIihmdW5jdGlvbiBlKHQsbixyKXtmdW5jdGlvbiBzKG8sdSl7aWYoIW5bb10pe2lmKCF0W29dKXt2YXIgYT10eXBlb2YgcmVxdWlyZT09XCJmdW5jdGlvblwiJiZyZXF1aXJlO2lmKCF1JiZhKXJldHVybiBhKG8sITApO2lmKGkpcmV0dXJuIGkobywhMCk7dmFyIGY9bmV3IEVycm9yKFwiQ2Fubm90IGZpbmQgbW9kdWxlICdcIitvK1wiJ1wiKTt0aHJvdyBmLmNvZGU9XCJNT0RVTEVfTk9UX0ZPVU5EXCIsZn12YXIgbD1uW29dPXtleHBvcnRzOnt9fTt0W29dWzBdLmNhbGwobC5leHBvcnRzLGZ1bmN0aW9uKGUpe3ZhciBuPXRbb11bMV1bZV07cmV0dXJuIHMobj9uOmUpfSxsLGwuZXhwb3J0cyxlLHQsbixyKX1yZXR1cm4gbltvXS5leHBvcnRzfXZhciBpPXR5cGVvZiByZXF1aXJlPT1cImZ1bmN0aW9uXCImJnJlcXVpcmU7Zm9yKHZhciBvPTA7bzxyLmxlbmd0aDtvKyspcyhyW29dKTtyZXR1cm4gc30pIiwiLy8gPT1DbG9zdXJlQ29tcGlsZXI9PVxyXG4vLyBAb3V0cHV0X2ZpbGVfbmFtZSBmaXQtY3VydmUubWluLmpzXHJcbi8vIEBjb21waWxhdGlvbl9sZXZlbCBTSU1QTEVfT1BUSU1JWkFUSU9OU1xyXG4vLyA9PS9DbG9zdXJlQ29tcGlsZXI9PVxyXG5cclxuLyoqXHJcbiAqICBAcHJlc2VydmUgIEphdmFTY3JpcHQgaW1wbGVtZW50YXRpb24gb2ZcclxuICogIEFsZ29yaXRobSBmb3IgQXV0b21hdGljYWxseSBGaXR0aW5nIERpZ2l0aXplZCBDdXJ2ZXNcclxuICogIGJ5IFBoaWxpcCBKLiBTY2huZWlkZXJcclxuICogIFwiR3JhcGhpY3MgR2Vtc1wiLCBBY2FkZW1pYyBQcmVzcywgMTk5MFxyXG4gKlxyXG4gKiAgVGhlIE1JVCBMaWNlbnNlIChNSVQpXHJcbiAqXHJcbiAqICBodHRwczovL2dpdGh1Yi5jb20vc29zd293L2ZpdC1jdXJ2ZXNcclxuICovXHJcblxyXG4vKipcclxuICogRml0IG9uZSBvciBtb3JlIEJlemllciBjdXJ2ZXMgdG8gYSBzZXQgb2YgcG9pbnRzLlxyXG4gKlxyXG4gKiBAcGFyYW0ge0FycmF5PEFycmF5PE51bWJlcj4+fSBwb2ludHMgLSBBcnJheSBvZiBkaWdpdGl6ZWQgcG9pbnRzLCBlLmcuIFtbNSw1XSxbNSw1MF0sWzExMCwxNDBdLFsyMTAsMTYwXSxbMzIwLDExMF1dXHJcbiAqIEBwYXJhbSB7TnVtYmVyfSBtYXhFcnJvciAtIFRvbGVyYW5jZSwgc3F1YXJlZCBlcnJvciBiZXR3ZWVuIHBvaW50cyBhbmQgZml0dGVkIGN1cnZlXHJcbiAqIEByZXR1cm5zIHtBcnJheTxBcnJheTxBcnJheTxOdW1iZXI+Pj59IEFycmF5IG9mIEJlemllciBjdXJ2ZXMsIHdoZXJlIGVhY2ggZWxlbWVudCBpcyBbZmlyc3QtcG9pbnQsIGNvbnRyb2wtcG9pbnQtMSwgY29udHJvbC1wb2ludC0yLCBzZWNvbmQtcG9pbnRdIGFuZCBwb2ludHMgYXJlIFt4LCB5XVxyXG4gKi9cclxuZnVuY3Rpb24gZml0Q3VydmUocG9pbnRzLCBtYXhFcnJvciwgcHJvZ3Jlc3NDYWxsYmFjaykge1xyXG4gICAgaWYgKCFBcnJheS5pc0FycmF5KHBvaW50cykpIHtcclxuICAgICAgICB0aHJvdyBuZXcgVHlwZUVycm9yKFwiRmlyc3QgYXJndW1lbnQgc2hvdWxkIGJlIGFuIGFycmF5XCIpO1xyXG4gICAgfVxyXG4gICAgcG9pbnRzLmZvckVhY2goKHBvaW50KSA9PiB7XHJcbiAgICAgICAgaWYoIUFycmF5LmlzQXJyYXkocG9pbnQpIHx8IHBvaW50Lmxlbmd0aCAhPT0gMlxyXG4gICAgICAgIHx8IHR5cGVvZiBwb2ludFswXSAhPT0gJ251bWJlcicgfHwgdHlwZW9mIHBvaW50WzFdICE9PSAnbnVtYmVyJyl7XHJcbiAgICAgICAgICAgIHRocm93IEVycm9yKFwiRWFjaCBwb2ludCBzaG91bGQgYmUgYW4gYXJyYXkgb2YgdHdvIG51bWJlcnNcIilcclxuICAgICAgICB9XHJcbiAgICB9KTtcclxuICAgIC8vIFJlbW92ZSBkdXBsaWNhdGUgcG9pbnRzXHJcbiAgICBwb2ludHMgPSBwb2ludHMuZmlsdGVyKChwb2ludCwgaSkgPT5cclxuICAgICAgICBpID09PSAwIHx8ICEocG9pbnRbMF0gPT09IHBvaW50c1tpLTFdWzBdICYmIHBvaW50WzFdID09PSBwb2ludHNbaS0xXVsxXSlcclxuICAgICk7XHJcblxyXG4gICAgaWYgKHBvaW50cy5sZW5ndGggPCAyKSB7XHJcbiAgICAgICAgcmV0dXJuIFtdO1xyXG4gICAgfVxyXG5cclxuICAgIGNvbnN0IGxlbiA9IHBvaW50cy5sZW5ndGg7XHJcbiAgICBjb25zdCBsZWZ0VGFuZ2VudCA9IGNyZWF0ZVRhbmdlbnQocG9pbnRzWzFdLCBwb2ludHNbMF0pO1xyXG4gICAgY29uc3QgcmlnaHRUYW5nZW50ID0gY3JlYXRlVGFuZ2VudChwb2ludHNbbGVuIC0gMl0sIHBvaW50c1tsZW4gLSAxXSk7XHJcblxyXG4gICAgcmV0dXJuIGZpdEN1YmljKHBvaW50cywgbGVmdFRhbmdlbnQsIHJpZ2h0VGFuZ2VudCwgbWF4RXJyb3IsIHByb2dyZXNzQ2FsbGJhY2spO1xyXG59XHJcblxyXG4vKipcclxuICogRml0IGEgQmV6aWVyIGN1cnZlIHRvIGEgKHN1YilzZXQgb2YgZGlnaXRpemVkIHBvaW50cy5cclxuICogWW91ciBjb2RlIHNob3VsZCBub3QgY2FsbCB0aGlzIGZ1bmN0aW9uIGRpcmVjdGx5LiBVc2Uge0BsaW5rIGZpdEN1cnZlfSBpbnN0ZWFkLlxyXG4gKlxyXG4gKiBAcGFyYW0ge0FycmF5PEFycmF5PE51bWJlcj4+fSBwb2ludHMgLSBBcnJheSBvZiBkaWdpdGl6ZWQgcG9pbnRzLCBlLmcuIFtbNSw1XSxbNSw1MF0sWzExMCwxNDBdLFsyMTAsMTYwXSxbMzIwLDExMF1dXHJcbiAqIEBwYXJhbSB7QXJyYXk8TnVtYmVyPn0gbGVmdFRhbmdlbnQgLSBVbml0IHRhbmdlbnQgdmVjdG9yIGF0IHN0YXJ0IHBvaW50XHJcbiAqIEBwYXJhbSB7QXJyYXk8TnVtYmVyPn0gcmlnaHRUYW5nZW50IC0gVW5pdCB0YW5nZW50IHZlY3RvciBhdCBlbmQgcG9pbnRcclxuICogQHBhcmFtIHtOdW1iZXJ9IGVycm9yIC0gVG9sZXJhbmNlLCBzcXVhcmVkIGVycm9yIGJldHdlZW4gcG9pbnRzIGFuZCBmaXR0ZWQgY3VydmVcclxuICogQHJldHVybnMge0FycmF5PEFycmF5PEFycmF5PE51bWJlcj4+Pn0gQXJyYXkgb2YgQmV6aWVyIGN1cnZlcywgd2hlcmUgZWFjaCBlbGVtZW50IGlzIFtmaXJzdC1wb2ludCwgY29udHJvbC1wb2ludC0xLCBjb250cm9sLXBvaW50LTIsIHNlY29uZC1wb2ludF0gYW5kIHBvaW50cyBhcmUgW3gsIHldXHJcbiAqL1xyXG5mdW5jdGlvbiBmaXRDdWJpYyhwb2ludHMsIGxlZnRUYW5nZW50LCByaWdodFRhbmdlbnQsIGVycm9yLCBwcm9ncmVzc0NhbGxiYWNrKSB7XHJcbiAgICBjb25zdCBNYXhJdGVyYXRpb25zID0gMjA7ICAgLy9NYXggdGltZXMgdG8gdHJ5IGl0ZXJhdGluZyAodG8gZmluZCBhbiBhY2NlcHRhYmxlIGN1cnZlKVxyXG5cclxuICAgIHZhciBiZXpDdXJ2ZSwgICAgICAgICAgICAgICAvL0NvbnRyb2wgcG9pbnRzIG9mIGZpdHRlZCBCZXppZXIgY3VydmVcclxuICAgICAgICB1LCAgICAgICAgICAgICAgICAgICAgICAvL1BhcmFtZXRlciB2YWx1ZXMgZm9yIHBvaW50XHJcbiAgICAgICAgdVByaW1lLCAgICAgICAgICAgICAgICAgLy9JbXByb3ZlZCBwYXJhbWV0ZXIgdmFsdWVzXHJcbiAgICAgICAgbWF4RXJyb3IsIHByZXZFcnIsICAgICAgLy9NYXhpbXVtIGZpdHRpbmcgZXJyb3JcclxuICAgICAgICBzcGxpdFBvaW50LCBwcmV2U3BsaXQsICAvL1BvaW50IHRvIHNwbGl0IHBvaW50IHNldCBhdCBpZiB3ZSBuZWVkIG1vcmUgdGhhbiBvbmUgY3VydmVcclxuICAgICAgICBjZW50ZXJWZWN0b3IsIHRvQ2VudGVyVGFuZ2VudCwgZnJvbUNlbnRlclRhbmdlbnQsICAvL1VuaXQgdGFuZ2VudCB2ZWN0b3IocykgYXQgc3BsaXRQb2ludFxyXG4gICAgICAgIGJlemllcnMsICAgICAgICAgICAgICAgIC8vQXJyYXkgb2YgZml0dGVkIEJlemllciBjdXJ2ZXMgaWYgd2UgbmVlZCBtb3JlIHRoYW4gb25lIGN1cnZlXHJcbiAgICAgICAgZGlzdCwgaTtcclxuXHJcbiAgICAvL2NvbnNvbGUubG9nKCdmaXRDdWJpYywgJywgcG9pbnRzLmxlbmd0aCk7XHJcblxyXG4gICAgLy9Vc2UgaGV1cmlzdGljIGlmIHJlZ2lvbiBvbmx5IGhhcyB0d28gcG9pbnRzIGluIGl0XHJcbiAgICBpZiAocG9pbnRzLmxlbmd0aCA9PT0gMikge1xyXG4gICAgICAgIGRpc3QgPSBtYXRocy52ZWN0b3JMZW4obWF0aHMuc3VidHJhY3QocG9pbnRzWzBdLCBwb2ludHNbMV0pKSAvIDMuMDtcclxuICAgICAgICBiZXpDdXJ2ZSA9IFtcclxuICAgICAgICAgICAgcG9pbnRzWzBdLFxyXG4gICAgICAgICAgICBtYXRocy5hZGRBcnJheXMocG9pbnRzWzBdLCBtYXRocy5tdWxJdGVtcyhsZWZ0VGFuZ2VudCwgIGRpc3QpKSxcclxuICAgICAgICAgICAgbWF0aHMuYWRkQXJyYXlzKHBvaW50c1sxXSwgbWF0aHMubXVsSXRlbXMocmlnaHRUYW5nZW50LCBkaXN0KSksXHJcbiAgICAgICAgICAgIHBvaW50c1sxXVxyXG4gICAgICAgIF07XHJcbiAgICAgICAgcmV0dXJuIFtiZXpDdXJ2ZV07XHJcbiAgICB9XHJcblxyXG4gICAgLy9QYXJhbWV0ZXJpemUgcG9pbnRzLCBhbmQgYXR0ZW1wdCB0byBmaXQgY3VydmVcclxuICAgIHUgPSBjaG9yZExlbmd0aFBhcmFtZXRlcml6ZShwb2ludHMpO1xyXG4gICAgW2JlekN1cnZlLCBtYXhFcnJvciwgc3BsaXRQb2ludF0gPSBnZW5lcmF0ZUFuZFJlcG9ydChwb2ludHMsIHUsIHUsIGxlZnRUYW5nZW50LCByaWdodFRhbmdlbnQsIHByb2dyZXNzQ2FsbGJhY2spXHJcblxyXG4gICAgaWYgKG1heEVycm9yIDwgZXJyb3IpIHtcclxuICAgICAgICByZXR1cm4gW2JlekN1cnZlXTtcclxuICAgIH1cclxuICAgIC8vSWYgZXJyb3Igbm90IHRvbyBsYXJnZSwgdHJ5IHNvbWUgcmVwYXJhbWV0ZXJpemF0aW9uIGFuZCBpdGVyYXRpb25cclxuICAgIGlmIChtYXhFcnJvciA8IChlcnJvciplcnJvcikpIHtcclxuXHJcbiAgICAgICAgdVByaW1lID0gdTtcclxuICAgICAgICBwcmV2RXJyID0gbWF4RXJyb3I7XHJcbiAgICAgICAgcHJldlNwbGl0ID0gc3BsaXRQb2ludDtcclxuXHJcbiAgICAgICAgZm9yIChpID0gMDsgaSA8IE1heEl0ZXJhdGlvbnM7IGkrKykge1xyXG5cclxuICAgICAgICAgICAgdVByaW1lID0gcmVwYXJhbWV0ZXJpemUoYmV6Q3VydmUsIHBvaW50cywgdVByaW1lKTtcclxuICAgICAgICAgICAgW2JlekN1cnZlLCBtYXhFcnJvciwgc3BsaXRQb2ludF0gPSBnZW5lcmF0ZUFuZFJlcG9ydChwb2ludHMsIHUsIHVQcmltZSwgbGVmdFRhbmdlbnQsIHJpZ2h0VGFuZ2VudCwgcHJvZ3Jlc3NDYWxsYmFjayk7XHJcblxyXG4gICAgICAgICAgICBpZiAobWF4RXJyb3IgPCBlcnJvcikge1xyXG4gICAgICAgICAgICAgICAgcmV0dXJuIFtiZXpDdXJ2ZV07XHJcbiAgICAgICAgICAgIH1cclxuICAgICAgICAgICAgLy9JZiB0aGUgZGV2ZWxvcG1lbnQgb2YgdGhlIGZpdHRlZCBjdXJ2ZSBncmluZHMgdG8gYSBoYWx0LFxyXG4gICAgICAgICAgICAvL3dlIGFib3J0IHRoaXMgYXR0ZW1wdCAoYW5kIHRyeSBhIHNob3J0ZXIgY3VydmUpOlxyXG4gICAgICAgICAgICBlbHNlIGlmKHNwbGl0UG9pbnQgPT09IHByZXZTcGxpdCkge1xyXG4gICAgICAgICAgICAgICAgbGV0IGVyckNoYW5nZSA9IG1heEVycm9yL3ByZXZFcnI7XHJcbiAgICAgICAgICAgICAgICBpZigoZXJyQ2hhbmdlID4gLjk5OTkpICYmIChlcnJDaGFuZ2UgPCAxLjAwMDEpKSB7XHJcbiAgICAgICAgICAgICAgICAgICAgYnJlYWs7XHJcbiAgICAgICAgICAgICAgICB9XHJcbiAgICAgICAgICAgIH1cclxuXHJcbiAgICAgICAgICAgIHByZXZFcnIgPSBtYXhFcnJvcjtcclxuICAgICAgICAgICAgcHJldlNwbGl0ID0gc3BsaXRQb2ludDtcclxuICAgICAgICB9XHJcbiAgICB9XHJcblxyXG4gICAgLy9GaXR0aW5nIGZhaWxlZCAtLSBzcGxpdCBhdCBtYXggZXJyb3IgcG9pbnQgYW5kIGZpdCByZWN1cnNpdmVseVxyXG4gICAgYmV6aWVycyA9IFtdO1xyXG5cclxuICAgIC8vVG8gY3JlYXRlIGEgc21vb3RoIHRyYW5zaXRpb24gZnJvbSBvbmUgY3VydmUgc2VnbWVudCB0byB0aGUgbmV4dCxcclxuICAgIC8vd2UgY2FsY3VsYXRlIHRoZSB0YW5nZW50IG9mIHRoZSBwb2ludHMgZGlyZWN0bHkgYmVmb3JlIGFuZCBhZnRlciB0aGUgY2VudGVyLFxyXG4gICAgLy9hbmQgdXNlIHRoYXQgc2FtZSB0YW5nZW50IGJvdGggdG8gYW5kIGZyb20gdGhlIGNlbnRlciBwb2ludC5cclxuICAgIGNlbnRlclZlY3RvciA9IG1hdGhzLnN1YnRyYWN0KHBvaW50c1tzcGxpdFBvaW50IC0gMV0sIHBvaW50c1tzcGxpdFBvaW50ICsgMV0pO1xyXG4gICAgLy9Ib3dldmVyLCBzaG91bGQgdGhvc2UgdHdvIHBvaW50cyBiZSBlcXVhbCwgdGhlIG5vcm1hbCB0YW5nZW50IGNhbGN1bGF0aW9uIHdpbGwgZmFpbC5cclxuICAgIC8vSW5zdGVhZCwgd2UgY2FsY3VsYXRlIHRoZSB0YW5nZW50IGZyb20gdGhhdCBcImRvdWJsZS1wb2ludFwiIHRvIHRoZSBjZW50ZXIgcG9pbnQsIGFuZCByb3RhdGUgOTBkZWcuXHJcbiAgICBpZigoY2VudGVyVmVjdG9yWzBdID09PSAwKSAmJiAoY2VudGVyVmVjdG9yWzFdID09PSAwKSkge1xyXG4gICAgICAgIC8vdG9DZW50ZXJUYW5nZW50ID0gY3JlYXRlVGFuZ2VudChwb2ludHNbc3BsaXRQb2ludCAtIDFdLCBwb2ludHNbc3BsaXRQb2ludF0pO1xyXG4gICAgICAgIC8vZnJvbUNlbnRlclRhbmdlbnQgPSBjcmVhdGVUYW5nZW50KHBvaW50c1tzcGxpdFBvaW50ICsgMV0sIHBvaW50c1tzcGxpdFBvaW50XSk7XHJcblxyXG4gICAgICAgIC8vW3gseV0gLT4gWy15LHhdOiBodHRwOi8vc3RhY2tvdmVyZmxvdy5jb20vYS80NzgwMTQxLzE4Njk2NjBcclxuICAgICAgICBjZW50ZXJWZWN0b3IgPSBtYXRocy5zdWJ0cmFjdChwb2ludHNbc3BsaXRQb2ludCAtIDFdLCBwb2ludHNbc3BsaXRQb2ludF0pXHJcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAucmV2ZXJzZSgpO1xyXG4gICAgICAgIGNlbnRlclZlY3RvclswXSA9IC1jZW50ZXJWZWN0b3JbMF07XHJcbiAgICB9XHJcbiAgICB0b0NlbnRlclRhbmdlbnQgPSBtYXRocy5ub3JtYWxpemUoY2VudGVyVmVjdG9yKTtcclxuICAgIC8vVG8gYW5kIGZyb20gbmVlZCB0byBwb2ludCBpbiBvcHBvc2l0ZSBkaXJlY3Rpb25zOlxyXG4gICAgZnJvbUNlbnRlclRhbmdlbnQgPSBtYXRocy5tdWxJdGVtcyh0b0NlbnRlclRhbmdlbnQsIC0xKTtcclxuXHJcbiAgICAvKlxyXG4gICAgTm90ZTogQW4gYWx0ZXJuYXRpdmUgdG8gdGhpcyBcImRpdmlkZSBhbmQgY29ucXVlclwiIHJlY3Vyc2lvbiBjb3VsZCBiZSB0byBhbHdheXNcclxuICAgICAgICAgIGxldCBuZXcgY3VydmUgc2VnbWVudHMgc3RhcnQgYnkgdHJ5aW5nIHRvIGdvIGFsbCB0aGUgd2F5IHRvIHRoZSBlbmQsXHJcbiAgICAgICAgICBpbnN0ZWFkIG9mIG9ubHkgdG8gdGhlIGVuZCBvZiB0aGUgY3VycmVudCBzdWJkaXZpZGVkIHBvbHlsaW5lLlxyXG4gICAgICAgICAgVGhhdCBtaWdodCBsZXQgbWFueSBzZWdtZW50cyBmaXQgYSBmZXcgcG9pbnRzIG1vcmUsIHJlZHVjaW5nIHRoZSBudW1iZXIgb2YgdG90YWwgc2VnbWVudHMuXHJcblxyXG4gICAgICAgICAgSG93ZXZlciwgYSBmZXcgdGVzdHMgaGF2ZSBzaG93biB0aGF0IHRoZSBzZWdtZW50IHJlZHVjdGlvbiBpcyBpbnNpZ25pZmljYW50XHJcbiAgICAgICAgICAoMjQwIHB0cywgMTAwIGVycjogMjUgY3VydmVzIHZzIDI3IGN1cnZlcy4gMTQwIHB0cywgMTAwIGVycjogMTcgY3VydmVzIG9uIGJvdGgpLFxyXG4gICAgICAgICAgYW5kIHRoZSByZXN1bHRzIHRha2UgdHdpY2UgYXMgbWFueSBzdGVwcyBhbmQgbWlsbGlzZWNvbmRzIHRvIGZpbmlzaCxcclxuICAgICAgICAgIHdpdGhvdXQgbG9va2luZyBhbnkgYmV0dGVyIHRoYW4gd2hhdCB3ZSBhbHJlYWR5IGhhdmUuXHJcbiAgICAqL1xyXG4gICAgYmV6aWVycyA9IGJlemllcnMuY29uY2F0KGZpdEN1YmljKHBvaW50cy5zbGljZSgwLCBzcGxpdFBvaW50ICsgMSksIGxlZnRUYW5nZW50LCB0b0NlbnRlclRhbmdlbnQsICAgIGVycm9yLCBwcm9ncmVzc0NhbGxiYWNrKSk7XHJcbiAgICBiZXppZXJzID0gYmV6aWVycy5jb25jYXQoZml0Q3ViaWMocG9pbnRzLnNsaWNlKHNwbGl0UG9pbnQpLCAgICAgICAgZnJvbUNlbnRlclRhbmdlbnQsIHJpZ2h0VGFuZ2VudCwgZXJyb3IsIHByb2dyZXNzQ2FsbGJhY2spKTtcclxuICAgIHJldHVybiBiZXppZXJzO1xyXG59O1xyXG5cclxuZnVuY3Rpb24gZ2VuZXJhdGVBbmRSZXBvcnQocG9pbnRzLCBwYXJhbXNPcmlnLCBwYXJhbXNQcmltZSwgbGVmdFRhbmdlbnQsIHJpZ2h0VGFuZ2VudCwgcHJvZ3Jlc3NDYWxsYmFjaykge1xyXG4gICAgdmFyIGJlekN1cnZlLCBtYXhFcnJvciwgc3BsaXRQb2ludDtcclxuXHJcbiAgICBiZXpDdXJ2ZSA9IGdlbmVyYXRlQmV6aWVyKHBvaW50cywgcGFyYW1zUHJpbWUsIGxlZnRUYW5nZW50LCByaWdodFRhbmdlbnQsIHByb2dyZXNzQ2FsbGJhY2spO1xyXG4gICAgLy9GaW5kIG1heCBkZXZpYXRpb24gb2YgcG9pbnRzIHRvIGZpdHRlZCBjdXJ2ZS5cclxuICAgIC8vSGVyZSB3ZSBhbHdheXMgdXNlIHRoZSBvcmlnaW5hbCBwYXJhbWV0ZXJzIChmcm9tIGNob3JkTGVuZ3RoUGFyYW1ldGVyaXplKCkpLFxyXG4gICAgLy9iZWNhdXNlIHdlIG5lZWQgdG8gY29tcGFyZSB0aGUgY3VycmVudCBjdXJ2ZSB0byB0aGUgYWN0dWFsIHNvdXJjZSBwb2x5bGluZSxcclxuICAgIC8vYW5kIG5vdCB0aGUgY3VycmVudGx5IGl0ZXJhdGVkIHBhcmFtZXRlcnMgd2hpY2ggcmVwYXJhbWV0ZXJpemUoKSAmIGdlbmVyYXRlQmV6aWVyKCkgdXNlLFxyXG4gICAgLy9hcyB0aG9zZSBoYXZlIHByb2JhYmx5IGRyaWZ0ZWQgZmFyIGF3YXkgYW5kIG1heSBubyBsb25nZXIgYmUgaW4gYXNjZW5kaW5nIG9yZGVyLlxyXG4gICAgW21heEVycm9yLCBzcGxpdFBvaW50XSA9IGNvbXB1dGVNYXhFcnJvcihwb2ludHMsIGJlekN1cnZlLCBwYXJhbXNPcmlnKTtcclxuXHJcbiAgICBpZihwcm9ncmVzc0NhbGxiYWNrKSB7XHJcbiAgICAgICAgcHJvZ3Jlc3NDYWxsYmFjayh7XHJcbiAgICAgICAgICAgIGJlejogYmV6Q3VydmUsXHJcbiAgICAgICAgICAgIHBvaW50czogcG9pbnRzLFxyXG4gICAgICAgICAgICBwYXJhbXM6IHBhcmFtc09yaWcsXHJcbiAgICAgICAgICAgIG1heEVycjogbWF4RXJyb3IsXHJcbiAgICAgICAgICAgIG1heFBvaW50OiBzcGxpdFBvaW50LFxyXG4gICAgICAgIH0pO1xyXG4gICAgfVxyXG5cclxuICAgIHJldHVybiBbYmV6Q3VydmUsIG1heEVycm9yLCBzcGxpdFBvaW50XTtcclxufVxyXG5cclxuLyoqXHJcbiAqIFVzZSBsZWFzdC1zcXVhcmVzIG1ldGhvZCB0byBmaW5kIEJlemllciBjb250cm9sIHBvaW50cyBmb3IgcmVnaW9uLlxyXG4gKlxyXG4gKiBAcGFyYW0ge0FycmF5PEFycmF5PE51bWJlcj4+fSBwb2ludHMgLSBBcnJheSBvZiBkaWdpdGl6ZWQgcG9pbnRzXHJcbiAqIEBwYXJhbSB7QXJyYXk8TnVtYmVyPn0gcGFyYW1ldGVycyAtIFBhcmFtZXRlciB2YWx1ZXMgZm9yIHJlZ2lvblxyXG4gKiBAcGFyYW0ge0FycmF5PE51bWJlcj59IGxlZnRUYW5nZW50IC0gVW5pdCB0YW5nZW50IHZlY3RvciBhdCBzdGFydCBwb2ludFxyXG4gKiBAcGFyYW0ge0FycmF5PE51bWJlcj59IHJpZ2h0VGFuZ2VudCAtIFVuaXQgdGFuZ2VudCB2ZWN0b3IgYXQgZW5kIHBvaW50XHJcbiAqIEByZXR1cm5zIHtBcnJheTxBcnJheTxOdW1iZXI+Pn0gQXBwcm94aW1hdGVkIEJlemllciBjdXJ2ZTogW2ZpcnN0LXBvaW50LCBjb250cm9sLXBvaW50LTEsIGNvbnRyb2wtcG9pbnQtMiwgc2Vjb25kLXBvaW50XSB3aGVyZSBwb2ludHMgYXJlIFt4LCB5XVxyXG4gKi9cclxuZnVuY3Rpb24gZ2VuZXJhdGVCZXppZXIocG9pbnRzLCBwYXJhbWV0ZXJzLCBsZWZ0VGFuZ2VudCwgcmlnaHRUYW5nZW50KSB7XHJcbiAgICB2YXIgYmV6Q3VydmUsICAgICAgICAgICAgICAgICAgICAgICAvL0JlemllciBjdXJ2ZSBjdGwgcHRzXHJcbiAgICAgICAgQSwgYSwgICAgICAgICAgICAgICAgICAgICAgICAgICAvL1ByZWNvbXB1dGVkIHJocyBmb3IgZXFuXHJcbiAgICAgICAgQywgWCwgICAgICAgICAgICAgICAgICAgICAgICAgICAvL01hdHJpY2VzIEMgJiBYXHJcbiAgICAgICAgZGV0X0MwX0MxLCBkZXRfQzBfWCwgZGV0X1hfQzEsICAvL0RldGVybWluYW50cyBvZiBtYXRyaWNlc1xyXG4gICAgICAgIGFscGhhX2wsIGFscGhhX3IsICAgICAgICAgICAgICAgLy9BbHBoYSB2YWx1ZXMsIGxlZnQgYW5kIHJpZ2h0XHJcblxyXG4gICAgICAgIGVwc2lsb24sIHNlZ0xlbmd0aCxcclxuICAgICAgICBpLCBsZW4sIHRtcCwgdSwgdXgsXHJcbiAgICAgICAgZmlyc3RQb2ludCA9IHBvaW50c1swXSxcclxuICAgICAgICBsYXN0UG9pbnQgPSBwb2ludHNbcG9pbnRzLmxlbmd0aC0xXTtcclxuXHJcbiAgICBiZXpDdXJ2ZSA9IFtmaXJzdFBvaW50LCBudWxsLCBudWxsLCBsYXN0UG9pbnRdO1xyXG4gICAgLy9jb25zb2xlLmxvZygnZ2InLCBwYXJhbWV0ZXJzLmxlbmd0aCk7XHJcblxyXG4gICAgLy9Db21wdXRlIHRoZSBBJ3NcclxuICAgIEEgPSBtYXRocy56ZXJvc19YeDJ4MihwYXJhbWV0ZXJzLmxlbmd0aCk7XHJcbiAgICBmb3IgKGkgPSAwLCBsZW4gPSBwYXJhbWV0ZXJzLmxlbmd0aDsgaSA8IGxlbjsgaSsrKSB7XHJcbiAgICAgICAgdSA9IHBhcmFtZXRlcnNbaV07XHJcbiAgICAgICAgdXggPSAxIC0gdTtcclxuICAgICAgICBhID0gQVtpXTtcclxuXHJcbiAgICAgICAgYVswXSA9IG1hdGhzLm11bEl0ZW1zKGxlZnRUYW5nZW50LCAgMyAqIHUgICogKHV4KnV4KSk7XHJcbiAgICAgICAgYVsxXSA9IG1hdGhzLm11bEl0ZW1zKHJpZ2h0VGFuZ2VudCwgMyAqIHV4ICogKHUqdSkpO1xyXG4gICAgfVxyXG5cclxuICAgIC8vQ3JlYXRlIHRoZSBDIGFuZCBYIG1hdHJpY2VzXHJcbiAgICBDID0gW1swLDBdLCBbMCwwXV07XHJcbiAgICBYID0gWzAsMF07XHJcbiAgICBmb3IgKGkgPSAwLCBsZW4gPSBwb2ludHMubGVuZ3RoOyBpIDwgbGVuOyBpKyspIHtcclxuICAgICAgICB1ID0gcGFyYW1ldGVyc1tpXTtcclxuICAgICAgICBhID0gQVtpXTtcclxuXHJcbiAgICAgICAgQ1swXVswXSArPSBtYXRocy5kb3QoYVswXSwgYVswXSk7XHJcbiAgICAgICAgQ1swXVsxXSArPSBtYXRocy5kb3QoYVswXSwgYVsxXSk7XHJcbiAgICAgICAgQ1sxXVswXSArPSBtYXRocy5kb3QoYVswXSwgYVsxXSk7XHJcbiAgICAgICAgQ1sxXVsxXSArPSBtYXRocy5kb3QoYVsxXSwgYVsxXSk7XHJcblxyXG4gICAgICAgIHRtcCA9IG1hdGhzLnN1YnRyYWN0KHBvaW50c1tpXSwgYmV6aWVyLnEoW2ZpcnN0UG9pbnQsIGZpcnN0UG9pbnQsIGxhc3RQb2ludCwgbGFzdFBvaW50XSwgdSkpO1xyXG5cclxuICAgICAgICBYWzBdICs9IG1hdGhzLmRvdChhWzBdLCB0bXApO1xyXG4gICAgICAgIFhbMV0gKz0gbWF0aHMuZG90KGFbMV0sIHRtcCk7XHJcbiAgICB9XHJcblxyXG4gICAgLy9Db21wdXRlIHRoZSBkZXRlcm1pbmFudHMgb2YgQyBhbmQgWFxyXG4gICAgZGV0X0MwX0MxID0gKENbMF1bMF0gKiBDWzFdWzFdKSAtIChDWzFdWzBdICogQ1swXVsxXSk7XHJcbiAgICBkZXRfQzBfWCAgPSAoQ1swXVswXSAqIFhbMV0gICApIC0gKENbMV1bMF0gKiBYWzBdICAgKTtcclxuICAgIGRldF9YX0MxICA9IChYWzBdICAgICogQ1sxXVsxXSkgLSAoWFsxXSAgICAqIENbMF1bMV0pO1xyXG5cclxuICAgIC8vRmluYWxseSwgZGVyaXZlIGFscGhhIHZhbHVlc1xyXG4gICAgYWxwaGFfbCA9IGRldF9DMF9DMSA9PT0gMCA/IDAgOiBkZXRfWF9DMSAvIGRldF9DMF9DMTtcclxuICAgIGFscGhhX3IgPSBkZXRfQzBfQzEgPT09IDAgPyAwIDogZGV0X0MwX1ggLyBkZXRfQzBfQzE7XHJcblxyXG4gICAgLy9JZiBhbHBoYSBuZWdhdGl2ZSwgdXNlIHRoZSBXdS9CYXJza3kgaGV1cmlzdGljIChzZWUgdGV4dCkuXHJcbiAgICAvL0lmIGFscGhhIGlzIDAsIHlvdSBnZXQgY29pbmNpZGVudCBjb250cm9sIHBvaW50cyB0aGF0IGxlYWQgdG9cclxuICAgIC8vZGl2aWRlIGJ5IHplcm8gaW4gYW55IHN1YnNlcXVlbnQgTmV3dG9uUmFwaHNvblJvb3RGaW5kKCkgY2FsbC5cclxuICAgIHNlZ0xlbmd0aCA9IG1hdGhzLnZlY3RvckxlbihtYXRocy5zdWJ0cmFjdChmaXJzdFBvaW50LCBsYXN0UG9pbnQpKTtcclxuICAgIGVwc2lsb24gPSAxLjBlLTYgKiBzZWdMZW5ndGg7XHJcbiAgICBpZiAoYWxwaGFfbCA8IGVwc2lsb24gfHwgYWxwaGFfciA8IGVwc2lsb24pIHtcclxuICAgICAgICAvL0ZhbGwgYmFjayBvbiBzdGFuZGFyZCAocHJvYmFibHkgaW5hY2N1cmF0ZSkgZm9ybXVsYSwgYW5kIHN1YmRpdmlkZSBmdXJ0aGVyIGlmIG5lZWRlZC5cclxuICAgICAgICBiZXpDdXJ2ZVsxXSA9IG1hdGhzLmFkZEFycmF5cyhmaXJzdFBvaW50LCBtYXRocy5tdWxJdGVtcyhsZWZ0VGFuZ2VudCwgIHNlZ0xlbmd0aCAvIDMuMCkpO1xyXG4gICAgICAgIGJlekN1cnZlWzJdID0gbWF0aHMuYWRkQXJyYXlzKGxhc3RQb2ludCwgIG1hdGhzLm11bEl0ZW1zKHJpZ2h0VGFuZ2VudCwgc2VnTGVuZ3RoIC8gMy4wKSk7XHJcbiAgICB9IGVsc2Uge1xyXG4gICAgICAgIC8vRmlyc3QgYW5kIGxhc3QgY29udHJvbCBwb2ludHMgb2YgdGhlIEJlemllciBjdXJ2ZSBhcmVcclxuICAgICAgICAvL3Bvc2l0aW9uZWQgZXhhY3RseSBhdCB0aGUgZmlyc3QgYW5kIGxhc3QgZGF0YSBwb2ludHNcclxuICAgICAgICAvL0NvbnRyb2wgcG9pbnRzIDEgYW5kIDIgYXJlIHBvc2l0aW9uZWQgYW4gYWxwaGEgZGlzdGFuY2Ugb3V0XHJcbiAgICAgICAgLy9vbiB0aGUgdGFuZ2VudCB2ZWN0b3JzLCBsZWZ0IGFuZCByaWdodCwgcmVzcGVjdGl2ZWx5XHJcbiAgICAgICAgYmV6Q3VydmVbMV0gPSBtYXRocy5hZGRBcnJheXMoZmlyc3RQb2ludCwgbWF0aHMubXVsSXRlbXMobGVmdFRhbmdlbnQsICBhbHBoYV9sKSk7XHJcbiAgICAgICAgYmV6Q3VydmVbMl0gPSBtYXRocy5hZGRBcnJheXMobGFzdFBvaW50LCAgbWF0aHMubXVsSXRlbXMocmlnaHRUYW5nZW50LCBhbHBoYV9yKSk7XHJcbiAgICB9XHJcblxyXG4gICAgcmV0dXJuIGJlekN1cnZlO1xyXG59O1xyXG5cclxuLyoqXHJcbiAqIEdpdmVuIHNldCBvZiBwb2ludHMgYW5kIHRoZWlyIHBhcmFtZXRlcml6YXRpb24sIHRyeSB0byBmaW5kIGEgYmV0dGVyIHBhcmFtZXRlcml6YXRpb24uXHJcbiAqXHJcbiAqIEBwYXJhbSB7QXJyYXk8QXJyYXk8TnVtYmVyPj59IGJlemllciAtIEN1cnJlbnQgZml0dGVkIGN1cnZlXHJcbiAqIEBwYXJhbSB7QXJyYXk8QXJyYXk8TnVtYmVyPj59IHBvaW50cyAtIEFycmF5IG9mIGRpZ2l0aXplZCBwb2ludHNcclxuICogQHBhcmFtIHtBcnJheTxOdW1iZXI+fSBwYXJhbWV0ZXJzIC0gQ3VycmVudCBwYXJhbWV0ZXIgdmFsdWVzXHJcbiAqIEByZXR1cm5zIHtBcnJheTxOdW1iZXI+fSBOZXcgcGFyYW1ldGVyIHZhbHVlc1xyXG4gKi9cclxuZnVuY3Rpb24gcmVwYXJhbWV0ZXJpemUoYmV6aWVyLCBwb2ludHMsIHBhcmFtZXRlcnMpIHtcclxuICAgIC8qXHJcbiAgICB2YXIgaiwgbGVuLCBwb2ludCwgcmVzdWx0cywgdTtcclxuICAgIHJlc3VsdHMgPSBbXTtcclxuICAgIGZvciAoaiA9IDAsIGxlbiA9IHBvaW50cy5sZW5ndGg7IGogPCBsZW47IGorKykge1xyXG4gICAgICAgIHBvaW50ID0gcG9pbnRzW2pdLCB1ID0gcGFyYW1ldGVyc1tqXTtcclxuXHJcbiAgICAgICAgcmVzdWx0cy5wdXNoKG5ld3RvblJhcGhzb25Sb290RmluZChiZXppZXIsIHBvaW50LCB1KSk7XHJcbiAgICB9XHJcbiAgICByZXR1cm4gcmVzdWx0cztcclxuICAgIC8vKi9cclxuICAgIHJldHVybiBwYXJhbWV0ZXJzLm1hcCgocCwgaSkgPT4gbmV3dG9uUmFwaHNvblJvb3RGaW5kKGJlemllciwgcG9pbnRzW2ldLCBwKSk7XHJcbn07XHJcblxyXG4vKipcclxuICogVXNlIE5ld3Rvbi1SYXBoc29uIGl0ZXJhdGlvbiB0byBmaW5kIGJldHRlciByb290LlxyXG4gKlxyXG4gKiBAcGFyYW0ge0FycmF5PEFycmF5PE51bWJlcj4+fSBiZXogLSBDdXJyZW50IGZpdHRlZCBjdXJ2ZVxyXG4gKiBAcGFyYW0ge0FycmF5PE51bWJlcj59IHBvaW50IC0gRGlnaXRpemVkIHBvaW50XHJcbiAqIEBwYXJhbSB7TnVtYmVyfSB1IC0gUGFyYW1ldGVyIHZhbHVlIGZvciBcIlBcIlxyXG4gKiBAcmV0dXJucyB7TnVtYmVyfSBOZXcgdVxyXG4gKi9cclxuZnVuY3Rpb24gbmV3dG9uUmFwaHNvblJvb3RGaW5kKGJleiwgcG9pbnQsIHUpIHtcclxuICAgIC8qXHJcbiAgICAgICAgTmV3dG9uJ3Mgcm9vdCBmaW5kaW5nIGFsZ29yaXRobSBjYWxjdWxhdGVzIGYoeCk9MCBieSByZWl0ZXJhdGluZ1xyXG4gICAgICAgIHhfbisxID0geF9uIC0gZih4X24pL2YnKHhfbilcclxuICAgICAgICBXZSBhcmUgdHJ5aW5nIHRvIGZpbmQgY3VydmUgcGFyYW1ldGVyIHUgZm9yIHNvbWUgcG9pbnQgcCB0aGF0IG1pbmltaXplc1xyXG4gICAgICAgIHRoZSBkaXN0YW5jZSBmcm9tIHRoYXQgcG9pbnQgdG8gdGhlIGN1cnZlLiBEaXN0YW5jZSBwb2ludCB0byBjdXJ2ZSBpcyBkPXEodSktcC5cclxuICAgICAgICBBdCBtaW5pbXVtIGRpc3RhbmNlIHRoZSBwb2ludCBpcyBwZXJwZW5kaWN1bGFyIHRvIHRoZSBjdXJ2ZS5cclxuICAgICAgICBXZSBhcmUgc29sdmluZ1xyXG4gICAgICAgIGYgPSBxKHUpLXAgKiBxJyh1KSA9IDBcclxuICAgICAgICB3aXRoXHJcbiAgICAgICAgZicgPSBxJyh1KSAqIHEnKHUpICsgcSh1KS1wICogcScnKHUpXHJcbiAgICAgICAgZ2l2ZXNcclxuICAgICAgICB1X24rMSA9IHVfbiAtIHxxKHVfbiktcCAqIHEnKHVfbil8IC8gfHEnKHVfbikqKjIgKyBxKHVfbiktcCAqIHEnJyh1X24pfFxyXG4gICAgKi9cclxuXHJcbiAgICB2YXIgZCA9IG1hdGhzLnN1YnRyYWN0KGJlemllci5xKGJleiwgdSksIHBvaW50KSxcclxuICAgICAgICBxcHJpbWUgPSBiZXppZXIucXByaW1lKGJleiwgdSksXHJcbiAgICAgICAgbnVtZXJhdG9yID0gLypzdW0oKi9tYXRocy5tdWxNYXRyaXgoZCwgcXByaW1lKS8qKSovLFxyXG4gICAgICAgIGRlbm9taW5hdG9yID0gbWF0aHMuc3VtKG1hdGhzLmFkZEl0ZW1zKCBtYXRocy5zcXVhcmVJdGVtcyhxcHJpbWUpLCBtYXRocy5tdWxNYXRyaXgoZCwgYmV6aWVyLnFwcmltZXByaW1lKGJleiwgdSkpICkpO1xyXG5cclxuICAgIGlmIChkZW5vbWluYXRvciA9PT0gMCkge1xyXG4gICAgICAgIHJldHVybiB1O1xyXG4gICAgfSBlbHNlIHtcclxuICAgICAgICByZXR1cm4gdSAtIChudW1lcmF0b3IvZGVub21pbmF0b3IpO1xyXG4gICAgfVxyXG59O1xyXG5cclxuLyoqXHJcbiAqIEFzc2lnbiBwYXJhbWV0ZXIgdmFsdWVzIHRvIGRpZ2l0aXplZCBwb2ludHMgdXNpbmcgcmVsYXRpdmUgZGlzdGFuY2VzIGJldHdlZW4gcG9pbnRzLlxyXG4gKlxyXG4gKiBAcGFyYW0ge0FycmF5PEFycmF5PE51bWJlcj4+fSBwb2ludHMgLSBBcnJheSBvZiBkaWdpdGl6ZWQgcG9pbnRzXHJcbiAqIEByZXR1cm5zIHtBcnJheTxOdW1iZXI+fSBQYXJhbWV0ZXIgdmFsdWVzXHJcbiAqL1xyXG5mdW5jdGlvbiBjaG9yZExlbmd0aFBhcmFtZXRlcml6ZShwb2ludHMpIHtcclxuICAgIHZhciB1ID0gW10sIGN1cnJVLCBwcmV2VSwgcHJldlA7XHJcblxyXG4gICAgcG9pbnRzLmZvckVhY2goKHAsIGkpID0+IHtcclxuICAgICAgICBjdXJyVSA9IGkgPyBwcmV2VSArIG1hdGhzLnZlY3RvckxlbihtYXRocy5zdWJ0cmFjdChwLCBwcmV2UCkpXHJcbiAgICAgICAgICAgICAgICAgIDogMDtcclxuICAgICAgICB1LnB1c2goY3VyclUpO1xyXG5cclxuICAgICAgICBwcmV2VSA9IGN1cnJVO1xyXG4gICAgICAgIHByZXZQID0gcDtcclxuICAgIH0pXHJcbiAgICB1ID0gdS5tYXAoeCA9PiB4L3ByZXZVKTtcclxuXHJcbiAgICByZXR1cm4gdTtcclxufTtcclxuXHJcbi8qKlxyXG4gKiBGaW5kIHRoZSBtYXhpbXVtIHNxdWFyZWQgZGlzdGFuY2Ugb2YgZGlnaXRpemVkIHBvaW50cyB0byBmaXR0ZWQgY3VydmUuXHJcbiAqXHJcbiAqIEBwYXJhbSB7QXJyYXk8QXJyYXk8TnVtYmVyPj59IHBvaW50cyAtIEFycmF5IG9mIGRpZ2l0aXplZCBwb2ludHNcclxuICogQHBhcmFtIHtBcnJheTxBcnJheTxOdW1iZXI+Pn0gYmV6IC0gRml0dGVkIGN1cnZlXHJcbiAqIEBwYXJhbSB7QXJyYXk8TnVtYmVyPn0gcGFyYW1ldGVycyAtIFBhcmFtZXRlcml6YXRpb24gb2YgcG9pbnRzXHJcbiAqIEByZXR1cm5zIHtBcnJheTxOdW1iZXI+fSBNYXhpbXVtIGVycm9yIChzcXVhcmVkKSBhbmQgcG9pbnQgb2YgbWF4IGVycm9yXHJcbiAqL1xyXG5mdW5jdGlvbiBjb21wdXRlTWF4RXJyb3IocG9pbnRzLCBiZXosIHBhcmFtZXRlcnMpIHtcclxuICAgIHZhciBkaXN0LCAgICAgICAvL0N1cnJlbnQgZXJyb3JcclxuICAgICAgICBtYXhEaXN0LCAgICAvL01heGltdW0gZXJyb3JcclxuICAgICAgICBzcGxpdFBvaW50LCAvL1BvaW50IG9mIG1heGltdW0gZXJyb3JcclxuICAgICAgICB2LCAgICAgICAgICAvL1ZlY3RvciBmcm9tIHBvaW50IHRvIGN1cnZlXHJcbiAgICAgICAgaSwgY291bnQsIHBvaW50LCB0O1xyXG5cclxuICAgIG1heERpc3QgPSAwO1xyXG4gICAgc3BsaXRQb2ludCA9IHBvaW50cy5sZW5ndGggLyAyO1xyXG5cclxuICAgIGNvbnN0IHRfZGlzdE1hcCA9IG1hcFR0b1JlbGF0aXZlRGlzdGFuY2VzKGJleiwgMTApO1xyXG5cclxuICAgIGZvciAoaSA9IDAsIGNvdW50ID0gcG9pbnRzLmxlbmd0aDsgaSA8IGNvdW50OyBpKyspIHtcclxuICAgICAgICBwb2ludCA9IHBvaW50c1tpXTtcclxuICAgICAgICAvL0ZpbmQgJ3QnIGZvciBhIHBvaW50IG9uIHRoZSBiZXogY3VydmUgdGhhdCdzIGFzIGNsb3NlIHRvICdwb2ludCcgYXMgcG9zc2libGU6XHJcbiAgICAgICAgdCA9IGZpbmRfdChiZXosIHBhcmFtZXRlcnNbaV0sIHRfZGlzdE1hcCwgMTApO1xyXG5cclxuICAgICAgICB2ID0gbWF0aHMuc3VidHJhY3QoYmV6aWVyLnEoYmV6LCB0KSwgcG9pbnQpO1xyXG4gICAgICAgIGRpc3QgPSB2WzBdKnZbMF0gKyB2WzFdKnZbMV07XHJcblxyXG4gICAgICAgIGlmIChkaXN0ID4gbWF4RGlzdCkge1xyXG4gICAgICAgICAgICBtYXhEaXN0ID0gZGlzdDtcclxuICAgICAgICAgICAgc3BsaXRQb2ludCA9IGk7XHJcbiAgICAgICAgfVxyXG4gICAgfVxyXG5cclxuICAgIHJldHVybiBbbWF4RGlzdCwgc3BsaXRQb2ludF07XHJcbn07XHJcblxyXG4vL1NhbXBsZSAndCdzIGFuZCBtYXAgdGhlbSB0byByZWxhdGl2ZSBkaXN0YW5jZXMgYWxvbmcgdGhlIGN1cnZlOlxyXG52YXIgbWFwVHRvUmVsYXRpdmVEaXN0YW5jZXMgPSBmdW5jdGlvbiAoYmV6LCBCX3BhcnRzKSB7XHJcbiAgICB2YXIgQl90X2N1cnI7XHJcbiAgICB2YXIgQl90X2Rpc3QgPSBbMF07XHJcbiAgICB2YXIgQl90X3ByZXYgPSBiZXpbMF07XHJcbiAgICB2YXIgc3VtTGVuID0gMDtcclxuXHJcbiAgICBmb3IgKHZhciBpPTE7IGk8PUJfcGFydHM7IGkrKykge1xyXG4gICAgICBCX3RfY3VyciA9IGJlemllci5xKGJleiwgaS9CX3BhcnRzKTtcclxuXHJcbiAgICAgIHN1bUxlbiArPSBtYXRocy52ZWN0b3JMZW4obWF0aHMuc3VidHJhY3QoQl90X2N1cnIsIEJfdF9wcmV2KSk7XHJcblxyXG4gICAgICBCX3RfZGlzdC5wdXNoKHN1bUxlbik7XHJcbiAgICAgIEJfdF9wcmV2ID0gQl90X2N1cnI7XHJcbiAgICB9XHJcblxyXG4gICAgLy9Ob3JtYWxpemUgQl9sZW5ndGggdG8gdGhlIHNhbWUgaW50ZXJ2YWwgYXMgdGhlIHBhcmFtZXRlciBkaXN0YW5jZXM7IDAgdG8gMTpcclxuICAgIEJfdF9kaXN0ID0gQl90X2Rpc3QubWFwKHggPT4geC9zdW1MZW4pO1xyXG4gICAgcmV0dXJuIEJfdF9kaXN0O1xyXG59O1xyXG5cclxuZnVuY3Rpb24gZmluZF90KGJleiwgcGFyYW0sIHRfZGlzdE1hcCwgQl9wYXJ0cykge1xyXG4gICAgaWYocGFyYW0gPCAwKSB7IHJldHVybiAwOyB9XHJcbiAgICBpZihwYXJhbSA+IDEpIHsgcmV0dXJuIDE7IH1cclxuXHJcbiAgICAvKlxyXG4gICAgICAgICdwYXJhbScgaXMgYSB2YWx1ZSBiZXR3ZWVuIDAgYW5kIDEgdGVsbGluZyB1cyB0aGUgcmVsYXRpdmUgcG9zaXRpb25cclxuICAgICAgICBvZiBhIHBvaW50IG9uIHRoZSBzb3VyY2UgcG9seWxpbmUgKGxpbmVhcmx5IGZyb20gdGhlIHN0YXJ0ICgwKSB0byB0aGUgZW5kICgxKSkuXHJcbiAgICAgICAgVG8gc2VlIGlmIGEgZ2l2ZW4gY3VydmUgLSAnYmV6JyAtIGlzIGEgY2xvc2UgYXBwcm94aW1hdGlvbiBvZiB0aGUgcG9seWxpbmUsXHJcbiAgICAgICAgd2UgY29tcGFyZSBzdWNoIGEgcG9seS1wb2ludCB0byB0aGUgcG9pbnQgb24gdGhlIGN1cnZlIHRoYXQncyB0aGUgc2FtZVxyXG4gICAgICAgIHJlbGF0aXZlIGRpc3RhbmNlIGFsb25nIHRoZSBjdXJ2ZSdzIGxlbmd0aC5cclxuXHJcbiAgICAgICAgQnV0IGZpbmRpbmcgdGhhdCBjdXJ2ZS1wb2ludCB0YWtlcyBhIGxpdHRsZSB3b3JrOlxyXG4gICAgICAgIFRoZXJlIGlzIGEgZnVuY3Rpb24gXCJCKHQpXCIgdG8gZmluZCBwb2ludHMgYWxvbmcgYSBjdXJ2ZSBmcm9tIHRoZSBwYXJhbWV0cmljIHBhcmFtZXRlciAndCdcclxuICAgICAgICAoYWxzbyByZWxhdGl2ZSBmcm9tIDAgdG8gMTogaHR0cDovL3N0YWNrb3ZlcmZsb3cuY29tL2EvMzI4NDE3NjQvMTg2OTY2MFxyXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBodHRwOi8vcG9tYXguZ2l0aHViLmlvL2JlemllcmluZm8vI2V4cGxhbmF0aW9uKSxcclxuICAgICAgICBidXQgJ3QnIGlzbid0IGxpbmVhciBieSBsZW5ndGggKGh0dHA6Ly9nYW1lZGV2LnN0YWNrZXhjaGFuZ2UuY29tL3F1ZXN0aW9ucy8xMDUyMzApLlxyXG5cclxuICAgICAgICBTbywgd2Ugc2FtcGxlIHNvbWUgcG9pbnRzIGFsb25nIHRoZSBjdXJ2ZSB1c2luZyBhIGhhbmRmdWwgb2YgdmFsdWVzIGZvciAndCcuXHJcbiAgICAgICAgVGhlbiwgd2UgY2FsY3VsYXRlIHRoZSBsZW5ndGggYmV0d2VlbiB0aG9zZSBzYW1wbGVzIHZpYSBwbGFpbiBldWNsaWRlYW4gZGlzdGFuY2U7XHJcbiAgICAgICAgQih0KSBjb25jZW50cmF0ZXMgdGhlIHBvaW50cyBhcm91bmQgc2hhcnAgdHVybnMsIHNvIHRoaXMgc2hvdWxkIGdpdmUgdXMgYSBnb29kLWVub3VnaCBvdXRsaW5lIG9mIHRoZSBjdXJ2ZS5cclxuICAgICAgICBUaHVzLCBmb3IgYSBnaXZlbiByZWxhdGl2ZSBkaXN0YW5jZSAoJ3BhcmFtJyksIHdlIGNhbiBub3cgZmluZCBhbiB1cHBlciBhbmQgbG93ZXIgdmFsdWVcclxuICAgICAgICBmb3IgdGhlIGNvcnJlc3BvbmRpbmcgJ3QnIGJ5IHNlYXJjaGluZyB0aHJvdWdoIHRob3NlIHNhbXBsZWQgZGlzdGFuY2VzLlxyXG4gICAgICAgIEZpbmFsbHksIHdlIGp1c3QgdXNlIGxpbmVhciBpbnRlcnBvbGF0aW9uIHRvIGZpbmQgYSBiZXR0ZXIgdmFsdWUgZm9yIHRoZSBleGFjdCAndCcuXHJcblxyXG4gICAgICAgIE1vcmUgaW5mbzpcclxuICAgICAgICAgICAgaHR0cDovL2dhbWVkZXYuc3RhY2tleGNoYW5nZS5jb20vcXVlc3Rpb25zLzEwNTIzMC9wb2ludHMtZXZlbmx5LXNwYWNlZC1hbG9uZy1hLWJlemllci1jdXJ2ZVxyXG4gICAgICAgICAgICBodHRwOi8vc3RhY2tvdmVyZmxvdy5jb20vcXVlc3Rpb25zLzI5NDM4Mzk4L2NoZWFwLXdheS1vZi1jYWxjdWxhdGluZy1jdWJpYy1iZXppZXItbGVuZ3RoXHJcbiAgICAgICAgICAgIGh0dHA6Ly9zdGV2ZS5ob2xsYXNjaC5uZXQvY2dpbmRleC9jdXJ2ZXMvY2JlemFyY2xlbi5odG1sXHJcbiAgICAgICAgICAgIGh0dHBzOi8vZ2l0aHViLmNvbS9yZXR1eHgvdGlueXNwbGluZVxyXG4gICAgKi9cclxuICAgIHZhciBsZW5NYXgsIGxlbk1pbiwgdE1heCwgdE1pbiwgdDtcclxuXHJcbiAgICAvL0ZpbmQgdGhlIHR3byB0LXMgdGhhdCB0aGUgY3VycmVudCBwYXJhbSBkaXN0YW5jZSBsaWVzIGJldHdlZW4sXHJcbiAgICAvL2FuZCB0aGVuIGludGVycG9sYXRlIGEgc29tZXdoYXQgYWNjdXJhdGUgdmFsdWUgZm9yIHRoZSBleGFjdCB0OlxyXG4gICAgZm9yKHZhciBpID0gMTsgaSA8PSBCX3BhcnRzOyBpKyspIHtcclxuXHJcbiAgICAgICAgaWYocGFyYW0gPD0gdF9kaXN0TWFwW2ldKSB7XHJcbiAgICAgICAgICAgIHRNaW4gICA9IChpLTEpIC8gQl9wYXJ0cztcclxuICAgICAgICAgICAgdE1heCAgID0gaSAvIEJfcGFydHM7XHJcbiAgICAgICAgICAgIGxlbk1pbiA9IHRfZGlzdE1hcFtpLTFdO1xyXG4gICAgICAgICAgICBsZW5NYXggPSB0X2Rpc3RNYXBbaV07XHJcblxyXG4gICAgICAgICAgICB0ID0gKHBhcmFtLWxlbk1pbikvKGxlbk1heC1sZW5NaW4pICogKHRNYXgtdE1pbikgKyB0TWluO1xyXG4gICAgICAgICAgICBicmVhaztcclxuICAgICAgICB9XHJcbiAgICB9XHJcbiAgICByZXR1cm4gdDtcclxufVxyXG5cclxuLyoqXHJcbiAqIENyZWF0ZXMgYSB2ZWN0b3Igb2YgbGVuZ3RoIDEgd2hpY2ggc2hvd3MgdGhlIGRpcmVjdGlvbiBmcm9tIEIgdG8gQVxyXG4gKi9cclxuZnVuY3Rpb24gY3JlYXRlVGFuZ2VudChwb2ludEEsIHBvaW50Qikge1xyXG4gICAgcmV0dXJuIG1hdGhzLm5vcm1hbGl6ZShtYXRocy5zdWJ0cmFjdChwb2ludEEsIHBvaW50QikpO1xyXG59XHJcblxyXG4vKlxyXG4gICAgU2ltcGxpZmllZCB2ZXJzaW9ucyBvZiB3aGF0IHdlIG5lZWQgZnJvbSBtYXRoLmpzXHJcbiAgICBPcHRpbWl6ZWQgZm9yIG91ciBpbnB1dCwgd2hpY2ggaXMgb25seSBudW1iZXJzIGFuZCAxeDIgYXJyYXlzIChpLmUuIFt4LCB5XSBjb29yZGluYXRlcykuXHJcbiovXHJcbmNsYXNzIG1hdGhzIHtcclxuICAgIC8vemVyb3MgPSBsb2dBbmRSdW4obWF0aC56ZXJvcyk7XHJcbiAgICBzdGF0aWMgemVyb3NfWHgyeDIoeCkge1xyXG4gICAgICAgIHZhciB6cyA9IFtdO1xyXG4gICAgICAgIHdoaWxlKHgtLSkgeyB6cy5wdXNoKFswLDBdKTsgfVxyXG4gICAgICAgIHJldHVybiB6c1xyXG4gICAgfVxyXG5cclxuICAgIC8vbXVsdGlwbHkgPSBsb2dBbmRSdW4obWF0aC5tdWx0aXBseSk7XHJcbiAgICBzdGF0aWMgbXVsSXRlbXMoaXRlbXMsIG11bHRpcGxpZXIpIHtcclxuICAgICAgICAvL3JldHVybiBpdGVtcy5tYXAoeCA9PiB4Km11bHRpcGxpZXIpO1xyXG4gICAgICAgIHJldHVybiBbaXRlbXNbMF0qbXVsdGlwbGllciwgaXRlbXNbMV0qbXVsdGlwbGllcl07XHJcbiAgICB9XHJcbiAgICBzdGF0aWMgbXVsTWF0cml4KG0xLCBtMikge1xyXG4gICAgICAgIC8vaHR0cHM6Ly9lbi53aWtpcGVkaWEub3JnL3dpa2kvTWF0cml4X211bHRpcGxpY2F0aW9uI01hdHJpeF9wcm9kdWN0Xy4yOHR3b19tYXRyaWNlcy4yOVxyXG4gICAgICAgIC8vU2ltcGxpZmllZCB0byBvbmx5IGhhbmRsZSAxLWRpbWVuc2lvbmFsIG1hdHJpY2VzIChpLmUuIGFycmF5cykgb2YgZXF1YWwgbGVuZ3RoOlxyXG4gICAgICAgIC8vICByZXR1cm4gbTEucmVkdWNlKChzdW0seDEsaSkgPT4gc3VtICsgKHgxKm0yW2ldKSxcclxuICAgICAgICAvLyAgICAgICAgICAgICAgICAgICAwKTtcclxuICAgICAgICByZXR1cm4gKG0xWzBdKm0yWzBdKSArIChtMVsxXSptMlsxXSk7XHJcbiAgICB9XHJcblxyXG4gICAgLy9Pbmx5IHVzZWQgdG8gc3VicmFjdCB0byBwb2ludHMgKG9yIGF0IGxlYXN0IGFycmF5cyk6XHJcbiAgICAvLyAgc3VidHJhY3QgPSBsb2dBbmRSdW4obWF0aC5zdWJ0cmFjdCk7XHJcbiAgICBzdGF0aWMgc3VidHJhY3QoYXJyMSwgYXJyMikge1xyXG4gICAgICAgIC8vcmV0dXJuIGFycjEubWFwKCh4MSwgaSkgPT4geDEgLSBhcnIyW2ldKTtcclxuICAgICAgICByZXR1cm4gW2FycjFbMF0tYXJyMlswXSwgYXJyMVsxXS1hcnIyWzFdXTtcclxuICAgIH1cclxuXHJcbiAgICAvL2FkZCA9IGxvZ0FuZFJ1bihtYXRoLmFkZCk7XHJcbiAgICBzdGF0aWMgYWRkQXJyYXlzKGFycjEsIGFycjIpIHtcclxuICAgICAgICAvL3JldHVybiBhcnIxLm1hcCgoeDEsIGkpID0+IHgxICsgYXJyMltpXSk7XHJcbiAgICAgICAgcmV0dXJuIFthcnIxWzBdK2FycjJbMF0sIGFycjFbMV0rYXJyMlsxXV07XHJcbiAgICB9XHJcbiAgICBzdGF0aWMgYWRkSXRlbXMoaXRlbXMsIGFkZGl0aW9uKSB7XHJcbiAgICAgICAgLy9yZXR1cm4gaXRlbXMubWFwKHggPT4geCthZGRpdGlvbik7XHJcbiAgICAgICAgcmV0dXJuIFtpdGVtc1swXSthZGRpdGlvbiwgaXRlbXNbMV0rYWRkaXRpb25dO1xyXG4gICAgfVxyXG5cclxuICAgIC8vdmFyIHN1bSA9IGxvZ0FuZFJ1bihtYXRoLnN1bSk7XHJcbiAgICBzdGF0aWMgc3VtKGl0ZW1zKSB7XHJcbiAgICAgICAgcmV0dXJuIGl0ZW1zLnJlZHVjZSgoc3VtLHgpID0+IHN1bSArIHgpO1xyXG4gICAgfVxyXG5cclxuICAgIC8vY2hhaW4gPSBtYXRoLmNoYWluO1xyXG5cclxuICAgIC8vT25seSB1c2VkIG9uIHR3byBhcnJheXMuIFRoZSBkb3QgcHJvZHVjdCBpcyBlcXVhbCB0byB0aGUgbWF0cml4IHByb2R1Y3QgaW4gdGhpcyBjYXNlOlxyXG4gICAgLy8gIGRvdCA9IGxvZ0FuZFJ1bihtYXRoLmRvdCk7XHJcbiAgICBzdGF0aWMgZG90KG0xLCBtMikge1xyXG4gICAgICAgIHJldHVybiBtYXRocy5tdWxNYXRyaXgobTEsIG0yKTtcclxuICAgIH1cclxuXHJcbiAgICAvL2h0dHBzOi8vZW4ud2lraXBlZGlhLm9yZy93aWtpL05vcm1fKG1hdGhlbWF0aWNzKSNFdWNsaWRlYW5fbm9ybVxyXG4gICAgLy8gIHZhciBub3JtID0gbG9nQW5kUnVuKG1hdGgubm9ybSk7XHJcbiAgICBzdGF0aWMgdmVjdG9yTGVuKHYpIHtcclxuICAgICAgICB2YXIgYSA9IHZbMF0sIGIgPSB2WzFdO1xyXG4gICAgICAgIHJldHVybiBNYXRoLnNxcnQoYSphICsgYipiKTtcclxuICAgIH1cclxuXHJcbiAgICAvL21hdGguZGl2aWRlID0gbG9nQW5kUnVuKG1hdGguZGl2aWRlKTtcclxuICAgIHN0YXRpYyBkaXZJdGVtcyhpdGVtcywgZGl2aXNvcikge1xyXG4gICAgICAgIC8vcmV0dXJuIGl0ZW1zLm1hcCh4ID0+IHgvZGl2aXNvcik7XHJcbiAgICAgICAgcmV0dXJuIFtpdGVtc1swXS9kaXZpc29yLCBpdGVtc1sxXS9kaXZpc29yXTtcclxuICAgIH1cclxuXHJcbiAgICAvL3ZhciBkb3RQb3cgPSBsb2dBbmRSdW4obWF0aC5kb3RQb3cpO1xyXG4gICAgc3RhdGljIHNxdWFyZUl0ZW1zKGl0ZW1zKSB7XHJcbiAgICAgICAgLy9yZXR1cm4gaXRlbXMubWFwKHggPT4geCp4KTtcclxuICAgICAgICB2YXIgYSA9IGl0ZW1zWzBdLCBiID0gaXRlbXNbMV07XHJcbiAgICAgICAgcmV0dXJuIFthKmEsIGIqYl07XHJcbiAgICB9XHJcblxyXG4gICAgc3RhdGljIG5vcm1hbGl6ZSh2KSB7XHJcbiAgICAgICAgcmV0dXJuIHRoaXMuZGl2SXRlbXModiwgdGhpcy52ZWN0b3JMZW4odikpO1xyXG4gICAgfVxyXG5cclxuICAgIC8vTWF0aC5wb3cgPSBsb2dBbmRSdW4oTWF0aC5wb3cpO1xyXG59XHJcblxyXG5cclxuY2xhc3MgYmV6aWVyIHtcclxuICAgIC8vRXZhbHVhdGVzIGN1YmljIGJlemllciBhdCB0LCByZXR1cm4gcG9pbnRcclxuICAgIHN0YXRpYyBxKGN0cmxQb2x5LCB0KSB7XHJcbiAgICAgICAgdmFyIHR4ID0gMS4wIC0gdDtcclxuICAgICAgICB2YXIgcEEgPSBtYXRocy5tdWxJdGVtcyggY3RybFBvbHlbMF0sICAgICAgdHggKiB0eCAqIHR4ICksXHJcbiAgICAgICAgICAgIHBCID0gbWF0aHMubXVsSXRlbXMoIGN0cmxQb2x5WzFdLCAgMyAqIHR4ICogdHggKiAgdCApLFxyXG4gICAgICAgICAgICBwQyA9IG1hdGhzLm11bEl0ZW1zKCBjdHJsUG9seVsyXSwgIDMgKiB0eCAqICB0ICogIHQgKSxcclxuICAgICAgICAgICAgcEQgPSBtYXRocy5tdWxJdGVtcyggY3RybFBvbHlbM10sICAgICAgIHQgKiAgdCAqICB0ICk7XHJcbiAgICAgICAgcmV0dXJuIG1hdGhzLmFkZEFycmF5cyhtYXRocy5hZGRBcnJheXMocEEsIHBCKSwgbWF0aHMuYWRkQXJyYXlzKHBDLCBwRCkpO1xyXG4gICAgfVxyXG5cclxuICAgIC8vRXZhbHVhdGVzIGN1YmljIGJlemllciBmaXJzdCBkZXJpdmF0aXZlIGF0IHQsIHJldHVybiBwb2ludFxyXG4gICAgc3RhdGljIHFwcmltZShjdHJsUG9seSwgdCkge1xyXG4gICAgICAgIHZhciB0eCA9IDEuMCAtIHQ7XHJcbiAgICAgICAgdmFyIHBBID0gbWF0aHMubXVsSXRlbXMoIG1hdGhzLnN1YnRyYWN0KGN0cmxQb2x5WzFdLCBjdHJsUG9seVswXSksICAzICogdHggKiB0eCApLFxyXG4gICAgICAgICAgICBwQiA9IG1hdGhzLm11bEl0ZW1zKCBtYXRocy5zdWJ0cmFjdChjdHJsUG9seVsyXSwgY3RybFBvbHlbMV0pLCAgNiAqIHR4ICogIHQgKSxcclxuICAgICAgICAgICAgcEMgPSBtYXRocy5tdWxJdGVtcyggbWF0aHMuc3VidHJhY3QoY3RybFBvbHlbM10sIGN0cmxQb2x5WzJdKSwgIDMgKiAgdCAqICB0ICk7XHJcbiAgICAgICAgcmV0dXJuIG1hdGhzLmFkZEFycmF5cyhtYXRocy5hZGRBcnJheXMocEEsIHBCKSwgcEMpO1xyXG4gICAgfVxyXG5cclxuICAgIC8vRXZhbHVhdGVzIGN1YmljIGJlemllciBzZWNvbmQgZGVyaXZhdGl2ZSBhdCB0LCByZXR1cm4gcG9pbnRcclxuICAgIHN0YXRpYyBxcHJpbWVwcmltZShjdHJsUG9seSwgdCkge1xyXG4gICAgICAgIHJldHVybiBtYXRocy5hZGRBcnJheXMobWF0aHMubXVsSXRlbXMoIG1hdGhzLmFkZEFycmF5cyhtYXRocy5zdWJ0cmFjdChjdHJsUG9seVsyXSwgbWF0aHMubXVsSXRlbXMoY3RybFBvbHlbMV0sIDIpKSwgY3RybFBvbHlbMF0pLCAgNiAqICgxLjAgLSB0KSApLFxyXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgbWF0aHMubXVsSXRlbXMoIG1hdGhzLmFkZEFycmF5cyhtYXRocy5zdWJ0cmFjdChjdHJsUG9seVszXSwgbWF0aHMubXVsSXRlbXMoY3RybFBvbHlbMl0sIDIpKSwgY3RybFBvbHlbMV0pLCAgNiAqICAgICAgICB0ICApKTtcclxuICAgIH1cclxufVxyXG5cclxubW9kdWxlLmV4cG9ydHMgPSBmaXRDdXJ2ZTtcclxuIiwiLyohXG4qIHN2Zy5qcyAtIEEgbGlnaHR3ZWlnaHQgbGlicmFyeSBmb3IgbWFuaXB1bGF0aW5nIGFuZCBhbmltYXRpbmcgU1ZHLlxuKiBAdmVyc2lvbiAyLjMuN1xuKiBodHRwczovL3N2Z2RvdGpzLmdpdGh1Yi5pby9cbipcbiogQGNvcHlyaWdodCBXb3V0IEZpZXJlbnMgPHdvdXRAbWljay13b3V0LmNvbT5cbiogQGxpY2Vuc2UgTUlUXG4qXG4qIEJVSUxUOiBTYXQgSmFuIDE0IDIwMTcgMDc6MjM6MTggR01UKzAxMDAgKENFVClcbiovO1xuKGZ1bmN0aW9uKHJvb3QsIGZhY3RvcnkpIHtcbiAgaWYgKHR5cGVvZiBkZWZpbmUgPT09ICdmdW5jdGlvbicgJiYgZGVmaW5lLmFtZCkge1xuICAgIGRlZmluZShmdW5jdGlvbigpe1xuICAgICAgcmV0dXJuIGZhY3Rvcnkocm9vdCwgcm9vdC5kb2N1bWVudClcbiAgICB9KVxuICB9IGVsc2UgaWYgKHR5cGVvZiBleHBvcnRzID09PSAnb2JqZWN0Jykge1xuICAgIG1vZHVsZS5leHBvcnRzID0gcm9vdC5kb2N1bWVudCA/IGZhY3Rvcnkocm9vdCwgcm9vdC5kb2N1bWVudCkgOiBmdW5jdGlvbih3KXsgcmV0dXJuIGZhY3Rvcnkodywgdy5kb2N1bWVudCkgfVxuICB9IGVsc2Uge1xuICAgIHJvb3QuU1ZHID0gZmFjdG9yeShyb290LCByb290LmRvY3VtZW50KVxuICB9XG59KHR5cGVvZiB3aW5kb3cgIT09IFwidW5kZWZpbmVkXCIgPyB3aW5kb3cgOiB0aGlzLCBmdW5jdGlvbih3aW5kb3csIGRvY3VtZW50KSB7XG5cbi8vIFRoZSBtYWluIHdyYXBwaW5nIGVsZW1lbnRcbnZhciBTVkcgPSB0aGlzLlNWRyA9IGZ1bmN0aW9uKGVsZW1lbnQpIHtcbiAgaWYgKFNWRy5zdXBwb3J0ZWQpIHtcbiAgICBlbGVtZW50ID0gbmV3IFNWRy5Eb2MoZWxlbWVudClcblxuICAgIGlmKCFTVkcucGFyc2VyLmRyYXcpXG4gICAgICBTVkcucHJlcGFyZSgpXG5cbiAgICByZXR1cm4gZWxlbWVudFxuICB9XG59XG5cbi8vIERlZmF1bHQgbmFtZXNwYWNlc1xuU1ZHLm5zICAgID0gJ2h0dHA6Ly93d3cudzMub3JnLzIwMDAvc3ZnJ1xuU1ZHLnhtbG5zID0gJ2h0dHA6Ly93d3cudzMub3JnLzIwMDAveG1sbnMvJ1xuU1ZHLnhsaW5rID0gJ2h0dHA6Ly93d3cudzMub3JnLzE5OTkveGxpbmsnXG5TVkcuc3ZnanMgPSAnaHR0cDovL3N2Z2pzLmNvbS9zdmdqcydcblxuLy8gU3ZnIHN1cHBvcnQgdGVzdFxuU1ZHLnN1cHBvcnRlZCA9IChmdW5jdGlvbigpIHtcbiAgcmV0dXJuICEhIGRvY3VtZW50LmNyZWF0ZUVsZW1lbnROUyAmJlxuICAgICAgICAgISEgZG9jdW1lbnQuY3JlYXRlRWxlbWVudE5TKFNWRy5ucywnc3ZnJykuY3JlYXRlU1ZHUmVjdFxufSkoKVxuXG4vLyBEb24ndCBib3RoZXIgdG8gY29udGludWUgaWYgU1ZHIGlzIG5vdCBzdXBwb3J0ZWRcbmlmICghU1ZHLnN1cHBvcnRlZCkgcmV0dXJuIGZhbHNlXG5cbi8vIEVsZW1lbnQgaWQgc2VxdWVuY2VcblNWRy5kaWQgID0gMTAwMFxuXG4vLyBHZXQgbmV4dCBuYW1lZCBlbGVtZW50IGlkXG5TVkcuZWlkID0gZnVuY3Rpb24obmFtZSkge1xuICByZXR1cm4gJ1N2Z2pzJyArIGNhcGl0YWxpemUobmFtZSkgKyAoU1ZHLmRpZCsrKVxufVxuXG4vLyBNZXRob2QgZm9yIGVsZW1lbnQgY3JlYXRpb25cblNWRy5jcmVhdGUgPSBmdW5jdGlvbihuYW1lKSB7XG4gIC8vIGNyZWF0ZSBlbGVtZW50XG4gIHZhciBlbGVtZW50ID0gZG9jdW1lbnQuY3JlYXRlRWxlbWVudE5TKHRoaXMubnMsIG5hbWUpXG5cbiAgLy8gYXBwbHkgdW5pcXVlIGlkXG4gIGVsZW1lbnQuc2V0QXR0cmlidXRlKCdpZCcsIHRoaXMuZWlkKG5hbWUpKVxuXG4gIHJldHVybiBlbGVtZW50XG59XG5cbi8vIE1ldGhvZCBmb3IgZXh0ZW5kaW5nIG9iamVjdHNcblNWRy5leHRlbmQgPSBmdW5jdGlvbigpIHtcbiAgdmFyIG1vZHVsZXMsIG1ldGhvZHMsIGtleSwgaVxuXG4gIC8vIEdldCBsaXN0IG9mIG1vZHVsZXNcbiAgbW9kdWxlcyA9IFtdLnNsaWNlLmNhbGwoYXJndW1lbnRzKVxuXG4gIC8vIEdldCBvYmplY3Qgd2l0aCBleHRlbnNpb25zXG4gIG1ldGhvZHMgPSBtb2R1bGVzLnBvcCgpXG5cbiAgZm9yIChpID0gbW9kdWxlcy5sZW5ndGggLSAxOyBpID49IDA7IGktLSlcbiAgICBpZiAobW9kdWxlc1tpXSlcbiAgICAgIGZvciAoa2V5IGluIG1ldGhvZHMpXG4gICAgICAgIG1vZHVsZXNbaV0ucHJvdG90eXBlW2tleV0gPSBtZXRob2RzW2tleV1cblxuICAvLyBNYWtlIHN1cmUgU1ZHLlNldCBpbmhlcml0cyBhbnkgbmV3bHkgYWRkZWQgbWV0aG9kc1xuICBpZiAoU1ZHLlNldCAmJiBTVkcuU2V0LmluaGVyaXQpXG4gICAgU1ZHLlNldC5pbmhlcml0KClcbn1cblxuLy8gSW52ZW50IG5ldyBlbGVtZW50XG5TVkcuaW52ZW50ID0gZnVuY3Rpb24oY29uZmlnKSB7XG4gIC8vIENyZWF0ZSBlbGVtZW50IGluaXRpYWxpemVyXG4gIHZhciBpbml0aWFsaXplciA9IHR5cGVvZiBjb25maWcuY3JlYXRlID09ICdmdW5jdGlvbicgP1xuICAgIGNvbmZpZy5jcmVhdGUgOlxuICAgIGZ1bmN0aW9uKCkge1xuICAgICAgdGhpcy5jb25zdHJ1Y3Rvci5jYWxsKHRoaXMsIFNWRy5jcmVhdGUoY29uZmlnLmNyZWF0ZSkpXG4gICAgfVxuXG4gIC8vIEluaGVyaXQgcHJvdG90eXBlXG4gIGlmIChjb25maWcuaW5oZXJpdClcbiAgICBpbml0aWFsaXplci5wcm90b3R5cGUgPSBuZXcgY29uZmlnLmluaGVyaXRcblxuICAvLyBFeHRlbmQgd2l0aCBtZXRob2RzXG4gIGlmIChjb25maWcuZXh0ZW5kKVxuICAgIFNWRy5leHRlbmQoaW5pdGlhbGl6ZXIsIGNvbmZpZy5leHRlbmQpXG5cbiAgLy8gQXR0YWNoIGNvbnN0cnVjdCBtZXRob2QgdG8gcGFyZW50XG4gIGlmIChjb25maWcuY29uc3RydWN0KVxuICAgIFNWRy5leHRlbmQoY29uZmlnLnBhcmVudCB8fCBTVkcuQ29udGFpbmVyLCBjb25maWcuY29uc3RydWN0KVxuXG4gIHJldHVybiBpbml0aWFsaXplclxufVxuXG4vLyBBZG9wdCBleGlzdGluZyBzdmcgZWxlbWVudHNcblNWRy5hZG9wdCA9IGZ1bmN0aW9uKG5vZGUpIHtcbiAgLy8gY2hlY2sgZm9yIHByZXNlbmNlIG9mIG5vZGVcbiAgaWYgKCFub2RlKSByZXR1cm4gbnVsbFxuXG4gIC8vIG1ha2Ugc3VyZSBhIG5vZGUgaXNuJ3QgYWxyZWFkeSBhZG9wdGVkXG4gIGlmIChub2RlLmluc3RhbmNlKSByZXR1cm4gbm9kZS5pbnN0YW5jZVxuXG4gIC8vIGluaXRpYWxpemUgdmFyaWFibGVzXG4gIHZhciBlbGVtZW50XG5cbiAgLy8gYWRvcHQgd2l0aCBlbGVtZW50LXNwZWNpZmljIHNldHRpbmdzXG4gIGlmIChub2RlLm5vZGVOYW1lID09ICdzdmcnKVxuICAgIGVsZW1lbnQgPSBub2RlLnBhcmVudE5vZGUgaW5zdGFuY2VvZiBTVkdFbGVtZW50ID8gbmV3IFNWRy5OZXN0ZWQgOiBuZXcgU1ZHLkRvY1xuICBlbHNlIGlmIChub2RlLm5vZGVOYW1lID09ICdsaW5lYXJHcmFkaWVudCcpXG4gICAgZWxlbWVudCA9IG5ldyBTVkcuR3JhZGllbnQoJ2xpbmVhcicpXG4gIGVsc2UgaWYgKG5vZGUubm9kZU5hbWUgPT0gJ3JhZGlhbEdyYWRpZW50JylcbiAgICBlbGVtZW50ID0gbmV3IFNWRy5HcmFkaWVudCgncmFkaWFsJylcbiAgZWxzZSBpZiAoU1ZHW2NhcGl0YWxpemUobm9kZS5ub2RlTmFtZSldKVxuICAgIGVsZW1lbnQgPSBuZXcgU1ZHW2NhcGl0YWxpemUobm9kZS5ub2RlTmFtZSldXG4gIGVsc2VcbiAgICBlbGVtZW50ID0gbmV3IFNWRy5FbGVtZW50KG5vZGUpXG5cbiAgLy8gZW5zdXJlIHJlZmVyZW5jZXNcbiAgZWxlbWVudC50eXBlICA9IG5vZGUubm9kZU5hbWVcbiAgZWxlbWVudC5ub2RlICA9IG5vZGVcbiAgbm9kZS5pbnN0YW5jZSA9IGVsZW1lbnRcblxuICAvLyBTVkcuQ2xhc3Mgc3BlY2lmaWMgcHJlcGFyYXRpb25zXG4gIGlmIChlbGVtZW50IGluc3RhbmNlb2YgU1ZHLkRvYylcbiAgICBlbGVtZW50Lm5hbWVzcGFjZSgpLmRlZnMoKVxuXG4gIC8vIHB1bGwgc3ZnanMgZGF0YSBmcm9tIHRoZSBkb20gKGdldEF0dHJpYnV0ZU5TIGRvZXNuJ3Qgd29yayBpbiBodG1sNSlcbiAgZWxlbWVudC5zZXREYXRhKEpTT04ucGFyc2Uobm9kZS5nZXRBdHRyaWJ1dGUoJ3N2Z2pzOmRhdGEnKSkgfHwge30pXG5cbiAgcmV0dXJuIGVsZW1lbnRcbn1cblxuLy8gSW5pdGlhbGl6ZSBwYXJzaW5nIGVsZW1lbnRcblNWRy5wcmVwYXJlID0gZnVuY3Rpb24oKSB7XG4gIC8vIFNlbGVjdCBkb2N1bWVudCBib2R5IGFuZCBjcmVhdGUgaW52aXNpYmxlIHN2ZyBlbGVtZW50XG4gIHZhciBib2R5ID0gZG9jdW1lbnQuZ2V0RWxlbWVudHNCeVRhZ05hbWUoJ2JvZHknKVswXVxuICAgICwgZHJhdyA9IChib2R5ID8gbmV3IFNWRy5Eb2MoYm9keSkgOiAgbmV3IFNWRy5Eb2MoZG9jdW1lbnQuZG9jdW1lbnRFbGVtZW50KS5uZXN0ZWQoKSkuc2l6ZSgyLCAwKVxuXG4gIC8vIENyZWF0ZSBwYXJzZXIgb2JqZWN0XG4gIFNWRy5wYXJzZXIgPSB7XG4gICAgYm9keTogYm9keSB8fCBkb2N1bWVudC5kb2N1bWVudEVsZW1lbnRcbiAgLCBkcmF3OiBkcmF3LnN0eWxlKCdvcGFjaXR5OjA7cG9zaXRpb246Zml4ZWQ7bGVmdDoxMDAlO3RvcDoxMDAlO292ZXJmbG93OmhpZGRlbicpXG4gICwgcG9seTogZHJhdy5wb2x5bGluZSgpLm5vZGVcbiAgLCBwYXRoOiBkcmF3LnBhdGgoKS5ub2RlXG4gICwgbmF0aXZlOiBTVkcuY3JlYXRlKCdzdmcnKVxuICB9XG59XG5cblNWRy5wYXJzZXIgPSB7XG4gIG5hdGl2ZTogU1ZHLmNyZWF0ZSgnc3ZnJylcbn1cblxuZG9jdW1lbnQuYWRkRXZlbnRMaXN0ZW5lcignRE9NQ29udGVudExvYWRlZCcsIGZ1bmN0aW9uKCkge1xuICBpZighU1ZHLnBhcnNlci5kcmF3KVxuICAgIFNWRy5wcmVwYXJlKClcbn0sIGZhbHNlKVxuXG4vLyBTdG9yYWdlIGZvciByZWd1bGFyIGV4cHJlc3Npb25zXG5TVkcucmVnZXggPSB7XG4gIC8vIFBhcnNlIHVuaXQgdmFsdWVcbiAgbnVtYmVyQW5kVW5pdDogICAgL14oWystXT8oXFxkKyhcXC5cXGQqKT98XFwuXFxkKykoZVsrLV0/XFxkKyk/KShbYS16JV0qKSQvaVxuXG4gIC8vIFBhcnNlIGhleCB2YWx1ZVxuLCBoZXg6ICAgICAgICAgICAgICAvXiM/KFthLWZcXGRdezJ9KShbYS1mXFxkXXsyfSkoW2EtZlxcZF17Mn0pJC9pXG5cbiAgLy8gUGFyc2UgcmdiIHZhbHVlXG4sIHJnYjogICAgICAgICAgICAgIC9yZ2JcXCgoXFxkKyksKFxcZCspLChcXGQrKVxcKS9cblxuICAvLyBQYXJzZSByZWZlcmVuY2UgaWRcbiwgcmVmZXJlbmNlOiAgICAgICAgLyMoW2EtejAtOVxcLV9dKykvaVxuXG4gIC8vIFBhcnNlIG1hdHJpeCB3cmFwcGVyXG4sIG1hdHJpeDogICAgICAgICAgIC9tYXRyaXhcXCh8XFwpL2dcblxuICAvLyBFbGVtZW50cyBvZiBhIG1hdHJpeFxuLCBtYXRyaXhFbGVtZW50czogICAvLCpcXHMrfCwvXG5cbiAgLy8gV2hpdGVzcGFjZVxuLCB3aGl0ZXNwYWNlOiAgICAgICAvXFxzL2dcblxuICAvLyBUZXN0IGhleCB2YWx1ZVxuLCBpc0hleDogICAgICAgICAgICAvXiNbYS1mMC05XXszLDZ9JC9pXG5cbiAgLy8gVGVzdCByZ2IgdmFsdWVcbiwgaXNSZ2I6ICAgICAgICAgICAgL15yZ2JcXCgvXG5cbiAgLy8gVGVzdCBjc3MgZGVjbGFyYXRpb25cbiwgaXNDc3M6ICAgICAgICAgICAgL1teOl0rOlteO10rOz8vXG5cbiAgLy8gVGVzdCBmb3IgYmxhbmsgc3RyaW5nXG4sIGlzQmxhbms6ICAgICAgICAgIC9eKFxccyspPyQvXG5cbiAgLy8gVGVzdCBmb3IgbnVtZXJpYyBzdHJpbmdcbiwgaXNOdW1iZXI6ICAgICAgICAgL15bKy1dPyhcXGQrKFxcLlxcZCopP3xcXC5cXGQrKShlWystXT9cXGQrKT8kL2lcblxuICAvLyBUZXN0IGZvciBwZXJjZW50IHZhbHVlXG4sIGlzUGVyY2VudDogICAgICAgIC9eLT9bXFxkXFwuXSslJC9cblxuICAvLyBUZXN0IGZvciBpbWFnZSB1cmxcbiwgaXNJbWFnZTogICAgICAgICAgL1xcLihqcGd8anBlZ3xwbmd8Z2lmfHN2ZykoXFw/W149XSsuKik/L2lcblxuICAvLyBUaGUgZm9sbG93aW5nIHJlZ2V4IGFyZSB1c2VkIHRvIHBhcnNlIHRoZSBkIGF0dHJpYnV0ZSBvZiBhIHBhdGhcblxuICAvLyBSZXBsYWNlcyBhbGwgbmVnYXRpdmUgZXhwb25lbnRzXG4sIG5lZ0V4cDogICAgICAgICAgIC9lXFwtL2dpXG5cbiAgLy8gUmVwbGFjZXMgYWxsIGNvbW1hXG4sIGNvbW1hOiAgICAgICAgICAgIC8sL2dcblxuICAvLyBSZXBsYWNlcyBhbGwgaHlwaGVuc1xuLCBoeXBoZW46ICAgICAgICAgICAvXFwtL2dcblxuICAvLyBSZXBsYWNlcyBhbmQgdGVzdHMgZm9yIGFsbCBwYXRoIGxldHRlcnNcbiwgcGF0aExldHRlcnM6ICAgICAgL1tNTEhWQ1NRVEFaXS9naVxuXG4gIC8vIHllcyB3ZSBuZWVkIHRoaXMgb25lLCB0b29cbiwgaXNQYXRoTGV0dGVyOiAgICAgL1tNTEhWQ1NRVEFaXS9pXG5cbiAgLy8gc3BsaXQgYXQgd2hpdGVzcGFjZXNcbiwgd2hpdGVzcGFjZXM6ICAgICAgL1xccysvXG5cbiAgLy8gbWF0Y2hlcyBYXG4sIFg6ICAgICAgICAgICAgICAgIC9YL2dcbn1cblxuU1ZHLnV0aWxzID0ge1xuICAvLyBNYXAgZnVuY3Rpb25cbiAgbWFwOiBmdW5jdGlvbihhcnJheSwgYmxvY2spIHtcbiAgICB2YXIgaVxuICAgICAgLCBpbCA9IGFycmF5Lmxlbmd0aFxuICAgICAgLCByZXN1bHQgPSBbXVxuXG4gICAgZm9yIChpID0gMDsgaSA8IGlsOyBpKyspXG4gICAgICByZXN1bHQucHVzaChibG9jayhhcnJheVtpXSkpXG5cbiAgICByZXR1cm4gcmVzdWx0XG4gIH1cblxuICAvLyBGaWx0ZXIgZnVuY3Rpb25cbiwgZmlsdGVyOiBmdW5jdGlvbihhcnJheSwgYmxvY2spIHtcbiAgICB2YXIgaVxuICAgICAgLCBpbCA9IGFycmF5Lmxlbmd0aFxuICAgICAgLCByZXN1bHQgPSBbXVxuXG4gICAgZm9yIChpID0gMDsgaSA8IGlsOyBpKyspXG4gICAgICBpZiAoYmxvY2soYXJyYXlbaV0pKVxuICAgICAgICByZXN1bHQucHVzaChhcnJheVtpXSlcblxuICAgIHJldHVybiByZXN1bHRcbiAgfVxuXG4gIC8vIERlZ3JlZXMgdG8gcmFkaWFuc1xuLCByYWRpYW5zOiBmdW5jdGlvbihkKSB7XG4gICAgcmV0dXJuIGQgJSAzNjAgKiBNYXRoLlBJIC8gMTgwXG4gIH1cblxuICAvLyBSYWRpYW5zIHRvIGRlZ3JlZXNcbiwgZGVncmVlczogZnVuY3Rpb24ocikge1xuICAgIHJldHVybiByICogMTgwIC8gTWF0aC5QSSAlIDM2MFxuICB9XG5cbiwgZmlsdGVyU1ZHRWxlbWVudHM6IGZ1bmN0aW9uKG5vZGVzKSB7XG4gICAgcmV0dXJuIHRoaXMuZmlsdGVyKCBub2RlcywgZnVuY3Rpb24oZWwpIHsgcmV0dXJuIGVsIGluc3RhbmNlb2YgU1ZHRWxlbWVudCB9KVxuICB9XG5cbn1cblxuU1ZHLmRlZmF1bHRzID0ge1xuICAvLyBEZWZhdWx0IGF0dHJpYnV0ZSB2YWx1ZXNcbiAgYXR0cnM6IHtcbiAgICAvLyBmaWxsIGFuZCBzdHJva2VcbiAgICAnZmlsbC1vcGFjaXR5JzogICAgIDFcbiAgLCAnc3Ryb2tlLW9wYWNpdHknOiAgIDFcbiAgLCAnc3Ryb2tlLXdpZHRoJzogICAgIDBcbiAgLCAnc3Ryb2tlLWxpbmVqb2luJzogICdtaXRlcidcbiAgLCAnc3Ryb2tlLWxpbmVjYXAnOiAgICdidXR0J1xuICAsIGZpbGw6ICAgICAgICAgICAgICAgJyMwMDAwMDAnXG4gICwgc3Ryb2tlOiAgICAgICAgICAgICAnIzAwMDAwMCdcbiAgLCBvcGFjaXR5OiAgICAgICAgICAgIDFcbiAgICAvLyBwb3NpdGlvblxuICAsIHg6ICAgICAgICAgICAgICAgICAgMFxuICAsIHk6ICAgICAgICAgICAgICAgICAgMFxuICAsIGN4OiAgICAgICAgICAgICAgICAgMFxuICAsIGN5OiAgICAgICAgICAgICAgICAgMFxuICAgIC8vIHNpemVcbiAgLCB3aWR0aDogICAgICAgICAgICAgIDBcbiAgLCBoZWlnaHQ6ICAgICAgICAgICAgIDBcbiAgICAvLyByYWRpdXNcbiAgLCByOiAgICAgICAgICAgICAgICAgIDBcbiAgLCByeDogICAgICAgICAgICAgICAgIDBcbiAgLCByeTogICAgICAgICAgICAgICAgIDBcbiAgICAvLyBncmFkaWVudFxuICAsIG9mZnNldDogICAgICAgICAgICAgMFxuICAsICdzdG9wLW9wYWNpdHknOiAgICAgMVxuICAsICdzdG9wLWNvbG9yJzogICAgICAgJyMwMDAwMDAnXG4gICAgLy8gdGV4dFxuICAsICdmb250LXNpemUnOiAgICAgICAgMTZcbiAgLCAnZm9udC1mYW1pbHknOiAgICAgICdIZWx2ZXRpY2EsIEFyaWFsLCBzYW5zLXNlcmlmJ1xuICAsICd0ZXh0LWFuY2hvcic6ICAgICAgJ3N0YXJ0J1xuICB9XG5cbn1cbi8vIE1vZHVsZSBmb3IgY29sb3IgY29udmVydGlvbnNcblNWRy5Db2xvciA9IGZ1bmN0aW9uKGNvbG9yKSB7XG4gIHZhciBtYXRjaFxuXG4gIC8vIGluaXRpYWxpemUgZGVmYXVsdHNcbiAgdGhpcy5yID0gMFxuICB0aGlzLmcgPSAwXG4gIHRoaXMuYiA9IDBcblxuICBpZighY29sb3IpIHJldHVyblxuXG4gIC8vIHBhcnNlIGNvbG9yXG4gIGlmICh0eXBlb2YgY29sb3IgPT09ICdzdHJpbmcnKSB7XG4gICAgaWYgKFNWRy5yZWdleC5pc1JnYi50ZXN0KGNvbG9yKSkge1xuICAgICAgLy8gZ2V0IHJnYiB2YWx1ZXNcbiAgICAgIG1hdGNoID0gU1ZHLnJlZ2V4LnJnYi5leGVjKGNvbG9yLnJlcGxhY2UoL1xccy9nLCcnKSlcblxuICAgICAgLy8gcGFyc2UgbnVtZXJpYyB2YWx1ZXNcbiAgICAgIHRoaXMuciA9IHBhcnNlSW50KG1hdGNoWzFdKVxuICAgICAgdGhpcy5nID0gcGFyc2VJbnQobWF0Y2hbMl0pXG4gICAgICB0aGlzLmIgPSBwYXJzZUludChtYXRjaFszXSlcblxuICAgIH0gZWxzZSBpZiAoU1ZHLnJlZ2V4LmlzSGV4LnRlc3QoY29sb3IpKSB7XG4gICAgICAvLyBnZXQgaGV4IHZhbHVlc1xuICAgICAgbWF0Y2ggPSBTVkcucmVnZXguaGV4LmV4ZWMoZnVsbEhleChjb2xvcikpXG5cbiAgICAgIC8vIHBhcnNlIG51bWVyaWMgdmFsdWVzXG4gICAgICB0aGlzLnIgPSBwYXJzZUludChtYXRjaFsxXSwgMTYpXG4gICAgICB0aGlzLmcgPSBwYXJzZUludChtYXRjaFsyXSwgMTYpXG4gICAgICB0aGlzLmIgPSBwYXJzZUludChtYXRjaFszXSwgMTYpXG5cbiAgICB9XG5cbiAgfSBlbHNlIGlmICh0eXBlb2YgY29sb3IgPT09ICdvYmplY3QnKSB7XG4gICAgdGhpcy5yID0gY29sb3IuclxuICAgIHRoaXMuZyA9IGNvbG9yLmdcbiAgICB0aGlzLmIgPSBjb2xvci5iXG5cbiAgfVxuXG59XG5cblNWRy5leHRlbmQoU1ZHLkNvbG9yLCB7XG4gIC8vIERlZmF1bHQgdG8gaGV4IGNvbnZlcnNpb25cbiAgdG9TdHJpbmc6IGZ1bmN0aW9uKCkge1xuICAgIHJldHVybiB0aGlzLnRvSGV4KClcbiAgfVxuICAvLyBCdWlsZCBoZXggdmFsdWVcbiwgdG9IZXg6IGZ1bmN0aW9uKCkge1xuICAgIHJldHVybiAnIydcbiAgICAgICsgY29tcFRvSGV4KHRoaXMucilcbiAgICAgICsgY29tcFRvSGV4KHRoaXMuZylcbiAgICAgICsgY29tcFRvSGV4KHRoaXMuYilcbiAgfVxuICAvLyBCdWlsZCByZ2IgdmFsdWVcbiwgdG9SZ2I6IGZ1bmN0aW9uKCkge1xuICAgIHJldHVybiAncmdiKCcgKyBbdGhpcy5yLCB0aGlzLmcsIHRoaXMuYl0uam9pbigpICsgJyknXG4gIH1cbiAgLy8gQ2FsY3VsYXRlIHRydWUgYnJpZ2h0bmVzc1xuLCBicmlnaHRuZXNzOiBmdW5jdGlvbigpIHtcbiAgICByZXR1cm4gKHRoaXMuciAvIDI1NSAqIDAuMzApXG4gICAgICAgICArICh0aGlzLmcgLyAyNTUgKiAwLjU5KVxuICAgICAgICAgKyAodGhpcy5iIC8gMjU1ICogMC4xMSlcbiAgfVxuICAvLyBNYWtlIGNvbG9yIG1vcnBoYWJsZVxuLCBtb3JwaDogZnVuY3Rpb24oY29sb3IpIHtcbiAgICB0aGlzLmRlc3RpbmF0aW9uID0gbmV3IFNWRy5Db2xvcihjb2xvcilcblxuICAgIHJldHVybiB0aGlzXG4gIH1cbiAgLy8gR2V0IG1vcnBoZWQgY29sb3IgYXQgZ2l2ZW4gcG9zaXRpb25cbiwgYXQ6IGZ1bmN0aW9uKHBvcykge1xuICAgIC8vIG1ha2Ugc3VyZSBhIGRlc3RpbmF0aW9uIGlzIGRlZmluZWRcbiAgICBpZiAoIXRoaXMuZGVzdGluYXRpb24pIHJldHVybiB0aGlzXG5cbiAgICAvLyBub3JtYWxpc2UgcG9zXG4gICAgcG9zID0gcG9zIDwgMCA/IDAgOiBwb3MgPiAxID8gMSA6IHBvc1xuXG4gICAgLy8gZ2VuZXJhdGUgbW9ycGhlZCBjb2xvclxuICAgIHJldHVybiBuZXcgU1ZHLkNvbG9yKHtcbiAgICAgIHI6IH5+KHRoaXMuciArICh0aGlzLmRlc3RpbmF0aW9uLnIgLSB0aGlzLnIpICogcG9zKVxuICAgICwgZzogfn4odGhpcy5nICsgKHRoaXMuZGVzdGluYXRpb24uZyAtIHRoaXMuZykgKiBwb3MpXG4gICAgLCBiOiB+fih0aGlzLmIgKyAodGhpcy5kZXN0aW5hdGlvbi5iIC0gdGhpcy5iKSAqIHBvcylcbiAgICB9KVxuICB9XG5cbn0pXG5cbi8vIFRlc3RlcnNcblxuLy8gVGVzdCBpZiBnaXZlbiB2YWx1ZSBpcyBhIGNvbG9yIHN0cmluZ1xuU1ZHLkNvbG9yLnRlc3QgPSBmdW5jdGlvbihjb2xvcikge1xuICBjb2xvciArPSAnJ1xuICByZXR1cm4gU1ZHLnJlZ2V4LmlzSGV4LnRlc3QoY29sb3IpXG4gICAgICB8fCBTVkcucmVnZXguaXNSZ2IudGVzdChjb2xvcilcbn1cblxuLy8gVGVzdCBpZiBnaXZlbiB2YWx1ZSBpcyBhIHJnYiBvYmplY3RcblNWRy5Db2xvci5pc1JnYiA9IGZ1bmN0aW9uKGNvbG9yKSB7XG4gIHJldHVybiBjb2xvciAmJiB0eXBlb2YgY29sb3IuciA9PSAnbnVtYmVyJ1xuICAgICAgICAgICAgICAgJiYgdHlwZW9mIGNvbG9yLmcgPT0gJ251bWJlcidcbiAgICAgICAgICAgICAgICYmIHR5cGVvZiBjb2xvci5iID09ICdudW1iZXInXG59XG5cbi8vIFRlc3QgaWYgZ2l2ZW4gdmFsdWUgaXMgYSBjb2xvclxuU1ZHLkNvbG9yLmlzQ29sb3IgPSBmdW5jdGlvbihjb2xvcikge1xuICByZXR1cm4gU1ZHLkNvbG9yLmlzUmdiKGNvbG9yKSB8fCBTVkcuQ29sb3IudGVzdChjb2xvcilcbn1cbi8vIE1vZHVsZSBmb3IgYXJyYXkgY29udmVyc2lvblxuU1ZHLkFycmF5ID0gZnVuY3Rpb24oYXJyYXksIGZhbGxiYWNrKSB7XG4gIGFycmF5ID0gKGFycmF5IHx8IFtdKS52YWx1ZU9mKClcblxuICAvLyBpZiBhcnJheSBpcyBlbXB0eSBhbmQgZmFsbGJhY2sgaXMgcHJvdmlkZWQsIHVzZSBmYWxsYmFja1xuICBpZiAoYXJyYXkubGVuZ3RoID09IDAgJiYgZmFsbGJhY2spXG4gICAgYXJyYXkgPSBmYWxsYmFjay52YWx1ZU9mKClcblxuICAvLyBwYXJzZSBhcnJheVxuICB0aGlzLnZhbHVlID0gdGhpcy5wYXJzZShhcnJheSlcbn1cblxuU1ZHLmV4dGVuZChTVkcuQXJyYXksIHtcbiAgLy8gTWFrZSBhcnJheSBtb3JwaGFibGVcbiAgbW9ycGg6IGZ1bmN0aW9uKGFycmF5KSB7XG4gICAgdGhpcy5kZXN0aW5hdGlvbiA9IHRoaXMucGFyc2UoYXJyYXkpXG5cbiAgICAvLyBub3JtYWxpemUgbGVuZ3RoIG9mIGFycmF5c1xuICAgIGlmICh0aGlzLnZhbHVlLmxlbmd0aCAhPSB0aGlzLmRlc3RpbmF0aW9uLmxlbmd0aCkge1xuICAgICAgdmFyIGxhc3RWYWx1ZSAgICAgICA9IHRoaXMudmFsdWVbdGhpcy52YWx1ZS5sZW5ndGggLSAxXVxuICAgICAgICAsIGxhc3REZXN0aW5hdGlvbiA9IHRoaXMuZGVzdGluYXRpb25bdGhpcy5kZXN0aW5hdGlvbi5sZW5ndGggLSAxXVxuXG4gICAgICB3aGlsZSh0aGlzLnZhbHVlLmxlbmd0aCA+IHRoaXMuZGVzdGluYXRpb24ubGVuZ3RoKVxuICAgICAgICB0aGlzLmRlc3RpbmF0aW9uLnB1c2gobGFzdERlc3RpbmF0aW9uKVxuICAgICAgd2hpbGUodGhpcy52YWx1ZS5sZW5ndGggPCB0aGlzLmRlc3RpbmF0aW9uLmxlbmd0aClcbiAgICAgICAgdGhpcy52YWx1ZS5wdXNoKGxhc3RWYWx1ZSlcbiAgICB9XG5cbiAgICByZXR1cm4gdGhpc1xuICB9XG4gIC8vIENsZWFuIHVwIGFueSBkdXBsaWNhdGUgcG9pbnRzXG4sIHNldHRsZTogZnVuY3Rpb24oKSB7XG4gICAgLy8gZmluZCBhbGwgdW5pcXVlIHZhbHVlc1xuICAgIGZvciAodmFyIGkgPSAwLCBpbCA9IHRoaXMudmFsdWUubGVuZ3RoLCBzZWVuID0gW107IGkgPCBpbDsgaSsrKVxuICAgICAgaWYgKHNlZW4uaW5kZXhPZih0aGlzLnZhbHVlW2ldKSA9PSAtMSlcbiAgICAgICAgc2Vlbi5wdXNoKHRoaXMudmFsdWVbaV0pXG5cbiAgICAvLyBzZXQgbmV3IHZhbHVlXG4gICAgcmV0dXJuIHRoaXMudmFsdWUgPSBzZWVuXG4gIH1cbiAgLy8gR2V0IG1vcnBoZWQgYXJyYXkgYXQgZ2l2ZW4gcG9zaXRpb25cbiwgYXQ6IGZ1bmN0aW9uKHBvcykge1xuICAgIC8vIG1ha2Ugc3VyZSBhIGRlc3RpbmF0aW9uIGlzIGRlZmluZWRcbiAgICBpZiAoIXRoaXMuZGVzdGluYXRpb24pIHJldHVybiB0aGlzXG5cbiAgICAvLyBnZW5lcmF0ZSBtb3JwaGVkIGFycmF5XG4gICAgZm9yICh2YXIgaSA9IDAsIGlsID0gdGhpcy52YWx1ZS5sZW5ndGgsIGFycmF5ID0gW107IGkgPCBpbDsgaSsrKVxuICAgICAgYXJyYXkucHVzaCh0aGlzLnZhbHVlW2ldICsgKHRoaXMuZGVzdGluYXRpb25baV0gLSB0aGlzLnZhbHVlW2ldKSAqIHBvcylcblxuICAgIHJldHVybiBuZXcgU1ZHLkFycmF5KGFycmF5KVxuICB9XG4gIC8vIENvbnZlcnQgYXJyYXkgdG8gc3RyaW5nXG4sIHRvU3RyaW5nOiBmdW5jdGlvbigpIHtcbiAgICByZXR1cm4gdGhpcy52YWx1ZS5qb2luKCcgJylcbiAgfVxuICAvLyBSZWFsIHZhbHVlXG4sIHZhbHVlT2Y6IGZ1bmN0aW9uKCkge1xuICAgIHJldHVybiB0aGlzLnZhbHVlXG4gIH1cbiAgLy8gUGFyc2Ugd2hpdGVzcGFjZSBzZXBhcmF0ZWQgc3RyaW5nXG4sIHBhcnNlOiBmdW5jdGlvbihhcnJheSkge1xuICAgIGFycmF5ID0gYXJyYXkudmFsdWVPZigpXG5cbiAgICAvLyBpZiBhbHJlYWR5IGlzIGFuIGFycmF5LCBubyBuZWVkIHRvIHBhcnNlIGl0XG4gICAgaWYgKEFycmF5LmlzQXJyYXkoYXJyYXkpKSByZXR1cm4gYXJyYXlcblxuICAgIHJldHVybiB0aGlzLnNwbGl0KGFycmF5KVxuICB9XG4gIC8vIFN0cmlwIHVubmVjZXNzYXJ5IHdoaXRlc3BhY2Vcbiwgc3BsaXQ6IGZ1bmN0aW9uKHN0cmluZykge1xuICAgIHJldHVybiBzdHJpbmcudHJpbSgpLnNwbGl0KC9cXHMrLylcbiAgfVxuICAvLyBSZXZlcnNlIGFycmF5XG4sIHJldmVyc2U6IGZ1bmN0aW9uKCkge1xuICAgIHRoaXMudmFsdWUucmV2ZXJzZSgpXG5cbiAgICByZXR1cm4gdGhpc1xuICB9XG5cbn0pXG4vLyBQb2x5IHBvaW50cyBhcnJheVxuU1ZHLlBvaW50QXJyYXkgPSBmdW5jdGlvbihhcnJheSwgZmFsbGJhY2spIHtcbiAgdGhpcy5jb25zdHJ1Y3Rvci5jYWxsKHRoaXMsIGFycmF5LCBmYWxsYmFjayB8fCBbWzAsMF1dKVxufVxuXG4vLyBJbmhlcml0IGZyb20gU1ZHLkFycmF5XG5TVkcuUG9pbnRBcnJheS5wcm90b3R5cGUgPSBuZXcgU1ZHLkFycmF5XG5cblNWRy5leHRlbmQoU1ZHLlBvaW50QXJyYXksIHtcbiAgLy8gQ29udmVydCBhcnJheSB0byBzdHJpbmdcbiAgdG9TdHJpbmc6IGZ1bmN0aW9uKCkge1xuICAgIC8vIGNvbnZlcnQgdG8gYSBwb2x5IHBvaW50IHN0cmluZ1xuICAgIGZvciAodmFyIGkgPSAwLCBpbCA9IHRoaXMudmFsdWUubGVuZ3RoLCBhcnJheSA9IFtdOyBpIDwgaWw7IGkrKylcbiAgICAgIGFycmF5LnB1c2godGhpcy52YWx1ZVtpXS5qb2luKCcsJykpXG5cbiAgICByZXR1cm4gYXJyYXkuam9pbignICcpXG4gIH1cbiAgLy8gQ29udmVydCBhcnJheSB0byBsaW5lIG9iamVjdFxuLCB0b0xpbmU6IGZ1bmN0aW9uKCkge1xuICAgIHJldHVybiB7XG4gICAgICB4MTogdGhpcy52YWx1ZVswXVswXVxuICAgICwgeTE6IHRoaXMudmFsdWVbMF1bMV1cbiAgICAsIHgyOiB0aGlzLnZhbHVlWzFdWzBdXG4gICAgLCB5MjogdGhpcy52YWx1ZVsxXVsxXVxuICAgIH1cbiAgfVxuICAvLyBHZXQgbW9ycGhlZCBhcnJheSBhdCBnaXZlbiBwb3NpdGlvblxuLCBhdDogZnVuY3Rpb24ocG9zKSB7XG4gICAgLy8gbWFrZSBzdXJlIGEgZGVzdGluYXRpb24gaXMgZGVmaW5lZFxuICAgIGlmICghdGhpcy5kZXN0aW5hdGlvbikgcmV0dXJuIHRoaXNcblxuICAgIC8vIGdlbmVyYXRlIG1vcnBoZWQgcG9pbnQgc3RyaW5nXG4gICAgZm9yICh2YXIgaSA9IDAsIGlsID0gdGhpcy52YWx1ZS5sZW5ndGgsIGFycmF5ID0gW107IGkgPCBpbDsgaSsrKVxuICAgICAgYXJyYXkucHVzaChbXG4gICAgICAgIHRoaXMudmFsdWVbaV1bMF0gKyAodGhpcy5kZXN0aW5hdGlvbltpXVswXSAtIHRoaXMudmFsdWVbaV1bMF0pICogcG9zXG4gICAgICAsIHRoaXMudmFsdWVbaV1bMV0gKyAodGhpcy5kZXN0aW5hdGlvbltpXVsxXSAtIHRoaXMudmFsdWVbaV1bMV0pICogcG9zXG4gICAgICBdKVxuXG4gICAgcmV0dXJuIG5ldyBTVkcuUG9pbnRBcnJheShhcnJheSlcbiAgfVxuICAvLyBQYXJzZSBwb2ludCBzdHJpbmdcbiwgcGFyc2U6IGZ1bmN0aW9uKGFycmF5KSB7XG4gICAgdmFyIHBvaW50cyA9IFtdXG5cbiAgICBhcnJheSA9IGFycmF5LnZhbHVlT2YoKVxuXG4gICAgLy8gaWYgYWxyZWFkeSBpcyBhbiBhcnJheSwgbm8gbmVlZCB0byBwYXJzZSBpdFxuICAgIGlmIChBcnJheS5pc0FycmF5KGFycmF5KSkgcmV0dXJuIGFycmF5XG5cbiAgICAvLyBwYXJzZSBwb2ludHNcbiAgICBhcnJheSA9IGFycmF5LnRyaW0oKS5zcGxpdCgvXFxzK3wsLylcblxuICAgIC8vIHZhbGlkYXRlIHBvaW50cyAtIGh0dHBzOi8vc3Znd2cub3JnL3N2ZzItZHJhZnQvc2hhcGVzLmh0bWwjRGF0YVR5cGVQb2ludHNcbiAgICAvLyBPZGQgbnVtYmVyIG9mIGNvb3JkaW5hdGVzIGlzIGFuIGVycm9yLiBJbiBzdWNoIGNhc2VzLCBkcm9wIHRoZSBsYXN0IG9kZCBjb29yZGluYXRlLlxuICAgIGlmIChhcnJheS5sZW5ndGggJSAyICE9PSAwKSBhcnJheS5wb3AoKVxuXG4gICAgLy8gd3JhcCBwb2ludHMgaW4gdHdvLXR1cGxlcyBhbmQgcGFyc2UgcG9pbnRzIGFzIGZsb2F0c1xuICAgIGZvcih2YXIgaSA9IDAsIGxlbiA9IGFycmF5Lmxlbmd0aDsgaSA8IGxlbjsgaSA9IGkgKyAyKVxuICAgICAgcG9pbnRzLnB1c2goWyBwYXJzZUZsb2F0KGFycmF5W2ldKSwgcGFyc2VGbG9hdChhcnJheVtpKzFdKSBdKVxuXG4gICAgcmV0dXJuIHBvaW50c1xuICB9XG4gIC8vIE1vdmUgcG9pbnQgc3RyaW5nXG4sIG1vdmU6IGZ1bmN0aW9uKHgsIHkpIHtcbiAgICB2YXIgYm94ID0gdGhpcy5iYm94KClcblxuICAgIC8vIGdldCByZWxhdGl2ZSBvZmZzZXRcbiAgICB4IC09IGJveC54XG4gICAgeSAtPSBib3gueVxuXG4gICAgLy8gbW92ZSBldmVyeSBwb2ludFxuICAgIGlmICghaXNOYU4oeCkgJiYgIWlzTmFOKHkpKVxuICAgICAgZm9yICh2YXIgaSA9IHRoaXMudmFsdWUubGVuZ3RoIC0gMTsgaSA+PSAwOyBpLS0pXG4gICAgICAgIHRoaXMudmFsdWVbaV0gPSBbdGhpcy52YWx1ZVtpXVswXSArIHgsIHRoaXMudmFsdWVbaV1bMV0gKyB5XVxuXG4gICAgcmV0dXJuIHRoaXNcbiAgfVxuICAvLyBSZXNpemUgcG9seSBzdHJpbmdcbiwgc2l6ZTogZnVuY3Rpb24od2lkdGgsIGhlaWdodCkge1xuICAgIHZhciBpLCBib3ggPSB0aGlzLmJib3goKVxuXG4gICAgLy8gcmVjYWxjdWxhdGUgcG9zaXRpb24gb2YgYWxsIHBvaW50cyBhY2NvcmRpbmcgdG8gbmV3IHNpemVcbiAgICBmb3IgKGkgPSB0aGlzLnZhbHVlLmxlbmd0aCAtIDE7IGkgPj0gMDsgaS0tKSB7XG4gICAgICB0aGlzLnZhbHVlW2ldWzBdID0gKCh0aGlzLnZhbHVlW2ldWzBdIC0gYm94LngpICogd2lkdGgpICAvIGJveC53aWR0aCAgKyBib3gueFxuICAgICAgdGhpcy52YWx1ZVtpXVsxXSA9ICgodGhpcy52YWx1ZVtpXVsxXSAtIGJveC55KSAqIGhlaWdodCkgLyBib3guaGVpZ2h0ICsgYm94LnlcbiAgICB9XG5cbiAgICByZXR1cm4gdGhpc1xuICB9XG4gIC8vIEdldCBib3VuZGluZyBib3ggb2YgcG9pbnRzXG4sIGJib3g6IGZ1bmN0aW9uKCkge1xuICAgIFNWRy5wYXJzZXIucG9seS5zZXRBdHRyaWJ1dGUoJ3BvaW50cycsIHRoaXMudG9TdHJpbmcoKSlcblxuICAgIHJldHVybiBTVkcucGFyc2VyLnBvbHkuZ2V0QkJveCgpXG4gIH1cblxufSlcbi8vIFBhdGggcG9pbnRzIGFycmF5XG5TVkcuUGF0aEFycmF5ID0gZnVuY3Rpb24oYXJyYXksIGZhbGxiYWNrKSB7XG4gIHRoaXMuY29uc3RydWN0b3IuY2FsbCh0aGlzLCBhcnJheSwgZmFsbGJhY2sgfHwgW1snTScsIDAsIDBdXSlcbn1cblxuLy8gSW5oZXJpdCBmcm9tIFNWRy5BcnJheVxuU1ZHLlBhdGhBcnJheS5wcm90b3R5cGUgPSBuZXcgU1ZHLkFycmF5XG5cblNWRy5leHRlbmQoU1ZHLlBhdGhBcnJheSwge1xuICAvLyBDb252ZXJ0IGFycmF5IHRvIHN0cmluZ1xuICB0b1N0cmluZzogZnVuY3Rpb24oKSB7XG4gICAgcmV0dXJuIGFycmF5VG9TdHJpbmcodGhpcy52YWx1ZSlcbiAgfVxuICAvLyBNb3ZlIHBhdGggc3RyaW5nXG4sIG1vdmU6IGZ1bmN0aW9uKHgsIHkpIHtcbiAgICAvLyBnZXQgYm91bmRpbmcgYm94IG9mIGN1cnJlbnQgc2l0dWF0aW9uXG4gICAgdmFyIGJveCA9IHRoaXMuYmJveCgpXG5cbiAgICAvLyBnZXQgcmVsYXRpdmUgb2Zmc2V0XG4gICAgeCAtPSBib3gueFxuICAgIHkgLT0gYm94LnlcblxuICAgIGlmICghaXNOYU4oeCkgJiYgIWlzTmFOKHkpKSB7XG4gICAgICAvLyBtb3ZlIGV2ZXJ5IHBvaW50XG4gICAgICBmb3IgKHZhciBsLCBpID0gdGhpcy52YWx1ZS5sZW5ndGggLSAxOyBpID49IDA7IGktLSkge1xuICAgICAgICBsID0gdGhpcy52YWx1ZVtpXVswXVxuXG4gICAgICAgIGlmIChsID09ICdNJyB8fCBsID09ICdMJyB8fCBsID09ICdUJykgIHtcbiAgICAgICAgICB0aGlzLnZhbHVlW2ldWzFdICs9IHhcbiAgICAgICAgICB0aGlzLnZhbHVlW2ldWzJdICs9IHlcblxuICAgICAgICB9IGVsc2UgaWYgKGwgPT0gJ0gnKSAge1xuICAgICAgICAgIHRoaXMudmFsdWVbaV1bMV0gKz0geFxuXG4gICAgICAgIH0gZWxzZSBpZiAobCA9PSAnVicpICB7XG4gICAgICAgICAgdGhpcy52YWx1ZVtpXVsxXSArPSB5XG5cbiAgICAgICAgfSBlbHNlIGlmIChsID09ICdDJyB8fCBsID09ICdTJyB8fCBsID09ICdRJykgIHtcbiAgICAgICAgICB0aGlzLnZhbHVlW2ldWzFdICs9IHhcbiAgICAgICAgICB0aGlzLnZhbHVlW2ldWzJdICs9IHlcbiAgICAgICAgICB0aGlzLnZhbHVlW2ldWzNdICs9IHhcbiAgICAgICAgICB0aGlzLnZhbHVlW2ldWzRdICs9IHlcblxuICAgICAgICAgIGlmIChsID09ICdDJykgIHtcbiAgICAgICAgICAgIHRoaXMudmFsdWVbaV1bNV0gKz0geFxuICAgICAgICAgICAgdGhpcy52YWx1ZVtpXVs2XSArPSB5XG4gICAgICAgICAgfVxuXG4gICAgICAgIH0gZWxzZSBpZiAobCA9PSAnQScpICB7XG4gICAgICAgICAgdGhpcy52YWx1ZVtpXVs2XSArPSB4XG4gICAgICAgICAgdGhpcy52YWx1ZVtpXVs3XSArPSB5XG4gICAgICAgIH1cblxuICAgICAgfVxuICAgIH1cblxuICAgIHJldHVybiB0aGlzXG4gIH1cbiAgLy8gUmVzaXplIHBhdGggc3RyaW5nXG4sIHNpemU6IGZ1bmN0aW9uKHdpZHRoLCBoZWlnaHQpIHtcbiAgICAvLyBnZXQgYm91bmRpbmcgYm94IG9mIGN1cnJlbnQgc2l0dWF0aW9uXG4gICAgdmFyIGksIGwsIGJveCA9IHRoaXMuYmJveCgpXG5cbiAgICAvLyByZWNhbGN1bGF0ZSBwb3NpdGlvbiBvZiBhbGwgcG9pbnRzIGFjY29yZGluZyB0byBuZXcgc2l6ZVxuICAgIGZvciAoaSA9IHRoaXMudmFsdWUubGVuZ3RoIC0gMTsgaSA+PSAwOyBpLS0pIHtcbiAgICAgIGwgPSB0aGlzLnZhbHVlW2ldWzBdXG5cbiAgICAgIGlmIChsID09ICdNJyB8fCBsID09ICdMJyB8fCBsID09ICdUJykgIHtcbiAgICAgICAgdGhpcy52YWx1ZVtpXVsxXSA9ICgodGhpcy52YWx1ZVtpXVsxXSAtIGJveC54KSAqIHdpZHRoKSAgLyBib3gud2lkdGggICsgYm94LnhcbiAgICAgICAgdGhpcy52YWx1ZVtpXVsyXSA9ICgodGhpcy52YWx1ZVtpXVsyXSAtIGJveC55KSAqIGhlaWdodCkgLyBib3guaGVpZ2h0ICsgYm94LnlcblxuICAgICAgfSBlbHNlIGlmIChsID09ICdIJykgIHtcbiAgICAgICAgdGhpcy52YWx1ZVtpXVsxXSA9ICgodGhpcy52YWx1ZVtpXVsxXSAtIGJveC54KSAqIHdpZHRoKSAgLyBib3gud2lkdGggICsgYm94LnhcblxuICAgICAgfSBlbHNlIGlmIChsID09ICdWJykgIHtcbiAgICAgICAgdGhpcy52YWx1ZVtpXVsxXSA9ICgodGhpcy52YWx1ZVtpXVsxXSAtIGJveC55KSAqIGhlaWdodCkgLyBib3guaGVpZ2h0ICsgYm94LnlcblxuICAgICAgfSBlbHNlIGlmIChsID09ICdDJyB8fCBsID09ICdTJyB8fCBsID09ICdRJykgIHtcbiAgICAgICAgdGhpcy52YWx1ZVtpXVsxXSA9ICgodGhpcy52YWx1ZVtpXVsxXSAtIGJveC54KSAqIHdpZHRoKSAgLyBib3gud2lkdGggICsgYm94LnhcbiAgICAgICAgdGhpcy52YWx1ZVtpXVsyXSA9ICgodGhpcy52YWx1ZVtpXVsyXSAtIGJveC55KSAqIGhlaWdodCkgLyBib3guaGVpZ2h0ICsgYm94LnlcbiAgICAgICAgdGhpcy52YWx1ZVtpXVszXSA9ICgodGhpcy52YWx1ZVtpXVszXSAtIGJveC54KSAqIHdpZHRoKSAgLyBib3gud2lkdGggICsgYm94LnhcbiAgICAgICAgdGhpcy52YWx1ZVtpXVs0XSA9ICgodGhpcy52YWx1ZVtpXVs0XSAtIGJveC55KSAqIGhlaWdodCkgLyBib3guaGVpZ2h0ICsgYm94LnlcblxuICAgICAgICBpZiAobCA9PSAnQycpICB7XG4gICAgICAgICAgdGhpcy52YWx1ZVtpXVs1XSA9ICgodGhpcy52YWx1ZVtpXVs1XSAtIGJveC54KSAqIHdpZHRoKSAgLyBib3gud2lkdGggICsgYm94LnhcbiAgICAgICAgICB0aGlzLnZhbHVlW2ldWzZdID0gKCh0aGlzLnZhbHVlW2ldWzZdIC0gYm94LnkpICogaGVpZ2h0KSAvIGJveC5oZWlnaHQgKyBib3gueVxuICAgICAgICB9XG5cbiAgICAgIH0gZWxzZSBpZiAobCA9PSAnQScpICB7XG4gICAgICAgIC8vIHJlc2l6ZSByYWRpaVxuICAgICAgICB0aGlzLnZhbHVlW2ldWzFdID0gKHRoaXMudmFsdWVbaV1bMV0gKiB3aWR0aCkgIC8gYm94LndpZHRoXG4gICAgICAgIHRoaXMudmFsdWVbaV1bMl0gPSAodGhpcy52YWx1ZVtpXVsyXSAqIGhlaWdodCkgLyBib3guaGVpZ2h0XG5cbiAgICAgICAgLy8gbW92ZSBwb3NpdGlvbiB2YWx1ZXNcbiAgICAgICAgdGhpcy52YWx1ZVtpXVs2XSA9ICgodGhpcy52YWx1ZVtpXVs2XSAtIGJveC54KSAqIHdpZHRoKSAgLyBib3gud2lkdGggICsgYm94LnhcbiAgICAgICAgdGhpcy52YWx1ZVtpXVs3XSA9ICgodGhpcy52YWx1ZVtpXVs3XSAtIGJveC55KSAqIGhlaWdodCkgLyBib3guaGVpZ2h0ICsgYm94LnlcbiAgICAgIH1cblxuICAgIH1cblxuICAgIHJldHVybiB0aGlzXG4gIH1cbiAgLy8gVGVzdCBpZiB0aGUgcGFzc2VkIHBhdGggYXJyYXkgdXNlIHRoZSBzYW1lIHBhdGggZGF0YSBjb21tYW5kcyBhcyB0aGlzIHBhdGggYXJyYXlcbiwgZXF1YWxDb21tYW5kczogZnVuY3Rpb24ocGF0aEFycmF5KSB7XG4gICAgdmFyIGksIGlsLCBlcXVhbENvbW1hbmRzXG5cbiAgICBwYXRoQXJyYXkgPSBuZXcgU1ZHLlBhdGhBcnJheShwYXRoQXJyYXkpXG5cbiAgICBlcXVhbENvbW1hbmRzID0gdGhpcy52YWx1ZS5sZW5ndGggPT09IHBhdGhBcnJheS52YWx1ZS5sZW5ndGhcbiAgICBmb3IoaSA9IDAsIGlsID0gdGhpcy52YWx1ZS5sZW5ndGg7IGVxdWFsQ29tbWFuZHMgJiYgaSA8IGlsOyBpKyspIHtcbiAgICAgIGVxdWFsQ29tbWFuZHMgPSB0aGlzLnZhbHVlW2ldWzBdID09PSBwYXRoQXJyYXkudmFsdWVbaV1bMF1cbiAgICB9XG5cbiAgICByZXR1cm4gZXF1YWxDb21tYW5kc1xuICB9XG4gIC8vIE1ha2UgcGF0aCBhcnJheSBtb3JwaGFibGVcbiwgbW9ycGg6IGZ1bmN0aW9uKHBhdGhBcnJheSkge1xuICAgIHBhdGhBcnJheSA9IG5ldyBTVkcuUGF0aEFycmF5KHBhdGhBcnJheSlcblxuICAgIGlmKHRoaXMuZXF1YWxDb21tYW5kcyhwYXRoQXJyYXkpKSB7XG4gICAgICB0aGlzLmRlc3RpbmF0aW9uID0gcGF0aEFycmF5XG4gICAgfSBlbHNlIHtcbiAgICAgIHRoaXMuZGVzdGluYXRpb24gPSBudWxsXG4gICAgfVxuXG4gICAgcmV0dXJuIHRoaXNcbiAgfVxuICAvLyBHZXQgbW9ycGhlZCBwYXRoIGFycmF5IGF0IGdpdmVuIHBvc2l0aW9uXG4sIGF0OiBmdW5jdGlvbihwb3MpIHtcbiAgICAvLyBtYWtlIHN1cmUgYSBkZXN0aW5hdGlvbiBpcyBkZWZpbmVkXG4gICAgaWYgKCF0aGlzLmRlc3RpbmF0aW9uKSByZXR1cm4gdGhpc1xuXG4gICAgdmFyIHNvdXJjZUFycmF5ID0gdGhpcy52YWx1ZVxuICAgICAgLCBkZXN0aW5hdGlvbkFycmF5ID0gdGhpcy5kZXN0aW5hdGlvbi52YWx1ZVxuICAgICAgLCBhcnJheSA9IFtdLCBwYXRoQXJyYXkgPSBuZXcgU1ZHLlBhdGhBcnJheSgpXG4gICAgICAsIGksIGlsLCBqLCBqbFxuXG4gICAgLy8gQW5pbWF0ZSBoYXMgc3BlY2lmaWVkIGluIHRoZSBTVkcgc3BlY1xuICAgIC8vIFNlZTogaHR0cHM6Ly93d3cudzMub3JnL1RSL1NWRzExL3BhdGhzLmh0bWwjUGF0aEVsZW1lbnRcbiAgICBmb3IgKGkgPSAwLCBpbCA9IHNvdXJjZUFycmF5Lmxlbmd0aDsgaSA8IGlsOyBpKyspIHtcbiAgICAgIGFycmF5W2ldID0gW3NvdXJjZUFycmF5W2ldWzBdXVxuICAgICAgZm9yKGogPSAxLCBqbCA9IHNvdXJjZUFycmF5W2ldLmxlbmd0aDsgaiA8IGpsOyBqKyspIHtcbiAgICAgICAgYXJyYXlbaV1bal0gPSBzb3VyY2VBcnJheVtpXVtqXSArIChkZXN0aW5hdGlvbkFycmF5W2ldW2pdIC0gc291cmNlQXJyYXlbaV1bal0pICogcG9zXG4gICAgICB9XG4gICAgICAvLyBGb3IgdGhlIHR3byBmbGFncyBvZiB0aGUgZWxsaXB0aWNhbCBhcmMgY29tbWFuZCwgdGhlIFNWRyBzcGVjIHNheTpcbiAgICAgIC8vIEZsYWdzIGFuZCBib29sZWFucyBhcmUgaW50ZXJwb2xhdGVkIGFzIGZyYWN0aW9ucyBiZXR3ZWVuIHplcm8gYW5kIG9uZSwgd2l0aCBhbnkgbm9uLXplcm8gdmFsdWUgY29uc2lkZXJlZCB0byBiZSBhIHZhbHVlIG9mIG9uZS90cnVlXG4gICAgICAvLyBFbGxpcHRpY2FsIGFyYyBjb21tYW5kIGFzIGFuIGFycmF5IGZvbGxvd2VkIGJ5IGNvcnJlc3BvbmRpbmcgaW5kZXhlczpcbiAgICAgIC8vIFsnQScsIHJ4LCByeSwgeC1heGlzLXJvdGF0aW9uLCBsYXJnZS1hcmMtZmxhZywgc3dlZXAtZmxhZywgeCwgeV1cbiAgICAgIC8vICAgMCAgICAxICAgMiAgICAgICAgMyAgICAgICAgICAgICAgICAgNCAgICAgICAgICAgICA1ICAgICAgNiAgN1xuICAgICAgaWYoYXJyYXlbaV1bMF0gPT09ICdBJykge1xuICAgICAgICBhcnJheVtpXVs0XSA9ICsoYXJyYXlbaV1bNF0gIT0gMClcbiAgICAgICAgYXJyYXlbaV1bNV0gPSArKGFycmF5W2ldWzVdICE9IDApXG4gICAgICB9XG4gICAgfVxuXG4gICAgLy8gRGlyZWN0bHkgbW9kaWZ5IHRoZSB2YWx1ZSBvZiBhIHBhdGggYXJyYXksIHRoaXMgaXMgZG9uZSB0aGlzIHdheSBmb3IgcGVyZm9ybWFuY2VcbiAgICBwYXRoQXJyYXkudmFsdWUgPSBhcnJheVxuICAgIHJldHVybiBwYXRoQXJyYXlcbiAgfVxuICAvLyBBYnNvbHV0aXplIGFuZCBwYXJzZSBwYXRoIHRvIGFycmF5XG4sIHBhcnNlOiBmdW5jdGlvbihhcnJheSkge1xuICAgIC8vIGlmIGl0J3MgYWxyZWFkeSBhIHBhdGhhcnJheSwgbm8gbmVlZCB0byBwYXJzZSBpdFxuICAgIGlmIChhcnJheSBpbnN0YW5jZW9mIFNWRy5QYXRoQXJyYXkpIHJldHVybiBhcnJheS52YWx1ZU9mKClcblxuICAgIC8vIHByZXBhcmUgZm9yIHBhcnNpbmdcbiAgICB2YXIgaSwgeDAsIHkwLCBzLCBzZWcsIGFyclxuICAgICAgLCB4ID0gMFxuICAgICAgLCB5ID0gMFxuICAgICAgLCBwYXJhbUNudCA9IHsgJ00nOjIsICdMJzoyLCAnSCc6MSwgJ1YnOjEsICdDJzo2LCAnUyc6NCwgJ1EnOjQsICdUJzoyLCAnQSc6NyB9XG5cbiAgICBpZih0eXBlb2YgYXJyYXkgPT0gJ3N0cmluZycpe1xuXG4gICAgICBhcnJheSA9IGFycmF5XG4gICAgICAgIC5yZXBsYWNlKFNWRy5yZWdleC5uZWdFeHAsICdYJykgICAgICAgICAvLyByZXBsYWNlIGFsbCBuZWdhdGl2ZSBleHBvbmVudHMgd2l0aCBjZXJ0YWluIGNoYXJcbiAgICAgICAgLnJlcGxhY2UoU1ZHLnJlZ2V4LnBhdGhMZXR0ZXJzLCAnICQmICcpIC8vIHB1dCBzb21lIHJvb20gYmV0d2VlbiBsZXR0ZXJzIGFuZCBudW1iZXJzXG4gICAgICAgIC5yZXBsYWNlKFNWRy5yZWdleC5oeXBoZW4sICcgLScpICAgICAgICAvLyBhZGQgc3BhY2UgYmVmb3JlIGh5cGhlblxuICAgICAgICAucmVwbGFjZShTVkcucmVnZXguY29tbWEsICcgJykgICAgICAgICAgLy8gdW5pZnkgYWxsIHNwYWNlc1xuICAgICAgICAucmVwbGFjZShTVkcucmVnZXguWCwgJ2UtJykgICAgICAgICAgICAgLy8gYWRkIGJhY2sgdGhlIGV4cG9lbnRcbiAgICAgICAgLnRyaW0oKSAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIC8vIHRyaW1cbiAgICAgICAgLnNwbGl0KFNWRy5yZWdleC53aGl0ZXNwYWNlcykgICAgICAgICAgIC8vIHNwbGl0IGludG8gYXJyYXlcblxuICAgICAgLy8gYXQgdGhpcyBwbGFjZSB0aGVyZSBjb3VsZCBiZSBwYXJ0cyBsaWtlIFsnMy4xMjQuODU0LjMyJ10gYmVjYXVzZSB3ZSBjb3VsZCBub3QgZGV0ZXJtaW5lIHRoZSBwb2ludCBhcyBzZXBlcmF0b3IgdGlsbCBub3dcbiAgICAgIC8vIHdlIGZpeCB0aGlzIGVsZW1lbnRzIGluIHRoZSBuZXh0IGxvb3BcbiAgICAgIGZvcihpID0gYXJyYXkubGVuZ3RoOyAtLWk7KXtcbiAgICAgICAgaWYoYXJyYXlbaV0uaW5kZXhPZignLicpICE9IGFycmF5W2ldLmxhc3RJbmRleE9mKCcuJykpe1xuICAgICAgICAgIHZhciBzcGxpdCA9IGFycmF5W2ldLnNwbGl0KCcuJykgLy8gc3BsaXQgYXQgdGhlIHBvaW50XG4gICAgICAgICAgdmFyIGZpcnN0ID0gW3NwbGl0LnNoaWZ0KCksIHNwbGl0LnNoaWZ0KCldLmpvaW4oJy4nKSAvLyBqb2luIHRoZSBmaXJzdCBudW1iZXIgdG9nZXRoZXJcbiAgICAgICAgICBhcnJheS5zcGxpY2UuYXBwbHkoYXJyYXksIFtpLCAxXS5jb25jYXQoZmlyc3QsIHNwbGl0Lm1hcChmdW5jdGlvbihlbCl7IHJldHVybiAnLicrZWwgfSkpKSAvLyBhZGQgZmlyc3QgYW5kIGFsbCBvdGhlciBlbnRyaWVzIGJhY2sgdG8gYXJyYXlcbiAgICAgICAgfVxuICAgICAgfVxuXG4gICAgfWVsc2V7XG4gICAgICBhcnJheSA9IGFycmF5LnJlZHVjZShmdW5jdGlvbihwcmV2LCBjdXJyKXtcbiAgICAgICAgcmV0dXJuIFtdLmNvbmNhdC5hcHBseShwcmV2LCBjdXJyKVxuICAgICAgfSwgW10pXG4gICAgfVxuXG4gICAgLy8gYXJyYXkgbm93IGlzIGFuIGFycmF5IGNvbnRhaW5pbmcgYWxsIHBhcnRzIG9mIGEgcGF0aCBlLmcuIFsnTScsICcwJywgJzAnLCAnTCcsICczMCcsICczMCcgLi4uXVxuXG4gICAgdmFyIGFyciA9IFtdXG5cbiAgICBkb3tcblxuICAgICAgLy8gVGVzdCBpZiB3ZSBoYXZlIGEgcGF0aCBsZXR0ZXJcbiAgICAgIGlmKFNWRy5yZWdleC5pc1BhdGhMZXR0ZXIudGVzdChhcnJheVswXSkpe1xuICAgICAgICBzID0gYXJyYXlbMF1cbiAgICAgICAgYXJyYXkuc2hpZnQoKVxuICAgICAgLy8gSWYgbGFzdCBsZXR0ZXIgd2FzIGEgbW92ZSBjb21tYW5kIGFuZCB3ZSBnb3Qgbm8gbmV3LCBpdCBkZWZhdWx0cyB0byBbTF1pbmVcbiAgICAgIH1lbHNlIGlmKHMgPT0gJ00nKXtcbiAgICAgICAgcyA9ICdMJ1xuICAgICAgfWVsc2UgaWYocyA9PSAnbScpe1xuICAgICAgICBzID0gJ2wnXG4gICAgICB9XG5cbiAgICAgIC8vIGFkZCBwYXRoIGxldHRlciBhcyBmaXJzdCBlbGVtZW50XG4gICAgICBzZWcgPSBbcy50b1VwcGVyQ2FzZSgpXVxuXG4gICAgICAvLyBwdXNoIGFsbCBuZWNlc3NhcnkgcGFyYW1ldGVycyB0byBzZWdtZW50XG4gICAgICBmb3IoaSA9IDA7IGkgPCBwYXJhbUNudFtzZWdbMF1dOyArK2kpe1xuICAgICAgICBzZWcucHVzaChwYXJzZUZsb2F0KGFycmF5LnNoaWZ0KCkpKVxuICAgICAgfVxuXG4gICAgICAvLyB1cHBlciBjYXNlXG4gICAgICBpZihzID09IHNlZ1swXSl7XG5cbiAgICAgICAgaWYocyA9PSAnTScgfHwgcyA9PSAnTCcgfHwgcyA9PSAnQycgfHwgcyA9PSAnUScgfHwgcyA9PSAnUycgfHwgcyA9PSAnVCcpe1xuICAgICAgICAgIHggPSBzZWdbcGFyYW1DbnRbc2VnWzBdXS0xXVxuICAgICAgICAgIHkgPSBzZWdbcGFyYW1DbnRbc2VnWzBdXV1cbiAgICAgICAgfWVsc2UgaWYocyA9PSAnVicpe1xuICAgICAgICAgIHkgPSBzZWdbMV1cbiAgICAgICAgfWVsc2UgaWYocyA9PSAnSCcpe1xuICAgICAgICAgIHggPSBzZWdbMV1cbiAgICAgICAgfWVsc2UgaWYocyA9PSAnQScpe1xuICAgICAgICAgIHggPSBzZWdbNl1cbiAgICAgICAgICB5ID0gc2VnWzddXG4gICAgICAgIH1cblxuICAgICAgLy8gbG93ZXIgY2FzZVxuICAgICAgfWVsc2V7XG5cbiAgICAgICAgLy8gY29udmVydCByZWxhdGl2ZSB0byBhYnNvbHV0ZSB2YWx1ZXNcbiAgICAgICAgaWYocyA9PSAnbScgfHwgcyA9PSAnbCcgfHwgcyA9PSAnYycgfHwgcyA9PSAncycgfHwgcyA9PSAncScgfHwgcyA9PSAndCcpe1xuXG4gICAgICAgICAgc2VnWzFdICs9IHhcbiAgICAgICAgICBzZWdbMl0gKz0geVxuXG4gICAgICAgICAgaWYoc2VnWzNdICE9IG51bGwpe1xuICAgICAgICAgICAgc2VnWzNdICs9IHhcbiAgICAgICAgICAgIHNlZ1s0XSArPSB5XG4gICAgICAgICAgfVxuXG4gICAgICAgICAgaWYoc2VnWzVdICE9IG51bGwpe1xuICAgICAgICAgICAgc2VnWzVdICs9IHhcbiAgICAgICAgICAgIHNlZ1s2XSArPSB5XG4gICAgICAgICAgfVxuXG4gICAgICAgICAgLy8gbW92ZSBwb2ludGVyXG4gICAgICAgICAgeCA9IHNlZ1twYXJhbUNudFtzZWdbMF1dLTFdXG4gICAgICAgICAgeSA9IHNlZ1twYXJhbUNudFtzZWdbMF1dXVxuXG4gICAgICAgIH1lbHNlIGlmKHMgPT0gJ3YnKXtcbiAgICAgICAgICBzZWdbMV0gKz0geVxuICAgICAgICAgIHkgPSBzZWdbMV1cbiAgICAgICAgfWVsc2UgaWYocyA9PSAnaCcpe1xuICAgICAgICAgIHNlZ1sxXSArPSB4XG4gICAgICAgICAgeCA9IHNlZ1sxXVxuICAgICAgICB9ZWxzZSBpZihzID09ICdhJyl7XG4gICAgICAgICAgc2VnWzZdICs9IHhcbiAgICAgICAgICBzZWdbN10gKz0geVxuICAgICAgICAgIHggPSBzZWdbNl1cbiAgICAgICAgICB5ID0gc2VnWzddXG4gICAgICAgIH1cblxuICAgICAgfVxuXG4gICAgICBpZihzZWdbMF0gPT0gJ00nKXtcbiAgICAgICAgeDAgPSB4XG4gICAgICAgIHkwID0geVxuICAgICAgfVxuXG4gICAgICBpZihzZWdbMF0gPT0gJ1onKXtcbiAgICAgICAgeCA9IHgwXG4gICAgICAgIHkgPSB5MFxuICAgICAgfVxuXG4gICAgICBhcnIucHVzaChzZWcpXG5cbiAgICB9d2hpbGUoYXJyYXkubGVuZ3RoKVxuXG4gICAgcmV0dXJuIGFyclxuXG4gIH1cbiAgLy8gR2V0IGJvdW5kaW5nIGJveCBvZiBwYXRoXG4sIGJib3g6IGZ1bmN0aW9uKCkge1xuICAgIFNWRy5wYXJzZXIucGF0aC5zZXRBdHRyaWJ1dGUoJ2QnLCB0aGlzLnRvU3RyaW5nKCkpXG5cbiAgICByZXR1cm4gU1ZHLnBhcnNlci5wYXRoLmdldEJCb3goKVxuICB9XG5cbn0pXG5cbi8vIE1vZHVsZSBmb3IgdW5pdCBjb252ZXJ0aW9uc1xuU1ZHLk51bWJlciA9IFNWRy5pbnZlbnQoe1xuICAvLyBJbml0aWFsaXplXG4gIGNyZWF0ZTogZnVuY3Rpb24odmFsdWUsIHVuaXQpIHtcbiAgICAvLyBpbml0aWFsaXplIGRlZmF1bHRzXG4gICAgdGhpcy52YWx1ZSA9IDBcbiAgICB0aGlzLnVuaXQgID0gdW5pdCB8fCAnJ1xuXG4gICAgLy8gcGFyc2UgdmFsdWVcbiAgICBpZiAodHlwZW9mIHZhbHVlID09PSAnbnVtYmVyJykge1xuICAgICAgLy8gZW5zdXJlIGEgdmFsaWQgbnVtZXJpYyB2YWx1ZVxuICAgICAgdGhpcy52YWx1ZSA9IGlzTmFOKHZhbHVlKSA/IDAgOiAhaXNGaW5pdGUodmFsdWUpID8gKHZhbHVlIDwgMCA/IC0zLjRlKzM4IDogKzMuNGUrMzgpIDogdmFsdWVcblxuICAgIH0gZWxzZSBpZiAodHlwZW9mIHZhbHVlID09PSAnc3RyaW5nJykge1xuICAgICAgdW5pdCA9IHZhbHVlLm1hdGNoKFNWRy5yZWdleC5udW1iZXJBbmRVbml0KVxuXG4gICAgICBpZiAodW5pdCkge1xuICAgICAgICAvLyBtYWtlIHZhbHVlIG51bWVyaWNcbiAgICAgICAgdGhpcy52YWx1ZSA9IHBhcnNlRmxvYXQodW5pdFsxXSlcblxuICAgICAgICAvLyBub3JtYWxpemVcbiAgICAgICAgaWYgKHVuaXRbNV0gPT0gJyUnKVxuICAgICAgICAgIHRoaXMudmFsdWUgLz0gMTAwXG4gICAgICAgIGVsc2UgaWYgKHVuaXRbNV0gPT0gJ3MnKVxuICAgICAgICAgIHRoaXMudmFsdWUgKj0gMTAwMFxuXG4gICAgICAgIC8vIHN0b3JlIHVuaXRcbiAgICAgICAgdGhpcy51bml0ID0gdW5pdFs1XVxuICAgICAgfVxuXG4gICAgfSBlbHNlIHtcbiAgICAgIGlmICh2YWx1ZSBpbnN0YW5jZW9mIFNWRy5OdW1iZXIpIHtcbiAgICAgICAgdGhpcy52YWx1ZSA9IHZhbHVlLnZhbHVlT2YoKVxuICAgICAgICB0aGlzLnVuaXQgID0gdmFsdWUudW5pdFxuICAgICAgfVxuICAgIH1cblxuICB9XG4gIC8vIEFkZCBtZXRob2RzXG4sIGV4dGVuZDoge1xuICAgIC8vIFN0cmluZ2FsaXplXG4gICAgdG9TdHJpbmc6IGZ1bmN0aW9uKCkge1xuICAgICAgcmV0dXJuIChcbiAgICAgICAgdGhpcy51bml0ID09ICclJyA/XG4gICAgICAgICAgfn4odGhpcy52YWx1ZSAqIDFlOCkgLyAxZTY6XG4gICAgICAgIHRoaXMudW5pdCA9PSAncycgP1xuICAgICAgICAgIHRoaXMudmFsdWUgLyAxZTMgOlxuICAgICAgICAgIHRoaXMudmFsdWVcbiAgICAgICkgKyB0aGlzLnVuaXRcbiAgICB9XG4gICwgdG9KU09OOiBmdW5jdGlvbigpIHtcbiAgICAgIHJldHVybiB0aGlzLnRvU3RyaW5nKClcbiAgICB9XG4gICwgLy8gQ29udmVydCB0byBwcmltaXRpdmVcbiAgICB2YWx1ZU9mOiBmdW5jdGlvbigpIHtcbiAgICAgIHJldHVybiB0aGlzLnZhbHVlXG4gICAgfVxuICAgIC8vIEFkZCBudW1iZXJcbiAgLCBwbHVzOiBmdW5jdGlvbihudW1iZXIpIHtcbiAgICAgIHJldHVybiBuZXcgU1ZHLk51bWJlcih0aGlzICsgbmV3IFNWRy5OdW1iZXIobnVtYmVyKSwgdGhpcy51bml0KVxuICAgIH1cbiAgICAvLyBTdWJ0cmFjdCBudW1iZXJcbiAgLCBtaW51czogZnVuY3Rpb24obnVtYmVyKSB7XG4gICAgICByZXR1cm4gdGhpcy5wbHVzKC1uZXcgU1ZHLk51bWJlcihudW1iZXIpKVxuICAgIH1cbiAgICAvLyBNdWx0aXBseSBudW1iZXJcbiAgLCB0aW1lczogZnVuY3Rpb24obnVtYmVyKSB7XG4gICAgICByZXR1cm4gbmV3IFNWRy5OdW1iZXIodGhpcyAqIG5ldyBTVkcuTnVtYmVyKG51bWJlciksIHRoaXMudW5pdClcbiAgICB9XG4gICAgLy8gRGl2aWRlIG51bWJlclxuICAsIGRpdmlkZTogZnVuY3Rpb24obnVtYmVyKSB7XG4gICAgICByZXR1cm4gbmV3IFNWRy5OdW1iZXIodGhpcyAvIG5ldyBTVkcuTnVtYmVyKG51bWJlciksIHRoaXMudW5pdClcbiAgICB9XG4gICAgLy8gQ29udmVydCB0byBkaWZmZXJlbnQgdW5pdFxuICAsIHRvOiBmdW5jdGlvbih1bml0KSB7XG4gICAgICB2YXIgbnVtYmVyID0gbmV3IFNWRy5OdW1iZXIodGhpcylcblxuICAgICAgaWYgKHR5cGVvZiB1bml0ID09PSAnc3RyaW5nJylcbiAgICAgICAgbnVtYmVyLnVuaXQgPSB1bml0XG5cbiAgICAgIHJldHVybiBudW1iZXJcbiAgICB9XG4gICAgLy8gTWFrZSBudW1iZXIgbW9ycGhhYmxlXG4gICwgbW9ycGg6IGZ1bmN0aW9uKG51bWJlcikge1xuICAgICAgdGhpcy5kZXN0aW5hdGlvbiA9IG5ldyBTVkcuTnVtYmVyKG51bWJlcilcblxuICAgICAgcmV0dXJuIHRoaXNcbiAgICB9XG4gICAgLy8gR2V0IG1vcnBoZWQgbnVtYmVyIGF0IGdpdmVuIHBvc2l0aW9uXG4gICwgYXQ6IGZ1bmN0aW9uKHBvcykge1xuICAgICAgLy8gTWFrZSBzdXJlIGEgZGVzdGluYXRpb24gaXMgZGVmaW5lZFxuICAgICAgaWYgKCF0aGlzLmRlc3RpbmF0aW9uKSByZXR1cm4gdGhpc1xuXG4gICAgICAvLyBHZW5lcmF0ZSBuZXcgbW9ycGhlZCBudW1iZXJcbiAgICAgIHJldHVybiBuZXcgU1ZHLk51bWJlcih0aGlzLmRlc3RpbmF0aW9uKVxuICAgICAgICAgIC5taW51cyh0aGlzKVxuICAgICAgICAgIC50aW1lcyhwb3MpXG4gICAgICAgICAgLnBsdXModGhpcylcbiAgICB9XG5cbiAgfVxufSlcblxuU1ZHLkVsZW1lbnQgPSBTVkcuaW52ZW50KHtcbiAgLy8gSW5pdGlhbGl6ZSBub2RlXG4gIGNyZWF0ZTogZnVuY3Rpb24obm9kZSkge1xuICAgIC8vIG1ha2Ugc3Ryb2tlIHZhbHVlIGFjY2Vzc2libGUgZHluYW1pY2FsbHlcbiAgICB0aGlzLl9zdHJva2UgPSBTVkcuZGVmYXVsdHMuYXR0cnMuc3Ryb2tlXG5cbiAgICAvLyBpbml0aWFsaXplIGRhdGEgb2JqZWN0XG4gICAgdGhpcy5kb20gPSB7fVxuXG4gICAgLy8gY3JlYXRlIGNpcmN1bGFyIHJlZmVyZW5jZVxuICAgIGlmICh0aGlzLm5vZGUgPSBub2RlKSB7XG4gICAgICB0aGlzLnR5cGUgPSBub2RlLm5vZGVOYW1lXG4gICAgICB0aGlzLm5vZGUuaW5zdGFuY2UgPSB0aGlzXG5cbiAgICAgIC8vIHN0b3JlIGN1cnJlbnQgYXR0cmlidXRlIHZhbHVlXG4gICAgICB0aGlzLl9zdHJva2UgPSBub2RlLmdldEF0dHJpYnV0ZSgnc3Ryb2tlJykgfHwgdGhpcy5fc3Ryb2tlXG4gICAgfVxuICB9XG5cbiAgLy8gQWRkIGNsYXNzIG1ldGhvZHNcbiwgZXh0ZW5kOiB7XG4gICAgLy8gTW92ZSBvdmVyIHgtYXhpc1xuICAgIHg6IGZ1bmN0aW9uKHgpIHtcbiAgICAgIHJldHVybiB0aGlzLmF0dHIoJ3gnLCB4KVxuICAgIH1cbiAgICAvLyBNb3ZlIG92ZXIgeS1heGlzXG4gICwgeTogZnVuY3Rpb24oeSkge1xuICAgICAgcmV0dXJuIHRoaXMuYXR0cigneScsIHkpXG4gICAgfVxuICAgIC8vIE1vdmUgYnkgY2VudGVyIG92ZXIgeC1heGlzXG4gICwgY3g6IGZ1bmN0aW9uKHgpIHtcbiAgICAgIHJldHVybiB4ID09IG51bGwgPyB0aGlzLngoKSArIHRoaXMud2lkdGgoKSAvIDIgOiB0aGlzLngoeCAtIHRoaXMud2lkdGgoKSAvIDIpXG4gICAgfVxuICAgIC8vIE1vdmUgYnkgY2VudGVyIG92ZXIgeS1heGlzXG4gICwgY3k6IGZ1bmN0aW9uKHkpIHtcbiAgICAgIHJldHVybiB5ID09IG51bGwgPyB0aGlzLnkoKSArIHRoaXMuaGVpZ2h0KCkgLyAyIDogdGhpcy55KHkgLSB0aGlzLmhlaWdodCgpIC8gMilcbiAgICB9XG4gICAgLy8gTW92ZSBlbGVtZW50IHRvIGdpdmVuIHggYW5kIHkgdmFsdWVzXG4gICwgbW92ZTogZnVuY3Rpb24oeCwgeSkge1xuICAgICAgcmV0dXJuIHRoaXMueCh4KS55KHkpXG4gICAgfVxuICAgIC8vIE1vdmUgZWxlbWVudCBieSBpdHMgY2VudGVyXG4gICwgY2VudGVyOiBmdW5jdGlvbih4LCB5KSB7XG4gICAgICByZXR1cm4gdGhpcy5jeCh4KS5jeSh5KVxuICAgIH1cbiAgICAvLyBTZXQgd2lkdGggb2YgZWxlbWVudFxuICAsIHdpZHRoOiBmdW5jdGlvbih3aWR0aCkge1xuICAgICAgcmV0dXJuIHRoaXMuYXR0cignd2lkdGgnLCB3aWR0aClcbiAgICB9XG4gICAgLy8gU2V0IGhlaWdodCBvZiBlbGVtZW50XG4gICwgaGVpZ2h0OiBmdW5jdGlvbihoZWlnaHQpIHtcbiAgICAgIHJldHVybiB0aGlzLmF0dHIoJ2hlaWdodCcsIGhlaWdodClcbiAgICB9XG4gICAgLy8gU2V0IGVsZW1lbnQgc2l6ZSB0byBnaXZlbiB3aWR0aCBhbmQgaGVpZ2h0XG4gICwgc2l6ZTogZnVuY3Rpb24od2lkdGgsIGhlaWdodCkge1xuICAgICAgdmFyIHAgPSBwcm9wb3J0aW9uYWxTaXplKHRoaXMsIHdpZHRoLCBoZWlnaHQpXG5cbiAgICAgIHJldHVybiB0aGlzXG4gICAgICAgIC53aWR0aChuZXcgU1ZHLk51bWJlcihwLndpZHRoKSlcbiAgICAgICAgLmhlaWdodChuZXcgU1ZHLk51bWJlcihwLmhlaWdodCkpXG4gICAgfVxuICAgIC8vIENsb25lIGVsZW1lbnRcbiAgLCBjbG9uZTogZnVuY3Rpb24ocGFyZW50KSB7XG4gICAgICAvLyBjbG9uZSBlbGVtZW50IGFuZCBhc3NpZ24gbmV3IGlkXG4gICAgICB2YXIgY2xvbmUgPSBhc3NpZ25OZXdJZCh0aGlzLm5vZGUuY2xvbmVOb2RlKHRydWUpKVxuXG4gICAgICAvLyBpbnNlcnQgdGhlIGNsb25lIGluIHRoZSBnaXZlbiBwYXJlbnQgb3IgYWZ0ZXIgbXlzZWxmXG4gICAgICBpZihwYXJlbnQpIHBhcmVudC5hZGQoY2xvbmUpXG4gICAgICBlbHNlIHRoaXMuYWZ0ZXIoY2xvbmUpXG5cbiAgICAgIHJldHVybiBjbG9uZVxuICAgIH1cbiAgICAvLyBSZW1vdmUgZWxlbWVudFxuICAsIHJlbW92ZTogZnVuY3Rpb24oKSB7XG4gICAgICBpZiAodGhpcy5wYXJlbnQoKSlcbiAgICAgICAgdGhpcy5wYXJlbnQoKS5yZW1vdmVFbGVtZW50KHRoaXMpXG5cbiAgICAgIHJldHVybiB0aGlzXG4gICAgfVxuICAgIC8vIFJlcGxhY2UgZWxlbWVudFxuICAsIHJlcGxhY2U6IGZ1bmN0aW9uKGVsZW1lbnQpIHtcbiAgICAgIHRoaXMuYWZ0ZXIoZWxlbWVudCkucmVtb3ZlKClcblxuICAgICAgcmV0dXJuIGVsZW1lbnRcbiAgICB9XG4gICAgLy8gQWRkIGVsZW1lbnQgdG8gZ2l2ZW4gY29udGFpbmVyIGFuZCByZXR1cm4gc2VsZlxuICAsIGFkZFRvOiBmdW5jdGlvbihwYXJlbnQpIHtcbiAgICAgIHJldHVybiBwYXJlbnQucHV0KHRoaXMpXG4gICAgfVxuICAgIC8vIEFkZCBlbGVtZW50IHRvIGdpdmVuIGNvbnRhaW5lciBhbmQgcmV0dXJuIGNvbnRhaW5lclxuICAsIHB1dEluOiBmdW5jdGlvbihwYXJlbnQpIHtcbiAgICAgIHJldHVybiBwYXJlbnQuYWRkKHRoaXMpXG4gICAgfVxuICAgIC8vIEdldCAvIHNldCBpZFxuICAsIGlkOiBmdW5jdGlvbihpZCkge1xuICAgICAgcmV0dXJuIHRoaXMuYXR0cignaWQnLCBpZClcbiAgICB9XG4gICAgLy8gQ2hlY2tzIHdoZXRoZXIgdGhlIGdpdmVuIHBvaW50IGluc2lkZSB0aGUgYm91bmRpbmcgYm94IG9mIHRoZSBlbGVtZW50XG4gICwgaW5zaWRlOiBmdW5jdGlvbih4LCB5KSB7XG4gICAgICB2YXIgYm94ID0gdGhpcy5iYm94KClcblxuICAgICAgcmV0dXJuIHggPiBib3gueFxuICAgICAgICAgICYmIHkgPiBib3gueVxuICAgICAgICAgICYmIHggPCBib3gueCArIGJveC53aWR0aFxuICAgICAgICAgICYmIHkgPCBib3gueSArIGJveC5oZWlnaHRcbiAgICB9XG4gICAgLy8gU2hvdyBlbGVtZW50XG4gICwgc2hvdzogZnVuY3Rpb24oKSB7XG4gICAgICByZXR1cm4gdGhpcy5zdHlsZSgnZGlzcGxheScsICcnKVxuICAgIH1cbiAgICAvLyBIaWRlIGVsZW1lbnRcbiAgLCBoaWRlOiBmdW5jdGlvbigpIHtcbiAgICAgIHJldHVybiB0aGlzLnN0eWxlKCdkaXNwbGF5JywgJ25vbmUnKVxuICAgIH1cbiAgICAvLyBJcyBlbGVtZW50IHZpc2libGU/XG4gICwgdmlzaWJsZTogZnVuY3Rpb24oKSB7XG4gICAgICByZXR1cm4gdGhpcy5zdHlsZSgnZGlzcGxheScpICE9ICdub25lJ1xuICAgIH1cbiAgICAvLyBSZXR1cm4gaWQgb24gc3RyaW5nIGNvbnZlcnNpb25cbiAgLCB0b1N0cmluZzogZnVuY3Rpb24oKSB7XG4gICAgICByZXR1cm4gdGhpcy5hdHRyKCdpZCcpXG4gICAgfVxuICAgIC8vIFJldHVybiBhcnJheSBvZiBjbGFzc2VzIG9uIHRoZSBub2RlXG4gICwgY2xhc3NlczogZnVuY3Rpb24oKSB7XG4gICAgICB2YXIgYXR0ciA9IHRoaXMuYXR0cignY2xhc3MnKVxuXG4gICAgICByZXR1cm4gYXR0ciA9PSBudWxsID8gW10gOiBhdHRyLnRyaW0oKS5zcGxpdCgvXFxzKy8pXG4gICAgfVxuICAgIC8vIFJldHVybiB0cnVlIGlmIGNsYXNzIGV4aXN0cyBvbiB0aGUgbm9kZSwgZmFsc2Ugb3RoZXJ3aXNlXG4gICwgaGFzQ2xhc3M6IGZ1bmN0aW9uKG5hbWUpIHtcbiAgICAgIHJldHVybiB0aGlzLmNsYXNzZXMoKS5pbmRleE9mKG5hbWUpICE9IC0xXG4gICAgfVxuICAgIC8vIEFkZCBjbGFzcyB0byB0aGUgbm9kZVxuICAsIGFkZENsYXNzOiBmdW5jdGlvbihuYW1lKSB7XG4gICAgICBpZiAoIXRoaXMuaGFzQ2xhc3MobmFtZSkpIHtcbiAgICAgICAgdmFyIGFycmF5ID0gdGhpcy5jbGFzc2VzKClcbiAgICAgICAgYXJyYXkucHVzaChuYW1lKVxuICAgICAgICB0aGlzLmF0dHIoJ2NsYXNzJywgYXJyYXkuam9pbignICcpKVxuICAgICAgfVxuXG4gICAgICByZXR1cm4gdGhpc1xuICAgIH1cbiAgICAvLyBSZW1vdmUgY2xhc3MgZnJvbSB0aGUgbm9kZVxuICAsIHJlbW92ZUNsYXNzOiBmdW5jdGlvbihuYW1lKSB7XG4gICAgICBpZiAodGhpcy5oYXNDbGFzcyhuYW1lKSkge1xuICAgICAgICB0aGlzLmF0dHIoJ2NsYXNzJywgdGhpcy5jbGFzc2VzKCkuZmlsdGVyKGZ1bmN0aW9uKGMpIHtcbiAgICAgICAgICByZXR1cm4gYyAhPSBuYW1lXG4gICAgICAgIH0pLmpvaW4oJyAnKSlcbiAgICAgIH1cblxuICAgICAgcmV0dXJuIHRoaXNcbiAgICB9XG4gICAgLy8gVG9nZ2xlIHRoZSBwcmVzZW5jZSBvZiBhIGNsYXNzIG9uIHRoZSBub2RlXG4gICwgdG9nZ2xlQ2xhc3M6IGZ1bmN0aW9uKG5hbWUpIHtcbiAgICAgIHJldHVybiB0aGlzLmhhc0NsYXNzKG5hbWUpID8gdGhpcy5yZW1vdmVDbGFzcyhuYW1lKSA6IHRoaXMuYWRkQ2xhc3MobmFtZSlcbiAgICB9XG4gICAgLy8gR2V0IHJlZmVyZW5jZWQgZWxlbWVudCBmb3JtIGF0dHJpYnV0ZSB2YWx1ZVxuICAsIHJlZmVyZW5jZTogZnVuY3Rpb24oYXR0cikge1xuICAgICAgcmV0dXJuIFNWRy5nZXQodGhpcy5hdHRyKGF0dHIpKVxuICAgIH1cbiAgICAvLyBSZXR1cm5zIHRoZSBwYXJlbnQgZWxlbWVudCBpbnN0YW5jZVxuICAsIHBhcmVudDogZnVuY3Rpb24odHlwZSkge1xuICAgICAgdmFyIHBhcmVudCA9IHRoaXNcblxuICAgICAgLy8gY2hlY2sgZm9yIHBhcmVudFxuICAgICAgaWYoIXBhcmVudC5ub2RlLnBhcmVudE5vZGUpIHJldHVybiBudWxsXG5cbiAgICAgIC8vIGdldCBwYXJlbnQgZWxlbWVudFxuICAgICAgcGFyZW50ID0gU1ZHLmFkb3B0KHBhcmVudC5ub2RlLnBhcmVudE5vZGUpXG5cbiAgICAgIGlmKCF0eXBlKSByZXR1cm4gcGFyZW50XG5cbiAgICAgIC8vIGxvb3AgdHJvdWdoIGFuY2VzdG9ycyBpZiB0eXBlIGlzIGdpdmVuXG4gICAgICB3aGlsZShwYXJlbnQgJiYgcGFyZW50Lm5vZGUgaW5zdGFuY2VvZiBTVkdFbGVtZW50KXtcbiAgICAgICAgaWYodHlwZW9mIHR5cGUgPT09ICdzdHJpbmcnID8gcGFyZW50Lm1hdGNoZXModHlwZSkgOiBwYXJlbnQgaW5zdGFuY2VvZiB0eXBlKSByZXR1cm4gcGFyZW50XG4gICAgICAgIHBhcmVudCA9IFNWRy5hZG9wdChwYXJlbnQubm9kZS5wYXJlbnROb2RlKVxuICAgICAgfVxuICAgIH1cbiAgICAvLyBHZXQgcGFyZW50IGRvY3VtZW50XG4gICwgZG9jOiBmdW5jdGlvbigpIHtcbiAgICAgIHJldHVybiB0aGlzIGluc3RhbmNlb2YgU1ZHLkRvYyA/IHRoaXMgOiB0aGlzLnBhcmVudChTVkcuRG9jKVxuICAgIH1cbiAgICAvLyByZXR1cm4gYXJyYXkgb2YgYWxsIGFuY2VzdG9ycyBvZiBnaXZlbiB0eXBlIHVwIHRvIHRoZSByb290IHN2Z1xuICAsIHBhcmVudHM6IGZ1bmN0aW9uKHR5cGUpIHtcbiAgICAgIHZhciBwYXJlbnRzID0gW10sIHBhcmVudCA9IHRoaXNcblxuICAgICAgZG97XG4gICAgICAgIHBhcmVudCA9IHBhcmVudC5wYXJlbnQodHlwZSlcbiAgICAgICAgaWYoIXBhcmVudCB8fCAhcGFyZW50Lm5vZGUpIGJyZWFrXG5cbiAgICAgICAgcGFyZW50cy5wdXNoKHBhcmVudClcbiAgICAgIH0gd2hpbGUocGFyZW50LnBhcmVudClcblxuICAgICAgcmV0dXJuIHBhcmVudHNcbiAgICB9XG4gICAgLy8gbWF0Y2hlcyB0aGUgZWxlbWVudCB2cyBhIGNzcyBzZWxlY3RvclxuICAsIG1hdGNoZXM6IGZ1bmN0aW9uKHNlbGVjdG9yKXtcbiAgICAgIHJldHVybiBtYXRjaGVzKHRoaXMubm9kZSwgc2VsZWN0b3IpXG4gICAgfVxuICAgIC8vIFJldHVybnMgdGhlIHN2ZyBub2RlIHRvIGNhbGwgbmF0aXZlIHN2ZyBtZXRob2RzIG9uIGl0XG4gICwgbmF0aXZlOiBmdW5jdGlvbigpIHtcbiAgICAgIHJldHVybiB0aGlzLm5vZGVcbiAgICB9XG4gICAgLy8gSW1wb3J0IHJhdyBzdmdcbiAgLCBzdmc6IGZ1bmN0aW9uKHN2Zykge1xuICAgICAgLy8gY3JlYXRlIHRlbXBvcmFyeSBob2xkZXJcbiAgICAgIHZhciB3ZWxsID0gZG9jdW1lbnQuY3JlYXRlRWxlbWVudCgnc3ZnJylcblxuICAgICAgLy8gYWN0IGFzIGEgc2V0dGVyIGlmIHN2ZyBpcyBnaXZlblxuICAgICAgaWYgKHN2ZyAmJiB0aGlzIGluc3RhbmNlb2YgU1ZHLlBhcmVudCkge1xuICAgICAgICAvLyBkdW1wIHJhdyBzdmdcbiAgICAgICAgd2VsbC5pbm5lckhUTUwgPSAnPHN2Zz4nICsgc3ZnLnJlcGxhY2UoL1xcbi8sICcnKS5yZXBsYWNlKC88KFxcdyspKFtePF0rPylcXC8+L2csICc8JDEkMj48LyQxPicpICsgJzwvc3ZnPidcblxuICAgICAgICAvLyB0cmFuc3BsYW50IG5vZGVzXG4gICAgICAgIGZvciAodmFyIGkgPSAwLCBpbCA9IHdlbGwuZmlyc3RDaGlsZC5jaGlsZE5vZGVzLmxlbmd0aDsgaSA8IGlsOyBpKyspXG4gICAgICAgICAgdGhpcy5ub2RlLmFwcGVuZENoaWxkKHdlbGwuZmlyc3RDaGlsZC5maXJzdENoaWxkKVxuXG4gICAgICAvLyBvdGhlcndpc2UgYWN0IGFzIGEgZ2V0dGVyXG4gICAgICB9IGVsc2Uge1xuICAgICAgICAvLyBjcmVhdGUgYSB3cmFwcGluZyBzdmcgZWxlbWVudCBpbiBjYXNlIG9mIHBhcnRpYWwgY29udGVudFxuICAgICAgICB3ZWxsLmFwcGVuZENoaWxkKHN2ZyA9IGRvY3VtZW50LmNyZWF0ZUVsZW1lbnQoJ3N2ZycpKVxuXG4gICAgICAgIC8vIHdyaXRlIHN2Z2pzIGRhdGEgdG8gdGhlIGRvbVxuICAgICAgICB0aGlzLndyaXRlRGF0YVRvRG9tKClcblxuICAgICAgICAvLyBpbnNlcnQgYSBjb3B5IG9mIHRoaXMgbm9kZVxuICAgICAgICBzdmcuYXBwZW5kQ2hpbGQodGhpcy5ub2RlLmNsb25lTm9kZSh0cnVlKSlcblxuICAgICAgICAvLyByZXR1cm4gdGFyZ2V0IGVsZW1lbnRcbiAgICAgICAgcmV0dXJuIHdlbGwuaW5uZXJIVE1MLnJlcGxhY2UoL148c3ZnPi8sICcnKS5yZXBsYWNlKC88XFwvc3ZnPiQvLCAnJylcbiAgICAgIH1cblxuICAgICAgcmV0dXJuIHRoaXNcbiAgICB9XG4gIC8vIHdyaXRlIHN2Z2pzIGRhdGEgdG8gdGhlIGRvbVxuICAsIHdyaXRlRGF0YVRvRG9tOiBmdW5jdGlvbigpIHtcblxuICAgICAgLy8gZHVtcCB2YXJpYWJsZXMgcmVjdXJzaXZlbHlcbiAgICAgIGlmKHRoaXMuZWFjaCB8fCB0aGlzLmxpbmVzKXtcbiAgICAgICAgdmFyIGZuID0gdGhpcy5lYWNoID8gdGhpcyA6IHRoaXMubGluZXMoKTtcbiAgICAgICAgZm4uZWFjaChmdW5jdGlvbigpe1xuICAgICAgICAgIHRoaXMud3JpdGVEYXRhVG9Eb20oKVxuICAgICAgICB9KVxuICAgICAgfVxuXG4gICAgICAvLyByZW1vdmUgcHJldmlvdXNseSBzZXQgZGF0YVxuICAgICAgdGhpcy5ub2RlLnJlbW92ZUF0dHJpYnV0ZSgnc3ZnanM6ZGF0YScpXG5cbiAgICAgIGlmKE9iamVjdC5rZXlzKHRoaXMuZG9tKS5sZW5ndGgpXG4gICAgICAgIHRoaXMubm9kZS5zZXRBdHRyaWJ1dGUoJ3N2Z2pzOmRhdGEnLCBKU09OLnN0cmluZ2lmeSh0aGlzLmRvbSkpIC8vIHNlZSAjNDI4XG5cbiAgICAgIHJldHVybiB0aGlzXG4gICAgfVxuICAvLyBzZXQgZ2l2ZW4gZGF0YSB0byB0aGUgZWxlbWVudHMgZGF0YSBwcm9wZXJ0eVxuICAsIHNldERhdGE6IGZ1bmN0aW9uKG8pe1xuICAgICAgdGhpcy5kb20gPSBvXG4gICAgICByZXR1cm4gdGhpc1xuICAgIH1cbiAgLCBpczogZnVuY3Rpb24ob2JqKXtcbiAgICAgIHJldHVybiBpcyh0aGlzLCBvYmopXG4gICAgfVxuICB9XG59KVxuXG5TVkcuZWFzaW5nID0ge1xuICAnLSc6IGZ1bmN0aW9uKHBvcyl7cmV0dXJuIHBvc31cbiwgJzw+JzpmdW5jdGlvbihwb3Mpe3JldHVybiAtTWF0aC5jb3MocG9zICogTWF0aC5QSSkgLyAyICsgMC41fVxuLCAnPic6IGZ1bmN0aW9uKHBvcyl7cmV0dXJuICBNYXRoLnNpbihwb3MgKiBNYXRoLlBJIC8gMil9XG4sICc8JzogZnVuY3Rpb24ocG9zKXtyZXR1cm4gLU1hdGguY29zKHBvcyAqIE1hdGguUEkgLyAyKSArIDF9XG59XG5cblNWRy5tb3JwaCA9IGZ1bmN0aW9uKHBvcyl7XG4gIHJldHVybiBmdW5jdGlvbihmcm9tLCB0bykge1xuICAgIHJldHVybiBuZXcgU1ZHLk1vcnBoT2JqKGZyb20sIHRvKS5hdChwb3MpXG4gIH1cbn1cblxuU1ZHLlNpdHVhdGlvbiA9IFNWRy5pbnZlbnQoe1xuXG4gIGNyZWF0ZTogZnVuY3Rpb24obyl7XG4gICAgdGhpcy5pbml0ID0gZmFsc2VcbiAgICB0aGlzLnJldmVyc2VkID0gZmFsc2VcbiAgICB0aGlzLnJldmVyc2luZyA9IGZhbHNlXG5cbiAgICB0aGlzLmR1cmF0aW9uID0gbmV3IFNWRy5OdW1iZXIoby5kdXJhdGlvbikudmFsdWVPZigpXG4gICAgdGhpcy5kZWxheSA9IG5ldyBTVkcuTnVtYmVyKG8uZGVsYXkpLnZhbHVlT2YoKVxuXG4gICAgdGhpcy5zdGFydCA9ICtuZXcgRGF0ZSgpICsgdGhpcy5kZWxheVxuICAgIHRoaXMuZmluaXNoID0gdGhpcy5zdGFydCArIHRoaXMuZHVyYXRpb25cbiAgICB0aGlzLmVhc2UgPSBvLmVhc2VcblxuICAgIC8vIHRoaXMubG9vcCBpcyBpbmNyZW1lbnRlZCBmcm9tIDAgdG8gdGhpcy5sb29wc1xuICAgIC8vIGl0IGlzIGFsc28gaW5jcmVtZW50ZWQgd2hlbiBpbiBhbiBpbmZpbml0ZSBsb29wICh3aGVuIHRoaXMubG9vcHMgaXMgdHJ1ZSlcbiAgICB0aGlzLmxvb3AgPSAwXG4gICAgdGhpcy5sb29wcyA9IGZhbHNlXG5cbiAgICB0aGlzLmFuaW1hdGlvbnMgPSB7XG4gICAgICAvLyBmdW5jdGlvblRvQ2FsbDogW2xpc3Qgb2YgbW9ycGhhYmxlIG9iamVjdHNdXG4gICAgICAvLyBlLmcuIG1vdmU6IFtTVkcuTnVtYmVyLCBTVkcuTnVtYmVyXVxuICAgIH1cblxuICAgIHRoaXMuYXR0cnMgPSB7XG4gICAgICAvLyBob2xkcyBhbGwgYXR0cmlidXRlcyB3aGljaCBhcmUgbm90IHJlcHJlc2VudGVkIGZyb20gYSBmdW5jdGlvbiBzdmcuanMgcHJvdmlkZXNcbiAgICAgIC8vIGUuZy4gc29tZUF0dHI6IFNWRy5OdW1iZXJcbiAgICB9XG5cbiAgICB0aGlzLnN0eWxlcyA9IHtcbiAgICAgIC8vIGhvbGRzIGFsbCBzdHlsZXMgd2hpY2ggc2hvdWxkIGJlIGFuaW1hdGVkXG4gICAgICAvLyBlLmcuIGZpbGwtY29sb3I6IFNWRy5Db2xvclxuICAgIH1cblxuICAgIHRoaXMudHJhbnNmb3JtcyA9IFtcbiAgICAgIC8vIGhvbGRzIGFsbCB0cmFuc2Zvcm1hdGlvbnMgYXMgdHJhbnNmb3JtYXRpb24gb2JqZWN0c1xuICAgICAgLy8gZS5nLiBbU1ZHLlJvdGF0ZSwgU1ZHLlRyYW5zbGF0ZSwgU1ZHLk1hdHJpeF1cbiAgICBdXG5cbiAgICB0aGlzLm9uY2UgPSB7XG4gICAgICAvLyBmdW5jdGlvbnMgdG8gZmlyZSBhdCBhIHNwZWNpZmljIHBvc2l0aW9uXG4gICAgICAvLyBlLmcuIFwiMC41XCI6IGZ1bmN0aW9uIGZvbygpe31cbiAgICB9XG5cbiAgfVxuXG59KVxuXG5cblNWRy5GWCA9IFNWRy5pbnZlbnQoe1xuXG4gIGNyZWF0ZTogZnVuY3Rpb24oZWxlbWVudCkge1xuICAgIHRoaXMuX3RhcmdldCA9IGVsZW1lbnRcbiAgICB0aGlzLnNpdHVhdGlvbnMgPSBbXVxuICAgIHRoaXMuYWN0aXZlID0gZmFsc2VcbiAgICB0aGlzLnNpdHVhdGlvbiA9IG51bGxcbiAgICB0aGlzLnBhdXNlZCA9IGZhbHNlXG4gICAgdGhpcy5sYXN0UG9zID0gMFxuICAgIHRoaXMucG9zID0gMFxuICAgIC8vIFRoZSBhYnNvbHV0ZSBwb3NpdGlvbiBvZiBhbiBhbmltYXRpb24gaXMgaXRzIHBvc2l0aW9uIGluIHRoZSBjb250ZXh0IG9mIGl0cyBjb21wbGV0ZSBkdXJhdGlvbiAoaW5jbHVkaW5nIGRlbGF5IGFuZCBsb29wcylcbiAgICAvLyBXaGVuIHBlcmZvcm1pbmcgYSBkZWxheSwgYWJzUG9zIGlzIGJlbG93IDAgYW5kIHdoZW4gcGVyZm9ybWluZyBhIGxvb3AsIGl0cyB2YWx1ZSBpcyBhYm92ZSAxXG4gICAgdGhpcy5hYnNQb3MgPSAwXG4gICAgdGhpcy5fc3BlZWQgPSAxXG4gIH1cblxuLCBleHRlbmQ6IHtcblxuICAgIC8qKlxuICAgICAqIHNldHMgb3IgcmV0dXJucyB0aGUgdGFyZ2V0IG9mIHRoaXMgYW5pbWF0aW9uXG4gICAgICogQHBhcmFtIG8gb2JqZWN0IHx8IG51bWJlciBJbiBjYXNlIG9mIE9iamVjdCBpdCBob2xkcyBhbGwgcGFyYW1ldGVycy4gSW4gY2FzZSBvZiBudW1iZXIgaXRzIHRoZSBkdXJhdGlvbiBvZiB0aGUgYW5pbWF0aW9uXG4gICAgICogQHBhcmFtIGVhc2UgZnVuY3Rpb24gfHwgc3RyaW5nIEZ1bmN0aW9uIHdoaWNoIHNob3VsZCBiZSB1c2VkIGZvciBlYXNpbmcgb3IgZWFzaW5nIGtleXdvcmRcbiAgICAgKiBAcGFyYW0gZGVsYXkgTnVtYmVyIGluZGljYXRpbmcgdGhlIGRlbGF5IGJlZm9yZSB0aGUgYW5pbWF0aW9uIHN0YXJ0c1xuICAgICAqIEByZXR1cm4gdGFyZ2V0IHx8IHRoaXNcbiAgICAgKi9cbiAgICBhbmltYXRlOiBmdW5jdGlvbihvLCBlYXNlLCBkZWxheSl7XG5cbiAgICAgIGlmKHR5cGVvZiBvID09ICdvYmplY3QnKXtcbiAgICAgICAgZWFzZSA9IG8uZWFzZVxuICAgICAgICBkZWxheSA9IG8uZGVsYXlcbiAgICAgICAgbyA9IG8uZHVyYXRpb25cbiAgICAgIH1cblxuICAgICAgdmFyIHNpdHVhdGlvbiA9IG5ldyBTVkcuU2l0dWF0aW9uKHtcbiAgICAgICAgZHVyYXRpb246IG8gfHwgMTAwMCxcbiAgICAgICAgZGVsYXk6IGRlbGF5IHx8IDAsXG4gICAgICAgIGVhc2U6IFNWRy5lYXNpbmdbZWFzZSB8fCAnLSddIHx8IGVhc2VcbiAgICAgIH0pXG5cbiAgICAgIHRoaXMucXVldWUoc2l0dWF0aW9uKVxuXG4gICAgICByZXR1cm4gdGhpc1xuICAgIH1cblxuICAgIC8qKlxuICAgICAqIHNldHMgYSBkZWxheSBiZWZvcmUgdGhlIG5leHQgZWxlbWVudCBvZiB0aGUgcXVldWUgaXMgY2FsbGVkXG4gICAgICogQHBhcmFtIGRlbGF5IER1cmF0aW9uIG9mIGRlbGF5IGluIG1pbGxpc2Vjb25kc1xuICAgICAqIEByZXR1cm4gdGhpcy50YXJnZXQoKVxuICAgICAqL1xuICAsIGRlbGF5OiBmdW5jdGlvbihkZWxheSl7XG4gICAgICAvLyBUaGUgZGVsYXkgaXMgcGVyZm9ybWVkIGJ5IGFuIGVtcHR5IHNpdHVhdGlvbiB3aXRoIGl0cyBkdXJhdGlvblxuICAgICAgLy8gYXR0cmlidXRlIHNldCB0byB0aGUgZHVyYXRpb24gb2YgdGhlIGRlbGF5XG4gICAgICB2YXIgc2l0dWF0aW9uID0gbmV3IFNWRy5TaXR1YXRpb24oe1xuICAgICAgICBkdXJhdGlvbjogZGVsYXksXG4gICAgICAgIGRlbGF5OiAwLFxuICAgICAgICBlYXNlOiBTVkcuZWFzaW5nWyctJ11cbiAgICAgIH0pXG5cbiAgICAgIHJldHVybiB0aGlzLnF1ZXVlKHNpdHVhdGlvbilcbiAgICB9XG5cbiAgICAvKipcbiAgICAgKiBzZXRzIG9yIHJldHVybnMgdGhlIHRhcmdldCBvZiB0aGlzIGFuaW1hdGlvblxuICAgICAqIEBwYXJhbSBudWxsIHx8IHRhcmdldCBTVkcuRWxlbWVudCB3aGljaCBzaG91bGQgYmUgc2V0IGFzIG5ldyB0YXJnZXRcbiAgICAgKiBAcmV0dXJuIHRhcmdldCB8fCB0aGlzXG4gICAgICovXG4gICwgdGFyZ2V0OiBmdW5jdGlvbih0YXJnZXQpe1xuICAgICAgaWYodGFyZ2V0ICYmIHRhcmdldCBpbnN0YW5jZW9mIFNWRy5FbGVtZW50KXtcbiAgICAgICAgdGhpcy5fdGFyZ2V0ID0gdGFyZ2V0XG4gICAgICAgIHJldHVybiB0aGlzXG4gICAgICB9XG5cbiAgICAgIHJldHVybiB0aGlzLl90YXJnZXRcbiAgICB9XG5cbiAgICAvLyByZXR1cm5zIHRoZSBhYnNvbHV0ZSBwb3NpdGlvbiBhdCBhIGdpdmVuIHRpbWVcbiAgLCB0aW1lVG9BYnNQb3M6IGZ1bmN0aW9uKHRpbWVzdGFtcCl7XG4gICAgICByZXR1cm4gKHRpbWVzdGFtcCAtIHRoaXMuc2l0dWF0aW9uLnN0YXJ0KSAvICh0aGlzLnNpdHVhdGlvbi5kdXJhdGlvbi90aGlzLl9zcGVlZClcbiAgICB9XG5cbiAgICAvLyByZXR1cm5zIHRoZSB0aW1lc3RhbXAgZnJvbSBhIGdpdmVuIGFic29sdXRlIHBvc2l0b25cbiAgLCBhYnNQb3NUb1RpbWU6IGZ1bmN0aW9uKGFic1Bvcyl7XG4gICAgICByZXR1cm4gdGhpcy5zaXR1YXRpb24uZHVyYXRpb24vdGhpcy5fc3BlZWQgKiBhYnNQb3MgKyB0aGlzLnNpdHVhdGlvbi5zdGFydFxuICAgIH1cblxuICAgIC8vIHN0YXJ0cyB0aGUgYW5pbWF0aW9ubG9vcFxuICAsIHN0YXJ0QW5pbUZyYW1lOiBmdW5jdGlvbigpe1xuICAgICAgdGhpcy5zdG9wQW5pbUZyYW1lKClcbiAgICAgIHRoaXMuYW5pbWF0aW9uRnJhbWUgPSByZXF1ZXN0QW5pbWF0aW9uRnJhbWUoZnVuY3Rpb24oKXsgdGhpcy5zdGVwKCkgfS5iaW5kKHRoaXMpKVxuICAgIH1cblxuICAgIC8vIGNhbmNlbHMgdGhlIGFuaW1hdGlvbmZyYW1lXG4gICwgc3RvcEFuaW1GcmFtZTogZnVuY3Rpb24oKXtcbiAgICAgIGNhbmNlbEFuaW1hdGlvbkZyYW1lKHRoaXMuYW5pbWF0aW9uRnJhbWUpXG4gICAgfVxuXG4gICAgLy8ga2lja3Mgb2ZmIHRoZSBhbmltYXRpb24gLSBvbmx5IGRvZXMgc29tZXRoaW5nIHdoZW4gdGhlIHF1ZXVlIGlzIGN1cnJlbnRseSBub3QgYWN0aXZlIGFuZCBhdCBsZWFzdCBvbmUgc2l0dWF0aW9uIGlzIHNldFxuICAsIHN0YXJ0OiBmdW5jdGlvbigpe1xuICAgICAgLy8gZG9udCBzdGFydCBpZiBhbHJlYWR5IHN0YXJ0ZWRcbiAgICAgIGlmKCF0aGlzLmFjdGl2ZSAmJiB0aGlzLnNpdHVhdGlvbil7XG4gICAgICAgIHRoaXMuYWN0aXZlID0gdHJ1ZVxuICAgICAgICB0aGlzLnN0YXJ0Q3VycmVudCgpXG4gICAgICB9XG5cbiAgICAgIHJldHVybiB0aGlzXG4gICAgfVxuXG4gICAgLy8gc3RhcnQgdGhlIGN1cnJlbnQgc2l0dWF0aW9uXG4gICwgc3RhcnRDdXJyZW50OiBmdW5jdGlvbigpe1xuICAgICAgdGhpcy5zaXR1YXRpb24uc3RhcnQgPSArbmV3IERhdGUgKyB0aGlzLnNpdHVhdGlvbi5kZWxheS90aGlzLl9zcGVlZFxuICAgICAgdGhpcy5zaXR1YXRpb24uZmluaXNoID0gdGhpcy5zaXR1YXRpb24uc3RhcnQgKyB0aGlzLnNpdHVhdGlvbi5kdXJhdGlvbi90aGlzLl9zcGVlZFxuICAgICAgcmV0dXJuIHRoaXMuaW5pdEFuaW1hdGlvbnMoKS5zdGVwKClcbiAgICB9XG5cbiAgICAvKipcbiAgICAgKiBhZGRzIGEgZnVuY3Rpb24gLyBTaXR1YXRpb24gdG8gdGhlIGFuaW1hdGlvbiBxdWV1ZVxuICAgICAqIEBwYXJhbSBmbiBmdW5jdGlvbiAvIHNpdHVhdGlvbiB0byBhZGRcbiAgICAgKiBAcmV0dXJuIHRoaXNcbiAgICAgKi9cbiAgLCBxdWV1ZTogZnVuY3Rpb24oZm4pe1xuICAgICAgaWYodHlwZW9mIGZuID09ICdmdW5jdGlvbicgfHwgZm4gaW5zdGFuY2VvZiBTVkcuU2l0dWF0aW9uKVxuICAgICAgICB0aGlzLnNpdHVhdGlvbnMucHVzaChmbilcblxuICAgICAgaWYoIXRoaXMuc2l0dWF0aW9uKSB0aGlzLnNpdHVhdGlvbiA9IHRoaXMuc2l0dWF0aW9ucy5zaGlmdCgpXG5cbiAgICAgIHJldHVybiB0aGlzXG4gICAgfVxuXG4gICAgLyoqXG4gICAgICogcHVsbHMgbmV4dCBlbGVtZW50IGZyb20gdGhlIHF1ZXVlIGFuZCBleGVjdXRlIGl0XG4gICAgICogQHJldHVybiB0aGlzXG4gICAgICovXG4gICwgZGVxdWV1ZTogZnVuY3Rpb24oKXtcbiAgICAgIC8vIHN0b3AgY3VycmVudCBhbmltYXRpb25cbiAgICAgIHRoaXMuc2l0dWF0aW9uICYmIHRoaXMuc2l0dWF0aW9uLnN0b3AgJiYgdGhpcy5zaXR1YXRpb24uc3RvcCgpXG5cbiAgICAgIC8vIGdldCBuZXh0IGFuaW1hdGlvbiBmcm9tIHF1ZXVlXG4gICAgICB0aGlzLnNpdHVhdGlvbiA9IHRoaXMuc2l0dWF0aW9ucy5zaGlmdCgpXG5cbiAgICAgIGlmKHRoaXMuc2l0dWF0aW9uKXtcbiAgICAgICAgaWYodGhpcy5zaXR1YXRpb24gaW5zdGFuY2VvZiBTVkcuU2l0dWF0aW9uKSB7XG4gICAgICAgICAgdGhpcy5zdGFydEN1cnJlbnQoKVxuICAgICAgICB9IGVsc2Uge1xuICAgICAgICAgIC8vIElmIGl0IGlzIG5vdCBhIFNWRy5TaXR1YXRpb24sIHRoZW4gaXQgaXMgYSBmdW5jdGlvbiwgd2UgZXhlY3V0ZSBpdFxuICAgICAgICAgIHRoaXMuc2l0dWF0aW9uLmNhbGwodGhpcylcbiAgICAgICAgfVxuICAgICAgfVxuXG4gICAgICByZXR1cm4gdGhpc1xuICAgIH1cblxuICAgIC8vIHVwZGF0ZXMgYWxsIGFuaW1hdGlvbnMgdG8gdGhlIGN1cnJlbnQgc3RhdGUgb2YgdGhlIGVsZW1lbnRcbiAgICAvLyB0aGlzIGlzIGltcG9ydGFudCB3aGVuIG9uZSBwcm9wZXJ0eSBjb3VsZCBiZSBjaGFuZ2VkIGZyb20gYW5vdGhlciBwcm9wZXJ0eVxuICAsIGluaXRBbmltYXRpb25zOiBmdW5jdGlvbigpIHtcbiAgICAgIHZhciBpXG4gICAgICB2YXIgcyA9IHRoaXMuc2l0dWF0aW9uXG5cbiAgICAgIGlmKHMuaW5pdCkgcmV0dXJuIHRoaXNcblxuICAgICAgZm9yKGkgaW4gcy5hbmltYXRpb25zKXtcblxuICAgICAgICBpZihpID09ICd2aWV3Ym94Jyl7XG4gICAgICAgICAgcy5hbmltYXRpb25zW2ldID0gdGhpcy50YXJnZXQoKS52aWV3Ym94KCkubW9ycGgocy5hbmltYXRpb25zW2ldKVxuICAgICAgICB9ZWxzZXtcblxuICAgICAgICAgIC8vIFRPRE86IHRoaXMgaXMgbm90IGEgY2xlYW4gY2xvbmUgb2YgdGhlIGFycmF5LiBXZSBtYXkgaGF2ZSBzb21lIHVuY2hlY2tlZCByZWZlcmVuY2VzXG4gICAgICAgICAgcy5hbmltYXRpb25zW2ldLnZhbHVlID0gKGkgPT0gJ3Bsb3QnID8gdGhpcy50YXJnZXQoKS5hcnJheSgpLnZhbHVlIDogdGhpcy50YXJnZXQoKVtpXSgpKVxuXG4gICAgICAgICAgLy8gc29tZXRpbWVzIHdlIGdldCBiYWNrIGFuIG9iamVjdCBhbmQgbm90IHRoZSByZWFsIHZhbHVlLCBmaXggdGhpc1xuICAgICAgICAgIGlmKHMuYW5pbWF0aW9uc1tpXS52YWx1ZS52YWx1ZSl7XG4gICAgICAgICAgICBzLmFuaW1hdGlvbnNbaV0udmFsdWUgPSBzLmFuaW1hdGlvbnNbaV0udmFsdWUudmFsdWVcbiAgICAgICAgICB9XG5cbiAgICAgICAgICBpZihzLmFuaW1hdGlvbnNbaV0ucmVsYXRpdmUpXG4gICAgICAgICAgICBzLmFuaW1hdGlvbnNbaV0uZGVzdGluYXRpb24udmFsdWUgPSBzLmFuaW1hdGlvbnNbaV0uZGVzdGluYXRpb24udmFsdWUgKyBzLmFuaW1hdGlvbnNbaV0udmFsdWVcblxuICAgICAgICB9XG5cbiAgICAgIH1cblxuICAgICAgZm9yKGkgaW4gcy5hdHRycyl7XG4gICAgICAgIGlmKHMuYXR0cnNbaV0gaW5zdGFuY2VvZiBTVkcuQ29sb3Ipe1xuICAgICAgICAgIHZhciBjb2xvciA9IG5ldyBTVkcuQ29sb3IodGhpcy50YXJnZXQoKS5hdHRyKGkpKVxuICAgICAgICAgIHMuYXR0cnNbaV0uciA9IGNvbG9yLnJcbiAgICAgICAgICBzLmF0dHJzW2ldLmcgPSBjb2xvci5nXG4gICAgICAgICAgcy5hdHRyc1tpXS5iID0gY29sb3IuYlxuICAgICAgICB9ZWxzZXtcbiAgICAgICAgICBzLmF0dHJzW2ldLnZhbHVlID0gdGhpcy50YXJnZXQoKS5hdHRyKGkpLy8gKyBzLmF0dHJzW2ldLnZhbHVlXG4gICAgICAgIH1cbiAgICAgIH1cblxuICAgICAgZm9yKGkgaW4gcy5zdHlsZXMpe1xuICAgICAgICBzLnN0eWxlc1tpXS52YWx1ZSA9IHRoaXMudGFyZ2V0KCkuc3R5bGUoaSlcbiAgICAgIH1cblxuICAgICAgcy5pbml0aWFsVHJhbnNmb3JtYXRpb24gPSB0aGlzLnRhcmdldCgpLm1hdHJpeGlmeSgpXG5cbiAgICAgIHMuaW5pdCA9IHRydWVcbiAgICAgIHJldHVybiB0aGlzXG4gICAgfVxuICAsIGNsZWFyUXVldWU6IGZ1bmN0aW9uKCl7XG4gICAgICB0aGlzLnNpdHVhdGlvbnMgPSBbXVxuICAgICAgcmV0dXJuIHRoaXNcbiAgICB9XG4gICwgY2xlYXJDdXJyZW50OiBmdW5jdGlvbigpe1xuICAgICAgdGhpcy5zaXR1YXRpb24gPSBudWxsXG4gICAgICByZXR1cm4gdGhpc1xuICAgIH1cbiAgICAvKiogc3RvcHMgdGhlIGFuaW1hdGlvbiBpbW1lZGlhdGVseVxuICAgICAqIEBwYXJhbSBqdW1wVG9FbmQgQSBCb29sZWFuIGluZGljYXRpbmcgd2hldGhlciB0byBjb21wbGV0ZSB0aGUgY3VycmVudCBhbmltYXRpb24gaW1tZWRpYXRlbHkuXG4gICAgICogQHBhcmFtIGNsZWFyUXVldWUgQSBCb29sZWFuIGluZGljYXRpbmcgd2hldGhlciB0byByZW1vdmUgcXVldWVkIGFuaW1hdGlvbiBhcyB3ZWxsLlxuICAgICAqIEByZXR1cm4gdGhpc1xuICAgICAqL1xuICAsIHN0b3A6IGZ1bmN0aW9uKGp1bXBUb0VuZCwgY2xlYXJRdWV1ZSl7XG4gICAgICBpZighdGhpcy5hY3RpdmUpIHRoaXMuc3RhcnQoKVxuXG4gICAgICBpZihjbGVhclF1ZXVlKXtcbiAgICAgICAgdGhpcy5jbGVhclF1ZXVlKClcbiAgICAgIH1cblxuICAgICAgdGhpcy5hY3RpdmUgPSBmYWxzZVxuXG4gICAgICBpZihqdW1wVG9FbmQgJiYgdGhpcy5zaXR1YXRpb24pe1xuICAgICAgICB0aGlzLmF0RW5kKClcbiAgICAgIH1cblxuICAgICAgdGhpcy5zdG9wQW5pbUZyYW1lKClcblxuICAgICAgcmV0dXJuIHRoaXMuY2xlYXJDdXJyZW50KClcbiAgICB9XG5cbiAgICAvKiogcmVzZXRzIHRoZSBlbGVtZW50IHRvIHRoZSBzdGF0ZSB3aGVyZSB0aGUgY3VycmVudCBlbGVtZW50IGhhcyBzdGFydGVkXG4gICAgICogQHJldHVybiB0aGlzXG4gICAgICovXG4gICwgcmVzZXQ6IGZ1bmN0aW9uKCl7XG4gICAgICBpZih0aGlzLnNpdHVhdGlvbil7XG4gICAgICAgIHZhciB0ZW1wID0gdGhpcy5zaXR1YXRpb25cbiAgICAgICAgdGhpcy5zdG9wKClcbiAgICAgICAgdGhpcy5zaXR1YXRpb24gPSB0ZW1wXG4gICAgICAgIHRoaXMuYXRTdGFydCgpXG4gICAgICB9XG4gICAgICByZXR1cm4gdGhpc1xuICAgIH1cblxuICAgIC8vIFN0b3AgdGhlIGN1cnJlbnRseS1ydW5uaW5nIGFuaW1hdGlvbiwgcmVtb3ZlIGFsbCBxdWV1ZWQgYW5pbWF0aW9ucywgYW5kIGNvbXBsZXRlIGFsbCBhbmltYXRpb25zIGZvciB0aGUgZWxlbWVudC5cbiAgLCBmaW5pc2g6IGZ1bmN0aW9uKCl7XG5cbiAgICAgIHRoaXMuc3RvcCh0cnVlLCBmYWxzZSlcblxuICAgICAgd2hpbGUodGhpcy5kZXF1ZXVlKCkuc2l0dWF0aW9uICYmIHRoaXMuc3RvcCh0cnVlLCBmYWxzZSkpO1xuXG4gICAgICB0aGlzLmNsZWFyUXVldWUoKS5jbGVhckN1cnJlbnQoKVxuXG4gICAgICByZXR1cm4gdGhpc1xuICAgIH1cblxuICAgIC8vIHNldCB0aGUgaW50ZXJuYWwgYW5pbWF0aW9uIHBvaW50ZXIgYXQgdGhlIHN0YXJ0IHBvc2l0aW9uLCBiZWZvcmUgYW55IGxvb3BzLCBhbmQgdXBkYXRlcyB0aGUgdmlzdWFsaXNhdGlvblxuICAsIGF0U3RhcnQ6IGZ1bmN0aW9uKCkge1xuICAgIHJldHVybiB0aGlzLmF0KDAsIHRydWUpXG4gIH1cblxuICAgIC8vIHNldCB0aGUgaW50ZXJuYWwgYW5pbWF0aW9uIHBvaW50ZXIgYXQgdGhlIGVuZCBwb3NpdGlvbiwgYWZ0ZXIgYWxsIHRoZSBsb29wcywgYW5kIHVwZGF0ZXMgdGhlIHZpc3VhbGlzYXRpb25cbiAgLCBhdEVuZDogZnVuY3Rpb24oKSB7XG4gICAgaWYgKHRoaXMuc2l0dWF0aW9uLmxvb3BzID09PSB0cnVlKSB7XG4gICAgICAvLyBJZiBpbiBhIGluZmluaXRlIGxvb3AsIHdlIGVuZCB0aGUgY3VycmVudCBpdGVyYXRpb25cbiAgICAgIHJldHVybiB0aGlzLmF0KHRoaXMuc2l0dWF0aW9uLmxvb3ArMSwgdHJ1ZSlcbiAgICB9IGVsc2UgaWYodHlwZW9mIHRoaXMuc2l0dWF0aW9uLmxvb3BzID09ICdudW1iZXInKSB7XG4gICAgICAvLyBJZiBwZXJmb3JtaW5nIGEgZmluaXRlIG51bWJlciBvZiBsb29wcywgd2UgZ28gYWZ0ZXIgYWxsIHRoZSBsb29wc1xuICAgICAgcmV0dXJuIHRoaXMuYXQodGhpcy5zaXR1YXRpb24ubG9vcHMsIHRydWUpXG4gICAgfSBlbHNlIHtcbiAgICAgIC8vIElmIG5vIGxvb3BzLCB3ZSBqdXN0IGdvIGF0IHRoZSBlbmRcbiAgICAgIHJldHVybiB0aGlzLmF0KDEsIHRydWUpXG4gICAgfVxuICB9XG5cbiAgICAvLyBzZXQgdGhlIGludGVybmFsIGFuaW1hdGlvbiBwb2ludGVyIHRvIHRoZSBzcGVjaWZpZWQgcG9zaXRpb24gYW5kIHVwZGF0ZXMgdGhlIHZpc3VhbGlzYXRpb25cbiAgICAvLyBpZiBpc0Fic1BvcyBpcyB0cnVlLCBwb3MgaXMgdHJlYXRlZCBhcyBhbiBhYnNvbHV0ZSBwb3NpdGlvblxuICAsIGF0OiBmdW5jdGlvbihwb3MsIGlzQWJzUG9zKXtcbiAgICAgIHZhciBkdXJEaXZTcGQgPSB0aGlzLnNpdHVhdGlvbi5kdXJhdGlvbi90aGlzLl9zcGVlZFxuXG4gICAgICB0aGlzLmFic1BvcyA9IHBvc1xuICAgICAgLy8gSWYgcG9zIGlzIG5vdCBhbiBhYnNvbHV0ZSBwb3NpdGlvbiwgd2UgY29udmVydCBpdCBpbnRvIG9uZVxuICAgICAgaWYgKCFpc0Fic1Bvcykge1xuICAgICAgICBpZiAodGhpcy5zaXR1YXRpb24ucmV2ZXJzZWQpIHRoaXMuYWJzUG9zID0gMSAtIHRoaXMuYWJzUG9zXG4gICAgICAgIHRoaXMuYWJzUG9zICs9IHRoaXMuc2l0dWF0aW9uLmxvb3BcbiAgICAgIH1cblxuICAgICAgdGhpcy5zaXR1YXRpb24uc3RhcnQgPSArbmV3IERhdGUgLSB0aGlzLmFic1BvcyAqIGR1ckRpdlNwZFxuICAgICAgdGhpcy5zaXR1YXRpb24uZmluaXNoID0gdGhpcy5zaXR1YXRpb24uc3RhcnQgKyBkdXJEaXZTcGRcblxuICAgICAgcmV0dXJuIHRoaXMuc3RlcCh0cnVlKVxuICAgIH1cblxuICAgIC8qKlxuICAgICAqIHNldHMgb3IgcmV0dXJucyB0aGUgc3BlZWQgb2YgdGhlIGFuaW1hdGlvbnNcbiAgICAgKiBAcGFyYW0gc3BlZWQgbnVsbCB8fCBOdW1iZXIgVGhlIG5ldyBzcGVlZCBvZiB0aGUgYW5pbWF0aW9uc1xuICAgICAqIEByZXR1cm4gTnVtYmVyIHx8IHRoaXNcbiAgICAgKi9cbiAgLCBzcGVlZDogZnVuY3Rpb24oc3BlZWQpe1xuICAgICAgaWYgKHNwZWVkID09PSAwKSByZXR1cm4gdGhpcy5wYXVzZSgpXG5cbiAgICAgIGlmIChzcGVlZCkge1xuICAgICAgICB0aGlzLl9zcGVlZCA9IHNwZWVkXG4gICAgICAgIC8vIFdlIHVzZSBhbiBhYnNvbHV0ZSBwb3NpdGlvbiBoZXJlIHNvIHRoYXQgc3BlZWQgY2FuIGFmZmVjdCB0aGUgZGVsYXkgYmVmb3JlIHRoZSBhbmltYXRpb25cbiAgICAgICAgcmV0dXJuIHRoaXMuYXQodGhpcy5hYnNQb3MsIHRydWUpXG4gICAgICB9IGVsc2UgcmV0dXJuIHRoaXMuX3NwZWVkXG4gICAgfVxuXG4gICAgLy8gTWFrZSBsb29wYWJsZVxuICAsIGxvb3A6IGZ1bmN0aW9uKHRpbWVzLCByZXZlcnNlKSB7XG4gICAgICB2YXIgYyA9IHRoaXMubGFzdCgpXG5cbiAgICAgIC8vIHN0b3JlIHRvdGFsIGxvb3BzXG4gICAgICBjLmxvb3BzID0gKHRpbWVzICE9IG51bGwpID8gdGltZXMgOiB0cnVlXG4gICAgICBjLmxvb3AgPSAwXG5cbiAgICAgIGlmKHJldmVyc2UpIGMucmV2ZXJzaW5nID0gdHJ1ZVxuICAgICAgcmV0dXJuIHRoaXNcbiAgICB9XG5cbiAgICAvLyBwYXVzZXMgdGhlIGFuaW1hdGlvblxuICAsIHBhdXNlOiBmdW5jdGlvbigpe1xuICAgICAgdGhpcy5wYXVzZWQgPSB0cnVlXG4gICAgICB0aGlzLnN0b3BBbmltRnJhbWUoKVxuXG4gICAgICByZXR1cm4gdGhpc1xuICAgIH1cblxuICAgIC8vIHVucGF1c2UgdGhlIGFuaW1hdGlvblxuICAsIHBsYXk6IGZ1bmN0aW9uKCl7XG4gICAgICBpZighdGhpcy5wYXVzZWQpIHJldHVybiB0aGlzXG4gICAgICB0aGlzLnBhdXNlZCA9IGZhbHNlXG4gICAgICAvLyBXZSB1c2UgYW4gYWJzb2x1dGUgcG9zaXRpb24gaGVyZSBzbyB0aGF0IHRoZSBkZWxheSBiZWZvcmUgdGhlIGFuaW1hdGlvbiBjYW4gYmUgcGF1c2VkXG4gICAgICByZXR1cm4gdGhpcy5hdCh0aGlzLmFic1BvcywgdHJ1ZSlcbiAgICB9XG5cbiAgICAvKipcbiAgICAgKiB0b2dnbGUgb3Igc2V0IHRoZSBkaXJlY3Rpb24gb2YgdGhlIGFuaW1hdGlvblxuICAgICAqIHRydWUgc2V0cyBkaXJlY3Rpb24gdG8gYmFja3dhcmRzIHdoaWxlIGZhbHNlIHNldHMgaXQgdG8gZm9yd2FyZHNcbiAgICAgKiBAcGFyYW0gcmV2ZXJzZWQgQm9vbGVhbiBpbmRpY2F0aW5nIHdoZXRoZXIgdG8gcmV2ZXJzZSB0aGUgYW5pbWF0aW9uIG9yIG5vdCAoZGVmYXVsdDogdG9nZ2xlIHRoZSByZXZlcnNlIHN0YXR1cylcbiAgICAgKiBAcmV0dXJuIHRoaXNcbiAgICAgKi9cbiAgLCByZXZlcnNlOiBmdW5jdGlvbihyZXZlcnNlZCl7XG4gICAgICB2YXIgYyA9IHRoaXMubGFzdCgpXG5cbiAgICAgIGlmKHR5cGVvZiByZXZlcnNlZCA9PSAndW5kZWZpbmVkJykgYy5yZXZlcnNlZCA9ICFjLnJldmVyc2VkXG4gICAgICBlbHNlIGMucmV2ZXJzZWQgPSByZXZlcnNlZFxuXG4gICAgICByZXR1cm4gdGhpc1xuICAgIH1cblxuXG4gICAgLyoqXG4gICAgICogcmV0dXJucyBhIGZsb2F0IGZyb20gMC0xIGluZGljYXRpbmcgdGhlIHByb2dyZXNzIG9mIHRoZSBjdXJyZW50IGFuaW1hdGlvblxuICAgICAqIEBwYXJhbSBlYXNlZCBCb29sZWFuIGluZGljYXRpbmcgd2hldGhlciB0aGUgcmV0dXJuZWQgcG9zaXRpb24gc2hvdWxkIGJlIGVhc2VkIG9yIG5vdFxuICAgICAqIEByZXR1cm4gbnVtYmVyXG4gICAgICovXG4gICwgcHJvZ3Jlc3M6IGZ1bmN0aW9uKGVhc2VJdCl7XG4gICAgICByZXR1cm4gZWFzZUl0ID8gdGhpcy5zaXR1YXRpb24uZWFzZSh0aGlzLnBvcykgOiB0aGlzLnBvc1xuICAgIH1cblxuICAgIC8qKlxuICAgICAqIGFkZHMgYSBjYWxsYmFjayBmdW5jdGlvbiB3aGljaCBpcyBjYWxsZWQgd2hlbiB0aGUgY3VycmVudCBhbmltYXRpb24gaXMgZmluaXNoZWRcbiAgICAgKiBAcGFyYW0gZm4gRnVuY3Rpb24gd2hpY2ggc2hvdWxkIGJlIGV4ZWN1dGVkIGFzIGNhbGxiYWNrXG4gICAgICogQHJldHVybiBudW1iZXJcbiAgICAgKi9cbiAgLCBhZnRlcjogZnVuY3Rpb24oZm4pe1xuICAgICAgdmFyIGMgPSB0aGlzLmxhc3QoKVxuICAgICAgICAsIHdyYXBwZXIgPSBmdW5jdGlvbiB3cmFwcGVyKGUpe1xuICAgICAgICAgICAgaWYoZS5kZXRhaWwuc2l0dWF0aW9uID09IGMpe1xuICAgICAgICAgICAgICBmbi5jYWxsKHRoaXMsIGMpXG4gICAgICAgICAgICAgIHRoaXMub2ZmKCdmaW5pc2hlZC5meCcsIHdyYXBwZXIpIC8vIHByZXZlbnQgbWVtb3J5IGxlYWtcbiAgICAgICAgICAgIH1cbiAgICAgICAgICB9XG5cbiAgICAgIHRoaXMudGFyZ2V0KCkub24oJ2ZpbmlzaGVkLmZ4Jywgd3JhcHBlcilcbiAgICAgIHJldHVybiB0aGlzXG4gICAgfVxuXG4gICAgLy8gYWRkcyBhIGNhbGxiYWNrIHdoaWNoIGlzIGNhbGxlZCB3aGVuZXZlciBvbmUgYW5pbWF0aW9uIHN0ZXAgaXMgcGVyZm9ybWVkXG4gICwgZHVyaW5nOiBmdW5jdGlvbihmbil7XG4gICAgICB2YXIgYyA9IHRoaXMubGFzdCgpXG4gICAgICAgICwgd3JhcHBlciA9IGZ1bmN0aW9uKGUpe1xuICAgICAgICAgICAgaWYoZS5kZXRhaWwuc2l0dWF0aW9uID09IGMpe1xuICAgICAgICAgICAgICBmbi5jYWxsKHRoaXMsIGUuZGV0YWlsLnBvcywgU1ZHLm1vcnBoKGUuZGV0YWlsLnBvcyksIGUuZGV0YWlsLmVhc2VkLCBjKVxuICAgICAgICAgICAgfVxuICAgICAgICAgIH1cblxuICAgICAgLy8gc2VlIGFib3ZlXG4gICAgICB0aGlzLnRhcmdldCgpLm9mZignZHVyaW5nLmZ4Jywgd3JhcHBlcikub24oJ2R1cmluZy5meCcsIHdyYXBwZXIpXG5cbiAgICAgIHJldHVybiB0aGlzLmFmdGVyKGZ1bmN0aW9uKCl7XG4gICAgICAgIHRoaXMub2ZmKCdkdXJpbmcuZngnLCB3cmFwcGVyKVxuICAgICAgfSlcbiAgICB9XG5cbiAgICAvLyBjYWxscyBhZnRlciBBTEwgYW5pbWF0aW9ucyBpbiB0aGUgcXVldWUgYXJlIGZpbmlzaGVkXG4gICwgYWZ0ZXJBbGw6IGZ1bmN0aW9uKGZuKXtcbiAgICAgIHZhciB3cmFwcGVyID0gZnVuY3Rpb24gd3JhcHBlcihlKXtcbiAgICAgICAgICAgIGZuLmNhbGwodGhpcylcbiAgICAgICAgICAgIHRoaXMub2ZmKCdhbGxmaW5pc2hlZC5meCcsIHdyYXBwZXIpXG4gICAgICAgICAgfVxuXG4gICAgICAvLyBzZWUgYWJvdmVcbiAgICAgIHRoaXMudGFyZ2V0KCkub2ZmKCdhbGxmaW5pc2hlZC5meCcsIHdyYXBwZXIpLm9uKCdhbGxmaW5pc2hlZC5meCcsIHdyYXBwZXIpXG4gICAgICByZXR1cm4gdGhpc1xuICAgIH1cblxuICAgIC8vIGNhbGxzIG9uIGV2ZXJ5IGFuaW1hdGlvbiBzdGVwIGZvciBhbGwgYW5pbWF0aW9uc1xuICAsIGR1cmluZ0FsbDogZnVuY3Rpb24oZm4pe1xuICAgICAgdmFyIHdyYXBwZXIgPSBmdW5jdGlvbihlKXtcbiAgICAgICAgICAgIGZuLmNhbGwodGhpcywgZS5kZXRhaWwucG9zLCBTVkcubW9ycGgoZS5kZXRhaWwucG9zKSwgZS5kZXRhaWwuZWFzZWQsIGUuZGV0YWlsLnNpdHVhdGlvbilcbiAgICAgICAgICB9XG5cbiAgICAgIHRoaXMudGFyZ2V0KCkub2ZmKCdkdXJpbmcuZngnLCB3cmFwcGVyKS5vbignZHVyaW5nLmZ4Jywgd3JhcHBlcilcblxuICAgICAgcmV0dXJuIHRoaXMuYWZ0ZXJBbGwoZnVuY3Rpb24oKXtcbiAgICAgICAgdGhpcy5vZmYoJ2R1cmluZy5meCcsIHdyYXBwZXIpXG4gICAgICB9KVxuICAgIH1cblxuICAsIGxhc3Q6IGZ1bmN0aW9uKCl7XG4gICAgICByZXR1cm4gdGhpcy5zaXR1YXRpb25zLmxlbmd0aCA/IHRoaXMuc2l0dWF0aW9uc1t0aGlzLnNpdHVhdGlvbnMubGVuZ3RoLTFdIDogdGhpcy5zaXR1YXRpb25cbiAgICB9XG5cbiAgICAvLyBhZGRzIG9uZSBwcm9wZXJ0eSB0byB0aGUgYW5pbWF0aW9uc1xuICAsIGFkZDogZnVuY3Rpb24obWV0aG9kLCBhcmdzLCB0eXBlKXtcbiAgICAgIHRoaXMubGFzdCgpW3R5cGUgfHwgJ2FuaW1hdGlvbnMnXVttZXRob2RdID0gYXJnc1xuICAgICAgc2V0VGltZW91dChmdW5jdGlvbigpe3RoaXMuc3RhcnQoKX0uYmluZCh0aGlzKSwgMClcbiAgICAgIHJldHVybiB0aGlzXG4gICAgfVxuXG4gICAgLyoqIHBlcmZvcm0gb25lIHN0ZXAgb2YgdGhlIGFuaW1hdGlvblxuICAgICAqICBAcGFyYW0gaWdub3JlVGltZSBCb29sZWFuIGluZGljYXRpbmcgd2hldGhlciB0byBpZ25vcmUgdGltZSBhbmQgdXNlIHBvc2l0aW9uIGRpcmVjdGx5IG9yIHJlY2FsY3VsYXRlIHBvc2l0aW9uIGJhc2VkIG9uIHRpbWVcbiAgICAgKiAgQHJldHVybiB0aGlzXG4gICAgICovXG4gICwgc3RlcDogZnVuY3Rpb24oaWdub3JlVGltZSl7XG5cbiAgICAgIC8vIGNvbnZlcnQgY3VycmVudCB0aW1lIHRvIGFuIGFic29sdXRlIHBvc2l0aW9uXG4gICAgICBpZighaWdub3JlVGltZSkgdGhpcy5hYnNQb3MgPSB0aGlzLnRpbWVUb0Fic1BvcygrbmV3IERhdGUpXG5cbiAgICAgIC8vIFRoaXMgcGFydCBjb252ZXJ0IGFuIGFic29sdXRlIHBvc2l0aW9uIHRvIGEgcG9zaXRpb25cbiAgICAgIGlmKHRoaXMuc2l0dWF0aW9uLmxvb3BzICE9PSBmYWxzZSkge1xuICAgICAgICB2YXIgYWJzUG9zLCBhYnNQb3NJbnQsIGxhc3RMb29wXG5cbiAgICAgICAgLy8gSWYgdGhlIGFic29sdXRlIHBvc2l0aW9uIGlzIGJlbG93IDAsIHdlIGp1c3QgdHJlYXQgaXQgYXMgaWYgaXQgd2FzIDBcbiAgICAgICAgYWJzUG9zID0gTWF0aC5tYXgodGhpcy5hYnNQb3MsIDApXG4gICAgICAgIGFic1Bvc0ludCA9IE1hdGguZmxvb3IoYWJzUG9zKVxuXG4gICAgICAgIGlmKHRoaXMuc2l0dWF0aW9uLmxvb3BzID09PSB0cnVlIHx8IGFic1Bvc0ludCA8IHRoaXMuc2l0dWF0aW9uLmxvb3BzKSB7XG4gICAgICAgICAgdGhpcy5wb3MgPSBhYnNQb3MgLSBhYnNQb3NJbnRcbiAgICAgICAgICBsYXN0TG9vcCA9IHRoaXMuc2l0dWF0aW9uLmxvb3BcbiAgICAgICAgICB0aGlzLnNpdHVhdGlvbi5sb29wID0gYWJzUG9zSW50XG4gICAgICAgIH0gZWxzZSB7XG4gICAgICAgICAgdGhpcy5hYnNQb3MgPSB0aGlzLnNpdHVhdGlvbi5sb29wc1xuICAgICAgICAgIHRoaXMucG9zID0gMVxuICAgICAgICAgIC8vIFRoZSAtMSBoZXJlIGlzIGJlY2F1c2Ugd2UgZG9uJ3Qgd2FudCB0byB0b2dnbGUgcmV2ZXJzZWQgd2hlbiBhbGwgdGhlIGxvb3BzIGhhdmUgYmVlbiBjb21wbGV0ZWRcbiAgICAgICAgICBsYXN0TG9vcCA9IHRoaXMuc2l0dWF0aW9uLmxvb3AgLSAxXG4gICAgICAgICAgdGhpcy5zaXR1YXRpb24ubG9vcCA9IHRoaXMuc2l0dWF0aW9uLmxvb3BzXG4gICAgICAgIH1cblxuICAgICAgICBpZih0aGlzLnNpdHVhdGlvbi5yZXZlcnNpbmcpIHtcbiAgICAgICAgICAvLyBUb2dnbGUgcmV2ZXJzZWQgaWYgYW4gb2RkIG51bWJlciBvZiBsb29wcyBhcyBvY2N1cmVkIHNpbmNlIHRoZSBsYXN0IGNhbGwgb2Ygc3RlcFxuICAgICAgICAgIHRoaXMuc2l0dWF0aW9uLnJldmVyc2VkID0gdGhpcy5zaXR1YXRpb24ucmV2ZXJzZWQgIT0gQm9vbGVhbigodGhpcy5zaXR1YXRpb24ubG9vcCAtIGxhc3RMb29wKSAlIDIpXG4gICAgICAgIH1cblxuICAgICAgfSBlbHNlIHtcbiAgICAgICAgLy8gSWYgdGhlcmUgYXJlIG5vIGxvb3AsIHRoZSBhYnNvbHV0ZSBwb3NpdGlvbiBtdXN0IG5vdCBiZSBhYm92ZSAxXG4gICAgICAgIHRoaXMuYWJzUG9zID0gTWF0aC5taW4odGhpcy5hYnNQb3MsIDEpXG4gICAgICAgIHRoaXMucG9zID0gdGhpcy5hYnNQb3NcbiAgICAgIH1cblxuICAgICAgLy8gd2hpbGUgdGhlIGFic29sdXRlIHBvc2l0aW9uIGNhbiBiZSBiZWxvdyAwLCB0aGUgcG9zaXRpb24gbXVzdCBub3QgYmUgYmVsb3cgMFxuICAgICAgaWYodGhpcy5wb3MgPCAwKSB0aGlzLnBvcyA9IDBcblxuICAgICAgaWYodGhpcy5zaXR1YXRpb24ucmV2ZXJzZWQpIHRoaXMucG9zID0gMSAtIHRoaXMucG9zXG5cblxuICAgICAgLy8gYXBwbHkgZWFzaW5nXG4gICAgICB2YXIgZWFzZWQgPSB0aGlzLnNpdHVhdGlvbi5lYXNlKHRoaXMucG9zKVxuXG4gICAgICAvLyBjYWxsIG9uY2UtY2FsbGJhY2tzXG4gICAgICBmb3IodmFyIGkgaW4gdGhpcy5zaXR1YXRpb24ub25jZSl7XG4gICAgICAgIGlmKGkgPiB0aGlzLmxhc3RQb3MgJiYgaSA8PSBlYXNlZCl7XG4gICAgICAgICAgdGhpcy5zaXR1YXRpb24ub25jZVtpXS5jYWxsKHRoaXMudGFyZ2V0KCksIHRoaXMucG9zLCBlYXNlZClcbiAgICAgICAgICBkZWxldGUgdGhpcy5zaXR1YXRpb24ub25jZVtpXVxuICAgICAgICB9XG4gICAgICB9XG5cbiAgICAgIC8vIGZpcmUgZHVyaW5nIGNhbGxiYWNrIHdpdGggcG9zaXRpb24sIGVhc2VkIHBvc2l0aW9uIGFuZCBjdXJyZW50IHNpdHVhdGlvbiBhcyBwYXJhbWV0ZXJcbiAgICAgIGlmKHRoaXMuYWN0aXZlKSB0aGlzLnRhcmdldCgpLmZpcmUoJ2R1cmluZycsIHtwb3M6IHRoaXMucG9zLCBlYXNlZDogZWFzZWQsIGZ4OiB0aGlzLCBzaXR1YXRpb246IHRoaXMuc2l0dWF0aW9ufSlcblxuICAgICAgLy8gdGhlIHVzZXIgbWF5IGNhbGwgc3RvcCBvciBmaW5pc2ggaW4gdGhlIGR1cmluZyBjYWxsYmFja1xuICAgICAgLy8gc28gbWFrZSBzdXJlIHRoYXQgd2Ugc3RpbGwgaGF2ZSBhIHZhbGlkIHNpdHVhdGlvblxuICAgICAgaWYoIXRoaXMuc2l0dWF0aW9uKXtcbiAgICAgICAgcmV0dXJuIHRoaXNcbiAgICAgIH1cblxuICAgICAgLy8gYXBwbHkgdGhlIGFjdHVhbCBhbmltYXRpb24gdG8gZXZlcnkgcHJvcGVydHlcbiAgICAgIHRoaXMuZWFjaEF0KClcblxuICAgICAgLy8gZG8gZmluYWwgY29kZSB3aGVuIHNpdHVhdGlvbiBpcyBmaW5pc2hlZFxuICAgICAgaWYoKHRoaXMucG9zID09IDEgJiYgIXRoaXMuc2l0dWF0aW9uLnJldmVyc2VkKSB8fCAodGhpcy5zaXR1YXRpb24ucmV2ZXJzZWQgJiYgdGhpcy5wb3MgPT0gMCkpe1xuXG4gICAgICAgIC8vIHN0b3AgYW5pbWF0aW9uIGNhbGxiYWNrXG4gICAgICAgIHRoaXMuc3RvcEFuaW1GcmFtZSgpXG5cbiAgICAgICAgLy8gZmlyZSBmaW5pc2hlZCBjYWxsYmFjayB3aXRoIGN1cnJlbnQgc2l0dWF0aW9uIGFzIHBhcmFtZXRlclxuICAgICAgICB0aGlzLnRhcmdldCgpLmZpcmUoJ2ZpbmlzaGVkJywge2Z4OnRoaXMsIHNpdHVhdGlvbjogdGhpcy5zaXR1YXRpb259KVxuXG4gICAgICAgIGlmKCF0aGlzLnNpdHVhdGlvbnMubGVuZ3RoKXtcbiAgICAgICAgICB0aGlzLnRhcmdldCgpLmZpcmUoJ2FsbGZpbmlzaGVkJylcbiAgICAgICAgICB0aGlzLnRhcmdldCgpLm9mZignLmZ4JykgLy8gdGhlcmUgc2hvdWxkbnQgYmUgYW55IGJpbmRpbmcgbGVmdCwgYnV0IHRvIG1ha2Ugc3VyZS4uLlxuICAgICAgICAgIHRoaXMuYWN0aXZlID0gZmFsc2VcbiAgICAgICAgfVxuXG4gICAgICAgIC8vIHN0YXJ0IG5leHQgYW5pbWF0aW9uXG4gICAgICAgIGlmKHRoaXMuYWN0aXZlKSB0aGlzLmRlcXVldWUoKVxuICAgICAgICBlbHNlIHRoaXMuY2xlYXJDdXJyZW50KClcblxuICAgICAgfWVsc2UgaWYoIXRoaXMucGF1c2VkICYmIHRoaXMuYWN0aXZlKXtcbiAgICAgICAgLy8gd2UgY29udGludWUgYW5pbWF0aW5nIHdoZW4gd2UgYXJlIG5vdCBhdCB0aGUgZW5kXG4gICAgICAgIHRoaXMuc3RhcnRBbmltRnJhbWUoKVxuICAgICAgfVxuXG4gICAgICAvLyBzYXZlIGxhc3QgZWFzZWQgcG9zaXRpb24gZm9yIG9uY2UgY2FsbGJhY2sgdHJpZ2dlcmluZ1xuICAgICAgdGhpcy5sYXN0UG9zID0gZWFzZWRcbiAgICAgIHJldHVybiB0aGlzXG5cbiAgICB9XG5cbiAgICAvLyBjYWxjdWxhdGVzIHRoZSBzdGVwIGZvciBldmVyeSBwcm9wZXJ0eSBhbmQgY2FsbHMgYmxvY2sgd2l0aCBpdFxuICAsIGVhY2hBdDogZnVuY3Rpb24oKXtcbiAgICAgIHZhciBpLCBhdCwgc2VsZiA9IHRoaXMsIHRhcmdldCA9IHRoaXMudGFyZ2V0KCksIHMgPSB0aGlzLnNpdHVhdGlvblxuXG4gICAgICAvLyBhcHBseSBhbmltYXRpb25zIHdoaWNoIGNhbiBiZSBjYWxsZWQgdHJvdWdoIGEgbWV0aG9kXG4gICAgICBmb3IoaSBpbiBzLmFuaW1hdGlvbnMpe1xuXG4gICAgICAgIGF0ID0gW10uY29uY2F0KHMuYW5pbWF0aW9uc1tpXSkubWFwKGZ1bmN0aW9uKGVsKXtcbiAgICAgICAgICByZXR1cm4gdHlwZW9mIGVsICE9PSAnc3RyaW5nJyAmJiBlbC5hdCA/IGVsLmF0KHMuZWFzZShzZWxmLnBvcyksIHNlbGYucG9zKSA6IGVsXG4gICAgICAgIH0pXG5cbiAgICAgICAgdGFyZ2V0W2ldLmFwcGx5KHRhcmdldCwgYXQpXG5cbiAgICAgIH1cblxuICAgICAgLy8gYXBwbHkgYW5pbWF0aW9uIHdoaWNoIGhhcyB0byBiZSBhcHBsaWVkIHdpdGggYXR0cigpXG4gICAgICBmb3IoaSBpbiBzLmF0dHJzKXtcblxuICAgICAgICBhdCA9IFtpXS5jb25jYXQocy5hdHRyc1tpXSkubWFwKGZ1bmN0aW9uKGVsKXtcbiAgICAgICAgICByZXR1cm4gdHlwZW9mIGVsICE9PSAnc3RyaW5nJyAmJiBlbC5hdCA/IGVsLmF0KHMuZWFzZShzZWxmLnBvcyksIHNlbGYucG9zKSA6IGVsXG4gICAgICAgIH0pXG5cbiAgICAgICAgdGFyZ2V0LmF0dHIuYXBwbHkodGFyZ2V0LCBhdClcblxuICAgICAgfVxuXG4gICAgICAvLyBhcHBseSBhbmltYXRpb24gd2hpY2ggaGFzIHRvIGJlIGFwcGxpZWQgd2l0aCBzdHlsZSgpXG4gICAgICBmb3IoaSBpbiBzLnN0eWxlcyl7XG5cbiAgICAgICAgYXQgPSBbaV0uY29uY2F0KHMuc3R5bGVzW2ldKS5tYXAoZnVuY3Rpb24oZWwpe1xuICAgICAgICAgIHJldHVybiB0eXBlb2YgZWwgIT09ICdzdHJpbmcnICYmIGVsLmF0ID8gZWwuYXQocy5lYXNlKHNlbGYucG9zKSwgc2VsZi5wb3MpIDogZWxcbiAgICAgICAgfSlcblxuICAgICAgICB0YXJnZXQuc3R5bGUuYXBwbHkodGFyZ2V0LCBhdClcblxuICAgICAgfVxuXG4gICAgICAvLyBhbmltYXRlIGluaXRpYWxUcmFuc2Zvcm1hdGlvbiB3aGljaCBoYXMgdG8gYmUgY2hhaW5lZFxuICAgICAgaWYocy50cmFuc2Zvcm1zLmxlbmd0aCl7XG5cbiAgICAgICAgLy8gZ2V0IGluaXRpYWwgaW5pdGlhbFRyYW5zZm9ybWF0aW9uXG4gICAgICAgIGF0ID0gcy5pbml0aWFsVHJhbnNmb3JtYXRpb25cbiAgICAgICAgZm9yKGkgPSAwLCBsZW4gPSBzLnRyYW5zZm9ybXMubGVuZ3RoOyBpIDwgbGVuOyBpKyspe1xuXG4gICAgICAgICAgLy8gZ2V0IG5leHQgdHJhbnNmb3JtYXRpb24gaW4gY2hhaW5cbiAgICAgICAgICB2YXIgYSA9IHMudHJhbnNmb3Jtc1tpXVxuXG4gICAgICAgICAgLy8gbXVsdGlwbHkgbWF0cml4IGRpcmVjdGx5XG4gICAgICAgICAgaWYoYSBpbnN0YW5jZW9mIFNWRy5NYXRyaXgpe1xuXG4gICAgICAgICAgICBpZihhLnJlbGF0aXZlKXtcbiAgICAgICAgICAgICAgYXQgPSBhdC5tdWx0aXBseShuZXcgU1ZHLk1hdHJpeCgpLm1vcnBoKGEpLmF0KHMuZWFzZSh0aGlzLnBvcykpKVxuICAgICAgICAgICAgfWVsc2V7XG4gICAgICAgICAgICAgIGF0ID0gYXQubW9ycGgoYSkuYXQocy5lYXNlKHRoaXMucG9zKSlcbiAgICAgICAgICAgIH1cbiAgICAgICAgICAgIGNvbnRpbnVlXG4gICAgICAgICAgfVxuXG4gICAgICAgICAgLy8gd2hlbiB0cmFuc2Zvcm1hdGlvbiBpcyBhYnNvbHV0ZSB3ZSBoYXZlIHRvIHJlc2V0IHRoZSBuZWVkZWQgdHJhbnNmb3JtYXRpb24gZmlyc3RcbiAgICAgICAgICBpZighYS5yZWxhdGl2ZSlcbiAgICAgICAgICAgIGEudW5kbyhhdC5leHRyYWN0KCkpXG5cbiAgICAgICAgICAvLyBhbmQgcmVhcHBseSBpdCBhZnRlclxuICAgICAgICAgIGF0ID0gYXQubXVsdGlwbHkoYS5hdChzLmVhc2UodGhpcy5wb3MpKSlcblxuICAgICAgICB9XG5cbiAgICAgICAgLy8gc2V0IG5ldyBtYXRyaXggb24gZWxlbWVudFxuICAgICAgICB0YXJnZXQubWF0cml4KGF0KVxuICAgICAgfVxuXG4gICAgICByZXR1cm4gdGhpc1xuXG4gICAgfVxuXG5cbiAgICAvLyBhZGRzIGFuIG9uY2UtY2FsbGJhY2sgd2hpY2ggaXMgY2FsbGVkIGF0IGEgc3BlY2lmaWMgcG9zaXRpb24gYW5kIG5ldmVyIGFnYWluXG4gICwgb25jZTogZnVuY3Rpb24ocG9zLCBmbiwgaXNFYXNlZCl7XG5cbiAgICAgIGlmKCFpc0Vhc2VkKXBvcyA9IHRoaXMuc2l0dWF0aW9uLmVhc2UocG9zKVxuXG4gICAgICB0aGlzLnNpdHVhdGlvbi5vbmNlW3Bvc10gPSBmblxuXG4gICAgICByZXR1cm4gdGhpc1xuICAgIH1cblxuICB9XG5cbiwgcGFyZW50OiBTVkcuRWxlbWVudFxuXG4gIC8vIEFkZCBtZXRob2QgdG8gcGFyZW50IGVsZW1lbnRzXG4sIGNvbnN0cnVjdDoge1xuICAgIC8vIEdldCBmeCBtb2R1bGUgb3IgY3JlYXRlIGEgbmV3IG9uZSwgdGhlbiBhbmltYXRlIHdpdGggZ2l2ZW4gZHVyYXRpb24gYW5kIGVhc2VcbiAgICBhbmltYXRlOiBmdW5jdGlvbihvLCBlYXNlLCBkZWxheSkge1xuICAgICAgcmV0dXJuICh0aGlzLmZ4IHx8ICh0aGlzLmZ4ID0gbmV3IFNWRy5GWCh0aGlzKSkpLmFuaW1hdGUobywgZWFzZSwgZGVsYXkpXG4gICAgfVxuICAsIGRlbGF5OiBmdW5jdGlvbihkZWxheSl7XG4gICAgICByZXR1cm4gKHRoaXMuZnggfHwgKHRoaXMuZnggPSBuZXcgU1ZHLkZYKHRoaXMpKSkuZGVsYXkoZGVsYXkpXG4gICAgfVxuICAsIHN0b3A6IGZ1bmN0aW9uKGp1bXBUb0VuZCwgY2xlYXJRdWV1ZSkge1xuICAgICAgaWYgKHRoaXMuZngpXG4gICAgICAgIHRoaXMuZnguc3RvcChqdW1wVG9FbmQsIGNsZWFyUXVldWUpXG5cbiAgICAgIHJldHVybiB0aGlzXG4gICAgfVxuICAsIGZpbmlzaDogZnVuY3Rpb24oKSB7XG4gICAgICBpZiAodGhpcy5meClcbiAgICAgICAgdGhpcy5meC5maW5pc2goKVxuXG4gICAgICByZXR1cm4gdGhpc1xuICAgIH1cbiAgICAvLyBQYXVzZSBjdXJyZW50IGFuaW1hdGlvblxuICAsIHBhdXNlOiBmdW5jdGlvbigpIHtcbiAgICAgIGlmICh0aGlzLmZ4KVxuICAgICAgICB0aGlzLmZ4LnBhdXNlKClcblxuICAgICAgcmV0dXJuIHRoaXNcbiAgICB9XG4gICAgLy8gUGxheSBwYXVzZWQgY3VycmVudCBhbmltYXRpb25cbiAgLCBwbGF5OiBmdW5jdGlvbigpIHtcbiAgICAgIGlmICh0aGlzLmZ4KVxuICAgICAgICB0aGlzLmZ4LnBsYXkoKVxuXG4gICAgICByZXR1cm4gdGhpc1xuICAgIH1cbiAgICAvLyBTZXQvR2V0IHRoZSBzcGVlZCBvZiB0aGUgYW5pbWF0aW9uc1xuICAsIHNwZWVkOiBmdW5jdGlvbihzcGVlZCkge1xuICAgICAgaWYgKHRoaXMuZngpXG4gICAgICAgIGlmIChzcGVlZCA9PSBudWxsKVxuICAgICAgICAgIHJldHVybiB0aGlzLmZ4LnNwZWVkKClcbiAgICAgICAgZWxzZVxuICAgICAgICAgIHRoaXMuZnguc3BlZWQoc3BlZWQpXG5cbiAgICAgIHJldHVybiB0aGlzXG4gICAgfVxuICB9XG5cbn0pXG5cbi8vIE1vcnBoT2JqIGlzIHVzZWQgd2hlbmV2ZXIgbm8gbW9ycGhhYmxlIG9iamVjdCBpcyBnaXZlblxuU1ZHLk1vcnBoT2JqID0gU1ZHLmludmVudCh7XG5cbiAgY3JlYXRlOiBmdW5jdGlvbihmcm9tLCB0byl7XG4gICAgLy8gcHJlcGFyZSBjb2xvciBmb3IgbW9ycGhpbmdcbiAgICBpZihTVkcuQ29sb3IuaXNDb2xvcih0bykpIHJldHVybiBuZXcgU1ZHLkNvbG9yKGZyb20pLm1vcnBoKHRvKVxuICAgIC8vIHByZXBhcmUgbnVtYmVyIGZvciBtb3JwaGluZ1xuICAgIGlmKFNWRy5yZWdleC5udW1iZXJBbmRVbml0LnRlc3QodG8pKSByZXR1cm4gbmV3IFNWRy5OdW1iZXIoZnJvbSkubW9ycGgodG8pXG5cbiAgICAvLyBwcmVwYXJlIGZvciBwbGFpbiBtb3JwaGluZ1xuICAgIHRoaXMudmFsdWUgPSAwXG4gICAgdGhpcy5kZXN0aW5hdGlvbiA9IHRvXG4gIH1cblxuLCBleHRlbmQ6IHtcbiAgICBhdDogZnVuY3Rpb24ocG9zLCByZWFsKXtcbiAgICAgIHJldHVybiByZWFsIDwgMSA/IHRoaXMudmFsdWUgOiB0aGlzLmRlc3RpbmF0aW9uXG4gICAgfSxcblxuICAgIHZhbHVlT2Y6IGZ1bmN0aW9uKCl7XG4gICAgICByZXR1cm4gdGhpcy52YWx1ZVxuICAgIH1cbiAgfVxuXG59KVxuXG5TVkcuZXh0ZW5kKFNWRy5GWCwge1xuICAvLyBBZGQgYW5pbWF0YWJsZSBhdHRyaWJ1dGVzXG4gIGF0dHI6IGZ1bmN0aW9uKGEsIHYsIHJlbGF0aXZlKSB7XG4gICAgLy8gYXBwbHkgYXR0cmlidXRlcyBpbmRpdmlkdWFsbHlcbiAgICBpZiAodHlwZW9mIGEgPT0gJ29iamVjdCcpIHtcbiAgICAgIGZvciAodmFyIGtleSBpbiBhKVxuICAgICAgICB0aGlzLmF0dHIoa2V5LCBhW2tleV0pXG5cbiAgICB9IGVsc2Uge1xuICAgICAgLy8gdGhlIE1vcnBoT2JqIHRha2VzIGNhcmUgYWJvdXQgdGhlIHJpZ2h0IGZ1bmN0aW9uIHVzZWRcbiAgICAgIHRoaXMuYWRkKGEsIG5ldyBTVkcuTW9ycGhPYmoobnVsbCwgdiksICdhdHRycycpXG4gICAgfVxuXG4gICAgcmV0dXJuIHRoaXNcbiAgfVxuICAvLyBBZGQgYW5pbWF0YWJsZSBzdHlsZXNcbiwgc3R5bGU6IGZ1bmN0aW9uKHMsIHYpIHtcbiAgICBpZiAodHlwZW9mIHMgPT0gJ29iamVjdCcpXG4gICAgICBmb3IgKHZhciBrZXkgaW4gcylcbiAgICAgICAgdGhpcy5zdHlsZShrZXksIHNba2V5XSlcblxuICAgIGVsc2VcbiAgICAgIHRoaXMuYWRkKHMsIG5ldyBTVkcuTW9ycGhPYmoobnVsbCwgdiksICdzdHlsZXMnKVxuXG4gICAgcmV0dXJuIHRoaXNcbiAgfVxuICAvLyBBbmltYXRhYmxlIHgtYXhpc1xuLCB4OiBmdW5jdGlvbih4LCByZWxhdGl2ZSkge1xuICAgIGlmKHRoaXMudGFyZ2V0KCkgaW5zdGFuY2VvZiBTVkcuRyl7XG4gICAgICB0aGlzLnRyYW5zZm9ybSh7eDp4fSwgcmVsYXRpdmUpXG4gICAgICByZXR1cm4gdGhpc1xuICAgIH1cblxuICAgIHZhciBudW0gPSBuZXcgU1ZHLk51bWJlcigpLm1vcnBoKHgpXG4gICAgbnVtLnJlbGF0aXZlID0gcmVsYXRpdmVcbiAgICByZXR1cm4gdGhpcy5hZGQoJ3gnLCBudW0pXG4gIH1cbiAgLy8gQW5pbWF0YWJsZSB5LWF4aXNcbiwgeTogZnVuY3Rpb24oeSwgcmVsYXRpdmUpIHtcbiAgICBpZih0aGlzLnRhcmdldCgpIGluc3RhbmNlb2YgU1ZHLkcpe1xuICAgICAgdGhpcy50cmFuc2Zvcm0oe3k6eX0sIHJlbGF0aXZlKVxuICAgICAgcmV0dXJuIHRoaXNcbiAgICB9XG5cbiAgICB2YXIgbnVtID0gbmV3IFNWRy5OdW1iZXIoKS5tb3JwaCh5KVxuICAgIG51bS5yZWxhdGl2ZSA9IHJlbGF0aXZlXG4gICAgcmV0dXJuIHRoaXMuYWRkKCd5JywgbnVtKVxuICB9XG4gIC8vIEFuaW1hdGFibGUgY2VudGVyIHgtYXhpc1xuLCBjeDogZnVuY3Rpb24oeCkge1xuICAgIHJldHVybiB0aGlzLmFkZCgnY3gnLCBuZXcgU1ZHLk51bWJlcigpLm1vcnBoKHgpKVxuICB9XG4gIC8vIEFuaW1hdGFibGUgY2VudGVyIHktYXhpc1xuLCBjeTogZnVuY3Rpb24oeSkge1xuICAgIHJldHVybiB0aGlzLmFkZCgnY3knLCBuZXcgU1ZHLk51bWJlcigpLm1vcnBoKHkpKVxuICB9XG4gIC8vIEFkZCBhbmltYXRhYmxlIG1vdmVcbiwgbW92ZTogZnVuY3Rpb24oeCwgeSkge1xuICAgIHJldHVybiB0aGlzLngoeCkueSh5KVxuICB9XG4gIC8vIEFkZCBhbmltYXRhYmxlIGNlbnRlclxuLCBjZW50ZXI6IGZ1bmN0aW9uKHgsIHkpIHtcbiAgICByZXR1cm4gdGhpcy5jeCh4KS5jeSh5KVxuICB9XG4gIC8vIEFkZCBhbmltYXRhYmxlIHNpemVcbiwgc2l6ZTogZnVuY3Rpb24od2lkdGgsIGhlaWdodCkge1xuICAgIGlmICh0aGlzLnRhcmdldCgpIGluc3RhbmNlb2YgU1ZHLlRleHQpIHtcbiAgICAgIC8vIGFuaW1hdGUgZm9udCBzaXplIGZvciBUZXh0IGVsZW1lbnRzXG4gICAgICB0aGlzLmF0dHIoJ2ZvbnQtc2l6ZScsIHdpZHRoKVxuXG4gICAgfSBlbHNlIHtcbiAgICAgIC8vIGFuaW1hdGUgYmJveCBiYXNlZCBzaXplIGZvciBhbGwgb3RoZXIgZWxlbWVudHNcbiAgICAgIHZhciBib3hcblxuICAgICAgaWYoIXdpZHRoIHx8ICFoZWlnaHQpe1xuICAgICAgICBib3ggPSB0aGlzLnRhcmdldCgpLmJib3goKVxuICAgICAgfVxuXG4gICAgICBpZighd2lkdGgpe1xuICAgICAgICB3aWR0aCA9IGJveC53aWR0aCAvIGJveC5oZWlnaHQgICogaGVpZ2h0XG4gICAgICB9XG5cbiAgICAgIGlmKCFoZWlnaHQpe1xuICAgICAgICBoZWlnaHQgPSBib3guaGVpZ2h0IC8gYm94LndpZHRoICAqIHdpZHRoXG4gICAgICB9XG5cbiAgICAgIHRoaXMuYWRkKCd3aWR0aCcgLCBuZXcgU1ZHLk51bWJlcigpLm1vcnBoKHdpZHRoKSlcbiAgICAgICAgICAuYWRkKCdoZWlnaHQnLCBuZXcgU1ZHLk51bWJlcigpLm1vcnBoKGhlaWdodCkpXG5cbiAgICB9XG5cbiAgICByZXR1cm4gdGhpc1xuICB9XG4gIC8vIEFkZCBhbmltYXRhYmxlIHBsb3RcbiwgcGxvdDogZnVuY3Rpb24ocCkge1xuICAgIHJldHVybiB0aGlzLmFkZCgncGxvdCcsIHRoaXMudGFyZ2V0KCkuYXJyYXkoKS5tb3JwaChwKSlcbiAgfVxuICAvLyBBZGQgbGVhZGluZyBtZXRob2RcbiwgbGVhZGluZzogZnVuY3Rpb24odmFsdWUpIHtcbiAgICByZXR1cm4gdGhpcy50YXJnZXQoKS5sZWFkaW5nID9cbiAgICAgIHRoaXMuYWRkKCdsZWFkaW5nJywgbmV3IFNWRy5OdW1iZXIoKS5tb3JwaCh2YWx1ZSkpIDpcbiAgICAgIHRoaXNcbiAgfVxuICAvLyBBZGQgYW5pbWF0YWJsZSB2aWV3Ym94XG4sIHZpZXdib3g6IGZ1bmN0aW9uKHgsIHksIHdpZHRoLCBoZWlnaHQpIHtcbiAgICBpZiAodGhpcy50YXJnZXQoKSBpbnN0YW5jZW9mIFNWRy5Db250YWluZXIpIHtcbiAgICAgIHRoaXMuYWRkKCd2aWV3Ym94JywgbmV3IFNWRy5WaWV3Qm94KHgsIHksIHdpZHRoLCBoZWlnaHQpKVxuICAgIH1cblxuICAgIHJldHVybiB0aGlzXG4gIH1cbiwgdXBkYXRlOiBmdW5jdGlvbihvKSB7XG4gICAgaWYgKHRoaXMudGFyZ2V0KCkgaW5zdGFuY2VvZiBTVkcuU3RvcCkge1xuICAgICAgaWYgKHR5cGVvZiBvID09ICdudW1iZXInIHx8IG8gaW5zdGFuY2VvZiBTVkcuTnVtYmVyKSB7XG4gICAgICAgIHJldHVybiB0aGlzLnVwZGF0ZSh7XG4gICAgICAgICAgb2Zmc2V0OiAgYXJndW1lbnRzWzBdXG4gICAgICAgICwgY29sb3I6ICAgYXJndW1lbnRzWzFdXG4gICAgICAgICwgb3BhY2l0eTogYXJndW1lbnRzWzJdXG4gICAgICAgIH0pXG4gICAgICB9XG5cbiAgICAgIGlmIChvLm9wYWNpdHkgIT0gbnVsbCkgdGhpcy5hdHRyKCdzdG9wLW9wYWNpdHknLCBvLm9wYWNpdHkpXG4gICAgICBpZiAoby5jb2xvciAgICE9IG51bGwpIHRoaXMuYXR0cignc3RvcC1jb2xvcicsIG8uY29sb3IpXG4gICAgICBpZiAoby5vZmZzZXQgICE9IG51bGwpIHRoaXMuYXR0cignb2Zmc2V0Jywgby5vZmZzZXQpXG4gICAgfVxuXG4gICAgcmV0dXJuIHRoaXNcbiAgfVxufSlcblxuU1ZHLkJCb3ggPSBTVkcuaW52ZW50KHtcbiAgLy8gSW5pdGlhbGl6ZVxuICBjcmVhdGU6IGZ1bmN0aW9uKGVsZW1lbnQpIHtcbiAgICAvLyBnZXQgdmFsdWVzIGlmIGVsZW1lbnQgaXMgZ2l2ZW5cbiAgICBpZiAoZWxlbWVudCkge1xuICAgICAgdmFyIGJveFxuXG4gICAgICAvLyB5ZXMgdGhpcyBpcyB1Z2x5LCBidXQgRmlyZWZveCBjYW4gYmUgYSBiaXRjaCB3aGVuIGl0IGNvbWVzIHRvIGVsZW1lbnRzIHRoYXQgYXJlIG5vdCB5ZXQgcmVuZGVyZWRcbiAgICAgIHRyeSB7XG5cbiAgICAgICAgLy8gdGhlIGVsZW1lbnQgaXMgTk9UIGluIHRoZSBkb20sIHRocm93IGVycm9yXG4gICAgICAgIGlmKCFkb2N1bWVudC5kb2N1bWVudEVsZW1lbnQuY29udGFpbnMoZWxlbWVudC5ub2RlKSkgdGhyb3cgbmV3IEV4Y2VwdGlvbignRWxlbWVudCBub3QgaW4gdGhlIGRvbScpXG5cbiAgICAgICAgLy8gZmluZCBuYXRpdmUgYmJveFxuICAgICAgICBib3ggPSBlbGVtZW50Lm5vZGUuZ2V0QkJveCgpXG4gICAgICB9IGNhdGNoKGUpIHtcbiAgICAgICAgaWYoZWxlbWVudCBpbnN0YW5jZW9mIFNWRy5TaGFwZSl7XG4gICAgICAgICAgdmFyIGNsb25lID0gZWxlbWVudC5jbG9uZShTVkcucGFyc2VyLmRyYXcpLnNob3coKVxuICAgICAgICAgIGJveCA9IGNsb25lLmJib3goKVxuICAgICAgICAgIGNsb25lLnJlbW92ZSgpXG4gICAgICAgIH1lbHNle1xuICAgICAgICAgIGJveCA9IHtcbiAgICAgICAgICAgIHg6ICAgICAgZWxlbWVudC5ub2RlLmNsaWVudExlZnRcbiAgICAgICAgICAsIHk6ICAgICAgZWxlbWVudC5ub2RlLmNsaWVudFRvcFxuICAgICAgICAgICwgd2lkdGg6ICBlbGVtZW50Lm5vZGUuY2xpZW50V2lkdGhcbiAgICAgICAgICAsIGhlaWdodDogZWxlbWVudC5ub2RlLmNsaWVudEhlaWdodFxuICAgICAgICAgIH1cbiAgICAgICAgfVxuICAgICAgfVxuXG4gICAgICAvLyBwbGFpbiB4IGFuZCB5XG4gICAgICB0aGlzLnggPSBib3gueFxuICAgICAgdGhpcy55ID0gYm94LnlcblxuICAgICAgLy8gcGxhaW4gd2lkdGggYW5kIGhlaWdodFxuICAgICAgdGhpcy53aWR0aCAgPSBib3gud2lkdGhcbiAgICAgIHRoaXMuaGVpZ2h0ID0gYm94LmhlaWdodFxuICAgIH1cblxuICAgIC8vIGFkZCBjZW50ZXIsIHJpZ2h0IGFuZCBib3R0b21cbiAgICBmdWxsQm94KHRoaXMpXG4gIH1cblxuICAvLyBEZWZpbmUgUGFyZW50XG4sIHBhcmVudDogU1ZHLkVsZW1lbnRcblxuICAvLyBDb25zdHJ1Y3RvclxuLCBjb25zdHJ1Y3Q6IHtcbiAgICAvLyBHZXQgYm91bmRpbmcgYm94XG4gICAgYmJveDogZnVuY3Rpb24oKSB7XG4gICAgICByZXR1cm4gbmV3IFNWRy5CQm94KHRoaXMpXG4gICAgfVxuICB9XG5cbn0pXG5cblNWRy5UQm94ID0gU1ZHLmludmVudCh7XG4gIC8vIEluaXRpYWxpemVcbiAgY3JlYXRlOiBmdW5jdGlvbihlbGVtZW50KSB7XG4gICAgLy8gZ2V0IHZhbHVlcyBpZiBlbGVtZW50IGlzIGdpdmVuXG4gICAgaWYgKGVsZW1lbnQpIHtcbiAgICAgIHZhciB0ICAgPSBlbGVtZW50LmN0bSgpLmV4dHJhY3QoKVxuICAgICAgICAsIGJveCA9IGVsZW1lbnQuYmJveCgpXG5cbiAgICAgIC8vIHdpZHRoIGFuZCBoZWlnaHQgaW5jbHVkaW5nIHRyYW5zZm9ybWF0aW9uc1xuICAgICAgdGhpcy53aWR0aCAgPSBib3gud2lkdGggICogdC5zY2FsZVhcbiAgICAgIHRoaXMuaGVpZ2h0ID0gYm94LmhlaWdodCAqIHQuc2NhbGVZXG5cbiAgICAgIC8vIHggYW5kIHkgaW5jbHVkaW5nIHRyYW5zZm9ybWF0aW9uc1xuICAgICAgdGhpcy54ID0gYm94LnggKyB0LnhcbiAgICAgIHRoaXMueSA9IGJveC55ICsgdC55XG4gICAgfVxuXG4gICAgLy8gYWRkIGNlbnRlciwgcmlnaHQgYW5kIGJvdHRvbVxuICAgIGZ1bGxCb3godGhpcylcbiAgfVxuXG4gIC8vIERlZmluZSBQYXJlbnRcbiwgcGFyZW50OiBTVkcuRWxlbWVudFxuXG4gIC8vIENvbnN0cnVjdG9yXG4sIGNvbnN0cnVjdDoge1xuICAgIC8vIEdldCB0cmFuc2Zvcm1lZCBib3VuZGluZyBib3hcbiAgICB0Ym94OiBmdW5jdGlvbigpIHtcbiAgICAgIHJldHVybiBuZXcgU1ZHLlRCb3godGhpcylcbiAgICB9XG4gIH1cblxufSlcblxuXG5TVkcuUkJveCA9IFNWRy5pbnZlbnQoe1xuICAvLyBJbml0aWFsaXplXG4gIGNyZWF0ZTogZnVuY3Rpb24oZWxlbWVudCkge1xuICAgIGlmIChlbGVtZW50KSB7XG4gICAgICB2YXIgZSAgICA9IGVsZW1lbnQuZG9jKCkucGFyZW50KClcbiAgICAgICAgLCBib3ggID0gZWxlbWVudC5ub2RlLmdldEJvdW5kaW5nQ2xpZW50UmVjdCgpXG4gICAgICAgICwgem9vbSA9IDFcblxuICAgICAgLy8gZ2V0IHNjcmVlbiBvZmZzZXRcbiAgICAgIHRoaXMueCA9IGJveC5sZWZ0XG4gICAgICB0aGlzLnkgPSBib3gudG9wXG5cbiAgICAgIC8vIHN1YnRyYWN0IHBhcmVudCBvZmZzZXRcbiAgICAgIHRoaXMueCAtPSBlLm9mZnNldExlZnRcbiAgICAgIHRoaXMueSAtPSBlLm9mZnNldFRvcFxuXG4gICAgICB3aGlsZSAoZSA9IGUub2Zmc2V0UGFyZW50KSB7XG4gICAgICAgIHRoaXMueCAtPSBlLm9mZnNldExlZnRcbiAgICAgICAgdGhpcy55IC09IGUub2Zmc2V0VG9wXG4gICAgICB9XG5cbiAgICAgIC8vIGNhbGN1bGF0ZSBjdW11bGF0aXZlIHpvb20gZnJvbSBzdmcgZG9jdW1lbnRzXG4gICAgICBlID0gZWxlbWVudFxuICAgICAgd2hpbGUgKGUucGFyZW50ICYmIChlID0gZS5wYXJlbnQoKSkpIHtcbiAgICAgICAgaWYgKGUudmlld2JveCkge1xuICAgICAgICAgIHpvb20gKj0gZS52aWV3Ym94KCkuem9vbVxuICAgICAgICAgIHRoaXMueCAtPSBlLngoKSB8fCAwXG4gICAgICAgICAgdGhpcy55IC09IGUueSgpIHx8IDBcbiAgICAgICAgfVxuICAgICAgfVxuXG4gICAgICAvLyByZWNhbGN1bGF0ZSB2aWV3Ym94IGRpc3RvcnRpb25cbiAgICAgIHRoaXMud2lkdGggID0gYm94LndpZHRoICAvPSB6b29tXG4gICAgICB0aGlzLmhlaWdodCA9IGJveC5oZWlnaHQgLz0gem9vbVxuICAgIH1cblxuICAgIC8vIGFkZCBjZW50ZXIsIHJpZ2h0IGFuZCBib3R0b21cbiAgICBmdWxsQm94KHRoaXMpXG5cbiAgICAvLyBvZmZzZXQgYnkgd2luZG93IHNjcm9sbCBwb3NpdGlvbiwgYmVjYXVzZSBnZXRCb3VuZGluZ0NsaWVudFJlY3QgY2hhbmdlcyB3aGVuIHdpbmRvdyBpcyBzY3JvbGxlZFxuICAgIHRoaXMueCArPSB3aW5kb3cucGFnZVhPZmZzZXRcbiAgICB0aGlzLnkgKz0gd2luZG93LnBhZ2VZT2Zmc2V0XG4gIH1cblxuICAvLyBkZWZpbmUgUGFyZW50XG4sIHBhcmVudDogU1ZHLkVsZW1lbnRcblxuICAvLyBDb25zdHJ1Y3RvclxuLCBjb25zdHJ1Y3Q6IHtcbiAgICAvLyBHZXQgcmVjdCBib3hcbiAgICByYm94OiBmdW5jdGlvbigpIHtcbiAgICAgIHJldHVybiBuZXcgU1ZHLlJCb3godGhpcylcbiAgICB9XG4gIH1cblxufSlcblxuLy8gQWRkIHVuaXZlcnNhbCBtZXJnZSBtZXRob2RcbjtbU1ZHLkJCb3gsIFNWRy5UQm94LCBTVkcuUkJveF0uZm9yRWFjaChmdW5jdGlvbihjKSB7XG5cbiAgU1ZHLmV4dGVuZChjLCB7XG4gICAgLy8gTWVyZ2UgcmVjdCBib3ggd2l0aCBhbm90aGVyLCByZXR1cm4gYSBuZXcgaW5zdGFuY2VcbiAgICBtZXJnZTogZnVuY3Rpb24oYm94KSB7XG4gICAgICB2YXIgYiA9IG5ldyBjKClcblxuICAgICAgLy8gbWVyZ2UgYm94ZXNcbiAgICAgIGIueCAgICAgID0gTWF0aC5taW4odGhpcy54LCBib3gueClcbiAgICAgIGIueSAgICAgID0gTWF0aC5taW4odGhpcy55LCBib3gueSlcbiAgICAgIGIud2lkdGggID0gTWF0aC5tYXgodGhpcy54ICsgdGhpcy53aWR0aCwgIGJveC54ICsgYm94LndpZHRoKSAgLSBiLnhcbiAgICAgIGIuaGVpZ2h0ID0gTWF0aC5tYXgodGhpcy55ICsgdGhpcy5oZWlnaHQsIGJveC55ICsgYm94LmhlaWdodCkgLSBiLnlcblxuICAgICAgcmV0dXJuIGZ1bGxCb3goYilcbiAgICB9XG5cbiAgfSlcblxufSlcblxuU1ZHLk1hdHJpeCA9IFNWRy5pbnZlbnQoe1xuICAvLyBJbml0aWFsaXplXG4gIGNyZWF0ZTogZnVuY3Rpb24oc291cmNlKSB7XG4gICAgdmFyIGksIGJhc2UgPSBhcnJheVRvTWF0cml4KFsxLCAwLCAwLCAxLCAwLCAwXSlcblxuICAgIC8vIGVuc3VyZSBzb3VyY2UgYXMgb2JqZWN0XG4gICAgc291cmNlID0gc291cmNlIGluc3RhbmNlb2YgU1ZHLkVsZW1lbnQgP1xuICAgICAgc291cmNlLm1hdHJpeGlmeSgpIDpcbiAgICB0eXBlb2Ygc291cmNlID09PSAnc3RyaW5nJyA/XG4gICAgICBzdHJpbmdUb01hdHJpeChzb3VyY2UpIDpcbiAgICBhcmd1bWVudHMubGVuZ3RoID09IDYgP1xuICAgICAgYXJyYXlUb01hdHJpeChbXS5zbGljZS5jYWxsKGFyZ3VtZW50cykpIDpcbiAgICB0eXBlb2Ygc291cmNlID09PSAnb2JqZWN0JyA/XG4gICAgICBzb3VyY2UgOiBiYXNlXG5cbiAgICAvLyBtZXJnZSBzb3VyY2VcbiAgICBmb3IgKGkgPSBhYmNkZWYubGVuZ3RoIC0gMTsgaSA+PSAwOyAtLWkpXG4gICAgICB0aGlzW2FiY2RlZltpXV0gPSBzb3VyY2UgJiYgdHlwZW9mIHNvdXJjZVthYmNkZWZbaV1dID09PSAnbnVtYmVyJyA/XG4gICAgICAgIHNvdXJjZVthYmNkZWZbaV1dIDogYmFzZVthYmNkZWZbaV1dXG4gIH1cblxuICAvLyBBZGQgbWV0aG9kc1xuLCBleHRlbmQ6IHtcbiAgICAvLyBFeHRyYWN0IGluZGl2aWR1YWwgdHJhbnNmb3JtYXRpb25zXG4gICAgZXh0cmFjdDogZnVuY3Rpb24oKSB7XG4gICAgICAvLyBmaW5kIGRlbHRhIHRyYW5zZm9ybSBwb2ludHNcbiAgICAgIHZhciBweCAgICA9IGRlbHRhVHJhbnNmb3JtUG9pbnQodGhpcywgMCwgMSlcbiAgICAgICAgLCBweSAgICA9IGRlbHRhVHJhbnNmb3JtUG9pbnQodGhpcywgMSwgMClcbiAgICAgICAgLCBza2V3WCA9IDE4MCAvIE1hdGguUEkgKiBNYXRoLmF0YW4yKHB4LnksIHB4LngpIC0gOTBcblxuICAgICAgcmV0dXJuIHtcbiAgICAgICAgLy8gdHJhbnNsYXRpb25cbiAgICAgICAgeDogICAgICAgIHRoaXMuZVxuICAgICAgLCB5OiAgICAgICAgdGhpcy5mXG4gICAgICAsIHRyYW5zZm9ybWVkWDoodGhpcy5lICogTWF0aC5jb3Moc2tld1ggKiBNYXRoLlBJIC8gMTgwKSArIHRoaXMuZiAqIE1hdGguc2luKHNrZXdYICogTWF0aC5QSSAvIDE4MCkpIC8gTWF0aC5zcXJ0KHRoaXMuYSAqIHRoaXMuYSArIHRoaXMuYiAqIHRoaXMuYilcbiAgICAgICwgdHJhbnNmb3JtZWRZOih0aGlzLmYgKiBNYXRoLmNvcyhza2V3WCAqIE1hdGguUEkgLyAxODApICsgdGhpcy5lICogTWF0aC5zaW4oLXNrZXdYICogTWF0aC5QSSAvIDE4MCkpIC8gTWF0aC5zcXJ0KHRoaXMuYyAqIHRoaXMuYyArIHRoaXMuZCAqIHRoaXMuZClcbiAgICAgICAgLy8gc2tld1xuICAgICAgLCBza2V3WDogICAgLXNrZXdYXG4gICAgICAsIHNrZXdZOiAgICAxODAgLyBNYXRoLlBJICogTWF0aC5hdGFuMihweS55LCBweS54KVxuICAgICAgICAvLyBzY2FsZVxuICAgICAgLCBzY2FsZVg6ICAgTWF0aC5zcXJ0KHRoaXMuYSAqIHRoaXMuYSArIHRoaXMuYiAqIHRoaXMuYilcbiAgICAgICwgc2NhbGVZOiAgIE1hdGguc3FydCh0aGlzLmMgKiB0aGlzLmMgKyB0aGlzLmQgKiB0aGlzLmQpXG4gICAgICAgIC8vIHJvdGF0aW9uXG4gICAgICAsIHJvdGF0aW9uOiBza2V3WFxuICAgICAgLCBhOiB0aGlzLmFcbiAgICAgICwgYjogdGhpcy5iXG4gICAgICAsIGM6IHRoaXMuY1xuICAgICAgLCBkOiB0aGlzLmRcbiAgICAgICwgZTogdGhpcy5lXG4gICAgICAsIGY6IHRoaXMuZlxuICAgICAgLCBtYXRyaXg6IG5ldyBTVkcuTWF0cml4KHRoaXMpXG4gICAgICB9XG4gICAgfVxuICAgIC8vIENsb25lIG1hdHJpeFxuICAsIGNsb25lOiBmdW5jdGlvbigpIHtcbiAgICAgIHJldHVybiBuZXcgU1ZHLk1hdHJpeCh0aGlzKVxuICAgIH1cbiAgICAvLyBNb3JwaCBvbmUgbWF0cml4IGludG8gYW5vdGhlclxuICAsIG1vcnBoOiBmdW5jdGlvbihtYXRyaXgpIHtcbiAgICAgIC8vIHN0b3JlIG5ldyBkZXN0aW5hdGlvblxuICAgICAgdGhpcy5kZXN0aW5hdGlvbiA9IG5ldyBTVkcuTWF0cml4KG1hdHJpeClcblxuICAgICAgcmV0dXJuIHRoaXNcbiAgICB9XG4gICAgLy8gR2V0IG1vcnBoZWQgbWF0cml4IGF0IGEgZ2l2ZW4gcG9zaXRpb25cbiAgLCBhdDogZnVuY3Rpb24ocG9zKSB7XG4gICAgICAvLyBtYWtlIHN1cmUgYSBkZXN0aW5hdGlvbiBpcyBkZWZpbmVkXG4gICAgICBpZiAoIXRoaXMuZGVzdGluYXRpb24pIHJldHVybiB0aGlzXG5cbiAgICAgIC8vIGNhbGN1bGF0ZSBtb3JwaGVkIG1hdHJpeCBhdCBhIGdpdmVuIHBvc2l0aW9uXG4gICAgICB2YXIgbWF0cml4ID0gbmV3IFNWRy5NYXRyaXgoe1xuICAgICAgICBhOiB0aGlzLmEgKyAodGhpcy5kZXN0aW5hdGlvbi5hIC0gdGhpcy5hKSAqIHBvc1xuICAgICAgLCBiOiB0aGlzLmIgKyAodGhpcy5kZXN0aW5hdGlvbi5iIC0gdGhpcy5iKSAqIHBvc1xuICAgICAgLCBjOiB0aGlzLmMgKyAodGhpcy5kZXN0aW5hdGlvbi5jIC0gdGhpcy5jKSAqIHBvc1xuICAgICAgLCBkOiB0aGlzLmQgKyAodGhpcy5kZXN0aW5hdGlvbi5kIC0gdGhpcy5kKSAqIHBvc1xuICAgICAgLCBlOiB0aGlzLmUgKyAodGhpcy5kZXN0aW5hdGlvbi5lIC0gdGhpcy5lKSAqIHBvc1xuICAgICAgLCBmOiB0aGlzLmYgKyAodGhpcy5kZXN0aW5hdGlvbi5mIC0gdGhpcy5mKSAqIHBvc1xuICAgICAgfSlcblxuICAgICAgLy8gcHJvY2VzcyBwYXJhbWV0cmljIHJvdGF0aW9uIGlmIHByZXNlbnRcbiAgICAgIGlmICh0aGlzLnBhcmFtICYmIHRoaXMucGFyYW0udG8pIHtcbiAgICAgICAgLy8gY2FsY3VsYXRlIGN1cnJlbnQgcGFyYW1ldHJpYyBwb3NpdGlvblxuICAgICAgICB2YXIgcGFyYW0gPSB7XG4gICAgICAgICAgcm90YXRpb246IHRoaXMucGFyYW0uZnJvbS5yb3RhdGlvbiArICh0aGlzLnBhcmFtLnRvLnJvdGF0aW9uIC0gdGhpcy5wYXJhbS5mcm9tLnJvdGF0aW9uKSAqIHBvc1xuICAgICAgICAsIGN4OiAgICAgICB0aGlzLnBhcmFtLmZyb20uY3hcbiAgICAgICAgLCBjeTogICAgICAgdGhpcy5wYXJhbS5mcm9tLmN5XG4gICAgICAgIH1cblxuICAgICAgICAvLyByb3RhdGUgbWF0cml4XG4gICAgICAgIG1hdHJpeCA9IG1hdHJpeC5yb3RhdGUoXG4gICAgICAgICAgKHRoaXMucGFyYW0udG8ucm90YXRpb24gLSB0aGlzLnBhcmFtLmZyb20ucm90YXRpb24gKiAyKSAqIHBvc1xuICAgICAgICAsIHBhcmFtLmN4XG4gICAgICAgICwgcGFyYW0uY3lcbiAgICAgICAgKVxuXG4gICAgICAgIC8vIHN0b3JlIGN1cnJlbnQgcGFyYW1ldHJpYyB2YWx1ZXNcbiAgICAgICAgbWF0cml4LnBhcmFtID0gcGFyYW1cbiAgICAgIH1cblxuICAgICAgcmV0dXJuIG1hdHJpeFxuICAgIH1cbiAgICAvLyBNdWx0aXBsaWVzIGJ5IGdpdmVuIG1hdHJpeFxuICAsIG11bHRpcGx5OiBmdW5jdGlvbihtYXRyaXgpIHtcbiAgICAgIHJldHVybiBuZXcgU1ZHLk1hdHJpeCh0aGlzLm5hdGl2ZSgpLm11bHRpcGx5KHBhcnNlTWF0cml4KG1hdHJpeCkubmF0aXZlKCkpKVxuICAgIH1cbiAgICAvLyBJbnZlcnNlcyBtYXRyaXhcbiAgLCBpbnZlcnNlOiBmdW5jdGlvbigpIHtcbiAgICAgIHJldHVybiBuZXcgU1ZHLk1hdHJpeCh0aGlzLm5hdGl2ZSgpLmludmVyc2UoKSlcbiAgICB9XG4gICAgLy8gVHJhbnNsYXRlIG1hdHJpeFxuICAsIHRyYW5zbGF0ZTogZnVuY3Rpb24oeCwgeSkge1xuICAgICAgcmV0dXJuIG5ldyBTVkcuTWF0cml4KHRoaXMubmF0aXZlKCkudHJhbnNsYXRlKHggfHwgMCwgeSB8fCAwKSlcbiAgICB9XG4gICAgLy8gU2NhbGUgbWF0cml4XG4gICwgc2NhbGU6IGZ1bmN0aW9uKHgsIHksIGN4LCBjeSkge1xuICAgICAgLy8gc3VwcG9ydCB1bmlmb3JtYWwgc2NhbGVcbiAgICAgIGlmIChhcmd1bWVudHMubGVuZ3RoID09IDEpIHtcbiAgICAgICAgeSA9IHhcbiAgICAgIH0gZWxzZSBpZiAoYXJndW1lbnRzLmxlbmd0aCA9PSAzKSB7XG4gICAgICAgIGN5ID0gY3hcbiAgICAgICAgY3ggPSB5XG4gICAgICAgIHkgPSB4XG4gICAgICB9XG5cbiAgICAgIHJldHVybiB0aGlzLmFyb3VuZChjeCwgY3ksIG5ldyBTVkcuTWF0cml4KHgsIDAsIDAsIHksIDAsIDApKVxuICAgIH1cbiAgICAvLyBSb3RhdGUgbWF0cml4XG4gICwgcm90YXRlOiBmdW5jdGlvbihyLCBjeCwgY3kpIHtcbiAgICAgIC8vIGNvbnZlcnQgZGVncmVlcyB0byByYWRpYW5zXG4gICAgICByID0gU1ZHLnV0aWxzLnJhZGlhbnMocilcblxuICAgICAgcmV0dXJuIHRoaXMuYXJvdW5kKGN4LCBjeSwgbmV3IFNWRy5NYXRyaXgoTWF0aC5jb3MociksIE1hdGguc2luKHIpLCAtTWF0aC5zaW4ociksIE1hdGguY29zKHIpLCAwLCAwKSlcbiAgICB9XG4gICAgLy8gRmxpcCBtYXRyaXggb24geCBvciB5LCBhdCBhIGdpdmVuIG9mZnNldFxuICAsIGZsaXA6IGZ1bmN0aW9uKGEsIG8pIHtcbiAgICAgIHJldHVybiBhID09ICd4JyA/IHRoaXMuc2NhbGUoLTEsIDEsIG8sIDApIDogdGhpcy5zY2FsZSgxLCAtMSwgMCwgbylcbiAgICB9XG4gICAgLy8gU2tld1xuICAsIHNrZXc6IGZ1bmN0aW9uKHgsIHksIGN4LCBjeSkge1xuICAgICAgLy8gc3VwcG9ydCB1bmlmb3JtYWwgc2tld1xuICAgICAgaWYgKGFyZ3VtZW50cy5sZW5ndGggPT0gMSkge1xuICAgICAgICB5ID0geFxuICAgICAgfSBlbHNlIGlmIChhcmd1bWVudHMubGVuZ3RoID09IDMpIHtcbiAgICAgICAgY3kgPSBjeFxuICAgICAgICBjeCA9IHlcbiAgICAgICAgeSA9IHhcbiAgICAgIH1cblxuICAgICAgLy8gY29udmVydCBkZWdyZWVzIHRvIHJhZGlhbnNcbiAgICAgIHggPSBTVkcudXRpbHMucmFkaWFucyh4KVxuICAgICAgeSA9IFNWRy51dGlscy5yYWRpYW5zKHkpXG5cbiAgICAgIHJldHVybiB0aGlzLmFyb3VuZChjeCwgY3ksIG5ldyBTVkcuTWF0cml4KDEsIE1hdGgudGFuKHkpLCBNYXRoLnRhbih4KSwgMSwgMCwgMCkpXG4gICAgfVxuICAgIC8vIFNrZXdYXG4gICwgc2tld1g6IGZ1bmN0aW9uKHgsIGN4LCBjeSkge1xuICAgICAgcmV0dXJuIHRoaXMuc2tldyh4LCAwLCBjeCwgY3kpXG4gICAgfVxuICAgIC8vIFNrZXdZXG4gICwgc2tld1k6IGZ1bmN0aW9uKHksIGN4LCBjeSkge1xuICAgICAgcmV0dXJuIHRoaXMuc2tldygwLCB5LCBjeCwgY3kpXG4gICAgfVxuICAgIC8vIFRyYW5zZm9ybSBhcm91bmQgYSBjZW50ZXIgcG9pbnRcbiAgLCBhcm91bmQ6IGZ1bmN0aW9uKGN4LCBjeSwgbWF0cml4KSB7XG4gICAgICByZXR1cm4gdGhpc1xuICAgICAgICAubXVsdGlwbHkobmV3IFNWRy5NYXRyaXgoMSwgMCwgMCwgMSwgY3ggfHwgMCwgY3kgfHwgMCkpXG4gICAgICAgIC5tdWx0aXBseShtYXRyaXgpXG4gICAgICAgIC5tdWx0aXBseShuZXcgU1ZHLk1hdHJpeCgxLCAwLCAwLCAxLCAtY3ggfHwgMCwgLWN5IHx8IDApKVxuICAgIH1cbiAgICAvLyBDb252ZXJ0IHRvIG5hdGl2ZSBTVkdNYXRyaXhcbiAgLCBuYXRpdmU6IGZ1bmN0aW9uKCkge1xuICAgICAgLy8gY3JlYXRlIG5ldyBtYXRyaXhcbiAgICAgIHZhciBtYXRyaXggPSBTVkcucGFyc2VyLm5hdGl2ZS5jcmVhdGVTVkdNYXRyaXgoKVxuXG4gICAgICAvLyB1cGRhdGUgd2l0aCBjdXJyZW50IHZhbHVlc1xuICAgICAgZm9yICh2YXIgaSA9IGFiY2RlZi5sZW5ndGggLSAxOyBpID49IDA7IGktLSlcbiAgICAgICAgbWF0cml4W2FiY2RlZltpXV0gPSB0aGlzW2FiY2RlZltpXV1cblxuICAgICAgcmV0dXJuIG1hdHJpeFxuICAgIH1cbiAgICAvLyBDb252ZXJ0IG1hdHJpeCB0byBzdHJpbmdcbiAgLCB0b1N0cmluZzogZnVuY3Rpb24oKSB7XG4gICAgICByZXR1cm4gJ21hdHJpeCgnICsgdGhpcy5hICsgJywnICsgdGhpcy5iICsgJywnICsgdGhpcy5jICsgJywnICsgdGhpcy5kICsgJywnICsgdGhpcy5lICsgJywnICsgdGhpcy5mICsgJyknXG4gICAgfVxuICB9XG5cbiAgLy8gRGVmaW5lIHBhcmVudFxuLCBwYXJlbnQ6IFNWRy5FbGVtZW50XG5cbiAgLy8gQWRkIHBhcmVudCBtZXRob2RcbiwgY29uc3RydWN0OiB7XG4gICAgLy8gR2V0IGN1cnJlbnQgbWF0cml4XG4gICAgY3RtOiBmdW5jdGlvbigpIHtcbiAgICAgIHJldHVybiBuZXcgU1ZHLk1hdHJpeCh0aGlzLm5vZGUuZ2V0Q1RNKCkpXG4gICAgfSxcbiAgICAvLyBHZXQgY3VycmVudCBzY3JlZW4gbWF0cml4XG4gICAgc2NyZWVuQ1RNOiBmdW5jdGlvbigpIHtcbiAgICAgIHJldHVybiBuZXcgU1ZHLk1hdHJpeCh0aGlzLm5vZGUuZ2V0U2NyZWVuQ1RNKCkpXG4gICAgfVxuXG4gIH1cblxufSlcblxuU1ZHLlBvaW50ID0gU1ZHLmludmVudCh7XG4gIC8vIEluaXRpYWxpemVcbiAgY3JlYXRlOiBmdW5jdGlvbih4LHkpIHtcbiAgICB2YXIgaSwgc291cmNlLCBiYXNlID0ge3g6MCwgeTowfVxuXG4gICAgLy8gZW5zdXJlIHNvdXJjZSBhcyBvYmplY3RcbiAgICBzb3VyY2UgPSBBcnJheS5pc0FycmF5KHgpID9cbiAgICAgIHt4OnhbMF0sIHk6eFsxXX0gOlxuICAgIHR5cGVvZiB4ID09PSAnb2JqZWN0JyA/XG4gICAgICB7eDp4LngsIHk6eC55fSA6XG4gICAgeCAhPSBudWxsID9cbiAgICAgIHt4OngsIHk6KHkgIT0gbnVsbCA/IHkgOiB4KX0gOiBiYXNlIC8vIElmIHkgaGFzIG5vIHZhbHVlLCB0aGVuIHggaXMgdXNlZCBoYXMgaXRzIHZhbHVlXG5cbiAgICAvLyBtZXJnZSBzb3VyY2VcbiAgICB0aGlzLnggPSBzb3VyY2UueFxuICAgIHRoaXMueSA9IHNvdXJjZS55XG4gIH1cblxuICAvLyBBZGQgbWV0aG9kc1xuLCBleHRlbmQ6IHtcbiAgICAvLyBDbG9uZSBwb2ludFxuICAgIGNsb25lOiBmdW5jdGlvbigpIHtcbiAgICAgIHJldHVybiBuZXcgU1ZHLlBvaW50KHRoaXMpXG4gICAgfVxuICAgIC8vIE1vcnBoIG9uZSBwb2ludCBpbnRvIGFub3RoZXJcbiAgLCBtb3JwaDogZnVuY3Rpb24oeCwgeSkge1xuICAgICAgLy8gc3RvcmUgbmV3IGRlc3RpbmF0aW9uXG4gICAgICB0aGlzLmRlc3RpbmF0aW9uID0gbmV3IFNWRy5Qb2ludCh4LCB5KVxuXG4gICAgICByZXR1cm4gdGhpc1xuICAgIH1cbiAgICAvLyBHZXQgbW9ycGhlZCBwb2ludCBhdCBhIGdpdmVuIHBvc2l0aW9uXG4gICwgYXQ6IGZ1bmN0aW9uKHBvcykge1xuICAgICAgLy8gbWFrZSBzdXJlIGEgZGVzdGluYXRpb24gaXMgZGVmaW5lZFxuICAgICAgaWYgKCF0aGlzLmRlc3RpbmF0aW9uKSByZXR1cm4gdGhpc1xuXG4gICAgICAvLyBjYWxjdWxhdGUgbW9ycGhlZCBtYXRyaXggYXQgYSBnaXZlbiBwb3NpdGlvblxuICAgICAgdmFyIHBvaW50ID0gbmV3IFNWRy5Qb2ludCh7XG4gICAgICAgIHg6IHRoaXMueCArICh0aGlzLmRlc3RpbmF0aW9uLnggLSB0aGlzLngpICogcG9zXG4gICAgICAsIHk6IHRoaXMueSArICh0aGlzLmRlc3RpbmF0aW9uLnkgLSB0aGlzLnkpICogcG9zXG4gICAgICB9KVxuXG4gICAgICByZXR1cm4gcG9pbnRcbiAgICB9XG4gICAgLy8gQ29udmVydCB0byBuYXRpdmUgU1ZHUG9pbnRcbiAgLCBuYXRpdmU6IGZ1bmN0aW9uKCkge1xuICAgICAgLy8gY3JlYXRlIG5ldyBwb2ludFxuICAgICAgdmFyIHBvaW50ID0gU1ZHLnBhcnNlci5uYXRpdmUuY3JlYXRlU1ZHUG9pbnQoKVxuXG4gICAgICAvLyB1cGRhdGUgd2l0aCBjdXJyZW50IHZhbHVlc1xuICAgICAgcG9pbnQueCA9IHRoaXMueFxuICAgICAgcG9pbnQueSA9IHRoaXMueVxuXG4gICAgICByZXR1cm4gcG9pbnRcbiAgICB9XG4gICAgLy8gdHJhbnNmb3JtIHBvaW50IHdpdGggbWF0cml4XG4gICwgdHJhbnNmb3JtOiBmdW5jdGlvbihtYXRyaXgpIHtcbiAgICAgIHJldHVybiBuZXcgU1ZHLlBvaW50KHRoaXMubmF0aXZlKCkubWF0cml4VHJhbnNmb3JtKG1hdHJpeC5uYXRpdmUoKSkpXG4gICAgfVxuXG4gIH1cblxufSlcblxuU1ZHLmV4dGVuZChTVkcuRWxlbWVudCwge1xuXG4gIC8vIEdldCBwb2ludFxuICBwb2ludDogZnVuY3Rpb24oeCwgeSkge1xuICAgIHJldHVybiBuZXcgU1ZHLlBvaW50KHgseSkudHJhbnNmb3JtKHRoaXMuc2NyZWVuQ1RNKCkuaW52ZXJzZSgpKTtcbiAgfVxuXG59KVxuXG5TVkcuZXh0ZW5kKFNWRy5FbGVtZW50LCB7XG4gIC8vIFNldCBzdmcgZWxlbWVudCBhdHRyaWJ1dGVcbiAgYXR0cjogZnVuY3Rpb24oYSwgdiwgbikge1xuICAgIC8vIGFjdCBhcyBmdWxsIGdldHRlclxuICAgIGlmIChhID09IG51bGwpIHtcbiAgICAgIC8vIGdldCBhbiBvYmplY3Qgb2YgYXR0cmlidXRlc1xuICAgICAgYSA9IHt9XG4gICAgICB2ID0gdGhpcy5ub2RlLmF0dHJpYnV0ZXNcbiAgICAgIGZvciAobiA9IHYubGVuZ3RoIC0gMTsgbiA+PSAwOyBuLS0pXG4gICAgICAgIGFbdltuXS5ub2RlTmFtZV0gPSBTVkcucmVnZXguaXNOdW1iZXIudGVzdCh2W25dLm5vZGVWYWx1ZSkgPyBwYXJzZUZsb2F0KHZbbl0ubm9kZVZhbHVlKSA6IHZbbl0ubm9kZVZhbHVlXG5cbiAgICAgIHJldHVybiBhXG5cbiAgICB9IGVsc2UgaWYgKHR5cGVvZiBhID09ICdvYmplY3QnKSB7XG4gICAgICAvLyBhcHBseSBldmVyeSBhdHRyaWJ1dGUgaW5kaXZpZHVhbGx5IGlmIGFuIG9iamVjdCBpcyBwYXNzZWRcbiAgICAgIGZvciAodiBpbiBhKSB0aGlzLmF0dHIodiwgYVt2XSlcblxuICAgIH0gZWxzZSBpZiAodiA9PT0gbnVsbCkge1xuICAgICAgICAvLyByZW1vdmUgdmFsdWVcbiAgICAgICAgdGhpcy5ub2RlLnJlbW92ZUF0dHJpYnV0ZShhKVxuXG4gICAgfSBlbHNlIGlmICh2ID09IG51bGwpIHtcbiAgICAgIC8vIGFjdCBhcyBhIGdldHRlciBpZiB0aGUgZmlyc3QgYW5kIG9ubHkgYXJndW1lbnQgaXMgbm90IGFuIG9iamVjdFxuICAgICAgdiA9IHRoaXMubm9kZS5nZXRBdHRyaWJ1dGUoYSlcbiAgICAgIHJldHVybiB2ID09IG51bGwgP1xuICAgICAgICBTVkcuZGVmYXVsdHMuYXR0cnNbYV0gOlxuICAgICAgU1ZHLnJlZ2V4LmlzTnVtYmVyLnRlc3QodikgP1xuICAgICAgICBwYXJzZUZsb2F0KHYpIDogdlxuXG4gICAgfSBlbHNlIHtcbiAgICAgIC8vIEJVRyBGSVg6IHNvbWUgYnJvd3NlcnMgd2lsbCByZW5kZXIgYSBzdHJva2UgaWYgYSBjb2xvciBpcyBnaXZlbiBldmVuIHRob3VnaCBzdHJva2Ugd2lkdGggaXMgMFxuICAgICAgaWYgKGEgPT0gJ3N0cm9rZS13aWR0aCcpXG4gICAgICAgIHRoaXMuYXR0cignc3Ryb2tlJywgcGFyc2VGbG9hdCh2KSA+IDAgPyB0aGlzLl9zdHJva2UgOiBudWxsKVxuICAgICAgZWxzZSBpZiAoYSA9PSAnc3Ryb2tlJylcbiAgICAgICAgdGhpcy5fc3Ryb2tlID0gdlxuXG4gICAgICAvLyBjb252ZXJ0IGltYWdlIGZpbGwgYW5kIHN0cm9rZSB0byBwYXR0ZXJuc1xuICAgICAgaWYgKGEgPT0gJ2ZpbGwnIHx8IGEgPT0gJ3N0cm9rZScpIHtcbiAgICAgICAgaWYgKFNWRy5yZWdleC5pc0ltYWdlLnRlc3QodikpXG4gICAgICAgICAgdiA9IHRoaXMuZG9jKCkuZGVmcygpLmltYWdlKHYsIDAsIDApXG5cbiAgICAgICAgaWYgKHYgaW5zdGFuY2VvZiBTVkcuSW1hZ2UpXG4gICAgICAgICAgdiA9IHRoaXMuZG9jKCkuZGVmcygpLnBhdHRlcm4oMCwgMCwgZnVuY3Rpb24oKSB7XG4gICAgICAgICAgICB0aGlzLmFkZCh2KVxuICAgICAgICAgIH0pXG4gICAgICB9XG5cbiAgICAgIC8vIGVuc3VyZSBjb3JyZWN0IG51bWVyaWMgdmFsdWVzIChhbHNvIGFjY2VwdHMgTmFOIGFuZCBJbmZpbml0eSlcbiAgICAgIGlmICh0eXBlb2YgdiA9PT0gJ251bWJlcicpXG4gICAgICAgIHYgPSBuZXcgU1ZHLk51bWJlcih2KVxuXG4gICAgICAvLyBlbnN1cmUgZnVsbCBoZXggY29sb3JcbiAgICAgIGVsc2UgaWYgKFNWRy5Db2xvci5pc0NvbG9yKHYpKVxuICAgICAgICB2ID0gbmV3IFNWRy5Db2xvcih2KVxuXG4gICAgICAvLyBwYXJzZSBhcnJheSB2YWx1ZXNcbiAgICAgIGVsc2UgaWYgKEFycmF5LmlzQXJyYXkodikpXG4gICAgICAgIHYgPSBuZXcgU1ZHLkFycmF5KHYpXG5cbiAgICAgIC8vIHN0b3JlIHBhcmFtZXRyaWMgdHJhbnNmb3JtYXRpb24gdmFsdWVzIGxvY2FsbHlcbiAgICAgIGVsc2UgaWYgKHYgaW5zdGFuY2VvZiBTVkcuTWF0cml4ICYmIHYucGFyYW0pXG4gICAgICAgIHRoaXMucGFyYW0gPSB2LnBhcmFtXG5cbiAgICAgIC8vIGlmIHRoZSBwYXNzZWQgYXR0cmlidXRlIGlzIGxlYWRpbmcuLi5cbiAgICAgIGlmIChhID09ICdsZWFkaW5nJykge1xuICAgICAgICAvLyAuLi4gY2FsbCB0aGUgbGVhZGluZyBtZXRob2QgaW5zdGVhZFxuICAgICAgICBpZiAodGhpcy5sZWFkaW5nKVxuICAgICAgICAgIHRoaXMubGVhZGluZyh2KVxuICAgICAgfSBlbHNlIHtcbiAgICAgICAgLy8gc2V0IGdpdmVuIGF0dHJpYnV0ZSBvbiBub2RlXG4gICAgICAgIHR5cGVvZiBuID09PSAnc3RyaW5nJyA/XG4gICAgICAgICAgdGhpcy5ub2RlLnNldEF0dHJpYnV0ZU5TKG4sIGEsIHYudG9TdHJpbmcoKSkgOlxuICAgICAgICAgIHRoaXMubm9kZS5zZXRBdHRyaWJ1dGUoYSwgdi50b1N0cmluZygpKVxuICAgICAgfVxuXG4gICAgICAvLyByZWJ1aWxkIGlmIHJlcXVpcmVkXG4gICAgICBpZiAodGhpcy5yZWJ1aWxkICYmIChhID09ICdmb250LXNpemUnIHx8IGEgPT0gJ3gnKSlcbiAgICAgICAgdGhpcy5yZWJ1aWxkKGEsIHYpXG4gICAgfVxuXG4gICAgcmV0dXJuIHRoaXNcbiAgfVxufSlcblNWRy5leHRlbmQoU1ZHLkVsZW1lbnQsIHtcbiAgLy8gQWRkIHRyYW5zZm9ybWF0aW9uc1xuICB0cmFuc2Zvcm06IGZ1bmN0aW9uKG8sIHJlbGF0aXZlKSB7XG4gICAgLy8gZ2V0IHRhcmdldCBpbiBjYXNlIG9mIHRoZSBmeCBtb2R1bGUsIG90aGVyd2lzZSByZWZlcmVuY2UgdGhpc1xuICAgIHZhciB0YXJnZXQgPSB0aGlzXG4gICAgICAsIG1hdHJpeFxuXG4gICAgLy8gYWN0IGFzIGEgZ2V0dGVyXG4gICAgaWYgKHR5cGVvZiBvICE9PSAnb2JqZWN0Jykge1xuICAgICAgLy8gZ2V0IGN1cnJlbnQgbWF0cml4XG4gICAgICBtYXRyaXggPSBuZXcgU1ZHLk1hdHJpeCh0YXJnZXQpLmV4dHJhY3QoKVxuXG4gICAgICByZXR1cm4gdHlwZW9mIG8gPT09ICdzdHJpbmcnID8gbWF0cml4W29dIDogbWF0cml4XG4gICAgfVxuXG4gICAgLy8gZ2V0IGN1cnJlbnQgbWF0cml4XG4gICAgbWF0cml4ID0gbmV3IFNWRy5NYXRyaXgodGFyZ2V0KVxuXG4gICAgLy8gZW5zdXJlIHJlbGF0aXZlIGZsYWdcbiAgICByZWxhdGl2ZSA9ICEhcmVsYXRpdmUgfHwgISFvLnJlbGF0aXZlXG5cbiAgICAvLyBhY3Qgb24gbWF0cml4XG4gICAgaWYgKG8uYSAhPSBudWxsKSB7XG4gICAgICBtYXRyaXggPSByZWxhdGl2ZSA/XG4gICAgICAgIC8vIHJlbGF0aXZlXG4gICAgICAgIG1hdHJpeC5tdWx0aXBseShuZXcgU1ZHLk1hdHJpeChvKSkgOlxuICAgICAgICAvLyBhYnNvbHV0ZVxuICAgICAgICBuZXcgU1ZHLk1hdHJpeChvKVxuXG4gICAgLy8gYWN0IG9uIHJvdGF0aW9uXG4gICAgfSBlbHNlIGlmIChvLnJvdGF0aW9uICE9IG51bGwpIHtcbiAgICAgIC8vIGVuc3VyZSBjZW50cmUgcG9pbnRcbiAgICAgIGVuc3VyZUNlbnRyZShvLCB0YXJnZXQpXG5cbiAgICAgIC8vIGFwcGx5IHRyYW5zZm9ybWF0aW9uXG4gICAgICBtYXRyaXggPSByZWxhdGl2ZSA/XG4gICAgICAgIC8vIHJlbGF0aXZlXG4gICAgICAgIG1hdHJpeC5yb3RhdGUoby5yb3RhdGlvbiwgby5jeCwgby5jeSkgOlxuICAgICAgICAvLyBhYnNvbHV0ZVxuICAgICAgICBtYXRyaXgucm90YXRlKG8ucm90YXRpb24gLSBtYXRyaXguZXh0cmFjdCgpLnJvdGF0aW9uLCBvLmN4LCBvLmN5KVxuXG4gICAgLy8gYWN0IG9uIHNjYWxlXG4gICAgfSBlbHNlIGlmIChvLnNjYWxlICE9IG51bGwgfHwgby5zY2FsZVggIT0gbnVsbCB8fCBvLnNjYWxlWSAhPSBudWxsKSB7XG4gICAgICAvLyBlbnN1cmUgY2VudHJlIHBvaW50XG4gICAgICBlbnN1cmVDZW50cmUobywgdGFyZ2V0KVxuXG4gICAgICAvLyBlbnN1cmUgc2NhbGUgdmFsdWVzIG9uIGJvdGggYXhlc1xuICAgICAgby5zY2FsZVggPSBvLnNjYWxlICE9IG51bGwgPyBvLnNjYWxlIDogby5zY2FsZVggIT0gbnVsbCA/IG8uc2NhbGVYIDogMVxuICAgICAgby5zY2FsZVkgPSBvLnNjYWxlICE9IG51bGwgPyBvLnNjYWxlIDogby5zY2FsZVkgIT0gbnVsbCA/IG8uc2NhbGVZIDogMVxuXG4gICAgICBpZiAoIXJlbGF0aXZlKSB7XG4gICAgICAgIC8vIGFic29sdXRlOyBtdWx0aXBseSBpbnZlcnNlZCB2YWx1ZXNcbiAgICAgICAgdmFyIGUgPSBtYXRyaXguZXh0cmFjdCgpXG4gICAgICAgIG8uc2NhbGVYID0gby5zY2FsZVggKiAxIC8gZS5zY2FsZVhcbiAgICAgICAgby5zY2FsZVkgPSBvLnNjYWxlWSAqIDEgLyBlLnNjYWxlWVxuICAgICAgfVxuXG4gICAgICBtYXRyaXggPSBtYXRyaXguc2NhbGUoby5zY2FsZVgsIG8uc2NhbGVZLCBvLmN4LCBvLmN5KVxuXG4gICAgLy8gYWN0IG9uIHNrZXdcbiAgICB9IGVsc2UgaWYgKG8uc2tldyAhPSBudWxsIHx8IG8uc2tld1ggIT0gbnVsbCB8fCBvLnNrZXdZICE9IG51bGwpIHtcbiAgICAgIC8vIGVuc3VyZSBjZW50cmUgcG9pbnRcbiAgICAgIGVuc3VyZUNlbnRyZShvLCB0YXJnZXQpXG5cbiAgICAgIC8vIGVuc3VyZSBza2V3IHZhbHVlcyBvbiBib3RoIGF4ZXNcbiAgICAgIG8uc2tld1ggPSBvLnNrZXcgIT0gbnVsbCA/IG8uc2tldyA6IG8uc2tld1ggIT0gbnVsbCA/IG8uc2tld1ggOiAwXG4gICAgICBvLnNrZXdZID0gby5za2V3ICE9IG51bGwgPyBvLnNrZXcgOiBvLnNrZXdZICE9IG51bGwgPyBvLnNrZXdZIDogMFxuXG4gICAgICBpZiAoIXJlbGF0aXZlKSB7XG4gICAgICAgIC8vIGFic29sdXRlOyByZXNldCBza2V3IHZhbHVlc1xuICAgICAgICB2YXIgZSA9IG1hdHJpeC5leHRyYWN0KClcbiAgICAgICAgbWF0cml4ID0gbWF0cml4Lm11bHRpcGx5KG5ldyBTVkcuTWF0cml4KCkuc2tldyhlLnNrZXdYLCBlLnNrZXdZLCBvLmN4LCBvLmN5KS5pbnZlcnNlKCkpXG4gICAgICB9XG5cbiAgICAgIG1hdHJpeCA9IG1hdHJpeC5za2V3KG8uc2tld1gsIG8uc2tld1ksIG8uY3gsIG8uY3kpXG5cbiAgICAvLyBhY3Qgb24gZmxpcFxuICAgIH0gZWxzZSBpZiAoby5mbGlwKSB7XG4gICAgICBtYXRyaXggPSBtYXRyaXguZmxpcChcbiAgICAgICAgby5mbGlwXG4gICAgICAsIG8ub2Zmc2V0ID09IG51bGwgPyB0YXJnZXQuYmJveCgpWydjJyArIG8uZmxpcF0gOiBvLm9mZnNldFxuICAgICAgKVxuXG4gICAgLy8gYWN0IG9uIHRyYW5zbGF0ZVxuICAgIH0gZWxzZSBpZiAoby54ICE9IG51bGwgfHwgby55ICE9IG51bGwpIHtcbiAgICAgIGlmIChyZWxhdGl2ZSkge1xuICAgICAgICAvLyByZWxhdGl2ZVxuICAgICAgICBtYXRyaXggPSBtYXRyaXgudHJhbnNsYXRlKG8ueCwgby55KVxuICAgICAgfSBlbHNlIHtcbiAgICAgICAgLy8gYWJzb2x1dGVcbiAgICAgICAgaWYgKG8ueCAhPSBudWxsKSBtYXRyaXguZSA9IG8ueFxuICAgICAgICBpZiAoby55ICE9IG51bGwpIG1hdHJpeC5mID0gby55XG4gICAgICB9XG4gICAgfVxuXG4gICAgcmV0dXJuIHRoaXMuYXR0cigndHJhbnNmb3JtJywgbWF0cml4KVxuICB9XG59KVxuXG5TVkcuZXh0ZW5kKFNWRy5GWCwge1xuICB0cmFuc2Zvcm06IGZ1bmN0aW9uKG8sIHJlbGF0aXZlKSB7XG4gICAgLy8gZ2V0IHRhcmdldCBpbiBjYXNlIG9mIHRoZSBmeCBtb2R1bGUsIG90aGVyd2lzZSByZWZlcmVuY2UgdGhpc1xuICAgIHZhciB0YXJnZXQgPSB0aGlzLnRhcmdldCgpXG4gICAgICAsIG1hdHJpeFxuXG4gICAgLy8gYWN0IGFzIGEgZ2V0dGVyXG4gICAgaWYgKHR5cGVvZiBvICE9PSAnb2JqZWN0Jykge1xuICAgICAgLy8gZ2V0IGN1cnJlbnQgbWF0cml4XG4gICAgICBtYXRyaXggPSBuZXcgU1ZHLk1hdHJpeCh0YXJnZXQpLmV4dHJhY3QoKVxuXG4gICAgICByZXR1cm4gdHlwZW9mIG8gPT09ICdzdHJpbmcnID8gbWF0cml4W29dIDogbWF0cml4XG4gICAgfVxuXG4gICAgLy8gZW5zdXJlIHJlbGF0aXZlIGZsYWdcbiAgICByZWxhdGl2ZSA9ICEhcmVsYXRpdmUgfHwgISFvLnJlbGF0aXZlXG5cbiAgICAvLyBhY3Qgb24gbWF0cml4XG4gICAgaWYgKG8uYSAhPSBudWxsKSB7XG4gICAgICBtYXRyaXggPSBuZXcgU1ZHLk1hdHJpeChvKVxuXG4gICAgLy8gYWN0IG9uIHJvdGF0aW9uXG4gICAgfSBlbHNlIGlmIChvLnJvdGF0aW9uICE9IG51bGwpIHtcbiAgICAgIC8vIGVuc3VyZSBjZW50cmUgcG9pbnRcbiAgICAgIGVuc3VyZUNlbnRyZShvLCB0YXJnZXQpXG5cbiAgICAgIC8vIGFwcGx5IHRyYW5zZm9ybWF0aW9uXG4gICAgICBtYXRyaXggPSBuZXcgU1ZHLlJvdGF0ZShvLnJvdGF0aW9uLCBvLmN4LCBvLmN5KVxuXG4gICAgLy8gYWN0IG9uIHNjYWxlXG4gICAgfSBlbHNlIGlmIChvLnNjYWxlICE9IG51bGwgfHwgby5zY2FsZVggIT0gbnVsbCB8fCBvLnNjYWxlWSAhPSBudWxsKSB7XG4gICAgICAvLyBlbnN1cmUgY2VudHJlIHBvaW50XG4gICAgICBlbnN1cmVDZW50cmUobywgdGFyZ2V0KVxuXG4gICAgICAvLyBlbnN1cmUgc2NhbGUgdmFsdWVzIG9uIGJvdGggYXhlc1xuICAgICAgby5zY2FsZVggPSBvLnNjYWxlICE9IG51bGwgPyBvLnNjYWxlIDogby5zY2FsZVggIT0gbnVsbCA/IG8uc2NhbGVYIDogMVxuICAgICAgby5zY2FsZVkgPSBvLnNjYWxlICE9IG51bGwgPyBvLnNjYWxlIDogby5zY2FsZVkgIT0gbnVsbCA/IG8uc2NhbGVZIDogMVxuXG4gICAgICBtYXRyaXggPSBuZXcgU1ZHLlNjYWxlKG8uc2NhbGVYLCBvLnNjYWxlWSwgby5jeCwgby5jeSlcblxuICAgIC8vIGFjdCBvbiBza2V3XG4gICAgfSBlbHNlIGlmIChvLnNrZXdYICE9IG51bGwgfHwgby5za2V3WSAhPSBudWxsKSB7XG4gICAgICAvLyBlbnN1cmUgY2VudHJlIHBvaW50XG4gICAgICBlbnN1cmVDZW50cmUobywgdGFyZ2V0KVxuXG4gICAgICAvLyBlbnN1cmUgc2tldyB2YWx1ZXMgb24gYm90aCBheGVzXG4gICAgICBvLnNrZXdYID0gby5za2V3WCAhPSBudWxsID8gby5za2V3WCA6IDBcbiAgICAgIG8uc2tld1kgPSBvLnNrZXdZICE9IG51bGwgPyBvLnNrZXdZIDogMFxuXG4gICAgICBtYXRyaXggPSBuZXcgU1ZHLlNrZXcoby5za2V3WCwgby5za2V3WSwgby5jeCwgby5jeSlcblxuICAgIC8vIGFjdCBvbiBmbGlwXG4gICAgfSBlbHNlIGlmIChvLmZsaXApIHtcbiAgICAgIG1hdHJpeCA9IG5ldyBTVkcuTWF0cml4KCkubW9ycGgobmV3IFNWRy5NYXRyaXgoKS5mbGlwKFxuICAgICAgICBvLmZsaXBcbiAgICAgICwgby5vZmZzZXQgPT0gbnVsbCA/IHRhcmdldC5iYm94KClbJ2MnICsgby5mbGlwXSA6IG8ub2Zmc2V0XG4gICAgICApKVxuXG4gICAgLy8gYWN0IG9uIHRyYW5zbGF0ZVxuICAgIH0gZWxzZSBpZiAoby54ICE9IG51bGwgfHwgby55ICE9IG51bGwpIHtcbiAgICAgIG1hdHJpeCA9IG5ldyBTVkcuVHJhbnNsYXRlKG8ueCwgby55KVxuICAgIH1cblxuICAgIGlmKCFtYXRyaXgpIHJldHVybiB0aGlzXG5cbiAgICBtYXRyaXgucmVsYXRpdmUgPSByZWxhdGl2ZVxuXG4gICAgdGhpcy5sYXN0KCkudHJhbnNmb3Jtcy5wdXNoKG1hdHJpeClcblxuICAgIHNldFRpbWVvdXQoZnVuY3Rpb24oKXt0aGlzLnN0YXJ0KCl9LmJpbmQodGhpcyksIDApXG5cbiAgICByZXR1cm4gdGhpc1xuICB9XG59KVxuXG5TVkcuZXh0ZW5kKFNWRy5FbGVtZW50LCB7XG4gIC8vIFJlc2V0IGFsbCB0cmFuc2Zvcm1hdGlvbnNcbiAgdW50cmFuc2Zvcm06IGZ1bmN0aW9uKCkge1xuICAgIHJldHVybiB0aGlzLmF0dHIoJ3RyYW5zZm9ybScsIG51bGwpXG4gIH0sXG4gIC8vIG1lcmdlIHRoZSB3aG9sZSB0cmFuc2Zvcm1hdGlvbiBjaGFpbiBpbnRvIG9uZSBtYXRyaXggYW5kIHJldHVybnMgaXRcbiAgbWF0cml4aWZ5OiBmdW5jdGlvbigpIHtcblxuICAgIHZhciBtYXRyaXggPSAodGhpcy5hdHRyKCd0cmFuc2Zvcm0nKSB8fCAnJylcbiAgICAgIC8vIHNwbGl0IHRyYW5zZm9ybWF0aW9uc1xuICAgICAgLnNwbGl0KC9cXClcXHMqLD9cXHMqLykuc2xpY2UoMCwtMSkubWFwKGZ1bmN0aW9uKHN0cil7XG4gICAgICAgIC8vIGdlbmVyYXRlIGtleSA9PiB2YWx1ZSBwYWlyc1xuICAgICAgICB2YXIga3YgPSBzdHIudHJpbSgpLnNwbGl0KCcoJylcbiAgICAgICAgcmV0dXJuIFtrdlswXSwga3ZbMV0uc3BsaXQoU1ZHLnJlZ2V4Lm1hdHJpeEVsZW1lbnRzKS5tYXAoZnVuY3Rpb24oc3RyKXsgcmV0dXJuIHBhcnNlRmxvYXQoc3RyKSB9KV1cbiAgICAgIH0pXG4gICAgICAvLyBjYWxjdWxhdGUgZXZlcnkgdHJhbnNmb3JtYXRpb24gaW50byBvbmUgbWF0cml4XG4gICAgICAucmVkdWNlKGZ1bmN0aW9uKG1hdHJpeCwgdHJhbnNmb3JtKXtcblxuICAgICAgICBpZih0cmFuc2Zvcm1bMF0gPT0gJ21hdHJpeCcpIHJldHVybiBtYXRyaXgubXVsdGlwbHkoYXJyYXlUb01hdHJpeCh0cmFuc2Zvcm1bMV0pKVxuICAgICAgICByZXR1cm4gbWF0cml4W3RyYW5zZm9ybVswXV0uYXBwbHkobWF0cml4LCB0cmFuc2Zvcm1bMV0pXG5cbiAgICAgIH0sIG5ldyBTVkcuTWF0cml4KCkpXG5cbiAgICByZXR1cm4gbWF0cml4XG4gIH0sXG4gIC8vIGFkZCBhbiBlbGVtZW50IHRvIGFub3RoZXIgcGFyZW50IHdpdGhvdXQgY2hhbmdpbmcgdGhlIHZpc3VhbCByZXByZXNlbnRhdGlvbiBvbiB0aGUgc2NyZWVuXG4gIHRvUGFyZW50OiBmdW5jdGlvbihwYXJlbnQpIHtcbiAgICBpZih0aGlzID09IHBhcmVudCkgcmV0dXJuIHRoaXNcbiAgICB2YXIgY3RtID0gdGhpcy5zY3JlZW5DVE0oKVxuICAgIHZhciB0ZW1wID0gcGFyZW50LnJlY3QoMSwxKVxuICAgIHZhciBwQ3RtID0gdGVtcC5zY3JlZW5DVE0oKS5pbnZlcnNlKClcbiAgICB0ZW1wLnJlbW92ZSgpXG5cbiAgICB0aGlzLmFkZFRvKHBhcmVudCkudW50cmFuc2Zvcm0oKS50cmFuc2Zvcm0ocEN0bS5tdWx0aXBseShjdG0pKVxuXG4gICAgcmV0dXJuIHRoaXNcbiAgfSxcbiAgLy8gc2FtZSBhcyBhYm92ZSB3aXRoIHBhcmVudCBlcXVhbHMgcm9vdC1zdmdcbiAgdG9Eb2M6IGZ1bmN0aW9uKCkge1xuICAgIHJldHVybiB0aGlzLnRvUGFyZW50KHRoaXMuZG9jKCkpXG4gIH1cblxufSlcblxuU1ZHLlRyYW5zZm9ybWF0aW9uID0gU1ZHLmludmVudCh7XG5cbiAgY3JlYXRlOiBmdW5jdGlvbihzb3VyY2UsIGludmVyc2VkKXtcblxuICAgIGlmKGFyZ3VtZW50cy5sZW5ndGggPiAxICYmIHR5cGVvZiBpbnZlcnNlZCAhPSAnYm9vbGVhbicpe1xuICAgICAgcmV0dXJuIHRoaXMuY3JlYXRlKFtdLnNsaWNlLmNhbGwoYXJndW1lbnRzKSlcbiAgICB9XG5cbiAgICBpZih0eXBlb2Ygc291cmNlID09ICdvYmplY3QnKXtcbiAgICAgIGZvcih2YXIgaSA9IDAsIGxlbiA9IHRoaXMuYXJndW1lbnRzLmxlbmd0aDsgaSA8IGxlbjsgKytpKXtcbiAgICAgICAgdGhpc1t0aGlzLmFyZ3VtZW50c1tpXV0gPSBzb3VyY2VbdGhpcy5hcmd1bWVudHNbaV1dXG4gICAgICB9XG4gICAgfVxuXG4gICAgaWYoQXJyYXkuaXNBcnJheShzb3VyY2UpKXtcbiAgICAgIGZvcih2YXIgaSA9IDAsIGxlbiA9IHRoaXMuYXJndW1lbnRzLmxlbmd0aDsgaSA8IGxlbjsgKytpKXtcbiAgICAgICAgdGhpc1t0aGlzLmFyZ3VtZW50c1tpXV0gPSBzb3VyY2VbaV1cbiAgICAgIH1cbiAgICB9XG5cbiAgICB0aGlzLmludmVyc2VkID0gZmFsc2VcblxuICAgIGlmKGludmVyc2VkID09PSB0cnVlKXtcbiAgICAgIHRoaXMuaW52ZXJzZWQgPSB0cnVlXG4gICAgfVxuXG4gIH1cblxuLCBleHRlbmQ6IHtcblxuICAgIGF0OiBmdW5jdGlvbihwb3Mpe1xuXG4gICAgICB2YXIgcGFyYW1zID0gW11cblxuICAgICAgZm9yKHZhciBpID0gMCwgbGVuID0gdGhpcy5hcmd1bWVudHMubGVuZ3RoOyBpIDwgbGVuOyArK2kpe1xuICAgICAgICBwYXJhbXMucHVzaCh0aGlzW3RoaXMuYXJndW1lbnRzW2ldXSlcbiAgICAgIH1cblxuICAgICAgdmFyIG0gPSB0aGlzLl91bmRvIHx8IG5ldyBTVkcuTWF0cml4KClcblxuICAgICAgbSA9IG5ldyBTVkcuTWF0cml4KCkubW9ycGgoU1ZHLk1hdHJpeC5wcm90b3R5cGVbdGhpcy5tZXRob2RdLmFwcGx5KG0sIHBhcmFtcykpLmF0KHBvcylcblxuICAgICAgcmV0dXJuIHRoaXMuaW52ZXJzZWQgPyBtLmludmVyc2UoKSA6IG1cblxuICAgIH1cblxuICAsIHVuZG86IGZ1bmN0aW9uKG8pe1xuICAgICAgZm9yKHZhciBpID0gMCwgbGVuID0gdGhpcy5hcmd1bWVudHMubGVuZ3RoOyBpIDwgbGVuOyArK2kpe1xuICAgICAgICBvW3RoaXMuYXJndW1lbnRzW2ldXSA9IHR5cGVvZiB0aGlzW3RoaXMuYXJndW1lbnRzW2ldXSA9PSAndW5kZWZpbmVkJyA/IDAgOiBvW3RoaXMuYXJndW1lbnRzW2ldXVxuICAgICAgfVxuXG4gICAgICAvLyBUaGUgbWV0aG9kIFNWRy5NYXRyaXguZXh0cmFjdCB3aGljaCB3YXMgdXNlZCBiZWZvcmUgY2FsbGluZyB0aGlzXG4gICAgICAvLyBtZXRob2QgdG8gb2J0YWluIGEgdmFsdWUgZm9yIHRoZSBwYXJhbWV0ZXIgbyBkb2Vzbid0IHJldHVybiBhIGN4IGFuZFxuICAgICAgLy8gYSBjeSBzbyB3ZSB1c2UgdGhlIG9uZXMgdGhhdCB3ZXJlIHByb3ZpZGVkIHRvIHRoaXMgb2JqZWN0IGF0IGl0cyBjcmVhdGlvblxuICAgICAgby5jeCA9IHRoaXMuY3hcbiAgICAgIG8uY3kgPSB0aGlzLmN5XG5cbiAgICAgIHRoaXMuX3VuZG8gPSBuZXcgU1ZHW2NhcGl0YWxpemUodGhpcy5tZXRob2QpXShvLCB0cnVlKS5hdCgxKVxuXG4gICAgICByZXR1cm4gdGhpc1xuICAgIH1cblxuICB9XG5cbn0pXG5cblNWRy5UcmFuc2xhdGUgPSBTVkcuaW52ZW50KHtcblxuICBwYXJlbnQ6IFNWRy5NYXRyaXhcbiwgaW5oZXJpdDogU1ZHLlRyYW5zZm9ybWF0aW9uXG5cbiwgY3JlYXRlOiBmdW5jdGlvbihzb3VyY2UsIGludmVyc2VkKXtcbiAgICBpZih0eXBlb2Ygc291cmNlID09ICdvYmplY3QnKSB0aGlzLmNvbnN0cnVjdG9yLmNhbGwodGhpcywgc291cmNlLCBpbnZlcnNlZClcbiAgICBlbHNlIHRoaXMuY29uc3RydWN0b3IuY2FsbCh0aGlzLCBbXS5zbGljZS5jYWxsKGFyZ3VtZW50cykpXG4gIH1cblxuLCBleHRlbmQ6IHtcbiAgICBhcmd1bWVudHM6IFsndHJhbnNmb3JtZWRYJywgJ3RyYW5zZm9ybWVkWSddXG4gICwgbWV0aG9kOiAndHJhbnNsYXRlJ1xuICB9XG5cbn0pXG5cblNWRy5Sb3RhdGUgPSBTVkcuaW52ZW50KHtcblxuICBwYXJlbnQ6IFNWRy5NYXRyaXhcbiwgaW5oZXJpdDogU1ZHLlRyYW5zZm9ybWF0aW9uXG5cbiwgY3JlYXRlOiBmdW5jdGlvbihzb3VyY2UsIGludmVyc2VkKXtcbiAgICBpZih0eXBlb2Ygc291cmNlID09ICdvYmplY3QnKSB0aGlzLmNvbnN0cnVjdG9yLmNhbGwodGhpcywgc291cmNlLCBpbnZlcnNlZClcbiAgICBlbHNlIHRoaXMuY29uc3RydWN0b3IuY2FsbCh0aGlzLCBbXS5zbGljZS5jYWxsKGFyZ3VtZW50cykpXG4gIH1cblxuLCBleHRlbmQ6IHtcbiAgICBhcmd1bWVudHM6IFsncm90YXRpb24nLCAnY3gnLCAnY3knXVxuICAsIG1ldGhvZDogJ3JvdGF0ZSdcbiAgLCBhdDogZnVuY3Rpb24ocG9zKXtcbiAgICAgIHZhciBtID0gbmV3IFNWRy5NYXRyaXgoKS5yb3RhdGUobmV3IFNWRy5OdW1iZXIoKS5tb3JwaCh0aGlzLnJvdGF0aW9uIC0gKHRoaXMuX3VuZG8gPyB0aGlzLl91bmRvLnJvdGF0aW9uIDogMCkpLmF0KHBvcyksIHRoaXMuY3gsIHRoaXMuY3kpXG4gICAgICByZXR1cm4gdGhpcy5pbnZlcnNlZCA/IG0uaW52ZXJzZSgpIDogbVxuICAgIH1cbiAgLCB1bmRvOiBmdW5jdGlvbihvKXtcbiAgICAgIHRoaXMuX3VuZG8gPSBvXG4gICAgfVxuICB9XG5cbn0pXG5cblNWRy5TY2FsZSA9IFNWRy5pbnZlbnQoe1xuXG4gIHBhcmVudDogU1ZHLk1hdHJpeFxuLCBpbmhlcml0OiBTVkcuVHJhbnNmb3JtYXRpb25cblxuLCBjcmVhdGU6IGZ1bmN0aW9uKHNvdXJjZSwgaW52ZXJzZWQpe1xuICAgIGlmKHR5cGVvZiBzb3VyY2UgPT0gJ29iamVjdCcpIHRoaXMuY29uc3RydWN0b3IuY2FsbCh0aGlzLCBzb3VyY2UsIGludmVyc2VkKVxuICAgIGVsc2UgdGhpcy5jb25zdHJ1Y3Rvci5jYWxsKHRoaXMsIFtdLnNsaWNlLmNhbGwoYXJndW1lbnRzKSlcbiAgfVxuXG4sIGV4dGVuZDoge1xuICAgIGFyZ3VtZW50czogWydzY2FsZVgnLCAnc2NhbGVZJywgJ2N4JywgJ2N5J11cbiAgLCBtZXRob2Q6ICdzY2FsZSdcbiAgfVxuXG59KVxuXG5TVkcuU2tldyA9IFNWRy5pbnZlbnQoe1xuXG4gIHBhcmVudDogU1ZHLk1hdHJpeFxuLCBpbmhlcml0OiBTVkcuVHJhbnNmb3JtYXRpb25cblxuLCBjcmVhdGU6IGZ1bmN0aW9uKHNvdXJjZSwgaW52ZXJzZWQpe1xuICAgIGlmKHR5cGVvZiBzb3VyY2UgPT0gJ29iamVjdCcpIHRoaXMuY29uc3RydWN0b3IuY2FsbCh0aGlzLCBzb3VyY2UsIGludmVyc2VkKVxuICAgIGVsc2UgdGhpcy5jb25zdHJ1Y3Rvci5jYWxsKHRoaXMsIFtdLnNsaWNlLmNhbGwoYXJndW1lbnRzKSlcbiAgfVxuXG4sIGV4dGVuZDoge1xuICAgIGFyZ3VtZW50czogWydza2V3WCcsICdza2V3WScsICdjeCcsICdjeSddXG4gICwgbWV0aG9kOiAnc2tldydcbiAgfVxuXG59KVxuXG5TVkcuZXh0ZW5kKFNWRy5FbGVtZW50LCB7XG4gIC8vIER5bmFtaWMgc3R5bGUgZ2VuZXJhdG9yXG4gIHN0eWxlOiBmdW5jdGlvbihzLCB2KSB7XG4gICAgaWYgKGFyZ3VtZW50cy5sZW5ndGggPT0gMCkge1xuICAgICAgLy8gZ2V0IGZ1bGwgc3R5bGVcbiAgICAgIHJldHVybiB0aGlzLm5vZGUuc3R5bGUuY3NzVGV4dCB8fCAnJ1xuXG4gICAgfSBlbHNlIGlmIChhcmd1bWVudHMubGVuZ3RoIDwgMikge1xuICAgICAgLy8gYXBwbHkgZXZlcnkgc3R5bGUgaW5kaXZpZHVhbGx5IGlmIGFuIG9iamVjdCBpcyBwYXNzZWRcbiAgICAgIGlmICh0eXBlb2YgcyA9PSAnb2JqZWN0Jykge1xuICAgICAgICBmb3IgKHYgaW4gcykgdGhpcy5zdHlsZSh2LCBzW3ZdKVxuXG4gICAgICB9IGVsc2UgaWYgKFNWRy5yZWdleC5pc0Nzcy50ZXN0KHMpKSB7XG4gICAgICAgIC8vIHBhcnNlIGNzcyBzdHJpbmdcbiAgICAgICAgcyA9IHMuc3BsaXQoJzsnKVxuXG4gICAgICAgIC8vIGFwcGx5IGV2ZXJ5IGRlZmluaXRpb24gaW5kaXZpZHVhbGx5XG4gICAgICAgIGZvciAodmFyIGkgPSAwOyBpIDwgcy5sZW5ndGg7IGkrKykge1xuICAgICAgICAgIHYgPSBzW2ldLnNwbGl0KCc6JylcbiAgICAgICAgICB0aGlzLnN0eWxlKHZbMF0ucmVwbGFjZSgvXFxzKy9nLCAnJyksIHZbMV0pXG4gICAgICAgIH1cbiAgICAgIH0gZWxzZSB7XG4gICAgICAgIC8vIGFjdCBhcyBhIGdldHRlciBpZiB0aGUgZmlyc3QgYW5kIG9ubHkgYXJndW1lbnQgaXMgbm90IGFuIG9iamVjdFxuICAgICAgICByZXR1cm4gdGhpcy5ub2RlLnN0eWxlW2NhbWVsQ2FzZShzKV1cbiAgICAgIH1cblxuICAgIH0gZWxzZSB7XG4gICAgICB0aGlzLm5vZGUuc3R5bGVbY2FtZWxDYXNlKHMpXSA9IHYgPT09IG51bGwgfHwgU1ZHLnJlZ2V4LmlzQmxhbmsudGVzdCh2KSA/ICcnIDogdlxuICAgIH1cblxuICAgIHJldHVybiB0aGlzXG4gIH1cbn0pXG5TVkcuUGFyZW50ID0gU1ZHLmludmVudCh7XG4gIC8vIEluaXRpYWxpemUgbm9kZVxuICBjcmVhdGU6IGZ1bmN0aW9uKGVsZW1lbnQpIHtcbiAgICB0aGlzLmNvbnN0cnVjdG9yLmNhbGwodGhpcywgZWxlbWVudClcbiAgfVxuXG4gIC8vIEluaGVyaXQgZnJvbVxuLCBpbmhlcml0OiBTVkcuRWxlbWVudFxuXG4gIC8vIEFkZCBjbGFzcyBtZXRob2RzXG4sIGV4dGVuZDoge1xuICAgIC8vIFJldHVybnMgYWxsIGNoaWxkIGVsZW1lbnRzXG4gICAgY2hpbGRyZW46IGZ1bmN0aW9uKCkge1xuICAgICAgcmV0dXJuIFNWRy51dGlscy5tYXAoU1ZHLnV0aWxzLmZpbHRlclNWR0VsZW1lbnRzKHRoaXMubm9kZS5jaGlsZE5vZGVzKSwgZnVuY3Rpb24obm9kZSkge1xuICAgICAgICByZXR1cm4gU1ZHLmFkb3B0KG5vZGUpXG4gICAgICB9KVxuICAgIH1cbiAgICAvLyBBZGQgZ2l2ZW4gZWxlbWVudCBhdCBhIHBvc2l0aW9uXG4gICwgYWRkOiBmdW5jdGlvbihlbGVtZW50LCBpKSB7XG4gICAgICBpZiAoaSA9PSBudWxsKVxuICAgICAgICB0aGlzLm5vZGUuYXBwZW5kQ2hpbGQoZWxlbWVudC5ub2RlKVxuICAgICAgZWxzZSBpZiAoZWxlbWVudC5ub2RlICE9IHRoaXMubm9kZS5jaGlsZE5vZGVzW2ldKVxuICAgICAgICB0aGlzLm5vZGUuaW5zZXJ0QmVmb3JlKGVsZW1lbnQubm9kZSwgdGhpcy5ub2RlLmNoaWxkTm9kZXNbaV0pXG5cbiAgICAgIHJldHVybiB0aGlzXG4gICAgfVxuICAgIC8vIEJhc2ljYWxseSBkb2VzIHRoZSBzYW1lIGFzIGBhZGQoKWAgYnV0IHJldHVybnMgdGhlIGFkZGVkIGVsZW1lbnQgaW5zdGVhZFxuICAsIHB1dDogZnVuY3Rpb24oZWxlbWVudCwgaSkge1xuICAgICAgdGhpcy5hZGQoZWxlbWVudCwgaSlcbiAgICAgIHJldHVybiBlbGVtZW50XG4gICAgfVxuICAgIC8vIENoZWNrcyBpZiB0aGUgZ2l2ZW4gZWxlbWVudCBpcyBhIGNoaWxkXG4gICwgaGFzOiBmdW5jdGlvbihlbGVtZW50KSB7XG4gICAgICByZXR1cm4gdGhpcy5pbmRleChlbGVtZW50KSA+PSAwXG4gICAgfVxuICAgIC8vIEdldHMgaW5kZXggb2YgZ2l2ZW4gZWxlbWVudFxuICAsIGluZGV4OiBmdW5jdGlvbihlbGVtZW50KSB7XG4gICAgICByZXR1cm4gW10uc2xpY2UuY2FsbCh0aGlzLm5vZGUuY2hpbGROb2RlcykuaW5kZXhPZihlbGVtZW50Lm5vZGUpXG4gICAgfVxuICAgIC8vIEdldCBhIGVsZW1lbnQgYXQgdGhlIGdpdmVuIGluZGV4XG4gICwgZ2V0OiBmdW5jdGlvbihpKSB7XG4gICAgICByZXR1cm4gU1ZHLmFkb3B0KHRoaXMubm9kZS5jaGlsZE5vZGVzW2ldKVxuICAgIH1cbiAgICAvLyBHZXQgZmlyc3QgY2hpbGRcbiAgLCBmaXJzdDogZnVuY3Rpb24oKSB7XG4gICAgICByZXR1cm4gdGhpcy5nZXQoMClcbiAgICB9XG4gICAgLy8gR2V0IHRoZSBsYXN0IGNoaWxkXG4gICwgbGFzdDogZnVuY3Rpb24oKSB7XG4gICAgICByZXR1cm4gdGhpcy5nZXQodGhpcy5ub2RlLmNoaWxkTm9kZXMubGVuZ3RoIC0gMSlcbiAgICB9XG4gICAgLy8gSXRlcmF0ZXMgb3ZlciBhbGwgY2hpbGRyZW4gYW5kIGludm9rZXMgYSBnaXZlbiBibG9ja1xuICAsIGVhY2g6IGZ1bmN0aW9uKGJsb2NrLCBkZWVwKSB7XG4gICAgICB2YXIgaSwgaWxcbiAgICAgICAgLCBjaGlsZHJlbiA9IHRoaXMuY2hpbGRyZW4oKVxuXG4gICAgICBmb3IgKGkgPSAwLCBpbCA9IGNoaWxkcmVuLmxlbmd0aDsgaSA8IGlsOyBpKyspIHtcbiAgICAgICAgaWYgKGNoaWxkcmVuW2ldIGluc3RhbmNlb2YgU1ZHLkVsZW1lbnQpXG4gICAgICAgICAgYmxvY2suYXBwbHkoY2hpbGRyZW5baV0sIFtpLCBjaGlsZHJlbl0pXG5cbiAgICAgICAgaWYgKGRlZXAgJiYgKGNoaWxkcmVuW2ldIGluc3RhbmNlb2YgU1ZHLkNvbnRhaW5lcikpXG4gICAgICAgICAgY2hpbGRyZW5baV0uZWFjaChibG9jaywgZGVlcClcbiAgICAgIH1cblxuICAgICAgcmV0dXJuIHRoaXNcbiAgICB9XG4gICAgLy8gUmVtb3ZlIGEgZ2l2ZW4gY2hpbGRcbiAgLCByZW1vdmVFbGVtZW50OiBmdW5jdGlvbihlbGVtZW50KSB7XG4gICAgICB0aGlzLm5vZGUucmVtb3ZlQ2hpbGQoZWxlbWVudC5ub2RlKVxuXG4gICAgICByZXR1cm4gdGhpc1xuICAgIH1cbiAgICAvLyBSZW1vdmUgYWxsIGVsZW1lbnRzIGluIHRoaXMgY29udGFpbmVyXG4gICwgY2xlYXI6IGZ1bmN0aW9uKCkge1xuICAgICAgLy8gcmVtb3ZlIGNoaWxkcmVuXG4gICAgICB3aGlsZSh0aGlzLm5vZGUuaGFzQ2hpbGROb2RlcygpKVxuICAgICAgICB0aGlzLm5vZGUucmVtb3ZlQ2hpbGQodGhpcy5ub2RlLmxhc3RDaGlsZClcblxuICAgICAgLy8gcmVtb3ZlIGRlZnMgcmVmZXJlbmNlXG4gICAgICBkZWxldGUgdGhpcy5fZGVmc1xuXG4gICAgICByZXR1cm4gdGhpc1xuICAgIH1cbiAgLCAvLyBHZXQgZGVmc1xuICAgIGRlZnM6IGZ1bmN0aW9uKCkge1xuICAgICAgcmV0dXJuIHRoaXMuZG9jKCkuZGVmcygpXG4gICAgfVxuICB9XG5cbn0pXG5cblNWRy5leHRlbmQoU1ZHLlBhcmVudCwge1xuXG4gIHVuZ3JvdXA6IGZ1bmN0aW9uKHBhcmVudCwgZGVwdGgpIHtcbiAgICBpZihkZXB0aCA9PT0gMCB8fCB0aGlzIGluc3RhbmNlb2YgU1ZHLkRlZnMpIHJldHVybiB0aGlzXG5cbiAgICBwYXJlbnQgPSBwYXJlbnQgfHwgKHRoaXMgaW5zdGFuY2VvZiBTVkcuRG9jID8gdGhpcyA6IHRoaXMucGFyZW50KFNWRy5QYXJlbnQpKVxuICAgIGRlcHRoID0gZGVwdGggfHwgSW5maW5pdHlcblxuICAgIHRoaXMuZWFjaChmdW5jdGlvbigpe1xuICAgICAgaWYodGhpcyBpbnN0YW5jZW9mIFNWRy5EZWZzKSByZXR1cm4gdGhpc1xuICAgICAgaWYodGhpcyBpbnN0YW5jZW9mIFNWRy5QYXJlbnQpIHJldHVybiB0aGlzLnVuZ3JvdXAocGFyZW50LCBkZXB0aC0xKVxuICAgICAgcmV0dXJuIHRoaXMudG9QYXJlbnQocGFyZW50KVxuICAgIH0pXG5cbiAgICB0aGlzLm5vZGUuZmlyc3RDaGlsZCB8fCB0aGlzLnJlbW92ZSgpXG5cbiAgICByZXR1cm4gdGhpc1xuICB9LFxuXG4gIGZsYXR0ZW46IGZ1bmN0aW9uKHBhcmVudCwgZGVwdGgpIHtcbiAgICByZXR1cm4gdGhpcy51bmdyb3VwKHBhcmVudCwgZGVwdGgpXG4gIH1cblxufSlcblNWRy5Db250YWluZXIgPSBTVkcuaW52ZW50KHtcbiAgLy8gSW5pdGlhbGl6ZSBub2RlXG4gIGNyZWF0ZTogZnVuY3Rpb24oZWxlbWVudCkge1xuICAgIHRoaXMuY29uc3RydWN0b3IuY2FsbCh0aGlzLCBlbGVtZW50KVxuICB9XG5cbiAgLy8gSW5oZXJpdCBmcm9tXG4sIGluaGVyaXQ6IFNWRy5QYXJlbnRcblxufSlcblxuU1ZHLlZpZXdCb3ggPSBTVkcuaW52ZW50KHtcblxuICBjcmVhdGU6IGZ1bmN0aW9uKHNvdXJjZSkge1xuICAgIHZhciBpLCBiYXNlID0gWzAsIDAsIDAsIDBdXG5cbiAgICB2YXIgeCwgeSwgd2lkdGgsIGhlaWdodCwgYm94LCB2aWV3LCB3ZSwgaGVcbiAgICAgICwgd20gICA9IDEgLy8gd2lkdGggbXVsdGlwbGllclxuICAgICAgLCBobSAgID0gMSAvLyBoZWlnaHQgbXVsdGlwbGllclxuICAgICAgLCByZWcgID0gL1srLV0/KD86XFxkKyg/OlxcLlxcZCopP3xcXC5cXGQrKSg/OmVbKy1dP1xcZCspPy9naVxuXG4gICAgaWYoc291cmNlIGluc3RhbmNlb2YgU1ZHLkVsZW1lbnQpe1xuXG4gICAgICB3ZSA9IHNvdXJjZVxuICAgICAgaGUgPSBzb3VyY2VcbiAgICAgIHZpZXcgPSAoc291cmNlLmF0dHIoJ3ZpZXdCb3gnKSB8fCAnJykubWF0Y2gocmVnKVxuICAgICAgYm94ID0gc291cmNlLmJib3hcblxuICAgICAgLy8gZ2V0IGRpbWVuc2lvbnMgb2YgY3VycmVudCBub2RlXG4gICAgICB3aWR0aCAgPSBuZXcgU1ZHLk51bWJlcihzb3VyY2Uud2lkdGgoKSlcbiAgICAgIGhlaWdodCA9IG5ldyBTVkcuTnVtYmVyKHNvdXJjZS5oZWlnaHQoKSlcblxuICAgICAgLy8gZmluZCBuZWFyZXN0IG5vbi1wZXJjZW50dWFsIGRpbWVuc2lvbnNcbiAgICAgIHdoaWxlICh3aWR0aC51bml0ID09ICclJykge1xuICAgICAgICB3bSAqPSB3aWR0aC52YWx1ZVxuICAgICAgICB3aWR0aCA9IG5ldyBTVkcuTnVtYmVyKHdlIGluc3RhbmNlb2YgU1ZHLkRvYyA/IHdlLnBhcmVudCgpLm9mZnNldFdpZHRoIDogd2UucGFyZW50KCkud2lkdGgoKSlcbiAgICAgICAgd2UgPSB3ZS5wYXJlbnQoKVxuICAgICAgfVxuICAgICAgd2hpbGUgKGhlaWdodC51bml0ID09ICclJykge1xuICAgICAgICBobSAqPSBoZWlnaHQudmFsdWVcbiAgICAgICAgaGVpZ2h0ID0gbmV3IFNWRy5OdW1iZXIoaGUgaW5zdGFuY2VvZiBTVkcuRG9jID8gaGUucGFyZW50KCkub2Zmc2V0SGVpZ2h0IDogaGUucGFyZW50KCkuaGVpZ2h0KCkpXG4gICAgICAgIGhlID0gaGUucGFyZW50KClcbiAgICAgIH1cblxuICAgICAgLy8gZW5zdXJlIGRlZmF1bHRzXG4gICAgICB0aGlzLnggICAgICA9IDBcbiAgICAgIHRoaXMueSAgICAgID0gMFxuICAgICAgdGhpcy53aWR0aCAgPSB3aWR0aCAgKiB3bVxuICAgICAgdGhpcy5oZWlnaHQgPSBoZWlnaHQgKiBobVxuICAgICAgdGhpcy56b29tICAgPSAxXG5cbiAgICAgIGlmICh2aWV3KSB7XG4gICAgICAgIC8vIGdldCB3aWR0aCBhbmQgaGVpZ2h0IGZyb20gdmlld2JveFxuICAgICAgICB4ICAgICAgPSBwYXJzZUZsb2F0KHZpZXdbMF0pXG4gICAgICAgIHkgICAgICA9IHBhcnNlRmxvYXQodmlld1sxXSlcbiAgICAgICAgd2lkdGggID0gcGFyc2VGbG9hdCh2aWV3WzJdKVxuICAgICAgICBoZWlnaHQgPSBwYXJzZUZsb2F0KHZpZXdbM10pXG5cbiAgICAgICAgLy8gY2FsY3VsYXRlIHpvb20gYWNjb3JpbmcgdG8gdmlld2JveFxuICAgICAgICB0aGlzLnpvb20gPSAoKHRoaXMud2lkdGggLyB0aGlzLmhlaWdodCkgPiAod2lkdGggLyBoZWlnaHQpKSA/XG4gICAgICAgICAgdGhpcy5oZWlnaHQgLyBoZWlnaHQgOlxuICAgICAgICAgIHRoaXMud2lkdGggIC8gd2lkdGhcblxuICAgICAgICAvLyBjYWxjdWxhdGUgcmVhbCBwaXhlbCBkaW1lbnNpb25zIG9uIHBhcmVudCBTVkcuRG9jIGVsZW1lbnRcbiAgICAgICAgdGhpcy54ICAgICAgPSB4XG4gICAgICAgIHRoaXMueSAgICAgID0geVxuICAgICAgICB0aGlzLndpZHRoICA9IHdpZHRoXG4gICAgICAgIHRoaXMuaGVpZ2h0ID0gaGVpZ2h0XG5cbiAgICAgIH1cblxuICAgIH1lbHNle1xuXG4gICAgICAvLyBlbnN1cmUgc291cmNlIGFzIG9iamVjdFxuICAgICAgc291cmNlID0gdHlwZW9mIHNvdXJjZSA9PT0gJ3N0cmluZycgP1xuICAgICAgICBzb3VyY2UubWF0Y2gocmVnKS5tYXAoZnVuY3Rpb24oZWwpeyByZXR1cm4gcGFyc2VGbG9hdChlbCkgfSkgOlxuICAgICAgQXJyYXkuaXNBcnJheShzb3VyY2UpID9cbiAgICAgICAgc291cmNlIDpcbiAgICAgIHR5cGVvZiBzb3VyY2UgPT0gJ29iamVjdCcgP1xuICAgICAgICBbc291cmNlLngsIHNvdXJjZS55LCBzb3VyY2Uud2lkdGgsIHNvdXJjZS5oZWlnaHRdIDpcbiAgICAgIGFyZ3VtZW50cy5sZW5ndGggPT0gNCA/XG4gICAgICAgIFtdLnNsaWNlLmNhbGwoYXJndW1lbnRzKSA6XG4gICAgICAgIGJhc2VcblxuICAgICAgdGhpcy54ID0gc291cmNlWzBdXG4gICAgICB0aGlzLnkgPSBzb3VyY2VbMV1cbiAgICAgIHRoaXMud2lkdGggPSBzb3VyY2VbMl1cbiAgICAgIHRoaXMuaGVpZ2h0ID0gc291cmNlWzNdXG4gICAgfVxuXG5cbiAgfVxuXG4sIGV4dGVuZDoge1xuXG4gICAgdG9TdHJpbmc6IGZ1bmN0aW9uKCkge1xuICAgICAgcmV0dXJuIHRoaXMueCArICcgJyArIHRoaXMueSArICcgJyArIHRoaXMud2lkdGggKyAnICcgKyB0aGlzLmhlaWdodFxuICAgIH1cbiAgLCBtb3JwaDogZnVuY3Rpb24odil7XG5cbiAgICAgIHZhciB2ID0gYXJndW1lbnRzLmxlbmd0aCA9PSAxID9cbiAgICAgICAgW3YueCwgdi55LCB2LndpZHRoLCB2LmhlaWdodF0gOlxuICAgICAgICBbXS5zbGljZS5jYWxsKGFyZ3VtZW50cylcblxuICAgICAgdGhpcy5kZXN0aW5hdGlvbiA9IG5ldyBTVkcuVmlld0JveCh2KVxuXG4gICAgICByZXR1cm4gdGhpc1xuXG4gICAgfVxuXG4gICwgYXQ6IGZ1bmN0aW9uKHBvcykge1xuXG4gICAgaWYoIXRoaXMuZGVzdGluYXRpb24pIHJldHVybiB0aGlzXG5cbiAgICByZXR1cm4gbmV3IFNWRy5WaWV3Qm94KFtcbiAgICAgICAgdGhpcy54ICsgKHRoaXMuZGVzdGluYXRpb24ueCAtIHRoaXMueCkgKiBwb3NcbiAgICAgICwgdGhpcy55ICsgKHRoaXMuZGVzdGluYXRpb24ueSAtIHRoaXMueSkgKiBwb3NcbiAgICAgICwgdGhpcy53aWR0aCArICh0aGlzLmRlc3RpbmF0aW9uLndpZHRoIC0gdGhpcy53aWR0aCkgKiBwb3NcbiAgICAgICwgdGhpcy5oZWlnaHQgKyAodGhpcy5kZXN0aW5hdGlvbi5oZWlnaHQgLSB0aGlzLmhlaWdodCkgKiBwb3NcbiAgICBdKVxuXG4gICAgfVxuXG4gIH1cblxuICAvLyBEZWZpbmUgcGFyZW50XG4sIHBhcmVudDogU1ZHLkNvbnRhaW5lclxuXG4gIC8vIEFkZCBwYXJlbnQgbWV0aG9kXG4sIGNvbnN0cnVjdDoge1xuXG4gICAgLy8gZ2V0L3NldCB2aWV3Ym94XG4gICAgdmlld2JveDogZnVuY3Rpb24odikge1xuICAgICAgaWYgKGFyZ3VtZW50cy5sZW5ndGggPT0gMClcbiAgICAgICAgLy8gYWN0IGFzIGEgZ2V0dGVyIGlmIHRoZXJlIGFyZSBubyBhcmd1bWVudHNcbiAgICAgICAgcmV0dXJuIG5ldyBTVkcuVmlld0JveCh0aGlzKVxuXG4gICAgICAvLyBvdGhlcndpc2UgYWN0IGFzIGEgc2V0dGVyXG4gICAgICB2ID0gYXJndW1lbnRzLmxlbmd0aCA9PSAxID9cbiAgICAgICAgW3YueCwgdi55LCB2LndpZHRoLCB2LmhlaWdodF0gOlxuICAgICAgICBbXS5zbGljZS5jYWxsKGFyZ3VtZW50cylcblxuICAgICAgcmV0dXJuIHRoaXMuYXR0cigndmlld0JveCcsIHYpXG4gICAgfVxuXG4gIH1cblxufSlcbi8vIEFkZCBldmVudHMgdG8gZWxlbWVudHNcbjtbICAnY2xpY2snXG4gICwgJ2RibGNsaWNrJ1xuICAsICdtb3VzZWRvd24nXG4gICwgJ21vdXNldXAnXG4gICwgJ21vdXNlb3ZlcidcbiAgLCAnbW91c2VvdXQnXG4gICwgJ21vdXNlbW92ZSdcbiAgLy8gLCAnbW91c2VlbnRlcicgLT4gbm90IHN1cHBvcnRlZCBieSBJRVxuICAvLyAsICdtb3VzZWxlYXZlJyAtPiBub3Qgc3VwcG9ydGVkIGJ5IElFXG4gICwgJ3RvdWNoc3RhcnQnXG4gICwgJ3RvdWNobW92ZSdcbiAgLCAndG91Y2hsZWF2ZSdcbiAgLCAndG91Y2hlbmQnXG4gICwgJ3RvdWNoY2FuY2VsJyBdLmZvckVhY2goZnVuY3Rpb24oZXZlbnQpIHtcblxuICAvLyBhZGQgZXZlbnQgdG8gU1ZHLkVsZW1lbnRcbiAgU1ZHLkVsZW1lbnQucHJvdG90eXBlW2V2ZW50XSA9IGZ1bmN0aW9uKGYpIHtcbiAgICB2YXIgc2VsZiA9IHRoaXNcblxuICAgIC8vIGJpbmQgZXZlbnQgdG8gZWxlbWVudCByYXRoZXIgdGhhbiBlbGVtZW50IG5vZGVcbiAgICB0aGlzLm5vZGVbJ29uJyArIGV2ZW50XSA9IHR5cGVvZiBmID09ICdmdW5jdGlvbicgP1xuICAgICAgZnVuY3Rpb24oKSB7IHJldHVybiBmLmFwcGx5KHNlbGYsIGFyZ3VtZW50cykgfSA6IG51bGxcblxuICAgIHJldHVybiB0aGlzXG4gIH1cblxufSlcblxuLy8gSW5pdGlhbGl6ZSBsaXN0ZW5lcnMgc3RhY2tcblNWRy5saXN0ZW5lcnMgPSBbXVxuU1ZHLmhhbmRsZXJNYXAgPSBbXVxuU1ZHLmxpc3RlbmVySWQgPSAwXG5cbi8vIEFkZCBldmVudCBiaW5kZXIgaW4gdGhlIFNWRyBuYW1lc3BhY2VcblNWRy5vbiA9IGZ1bmN0aW9uKG5vZGUsIGV2ZW50LCBsaXN0ZW5lciwgYmluZGluZykge1xuICAvLyBjcmVhdGUgbGlzdGVuZXIsIGdldCBvYmplY3QtaW5kZXhcbiAgdmFyIGwgICAgID0gbGlzdGVuZXIuYmluZChiaW5kaW5nIHx8IG5vZGUuaW5zdGFuY2UgfHwgbm9kZSlcbiAgICAsIGluZGV4ID0gKFNWRy5oYW5kbGVyTWFwLmluZGV4T2Yobm9kZSkgKyAxIHx8IFNWRy5oYW5kbGVyTWFwLnB1c2gobm9kZSkpIC0gMVxuICAgICwgZXYgICAgPSBldmVudC5zcGxpdCgnLicpWzBdXG4gICAgLCBucyAgICA9IGV2ZW50LnNwbGl0KCcuJylbMV0gfHwgJyonXG5cblxuICAvLyBlbnN1cmUgdmFsaWQgb2JqZWN0XG4gIFNWRy5saXN0ZW5lcnNbaW5kZXhdICAgICAgICAgPSBTVkcubGlzdGVuZXJzW2luZGV4XSAgICAgICAgIHx8IHt9XG4gIFNWRy5saXN0ZW5lcnNbaW5kZXhdW2V2XSAgICAgPSBTVkcubGlzdGVuZXJzW2luZGV4XVtldl0gICAgIHx8IHt9XG4gIFNWRy5saXN0ZW5lcnNbaW5kZXhdW2V2XVtuc10gPSBTVkcubGlzdGVuZXJzW2luZGV4XVtldl1bbnNdIHx8IHt9XG5cbiAgaWYoIWxpc3RlbmVyLl9zdmdqc0xpc3RlbmVySWQpXG4gICAgbGlzdGVuZXIuX3N2Z2pzTGlzdGVuZXJJZCA9ICsrU1ZHLmxpc3RlbmVySWRcblxuICAvLyByZWZlcmVuY2UgbGlzdGVuZXJcbiAgU1ZHLmxpc3RlbmVyc1tpbmRleF1bZXZdW25zXVtsaXN0ZW5lci5fc3ZnanNMaXN0ZW5lcklkXSA9IGxcblxuICAvLyBhZGQgbGlzdGVuZXJcbiAgbm9kZS5hZGRFdmVudExpc3RlbmVyKGV2LCBsLCBmYWxzZSlcbn1cblxuLy8gQWRkIGV2ZW50IHVuYmluZGVyIGluIHRoZSBTVkcgbmFtZXNwYWNlXG5TVkcub2ZmID0gZnVuY3Rpb24obm9kZSwgZXZlbnQsIGxpc3RlbmVyKSB7XG4gIHZhciBpbmRleCA9IFNWRy5oYW5kbGVyTWFwLmluZGV4T2Yobm9kZSlcbiAgICAsIGV2ICAgID0gZXZlbnQgJiYgZXZlbnQuc3BsaXQoJy4nKVswXVxuICAgICwgbnMgICAgPSBldmVudCAmJiBldmVudC5zcGxpdCgnLicpWzFdXG5cbiAgaWYoaW5kZXggPT0gLTEpIHJldHVyblxuXG4gIGlmIChsaXN0ZW5lcikge1xuICAgIGlmKHR5cGVvZiBsaXN0ZW5lciA9PSAnZnVuY3Rpb24nKSBsaXN0ZW5lciA9IGxpc3RlbmVyLl9zdmdqc0xpc3RlbmVySWRcbiAgICBpZighbGlzdGVuZXIpIHJldHVyblxuXG4gICAgLy8gcmVtb3ZlIGxpc3RlbmVyIHJlZmVyZW5jZVxuICAgIGlmIChTVkcubGlzdGVuZXJzW2luZGV4XVtldl0gJiYgU1ZHLmxpc3RlbmVyc1tpbmRleF1bZXZdW25zIHx8ICcqJ10pIHtcbiAgICAgIC8vIHJlbW92ZSBsaXN0ZW5lclxuICAgICAgbm9kZS5yZW1vdmVFdmVudExpc3RlbmVyKGV2LCBTVkcubGlzdGVuZXJzW2luZGV4XVtldl1bbnMgfHwgJyonXVtsaXN0ZW5lcl0sIGZhbHNlKVxuXG4gICAgICBkZWxldGUgU1ZHLmxpc3RlbmVyc1tpbmRleF1bZXZdW25zIHx8ICcqJ11bbGlzdGVuZXJdXG4gICAgfVxuXG4gIH0gZWxzZSBpZiAobnMgJiYgZXYpIHtcbiAgICAvLyByZW1vdmUgYWxsIGxpc3RlbmVycyBmb3IgYSBuYW1lc3BhY2VkIGV2ZW50XG4gICAgaWYgKFNWRy5saXN0ZW5lcnNbaW5kZXhdW2V2XSAmJiBTVkcubGlzdGVuZXJzW2luZGV4XVtldl1bbnNdKSB7XG4gICAgICBmb3IgKGxpc3RlbmVyIGluIFNWRy5saXN0ZW5lcnNbaW5kZXhdW2V2XVtuc10pXG4gICAgICAgIFNWRy5vZmYobm9kZSwgW2V2LCBuc10uam9pbignLicpLCBsaXN0ZW5lcilcblxuICAgICAgZGVsZXRlIFNWRy5saXN0ZW5lcnNbaW5kZXhdW2V2XVtuc11cbiAgICB9XG5cbiAgfSBlbHNlIGlmIChucyl7XG4gICAgLy8gcmVtb3ZlIGFsbCBsaXN0ZW5lcnMgZm9yIGEgc3BlY2lmaWMgbmFtZXNwYWNlXG4gICAgZm9yKGV2ZW50IGluIFNWRy5saXN0ZW5lcnNbaW5kZXhdKXtcbiAgICAgICAgZm9yKG5hbWVzcGFjZSBpbiBTVkcubGlzdGVuZXJzW2luZGV4XVtldmVudF0pe1xuICAgICAgICAgICAgaWYobnMgPT09IG5hbWVzcGFjZSl7XG4gICAgICAgICAgICAgICAgU1ZHLm9mZihub2RlLCBbZXZlbnQsIG5zXS5qb2luKCcuJykpXG4gICAgICAgICAgICB9XG4gICAgICAgIH1cbiAgICB9XG5cbiAgfSBlbHNlIGlmIChldikge1xuICAgIC8vIHJlbW92ZSBhbGwgbGlzdGVuZXJzIGZvciB0aGUgZXZlbnRcbiAgICBpZiAoU1ZHLmxpc3RlbmVyc1tpbmRleF1bZXZdKSB7XG4gICAgICBmb3IgKG5hbWVzcGFjZSBpbiBTVkcubGlzdGVuZXJzW2luZGV4XVtldl0pXG4gICAgICAgIFNWRy5vZmYobm9kZSwgW2V2LCBuYW1lc3BhY2VdLmpvaW4oJy4nKSlcblxuICAgICAgZGVsZXRlIFNWRy5saXN0ZW5lcnNbaW5kZXhdW2V2XVxuICAgIH1cblxuICB9IGVsc2Uge1xuICAgIC8vIHJlbW92ZSBhbGwgbGlzdGVuZXJzIG9uIGEgZ2l2ZW4gbm9kZVxuICAgIGZvciAoZXZlbnQgaW4gU1ZHLmxpc3RlbmVyc1tpbmRleF0pXG4gICAgICBTVkcub2ZmKG5vZGUsIGV2ZW50KVxuXG4gICAgZGVsZXRlIFNWRy5saXN0ZW5lcnNbaW5kZXhdXG5cbiAgfVxufVxuXG4vL1xuU1ZHLmV4dGVuZChTVkcuRWxlbWVudCwge1xuICAvLyBCaW5kIGdpdmVuIGV2ZW50IHRvIGxpc3RlbmVyXG4gIG9uOiBmdW5jdGlvbihldmVudCwgbGlzdGVuZXIsIGJpbmRpbmcpIHtcbiAgICBTVkcub24odGhpcy5ub2RlLCBldmVudCwgbGlzdGVuZXIsIGJpbmRpbmcpXG5cbiAgICByZXR1cm4gdGhpc1xuICB9XG4gIC8vIFVuYmluZCBldmVudCBmcm9tIGxpc3RlbmVyXG4sIG9mZjogZnVuY3Rpb24oZXZlbnQsIGxpc3RlbmVyKSB7XG4gICAgU1ZHLm9mZih0aGlzLm5vZGUsIGV2ZW50LCBsaXN0ZW5lcilcblxuICAgIHJldHVybiB0aGlzXG4gIH1cbiAgLy8gRmlyZSBnaXZlbiBldmVudFxuLCBmaXJlOiBmdW5jdGlvbihldmVudCwgZGF0YSkge1xuXG4gICAgLy8gRGlzcGF0Y2ggZXZlbnRcbiAgICBpZihldmVudCBpbnN0YW5jZW9mIEV2ZW50KXtcbiAgICAgICAgdGhpcy5ub2RlLmRpc3BhdGNoRXZlbnQoZXZlbnQpXG4gICAgfWVsc2V7XG4gICAgICAgIHRoaXMubm9kZS5kaXNwYXRjaEV2ZW50KG5ldyBDdXN0b21FdmVudChldmVudCwge2RldGFpbDpkYXRhfSkpXG4gICAgfVxuXG4gICAgcmV0dXJuIHRoaXNcbiAgfVxufSlcblxuU1ZHLkRlZnMgPSBTVkcuaW52ZW50KHtcbiAgLy8gSW5pdGlhbGl6ZSBub2RlXG4gIGNyZWF0ZTogJ2RlZnMnXG5cbiAgLy8gSW5oZXJpdCBmcm9tXG4sIGluaGVyaXQ6IFNWRy5Db250YWluZXJcblxufSlcblNWRy5HID0gU1ZHLmludmVudCh7XG4gIC8vIEluaXRpYWxpemUgbm9kZVxuICBjcmVhdGU6ICdnJ1xuXG4gIC8vIEluaGVyaXQgZnJvbVxuLCBpbmhlcml0OiBTVkcuQ29udGFpbmVyXG5cbiAgLy8gQWRkIGNsYXNzIG1ldGhvZHNcbiwgZXh0ZW5kOiB7XG4gICAgLy8gTW92ZSBvdmVyIHgtYXhpc1xuICAgIHg6IGZ1bmN0aW9uKHgpIHtcbiAgICAgIHJldHVybiB4ID09IG51bGwgPyB0aGlzLnRyYW5zZm9ybSgneCcpIDogdGhpcy50cmFuc2Zvcm0oeyB4OiB4IC0gdGhpcy54KCkgfSwgdHJ1ZSlcbiAgICB9XG4gICAgLy8gTW92ZSBvdmVyIHktYXhpc1xuICAsIHk6IGZ1bmN0aW9uKHkpIHtcbiAgICAgIHJldHVybiB5ID09IG51bGwgPyB0aGlzLnRyYW5zZm9ybSgneScpIDogdGhpcy50cmFuc2Zvcm0oeyB5OiB5IC0gdGhpcy55KCkgfSwgdHJ1ZSlcbiAgICB9XG4gICAgLy8gTW92ZSBieSBjZW50ZXIgb3ZlciB4LWF4aXNcbiAgLCBjeDogZnVuY3Rpb24oeCkge1xuICAgICAgcmV0dXJuIHggPT0gbnVsbCA/IHRoaXMuZ2JveCgpLmN4IDogdGhpcy54KHggLSB0aGlzLmdib3goKS53aWR0aCAvIDIpXG4gICAgfVxuICAgIC8vIE1vdmUgYnkgY2VudGVyIG92ZXIgeS1heGlzXG4gICwgY3k6IGZ1bmN0aW9uKHkpIHtcbiAgICAgIHJldHVybiB5ID09IG51bGwgPyB0aGlzLmdib3goKS5jeSA6IHRoaXMueSh5IC0gdGhpcy5nYm94KCkuaGVpZ2h0IC8gMilcbiAgICB9XG4gICwgZ2JveDogZnVuY3Rpb24oKSB7XG5cbiAgICAgIHZhciBiYm94ICA9IHRoaXMuYmJveCgpXG4gICAgICAgICwgdHJhbnMgPSB0aGlzLnRyYW5zZm9ybSgpXG5cbiAgICAgIGJib3gueCAgKz0gdHJhbnMueFxuICAgICAgYmJveC54MiArPSB0cmFucy54XG4gICAgICBiYm94LmN4ICs9IHRyYW5zLnhcblxuICAgICAgYmJveC55ICArPSB0cmFucy55XG4gICAgICBiYm94LnkyICs9IHRyYW5zLnlcbiAgICAgIGJib3guY3kgKz0gdHJhbnMueVxuXG4gICAgICByZXR1cm4gYmJveFxuICAgIH1cbiAgfVxuXG4gIC8vIEFkZCBwYXJlbnQgbWV0aG9kXG4sIGNvbnN0cnVjdDoge1xuICAgIC8vIENyZWF0ZSBhIGdyb3VwIGVsZW1lbnRcbiAgICBncm91cDogZnVuY3Rpb24oKSB7XG4gICAgICByZXR1cm4gdGhpcy5wdXQobmV3IFNWRy5HKVxuICAgIH1cbiAgfVxufSlcblxuLy8gIyMjIFRoaXMgbW9kdWxlIGFkZHMgYmFja3dhcmQgLyBmb3J3YXJkIGZ1bmN0aW9uYWxpdHkgdG8gZWxlbWVudHMuXG5cbi8vXG5TVkcuZXh0ZW5kKFNWRy5FbGVtZW50LCB7XG4gIC8vIEdldCBhbGwgc2libGluZ3MsIGluY2x1ZGluZyBteXNlbGZcbiAgc2libGluZ3M6IGZ1bmN0aW9uKCkge1xuICAgIHJldHVybiB0aGlzLnBhcmVudCgpLmNoaWxkcmVuKClcbiAgfVxuICAvLyBHZXQgdGhlIGN1cmVudCBwb3NpdGlvbiBzaWJsaW5nc1xuLCBwb3NpdGlvbjogZnVuY3Rpb24oKSB7XG4gICAgcmV0dXJuIHRoaXMucGFyZW50KCkuaW5kZXgodGhpcylcbiAgfVxuICAvLyBHZXQgdGhlIG5leHQgZWxlbWVudCAod2lsbCByZXR1cm4gbnVsbCBpZiB0aGVyZSBpcyBub25lKVxuLCBuZXh0OiBmdW5jdGlvbigpIHtcbiAgICByZXR1cm4gdGhpcy5zaWJsaW5ncygpW3RoaXMucG9zaXRpb24oKSArIDFdXG4gIH1cbiAgLy8gR2V0IHRoZSBuZXh0IGVsZW1lbnQgKHdpbGwgcmV0dXJuIG51bGwgaWYgdGhlcmUgaXMgbm9uZSlcbiwgcHJldmlvdXM6IGZ1bmN0aW9uKCkge1xuICAgIHJldHVybiB0aGlzLnNpYmxpbmdzKClbdGhpcy5wb3NpdGlvbigpIC0gMV1cbiAgfVxuICAvLyBTZW5kIGdpdmVuIGVsZW1lbnQgb25lIHN0ZXAgZm9yd2FyZFxuLCBmb3J3YXJkOiBmdW5jdGlvbigpIHtcbiAgICB2YXIgaSA9IHRoaXMucG9zaXRpb24oKSArIDFcbiAgICAgICwgcCA9IHRoaXMucGFyZW50KClcblxuICAgIC8vIG1vdmUgbm9kZSBvbmUgc3RlcCBmb3J3YXJkXG4gICAgcC5yZW1vdmVFbGVtZW50KHRoaXMpLmFkZCh0aGlzLCBpKVxuXG4gICAgLy8gbWFrZSBzdXJlIGRlZnMgbm9kZSBpcyBhbHdheXMgYXQgdGhlIHRvcFxuICAgIGlmIChwIGluc3RhbmNlb2YgU1ZHLkRvYylcbiAgICAgIHAubm9kZS5hcHBlbmRDaGlsZChwLmRlZnMoKS5ub2RlKVxuXG4gICAgcmV0dXJuIHRoaXNcbiAgfVxuICAvLyBTZW5kIGdpdmVuIGVsZW1lbnQgb25lIHN0ZXAgYmFja3dhcmRcbiwgYmFja3dhcmQ6IGZ1bmN0aW9uKCkge1xuICAgIHZhciBpID0gdGhpcy5wb3NpdGlvbigpXG5cbiAgICBpZiAoaSA+IDApXG4gICAgICB0aGlzLnBhcmVudCgpLnJlbW92ZUVsZW1lbnQodGhpcykuYWRkKHRoaXMsIGkgLSAxKVxuXG4gICAgcmV0dXJuIHRoaXNcbiAgfVxuICAvLyBTZW5kIGdpdmVuIGVsZW1lbnQgYWxsIHRoZSB3YXkgdG8gdGhlIGZyb250XG4sIGZyb250OiBmdW5jdGlvbigpIHtcbiAgICB2YXIgcCA9IHRoaXMucGFyZW50KClcblxuICAgIC8vIE1vdmUgbm9kZSBmb3J3YXJkXG4gICAgcC5ub2RlLmFwcGVuZENoaWxkKHRoaXMubm9kZSlcblxuICAgIC8vIE1ha2Ugc3VyZSBkZWZzIG5vZGUgaXMgYWx3YXlzIGF0IHRoZSB0b3BcbiAgICBpZiAocCBpbnN0YW5jZW9mIFNWRy5Eb2MpXG4gICAgICBwLm5vZGUuYXBwZW5kQ2hpbGQocC5kZWZzKCkubm9kZSlcblxuICAgIHJldHVybiB0aGlzXG4gIH1cbiAgLy8gU2VuZCBnaXZlbiBlbGVtZW50IGFsbCB0aGUgd2F5IHRvIHRoZSBiYWNrXG4sIGJhY2s6IGZ1bmN0aW9uKCkge1xuICAgIGlmICh0aGlzLnBvc2l0aW9uKCkgPiAwKVxuICAgICAgdGhpcy5wYXJlbnQoKS5yZW1vdmVFbGVtZW50KHRoaXMpLmFkZCh0aGlzLCAwKVxuXG4gICAgcmV0dXJuIHRoaXNcbiAgfVxuICAvLyBJbnNlcnRzIGEgZ2l2ZW4gZWxlbWVudCBiZWZvcmUgdGhlIHRhcmdldGVkIGVsZW1lbnRcbiwgYmVmb3JlOiBmdW5jdGlvbihlbGVtZW50KSB7XG4gICAgZWxlbWVudC5yZW1vdmUoKVxuXG4gICAgdmFyIGkgPSB0aGlzLnBvc2l0aW9uKClcblxuICAgIHRoaXMucGFyZW50KCkuYWRkKGVsZW1lbnQsIGkpXG5cbiAgICByZXR1cm4gdGhpc1xuICB9XG4gIC8vIEluc3RlcnMgYSBnaXZlbiBlbGVtZW50IGFmdGVyIHRoZSB0YXJnZXRlZCBlbGVtZW50XG4sIGFmdGVyOiBmdW5jdGlvbihlbGVtZW50KSB7XG4gICAgZWxlbWVudC5yZW1vdmUoKVxuXG4gICAgdmFyIGkgPSB0aGlzLnBvc2l0aW9uKClcblxuICAgIHRoaXMucGFyZW50KCkuYWRkKGVsZW1lbnQsIGkgKyAxKVxuXG4gICAgcmV0dXJuIHRoaXNcbiAgfVxuXG59KVxuU1ZHLk1hc2sgPSBTVkcuaW52ZW50KHtcbiAgLy8gSW5pdGlhbGl6ZSBub2RlXG4gIGNyZWF0ZTogZnVuY3Rpb24oKSB7XG4gICAgdGhpcy5jb25zdHJ1Y3Rvci5jYWxsKHRoaXMsIFNWRy5jcmVhdGUoJ21hc2snKSlcblxuICAgIC8vIGtlZXAgcmVmZXJlbmNlcyB0byBtYXNrZWQgZWxlbWVudHNcbiAgICB0aGlzLnRhcmdldHMgPSBbXVxuICB9XG5cbiAgLy8gSW5oZXJpdCBmcm9tXG4sIGluaGVyaXQ6IFNWRy5Db250YWluZXJcblxuICAvLyBBZGQgY2xhc3MgbWV0aG9kc1xuLCBleHRlbmQ6IHtcbiAgICAvLyBVbm1hc2sgYWxsIG1hc2tlZCBlbGVtZW50cyBhbmQgcmVtb3ZlIGl0c2VsZlxuICAgIHJlbW92ZTogZnVuY3Rpb24oKSB7XG4gICAgICAvLyB1bm1hc2sgYWxsIHRhcmdldHNcbiAgICAgIGZvciAodmFyIGkgPSB0aGlzLnRhcmdldHMubGVuZ3RoIC0gMTsgaSA+PSAwOyBpLS0pXG4gICAgICAgIGlmICh0aGlzLnRhcmdldHNbaV0pXG4gICAgICAgICAgdGhpcy50YXJnZXRzW2ldLnVubWFzaygpXG4gICAgICB0aGlzLnRhcmdldHMgPSBbXVxuXG4gICAgICAvLyByZW1vdmUgbWFzayBmcm9tIHBhcmVudFxuICAgICAgdGhpcy5wYXJlbnQoKS5yZW1vdmVFbGVtZW50KHRoaXMpXG5cbiAgICAgIHJldHVybiB0aGlzXG4gICAgfVxuICB9XG5cbiAgLy8gQWRkIHBhcmVudCBtZXRob2RcbiwgY29uc3RydWN0OiB7XG4gICAgLy8gQ3JlYXRlIG1hc2tpbmcgZWxlbWVudFxuICAgIG1hc2s6IGZ1bmN0aW9uKCkge1xuICAgICAgcmV0dXJuIHRoaXMuZGVmcygpLnB1dChuZXcgU1ZHLk1hc2spXG4gICAgfVxuICB9XG59KVxuXG5cblNWRy5leHRlbmQoU1ZHLkVsZW1lbnQsIHtcbiAgLy8gRGlzdHJpYnV0ZSBtYXNrIHRvIHN2ZyBlbGVtZW50XG4gIG1hc2tXaXRoOiBmdW5jdGlvbihlbGVtZW50KSB7XG4gICAgLy8gdXNlIGdpdmVuIG1hc2sgb3IgY3JlYXRlIGEgbmV3IG9uZVxuICAgIHRoaXMubWFza2VyID0gZWxlbWVudCBpbnN0YW5jZW9mIFNWRy5NYXNrID8gZWxlbWVudCA6IHRoaXMucGFyZW50KCkubWFzaygpLmFkZChlbGVtZW50KVxuXG4gICAgLy8gc3RvcmUgcmV2ZXJlbmNlIG9uIHNlbGYgaW4gbWFza1xuICAgIHRoaXMubWFza2VyLnRhcmdldHMucHVzaCh0aGlzKVxuXG4gICAgLy8gYXBwbHkgbWFza1xuICAgIHJldHVybiB0aGlzLmF0dHIoJ21hc2snLCAndXJsKFwiIycgKyB0aGlzLm1hc2tlci5hdHRyKCdpZCcpICsgJ1wiKScpXG4gIH1cbiAgLy8gVW5tYXNrIGVsZW1lbnRcbiwgdW5tYXNrOiBmdW5jdGlvbigpIHtcbiAgICBkZWxldGUgdGhpcy5tYXNrZXJcbiAgICByZXR1cm4gdGhpcy5hdHRyKCdtYXNrJywgbnVsbClcbiAgfVxuXG59KVxuXG5TVkcuQ2xpcFBhdGggPSBTVkcuaW52ZW50KHtcbiAgLy8gSW5pdGlhbGl6ZSBub2RlXG4gIGNyZWF0ZTogZnVuY3Rpb24oKSB7XG4gICAgdGhpcy5jb25zdHJ1Y3Rvci5jYWxsKHRoaXMsIFNWRy5jcmVhdGUoJ2NsaXBQYXRoJykpXG5cbiAgICAvLyBrZWVwIHJlZmVyZW5jZXMgdG8gY2xpcHBlZCBlbGVtZW50c1xuICAgIHRoaXMudGFyZ2V0cyA9IFtdXG4gIH1cblxuICAvLyBJbmhlcml0IGZyb21cbiwgaW5oZXJpdDogU1ZHLkNvbnRhaW5lclxuXG4gIC8vIEFkZCBjbGFzcyBtZXRob2RzXG4sIGV4dGVuZDoge1xuICAgIC8vIFVuY2xpcCBhbGwgY2xpcHBlZCBlbGVtZW50cyBhbmQgcmVtb3ZlIGl0c2VsZlxuICAgIHJlbW92ZTogZnVuY3Rpb24oKSB7XG4gICAgICAvLyB1bmNsaXAgYWxsIHRhcmdldHNcbiAgICAgIGZvciAodmFyIGkgPSB0aGlzLnRhcmdldHMubGVuZ3RoIC0gMTsgaSA+PSAwOyBpLS0pXG4gICAgICAgIGlmICh0aGlzLnRhcmdldHNbaV0pXG4gICAgICAgICAgdGhpcy50YXJnZXRzW2ldLnVuY2xpcCgpXG4gICAgICB0aGlzLnRhcmdldHMgPSBbXVxuXG4gICAgICAvLyByZW1vdmUgY2xpcFBhdGggZnJvbSBwYXJlbnRcbiAgICAgIHRoaXMucGFyZW50KCkucmVtb3ZlRWxlbWVudCh0aGlzKVxuXG4gICAgICByZXR1cm4gdGhpc1xuICAgIH1cbiAgfVxuXG4gIC8vIEFkZCBwYXJlbnQgbWV0aG9kXG4sIGNvbnN0cnVjdDoge1xuICAgIC8vIENyZWF0ZSBjbGlwcGluZyBlbGVtZW50XG4gICAgY2xpcDogZnVuY3Rpb24oKSB7XG4gICAgICByZXR1cm4gdGhpcy5kZWZzKCkucHV0KG5ldyBTVkcuQ2xpcFBhdGgpXG4gICAgfVxuICB9XG59KVxuXG4vL1xuU1ZHLmV4dGVuZChTVkcuRWxlbWVudCwge1xuICAvLyBEaXN0cmlidXRlIGNsaXBQYXRoIHRvIHN2ZyBlbGVtZW50XG4gIGNsaXBXaXRoOiBmdW5jdGlvbihlbGVtZW50KSB7XG4gICAgLy8gdXNlIGdpdmVuIGNsaXAgb3IgY3JlYXRlIGEgbmV3IG9uZVxuICAgIHRoaXMuY2xpcHBlciA9IGVsZW1lbnQgaW5zdGFuY2VvZiBTVkcuQ2xpcFBhdGggPyBlbGVtZW50IDogdGhpcy5wYXJlbnQoKS5jbGlwKCkuYWRkKGVsZW1lbnQpXG5cbiAgICAvLyBzdG9yZSByZXZlcmVuY2Ugb24gc2VsZiBpbiBtYXNrXG4gICAgdGhpcy5jbGlwcGVyLnRhcmdldHMucHVzaCh0aGlzKVxuXG4gICAgLy8gYXBwbHkgbWFza1xuICAgIHJldHVybiB0aGlzLmF0dHIoJ2NsaXAtcGF0aCcsICd1cmwoXCIjJyArIHRoaXMuY2xpcHBlci5hdHRyKCdpZCcpICsgJ1wiKScpXG4gIH1cbiAgLy8gVW5jbGlwIGVsZW1lbnRcbiwgdW5jbGlwOiBmdW5jdGlvbigpIHtcbiAgICBkZWxldGUgdGhpcy5jbGlwcGVyXG4gICAgcmV0dXJuIHRoaXMuYXR0cignY2xpcC1wYXRoJywgbnVsbClcbiAgfVxuXG59KVxuU1ZHLkdyYWRpZW50ID0gU1ZHLmludmVudCh7XG4gIC8vIEluaXRpYWxpemUgbm9kZVxuICBjcmVhdGU6IGZ1bmN0aW9uKHR5cGUpIHtcbiAgICB0aGlzLmNvbnN0cnVjdG9yLmNhbGwodGhpcywgU1ZHLmNyZWF0ZSh0eXBlICsgJ0dyYWRpZW50JykpXG5cbiAgICAvLyBzdG9yZSB0eXBlXG4gICAgdGhpcy50eXBlID0gdHlwZVxuICB9XG5cbiAgLy8gSW5oZXJpdCBmcm9tXG4sIGluaGVyaXQ6IFNWRy5Db250YWluZXJcblxuICAvLyBBZGQgY2xhc3MgbWV0aG9kc1xuLCBleHRlbmQ6IHtcbiAgICAvLyBBZGQgYSBjb2xvciBzdG9wXG4gICAgYXQ6IGZ1bmN0aW9uKG9mZnNldCwgY29sb3IsIG9wYWNpdHkpIHtcbiAgICAgIHJldHVybiB0aGlzLnB1dChuZXcgU1ZHLlN0b3ApLnVwZGF0ZShvZmZzZXQsIGNvbG9yLCBvcGFjaXR5KVxuICAgIH1cbiAgICAvLyBVcGRhdGUgZ3JhZGllbnRcbiAgLCB1cGRhdGU6IGZ1bmN0aW9uKGJsb2NrKSB7XG4gICAgICAvLyByZW1vdmUgYWxsIHN0b3BzXG4gICAgICB0aGlzLmNsZWFyKClcblxuICAgICAgLy8gaW52b2tlIHBhc3NlZCBibG9ja1xuICAgICAgaWYgKHR5cGVvZiBibG9jayA9PSAnZnVuY3Rpb24nKVxuICAgICAgICBibG9jay5jYWxsKHRoaXMsIHRoaXMpXG5cbiAgICAgIHJldHVybiB0aGlzXG4gICAgfVxuICAgIC8vIFJldHVybiB0aGUgZmlsbCBpZFxuICAsIGZpbGw6IGZ1bmN0aW9uKCkge1xuICAgICAgcmV0dXJuICd1cmwoIycgKyB0aGlzLmlkKCkgKyAnKSdcbiAgICB9XG4gICAgLy8gQWxpYXMgc3RyaW5nIGNvbnZlcnRpb24gdG8gZmlsbFxuICAsIHRvU3RyaW5nOiBmdW5jdGlvbigpIHtcbiAgICAgIHJldHVybiB0aGlzLmZpbGwoKVxuICAgIH1cbiAgICAvLyBjdXN0b20gYXR0ciB0byBoYW5kbGUgdHJhbnNmb3JtXG4gICwgYXR0cjogZnVuY3Rpb24oYSwgYiwgYykge1xuICAgICAgaWYoYSA9PSAndHJhbnNmb3JtJykgYSA9ICdncmFkaWVudFRyYW5zZm9ybSdcbiAgICAgIHJldHVybiBTVkcuQ29udGFpbmVyLnByb3RvdHlwZS5hdHRyLmNhbGwodGhpcywgYSwgYiwgYylcbiAgICB9XG4gIH1cblxuICAvLyBBZGQgcGFyZW50IG1ldGhvZFxuLCBjb25zdHJ1Y3Q6IHtcbiAgICAvLyBDcmVhdGUgZ3JhZGllbnQgZWxlbWVudCBpbiBkZWZzXG4gICAgZ3JhZGllbnQ6IGZ1bmN0aW9uKHR5cGUsIGJsb2NrKSB7XG4gICAgICByZXR1cm4gdGhpcy5kZWZzKCkuZ3JhZGllbnQodHlwZSwgYmxvY2spXG4gICAgfVxuICB9XG59KVxuXG4vLyBBZGQgYW5pbWF0YWJsZSBtZXRob2RzIHRvIGJvdGggZ3JhZGllbnQgYW5kIGZ4IG1vZHVsZVxuU1ZHLmV4dGVuZChTVkcuR3JhZGllbnQsIFNWRy5GWCwge1xuICAvLyBGcm9tIHBvc2l0aW9uXG4gIGZyb206IGZ1bmN0aW9uKHgsIHkpIHtcbiAgICByZXR1cm4gKHRoaXMuX3RhcmdldCB8fCB0aGlzKS50eXBlID09ICdyYWRpYWwnID9cbiAgICAgIHRoaXMuYXR0cih7IGZ4OiBuZXcgU1ZHLk51bWJlcih4KSwgZnk6IG5ldyBTVkcuTnVtYmVyKHkpIH0pIDpcbiAgICAgIHRoaXMuYXR0cih7IHgxOiBuZXcgU1ZHLk51bWJlcih4KSwgeTE6IG5ldyBTVkcuTnVtYmVyKHkpIH0pXG4gIH1cbiAgLy8gVG8gcG9zaXRpb25cbiwgdG86IGZ1bmN0aW9uKHgsIHkpIHtcbiAgICByZXR1cm4gKHRoaXMuX3RhcmdldCB8fCB0aGlzKS50eXBlID09ICdyYWRpYWwnID9cbiAgICAgIHRoaXMuYXR0cih7IGN4OiBuZXcgU1ZHLk51bWJlcih4KSwgY3k6IG5ldyBTVkcuTnVtYmVyKHkpIH0pIDpcbiAgICAgIHRoaXMuYXR0cih7IHgyOiBuZXcgU1ZHLk51bWJlcih4KSwgeTI6IG5ldyBTVkcuTnVtYmVyKHkpIH0pXG4gIH1cbn0pXG5cbi8vIEJhc2UgZ3JhZGllbnQgZ2VuZXJhdGlvblxuU1ZHLmV4dGVuZChTVkcuRGVmcywge1xuICAvLyBkZWZpbmUgZ3JhZGllbnRcbiAgZ3JhZGllbnQ6IGZ1bmN0aW9uKHR5cGUsIGJsb2NrKSB7XG4gICAgcmV0dXJuIHRoaXMucHV0KG5ldyBTVkcuR3JhZGllbnQodHlwZSkpLnVwZGF0ZShibG9jaylcbiAgfVxuXG59KVxuXG5TVkcuU3RvcCA9IFNWRy5pbnZlbnQoe1xuICAvLyBJbml0aWFsaXplIG5vZGVcbiAgY3JlYXRlOiAnc3RvcCdcblxuICAvLyBJbmhlcml0IGZyb21cbiwgaW5oZXJpdDogU1ZHLkVsZW1lbnRcblxuICAvLyBBZGQgY2xhc3MgbWV0aG9kc1xuLCBleHRlbmQ6IHtcbiAgICAvLyBhZGQgY29sb3Igc3RvcHNcbiAgICB1cGRhdGU6IGZ1bmN0aW9uKG8pIHtcbiAgICAgIGlmICh0eXBlb2YgbyA9PSAnbnVtYmVyJyB8fCBvIGluc3RhbmNlb2YgU1ZHLk51bWJlcikge1xuICAgICAgICBvID0ge1xuICAgICAgICAgIG9mZnNldDogIGFyZ3VtZW50c1swXVxuICAgICAgICAsIGNvbG9yOiAgIGFyZ3VtZW50c1sxXVxuICAgICAgICAsIG9wYWNpdHk6IGFyZ3VtZW50c1syXVxuICAgICAgICB9XG4gICAgICB9XG5cbiAgICAgIC8vIHNldCBhdHRyaWJ1dGVzXG4gICAgICBpZiAoby5vcGFjaXR5ICE9IG51bGwpIHRoaXMuYXR0cignc3RvcC1vcGFjaXR5Jywgby5vcGFjaXR5KVxuICAgICAgaWYgKG8uY29sb3IgICAhPSBudWxsKSB0aGlzLmF0dHIoJ3N0b3AtY29sb3InLCBvLmNvbG9yKVxuICAgICAgaWYgKG8ub2Zmc2V0ICAhPSBudWxsKSB0aGlzLmF0dHIoJ29mZnNldCcsIG5ldyBTVkcuTnVtYmVyKG8ub2Zmc2V0KSlcblxuICAgICAgcmV0dXJuIHRoaXNcbiAgICB9XG4gIH1cblxufSlcblxuU1ZHLlBhdHRlcm4gPSBTVkcuaW52ZW50KHtcbiAgLy8gSW5pdGlhbGl6ZSBub2RlXG4gIGNyZWF0ZTogJ3BhdHRlcm4nXG5cbiAgLy8gSW5oZXJpdCBmcm9tXG4sIGluaGVyaXQ6IFNWRy5Db250YWluZXJcblxuICAvLyBBZGQgY2xhc3MgbWV0aG9kc1xuLCBleHRlbmQ6IHtcbiAgICAvLyBSZXR1cm4gdGhlIGZpbGwgaWRcbiAgICBmaWxsOiBmdW5jdGlvbigpIHtcbiAgICAgIHJldHVybiAndXJsKCMnICsgdGhpcy5pZCgpICsgJyknXG4gICAgfVxuICAgIC8vIFVwZGF0ZSBwYXR0ZXJuIGJ5IHJlYnVpbGRpbmdcbiAgLCB1cGRhdGU6IGZ1bmN0aW9uKGJsb2NrKSB7XG4gICAgICAvLyByZW1vdmUgY29udGVudFxuICAgICAgdGhpcy5jbGVhcigpXG5cbiAgICAgIC8vIGludm9rZSBwYXNzZWQgYmxvY2tcbiAgICAgIGlmICh0eXBlb2YgYmxvY2sgPT0gJ2Z1bmN0aW9uJylcbiAgICAgICAgYmxvY2suY2FsbCh0aGlzLCB0aGlzKVxuXG4gICAgICByZXR1cm4gdGhpc1xuICAgIH1cbiAgICAvLyBBbGlhcyBzdHJpbmcgY29udmVydGlvbiB0byBmaWxsXG4gICwgdG9TdHJpbmc6IGZ1bmN0aW9uKCkge1xuICAgICAgcmV0dXJuIHRoaXMuZmlsbCgpXG4gICAgfVxuICAgIC8vIGN1c3RvbSBhdHRyIHRvIGhhbmRsZSB0cmFuc2Zvcm1cbiAgLCBhdHRyOiBmdW5jdGlvbihhLCBiLCBjKSB7XG4gICAgICBpZihhID09ICd0cmFuc2Zvcm0nKSBhID0gJ3BhdHRlcm5UcmFuc2Zvcm0nXG4gICAgICByZXR1cm4gU1ZHLkNvbnRhaW5lci5wcm90b3R5cGUuYXR0ci5jYWxsKHRoaXMsIGEsIGIsIGMpXG4gICAgfVxuXG4gIH1cblxuICAvLyBBZGQgcGFyZW50IG1ldGhvZFxuLCBjb25zdHJ1Y3Q6IHtcbiAgICAvLyBDcmVhdGUgcGF0dGVybiBlbGVtZW50IGluIGRlZnNcbiAgICBwYXR0ZXJuOiBmdW5jdGlvbih3aWR0aCwgaGVpZ2h0LCBibG9jaykge1xuICAgICAgcmV0dXJuIHRoaXMuZGVmcygpLnBhdHRlcm4od2lkdGgsIGhlaWdodCwgYmxvY2spXG4gICAgfVxuICB9XG59KVxuXG5TVkcuZXh0ZW5kKFNWRy5EZWZzLCB7XG4gIC8vIERlZmluZSBncmFkaWVudFxuICBwYXR0ZXJuOiBmdW5jdGlvbih3aWR0aCwgaGVpZ2h0LCBibG9jaykge1xuICAgIHJldHVybiB0aGlzLnB1dChuZXcgU1ZHLlBhdHRlcm4pLnVwZGF0ZShibG9jaykuYXR0cih7XG4gICAgICB4OiAgICAgICAgICAgIDBcbiAgICAsIHk6ICAgICAgICAgICAgMFxuICAgICwgd2lkdGg6ICAgICAgICB3aWR0aFxuICAgICwgaGVpZ2h0OiAgICAgICBoZWlnaHRcbiAgICAsIHBhdHRlcm5Vbml0czogJ3VzZXJTcGFjZU9uVXNlJ1xuICAgIH0pXG4gIH1cblxufSlcblNWRy5Eb2MgPSBTVkcuaW52ZW50KHtcbiAgLy8gSW5pdGlhbGl6ZSBub2RlXG4gIGNyZWF0ZTogZnVuY3Rpb24oZWxlbWVudCkge1xuICAgIGlmIChlbGVtZW50KSB7XG4gICAgICAvLyBlbnN1cmUgdGhlIHByZXNlbmNlIG9mIGEgZG9tIGVsZW1lbnRcbiAgICAgIGVsZW1lbnQgPSB0eXBlb2YgZWxlbWVudCA9PSAnc3RyaW5nJyA/XG4gICAgICAgIGRvY3VtZW50LmdldEVsZW1lbnRCeUlkKGVsZW1lbnQpIDpcbiAgICAgICAgZWxlbWVudFxuXG4gICAgICAvLyBJZiB0aGUgdGFyZ2V0IGlzIGFuIHN2ZyBlbGVtZW50LCB1c2UgdGhhdCBlbGVtZW50IGFzIHRoZSBtYWluIHdyYXBwZXIuXG4gICAgICAvLyBUaGlzIGFsbG93cyBzdmcuanMgdG8gd29yayB3aXRoIHN2ZyBkb2N1bWVudHMgYXMgd2VsbC5cbiAgICAgIGlmIChlbGVtZW50Lm5vZGVOYW1lID09ICdzdmcnKSB7XG4gICAgICAgIHRoaXMuY29uc3RydWN0b3IuY2FsbCh0aGlzLCBlbGVtZW50KVxuICAgICAgfSBlbHNlIHtcbiAgICAgICAgdGhpcy5jb25zdHJ1Y3Rvci5jYWxsKHRoaXMsIFNWRy5jcmVhdGUoJ3N2ZycpKVxuICAgICAgICBlbGVtZW50LmFwcGVuZENoaWxkKHRoaXMubm9kZSlcbiAgICAgICAgdGhpcy5zaXplKCcxMDAlJywgJzEwMCUnKVxuICAgICAgfVxuXG4gICAgICAvLyBzZXQgc3ZnIGVsZW1lbnQgYXR0cmlidXRlcyBhbmQgZW5zdXJlIGRlZnMgbm9kZVxuICAgICAgdGhpcy5uYW1lc3BhY2UoKS5kZWZzKClcbiAgICB9XG4gIH1cblxuICAvLyBJbmhlcml0IGZyb21cbiwgaW5oZXJpdDogU1ZHLkNvbnRhaW5lclxuXG4gIC8vIEFkZCBjbGFzcyBtZXRob2RzXG4sIGV4dGVuZDoge1xuICAgIC8vIEFkZCBuYW1lc3BhY2VzXG4gICAgbmFtZXNwYWNlOiBmdW5jdGlvbigpIHtcbiAgICAgIHJldHVybiB0aGlzXG4gICAgICAgIC5hdHRyKHsgeG1sbnM6IFNWRy5ucywgdmVyc2lvbjogJzEuMScgfSlcbiAgICAgICAgLmF0dHIoJ3htbG5zOnhsaW5rJywgU1ZHLnhsaW5rLCBTVkcueG1sbnMpXG4gICAgICAgIC5hdHRyKCd4bWxuczpzdmdqcycsIFNWRy5zdmdqcywgU1ZHLnhtbG5zKVxuICAgIH1cbiAgICAvLyBDcmVhdGVzIGFuZCByZXR1cm5zIGRlZnMgZWxlbWVudFxuICAsIGRlZnM6IGZ1bmN0aW9uKCkge1xuICAgICAgaWYgKCF0aGlzLl9kZWZzKSB7XG4gICAgICAgIHZhciBkZWZzXG5cbiAgICAgICAgLy8gRmluZCBvciBjcmVhdGUgYSBkZWZzIGVsZW1lbnQgaW4gdGhpcyBpbnN0YW5jZVxuICAgICAgICBpZiAoZGVmcyA9IHRoaXMubm9kZS5nZXRFbGVtZW50c0J5VGFnTmFtZSgnZGVmcycpWzBdKVxuICAgICAgICAgIHRoaXMuX2RlZnMgPSBTVkcuYWRvcHQoZGVmcylcbiAgICAgICAgZWxzZVxuICAgICAgICAgIHRoaXMuX2RlZnMgPSBuZXcgU1ZHLkRlZnNcblxuICAgICAgICAvLyBNYWtlIHN1cmUgdGhlIGRlZnMgbm9kZSBpcyBhdCB0aGUgZW5kIG9mIHRoZSBzdGFja1xuICAgICAgICB0aGlzLm5vZGUuYXBwZW5kQ2hpbGQodGhpcy5fZGVmcy5ub2RlKVxuICAgICAgfVxuXG4gICAgICByZXR1cm4gdGhpcy5fZGVmc1xuICAgIH1cbiAgICAvLyBjdXN0b20gcGFyZW50IG1ldGhvZFxuICAsIHBhcmVudDogZnVuY3Rpb24oKSB7XG4gICAgICByZXR1cm4gdGhpcy5ub2RlLnBhcmVudE5vZGUubm9kZU5hbWUgPT0gJyNkb2N1bWVudCcgPyBudWxsIDogdGhpcy5ub2RlLnBhcmVudE5vZGVcbiAgICB9XG4gICAgLy8gRml4IGZvciBwb3NzaWJsZSBzdWItcGl4ZWwgb2Zmc2V0LiBTZWU6XG4gICAgLy8gaHR0cHM6Ly9idWd6aWxsYS5tb3ppbGxhLm9yZy9zaG93X2J1Zy5jZ2k/aWQ9NjA4ODEyXG4gICwgc3BvZjogZnVuY3Rpb24oc3BvZikge1xuICAgICAgdmFyIHBvcyA9IHRoaXMubm9kZS5nZXRTY3JlZW5DVE0oKVxuXG4gICAgICBpZiAocG9zKVxuICAgICAgICB0aGlzXG4gICAgICAgICAgLnN0eWxlKCdsZWZ0JywgKC1wb3MuZSAlIDEpICsgJ3B4JylcbiAgICAgICAgICAuc3R5bGUoJ3RvcCcsICAoLXBvcy5mICUgMSkgKyAncHgnKVxuXG4gICAgICByZXR1cm4gdGhpc1xuICAgIH1cblxuICAgICAgLy8gUmVtb3ZlcyB0aGUgZG9jIGZyb20gdGhlIERPTVxuICAsIHJlbW92ZTogZnVuY3Rpb24oKSB7XG4gICAgICBpZih0aGlzLnBhcmVudCgpKSB7XG4gICAgICAgIHRoaXMucGFyZW50KCkucmVtb3ZlQ2hpbGQodGhpcy5ub2RlKTtcbiAgICAgIH1cblxuICAgICAgcmV0dXJuIHRoaXM7XG4gICAgfVxuICB9XG5cbn0pXG5cblNWRy5TaGFwZSA9IFNWRy5pbnZlbnQoe1xuICAvLyBJbml0aWFsaXplIG5vZGVcbiAgY3JlYXRlOiBmdW5jdGlvbihlbGVtZW50KSB7XG4gICAgdGhpcy5jb25zdHJ1Y3Rvci5jYWxsKHRoaXMsIGVsZW1lbnQpXG4gIH1cblxuICAvLyBJbmhlcml0IGZyb21cbiwgaW5oZXJpdDogU1ZHLkVsZW1lbnRcblxufSlcblxuU1ZHLkJhcmUgPSBTVkcuaW52ZW50KHtcbiAgLy8gSW5pdGlhbGl6ZVxuICBjcmVhdGU6IGZ1bmN0aW9uKGVsZW1lbnQsIGluaGVyaXQpIHtcbiAgICAvLyBjb25zdHJ1Y3QgZWxlbWVudFxuICAgIHRoaXMuY29uc3RydWN0b3IuY2FsbCh0aGlzLCBTVkcuY3JlYXRlKGVsZW1lbnQpKVxuXG4gICAgLy8gaW5oZXJpdCBjdXN0b20gbWV0aG9kc1xuICAgIGlmIChpbmhlcml0KVxuICAgICAgZm9yICh2YXIgbWV0aG9kIGluIGluaGVyaXQucHJvdG90eXBlKVxuICAgICAgICBpZiAodHlwZW9mIGluaGVyaXQucHJvdG90eXBlW21ldGhvZF0gPT09ICdmdW5jdGlvbicpXG4gICAgICAgICAgdGhpc1ttZXRob2RdID0gaW5oZXJpdC5wcm90b3R5cGVbbWV0aG9kXVxuICB9XG5cbiAgLy8gSW5oZXJpdCBmcm9tXG4sIGluaGVyaXQ6IFNWRy5FbGVtZW50XG5cbiAgLy8gQWRkIG1ldGhvZHNcbiwgZXh0ZW5kOiB7XG4gICAgLy8gSW5zZXJ0IHNvbWUgcGxhaW4gdGV4dFxuICAgIHdvcmRzOiBmdW5jdGlvbih0ZXh0KSB7XG4gICAgICAvLyByZW1vdmUgY29udGVudHNcbiAgICAgIHdoaWxlICh0aGlzLm5vZGUuaGFzQ2hpbGROb2RlcygpKVxuICAgICAgICB0aGlzLm5vZGUucmVtb3ZlQ2hpbGQodGhpcy5ub2RlLmxhc3RDaGlsZClcblxuICAgICAgLy8gY3JlYXRlIHRleHQgbm9kZVxuICAgICAgdGhpcy5ub2RlLmFwcGVuZENoaWxkKGRvY3VtZW50LmNyZWF0ZVRleHROb2RlKHRleHQpKVxuXG4gICAgICByZXR1cm4gdGhpc1xuICAgIH1cbiAgfVxufSlcblxuXG5TVkcuZXh0ZW5kKFNWRy5QYXJlbnQsIHtcbiAgLy8gQ3JlYXRlIGFuIGVsZW1lbnQgdGhhdCBpcyBub3QgZGVzY3JpYmVkIGJ5IFNWRy5qc1xuICBlbGVtZW50OiBmdW5jdGlvbihlbGVtZW50LCBpbmhlcml0KSB7XG4gICAgcmV0dXJuIHRoaXMucHV0KG5ldyBTVkcuQmFyZShlbGVtZW50LCBpbmhlcml0KSlcbiAgfVxuICAvLyBBZGQgc3ltYm9sIGVsZW1lbnRcbiwgc3ltYm9sOiBmdW5jdGlvbigpIHtcbiAgICByZXR1cm4gdGhpcy5kZWZzKCkuZWxlbWVudCgnc3ltYm9sJywgU1ZHLkNvbnRhaW5lcilcbiAgfVxuXG59KVxuU1ZHLlVzZSA9IFNWRy5pbnZlbnQoe1xuICAvLyBJbml0aWFsaXplIG5vZGVcbiAgY3JlYXRlOiAndXNlJ1xuXG4gIC8vIEluaGVyaXQgZnJvbVxuLCBpbmhlcml0OiBTVkcuU2hhcGVcblxuICAvLyBBZGQgY2xhc3MgbWV0aG9kc1xuLCBleHRlbmQ6IHtcbiAgICAvLyBVc2UgZWxlbWVudCBhcyBhIHJlZmVyZW5jZVxuICAgIGVsZW1lbnQ6IGZ1bmN0aW9uKGVsZW1lbnQsIGZpbGUpIHtcbiAgICAgIC8vIFNldCBsaW5lZCBlbGVtZW50XG4gICAgICByZXR1cm4gdGhpcy5hdHRyKCdocmVmJywgKGZpbGUgfHwgJycpICsgJyMnICsgZWxlbWVudCwgU1ZHLnhsaW5rKVxuICAgIH1cbiAgfVxuXG4gIC8vIEFkZCBwYXJlbnQgbWV0aG9kXG4sIGNvbnN0cnVjdDoge1xuICAgIC8vIENyZWF0ZSBhIHVzZSBlbGVtZW50XG4gICAgdXNlOiBmdW5jdGlvbihlbGVtZW50LCBmaWxlKSB7XG4gICAgICByZXR1cm4gdGhpcy5wdXQobmV3IFNWRy5Vc2UpLmVsZW1lbnQoZWxlbWVudCwgZmlsZSlcbiAgICB9XG4gIH1cbn0pXG5TVkcuUmVjdCA9IFNWRy5pbnZlbnQoe1xuICAvLyBJbml0aWFsaXplIG5vZGVcbiAgY3JlYXRlOiAncmVjdCdcblxuICAvLyBJbmhlcml0IGZyb21cbiwgaW5oZXJpdDogU1ZHLlNoYXBlXG5cbiAgLy8gQWRkIHBhcmVudCBtZXRob2RcbiwgY29uc3RydWN0OiB7XG4gICAgLy8gQ3JlYXRlIGEgcmVjdCBlbGVtZW50XG4gICAgcmVjdDogZnVuY3Rpb24od2lkdGgsIGhlaWdodCkge1xuICAgICAgcmV0dXJuIHRoaXMucHV0KG5ldyBTVkcuUmVjdCgpKS5zaXplKHdpZHRoLCBoZWlnaHQpXG4gICAgfVxuICB9XG59KVxuU1ZHLkNpcmNsZSA9IFNWRy5pbnZlbnQoe1xuICAvLyBJbml0aWFsaXplIG5vZGVcbiAgY3JlYXRlOiAnY2lyY2xlJ1xuXG4gIC8vIEluaGVyaXQgZnJvbVxuLCBpbmhlcml0OiBTVkcuU2hhcGVcblxuICAvLyBBZGQgcGFyZW50IG1ldGhvZFxuLCBjb25zdHJ1Y3Q6IHtcbiAgICAvLyBDcmVhdGUgY2lyY2xlIGVsZW1lbnQsIGJhc2VkIG9uIGVsbGlwc2VcbiAgICBjaXJjbGU6IGZ1bmN0aW9uKHNpemUpIHtcbiAgICAgIHJldHVybiB0aGlzLnB1dChuZXcgU1ZHLkNpcmNsZSkucngobmV3IFNWRy5OdW1iZXIoc2l6ZSkuZGl2aWRlKDIpKS5tb3ZlKDAsIDApXG4gICAgfVxuICB9XG59KVxuXG5TVkcuZXh0ZW5kKFNWRy5DaXJjbGUsIFNWRy5GWCwge1xuICAvLyBSYWRpdXMgeCB2YWx1ZVxuICByeDogZnVuY3Rpb24ocngpIHtcbiAgICByZXR1cm4gdGhpcy5hdHRyKCdyJywgcngpXG4gIH1cbiAgLy8gQWxpYXMgcmFkaXVzIHggdmFsdWVcbiwgcnk6IGZ1bmN0aW9uKHJ5KSB7XG4gICAgcmV0dXJuIHRoaXMucngocnkpXG4gIH1cbn0pXG5cblNWRy5FbGxpcHNlID0gU1ZHLmludmVudCh7XG4gIC8vIEluaXRpYWxpemUgbm9kZVxuICBjcmVhdGU6ICdlbGxpcHNlJ1xuXG4gIC8vIEluaGVyaXQgZnJvbVxuLCBpbmhlcml0OiBTVkcuU2hhcGVcblxuICAvLyBBZGQgcGFyZW50IG1ldGhvZFxuLCBjb25zdHJ1Y3Q6IHtcbiAgICAvLyBDcmVhdGUgYW4gZWxsaXBzZVxuICAgIGVsbGlwc2U6IGZ1bmN0aW9uKHdpZHRoLCBoZWlnaHQpIHtcbiAgICAgIHJldHVybiB0aGlzLnB1dChuZXcgU1ZHLkVsbGlwc2UpLnNpemUod2lkdGgsIGhlaWdodCkubW92ZSgwLCAwKVxuICAgIH1cbiAgfVxufSlcblxuU1ZHLmV4dGVuZChTVkcuRWxsaXBzZSwgU1ZHLlJlY3QsIFNWRy5GWCwge1xuICAvLyBSYWRpdXMgeCB2YWx1ZVxuICByeDogZnVuY3Rpb24ocngpIHtcbiAgICByZXR1cm4gdGhpcy5hdHRyKCdyeCcsIHJ4KVxuICB9XG4gIC8vIFJhZGl1cyB5IHZhbHVlXG4sIHJ5OiBmdW5jdGlvbihyeSkge1xuICAgIHJldHVybiB0aGlzLmF0dHIoJ3J5JywgcnkpXG4gIH1cbn0pXG5cbi8vIEFkZCBjb21tb24gbWV0aG9kXG5TVkcuZXh0ZW5kKFNWRy5DaXJjbGUsIFNWRy5FbGxpcHNlLCB7XG4gICAgLy8gTW92ZSBvdmVyIHgtYXhpc1xuICAgIHg6IGZ1bmN0aW9uKHgpIHtcbiAgICAgIHJldHVybiB4ID09IG51bGwgPyB0aGlzLmN4KCkgLSB0aGlzLnJ4KCkgOiB0aGlzLmN4KHggKyB0aGlzLnJ4KCkpXG4gICAgfVxuICAgIC8vIE1vdmUgb3ZlciB5LWF4aXNcbiAgLCB5OiBmdW5jdGlvbih5KSB7XG4gICAgICByZXR1cm4geSA9PSBudWxsID8gdGhpcy5jeSgpIC0gdGhpcy5yeSgpIDogdGhpcy5jeSh5ICsgdGhpcy5yeSgpKVxuICAgIH1cbiAgICAvLyBNb3ZlIGJ5IGNlbnRlciBvdmVyIHgtYXhpc1xuICAsIGN4OiBmdW5jdGlvbih4KSB7XG4gICAgICByZXR1cm4geCA9PSBudWxsID8gdGhpcy5hdHRyKCdjeCcpIDogdGhpcy5hdHRyKCdjeCcsIHgpXG4gICAgfVxuICAgIC8vIE1vdmUgYnkgY2VudGVyIG92ZXIgeS1heGlzXG4gICwgY3k6IGZ1bmN0aW9uKHkpIHtcbiAgICAgIHJldHVybiB5ID09IG51bGwgPyB0aGlzLmF0dHIoJ2N5JykgOiB0aGlzLmF0dHIoJ2N5JywgeSlcbiAgICB9XG4gICAgLy8gU2V0IHdpZHRoIG9mIGVsZW1lbnRcbiAgLCB3aWR0aDogZnVuY3Rpb24od2lkdGgpIHtcbiAgICAgIHJldHVybiB3aWR0aCA9PSBudWxsID8gdGhpcy5yeCgpICogMiA6IHRoaXMucngobmV3IFNWRy5OdW1iZXIod2lkdGgpLmRpdmlkZSgyKSlcbiAgICB9XG4gICAgLy8gU2V0IGhlaWdodCBvZiBlbGVtZW50XG4gICwgaGVpZ2h0OiBmdW5jdGlvbihoZWlnaHQpIHtcbiAgICAgIHJldHVybiBoZWlnaHQgPT0gbnVsbCA/IHRoaXMucnkoKSAqIDIgOiB0aGlzLnJ5KG5ldyBTVkcuTnVtYmVyKGhlaWdodCkuZGl2aWRlKDIpKVxuICAgIH1cbiAgICAvLyBDdXN0b20gc2l6ZSBmdW5jdGlvblxuICAsIHNpemU6IGZ1bmN0aW9uKHdpZHRoLCBoZWlnaHQpIHtcbiAgICAgIHZhciBwID0gcHJvcG9ydGlvbmFsU2l6ZSh0aGlzLCB3aWR0aCwgaGVpZ2h0KVxuXG4gICAgICByZXR1cm4gdGhpc1xuICAgICAgICAucngobmV3IFNWRy5OdW1iZXIocC53aWR0aCkuZGl2aWRlKDIpKVxuICAgICAgICAucnkobmV3IFNWRy5OdW1iZXIocC5oZWlnaHQpLmRpdmlkZSgyKSlcbiAgICB9XG59KVxuU1ZHLkxpbmUgPSBTVkcuaW52ZW50KHtcbiAgLy8gSW5pdGlhbGl6ZSBub2RlXG4gIGNyZWF0ZTogJ2xpbmUnXG5cbiAgLy8gSW5oZXJpdCBmcm9tXG4sIGluaGVyaXQ6IFNWRy5TaGFwZVxuXG4gIC8vIEFkZCBjbGFzcyBtZXRob2RzXG4sIGV4dGVuZDoge1xuICAgIC8vIEdldCBhcnJheVxuICAgIGFycmF5OiBmdW5jdGlvbigpIHtcbiAgICAgIHJldHVybiBuZXcgU1ZHLlBvaW50QXJyYXkoW1xuICAgICAgICBbIHRoaXMuYXR0cigneDEnKSwgdGhpcy5hdHRyKCd5MScpIF1cbiAgICAgICwgWyB0aGlzLmF0dHIoJ3gyJyksIHRoaXMuYXR0cigneTInKSBdXG4gICAgICBdKVxuICAgIH1cbiAgICAvLyBPdmVyd3JpdGUgbmF0aXZlIHBsb3QoKSBtZXRob2RcbiAgLCBwbG90OiBmdW5jdGlvbih4MSwgeTEsIHgyLCB5Mikge1xuICAgICAgaWYgKHR5cGVvZiB5MSAhPT0gJ3VuZGVmaW5lZCcpXG4gICAgICAgIHgxID0geyB4MTogeDEsIHkxOiB5MSwgeDI6IHgyLCB5MjogeTIgfVxuICAgICAgZWxzZVxuICAgICAgICB4MSA9IG5ldyBTVkcuUG9pbnRBcnJheSh4MSkudG9MaW5lKClcblxuICAgICAgcmV0dXJuIHRoaXMuYXR0cih4MSlcbiAgICB9XG4gICAgLy8gTW92ZSBieSBsZWZ0IHRvcCBjb3JuZXJcbiAgLCBtb3ZlOiBmdW5jdGlvbih4LCB5KSB7XG4gICAgICByZXR1cm4gdGhpcy5hdHRyKHRoaXMuYXJyYXkoKS5tb3ZlKHgsIHkpLnRvTGluZSgpKVxuICAgIH1cbiAgICAvLyBTZXQgZWxlbWVudCBzaXplIHRvIGdpdmVuIHdpZHRoIGFuZCBoZWlnaHRcbiAgLCBzaXplOiBmdW5jdGlvbih3aWR0aCwgaGVpZ2h0KSB7XG4gICAgICB2YXIgcCA9IHByb3BvcnRpb25hbFNpemUodGhpcywgd2lkdGgsIGhlaWdodClcblxuICAgICAgcmV0dXJuIHRoaXMuYXR0cih0aGlzLmFycmF5KCkuc2l6ZShwLndpZHRoLCBwLmhlaWdodCkudG9MaW5lKCkpXG4gICAgfVxuICB9XG5cbiAgLy8gQWRkIHBhcmVudCBtZXRob2RcbiwgY29uc3RydWN0OiB7XG4gICAgLy8gQ3JlYXRlIGEgbGluZSBlbGVtZW50XG4gICAgbGluZTogZnVuY3Rpb24oeDEsIHkxLCB4MiwgeTIpIHtcbiAgICAgIHJldHVybiB0aGlzLnB1dChuZXcgU1ZHLkxpbmUpLnBsb3QoeDEsIHkxLCB4MiwgeTIpXG4gICAgfVxuICB9XG59KVxuXG5TVkcuUG9seWxpbmUgPSBTVkcuaW52ZW50KHtcbiAgLy8gSW5pdGlhbGl6ZSBub2RlXG4gIGNyZWF0ZTogJ3BvbHlsaW5lJ1xuXG4gIC8vIEluaGVyaXQgZnJvbVxuLCBpbmhlcml0OiBTVkcuU2hhcGVcblxuICAvLyBBZGQgcGFyZW50IG1ldGhvZFxuLCBjb25zdHJ1Y3Q6IHtcbiAgICAvLyBDcmVhdGUgYSB3cmFwcGVkIHBvbHlsaW5lIGVsZW1lbnRcbiAgICBwb2x5bGluZTogZnVuY3Rpb24ocCkge1xuICAgICAgcmV0dXJuIHRoaXMucHV0KG5ldyBTVkcuUG9seWxpbmUpLnBsb3QocClcbiAgICB9XG4gIH1cbn0pXG5cblNWRy5Qb2x5Z29uID0gU1ZHLmludmVudCh7XG4gIC8vIEluaXRpYWxpemUgbm9kZVxuICBjcmVhdGU6ICdwb2x5Z29uJ1xuXG4gIC8vIEluaGVyaXQgZnJvbVxuLCBpbmhlcml0OiBTVkcuU2hhcGVcblxuICAvLyBBZGQgcGFyZW50IG1ldGhvZFxuLCBjb25zdHJ1Y3Q6IHtcbiAgICAvLyBDcmVhdGUgYSB3cmFwcGVkIHBvbHlnb24gZWxlbWVudFxuICAgIHBvbHlnb246IGZ1bmN0aW9uKHApIHtcbiAgICAgIHJldHVybiB0aGlzLnB1dChuZXcgU1ZHLlBvbHlnb24pLnBsb3QocClcbiAgICB9XG4gIH1cbn0pXG5cbi8vIEFkZCBwb2x5Z29uLXNwZWNpZmljIGZ1bmN0aW9uc1xuU1ZHLmV4dGVuZChTVkcuUG9seWxpbmUsIFNWRy5Qb2x5Z29uLCB7XG4gIC8vIEdldCBhcnJheVxuICBhcnJheTogZnVuY3Rpb24oKSB7XG4gICAgcmV0dXJuIHRoaXMuX2FycmF5IHx8ICh0aGlzLl9hcnJheSA9IG5ldyBTVkcuUG9pbnRBcnJheSh0aGlzLmF0dHIoJ3BvaW50cycpKSlcbiAgfVxuICAvLyBQbG90IG5ldyBwYXRoXG4sIHBsb3Q6IGZ1bmN0aW9uKHApIHtcbiAgICByZXR1cm4gdGhpcy5hdHRyKCdwb2ludHMnLCAodGhpcy5fYXJyYXkgPSBuZXcgU1ZHLlBvaW50QXJyYXkocCkpKVxuICB9XG4gIC8vIE1vdmUgYnkgbGVmdCB0b3AgY29ybmVyXG4sIG1vdmU6IGZ1bmN0aW9uKHgsIHkpIHtcbiAgICByZXR1cm4gdGhpcy5hdHRyKCdwb2ludHMnLCB0aGlzLmFycmF5KCkubW92ZSh4LCB5KSlcbiAgfVxuICAvLyBTZXQgZWxlbWVudCBzaXplIHRvIGdpdmVuIHdpZHRoIGFuZCBoZWlnaHRcbiwgc2l6ZTogZnVuY3Rpb24od2lkdGgsIGhlaWdodCkge1xuICAgIHZhciBwID0gcHJvcG9ydGlvbmFsU2l6ZSh0aGlzLCB3aWR0aCwgaGVpZ2h0KVxuXG4gICAgcmV0dXJuIHRoaXMuYXR0cigncG9pbnRzJywgdGhpcy5hcnJheSgpLnNpemUocC53aWR0aCwgcC5oZWlnaHQpKVxuICB9XG5cbn0pXG4vLyB1bmlmeSBhbGwgcG9pbnQgdG8gcG9pbnQgZWxlbWVudHNcblNWRy5leHRlbmQoU1ZHLkxpbmUsIFNWRy5Qb2x5bGluZSwgU1ZHLlBvbHlnb24sIHtcbiAgLy8gRGVmaW5lIG1vcnBoYWJsZSBhcnJheVxuICBtb3JwaEFycmF5OiAgU1ZHLlBvaW50QXJyYXlcbiAgLy8gTW92ZSBieSBsZWZ0IHRvcCBjb3JuZXIgb3ZlciB4LWF4aXNcbiwgeDogZnVuY3Rpb24oeCkge1xuICAgIHJldHVybiB4ID09IG51bGwgPyB0aGlzLmJib3goKS54IDogdGhpcy5tb3ZlKHgsIHRoaXMuYmJveCgpLnkpXG4gIH1cbiAgLy8gTW92ZSBieSBsZWZ0IHRvcCBjb3JuZXIgb3ZlciB5LWF4aXNcbiwgeTogZnVuY3Rpb24oeSkge1xuICAgIHJldHVybiB5ID09IG51bGwgPyB0aGlzLmJib3goKS55IDogdGhpcy5tb3ZlKHRoaXMuYmJveCgpLngsIHkpXG4gIH1cbiAgLy8gU2V0IHdpZHRoIG9mIGVsZW1lbnRcbiwgd2lkdGg6IGZ1bmN0aW9uKHdpZHRoKSB7XG4gICAgdmFyIGIgPSB0aGlzLmJib3goKVxuXG4gICAgcmV0dXJuIHdpZHRoID09IG51bGwgPyBiLndpZHRoIDogdGhpcy5zaXplKHdpZHRoLCBiLmhlaWdodClcbiAgfVxuICAvLyBTZXQgaGVpZ2h0IG9mIGVsZW1lbnRcbiwgaGVpZ2h0OiBmdW5jdGlvbihoZWlnaHQpIHtcbiAgICB2YXIgYiA9IHRoaXMuYmJveCgpXG5cbiAgICByZXR1cm4gaGVpZ2h0ID09IG51bGwgPyBiLmhlaWdodCA6IHRoaXMuc2l6ZShiLndpZHRoLCBoZWlnaHQpXG4gIH1cbn0pXG5TVkcuUGF0aCA9IFNWRy5pbnZlbnQoe1xuICAvLyBJbml0aWFsaXplIG5vZGVcbiAgY3JlYXRlOiAncGF0aCdcblxuICAvLyBJbmhlcml0IGZyb21cbiwgaW5oZXJpdDogU1ZHLlNoYXBlXG5cbiAgLy8gQWRkIGNsYXNzIG1ldGhvZHNcbiwgZXh0ZW5kOiB7XG4gICAgLy8gRGVmaW5lIG1vcnBoYWJsZSBhcnJheVxuICAgIG1vcnBoQXJyYXk6ICBTVkcuUGF0aEFycmF5XG4gICAgLy8gR2V0IGFycmF5XG4gICwgYXJyYXk6IGZ1bmN0aW9uKCkge1xuICAgICAgcmV0dXJuIHRoaXMuX2FycmF5IHx8ICh0aGlzLl9hcnJheSA9IG5ldyBTVkcuUGF0aEFycmF5KHRoaXMuYXR0cignZCcpKSlcbiAgICB9XG4gICAgLy8gUGxvdCBuZXcgcG9seSBwb2ludHNcbiAgLCBwbG90OiBmdW5jdGlvbihwKSB7XG4gICAgICByZXR1cm4gdGhpcy5hdHRyKCdkJywgKHRoaXMuX2FycmF5ID0gbmV3IFNWRy5QYXRoQXJyYXkocCkpKVxuICAgIH1cbiAgICAvLyBNb3ZlIGJ5IGxlZnQgdG9wIGNvcm5lclxuICAsIG1vdmU6IGZ1bmN0aW9uKHgsIHkpIHtcbiAgICAgIHJldHVybiB0aGlzLmF0dHIoJ2QnLCB0aGlzLmFycmF5KCkubW92ZSh4LCB5KSlcbiAgICB9XG4gICAgLy8gTW92ZSBieSBsZWZ0IHRvcCBjb3JuZXIgb3ZlciB4LWF4aXNcbiAgLCB4OiBmdW5jdGlvbih4KSB7XG4gICAgICByZXR1cm4geCA9PSBudWxsID8gdGhpcy5iYm94KCkueCA6IHRoaXMubW92ZSh4LCB0aGlzLmJib3goKS55KVxuICAgIH1cbiAgICAvLyBNb3ZlIGJ5IGxlZnQgdG9wIGNvcm5lciBvdmVyIHktYXhpc1xuICAsIHk6IGZ1bmN0aW9uKHkpIHtcbiAgICAgIHJldHVybiB5ID09IG51bGwgPyB0aGlzLmJib3goKS55IDogdGhpcy5tb3ZlKHRoaXMuYmJveCgpLngsIHkpXG4gICAgfVxuICAgIC8vIFNldCBlbGVtZW50IHNpemUgdG8gZ2l2ZW4gd2lkdGggYW5kIGhlaWdodFxuICAsIHNpemU6IGZ1bmN0aW9uKHdpZHRoLCBoZWlnaHQpIHtcbiAgICAgIHZhciBwID0gcHJvcG9ydGlvbmFsU2l6ZSh0aGlzLCB3aWR0aCwgaGVpZ2h0KVxuXG4gICAgICByZXR1cm4gdGhpcy5hdHRyKCdkJywgdGhpcy5hcnJheSgpLnNpemUocC53aWR0aCwgcC5oZWlnaHQpKVxuICAgIH1cbiAgICAvLyBTZXQgd2lkdGggb2YgZWxlbWVudFxuICAsIHdpZHRoOiBmdW5jdGlvbih3aWR0aCkge1xuICAgICAgcmV0dXJuIHdpZHRoID09IG51bGwgPyB0aGlzLmJib3goKS53aWR0aCA6IHRoaXMuc2l6ZSh3aWR0aCwgdGhpcy5iYm94KCkuaGVpZ2h0KVxuICAgIH1cbiAgICAvLyBTZXQgaGVpZ2h0IG9mIGVsZW1lbnRcbiAgLCBoZWlnaHQ6IGZ1bmN0aW9uKGhlaWdodCkge1xuICAgICAgcmV0dXJuIGhlaWdodCA9PSBudWxsID8gdGhpcy5iYm94KCkuaGVpZ2h0IDogdGhpcy5zaXplKHRoaXMuYmJveCgpLndpZHRoLCBoZWlnaHQpXG4gICAgfVxuXG4gIH1cblxuICAvLyBBZGQgcGFyZW50IG1ldGhvZFxuLCBjb25zdHJ1Y3Q6IHtcbiAgICAvLyBDcmVhdGUgYSB3cmFwcGVkIHBhdGggZWxlbWVudFxuICAgIHBhdGg6IGZ1bmN0aW9uKGQpIHtcbiAgICAgIHJldHVybiB0aGlzLnB1dChuZXcgU1ZHLlBhdGgpLnBsb3QoZClcbiAgICB9XG4gIH1cbn0pXG5TVkcuSW1hZ2UgPSBTVkcuaW52ZW50KHtcbiAgLy8gSW5pdGlhbGl6ZSBub2RlXG4gIGNyZWF0ZTogJ2ltYWdlJ1xuXG4gIC8vIEluaGVyaXQgZnJvbVxuLCBpbmhlcml0OiBTVkcuU2hhcGVcblxuICAvLyBBZGQgY2xhc3MgbWV0aG9kc1xuLCBleHRlbmQ6IHtcbiAgICAvLyAocmUpbG9hZCBpbWFnZVxuICAgIGxvYWQ6IGZ1bmN0aW9uKHVybCkge1xuICAgICAgaWYgKCF1cmwpIHJldHVybiB0aGlzXG5cbiAgICAgIHZhciBzZWxmID0gdGhpc1xuICAgICAgICAsIGltZyAgPSBkb2N1bWVudC5jcmVhdGVFbGVtZW50KCdpbWcnKVxuXG4gICAgICAvLyBwcmVsb2FkIGltYWdlXG4gICAgICBpbWcub25sb2FkID0gZnVuY3Rpb24oKSB7XG4gICAgICAgIHZhciBwID0gc2VsZi5wYXJlbnQoU1ZHLlBhdHRlcm4pXG5cbiAgICAgICAgaWYocCA9PT0gbnVsbCkgcmV0dXJuXG5cbiAgICAgICAgLy8gZW5zdXJlIGltYWdlIHNpemVcbiAgICAgICAgaWYgKHNlbGYud2lkdGgoKSA9PSAwICYmIHNlbGYuaGVpZ2h0KCkgPT0gMClcbiAgICAgICAgICBzZWxmLnNpemUoaW1nLndpZHRoLCBpbWcuaGVpZ2h0KVxuXG4gICAgICAgIC8vIGVuc3VyZSBwYXR0ZXJuIHNpemUgaWYgbm90IHNldFxuICAgICAgICBpZiAocCAmJiBwLndpZHRoKCkgPT0gMCAmJiBwLmhlaWdodCgpID09IDApXG4gICAgICAgICAgcC5zaXplKHNlbGYud2lkdGgoKSwgc2VsZi5oZWlnaHQoKSlcblxuICAgICAgICAvLyBjYWxsYmFja1xuICAgICAgICBpZiAodHlwZW9mIHNlbGYuX2xvYWRlZCA9PT0gJ2Z1bmN0aW9uJylcbiAgICAgICAgICBzZWxmLl9sb2FkZWQuY2FsbChzZWxmLCB7XG4gICAgICAgICAgICB3aWR0aDogIGltZy53aWR0aFxuICAgICAgICAgICwgaGVpZ2h0OiBpbWcuaGVpZ2h0XG4gICAgICAgICAgLCByYXRpbzogIGltZy53aWR0aCAvIGltZy5oZWlnaHRcbiAgICAgICAgICAsIHVybDogICAgdXJsXG4gICAgICAgICAgfSlcbiAgICAgIH1cblxuICAgICAgaW1nLm9uZXJyb3IgPSBmdW5jdGlvbihlKXtcbiAgICAgICAgaWYgKHR5cGVvZiBzZWxmLl9lcnJvciA9PT0gJ2Z1bmN0aW9uJyl7XG4gICAgICAgICAgICBzZWxmLl9lcnJvci5jYWxsKHNlbGYsIGUpXG4gICAgICAgIH1cbiAgICAgIH1cblxuICAgICAgcmV0dXJuIHRoaXMuYXR0cignaHJlZicsIChpbWcuc3JjID0gdGhpcy5zcmMgPSB1cmwpLCBTVkcueGxpbmspXG4gICAgfVxuICAgIC8vIEFkZCBsb2FkZWQgY2FsbGJhY2tcbiAgLCBsb2FkZWQ6IGZ1bmN0aW9uKGxvYWRlZCkge1xuICAgICAgdGhpcy5fbG9hZGVkID0gbG9hZGVkXG4gICAgICByZXR1cm4gdGhpc1xuICAgIH1cblxuICAsIGVycm9yOiBmdW5jdGlvbihlcnJvcikge1xuICAgICAgdGhpcy5fZXJyb3IgPSBlcnJvclxuICAgICAgcmV0dXJuIHRoaXNcbiAgICB9XG4gIH1cblxuICAvLyBBZGQgcGFyZW50IG1ldGhvZFxuLCBjb25zdHJ1Y3Q6IHtcbiAgICAvLyBjcmVhdGUgaW1hZ2UgZWxlbWVudCwgbG9hZCBpbWFnZSBhbmQgc2V0IGl0cyBzaXplXG4gICAgaW1hZ2U6IGZ1bmN0aW9uKHNvdXJjZSwgd2lkdGgsIGhlaWdodCkge1xuICAgICAgcmV0dXJuIHRoaXMucHV0KG5ldyBTVkcuSW1hZ2UpLmxvYWQoc291cmNlKS5zaXplKHdpZHRoIHx8IDAsIGhlaWdodCB8fCB3aWR0aCB8fCAwKVxuICAgIH1cbiAgfVxuXG59KVxuU1ZHLlRleHQgPSBTVkcuaW52ZW50KHtcbiAgLy8gSW5pdGlhbGl6ZSBub2RlXG4gIGNyZWF0ZTogZnVuY3Rpb24oKSB7XG4gICAgdGhpcy5jb25zdHJ1Y3Rvci5jYWxsKHRoaXMsIFNWRy5jcmVhdGUoJ3RleHQnKSlcblxuICAgIHRoaXMuZG9tLmxlYWRpbmcgPSBuZXcgU1ZHLk51bWJlcigxLjMpICAgIC8vIHN0b3JlIGxlYWRpbmcgdmFsdWUgZm9yIHJlYnVpbGRpbmdcbiAgICB0aGlzLl9yZWJ1aWxkID0gdHJ1ZSAgICAgICAgICAgICAgICAgICAgICAvLyBlbmFibGUgYXV0b21hdGljIHVwZGF0aW5nIG9mIGR5IHZhbHVlc1xuICAgIHRoaXMuX2J1aWxkICAgPSBmYWxzZSAgICAgICAgICAgICAgICAgICAgIC8vIGRpc2FibGUgYnVpbGQgbW9kZSBmb3IgYWRkaW5nIG11bHRpcGxlIGxpbmVzXG5cbiAgICAvLyBzZXQgZGVmYXVsdCBmb250XG4gICAgdGhpcy5hdHRyKCdmb250LWZhbWlseScsIFNWRy5kZWZhdWx0cy5hdHRyc1snZm9udC1mYW1pbHknXSlcbiAgfVxuXG4gIC8vIEluaGVyaXQgZnJvbVxuLCBpbmhlcml0OiBTVkcuU2hhcGVcblxuICAvLyBBZGQgY2xhc3MgbWV0aG9kc1xuLCBleHRlbmQ6IHtcbiAgICAvLyBNb3ZlIG92ZXIgeC1heGlzXG4gICAgeDogZnVuY3Rpb24oeCkge1xuICAgICAgLy8gYWN0IGFzIGdldHRlclxuICAgICAgaWYgKHggPT0gbnVsbClcbiAgICAgICAgcmV0dXJuIHRoaXMuYXR0cigneCcpXG5cbiAgICAgIC8vIG1vdmUgbGluZXMgYXMgd2VsbCBpZiBubyB0ZXh0UGF0aCBpcyBwcmVzZW50XG4gICAgICBpZiAoIXRoaXMudGV4dFBhdGgpXG4gICAgICAgIHRoaXMubGluZXMoKS5lYWNoKGZ1bmN0aW9uKCkgeyBpZiAodGhpcy5kb20ubmV3TGluZWQpIHRoaXMueCh4KSB9KVxuXG4gICAgICByZXR1cm4gdGhpcy5hdHRyKCd4JywgeClcbiAgICB9XG4gICAgLy8gTW92ZSBvdmVyIHktYXhpc1xuICAsIHk6IGZ1bmN0aW9uKHkpIHtcbiAgICAgIHZhciBveSA9IHRoaXMuYXR0cigneScpXG4gICAgICAgICwgbyAgPSB0eXBlb2Ygb3kgPT09ICdudW1iZXInID8gb3kgLSB0aGlzLmJib3goKS55IDogMFxuXG4gICAgICAvLyBhY3QgYXMgZ2V0dGVyXG4gICAgICBpZiAoeSA9PSBudWxsKVxuICAgICAgICByZXR1cm4gdHlwZW9mIG95ID09PSAnbnVtYmVyJyA/IG95IC0gbyA6IG95XG5cbiAgICAgIHJldHVybiB0aGlzLmF0dHIoJ3knLCB0eXBlb2YgeSA9PT0gJ251bWJlcicgPyB5ICsgbyA6IHkpXG4gICAgfVxuICAgIC8vIE1vdmUgY2VudGVyIG92ZXIgeC1heGlzXG4gICwgY3g6IGZ1bmN0aW9uKHgpIHtcbiAgICAgIHJldHVybiB4ID09IG51bGwgPyB0aGlzLmJib3goKS5jeCA6IHRoaXMueCh4IC0gdGhpcy5iYm94KCkud2lkdGggLyAyKVxuICAgIH1cbiAgICAvLyBNb3ZlIGNlbnRlciBvdmVyIHktYXhpc1xuICAsIGN5OiBmdW5jdGlvbih5KSB7XG4gICAgICByZXR1cm4geSA9PSBudWxsID8gdGhpcy5iYm94KCkuY3kgOiB0aGlzLnkoeSAtIHRoaXMuYmJveCgpLmhlaWdodCAvIDIpXG4gICAgfVxuICAgIC8vIFNldCB0aGUgdGV4dCBjb250ZW50XG4gICwgdGV4dDogZnVuY3Rpb24odGV4dCkge1xuICAgICAgLy8gYWN0IGFzIGdldHRlclxuICAgICAgaWYgKHR5cGVvZiB0ZXh0ID09PSAndW5kZWZpbmVkJyl7XG4gICAgICAgIHZhciB0ZXh0ID0gJydcbiAgICAgICAgdmFyIGNoaWxkcmVuID0gdGhpcy5ub2RlLmNoaWxkTm9kZXNcbiAgICAgICAgZm9yKHZhciBpID0gMCwgbGVuID0gY2hpbGRyZW4ubGVuZ3RoOyBpIDwgbGVuOyArK2kpe1xuXG4gICAgICAgICAgLy8gYWRkIG5ld2xpbmUgaWYgaXRzIG5vdCB0aGUgZmlyc3QgY2hpbGQgYW5kIG5ld0xpbmVkIGlzIHNldCB0byB0cnVlXG4gICAgICAgICAgaWYoaSAhPSAwICYmIGNoaWxkcmVuW2ldLm5vZGVUeXBlICE9IDMgJiYgU1ZHLmFkb3B0KGNoaWxkcmVuW2ldKS5kb20ubmV3TGluZWQgPT0gdHJ1ZSl7XG4gICAgICAgICAgICB0ZXh0ICs9ICdcXG4nXG4gICAgICAgICAgfVxuXG4gICAgICAgICAgLy8gYWRkIGNvbnRlbnQgb2YgdGhpcyBub2RlXG4gICAgICAgICAgdGV4dCArPSBjaGlsZHJlbltpXS50ZXh0Q29udGVudFxuICAgICAgICB9XG5cbiAgICAgICAgcmV0dXJuIHRleHRcbiAgICAgIH1cblxuICAgICAgLy8gcmVtb3ZlIGV4aXN0aW5nIGNvbnRlbnRcbiAgICAgIHRoaXMuY2xlYXIoKS5idWlsZCh0cnVlKVxuXG4gICAgICBpZiAodHlwZW9mIHRleHQgPT09ICdmdW5jdGlvbicpIHtcbiAgICAgICAgLy8gY2FsbCBibG9ja1xuICAgICAgICB0ZXh0LmNhbGwodGhpcywgdGhpcylcblxuICAgICAgfSBlbHNlIHtcbiAgICAgICAgLy8gc3RvcmUgdGV4dCBhbmQgbWFrZSBzdXJlIHRleHQgaXMgbm90IGJsYW5rXG4gICAgICAgIHRleHQgPSB0ZXh0LnNwbGl0KCdcXG4nKVxuXG4gICAgICAgIC8vIGJ1aWxkIG5ldyBsaW5lc1xuICAgICAgICBmb3IgKHZhciBpID0gMCwgaWwgPSB0ZXh0Lmxlbmd0aDsgaSA8IGlsOyBpKyspXG4gICAgICAgICAgdGhpcy50c3Bhbih0ZXh0W2ldKS5uZXdMaW5lKClcbiAgICAgIH1cblxuICAgICAgLy8gZGlzYWJsZSBidWlsZCBtb2RlIGFuZCByZWJ1aWxkIGxpbmVzXG4gICAgICByZXR1cm4gdGhpcy5idWlsZChmYWxzZSkucmVidWlsZCgpXG4gICAgfVxuICAgIC8vIFNldCBmb250IHNpemVcbiAgLCBzaXplOiBmdW5jdGlvbihzaXplKSB7XG4gICAgICByZXR1cm4gdGhpcy5hdHRyKCdmb250LXNpemUnLCBzaXplKS5yZWJ1aWxkKClcbiAgICB9XG4gICAgLy8gU2V0IC8gZ2V0IGxlYWRpbmdcbiAgLCBsZWFkaW5nOiBmdW5jdGlvbih2YWx1ZSkge1xuICAgICAgLy8gYWN0IGFzIGdldHRlclxuICAgICAgaWYgKHZhbHVlID09IG51bGwpXG4gICAgICAgIHJldHVybiB0aGlzLmRvbS5sZWFkaW5nXG5cbiAgICAgIC8vIGFjdCBhcyBzZXR0ZXJcbiAgICAgIHRoaXMuZG9tLmxlYWRpbmcgPSBuZXcgU1ZHLk51bWJlcih2YWx1ZSlcblxuICAgICAgcmV0dXJuIHRoaXMucmVidWlsZCgpXG4gICAgfVxuICAgIC8vIEdldCBhbGwgdGhlIGZpcnN0IGxldmVsIGxpbmVzXG4gICwgbGluZXM6IGZ1bmN0aW9uKCkge1xuICAgICAgdmFyIG5vZGUgPSAodGhpcy50ZXh0UGF0aCAmJiB0aGlzLnRleHRQYXRoKCkgfHwgdGhpcykubm9kZVxuXG4gICAgICAvLyBmaWx0ZXIgdHNwYW5zIGFuZCBtYXAgdGhlbSB0byBTVkcuanMgaW5zdGFuY2VzXG4gICAgICB2YXIgbGluZXMgPSBTVkcudXRpbHMubWFwKFNWRy51dGlscy5maWx0ZXJTVkdFbGVtZW50cyhub2RlLmNoaWxkTm9kZXMpLCBmdW5jdGlvbihlbCl7XG4gICAgICAgIHJldHVybiBTVkcuYWRvcHQoZWwpXG4gICAgICB9KVxuXG4gICAgICAvLyByZXR1cm4gYW4gaW5zdGFuY2Ugb2YgU1ZHLnNldFxuICAgICAgcmV0dXJuIG5ldyBTVkcuU2V0KGxpbmVzKVxuICAgIH1cbiAgICAvLyBSZWJ1aWxkIGFwcGVhcmFuY2UgdHlwZVxuICAsIHJlYnVpbGQ6IGZ1bmN0aW9uKHJlYnVpbGQpIHtcbiAgICAgIC8vIHN0b3JlIG5ldyByZWJ1aWxkIGZsYWcgaWYgZ2l2ZW5cbiAgICAgIGlmICh0eXBlb2YgcmVidWlsZCA9PSAnYm9vbGVhbicpXG4gICAgICAgIHRoaXMuX3JlYnVpbGQgPSByZWJ1aWxkXG5cbiAgICAgIC8vIGRlZmluZSBwb3NpdGlvbiBvZiBhbGwgbGluZXNcbiAgICAgIGlmICh0aGlzLl9yZWJ1aWxkKSB7XG4gICAgICAgIHZhciBzZWxmID0gdGhpc1xuICAgICAgICAgICwgYmxhbmtMaW5lT2Zmc2V0ID0gMFxuICAgICAgICAgICwgZHkgPSB0aGlzLmRvbS5sZWFkaW5nICogbmV3IFNWRy5OdW1iZXIodGhpcy5hdHRyKCdmb250LXNpemUnKSlcblxuICAgICAgICB0aGlzLmxpbmVzKCkuZWFjaChmdW5jdGlvbigpIHtcbiAgICAgICAgICBpZiAodGhpcy5kb20ubmV3TGluZWQpIHtcbiAgICAgICAgICAgIGlmICghdGhpcy50ZXh0UGF0aClcbiAgICAgICAgICAgICAgdGhpcy5hdHRyKCd4Jywgc2VsZi5hdHRyKCd4JykpXG5cbiAgICAgICAgICAgIGlmKHRoaXMudGV4dCgpID09ICdcXG4nKSB7XG4gICAgICAgICAgICAgIGJsYW5rTGluZU9mZnNldCArPSBkeVxuICAgICAgICAgICAgfWVsc2V7XG4gICAgICAgICAgICAgIHRoaXMuYXR0cignZHknLCBkeSArIGJsYW5rTGluZU9mZnNldClcbiAgICAgICAgICAgICAgYmxhbmtMaW5lT2Zmc2V0ID0gMFxuICAgICAgICAgICAgfVxuICAgICAgICAgIH1cbiAgICAgICAgfSlcblxuICAgICAgICB0aGlzLmZpcmUoJ3JlYnVpbGQnKVxuICAgICAgfVxuXG4gICAgICByZXR1cm4gdGhpc1xuICAgIH1cbiAgICAvLyBFbmFibGUgLyBkaXNhYmxlIGJ1aWxkIG1vZGVcbiAgLCBidWlsZDogZnVuY3Rpb24oYnVpbGQpIHtcbiAgICAgIHRoaXMuX2J1aWxkID0gISFidWlsZFxuICAgICAgcmV0dXJuIHRoaXNcbiAgICB9XG4gICAgLy8gb3ZlcndyaXRlIG1ldGhvZCBmcm9tIHBhcmVudCB0byBzZXQgZGF0YSBwcm9wZXJseVxuICAsIHNldERhdGE6IGZ1bmN0aW9uKG8pe1xuICAgICAgdGhpcy5kb20gPSBvXG4gICAgICB0aGlzLmRvbS5sZWFkaW5nID0gbmV3IFNWRy5OdW1iZXIoby5sZWFkaW5nIHx8IDEuMylcbiAgICAgIHJldHVybiB0aGlzXG4gICAgfVxuICB9XG5cbiAgLy8gQWRkIHBhcmVudCBtZXRob2RcbiwgY29uc3RydWN0OiB7XG4gICAgLy8gQ3JlYXRlIHRleHQgZWxlbWVudFxuICAgIHRleHQ6IGZ1bmN0aW9uKHRleHQpIHtcbiAgICAgIHJldHVybiB0aGlzLnB1dChuZXcgU1ZHLlRleHQpLnRleHQodGV4dClcbiAgICB9XG4gICAgLy8gQ3JlYXRlIHBsYWluIHRleHQgZWxlbWVudFxuICAsIHBsYWluOiBmdW5jdGlvbih0ZXh0KSB7XG4gICAgICByZXR1cm4gdGhpcy5wdXQobmV3IFNWRy5UZXh0KS5wbGFpbih0ZXh0KVxuICAgIH1cbiAgfVxuXG59KVxuXG5TVkcuVHNwYW4gPSBTVkcuaW52ZW50KHtcbiAgLy8gSW5pdGlhbGl6ZSBub2RlXG4gIGNyZWF0ZTogJ3RzcGFuJ1xuXG4gIC8vIEluaGVyaXQgZnJvbVxuLCBpbmhlcml0OiBTVkcuU2hhcGVcblxuICAvLyBBZGQgY2xhc3MgbWV0aG9kc1xuLCBleHRlbmQ6IHtcbiAgICAvLyBTZXQgdGV4dCBjb250ZW50XG4gICAgdGV4dDogZnVuY3Rpb24odGV4dCkge1xuICAgICAgaWYodGV4dCA9PSBudWxsKSByZXR1cm4gdGhpcy5ub2RlLnRleHRDb250ZW50ICsgKHRoaXMuZG9tLm5ld0xpbmVkID8gJ1xcbicgOiAnJylcblxuICAgICAgdHlwZW9mIHRleHQgPT09ICdmdW5jdGlvbicgPyB0ZXh0LmNhbGwodGhpcywgdGhpcykgOiB0aGlzLnBsYWluKHRleHQpXG5cbiAgICAgIHJldHVybiB0aGlzXG4gICAgfVxuICAgIC8vIFNob3J0Y3V0IGR4XG4gICwgZHg6IGZ1bmN0aW9uKGR4KSB7XG4gICAgICByZXR1cm4gdGhpcy5hdHRyKCdkeCcsIGR4KVxuICAgIH1cbiAgICAvLyBTaG9ydGN1dCBkeVxuICAsIGR5OiBmdW5jdGlvbihkeSkge1xuICAgICAgcmV0dXJuIHRoaXMuYXR0cignZHknLCBkeSlcbiAgICB9XG4gICAgLy8gQ3JlYXRlIG5ldyBsaW5lXG4gICwgbmV3TGluZTogZnVuY3Rpb24oKSB7XG4gICAgICAvLyBmZXRjaCB0ZXh0IHBhcmVudFxuICAgICAgdmFyIHQgPSB0aGlzLnBhcmVudChTVkcuVGV4dClcblxuICAgICAgLy8gbWFyayBuZXcgbGluZVxuICAgICAgdGhpcy5kb20ubmV3TGluZWQgPSB0cnVlXG5cbiAgICAgIC8vIGFwcGx5IG5ldyBoecKhblxuICAgICAgcmV0dXJuIHRoaXMuZHkodC5kb20ubGVhZGluZyAqIHQuYXR0cignZm9udC1zaXplJykpLmF0dHIoJ3gnLCB0LngoKSlcbiAgICB9XG4gIH1cblxufSlcblxuU1ZHLmV4dGVuZChTVkcuVGV4dCwgU1ZHLlRzcGFuLCB7XG4gIC8vIENyZWF0ZSBwbGFpbiB0ZXh0IG5vZGVcbiAgcGxhaW46IGZ1bmN0aW9uKHRleHQpIHtcbiAgICAvLyBjbGVhciBpZiBidWlsZCBtb2RlIGlzIGRpc2FibGVkXG4gICAgaWYgKHRoaXMuX2J1aWxkID09PSBmYWxzZSlcbiAgICAgIHRoaXMuY2xlYXIoKVxuXG4gICAgLy8gY3JlYXRlIHRleHQgbm9kZVxuICAgIHRoaXMubm9kZS5hcHBlbmRDaGlsZChkb2N1bWVudC5jcmVhdGVUZXh0Tm9kZSh0ZXh0KSlcblxuICAgIHJldHVybiB0aGlzXG4gIH1cbiAgLy8gQ3JlYXRlIGEgdHNwYW5cbiwgdHNwYW46IGZ1bmN0aW9uKHRleHQpIHtcbiAgICB2YXIgbm9kZSAgPSAodGhpcy50ZXh0UGF0aCAmJiB0aGlzLnRleHRQYXRoKCkgfHwgdGhpcykubm9kZVxuICAgICAgLCB0c3BhbiA9IG5ldyBTVkcuVHNwYW5cblxuICAgIC8vIGNsZWFyIGlmIGJ1aWxkIG1vZGUgaXMgZGlzYWJsZWRcbiAgICBpZiAodGhpcy5fYnVpbGQgPT09IGZhbHNlKVxuICAgICAgdGhpcy5jbGVhcigpXG5cbiAgICAvLyBhZGQgbmV3IHRzcGFuXG4gICAgbm9kZS5hcHBlbmRDaGlsZCh0c3Bhbi5ub2RlKVxuXG4gICAgcmV0dXJuIHRzcGFuLnRleHQodGV4dClcbiAgfVxuICAvLyBDbGVhciBhbGwgbGluZXNcbiwgY2xlYXI6IGZ1bmN0aW9uKCkge1xuICAgIHZhciBub2RlID0gKHRoaXMudGV4dFBhdGggJiYgdGhpcy50ZXh0UGF0aCgpIHx8IHRoaXMpLm5vZGVcblxuICAgIC8vIHJlbW92ZSBleGlzdGluZyBjaGlsZCBub2Rlc1xuICAgIHdoaWxlIChub2RlLmhhc0NoaWxkTm9kZXMoKSlcbiAgICAgIG5vZGUucmVtb3ZlQ2hpbGQobm9kZS5sYXN0Q2hpbGQpXG5cbiAgICByZXR1cm4gdGhpc1xuICB9XG4gIC8vIEdldCBsZW5ndGggb2YgdGV4dCBlbGVtZW50XG4sIGxlbmd0aDogZnVuY3Rpb24oKSB7XG4gICAgcmV0dXJuIHRoaXMubm9kZS5nZXRDb21wdXRlZFRleHRMZW5ndGgoKVxuICB9XG59KVxuXG5TVkcuVGV4dFBhdGggPSBTVkcuaW52ZW50KHtcbiAgLy8gSW5pdGlhbGl6ZSBub2RlXG4gIGNyZWF0ZTogJ3RleHRQYXRoJ1xuXG4gIC8vIEluaGVyaXQgZnJvbVxuLCBpbmhlcml0OiBTVkcuUGFyZW50XG5cbiAgLy8gRGVmaW5lIHBhcmVudCBjbGFzc1xuLCBwYXJlbnQ6IFNWRy5UZXh0XG5cbiAgLy8gQWRkIHBhcmVudCBtZXRob2RcbiwgY29uc3RydWN0OiB7XG4gICAgLy8gQ3JlYXRlIHBhdGggZm9yIHRleHQgdG8gcnVuIG9uXG4gICAgcGF0aDogZnVuY3Rpb24oZCkge1xuICAgICAgLy8gY3JlYXRlIHRleHRQYXRoIGVsZW1lbnRcbiAgICAgIHZhciBwYXRoICA9IG5ldyBTVkcuVGV4dFBhdGhcbiAgICAgICAgLCB0cmFjayA9IHRoaXMuZG9jKCkuZGVmcygpLnBhdGgoZClcblxuICAgICAgLy8gbW92ZSBsaW5lcyB0byB0ZXh0cGF0aFxuICAgICAgd2hpbGUgKHRoaXMubm9kZS5oYXNDaGlsZE5vZGVzKCkpXG4gICAgICAgIHBhdGgubm9kZS5hcHBlbmRDaGlsZCh0aGlzLm5vZGUuZmlyc3RDaGlsZClcblxuICAgICAgLy8gYWRkIHRleHRQYXRoIGVsZW1lbnQgYXMgY2hpbGQgbm9kZVxuICAgICAgdGhpcy5ub2RlLmFwcGVuZENoaWxkKHBhdGgubm9kZSlcblxuICAgICAgLy8gbGluayB0ZXh0UGF0aCB0byBwYXRoIGFuZCBhZGQgY29udGVudFxuICAgICAgcGF0aC5hdHRyKCdocmVmJywgJyMnICsgdHJhY2ssIFNWRy54bGluaylcblxuICAgICAgcmV0dXJuIHRoaXNcbiAgICB9XG4gICAgLy8gUGxvdCBwYXRoIGlmIGFueVxuICAsIHBsb3Q6IGZ1bmN0aW9uKGQpIHtcbiAgICAgIHZhciB0cmFjayA9IHRoaXMudHJhY2soKVxuXG4gICAgICBpZiAodHJhY2spXG4gICAgICAgIHRyYWNrLnBsb3QoZClcblxuICAgICAgcmV0dXJuIHRoaXNcbiAgICB9XG4gICAgLy8gR2V0IHRoZSBwYXRoIHRyYWNrIGVsZW1lbnRcbiAgLCB0cmFjazogZnVuY3Rpb24oKSB7XG4gICAgICB2YXIgcGF0aCA9IHRoaXMudGV4dFBhdGgoKVxuXG4gICAgICBpZiAocGF0aClcbiAgICAgICAgcmV0dXJuIHBhdGgucmVmZXJlbmNlKCdocmVmJylcbiAgICB9XG4gICAgLy8gR2V0IHRoZSB0ZXh0UGF0aCBjaGlsZFxuICAsIHRleHRQYXRoOiBmdW5jdGlvbigpIHtcbiAgICAgIGlmICh0aGlzLm5vZGUuZmlyc3RDaGlsZCAmJiB0aGlzLm5vZGUuZmlyc3RDaGlsZC5ub2RlTmFtZSA9PSAndGV4dFBhdGgnKVxuICAgICAgICByZXR1cm4gU1ZHLmFkb3B0KHRoaXMubm9kZS5maXJzdENoaWxkKVxuICAgIH1cbiAgfVxufSlcblNWRy5OZXN0ZWQgPSBTVkcuaW52ZW50KHtcbiAgLy8gSW5pdGlhbGl6ZSBub2RlXG4gIGNyZWF0ZTogZnVuY3Rpb24oKSB7XG4gICAgdGhpcy5jb25zdHJ1Y3Rvci5jYWxsKHRoaXMsIFNWRy5jcmVhdGUoJ3N2ZycpKVxuXG4gICAgdGhpcy5zdHlsZSgnb3ZlcmZsb3cnLCAndmlzaWJsZScpXG4gIH1cblxuICAvLyBJbmhlcml0IGZyb21cbiwgaW5oZXJpdDogU1ZHLkNvbnRhaW5lclxuXG4gIC8vIEFkZCBwYXJlbnQgbWV0aG9kXG4sIGNvbnN0cnVjdDoge1xuICAgIC8vIENyZWF0ZSBuZXN0ZWQgc3ZnIGRvY3VtZW50XG4gICAgbmVzdGVkOiBmdW5jdGlvbigpIHtcbiAgICAgIHJldHVybiB0aGlzLnB1dChuZXcgU1ZHLk5lc3RlZClcbiAgICB9XG4gIH1cbn0pXG5TVkcuQSA9IFNWRy5pbnZlbnQoe1xuICAvLyBJbml0aWFsaXplIG5vZGVcbiAgY3JlYXRlOiAnYSdcblxuICAvLyBJbmhlcml0IGZyb21cbiwgaW5oZXJpdDogU1ZHLkNvbnRhaW5lclxuXG4gIC8vIEFkZCBjbGFzcyBtZXRob2RzXG4sIGV4dGVuZDoge1xuICAgIC8vIExpbmsgdXJsXG4gICAgdG86IGZ1bmN0aW9uKHVybCkge1xuICAgICAgcmV0dXJuIHRoaXMuYXR0cignaHJlZicsIHVybCwgU1ZHLnhsaW5rKVxuICAgIH1cbiAgICAvLyBMaW5rIHNob3cgYXR0cmlidXRlXG4gICwgc2hvdzogZnVuY3Rpb24odGFyZ2V0KSB7XG4gICAgICByZXR1cm4gdGhpcy5hdHRyKCdzaG93JywgdGFyZ2V0LCBTVkcueGxpbmspXG4gICAgfVxuICAgIC8vIExpbmsgdGFyZ2V0IGF0dHJpYnV0ZVxuICAsIHRhcmdldDogZnVuY3Rpb24odGFyZ2V0KSB7XG4gICAgICByZXR1cm4gdGhpcy5hdHRyKCd0YXJnZXQnLCB0YXJnZXQpXG4gICAgfVxuICB9XG5cbiAgLy8gQWRkIHBhcmVudCBtZXRob2RcbiwgY29uc3RydWN0OiB7XG4gICAgLy8gQ3JlYXRlIGEgaHlwZXJsaW5rIGVsZW1lbnRcbiAgICBsaW5rOiBmdW5jdGlvbih1cmwpIHtcbiAgICAgIHJldHVybiB0aGlzLnB1dChuZXcgU1ZHLkEpLnRvKHVybClcbiAgICB9XG4gIH1cbn0pXG5cblNWRy5leHRlbmQoU1ZHLkVsZW1lbnQsIHtcbiAgLy8gQ3JlYXRlIGEgaHlwZXJsaW5rIGVsZW1lbnRcbiAgbGlua1RvOiBmdW5jdGlvbih1cmwpIHtcbiAgICB2YXIgbGluayA9IG5ldyBTVkcuQVxuXG4gICAgaWYgKHR5cGVvZiB1cmwgPT0gJ2Z1bmN0aW9uJylcbiAgICAgIHVybC5jYWxsKGxpbmssIGxpbmspXG4gICAgZWxzZVxuICAgICAgbGluay50byh1cmwpXG5cbiAgICByZXR1cm4gdGhpcy5wYXJlbnQoKS5wdXQobGluaykucHV0KHRoaXMpXG4gIH1cblxufSlcblNWRy5NYXJrZXIgPSBTVkcuaW52ZW50KHtcbiAgLy8gSW5pdGlhbGl6ZSBub2RlXG4gIGNyZWF0ZTogJ21hcmtlcidcblxuICAvLyBJbmhlcml0IGZyb21cbiwgaW5oZXJpdDogU1ZHLkNvbnRhaW5lclxuXG4gIC8vIEFkZCBjbGFzcyBtZXRob2RzXG4sIGV4dGVuZDoge1xuICAgIC8vIFNldCB3aWR0aCBvZiBlbGVtZW50XG4gICAgd2lkdGg6IGZ1bmN0aW9uKHdpZHRoKSB7XG4gICAgICByZXR1cm4gdGhpcy5hdHRyKCdtYXJrZXJXaWR0aCcsIHdpZHRoKVxuICAgIH1cbiAgICAvLyBTZXQgaGVpZ2h0IG9mIGVsZW1lbnRcbiAgLCBoZWlnaHQ6IGZ1bmN0aW9uKGhlaWdodCkge1xuICAgICAgcmV0dXJuIHRoaXMuYXR0cignbWFya2VySGVpZ2h0JywgaGVpZ2h0KVxuICAgIH1cbiAgICAvLyBTZXQgbWFya2VyIHJlZlggYW5kIHJlZllcbiAgLCByZWY6IGZ1bmN0aW9uKHgsIHkpIHtcbiAgICAgIHJldHVybiB0aGlzLmF0dHIoJ3JlZlgnLCB4KS5hdHRyKCdyZWZZJywgeSlcbiAgICB9XG4gICAgLy8gVXBkYXRlIG1hcmtlclxuICAsIHVwZGF0ZTogZnVuY3Rpb24oYmxvY2spIHtcbiAgICAgIC8vIHJlbW92ZSBhbGwgY29udGVudFxuICAgICAgdGhpcy5jbGVhcigpXG5cbiAgICAgIC8vIGludm9rZSBwYXNzZWQgYmxvY2tcbiAgICAgIGlmICh0eXBlb2YgYmxvY2sgPT0gJ2Z1bmN0aW9uJylcbiAgICAgICAgYmxvY2suY2FsbCh0aGlzLCB0aGlzKVxuXG4gICAgICByZXR1cm4gdGhpc1xuICAgIH1cbiAgICAvLyBSZXR1cm4gdGhlIGZpbGwgaWRcbiAgLCB0b1N0cmluZzogZnVuY3Rpb24oKSB7XG4gICAgICByZXR1cm4gJ3VybCgjJyArIHRoaXMuaWQoKSArICcpJ1xuICAgIH1cbiAgfVxuXG4gIC8vIEFkZCBwYXJlbnQgbWV0aG9kXG4sIGNvbnN0cnVjdDoge1xuICAgIG1hcmtlcjogZnVuY3Rpb24od2lkdGgsIGhlaWdodCwgYmxvY2spIHtcbiAgICAgIC8vIENyZWF0ZSBtYXJrZXIgZWxlbWVudCBpbiBkZWZzXG4gICAgICByZXR1cm4gdGhpcy5kZWZzKCkubWFya2VyKHdpZHRoLCBoZWlnaHQsIGJsb2NrKVxuICAgIH1cbiAgfVxuXG59KVxuXG5TVkcuZXh0ZW5kKFNWRy5EZWZzLCB7XG4gIC8vIENyZWF0ZSBtYXJrZXJcbiAgbWFya2VyOiBmdW5jdGlvbih3aWR0aCwgaGVpZ2h0LCBibG9jaykge1xuICAgIC8vIFNldCBkZWZhdWx0IHZpZXdib3ggdG8gbWF0Y2ggdGhlIHdpZHRoIGFuZCBoZWlnaHQsIHNldCByZWYgdG8gY3ggYW5kIGN5IGFuZCBzZXQgb3JpZW50IHRvIGF1dG9cbiAgICByZXR1cm4gdGhpcy5wdXQobmV3IFNWRy5NYXJrZXIpXG4gICAgICAuc2l6ZSh3aWR0aCwgaGVpZ2h0KVxuICAgICAgLnJlZih3aWR0aCAvIDIsIGhlaWdodCAvIDIpXG4gICAgICAudmlld2JveCgwLCAwLCB3aWR0aCwgaGVpZ2h0KVxuICAgICAgLmF0dHIoJ29yaWVudCcsICdhdXRvJylcbiAgICAgIC51cGRhdGUoYmxvY2spXG4gIH1cblxufSlcblxuU1ZHLmV4dGVuZChTVkcuTGluZSwgU1ZHLlBvbHlsaW5lLCBTVkcuUG9seWdvbiwgU1ZHLlBhdGgsIHtcbiAgLy8gQ3JlYXRlIGFuZCBhdHRhY2ggbWFya2Vyc1xuICBtYXJrZXI6IGZ1bmN0aW9uKG1hcmtlciwgd2lkdGgsIGhlaWdodCwgYmxvY2spIHtcbiAgICB2YXIgYXR0ciA9IFsnbWFya2VyJ11cblxuICAgIC8vIEJ1aWxkIGF0dHJpYnV0ZSBuYW1lXG4gICAgaWYgKG1hcmtlciAhPSAnYWxsJykgYXR0ci5wdXNoKG1hcmtlcilcbiAgICBhdHRyID0gYXR0ci5qb2luKCctJylcblxuICAgIC8vIFNldCBtYXJrZXIgYXR0cmlidXRlXG4gICAgbWFya2VyID0gYXJndW1lbnRzWzFdIGluc3RhbmNlb2YgU1ZHLk1hcmtlciA/XG4gICAgICBhcmd1bWVudHNbMV0gOlxuICAgICAgdGhpcy5kb2MoKS5tYXJrZXIod2lkdGgsIGhlaWdodCwgYmxvY2spXG5cbiAgICByZXR1cm4gdGhpcy5hdHRyKGF0dHIsIG1hcmtlcilcbiAgfVxuXG59KVxuLy8gRGVmaW5lIGxpc3Qgb2YgYXZhaWxhYmxlIGF0dHJpYnV0ZXMgZm9yIHN0cm9rZSBhbmQgZmlsbFxudmFyIHN1Z2FyID0ge1xuICBzdHJva2U6IFsnY29sb3InLCAnd2lkdGgnLCAnb3BhY2l0eScsICdsaW5lY2FwJywgJ2xpbmVqb2luJywgJ21pdGVybGltaXQnLCAnZGFzaGFycmF5JywgJ2Rhc2hvZmZzZXQnXVxuLCBmaWxsOiAgIFsnY29sb3InLCAnb3BhY2l0eScsICdydWxlJ11cbiwgcHJlZml4OiBmdW5jdGlvbih0LCBhKSB7XG4gICAgcmV0dXJuIGEgPT0gJ2NvbG9yJyA/IHQgOiB0ICsgJy0nICsgYVxuICB9XG59XG5cbi8vIEFkZCBzdWdhciBmb3IgZmlsbCBhbmQgc3Ryb2tlXG47WydmaWxsJywgJ3N0cm9rZSddLmZvckVhY2goZnVuY3Rpb24obSkge1xuICB2YXIgaSwgZXh0ZW5zaW9uID0ge31cblxuICBleHRlbnNpb25bbV0gPSBmdW5jdGlvbihvKSB7XG4gICAgaWYgKHR5cGVvZiBvID09ICd1bmRlZmluZWQnKVxuICAgICAgcmV0dXJuIHRoaXNcbiAgICBpZiAodHlwZW9mIG8gPT0gJ3N0cmluZycgfHwgU1ZHLkNvbG9yLmlzUmdiKG8pIHx8IChvICYmIHR5cGVvZiBvLmZpbGwgPT09ICdmdW5jdGlvbicpKVxuICAgICAgdGhpcy5hdHRyKG0sIG8pXG5cbiAgICBlbHNlXG4gICAgICAvLyBzZXQgYWxsIGF0dHJpYnV0ZXMgZnJvbSBzdWdhci5maWxsIGFuZCBzdWdhci5zdHJva2UgbGlzdFxuICAgICAgZm9yIChpID0gc3VnYXJbbV0ubGVuZ3RoIC0gMTsgaSA+PSAwOyBpLS0pXG4gICAgICAgIGlmIChvW3N1Z2FyW21dW2ldXSAhPSBudWxsKVxuICAgICAgICAgIHRoaXMuYXR0cihzdWdhci5wcmVmaXgobSwgc3VnYXJbbV1baV0pLCBvW3N1Z2FyW21dW2ldXSlcblxuICAgIHJldHVybiB0aGlzXG4gIH1cblxuICBTVkcuZXh0ZW5kKFNWRy5FbGVtZW50LCBTVkcuRlgsIGV4dGVuc2lvbilcblxufSlcblxuU1ZHLmV4dGVuZChTVkcuRWxlbWVudCwgU1ZHLkZYLCB7XG4gIC8vIE1hcCByb3RhdGlvbiB0byB0cmFuc2Zvcm1cbiAgcm90YXRlOiBmdW5jdGlvbihkLCBjeCwgY3kpIHtcbiAgICByZXR1cm4gdGhpcy50cmFuc2Zvcm0oeyByb3RhdGlvbjogZCwgY3g6IGN4LCBjeTogY3kgfSlcbiAgfVxuICAvLyBNYXAgc2tldyB0byB0cmFuc2Zvcm1cbiwgc2tldzogZnVuY3Rpb24oeCwgeSwgY3gsIGN5KSB7XG4gICAgcmV0dXJuIGFyZ3VtZW50cy5sZW5ndGggPT0gMSAgfHwgYXJndW1lbnRzLmxlbmd0aCA9PSAzID9cbiAgICAgIHRoaXMudHJhbnNmb3JtKHsgc2tldzogeCwgY3g6IHksIGN5OiBjeCB9KSA6XG4gICAgICB0aGlzLnRyYW5zZm9ybSh7IHNrZXdYOiB4LCBza2V3WTogeSwgY3g6IGN4LCBjeTogY3kgfSlcbiAgfVxuICAvLyBNYXAgc2NhbGUgdG8gdHJhbnNmb3JtXG4sIHNjYWxlOiBmdW5jdGlvbih4LCB5LCBjeCwgY3kpIHtcbiAgICByZXR1cm4gYXJndW1lbnRzLmxlbmd0aCA9PSAxICB8fCBhcmd1bWVudHMubGVuZ3RoID09IDMgP1xuICAgICAgdGhpcy50cmFuc2Zvcm0oeyBzY2FsZTogeCwgY3g6IHksIGN5OiBjeCB9KSA6XG4gICAgICB0aGlzLnRyYW5zZm9ybSh7IHNjYWxlWDogeCwgc2NhbGVZOiB5LCBjeDogY3gsIGN5OiBjeSB9KVxuICB9XG4gIC8vIE1hcCB0cmFuc2xhdGUgdG8gdHJhbnNmb3JtXG4sIHRyYW5zbGF0ZTogZnVuY3Rpb24oeCwgeSkge1xuICAgIHJldHVybiB0aGlzLnRyYW5zZm9ybSh7IHg6IHgsIHk6IHkgfSlcbiAgfVxuICAvLyBNYXAgZmxpcCB0byB0cmFuc2Zvcm1cbiwgZmxpcDogZnVuY3Rpb24oYSwgbykge1xuICAgIHJldHVybiB0aGlzLnRyYW5zZm9ybSh7IGZsaXA6IGEsIG9mZnNldDogbyB9KVxuICB9XG4gIC8vIE1hcCBtYXRyaXggdG8gdHJhbnNmb3JtXG4sIG1hdHJpeDogZnVuY3Rpb24obSkge1xuICAgIHJldHVybiB0aGlzLmF0dHIoJ3RyYW5zZm9ybScsIG5ldyBTVkcuTWF0cml4KG0pKVxuICB9XG4gIC8vIE9wYWNpdHlcbiwgb3BhY2l0eTogZnVuY3Rpb24odmFsdWUpIHtcbiAgICByZXR1cm4gdGhpcy5hdHRyKCdvcGFjaXR5JywgdmFsdWUpXG4gIH1cbiAgLy8gUmVsYXRpdmUgbW92ZSBvdmVyIHggYXhpc1xuLCBkeDogZnVuY3Rpb24oeCkge1xuICAgIHJldHVybiB0aGlzLngoKHRoaXMgaW5zdGFuY2VvZiBTVkcuRlggPyAwIDogdGhpcy54KCkpICsgeCwgdHJ1ZSlcbiAgfVxuICAvLyBSZWxhdGl2ZSBtb3ZlIG92ZXIgeSBheGlzXG4sIGR5OiBmdW5jdGlvbih5KSB7XG4gICAgcmV0dXJuIHRoaXMueSgodGhpcyBpbnN0YW5jZW9mIFNWRy5GWCA/IDAgOiB0aGlzLnkoKSkgKyB5LCB0cnVlKVxuICB9XG4gIC8vIFJlbGF0aXZlIG1vdmUgb3ZlciB4IGFuZCB5IGF4ZXNcbiwgZG1vdmU6IGZ1bmN0aW9uKHgsIHkpIHtcbiAgICByZXR1cm4gdGhpcy5keCh4KS5keSh5KVxuICB9XG59KVxuXG5TVkcuZXh0ZW5kKFNWRy5SZWN0LCBTVkcuRWxsaXBzZSwgU1ZHLkNpcmNsZSwgU1ZHLkdyYWRpZW50LCBTVkcuRlgsIHtcbiAgLy8gQWRkIHggYW5kIHkgcmFkaXVzXG4gIHJhZGl1czogZnVuY3Rpb24oeCwgeSkge1xuICAgIHZhciB0eXBlID0gKHRoaXMuX3RhcmdldCB8fCB0aGlzKS50eXBlO1xuICAgIHJldHVybiB0eXBlID09ICdyYWRpYWwnIHx8IHR5cGUgPT0gJ2NpcmNsZScgP1xuICAgICAgdGhpcy5hdHRyKCdyJywgbmV3IFNWRy5OdW1iZXIoeCkpIDpcbiAgICAgIHRoaXMucngoeCkucnkoeSA9PSBudWxsID8geCA6IHkpXG4gIH1cbn0pXG5cblNWRy5leHRlbmQoU1ZHLlBhdGgsIHtcbiAgLy8gR2V0IHBhdGggbGVuZ3RoXG4gIGxlbmd0aDogZnVuY3Rpb24oKSB7XG4gICAgcmV0dXJuIHRoaXMubm9kZS5nZXRUb3RhbExlbmd0aCgpXG4gIH1cbiAgLy8gR2V0IHBvaW50IGF0IGxlbmd0aFxuLCBwb2ludEF0OiBmdW5jdGlvbihsZW5ndGgpIHtcbiAgICByZXR1cm4gdGhpcy5ub2RlLmdldFBvaW50QXRMZW5ndGgobGVuZ3RoKVxuICB9XG59KVxuXG5TVkcuZXh0ZW5kKFNWRy5QYXJlbnQsIFNWRy5UZXh0LCBTVkcuRlgsIHtcbiAgLy8gU2V0IGZvbnRcbiAgZm9udDogZnVuY3Rpb24obykge1xuICAgIGZvciAodmFyIGsgaW4gbylcbiAgICAgIGsgPT0gJ2xlYWRpbmcnID9cbiAgICAgICAgdGhpcy5sZWFkaW5nKG9ba10pIDpcbiAgICAgIGsgPT0gJ2FuY2hvcicgP1xuICAgICAgICB0aGlzLmF0dHIoJ3RleHQtYW5jaG9yJywgb1trXSkgOlxuICAgICAgayA9PSAnc2l6ZScgfHwgayA9PSAnZmFtaWx5JyB8fCBrID09ICd3ZWlnaHQnIHx8IGsgPT0gJ3N0cmV0Y2gnIHx8IGsgPT0gJ3ZhcmlhbnQnIHx8IGsgPT0gJ3N0eWxlJyA/XG4gICAgICAgIHRoaXMuYXR0cignZm9udC0nKyBrLCBvW2tdKSA6XG4gICAgICAgIHRoaXMuYXR0cihrLCBvW2tdKVxuXG4gICAgcmV0dXJuIHRoaXNcbiAgfVxufSlcblxuU1ZHLlNldCA9IFNWRy5pbnZlbnQoe1xuICAvLyBJbml0aWFsaXplXG4gIGNyZWF0ZTogZnVuY3Rpb24obWVtYmVycykge1xuICAgIC8vIFNldCBpbml0aWFsIHN0YXRlXG4gICAgQXJyYXkuaXNBcnJheShtZW1iZXJzKSA/IHRoaXMubWVtYmVycyA9IG1lbWJlcnMgOiB0aGlzLmNsZWFyKClcbiAgfVxuXG4gIC8vIEFkZCBjbGFzcyBtZXRob2RzXG4sIGV4dGVuZDoge1xuICAgIC8vIEFkZCBlbGVtZW50IHRvIHNldFxuICAgIGFkZDogZnVuY3Rpb24oKSB7XG4gICAgICB2YXIgaSwgaWwsIGVsZW1lbnRzID0gW10uc2xpY2UuY2FsbChhcmd1bWVudHMpXG5cbiAgICAgIGZvciAoaSA9IDAsIGlsID0gZWxlbWVudHMubGVuZ3RoOyBpIDwgaWw7IGkrKylcbiAgICAgICAgdGhpcy5tZW1iZXJzLnB1c2goZWxlbWVudHNbaV0pXG5cbiAgICAgIHJldHVybiB0aGlzXG4gICAgfVxuICAgIC8vIFJlbW92ZSBlbGVtZW50IGZyb20gc2V0XG4gICwgcmVtb3ZlOiBmdW5jdGlvbihlbGVtZW50KSB7XG4gICAgICB2YXIgaSA9IHRoaXMuaW5kZXgoZWxlbWVudClcblxuICAgICAgLy8gcmVtb3ZlIGdpdmVuIGNoaWxkXG4gICAgICBpZiAoaSA+IC0xKVxuICAgICAgICB0aGlzLm1lbWJlcnMuc3BsaWNlKGksIDEpXG5cbiAgICAgIHJldHVybiB0aGlzXG4gICAgfVxuICAgIC8vIEl0ZXJhdGUgb3ZlciBhbGwgbWVtYmVyc1xuICAsIGVhY2g6IGZ1bmN0aW9uKGJsb2NrKSB7XG4gICAgICBmb3IgKHZhciBpID0gMCwgaWwgPSB0aGlzLm1lbWJlcnMubGVuZ3RoOyBpIDwgaWw7IGkrKylcbiAgICAgICAgYmxvY2suYXBwbHkodGhpcy5tZW1iZXJzW2ldLCBbaSwgdGhpcy5tZW1iZXJzXSlcblxuICAgICAgcmV0dXJuIHRoaXNcbiAgICB9XG4gICAgLy8gUmVzdG9yZSB0byBkZWZhdWx0c1xuICAsIGNsZWFyOiBmdW5jdGlvbigpIHtcbiAgICAgIC8vIGluaXRpYWxpemUgc3RvcmVcbiAgICAgIHRoaXMubWVtYmVycyA9IFtdXG5cbiAgICAgIHJldHVybiB0aGlzXG4gICAgfVxuICAgIC8vIEdldCB0aGUgbGVuZ3RoIG9mIGEgc2V0XG4gICwgbGVuZ3RoOiBmdW5jdGlvbigpIHtcbiAgICAgIHJldHVybiB0aGlzLm1lbWJlcnMubGVuZ3RoXG4gICAgfVxuICAgIC8vIENoZWNrcyBpZiBhIGdpdmVuIGVsZW1lbnQgaXMgcHJlc2VudCBpbiBzZXRcbiAgLCBoYXM6IGZ1bmN0aW9uKGVsZW1lbnQpIHtcbiAgICAgIHJldHVybiB0aGlzLmluZGV4KGVsZW1lbnQpID49IDBcbiAgICB9XG4gICAgLy8gcmV0dW5zIGluZGV4IG9mIGdpdmVuIGVsZW1lbnQgaW4gc2V0XG4gICwgaW5kZXg6IGZ1bmN0aW9uKGVsZW1lbnQpIHtcbiAgICAgIHJldHVybiB0aGlzLm1lbWJlcnMuaW5kZXhPZihlbGVtZW50KVxuICAgIH1cbiAgICAvLyBHZXQgbWVtYmVyIGF0IGdpdmVuIGluZGV4XG4gICwgZ2V0OiBmdW5jdGlvbihpKSB7XG4gICAgICByZXR1cm4gdGhpcy5tZW1iZXJzW2ldXG4gICAgfVxuICAgIC8vIEdldCBmaXJzdCBtZW1iZXJcbiAgLCBmaXJzdDogZnVuY3Rpb24oKSB7XG4gICAgICByZXR1cm4gdGhpcy5nZXQoMClcbiAgICB9XG4gICAgLy8gR2V0IGxhc3QgbWVtYmVyXG4gICwgbGFzdDogZnVuY3Rpb24oKSB7XG4gICAgICByZXR1cm4gdGhpcy5nZXQodGhpcy5tZW1iZXJzLmxlbmd0aCAtIDEpXG4gICAgfVxuICAgIC8vIERlZmF1bHQgdmFsdWVcbiAgLCB2YWx1ZU9mOiBmdW5jdGlvbigpIHtcbiAgICAgIHJldHVybiB0aGlzLm1lbWJlcnNcbiAgICB9XG4gICAgLy8gR2V0IHRoZSBib3VuZGluZyBib3ggb2YgYWxsIG1lbWJlcnMgaW5jbHVkZWQgb3IgZW1wdHkgYm94IGlmIHNldCBoYXMgbm8gaXRlbXNcbiAgLCBiYm94OiBmdW5jdGlvbigpe1xuICAgICAgdmFyIGJveCA9IG5ldyBTVkcuQkJveCgpXG5cbiAgICAgIC8vIHJldHVybiBhbiBlbXB0eSBib3ggb2YgdGhlcmUgYXJlIG5vIG1lbWJlcnNcbiAgICAgIGlmICh0aGlzLm1lbWJlcnMubGVuZ3RoID09IDApXG4gICAgICAgIHJldHVybiBib3hcblxuICAgICAgLy8gZ2V0IHRoZSBmaXJzdCByYm94IGFuZCB1cGRhdGUgdGhlIHRhcmdldCBiYm94XG4gICAgICB2YXIgcmJveCA9IHRoaXMubWVtYmVyc1swXS5yYm94KClcbiAgICAgIGJveC54ICAgICAgPSByYm94LnhcbiAgICAgIGJveC55ICAgICAgPSByYm94LnlcbiAgICAgIGJveC53aWR0aCAgPSByYm94LndpZHRoXG4gICAgICBib3guaGVpZ2h0ID0gcmJveC5oZWlnaHRcblxuICAgICAgdGhpcy5lYWNoKGZ1bmN0aW9uKCkge1xuICAgICAgICAvLyB1c2VyIHJib3ggZm9yIGNvcnJlY3QgcG9zaXRpb24gYW5kIHZpc3VhbCByZXByZXNlbnRhdGlvblxuICAgICAgICBib3ggPSBib3gubWVyZ2UodGhpcy5yYm94KCkpXG4gICAgICB9KVxuXG4gICAgICByZXR1cm4gYm94XG4gICAgfVxuICB9XG5cbiAgLy8gQWRkIHBhcmVudCBtZXRob2RcbiwgY29uc3RydWN0OiB7XG4gICAgLy8gQ3JlYXRlIGEgbmV3IHNldFxuICAgIHNldDogZnVuY3Rpb24obWVtYmVycykge1xuICAgICAgcmV0dXJuIG5ldyBTVkcuU2V0KG1lbWJlcnMpXG4gICAgfVxuICB9XG59KVxuXG5TVkcuRlguU2V0ID0gU1ZHLmludmVudCh7XG4gIC8vIEluaXRpYWxpemUgbm9kZVxuICBjcmVhdGU6IGZ1bmN0aW9uKHNldCkge1xuICAgIC8vIHN0b3JlIHJlZmVyZW5jZSB0byBzZXRcbiAgICB0aGlzLnNldCA9IHNldFxuICB9XG5cbn0pXG5cbi8vIEFsaWFzIG1ldGhvZHNcblNWRy5TZXQuaW5oZXJpdCA9IGZ1bmN0aW9uKCkge1xuICB2YXIgbVxuICAgICwgbWV0aG9kcyA9IFtdXG5cbiAgLy8gZ2F0aGVyIHNoYXBlIG1ldGhvZHNcbiAgZm9yKHZhciBtIGluIFNWRy5TaGFwZS5wcm90b3R5cGUpXG4gICAgaWYgKHR5cGVvZiBTVkcuU2hhcGUucHJvdG90eXBlW21dID09ICdmdW5jdGlvbicgJiYgdHlwZW9mIFNWRy5TZXQucHJvdG90eXBlW21dICE9ICdmdW5jdGlvbicpXG4gICAgICBtZXRob2RzLnB1c2gobSlcblxuICAvLyBhcHBseSBzaGFwZSBhbGlhc3Nlc1xuICBtZXRob2RzLmZvckVhY2goZnVuY3Rpb24obWV0aG9kKSB7XG4gICAgU1ZHLlNldC5wcm90b3R5cGVbbWV0aG9kXSA9IGZ1bmN0aW9uKCkge1xuICAgICAgZm9yICh2YXIgaSA9IDAsIGlsID0gdGhpcy5tZW1iZXJzLmxlbmd0aDsgaSA8IGlsOyBpKyspXG4gICAgICAgIGlmICh0aGlzLm1lbWJlcnNbaV0gJiYgdHlwZW9mIHRoaXMubWVtYmVyc1tpXVttZXRob2RdID09ICdmdW5jdGlvbicpXG4gICAgICAgICAgdGhpcy5tZW1iZXJzW2ldW21ldGhvZF0uYXBwbHkodGhpcy5tZW1iZXJzW2ldLCBhcmd1bWVudHMpXG5cbiAgICAgIHJldHVybiBtZXRob2QgPT0gJ2FuaW1hdGUnID8gKHRoaXMuZnggfHwgKHRoaXMuZnggPSBuZXcgU1ZHLkZYLlNldCh0aGlzKSkpIDogdGhpc1xuICAgIH1cbiAgfSlcblxuICAvLyBjbGVhciBtZXRob2RzIGZvciB0aGUgbmV4dCByb3VuZFxuICBtZXRob2RzID0gW11cblxuICAvLyBnYXRoZXIgZnggbWV0aG9kc1xuICBmb3IodmFyIG0gaW4gU1ZHLkZYLnByb3RvdHlwZSlcbiAgICBpZiAodHlwZW9mIFNWRy5GWC5wcm90b3R5cGVbbV0gPT0gJ2Z1bmN0aW9uJyAmJiB0eXBlb2YgU1ZHLkZYLlNldC5wcm90b3R5cGVbbV0gIT0gJ2Z1bmN0aW9uJylcbiAgICAgIG1ldGhvZHMucHVzaChtKVxuXG4gIC8vIGFwcGx5IGZ4IGFsaWFzc2VzXG4gIG1ldGhvZHMuZm9yRWFjaChmdW5jdGlvbihtZXRob2QpIHtcbiAgICBTVkcuRlguU2V0LnByb3RvdHlwZVttZXRob2RdID0gZnVuY3Rpb24oKSB7XG4gICAgICBmb3IgKHZhciBpID0gMCwgaWwgPSB0aGlzLnNldC5tZW1iZXJzLmxlbmd0aDsgaSA8IGlsOyBpKyspXG4gICAgICAgIHRoaXMuc2V0Lm1lbWJlcnNbaV0uZnhbbWV0aG9kXS5hcHBseSh0aGlzLnNldC5tZW1iZXJzW2ldLmZ4LCBhcmd1bWVudHMpXG5cbiAgICAgIHJldHVybiB0aGlzXG4gICAgfVxuICB9KVxufVxuXG5cblxuXG5TVkcuZXh0ZW5kKFNWRy5FbGVtZW50LCB7XG4gIC8vIFN0b3JlIGRhdGEgdmFsdWVzIG9uIHN2ZyBub2Rlc1xuICBkYXRhOiBmdW5jdGlvbihhLCB2LCByKSB7XG4gICAgaWYgKHR5cGVvZiBhID09ICdvYmplY3QnKSB7XG4gICAgICBmb3IgKHYgaW4gYSlcbiAgICAgICAgdGhpcy5kYXRhKHYsIGFbdl0pXG5cbiAgICB9IGVsc2UgaWYgKGFyZ3VtZW50cy5sZW5ndGggPCAyKSB7XG4gICAgICB0cnkge1xuICAgICAgICByZXR1cm4gSlNPTi5wYXJzZSh0aGlzLmF0dHIoJ2RhdGEtJyArIGEpKVxuICAgICAgfSBjYXRjaChlKSB7XG4gICAgICAgIHJldHVybiB0aGlzLmF0dHIoJ2RhdGEtJyArIGEpXG4gICAgICB9XG5cbiAgICB9IGVsc2Uge1xuICAgICAgdGhpcy5hdHRyKFxuICAgICAgICAnZGF0YS0nICsgYVxuICAgICAgLCB2ID09PSBudWxsID9cbiAgICAgICAgICBudWxsIDpcbiAgICAgICAgciA9PT0gdHJ1ZSB8fCB0eXBlb2YgdiA9PT0gJ3N0cmluZycgfHwgdHlwZW9mIHYgPT09ICdudW1iZXInID9cbiAgICAgICAgICB2IDpcbiAgICAgICAgICBKU09OLnN0cmluZ2lmeSh2KVxuICAgICAgKVxuICAgIH1cblxuICAgIHJldHVybiB0aGlzXG4gIH1cbn0pXG5TVkcuZXh0ZW5kKFNWRy5FbGVtZW50LCB7XG4gIC8vIFJlbWVtYmVyIGFyYml0cmFyeSBkYXRhXG4gIHJlbWVtYmVyOiBmdW5jdGlvbihrLCB2KSB7XG4gICAgLy8gcmVtZW1iZXIgZXZlcnkgaXRlbSBpbiBhbiBvYmplY3QgaW5kaXZpZHVhbGx5XG4gICAgaWYgKHR5cGVvZiBhcmd1bWVudHNbMF0gPT0gJ29iamVjdCcpXG4gICAgICBmb3IgKHZhciB2IGluIGspXG4gICAgICAgIHRoaXMucmVtZW1iZXIodiwga1t2XSlcblxuICAgIC8vIHJldHJpZXZlIG1lbW9yeVxuICAgIGVsc2UgaWYgKGFyZ3VtZW50cy5sZW5ndGggPT0gMSlcbiAgICAgIHJldHVybiB0aGlzLm1lbW9yeSgpW2tdXG5cbiAgICAvLyBzdG9yZSBtZW1vcnlcbiAgICBlbHNlXG4gICAgICB0aGlzLm1lbW9yeSgpW2tdID0gdlxuXG4gICAgcmV0dXJuIHRoaXNcbiAgfVxuXG4gIC8vIEVyYXNlIGEgZ2l2ZW4gbWVtb3J5XG4sIGZvcmdldDogZnVuY3Rpb24oKSB7XG4gICAgaWYgKGFyZ3VtZW50cy5sZW5ndGggPT0gMClcbiAgICAgIHRoaXMuX21lbW9yeSA9IHt9XG4gICAgZWxzZVxuICAgICAgZm9yICh2YXIgaSA9IGFyZ3VtZW50cy5sZW5ndGggLSAxOyBpID49IDA7IGktLSlcbiAgICAgICAgZGVsZXRlIHRoaXMubWVtb3J5KClbYXJndW1lbnRzW2ldXVxuXG4gICAgcmV0dXJuIHRoaXNcbiAgfVxuXG4gIC8vIEluaXRpYWxpemUgb3IgcmV0dXJuIGxvY2FsIG1lbW9yeSBvYmplY3RcbiwgbWVtb3J5OiBmdW5jdGlvbigpIHtcbiAgICByZXR1cm4gdGhpcy5fbWVtb3J5IHx8ICh0aGlzLl9tZW1vcnkgPSB7fSlcbiAgfVxuXG59KVxuLy8gTWV0aG9kIGZvciBnZXR0aW5nIGFuIGVsZW1lbnQgYnkgaWRcblNWRy5nZXQgPSBmdW5jdGlvbihpZCkge1xuICB2YXIgbm9kZSA9IGRvY3VtZW50LmdldEVsZW1lbnRCeUlkKGlkRnJvbVJlZmVyZW5jZShpZCkgfHwgaWQpXG4gIHJldHVybiBTVkcuYWRvcHQobm9kZSlcbn1cblxuLy8gU2VsZWN0IGVsZW1lbnRzIGJ5IHF1ZXJ5IHN0cmluZ1xuU1ZHLnNlbGVjdCA9IGZ1bmN0aW9uKHF1ZXJ5LCBwYXJlbnQpIHtcbiAgcmV0dXJuIG5ldyBTVkcuU2V0KFxuICAgIFNWRy51dGlscy5tYXAoKHBhcmVudCB8fCBkb2N1bWVudCkucXVlcnlTZWxlY3RvckFsbChxdWVyeSksIGZ1bmN0aW9uKG5vZGUpIHtcbiAgICAgIHJldHVybiBTVkcuYWRvcHQobm9kZSlcbiAgICB9KVxuICApXG59XG5cblNWRy5leHRlbmQoU1ZHLlBhcmVudCwge1xuICAvLyBTY29wZWQgc2VsZWN0IG1ldGhvZFxuICBzZWxlY3Q6IGZ1bmN0aW9uKHF1ZXJ5KSB7XG4gICAgcmV0dXJuIFNWRy5zZWxlY3QocXVlcnksIHRoaXMubm9kZSlcbiAgfVxuXG59KVxuZnVuY3Rpb24gaXMoZWwsIG9iail7XG4gIHJldHVybiBlbCBpbnN0YW5jZW9mIG9ialxufVxuXG4vLyB0ZXN0cyBpZiBhIGdpdmVuIHNlbGVjdG9yIG1hdGNoZXMgYW4gZWxlbWVudFxuZnVuY3Rpb24gbWF0Y2hlcyhlbCwgc2VsZWN0b3IpIHtcbiAgcmV0dXJuIChlbC5tYXRjaGVzIHx8IGVsLm1hdGNoZXNTZWxlY3RvciB8fCBlbC5tc01hdGNoZXNTZWxlY3RvciB8fCBlbC5tb3pNYXRjaGVzU2VsZWN0b3IgfHwgZWwud2Via2l0TWF0Y2hlc1NlbGVjdG9yIHx8IGVsLm9NYXRjaGVzU2VsZWN0b3IpLmNhbGwoZWwsIHNlbGVjdG9yKTtcbn1cblxuLy8gQ29udmVydCBkYXNoLXNlcGFyYXRlZC1zdHJpbmcgdG8gY2FtZWxDYXNlXG5mdW5jdGlvbiBjYW1lbENhc2Uocykge1xuICByZXR1cm4gcy50b0xvd2VyQ2FzZSgpLnJlcGxhY2UoLy0oLikvZywgZnVuY3Rpb24obSwgZykge1xuICAgIHJldHVybiBnLnRvVXBwZXJDYXNlKClcbiAgfSlcbn1cblxuLy8gQ2FwaXRhbGl6ZSBmaXJzdCBsZXR0ZXIgb2YgYSBzdHJpbmdcbmZ1bmN0aW9uIGNhcGl0YWxpemUocykge1xuICByZXR1cm4gcy5jaGFyQXQoMCkudG9VcHBlckNhc2UoKSArIHMuc2xpY2UoMSlcbn1cblxuLy8gRW5zdXJlIHRvIHNpeC1iYXNlZCBoZXhcbmZ1bmN0aW9uIGZ1bGxIZXgoaGV4KSB7XG4gIHJldHVybiBoZXgubGVuZ3RoID09IDQgP1xuICAgIFsgJyMnLFxuICAgICAgaGV4LnN1YnN0cmluZygxLCAyKSwgaGV4LnN1YnN0cmluZygxLCAyKVxuICAgICwgaGV4LnN1YnN0cmluZygyLCAzKSwgaGV4LnN1YnN0cmluZygyLCAzKVxuICAgICwgaGV4LnN1YnN0cmluZygzLCA0KSwgaGV4LnN1YnN0cmluZygzLCA0KVxuICAgIF0uam9pbignJykgOiBoZXhcbn1cblxuLy8gQ29tcG9uZW50IHRvIGhleCB2YWx1ZVxuZnVuY3Rpb24gY29tcFRvSGV4KGNvbXApIHtcbiAgdmFyIGhleCA9IGNvbXAudG9TdHJpbmcoMTYpXG4gIHJldHVybiBoZXgubGVuZ3RoID09IDEgPyAnMCcgKyBoZXggOiBoZXhcbn1cblxuLy8gQ2FsY3VsYXRlIHByb3BvcnRpb25hbCB3aWR0aCBhbmQgaGVpZ2h0IHZhbHVlcyB3aGVuIG5lY2Vzc2FyeVxuZnVuY3Rpb24gcHJvcG9ydGlvbmFsU2l6ZShlbGVtZW50LCB3aWR0aCwgaGVpZ2h0KSB7XG4gIGlmICh3aWR0aCA9PSBudWxsIHx8IGhlaWdodCA9PSBudWxsKSB7XG4gICAgdmFyIGJveCA9IGVsZW1lbnQuYmJveCgpXG5cbiAgICBpZiAod2lkdGggPT0gbnVsbClcbiAgICAgIHdpZHRoID0gYm94LndpZHRoIC8gYm94LmhlaWdodCAqIGhlaWdodFxuICAgIGVsc2UgaWYgKGhlaWdodCA9PSBudWxsKVxuICAgICAgaGVpZ2h0ID0gYm94LmhlaWdodCAvIGJveC53aWR0aCAqIHdpZHRoXG4gIH1cblxuICByZXR1cm4ge1xuICAgIHdpZHRoOiAgd2lkdGhcbiAgLCBoZWlnaHQ6IGhlaWdodFxuICB9XG59XG5cbi8vIERlbHRhIHRyYW5zZm9ybSBwb2ludFxuZnVuY3Rpb24gZGVsdGFUcmFuc2Zvcm1Qb2ludChtYXRyaXgsIHgsIHkpIHtcbiAgcmV0dXJuIHtcbiAgICB4OiB4ICogbWF0cml4LmEgKyB5ICogbWF0cml4LmMgKyAwXG4gICwgeTogeCAqIG1hdHJpeC5iICsgeSAqIG1hdHJpeC5kICsgMFxuICB9XG59XG5cbi8vIE1hcCBtYXRyaXggYXJyYXkgdG8gb2JqZWN0XG5mdW5jdGlvbiBhcnJheVRvTWF0cml4KGEpIHtcbiAgcmV0dXJuIHsgYTogYVswXSwgYjogYVsxXSwgYzogYVsyXSwgZDogYVszXSwgZTogYVs0XSwgZjogYVs1XSB9XG59XG5cbi8vIFBhcnNlIG1hdHJpeCBpZiByZXF1aXJlZFxuZnVuY3Rpb24gcGFyc2VNYXRyaXgobWF0cml4KSB7XG4gIGlmICghKG1hdHJpeCBpbnN0YW5jZW9mIFNWRy5NYXRyaXgpKVxuICAgIG1hdHJpeCA9IG5ldyBTVkcuTWF0cml4KG1hdHJpeClcblxuICByZXR1cm4gbWF0cml4XG59XG5cbi8vIEFkZCBjZW50cmUgcG9pbnQgdG8gdHJhbnNmb3JtIG9iamVjdFxuZnVuY3Rpb24gZW5zdXJlQ2VudHJlKG8sIHRhcmdldCkge1xuICBvLmN4ID0gby5jeCA9PSBudWxsID8gdGFyZ2V0LmJib3goKS5jeCA6IG8uY3hcbiAgby5jeSA9IG8uY3kgPT0gbnVsbCA/IHRhcmdldC5iYm94KCkuY3kgOiBvLmN5XG59XG5cbi8vIENvbnZlcnQgc3RyaW5nIHRvIG1hdHJpeFxuZnVuY3Rpb24gc3RyaW5nVG9NYXRyaXgoc291cmNlKSB7XG4gIC8vIHJlbW92ZSBtYXRyaXggd3JhcHBlciBhbmQgc3BsaXQgdG8gaW5kaXZpZHVhbCBudW1iZXJzXG4gIHNvdXJjZSA9IHNvdXJjZVxuICAgIC5yZXBsYWNlKFNWRy5yZWdleC53aGl0ZXNwYWNlLCAnJylcbiAgICAucmVwbGFjZShTVkcucmVnZXgubWF0cml4LCAnJylcbiAgICAuc3BsaXQoU1ZHLnJlZ2V4Lm1hdHJpeEVsZW1lbnRzKVxuXG4gIC8vIGNvbnZlcnQgc3RyaW5nIHZhbHVlcyB0byBmbG9hdHMgYW5kIGNvbnZlcnQgdG8gYSBtYXRyaXgtZm9ybWF0dGVkIG9iamVjdFxuICByZXR1cm4gYXJyYXlUb01hdHJpeChcbiAgICBTVkcudXRpbHMubWFwKHNvdXJjZSwgZnVuY3Rpb24obikge1xuICAgICAgcmV0dXJuIHBhcnNlRmxvYXQobilcbiAgICB9KVxuICApXG59XG5cbi8vIENhbGN1bGF0ZSBwb3NpdGlvbiBhY2NvcmRpbmcgdG8gZnJvbSBhbmQgdG9cbmZ1bmN0aW9uIGF0KG8sIHBvcykge1xuICAvLyBudW1iZXIgcmVjYWxjdWxhdGlvbiAoZG9uJ3QgYm90aGVyIGNvbnZlcnRpbmcgdG8gU1ZHLk51bWJlciBmb3IgcGVyZm9ybWFuY2UgcmVhc29ucylcbiAgcmV0dXJuIHR5cGVvZiBvLmZyb20gPT0gJ251bWJlcicgP1xuICAgIG8uZnJvbSArIChvLnRvIC0gby5mcm9tKSAqIHBvcyA6XG5cbiAgLy8gaW5zdGFuY2UgcmVjYWxjdWxhdGlvblxuICBvIGluc3RhbmNlb2YgU1ZHLkNvbG9yIHx8IG8gaW5zdGFuY2VvZiBTVkcuTnVtYmVyIHx8IG8gaW5zdGFuY2VvZiBTVkcuTWF0cml4ID8gby5hdChwb3MpIDpcblxuICAvLyBmb3IgYWxsIG90aGVyIHZhbHVlcyB3YWl0IHVudGlsIHBvcyBoYXMgcmVhY2hlZCAxIHRvIHJldHVybiB0aGUgZmluYWwgdmFsdWVcbiAgcG9zIDwgMSA/IG8uZnJvbSA6IG8udG9cbn1cblxuLy8gUGF0aEFycmF5IEhlbHBlcnNcbmZ1bmN0aW9uIGFycmF5VG9TdHJpbmcoYSkge1xuICBmb3IgKHZhciBpID0gMCwgaWwgPSBhLmxlbmd0aCwgcyA9ICcnOyBpIDwgaWw7IGkrKykge1xuICAgIHMgKz0gYVtpXVswXVxuXG4gICAgaWYgKGFbaV1bMV0gIT0gbnVsbCkge1xuICAgICAgcyArPSBhW2ldWzFdXG5cbiAgICAgIGlmIChhW2ldWzJdICE9IG51bGwpIHtcbiAgICAgICAgcyArPSAnICdcbiAgICAgICAgcyArPSBhW2ldWzJdXG5cbiAgICAgICAgaWYgKGFbaV1bM10gIT0gbnVsbCkge1xuICAgICAgICAgIHMgKz0gJyAnXG4gICAgICAgICAgcyArPSBhW2ldWzNdXG4gICAgICAgICAgcyArPSAnICdcbiAgICAgICAgICBzICs9IGFbaV1bNF1cblxuICAgICAgICAgIGlmIChhW2ldWzVdICE9IG51bGwpIHtcbiAgICAgICAgICAgIHMgKz0gJyAnXG4gICAgICAgICAgICBzICs9IGFbaV1bNV1cbiAgICAgICAgICAgIHMgKz0gJyAnXG4gICAgICAgICAgICBzICs9IGFbaV1bNl1cblxuICAgICAgICAgICAgaWYgKGFbaV1bN10gIT0gbnVsbCkge1xuICAgICAgICAgICAgICBzICs9ICcgJ1xuICAgICAgICAgICAgICBzICs9IGFbaV1bN11cbiAgICAgICAgICAgIH1cbiAgICAgICAgICB9XG4gICAgICAgIH1cbiAgICAgIH1cbiAgICB9XG4gIH1cblxuICByZXR1cm4gcyArICcgJ1xufVxuXG4vLyBEZWVwIG5ldyBpZCBhc3NpZ25tZW50XG5mdW5jdGlvbiBhc3NpZ25OZXdJZChub2RlKSB7XG4gIC8vIGRvIHRoZSBzYW1lIGZvciBTVkcgY2hpbGQgbm9kZXMgYXMgd2VsbFxuICBmb3IgKHZhciBpID0gbm9kZS5jaGlsZE5vZGVzLmxlbmd0aCAtIDE7IGkgPj0gMDsgaS0tKVxuICAgIGlmIChub2RlLmNoaWxkTm9kZXNbaV0gaW5zdGFuY2VvZiBTVkdFbGVtZW50KVxuICAgICAgYXNzaWduTmV3SWQobm9kZS5jaGlsZE5vZGVzW2ldKVxuXG4gIHJldHVybiBTVkcuYWRvcHQobm9kZSkuaWQoU1ZHLmVpZChub2RlLm5vZGVOYW1lKSlcbn1cblxuLy8gQWRkIG1vcmUgYm91bmRpbmcgYm94IHByb3BlcnRpZXNcbmZ1bmN0aW9uIGZ1bGxCb3goYikge1xuICBpZiAoYi54ID09IG51bGwpIHtcbiAgICBiLnggICAgICA9IDBcbiAgICBiLnkgICAgICA9IDBcbiAgICBiLndpZHRoICA9IDBcbiAgICBiLmhlaWdodCA9IDBcbiAgfVxuXG4gIGIudyAgPSBiLndpZHRoXG4gIGIuaCAgPSBiLmhlaWdodFxuICBiLngyID0gYi54ICsgYi53aWR0aFxuICBiLnkyID0gYi55ICsgYi5oZWlnaHRcbiAgYi5jeCA9IGIueCArIGIud2lkdGggLyAyXG4gIGIuY3kgPSBiLnkgKyBiLmhlaWdodCAvIDJcblxuICByZXR1cm4gYlxufVxuXG4vLyBHZXQgaWQgZnJvbSByZWZlcmVuY2Ugc3RyaW5nXG5mdW5jdGlvbiBpZEZyb21SZWZlcmVuY2UodXJsKSB7XG4gIHZhciBtID0gdXJsLnRvU3RyaW5nKCkubWF0Y2goU1ZHLnJlZ2V4LnJlZmVyZW5jZSlcblxuICBpZiAobSkgcmV0dXJuIG1bMV1cbn1cblxuLy8gQ3JlYXRlIG1hdHJpeCBhcnJheSBmb3IgbG9vcGluZ1xudmFyIGFiY2RlZiA9ICdhYmNkZWYnLnNwbGl0KCcnKVxuLy8gQWRkIEN1c3RvbUV2ZW50IHRvIElFOSBhbmQgSUUxMFxuaWYgKHR5cGVvZiBDdXN0b21FdmVudCAhPT0gJ2Z1bmN0aW9uJykge1xuICAvLyBDb2RlIGZyb206IGh0dHBzOi8vZGV2ZWxvcGVyLm1vemlsbGEub3JnL2VuLVVTL2RvY3MvV2ViL0FQSS9DdXN0b21FdmVudFxuICB2YXIgQ3VzdG9tRXZlbnQgPSBmdW5jdGlvbihldmVudCwgb3B0aW9ucykge1xuICAgIG9wdGlvbnMgPSBvcHRpb25zIHx8IHsgYnViYmxlczogZmFsc2UsIGNhbmNlbGFibGU6IGZhbHNlLCBkZXRhaWw6IHVuZGVmaW5lZCB9XG4gICAgdmFyIGUgPSBkb2N1bWVudC5jcmVhdGVFdmVudCgnQ3VzdG9tRXZlbnQnKVxuICAgIGUuaW5pdEN1c3RvbUV2ZW50KGV2ZW50LCBvcHRpb25zLmJ1YmJsZXMsIG9wdGlvbnMuY2FuY2VsYWJsZSwgb3B0aW9ucy5kZXRhaWwpXG4gICAgcmV0dXJuIGVcbiAgfVxuXG4gIEN1c3RvbUV2ZW50LnByb3RvdHlwZSA9IHdpbmRvdy5FdmVudC5wcm90b3R5cGVcblxuICB3aW5kb3cuQ3VzdG9tRXZlbnQgPSBDdXN0b21FdmVudFxufVxuXG4vLyByZXF1ZXN0QW5pbWF0aW9uRnJhbWUgLyBjYW5jZWxBbmltYXRpb25GcmFtZSBQb2x5ZmlsbCB3aXRoIGZhbGxiYWNrIGJhc2VkIG9uIFBhdWwgSXJpc2hcbihmdW5jdGlvbih3KSB7XG4gIHZhciBsYXN0VGltZSA9IDBcbiAgdmFyIHZlbmRvcnMgPSBbJ21veicsICd3ZWJraXQnXVxuXG4gIGZvcih2YXIgeCA9IDA7IHggPCB2ZW5kb3JzLmxlbmd0aCAmJiAhd2luZG93LnJlcXVlc3RBbmltYXRpb25GcmFtZTsgKyt4KSB7XG4gICAgdy5yZXF1ZXN0QW5pbWF0aW9uRnJhbWUgPSB3W3ZlbmRvcnNbeF0gKyAnUmVxdWVzdEFuaW1hdGlvbkZyYW1lJ11cbiAgICB3LmNhbmNlbEFuaW1hdGlvbkZyYW1lICA9IHdbdmVuZG9yc1t4XSArICdDYW5jZWxBbmltYXRpb25GcmFtZSddIHx8XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICB3W3ZlbmRvcnNbeF0gKyAnQ2FuY2VsUmVxdWVzdEFuaW1hdGlvbkZyYW1lJ11cbiAgfVxuXG4gIHcucmVxdWVzdEFuaW1hdGlvbkZyYW1lID0gdy5yZXF1ZXN0QW5pbWF0aW9uRnJhbWUgfHxcbiAgICBmdW5jdGlvbihjYWxsYmFjaykge1xuICAgICAgdmFyIGN1cnJUaW1lID0gbmV3IERhdGUoKS5nZXRUaW1lKClcbiAgICAgIHZhciB0aW1lVG9DYWxsID0gTWF0aC5tYXgoMCwgMTYgLSAoY3VyclRpbWUgLSBsYXN0VGltZSkpXG5cbiAgICAgIHZhciBpZCA9IHcuc2V0VGltZW91dChmdW5jdGlvbigpIHtcbiAgICAgICAgY2FsbGJhY2soY3VyclRpbWUgKyB0aW1lVG9DYWxsKVxuICAgICAgfSwgdGltZVRvQ2FsbClcblxuICAgICAgbGFzdFRpbWUgPSBjdXJyVGltZSArIHRpbWVUb0NhbGxcbiAgICAgIHJldHVybiBpZFxuICAgIH1cblxuICB3LmNhbmNlbEFuaW1hdGlvbkZyYW1lID0gdy5jYW5jZWxBbmltYXRpb25GcmFtZSB8fCB3LmNsZWFyVGltZW91dDtcblxufSh3aW5kb3cpKVxuXG5yZXR1cm4gU1ZHXG5cbn0pKTsiLCJpbXBvcnQgZml0Q3VydmUgZnJvbSAnZml0LWN1cnZlJztcblxuY29uc3QgZXJyb3IgPSAxMDA7XG5cbmZ1bmN0aW9uIFBhaW50Q29udHJvbChwYW5uZWwpIHtcblx0bGV0IHJhd1BvaW50RGF0YSA9IFtdO1xuXHRsZXQgcGFpbnRpbmdQb2x5TGluZSA9IHVuZGVmaW5lZDtcblxuXHR0aGlzLnN0YXJ0ID0gZnVuY3Rpb24oIHBvaW50ICkge1xuXHRcdHJhd1BvaW50RGF0YS5wdXNoKCBwb2ludCApO1xuXHRcdHBhaW50aW5nUG9seUxpbmUgPSBwYW5uZWwucG9seWxpbmUoKS5maWxsKCdub25lJykuc3Ryb2tlKHsgd2lkdGg6IDEgfSk7XG5cblx0fTtcblx0dGhpcy51cGRhdGUgPSBmdW5jdGlvbiggcG9pbnQgKSB7XG5cdFx0cmF3UG9pbnREYXRhLnB1c2goIHBvaW50ICk7XG5cdFx0dXBkYXRlTGluZXMoIHBhaW50aW5nUG9seUxpbmUsIHJhd1BvaW50RGF0YSk7XG5cdH07XG5cblx0dGhpcy5lbmQgPSBmdW5jdGlvbigpIHtcblx0XHRsZXQgc21vb3RoQml6ZXIgPSBmaXRDdXJ2ZSggcmF3UG9pbnREYXRhLCBlcnJvciApO1xuXHRcdGxldCBwYXRoU3RyaW5nID0gZml0dGVkQ3VydmVUb1BhdGhTdHJpbmcoc21vb3RoQml6ZXIpO1xuXG5cdFx0Ly8gZHJhdyBtYWduZXRpYyBjdXJ2ZVxuXHRcdGRyYXdPblBhbm5lbChwYW5uZWwsIHBhdGhTdHJpbmcpO1xuXHRcdGNsZWFyUmF3RGF0YSgpO1xuXHR9O1xuXG5cdGZ1bmN0aW9uIHVwZGF0ZUxpbmVzKHBhaW50aW5nUG9seUxpbmUsIHJhd1BvaW50RGF0YSkge1xuXHRcdHBhaW50aW5nUG9seUxpbmUucGxvdCggcmF3UG9pbnREYXRhICk7XG5cdH1cblx0ZnVuY3Rpb24gZml0dGVkQ3VydmVUb1BhdGhTdHJpbmcoZml0dGVkTGluZURhdGEpIHtcblx0XHR2YXIgc3RyID0gJyc7XG5cdFx0Ly9iZXppZXIgOiBbIFtjMF0sIFtjMV0sIFtjMl0sIFtjM10gXVxuXHRcdGZpdHRlZExpbmVEYXRhLm1hcChmdW5jdGlvbiAoYmV6aWVyLCBpKSB7XG5cdFx0XHRpZiAoaSA9PSAwKSB7XG5cdFx0XHRcdHN0ciArPSAnTSAnICsgYmV6aWVyWzBdWzBdICsgJyAnICsgYmV6aWVyWzBdWzFdO1xuXHRcdFx0fVxuXG5cdFx0XHRzdHIgKz0gJ0MgJyArIGJlemllclsxXVswXSArICcgJyArIGJlemllclsxXVsxXSArICcsICcgK1xuXHRcdFx0YmV6aWVyWzJdWzBdICsgJyAnICsgYmV6aWVyWzJdWzFdICsgJywgJyArXG5cdFx0XHRiZXppZXJbM11bMF0gKyAnICcgKyBiZXppZXJbM11bMV0gKyAnICc7XHRcblx0XHRcdFx0XHRcdFxuXHRcdH0pO1xuXG5cdFx0cmV0dXJuIHN0cjtcblx0fVxuXHRmdW5jdGlvbiBkcmF3T25QYW5uZWwocGFubmVsLCBwYXRoU3RyaW5nKXtcblx0XHRwYW5uZWwucGF0aCggcGF0aFN0cmluZyApLmZpbGwoJ25vbmUnKS5zdHJva2UoeyB3aWR0aDogMyB9KS5zdHJva2UoJyNmMDYnKTtcblx0fVxuXHRmdW5jdGlvbiBjbGVhclJhd0RhdGEoKXtcblx0XHRyYXdQb2ludERhdGEubGVuZ3RoID0gMDtcblx0XHRwYWludGluZ1BvbHlMaW5lLnJlbW92ZSgpO1xuXHR9XG59XG5cbmV4cG9ydCBkZWZhdWx0IFBhaW50Q29udHJvbDtcblxuIiwiaW1wb3J0IFBhaW50Q29udHJvbCBmcm9tICcuL0NvbnRyb2xzL1BhaW50Q29udHJvbCc7XG5pbXBvcnQgU1ZHIGZyb20gJ3N2Zy5qcyc7XG5cbnZhciBkcmF3ID0gU1ZHKCdkcmF3aW5nJykuc2l6ZSgzMDAsIDMwMCk7XG5cblxuc2V0Q29udHJvbChkcmF3KTtcblxuZnVuY3Rpb24gc2V0Q29udHJvbChfY29udGFpbmVyKSB7XG5cdGxldCBpc01vdXNlRG93biA9IGZhbHNlO1xuXHRsZXQgY3Vycm5ldENvbnRyb2wgPSBuZXcgUGFpbnRDb250cm9sKGRyYXcpO1xuXHRjb25zdCB0b3AgPSBkcmF3Lm5vZGUuZ2V0Qm91bmRpbmdDbGllbnRSZWN0KCkudG9wO1xuXHRjb25zdCBsZWZ0ID0gZHJhdy5ub2RlLmdldEJvdW5kaW5nQ2xpZW50UmVjdCgpLmxlZnQ7XG5cdFxuXHRfY29udGFpbmVyLm9uKCdtb3VzZWRvd24nLCBmdW5jdGlvbiAoZSkge1xuXHRcdGNvbnN0IHBvaW50ID0gW1xuXHRcdFx0ZS5jbGllbnRYIC0gdG9wLFxuXHRcdFx0ZS5jbGllbnRZIC0gbGVmdFxuXHRcdF07XG5cdFx0aXNNb3VzZURvd24gPSB0cnVlO1xuXHRcdGN1cnJuZXRDb250cm9sLnN0YXJ0KHBvaW50KTtcblxuXHR9KTtcblx0X2NvbnRhaW5lci5vbignbW91c2V1cCcsIGZ1bmN0aW9uICgpIHtcblx0XHRpc01vdXNlRG93biA9IGZhbHNlO1xuXHRcdGN1cnJuZXRDb250cm9sLmVuZCgpO1xuXHR9KTtcblx0X2NvbnRhaW5lci5vbignbW91c2Vtb3ZlJywgZnVuY3Rpb24gKGUpIHtcblx0XHR2YXIgeCA9IGUub2Zmc2V0WDtcblx0XHR2YXIgeSA9IGUub2Zmc2V0WTtcblx0XHRpZiAoaXNNb3VzZURvd24pIHtcblx0XHRcdGN1cnJuZXRDb250cm9sLnVwZGF0ZShbeCwgeV0pO1xuXHRcdH1cblx0fSk7XG59XG5cbiJdfQ==
