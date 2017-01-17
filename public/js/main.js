(function e(t,n,r){function s(o,u){if(!n[o]){if(!t[o]){var a=typeof require=='function'&&require;if(!u&&a)return a(o,!0);if(i)return i(o,!0);var f=new Error('Cannot find module \''+o+'\'');throw f.code='MODULE_NOT_FOUND',f;}var l=n[o]={exports:{}};t[o][0].call(l.exports,function(e){var n=t[o][1][e];return s(n?n:e);},l,l.exports,e,t,n,r);}return n[o].exports;}var i=typeof require=='function'&&require;for(var o=0;o<r.length;o++)s(r[o]);return s;})({1:[function(require,module,exports){
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
			throw new TypeError('First argument should be an array');
		}
		points.forEach((point) => {
			if(!Array.isArray(point) || point.length !== 2
        || typeof point[0] !== 'number' || typeof point[1] !== 'number'){
				throw Error('Each point should be an array of two numbers');
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
		[bezCurve, maxError, splitPoint] = generateAndReport(points, u, u, leftTangent, rightTangent, progressCallback);

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
	}

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
	}

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
	}

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
	}

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
		});
		u = u.map(x => x/prevU);

		return u;
	}

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
	}

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
			return zs;
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
*/
	(function(root, factory) {
		if (typeof define === 'function' && define.amd) {
			define(function(){
				return factory(root, root.document);
			});
		} else if (typeof exports === 'object') {
			module.exports = root.document ? factory(root, root.document) : function(w){ return factory(w, w.document); };
		} else {
			root.SVG = factory(root, root.document);
		}
	}(typeof window !== 'undefined' ? window : this, function(window, document) {

// The main wrapping element
		var SVG = this.SVG = function(element) {
			if (SVG.supported) {
				element = new SVG.Doc(element);

				if(!SVG.parser.draw)
					SVG.prepare();

				return element;
			}
		};

// Default namespaces
		SVG.ns    = 'http://www.w3.org/2000/svg';
		SVG.xmlns = 'http://www.w3.org/2000/xmlns/';
		SVG.xlink = 'http://www.w3.org/1999/xlink';
		SVG.svgjs = 'http://svgjs.com/svgjs';

// Svg support test
		SVG.supported = (function() {
			return !! document.createElementNS &&
         !! document.createElementNS(SVG.ns,'svg').createSVGRect;
		})();

// Don't bother to continue if SVG is not supported
		if (!SVG.supported) return false;

// Element id sequence
		SVG.did  = 1000;

// Get next named element id
		SVG.eid = function(name) {
			return 'Svgjs' + capitalize(name) + (SVG.did++);
		};

// Method for element creation
		SVG.create = function(name) {
  // create element
			var element = document.createElementNS(this.ns, name);

  // apply unique id
			element.setAttribute('id', this.eid(name));

			return element;
		};

// Method for extending objects
		SVG.extend = function() {
			var modules, methods, key, i;

  // Get list of modules
			modules = [].slice.call(arguments);

  // Get object with extensions
			methods = modules.pop();

			for (i = modules.length - 1; i >= 0; i--)
				if (modules[i])
					for (key in methods)
						modules[i].prototype[key] = methods[key];

  // Make sure SVG.Set inherits any newly added methods
			if (SVG.Set && SVG.Set.inherit)
				SVG.Set.inherit();
		};

// Invent new element
		SVG.invent = function(config) {
  // Create element initializer
			var initializer = typeof config.create == 'function' ?
    config.create :
    function() {
	this.constructor.call(this, SVG.create(config.create));
};

  // Inherit prototype
			if (config.inherit)
				initializer.prototype = new config.inherit;

  // Extend with methods
			if (config.extend)
				SVG.extend(initializer, config.extend);

  // Attach construct method to parent
			if (config.construct)
				SVG.extend(config.parent || SVG.Container, config.construct);

			return initializer;
		};

// Adopt existing svg elements
		SVG.adopt = function(node) {
  // check for presence of node
			if (!node) return null;

  // make sure a node isn't already adopted
			if (node.instance) return node.instance;

  // initialize variables
			var element;

  // adopt with element-specific settings
			if (node.nodeName == 'svg')
				element = node.parentNode instanceof SVGElement ? new SVG.Nested : new SVG.Doc;
			else if (node.nodeName == 'linearGradient')
				element = new SVG.Gradient('linear');
			else if (node.nodeName == 'radialGradient')
				element = new SVG.Gradient('radial');
			else if (SVG[capitalize(node.nodeName)])
				element = new SVG[capitalize(node.nodeName)];
			else
    element = new SVG.Element(node);

  // ensure references
			element.type  = node.nodeName;
			element.node  = node;
			node.instance = element;

  // SVG.Class specific preparations
			if (element instanceof SVG.Doc)
				element.namespace().defs();

  // pull svgjs data from the dom (getAttributeNS doesn't work in html5)
			element.setData(JSON.parse(node.getAttribute('svgjs:data')) || {});

			return element;
		};

// Initialize parsing element
		SVG.prepare = function() {
  // Select document body and create invisible svg element
			var body = document.getElementsByTagName('body')[0]
    , draw = (body ? new SVG.Doc(body) :  new SVG.Doc(document.documentElement).nested()).size(2, 0);

  // Create parser object
			SVG.parser = {
				body: body || document.documentElement
  , draw: draw.style('opacity:0;position:fixed;left:100%;top:100%;overflow:hidden')
  , poly: draw.polyline().node
  , path: draw.path().node
  , native: SVG.create('svg')
			};
		};

		SVG.parser = {
			native: SVG.create('svg')
		};

		document.addEventListener('DOMContentLoaded', function() {
			if(!SVG.parser.draw)
				SVG.prepare();
		}, false);

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
		};

		SVG.utils = {
  // Map function
			map: function(array, block) {
				var i
      , il = array.length
      , result = [];

				for (i = 0; i < il; i++)
					result.push(block(array[i]));

				return result;
			}

  // Filter function
, filter: function(array, block) {
	var i
      , il = array.length
      , result = [];

	for (i = 0; i < il; i++)
		if (block(array[i]))
			result.push(array[i]);

	return result;
}

  // Degrees to radians
, radians: function(d) {
	return d % 360 * Math.PI / 180;
}

  // Radians to degrees
, degrees: function(r) {
	return r * 180 / Math.PI % 360;
}

, filterSVGElements: function(nodes) {
	return this.filter( nodes, function(el) { return el instanceof SVGElement; });
}

		};

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

		};
// Module for color convertions
		SVG.Color = function(color) {
			var match;

  // initialize defaults
			this.r = 0;
			this.g = 0;
			this.b = 0;

			if(!color) return;

  // parse color
			if (typeof color === 'string') {
				if (SVG.regex.isRgb.test(color)) {
      // get rgb values
					match = SVG.regex.rgb.exec(color.replace(/\s/g,''));

      // parse numeric values
					this.r = parseInt(match[1]);
					this.g = parseInt(match[2]);
					this.b = parseInt(match[3]);

				} else if (SVG.regex.isHex.test(color)) {
      // get hex values
					match = SVG.regex.hex.exec(fullHex(color));

      // parse numeric values
					this.r = parseInt(match[1], 16);
					this.g = parseInt(match[2], 16);
					this.b = parseInt(match[3], 16);

				}

			} else if (typeof color === 'object') {
				this.r = color.r;
				this.g = color.g;
				this.b = color.b;

			}

		};

		SVG.extend(SVG.Color, {
  // Default to hex conversion
			toString: function() {
				return this.toHex();
			}
  // Build hex value
, toHex: function() {
	return '#'
      + compToHex(this.r)
      + compToHex(this.g)
      + compToHex(this.b);
}
  // Build rgb value
, toRgb: function() {
	return 'rgb(' + [this.r, this.g, this.b].join() + ')';
}
  // Calculate true brightness
, brightness: function() {
	return (this.r / 255 * 0.30)
         + (this.g / 255 * 0.59)
         + (this.b / 255 * 0.11);
}
  // Make color morphable
, morph: function(color) {
	this.destination = new SVG.Color(color);

	return this;
}
  // Get morphed color at given position
, at: function(pos) {
    // make sure a destination is defined
	if (!this.destination) return this;

    // normalise pos
	pos = pos < 0 ? 0 : pos > 1 ? 1 : pos;

    // generate morphed color
	return new SVG.Color({
		r: ~~(this.r + (this.destination.r - this.r) * pos)
    , g: ~~(this.g + (this.destination.g - this.g) * pos)
    , b: ~~(this.b + (this.destination.b - this.b) * pos)
	});
}

		});

// Testers

// Test if given value is a color string
		SVG.Color.test = function(color) {
			color += '';
			return SVG.regex.isHex.test(color)
      || SVG.regex.isRgb.test(color);
		};

// Test if given value is a rgb object
		SVG.Color.isRgb = function(color) {
			return color && typeof color.r == 'number'
               && typeof color.g == 'number'
               && typeof color.b == 'number';
		};

// Test if given value is a color
		SVG.Color.isColor = function(color) {
			return SVG.Color.isRgb(color) || SVG.Color.test(color);
		};
// Module for array conversion
		SVG.Array = function(array, fallback) {
			array = (array || []).valueOf();

  // if array is empty and fallback is provided, use fallback
			if (array.length == 0 && fallback)
				array = fallback.valueOf();

  // parse array
			this.value = this.parse(array);
		};

		SVG.extend(SVG.Array, {
  // Make array morphable
			morph: function(array) {
				this.destination = this.parse(array);

    // normalize length of arrays
				if (this.value.length != this.destination.length) {
					var lastValue       = this.value[this.value.length - 1]
        , lastDestination = this.destination[this.destination.length - 1];

					while(this.value.length > this.destination.length)
						this.destination.push(lastDestination);
					while(this.value.length < this.destination.length)
						this.value.push(lastValue);
				}

				return this;
			}
  // Clean up any duplicate points
, settle: function() {
    // find all unique values
	for (var i = 0, il = this.value.length, seen = []; i < il; i++)
		if (seen.indexOf(this.value[i]) == -1)
			seen.push(this.value[i]);

    // set new value
	return this.value = seen;
}
  // Get morphed array at given position
, at: function(pos) {
    // make sure a destination is defined
	if (!this.destination) return this;

    // generate morphed array
	for (var i = 0, il = this.value.length, array = []; i < il; i++)
		array.push(this.value[i] + (this.destination[i] - this.value[i]) * pos);

	return new SVG.Array(array);
}
  // Convert array to string
, toString: function() {
	return this.value.join(' ');
}
  // Real value
, valueOf: function() {
	return this.value;
}
  // Parse whitespace separated string
, parse: function(array) {
	array = array.valueOf();

    // if already is an array, no need to parse it
	if (Array.isArray(array)) return array;

	return this.split(array);
}
  // Strip unnecessary whitespace
, split: function(string) {
	return string.trim().split(/\s+/);
}
  // Reverse array
, reverse: function() {
	this.value.reverse();

	return this;
}

		});
// Poly points array
		SVG.PointArray = function(array, fallback) {
			this.constructor.call(this, array, fallback || [[0,0]]);
		};

// Inherit from SVG.Array
		SVG.PointArray.prototype = new SVG.Array;

		SVG.extend(SVG.PointArray, {
  // Convert array to string
			toString: function() {
    // convert to a poly point string
				for (var i = 0, il = this.value.length, array = []; i < il; i++)
					array.push(this.value[i].join(','));

				return array.join(' ');
			}
  // Convert array to line object
, toLine: function() {
	return {
		x1: this.value[0][0]
    , y1: this.value[0][1]
    , x2: this.value[1][0]
    , y2: this.value[1][1]
	};
}
  // Get morphed array at given position
, at: function(pos) {
    // make sure a destination is defined
	if (!this.destination) return this;

    // generate morphed point string
	for (var i = 0, il = this.value.length, array = []; i < il; i++)
		array.push([
			this.value[i][0] + (this.destination[i][0] - this.value[i][0]) * pos
      , this.value[i][1] + (this.destination[i][1] - this.value[i][1]) * pos
		]);

	return new SVG.PointArray(array);
}
  // Parse point string
, parse: function(array) {
	var points = [];

	array = array.valueOf();

    // if already is an array, no need to parse it
	if (Array.isArray(array)) return array;

    // parse points
	array = array.trim().split(/\s+|,/);

    // validate points - https://svgwg.org/svg2-draft/shapes.html#DataTypePoints
    // Odd number of coordinates is an error. In such cases, drop the last odd coordinate.
	if (array.length % 2 !== 0) array.pop();

    // wrap points in two-tuples and parse points as floats
	for(var i = 0, len = array.length; i < len; i = i + 2)
		points.push([ parseFloat(array[i]), parseFloat(array[i+1]) ]);

	return points;
}
  // Move point string
, move: function(x, y) {
	var box = this.bbox();

    // get relative offset
	x -= box.x;
	y -= box.y;

    // move every point
	if (!isNaN(x) && !isNaN(y))
		for (var i = this.value.length - 1; i >= 0; i--)
			this.value[i] = [this.value[i][0] + x, this.value[i][1] + y];

	return this;
}
  // Resize poly string
, size: function(width, height) {
	var i, box = this.bbox();

    // recalculate position of all points according to new size
	for (i = this.value.length - 1; i >= 0; i--) {
		this.value[i][0] = ((this.value[i][0] - box.x) * width)  / box.width  + box.x;
		this.value[i][1] = ((this.value[i][1] - box.y) * height) / box.height + box.y;
	}

	return this;
}
  // Get bounding box of points
, bbox: function() {
	SVG.parser.poly.setAttribute('points', this.toString());

	return SVG.parser.poly.getBBox();
}

		});
// Path points array
		SVG.PathArray = function(array, fallback) {
			this.constructor.call(this, array, fallback || [['M', 0, 0]]);
		};

// Inherit from SVG.Array
		SVG.PathArray.prototype = new SVG.Array;

		SVG.extend(SVG.PathArray, {
  // Convert array to string
			toString: function() {
				return arrayToString(this.value);
			}
  // Move path string
, move: function(x, y) {
    // get bounding box of current situation
	var box = this.bbox();

    // get relative offset
	x -= box.x;
	y -= box.y;

	if (!isNaN(x) && !isNaN(y)) {
      // move every point
		for (var l, i = this.value.length - 1; i >= 0; i--) {
			l = this.value[i][0];

			if (l == 'M' || l == 'L' || l == 'T')  {
				this.value[i][1] += x;
				this.value[i][2] += y;

			} else if (l == 'H')  {
				this.value[i][1] += x;

			} else if (l == 'V')  {
				this.value[i][1] += y;

			} else if (l == 'C' || l == 'S' || l == 'Q')  {
				this.value[i][1] += x;
				this.value[i][2] += y;
				this.value[i][3] += x;
				this.value[i][4] += y;

				if (l == 'C')  {
					this.value[i][5] += x;
					this.value[i][6] += y;
				}

			} else if (l == 'A')  {
				this.value[i][6] += x;
				this.value[i][7] += y;
			}

		}
	}

	return this;
}
  // Resize path string
, size: function(width, height) {
    // get bounding box of current situation
	var i, l, box = this.bbox();

    // recalculate position of all points according to new size
	for (i = this.value.length - 1; i >= 0; i--) {
		l = this.value[i][0];

		if (l == 'M' || l == 'L' || l == 'T')  {
			this.value[i][1] = ((this.value[i][1] - box.x) * width)  / box.width  + box.x;
			this.value[i][2] = ((this.value[i][2] - box.y) * height) / box.height + box.y;

		} else if (l == 'H')  {
			this.value[i][1] = ((this.value[i][1] - box.x) * width)  / box.width  + box.x;

		} else if (l == 'V')  {
			this.value[i][1] = ((this.value[i][1] - box.y) * height) / box.height + box.y;

		} else if (l == 'C' || l == 'S' || l == 'Q')  {
			this.value[i][1] = ((this.value[i][1] - box.x) * width)  / box.width  + box.x;
			this.value[i][2] = ((this.value[i][2] - box.y) * height) / box.height + box.y;
			this.value[i][3] = ((this.value[i][3] - box.x) * width)  / box.width  + box.x;
			this.value[i][4] = ((this.value[i][4] - box.y) * height) / box.height + box.y;

			if (l == 'C')  {
				this.value[i][5] = ((this.value[i][5] - box.x) * width)  / box.width  + box.x;
				this.value[i][6] = ((this.value[i][6] - box.y) * height) / box.height + box.y;
			}

		} else if (l == 'A')  {
        // resize radii
			this.value[i][1] = (this.value[i][1] * width)  / box.width;
			this.value[i][2] = (this.value[i][2] * height) / box.height;

        // move position values
			this.value[i][6] = ((this.value[i][6] - box.x) * width)  / box.width  + box.x;
			this.value[i][7] = ((this.value[i][7] - box.y) * height) / box.height + box.y;
		}

	}

	return this;
}
  // Test if the passed path array use the same path data commands as this path array
, equalCommands: function(pathArray) {
	var i, il, equalCommands;

	pathArray = new SVG.PathArray(pathArray);

	equalCommands = this.value.length === pathArray.value.length;
	for(i = 0, il = this.value.length; equalCommands && i < il; i++) {
		equalCommands = this.value[i][0] === pathArray.value[i][0];
	}

	return equalCommands;
}
  // Make path array morphable
, morph: function(pathArray) {
	pathArray = new SVG.PathArray(pathArray);

	if(this.equalCommands(pathArray)) {
		this.destination = pathArray;
	} else {
		this.destination = null;
	}

	return this;
}
  // Get morphed path array at given position
, at: function(pos) {
    // make sure a destination is defined
	if (!this.destination) return this;

	var sourceArray = this.value
      , destinationArray = this.destination.value
      , array = [], pathArray = new SVG.PathArray()
      , i, il, j, jl;

    // Animate has specified in the SVG spec
    // See: https://www.w3.org/TR/SVG11/paths.html#PathElement
	for (i = 0, il = sourceArray.length; i < il; i++) {
		array[i] = [sourceArray[i][0]];
		for(j = 1, jl = sourceArray[i].length; j < jl; j++) {
			array[i][j] = sourceArray[i][j] + (destinationArray[i][j] - sourceArray[i][j]) * pos;
		}
      // For the two flags of the elliptical arc command, the SVG spec say:
      // Flags and booleans are interpolated as fractions between zero and one, with any non-zero value considered to be a value of one/true
      // Elliptical arc command as an array followed by corresponding indexes:
      // ['A', rx, ry, x-axis-rotation, large-arc-flag, sweep-flag, x, y]
      //   0    1   2        3                 4             5      6  7
		if(array[i][0] === 'A') {
			array[i][4] = +(array[i][4] != 0);
			array[i][5] = +(array[i][5] != 0);
		}
	}

    // Directly modify the value of a path array, this is done this way for performance
	pathArray.value = array;
	return pathArray;
}
  // Absolutize and parse path to array
, parse: function(array) {
    // if it's already a patharray, no need to parse it
	if (array instanceof SVG.PathArray) return array.valueOf();

    // prepare for parsing
	var i, x0, y0, s, seg, arr
      , x = 0
      , y = 0
      , paramCnt = { 'M':2, 'L':2, 'H':1, 'V':1, 'C':6, 'S':4, 'Q':4, 'T':2, 'A':7 };

	if(typeof array == 'string'){

		array = array
        .replace(SVG.regex.negExp, 'X')         // replace all negative exponents with certain char
        .replace(SVG.regex.pathLetters, ' $& ') // put some room between letters and numbers
        .replace(SVG.regex.hyphen, ' -')        // add space before hyphen
        .replace(SVG.regex.comma, ' ')          // unify all spaces
        .replace(SVG.regex.X, 'e-')             // add back the expoent
        .trim()                                 // trim
        .split(SVG.regex.whitespaces);           // split into array

      // at this place there could be parts like ['3.124.854.32'] because we could not determine the point as seperator till now
      // we fix this elements in the next loop
		for(i = array.length; --i;){
			if(array[i].indexOf('.') != array[i].lastIndexOf('.')){
				var split = array[i].split('.'); // split at the point
				var first = [split.shift(), split.shift()].join('.'); // join the first number together
				array.splice.apply(array, [i, 1].concat(first, split.map(function(el){ return '.'+el; }))); // add first and all other entries back to array
			}
		}

	}else{
		array = array.reduce(function(prev, curr){
			return [].concat.apply(prev, curr);
		}, []);
	}

    // array now is an array containing all parts of a path e.g. ['M', '0', '0', 'L', '30', '30' ...]

	var arr = [];

	do{

      // Test if we have a path letter
		if(SVG.regex.isPathLetter.test(array[0])){
			s = array[0];
			array.shift();
      // If last letter was a move command and we got no new, it defaults to [L]ine
		}else if(s == 'M'){
			s = 'L';
		}else if(s == 'm'){
			s = 'l';
		}

      // add path letter as first element
		seg = [s.toUpperCase()];

      // push all necessary parameters to segment
		for(i = 0; i < paramCnt[seg[0]]; ++i){
			seg.push(parseFloat(array.shift()));
		}

      // upper case
		if(s == seg[0]){

			if(s == 'M' || s == 'L' || s == 'C' || s == 'Q' || s == 'S' || s == 'T'){
				x = seg[paramCnt[seg[0]]-1];
				y = seg[paramCnt[seg[0]]];
			}else if(s == 'V'){
				y = seg[1];
			}else if(s == 'H'){
				x = seg[1];
			}else if(s == 'A'){
				x = seg[6];
				y = seg[7];
			}

      // lower case
		}else{

        // convert relative to absolute values
			if(s == 'm' || s == 'l' || s == 'c' || s == 's' || s == 'q' || s == 't'){

				seg[1] += x;
				seg[2] += y;

				if(seg[3] != null){
					seg[3] += x;
					seg[4] += y;
				}

				if(seg[5] != null){
					seg[5] += x;
					seg[6] += y;
				}

          // move pointer
				x = seg[paramCnt[seg[0]]-1];
				y = seg[paramCnt[seg[0]]];

			}else if(s == 'v'){
				seg[1] += y;
				y = seg[1];
			}else if(s == 'h'){
				seg[1] += x;
				x = seg[1];
			}else if(s == 'a'){
				seg[6] += x;
				seg[7] += y;
				x = seg[6];
				y = seg[7];
			}

		}

		if(seg[0] == 'M'){
			x0 = x;
			y0 = y;
		}

		if(seg[0] == 'Z'){
			x = x0;
			y = y0;
		}

		arr.push(seg);

	}while(array.length);

	return arr;

}
  // Get bounding box of path
, bbox: function() {
	SVG.parser.path.setAttribute('d', this.toString());

	return SVG.parser.path.getBBox();
}

		});

// Module for unit convertions
		SVG.Number = SVG.invent({
  // Initialize
			create: function(value, unit) {
    // initialize defaults
				this.value = 0;
				this.unit  = unit || '';

    // parse value
				if (typeof value === 'number') {
      // ensure a valid numeric value
					this.value = isNaN(value) ? 0 : !isFinite(value) ? (value < 0 ? -3.4e+38 : +3.4e+38) : value;

				} else if (typeof value === 'string') {
					unit = value.match(SVG.regex.numberAndUnit);

					if (unit) {
        // make value numeric
						this.value = parseFloat(unit[1]);

        // normalize
						if (unit[5] == '%')
							this.value /= 100;
						else if (unit[5] == 's')
							this.value *= 1000;

        // store unit
						this.unit = unit[5];
					}

				} else {
					if (value instanceof SVG.Number) {
						this.value = value.valueOf();
						this.unit  = value.unit;
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
      ) + this.unit;
	}
  , toJSON: function() {
	return this.toString();
}
  , // Convert to primitive
	valueOf: function() {
		return this.value;
	}
    // Add number
  , plus: function(number) {
	return new SVG.Number(this + new SVG.Number(number), this.unit);
}
    // Subtract number
  , minus: function(number) {
	return this.plus(-new SVG.Number(number));
}
    // Multiply number
  , times: function(number) {
	return new SVG.Number(this * new SVG.Number(number), this.unit);
}
    // Divide number
  , divide: function(number) {
	return new SVG.Number(this / new SVG.Number(number), this.unit);
}
    // Convert to different unit
  , to: function(unit) {
	var number = new SVG.Number(this);

	if (typeof unit === 'string')
		number.unit = unit;

	return number;
}
    // Make number morphable
  , morph: function(number) {
	this.destination = new SVG.Number(number);

	return this;
}
    // Get morphed number at given position
  , at: function(pos) {
      // Make sure a destination is defined
	if (!this.destination) return this;

      // Generate new morphed number
	return new SVG.Number(this.destination)
          .minus(this)
          .times(pos)
          .plus(this);
}

}
		});

		SVG.Element = SVG.invent({
  // Initialize node
			create: function(node) {
    // make stroke value accessible dynamically
				this._stroke = SVG.defaults.attrs.stroke;

    // initialize data object
				this.dom = {};

    // create circular reference
				if (this.node = node) {
					this.type = node.nodeName;
					this.node.instance = this;

      // store current attribute value
					this._stroke = node.getAttribute('stroke') || this._stroke;
				}
			}

  // Add class methods
, extend: {
    // Move over x-axis
	x: function(x) {
		return this.attr('x', x);
	}
    // Move over y-axis
  , y: function(y) {
	return this.attr('y', y);
}
    // Move by center over x-axis
  , cx: function(x) {
	return x == null ? this.x() + this.width() / 2 : this.x(x - this.width() / 2);
}
    // Move by center over y-axis
  , cy: function(y) {
	return y == null ? this.y() + this.height() / 2 : this.y(y - this.height() / 2);
}
    // Move element to given x and y values
  , move: function(x, y) {
	return this.x(x).y(y);
}
    // Move element by its center
  , center: function(x, y) {
	return this.cx(x).cy(y);
}
    // Set width of element
  , width: function(width) {
	return this.attr('width', width);
}
    // Set height of element
  , height: function(height) {
	return this.attr('height', height);
}
    // Set element size to given width and height
  , size: function(width, height) {
	var p = proportionalSize(this, width, height);

	return this
        .width(new SVG.Number(p.width))
        .height(new SVG.Number(p.height));
}
    // Clone element
  , clone: function(parent) {
      // clone element and assign new id
	var clone = assignNewId(this.node.cloneNode(true));

      // insert the clone in the given parent or after myself
	if(parent) parent.add(clone);
	else this.after(clone);

	return clone;
}
    // Remove element
  , remove: function() {
	if (this.parent())
		this.parent().removeElement(this);

	return this;
}
    // Replace element
  , replace: function(element) {
	this.after(element).remove();

	return element;
}
    // Add element to given container and return self
  , addTo: function(parent) {
	return parent.put(this);
}
    // Add element to given container and return container
  , putIn: function(parent) {
	return parent.add(this);
}
    // Get / set id
  , id: function(id) {
	return this.attr('id', id);
}
    // Checks whether the given point inside the bounding box of the element
  , inside: function(x, y) {
	var box = this.bbox();

	return x > box.x
          && y > box.y
          && x < box.x + box.width
          && y < box.y + box.height;
}
    // Show element
  , show: function() {
	return this.style('display', '');
}
    // Hide element
  , hide: function() {
	return this.style('display', 'none');
}
    // Is element visible?
  , visible: function() {
	return this.style('display') != 'none';
}
    // Return id on string conversion
  , toString: function() {
	return this.attr('id');
}
    // Return array of classes on the node
  , classes: function() {
	var attr = this.attr('class');

	return attr == null ? [] : attr.trim().split(/\s+/);
}
    // Return true if class exists on the node, false otherwise
  , hasClass: function(name) {
	return this.classes().indexOf(name) != -1;
}
    // Add class to the node
  , addClass: function(name) {
	if (!this.hasClass(name)) {
		var array = this.classes();
		array.push(name);
		this.attr('class', array.join(' '));
	}

	return this;
}
    // Remove class from the node
  , removeClass: function(name) {
	if (this.hasClass(name)) {
		this.attr('class', this.classes().filter(function(c) {
			return c != name;
		}).join(' '));
	}

	return this;
}
    // Toggle the presence of a class on the node
  , toggleClass: function(name) {
	return this.hasClass(name) ? this.removeClass(name) : this.addClass(name);
}
    // Get referenced element form attribute value
  , reference: function(attr) {
	return SVG.get(this.attr(attr));
}
    // Returns the parent element instance
  , parent: function(type) {
	var parent = this;

      // check for parent
	if(!parent.node.parentNode) return null;

      // get parent element
	parent = SVG.adopt(parent.node.parentNode);

	if(!type) return parent;

      // loop trough ancestors if type is given
	while(parent && parent.node instanceof SVGElement){
		if(typeof type === 'string' ? parent.matches(type) : parent instanceof type) return parent;
		parent = SVG.adopt(parent.node.parentNode);
	}
}
    // Get parent document
  , doc: function() {
	return this instanceof SVG.Doc ? this : this.parent(SVG.Doc);
}
    // return array of all ancestors of given type up to the root svg
  , parents: function(type) {
	var parents = [], parent = this;

	do{
		parent = parent.parent(type);
		if(!parent || !parent.node) break;

		parents.push(parent);
	} while(parent.parent);

	return parents;
}
    // matches the element vs a css selector
  , matches: function(selector){
	return matches(this.node, selector);
}
    // Returns the svg node to call native svg methods on it
  , native: function() {
	return this.node;
}
    // Import raw svg
  , svg: function(svg) {
      // create temporary holder
	var well = document.createElement('svg');

      // act as a setter if svg is given
	if (svg && this instanceof SVG.Parent) {
        // dump raw svg
		well.innerHTML = '<svg>' + svg.replace(/\n/, '').replace(/<(\w+)([^<]+?)\/>/g, '<$1$2></$1>') + '</svg>';

        // transplant nodes
		for (var i = 0, il = well.firstChild.childNodes.length; i < il; i++)
			this.node.appendChild(well.firstChild.firstChild);

      // otherwise act as a getter
	} else {
        // create a wrapping svg element in case of partial content
		well.appendChild(svg = document.createElement('svg'));

        // write svgjs data to the dom
		this.writeDataToDom();

        // insert a copy of this node
		svg.appendChild(this.node.cloneNode(true));

        // return target element
		return well.innerHTML.replace(/^<svg>/, '').replace(/<\/svg>$/, '');
	}

	return this;
}
  // write svgjs data to the dom
  , writeDataToDom: function() {

      // dump variables recursively
	if(this.each || this.lines){
		var fn = this.each ? this : this.lines();
		fn.each(function(){
			this.writeDataToDom();
		});
	}

      // remove previously set data
	this.node.removeAttribute('svgjs:data');

	if(Object.keys(this.dom).length)
		this.node.setAttribute('svgjs:data', JSON.stringify(this.dom)); // see #428

	return this;
}
  // set given data to the elements data property
  , setData: function(o){
	this.dom = o;
	return this;
}
  , is: function(obj){
	return is(this, obj);
}
}
		});

		SVG.easing = {
			'-': function(pos){return pos;}
, '<>':function(pos){return -Math.cos(pos * Math.PI) / 2 + 0.5;}
, '>': function(pos){return  Math.sin(pos * Math.PI / 2);}
, '<': function(pos){return -Math.cos(pos * Math.PI / 2) + 1;}
		};

		SVG.morph = function(pos){
			return function(from, to) {
				return new SVG.MorphObj(from, to).at(pos);
			};
		};

		SVG.Situation = SVG.invent({

			create: function(o){
				this.init = false;
				this.reversed = false;
				this.reversing = false;

				this.duration = new SVG.Number(o.duration).valueOf();
				this.delay = new SVG.Number(o.delay).valueOf();

				this.start = +new Date() + this.delay;
				this.finish = this.start + this.duration;
				this.ease = o.ease;

    // this.loop is incremented from 0 to this.loops
    // it is also incremented when in an infinite loop (when this.loops is true)
				this.loop = 0;
				this.loops = false;

				this.animations = {
      // functionToCall: [list of morphable objects]
      // e.g. move: [SVG.Number, SVG.Number]
				};

				this.attrs = {
      // holds all attributes which are not represented from a function svg.js provides
      // e.g. someAttr: SVG.Number
				};

				this.styles = {
      // holds all styles which should be animated
      // e.g. fill-color: SVG.Color
				};

				this.transforms = [
      // holds all transformations as transformation objects
      // e.g. [SVG.Rotate, SVG.Translate, SVG.Matrix]
				];

				this.once = {
      // functions to fire at a specific position
      // e.g. "0.5": function foo(){}
				};

			}

		});


		SVG.FX = SVG.invent({

			create: function(element) {
				this._target = element;
				this.situations = [];
				this.active = false;
				this.situation = null;
				this.paused = false;
				this.lastPos = 0;
				this.pos = 0;
    // The absolute position of an animation is its position in the context of its complete duration (including delay and loops)
    // When performing a delay, absPos is below 0 and when performing a loop, its value is above 1
				this.absPos = 0;
				this._speed = 1;
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
			ease = o.ease;
			delay = o.delay;
			o = o.duration;
		}

		var situation = new SVG.Situation({
			duration: o || 1000,
			delay: delay || 0,
			ease: SVG.easing[ease || '-'] || ease
		});

		this.queue(situation);

		return this;
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
	});

	return this.queue(situation);
}

    /**
     * sets or returns the target of this animation
     * @param null || target SVG.Element which should be set as new target
     * @return target || this
     */
  , target: function(target){
	if(target && target instanceof SVG.Element){
		this._target = target;
		return this;
	}

	return this._target;
}

    // returns the absolute position at a given time
  , timeToAbsPos: function(timestamp){
	return (timestamp - this.situation.start) / (this.situation.duration/this._speed);
}

    // returns the timestamp from a given absolute positon
  , absPosToTime: function(absPos){
	return this.situation.duration/this._speed * absPos + this.situation.start;
}

    // starts the animationloop
  , startAnimFrame: function(){
	this.stopAnimFrame();
	this.animationFrame = requestAnimationFrame(function(){ this.step(); }.bind(this));
}

    // cancels the animationframe
  , stopAnimFrame: function(){
	cancelAnimationFrame(this.animationFrame);
}

    // kicks off the animation - only does something when the queue is currently not active and at least one situation is set
  , start: function(){
      // dont start if already started
	if(!this.active && this.situation){
		this.active = true;
		this.startCurrent();
	}

	return this;
}

    // start the current situation
  , startCurrent: function(){
	this.situation.start = +new Date + this.situation.delay/this._speed;
	this.situation.finish = this.situation.start + this.situation.duration/this._speed;
	return this.initAnimations().step();
}

    /**
     * adds a function / Situation to the animation queue
     * @param fn function / situation to add
     * @return this
     */
  , queue: function(fn){
	if(typeof fn == 'function' || fn instanceof SVG.Situation)
		this.situations.push(fn);

	if(!this.situation) this.situation = this.situations.shift();

	return this;
}

    /**
     * pulls next element from the queue and execute it
     * @return this
     */
  , dequeue: function(){
      // stop current animation
	this.situation && this.situation.stop && this.situation.stop();

      // get next animation from queue
	this.situation = this.situations.shift();

	if(this.situation){
		if(this.situation instanceof SVG.Situation) {
			this.startCurrent();
		} else {
          // If it is not a SVG.Situation, then it is a function, we execute it
			this.situation.call(this);
		}
	}

	return this;
}

    // updates all animations to the current state of the element
    // this is important when one property could be changed from another property
  , initAnimations: function() {
	var i;
	var s = this.situation;

	if(s.init) return this;

	for(i in s.animations){

		if(i == 'viewbox'){
			s.animations[i] = this.target().viewbox().morph(s.animations[i]);
		}else{

          // TODO: this is not a clean clone of the array. We may have some unchecked references
			s.animations[i].value = (i == 'plot' ? this.target().array().value : this.target()[i]());

          // sometimes we get back an object and not the real value, fix this
			if(s.animations[i].value.value){
				s.animations[i].value = s.animations[i].value.value;
			}

			if(s.animations[i].relative)
				s.animations[i].destination.value = s.animations[i].destination.value + s.animations[i].value;

		}

	}

	for(i in s.attrs){
		if(s.attrs[i] instanceof SVG.Color){
			var color = new SVG.Color(this.target().attr(i));
			s.attrs[i].r = color.r;
			s.attrs[i].g = color.g;
			s.attrs[i].b = color.b;
		}else{
			s.attrs[i].value = this.target().attr(i);// + s.attrs[i].value
		}
	}

	for(i in s.styles){
		s.styles[i].value = this.target().style(i);
	}

	s.initialTransformation = this.target().matrixify();

	s.init = true;
	return this;
}
  , clearQueue: function(){
	this.situations = [];
	return this;
}
  , clearCurrent: function(){
	this.situation = null;
	return this;
}
    /** stops the animation immediately
     * @param jumpToEnd A Boolean indicating whether to complete the current animation immediately.
     * @param clearQueue A Boolean indicating whether to remove queued animation as well.
     * @return this
     */
  , stop: function(jumpToEnd, clearQueue){
	if(!this.active) this.start();

	if(clearQueue){
		this.clearQueue();
	}

	this.active = false;

	if(jumpToEnd && this.situation){
		this.atEnd();
	}

	this.stopAnimFrame();

	return this.clearCurrent();
}

    /** resets the element to the state where the current element has started
     * @return this
     */
  , reset: function(){
	if(this.situation){
		var temp = this.situation;
		this.stop();
		this.situation = temp;
		this.atStart();
	}
	return this;
}

    // Stop the currently-running animation, remove all queued animations, and complete all animations for the element.
  , finish: function(){

	this.stop(true, false);

	while(this.dequeue().situation && this.stop(true, false));

	this.clearQueue().clearCurrent();

	return this;
}

    // set the internal animation pointer at the start position, before any loops, and updates the visualisation
  , atStart: function() {
	return this.at(0, true);
}

    // set the internal animation pointer at the end position, after all the loops, and updates the visualisation
  , atEnd: function() {
	if (this.situation.loops === true) {
      // If in a infinite loop, we end the current iteration
		return this.at(this.situation.loop+1, true);
	} else if(typeof this.situation.loops == 'number') {
      // If performing a finite number of loops, we go after all the loops
		return this.at(this.situation.loops, true);
	} else {
      // If no loops, we just go at the end
		return this.at(1, true);
	}
}

    // set the internal animation pointer to the specified position and updates the visualisation
    // if isAbsPos is true, pos is treated as an absolute position
  , at: function(pos, isAbsPos){
	var durDivSpd = this.situation.duration/this._speed;

	this.absPos = pos;
      // If pos is not an absolute position, we convert it into one
	if (!isAbsPos) {
		if (this.situation.reversed) this.absPos = 1 - this.absPos;
		this.absPos += this.situation.loop;
	}

	this.situation.start = +new Date - this.absPos * durDivSpd;
	this.situation.finish = this.situation.start + durDivSpd;

	return this.step(true);
}

    /**
     * sets or returns the speed of the animations
     * @param speed null || Number The new speed of the animations
     * @return Number || this
     */
  , speed: function(speed){
	if (speed === 0) return this.pause();

	if (speed) {
		this._speed = speed;
        // We use an absolute position here so that speed can affect the delay before the animation
		return this.at(this.absPos, true);
	} else return this._speed;
}

    // Make loopable
  , loop: function(times, reverse) {
	var c = this.last();

      // store total loops
	c.loops = (times != null) ? times : true;
	c.loop = 0;

	if(reverse) c.reversing = true;
	return this;
}

    // pauses the animation
  , pause: function(){
	this.paused = true;
	this.stopAnimFrame();

	return this;
}

    // unpause the animation
  , play: function(){
	if(!this.paused) return this;
	this.paused = false;
      // We use an absolute position here so that the delay before the animation can be paused
	return this.at(this.absPos, true);
}

    /**
     * toggle or set the direction of the animation
     * true sets direction to backwards while false sets it to forwards
     * @param reversed Boolean indicating whether to reverse the animation or not (default: toggle the reverse status)
     * @return this
     */
  , reverse: function(reversed){
	var c = this.last();

	if(typeof reversed == 'undefined') c.reversed = !c.reversed;
	else c.reversed = reversed;

	return this;
}


    /**
     * returns a float from 0-1 indicating the progress of the current animation
     * @param eased Boolean indicating whether the returned position should be eased or not
     * @return number
     */
  , progress: function(easeIt){
	return easeIt ? this.situation.ease(this.pos) : this.pos;
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
		fn.call(this, c);
		this.off('finished.fx', wrapper); // prevent memory leak
	}
};

	this.target().on('finished.fx', wrapper);
	return this;
}

    // adds a callback which is called whenever one animation step is performed
  , during: function(fn){
	var c = this.last()
        , wrapper = function(e){
	if(e.detail.situation == c){
		fn.call(this, e.detail.pos, SVG.morph(e.detail.pos), e.detail.eased, c);
	}
};

      // see above
	this.target().off('during.fx', wrapper).on('during.fx', wrapper);

	return this.after(function(){
		this.off('during.fx', wrapper);
	});
}

    // calls after ALL animations in the queue are finished
  , afterAll: function(fn){
	var wrapper = function wrapper(e){
		fn.call(this);
		this.off('allfinished.fx', wrapper);
	};

      // see above
	this.target().off('allfinished.fx', wrapper).on('allfinished.fx', wrapper);
	return this;
}

    // calls on every animation step for all animations
  , duringAll: function(fn){
	var wrapper = function(e){
		fn.call(this, e.detail.pos, SVG.morph(e.detail.pos), e.detail.eased, e.detail.situation);
	};

	this.target().off('during.fx', wrapper).on('during.fx', wrapper);

	return this.afterAll(function(){
		this.off('during.fx', wrapper);
	});
}

  , last: function(){
	return this.situations.length ? this.situations[this.situations.length-1] : this.situation;
}

    // adds one property to the animations
  , add: function(method, args, type){
	this.last()[type || 'animations'][method] = args;
	setTimeout(function(){this.start();}.bind(this), 0);
	return this;
}

    /** perform one step of the animation
     *  @param ignoreTime Boolean indicating whether to ignore time and use position directly or recalculate position based on time
     *  @return this
     */
  , step: function(ignoreTime){

      // convert current time to an absolute position
	if(!ignoreTime) this.absPos = this.timeToAbsPos(+new Date);

      // This part convert an absolute position to a position
	if(this.situation.loops !== false) {
		var absPos, absPosInt, lastLoop;

        // If the absolute position is below 0, we just treat it as if it was 0
		absPos = Math.max(this.absPos, 0);
		absPosInt = Math.floor(absPos);

		if(this.situation.loops === true || absPosInt < this.situation.loops) {
			this.pos = absPos - absPosInt;
			lastLoop = this.situation.loop;
			this.situation.loop = absPosInt;
		} else {
			this.absPos = this.situation.loops;
			this.pos = 1;
          // The -1 here is because we don't want to toggle reversed when all the loops have been completed
			lastLoop = this.situation.loop - 1;
			this.situation.loop = this.situation.loops;
		}

		if(this.situation.reversing) {
          // Toggle reversed if an odd number of loops as occured since the last call of step
			this.situation.reversed = this.situation.reversed != Boolean((this.situation.loop - lastLoop) % 2);
		}

	} else {
        // If there are no loop, the absolute position must not be above 1
		this.absPos = Math.min(this.absPos, 1);
		this.pos = this.absPos;
	}

      // while the absolute position can be below 0, the position must not be below 0
	if(this.pos < 0) this.pos = 0;

	if(this.situation.reversed) this.pos = 1 - this.pos;


      // apply easing
	var eased = this.situation.ease(this.pos);

      // call once-callbacks
	for(var i in this.situation.once){
		if(i > this.lastPos && i <= eased){
			this.situation.once[i].call(this.target(), this.pos, eased);
			delete this.situation.once[i];
		}
	}

      // fire during callback with position, eased position and current situation as parameter
	if(this.active) this.target().fire('during', {pos: this.pos, eased: eased, fx: this, situation: this.situation});

      // the user may call stop or finish in the during callback
      // so make sure that we still have a valid situation
	if(!this.situation){
		return this;
	}

      // apply the actual animation to every property
	this.eachAt();

      // do final code when situation is finished
	if((this.pos == 1 && !this.situation.reversed) || (this.situation.reversed && this.pos == 0)){

        // stop animation callback
		this.stopAnimFrame();

        // fire finished callback with current situation as parameter
		this.target().fire('finished', {fx:this, situation: this.situation});

		if(!this.situations.length){
			this.target().fire('allfinished');
			this.target().off('.fx'); // there shouldnt be any binding left, but to make sure...
			this.active = false;
		}

        // start next animation
		if(this.active) this.dequeue();
		else this.clearCurrent();

	}else if(!this.paused && this.active){
        // we continue animating when we are not at the end
		this.startAnimFrame();
	}

      // save last eased position for once callback triggering
	this.lastPos = eased;
	return this;

}

    // calculates the step for every property and calls block with it
  , eachAt: function(){
	var i, at, self = this, target = this.target(), s = this.situation;

      // apply animations which can be called trough a method
	for(i in s.animations){

		at = [].concat(s.animations[i]).map(function(el){
			return typeof el !== 'string' && el.at ? el.at(s.ease(self.pos), self.pos) : el;
		});

		target[i].apply(target, at);

	}

      // apply animation which has to be applied with attr()
	for(i in s.attrs){

		at = [i].concat(s.attrs[i]).map(function(el){
			return typeof el !== 'string' && el.at ? el.at(s.ease(self.pos), self.pos) : el;
		});

		target.attr.apply(target, at);

	}

      // apply animation which has to be applied with style()
	for(i in s.styles){

		at = [i].concat(s.styles[i]).map(function(el){
			return typeof el !== 'string' && el.at ? el.at(s.ease(self.pos), self.pos) : el;
		});

		target.style.apply(target, at);

	}

      // animate initialTransformation which has to be chained
	if(s.transforms.length){

        // get initial initialTransformation
		at = s.initialTransformation;
		for(i = 0, len = s.transforms.length; i < len; i++){

          // get next transformation in chain
			var a = s.transforms[i];

          // multiply matrix directly
			if(a instanceof SVG.Matrix){

				if(a.relative){
					at = at.multiply(new SVG.Matrix().morph(a).at(s.ease(this.pos)));
				}else{
					at = at.morph(a).at(s.ease(this.pos));
				}
				continue;
			}

          // when transformation is absolute we have to reset the needed transformation first
			if(!a.relative)
				a.undo(at.extract());

          // and reapply it after
			at = at.multiply(a.at(s.ease(this.pos)));

		}

        // set new matrix on element
		target.matrix(at);
	}

	return this;

}


    // adds an once-callback which is called at a specific position and never again
  , once: function(pos, fn, isEased){

	if(!isEased)pos = this.situation.ease(pos);

	this.situation.once[pos] = fn;

	return this;
}

}

, parent: SVG.Element

  // Add method to parent elements
, construct: {
    // Get fx module or create a new one, then animate with given duration and ease
	animate: function(o, ease, delay) {
		return (this.fx || (this.fx = new SVG.FX(this))).animate(o, ease, delay);
	}
  , delay: function(delay){
	return (this.fx || (this.fx = new SVG.FX(this))).delay(delay);
}
  , stop: function(jumpToEnd, clearQueue) {
	if (this.fx)
		this.fx.stop(jumpToEnd, clearQueue);

	return this;
}
  , finish: function() {
	if (this.fx)
		this.fx.finish();

	return this;
}
    // Pause current animation
  , pause: function() {
	if (this.fx)
		this.fx.pause();

	return this;
}
    // Play paused current animation
  , play: function() {
	if (this.fx)
		this.fx.play();

	return this;
}
    // Set/Get the speed of the animations
  , speed: function(speed) {
	if (this.fx)
		if (speed == null)
			return this.fx.speed();
		else
          this.fx.speed(speed);

	return this;
}
}

		});

// MorphObj is used whenever no morphable object is given
		SVG.MorphObj = SVG.invent({

			create: function(from, to){
    // prepare color for morphing
				if(SVG.Color.isColor(to)) return new SVG.Color(from).morph(to);
    // prepare number for morphing
				if(SVG.regex.numberAndUnit.test(to)) return new SVG.Number(from).morph(to);

    // prepare for plain morphing
				this.value = 0;
				this.destination = to;
			}

, extend: {
	at: function(pos, real){
		return real < 1 ? this.value : this.destination;
	},

	valueOf: function(){
		return this.value;
	}
}

		});

		SVG.extend(SVG.FX, {
  // Add animatable attributes
			attr: function(a, v, relative) {
    // apply attributes individually
				if (typeof a == 'object') {
					for (var key in a)
						this.attr(key, a[key]);

				} else {
      // the MorphObj takes care about the right function used
					this.add(a, new SVG.MorphObj(null, v), 'attrs');
				}

				return this;
			}
  // Add animatable styles
, style: function(s, v) {
	if (typeof s == 'object')
		for (var key in s)
			this.style(key, s[key]);

	else
      this.add(s, new SVG.MorphObj(null, v), 'styles');

	return this;
}
  // Animatable x-axis
, x: function(x, relative) {
	if(this.target() instanceof SVG.G){
		this.transform({x:x}, relative);
		return this;
	}

	var num = new SVG.Number().morph(x);
	num.relative = relative;
	return this.add('x', num);
}
  // Animatable y-axis
, y: function(y, relative) {
	if(this.target() instanceof SVG.G){
		this.transform({y:y}, relative);
		return this;
	}

	var num = new SVG.Number().morph(y);
	num.relative = relative;
	return this.add('y', num);
}
  // Animatable center x-axis
, cx: function(x) {
	return this.add('cx', new SVG.Number().morph(x));
}
  // Animatable center y-axis
, cy: function(y) {
	return this.add('cy', new SVG.Number().morph(y));
}
  // Add animatable move
, move: function(x, y) {
	return this.x(x).y(y);
}
  // Add animatable center
, center: function(x, y) {
	return this.cx(x).cy(y);
}
  // Add animatable size
, size: function(width, height) {
	if (this.target() instanceof SVG.Text) {
      // animate font size for Text elements
		this.attr('font-size', width);

	} else {
      // animate bbox based size for all other elements
		var box;

		if(!width || !height){
			box = this.target().bbox();
		}

		if(!width){
			width = box.width / box.height  * height;
		}

		if(!height){
			height = box.height / box.width  * width;
		}

		this.add('width' , new SVG.Number().morph(width))
          .add('height', new SVG.Number().morph(height));

	}

	return this;
}
  // Add animatable plot
, plot: function(p) {
	return this.add('plot', this.target().array().morph(p));
}
  // Add leading method
, leading: function(value) {
	return this.target().leading ?
      this.add('leading', new SVG.Number().morph(value)) :
      this;
}
  // Add animatable viewbox
, viewbox: function(x, y, width, height) {
	if (this.target() instanceof SVG.Container) {
		this.add('viewbox', new SVG.ViewBox(x, y, width, height));
	}

	return this;
}
, update: function(o) {
	if (this.target() instanceof SVG.Stop) {
		if (typeof o == 'number' || o instanceof SVG.Number) {
			return this.update({
				offset:  arguments[0]
        , color:   arguments[1]
        , opacity: arguments[2]
			});
		}

		if (o.opacity != null) this.attr('stop-opacity', o.opacity);
		if (o.color   != null) this.attr('stop-color', o.color);
		if (o.offset  != null) this.attr('offset', o.offset);
	}

	return this;
}
		});

		SVG.BBox = SVG.invent({
  // Initialize
			create: function(element) {
    // get values if element is given
				if (element) {
					var box;

      // yes this is ugly, but Firefox can be a bitch when it comes to elements that are not yet rendered
					try {

        // the element is NOT in the dom, throw error
						if(!document.documentElement.contains(element.node)) throw new Exception('Element not in the dom');

        // find native bbox
						box = element.node.getBBox();
					} catch(e) {
						if(element instanceof SVG.Shape){
							var clone = element.clone(SVG.parser.draw).show();
							box = clone.bbox();
							clone.remove();
						}else{
							box = {
								x:      element.node.clientLeft
          , y:      element.node.clientTop
          , width:  element.node.clientWidth
          , height: element.node.clientHeight
							};
						}
					}

      // plain x and y
					this.x = box.x;
					this.y = box.y;

      // plain width and height
					this.width  = box.width;
					this.height = box.height;
				}

    // add center, right and bottom
				fullBox(this);
			}

  // Define Parent
, parent: SVG.Element

  // Constructor
, construct: {
    // Get bounding box
	bbox: function() {
		return new SVG.BBox(this);
	}
}

		});

		SVG.TBox = SVG.invent({
  // Initialize
			create: function(element) {
    // get values if element is given
				if (element) {
					var t   = element.ctm().extract()
        , box = element.bbox();

      // width and height including transformations
					this.width  = box.width  * t.scaleX;
					this.height = box.height * t.scaleY;

      // x and y including transformations
					this.x = box.x + t.x;
					this.y = box.y + t.y;
				}

    // add center, right and bottom
				fullBox(this);
			}

  // Define Parent
, parent: SVG.Element

  // Constructor
, construct: {
    // Get transformed bounding box
	tbox: function() {
		return new SVG.TBox(this);
	}
}

		});


		SVG.RBox = SVG.invent({
  // Initialize
			create: function(element) {
				if (element) {
					var e    = element.doc().parent()
        , box  = element.node.getBoundingClientRect()
        , zoom = 1;

      // get screen offset
					this.x = box.left;
					this.y = box.top;

      // subtract parent offset
					this.x -= e.offsetLeft;
					this.y -= e.offsetTop;

					while (e = e.offsetParent) {
						this.x -= e.offsetLeft;
						this.y -= e.offsetTop;
					}

      // calculate cumulative zoom from svg documents
					e = element;
					while (e.parent && (e = e.parent())) {
						if (e.viewbox) {
							zoom *= e.viewbox().zoom;
							this.x -= e.x() || 0;
							this.y -= e.y() || 0;
						}
					}

      // recalculate viewbox distortion
					this.width  = box.width  /= zoom;
					this.height = box.height /= zoom;
				}

    // add center, right and bottom
				fullBox(this);

    // offset by window scroll position, because getBoundingClientRect changes when window is scrolled
				this.x += window.pageXOffset;
				this.y += window.pageYOffset;
			}

  // define Parent
, parent: SVG.Element

  // Constructor
, construct: {
    // Get rect box
	rbox: function() {
		return new SVG.RBox(this);
	}
}

		})

// Add universal merge method
;[SVG.BBox, SVG.TBox, SVG.RBox].forEach(function(c) {

	SVG.extend(c, {
    // Merge rect box with another, return a new instance
		merge: function(box) {
			var b = new c();

      // merge boxes
			b.x      = Math.min(this.x, box.x);
			b.y      = Math.min(this.y, box.y);
			b.width  = Math.max(this.x + this.width,  box.x + box.width)  - b.x;
			b.height = Math.max(this.y + this.height, box.y + box.height) - b.y;

			return fullBox(b);
		}

	});

});

		SVG.Matrix = SVG.invent({
  // Initialize
			create: function(source) {
				var i, base = arrayToMatrix([1, 0, 0, 1, 0, 0]);

    // ensure source as object
				source = source instanceof SVG.Element ?
      source.matrixify() :
    typeof source === 'string' ?
      stringToMatrix(source) :
    arguments.length == 6 ?
      arrayToMatrix([].slice.call(arguments)) :
    typeof source === 'object' ?
      source : base;

    // merge source
				for (i = abcdef.length - 1; i >= 0; --i)
					this[abcdef[i]] = source && typeof source[abcdef[i]] === 'number' ?
        source[abcdef[i]] : base[abcdef[i]];
			}

  // Add methods
, extend: {
    // Extract individual transformations
	extract: function() {
      // find delta transform points
		var px    = deltaTransformPoint(this, 0, 1)
        , py    = deltaTransformPoint(this, 1, 0)
        , skewX = 180 / Math.PI * Math.atan2(px.y, px.x) - 90;

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
		};
	}
    // Clone matrix
  , clone: function() {
	return new SVG.Matrix(this);
}
    // Morph one matrix into another
  , morph: function(matrix) {
      // store new destination
	this.destination = new SVG.Matrix(matrix);

	return this;
}
    // Get morphed matrix at a given position
  , at: function(pos) {
      // make sure a destination is defined
	if (!this.destination) return this;

      // calculate morphed matrix at a given position
	var matrix = new SVG.Matrix({
		a: this.a + (this.destination.a - this.a) * pos
      , b: this.b + (this.destination.b - this.b) * pos
      , c: this.c + (this.destination.c - this.c) * pos
      , d: this.d + (this.destination.d - this.d) * pos
      , e: this.e + (this.destination.e - this.e) * pos
      , f: this.f + (this.destination.f - this.f) * pos
	});

      // process parametric rotation if present
	if (this.param && this.param.to) {
        // calculate current parametric position
		var param = {
			rotation: this.param.from.rotation + (this.param.to.rotation - this.param.from.rotation) * pos
        , cx:       this.param.from.cx
        , cy:       this.param.from.cy
		};

        // rotate matrix
		matrix = matrix.rotate(
          (this.param.to.rotation - this.param.from.rotation * 2) * pos
        , param.cx
        , param.cy
        );

        // store current parametric values
		matrix.param = param;
	}

	return matrix;
}
    // Multiplies by given matrix
  , multiply: function(matrix) {
	return new SVG.Matrix(this.native().multiply(parseMatrix(matrix).native()));
}
    // Inverses matrix
  , inverse: function() {
	return new SVG.Matrix(this.native().inverse());
}
    // Translate matrix
  , translate: function(x, y) {
	return new SVG.Matrix(this.native().translate(x || 0, y || 0));
}
    // Scale matrix
  , scale: function(x, y, cx, cy) {
      // support uniformal scale
	if (arguments.length == 1) {
		y = x;
	} else if (arguments.length == 3) {
		cy = cx;
		cx = y;
		y = x;
	}

	return this.around(cx, cy, new SVG.Matrix(x, 0, 0, y, 0, 0));
}
    // Rotate matrix
  , rotate: function(r, cx, cy) {
      // convert degrees to radians
	r = SVG.utils.radians(r);

	return this.around(cx, cy, new SVG.Matrix(Math.cos(r), Math.sin(r), -Math.sin(r), Math.cos(r), 0, 0));
}
    // Flip matrix on x or y, at a given offset
  , flip: function(a, o) {
	return a == 'x' ? this.scale(-1, 1, o, 0) : this.scale(1, -1, 0, o);
}
    // Skew
  , skew: function(x, y, cx, cy) {
      // support uniformal skew
	if (arguments.length == 1) {
		y = x;
	} else if (arguments.length == 3) {
		cy = cx;
		cx = y;
		y = x;
	}

      // convert degrees to radians
	x = SVG.utils.radians(x);
	y = SVG.utils.radians(y);

	return this.around(cx, cy, new SVG.Matrix(1, Math.tan(y), Math.tan(x), 1, 0, 0));
}
    // SkewX
  , skewX: function(x, cx, cy) {
	return this.skew(x, 0, cx, cy);
}
    // SkewY
  , skewY: function(y, cx, cy) {
	return this.skew(0, y, cx, cy);
}
    // Transform around a center point
  , around: function(cx, cy, matrix) {
	return this
        .multiply(new SVG.Matrix(1, 0, 0, 1, cx || 0, cy || 0))
        .multiply(matrix)
        .multiply(new SVG.Matrix(1, 0, 0, 1, -cx || 0, -cy || 0));
}
    // Convert to native SVGMatrix
  , native: function() {
      // create new matrix
	var matrix = SVG.parser.native.createSVGMatrix();

      // update with current values
	for (var i = abcdef.length - 1; i >= 0; i--)
		matrix[abcdef[i]] = this[abcdef[i]];

	return matrix;
}
    // Convert matrix to string
  , toString: function() {
	return 'matrix(' + this.a + ',' + this.b + ',' + this.c + ',' + this.d + ',' + this.e + ',' + this.f + ')';
}
}

  // Define parent
, parent: SVG.Element

  // Add parent method
, construct: {
    // Get current matrix
	ctm: function() {
		return new SVG.Matrix(this.node.getCTM());
	},
    // Get current screen matrix
	screenCTM: function() {
		return new SVG.Matrix(this.node.getScreenCTM());
	}

}

		});

		SVG.Point = SVG.invent({
  // Initialize
			create: function(x,y) {
				var i, source, base = {x:0, y:0};

    // ensure source as object
				source = Array.isArray(x) ?
      {x:x[0], y:x[1]} :
    typeof x === 'object' ?
      {x:x.x, y:x.y} :
    x != null ?
      {x:x, y:(y != null ? y : x)} : base; // If y has no value, then x is used has its value

    // merge source
				this.x = source.x;
				this.y = source.y;
			}

  // Add methods
, extend: {
    // Clone point
	clone: function() {
		return new SVG.Point(this);
	}
    // Morph one point into another
  , morph: function(x, y) {
      // store new destination
	this.destination = new SVG.Point(x, y);

	return this;
}
    // Get morphed point at a given position
  , at: function(pos) {
      // make sure a destination is defined
	if (!this.destination) return this;

      // calculate morphed matrix at a given position
	var point = new SVG.Point({
		x: this.x + (this.destination.x - this.x) * pos
      , y: this.y + (this.destination.y - this.y) * pos
	});

	return point;
}
    // Convert to native SVGPoint
  , native: function() {
      // create new point
	var point = SVG.parser.native.createSVGPoint();

      // update with current values
	point.x = this.x;
	point.y = this.y;

	return point;
}
    // transform point with matrix
  , transform: function(matrix) {
	return new SVG.Point(this.native().matrixTransform(matrix.native()));
}

}

		});

		SVG.extend(SVG.Element, {

  // Get point
			point: function(x, y) {
				return new SVG.Point(x,y).transform(this.screenCTM().inverse());
			}

		});

		SVG.extend(SVG.Element, {
  // Set svg element attribute
			attr: function(a, v, n) {
    // act as full getter
				if (a == null) {
      // get an object of attributes
					a = {};
					v = this.node.attributes;
					for (n = v.length - 1; n >= 0; n--)
						a[v[n].nodeName] = SVG.regex.isNumber.test(v[n].nodeValue) ? parseFloat(v[n].nodeValue) : v[n].nodeValue;

					return a;

				} else if (typeof a == 'object') {
      // apply every attribute individually if an object is passed
					for (v in a) this.attr(v, a[v]);

				} else if (v === null) {
        // remove value
					this.node.removeAttribute(a);

				} else if (v == null) {
      // act as a getter if the first and only argument is not an object
					v = this.node.getAttribute(a);
					return v == null ?
        SVG.defaults.attrs[a] :
      SVG.regex.isNumber.test(v) ?
        parseFloat(v) : v;

				} else {
      // BUG FIX: some browsers will render a stroke if a color is given even though stroke width is 0
					if (a == 'stroke-width')
						this.attr('stroke', parseFloat(v) > 0 ? this._stroke : null);
					else if (a == 'stroke')
						this._stroke = v;

      // convert image fill and stroke to patterns
					if (a == 'fill' || a == 'stroke') {
						if (SVG.regex.isImage.test(v))
						v = this.doc().defs().image(v, 0, 0);

						if (v instanceof SVG.Image)
						v = this.doc().defs().pattern(0, 0, function() {
							this.add(v);
						});
					}

      // ensure correct numeric values (also accepts NaN and Infinity)
					if (typeof v === 'number')
						v = new SVG.Number(v);

      // ensure full hex color
					else if (SVG.Color.isColor(v))
						v = new SVG.Color(v);

      // parse array values
					else if (Array.isArray(v))
						v = new SVG.Array(v);

      // store parametric transformation values locally
					else if (v instanceof SVG.Matrix && v.param)
						this.param = v.param;

      // if the passed attribute is leading...
					if (a == 'leading') {
        // ... call the leading method instead
						if (this.leading)
						this.leading(v);
					} else {
        // set given attribute on node
						typeof n === 'string' ?
          this.node.setAttributeNS(n, a, v.toString()) :
          this.node.setAttribute(a, v.toString());
					}

      // rebuild if required
					if (this.rebuild && (a == 'font-size' || a == 'x'))
						this.rebuild(a, v);
				}

				return this;
			}
		});
		SVG.extend(SVG.Element, {
  // Add transformations
			transform: function(o, relative) {
    // get target in case of the fx module, otherwise reference this
				var target = this
      , matrix;

    // act as a getter
				if (typeof o !== 'object') {
      // get current matrix
					matrix = new SVG.Matrix(target).extract();

					return typeof o === 'string' ? matrix[o] : matrix;
				}

    // get current matrix
				matrix = new SVG.Matrix(target);

    // ensure relative flag
				relative = !!relative || !!o.relative;

    // act on matrix
				if (o.a != null) {
					matrix = relative ?
        // relative
        matrix.multiply(new SVG.Matrix(o)) :
        // absolute
        new SVG.Matrix(o);

    // act on rotation
				} else if (o.rotation != null) {
      // ensure centre point
					ensureCentre(o, target);

      // apply transformation
					matrix = relative ?
        // relative
        matrix.rotate(o.rotation, o.cx, o.cy) :
        // absolute
        matrix.rotate(o.rotation - matrix.extract().rotation, o.cx, o.cy);

    // act on scale
				} else if (o.scale != null || o.scaleX != null || o.scaleY != null) {
      // ensure centre point
					ensureCentre(o, target);

      // ensure scale values on both axes
					o.scaleX = o.scale != null ? o.scale : o.scaleX != null ? o.scaleX : 1;
					o.scaleY = o.scale != null ? o.scale : o.scaleY != null ? o.scaleY : 1;

					if (!relative) {
        // absolute; multiply inversed values
						var e = matrix.extract();
						o.scaleX = o.scaleX * 1 / e.scaleX;
						o.scaleY = o.scaleY * 1 / e.scaleY;
					}

					matrix = matrix.scale(o.scaleX, o.scaleY, o.cx, o.cy);

    // act on skew
				} else if (o.skew != null || o.skewX != null || o.skewY != null) {
      // ensure centre point
					ensureCentre(o, target);

      // ensure skew values on both axes
					o.skewX = o.skew != null ? o.skew : o.skewX != null ? o.skewX : 0;
					o.skewY = o.skew != null ? o.skew : o.skewY != null ? o.skewY : 0;

					if (!relative) {
        // absolute; reset skew values
						var e = matrix.extract();
						matrix = matrix.multiply(new SVG.Matrix().skew(e.skewX, e.skewY, o.cx, o.cy).inverse());
					}

					matrix = matrix.skew(o.skewX, o.skewY, o.cx, o.cy);

    // act on flip
				} else if (o.flip) {
					matrix = matrix.flip(
        o.flip
      , o.offset == null ? target.bbox()['c' + o.flip] : o.offset
      );

    // act on translate
				} else if (o.x != null || o.y != null) {
				if (relative) {
        // relative
					matrix = matrix.translate(o.x, o.y);
				} else {
        // absolute
					if (o.x != null) matrix.e = o.x;
					if (o.y != null) matrix.f = o.y;
				}
			}

				return this.attr('transform', matrix);
			}
		});

		SVG.extend(SVG.FX, {
			transform: function(o, relative) {
    // get target in case of the fx module, otherwise reference this
				var target = this.target()
      , matrix;

    // act as a getter
				if (typeof o !== 'object') {
      // get current matrix
					matrix = new SVG.Matrix(target).extract();

					return typeof o === 'string' ? matrix[o] : matrix;
				}

    // ensure relative flag
				relative = !!relative || !!o.relative;

    // act on matrix
				if (o.a != null) {
					matrix = new SVG.Matrix(o);

    // act on rotation
				} else if (o.rotation != null) {
      // ensure centre point
					ensureCentre(o, target);

      // apply transformation
					matrix = new SVG.Rotate(o.rotation, o.cx, o.cy);

    // act on scale
				} else if (o.scale != null || o.scaleX != null || o.scaleY != null) {
      // ensure centre point
					ensureCentre(o, target);

      // ensure scale values on both axes
					o.scaleX = o.scale != null ? o.scale : o.scaleX != null ? o.scaleX : 1;
					o.scaleY = o.scale != null ? o.scale : o.scaleY != null ? o.scaleY : 1;

					matrix = new SVG.Scale(o.scaleX, o.scaleY, o.cx, o.cy);

    // act on skew
				} else if (o.skewX != null || o.skewY != null) {
      // ensure centre point
					ensureCentre(o, target);

      // ensure skew values on both axes
					o.skewX = o.skewX != null ? o.skewX : 0;
					o.skewY = o.skewY != null ? o.skewY : 0;

					matrix = new SVG.Skew(o.skewX, o.skewY, o.cx, o.cy);

    // act on flip
				} else if (o.flip) {
					matrix = new SVG.Matrix().morph(new SVG.Matrix().flip(
        o.flip
      , o.offset == null ? target.bbox()['c' + o.flip] : o.offset
      ));

    // act on translate
				} else if (o.x != null || o.y != null) {
				matrix = new SVG.Translate(o.x, o.y);
			}

				if(!matrix) return this;

				matrix.relative = relative;

				this.last().transforms.push(matrix);

				setTimeout(function(){this.start();}.bind(this), 0);

				return this;
			}
		});

		SVG.extend(SVG.Element, {
  // Reset all transformations
			untransform: function() {
				return this.attr('transform', null);
			},
  // merge the whole transformation chain into one matrix and returns it
			matrixify: function() {

				var matrix = (this.attr('transform') || '')
      // split transformations
      .split(/\)\s*,?\s*/).slice(0,-1).map(function(str){
        // generate key => value pairs
	var kv = str.trim().split('(');
	return [kv[0], kv[1].split(SVG.regex.matrixElements).map(function(str){ return parseFloat(str); })];
})
      // calculate every transformation into one matrix
      .reduce(function(matrix, transform){

	if(transform[0] == 'matrix') return matrix.multiply(arrayToMatrix(transform[1]));
	return matrix[transform[0]].apply(matrix, transform[1]);

}, new SVG.Matrix());

				return matrix;
			},
  // add an element to another parent without changing the visual representation on the screen
			toParent: function(parent) {
				if(this == parent) return this;
				var ctm = this.screenCTM();
				var temp = parent.rect(1,1);
				var pCtm = temp.screenCTM().inverse();
				temp.remove();

				this.addTo(parent).untransform().transform(pCtm.multiply(ctm));

				return this;
			},
  // same as above with parent equals root-svg
			toDoc: function() {
				return this.toParent(this.doc());
			}

		});

		SVG.Transformation = SVG.invent({

			create: function(source, inversed){

				if(arguments.length > 1 && typeof inversed != 'boolean'){
					return this.create([].slice.call(arguments));
				}

				if(typeof source == 'object'){
					for(var i = 0, len = this.arguments.length; i < len; ++i){
						this[this.arguments[i]] = source[this.arguments[i]];
					}
				}

				if(Array.isArray(source)){
					for(var i = 0, len = this.arguments.length; i < len; ++i){
						this[this.arguments[i]] = source[i];
					}
				}

				this.inversed = false;

				if(inversed === true){
					this.inversed = true;
				}

			}

, extend: {

	at: function(pos){

		var params = [];

		for(var i = 0, len = this.arguments.length; i < len; ++i){
			params.push(this[this.arguments[i]]);
		}

		var m = this._undo || new SVG.Matrix();

		m = new SVG.Matrix().morph(SVG.Matrix.prototype[this.method].apply(m, params)).at(pos);

		return this.inversed ? m.inverse() : m;

	}

  , undo: function(o){
	for(var i = 0, len = this.arguments.length; i < len; ++i){
		o[this.arguments[i]] = typeof this[this.arguments[i]] == 'undefined' ? 0 : o[this.arguments[i]];
	}

      // The method SVG.Matrix.extract which was used before calling this
      // method to obtain a value for the parameter o doesn't return a cx and
      // a cy so we use the ones that were provided to this object at its creation
	o.cx = this.cx;
	o.cy = this.cy;

	this._undo = new SVG[capitalize(this.method)](o, true).at(1);

	return this;
}

}

		});

		SVG.Translate = SVG.invent({

			parent: SVG.Matrix
, inherit: SVG.Transformation

, create: function(source, inversed){
	if(typeof source == 'object') this.constructor.call(this, source, inversed);
	else this.constructor.call(this, [].slice.call(arguments));
}

, extend: {
	arguments: ['transformedX', 'transformedY']
  , method: 'translate'
}

		});

		SVG.Rotate = SVG.invent({

			parent: SVG.Matrix
, inherit: SVG.Transformation

, create: function(source, inversed){
	if(typeof source == 'object') this.constructor.call(this, source, inversed);
	else this.constructor.call(this, [].slice.call(arguments));
}

, extend: {
	arguments: ['rotation', 'cx', 'cy']
  , method: 'rotate'
  , at: function(pos){
	var m = new SVG.Matrix().rotate(new SVG.Number().morph(this.rotation - (this._undo ? this._undo.rotation : 0)).at(pos), this.cx, this.cy);
	return this.inversed ? m.inverse() : m;
}
  , undo: function(o){
	this._undo = o;
}
}

		});

		SVG.Scale = SVG.invent({

			parent: SVG.Matrix
, inherit: SVG.Transformation

, create: function(source, inversed){
	if(typeof source == 'object') this.constructor.call(this, source, inversed);
	else this.constructor.call(this, [].slice.call(arguments));
}

, extend: {
	arguments: ['scaleX', 'scaleY', 'cx', 'cy']
  , method: 'scale'
}

		});

		SVG.Skew = SVG.invent({

			parent: SVG.Matrix
, inherit: SVG.Transformation

, create: function(source, inversed){
	if(typeof source == 'object') this.constructor.call(this, source, inversed);
	else this.constructor.call(this, [].slice.call(arguments));
}

, extend: {
	arguments: ['skewX', 'skewY', 'cx', 'cy']
  , method: 'skew'
}

		});

		SVG.extend(SVG.Element, {
  // Dynamic style generator
			style: function(s, v) {
				if (arguments.length == 0) {
      // get full style
					return this.node.style.cssText || '';

				} else if (arguments.length < 2) {
      // apply every style individually if an object is passed
					if (typeof s == 'object') {
						for (v in s) this.style(v, s[v]);

					} else if (SVG.regex.isCss.test(s)) {
        // parse css string
						s = s.split(';');

        // apply every definition individually
						for (var i = 0; i < s.length; i++) {
							v = s[i].split(':');
							this.style(v[0].replace(/\s+/g, ''), v[1]);
						}
					} else {
        // act as a getter if the first and only argument is not an object
						return this.node.style[camelCase(s)];
					}

				} else {
					this.node.style[camelCase(s)] = v === null || SVG.regex.isBlank.test(v) ? '' : v;
				}

				return this;
			}
		});
		SVG.Parent = SVG.invent({
  // Initialize node
			create: function(element) {
				this.constructor.call(this, element);
			}

  // Inherit from
, inherit: SVG.Element

  // Add class methods
, extend: {
    // Returns all child elements
	children: function() {
		return SVG.utils.map(SVG.utils.filterSVGElements(this.node.childNodes), function(node) {
			return SVG.adopt(node);
		});
	}
    // Add given element at a position
  , add: function(element, i) {
	if (i == null)
		this.node.appendChild(element.node);
	else if (element.node != this.node.childNodes[i])
		this.node.insertBefore(element.node, this.node.childNodes[i]);

	return this;
}
    // Basically does the same as `add()` but returns the added element instead
  , put: function(element, i) {
	this.add(element, i);
	return element;
}
    // Checks if the given element is a child
  , has: function(element) {
	return this.index(element) >= 0;
}
    // Gets index of given element
  , index: function(element) {
	return [].slice.call(this.node.childNodes).indexOf(element.node);
}
    // Get a element at the given index
  , get: function(i) {
	return SVG.adopt(this.node.childNodes[i]);
}
    // Get first child
  , first: function() {
	return this.get(0);
}
    // Get the last child
  , last: function() {
	return this.get(this.node.childNodes.length - 1);
}
    // Iterates over all children and invokes a given block
  , each: function(block, deep) {
	var i, il
        , children = this.children();

	for (i = 0, il = children.length; i < il; i++) {
		if (children[i] instanceof SVG.Element)
			block.apply(children[i], [i, children]);

		if (deep && (children[i] instanceof SVG.Container))
			children[i].each(block, deep);
	}

	return this;
}
    // Remove a given child
  , removeElement: function(element) {
	this.node.removeChild(element.node);

	return this;
}
    // Remove all elements in this container
  , clear: function() {
      // remove children
	while(this.node.hasChildNodes())
		this.node.removeChild(this.node.lastChild);

      // remove defs reference
	delete this._defs;

	return this;
}
  , // Get defs
	defs: function() {
		return this.doc().defs();
	}
}

		});

		SVG.extend(SVG.Parent, {

			ungroup: function(parent, depth) {
				if(depth === 0 || this instanceof SVG.Defs) return this;

				parent = parent || (this instanceof SVG.Doc ? this : this.parent(SVG.Parent));
				depth = depth || Infinity;

				this.each(function(){
					if(this instanceof SVG.Defs) return this;
					if(this instanceof SVG.Parent) return this.ungroup(parent, depth-1);
					return this.toParent(parent);
				});

				this.node.firstChild || this.remove();

				return this;
			},

			flatten: function(parent, depth) {
				return this.ungroup(parent, depth);
			}

		});
		SVG.Container = SVG.invent({
  // Initialize node
			create: function(element) {
				this.constructor.call(this, element);
			}

  // Inherit from
, inherit: SVG.Parent

		});

		SVG.ViewBox = SVG.invent({

			create: function(source) {
				var i, base = [0, 0, 0, 0];

				var x, y, width, height, box, view, we, he
      , wm   = 1 // width multiplier
      , hm   = 1 // height multiplier
      , reg  = /[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:e[+-]?\d+)?/gi;

				if(source instanceof SVG.Element){

					we = source;
					he = source;
					view = (source.attr('viewBox') || '').match(reg);
					box = source.bbox;

      // get dimensions of current node
					width  = new SVG.Number(source.width());
					height = new SVG.Number(source.height());

      // find nearest non-percentual dimensions
					while (width.unit == '%') {
						wm *= width.value;
						width = new SVG.Number(we instanceof SVG.Doc ? we.parent().offsetWidth : we.parent().width());
						we = we.parent();
					}
					while (height.unit == '%') {
						hm *= height.value;
						height = new SVG.Number(he instanceof SVG.Doc ? he.parent().offsetHeight : he.parent().height());
						he = he.parent();
					}

      // ensure defaults
					this.x      = 0;
					this.y      = 0;
					this.width  = width  * wm;
					this.height = height * hm;
					this.zoom   = 1;

					if (view) {
        // get width and height from viewbox
						x      = parseFloat(view[0]);
						y      = parseFloat(view[1]);
						width  = parseFloat(view[2]);
						height = parseFloat(view[3]);

        // calculate zoom accoring to viewbox
						this.zoom = ((this.width / this.height) > (width / height)) ?
          this.height / height :
          this.width  / width;

        // calculate real pixel dimensions on parent SVG.Doc element
						this.x      = x;
						this.y      = y;
						this.width  = width;
						this.height = height;

					}

				}else{

      // ensure source as object
					source = typeof source === 'string' ?
        source.match(reg).map(function(el){ return parseFloat(el); }) :
      Array.isArray(source) ?
        source :
      typeof source == 'object' ?
        [source.x, source.y, source.width, source.height] :
      arguments.length == 4 ?
        [].slice.call(arguments) :
        base;

					this.x = source[0];
					this.y = source[1];
					this.width = source[2];
					this.height = source[3];
				}


			}

, extend: {

	toString: function() {
		return this.x + ' ' + this.y + ' ' + this.width + ' ' + this.height;
	}
  , morph: function(v){

	var v = arguments.length == 1 ?
        [v.x, v.y, v.width, v.height] :
        [].slice.call(arguments);

	this.destination = new SVG.ViewBox(v);

	return this;

}

  , at: function(pos) {

	if(!this.destination) return this;

	return new SVG.ViewBox([
		this.x + (this.destination.x - this.x) * pos
      , this.y + (this.destination.y - this.y) * pos
      , this.width + (this.destination.width - this.width) * pos
      , this.height + (this.destination.height - this.height) * pos
	]);

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
			return new SVG.ViewBox(this);

      // otherwise act as a setter
		v = arguments.length == 1 ?
        [v.x, v.y, v.width, v.height] :
        [].slice.call(arguments);

		return this.attr('viewBox', v);
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
		var self = this;

    // bind event to element rather than element node
		this.node['on' + event] = typeof f == 'function' ?
      function() { return f.apply(self, arguments); } : null;

		return this;
	};

});

// Initialize listeners stack
		SVG.listeners = [];
		SVG.handlerMap = [];
		SVG.listenerId = 0;

// Add event binder in the SVG namespace
		SVG.on = function(node, event, listener, binding) {
  // create listener, get object-index
			var l     = listener.bind(binding || node.instance || node)
    , index = (SVG.handlerMap.indexOf(node) + 1 || SVG.handlerMap.push(node)) - 1
    , ev    = event.split('.')[0]
    , ns    = event.split('.')[1] || '*';


  // ensure valid object
			SVG.listeners[index]         = SVG.listeners[index]         || {};
			SVG.listeners[index][ev]     = SVG.listeners[index][ev]     || {};
			SVG.listeners[index][ev][ns] = SVG.listeners[index][ev][ns] || {};

			if(!listener._svgjsListenerId)
				listener._svgjsListenerId = ++SVG.listenerId;

  // reference listener
			SVG.listeners[index][ev][ns][listener._svgjsListenerId] = l;

  // add listener
			node.addEventListener(ev, l, false);
		};

// Add event unbinder in the SVG namespace
		SVG.off = function(node, event, listener) {
			var index = SVG.handlerMap.indexOf(node)
    , ev    = event && event.split('.')[0]
    , ns    = event && event.split('.')[1];

			if(index == -1) return;

			if (listener) {
				if(typeof listener == 'function') listener = listener._svgjsListenerId;
				if(!listener) return;

    // remove listener reference
				if (SVG.listeners[index][ev] && SVG.listeners[index][ev][ns || '*']) {
      // remove listener
					node.removeEventListener(ev, SVG.listeners[index][ev][ns || '*'][listener], false);

					delete SVG.listeners[index][ev][ns || '*'][listener];
				}

			} else if (ns && ev) {
    // remove all listeners for a namespaced event
				if (SVG.listeners[index][ev] && SVG.listeners[index][ev][ns]) {
					for (listener in SVG.listeners[index][ev][ns])
						SVG.off(node, [ev, ns].join('.'), listener);

					delete SVG.listeners[index][ev][ns];
				}

			} else if (ns){
    // remove all listeners for a specific namespace
				for(event in SVG.listeners[index]){
					for(namespace in SVG.listeners[index][event]){
						if(ns === namespace){
							SVG.off(node, [event, ns].join('.'));
						}
					}
				}

			} else if (ev) {
    // remove all listeners for the event
				if (SVG.listeners[index][ev]) {
					for (namespace in SVG.listeners[index][ev])
						SVG.off(node, [ev, namespace].join('.'));

					delete SVG.listeners[index][ev];
				}

			} else {
    // remove all listeners on a given node
				for (event in SVG.listeners[index])
					SVG.off(node, event);

				delete SVG.listeners[index];

			}
		};

//
		SVG.extend(SVG.Element, {
  // Bind given event to listener
			on: function(event, listener, binding) {
				SVG.on(this.node, event, listener, binding);

				return this;
			}
  // Unbind event from listener
, off: function(event, listener) {
	SVG.off(this.node, event, listener);

	return this;
}
  // Fire given event
, fire: function(event, data) {

    // Dispatch event
	if(event instanceof Event){
		this.node.dispatchEvent(event);
	}else{
		this.node.dispatchEvent(new CustomEvent(event, {detail:data}));
	}

	return this;
}
		});

		SVG.Defs = SVG.invent({
  // Initialize node
			create: 'defs'

  // Inherit from
, inherit: SVG.Container

		});
		SVG.G = SVG.invent({
  // Initialize node
			create: 'g'

  // Inherit from
, inherit: SVG.Container

  // Add class methods
, extend: {
    // Move over x-axis
	x: function(x) {
		return x == null ? this.transform('x') : this.transform({ x: x - this.x() }, true);
	}
    // Move over y-axis
  , y: function(y) {
	return y == null ? this.transform('y') : this.transform({ y: y - this.y() }, true);
}
    // Move by center over x-axis
  , cx: function(x) {
	return x == null ? this.gbox().cx : this.x(x - this.gbox().width / 2);
}
    // Move by center over y-axis
  , cy: function(y) {
	return y == null ? this.gbox().cy : this.y(y - this.gbox().height / 2);
}
  , gbox: function() {

	var bbox  = this.bbox()
        , trans = this.transform();

	bbox.x  += trans.x;
	bbox.x2 += trans.x;
	bbox.cx += trans.x;

	bbox.y  += trans.y;
	bbox.y2 += trans.y;
	bbox.cy += trans.y;

	return bbox;
}
}

  // Add parent method
, construct: {
    // Create a group element
	group: function() {
		return this.put(new SVG.G);
	}
}
		});

// ### This module adds backward / forward functionality to elements.

//
		SVG.extend(SVG.Element, {
  // Get all siblings, including myself
			siblings: function() {
				return this.parent().children();
			}
  // Get the curent position siblings
, position: function() {
	return this.parent().index(this);
}
  // Get the next element (will return null if there is none)
, next: function() {
	return this.siblings()[this.position() + 1];
}
  // Get the next element (will return null if there is none)
, previous: function() {
	return this.siblings()[this.position() - 1];
}
  // Send given element one step forward
, forward: function() {
	var i = this.position() + 1
      , p = this.parent();

    // move node one step forward
	p.removeElement(this).add(this, i);

    // make sure defs node is always at the top
	if (p instanceof SVG.Doc)
		p.node.appendChild(p.defs().node);

	return this;
}
  // Send given element one step backward
, backward: function() {
	var i = this.position();

	if (i > 0)
		this.parent().removeElement(this).add(this, i - 1);

	return this;
}
  // Send given element all the way to the front
, front: function() {
	var p = this.parent();

    // Move node forward
	p.node.appendChild(this.node);

    // Make sure defs node is always at the top
	if (p instanceof SVG.Doc)
		p.node.appendChild(p.defs().node);

	return this;
}
  // Send given element all the way to the back
, back: function() {
	if (this.position() > 0)
		this.parent().removeElement(this).add(this, 0);

	return this;
}
  // Inserts a given element before the targeted element
, before: function(element) {
	element.remove();

	var i = this.position();

	this.parent().add(element, i);

	return this;
}
  // Insters a given element after the targeted element
, after: function(element) {
	element.remove();

	var i = this.position();

	this.parent().add(element, i + 1);

	return this;
}

		});
		SVG.Mask = SVG.invent({
  // Initialize node
			create: function() {
				this.constructor.call(this, SVG.create('mask'));

    // keep references to masked elements
				this.targets = [];
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
				this.targets[i].unmask();
		this.targets = [];

      // remove mask from parent
		this.parent().removeElement(this);

		return this;
	}
}

  // Add parent method
, construct: {
    // Create masking element
	mask: function() {
		return this.defs().put(new SVG.Mask);
	}
}
		});


		SVG.extend(SVG.Element, {
  // Distribute mask to svg element
			maskWith: function(element) {
    // use given mask or create a new one
				this.masker = element instanceof SVG.Mask ? element : this.parent().mask().add(element);

    // store reverence on self in mask
				this.masker.targets.push(this);

    // apply mask
				return this.attr('mask', 'url("#' + this.masker.attr('id') + '")');
			}
  // Unmask element
, unmask: function() {
	delete this.masker;
	return this.attr('mask', null);
}

		});

		SVG.ClipPath = SVG.invent({
  // Initialize node
			create: function() {
				this.constructor.call(this, SVG.create('clipPath'));

    // keep references to clipped elements
				this.targets = [];
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
				this.targets[i].unclip();
		this.targets = [];

      // remove clipPath from parent
		this.parent().removeElement(this);

		return this;
	}
}

  // Add parent method
, construct: {
    // Create clipping element
	clip: function() {
		return this.defs().put(new SVG.ClipPath);
	}
}
		});

//
		SVG.extend(SVG.Element, {
  // Distribute clipPath to svg element
			clipWith: function(element) {
    // use given clip or create a new one
				this.clipper = element instanceof SVG.ClipPath ? element : this.parent().clip().add(element);

    // store reverence on self in mask
				this.clipper.targets.push(this);

    // apply mask
				return this.attr('clip-path', 'url("#' + this.clipper.attr('id') + '")');
			}
  // Unclip element
, unclip: function() {
	delete this.clipper;
	return this.attr('clip-path', null);
}

		});
		SVG.Gradient = SVG.invent({
  // Initialize node
			create: function(type) {
				this.constructor.call(this, SVG.create(type + 'Gradient'));

    // store type
				this.type = type;
			}

  // Inherit from
, inherit: SVG.Container

  // Add class methods
, extend: {
    // Add a color stop
	at: function(offset, color, opacity) {
		return this.put(new SVG.Stop).update(offset, color, opacity);
	}
    // Update gradient
  , update: function(block) {
      // remove all stops
	this.clear();

      // invoke passed block
	if (typeof block == 'function')
		block.call(this, this);

	return this;
}
    // Return the fill id
  , fill: function() {
	return 'url(#' + this.id() + ')';
}
    // Alias string convertion to fill
  , toString: function() {
	return this.fill();
}
    // custom attr to handle transform
  , attr: function(a, b, c) {
	if(a == 'transform') a = 'gradientTransform';
	return SVG.Container.prototype.attr.call(this, a, b, c);
}
}

  // Add parent method
, construct: {
    // Create gradient element in defs
	gradient: function(type, block) {
		return this.defs().gradient(type, block);
	}
}
		});

// Add animatable methods to both gradient and fx module
		SVG.extend(SVG.Gradient, SVG.FX, {
  // From position
			from: function(x, y) {
				return (this._target || this).type == 'radial' ?
      this.attr({ fx: new SVG.Number(x), fy: new SVG.Number(y) }) :
      this.attr({ x1: new SVG.Number(x), y1: new SVG.Number(y) });
			}
  // To position
, to: function(x, y) {
	return (this._target || this).type == 'radial' ?
      this.attr({ cx: new SVG.Number(x), cy: new SVG.Number(y) }) :
      this.attr({ x2: new SVG.Number(x), y2: new SVG.Number(y) });
}
		});

// Base gradient generation
		SVG.extend(SVG.Defs, {
  // define gradient
			gradient: function(type, block) {
				return this.put(new SVG.Gradient(type)).update(block);
			}

		});

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
			};
		}

      // set attributes
		if (o.opacity != null) this.attr('stop-opacity', o.opacity);
		if (o.color   != null) this.attr('stop-color', o.color);
		if (o.offset  != null) this.attr('offset', new SVG.Number(o.offset));

		return this;
	}
}

		});

		SVG.Pattern = SVG.invent({
  // Initialize node
			create: 'pattern'

  // Inherit from
, inherit: SVG.Container

  // Add class methods
, extend: {
    // Return the fill id
	fill: function() {
		return 'url(#' + this.id() + ')';
	}
    // Update pattern by rebuilding
  , update: function(block) {
      // remove content
	this.clear();

      // invoke passed block
	if (typeof block == 'function')
		block.call(this, this);

	return this;
}
    // Alias string convertion to fill
  , toString: function() {
	return this.fill();
}
    // custom attr to handle transform
  , attr: function(a, b, c) {
	if(a == 'transform') a = 'patternTransform';
	return SVG.Container.prototype.attr.call(this, a, b, c);
}

}

  // Add parent method
, construct: {
    // Create pattern element in defs
	pattern: function(width, height, block) {
		return this.defs().pattern(width, height, block);
	}
}
		});

		SVG.extend(SVG.Defs, {
  // Define gradient
			pattern: function(width, height, block) {
				return this.put(new SVG.Pattern).update(block).attr({
					x:            0
    , y:            0
    , width:        width
    , height:       height
    , patternUnits: 'userSpaceOnUse'
				});
			}

		});
		SVG.Doc = SVG.invent({
  // Initialize node
			create: function(element) {
				if (element) {
      // ensure the presence of a dom element
					element = typeof element == 'string' ?
        document.getElementById(element) :
        element;

      // If the target is an svg element, use that element as the main wrapper.
      // This allows svg.js to work with svg documents as well.
					if (element.nodeName == 'svg') {
						this.constructor.call(this, element);
					} else {
						this.constructor.call(this, SVG.create('svg'));
						element.appendChild(this.node);
						this.size('100%', '100%');
					}

      // set svg element attributes and ensure defs node
					this.namespace().defs();
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
        .attr('xmlns:svgjs', SVG.svgjs, SVG.xmlns);
	}
    // Creates and returns defs element
  , defs: function() {
	if (!this._defs) {
		var defs;

        // Find or create a defs element in this instance
		if (defs = this.node.getElementsByTagName('defs')[0])
			this._defs = SVG.adopt(defs);
		else
          this._defs = new SVG.Defs;

        // Make sure the defs node is at the end of the stack
		this.node.appendChild(this._defs.node);
	}

	return this._defs;
}
    // custom parent method
  , parent: function() {
	return this.node.parentNode.nodeName == '#document' ? null : this.node.parentNode;
}
    // Fix for possible sub-pixel offset. See:
    // https://bugzilla.mozilla.org/show_bug.cgi?id=608812
  , spof: function(spof) {
	var pos = this.node.getScreenCTM();

	if (pos)
		this
          .style('left', (-pos.e % 1) + 'px')
          .style('top',  (-pos.f % 1) + 'px');

	return this;
}

      // Removes the doc from the DOM
  , remove: function() {
	if(this.parent()) {
		this.parent().removeChild(this.node);
	}

	return this;
}
}

		});

		SVG.Shape = SVG.invent({
  // Initialize node
			create: function(element) {
				this.constructor.call(this, element);
			}

  // Inherit from
, inherit: SVG.Element

		});

		SVG.Bare = SVG.invent({
  // Initialize
			create: function(element, inherit) {
    // construct element
				this.constructor.call(this, SVG.create(element));

    // inherit custom methods
				if (inherit)
					for (var method in inherit.prototype)
						if (typeof inherit.prototype[method] === 'function')
							this[method] = inherit.prototype[method];
			}

  // Inherit from
, inherit: SVG.Element

  // Add methods
, extend: {
    // Insert some plain text
	words: function(text) {
      // remove contents
		while (this.node.hasChildNodes())
			this.node.removeChild(this.node.lastChild);

      // create text node
		this.node.appendChild(document.createTextNode(text));

		return this;
	}
}
		});


		SVG.extend(SVG.Parent, {
  // Create an element that is not described by SVG.js
			element: function(element, inherit) {
				return this.put(new SVG.Bare(element, inherit));
			}
  // Add symbol element
, symbol: function() {
	return this.defs().element('symbol', SVG.Container);
}

		});
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
		return this.attr('href', (file || '') + '#' + element, SVG.xlink);
	}
}

  // Add parent method
, construct: {
    // Create a use element
	use: function(element, file) {
		return this.put(new SVG.Use).element(element, file);
	}
}
		});
		SVG.Rect = SVG.invent({
  // Initialize node
			create: 'rect'

  // Inherit from
, inherit: SVG.Shape

  // Add parent method
, construct: {
    // Create a rect element
	rect: function(width, height) {
		return this.put(new SVG.Rect()).size(width, height);
	}
}
		});
		SVG.Circle = SVG.invent({
  // Initialize node
			create: 'circle'

  // Inherit from
, inherit: SVG.Shape

  // Add parent method
, construct: {
    // Create circle element, based on ellipse
	circle: function(size) {
		return this.put(new SVG.Circle).rx(new SVG.Number(size).divide(2)).move(0, 0);
	}
}
		});

		SVG.extend(SVG.Circle, SVG.FX, {
  // Radius x value
			rx: function(rx) {
				return this.attr('r', rx);
			}
  // Alias radius x value
, ry: function(ry) {
	return this.rx(ry);
}
		});

		SVG.Ellipse = SVG.invent({
  // Initialize node
			create: 'ellipse'

  // Inherit from
, inherit: SVG.Shape

  // Add parent method
, construct: {
    // Create an ellipse
	ellipse: function(width, height) {
		return this.put(new SVG.Ellipse).size(width, height).move(0, 0);
	}
}
		});

		SVG.extend(SVG.Ellipse, SVG.Rect, SVG.FX, {
  // Radius x value
			rx: function(rx) {
				return this.attr('rx', rx);
			}
  // Radius y value
, ry: function(ry) {
	return this.attr('ry', ry);
}
		});

// Add common method
		SVG.extend(SVG.Circle, SVG.Ellipse, {
    // Move over x-axis
			x: function(x) {
				return x == null ? this.cx() - this.rx() : this.cx(x + this.rx());
			}
    // Move over y-axis
  , y: function(y) {
	return y == null ? this.cy() - this.ry() : this.cy(y + this.ry());
}
    // Move by center over x-axis
  , cx: function(x) {
	return x == null ? this.attr('cx') : this.attr('cx', x);
}
    // Move by center over y-axis
  , cy: function(y) {
	return y == null ? this.attr('cy') : this.attr('cy', y);
}
    // Set width of element
  , width: function(width) {
	return width == null ? this.rx() * 2 : this.rx(new SVG.Number(width).divide(2));
}
    // Set height of element
  , height: function(height) {
	return height == null ? this.ry() * 2 : this.ry(new SVG.Number(height).divide(2));
}
    // Custom size function
  , size: function(width, height) {
	var p = proportionalSize(this, width, height);

	return this
        .rx(new SVG.Number(p.width).divide(2))
        .ry(new SVG.Number(p.height).divide(2));
}
		});
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
		]);
	}
    // Overwrite native plot() method
  , plot: function(x1, y1, x2, y2) {
	if (typeof y1 !== 'undefined')
		x1 = { x1: x1, y1: y1, x2: x2, y2: y2 };
	else
        x1 = new SVG.PointArray(x1).toLine();

	return this.attr(x1);
}
    // Move by left top corner
  , move: function(x, y) {
	return this.attr(this.array().move(x, y).toLine());
}
    // Set element size to given width and height
  , size: function(width, height) {
	var p = proportionalSize(this, width, height);

	return this.attr(this.array().size(p.width, p.height).toLine());
}
}

  // Add parent method
, construct: {
    // Create a line element
	line: function(x1, y1, x2, y2) {
		return this.put(new SVG.Line).plot(x1, y1, x2, y2);
	}
}
		});

		SVG.Polyline = SVG.invent({
  // Initialize node
			create: 'polyline'

  // Inherit from
, inherit: SVG.Shape

  // Add parent method
, construct: {
    // Create a wrapped polyline element
	polyline: function(p) {
		return this.put(new SVG.Polyline).plot(p);
	}
}
		});

		SVG.Polygon = SVG.invent({
  // Initialize node
			create: 'polygon'

  // Inherit from
, inherit: SVG.Shape

  // Add parent method
, construct: {
    // Create a wrapped polygon element
	polygon: function(p) {
		return this.put(new SVG.Polygon).plot(p);
	}
}
		});

// Add polygon-specific functions
		SVG.extend(SVG.Polyline, SVG.Polygon, {
  // Get array
			array: function() {
				return this._array || (this._array = new SVG.PointArray(this.attr('points')));
			}
  // Plot new path
, plot: function(p) {
	return this.attr('points', (this._array = new SVG.PointArray(p)));
}
  // Move by left top corner
, move: function(x, y) {
	return this.attr('points', this.array().move(x, y));
}
  // Set element size to given width and height
, size: function(width, height) {
	var p = proportionalSize(this, width, height);

	return this.attr('points', this.array().size(p.width, p.height));
}

		});
// unify all point to point elements
		SVG.extend(SVG.Line, SVG.Polyline, SVG.Polygon, {
  // Define morphable array
			morphArray:  SVG.PointArray
  // Move by left top corner over x-axis
, x: function(x) {
	return x == null ? this.bbox().x : this.move(x, this.bbox().y);
}
  // Move by left top corner over y-axis
, y: function(y) {
	return y == null ? this.bbox().y : this.move(this.bbox().x, y);
}
  // Set width of element
, width: function(width) {
	var b = this.bbox();

	return width == null ? b.width : this.size(width, b.height);
}
  // Set height of element
, height: function(height) {
	var b = this.bbox();

	return height == null ? b.height : this.size(b.width, height);
}
		});
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
	return this._array || (this._array = new SVG.PathArray(this.attr('d')));
}
    // Plot new poly points
  , plot: function(p) {
	return this.attr('d', (this._array = new SVG.PathArray(p)));
}
    // Move by left top corner
  , move: function(x, y) {
	return this.attr('d', this.array().move(x, y));
}
    // Move by left top corner over x-axis
  , x: function(x) {
	return x == null ? this.bbox().x : this.move(x, this.bbox().y);
}
    // Move by left top corner over y-axis
  , y: function(y) {
	return y == null ? this.bbox().y : this.move(this.bbox().x, y);
}
    // Set element size to given width and height
  , size: function(width, height) {
	var p = proportionalSize(this, width, height);

	return this.attr('d', this.array().size(p.width, p.height));
}
    // Set width of element
  , width: function(width) {
	return width == null ? this.bbox().width : this.size(width, this.bbox().height);
}
    // Set height of element
  , height: function(height) {
	return height == null ? this.bbox().height : this.size(this.bbox().width, height);
}

}

  // Add parent method
, construct: {
    // Create a wrapped path element
	path: function(d) {
		return this.put(new SVG.Path).plot(d);
	}
}
		});
		SVG.Image = SVG.invent({
  // Initialize node
			create: 'image'

  // Inherit from
, inherit: SVG.Shape

  // Add class methods
, extend: {
    // (re)load image
	load: function(url) {
		if (!url) return this;

		var self = this
        , img  = document.createElement('img');

      // preload image
		img.onload = function() {
			var p = self.parent(SVG.Pattern);

			if(p === null) return;

        // ensure image size
			if (self.width() == 0 && self.height() == 0)
				self.size(img.width, img.height);

        // ensure pattern size if not set
			if (p && p.width() == 0 && p.height() == 0)
				p.size(self.width(), self.height());

        // callback
			if (typeof self._loaded === 'function')
				self._loaded.call(self, {
					width:  img.width
          , height: img.height
          , ratio:  img.width / img.height
          , url:    url
				});
		};

		img.onerror = function(e){
			if (typeof self._error === 'function'){
				self._error.call(self, e);
			}
		};

		return this.attr('href', (img.src = this.src = url), SVG.xlink);
	}
    // Add loaded callback
  , loaded: function(loaded) {
	this._loaded = loaded;
	return this;
}

  , error: function(error) {
	this._error = error;
	return this;
}
}

  // Add parent method
, construct: {
    // create image element, load image and set its size
	image: function(source, width, height) {
		return this.put(new SVG.Image).load(source).size(width || 0, height || width || 0);
	}
}

		});
		SVG.Text = SVG.invent({
  // Initialize node
			create: function() {
				this.constructor.call(this, SVG.create('text'));

				this.dom.leading = new SVG.Number(1.3);    // store leading value for rebuilding
				this._rebuild = true;                      // enable automatic updating of dy values
				this._build   = false;                     // disable build mode for adding multiple lines

    // set default font
				this.attr('font-family', SVG.defaults.attrs['font-family']);
			}

  // Inherit from
, inherit: SVG.Shape

  // Add class methods
, extend: {
    // Move over x-axis
	x: function(x) {
      // act as getter
		if (x == null)
			return this.attr('x');

      // move lines as well if no textPath is present
		if (!this.textPath)
			this.lines().each(function() { if (this.dom.newLined) this.x(x); });

		return this.attr('x', x);
	}
    // Move over y-axis
  , y: function(y) {
	var oy = this.attr('y')
        , o  = typeof oy === 'number' ? oy - this.bbox().y : 0;

      // act as getter
	if (y == null)
		return typeof oy === 'number' ? oy - o : oy;

	return this.attr('y', typeof y === 'number' ? y + o : y);
}
    // Move center over x-axis
  , cx: function(x) {
	return x == null ? this.bbox().cx : this.x(x - this.bbox().width / 2);
}
    // Move center over y-axis
  , cy: function(y) {
	return y == null ? this.bbox().cy : this.y(y - this.bbox().height / 2);
}
    // Set the text content
  , text: function(text) {
      // act as getter
	if (typeof text === 'undefined'){
		var text = '';
		var children = this.node.childNodes;
		for(var i = 0, len = children.length; i < len; ++i){

          // add newline if its not the first child and newLined is set to true
			if(i != 0 && children[i].nodeType != 3 && SVG.adopt(children[i]).dom.newLined == true){
				text += '\n';
			}

          // add content of this node
			text += children[i].textContent;
		}

		return text;
	}

      // remove existing content
	this.clear().build(true);

	if (typeof text === 'function') {
        // call block
		text.call(this, this);

	} else {
        // store text and make sure text is not blank
		text = text.split('\n');

        // build new lines
		for (var i = 0, il = text.length; i < il; i++)
			this.tspan(text[i]).newLine();
	}

      // disable build mode and rebuild lines
	return this.build(false).rebuild();
}
    // Set font size
  , size: function(size) {
	return this.attr('font-size', size).rebuild();
}
    // Set / get leading
  , leading: function(value) {
      // act as getter
	if (value == null)
		return this.dom.leading;

      // act as setter
	this.dom.leading = new SVG.Number(value);

	return this.rebuild();
}
    // Get all the first level lines
  , lines: function() {
	var node = (this.textPath && this.textPath() || this).node;

      // filter tspans and map them to SVG.js instances
	var lines = SVG.utils.map(SVG.utils.filterSVGElements(node.childNodes), function(el){
		return SVG.adopt(el);
	});

      // return an instance of SVG.set
	return new SVG.Set(lines);
}
    // Rebuild appearance type
  , rebuild: function(rebuild) {
      // store new rebuild flag if given
	if (typeof rebuild == 'boolean')
		this._rebuild = rebuild;

      // define position of all lines
	if (this._rebuild) {
		var self = this
          , blankLineOffset = 0
          , dy = this.dom.leading * new SVG.Number(this.attr('font-size'));

		this.lines().each(function() {
			if (this.dom.newLined) {
				if (!this.textPath)
					this.attr('x', self.attr('x'));

				if(this.text() == '\n') {
					blankLineOffset += dy;
				}else{
					this.attr('dy', dy + blankLineOffset);
					blankLineOffset = 0;
				}
			}
		});

		this.fire('rebuild');
	}

	return this;
}
    // Enable / disable build mode
  , build: function(build) {
	this._build = !!build;
	return this;
}
    // overwrite method from parent to set data properly
  , setData: function(o){
	this.dom = o;
	this.dom.leading = new SVG.Number(o.leading || 1.3);
	return this;
}
}

  // Add parent method
, construct: {
    // Create text element
	text: function(text) {
		return this.put(new SVG.Text).text(text);
	}
    // Create plain text element
  , plain: function(text) {
	return this.put(new SVG.Text).plain(text);
}
}

		});

		SVG.Tspan = SVG.invent({
  // Initialize node
			create: 'tspan'

  // Inherit from
, inherit: SVG.Shape

  // Add class methods
, extend: {
    // Set text content
	text: function(text) {
		if(text == null) return this.node.textContent + (this.dom.newLined ? '\n' : '');

		typeof text === 'function' ? text.call(this, this) : this.plain(text);

		return this;
	}
    // Shortcut dx
  , dx: function(dx) {
	return this.attr('dx', dx);
}
    // Shortcut dy
  , dy: function(dy) {
	return this.attr('dy', dy);
}
    // Create new line
  , newLine: function() {
      // fetch text parent
	var t = this.parent(SVG.Text);

      // mark new line
	this.dom.newLined = true;

      // apply new hy¡n
	return this.dy(t.dom.leading * t.attr('font-size')).attr('x', t.x());
}
}

		});

		SVG.extend(SVG.Text, SVG.Tspan, {
  // Create plain text node
			plain: function(text) {
    // clear if build mode is disabled
				if (this._build === false)
					this.clear();

    // create text node
				this.node.appendChild(document.createTextNode(text));

				return this;
			}
  // Create a tspan
, tspan: function(text) {
	var node  = (this.textPath && this.textPath() || this).node
      , tspan = new SVG.Tspan;

    // clear if build mode is disabled
	if (this._build === false)
		this.clear();

    // add new tspan
	node.appendChild(tspan.node);

	return tspan.text(text);
}
  // Clear all lines
, clear: function() {
	var node = (this.textPath && this.textPath() || this).node;

    // remove existing child nodes
	while (node.hasChildNodes())
		node.removeChild(node.lastChild);

	return this;
}
  // Get length of text element
, length: function() {
	return this.node.getComputedTextLength();
}
		});

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
        , track = this.doc().defs().path(d);

      // move lines to textpath
		while (this.node.hasChildNodes())
			path.node.appendChild(this.node.firstChild);

      // add textPath element as child node
		this.node.appendChild(path.node);

      // link textPath to path and add content
		path.attr('href', '#' + track, SVG.xlink);

		return this;
	}
    // Plot path if any
  , plot: function(d) {
	var track = this.track();

	if (track)
		track.plot(d);

	return this;
}
    // Get the path track element
  , track: function() {
	var path = this.textPath();

	if (path)
		return path.reference('href');
}
    // Get the textPath child
  , textPath: function() {
	if (this.node.firstChild && this.node.firstChild.nodeName == 'textPath')
		return SVG.adopt(this.node.firstChild);
}
}
		});
		SVG.Nested = SVG.invent({
  // Initialize node
			create: function() {
				this.constructor.call(this, SVG.create('svg'));

				this.style('overflow', 'visible');
			}

  // Inherit from
, inherit: SVG.Container

  // Add parent method
, construct: {
    // Create nested svg document
	nested: function() {
		return this.put(new SVG.Nested);
	}
}
		});
		SVG.A = SVG.invent({
  // Initialize node
			create: 'a'

  // Inherit from
, inherit: SVG.Container

  // Add class methods
, extend: {
    // Link url
	to: function(url) {
		return this.attr('href', url, SVG.xlink);
	}
    // Link show attribute
  , show: function(target) {
	return this.attr('show', target, SVG.xlink);
}
    // Link target attribute
  , target: function(target) {
	return this.attr('target', target);
}
}

  // Add parent method
, construct: {
    // Create a hyperlink element
	link: function(url) {
		return this.put(new SVG.A).to(url);
	}
}
		});

		SVG.extend(SVG.Element, {
  // Create a hyperlink element
			linkTo: function(url) {
				var link = new SVG.A;

				if (typeof url == 'function')
					url.call(link, link);
				else
      link.to(url);

				return this.parent().put(link).put(this);
			}

		});
		SVG.Marker = SVG.invent({
  // Initialize node
			create: 'marker'

  // Inherit from
, inherit: SVG.Container

  // Add class methods
, extend: {
    // Set width of element
	width: function(width) {
		return this.attr('markerWidth', width);
	}
    // Set height of element
  , height: function(height) {
	return this.attr('markerHeight', height);
}
    // Set marker refX and refY
  , ref: function(x, y) {
	return this.attr('refX', x).attr('refY', y);
}
    // Update marker
  , update: function(block) {
      // remove all content
	this.clear();

      // invoke passed block
	if (typeof block == 'function')
		block.call(this, this);

	return this;
}
    // Return the fill id
  , toString: function() {
	return 'url(#' + this.id() + ')';
}
}

  // Add parent method
, construct: {
	marker: function(width, height, block) {
      // Create marker element in defs
		return this.defs().marker(width, height, block);
	}
}

		});

		SVG.extend(SVG.Defs, {
  // Create marker
			marker: function(width, height, block) {
    // Set default viewbox to match the width and height, set ref to cx and cy and set orient to auto
				return this.put(new SVG.Marker)
      .size(width, height)
      .ref(width / 2, height / 2)
      .viewbox(0, 0, width, height)
      .attr('orient', 'auto')
      .update(block);
			}

		});

		SVG.extend(SVG.Line, SVG.Polyline, SVG.Polygon, SVG.Path, {
  // Create and attach markers
			marker: function(marker, width, height, block) {
				var attr = ['marker'];

    // Build attribute name
				if (marker != 'all') attr.push(marker);
				attr = attr.join('-');

    // Set marker attribute
				marker = arguments[1] instanceof SVG.Marker ?
      arguments[1] :
      this.doc().marker(width, height, block);

				return this.attr(attr, marker);
			}

		});
// Define list of available attributes for stroke and fill
		var sugar = {
			stroke: ['color', 'width', 'opacity', 'linecap', 'linejoin', 'miterlimit', 'dasharray', 'dashoffset']
, fill:   ['color', 'opacity', 'rule']
, prefix: function(t, a) {
	return a == 'color' ? t : t + '-' + a;
}
		}

// Add sugar for fill and stroke
;['fill', 'stroke'].forEach(function(m) {
	var i, extension = {};

	extension[m] = function(o) {
		if (typeof o == 'undefined')
			return this;
		if (typeof o == 'string' || SVG.Color.isRgb(o) || (o && typeof o.fill === 'function'))
			this.attr(m, o);

		else
      // set all attributes from sugar.fill and sugar.stroke list
      for (i = sugar[m].length - 1; i >= 0; i--)
	if (o[sugar[m][i]] != null)
		this.attr(sugar.prefix(m, sugar[m][i]), o[sugar[m][i]]);

		return this;
	};

	SVG.extend(SVG.Element, SVG.FX, extension);

});

		SVG.extend(SVG.Element, SVG.FX, {
  // Map rotation to transform
			rotate: function(d, cx, cy) {
				return this.transform({ rotation: d, cx: cx, cy: cy });
			}
  // Map skew to transform
, skew: function(x, y, cx, cy) {
	return arguments.length == 1  || arguments.length == 3 ?
      this.transform({ skew: x, cx: y, cy: cx }) :
      this.transform({ skewX: x, skewY: y, cx: cx, cy: cy });
}
  // Map scale to transform
, scale: function(x, y, cx, cy) {
	return arguments.length == 1  || arguments.length == 3 ?
      this.transform({ scale: x, cx: y, cy: cx }) :
      this.transform({ scaleX: x, scaleY: y, cx: cx, cy: cy });
}
  // Map translate to transform
, translate: function(x, y) {
	return this.transform({ x: x, y: y });
}
  // Map flip to transform
, flip: function(a, o) {
	return this.transform({ flip: a, offset: o });
}
  // Map matrix to transform
, matrix: function(m) {
	return this.attr('transform', new SVG.Matrix(m));
}
  // Opacity
, opacity: function(value) {
	return this.attr('opacity', value);
}
  // Relative move over x axis
, dx: function(x) {
	return this.x((this instanceof SVG.FX ? 0 : this.x()) + x, true);
}
  // Relative move over y axis
, dy: function(y) {
	return this.y((this instanceof SVG.FX ? 0 : this.y()) + y, true);
}
  // Relative move over x and y axes
, dmove: function(x, y) {
	return this.dx(x).dy(y);
}
		});

		SVG.extend(SVG.Rect, SVG.Ellipse, SVG.Circle, SVG.Gradient, SVG.FX, {
  // Add x and y radius
			radius: function(x, y) {
				var type = (this._target || this).type;
				return type == 'radial' || type == 'circle' ?
      this.attr('r', new SVG.Number(x)) :
      this.rx(x).ry(y == null ? x : y);
			}
		});

		SVG.extend(SVG.Path, {
  // Get path length
			length: function() {
				return this.node.getTotalLength();
			}
  // Get point at length
, pointAt: function(length) {
	return this.node.getPointAtLength(length);
}
		});

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
        this.attr(k, o[k]);

				return this;
			}
		});

		SVG.Set = SVG.invent({
  // Initialize
			create: function(members) {
    // Set initial state
				Array.isArray(members) ? this.members = members : this.clear();
			}

  // Add class methods
, extend: {
    // Add element to set
	add: function() {
		var i, il, elements = [].slice.call(arguments);

		for (i = 0, il = elements.length; i < il; i++)
			this.members.push(elements[i]);

		return this;
	}
    // Remove element from set
  , remove: function(element) {
	var i = this.index(element);

      // remove given child
	if (i > -1)
		this.members.splice(i, 1);

	return this;
}
    // Iterate over all members
  , each: function(block) {
	for (var i = 0, il = this.members.length; i < il; i++)
		block.apply(this.members[i], [i, this.members]);

	return this;
}
    // Restore to defaults
  , clear: function() {
      // initialize store
	this.members = [];

	return this;
}
    // Get the length of a set
  , length: function() {
	return this.members.length;
}
    // Checks if a given element is present in set
  , has: function(element) {
	return this.index(element) >= 0;
}
    // retuns index of given element in set
  , index: function(element) {
	return this.members.indexOf(element);
}
    // Get member at given index
  , get: function(i) {
	return this.members[i];
}
    // Get first member
  , first: function() {
	return this.get(0);
}
    // Get last member
  , last: function() {
	return this.get(this.members.length - 1);
}
    // Default value
  , valueOf: function() {
	return this.members;
}
    // Get the bounding box of all members included or empty box if set has no items
  , bbox: function(){
	var box = new SVG.BBox();

      // return an empty box of there are no members
	if (this.members.length == 0)
		return box;

      // get the first rbox and update the target bbox
	var rbox = this.members[0].rbox();
	box.x      = rbox.x;
	box.y      = rbox.y;
	box.width  = rbox.width;
	box.height = rbox.height;

	this.each(function() {
        // user rbox for correct position and visual representation
		box = box.merge(this.rbox());
	});

	return box;
}
}

  // Add parent method
, construct: {
    // Create a new set
	set: function(members) {
		return new SVG.Set(members);
	}
}
		});

		SVG.FX.Set = SVG.invent({
  // Initialize node
			create: function(set) {
    // store reference to set
				this.set = set;
			}

		});

// Alias methods
		SVG.Set.inherit = function() {
			var m
    , methods = [];

  // gather shape methods
			for(var m in SVG.Shape.prototype)
				if (typeof SVG.Shape.prototype[m] == 'function' && typeof SVG.Set.prototype[m] != 'function')
					methods.push(m);

  // apply shape aliasses
			methods.forEach(function(method) {
				SVG.Set.prototype[method] = function() {
					for (var i = 0, il = this.members.length; i < il; i++)
						if (this.members[i] && typeof this.members[i][method] == 'function')
							this.members[i][method].apply(this.members[i], arguments);

					return method == 'animate' ? (this.fx || (this.fx = new SVG.FX.Set(this))) : this;
				};
			});

  // clear methods for the next round
			methods = [];

  // gather fx methods
			for(var m in SVG.FX.prototype)
				if (typeof SVG.FX.prototype[m] == 'function' && typeof SVG.FX.Set.prototype[m] != 'function')
					methods.push(m);

  // apply fx aliasses
			methods.forEach(function(method) {
				SVG.FX.Set.prototype[method] = function() {
					for (var i = 0, il = this.set.members.length; i < il; i++)
						this.set.members[i].fx[method].apply(this.set.members[i].fx, arguments);

					return this;
				};
			});
		};




		SVG.extend(SVG.Element, {
  // Store data values on svg nodes
			data: function(a, v, r) {
				if (typeof a == 'object') {
					for (v in a)
						this.data(v, a[v]);

				} else if (arguments.length < 2) {
					try {
						return JSON.parse(this.attr('data-' + a));
					} catch(e) {
						return this.attr('data-' + a);
					}

				} else {
					this.attr(
        'data-' + a
      , v === null ?
          null :
        r === true || typeof v === 'string' || typeof v === 'number' ?
          v :
          JSON.stringify(v)
      );
				}

				return this;
			}
		});
		SVG.extend(SVG.Element, {
  // Remember arbitrary data
			remember: function(k, v) {
    // remember every item in an object individually
				if (typeof arguments[0] == 'object')
					for (var v in k)
						this.remember(v, k[v]);

    // retrieve memory
				else if (arguments.length == 1)
					return this.memory()[k];

    // store memory
				else
      this.memory()[k] = v;

				return this;
			}

  // Erase a given memory
, forget: function() {
	if (arguments.length == 0)
		this._memory = {};
	else
      for (var i = arguments.length - 1; i >= 0; i--)
	delete this.memory()[arguments[i]];

	return this;
}

  // Initialize or return local memory object
, memory: function() {
	return this._memory || (this._memory = {});
}

		});
// Method for getting an element by id
		SVG.get = function(id) {
			var node = document.getElementById(idFromReference(id) || id);
			return SVG.adopt(node);
		};

// Select elements by query string
		SVG.select = function(query, parent) {
			return new SVG.Set(
    SVG.utils.map((parent || document).querySelectorAll(query), function(node) {
	return SVG.adopt(node);
})
  );
		};

		SVG.extend(SVG.Parent, {
  // Scoped select method
			select: function(query) {
				return SVG.select(query, this.node);
			}

		});
		function is(el, obj){
			return el instanceof obj;
		}

// tests if a given selector matches an element
		function matches(el, selector) {
			return (el.matches || el.matchesSelector || el.msMatchesSelector || el.mozMatchesSelector || el.webkitMatchesSelector || el.oMatchesSelector).call(el, selector);
		}

// Convert dash-separated-string to camelCase
		function camelCase(s) {
			return s.toLowerCase().replace(/-(.)/g, function(m, g) {
				return g.toUpperCase();
			});
		}

// Capitalize first letter of a string
		function capitalize(s) {
			return s.charAt(0).toUpperCase() + s.slice(1);
		}

// Ensure to six-based hex
		function fullHex(hex) {
			return hex.length == 4 ?
    [ '#',
      hex.substring(1, 2), hex.substring(1, 2)
    , hex.substring(2, 3), hex.substring(2, 3)
    , hex.substring(3, 4), hex.substring(3, 4)
    ].join('') : hex;
		}

// Component to hex value
		function compToHex(comp) {
			var hex = comp.toString(16);
			return hex.length == 1 ? '0' + hex : hex;
		}

// Calculate proportional width and height values when necessary
		function proportionalSize(element, width, height) {
			if (width == null || height == null) {
				var box = element.bbox();

				if (width == null)
					width = box.width / box.height * height;
				else if (height == null)
					height = box.height / box.width * width;
			}

			return {
				width:  width
  , height: height
			};
		}

// Delta transform point
		function deltaTransformPoint(matrix, x, y) {
			return {
				x: x * matrix.a + y * matrix.c + 0
  , y: x * matrix.b + y * matrix.d + 0
			};
		}

// Map matrix array to object
		function arrayToMatrix(a) {
			return { a: a[0], b: a[1], c: a[2], d: a[3], e: a[4], f: a[5] };
		}

// Parse matrix if required
		function parseMatrix(matrix) {
			if (!(matrix instanceof SVG.Matrix))
				matrix = new SVG.Matrix(matrix);

			return matrix;
		}

// Add centre point to transform object
		function ensureCentre(o, target) {
			o.cx = o.cx == null ? target.bbox().cx : o.cx;
			o.cy = o.cy == null ? target.bbox().cy : o.cy;
		}

// Convert string to matrix
		function stringToMatrix(source) {
  // remove matrix wrapper and split to individual numbers
			source = source
    .replace(SVG.regex.whitespace, '')
    .replace(SVG.regex.matrix, '')
    .split(SVG.regex.matrixElements);

  // convert string values to floats and convert to a matrix-formatted object
			return arrayToMatrix(
    SVG.utils.map(source, function(n) {
	return parseFloat(n);
})
  );
		}

// Calculate position according to from and to
		function at(o, pos) {
  // number recalculation (don't bother converting to SVG.Number for performance reasons)
			return typeof o.from == 'number' ?
    o.from + (o.to - o.from) * pos :

  // instance recalculation
  o instanceof SVG.Color || o instanceof SVG.Number || o instanceof SVG.Matrix ? o.at(pos) :

  // for all other values wait until pos has reached 1 to return the final value
  pos < 1 ? o.from : o.to;
		}

// PathArray Helpers
		function arrayToString(a) {
			for (var i = 0, il = a.length, s = ''; i < il; i++) {
				s += a[i][0];

				if (a[i][1] != null) {
					s += a[i][1];

					if (a[i][2] != null) {
						s += ' ';
						s += a[i][2];

						if (a[i][3] != null) {
							s += ' ';
							s += a[i][3];
							s += ' ';
							s += a[i][4];

							if (a[i][5] != null) {
								s += ' ';
								s += a[i][5];
								s += ' ';
								s += a[i][6];

								if (a[i][7] != null) {
									s += ' ';
									s += a[i][7];
								}
							}
						}
					}
				}
			}

			return s + ' ';
		}

// Deep new id assignment
		function assignNewId(node) {
  // do the same for SVG child nodes as well
			for (var i = node.childNodes.length - 1; i >= 0; i--)
				if (node.childNodes[i] instanceof SVGElement)
					assignNewId(node.childNodes[i]);

			return SVG.adopt(node).id(SVG.eid(node.nodeName));
		}

// Add more bounding box properties
		function fullBox(b) {
			if (b.x == null) {
				b.x      = 0;
				b.y      = 0;
				b.width  = 0;
				b.height = 0;
			}

			b.w  = b.width;
			b.h  = b.height;
			b.x2 = b.x + b.width;
			b.y2 = b.y + b.height;
			b.cx = b.x + b.width / 2;
			b.cy = b.y + b.height / 2;

			return b;
		}

// Get id from reference string
		function idFromReference(url) {
			var m = url.toString().match(SVG.regex.reference);

			if (m) return m[1];
		}

// Create matrix array for looping
		var abcdef = 'abcdef'.split('');
// Add CustomEvent to IE9 and IE10
		if (typeof CustomEvent !== 'function') {
  // Code from: https://developer.mozilla.org/en-US/docs/Web/API/CustomEvent
			var CustomEvent = function(event, options) {
				options = options || { bubbles: false, cancelable: false, detail: undefined };
				var e = document.createEvent('CustomEvent');
				e.initCustomEvent(event, options.bubbles, options.cancelable, options.detail);
				return e;
			};

			CustomEvent.prototype = window.Event.prototype;

			window.CustomEvent = CustomEvent;
		}

// requestAnimationFrame / cancelAnimationFrame Polyfill with fallback based on Paul Irish
		(function(w) {
			var lastTime = 0;
			var vendors = ['moz', 'webkit'];

			for(var x = 0; x < vendors.length && !window.requestAnimationFrame; ++x) {
				w.requestAnimationFrame = w[vendors[x] + 'RequestAnimationFrame'];
				w.cancelAnimationFrame  = w[vendors[x] + 'CancelAnimationFrame'] ||
                              w[vendors[x] + 'CancelRequestAnimationFrame'];
			}

			w.requestAnimationFrame = w.requestAnimationFrame ||
    function(callback) {
	var currTime = new Date().getTime();
	var timeToCall = Math.max(0, 16 - (currTime - lastTime));

	var id = w.setTimeout(function() {
		callback(currTime + timeToCall);
	}, timeToCall);

	lastTime = currTime + timeToCall;
	return id;
};

			w.cancelAnimationFrame = w.cancelAnimationFrame || w.clearTimeout;

		}(window));

		return SVG;

	}));
},{}],3:[function(require,module,exports){
	'use strict';

	Object.defineProperty(exports, '__esModule', {
		value: true
	});

	var _svg = require('svg.js');

	var _svg2 = _interopRequireDefault(_svg);

	var _fitCurve = require('fit-curve');

	var _fitCurve2 = _interopRequireDefault(_fitCurve);

	function _interopRequireDefault(obj) { return obj && obj.__esModule ? obj : { default: obj }; }

	function GlyphsControl(pannel) {

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
			var pathString = fittedCurveDataToPathString(smoothBizer);
			pannel.path(pathString).fill('none').stroke({ width: 3 }).stroke('#f06');
			rawPointData.length = 0;
		};
	}

	exports.default = GlyphsControl;

},{'fit-curve':1,'svg.js':2}],4:[function(require,module,exports){
	'use strict';

	Object.defineProperty(exports, '__esModule', {
		value: true
	});

	var _svg = require('svg.js');

	var _svg2 = _interopRequireDefault(_svg);

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
			var pathString = fittedCurveDataToPathString(smoothBizer);
			pannel.path(pathString).fill('none').stroke({ width: 3 }).stroke('#f06');
			rawPointData.length = 0;
		};
	}

	exports.default = PaintControl;


	function updateLines(paintingPolyLine, rawPointData) {
		paintingPolyLine.plot(rawPointData);
	}
	function fittedCurveDataToPathString(fittedLineData) {
		var str = '';
		fittedLineData.map(function (bezier, i) {
			if (i == 0) {
				str += 'M ' + bezier[0][0] + ' ' + bezier[0][1];
			}
			str += 'C ' + bezier[1][0] + ' ' + bezier[1][1] + ', ' + bezier[2][0] + ' ' + bezier[2][1] + ', ' + bezier[3][0] + ' ' + bezier[3][1] + ' ';
		});

		return str;
	}

},{'fit-curve':1,'svg.js':2}],5:[function(require,module,exports){
	'use strict';

	var _GlyphsControl = require('./Controls/GlyphsControl');

	var _GlyphsControl2 = _interopRequireDefault(_GlyphsControl);

	var _PaintControl = require('./Controls/PaintControl');

	var _PaintControl2 = _interopRequireDefault(_PaintControl);

	var _svg = require('svg.js');

	var _svg2 = _interopRequireDefault(_svg);

	function _interopRequireDefault(obj) { return obj && obj.__esModule ? obj : { default: obj }; }

	var draw = (0, _svg2.default)('drawing').size(300, 300);

	setControl(draw);

	function setControl(_container) {
		var isMouseDown = false;
		var polyline = void 0;

		var currnetControl = new _PaintControl2.default(draw);

		_container.on('mousedown', function (e) {
			var point = [e.clientX, e.clientY];
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

},{'./Controls/GlyphsControl':3,'./Controls/PaintControl':4,'svg.js':2}]},{},[5]);
//# sourceMappingURL=data:application/json;charset=utf-8;base64,eyJ2ZXJzaW9uIjozLCJzb3VyY2VzIjpbIm5vZGVfbW9kdWxlcy9icm93c2VyLXBhY2svX3ByZWx1ZGUuanMiLCJub2RlX21vZHVsZXMvZml0LWN1cnZlL3NyYy9maXQtY3VydmUuanMiLCJub2RlX21vZHVsZXMvc3ZnLmpzL2Rpc3Qvc3ZnLmpzIiwic3JjL0NvbnRyb2xzL0dseXBoc0NvbnRyb2wuanMiLCJzcmMvQ29udHJvbHMvUGFpbnRDb250cm9sLmpzIiwic3JjL21haW4uanMiXSwibmFtZXMiOltdLCJtYXBwaW5ncyI6IkFBQUE7QUNBQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTs7QUMxakJBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTs7Ozs7Ozs7QUMvMUtBOzs7O0FBQ0E7Ozs7OztBQUVBLFNBQVMsYUFBVCxDQUF1QixNQUF2QixFQUErQjs7QUFFOUIsTUFBSyxLQUFMLEdBQWEsVUFBVSxLQUFWLEVBQWtCO0FBQzlCLGVBQWEsSUFBYixDQUFtQixLQUFuQjtBQUNNLHFCQUFtQixPQUFPLFFBQVAsR0FBa0IsSUFBbEIsQ0FBdUIsTUFBdkIsRUFBK0IsTUFBL0IsQ0FBc0MsRUFBRSxPQUFPLENBQVQsRUFBdEMsQ0FBbkI7QUFFTixFQUpEOztBQU1BLE1BQUssTUFBTCxHQUFjLFVBQVUsS0FBVixFQUFrQjtBQUMvQixlQUFhLElBQWIsQ0FBbUIsS0FBbkI7QUFDTSxjQUFhLGdCQUFiLEVBQStCLFlBQS9CO0FBQ04sRUFIRDs7QUFLQSxNQUFLLEdBQUwsR0FBVyxZQUFXO0FBQ3JCLE1BQUksY0FBYyx3QkFBVSxZQUFWLEVBQXdCLEtBQXhCLENBQWxCO0FBQ00sTUFBSSxhQUFhLDRCQUE0QixXQUE1QixDQUFqQjtBQUNBLFNBQU8sSUFBUCxDQUFhLFVBQWIsRUFBMEIsSUFBMUIsQ0FBK0IsTUFBL0IsRUFBdUMsTUFBdkMsQ0FBOEMsRUFBRSxPQUFPLENBQVQsRUFBOUMsRUFBNEQsTUFBNUQsQ0FBbUUsTUFBbkU7QUFDQSxlQUFhLE1BQWIsR0FBc0IsQ0FBdEI7QUFDTixFQUxEO0FBTUE7O2tCQUVjLGE7Ozs7Ozs7OztBQ3hCZjs7OztBQUNBOzs7Ozs7QUFFQSxJQUFNLFFBQVEsR0FBZDs7QUFFQSxTQUFTLFlBQVQsQ0FBc0IsTUFBdEIsRUFBOEI7QUFDN0IsTUFBSSxlQUFlLEVBQW5CO0FBQ0EsTUFBSSxtQkFBbUIsU0FBdkI7O0FBRUEsT0FBSyxLQUFMLEdBQWEsVUFBVSxLQUFWLEVBQWtCO0FBQzlCLGlCQUFhLElBQWIsQ0FBbUIsS0FBbkI7QUFDTSx1QkFBbUIsT0FBTyxRQUFQLEdBQWtCLElBQWxCLENBQXVCLE1BQXZCLEVBQStCLE1BQS9CLENBQXNDLEVBQUUsT0FBTyxDQUFULEVBQXRDLENBQW5CO0FBRU4sR0FKRDtBQUtBLE9BQUssTUFBTCxHQUFjLFVBQVUsS0FBVixFQUFrQjtBQUMvQixpQkFBYSxJQUFiLENBQW1CLEtBQW5CO0FBQ00sZ0JBQWEsZ0JBQWIsRUFBK0IsWUFBL0I7QUFDTixHQUhEOztBQUtBLE9BQUssR0FBTCxHQUFXLFlBQVc7QUFDckIsUUFBSSxjQUFjLHdCQUFVLFlBQVYsRUFBd0IsS0FBeEIsQ0FBbEI7QUFDTSxRQUFJLGFBQWEsNEJBQTRCLFdBQTVCLENBQWpCO0FBQ0EsV0FBTyxJQUFQLENBQWEsVUFBYixFQUEwQixJQUExQixDQUErQixNQUEvQixFQUF1QyxNQUF2QyxDQUE4QyxFQUFFLE9BQU8sQ0FBVCxFQUE5QyxFQUE0RCxNQUE1RCxDQUFtRSxNQUFuRTtBQUNBLGlCQUFhLE1BQWIsR0FBc0IsQ0FBdEI7QUFDTixHQUxEO0FBTUE7O2tCQUVjLFk7OztBQUVmLFNBQVMsV0FBVCxDQUFxQixnQkFBckIsRUFBdUMsWUFBdkMsRUFBcUQ7QUFDcEQsbUJBQWlCLElBQWpCLENBQXVCLFlBQXZCO0FBQ0E7QUFDRCxTQUFTLDJCQUFULENBQXFDLGNBQXJDLEVBQXFEO0FBQ2pELE1BQUksTUFBTSxFQUFWO0FBQ0EsaUJBQWUsR0FBZixDQUFtQixVQUFVLE1BQVYsRUFBa0IsQ0FBbEIsRUFBcUI7QUFDcEMsUUFBSSxLQUFLLENBQVQsRUFBWTtBQUNSLGFBQU8sT0FBTyxPQUFPLENBQVAsRUFBVSxDQUFWLENBQVAsR0FBc0IsR0FBdEIsR0FBNEIsT0FBTyxDQUFQLEVBQVUsQ0FBVixDQUFuQztBQUNIO0FBQ0QsV0FBTyxPQUFPLE9BQU8sQ0FBUCxFQUFVLENBQVYsQ0FBUCxHQUFzQixHQUF0QixHQUE0QixPQUFPLENBQVAsRUFBVSxDQUFWLENBQTVCLEdBQTJDLElBQTNDLEdBQ0gsT0FBTyxDQUFQLEVBQVUsQ0FBVixDQURHLEdBQ1ksR0FEWixHQUNrQixPQUFPLENBQVAsRUFBVSxDQUFWLENBRGxCLEdBQ2lDLElBRGpDLEdBRUgsT0FBTyxDQUFQLEVBQVUsQ0FBVixDQUZHLEdBRVksR0FGWixHQUVrQixPQUFPLENBQVAsRUFBVSxDQUFWLENBRmxCLEdBRWlDLEdBRnhDO0FBR0gsR0FQRDs7QUFTQSxTQUFPLEdBQVA7QUFDSDs7Ozs7QUM1Q0Q7Ozs7QUFDQTs7OztBQUNBOzs7Ozs7QUFFQSxJQUFJLE9BQU8sbUJBQUksU0FBSixFQUFlLElBQWYsQ0FBb0IsR0FBcEIsRUFBeUIsR0FBekIsQ0FBWDs7QUFHQSxXQUFXLElBQVg7O0FBSUEsU0FBUyxVQUFULENBQW9CLFVBQXBCLEVBQWdDO0FBQy9CLFFBQUksY0FBYyxLQUFsQjtBQUNBLFFBQUksaUJBQUo7O0FBRUEsUUFBSSxpQkFBaUIsMkJBQWlCLElBQWpCLENBQXJCOztBQUVHLGVBQVcsRUFBWCxDQUFjLFdBQWQsRUFBMkIsVUFBVSxDQUFWLEVBQWE7QUFDcEMsWUFBTSxRQUFRLENBQ2IsRUFBRSxPQURXLEVBRWIsRUFBRSxPQUZXLENBQWQ7QUFJQSxzQkFBYyxJQUFkO0FBQ0EsdUJBQWUsS0FBZixDQUFxQixLQUFyQjtBQUVILEtBUkQ7QUFTQSxlQUFXLEVBQVgsQ0FBYyxTQUFkLEVBQXlCLFlBQVk7QUFDakMsc0JBQWMsS0FBZDtBQUNBLHVCQUFlLEdBQWY7QUFDSCxLQUhEO0FBSUEsZUFBVyxFQUFYLENBQWMsV0FBZCxFQUEyQixVQUFVLENBQVYsRUFBYTtBQUNwQyxZQUFJLElBQUksRUFBRSxPQUFWO0FBQ0EsWUFBSSxJQUFJLEVBQUUsT0FBVjtBQUNBLFlBQUksV0FBSixFQUFpQjtBQUNiLDJCQUFlLE1BQWYsQ0FBc0IsQ0FBQyxDQUFELEVBQUksQ0FBSixDQUF0QjtBQUNIO0FBQ0osS0FORDtBQU9IIiwiZmlsZSI6ImdlbmVyYXRlZC5qcyIsInNvdXJjZVJvb3QiOiIiLCJzb3VyY2VzQ29udGVudCI6WyIoZnVuY3Rpb24gZSh0LG4scil7ZnVuY3Rpb24gcyhvLHUpe2lmKCFuW29dKXtpZighdFtvXSl7dmFyIGE9dHlwZW9mIHJlcXVpcmU9PVwiZnVuY3Rpb25cIiYmcmVxdWlyZTtpZighdSYmYSlyZXR1cm4gYShvLCEwKTtpZihpKXJldHVybiBpKG8sITApO3ZhciBmPW5ldyBFcnJvcihcIkNhbm5vdCBmaW5kIG1vZHVsZSAnXCIrbytcIidcIik7dGhyb3cgZi5jb2RlPVwiTU9EVUxFX05PVF9GT1VORFwiLGZ9dmFyIGw9bltvXT17ZXhwb3J0czp7fX07dFtvXVswXS5jYWxsKGwuZXhwb3J0cyxmdW5jdGlvbihlKXt2YXIgbj10W29dWzFdW2VdO3JldHVybiBzKG4/bjplKX0sbCxsLmV4cG9ydHMsZSx0LG4scil9cmV0dXJuIG5bb10uZXhwb3J0c312YXIgaT10eXBlb2YgcmVxdWlyZT09XCJmdW5jdGlvblwiJiZyZXF1aXJlO2Zvcih2YXIgbz0wO288ci5sZW5ndGg7bysrKXMocltvXSk7cmV0dXJuIHN9KSIsIi8vID09Q2xvc3VyZUNvbXBpbGVyPT1cclxuLy8gQG91dHB1dF9maWxlX25hbWUgZml0LWN1cnZlLm1pbi5qc1xyXG4vLyBAY29tcGlsYXRpb25fbGV2ZWwgU0lNUExFX09QVElNSVpBVElPTlNcclxuLy8gPT0vQ2xvc3VyZUNvbXBpbGVyPT1cclxuXHJcbi8qKlxyXG4gKiAgQHByZXNlcnZlICBKYXZhU2NyaXB0IGltcGxlbWVudGF0aW9uIG9mXHJcbiAqICBBbGdvcml0aG0gZm9yIEF1dG9tYXRpY2FsbHkgRml0dGluZyBEaWdpdGl6ZWQgQ3VydmVzXHJcbiAqICBieSBQaGlsaXAgSi4gU2NobmVpZGVyXHJcbiAqICBcIkdyYXBoaWNzIEdlbXNcIiwgQWNhZGVtaWMgUHJlc3MsIDE5OTBcclxuICpcclxuICogIFRoZSBNSVQgTGljZW5zZSAoTUlUKVxyXG4gKlxyXG4gKiAgaHR0cHM6Ly9naXRodWIuY29tL3Nvc3dvdy9maXQtY3VydmVzXHJcbiAqL1xyXG5cclxuLyoqXHJcbiAqIEZpdCBvbmUgb3IgbW9yZSBCZXppZXIgY3VydmVzIHRvIGEgc2V0IG9mIHBvaW50cy5cclxuICpcclxuICogQHBhcmFtIHtBcnJheTxBcnJheTxOdW1iZXI+Pn0gcG9pbnRzIC0gQXJyYXkgb2YgZGlnaXRpemVkIHBvaW50cywgZS5nLiBbWzUsNV0sWzUsNTBdLFsxMTAsMTQwXSxbMjEwLDE2MF0sWzMyMCwxMTBdXVxyXG4gKiBAcGFyYW0ge051bWJlcn0gbWF4RXJyb3IgLSBUb2xlcmFuY2UsIHNxdWFyZWQgZXJyb3IgYmV0d2VlbiBwb2ludHMgYW5kIGZpdHRlZCBjdXJ2ZVxyXG4gKiBAcmV0dXJucyB7QXJyYXk8QXJyYXk8QXJyYXk8TnVtYmVyPj4+fSBBcnJheSBvZiBCZXppZXIgY3VydmVzLCB3aGVyZSBlYWNoIGVsZW1lbnQgaXMgW2ZpcnN0LXBvaW50LCBjb250cm9sLXBvaW50LTEsIGNvbnRyb2wtcG9pbnQtMiwgc2Vjb25kLXBvaW50XSBhbmQgcG9pbnRzIGFyZSBbeCwgeV1cclxuICovXHJcbmZ1bmN0aW9uIGZpdEN1cnZlKHBvaW50cywgbWF4RXJyb3IsIHByb2dyZXNzQ2FsbGJhY2spIHtcclxuICAgIGlmICghQXJyYXkuaXNBcnJheShwb2ludHMpKSB7XHJcbiAgICAgICAgdGhyb3cgbmV3IFR5cGVFcnJvcihcIkZpcnN0IGFyZ3VtZW50IHNob3VsZCBiZSBhbiBhcnJheVwiKTtcclxuICAgIH1cclxuICAgIHBvaW50cy5mb3JFYWNoKChwb2ludCkgPT4ge1xyXG4gICAgICAgIGlmKCFBcnJheS5pc0FycmF5KHBvaW50KSB8fCBwb2ludC5sZW5ndGggIT09IDJcclxuICAgICAgICB8fCB0eXBlb2YgcG9pbnRbMF0gIT09ICdudW1iZXInIHx8IHR5cGVvZiBwb2ludFsxXSAhPT0gJ251bWJlcicpe1xyXG4gICAgICAgICAgICB0aHJvdyBFcnJvcihcIkVhY2ggcG9pbnQgc2hvdWxkIGJlIGFuIGFycmF5IG9mIHR3byBudW1iZXJzXCIpXHJcbiAgICAgICAgfVxyXG4gICAgfSk7XHJcbiAgICAvLyBSZW1vdmUgZHVwbGljYXRlIHBvaW50c1xyXG4gICAgcG9pbnRzID0gcG9pbnRzLmZpbHRlcigocG9pbnQsIGkpID0+XHJcbiAgICAgICAgaSA9PT0gMCB8fCAhKHBvaW50WzBdID09PSBwb2ludHNbaS0xXVswXSAmJiBwb2ludFsxXSA9PT0gcG9pbnRzW2ktMV1bMV0pXHJcbiAgICApO1xyXG5cclxuICAgIGlmIChwb2ludHMubGVuZ3RoIDwgMikge1xyXG4gICAgICAgIHJldHVybiBbXTtcclxuICAgIH1cclxuXHJcbiAgICBjb25zdCBsZW4gPSBwb2ludHMubGVuZ3RoO1xyXG4gICAgY29uc3QgbGVmdFRhbmdlbnQgPSBjcmVhdGVUYW5nZW50KHBvaW50c1sxXSwgcG9pbnRzWzBdKTtcclxuICAgIGNvbnN0IHJpZ2h0VGFuZ2VudCA9IGNyZWF0ZVRhbmdlbnQocG9pbnRzW2xlbiAtIDJdLCBwb2ludHNbbGVuIC0gMV0pO1xyXG5cclxuICAgIHJldHVybiBmaXRDdWJpYyhwb2ludHMsIGxlZnRUYW5nZW50LCByaWdodFRhbmdlbnQsIG1heEVycm9yLCBwcm9ncmVzc0NhbGxiYWNrKTtcclxufVxyXG5cclxuLyoqXHJcbiAqIEZpdCBhIEJlemllciBjdXJ2ZSB0byBhIChzdWIpc2V0IG9mIGRpZ2l0aXplZCBwb2ludHMuXHJcbiAqIFlvdXIgY29kZSBzaG91bGQgbm90IGNhbGwgdGhpcyBmdW5jdGlvbiBkaXJlY3RseS4gVXNlIHtAbGluayBmaXRDdXJ2ZX0gaW5zdGVhZC5cclxuICpcclxuICogQHBhcmFtIHtBcnJheTxBcnJheTxOdW1iZXI+Pn0gcG9pbnRzIC0gQXJyYXkgb2YgZGlnaXRpemVkIHBvaW50cywgZS5nLiBbWzUsNV0sWzUsNTBdLFsxMTAsMTQwXSxbMjEwLDE2MF0sWzMyMCwxMTBdXVxyXG4gKiBAcGFyYW0ge0FycmF5PE51bWJlcj59IGxlZnRUYW5nZW50IC0gVW5pdCB0YW5nZW50IHZlY3RvciBhdCBzdGFydCBwb2ludFxyXG4gKiBAcGFyYW0ge0FycmF5PE51bWJlcj59IHJpZ2h0VGFuZ2VudCAtIFVuaXQgdGFuZ2VudCB2ZWN0b3IgYXQgZW5kIHBvaW50XHJcbiAqIEBwYXJhbSB7TnVtYmVyfSBlcnJvciAtIFRvbGVyYW5jZSwgc3F1YXJlZCBlcnJvciBiZXR3ZWVuIHBvaW50cyBhbmQgZml0dGVkIGN1cnZlXHJcbiAqIEByZXR1cm5zIHtBcnJheTxBcnJheTxBcnJheTxOdW1iZXI+Pj59IEFycmF5IG9mIEJlemllciBjdXJ2ZXMsIHdoZXJlIGVhY2ggZWxlbWVudCBpcyBbZmlyc3QtcG9pbnQsIGNvbnRyb2wtcG9pbnQtMSwgY29udHJvbC1wb2ludC0yLCBzZWNvbmQtcG9pbnRdIGFuZCBwb2ludHMgYXJlIFt4LCB5XVxyXG4gKi9cclxuZnVuY3Rpb24gZml0Q3ViaWMocG9pbnRzLCBsZWZ0VGFuZ2VudCwgcmlnaHRUYW5nZW50LCBlcnJvciwgcHJvZ3Jlc3NDYWxsYmFjaykge1xyXG4gICAgY29uc3QgTWF4SXRlcmF0aW9ucyA9IDIwOyAgIC8vTWF4IHRpbWVzIHRvIHRyeSBpdGVyYXRpbmcgKHRvIGZpbmQgYW4gYWNjZXB0YWJsZSBjdXJ2ZSlcclxuXHJcbiAgICB2YXIgYmV6Q3VydmUsICAgICAgICAgICAgICAgLy9Db250cm9sIHBvaW50cyBvZiBmaXR0ZWQgQmV6aWVyIGN1cnZlXHJcbiAgICAgICAgdSwgICAgICAgICAgICAgICAgICAgICAgLy9QYXJhbWV0ZXIgdmFsdWVzIGZvciBwb2ludFxyXG4gICAgICAgIHVQcmltZSwgICAgICAgICAgICAgICAgIC8vSW1wcm92ZWQgcGFyYW1ldGVyIHZhbHVlc1xyXG4gICAgICAgIG1heEVycm9yLCBwcmV2RXJyLCAgICAgIC8vTWF4aW11bSBmaXR0aW5nIGVycm9yXHJcbiAgICAgICAgc3BsaXRQb2ludCwgcHJldlNwbGl0LCAgLy9Qb2ludCB0byBzcGxpdCBwb2ludCBzZXQgYXQgaWYgd2UgbmVlZCBtb3JlIHRoYW4gb25lIGN1cnZlXHJcbiAgICAgICAgY2VudGVyVmVjdG9yLCB0b0NlbnRlclRhbmdlbnQsIGZyb21DZW50ZXJUYW5nZW50LCAgLy9Vbml0IHRhbmdlbnQgdmVjdG9yKHMpIGF0IHNwbGl0UG9pbnRcclxuICAgICAgICBiZXppZXJzLCAgICAgICAgICAgICAgICAvL0FycmF5IG9mIGZpdHRlZCBCZXppZXIgY3VydmVzIGlmIHdlIG5lZWQgbW9yZSB0aGFuIG9uZSBjdXJ2ZVxyXG4gICAgICAgIGRpc3QsIGk7XHJcblxyXG4gICAgLy9jb25zb2xlLmxvZygnZml0Q3ViaWMsICcsIHBvaW50cy5sZW5ndGgpO1xyXG5cclxuICAgIC8vVXNlIGhldXJpc3RpYyBpZiByZWdpb24gb25seSBoYXMgdHdvIHBvaW50cyBpbiBpdFxyXG4gICAgaWYgKHBvaW50cy5sZW5ndGggPT09IDIpIHtcclxuICAgICAgICBkaXN0ID0gbWF0aHMudmVjdG9yTGVuKG1hdGhzLnN1YnRyYWN0KHBvaW50c1swXSwgcG9pbnRzWzFdKSkgLyAzLjA7XHJcbiAgICAgICAgYmV6Q3VydmUgPSBbXHJcbiAgICAgICAgICAgIHBvaW50c1swXSxcclxuICAgICAgICAgICAgbWF0aHMuYWRkQXJyYXlzKHBvaW50c1swXSwgbWF0aHMubXVsSXRlbXMobGVmdFRhbmdlbnQsICBkaXN0KSksXHJcbiAgICAgICAgICAgIG1hdGhzLmFkZEFycmF5cyhwb2ludHNbMV0sIG1hdGhzLm11bEl0ZW1zKHJpZ2h0VGFuZ2VudCwgZGlzdCkpLFxyXG4gICAgICAgICAgICBwb2ludHNbMV1cclxuICAgICAgICBdO1xyXG4gICAgICAgIHJldHVybiBbYmV6Q3VydmVdO1xyXG4gICAgfVxyXG5cclxuICAgIC8vUGFyYW1ldGVyaXplIHBvaW50cywgYW5kIGF0dGVtcHQgdG8gZml0IGN1cnZlXHJcbiAgICB1ID0gY2hvcmRMZW5ndGhQYXJhbWV0ZXJpemUocG9pbnRzKTtcclxuICAgIFtiZXpDdXJ2ZSwgbWF4RXJyb3IsIHNwbGl0UG9pbnRdID0gZ2VuZXJhdGVBbmRSZXBvcnQocG9pbnRzLCB1LCB1LCBsZWZ0VGFuZ2VudCwgcmlnaHRUYW5nZW50LCBwcm9ncmVzc0NhbGxiYWNrKVxyXG5cclxuICAgIGlmIChtYXhFcnJvciA8IGVycm9yKSB7XHJcbiAgICAgICAgcmV0dXJuIFtiZXpDdXJ2ZV07XHJcbiAgICB9XHJcbiAgICAvL0lmIGVycm9yIG5vdCB0b28gbGFyZ2UsIHRyeSBzb21lIHJlcGFyYW1ldGVyaXphdGlvbiBhbmQgaXRlcmF0aW9uXHJcbiAgICBpZiAobWF4RXJyb3IgPCAoZXJyb3IqZXJyb3IpKSB7XHJcblxyXG4gICAgICAgIHVQcmltZSA9IHU7XHJcbiAgICAgICAgcHJldkVyciA9IG1heEVycm9yO1xyXG4gICAgICAgIHByZXZTcGxpdCA9IHNwbGl0UG9pbnQ7XHJcblxyXG4gICAgICAgIGZvciAoaSA9IDA7IGkgPCBNYXhJdGVyYXRpb25zOyBpKyspIHtcclxuXHJcbiAgICAgICAgICAgIHVQcmltZSA9IHJlcGFyYW1ldGVyaXplKGJlekN1cnZlLCBwb2ludHMsIHVQcmltZSk7XHJcbiAgICAgICAgICAgIFtiZXpDdXJ2ZSwgbWF4RXJyb3IsIHNwbGl0UG9pbnRdID0gZ2VuZXJhdGVBbmRSZXBvcnQocG9pbnRzLCB1LCB1UHJpbWUsIGxlZnRUYW5nZW50LCByaWdodFRhbmdlbnQsIHByb2dyZXNzQ2FsbGJhY2spO1xyXG5cclxuICAgICAgICAgICAgaWYgKG1heEVycm9yIDwgZXJyb3IpIHtcclxuICAgICAgICAgICAgICAgIHJldHVybiBbYmV6Q3VydmVdO1xyXG4gICAgICAgICAgICB9XHJcbiAgICAgICAgICAgIC8vSWYgdGhlIGRldmVsb3BtZW50IG9mIHRoZSBmaXR0ZWQgY3VydmUgZ3JpbmRzIHRvIGEgaGFsdCxcclxuICAgICAgICAgICAgLy93ZSBhYm9ydCB0aGlzIGF0dGVtcHQgKGFuZCB0cnkgYSBzaG9ydGVyIGN1cnZlKTpcclxuICAgICAgICAgICAgZWxzZSBpZihzcGxpdFBvaW50ID09PSBwcmV2U3BsaXQpIHtcclxuICAgICAgICAgICAgICAgIGxldCBlcnJDaGFuZ2UgPSBtYXhFcnJvci9wcmV2RXJyO1xyXG4gICAgICAgICAgICAgICAgaWYoKGVyckNoYW5nZSA+IC45OTk5KSAmJiAoZXJyQ2hhbmdlIDwgMS4wMDAxKSkge1xyXG4gICAgICAgICAgICAgICAgICAgIGJyZWFrO1xyXG4gICAgICAgICAgICAgICAgfVxyXG4gICAgICAgICAgICB9XHJcblxyXG4gICAgICAgICAgICBwcmV2RXJyID0gbWF4RXJyb3I7XHJcbiAgICAgICAgICAgIHByZXZTcGxpdCA9IHNwbGl0UG9pbnQ7XHJcbiAgICAgICAgfVxyXG4gICAgfVxyXG5cclxuICAgIC8vRml0dGluZyBmYWlsZWQgLS0gc3BsaXQgYXQgbWF4IGVycm9yIHBvaW50IGFuZCBmaXQgcmVjdXJzaXZlbHlcclxuICAgIGJlemllcnMgPSBbXTtcclxuXHJcbiAgICAvL1RvIGNyZWF0ZSBhIHNtb290aCB0cmFuc2l0aW9uIGZyb20gb25lIGN1cnZlIHNlZ21lbnQgdG8gdGhlIG5leHQsXHJcbiAgICAvL3dlIGNhbGN1bGF0ZSB0aGUgdGFuZ2VudCBvZiB0aGUgcG9pbnRzIGRpcmVjdGx5IGJlZm9yZSBhbmQgYWZ0ZXIgdGhlIGNlbnRlcixcclxuICAgIC8vYW5kIHVzZSB0aGF0IHNhbWUgdGFuZ2VudCBib3RoIHRvIGFuZCBmcm9tIHRoZSBjZW50ZXIgcG9pbnQuXHJcbiAgICBjZW50ZXJWZWN0b3IgPSBtYXRocy5zdWJ0cmFjdChwb2ludHNbc3BsaXRQb2ludCAtIDFdLCBwb2ludHNbc3BsaXRQb2ludCArIDFdKTtcclxuICAgIC8vSG93ZXZlciwgc2hvdWxkIHRob3NlIHR3byBwb2ludHMgYmUgZXF1YWwsIHRoZSBub3JtYWwgdGFuZ2VudCBjYWxjdWxhdGlvbiB3aWxsIGZhaWwuXHJcbiAgICAvL0luc3RlYWQsIHdlIGNhbGN1bGF0ZSB0aGUgdGFuZ2VudCBmcm9tIHRoYXQgXCJkb3VibGUtcG9pbnRcIiB0byB0aGUgY2VudGVyIHBvaW50LCBhbmQgcm90YXRlIDkwZGVnLlxyXG4gICAgaWYoKGNlbnRlclZlY3RvclswXSA9PT0gMCkgJiYgKGNlbnRlclZlY3RvclsxXSA9PT0gMCkpIHtcclxuICAgICAgICAvL3RvQ2VudGVyVGFuZ2VudCA9IGNyZWF0ZVRhbmdlbnQocG9pbnRzW3NwbGl0UG9pbnQgLSAxXSwgcG9pbnRzW3NwbGl0UG9pbnRdKTtcclxuICAgICAgICAvL2Zyb21DZW50ZXJUYW5nZW50ID0gY3JlYXRlVGFuZ2VudChwb2ludHNbc3BsaXRQb2ludCArIDFdLCBwb2ludHNbc3BsaXRQb2ludF0pO1xyXG5cclxuICAgICAgICAvL1t4LHldIC0+IFsteSx4XTogaHR0cDovL3N0YWNrb3ZlcmZsb3cuY29tL2EvNDc4MDE0MS8xODY5NjYwXHJcbiAgICAgICAgY2VudGVyVmVjdG9yID0gbWF0aHMuc3VidHJhY3QocG9pbnRzW3NwbGl0UG9pbnQgLSAxXSwgcG9pbnRzW3NwbGl0UG9pbnRdKVxyXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgLnJldmVyc2UoKTtcclxuICAgICAgICBjZW50ZXJWZWN0b3JbMF0gPSAtY2VudGVyVmVjdG9yWzBdO1xyXG4gICAgfVxyXG4gICAgdG9DZW50ZXJUYW5nZW50ID0gbWF0aHMubm9ybWFsaXplKGNlbnRlclZlY3Rvcik7XHJcbiAgICAvL1RvIGFuZCBmcm9tIG5lZWQgdG8gcG9pbnQgaW4gb3Bwb3NpdGUgZGlyZWN0aW9uczpcclxuICAgIGZyb21DZW50ZXJUYW5nZW50ID0gbWF0aHMubXVsSXRlbXModG9DZW50ZXJUYW5nZW50LCAtMSk7XHJcblxyXG4gICAgLypcclxuICAgIE5vdGU6IEFuIGFsdGVybmF0aXZlIHRvIHRoaXMgXCJkaXZpZGUgYW5kIGNvbnF1ZXJcIiByZWN1cnNpb24gY291bGQgYmUgdG8gYWx3YXlzXHJcbiAgICAgICAgICBsZXQgbmV3IGN1cnZlIHNlZ21lbnRzIHN0YXJ0IGJ5IHRyeWluZyB0byBnbyBhbGwgdGhlIHdheSB0byB0aGUgZW5kLFxyXG4gICAgICAgICAgaW5zdGVhZCBvZiBvbmx5IHRvIHRoZSBlbmQgb2YgdGhlIGN1cnJlbnQgc3ViZGl2aWRlZCBwb2x5bGluZS5cclxuICAgICAgICAgIFRoYXQgbWlnaHQgbGV0IG1hbnkgc2VnbWVudHMgZml0IGEgZmV3IHBvaW50cyBtb3JlLCByZWR1Y2luZyB0aGUgbnVtYmVyIG9mIHRvdGFsIHNlZ21lbnRzLlxyXG5cclxuICAgICAgICAgIEhvd2V2ZXIsIGEgZmV3IHRlc3RzIGhhdmUgc2hvd24gdGhhdCB0aGUgc2VnbWVudCByZWR1Y3Rpb24gaXMgaW5zaWduaWZpY2FudFxyXG4gICAgICAgICAgKDI0MCBwdHMsIDEwMCBlcnI6IDI1IGN1cnZlcyB2cyAyNyBjdXJ2ZXMuIDE0MCBwdHMsIDEwMCBlcnI6IDE3IGN1cnZlcyBvbiBib3RoKSxcclxuICAgICAgICAgIGFuZCB0aGUgcmVzdWx0cyB0YWtlIHR3aWNlIGFzIG1hbnkgc3RlcHMgYW5kIG1pbGxpc2Vjb25kcyB0byBmaW5pc2gsXHJcbiAgICAgICAgICB3aXRob3V0IGxvb2tpbmcgYW55IGJldHRlciB0aGFuIHdoYXQgd2UgYWxyZWFkeSBoYXZlLlxyXG4gICAgKi9cclxuICAgIGJlemllcnMgPSBiZXppZXJzLmNvbmNhdChmaXRDdWJpYyhwb2ludHMuc2xpY2UoMCwgc3BsaXRQb2ludCArIDEpLCBsZWZ0VGFuZ2VudCwgdG9DZW50ZXJUYW5nZW50LCAgICBlcnJvciwgcHJvZ3Jlc3NDYWxsYmFjaykpO1xyXG4gICAgYmV6aWVycyA9IGJlemllcnMuY29uY2F0KGZpdEN1YmljKHBvaW50cy5zbGljZShzcGxpdFBvaW50KSwgICAgICAgIGZyb21DZW50ZXJUYW5nZW50LCByaWdodFRhbmdlbnQsIGVycm9yLCBwcm9ncmVzc0NhbGxiYWNrKSk7XHJcbiAgICByZXR1cm4gYmV6aWVycztcclxufTtcclxuXHJcbmZ1bmN0aW9uIGdlbmVyYXRlQW5kUmVwb3J0KHBvaW50cywgcGFyYW1zT3JpZywgcGFyYW1zUHJpbWUsIGxlZnRUYW5nZW50LCByaWdodFRhbmdlbnQsIHByb2dyZXNzQ2FsbGJhY2spIHtcclxuICAgIHZhciBiZXpDdXJ2ZSwgbWF4RXJyb3IsIHNwbGl0UG9pbnQ7XHJcblxyXG4gICAgYmV6Q3VydmUgPSBnZW5lcmF0ZUJlemllcihwb2ludHMsIHBhcmFtc1ByaW1lLCBsZWZ0VGFuZ2VudCwgcmlnaHRUYW5nZW50LCBwcm9ncmVzc0NhbGxiYWNrKTtcclxuICAgIC8vRmluZCBtYXggZGV2aWF0aW9uIG9mIHBvaW50cyB0byBmaXR0ZWQgY3VydmUuXHJcbiAgICAvL0hlcmUgd2UgYWx3YXlzIHVzZSB0aGUgb3JpZ2luYWwgcGFyYW1ldGVycyAoZnJvbSBjaG9yZExlbmd0aFBhcmFtZXRlcml6ZSgpKSxcclxuICAgIC8vYmVjYXVzZSB3ZSBuZWVkIHRvIGNvbXBhcmUgdGhlIGN1cnJlbnQgY3VydmUgdG8gdGhlIGFjdHVhbCBzb3VyY2UgcG9seWxpbmUsXHJcbiAgICAvL2FuZCBub3QgdGhlIGN1cnJlbnRseSBpdGVyYXRlZCBwYXJhbWV0ZXJzIHdoaWNoIHJlcGFyYW1ldGVyaXplKCkgJiBnZW5lcmF0ZUJlemllcigpIHVzZSxcclxuICAgIC8vYXMgdGhvc2UgaGF2ZSBwcm9iYWJseSBkcmlmdGVkIGZhciBhd2F5IGFuZCBtYXkgbm8gbG9uZ2VyIGJlIGluIGFzY2VuZGluZyBvcmRlci5cclxuICAgIFttYXhFcnJvciwgc3BsaXRQb2ludF0gPSBjb21wdXRlTWF4RXJyb3IocG9pbnRzLCBiZXpDdXJ2ZSwgcGFyYW1zT3JpZyk7XHJcblxyXG4gICAgaWYocHJvZ3Jlc3NDYWxsYmFjaykge1xyXG4gICAgICAgIHByb2dyZXNzQ2FsbGJhY2soe1xyXG4gICAgICAgICAgICBiZXo6IGJlekN1cnZlLFxyXG4gICAgICAgICAgICBwb2ludHM6IHBvaW50cyxcclxuICAgICAgICAgICAgcGFyYW1zOiBwYXJhbXNPcmlnLFxyXG4gICAgICAgICAgICBtYXhFcnI6IG1heEVycm9yLFxyXG4gICAgICAgICAgICBtYXhQb2ludDogc3BsaXRQb2ludCxcclxuICAgICAgICB9KTtcclxuICAgIH1cclxuXHJcbiAgICByZXR1cm4gW2JlekN1cnZlLCBtYXhFcnJvciwgc3BsaXRQb2ludF07XHJcbn1cclxuXHJcbi8qKlxyXG4gKiBVc2UgbGVhc3Qtc3F1YXJlcyBtZXRob2QgdG8gZmluZCBCZXppZXIgY29udHJvbCBwb2ludHMgZm9yIHJlZ2lvbi5cclxuICpcclxuICogQHBhcmFtIHtBcnJheTxBcnJheTxOdW1iZXI+Pn0gcG9pbnRzIC0gQXJyYXkgb2YgZGlnaXRpemVkIHBvaW50c1xyXG4gKiBAcGFyYW0ge0FycmF5PE51bWJlcj59IHBhcmFtZXRlcnMgLSBQYXJhbWV0ZXIgdmFsdWVzIGZvciByZWdpb25cclxuICogQHBhcmFtIHtBcnJheTxOdW1iZXI+fSBsZWZ0VGFuZ2VudCAtIFVuaXQgdGFuZ2VudCB2ZWN0b3IgYXQgc3RhcnQgcG9pbnRcclxuICogQHBhcmFtIHtBcnJheTxOdW1iZXI+fSByaWdodFRhbmdlbnQgLSBVbml0IHRhbmdlbnQgdmVjdG9yIGF0IGVuZCBwb2ludFxyXG4gKiBAcmV0dXJucyB7QXJyYXk8QXJyYXk8TnVtYmVyPj59IEFwcHJveGltYXRlZCBCZXppZXIgY3VydmU6IFtmaXJzdC1wb2ludCwgY29udHJvbC1wb2ludC0xLCBjb250cm9sLXBvaW50LTIsIHNlY29uZC1wb2ludF0gd2hlcmUgcG9pbnRzIGFyZSBbeCwgeV1cclxuICovXHJcbmZ1bmN0aW9uIGdlbmVyYXRlQmV6aWVyKHBvaW50cywgcGFyYW1ldGVycywgbGVmdFRhbmdlbnQsIHJpZ2h0VGFuZ2VudCkge1xyXG4gICAgdmFyIGJlekN1cnZlLCAgICAgICAgICAgICAgICAgICAgICAgLy9CZXppZXIgY3VydmUgY3RsIHB0c1xyXG4gICAgICAgIEEsIGEsICAgICAgICAgICAgICAgICAgICAgICAgICAgLy9QcmVjb21wdXRlZCByaHMgZm9yIGVxblxyXG4gICAgICAgIEMsIFgsICAgICAgICAgICAgICAgICAgICAgICAgICAgLy9NYXRyaWNlcyBDICYgWFxyXG4gICAgICAgIGRldF9DMF9DMSwgZGV0X0MwX1gsIGRldF9YX0MxLCAgLy9EZXRlcm1pbmFudHMgb2YgbWF0cmljZXNcclxuICAgICAgICBhbHBoYV9sLCBhbHBoYV9yLCAgICAgICAgICAgICAgIC8vQWxwaGEgdmFsdWVzLCBsZWZ0IGFuZCByaWdodFxyXG5cclxuICAgICAgICBlcHNpbG9uLCBzZWdMZW5ndGgsXHJcbiAgICAgICAgaSwgbGVuLCB0bXAsIHUsIHV4LFxyXG4gICAgICAgIGZpcnN0UG9pbnQgPSBwb2ludHNbMF0sXHJcbiAgICAgICAgbGFzdFBvaW50ID0gcG9pbnRzW3BvaW50cy5sZW5ndGgtMV07XHJcblxyXG4gICAgYmV6Q3VydmUgPSBbZmlyc3RQb2ludCwgbnVsbCwgbnVsbCwgbGFzdFBvaW50XTtcclxuICAgIC8vY29uc29sZS5sb2coJ2diJywgcGFyYW1ldGVycy5sZW5ndGgpO1xyXG5cclxuICAgIC8vQ29tcHV0ZSB0aGUgQSdzXHJcbiAgICBBID0gbWF0aHMuemVyb3NfWHgyeDIocGFyYW1ldGVycy5sZW5ndGgpO1xyXG4gICAgZm9yIChpID0gMCwgbGVuID0gcGFyYW1ldGVycy5sZW5ndGg7IGkgPCBsZW47IGkrKykge1xyXG4gICAgICAgIHUgPSBwYXJhbWV0ZXJzW2ldO1xyXG4gICAgICAgIHV4ID0gMSAtIHU7XHJcbiAgICAgICAgYSA9IEFbaV07XHJcblxyXG4gICAgICAgIGFbMF0gPSBtYXRocy5tdWxJdGVtcyhsZWZ0VGFuZ2VudCwgIDMgKiB1ICAqICh1eCp1eCkpO1xyXG4gICAgICAgIGFbMV0gPSBtYXRocy5tdWxJdGVtcyhyaWdodFRhbmdlbnQsIDMgKiB1eCAqICh1KnUpKTtcclxuICAgIH1cclxuXHJcbiAgICAvL0NyZWF0ZSB0aGUgQyBhbmQgWCBtYXRyaWNlc1xyXG4gICAgQyA9IFtbMCwwXSwgWzAsMF1dO1xyXG4gICAgWCA9IFswLDBdO1xyXG4gICAgZm9yIChpID0gMCwgbGVuID0gcG9pbnRzLmxlbmd0aDsgaSA8IGxlbjsgaSsrKSB7XHJcbiAgICAgICAgdSA9IHBhcmFtZXRlcnNbaV07XHJcbiAgICAgICAgYSA9IEFbaV07XHJcblxyXG4gICAgICAgIENbMF1bMF0gKz0gbWF0aHMuZG90KGFbMF0sIGFbMF0pO1xyXG4gICAgICAgIENbMF1bMV0gKz0gbWF0aHMuZG90KGFbMF0sIGFbMV0pO1xyXG4gICAgICAgIENbMV1bMF0gKz0gbWF0aHMuZG90KGFbMF0sIGFbMV0pO1xyXG4gICAgICAgIENbMV1bMV0gKz0gbWF0aHMuZG90KGFbMV0sIGFbMV0pO1xyXG5cclxuICAgICAgICB0bXAgPSBtYXRocy5zdWJ0cmFjdChwb2ludHNbaV0sIGJlemllci5xKFtmaXJzdFBvaW50LCBmaXJzdFBvaW50LCBsYXN0UG9pbnQsIGxhc3RQb2ludF0sIHUpKTtcclxuXHJcbiAgICAgICAgWFswXSArPSBtYXRocy5kb3QoYVswXSwgdG1wKTtcclxuICAgICAgICBYWzFdICs9IG1hdGhzLmRvdChhWzFdLCB0bXApO1xyXG4gICAgfVxyXG5cclxuICAgIC8vQ29tcHV0ZSB0aGUgZGV0ZXJtaW5hbnRzIG9mIEMgYW5kIFhcclxuICAgIGRldF9DMF9DMSA9IChDWzBdWzBdICogQ1sxXVsxXSkgLSAoQ1sxXVswXSAqIENbMF1bMV0pO1xyXG4gICAgZGV0X0MwX1ggID0gKENbMF1bMF0gKiBYWzFdICAgKSAtIChDWzFdWzBdICogWFswXSAgICk7XHJcbiAgICBkZXRfWF9DMSAgPSAoWFswXSAgICAqIENbMV1bMV0pIC0gKFhbMV0gICAgKiBDWzBdWzFdKTtcclxuXHJcbiAgICAvL0ZpbmFsbHksIGRlcml2ZSBhbHBoYSB2YWx1ZXNcclxuICAgIGFscGhhX2wgPSBkZXRfQzBfQzEgPT09IDAgPyAwIDogZGV0X1hfQzEgLyBkZXRfQzBfQzE7XHJcbiAgICBhbHBoYV9yID0gZGV0X0MwX0MxID09PSAwID8gMCA6IGRldF9DMF9YIC8gZGV0X0MwX0MxO1xyXG5cclxuICAgIC8vSWYgYWxwaGEgbmVnYXRpdmUsIHVzZSB0aGUgV3UvQmFyc2t5IGhldXJpc3RpYyAoc2VlIHRleHQpLlxyXG4gICAgLy9JZiBhbHBoYSBpcyAwLCB5b3UgZ2V0IGNvaW5jaWRlbnQgY29udHJvbCBwb2ludHMgdGhhdCBsZWFkIHRvXHJcbiAgICAvL2RpdmlkZSBieSB6ZXJvIGluIGFueSBzdWJzZXF1ZW50IE5ld3RvblJhcGhzb25Sb290RmluZCgpIGNhbGwuXHJcbiAgICBzZWdMZW5ndGggPSBtYXRocy52ZWN0b3JMZW4obWF0aHMuc3VidHJhY3QoZmlyc3RQb2ludCwgbGFzdFBvaW50KSk7XHJcbiAgICBlcHNpbG9uID0gMS4wZS02ICogc2VnTGVuZ3RoO1xyXG4gICAgaWYgKGFscGhhX2wgPCBlcHNpbG9uIHx8IGFscGhhX3IgPCBlcHNpbG9uKSB7XHJcbiAgICAgICAgLy9GYWxsIGJhY2sgb24gc3RhbmRhcmQgKHByb2JhYmx5IGluYWNjdXJhdGUpIGZvcm11bGEsIGFuZCBzdWJkaXZpZGUgZnVydGhlciBpZiBuZWVkZWQuXHJcbiAgICAgICAgYmV6Q3VydmVbMV0gPSBtYXRocy5hZGRBcnJheXMoZmlyc3RQb2ludCwgbWF0aHMubXVsSXRlbXMobGVmdFRhbmdlbnQsICBzZWdMZW5ndGggLyAzLjApKTtcclxuICAgICAgICBiZXpDdXJ2ZVsyXSA9IG1hdGhzLmFkZEFycmF5cyhsYXN0UG9pbnQsICBtYXRocy5tdWxJdGVtcyhyaWdodFRhbmdlbnQsIHNlZ0xlbmd0aCAvIDMuMCkpO1xyXG4gICAgfSBlbHNlIHtcclxuICAgICAgICAvL0ZpcnN0IGFuZCBsYXN0IGNvbnRyb2wgcG9pbnRzIG9mIHRoZSBCZXppZXIgY3VydmUgYXJlXHJcbiAgICAgICAgLy9wb3NpdGlvbmVkIGV4YWN0bHkgYXQgdGhlIGZpcnN0IGFuZCBsYXN0IGRhdGEgcG9pbnRzXHJcbiAgICAgICAgLy9Db250cm9sIHBvaW50cyAxIGFuZCAyIGFyZSBwb3NpdGlvbmVkIGFuIGFscGhhIGRpc3RhbmNlIG91dFxyXG4gICAgICAgIC8vb24gdGhlIHRhbmdlbnQgdmVjdG9ycywgbGVmdCBhbmQgcmlnaHQsIHJlc3BlY3RpdmVseVxyXG4gICAgICAgIGJlekN1cnZlWzFdID0gbWF0aHMuYWRkQXJyYXlzKGZpcnN0UG9pbnQsIG1hdGhzLm11bEl0ZW1zKGxlZnRUYW5nZW50LCAgYWxwaGFfbCkpO1xyXG4gICAgICAgIGJlekN1cnZlWzJdID0gbWF0aHMuYWRkQXJyYXlzKGxhc3RQb2ludCwgIG1hdGhzLm11bEl0ZW1zKHJpZ2h0VGFuZ2VudCwgYWxwaGFfcikpO1xyXG4gICAgfVxyXG5cclxuICAgIHJldHVybiBiZXpDdXJ2ZTtcclxufTtcclxuXHJcbi8qKlxyXG4gKiBHaXZlbiBzZXQgb2YgcG9pbnRzIGFuZCB0aGVpciBwYXJhbWV0ZXJpemF0aW9uLCB0cnkgdG8gZmluZCBhIGJldHRlciBwYXJhbWV0ZXJpemF0aW9uLlxyXG4gKlxyXG4gKiBAcGFyYW0ge0FycmF5PEFycmF5PE51bWJlcj4+fSBiZXppZXIgLSBDdXJyZW50IGZpdHRlZCBjdXJ2ZVxyXG4gKiBAcGFyYW0ge0FycmF5PEFycmF5PE51bWJlcj4+fSBwb2ludHMgLSBBcnJheSBvZiBkaWdpdGl6ZWQgcG9pbnRzXHJcbiAqIEBwYXJhbSB7QXJyYXk8TnVtYmVyPn0gcGFyYW1ldGVycyAtIEN1cnJlbnQgcGFyYW1ldGVyIHZhbHVlc1xyXG4gKiBAcmV0dXJucyB7QXJyYXk8TnVtYmVyPn0gTmV3IHBhcmFtZXRlciB2YWx1ZXNcclxuICovXHJcbmZ1bmN0aW9uIHJlcGFyYW1ldGVyaXplKGJlemllciwgcG9pbnRzLCBwYXJhbWV0ZXJzKSB7XHJcbiAgICAvKlxyXG4gICAgdmFyIGosIGxlbiwgcG9pbnQsIHJlc3VsdHMsIHU7XHJcbiAgICByZXN1bHRzID0gW107XHJcbiAgICBmb3IgKGogPSAwLCBsZW4gPSBwb2ludHMubGVuZ3RoOyBqIDwgbGVuOyBqKyspIHtcclxuICAgICAgICBwb2ludCA9IHBvaW50c1tqXSwgdSA9IHBhcmFtZXRlcnNbal07XHJcblxyXG4gICAgICAgIHJlc3VsdHMucHVzaChuZXd0b25SYXBoc29uUm9vdEZpbmQoYmV6aWVyLCBwb2ludCwgdSkpO1xyXG4gICAgfVxyXG4gICAgcmV0dXJuIHJlc3VsdHM7XHJcbiAgICAvLyovXHJcbiAgICByZXR1cm4gcGFyYW1ldGVycy5tYXAoKHAsIGkpID0+IG5ld3RvblJhcGhzb25Sb290RmluZChiZXppZXIsIHBvaW50c1tpXSwgcCkpO1xyXG59O1xyXG5cclxuLyoqXHJcbiAqIFVzZSBOZXd0b24tUmFwaHNvbiBpdGVyYXRpb24gdG8gZmluZCBiZXR0ZXIgcm9vdC5cclxuICpcclxuICogQHBhcmFtIHtBcnJheTxBcnJheTxOdW1iZXI+Pn0gYmV6IC0gQ3VycmVudCBmaXR0ZWQgY3VydmVcclxuICogQHBhcmFtIHtBcnJheTxOdW1iZXI+fSBwb2ludCAtIERpZ2l0aXplZCBwb2ludFxyXG4gKiBAcGFyYW0ge051bWJlcn0gdSAtIFBhcmFtZXRlciB2YWx1ZSBmb3IgXCJQXCJcclxuICogQHJldHVybnMge051bWJlcn0gTmV3IHVcclxuICovXHJcbmZ1bmN0aW9uIG5ld3RvblJhcGhzb25Sb290RmluZChiZXosIHBvaW50LCB1KSB7XHJcbiAgICAvKlxyXG4gICAgICAgIE5ld3RvbidzIHJvb3QgZmluZGluZyBhbGdvcml0aG0gY2FsY3VsYXRlcyBmKHgpPTAgYnkgcmVpdGVyYXRpbmdcclxuICAgICAgICB4X24rMSA9IHhfbiAtIGYoeF9uKS9mJyh4X24pXHJcbiAgICAgICAgV2UgYXJlIHRyeWluZyB0byBmaW5kIGN1cnZlIHBhcmFtZXRlciB1IGZvciBzb21lIHBvaW50IHAgdGhhdCBtaW5pbWl6ZXNcclxuICAgICAgICB0aGUgZGlzdGFuY2UgZnJvbSB0aGF0IHBvaW50IHRvIHRoZSBjdXJ2ZS4gRGlzdGFuY2UgcG9pbnQgdG8gY3VydmUgaXMgZD1xKHUpLXAuXHJcbiAgICAgICAgQXQgbWluaW11bSBkaXN0YW5jZSB0aGUgcG9pbnQgaXMgcGVycGVuZGljdWxhciB0byB0aGUgY3VydmUuXHJcbiAgICAgICAgV2UgYXJlIHNvbHZpbmdcclxuICAgICAgICBmID0gcSh1KS1wICogcScodSkgPSAwXHJcbiAgICAgICAgd2l0aFxyXG4gICAgICAgIGYnID0gcScodSkgKiBxJyh1KSArIHEodSktcCAqIHEnJyh1KVxyXG4gICAgICAgIGdpdmVzXHJcbiAgICAgICAgdV9uKzEgPSB1X24gLSB8cSh1X24pLXAgKiBxJyh1X24pfCAvIHxxJyh1X24pKioyICsgcSh1X24pLXAgKiBxJycodV9uKXxcclxuICAgICovXHJcblxyXG4gICAgdmFyIGQgPSBtYXRocy5zdWJ0cmFjdChiZXppZXIucShiZXosIHUpLCBwb2ludCksXHJcbiAgICAgICAgcXByaW1lID0gYmV6aWVyLnFwcmltZShiZXosIHUpLFxyXG4gICAgICAgIG51bWVyYXRvciA9IC8qc3VtKCovbWF0aHMubXVsTWF0cml4KGQsIHFwcmltZSkvKikqLyxcclxuICAgICAgICBkZW5vbWluYXRvciA9IG1hdGhzLnN1bShtYXRocy5hZGRJdGVtcyggbWF0aHMuc3F1YXJlSXRlbXMocXByaW1lKSwgbWF0aHMubXVsTWF0cml4KGQsIGJlemllci5xcHJpbWVwcmltZShiZXosIHUpKSApKTtcclxuXHJcbiAgICBpZiAoZGVub21pbmF0b3IgPT09IDApIHtcclxuICAgICAgICByZXR1cm4gdTtcclxuICAgIH0gZWxzZSB7XHJcbiAgICAgICAgcmV0dXJuIHUgLSAobnVtZXJhdG9yL2Rlbm9taW5hdG9yKTtcclxuICAgIH1cclxufTtcclxuXHJcbi8qKlxyXG4gKiBBc3NpZ24gcGFyYW1ldGVyIHZhbHVlcyB0byBkaWdpdGl6ZWQgcG9pbnRzIHVzaW5nIHJlbGF0aXZlIGRpc3RhbmNlcyBiZXR3ZWVuIHBvaW50cy5cclxuICpcclxuICogQHBhcmFtIHtBcnJheTxBcnJheTxOdW1iZXI+Pn0gcG9pbnRzIC0gQXJyYXkgb2YgZGlnaXRpemVkIHBvaW50c1xyXG4gKiBAcmV0dXJucyB7QXJyYXk8TnVtYmVyPn0gUGFyYW1ldGVyIHZhbHVlc1xyXG4gKi9cclxuZnVuY3Rpb24gY2hvcmRMZW5ndGhQYXJhbWV0ZXJpemUocG9pbnRzKSB7XHJcbiAgICB2YXIgdSA9IFtdLCBjdXJyVSwgcHJldlUsIHByZXZQO1xyXG5cclxuICAgIHBvaW50cy5mb3JFYWNoKChwLCBpKSA9PiB7XHJcbiAgICAgICAgY3VyclUgPSBpID8gcHJldlUgKyBtYXRocy52ZWN0b3JMZW4obWF0aHMuc3VidHJhY3QocCwgcHJldlApKVxyXG4gICAgICAgICAgICAgICAgICA6IDA7XHJcbiAgICAgICAgdS5wdXNoKGN1cnJVKTtcclxuXHJcbiAgICAgICAgcHJldlUgPSBjdXJyVTtcclxuICAgICAgICBwcmV2UCA9IHA7XHJcbiAgICB9KVxyXG4gICAgdSA9IHUubWFwKHggPT4geC9wcmV2VSk7XHJcblxyXG4gICAgcmV0dXJuIHU7XHJcbn07XHJcblxyXG4vKipcclxuICogRmluZCB0aGUgbWF4aW11bSBzcXVhcmVkIGRpc3RhbmNlIG9mIGRpZ2l0aXplZCBwb2ludHMgdG8gZml0dGVkIGN1cnZlLlxyXG4gKlxyXG4gKiBAcGFyYW0ge0FycmF5PEFycmF5PE51bWJlcj4+fSBwb2ludHMgLSBBcnJheSBvZiBkaWdpdGl6ZWQgcG9pbnRzXHJcbiAqIEBwYXJhbSB7QXJyYXk8QXJyYXk8TnVtYmVyPj59IGJleiAtIEZpdHRlZCBjdXJ2ZVxyXG4gKiBAcGFyYW0ge0FycmF5PE51bWJlcj59IHBhcmFtZXRlcnMgLSBQYXJhbWV0ZXJpemF0aW9uIG9mIHBvaW50c1xyXG4gKiBAcmV0dXJucyB7QXJyYXk8TnVtYmVyPn0gTWF4aW11bSBlcnJvciAoc3F1YXJlZCkgYW5kIHBvaW50IG9mIG1heCBlcnJvclxyXG4gKi9cclxuZnVuY3Rpb24gY29tcHV0ZU1heEVycm9yKHBvaW50cywgYmV6LCBwYXJhbWV0ZXJzKSB7XHJcbiAgICB2YXIgZGlzdCwgICAgICAgLy9DdXJyZW50IGVycm9yXHJcbiAgICAgICAgbWF4RGlzdCwgICAgLy9NYXhpbXVtIGVycm9yXHJcbiAgICAgICAgc3BsaXRQb2ludCwgLy9Qb2ludCBvZiBtYXhpbXVtIGVycm9yXHJcbiAgICAgICAgdiwgICAgICAgICAgLy9WZWN0b3IgZnJvbSBwb2ludCB0byBjdXJ2ZVxyXG4gICAgICAgIGksIGNvdW50LCBwb2ludCwgdDtcclxuXHJcbiAgICBtYXhEaXN0ID0gMDtcclxuICAgIHNwbGl0UG9pbnQgPSBwb2ludHMubGVuZ3RoIC8gMjtcclxuXHJcbiAgICBjb25zdCB0X2Rpc3RNYXAgPSBtYXBUdG9SZWxhdGl2ZURpc3RhbmNlcyhiZXosIDEwKTtcclxuXHJcbiAgICBmb3IgKGkgPSAwLCBjb3VudCA9IHBvaW50cy5sZW5ndGg7IGkgPCBjb3VudDsgaSsrKSB7XHJcbiAgICAgICAgcG9pbnQgPSBwb2ludHNbaV07XHJcbiAgICAgICAgLy9GaW5kICd0JyBmb3IgYSBwb2ludCBvbiB0aGUgYmV6IGN1cnZlIHRoYXQncyBhcyBjbG9zZSB0byAncG9pbnQnIGFzIHBvc3NpYmxlOlxyXG4gICAgICAgIHQgPSBmaW5kX3QoYmV6LCBwYXJhbWV0ZXJzW2ldLCB0X2Rpc3RNYXAsIDEwKTtcclxuXHJcbiAgICAgICAgdiA9IG1hdGhzLnN1YnRyYWN0KGJlemllci5xKGJleiwgdCksIHBvaW50KTtcclxuICAgICAgICBkaXN0ID0gdlswXSp2WzBdICsgdlsxXSp2WzFdO1xyXG5cclxuICAgICAgICBpZiAoZGlzdCA+IG1heERpc3QpIHtcclxuICAgICAgICAgICAgbWF4RGlzdCA9IGRpc3Q7XHJcbiAgICAgICAgICAgIHNwbGl0UG9pbnQgPSBpO1xyXG4gICAgICAgIH1cclxuICAgIH1cclxuXHJcbiAgICByZXR1cm4gW21heERpc3QsIHNwbGl0UG9pbnRdO1xyXG59O1xyXG5cclxuLy9TYW1wbGUgJ3QncyBhbmQgbWFwIHRoZW0gdG8gcmVsYXRpdmUgZGlzdGFuY2VzIGFsb25nIHRoZSBjdXJ2ZTpcclxudmFyIG1hcFR0b1JlbGF0aXZlRGlzdGFuY2VzID0gZnVuY3Rpb24gKGJleiwgQl9wYXJ0cykge1xyXG4gICAgdmFyIEJfdF9jdXJyO1xyXG4gICAgdmFyIEJfdF9kaXN0ID0gWzBdO1xyXG4gICAgdmFyIEJfdF9wcmV2ID0gYmV6WzBdO1xyXG4gICAgdmFyIHN1bUxlbiA9IDA7XHJcblxyXG4gICAgZm9yICh2YXIgaT0xOyBpPD1CX3BhcnRzOyBpKyspIHtcclxuICAgICAgQl90X2N1cnIgPSBiZXppZXIucShiZXosIGkvQl9wYXJ0cyk7XHJcblxyXG4gICAgICBzdW1MZW4gKz0gbWF0aHMudmVjdG9yTGVuKG1hdGhzLnN1YnRyYWN0KEJfdF9jdXJyLCBCX3RfcHJldikpO1xyXG5cclxuICAgICAgQl90X2Rpc3QucHVzaChzdW1MZW4pO1xyXG4gICAgICBCX3RfcHJldiA9IEJfdF9jdXJyO1xyXG4gICAgfVxyXG5cclxuICAgIC8vTm9ybWFsaXplIEJfbGVuZ3RoIHRvIHRoZSBzYW1lIGludGVydmFsIGFzIHRoZSBwYXJhbWV0ZXIgZGlzdGFuY2VzOyAwIHRvIDE6XHJcbiAgICBCX3RfZGlzdCA9IEJfdF9kaXN0Lm1hcCh4ID0+IHgvc3VtTGVuKTtcclxuICAgIHJldHVybiBCX3RfZGlzdDtcclxufTtcclxuXHJcbmZ1bmN0aW9uIGZpbmRfdChiZXosIHBhcmFtLCB0X2Rpc3RNYXAsIEJfcGFydHMpIHtcclxuICAgIGlmKHBhcmFtIDwgMCkgeyByZXR1cm4gMDsgfVxyXG4gICAgaWYocGFyYW0gPiAxKSB7IHJldHVybiAxOyB9XHJcblxyXG4gICAgLypcclxuICAgICAgICAncGFyYW0nIGlzIGEgdmFsdWUgYmV0d2VlbiAwIGFuZCAxIHRlbGxpbmcgdXMgdGhlIHJlbGF0aXZlIHBvc2l0aW9uXHJcbiAgICAgICAgb2YgYSBwb2ludCBvbiB0aGUgc291cmNlIHBvbHlsaW5lIChsaW5lYXJseSBmcm9tIHRoZSBzdGFydCAoMCkgdG8gdGhlIGVuZCAoMSkpLlxyXG4gICAgICAgIFRvIHNlZSBpZiBhIGdpdmVuIGN1cnZlIC0gJ2JleicgLSBpcyBhIGNsb3NlIGFwcHJveGltYXRpb24gb2YgdGhlIHBvbHlsaW5lLFxyXG4gICAgICAgIHdlIGNvbXBhcmUgc3VjaCBhIHBvbHktcG9pbnQgdG8gdGhlIHBvaW50IG9uIHRoZSBjdXJ2ZSB0aGF0J3MgdGhlIHNhbWVcclxuICAgICAgICByZWxhdGl2ZSBkaXN0YW5jZSBhbG9uZyB0aGUgY3VydmUncyBsZW5ndGguXHJcblxyXG4gICAgICAgIEJ1dCBmaW5kaW5nIHRoYXQgY3VydmUtcG9pbnQgdGFrZXMgYSBsaXR0bGUgd29yazpcclxuICAgICAgICBUaGVyZSBpcyBhIGZ1bmN0aW9uIFwiQih0KVwiIHRvIGZpbmQgcG9pbnRzIGFsb25nIGEgY3VydmUgZnJvbSB0aGUgcGFyYW1ldHJpYyBwYXJhbWV0ZXIgJ3QnXHJcbiAgICAgICAgKGFsc28gcmVsYXRpdmUgZnJvbSAwIHRvIDE6IGh0dHA6Ly9zdGFja292ZXJmbG93LmNvbS9hLzMyODQxNzY0LzE4Njk2NjBcclxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgaHR0cDovL3BvbWF4LmdpdGh1Yi5pby9iZXppZXJpbmZvLyNleHBsYW5hdGlvbiksXHJcbiAgICAgICAgYnV0ICd0JyBpc24ndCBsaW5lYXIgYnkgbGVuZ3RoIChodHRwOi8vZ2FtZWRldi5zdGFja2V4Y2hhbmdlLmNvbS9xdWVzdGlvbnMvMTA1MjMwKS5cclxuXHJcbiAgICAgICAgU28sIHdlIHNhbXBsZSBzb21lIHBvaW50cyBhbG9uZyB0aGUgY3VydmUgdXNpbmcgYSBoYW5kZnVsIG9mIHZhbHVlcyBmb3IgJ3QnLlxyXG4gICAgICAgIFRoZW4sIHdlIGNhbGN1bGF0ZSB0aGUgbGVuZ3RoIGJldHdlZW4gdGhvc2Ugc2FtcGxlcyB2aWEgcGxhaW4gZXVjbGlkZWFuIGRpc3RhbmNlO1xyXG4gICAgICAgIEIodCkgY29uY2VudHJhdGVzIHRoZSBwb2ludHMgYXJvdW5kIHNoYXJwIHR1cm5zLCBzbyB0aGlzIHNob3VsZCBnaXZlIHVzIGEgZ29vZC1lbm91Z2ggb3V0bGluZSBvZiB0aGUgY3VydmUuXHJcbiAgICAgICAgVGh1cywgZm9yIGEgZ2l2ZW4gcmVsYXRpdmUgZGlzdGFuY2UgKCdwYXJhbScpLCB3ZSBjYW4gbm93IGZpbmQgYW4gdXBwZXIgYW5kIGxvd2VyIHZhbHVlXHJcbiAgICAgICAgZm9yIHRoZSBjb3JyZXNwb25kaW5nICd0JyBieSBzZWFyY2hpbmcgdGhyb3VnaCB0aG9zZSBzYW1wbGVkIGRpc3RhbmNlcy5cclxuICAgICAgICBGaW5hbGx5LCB3ZSBqdXN0IHVzZSBsaW5lYXIgaW50ZXJwb2xhdGlvbiB0byBmaW5kIGEgYmV0dGVyIHZhbHVlIGZvciB0aGUgZXhhY3QgJ3QnLlxyXG5cclxuICAgICAgICBNb3JlIGluZm86XHJcbiAgICAgICAgICAgIGh0dHA6Ly9nYW1lZGV2LnN0YWNrZXhjaGFuZ2UuY29tL3F1ZXN0aW9ucy8xMDUyMzAvcG9pbnRzLWV2ZW5seS1zcGFjZWQtYWxvbmctYS1iZXppZXItY3VydmVcclxuICAgICAgICAgICAgaHR0cDovL3N0YWNrb3ZlcmZsb3cuY29tL3F1ZXN0aW9ucy8yOTQzODM5OC9jaGVhcC13YXktb2YtY2FsY3VsYXRpbmctY3ViaWMtYmV6aWVyLWxlbmd0aFxyXG4gICAgICAgICAgICBodHRwOi8vc3RldmUuaG9sbGFzY2gubmV0L2NnaW5kZXgvY3VydmVzL2NiZXphcmNsZW4uaHRtbFxyXG4gICAgICAgICAgICBodHRwczovL2dpdGh1Yi5jb20vcmV0dXh4L3RpbnlzcGxpbmVcclxuICAgICovXHJcbiAgICB2YXIgbGVuTWF4LCBsZW5NaW4sIHRNYXgsIHRNaW4sIHQ7XHJcblxyXG4gICAgLy9GaW5kIHRoZSB0d28gdC1zIHRoYXQgdGhlIGN1cnJlbnQgcGFyYW0gZGlzdGFuY2UgbGllcyBiZXR3ZWVuLFxyXG4gICAgLy9hbmQgdGhlbiBpbnRlcnBvbGF0ZSBhIHNvbWV3aGF0IGFjY3VyYXRlIHZhbHVlIGZvciB0aGUgZXhhY3QgdDpcclxuICAgIGZvcih2YXIgaSA9IDE7IGkgPD0gQl9wYXJ0czsgaSsrKSB7XHJcblxyXG4gICAgICAgIGlmKHBhcmFtIDw9IHRfZGlzdE1hcFtpXSkge1xyXG4gICAgICAgICAgICB0TWluICAgPSAoaS0xKSAvIEJfcGFydHM7XHJcbiAgICAgICAgICAgIHRNYXggICA9IGkgLyBCX3BhcnRzO1xyXG4gICAgICAgICAgICBsZW5NaW4gPSB0X2Rpc3RNYXBbaS0xXTtcclxuICAgICAgICAgICAgbGVuTWF4ID0gdF9kaXN0TWFwW2ldO1xyXG5cclxuICAgICAgICAgICAgdCA9IChwYXJhbS1sZW5NaW4pLyhsZW5NYXgtbGVuTWluKSAqICh0TWF4LXRNaW4pICsgdE1pbjtcclxuICAgICAgICAgICAgYnJlYWs7XHJcbiAgICAgICAgfVxyXG4gICAgfVxyXG4gICAgcmV0dXJuIHQ7XHJcbn1cclxuXHJcbi8qKlxyXG4gKiBDcmVhdGVzIGEgdmVjdG9yIG9mIGxlbmd0aCAxIHdoaWNoIHNob3dzIHRoZSBkaXJlY3Rpb24gZnJvbSBCIHRvIEFcclxuICovXHJcbmZ1bmN0aW9uIGNyZWF0ZVRhbmdlbnQocG9pbnRBLCBwb2ludEIpIHtcclxuICAgIHJldHVybiBtYXRocy5ub3JtYWxpemUobWF0aHMuc3VidHJhY3QocG9pbnRBLCBwb2ludEIpKTtcclxufVxyXG5cclxuLypcclxuICAgIFNpbXBsaWZpZWQgdmVyc2lvbnMgb2Ygd2hhdCB3ZSBuZWVkIGZyb20gbWF0aC5qc1xyXG4gICAgT3B0aW1pemVkIGZvciBvdXIgaW5wdXQsIHdoaWNoIGlzIG9ubHkgbnVtYmVycyBhbmQgMXgyIGFycmF5cyAoaS5lLiBbeCwgeV0gY29vcmRpbmF0ZXMpLlxyXG4qL1xyXG5jbGFzcyBtYXRocyB7XHJcbiAgICAvL3plcm9zID0gbG9nQW5kUnVuKG1hdGguemVyb3MpO1xyXG4gICAgc3RhdGljIHplcm9zX1h4MngyKHgpIHtcclxuICAgICAgICB2YXIgenMgPSBbXTtcclxuICAgICAgICB3aGlsZSh4LS0pIHsgenMucHVzaChbMCwwXSk7IH1cclxuICAgICAgICByZXR1cm4genNcclxuICAgIH1cclxuXHJcbiAgICAvL211bHRpcGx5ID0gbG9nQW5kUnVuKG1hdGgubXVsdGlwbHkpO1xyXG4gICAgc3RhdGljIG11bEl0ZW1zKGl0ZW1zLCBtdWx0aXBsaWVyKSB7XHJcbiAgICAgICAgLy9yZXR1cm4gaXRlbXMubWFwKHggPT4geCptdWx0aXBsaWVyKTtcclxuICAgICAgICByZXR1cm4gW2l0ZW1zWzBdKm11bHRpcGxpZXIsIGl0ZW1zWzFdKm11bHRpcGxpZXJdO1xyXG4gICAgfVxyXG4gICAgc3RhdGljIG11bE1hdHJpeChtMSwgbTIpIHtcclxuICAgICAgICAvL2h0dHBzOi8vZW4ud2lraXBlZGlhLm9yZy93aWtpL01hdHJpeF9tdWx0aXBsaWNhdGlvbiNNYXRyaXhfcHJvZHVjdF8uMjh0d29fbWF0cmljZXMuMjlcclxuICAgICAgICAvL1NpbXBsaWZpZWQgdG8gb25seSBoYW5kbGUgMS1kaW1lbnNpb25hbCBtYXRyaWNlcyAoaS5lLiBhcnJheXMpIG9mIGVxdWFsIGxlbmd0aDpcclxuICAgICAgICAvLyAgcmV0dXJuIG0xLnJlZHVjZSgoc3VtLHgxLGkpID0+IHN1bSArICh4MSptMltpXSksXHJcbiAgICAgICAgLy8gICAgICAgICAgICAgICAgICAgMCk7XHJcbiAgICAgICAgcmV0dXJuIChtMVswXSptMlswXSkgKyAobTFbMV0qbTJbMV0pO1xyXG4gICAgfVxyXG5cclxuICAgIC8vT25seSB1c2VkIHRvIHN1YnJhY3QgdG8gcG9pbnRzIChvciBhdCBsZWFzdCBhcnJheXMpOlxyXG4gICAgLy8gIHN1YnRyYWN0ID0gbG9nQW5kUnVuKG1hdGguc3VidHJhY3QpO1xyXG4gICAgc3RhdGljIHN1YnRyYWN0KGFycjEsIGFycjIpIHtcclxuICAgICAgICAvL3JldHVybiBhcnIxLm1hcCgoeDEsIGkpID0+IHgxIC0gYXJyMltpXSk7XHJcbiAgICAgICAgcmV0dXJuIFthcnIxWzBdLWFycjJbMF0sIGFycjFbMV0tYXJyMlsxXV07XHJcbiAgICB9XHJcblxyXG4gICAgLy9hZGQgPSBsb2dBbmRSdW4obWF0aC5hZGQpO1xyXG4gICAgc3RhdGljIGFkZEFycmF5cyhhcnIxLCBhcnIyKSB7XHJcbiAgICAgICAgLy9yZXR1cm4gYXJyMS5tYXAoKHgxLCBpKSA9PiB4MSArIGFycjJbaV0pO1xyXG4gICAgICAgIHJldHVybiBbYXJyMVswXSthcnIyWzBdLCBhcnIxWzFdK2FycjJbMV1dO1xyXG4gICAgfVxyXG4gICAgc3RhdGljIGFkZEl0ZW1zKGl0ZW1zLCBhZGRpdGlvbikge1xyXG4gICAgICAgIC8vcmV0dXJuIGl0ZW1zLm1hcCh4ID0+IHgrYWRkaXRpb24pO1xyXG4gICAgICAgIHJldHVybiBbaXRlbXNbMF0rYWRkaXRpb24sIGl0ZW1zWzFdK2FkZGl0aW9uXTtcclxuICAgIH1cclxuXHJcbiAgICAvL3ZhciBzdW0gPSBsb2dBbmRSdW4obWF0aC5zdW0pO1xyXG4gICAgc3RhdGljIHN1bShpdGVtcykge1xyXG4gICAgICAgIHJldHVybiBpdGVtcy5yZWR1Y2UoKHN1bSx4KSA9PiBzdW0gKyB4KTtcclxuICAgIH1cclxuXHJcbiAgICAvL2NoYWluID0gbWF0aC5jaGFpbjtcclxuXHJcbiAgICAvL09ubHkgdXNlZCBvbiB0d28gYXJyYXlzLiBUaGUgZG90IHByb2R1Y3QgaXMgZXF1YWwgdG8gdGhlIG1hdHJpeCBwcm9kdWN0IGluIHRoaXMgY2FzZTpcclxuICAgIC8vICBkb3QgPSBsb2dBbmRSdW4obWF0aC5kb3QpO1xyXG4gICAgc3RhdGljIGRvdChtMSwgbTIpIHtcclxuICAgICAgICByZXR1cm4gbWF0aHMubXVsTWF0cml4KG0xLCBtMik7XHJcbiAgICB9XHJcblxyXG4gICAgLy9odHRwczovL2VuLndpa2lwZWRpYS5vcmcvd2lraS9Ob3JtXyhtYXRoZW1hdGljcykjRXVjbGlkZWFuX25vcm1cclxuICAgIC8vICB2YXIgbm9ybSA9IGxvZ0FuZFJ1bihtYXRoLm5vcm0pO1xyXG4gICAgc3RhdGljIHZlY3Rvckxlbih2KSB7XHJcbiAgICAgICAgdmFyIGEgPSB2WzBdLCBiID0gdlsxXTtcclxuICAgICAgICByZXR1cm4gTWF0aC5zcXJ0KGEqYSArIGIqYik7XHJcbiAgICB9XHJcblxyXG4gICAgLy9tYXRoLmRpdmlkZSA9IGxvZ0FuZFJ1bihtYXRoLmRpdmlkZSk7XHJcbiAgICBzdGF0aWMgZGl2SXRlbXMoaXRlbXMsIGRpdmlzb3IpIHtcclxuICAgICAgICAvL3JldHVybiBpdGVtcy5tYXAoeCA9PiB4L2Rpdmlzb3IpO1xyXG4gICAgICAgIHJldHVybiBbaXRlbXNbMF0vZGl2aXNvciwgaXRlbXNbMV0vZGl2aXNvcl07XHJcbiAgICB9XHJcblxyXG4gICAgLy92YXIgZG90UG93ID0gbG9nQW5kUnVuKG1hdGguZG90UG93KTtcclxuICAgIHN0YXRpYyBzcXVhcmVJdGVtcyhpdGVtcykge1xyXG4gICAgICAgIC8vcmV0dXJuIGl0ZW1zLm1hcCh4ID0+IHgqeCk7XHJcbiAgICAgICAgdmFyIGEgPSBpdGVtc1swXSwgYiA9IGl0ZW1zWzFdO1xyXG4gICAgICAgIHJldHVybiBbYSphLCBiKmJdO1xyXG4gICAgfVxyXG5cclxuICAgIHN0YXRpYyBub3JtYWxpemUodikge1xyXG4gICAgICAgIHJldHVybiB0aGlzLmRpdkl0ZW1zKHYsIHRoaXMudmVjdG9yTGVuKHYpKTtcclxuICAgIH1cclxuXHJcbiAgICAvL01hdGgucG93ID0gbG9nQW5kUnVuKE1hdGgucG93KTtcclxufVxyXG5cclxuXHJcbmNsYXNzIGJlemllciB7XHJcbiAgICAvL0V2YWx1YXRlcyBjdWJpYyBiZXppZXIgYXQgdCwgcmV0dXJuIHBvaW50XHJcbiAgICBzdGF0aWMgcShjdHJsUG9seSwgdCkge1xyXG4gICAgICAgIHZhciB0eCA9IDEuMCAtIHQ7XHJcbiAgICAgICAgdmFyIHBBID0gbWF0aHMubXVsSXRlbXMoIGN0cmxQb2x5WzBdLCAgICAgIHR4ICogdHggKiB0eCApLFxyXG4gICAgICAgICAgICBwQiA9IG1hdGhzLm11bEl0ZW1zKCBjdHJsUG9seVsxXSwgIDMgKiB0eCAqIHR4ICogIHQgKSxcclxuICAgICAgICAgICAgcEMgPSBtYXRocy5tdWxJdGVtcyggY3RybFBvbHlbMl0sICAzICogdHggKiAgdCAqICB0ICksXHJcbiAgICAgICAgICAgIHBEID0gbWF0aHMubXVsSXRlbXMoIGN0cmxQb2x5WzNdLCAgICAgICB0ICogIHQgKiAgdCApO1xyXG4gICAgICAgIHJldHVybiBtYXRocy5hZGRBcnJheXMobWF0aHMuYWRkQXJyYXlzKHBBLCBwQiksIG1hdGhzLmFkZEFycmF5cyhwQywgcEQpKTtcclxuICAgIH1cclxuXHJcbiAgICAvL0V2YWx1YXRlcyBjdWJpYyBiZXppZXIgZmlyc3QgZGVyaXZhdGl2ZSBhdCB0LCByZXR1cm4gcG9pbnRcclxuICAgIHN0YXRpYyBxcHJpbWUoY3RybFBvbHksIHQpIHtcclxuICAgICAgICB2YXIgdHggPSAxLjAgLSB0O1xyXG4gICAgICAgIHZhciBwQSA9IG1hdGhzLm11bEl0ZW1zKCBtYXRocy5zdWJ0cmFjdChjdHJsUG9seVsxXSwgY3RybFBvbHlbMF0pLCAgMyAqIHR4ICogdHggKSxcclxuICAgICAgICAgICAgcEIgPSBtYXRocy5tdWxJdGVtcyggbWF0aHMuc3VidHJhY3QoY3RybFBvbHlbMl0sIGN0cmxQb2x5WzFdKSwgIDYgKiB0eCAqICB0ICksXHJcbiAgICAgICAgICAgIHBDID0gbWF0aHMubXVsSXRlbXMoIG1hdGhzLnN1YnRyYWN0KGN0cmxQb2x5WzNdLCBjdHJsUG9seVsyXSksICAzICogIHQgKiAgdCApO1xyXG4gICAgICAgIHJldHVybiBtYXRocy5hZGRBcnJheXMobWF0aHMuYWRkQXJyYXlzKHBBLCBwQiksIHBDKTtcclxuICAgIH1cclxuXHJcbiAgICAvL0V2YWx1YXRlcyBjdWJpYyBiZXppZXIgc2Vjb25kIGRlcml2YXRpdmUgYXQgdCwgcmV0dXJuIHBvaW50XHJcbiAgICBzdGF0aWMgcXByaW1lcHJpbWUoY3RybFBvbHksIHQpIHtcclxuICAgICAgICByZXR1cm4gbWF0aHMuYWRkQXJyYXlzKG1hdGhzLm11bEl0ZW1zKCBtYXRocy5hZGRBcnJheXMobWF0aHMuc3VidHJhY3QoY3RybFBvbHlbMl0sIG1hdGhzLm11bEl0ZW1zKGN0cmxQb2x5WzFdLCAyKSksIGN0cmxQb2x5WzBdKSwgIDYgKiAoMS4wIC0gdCkgKSxcclxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIG1hdGhzLm11bEl0ZW1zKCBtYXRocy5hZGRBcnJheXMobWF0aHMuc3VidHJhY3QoY3RybFBvbHlbM10sIG1hdGhzLm11bEl0ZW1zKGN0cmxQb2x5WzJdLCAyKSksIGN0cmxQb2x5WzFdKSwgIDYgKiAgICAgICAgdCAgKSk7XHJcbiAgICB9XHJcbn1cclxuXHJcbm1vZHVsZS5leHBvcnRzID0gZml0Q3VydmU7XHJcbiIsIi8qIVxuKiBzdmcuanMgLSBBIGxpZ2h0d2VpZ2h0IGxpYnJhcnkgZm9yIG1hbmlwdWxhdGluZyBhbmQgYW5pbWF0aW5nIFNWRy5cbiogQHZlcnNpb24gMi4zLjdcbiogaHR0cHM6Ly9zdmdkb3Rqcy5naXRodWIuaW8vXG4qXG4qIEBjb3B5cmlnaHQgV291dCBGaWVyZW5zIDx3b3V0QG1pY2std291dC5jb20+XG4qIEBsaWNlbnNlIE1JVFxuKlxuKiBCVUlMVDogU2F0IEphbiAxNCAyMDE3IDA3OjIzOjE4IEdNVCswMTAwIChDRVQpXG4qLztcbihmdW5jdGlvbihyb290LCBmYWN0b3J5KSB7XG4gIGlmICh0eXBlb2YgZGVmaW5lID09PSAnZnVuY3Rpb24nICYmIGRlZmluZS5hbWQpIHtcbiAgICBkZWZpbmUoZnVuY3Rpb24oKXtcbiAgICAgIHJldHVybiBmYWN0b3J5KHJvb3QsIHJvb3QuZG9jdW1lbnQpXG4gICAgfSlcbiAgfSBlbHNlIGlmICh0eXBlb2YgZXhwb3J0cyA9PT0gJ29iamVjdCcpIHtcbiAgICBtb2R1bGUuZXhwb3J0cyA9IHJvb3QuZG9jdW1lbnQgPyBmYWN0b3J5KHJvb3QsIHJvb3QuZG9jdW1lbnQpIDogZnVuY3Rpb24odyl7IHJldHVybiBmYWN0b3J5KHcsIHcuZG9jdW1lbnQpIH1cbiAgfSBlbHNlIHtcbiAgICByb290LlNWRyA9IGZhY3Rvcnkocm9vdCwgcm9vdC5kb2N1bWVudClcbiAgfVxufSh0eXBlb2Ygd2luZG93ICE9PSBcInVuZGVmaW5lZFwiID8gd2luZG93IDogdGhpcywgZnVuY3Rpb24od2luZG93LCBkb2N1bWVudCkge1xuXG4vLyBUaGUgbWFpbiB3cmFwcGluZyBlbGVtZW50XG52YXIgU1ZHID0gdGhpcy5TVkcgPSBmdW5jdGlvbihlbGVtZW50KSB7XG4gIGlmIChTVkcuc3VwcG9ydGVkKSB7XG4gICAgZWxlbWVudCA9IG5ldyBTVkcuRG9jKGVsZW1lbnQpXG5cbiAgICBpZighU1ZHLnBhcnNlci5kcmF3KVxuICAgICAgU1ZHLnByZXBhcmUoKVxuXG4gICAgcmV0dXJuIGVsZW1lbnRcbiAgfVxufVxuXG4vLyBEZWZhdWx0IG5hbWVzcGFjZXNcblNWRy5ucyAgICA9ICdodHRwOi8vd3d3LnczLm9yZy8yMDAwL3N2ZydcblNWRy54bWxucyA9ICdodHRwOi8vd3d3LnczLm9yZy8yMDAwL3htbG5zLydcblNWRy54bGluayA9ICdodHRwOi8vd3d3LnczLm9yZy8xOTk5L3hsaW5rJ1xuU1ZHLnN2Z2pzID0gJ2h0dHA6Ly9zdmdqcy5jb20vc3ZnanMnXG5cbi8vIFN2ZyBzdXBwb3J0IHRlc3RcblNWRy5zdXBwb3J0ZWQgPSAoZnVuY3Rpb24oKSB7XG4gIHJldHVybiAhISBkb2N1bWVudC5jcmVhdGVFbGVtZW50TlMgJiZcbiAgICAgICAgICEhIGRvY3VtZW50LmNyZWF0ZUVsZW1lbnROUyhTVkcubnMsJ3N2ZycpLmNyZWF0ZVNWR1JlY3Rcbn0pKClcblxuLy8gRG9uJ3QgYm90aGVyIHRvIGNvbnRpbnVlIGlmIFNWRyBpcyBub3Qgc3VwcG9ydGVkXG5pZiAoIVNWRy5zdXBwb3J0ZWQpIHJldHVybiBmYWxzZVxuXG4vLyBFbGVtZW50IGlkIHNlcXVlbmNlXG5TVkcuZGlkICA9IDEwMDBcblxuLy8gR2V0IG5leHQgbmFtZWQgZWxlbWVudCBpZFxuU1ZHLmVpZCA9IGZ1bmN0aW9uKG5hbWUpIHtcbiAgcmV0dXJuICdTdmdqcycgKyBjYXBpdGFsaXplKG5hbWUpICsgKFNWRy5kaWQrKylcbn1cblxuLy8gTWV0aG9kIGZvciBlbGVtZW50IGNyZWF0aW9uXG5TVkcuY3JlYXRlID0gZnVuY3Rpb24obmFtZSkge1xuICAvLyBjcmVhdGUgZWxlbWVudFxuICB2YXIgZWxlbWVudCA9IGRvY3VtZW50LmNyZWF0ZUVsZW1lbnROUyh0aGlzLm5zLCBuYW1lKVxuXG4gIC8vIGFwcGx5IHVuaXF1ZSBpZFxuICBlbGVtZW50LnNldEF0dHJpYnV0ZSgnaWQnLCB0aGlzLmVpZChuYW1lKSlcblxuICByZXR1cm4gZWxlbWVudFxufVxuXG4vLyBNZXRob2QgZm9yIGV4dGVuZGluZyBvYmplY3RzXG5TVkcuZXh0ZW5kID0gZnVuY3Rpb24oKSB7XG4gIHZhciBtb2R1bGVzLCBtZXRob2RzLCBrZXksIGlcblxuICAvLyBHZXQgbGlzdCBvZiBtb2R1bGVzXG4gIG1vZHVsZXMgPSBbXS5zbGljZS5jYWxsKGFyZ3VtZW50cylcblxuICAvLyBHZXQgb2JqZWN0IHdpdGggZXh0ZW5zaW9uc1xuICBtZXRob2RzID0gbW9kdWxlcy5wb3AoKVxuXG4gIGZvciAoaSA9IG1vZHVsZXMubGVuZ3RoIC0gMTsgaSA+PSAwOyBpLS0pXG4gICAgaWYgKG1vZHVsZXNbaV0pXG4gICAgICBmb3IgKGtleSBpbiBtZXRob2RzKVxuICAgICAgICBtb2R1bGVzW2ldLnByb3RvdHlwZVtrZXldID0gbWV0aG9kc1trZXldXG5cbiAgLy8gTWFrZSBzdXJlIFNWRy5TZXQgaW5oZXJpdHMgYW55IG5ld2x5IGFkZGVkIG1ldGhvZHNcbiAgaWYgKFNWRy5TZXQgJiYgU1ZHLlNldC5pbmhlcml0KVxuICAgIFNWRy5TZXQuaW5oZXJpdCgpXG59XG5cbi8vIEludmVudCBuZXcgZWxlbWVudFxuU1ZHLmludmVudCA9IGZ1bmN0aW9uKGNvbmZpZykge1xuICAvLyBDcmVhdGUgZWxlbWVudCBpbml0aWFsaXplclxuICB2YXIgaW5pdGlhbGl6ZXIgPSB0eXBlb2YgY29uZmlnLmNyZWF0ZSA9PSAnZnVuY3Rpb24nID9cbiAgICBjb25maWcuY3JlYXRlIDpcbiAgICBmdW5jdGlvbigpIHtcbiAgICAgIHRoaXMuY29uc3RydWN0b3IuY2FsbCh0aGlzLCBTVkcuY3JlYXRlKGNvbmZpZy5jcmVhdGUpKVxuICAgIH1cblxuICAvLyBJbmhlcml0IHByb3RvdHlwZVxuICBpZiAoY29uZmlnLmluaGVyaXQpXG4gICAgaW5pdGlhbGl6ZXIucHJvdG90eXBlID0gbmV3IGNvbmZpZy5pbmhlcml0XG5cbiAgLy8gRXh0ZW5kIHdpdGggbWV0aG9kc1xuICBpZiAoY29uZmlnLmV4dGVuZClcbiAgICBTVkcuZXh0ZW5kKGluaXRpYWxpemVyLCBjb25maWcuZXh0ZW5kKVxuXG4gIC8vIEF0dGFjaCBjb25zdHJ1Y3QgbWV0aG9kIHRvIHBhcmVudFxuICBpZiAoY29uZmlnLmNvbnN0cnVjdClcbiAgICBTVkcuZXh0ZW5kKGNvbmZpZy5wYXJlbnQgfHwgU1ZHLkNvbnRhaW5lciwgY29uZmlnLmNvbnN0cnVjdClcblxuICByZXR1cm4gaW5pdGlhbGl6ZXJcbn1cblxuLy8gQWRvcHQgZXhpc3Rpbmcgc3ZnIGVsZW1lbnRzXG5TVkcuYWRvcHQgPSBmdW5jdGlvbihub2RlKSB7XG4gIC8vIGNoZWNrIGZvciBwcmVzZW5jZSBvZiBub2RlXG4gIGlmICghbm9kZSkgcmV0dXJuIG51bGxcblxuICAvLyBtYWtlIHN1cmUgYSBub2RlIGlzbid0IGFscmVhZHkgYWRvcHRlZFxuICBpZiAobm9kZS5pbnN0YW5jZSkgcmV0dXJuIG5vZGUuaW5zdGFuY2VcblxuICAvLyBpbml0aWFsaXplIHZhcmlhYmxlc1xuICB2YXIgZWxlbWVudFxuXG4gIC8vIGFkb3B0IHdpdGggZWxlbWVudC1zcGVjaWZpYyBzZXR0aW5nc1xuICBpZiAobm9kZS5ub2RlTmFtZSA9PSAnc3ZnJylcbiAgICBlbGVtZW50ID0gbm9kZS5wYXJlbnROb2RlIGluc3RhbmNlb2YgU1ZHRWxlbWVudCA/IG5ldyBTVkcuTmVzdGVkIDogbmV3IFNWRy5Eb2NcbiAgZWxzZSBpZiAobm9kZS5ub2RlTmFtZSA9PSAnbGluZWFyR3JhZGllbnQnKVxuICAgIGVsZW1lbnQgPSBuZXcgU1ZHLkdyYWRpZW50KCdsaW5lYXInKVxuICBlbHNlIGlmIChub2RlLm5vZGVOYW1lID09ICdyYWRpYWxHcmFkaWVudCcpXG4gICAgZWxlbWVudCA9IG5ldyBTVkcuR3JhZGllbnQoJ3JhZGlhbCcpXG4gIGVsc2UgaWYgKFNWR1tjYXBpdGFsaXplKG5vZGUubm9kZU5hbWUpXSlcbiAgICBlbGVtZW50ID0gbmV3IFNWR1tjYXBpdGFsaXplKG5vZGUubm9kZU5hbWUpXVxuICBlbHNlXG4gICAgZWxlbWVudCA9IG5ldyBTVkcuRWxlbWVudChub2RlKVxuXG4gIC8vIGVuc3VyZSByZWZlcmVuY2VzXG4gIGVsZW1lbnQudHlwZSAgPSBub2RlLm5vZGVOYW1lXG4gIGVsZW1lbnQubm9kZSAgPSBub2RlXG4gIG5vZGUuaW5zdGFuY2UgPSBlbGVtZW50XG5cbiAgLy8gU1ZHLkNsYXNzIHNwZWNpZmljIHByZXBhcmF0aW9uc1xuICBpZiAoZWxlbWVudCBpbnN0YW5jZW9mIFNWRy5Eb2MpXG4gICAgZWxlbWVudC5uYW1lc3BhY2UoKS5kZWZzKClcblxuICAvLyBwdWxsIHN2Z2pzIGRhdGEgZnJvbSB0aGUgZG9tIChnZXRBdHRyaWJ1dGVOUyBkb2Vzbid0IHdvcmsgaW4gaHRtbDUpXG4gIGVsZW1lbnQuc2V0RGF0YShKU09OLnBhcnNlKG5vZGUuZ2V0QXR0cmlidXRlKCdzdmdqczpkYXRhJykpIHx8IHt9KVxuXG4gIHJldHVybiBlbGVtZW50XG59XG5cbi8vIEluaXRpYWxpemUgcGFyc2luZyBlbGVtZW50XG5TVkcucHJlcGFyZSA9IGZ1bmN0aW9uKCkge1xuICAvLyBTZWxlY3QgZG9jdW1lbnQgYm9keSBhbmQgY3JlYXRlIGludmlzaWJsZSBzdmcgZWxlbWVudFxuICB2YXIgYm9keSA9IGRvY3VtZW50LmdldEVsZW1lbnRzQnlUYWdOYW1lKCdib2R5JylbMF1cbiAgICAsIGRyYXcgPSAoYm9keSA/IG5ldyBTVkcuRG9jKGJvZHkpIDogIG5ldyBTVkcuRG9jKGRvY3VtZW50LmRvY3VtZW50RWxlbWVudCkubmVzdGVkKCkpLnNpemUoMiwgMClcblxuICAvLyBDcmVhdGUgcGFyc2VyIG9iamVjdFxuICBTVkcucGFyc2VyID0ge1xuICAgIGJvZHk6IGJvZHkgfHwgZG9jdW1lbnQuZG9jdW1lbnRFbGVtZW50XG4gICwgZHJhdzogZHJhdy5zdHlsZSgnb3BhY2l0eTowO3Bvc2l0aW9uOmZpeGVkO2xlZnQ6MTAwJTt0b3A6MTAwJTtvdmVyZmxvdzpoaWRkZW4nKVxuICAsIHBvbHk6IGRyYXcucG9seWxpbmUoKS5ub2RlXG4gICwgcGF0aDogZHJhdy5wYXRoKCkubm9kZVxuICAsIG5hdGl2ZTogU1ZHLmNyZWF0ZSgnc3ZnJylcbiAgfVxufVxuXG5TVkcucGFyc2VyID0ge1xuICBuYXRpdmU6IFNWRy5jcmVhdGUoJ3N2ZycpXG59XG5cbmRvY3VtZW50LmFkZEV2ZW50TGlzdGVuZXIoJ0RPTUNvbnRlbnRMb2FkZWQnLCBmdW5jdGlvbigpIHtcbiAgaWYoIVNWRy5wYXJzZXIuZHJhdylcbiAgICBTVkcucHJlcGFyZSgpXG59LCBmYWxzZSlcblxuLy8gU3RvcmFnZSBmb3IgcmVndWxhciBleHByZXNzaW9uc1xuU1ZHLnJlZ2V4ID0ge1xuICAvLyBQYXJzZSB1bml0IHZhbHVlXG4gIG51bWJlckFuZFVuaXQ6ICAgIC9eKFsrLV0/KFxcZCsoXFwuXFxkKik/fFxcLlxcZCspKGVbKy1dP1xcZCspPykoW2EteiVdKikkL2lcblxuICAvLyBQYXJzZSBoZXggdmFsdWVcbiwgaGV4OiAgICAgICAgICAgICAgL14jPyhbYS1mXFxkXXsyfSkoW2EtZlxcZF17Mn0pKFthLWZcXGRdezJ9KSQvaVxuXG4gIC8vIFBhcnNlIHJnYiB2YWx1ZVxuLCByZ2I6ICAgICAgICAgICAgICAvcmdiXFwoKFxcZCspLChcXGQrKSwoXFxkKylcXCkvXG5cbiAgLy8gUGFyc2UgcmVmZXJlbmNlIGlkXG4sIHJlZmVyZW5jZTogICAgICAgIC8jKFthLXowLTlcXC1fXSspL2lcblxuICAvLyBQYXJzZSBtYXRyaXggd3JhcHBlclxuLCBtYXRyaXg6ICAgICAgICAgICAvbWF0cml4XFwofFxcKS9nXG5cbiAgLy8gRWxlbWVudHMgb2YgYSBtYXRyaXhcbiwgbWF0cml4RWxlbWVudHM6ICAgLywqXFxzK3wsL1xuXG4gIC8vIFdoaXRlc3BhY2Vcbiwgd2hpdGVzcGFjZTogICAgICAgL1xccy9nXG5cbiAgLy8gVGVzdCBoZXggdmFsdWVcbiwgaXNIZXg6ICAgICAgICAgICAgL14jW2EtZjAtOV17Myw2fSQvaVxuXG4gIC8vIFRlc3QgcmdiIHZhbHVlXG4sIGlzUmdiOiAgICAgICAgICAgIC9ecmdiXFwoL1xuXG4gIC8vIFRlc3QgY3NzIGRlY2xhcmF0aW9uXG4sIGlzQ3NzOiAgICAgICAgICAgIC9bXjpdKzpbXjtdKzs/L1xuXG4gIC8vIFRlc3QgZm9yIGJsYW5rIHN0cmluZ1xuLCBpc0JsYW5rOiAgICAgICAgICAvXihcXHMrKT8kL1xuXG4gIC8vIFRlc3QgZm9yIG51bWVyaWMgc3RyaW5nXG4sIGlzTnVtYmVyOiAgICAgICAgIC9eWystXT8oXFxkKyhcXC5cXGQqKT98XFwuXFxkKykoZVsrLV0/XFxkKyk/JC9pXG5cbiAgLy8gVGVzdCBmb3IgcGVyY2VudCB2YWx1ZVxuLCBpc1BlcmNlbnQ6ICAgICAgICAvXi0/W1xcZFxcLl0rJSQvXG5cbiAgLy8gVGVzdCBmb3IgaW1hZ2UgdXJsXG4sIGlzSW1hZ2U6ICAgICAgICAgIC9cXC4oanBnfGpwZWd8cG5nfGdpZnxzdmcpKFxcP1tePV0rLiopPy9pXG5cbiAgLy8gVGhlIGZvbGxvd2luZyByZWdleCBhcmUgdXNlZCB0byBwYXJzZSB0aGUgZCBhdHRyaWJ1dGUgb2YgYSBwYXRoXG5cbiAgLy8gUmVwbGFjZXMgYWxsIG5lZ2F0aXZlIGV4cG9uZW50c1xuLCBuZWdFeHA6ICAgICAgICAgICAvZVxcLS9naVxuXG4gIC8vIFJlcGxhY2VzIGFsbCBjb21tYVxuLCBjb21tYTogICAgICAgICAgICAvLC9nXG5cbiAgLy8gUmVwbGFjZXMgYWxsIGh5cGhlbnNcbiwgaHlwaGVuOiAgICAgICAgICAgL1xcLS9nXG5cbiAgLy8gUmVwbGFjZXMgYW5kIHRlc3RzIGZvciBhbGwgcGF0aCBsZXR0ZXJzXG4sIHBhdGhMZXR0ZXJzOiAgICAgIC9bTUxIVkNTUVRBWl0vZ2lcblxuICAvLyB5ZXMgd2UgbmVlZCB0aGlzIG9uZSwgdG9vXG4sIGlzUGF0aExldHRlcjogICAgIC9bTUxIVkNTUVRBWl0vaVxuXG4gIC8vIHNwbGl0IGF0IHdoaXRlc3BhY2VzXG4sIHdoaXRlc3BhY2VzOiAgICAgIC9cXHMrL1xuXG4gIC8vIG1hdGNoZXMgWFxuLCBYOiAgICAgICAgICAgICAgICAvWC9nXG59XG5cblNWRy51dGlscyA9IHtcbiAgLy8gTWFwIGZ1bmN0aW9uXG4gIG1hcDogZnVuY3Rpb24oYXJyYXksIGJsb2NrKSB7XG4gICAgdmFyIGlcbiAgICAgICwgaWwgPSBhcnJheS5sZW5ndGhcbiAgICAgICwgcmVzdWx0ID0gW11cblxuICAgIGZvciAoaSA9IDA7IGkgPCBpbDsgaSsrKVxuICAgICAgcmVzdWx0LnB1c2goYmxvY2soYXJyYXlbaV0pKVxuXG4gICAgcmV0dXJuIHJlc3VsdFxuICB9XG5cbiAgLy8gRmlsdGVyIGZ1bmN0aW9uXG4sIGZpbHRlcjogZnVuY3Rpb24oYXJyYXksIGJsb2NrKSB7XG4gICAgdmFyIGlcbiAgICAgICwgaWwgPSBhcnJheS5sZW5ndGhcbiAgICAgICwgcmVzdWx0ID0gW11cblxuICAgIGZvciAoaSA9IDA7IGkgPCBpbDsgaSsrKVxuICAgICAgaWYgKGJsb2NrKGFycmF5W2ldKSlcbiAgICAgICAgcmVzdWx0LnB1c2goYXJyYXlbaV0pXG5cbiAgICByZXR1cm4gcmVzdWx0XG4gIH1cblxuICAvLyBEZWdyZWVzIHRvIHJhZGlhbnNcbiwgcmFkaWFuczogZnVuY3Rpb24oZCkge1xuICAgIHJldHVybiBkICUgMzYwICogTWF0aC5QSSAvIDE4MFxuICB9XG5cbiAgLy8gUmFkaWFucyB0byBkZWdyZWVzXG4sIGRlZ3JlZXM6IGZ1bmN0aW9uKHIpIHtcbiAgICByZXR1cm4gciAqIDE4MCAvIE1hdGguUEkgJSAzNjBcbiAgfVxuXG4sIGZpbHRlclNWR0VsZW1lbnRzOiBmdW5jdGlvbihub2Rlcykge1xuICAgIHJldHVybiB0aGlzLmZpbHRlciggbm9kZXMsIGZ1bmN0aW9uKGVsKSB7IHJldHVybiBlbCBpbnN0YW5jZW9mIFNWR0VsZW1lbnQgfSlcbiAgfVxuXG59XG5cblNWRy5kZWZhdWx0cyA9IHtcbiAgLy8gRGVmYXVsdCBhdHRyaWJ1dGUgdmFsdWVzXG4gIGF0dHJzOiB7XG4gICAgLy8gZmlsbCBhbmQgc3Ryb2tlXG4gICAgJ2ZpbGwtb3BhY2l0eSc6ICAgICAxXG4gICwgJ3N0cm9rZS1vcGFjaXR5JzogICAxXG4gICwgJ3N0cm9rZS13aWR0aCc6ICAgICAwXG4gICwgJ3N0cm9rZS1saW5lam9pbic6ICAnbWl0ZXInXG4gICwgJ3N0cm9rZS1saW5lY2FwJzogICAnYnV0dCdcbiAgLCBmaWxsOiAgICAgICAgICAgICAgICcjMDAwMDAwJ1xuICAsIHN0cm9rZTogICAgICAgICAgICAgJyMwMDAwMDAnXG4gICwgb3BhY2l0eTogICAgICAgICAgICAxXG4gICAgLy8gcG9zaXRpb25cbiAgLCB4OiAgICAgICAgICAgICAgICAgIDBcbiAgLCB5OiAgICAgICAgICAgICAgICAgIDBcbiAgLCBjeDogICAgICAgICAgICAgICAgIDBcbiAgLCBjeTogICAgICAgICAgICAgICAgIDBcbiAgICAvLyBzaXplXG4gICwgd2lkdGg6ICAgICAgICAgICAgICAwXG4gICwgaGVpZ2h0OiAgICAgICAgICAgICAwXG4gICAgLy8gcmFkaXVzXG4gICwgcjogICAgICAgICAgICAgICAgICAwXG4gICwgcng6ICAgICAgICAgICAgICAgICAwXG4gICwgcnk6ICAgICAgICAgICAgICAgICAwXG4gICAgLy8gZ3JhZGllbnRcbiAgLCBvZmZzZXQ6ICAgICAgICAgICAgIDBcbiAgLCAnc3RvcC1vcGFjaXR5JzogICAgIDFcbiAgLCAnc3RvcC1jb2xvcic6ICAgICAgICcjMDAwMDAwJ1xuICAgIC8vIHRleHRcbiAgLCAnZm9udC1zaXplJzogICAgICAgIDE2XG4gICwgJ2ZvbnQtZmFtaWx5JzogICAgICAnSGVsdmV0aWNhLCBBcmlhbCwgc2Fucy1zZXJpZidcbiAgLCAndGV4dC1hbmNob3InOiAgICAgICdzdGFydCdcbiAgfVxuXG59XG4vLyBNb2R1bGUgZm9yIGNvbG9yIGNvbnZlcnRpb25zXG5TVkcuQ29sb3IgPSBmdW5jdGlvbihjb2xvcikge1xuICB2YXIgbWF0Y2hcblxuICAvLyBpbml0aWFsaXplIGRlZmF1bHRzXG4gIHRoaXMuciA9IDBcbiAgdGhpcy5nID0gMFxuICB0aGlzLmIgPSAwXG5cbiAgaWYoIWNvbG9yKSByZXR1cm5cblxuICAvLyBwYXJzZSBjb2xvclxuICBpZiAodHlwZW9mIGNvbG9yID09PSAnc3RyaW5nJykge1xuICAgIGlmIChTVkcucmVnZXguaXNSZ2IudGVzdChjb2xvcikpIHtcbiAgICAgIC8vIGdldCByZ2IgdmFsdWVzXG4gICAgICBtYXRjaCA9IFNWRy5yZWdleC5yZ2IuZXhlYyhjb2xvci5yZXBsYWNlKC9cXHMvZywnJykpXG5cbiAgICAgIC8vIHBhcnNlIG51bWVyaWMgdmFsdWVzXG4gICAgICB0aGlzLnIgPSBwYXJzZUludChtYXRjaFsxXSlcbiAgICAgIHRoaXMuZyA9IHBhcnNlSW50KG1hdGNoWzJdKVxuICAgICAgdGhpcy5iID0gcGFyc2VJbnQobWF0Y2hbM10pXG5cbiAgICB9IGVsc2UgaWYgKFNWRy5yZWdleC5pc0hleC50ZXN0KGNvbG9yKSkge1xuICAgICAgLy8gZ2V0IGhleCB2YWx1ZXNcbiAgICAgIG1hdGNoID0gU1ZHLnJlZ2V4LmhleC5leGVjKGZ1bGxIZXgoY29sb3IpKVxuXG4gICAgICAvLyBwYXJzZSBudW1lcmljIHZhbHVlc1xuICAgICAgdGhpcy5yID0gcGFyc2VJbnQobWF0Y2hbMV0sIDE2KVxuICAgICAgdGhpcy5nID0gcGFyc2VJbnQobWF0Y2hbMl0sIDE2KVxuICAgICAgdGhpcy5iID0gcGFyc2VJbnQobWF0Y2hbM10sIDE2KVxuXG4gICAgfVxuXG4gIH0gZWxzZSBpZiAodHlwZW9mIGNvbG9yID09PSAnb2JqZWN0Jykge1xuICAgIHRoaXMuciA9IGNvbG9yLnJcbiAgICB0aGlzLmcgPSBjb2xvci5nXG4gICAgdGhpcy5iID0gY29sb3IuYlxuXG4gIH1cblxufVxuXG5TVkcuZXh0ZW5kKFNWRy5Db2xvciwge1xuICAvLyBEZWZhdWx0IHRvIGhleCBjb252ZXJzaW9uXG4gIHRvU3RyaW5nOiBmdW5jdGlvbigpIHtcbiAgICByZXR1cm4gdGhpcy50b0hleCgpXG4gIH1cbiAgLy8gQnVpbGQgaGV4IHZhbHVlXG4sIHRvSGV4OiBmdW5jdGlvbigpIHtcbiAgICByZXR1cm4gJyMnXG4gICAgICArIGNvbXBUb0hleCh0aGlzLnIpXG4gICAgICArIGNvbXBUb0hleCh0aGlzLmcpXG4gICAgICArIGNvbXBUb0hleCh0aGlzLmIpXG4gIH1cbiAgLy8gQnVpbGQgcmdiIHZhbHVlXG4sIHRvUmdiOiBmdW5jdGlvbigpIHtcbiAgICByZXR1cm4gJ3JnYignICsgW3RoaXMuciwgdGhpcy5nLCB0aGlzLmJdLmpvaW4oKSArICcpJ1xuICB9XG4gIC8vIENhbGN1bGF0ZSB0cnVlIGJyaWdodG5lc3NcbiwgYnJpZ2h0bmVzczogZnVuY3Rpb24oKSB7XG4gICAgcmV0dXJuICh0aGlzLnIgLyAyNTUgKiAwLjMwKVxuICAgICAgICAgKyAodGhpcy5nIC8gMjU1ICogMC41OSlcbiAgICAgICAgICsgKHRoaXMuYiAvIDI1NSAqIDAuMTEpXG4gIH1cbiAgLy8gTWFrZSBjb2xvciBtb3JwaGFibGVcbiwgbW9ycGg6IGZ1bmN0aW9uKGNvbG9yKSB7XG4gICAgdGhpcy5kZXN0aW5hdGlvbiA9IG5ldyBTVkcuQ29sb3IoY29sb3IpXG5cbiAgICByZXR1cm4gdGhpc1xuICB9XG4gIC8vIEdldCBtb3JwaGVkIGNvbG9yIGF0IGdpdmVuIHBvc2l0aW9uXG4sIGF0OiBmdW5jdGlvbihwb3MpIHtcbiAgICAvLyBtYWtlIHN1cmUgYSBkZXN0aW5hdGlvbiBpcyBkZWZpbmVkXG4gICAgaWYgKCF0aGlzLmRlc3RpbmF0aW9uKSByZXR1cm4gdGhpc1xuXG4gICAgLy8gbm9ybWFsaXNlIHBvc1xuICAgIHBvcyA9IHBvcyA8IDAgPyAwIDogcG9zID4gMSA/IDEgOiBwb3NcblxuICAgIC8vIGdlbmVyYXRlIG1vcnBoZWQgY29sb3JcbiAgICByZXR1cm4gbmV3IFNWRy5Db2xvcih7XG4gICAgICByOiB+fih0aGlzLnIgKyAodGhpcy5kZXN0aW5hdGlvbi5yIC0gdGhpcy5yKSAqIHBvcylcbiAgICAsIGc6IH5+KHRoaXMuZyArICh0aGlzLmRlc3RpbmF0aW9uLmcgLSB0aGlzLmcpICogcG9zKVxuICAgICwgYjogfn4odGhpcy5iICsgKHRoaXMuZGVzdGluYXRpb24uYiAtIHRoaXMuYikgKiBwb3MpXG4gICAgfSlcbiAgfVxuXG59KVxuXG4vLyBUZXN0ZXJzXG5cbi8vIFRlc3QgaWYgZ2l2ZW4gdmFsdWUgaXMgYSBjb2xvciBzdHJpbmdcblNWRy5Db2xvci50ZXN0ID0gZnVuY3Rpb24oY29sb3IpIHtcbiAgY29sb3IgKz0gJydcbiAgcmV0dXJuIFNWRy5yZWdleC5pc0hleC50ZXN0KGNvbG9yKVxuICAgICAgfHwgU1ZHLnJlZ2V4LmlzUmdiLnRlc3QoY29sb3IpXG59XG5cbi8vIFRlc3QgaWYgZ2l2ZW4gdmFsdWUgaXMgYSByZ2Igb2JqZWN0XG5TVkcuQ29sb3IuaXNSZ2IgPSBmdW5jdGlvbihjb2xvcikge1xuICByZXR1cm4gY29sb3IgJiYgdHlwZW9mIGNvbG9yLnIgPT0gJ251bWJlcidcbiAgICAgICAgICAgICAgICYmIHR5cGVvZiBjb2xvci5nID09ICdudW1iZXInXG4gICAgICAgICAgICAgICAmJiB0eXBlb2YgY29sb3IuYiA9PSAnbnVtYmVyJ1xufVxuXG4vLyBUZXN0IGlmIGdpdmVuIHZhbHVlIGlzIGEgY29sb3JcblNWRy5Db2xvci5pc0NvbG9yID0gZnVuY3Rpb24oY29sb3IpIHtcbiAgcmV0dXJuIFNWRy5Db2xvci5pc1JnYihjb2xvcikgfHwgU1ZHLkNvbG9yLnRlc3QoY29sb3IpXG59XG4vLyBNb2R1bGUgZm9yIGFycmF5IGNvbnZlcnNpb25cblNWRy5BcnJheSA9IGZ1bmN0aW9uKGFycmF5LCBmYWxsYmFjaykge1xuICBhcnJheSA9IChhcnJheSB8fCBbXSkudmFsdWVPZigpXG5cbiAgLy8gaWYgYXJyYXkgaXMgZW1wdHkgYW5kIGZhbGxiYWNrIGlzIHByb3ZpZGVkLCB1c2UgZmFsbGJhY2tcbiAgaWYgKGFycmF5Lmxlbmd0aCA9PSAwICYmIGZhbGxiYWNrKVxuICAgIGFycmF5ID0gZmFsbGJhY2sudmFsdWVPZigpXG5cbiAgLy8gcGFyc2UgYXJyYXlcbiAgdGhpcy52YWx1ZSA9IHRoaXMucGFyc2UoYXJyYXkpXG59XG5cblNWRy5leHRlbmQoU1ZHLkFycmF5LCB7XG4gIC8vIE1ha2UgYXJyYXkgbW9ycGhhYmxlXG4gIG1vcnBoOiBmdW5jdGlvbihhcnJheSkge1xuICAgIHRoaXMuZGVzdGluYXRpb24gPSB0aGlzLnBhcnNlKGFycmF5KVxuXG4gICAgLy8gbm9ybWFsaXplIGxlbmd0aCBvZiBhcnJheXNcbiAgICBpZiAodGhpcy52YWx1ZS5sZW5ndGggIT0gdGhpcy5kZXN0aW5hdGlvbi5sZW5ndGgpIHtcbiAgICAgIHZhciBsYXN0VmFsdWUgICAgICAgPSB0aGlzLnZhbHVlW3RoaXMudmFsdWUubGVuZ3RoIC0gMV1cbiAgICAgICAgLCBsYXN0RGVzdGluYXRpb24gPSB0aGlzLmRlc3RpbmF0aW9uW3RoaXMuZGVzdGluYXRpb24ubGVuZ3RoIC0gMV1cblxuICAgICAgd2hpbGUodGhpcy52YWx1ZS5sZW5ndGggPiB0aGlzLmRlc3RpbmF0aW9uLmxlbmd0aClcbiAgICAgICAgdGhpcy5kZXN0aW5hdGlvbi5wdXNoKGxhc3REZXN0aW5hdGlvbilcbiAgICAgIHdoaWxlKHRoaXMudmFsdWUubGVuZ3RoIDwgdGhpcy5kZXN0aW5hdGlvbi5sZW5ndGgpXG4gICAgICAgIHRoaXMudmFsdWUucHVzaChsYXN0VmFsdWUpXG4gICAgfVxuXG4gICAgcmV0dXJuIHRoaXNcbiAgfVxuICAvLyBDbGVhbiB1cCBhbnkgZHVwbGljYXRlIHBvaW50c1xuLCBzZXR0bGU6IGZ1bmN0aW9uKCkge1xuICAgIC8vIGZpbmQgYWxsIHVuaXF1ZSB2YWx1ZXNcbiAgICBmb3IgKHZhciBpID0gMCwgaWwgPSB0aGlzLnZhbHVlLmxlbmd0aCwgc2VlbiA9IFtdOyBpIDwgaWw7IGkrKylcbiAgICAgIGlmIChzZWVuLmluZGV4T2YodGhpcy52YWx1ZVtpXSkgPT0gLTEpXG4gICAgICAgIHNlZW4ucHVzaCh0aGlzLnZhbHVlW2ldKVxuXG4gICAgLy8gc2V0IG5ldyB2YWx1ZVxuICAgIHJldHVybiB0aGlzLnZhbHVlID0gc2VlblxuICB9XG4gIC8vIEdldCBtb3JwaGVkIGFycmF5IGF0IGdpdmVuIHBvc2l0aW9uXG4sIGF0OiBmdW5jdGlvbihwb3MpIHtcbiAgICAvLyBtYWtlIHN1cmUgYSBkZXN0aW5hdGlvbiBpcyBkZWZpbmVkXG4gICAgaWYgKCF0aGlzLmRlc3RpbmF0aW9uKSByZXR1cm4gdGhpc1xuXG4gICAgLy8gZ2VuZXJhdGUgbW9ycGhlZCBhcnJheVxuICAgIGZvciAodmFyIGkgPSAwLCBpbCA9IHRoaXMudmFsdWUubGVuZ3RoLCBhcnJheSA9IFtdOyBpIDwgaWw7IGkrKylcbiAgICAgIGFycmF5LnB1c2godGhpcy52YWx1ZVtpXSArICh0aGlzLmRlc3RpbmF0aW9uW2ldIC0gdGhpcy52YWx1ZVtpXSkgKiBwb3MpXG5cbiAgICByZXR1cm4gbmV3IFNWRy5BcnJheShhcnJheSlcbiAgfVxuICAvLyBDb252ZXJ0IGFycmF5IHRvIHN0cmluZ1xuLCB0b1N0cmluZzogZnVuY3Rpb24oKSB7XG4gICAgcmV0dXJuIHRoaXMudmFsdWUuam9pbignICcpXG4gIH1cbiAgLy8gUmVhbCB2YWx1ZVxuLCB2YWx1ZU9mOiBmdW5jdGlvbigpIHtcbiAgICByZXR1cm4gdGhpcy52YWx1ZVxuICB9XG4gIC8vIFBhcnNlIHdoaXRlc3BhY2Ugc2VwYXJhdGVkIHN0cmluZ1xuLCBwYXJzZTogZnVuY3Rpb24oYXJyYXkpIHtcbiAgICBhcnJheSA9IGFycmF5LnZhbHVlT2YoKVxuXG4gICAgLy8gaWYgYWxyZWFkeSBpcyBhbiBhcnJheSwgbm8gbmVlZCB0byBwYXJzZSBpdFxuICAgIGlmIChBcnJheS5pc0FycmF5KGFycmF5KSkgcmV0dXJuIGFycmF5XG5cbiAgICByZXR1cm4gdGhpcy5zcGxpdChhcnJheSlcbiAgfVxuICAvLyBTdHJpcCB1bm5lY2Vzc2FyeSB3aGl0ZXNwYWNlXG4sIHNwbGl0OiBmdW5jdGlvbihzdHJpbmcpIHtcbiAgICByZXR1cm4gc3RyaW5nLnRyaW0oKS5zcGxpdCgvXFxzKy8pXG4gIH1cbiAgLy8gUmV2ZXJzZSBhcnJheVxuLCByZXZlcnNlOiBmdW5jdGlvbigpIHtcbiAgICB0aGlzLnZhbHVlLnJldmVyc2UoKVxuXG4gICAgcmV0dXJuIHRoaXNcbiAgfVxuXG59KVxuLy8gUG9seSBwb2ludHMgYXJyYXlcblNWRy5Qb2ludEFycmF5ID0gZnVuY3Rpb24oYXJyYXksIGZhbGxiYWNrKSB7XG4gIHRoaXMuY29uc3RydWN0b3IuY2FsbCh0aGlzLCBhcnJheSwgZmFsbGJhY2sgfHwgW1swLDBdXSlcbn1cblxuLy8gSW5oZXJpdCBmcm9tIFNWRy5BcnJheVxuU1ZHLlBvaW50QXJyYXkucHJvdG90eXBlID0gbmV3IFNWRy5BcnJheVxuXG5TVkcuZXh0ZW5kKFNWRy5Qb2ludEFycmF5LCB7XG4gIC8vIENvbnZlcnQgYXJyYXkgdG8gc3RyaW5nXG4gIHRvU3RyaW5nOiBmdW5jdGlvbigpIHtcbiAgICAvLyBjb252ZXJ0IHRvIGEgcG9seSBwb2ludCBzdHJpbmdcbiAgICBmb3IgKHZhciBpID0gMCwgaWwgPSB0aGlzLnZhbHVlLmxlbmd0aCwgYXJyYXkgPSBbXTsgaSA8IGlsOyBpKyspXG4gICAgICBhcnJheS5wdXNoKHRoaXMudmFsdWVbaV0uam9pbignLCcpKVxuXG4gICAgcmV0dXJuIGFycmF5LmpvaW4oJyAnKVxuICB9XG4gIC8vIENvbnZlcnQgYXJyYXkgdG8gbGluZSBvYmplY3RcbiwgdG9MaW5lOiBmdW5jdGlvbigpIHtcbiAgICByZXR1cm4ge1xuICAgICAgeDE6IHRoaXMudmFsdWVbMF1bMF1cbiAgICAsIHkxOiB0aGlzLnZhbHVlWzBdWzFdXG4gICAgLCB4MjogdGhpcy52YWx1ZVsxXVswXVxuICAgICwgeTI6IHRoaXMudmFsdWVbMV1bMV1cbiAgICB9XG4gIH1cbiAgLy8gR2V0IG1vcnBoZWQgYXJyYXkgYXQgZ2l2ZW4gcG9zaXRpb25cbiwgYXQ6IGZ1bmN0aW9uKHBvcykge1xuICAgIC8vIG1ha2Ugc3VyZSBhIGRlc3RpbmF0aW9uIGlzIGRlZmluZWRcbiAgICBpZiAoIXRoaXMuZGVzdGluYXRpb24pIHJldHVybiB0aGlzXG5cbiAgICAvLyBnZW5lcmF0ZSBtb3JwaGVkIHBvaW50IHN0cmluZ1xuICAgIGZvciAodmFyIGkgPSAwLCBpbCA9IHRoaXMudmFsdWUubGVuZ3RoLCBhcnJheSA9IFtdOyBpIDwgaWw7IGkrKylcbiAgICAgIGFycmF5LnB1c2goW1xuICAgICAgICB0aGlzLnZhbHVlW2ldWzBdICsgKHRoaXMuZGVzdGluYXRpb25baV1bMF0gLSB0aGlzLnZhbHVlW2ldWzBdKSAqIHBvc1xuICAgICAgLCB0aGlzLnZhbHVlW2ldWzFdICsgKHRoaXMuZGVzdGluYXRpb25baV1bMV0gLSB0aGlzLnZhbHVlW2ldWzFdKSAqIHBvc1xuICAgICAgXSlcblxuICAgIHJldHVybiBuZXcgU1ZHLlBvaW50QXJyYXkoYXJyYXkpXG4gIH1cbiAgLy8gUGFyc2UgcG9pbnQgc3RyaW5nXG4sIHBhcnNlOiBmdW5jdGlvbihhcnJheSkge1xuICAgIHZhciBwb2ludHMgPSBbXVxuXG4gICAgYXJyYXkgPSBhcnJheS52YWx1ZU9mKClcblxuICAgIC8vIGlmIGFscmVhZHkgaXMgYW4gYXJyYXksIG5vIG5lZWQgdG8gcGFyc2UgaXRcbiAgICBpZiAoQXJyYXkuaXNBcnJheShhcnJheSkpIHJldHVybiBhcnJheVxuXG4gICAgLy8gcGFyc2UgcG9pbnRzXG4gICAgYXJyYXkgPSBhcnJheS50cmltKCkuc3BsaXQoL1xccyt8LC8pXG5cbiAgICAvLyB2YWxpZGF0ZSBwb2ludHMgLSBodHRwczovL3N2Z3dnLm9yZy9zdmcyLWRyYWZ0L3NoYXBlcy5odG1sI0RhdGFUeXBlUG9pbnRzXG4gICAgLy8gT2RkIG51bWJlciBvZiBjb29yZGluYXRlcyBpcyBhbiBlcnJvci4gSW4gc3VjaCBjYXNlcywgZHJvcCB0aGUgbGFzdCBvZGQgY29vcmRpbmF0ZS5cbiAgICBpZiAoYXJyYXkubGVuZ3RoICUgMiAhPT0gMCkgYXJyYXkucG9wKClcblxuICAgIC8vIHdyYXAgcG9pbnRzIGluIHR3by10dXBsZXMgYW5kIHBhcnNlIHBvaW50cyBhcyBmbG9hdHNcbiAgICBmb3IodmFyIGkgPSAwLCBsZW4gPSBhcnJheS5sZW5ndGg7IGkgPCBsZW47IGkgPSBpICsgMilcbiAgICAgIHBvaW50cy5wdXNoKFsgcGFyc2VGbG9hdChhcnJheVtpXSksIHBhcnNlRmxvYXQoYXJyYXlbaSsxXSkgXSlcblxuICAgIHJldHVybiBwb2ludHNcbiAgfVxuICAvLyBNb3ZlIHBvaW50IHN0cmluZ1xuLCBtb3ZlOiBmdW5jdGlvbih4LCB5KSB7XG4gICAgdmFyIGJveCA9IHRoaXMuYmJveCgpXG5cbiAgICAvLyBnZXQgcmVsYXRpdmUgb2Zmc2V0XG4gICAgeCAtPSBib3gueFxuICAgIHkgLT0gYm94LnlcblxuICAgIC8vIG1vdmUgZXZlcnkgcG9pbnRcbiAgICBpZiAoIWlzTmFOKHgpICYmICFpc05hTih5KSlcbiAgICAgIGZvciAodmFyIGkgPSB0aGlzLnZhbHVlLmxlbmd0aCAtIDE7IGkgPj0gMDsgaS0tKVxuICAgICAgICB0aGlzLnZhbHVlW2ldID0gW3RoaXMudmFsdWVbaV1bMF0gKyB4LCB0aGlzLnZhbHVlW2ldWzFdICsgeV1cblxuICAgIHJldHVybiB0aGlzXG4gIH1cbiAgLy8gUmVzaXplIHBvbHkgc3RyaW5nXG4sIHNpemU6IGZ1bmN0aW9uKHdpZHRoLCBoZWlnaHQpIHtcbiAgICB2YXIgaSwgYm94ID0gdGhpcy5iYm94KClcblxuICAgIC8vIHJlY2FsY3VsYXRlIHBvc2l0aW9uIG9mIGFsbCBwb2ludHMgYWNjb3JkaW5nIHRvIG5ldyBzaXplXG4gICAgZm9yIChpID0gdGhpcy52YWx1ZS5sZW5ndGggLSAxOyBpID49IDA7IGktLSkge1xuICAgICAgdGhpcy52YWx1ZVtpXVswXSA9ICgodGhpcy52YWx1ZVtpXVswXSAtIGJveC54KSAqIHdpZHRoKSAgLyBib3gud2lkdGggICsgYm94LnhcbiAgICAgIHRoaXMudmFsdWVbaV1bMV0gPSAoKHRoaXMudmFsdWVbaV1bMV0gLSBib3gueSkgKiBoZWlnaHQpIC8gYm94LmhlaWdodCArIGJveC55XG4gICAgfVxuXG4gICAgcmV0dXJuIHRoaXNcbiAgfVxuICAvLyBHZXQgYm91bmRpbmcgYm94IG9mIHBvaW50c1xuLCBiYm94OiBmdW5jdGlvbigpIHtcbiAgICBTVkcucGFyc2VyLnBvbHkuc2V0QXR0cmlidXRlKCdwb2ludHMnLCB0aGlzLnRvU3RyaW5nKCkpXG5cbiAgICByZXR1cm4gU1ZHLnBhcnNlci5wb2x5LmdldEJCb3goKVxuICB9XG5cbn0pXG4vLyBQYXRoIHBvaW50cyBhcnJheVxuU1ZHLlBhdGhBcnJheSA9IGZ1bmN0aW9uKGFycmF5LCBmYWxsYmFjaykge1xuICB0aGlzLmNvbnN0cnVjdG9yLmNhbGwodGhpcywgYXJyYXksIGZhbGxiYWNrIHx8IFtbJ00nLCAwLCAwXV0pXG59XG5cbi8vIEluaGVyaXQgZnJvbSBTVkcuQXJyYXlcblNWRy5QYXRoQXJyYXkucHJvdG90eXBlID0gbmV3IFNWRy5BcnJheVxuXG5TVkcuZXh0ZW5kKFNWRy5QYXRoQXJyYXksIHtcbiAgLy8gQ29udmVydCBhcnJheSB0byBzdHJpbmdcbiAgdG9TdHJpbmc6IGZ1bmN0aW9uKCkge1xuICAgIHJldHVybiBhcnJheVRvU3RyaW5nKHRoaXMudmFsdWUpXG4gIH1cbiAgLy8gTW92ZSBwYXRoIHN0cmluZ1xuLCBtb3ZlOiBmdW5jdGlvbih4LCB5KSB7XG4gICAgLy8gZ2V0IGJvdW5kaW5nIGJveCBvZiBjdXJyZW50IHNpdHVhdGlvblxuICAgIHZhciBib3ggPSB0aGlzLmJib3goKVxuXG4gICAgLy8gZ2V0IHJlbGF0aXZlIG9mZnNldFxuICAgIHggLT0gYm94LnhcbiAgICB5IC09IGJveC55XG5cbiAgICBpZiAoIWlzTmFOKHgpICYmICFpc05hTih5KSkge1xuICAgICAgLy8gbW92ZSBldmVyeSBwb2ludFxuICAgICAgZm9yICh2YXIgbCwgaSA9IHRoaXMudmFsdWUubGVuZ3RoIC0gMTsgaSA+PSAwOyBpLS0pIHtcbiAgICAgICAgbCA9IHRoaXMudmFsdWVbaV1bMF1cblxuICAgICAgICBpZiAobCA9PSAnTScgfHwgbCA9PSAnTCcgfHwgbCA9PSAnVCcpICB7XG4gICAgICAgICAgdGhpcy52YWx1ZVtpXVsxXSArPSB4XG4gICAgICAgICAgdGhpcy52YWx1ZVtpXVsyXSArPSB5XG5cbiAgICAgICAgfSBlbHNlIGlmIChsID09ICdIJykgIHtcbiAgICAgICAgICB0aGlzLnZhbHVlW2ldWzFdICs9IHhcblxuICAgICAgICB9IGVsc2UgaWYgKGwgPT0gJ1YnKSAge1xuICAgICAgICAgIHRoaXMudmFsdWVbaV1bMV0gKz0geVxuXG4gICAgICAgIH0gZWxzZSBpZiAobCA9PSAnQycgfHwgbCA9PSAnUycgfHwgbCA9PSAnUScpICB7XG4gICAgICAgICAgdGhpcy52YWx1ZVtpXVsxXSArPSB4XG4gICAgICAgICAgdGhpcy52YWx1ZVtpXVsyXSArPSB5XG4gICAgICAgICAgdGhpcy52YWx1ZVtpXVszXSArPSB4XG4gICAgICAgICAgdGhpcy52YWx1ZVtpXVs0XSArPSB5XG5cbiAgICAgICAgICBpZiAobCA9PSAnQycpICB7XG4gICAgICAgICAgICB0aGlzLnZhbHVlW2ldWzVdICs9IHhcbiAgICAgICAgICAgIHRoaXMudmFsdWVbaV1bNl0gKz0geVxuICAgICAgICAgIH1cblxuICAgICAgICB9IGVsc2UgaWYgKGwgPT0gJ0EnKSAge1xuICAgICAgICAgIHRoaXMudmFsdWVbaV1bNl0gKz0geFxuICAgICAgICAgIHRoaXMudmFsdWVbaV1bN10gKz0geVxuICAgICAgICB9XG5cbiAgICAgIH1cbiAgICB9XG5cbiAgICByZXR1cm4gdGhpc1xuICB9XG4gIC8vIFJlc2l6ZSBwYXRoIHN0cmluZ1xuLCBzaXplOiBmdW5jdGlvbih3aWR0aCwgaGVpZ2h0KSB7XG4gICAgLy8gZ2V0IGJvdW5kaW5nIGJveCBvZiBjdXJyZW50IHNpdHVhdGlvblxuICAgIHZhciBpLCBsLCBib3ggPSB0aGlzLmJib3goKVxuXG4gICAgLy8gcmVjYWxjdWxhdGUgcG9zaXRpb24gb2YgYWxsIHBvaW50cyBhY2NvcmRpbmcgdG8gbmV3IHNpemVcbiAgICBmb3IgKGkgPSB0aGlzLnZhbHVlLmxlbmd0aCAtIDE7IGkgPj0gMDsgaS0tKSB7XG4gICAgICBsID0gdGhpcy52YWx1ZVtpXVswXVxuXG4gICAgICBpZiAobCA9PSAnTScgfHwgbCA9PSAnTCcgfHwgbCA9PSAnVCcpICB7XG4gICAgICAgIHRoaXMudmFsdWVbaV1bMV0gPSAoKHRoaXMudmFsdWVbaV1bMV0gLSBib3gueCkgKiB3aWR0aCkgIC8gYm94LndpZHRoICArIGJveC54XG4gICAgICAgIHRoaXMudmFsdWVbaV1bMl0gPSAoKHRoaXMudmFsdWVbaV1bMl0gLSBib3gueSkgKiBoZWlnaHQpIC8gYm94LmhlaWdodCArIGJveC55XG5cbiAgICAgIH0gZWxzZSBpZiAobCA9PSAnSCcpICB7XG4gICAgICAgIHRoaXMudmFsdWVbaV1bMV0gPSAoKHRoaXMudmFsdWVbaV1bMV0gLSBib3gueCkgKiB3aWR0aCkgIC8gYm94LndpZHRoICArIGJveC54XG5cbiAgICAgIH0gZWxzZSBpZiAobCA9PSAnVicpICB7XG4gICAgICAgIHRoaXMudmFsdWVbaV1bMV0gPSAoKHRoaXMudmFsdWVbaV1bMV0gLSBib3gueSkgKiBoZWlnaHQpIC8gYm94LmhlaWdodCArIGJveC55XG5cbiAgICAgIH0gZWxzZSBpZiAobCA9PSAnQycgfHwgbCA9PSAnUycgfHwgbCA9PSAnUScpICB7XG4gICAgICAgIHRoaXMudmFsdWVbaV1bMV0gPSAoKHRoaXMudmFsdWVbaV1bMV0gLSBib3gueCkgKiB3aWR0aCkgIC8gYm94LndpZHRoICArIGJveC54XG4gICAgICAgIHRoaXMudmFsdWVbaV1bMl0gPSAoKHRoaXMudmFsdWVbaV1bMl0gLSBib3gueSkgKiBoZWlnaHQpIC8gYm94LmhlaWdodCArIGJveC55XG4gICAgICAgIHRoaXMudmFsdWVbaV1bM10gPSAoKHRoaXMudmFsdWVbaV1bM10gLSBib3gueCkgKiB3aWR0aCkgIC8gYm94LndpZHRoICArIGJveC54XG4gICAgICAgIHRoaXMudmFsdWVbaV1bNF0gPSAoKHRoaXMudmFsdWVbaV1bNF0gLSBib3gueSkgKiBoZWlnaHQpIC8gYm94LmhlaWdodCArIGJveC55XG5cbiAgICAgICAgaWYgKGwgPT0gJ0MnKSAge1xuICAgICAgICAgIHRoaXMudmFsdWVbaV1bNV0gPSAoKHRoaXMudmFsdWVbaV1bNV0gLSBib3gueCkgKiB3aWR0aCkgIC8gYm94LndpZHRoICArIGJveC54XG4gICAgICAgICAgdGhpcy52YWx1ZVtpXVs2XSA9ICgodGhpcy52YWx1ZVtpXVs2XSAtIGJveC55KSAqIGhlaWdodCkgLyBib3guaGVpZ2h0ICsgYm94LnlcbiAgICAgICAgfVxuXG4gICAgICB9IGVsc2UgaWYgKGwgPT0gJ0EnKSAge1xuICAgICAgICAvLyByZXNpemUgcmFkaWlcbiAgICAgICAgdGhpcy52YWx1ZVtpXVsxXSA9ICh0aGlzLnZhbHVlW2ldWzFdICogd2lkdGgpICAvIGJveC53aWR0aFxuICAgICAgICB0aGlzLnZhbHVlW2ldWzJdID0gKHRoaXMudmFsdWVbaV1bMl0gKiBoZWlnaHQpIC8gYm94LmhlaWdodFxuXG4gICAgICAgIC8vIG1vdmUgcG9zaXRpb24gdmFsdWVzXG4gICAgICAgIHRoaXMudmFsdWVbaV1bNl0gPSAoKHRoaXMudmFsdWVbaV1bNl0gLSBib3gueCkgKiB3aWR0aCkgIC8gYm94LndpZHRoICArIGJveC54XG4gICAgICAgIHRoaXMudmFsdWVbaV1bN10gPSAoKHRoaXMudmFsdWVbaV1bN10gLSBib3gueSkgKiBoZWlnaHQpIC8gYm94LmhlaWdodCArIGJveC55XG4gICAgICB9XG5cbiAgICB9XG5cbiAgICByZXR1cm4gdGhpc1xuICB9XG4gIC8vIFRlc3QgaWYgdGhlIHBhc3NlZCBwYXRoIGFycmF5IHVzZSB0aGUgc2FtZSBwYXRoIGRhdGEgY29tbWFuZHMgYXMgdGhpcyBwYXRoIGFycmF5XG4sIGVxdWFsQ29tbWFuZHM6IGZ1bmN0aW9uKHBhdGhBcnJheSkge1xuICAgIHZhciBpLCBpbCwgZXF1YWxDb21tYW5kc1xuXG4gICAgcGF0aEFycmF5ID0gbmV3IFNWRy5QYXRoQXJyYXkocGF0aEFycmF5KVxuXG4gICAgZXF1YWxDb21tYW5kcyA9IHRoaXMudmFsdWUubGVuZ3RoID09PSBwYXRoQXJyYXkudmFsdWUubGVuZ3RoXG4gICAgZm9yKGkgPSAwLCBpbCA9IHRoaXMudmFsdWUubGVuZ3RoOyBlcXVhbENvbW1hbmRzICYmIGkgPCBpbDsgaSsrKSB7XG4gICAgICBlcXVhbENvbW1hbmRzID0gdGhpcy52YWx1ZVtpXVswXSA9PT0gcGF0aEFycmF5LnZhbHVlW2ldWzBdXG4gICAgfVxuXG4gICAgcmV0dXJuIGVxdWFsQ29tbWFuZHNcbiAgfVxuICAvLyBNYWtlIHBhdGggYXJyYXkgbW9ycGhhYmxlXG4sIG1vcnBoOiBmdW5jdGlvbihwYXRoQXJyYXkpIHtcbiAgICBwYXRoQXJyYXkgPSBuZXcgU1ZHLlBhdGhBcnJheShwYXRoQXJyYXkpXG5cbiAgICBpZih0aGlzLmVxdWFsQ29tbWFuZHMocGF0aEFycmF5KSkge1xuICAgICAgdGhpcy5kZXN0aW5hdGlvbiA9IHBhdGhBcnJheVxuICAgIH0gZWxzZSB7XG4gICAgICB0aGlzLmRlc3RpbmF0aW9uID0gbnVsbFxuICAgIH1cblxuICAgIHJldHVybiB0aGlzXG4gIH1cbiAgLy8gR2V0IG1vcnBoZWQgcGF0aCBhcnJheSBhdCBnaXZlbiBwb3NpdGlvblxuLCBhdDogZnVuY3Rpb24ocG9zKSB7XG4gICAgLy8gbWFrZSBzdXJlIGEgZGVzdGluYXRpb24gaXMgZGVmaW5lZFxuICAgIGlmICghdGhpcy5kZXN0aW5hdGlvbikgcmV0dXJuIHRoaXNcblxuICAgIHZhciBzb3VyY2VBcnJheSA9IHRoaXMudmFsdWVcbiAgICAgICwgZGVzdGluYXRpb25BcnJheSA9IHRoaXMuZGVzdGluYXRpb24udmFsdWVcbiAgICAgICwgYXJyYXkgPSBbXSwgcGF0aEFycmF5ID0gbmV3IFNWRy5QYXRoQXJyYXkoKVxuICAgICAgLCBpLCBpbCwgaiwgamxcblxuICAgIC8vIEFuaW1hdGUgaGFzIHNwZWNpZmllZCBpbiB0aGUgU1ZHIHNwZWNcbiAgICAvLyBTZWU6IGh0dHBzOi8vd3d3LnczLm9yZy9UUi9TVkcxMS9wYXRocy5odG1sI1BhdGhFbGVtZW50XG4gICAgZm9yIChpID0gMCwgaWwgPSBzb3VyY2VBcnJheS5sZW5ndGg7IGkgPCBpbDsgaSsrKSB7XG4gICAgICBhcnJheVtpXSA9IFtzb3VyY2VBcnJheVtpXVswXV1cbiAgICAgIGZvcihqID0gMSwgamwgPSBzb3VyY2VBcnJheVtpXS5sZW5ndGg7IGogPCBqbDsgaisrKSB7XG4gICAgICAgIGFycmF5W2ldW2pdID0gc291cmNlQXJyYXlbaV1bal0gKyAoZGVzdGluYXRpb25BcnJheVtpXVtqXSAtIHNvdXJjZUFycmF5W2ldW2pdKSAqIHBvc1xuICAgICAgfVxuICAgICAgLy8gRm9yIHRoZSB0d28gZmxhZ3Mgb2YgdGhlIGVsbGlwdGljYWwgYXJjIGNvbW1hbmQsIHRoZSBTVkcgc3BlYyBzYXk6XG4gICAgICAvLyBGbGFncyBhbmQgYm9vbGVhbnMgYXJlIGludGVycG9sYXRlZCBhcyBmcmFjdGlvbnMgYmV0d2VlbiB6ZXJvIGFuZCBvbmUsIHdpdGggYW55IG5vbi16ZXJvIHZhbHVlIGNvbnNpZGVyZWQgdG8gYmUgYSB2YWx1ZSBvZiBvbmUvdHJ1ZVxuICAgICAgLy8gRWxsaXB0aWNhbCBhcmMgY29tbWFuZCBhcyBhbiBhcnJheSBmb2xsb3dlZCBieSBjb3JyZXNwb25kaW5nIGluZGV4ZXM6XG4gICAgICAvLyBbJ0EnLCByeCwgcnksIHgtYXhpcy1yb3RhdGlvbiwgbGFyZ2UtYXJjLWZsYWcsIHN3ZWVwLWZsYWcsIHgsIHldXG4gICAgICAvLyAgIDAgICAgMSAgIDIgICAgICAgIDMgICAgICAgICAgICAgICAgIDQgICAgICAgICAgICAgNSAgICAgIDYgIDdcbiAgICAgIGlmKGFycmF5W2ldWzBdID09PSAnQScpIHtcbiAgICAgICAgYXJyYXlbaV1bNF0gPSArKGFycmF5W2ldWzRdICE9IDApXG4gICAgICAgIGFycmF5W2ldWzVdID0gKyhhcnJheVtpXVs1XSAhPSAwKVxuICAgICAgfVxuICAgIH1cblxuICAgIC8vIERpcmVjdGx5IG1vZGlmeSB0aGUgdmFsdWUgb2YgYSBwYXRoIGFycmF5LCB0aGlzIGlzIGRvbmUgdGhpcyB3YXkgZm9yIHBlcmZvcm1hbmNlXG4gICAgcGF0aEFycmF5LnZhbHVlID0gYXJyYXlcbiAgICByZXR1cm4gcGF0aEFycmF5XG4gIH1cbiAgLy8gQWJzb2x1dGl6ZSBhbmQgcGFyc2UgcGF0aCB0byBhcnJheVxuLCBwYXJzZTogZnVuY3Rpb24oYXJyYXkpIHtcbiAgICAvLyBpZiBpdCdzIGFscmVhZHkgYSBwYXRoYXJyYXksIG5vIG5lZWQgdG8gcGFyc2UgaXRcbiAgICBpZiAoYXJyYXkgaW5zdGFuY2VvZiBTVkcuUGF0aEFycmF5KSByZXR1cm4gYXJyYXkudmFsdWVPZigpXG5cbiAgICAvLyBwcmVwYXJlIGZvciBwYXJzaW5nXG4gICAgdmFyIGksIHgwLCB5MCwgcywgc2VnLCBhcnJcbiAgICAgICwgeCA9IDBcbiAgICAgICwgeSA9IDBcbiAgICAgICwgcGFyYW1DbnQgPSB7ICdNJzoyLCAnTCc6MiwgJ0gnOjEsICdWJzoxLCAnQyc6NiwgJ1MnOjQsICdRJzo0LCAnVCc6MiwgJ0EnOjcgfVxuXG4gICAgaWYodHlwZW9mIGFycmF5ID09ICdzdHJpbmcnKXtcblxuICAgICAgYXJyYXkgPSBhcnJheVxuICAgICAgICAucmVwbGFjZShTVkcucmVnZXgubmVnRXhwLCAnWCcpICAgICAgICAgLy8gcmVwbGFjZSBhbGwgbmVnYXRpdmUgZXhwb25lbnRzIHdpdGggY2VydGFpbiBjaGFyXG4gICAgICAgIC5yZXBsYWNlKFNWRy5yZWdleC5wYXRoTGV0dGVycywgJyAkJiAnKSAvLyBwdXQgc29tZSByb29tIGJldHdlZW4gbGV0dGVycyBhbmQgbnVtYmVyc1xuICAgICAgICAucmVwbGFjZShTVkcucmVnZXguaHlwaGVuLCAnIC0nKSAgICAgICAgLy8gYWRkIHNwYWNlIGJlZm9yZSBoeXBoZW5cbiAgICAgICAgLnJlcGxhY2UoU1ZHLnJlZ2V4LmNvbW1hLCAnICcpICAgICAgICAgIC8vIHVuaWZ5IGFsbCBzcGFjZXNcbiAgICAgICAgLnJlcGxhY2UoU1ZHLnJlZ2V4LlgsICdlLScpICAgICAgICAgICAgIC8vIGFkZCBiYWNrIHRoZSBleHBvZW50XG4gICAgICAgIC50cmltKCkgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAvLyB0cmltXG4gICAgICAgIC5zcGxpdChTVkcucmVnZXgud2hpdGVzcGFjZXMpICAgICAgICAgICAvLyBzcGxpdCBpbnRvIGFycmF5XG5cbiAgICAgIC8vIGF0IHRoaXMgcGxhY2UgdGhlcmUgY291bGQgYmUgcGFydHMgbGlrZSBbJzMuMTI0Ljg1NC4zMiddIGJlY2F1c2Ugd2UgY291bGQgbm90IGRldGVybWluZSB0aGUgcG9pbnQgYXMgc2VwZXJhdG9yIHRpbGwgbm93XG4gICAgICAvLyB3ZSBmaXggdGhpcyBlbGVtZW50cyBpbiB0aGUgbmV4dCBsb29wXG4gICAgICBmb3IoaSA9IGFycmF5Lmxlbmd0aDsgLS1pOyl7XG4gICAgICAgIGlmKGFycmF5W2ldLmluZGV4T2YoJy4nKSAhPSBhcnJheVtpXS5sYXN0SW5kZXhPZignLicpKXtcbiAgICAgICAgICB2YXIgc3BsaXQgPSBhcnJheVtpXS5zcGxpdCgnLicpIC8vIHNwbGl0IGF0IHRoZSBwb2ludFxuICAgICAgICAgIHZhciBmaXJzdCA9IFtzcGxpdC5zaGlmdCgpLCBzcGxpdC5zaGlmdCgpXS5qb2luKCcuJykgLy8gam9pbiB0aGUgZmlyc3QgbnVtYmVyIHRvZ2V0aGVyXG4gICAgICAgICAgYXJyYXkuc3BsaWNlLmFwcGx5KGFycmF5LCBbaSwgMV0uY29uY2F0KGZpcnN0LCBzcGxpdC5tYXAoZnVuY3Rpb24oZWwpeyByZXR1cm4gJy4nK2VsIH0pKSkgLy8gYWRkIGZpcnN0IGFuZCBhbGwgb3RoZXIgZW50cmllcyBiYWNrIHRvIGFycmF5XG4gICAgICAgIH1cbiAgICAgIH1cblxuICAgIH1lbHNle1xuICAgICAgYXJyYXkgPSBhcnJheS5yZWR1Y2UoZnVuY3Rpb24ocHJldiwgY3Vycil7XG4gICAgICAgIHJldHVybiBbXS5jb25jYXQuYXBwbHkocHJldiwgY3VycilcbiAgICAgIH0sIFtdKVxuICAgIH1cblxuICAgIC8vIGFycmF5IG5vdyBpcyBhbiBhcnJheSBjb250YWluaW5nIGFsbCBwYXJ0cyBvZiBhIHBhdGggZS5nLiBbJ00nLCAnMCcsICcwJywgJ0wnLCAnMzAnLCAnMzAnIC4uLl1cblxuICAgIHZhciBhcnIgPSBbXVxuXG4gICAgZG97XG5cbiAgICAgIC8vIFRlc3QgaWYgd2UgaGF2ZSBhIHBhdGggbGV0dGVyXG4gICAgICBpZihTVkcucmVnZXguaXNQYXRoTGV0dGVyLnRlc3QoYXJyYXlbMF0pKXtcbiAgICAgICAgcyA9IGFycmF5WzBdXG4gICAgICAgIGFycmF5LnNoaWZ0KClcbiAgICAgIC8vIElmIGxhc3QgbGV0dGVyIHdhcyBhIG1vdmUgY29tbWFuZCBhbmQgd2UgZ290IG5vIG5ldywgaXQgZGVmYXVsdHMgdG8gW0xdaW5lXG4gICAgICB9ZWxzZSBpZihzID09ICdNJyl7XG4gICAgICAgIHMgPSAnTCdcbiAgICAgIH1lbHNlIGlmKHMgPT0gJ20nKXtcbiAgICAgICAgcyA9ICdsJ1xuICAgICAgfVxuXG4gICAgICAvLyBhZGQgcGF0aCBsZXR0ZXIgYXMgZmlyc3QgZWxlbWVudFxuICAgICAgc2VnID0gW3MudG9VcHBlckNhc2UoKV1cblxuICAgICAgLy8gcHVzaCBhbGwgbmVjZXNzYXJ5IHBhcmFtZXRlcnMgdG8gc2VnbWVudFxuICAgICAgZm9yKGkgPSAwOyBpIDwgcGFyYW1DbnRbc2VnWzBdXTsgKytpKXtcbiAgICAgICAgc2VnLnB1c2gocGFyc2VGbG9hdChhcnJheS5zaGlmdCgpKSlcbiAgICAgIH1cblxuICAgICAgLy8gdXBwZXIgY2FzZVxuICAgICAgaWYocyA9PSBzZWdbMF0pe1xuXG4gICAgICAgIGlmKHMgPT0gJ00nIHx8IHMgPT0gJ0wnIHx8IHMgPT0gJ0MnIHx8IHMgPT0gJ1EnIHx8IHMgPT0gJ1MnIHx8IHMgPT0gJ1QnKXtcbiAgICAgICAgICB4ID0gc2VnW3BhcmFtQ250W3NlZ1swXV0tMV1cbiAgICAgICAgICB5ID0gc2VnW3BhcmFtQ250W3NlZ1swXV1dXG4gICAgICAgIH1lbHNlIGlmKHMgPT0gJ1YnKXtcbiAgICAgICAgICB5ID0gc2VnWzFdXG4gICAgICAgIH1lbHNlIGlmKHMgPT0gJ0gnKXtcbiAgICAgICAgICB4ID0gc2VnWzFdXG4gICAgICAgIH1lbHNlIGlmKHMgPT0gJ0EnKXtcbiAgICAgICAgICB4ID0gc2VnWzZdXG4gICAgICAgICAgeSA9IHNlZ1s3XVxuICAgICAgICB9XG5cbiAgICAgIC8vIGxvd2VyIGNhc2VcbiAgICAgIH1lbHNle1xuXG4gICAgICAgIC8vIGNvbnZlcnQgcmVsYXRpdmUgdG8gYWJzb2x1dGUgdmFsdWVzXG4gICAgICAgIGlmKHMgPT0gJ20nIHx8IHMgPT0gJ2wnIHx8IHMgPT0gJ2MnIHx8IHMgPT0gJ3MnIHx8IHMgPT0gJ3EnIHx8IHMgPT0gJ3QnKXtcblxuICAgICAgICAgIHNlZ1sxXSArPSB4XG4gICAgICAgICAgc2VnWzJdICs9IHlcblxuICAgICAgICAgIGlmKHNlZ1szXSAhPSBudWxsKXtcbiAgICAgICAgICAgIHNlZ1szXSArPSB4XG4gICAgICAgICAgICBzZWdbNF0gKz0geVxuICAgICAgICAgIH1cblxuICAgICAgICAgIGlmKHNlZ1s1XSAhPSBudWxsKXtcbiAgICAgICAgICAgIHNlZ1s1XSArPSB4XG4gICAgICAgICAgICBzZWdbNl0gKz0geVxuICAgICAgICAgIH1cblxuICAgICAgICAgIC8vIG1vdmUgcG9pbnRlclxuICAgICAgICAgIHggPSBzZWdbcGFyYW1DbnRbc2VnWzBdXS0xXVxuICAgICAgICAgIHkgPSBzZWdbcGFyYW1DbnRbc2VnWzBdXV1cblxuICAgICAgICB9ZWxzZSBpZihzID09ICd2Jyl7XG4gICAgICAgICAgc2VnWzFdICs9IHlcbiAgICAgICAgICB5ID0gc2VnWzFdXG4gICAgICAgIH1lbHNlIGlmKHMgPT0gJ2gnKXtcbiAgICAgICAgICBzZWdbMV0gKz0geFxuICAgICAgICAgIHggPSBzZWdbMV1cbiAgICAgICAgfWVsc2UgaWYocyA9PSAnYScpe1xuICAgICAgICAgIHNlZ1s2XSArPSB4XG4gICAgICAgICAgc2VnWzddICs9IHlcbiAgICAgICAgICB4ID0gc2VnWzZdXG4gICAgICAgICAgeSA9IHNlZ1s3XVxuICAgICAgICB9XG5cbiAgICAgIH1cblxuICAgICAgaWYoc2VnWzBdID09ICdNJyl7XG4gICAgICAgIHgwID0geFxuICAgICAgICB5MCA9IHlcbiAgICAgIH1cblxuICAgICAgaWYoc2VnWzBdID09ICdaJyl7XG4gICAgICAgIHggPSB4MFxuICAgICAgICB5ID0geTBcbiAgICAgIH1cblxuICAgICAgYXJyLnB1c2goc2VnKVxuXG4gICAgfXdoaWxlKGFycmF5Lmxlbmd0aClcblxuICAgIHJldHVybiBhcnJcblxuICB9XG4gIC8vIEdldCBib3VuZGluZyBib3ggb2YgcGF0aFxuLCBiYm94OiBmdW5jdGlvbigpIHtcbiAgICBTVkcucGFyc2VyLnBhdGguc2V0QXR0cmlidXRlKCdkJywgdGhpcy50b1N0cmluZygpKVxuXG4gICAgcmV0dXJuIFNWRy5wYXJzZXIucGF0aC5nZXRCQm94KClcbiAgfVxuXG59KVxuXG4vLyBNb2R1bGUgZm9yIHVuaXQgY29udmVydGlvbnNcblNWRy5OdW1iZXIgPSBTVkcuaW52ZW50KHtcbiAgLy8gSW5pdGlhbGl6ZVxuICBjcmVhdGU6IGZ1bmN0aW9uKHZhbHVlLCB1bml0KSB7XG4gICAgLy8gaW5pdGlhbGl6ZSBkZWZhdWx0c1xuICAgIHRoaXMudmFsdWUgPSAwXG4gICAgdGhpcy51bml0ICA9IHVuaXQgfHwgJydcblxuICAgIC8vIHBhcnNlIHZhbHVlXG4gICAgaWYgKHR5cGVvZiB2YWx1ZSA9PT0gJ251bWJlcicpIHtcbiAgICAgIC8vIGVuc3VyZSBhIHZhbGlkIG51bWVyaWMgdmFsdWVcbiAgICAgIHRoaXMudmFsdWUgPSBpc05hTih2YWx1ZSkgPyAwIDogIWlzRmluaXRlKHZhbHVlKSA/ICh2YWx1ZSA8IDAgPyAtMy40ZSszOCA6ICszLjRlKzM4KSA6IHZhbHVlXG5cbiAgICB9IGVsc2UgaWYgKHR5cGVvZiB2YWx1ZSA9PT0gJ3N0cmluZycpIHtcbiAgICAgIHVuaXQgPSB2YWx1ZS5tYXRjaChTVkcucmVnZXgubnVtYmVyQW5kVW5pdClcblxuICAgICAgaWYgKHVuaXQpIHtcbiAgICAgICAgLy8gbWFrZSB2YWx1ZSBudW1lcmljXG4gICAgICAgIHRoaXMudmFsdWUgPSBwYXJzZUZsb2F0KHVuaXRbMV0pXG5cbiAgICAgICAgLy8gbm9ybWFsaXplXG4gICAgICAgIGlmICh1bml0WzVdID09ICclJylcbiAgICAgICAgICB0aGlzLnZhbHVlIC89IDEwMFxuICAgICAgICBlbHNlIGlmICh1bml0WzVdID09ICdzJylcbiAgICAgICAgICB0aGlzLnZhbHVlICo9IDEwMDBcblxuICAgICAgICAvLyBzdG9yZSB1bml0XG4gICAgICAgIHRoaXMudW5pdCA9IHVuaXRbNV1cbiAgICAgIH1cblxuICAgIH0gZWxzZSB7XG4gICAgICBpZiAodmFsdWUgaW5zdGFuY2VvZiBTVkcuTnVtYmVyKSB7XG4gICAgICAgIHRoaXMudmFsdWUgPSB2YWx1ZS52YWx1ZU9mKClcbiAgICAgICAgdGhpcy51bml0ICA9IHZhbHVlLnVuaXRcbiAgICAgIH1cbiAgICB9XG5cbiAgfVxuICAvLyBBZGQgbWV0aG9kc1xuLCBleHRlbmQ6IHtcbiAgICAvLyBTdHJpbmdhbGl6ZVxuICAgIHRvU3RyaW5nOiBmdW5jdGlvbigpIHtcbiAgICAgIHJldHVybiAoXG4gICAgICAgIHRoaXMudW5pdCA9PSAnJScgP1xuICAgICAgICAgIH5+KHRoaXMudmFsdWUgKiAxZTgpIC8gMWU2OlxuICAgICAgICB0aGlzLnVuaXQgPT0gJ3MnID9cbiAgICAgICAgICB0aGlzLnZhbHVlIC8gMWUzIDpcbiAgICAgICAgICB0aGlzLnZhbHVlXG4gICAgICApICsgdGhpcy51bml0XG4gICAgfVxuICAsIHRvSlNPTjogZnVuY3Rpb24oKSB7XG4gICAgICByZXR1cm4gdGhpcy50b1N0cmluZygpXG4gICAgfVxuICAsIC8vIENvbnZlcnQgdG8gcHJpbWl0aXZlXG4gICAgdmFsdWVPZjogZnVuY3Rpb24oKSB7XG4gICAgICByZXR1cm4gdGhpcy52YWx1ZVxuICAgIH1cbiAgICAvLyBBZGQgbnVtYmVyXG4gICwgcGx1czogZnVuY3Rpb24obnVtYmVyKSB7XG4gICAgICByZXR1cm4gbmV3IFNWRy5OdW1iZXIodGhpcyArIG5ldyBTVkcuTnVtYmVyKG51bWJlciksIHRoaXMudW5pdClcbiAgICB9XG4gICAgLy8gU3VidHJhY3QgbnVtYmVyXG4gICwgbWludXM6IGZ1bmN0aW9uKG51bWJlcikge1xuICAgICAgcmV0dXJuIHRoaXMucGx1cygtbmV3IFNWRy5OdW1iZXIobnVtYmVyKSlcbiAgICB9XG4gICAgLy8gTXVsdGlwbHkgbnVtYmVyXG4gICwgdGltZXM6IGZ1bmN0aW9uKG51bWJlcikge1xuICAgICAgcmV0dXJuIG5ldyBTVkcuTnVtYmVyKHRoaXMgKiBuZXcgU1ZHLk51bWJlcihudW1iZXIpLCB0aGlzLnVuaXQpXG4gICAgfVxuICAgIC8vIERpdmlkZSBudW1iZXJcbiAgLCBkaXZpZGU6IGZ1bmN0aW9uKG51bWJlcikge1xuICAgICAgcmV0dXJuIG5ldyBTVkcuTnVtYmVyKHRoaXMgLyBuZXcgU1ZHLk51bWJlcihudW1iZXIpLCB0aGlzLnVuaXQpXG4gICAgfVxuICAgIC8vIENvbnZlcnQgdG8gZGlmZmVyZW50IHVuaXRcbiAgLCB0bzogZnVuY3Rpb24odW5pdCkge1xuICAgICAgdmFyIG51bWJlciA9IG5ldyBTVkcuTnVtYmVyKHRoaXMpXG5cbiAgICAgIGlmICh0eXBlb2YgdW5pdCA9PT0gJ3N0cmluZycpXG4gICAgICAgIG51bWJlci51bml0ID0gdW5pdFxuXG4gICAgICByZXR1cm4gbnVtYmVyXG4gICAgfVxuICAgIC8vIE1ha2UgbnVtYmVyIG1vcnBoYWJsZVxuICAsIG1vcnBoOiBmdW5jdGlvbihudW1iZXIpIHtcbiAgICAgIHRoaXMuZGVzdGluYXRpb24gPSBuZXcgU1ZHLk51bWJlcihudW1iZXIpXG5cbiAgICAgIHJldHVybiB0aGlzXG4gICAgfVxuICAgIC8vIEdldCBtb3JwaGVkIG51bWJlciBhdCBnaXZlbiBwb3NpdGlvblxuICAsIGF0OiBmdW5jdGlvbihwb3MpIHtcbiAgICAgIC8vIE1ha2Ugc3VyZSBhIGRlc3RpbmF0aW9uIGlzIGRlZmluZWRcbiAgICAgIGlmICghdGhpcy5kZXN0aW5hdGlvbikgcmV0dXJuIHRoaXNcblxuICAgICAgLy8gR2VuZXJhdGUgbmV3IG1vcnBoZWQgbnVtYmVyXG4gICAgICByZXR1cm4gbmV3IFNWRy5OdW1iZXIodGhpcy5kZXN0aW5hdGlvbilcbiAgICAgICAgICAubWludXModGhpcylcbiAgICAgICAgICAudGltZXMocG9zKVxuICAgICAgICAgIC5wbHVzKHRoaXMpXG4gICAgfVxuXG4gIH1cbn0pXG5cblNWRy5FbGVtZW50ID0gU1ZHLmludmVudCh7XG4gIC8vIEluaXRpYWxpemUgbm9kZVxuICBjcmVhdGU6IGZ1bmN0aW9uKG5vZGUpIHtcbiAgICAvLyBtYWtlIHN0cm9rZSB2YWx1ZSBhY2Nlc3NpYmxlIGR5bmFtaWNhbGx5XG4gICAgdGhpcy5fc3Ryb2tlID0gU1ZHLmRlZmF1bHRzLmF0dHJzLnN0cm9rZVxuXG4gICAgLy8gaW5pdGlhbGl6ZSBkYXRhIG9iamVjdFxuICAgIHRoaXMuZG9tID0ge31cblxuICAgIC8vIGNyZWF0ZSBjaXJjdWxhciByZWZlcmVuY2VcbiAgICBpZiAodGhpcy5ub2RlID0gbm9kZSkge1xuICAgICAgdGhpcy50eXBlID0gbm9kZS5ub2RlTmFtZVxuICAgICAgdGhpcy5ub2RlLmluc3RhbmNlID0gdGhpc1xuXG4gICAgICAvLyBzdG9yZSBjdXJyZW50IGF0dHJpYnV0ZSB2YWx1ZVxuICAgICAgdGhpcy5fc3Ryb2tlID0gbm9kZS5nZXRBdHRyaWJ1dGUoJ3N0cm9rZScpIHx8IHRoaXMuX3N0cm9rZVxuICAgIH1cbiAgfVxuXG4gIC8vIEFkZCBjbGFzcyBtZXRob2RzXG4sIGV4dGVuZDoge1xuICAgIC8vIE1vdmUgb3ZlciB4LWF4aXNcbiAgICB4OiBmdW5jdGlvbih4KSB7XG4gICAgICByZXR1cm4gdGhpcy5hdHRyKCd4JywgeClcbiAgICB9XG4gICAgLy8gTW92ZSBvdmVyIHktYXhpc1xuICAsIHk6IGZ1bmN0aW9uKHkpIHtcbiAgICAgIHJldHVybiB0aGlzLmF0dHIoJ3knLCB5KVxuICAgIH1cbiAgICAvLyBNb3ZlIGJ5IGNlbnRlciBvdmVyIHgtYXhpc1xuICAsIGN4OiBmdW5jdGlvbih4KSB7XG4gICAgICByZXR1cm4geCA9PSBudWxsID8gdGhpcy54KCkgKyB0aGlzLndpZHRoKCkgLyAyIDogdGhpcy54KHggLSB0aGlzLndpZHRoKCkgLyAyKVxuICAgIH1cbiAgICAvLyBNb3ZlIGJ5IGNlbnRlciBvdmVyIHktYXhpc1xuICAsIGN5OiBmdW5jdGlvbih5KSB7XG4gICAgICByZXR1cm4geSA9PSBudWxsID8gdGhpcy55KCkgKyB0aGlzLmhlaWdodCgpIC8gMiA6IHRoaXMueSh5IC0gdGhpcy5oZWlnaHQoKSAvIDIpXG4gICAgfVxuICAgIC8vIE1vdmUgZWxlbWVudCB0byBnaXZlbiB4IGFuZCB5IHZhbHVlc1xuICAsIG1vdmU6IGZ1bmN0aW9uKHgsIHkpIHtcbiAgICAgIHJldHVybiB0aGlzLngoeCkueSh5KVxuICAgIH1cbiAgICAvLyBNb3ZlIGVsZW1lbnQgYnkgaXRzIGNlbnRlclxuICAsIGNlbnRlcjogZnVuY3Rpb24oeCwgeSkge1xuICAgICAgcmV0dXJuIHRoaXMuY3goeCkuY3koeSlcbiAgICB9XG4gICAgLy8gU2V0IHdpZHRoIG9mIGVsZW1lbnRcbiAgLCB3aWR0aDogZnVuY3Rpb24od2lkdGgpIHtcbiAgICAgIHJldHVybiB0aGlzLmF0dHIoJ3dpZHRoJywgd2lkdGgpXG4gICAgfVxuICAgIC8vIFNldCBoZWlnaHQgb2YgZWxlbWVudFxuICAsIGhlaWdodDogZnVuY3Rpb24oaGVpZ2h0KSB7XG4gICAgICByZXR1cm4gdGhpcy5hdHRyKCdoZWlnaHQnLCBoZWlnaHQpXG4gICAgfVxuICAgIC8vIFNldCBlbGVtZW50IHNpemUgdG8gZ2l2ZW4gd2lkdGggYW5kIGhlaWdodFxuICAsIHNpemU6IGZ1bmN0aW9uKHdpZHRoLCBoZWlnaHQpIHtcbiAgICAgIHZhciBwID0gcHJvcG9ydGlvbmFsU2l6ZSh0aGlzLCB3aWR0aCwgaGVpZ2h0KVxuXG4gICAgICByZXR1cm4gdGhpc1xuICAgICAgICAud2lkdGgobmV3IFNWRy5OdW1iZXIocC53aWR0aCkpXG4gICAgICAgIC5oZWlnaHQobmV3IFNWRy5OdW1iZXIocC5oZWlnaHQpKVxuICAgIH1cbiAgICAvLyBDbG9uZSBlbGVtZW50XG4gICwgY2xvbmU6IGZ1bmN0aW9uKHBhcmVudCkge1xuICAgICAgLy8gY2xvbmUgZWxlbWVudCBhbmQgYXNzaWduIG5ldyBpZFxuICAgICAgdmFyIGNsb25lID0gYXNzaWduTmV3SWQodGhpcy5ub2RlLmNsb25lTm9kZSh0cnVlKSlcblxuICAgICAgLy8gaW5zZXJ0IHRoZSBjbG9uZSBpbiB0aGUgZ2l2ZW4gcGFyZW50IG9yIGFmdGVyIG15c2VsZlxuICAgICAgaWYocGFyZW50KSBwYXJlbnQuYWRkKGNsb25lKVxuICAgICAgZWxzZSB0aGlzLmFmdGVyKGNsb25lKVxuXG4gICAgICByZXR1cm4gY2xvbmVcbiAgICB9XG4gICAgLy8gUmVtb3ZlIGVsZW1lbnRcbiAgLCByZW1vdmU6IGZ1bmN0aW9uKCkge1xuICAgICAgaWYgKHRoaXMucGFyZW50KCkpXG4gICAgICAgIHRoaXMucGFyZW50KCkucmVtb3ZlRWxlbWVudCh0aGlzKVxuXG4gICAgICByZXR1cm4gdGhpc1xuICAgIH1cbiAgICAvLyBSZXBsYWNlIGVsZW1lbnRcbiAgLCByZXBsYWNlOiBmdW5jdGlvbihlbGVtZW50KSB7XG4gICAgICB0aGlzLmFmdGVyKGVsZW1lbnQpLnJlbW92ZSgpXG5cbiAgICAgIHJldHVybiBlbGVtZW50XG4gICAgfVxuICAgIC8vIEFkZCBlbGVtZW50IHRvIGdpdmVuIGNvbnRhaW5lciBhbmQgcmV0dXJuIHNlbGZcbiAgLCBhZGRUbzogZnVuY3Rpb24ocGFyZW50KSB7XG4gICAgICByZXR1cm4gcGFyZW50LnB1dCh0aGlzKVxuICAgIH1cbiAgICAvLyBBZGQgZWxlbWVudCB0byBnaXZlbiBjb250YWluZXIgYW5kIHJldHVybiBjb250YWluZXJcbiAgLCBwdXRJbjogZnVuY3Rpb24ocGFyZW50KSB7XG4gICAgICByZXR1cm4gcGFyZW50LmFkZCh0aGlzKVxuICAgIH1cbiAgICAvLyBHZXQgLyBzZXQgaWRcbiAgLCBpZDogZnVuY3Rpb24oaWQpIHtcbiAgICAgIHJldHVybiB0aGlzLmF0dHIoJ2lkJywgaWQpXG4gICAgfVxuICAgIC8vIENoZWNrcyB3aGV0aGVyIHRoZSBnaXZlbiBwb2ludCBpbnNpZGUgdGhlIGJvdW5kaW5nIGJveCBvZiB0aGUgZWxlbWVudFxuICAsIGluc2lkZTogZnVuY3Rpb24oeCwgeSkge1xuICAgICAgdmFyIGJveCA9IHRoaXMuYmJveCgpXG5cbiAgICAgIHJldHVybiB4ID4gYm94LnhcbiAgICAgICAgICAmJiB5ID4gYm94LnlcbiAgICAgICAgICAmJiB4IDwgYm94LnggKyBib3gud2lkdGhcbiAgICAgICAgICAmJiB5IDwgYm94LnkgKyBib3guaGVpZ2h0XG4gICAgfVxuICAgIC8vIFNob3cgZWxlbWVudFxuICAsIHNob3c6IGZ1bmN0aW9uKCkge1xuICAgICAgcmV0dXJuIHRoaXMuc3R5bGUoJ2Rpc3BsYXknLCAnJylcbiAgICB9XG4gICAgLy8gSGlkZSBlbGVtZW50XG4gICwgaGlkZTogZnVuY3Rpb24oKSB7XG4gICAgICByZXR1cm4gdGhpcy5zdHlsZSgnZGlzcGxheScsICdub25lJylcbiAgICB9XG4gICAgLy8gSXMgZWxlbWVudCB2aXNpYmxlP1xuICAsIHZpc2libGU6IGZ1bmN0aW9uKCkge1xuICAgICAgcmV0dXJuIHRoaXMuc3R5bGUoJ2Rpc3BsYXknKSAhPSAnbm9uZSdcbiAgICB9XG4gICAgLy8gUmV0dXJuIGlkIG9uIHN0cmluZyBjb252ZXJzaW9uXG4gICwgdG9TdHJpbmc6IGZ1bmN0aW9uKCkge1xuICAgICAgcmV0dXJuIHRoaXMuYXR0cignaWQnKVxuICAgIH1cbiAgICAvLyBSZXR1cm4gYXJyYXkgb2YgY2xhc3NlcyBvbiB0aGUgbm9kZVxuICAsIGNsYXNzZXM6IGZ1bmN0aW9uKCkge1xuICAgICAgdmFyIGF0dHIgPSB0aGlzLmF0dHIoJ2NsYXNzJylcblxuICAgICAgcmV0dXJuIGF0dHIgPT0gbnVsbCA/IFtdIDogYXR0ci50cmltKCkuc3BsaXQoL1xccysvKVxuICAgIH1cbiAgICAvLyBSZXR1cm4gdHJ1ZSBpZiBjbGFzcyBleGlzdHMgb24gdGhlIG5vZGUsIGZhbHNlIG90aGVyd2lzZVxuICAsIGhhc0NsYXNzOiBmdW5jdGlvbihuYW1lKSB7XG4gICAgICByZXR1cm4gdGhpcy5jbGFzc2VzKCkuaW5kZXhPZihuYW1lKSAhPSAtMVxuICAgIH1cbiAgICAvLyBBZGQgY2xhc3MgdG8gdGhlIG5vZGVcbiAgLCBhZGRDbGFzczogZnVuY3Rpb24obmFtZSkge1xuICAgICAgaWYgKCF0aGlzLmhhc0NsYXNzKG5hbWUpKSB7XG4gICAgICAgIHZhciBhcnJheSA9IHRoaXMuY2xhc3NlcygpXG4gICAgICAgIGFycmF5LnB1c2gobmFtZSlcbiAgICAgICAgdGhpcy5hdHRyKCdjbGFzcycsIGFycmF5LmpvaW4oJyAnKSlcbiAgICAgIH1cblxuICAgICAgcmV0dXJuIHRoaXNcbiAgICB9XG4gICAgLy8gUmVtb3ZlIGNsYXNzIGZyb20gdGhlIG5vZGVcbiAgLCByZW1vdmVDbGFzczogZnVuY3Rpb24obmFtZSkge1xuICAgICAgaWYgKHRoaXMuaGFzQ2xhc3MobmFtZSkpIHtcbiAgICAgICAgdGhpcy5hdHRyKCdjbGFzcycsIHRoaXMuY2xhc3NlcygpLmZpbHRlcihmdW5jdGlvbihjKSB7XG4gICAgICAgICAgcmV0dXJuIGMgIT0gbmFtZVxuICAgICAgICB9KS5qb2luKCcgJykpXG4gICAgICB9XG5cbiAgICAgIHJldHVybiB0aGlzXG4gICAgfVxuICAgIC8vIFRvZ2dsZSB0aGUgcHJlc2VuY2Ugb2YgYSBjbGFzcyBvbiB0aGUgbm9kZVxuICAsIHRvZ2dsZUNsYXNzOiBmdW5jdGlvbihuYW1lKSB7XG4gICAgICByZXR1cm4gdGhpcy5oYXNDbGFzcyhuYW1lKSA/IHRoaXMucmVtb3ZlQ2xhc3MobmFtZSkgOiB0aGlzLmFkZENsYXNzKG5hbWUpXG4gICAgfVxuICAgIC8vIEdldCByZWZlcmVuY2VkIGVsZW1lbnQgZm9ybSBhdHRyaWJ1dGUgdmFsdWVcbiAgLCByZWZlcmVuY2U6IGZ1bmN0aW9uKGF0dHIpIHtcbiAgICAgIHJldHVybiBTVkcuZ2V0KHRoaXMuYXR0cihhdHRyKSlcbiAgICB9XG4gICAgLy8gUmV0dXJucyB0aGUgcGFyZW50IGVsZW1lbnQgaW5zdGFuY2VcbiAgLCBwYXJlbnQ6IGZ1bmN0aW9uKHR5cGUpIHtcbiAgICAgIHZhciBwYXJlbnQgPSB0aGlzXG5cbiAgICAgIC8vIGNoZWNrIGZvciBwYXJlbnRcbiAgICAgIGlmKCFwYXJlbnQubm9kZS5wYXJlbnROb2RlKSByZXR1cm4gbnVsbFxuXG4gICAgICAvLyBnZXQgcGFyZW50IGVsZW1lbnRcbiAgICAgIHBhcmVudCA9IFNWRy5hZG9wdChwYXJlbnQubm9kZS5wYXJlbnROb2RlKVxuXG4gICAgICBpZighdHlwZSkgcmV0dXJuIHBhcmVudFxuXG4gICAgICAvLyBsb29wIHRyb3VnaCBhbmNlc3RvcnMgaWYgdHlwZSBpcyBnaXZlblxuICAgICAgd2hpbGUocGFyZW50ICYmIHBhcmVudC5ub2RlIGluc3RhbmNlb2YgU1ZHRWxlbWVudCl7XG4gICAgICAgIGlmKHR5cGVvZiB0eXBlID09PSAnc3RyaW5nJyA/IHBhcmVudC5tYXRjaGVzKHR5cGUpIDogcGFyZW50IGluc3RhbmNlb2YgdHlwZSkgcmV0dXJuIHBhcmVudFxuICAgICAgICBwYXJlbnQgPSBTVkcuYWRvcHQocGFyZW50Lm5vZGUucGFyZW50Tm9kZSlcbiAgICAgIH1cbiAgICB9XG4gICAgLy8gR2V0IHBhcmVudCBkb2N1bWVudFxuICAsIGRvYzogZnVuY3Rpb24oKSB7XG4gICAgICByZXR1cm4gdGhpcyBpbnN0YW5jZW9mIFNWRy5Eb2MgPyB0aGlzIDogdGhpcy5wYXJlbnQoU1ZHLkRvYylcbiAgICB9XG4gICAgLy8gcmV0dXJuIGFycmF5IG9mIGFsbCBhbmNlc3RvcnMgb2YgZ2l2ZW4gdHlwZSB1cCB0byB0aGUgcm9vdCBzdmdcbiAgLCBwYXJlbnRzOiBmdW5jdGlvbih0eXBlKSB7XG4gICAgICB2YXIgcGFyZW50cyA9IFtdLCBwYXJlbnQgPSB0aGlzXG5cbiAgICAgIGRve1xuICAgICAgICBwYXJlbnQgPSBwYXJlbnQucGFyZW50KHR5cGUpXG4gICAgICAgIGlmKCFwYXJlbnQgfHwgIXBhcmVudC5ub2RlKSBicmVha1xuXG4gICAgICAgIHBhcmVudHMucHVzaChwYXJlbnQpXG4gICAgICB9IHdoaWxlKHBhcmVudC5wYXJlbnQpXG5cbiAgICAgIHJldHVybiBwYXJlbnRzXG4gICAgfVxuICAgIC8vIG1hdGNoZXMgdGhlIGVsZW1lbnQgdnMgYSBjc3Mgc2VsZWN0b3JcbiAgLCBtYXRjaGVzOiBmdW5jdGlvbihzZWxlY3Rvcil7XG4gICAgICByZXR1cm4gbWF0Y2hlcyh0aGlzLm5vZGUsIHNlbGVjdG9yKVxuICAgIH1cbiAgICAvLyBSZXR1cm5zIHRoZSBzdmcgbm9kZSB0byBjYWxsIG5hdGl2ZSBzdmcgbWV0aG9kcyBvbiBpdFxuICAsIG5hdGl2ZTogZnVuY3Rpb24oKSB7XG4gICAgICByZXR1cm4gdGhpcy5ub2RlXG4gICAgfVxuICAgIC8vIEltcG9ydCByYXcgc3ZnXG4gICwgc3ZnOiBmdW5jdGlvbihzdmcpIHtcbiAgICAgIC8vIGNyZWF0ZSB0ZW1wb3JhcnkgaG9sZGVyXG4gICAgICB2YXIgd2VsbCA9IGRvY3VtZW50LmNyZWF0ZUVsZW1lbnQoJ3N2ZycpXG5cbiAgICAgIC8vIGFjdCBhcyBhIHNldHRlciBpZiBzdmcgaXMgZ2l2ZW5cbiAgICAgIGlmIChzdmcgJiYgdGhpcyBpbnN0YW5jZW9mIFNWRy5QYXJlbnQpIHtcbiAgICAgICAgLy8gZHVtcCByYXcgc3ZnXG4gICAgICAgIHdlbGwuaW5uZXJIVE1MID0gJzxzdmc+JyArIHN2Zy5yZXBsYWNlKC9cXG4vLCAnJykucmVwbGFjZSgvPChcXHcrKShbXjxdKz8pXFwvPi9nLCAnPCQxJDI+PC8kMT4nKSArICc8L3N2Zz4nXG5cbiAgICAgICAgLy8gdHJhbnNwbGFudCBub2Rlc1xuICAgICAgICBmb3IgKHZhciBpID0gMCwgaWwgPSB3ZWxsLmZpcnN0Q2hpbGQuY2hpbGROb2Rlcy5sZW5ndGg7IGkgPCBpbDsgaSsrKVxuICAgICAgICAgIHRoaXMubm9kZS5hcHBlbmRDaGlsZCh3ZWxsLmZpcnN0Q2hpbGQuZmlyc3RDaGlsZClcblxuICAgICAgLy8gb3RoZXJ3aXNlIGFjdCBhcyBhIGdldHRlclxuICAgICAgfSBlbHNlIHtcbiAgICAgICAgLy8gY3JlYXRlIGEgd3JhcHBpbmcgc3ZnIGVsZW1lbnQgaW4gY2FzZSBvZiBwYXJ0aWFsIGNvbnRlbnRcbiAgICAgICAgd2VsbC5hcHBlbmRDaGlsZChzdmcgPSBkb2N1bWVudC5jcmVhdGVFbGVtZW50KCdzdmcnKSlcblxuICAgICAgICAvLyB3cml0ZSBzdmdqcyBkYXRhIHRvIHRoZSBkb21cbiAgICAgICAgdGhpcy53cml0ZURhdGFUb0RvbSgpXG5cbiAgICAgICAgLy8gaW5zZXJ0IGEgY29weSBvZiB0aGlzIG5vZGVcbiAgICAgICAgc3ZnLmFwcGVuZENoaWxkKHRoaXMubm9kZS5jbG9uZU5vZGUodHJ1ZSkpXG5cbiAgICAgICAgLy8gcmV0dXJuIHRhcmdldCBlbGVtZW50XG4gICAgICAgIHJldHVybiB3ZWxsLmlubmVySFRNTC5yZXBsYWNlKC9ePHN2Zz4vLCAnJykucmVwbGFjZSgvPFxcL3N2Zz4kLywgJycpXG4gICAgICB9XG5cbiAgICAgIHJldHVybiB0aGlzXG4gICAgfVxuICAvLyB3cml0ZSBzdmdqcyBkYXRhIHRvIHRoZSBkb21cbiAgLCB3cml0ZURhdGFUb0RvbTogZnVuY3Rpb24oKSB7XG5cbiAgICAgIC8vIGR1bXAgdmFyaWFibGVzIHJlY3Vyc2l2ZWx5XG4gICAgICBpZih0aGlzLmVhY2ggfHwgdGhpcy5saW5lcyl7XG4gICAgICAgIHZhciBmbiA9IHRoaXMuZWFjaCA/IHRoaXMgOiB0aGlzLmxpbmVzKCk7XG4gICAgICAgIGZuLmVhY2goZnVuY3Rpb24oKXtcbiAgICAgICAgICB0aGlzLndyaXRlRGF0YVRvRG9tKClcbiAgICAgICAgfSlcbiAgICAgIH1cblxuICAgICAgLy8gcmVtb3ZlIHByZXZpb3VzbHkgc2V0IGRhdGFcbiAgICAgIHRoaXMubm9kZS5yZW1vdmVBdHRyaWJ1dGUoJ3N2Z2pzOmRhdGEnKVxuXG4gICAgICBpZihPYmplY3Qua2V5cyh0aGlzLmRvbSkubGVuZ3RoKVxuICAgICAgICB0aGlzLm5vZGUuc2V0QXR0cmlidXRlKCdzdmdqczpkYXRhJywgSlNPTi5zdHJpbmdpZnkodGhpcy5kb20pKSAvLyBzZWUgIzQyOFxuXG4gICAgICByZXR1cm4gdGhpc1xuICAgIH1cbiAgLy8gc2V0IGdpdmVuIGRhdGEgdG8gdGhlIGVsZW1lbnRzIGRhdGEgcHJvcGVydHlcbiAgLCBzZXREYXRhOiBmdW5jdGlvbihvKXtcbiAgICAgIHRoaXMuZG9tID0gb1xuICAgICAgcmV0dXJuIHRoaXNcbiAgICB9XG4gICwgaXM6IGZ1bmN0aW9uKG9iail7XG4gICAgICByZXR1cm4gaXModGhpcywgb2JqKVxuICAgIH1cbiAgfVxufSlcblxuU1ZHLmVhc2luZyA9IHtcbiAgJy0nOiBmdW5jdGlvbihwb3Mpe3JldHVybiBwb3N9XG4sICc8Pic6ZnVuY3Rpb24ocG9zKXtyZXR1cm4gLU1hdGguY29zKHBvcyAqIE1hdGguUEkpIC8gMiArIDAuNX1cbiwgJz4nOiBmdW5jdGlvbihwb3Mpe3JldHVybiAgTWF0aC5zaW4ocG9zICogTWF0aC5QSSAvIDIpfVxuLCAnPCc6IGZ1bmN0aW9uKHBvcyl7cmV0dXJuIC1NYXRoLmNvcyhwb3MgKiBNYXRoLlBJIC8gMikgKyAxfVxufVxuXG5TVkcubW9ycGggPSBmdW5jdGlvbihwb3Mpe1xuICByZXR1cm4gZnVuY3Rpb24oZnJvbSwgdG8pIHtcbiAgICByZXR1cm4gbmV3IFNWRy5Nb3JwaE9iaihmcm9tLCB0bykuYXQocG9zKVxuICB9XG59XG5cblNWRy5TaXR1YXRpb24gPSBTVkcuaW52ZW50KHtcblxuICBjcmVhdGU6IGZ1bmN0aW9uKG8pe1xuICAgIHRoaXMuaW5pdCA9IGZhbHNlXG4gICAgdGhpcy5yZXZlcnNlZCA9IGZhbHNlXG4gICAgdGhpcy5yZXZlcnNpbmcgPSBmYWxzZVxuXG4gICAgdGhpcy5kdXJhdGlvbiA9IG5ldyBTVkcuTnVtYmVyKG8uZHVyYXRpb24pLnZhbHVlT2YoKVxuICAgIHRoaXMuZGVsYXkgPSBuZXcgU1ZHLk51bWJlcihvLmRlbGF5KS52YWx1ZU9mKClcblxuICAgIHRoaXMuc3RhcnQgPSArbmV3IERhdGUoKSArIHRoaXMuZGVsYXlcbiAgICB0aGlzLmZpbmlzaCA9IHRoaXMuc3RhcnQgKyB0aGlzLmR1cmF0aW9uXG4gICAgdGhpcy5lYXNlID0gby5lYXNlXG5cbiAgICAvLyB0aGlzLmxvb3AgaXMgaW5jcmVtZW50ZWQgZnJvbSAwIHRvIHRoaXMubG9vcHNcbiAgICAvLyBpdCBpcyBhbHNvIGluY3JlbWVudGVkIHdoZW4gaW4gYW4gaW5maW5pdGUgbG9vcCAod2hlbiB0aGlzLmxvb3BzIGlzIHRydWUpXG4gICAgdGhpcy5sb29wID0gMFxuICAgIHRoaXMubG9vcHMgPSBmYWxzZVxuXG4gICAgdGhpcy5hbmltYXRpb25zID0ge1xuICAgICAgLy8gZnVuY3Rpb25Ub0NhbGw6IFtsaXN0IG9mIG1vcnBoYWJsZSBvYmplY3RzXVxuICAgICAgLy8gZS5nLiBtb3ZlOiBbU1ZHLk51bWJlciwgU1ZHLk51bWJlcl1cbiAgICB9XG5cbiAgICB0aGlzLmF0dHJzID0ge1xuICAgICAgLy8gaG9sZHMgYWxsIGF0dHJpYnV0ZXMgd2hpY2ggYXJlIG5vdCByZXByZXNlbnRlZCBmcm9tIGEgZnVuY3Rpb24gc3ZnLmpzIHByb3ZpZGVzXG4gICAgICAvLyBlLmcuIHNvbWVBdHRyOiBTVkcuTnVtYmVyXG4gICAgfVxuXG4gICAgdGhpcy5zdHlsZXMgPSB7XG4gICAgICAvLyBob2xkcyBhbGwgc3R5bGVzIHdoaWNoIHNob3VsZCBiZSBhbmltYXRlZFxuICAgICAgLy8gZS5nLiBmaWxsLWNvbG9yOiBTVkcuQ29sb3JcbiAgICB9XG5cbiAgICB0aGlzLnRyYW5zZm9ybXMgPSBbXG4gICAgICAvLyBob2xkcyBhbGwgdHJhbnNmb3JtYXRpb25zIGFzIHRyYW5zZm9ybWF0aW9uIG9iamVjdHNcbiAgICAgIC8vIGUuZy4gW1NWRy5Sb3RhdGUsIFNWRy5UcmFuc2xhdGUsIFNWRy5NYXRyaXhdXG4gICAgXVxuXG4gICAgdGhpcy5vbmNlID0ge1xuICAgICAgLy8gZnVuY3Rpb25zIHRvIGZpcmUgYXQgYSBzcGVjaWZpYyBwb3NpdGlvblxuICAgICAgLy8gZS5nLiBcIjAuNVwiOiBmdW5jdGlvbiBmb28oKXt9XG4gICAgfVxuXG4gIH1cblxufSlcblxuXG5TVkcuRlggPSBTVkcuaW52ZW50KHtcblxuICBjcmVhdGU6IGZ1bmN0aW9uKGVsZW1lbnQpIHtcbiAgICB0aGlzLl90YXJnZXQgPSBlbGVtZW50XG4gICAgdGhpcy5zaXR1YXRpb25zID0gW11cbiAgICB0aGlzLmFjdGl2ZSA9IGZhbHNlXG4gICAgdGhpcy5zaXR1YXRpb24gPSBudWxsXG4gICAgdGhpcy5wYXVzZWQgPSBmYWxzZVxuICAgIHRoaXMubGFzdFBvcyA9IDBcbiAgICB0aGlzLnBvcyA9IDBcbiAgICAvLyBUaGUgYWJzb2x1dGUgcG9zaXRpb24gb2YgYW4gYW5pbWF0aW9uIGlzIGl0cyBwb3NpdGlvbiBpbiB0aGUgY29udGV4dCBvZiBpdHMgY29tcGxldGUgZHVyYXRpb24gKGluY2x1ZGluZyBkZWxheSBhbmQgbG9vcHMpXG4gICAgLy8gV2hlbiBwZXJmb3JtaW5nIGEgZGVsYXksIGFic1BvcyBpcyBiZWxvdyAwIGFuZCB3aGVuIHBlcmZvcm1pbmcgYSBsb29wLCBpdHMgdmFsdWUgaXMgYWJvdmUgMVxuICAgIHRoaXMuYWJzUG9zID0gMFxuICAgIHRoaXMuX3NwZWVkID0gMVxuICB9XG5cbiwgZXh0ZW5kOiB7XG5cbiAgICAvKipcbiAgICAgKiBzZXRzIG9yIHJldHVybnMgdGhlIHRhcmdldCBvZiB0aGlzIGFuaW1hdGlvblxuICAgICAqIEBwYXJhbSBvIG9iamVjdCB8fCBudW1iZXIgSW4gY2FzZSBvZiBPYmplY3QgaXQgaG9sZHMgYWxsIHBhcmFtZXRlcnMuIEluIGNhc2Ugb2YgbnVtYmVyIGl0cyB0aGUgZHVyYXRpb24gb2YgdGhlIGFuaW1hdGlvblxuICAgICAqIEBwYXJhbSBlYXNlIGZ1bmN0aW9uIHx8IHN0cmluZyBGdW5jdGlvbiB3aGljaCBzaG91bGQgYmUgdXNlZCBmb3IgZWFzaW5nIG9yIGVhc2luZyBrZXl3b3JkXG4gICAgICogQHBhcmFtIGRlbGF5IE51bWJlciBpbmRpY2F0aW5nIHRoZSBkZWxheSBiZWZvcmUgdGhlIGFuaW1hdGlvbiBzdGFydHNcbiAgICAgKiBAcmV0dXJuIHRhcmdldCB8fCB0aGlzXG4gICAgICovXG4gICAgYW5pbWF0ZTogZnVuY3Rpb24obywgZWFzZSwgZGVsYXkpe1xuXG4gICAgICBpZih0eXBlb2YgbyA9PSAnb2JqZWN0Jyl7XG4gICAgICAgIGVhc2UgPSBvLmVhc2VcbiAgICAgICAgZGVsYXkgPSBvLmRlbGF5XG4gICAgICAgIG8gPSBvLmR1cmF0aW9uXG4gICAgICB9XG5cbiAgICAgIHZhciBzaXR1YXRpb24gPSBuZXcgU1ZHLlNpdHVhdGlvbih7XG4gICAgICAgIGR1cmF0aW9uOiBvIHx8IDEwMDAsXG4gICAgICAgIGRlbGF5OiBkZWxheSB8fCAwLFxuICAgICAgICBlYXNlOiBTVkcuZWFzaW5nW2Vhc2UgfHwgJy0nXSB8fCBlYXNlXG4gICAgICB9KVxuXG4gICAgICB0aGlzLnF1ZXVlKHNpdHVhdGlvbilcblxuICAgICAgcmV0dXJuIHRoaXNcbiAgICB9XG5cbiAgICAvKipcbiAgICAgKiBzZXRzIGEgZGVsYXkgYmVmb3JlIHRoZSBuZXh0IGVsZW1lbnQgb2YgdGhlIHF1ZXVlIGlzIGNhbGxlZFxuICAgICAqIEBwYXJhbSBkZWxheSBEdXJhdGlvbiBvZiBkZWxheSBpbiBtaWxsaXNlY29uZHNcbiAgICAgKiBAcmV0dXJuIHRoaXMudGFyZ2V0KClcbiAgICAgKi9cbiAgLCBkZWxheTogZnVuY3Rpb24oZGVsYXkpe1xuICAgICAgLy8gVGhlIGRlbGF5IGlzIHBlcmZvcm1lZCBieSBhbiBlbXB0eSBzaXR1YXRpb24gd2l0aCBpdHMgZHVyYXRpb25cbiAgICAgIC8vIGF0dHJpYnV0ZSBzZXQgdG8gdGhlIGR1cmF0aW9uIG9mIHRoZSBkZWxheVxuICAgICAgdmFyIHNpdHVhdGlvbiA9IG5ldyBTVkcuU2l0dWF0aW9uKHtcbiAgICAgICAgZHVyYXRpb246IGRlbGF5LFxuICAgICAgICBkZWxheTogMCxcbiAgICAgICAgZWFzZTogU1ZHLmVhc2luZ1snLSddXG4gICAgICB9KVxuXG4gICAgICByZXR1cm4gdGhpcy5xdWV1ZShzaXR1YXRpb24pXG4gICAgfVxuXG4gICAgLyoqXG4gICAgICogc2V0cyBvciByZXR1cm5zIHRoZSB0YXJnZXQgb2YgdGhpcyBhbmltYXRpb25cbiAgICAgKiBAcGFyYW0gbnVsbCB8fCB0YXJnZXQgU1ZHLkVsZW1lbnQgd2hpY2ggc2hvdWxkIGJlIHNldCBhcyBuZXcgdGFyZ2V0XG4gICAgICogQHJldHVybiB0YXJnZXQgfHwgdGhpc1xuICAgICAqL1xuICAsIHRhcmdldDogZnVuY3Rpb24odGFyZ2V0KXtcbiAgICAgIGlmKHRhcmdldCAmJiB0YXJnZXQgaW5zdGFuY2VvZiBTVkcuRWxlbWVudCl7XG4gICAgICAgIHRoaXMuX3RhcmdldCA9IHRhcmdldFxuICAgICAgICByZXR1cm4gdGhpc1xuICAgICAgfVxuXG4gICAgICByZXR1cm4gdGhpcy5fdGFyZ2V0XG4gICAgfVxuXG4gICAgLy8gcmV0dXJucyB0aGUgYWJzb2x1dGUgcG9zaXRpb24gYXQgYSBnaXZlbiB0aW1lXG4gICwgdGltZVRvQWJzUG9zOiBmdW5jdGlvbih0aW1lc3RhbXApe1xuICAgICAgcmV0dXJuICh0aW1lc3RhbXAgLSB0aGlzLnNpdHVhdGlvbi5zdGFydCkgLyAodGhpcy5zaXR1YXRpb24uZHVyYXRpb24vdGhpcy5fc3BlZWQpXG4gICAgfVxuXG4gICAgLy8gcmV0dXJucyB0aGUgdGltZXN0YW1wIGZyb20gYSBnaXZlbiBhYnNvbHV0ZSBwb3NpdG9uXG4gICwgYWJzUG9zVG9UaW1lOiBmdW5jdGlvbihhYnNQb3Mpe1xuICAgICAgcmV0dXJuIHRoaXMuc2l0dWF0aW9uLmR1cmF0aW9uL3RoaXMuX3NwZWVkICogYWJzUG9zICsgdGhpcy5zaXR1YXRpb24uc3RhcnRcbiAgICB9XG5cbiAgICAvLyBzdGFydHMgdGhlIGFuaW1hdGlvbmxvb3BcbiAgLCBzdGFydEFuaW1GcmFtZTogZnVuY3Rpb24oKXtcbiAgICAgIHRoaXMuc3RvcEFuaW1GcmFtZSgpXG4gICAgICB0aGlzLmFuaW1hdGlvbkZyYW1lID0gcmVxdWVzdEFuaW1hdGlvbkZyYW1lKGZ1bmN0aW9uKCl7IHRoaXMuc3RlcCgpIH0uYmluZCh0aGlzKSlcbiAgICB9XG5cbiAgICAvLyBjYW5jZWxzIHRoZSBhbmltYXRpb25mcmFtZVxuICAsIHN0b3BBbmltRnJhbWU6IGZ1bmN0aW9uKCl7XG4gICAgICBjYW5jZWxBbmltYXRpb25GcmFtZSh0aGlzLmFuaW1hdGlvbkZyYW1lKVxuICAgIH1cblxuICAgIC8vIGtpY2tzIG9mZiB0aGUgYW5pbWF0aW9uIC0gb25seSBkb2VzIHNvbWV0aGluZyB3aGVuIHRoZSBxdWV1ZSBpcyBjdXJyZW50bHkgbm90IGFjdGl2ZSBhbmQgYXQgbGVhc3Qgb25lIHNpdHVhdGlvbiBpcyBzZXRcbiAgLCBzdGFydDogZnVuY3Rpb24oKXtcbiAgICAgIC8vIGRvbnQgc3RhcnQgaWYgYWxyZWFkeSBzdGFydGVkXG4gICAgICBpZighdGhpcy5hY3RpdmUgJiYgdGhpcy5zaXR1YXRpb24pe1xuICAgICAgICB0aGlzLmFjdGl2ZSA9IHRydWVcbiAgICAgICAgdGhpcy5zdGFydEN1cnJlbnQoKVxuICAgICAgfVxuXG4gICAgICByZXR1cm4gdGhpc1xuICAgIH1cblxuICAgIC8vIHN0YXJ0IHRoZSBjdXJyZW50IHNpdHVhdGlvblxuICAsIHN0YXJ0Q3VycmVudDogZnVuY3Rpb24oKXtcbiAgICAgIHRoaXMuc2l0dWF0aW9uLnN0YXJ0ID0gK25ldyBEYXRlICsgdGhpcy5zaXR1YXRpb24uZGVsYXkvdGhpcy5fc3BlZWRcbiAgICAgIHRoaXMuc2l0dWF0aW9uLmZpbmlzaCA9IHRoaXMuc2l0dWF0aW9uLnN0YXJ0ICsgdGhpcy5zaXR1YXRpb24uZHVyYXRpb24vdGhpcy5fc3BlZWRcbiAgICAgIHJldHVybiB0aGlzLmluaXRBbmltYXRpb25zKCkuc3RlcCgpXG4gICAgfVxuXG4gICAgLyoqXG4gICAgICogYWRkcyBhIGZ1bmN0aW9uIC8gU2l0dWF0aW9uIHRvIHRoZSBhbmltYXRpb24gcXVldWVcbiAgICAgKiBAcGFyYW0gZm4gZnVuY3Rpb24gLyBzaXR1YXRpb24gdG8gYWRkXG4gICAgICogQHJldHVybiB0aGlzXG4gICAgICovXG4gICwgcXVldWU6IGZ1bmN0aW9uKGZuKXtcbiAgICAgIGlmKHR5cGVvZiBmbiA9PSAnZnVuY3Rpb24nIHx8IGZuIGluc3RhbmNlb2YgU1ZHLlNpdHVhdGlvbilcbiAgICAgICAgdGhpcy5zaXR1YXRpb25zLnB1c2goZm4pXG5cbiAgICAgIGlmKCF0aGlzLnNpdHVhdGlvbikgdGhpcy5zaXR1YXRpb24gPSB0aGlzLnNpdHVhdGlvbnMuc2hpZnQoKVxuXG4gICAgICByZXR1cm4gdGhpc1xuICAgIH1cblxuICAgIC8qKlxuICAgICAqIHB1bGxzIG5leHQgZWxlbWVudCBmcm9tIHRoZSBxdWV1ZSBhbmQgZXhlY3V0ZSBpdFxuICAgICAqIEByZXR1cm4gdGhpc1xuICAgICAqL1xuICAsIGRlcXVldWU6IGZ1bmN0aW9uKCl7XG4gICAgICAvLyBzdG9wIGN1cnJlbnQgYW5pbWF0aW9uXG4gICAgICB0aGlzLnNpdHVhdGlvbiAmJiB0aGlzLnNpdHVhdGlvbi5zdG9wICYmIHRoaXMuc2l0dWF0aW9uLnN0b3AoKVxuXG4gICAgICAvLyBnZXQgbmV4dCBhbmltYXRpb24gZnJvbSBxdWV1ZVxuICAgICAgdGhpcy5zaXR1YXRpb24gPSB0aGlzLnNpdHVhdGlvbnMuc2hpZnQoKVxuXG4gICAgICBpZih0aGlzLnNpdHVhdGlvbil7XG4gICAgICAgIGlmKHRoaXMuc2l0dWF0aW9uIGluc3RhbmNlb2YgU1ZHLlNpdHVhdGlvbikge1xuICAgICAgICAgIHRoaXMuc3RhcnRDdXJyZW50KClcbiAgICAgICAgfSBlbHNlIHtcbiAgICAgICAgICAvLyBJZiBpdCBpcyBub3QgYSBTVkcuU2l0dWF0aW9uLCB0aGVuIGl0IGlzIGEgZnVuY3Rpb24sIHdlIGV4ZWN1dGUgaXRcbiAgICAgICAgICB0aGlzLnNpdHVhdGlvbi5jYWxsKHRoaXMpXG4gICAgICAgIH1cbiAgICAgIH1cblxuICAgICAgcmV0dXJuIHRoaXNcbiAgICB9XG5cbiAgICAvLyB1cGRhdGVzIGFsbCBhbmltYXRpb25zIHRvIHRoZSBjdXJyZW50IHN0YXRlIG9mIHRoZSBlbGVtZW50XG4gICAgLy8gdGhpcyBpcyBpbXBvcnRhbnQgd2hlbiBvbmUgcHJvcGVydHkgY291bGQgYmUgY2hhbmdlZCBmcm9tIGFub3RoZXIgcHJvcGVydHlcbiAgLCBpbml0QW5pbWF0aW9uczogZnVuY3Rpb24oKSB7XG4gICAgICB2YXIgaVxuICAgICAgdmFyIHMgPSB0aGlzLnNpdHVhdGlvblxuXG4gICAgICBpZihzLmluaXQpIHJldHVybiB0aGlzXG5cbiAgICAgIGZvcihpIGluIHMuYW5pbWF0aW9ucyl7XG5cbiAgICAgICAgaWYoaSA9PSAndmlld2JveCcpe1xuICAgICAgICAgIHMuYW5pbWF0aW9uc1tpXSA9IHRoaXMudGFyZ2V0KCkudmlld2JveCgpLm1vcnBoKHMuYW5pbWF0aW9uc1tpXSlcbiAgICAgICAgfWVsc2V7XG5cbiAgICAgICAgICAvLyBUT0RPOiB0aGlzIGlzIG5vdCBhIGNsZWFuIGNsb25lIG9mIHRoZSBhcnJheS4gV2UgbWF5IGhhdmUgc29tZSB1bmNoZWNrZWQgcmVmZXJlbmNlc1xuICAgICAgICAgIHMuYW5pbWF0aW9uc1tpXS52YWx1ZSA9IChpID09ICdwbG90JyA/IHRoaXMudGFyZ2V0KCkuYXJyYXkoKS52YWx1ZSA6IHRoaXMudGFyZ2V0KClbaV0oKSlcblxuICAgICAgICAgIC8vIHNvbWV0aW1lcyB3ZSBnZXQgYmFjayBhbiBvYmplY3QgYW5kIG5vdCB0aGUgcmVhbCB2YWx1ZSwgZml4IHRoaXNcbiAgICAgICAgICBpZihzLmFuaW1hdGlvbnNbaV0udmFsdWUudmFsdWUpe1xuICAgICAgICAgICAgcy5hbmltYXRpb25zW2ldLnZhbHVlID0gcy5hbmltYXRpb25zW2ldLnZhbHVlLnZhbHVlXG4gICAgICAgICAgfVxuXG4gICAgICAgICAgaWYocy5hbmltYXRpb25zW2ldLnJlbGF0aXZlKVxuICAgICAgICAgICAgcy5hbmltYXRpb25zW2ldLmRlc3RpbmF0aW9uLnZhbHVlID0gcy5hbmltYXRpb25zW2ldLmRlc3RpbmF0aW9uLnZhbHVlICsgcy5hbmltYXRpb25zW2ldLnZhbHVlXG5cbiAgICAgICAgfVxuXG4gICAgICB9XG5cbiAgICAgIGZvcihpIGluIHMuYXR0cnMpe1xuICAgICAgICBpZihzLmF0dHJzW2ldIGluc3RhbmNlb2YgU1ZHLkNvbG9yKXtcbiAgICAgICAgICB2YXIgY29sb3IgPSBuZXcgU1ZHLkNvbG9yKHRoaXMudGFyZ2V0KCkuYXR0cihpKSlcbiAgICAgICAgICBzLmF0dHJzW2ldLnIgPSBjb2xvci5yXG4gICAgICAgICAgcy5hdHRyc1tpXS5nID0gY29sb3IuZ1xuICAgICAgICAgIHMuYXR0cnNbaV0uYiA9IGNvbG9yLmJcbiAgICAgICAgfWVsc2V7XG4gICAgICAgICAgcy5hdHRyc1tpXS52YWx1ZSA9IHRoaXMudGFyZ2V0KCkuYXR0cihpKS8vICsgcy5hdHRyc1tpXS52YWx1ZVxuICAgICAgICB9XG4gICAgICB9XG5cbiAgICAgIGZvcihpIGluIHMuc3R5bGVzKXtcbiAgICAgICAgcy5zdHlsZXNbaV0udmFsdWUgPSB0aGlzLnRhcmdldCgpLnN0eWxlKGkpXG4gICAgICB9XG5cbiAgICAgIHMuaW5pdGlhbFRyYW5zZm9ybWF0aW9uID0gdGhpcy50YXJnZXQoKS5tYXRyaXhpZnkoKVxuXG4gICAgICBzLmluaXQgPSB0cnVlXG4gICAgICByZXR1cm4gdGhpc1xuICAgIH1cbiAgLCBjbGVhclF1ZXVlOiBmdW5jdGlvbigpe1xuICAgICAgdGhpcy5zaXR1YXRpb25zID0gW11cbiAgICAgIHJldHVybiB0aGlzXG4gICAgfVxuICAsIGNsZWFyQ3VycmVudDogZnVuY3Rpb24oKXtcbiAgICAgIHRoaXMuc2l0dWF0aW9uID0gbnVsbFxuICAgICAgcmV0dXJuIHRoaXNcbiAgICB9XG4gICAgLyoqIHN0b3BzIHRoZSBhbmltYXRpb24gaW1tZWRpYXRlbHlcbiAgICAgKiBAcGFyYW0ganVtcFRvRW5kIEEgQm9vbGVhbiBpbmRpY2F0aW5nIHdoZXRoZXIgdG8gY29tcGxldGUgdGhlIGN1cnJlbnQgYW5pbWF0aW9uIGltbWVkaWF0ZWx5LlxuICAgICAqIEBwYXJhbSBjbGVhclF1ZXVlIEEgQm9vbGVhbiBpbmRpY2F0aW5nIHdoZXRoZXIgdG8gcmVtb3ZlIHF1ZXVlZCBhbmltYXRpb24gYXMgd2VsbC5cbiAgICAgKiBAcmV0dXJuIHRoaXNcbiAgICAgKi9cbiAgLCBzdG9wOiBmdW5jdGlvbihqdW1wVG9FbmQsIGNsZWFyUXVldWUpe1xuICAgICAgaWYoIXRoaXMuYWN0aXZlKSB0aGlzLnN0YXJ0KClcblxuICAgICAgaWYoY2xlYXJRdWV1ZSl7XG4gICAgICAgIHRoaXMuY2xlYXJRdWV1ZSgpXG4gICAgICB9XG5cbiAgICAgIHRoaXMuYWN0aXZlID0gZmFsc2VcblxuICAgICAgaWYoanVtcFRvRW5kICYmIHRoaXMuc2l0dWF0aW9uKXtcbiAgICAgICAgdGhpcy5hdEVuZCgpXG4gICAgICB9XG5cbiAgICAgIHRoaXMuc3RvcEFuaW1GcmFtZSgpXG5cbiAgICAgIHJldHVybiB0aGlzLmNsZWFyQ3VycmVudCgpXG4gICAgfVxuXG4gICAgLyoqIHJlc2V0cyB0aGUgZWxlbWVudCB0byB0aGUgc3RhdGUgd2hlcmUgdGhlIGN1cnJlbnQgZWxlbWVudCBoYXMgc3RhcnRlZFxuICAgICAqIEByZXR1cm4gdGhpc1xuICAgICAqL1xuICAsIHJlc2V0OiBmdW5jdGlvbigpe1xuICAgICAgaWYodGhpcy5zaXR1YXRpb24pe1xuICAgICAgICB2YXIgdGVtcCA9IHRoaXMuc2l0dWF0aW9uXG4gICAgICAgIHRoaXMuc3RvcCgpXG4gICAgICAgIHRoaXMuc2l0dWF0aW9uID0gdGVtcFxuICAgICAgICB0aGlzLmF0U3RhcnQoKVxuICAgICAgfVxuICAgICAgcmV0dXJuIHRoaXNcbiAgICB9XG5cbiAgICAvLyBTdG9wIHRoZSBjdXJyZW50bHktcnVubmluZyBhbmltYXRpb24sIHJlbW92ZSBhbGwgcXVldWVkIGFuaW1hdGlvbnMsIGFuZCBjb21wbGV0ZSBhbGwgYW5pbWF0aW9ucyBmb3IgdGhlIGVsZW1lbnQuXG4gICwgZmluaXNoOiBmdW5jdGlvbigpe1xuXG4gICAgICB0aGlzLnN0b3AodHJ1ZSwgZmFsc2UpXG5cbiAgICAgIHdoaWxlKHRoaXMuZGVxdWV1ZSgpLnNpdHVhdGlvbiAmJiB0aGlzLnN0b3AodHJ1ZSwgZmFsc2UpKTtcblxuICAgICAgdGhpcy5jbGVhclF1ZXVlKCkuY2xlYXJDdXJyZW50KClcblxuICAgICAgcmV0dXJuIHRoaXNcbiAgICB9XG5cbiAgICAvLyBzZXQgdGhlIGludGVybmFsIGFuaW1hdGlvbiBwb2ludGVyIGF0IHRoZSBzdGFydCBwb3NpdGlvbiwgYmVmb3JlIGFueSBsb29wcywgYW5kIHVwZGF0ZXMgdGhlIHZpc3VhbGlzYXRpb25cbiAgLCBhdFN0YXJ0OiBmdW5jdGlvbigpIHtcbiAgICByZXR1cm4gdGhpcy5hdCgwLCB0cnVlKVxuICB9XG5cbiAgICAvLyBzZXQgdGhlIGludGVybmFsIGFuaW1hdGlvbiBwb2ludGVyIGF0IHRoZSBlbmQgcG9zaXRpb24sIGFmdGVyIGFsbCB0aGUgbG9vcHMsIGFuZCB1cGRhdGVzIHRoZSB2aXN1YWxpc2F0aW9uXG4gICwgYXRFbmQ6IGZ1bmN0aW9uKCkge1xuICAgIGlmICh0aGlzLnNpdHVhdGlvbi5sb29wcyA9PT0gdHJ1ZSkge1xuICAgICAgLy8gSWYgaW4gYSBpbmZpbml0ZSBsb29wLCB3ZSBlbmQgdGhlIGN1cnJlbnQgaXRlcmF0aW9uXG4gICAgICByZXR1cm4gdGhpcy5hdCh0aGlzLnNpdHVhdGlvbi5sb29wKzEsIHRydWUpXG4gICAgfSBlbHNlIGlmKHR5cGVvZiB0aGlzLnNpdHVhdGlvbi5sb29wcyA9PSAnbnVtYmVyJykge1xuICAgICAgLy8gSWYgcGVyZm9ybWluZyBhIGZpbml0ZSBudW1iZXIgb2YgbG9vcHMsIHdlIGdvIGFmdGVyIGFsbCB0aGUgbG9vcHNcbiAgICAgIHJldHVybiB0aGlzLmF0KHRoaXMuc2l0dWF0aW9uLmxvb3BzLCB0cnVlKVxuICAgIH0gZWxzZSB7XG4gICAgICAvLyBJZiBubyBsb29wcywgd2UganVzdCBnbyBhdCB0aGUgZW5kXG4gICAgICByZXR1cm4gdGhpcy5hdCgxLCB0cnVlKVxuICAgIH1cbiAgfVxuXG4gICAgLy8gc2V0IHRoZSBpbnRlcm5hbCBhbmltYXRpb24gcG9pbnRlciB0byB0aGUgc3BlY2lmaWVkIHBvc2l0aW9uIGFuZCB1cGRhdGVzIHRoZSB2aXN1YWxpc2F0aW9uXG4gICAgLy8gaWYgaXNBYnNQb3MgaXMgdHJ1ZSwgcG9zIGlzIHRyZWF0ZWQgYXMgYW4gYWJzb2x1dGUgcG9zaXRpb25cbiAgLCBhdDogZnVuY3Rpb24ocG9zLCBpc0Fic1Bvcyl7XG4gICAgICB2YXIgZHVyRGl2U3BkID0gdGhpcy5zaXR1YXRpb24uZHVyYXRpb24vdGhpcy5fc3BlZWRcblxuICAgICAgdGhpcy5hYnNQb3MgPSBwb3NcbiAgICAgIC8vIElmIHBvcyBpcyBub3QgYW4gYWJzb2x1dGUgcG9zaXRpb24sIHdlIGNvbnZlcnQgaXQgaW50byBvbmVcbiAgICAgIGlmICghaXNBYnNQb3MpIHtcbiAgICAgICAgaWYgKHRoaXMuc2l0dWF0aW9uLnJldmVyc2VkKSB0aGlzLmFic1BvcyA9IDEgLSB0aGlzLmFic1Bvc1xuICAgICAgICB0aGlzLmFic1BvcyArPSB0aGlzLnNpdHVhdGlvbi5sb29wXG4gICAgICB9XG5cbiAgICAgIHRoaXMuc2l0dWF0aW9uLnN0YXJ0ID0gK25ldyBEYXRlIC0gdGhpcy5hYnNQb3MgKiBkdXJEaXZTcGRcbiAgICAgIHRoaXMuc2l0dWF0aW9uLmZpbmlzaCA9IHRoaXMuc2l0dWF0aW9uLnN0YXJ0ICsgZHVyRGl2U3BkXG5cbiAgICAgIHJldHVybiB0aGlzLnN0ZXAodHJ1ZSlcbiAgICB9XG5cbiAgICAvKipcbiAgICAgKiBzZXRzIG9yIHJldHVybnMgdGhlIHNwZWVkIG9mIHRoZSBhbmltYXRpb25zXG4gICAgICogQHBhcmFtIHNwZWVkIG51bGwgfHwgTnVtYmVyIFRoZSBuZXcgc3BlZWQgb2YgdGhlIGFuaW1hdGlvbnNcbiAgICAgKiBAcmV0dXJuIE51bWJlciB8fCB0aGlzXG4gICAgICovXG4gICwgc3BlZWQ6IGZ1bmN0aW9uKHNwZWVkKXtcbiAgICAgIGlmIChzcGVlZCA9PT0gMCkgcmV0dXJuIHRoaXMucGF1c2UoKVxuXG4gICAgICBpZiAoc3BlZWQpIHtcbiAgICAgICAgdGhpcy5fc3BlZWQgPSBzcGVlZFxuICAgICAgICAvLyBXZSB1c2UgYW4gYWJzb2x1dGUgcG9zaXRpb24gaGVyZSBzbyB0aGF0IHNwZWVkIGNhbiBhZmZlY3QgdGhlIGRlbGF5IGJlZm9yZSB0aGUgYW5pbWF0aW9uXG4gICAgICAgIHJldHVybiB0aGlzLmF0KHRoaXMuYWJzUG9zLCB0cnVlKVxuICAgICAgfSBlbHNlIHJldHVybiB0aGlzLl9zcGVlZFxuICAgIH1cblxuICAgIC8vIE1ha2UgbG9vcGFibGVcbiAgLCBsb29wOiBmdW5jdGlvbih0aW1lcywgcmV2ZXJzZSkge1xuICAgICAgdmFyIGMgPSB0aGlzLmxhc3QoKVxuXG4gICAgICAvLyBzdG9yZSB0b3RhbCBsb29wc1xuICAgICAgYy5sb29wcyA9ICh0aW1lcyAhPSBudWxsKSA/IHRpbWVzIDogdHJ1ZVxuICAgICAgYy5sb29wID0gMFxuXG4gICAgICBpZihyZXZlcnNlKSBjLnJldmVyc2luZyA9IHRydWVcbiAgICAgIHJldHVybiB0aGlzXG4gICAgfVxuXG4gICAgLy8gcGF1c2VzIHRoZSBhbmltYXRpb25cbiAgLCBwYXVzZTogZnVuY3Rpb24oKXtcbiAgICAgIHRoaXMucGF1c2VkID0gdHJ1ZVxuICAgICAgdGhpcy5zdG9wQW5pbUZyYW1lKClcblxuICAgICAgcmV0dXJuIHRoaXNcbiAgICB9XG5cbiAgICAvLyB1bnBhdXNlIHRoZSBhbmltYXRpb25cbiAgLCBwbGF5OiBmdW5jdGlvbigpe1xuICAgICAgaWYoIXRoaXMucGF1c2VkKSByZXR1cm4gdGhpc1xuICAgICAgdGhpcy5wYXVzZWQgPSBmYWxzZVxuICAgICAgLy8gV2UgdXNlIGFuIGFic29sdXRlIHBvc2l0aW9uIGhlcmUgc28gdGhhdCB0aGUgZGVsYXkgYmVmb3JlIHRoZSBhbmltYXRpb24gY2FuIGJlIHBhdXNlZFxuICAgICAgcmV0dXJuIHRoaXMuYXQodGhpcy5hYnNQb3MsIHRydWUpXG4gICAgfVxuXG4gICAgLyoqXG4gICAgICogdG9nZ2xlIG9yIHNldCB0aGUgZGlyZWN0aW9uIG9mIHRoZSBhbmltYXRpb25cbiAgICAgKiB0cnVlIHNldHMgZGlyZWN0aW9uIHRvIGJhY2t3YXJkcyB3aGlsZSBmYWxzZSBzZXRzIGl0IHRvIGZvcndhcmRzXG4gICAgICogQHBhcmFtIHJldmVyc2VkIEJvb2xlYW4gaW5kaWNhdGluZyB3aGV0aGVyIHRvIHJldmVyc2UgdGhlIGFuaW1hdGlvbiBvciBub3QgKGRlZmF1bHQ6IHRvZ2dsZSB0aGUgcmV2ZXJzZSBzdGF0dXMpXG4gICAgICogQHJldHVybiB0aGlzXG4gICAgICovXG4gICwgcmV2ZXJzZTogZnVuY3Rpb24ocmV2ZXJzZWQpe1xuICAgICAgdmFyIGMgPSB0aGlzLmxhc3QoKVxuXG4gICAgICBpZih0eXBlb2YgcmV2ZXJzZWQgPT0gJ3VuZGVmaW5lZCcpIGMucmV2ZXJzZWQgPSAhYy5yZXZlcnNlZFxuICAgICAgZWxzZSBjLnJldmVyc2VkID0gcmV2ZXJzZWRcblxuICAgICAgcmV0dXJuIHRoaXNcbiAgICB9XG5cblxuICAgIC8qKlxuICAgICAqIHJldHVybnMgYSBmbG9hdCBmcm9tIDAtMSBpbmRpY2F0aW5nIHRoZSBwcm9ncmVzcyBvZiB0aGUgY3VycmVudCBhbmltYXRpb25cbiAgICAgKiBAcGFyYW0gZWFzZWQgQm9vbGVhbiBpbmRpY2F0aW5nIHdoZXRoZXIgdGhlIHJldHVybmVkIHBvc2l0aW9uIHNob3VsZCBiZSBlYXNlZCBvciBub3RcbiAgICAgKiBAcmV0dXJuIG51bWJlclxuICAgICAqL1xuICAsIHByb2dyZXNzOiBmdW5jdGlvbihlYXNlSXQpe1xuICAgICAgcmV0dXJuIGVhc2VJdCA/IHRoaXMuc2l0dWF0aW9uLmVhc2UodGhpcy5wb3MpIDogdGhpcy5wb3NcbiAgICB9XG5cbiAgICAvKipcbiAgICAgKiBhZGRzIGEgY2FsbGJhY2sgZnVuY3Rpb24gd2hpY2ggaXMgY2FsbGVkIHdoZW4gdGhlIGN1cnJlbnQgYW5pbWF0aW9uIGlzIGZpbmlzaGVkXG4gICAgICogQHBhcmFtIGZuIEZ1bmN0aW9uIHdoaWNoIHNob3VsZCBiZSBleGVjdXRlZCBhcyBjYWxsYmFja1xuICAgICAqIEByZXR1cm4gbnVtYmVyXG4gICAgICovXG4gICwgYWZ0ZXI6IGZ1bmN0aW9uKGZuKXtcbiAgICAgIHZhciBjID0gdGhpcy5sYXN0KClcbiAgICAgICAgLCB3cmFwcGVyID0gZnVuY3Rpb24gd3JhcHBlcihlKXtcbiAgICAgICAgICAgIGlmKGUuZGV0YWlsLnNpdHVhdGlvbiA9PSBjKXtcbiAgICAgICAgICAgICAgZm4uY2FsbCh0aGlzLCBjKVxuICAgICAgICAgICAgICB0aGlzLm9mZignZmluaXNoZWQuZngnLCB3cmFwcGVyKSAvLyBwcmV2ZW50IG1lbW9yeSBsZWFrXG4gICAgICAgICAgICB9XG4gICAgICAgICAgfVxuXG4gICAgICB0aGlzLnRhcmdldCgpLm9uKCdmaW5pc2hlZC5meCcsIHdyYXBwZXIpXG4gICAgICByZXR1cm4gdGhpc1xuICAgIH1cblxuICAgIC8vIGFkZHMgYSBjYWxsYmFjayB3aGljaCBpcyBjYWxsZWQgd2hlbmV2ZXIgb25lIGFuaW1hdGlvbiBzdGVwIGlzIHBlcmZvcm1lZFxuICAsIGR1cmluZzogZnVuY3Rpb24oZm4pe1xuICAgICAgdmFyIGMgPSB0aGlzLmxhc3QoKVxuICAgICAgICAsIHdyYXBwZXIgPSBmdW5jdGlvbihlKXtcbiAgICAgICAgICAgIGlmKGUuZGV0YWlsLnNpdHVhdGlvbiA9PSBjKXtcbiAgICAgICAgICAgICAgZm4uY2FsbCh0aGlzLCBlLmRldGFpbC5wb3MsIFNWRy5tb3JwaChlLmRldGFpbC5wb3MpLCBlLmRldGFpbC5lYXNlZCwgYylcbiAgICAgICAgICAgIH1cbiAgICAgICAgICB9XG5cbiAgICAgIC8vIHNlZSBhYm92ZVxuICAgICAgdGhpcy50YXJnZXQoKS5vZmYoJ2R1cmluZy5meCcsIHdyYXBwZXIpLm9uKCdkdXJpbmcuZngnLCB3cmFwcGVyKVxuXG4gICAgICByZXR1cm4gdGhpcy5hZnRlcihmdW5jdGlvbigpe1xuICAgICAgICB0aGlzLm9mZignZHVyaW5nLmZ4Jywgd3JhcHBlcilcbiAgICAgIH0pXG4gICAgfVxuXG4gICAgLy8gY2FsbHMgYWZ0ZXIgQUxMIGFuaW1hdGlvbnMgaW4gdGhlIHF1ZXVlIGFyZSBmaW5pc2hlZFxuICAsIGFmdGVyQWxsOiBmdW5jdGlvbihmbil7XG4gICAgICB2YXIgd3JhcHBlciA9IGZ1bmN0aW9uIHdyYXBwZXIoZSl7XG4gICAgICAgICAgICBmbi5jYWxsKHRoaXMpXG4gICAgICAgICAgICB0aGlzLm9mZignYWxsZmluaXNoZWQuZngnLCB3cmFwcGVyKVxuICAgICAgICAgIH1cblxuICAgICAgLy8gc2VlIGFib3ZlXG4gICAgICB0aGlzLnRhcmdldCgpLm9mZignYWxsZmluaXNoZWQuZngnLCB3cmFwcGVyKS5vbignYWxsZmluaXNoZWQuZngnLCB3cmFwcGVyKVxuICAgICAgcmV0dXJuIHRoaXNcbiAgICB9XG5cbiAgICAvLyBjYWxscyBvbiBldmVyeSBhbmltYXRpb24gc3RlcCBmb3IgYWxsIGFuaW1hdGlvbnNcbiAgLCBkdXJpbmdBbGw6IGZ1bmN0aW9uKGZuKXtcbiAgICAgIHZhciB3cmFwcGVyID0gZnVuY3Rpb24oZSl7XG4gICAgICAgICAgICBmbi5jYWxsKHRoaXMsIGUuZGV0YWlsLnBvcywgU1ZHLm1vcnBoKGUuZGV0YWlsLnBvcyksIGUuZGV0YWlsLmVhc2VkLCBlLmRldGFpbC5zaXR1YXRpb24pXG4gICAgICAgICAgfVxuXG4gICAgICB0aGlzLnRhcmdldCgpLm9mZignZHVyaW5nLmZ4Jywgd3JhcHBlcikub24oJ2R1cmluZy5meCcsIHdyYXBwZXIpXG5cbiAgICAgIHJldHVybiB0aGlzLmFmdGVyQWxsKGZ1bmN0aW9uKCl7XG4gICAgICAgIHRoaXMub2ZmKCdkdXJpbmcuZngnLCB3cmFwcGVyKVxuICAgICAgfSlcbiAgICB9XG5cbiAgLCBsYXN0OiBmdW5jdGlvbigpe1xuICAgICAgcmV0dXJuIHRoaXMuc2l0dWF0aW9ucy5sZW5ndGggPyB0aGlzLnNpdHVhdGlvbnNbdGhpcy5zaXR1YXRpb25zLmxlbmd0aC0xXSA6IHRoaXMuc2l0dWF0aW9uXG4gICAgfVxuXG4gICAgLy8gYWRkcyBvbmUgcHJvcGVydHkgdG8gdGhlIGFuaW1hdGlvbnNcbiAgLCBhZGQ6IGZ1bmN0aW9uKG1ldGhvZCwgYXJncywgdHlwZSl7XG4gICAgICB0aGlzLmxhc3QoKVt0eXBlIHx8ICdhbmltYXRpb25zJ11bbWV0aG9kXSA9IGFyZ3NcbiAgICAgIHNldFRpbWVvdXQoZnVuY3Rpb24oKXt0aGlzLnN0YXJ0KCl9LmJpbmQodGhpcyksIDApXG4gICAgICByZXR1cm4gdGhpc1xuICAgIH1cblxuICAgIC8qKiBwZXJmb3JtIG9uZSBzdGVwIG9mIHRoZSBhbmltYXRpb25cbiAgICAgKiAgQHBhcmFtIGlnbm9yZVRpbWUgQm9vbGVhbiBpbmRpY2F0aW5nIHdoZXRoZXIgdG8gaWdub3JlIHRpbWUgYW5kIHVzZSBwb3NpdGlvbiBkaXJlY3RseSBvciByZWNhbGN1bGF0ZSBwb3NpdGlvbiBiYXNlZCBvbiB0aW1lXG4gICAgICogIEByZXR1cm4gdGhpc1xuICAgICAqL1xuICAsIHN0ZXA6IGZ1bmN0aW9uKGlnbm9yZVRpbWUpe1xuXG4gICAgICAvLyBjb252ZXJ0IGN1cnJlbnQgdGltZSB0byBhbiBhYnNvbHV0ZSBwb3NpdGlvblxuICAgICAgaWYoIWlnbm9yZVRpbWUpIHRoaXMuYWJzUG9zID0gdGhpcy50aW1lVG9BYnNQb3MoK25ldyBEYXRlKVxuXG4gICAgICAvLyBUaGlzIHBhcnQgY29udmVydCBhbiBhYnNvbHV0ZSBwb3NpdGlvbiB0byBhIHBvc2l0aW9uXG4gICAgICBpZih0aGlzLnNpdHVhdGlvbi5sb29wcyAhPT0gZmFsc2UpIHtcbiAgICAgICAgdmFyIGFic1BvcywgYWJzUG9zSW50LCBsYXN0TG9vcFxuXG4gICAgICAgIC8vIElmIHRoZSBhYnNvbHV0ZSBwb3NpdGlvbiBpcyBiZWxvdyAwLCB3ZSBqdXN0IHRyZWF0IGl0IGFzIGlmIGl0IHdhcyAwXG4gICAgICAgIGFic1BvcyA9IE1hdGgubWF4KHRoaXMuYWJzUG9zLCAwKVxuICAgICAgICBhYnNQb3NJbnQgPSBNYXRoLmZsb29yKGFic1BvcylcblxuICAgICAgICBpZih0aGlzLnNpdHVhdGlvbi5sb29wcyA9PT0gdHJ1ZSB8fCBhYnNQb3NJbnQgPCB0aGlzLnNpdHVhdGlvbi5sb29wcykge1xuICAgICAgICAgIHRoaXMucG9zID0gYWJzUG9zIC0gYWJzUG9zSW50XG4gICAgICAgICAgbGFzdExvb3AgPSB0aGlzLnNpdHVhdGlvbi5sb29wXG4gICAgICAgICAgdGhpcy5zaXR1YXRpb24ubG9vcCA9IGFic1Bvc0ludFxuICAgICAgICB9IGVsc2Uge1xuICAgICAgICAgIHRoaXMuYWJzUG9zID0gdGhpcy5zaXR1YXRpb24ubG9vcHNcbiAgICAgICAgICB0aGlzLnBvcyA9IDFcbiAgICAgICAgICAvLyBUaGUgLTEgaGVyZSBpcyBiZWNhdXNlIHdlIGRvbid0IHdhbnQgdG8gdG9nZ2xlIHJldmVyc2VkIHdoZW4gYWxsIHRoZSBsb29wcyBoYXZlIGJlZW4gY29tcGxldGVkXG4gICAgICAgICAgbGFzdExvb3AgPSB0aGlzLnNpdHVhdGlvbi5sb29wIC0gMVxuICAgICAgICAgIHRoaXMuc2l0dWF0aW9uLmxvb3AgPSB0aGlzLnNpdHVhdGlvbi5sb29wc1xuICAgICAgICB9XG5cbiAgICAgICAgaWYodGhpcy5zaXR1YXRpb24ucmV2ZXJzaW5nKSB7XG4gICAgICAgICAgLy8gVG9nZ2xlIHJldmVyc2VkIGlmIGFuIG9kZCBudW1iZXIgb2YgbG9vcHMgYXMgb2NjdXJlZCBzaW5jZSB0aGUgbGFzdCBjYWxsIG9mIHN0ZXBcbiAgICAgICAgICB0aGlzLnNpdHVhdGlvbi5yZXZlcnNlZCA9IHRoaXMuc2l0dWF0aW9uLnJldmVyc2VkICE9IEJvb2xlYW4oKHRoaXMuc2l0dWF0aW9uLmxvb3AgLSBsYXN0TG9vcCkgJSAyKVxuICAgICAgICB9XG5cbiAgICAgIH0gZWxzZSB7XG4gICAgICAgIC8vIElmIHRoZXJlIGFyZSBubyBsb29wLCB0aGUgYWJzb2x1dGUgcG9zaXRpb24gbXVzdCBub3QgYmUgYWJvdmUgMVxuICAgICAgICB0aGlzLmFic1BvcyA9IE1hdGgubWluKHRoaXMuYWJzUG9zLCAxKVxuICAgICAgICB0aGlzLnBvcyA9IHRoaXMuYWJzUG9zXG4gICAgICB9XG5cbiAgICAgIC8vIHdoaWxlIHRoZSBhYnNvbHV0ZSBwb3NpdGlvbiBjYW4gYmUgYmVsb3cgMCwgdGhlIHBvc2l0aW9uIG11c3Qgbm90IGJlIGJlbG93IDBcbiAgICAgIGlmKHRoaXMucG9zIDwgMCkgdGhpcy5wb3MgPSAwXG5cbiAgICAgIGlmKHRoaXMuc2l0dWF0aW9uLnJldmVyc2VkKSB0aGlzLnBvcyA9IDEgLSB0aGlzLnBvc1xuXG5cbiAgICAgIC8vIGFwcGx5IGVhc2luZ1xuICAgICAgdmFyIGVhc2VkID0gdGhpcy5zaXR1YXRpb24uZWFzZSh0aGlzLnBvcylcblxuICAgICAgLy8gY2FsbCBvbmNlLWNhbGxiYWNrc1xuICAgICAgZm9yKHZhciBpIGluIHRoaXMuc2l0dWF0aW9uLm9uY2Upe1xuICAgICAgICBpZihpID4gdGhpcy5sYXN0UG9zICYmIGkgPD0gZWFzZWQpe1xuICAgICAgICAgIHRoaXMuc2l0dWF0aW9uLm9uY2VbaV0uY2FsbCh0aGlzLnRhcmdldCgpLCB0aGlzLnBvcywgZWFzZWQpXG4gICAgICAgICAgZGVsZXRlIHRoaXMuc2l0dWF0aW9uLm9uY2VbaV1cbiAgICAgICAgfVxuICAgICAgfVxuXG4gICAgICAvLyBmaXJlIGR1cmluZyBjYWxsYmFjayB3aXRoIHBvc2l0aW9uLCBlYXNlZCBwb3NpdGlvbiBhbmQgY3VycmVudCBzaXR1YXRpb24gYXMgcGFyYW1ldGVyXG4gICAgICBpZih0aGlzLmFjdGl2ZSkgdGhpcy50YXJnZXQoKS5maXJlKCdkdXJpbmcnLCB7cG9zOiB0aGlzLnBvcywgZWFzZWQ6IGVhc2VkLCBmeDogdGhpcywgc2l0dWF0aW9uOiB0aGlzLnNpdHVhdGlvbn0pXG5cbiAgICAgIC8vIHRoZSB1c2VyIG1heSBjYWxsIHN0b3Agb3IgZmluaXNoIGluIHRoZSBkdXJpbmcgY2FsbGJhY2tcbiAgICAgIC8vIHNvIG1ha2Ugc3VyZSB0aGF0IHdlIHN0aWxsIGhhdmUgYSB2YWxpZCBzaXR1YXRpb25cbiAgICAgIGlmKCF0aGlzLnNpdHVhdGlvbil7XG4gICAgICAgIHJldHVybiB0aGlzXG4gICAgICB9XG5cbiAgICAgIC8vIGFwcGx5IHRoZSBhY3R1YWwgYW5pbWF0aW9uIHRvIGV2ZXJ5IHByb3BlcnR5XG4gICAgICB0aGlzLmVhY2hBdCgpXG5cbiAgICAgIC8vIGRvIGZpbmFsIGNvZGUgd2hlbiBzaXR1YXRpb24gaXMgZmluaXNoZWRcbiAgICAgIGlmKCh0aGlzLnBvcyA9PSAxICYmICF0aGlzLnNpdHVhdGlvbi5yZXZlcnNlZCkgfHwgKHRoaXMuc2l0dWF0aW9uLnJldmVyc2VkICYmIHRoaXMucG9zID09IDApKXtcblxuICAgICAgICAvLyBzdG9wIGFuaW1hdGlvbiBjYWxsYmFja1xuICAgICAgICB0aGlzLnN0b3BBbmltRnJhbWUoKVxuXG4gICAgICAgIC8vIGZpcmUgZmluaXNoZWQgY2FsbGJhY2sgd2l0aCBjdXJyZW50IHNpdHVhdGlvbiBhcyBwYXJhbWV0ZXJcbiAgICAgICAgdGhpcy50YXJnZXQoKS5maXJlKCdmaW5pc2hlZCcsIHtmeDp0aGlzLCBzaXR1YXRpb246IHRoaXMuc2l0dWF0aW9ufSlcblxuICAgICAgICBpZighdGhpcy5zaXR1YXRpb25zLmxlbmd0aCl7XG4gICAgICAgICAgdGhpcy50YXJnZXQoKS5maXJlKCdhbGxmaW5pc2hlZCcpXG4gICAgICAgICAgdGhpcy50YXJnZXQoKS5vZmYoJy5meCcpIC8vIHRoZXJlIHNob3VsZG50IGJlIGFueSBiaW5kaW5nIGxlZnQsIGJ1dCB0byBtYWtlIHN1cmUuLi5cbiAgICAgICAgICB0aGlzLmFjdGl2ZSA9IGZhbHNlXG4gICAgICAgIH1cblxuICAgICAgICAvLyBzdGFydCBuZXh0IGFuaW1hdGlvblxuICAgICAgICBpZih0aGlzLmFjdGl2ZSkgdGhpcy5kZXF1ZXVlKClcbiAgICAgICAgZWxzZSB0aGlzLmNsZWFyQ3VycmVudCgpXG5cbiAgICAgIH1lbHNlIGlmKCF0aGlzLnBhdXNlZCAmJiB0aGlzLmFjdGl2ZSl7XG4gICAgICAgIC8vIHdlIGNvbnRpbnVlIGFuaW1hdGluZyB3aGVuIHdlIGFyZSBub3QgYXQgdGhlIGVuZFxuICAgICAgICB0aGlzLnN0YXJ0QW5pbUZyYW1lKClcbiAgICAgIH1cblxuICAgICAgLy8gc2F2ZSBsYXN0IGVhc2VkIHBvc2l0aW9uIGZvciBvbmNlIGNhbGxiYWNrIHRyaWdnZXJpbmdcbiAgICAgIHRoaXMubGFzdFBvcyA9IGVhc2VkXG4gICAgICByZXR1cm4gdGhpc1xuXG4gICAgfVxuXG4gICAgLy8gY2FsY3VsYXRlcyB0aGUgc3RlcCBmb3IgZXZlcnkgcHJvcGVydHkgYW5kIGNhbGxzIGJsb2NrIHdpdGggaXRcbiAgLCBlYWNoQXQ6IGZ1bmN0aW9uKCl7XG4gICAgICB2YXIgaSwgYXQsIHNlbGYgPSB0aGlzLCB0YXJnZXQgPSB0aGlzLnRhcmdldCgpLCBzID0gdGhpcy5zaXR1YXRpb25cblxuICAgICAgLy8gYXBwbHkgYW5pbWF0aW9ucyB3aGljaCBjYW4gYmUgY2FsbGVkIHRyb3VnaCBhIG1ldGhvZFxuICAgICAgZm9yKGkgaW4gcy5hbmltYXRpb25zKXtcblxuICAgICAgICBhdCA9IFtdLmNvbmNhdChzLmFuaW1hdGlvbnNbaV0pLm1hcChmdW5jdGlvbihlbCl7XG4gICAgICAgICAgcmV0dXJuIHR5cGVvZiBlbCAhPT0gJ3N0cmluZycgJiYgZWwuYXQgPyBlbC5hdChzLmVhc2Uoc2VsZi5wb3MpLCBzZWxmLnBvcykgOiBlbFxuICAgICAgICB9KVxuXG4gICAgICAgIHRhcmdldFtpXS5hcHBseSh0YXJnZXQsIGF0KVxuXG4gICAgICB9XG5cbiAgICAgIC8vIGFwcGx5IGFuaW1hdGlvbiB3aGljaCBoYXMgdG8gYmUgYXBwbGllZCB3aXRoIGF0dHIoKVxuICAgICAgZm9yKGkgaW4gcy5hdHRycyl7XG5cbiAgICAgICAgYXQgPSBbaV0uY29uY2F0KHMuYXR0cnNbaV0pLm1hcChmdW5jdGlvbihlbCl7XG4gICAgICAgICAgcmV0dXJuIHR5cGVvZiBlbCAhPT0gJ3N0cmluZycgJiYgZWwuYXQgPyBlbC5hdChzLmVhc2Uoc2VsZi5wb3MpLCBzZWxmLnBvcykgOiBlbFxuICAgICAgICB9KVxuXG4gICAgICAgIHRhcmdldC5hdHRyLmFwcGx5KHRhcmdldCwgYXQpXG5cbiAgICAgIH1cblxuICAgICAgLy8gYXBwbHkgYW5pbWF0aW9uIHdoaWNoIGhhcyB0byBiZSBhcHBsaWVkIHdpdGggc3R5bGUoKVxuICAgICAgZm9yKGkgaW4gcy5zdHlsZXMpe1xuXG4gICAgICAgIGF0ID0gW2ldLmNvbmNhdChzLnN0eWxlc1tpXSkubWFwKGZ1bmN0aW9uKGVsKXtcbiAgICAgICAgICByZXR1cm4gdHlwZW9mIGVsICE9PSAnc3RyaW5nJyAmJiBlbC5hdCA/IGVsLmF0KHMuZWFzZShzZWxmLnBvcyksIHNlbGYucG9zKSA6IGVsXG4gICAgICAgIH0pXG5cbiAgICAgICAgdGFyZ2V0LnN0eWxlLmFwcGx5KHRhcmdldCwgYXQpXG5cbiAgICAgIH1cblxuICAgICAgLy8gYW5pbWF0ZSBpbml0aWFsVHJhbnNmb3JtYXRpb24gd2hpY2ggaGFzIHRvIGJlIGNoYWluZWRcbiAgICAgIGlmKHMudHJhbnNmb3Jtcy5sZW5ndGgpe1xuXG4gICAgICAgIC8vIGdldCBpbml0aWFsIGluaXRpYWxUcmFuc2Zvcm1hdGlvblxuICAgICAgICBhdCA9IHMuaW5pdGlhbFRyYW5zZm9ybWF0aW9uXG4gICAgICAgIGZvcihpID0gMCwgbGVuID0gcy50cmFuc2Zvcm1zLmxlbmd0aDsgaSA8IGxlbjsgaSsrKXtcblxuICAgICAgICAgIC8vIGdldCBuZXh0IHRyYW5zZm9ybWF0aW9uIGluIGNoYWluXG4gICAgICAgICAgdmFyIGEgPSBzLnRyYW5zZm9ybXNbaV1cblxuICAgICAgICAgIC8vIG11bHRpcGx5IG1hdHJpeCBkaXJlY3RseVxuICAgICAgICAgIGlmKGEgaW5zdGFuY2VvZiBTVkcuTWF0cml4KXtcblxuICAgICAgICAgICAgaWYoYS5yZWxhdGl2ZSl7XG4gICAgICAgICAgICAgIGF0ID0gYXQubXVsdGlwbHkobmV3IFNWRy5NYXRyaXgoKS5tb3JwaChhKS5hdChzLmVhc2UodGhpcy5wb3MpKSlcbiAgICAgICAgICAgIH1lbHNle1xuICAgICAgICAgICAgICBhdCA9IGF0Lm1vcnBoKGEpLmF0KHMuZWFzZSh0aGlzLnBvcykpXG4gICAgICAgICAgICB9XG4gICAgICAgICAgICBjb250aW51ZVxuICAgICAgICAgIH1cblxuICAgICAgICAgIC8vIHdoZW4gdHJhbnNmb3JtYXRpb24gaXMgYWJzb2x1dGUgd2UgaGF2ZSB0byByZXNldCB0aGUgbmVlZGVkIHRyYW5zZm9ybWF0aW9uIGZpcnN0XG4gICAgICAgICAgaWYoIWEucmVsYXRpdmUpXG4gICAgICAgICAgICBhLnVuZG8oYXQuZXh0cmFjdCgpKVxuXG4gICAgICAgICAgLy8gYW5kIHJlYXBwbHkgaXQgYWZ0ZXJcbiAgICAgICAgICBhdCA9IGF0Lm11bHRpcGx5KGEuYXQocy5lYXNlKHRoaXMucG9zKSkpXG5cbiAgICAgICAgfVxuXG4gICAgICAgIC8vIHNldCBuZXcgbWF0cml4IG9uIGVsZW1lbnRcbiAgICAgICAgdGFyZ2V0Lm1hdHJpeChhdClcbiAgICAgIH1cblxuICAgICAgcmV0dXJuIHRoaXNcblxuICAgIH1cblxuXG4gICAgLy8gYWRkcyBhbiBvbmNlLWNhbGxiYWNrIHdoaWNoIGlzIGNhbGxlZCBhdCBhIHNwZWNpZmljIHBvc2l0aW9uIGFuZCBuZXZlciBhZ2FpblxuICAsIG9uY2U6IGZ1bmN0aW9uKHBvcywgZm4sIGlzRWFzZWQpe1xuXG4gICAgICBpZighaXNFYXNlZClwb3MgPSB0aGlzLnNpdHVhdGlvbi5lYXNlKHBvcylcblxuICAgICAgdGhpcy5zaXR1YXRpb24ub25jZVtwb3NdID0gZm5cblxuICAgICAgcmV0dXJuIHRoaXNcbiAgICB9XG5cbiAgfVxuXG4sIHBhcmVudDogU1ZHLkVsZW1lbnRcblxuICAvLyBBZGQgbWV0aG9kIHRvIHBhcmVudCBlbGVtZW50c1xuLCBjb25zdHJ1Y3Q6IHtcbiAgICAvLyBHZXQgZnggbW9kdWxlIG9yIGNyZWF0ZSBhIG5ldyBvbmUsIHRoZW4gYW5pbWF0ZSB3aXRoIGdpdmVuIGR1cmF0aW9uIGFuZCBlYXNlXG4gICAgYW5pbWF0ZTogZnVuY3Rpb24obywgZWFzZSwgZGVsYXkpIHtcbiAgICAgIHJldHVybiAodGhpcy5meCB8fCAodGhpcy5meCA9IG5ldyBTVkcuRlgodGhpcykpKS5hbmltYXRlKG8sIGVhc2UsIGRlbGF5KVxuICAgIH1cbiAgLCBkZWxheTogZnVuY3Rpb24oZGVsYXkpe1xuICAgICAgcmV0dXJuICh0aGlzLmZ4IHx8ICh0aGlzLmZ4ID0gbmV3IFNWRy5GWCh0aGlzKSkpLmRlbGF5KGRlbGF5KVxuICAgIH1cbiAgLCBzdG9wOiBmdW5jdGlvbihqdW1wVG9FbmQsIGNsZWFyUXVldWUpIHtcbiAgICAgIGlmICh0aGlzLmZ4KVxuICAgICAgICB0aGlzLmZ4LnN0b3AoanVtcFRvRW5kLCBjbGVhclF1ZXVlKVxuXG4gICAgICByZXR1cm4gdGhpc1xuICAgIH1cbiAgLCBmaW5pc2g6IGZ1bmN0aW9uKCkge1xuICAgICAgaWYgKHRoaXMuZngpXG4gICAgICAgIHRoaXMuZnguZmluaXNoKClcblxuICAgICAgcmV0dXJuIHRoaXNcbiAgICB9XG4gICAgLy8gUGF1c2UgY3VycmVudCBhbmltYXRpb25cbiAgLCBwYXVzZTogZnVuY3Rpb24oKSB7XG4gICAgICBpZiAodGhpcy5meClcbiAgICAgICAgdGhpcy5meC5wYXVzZSgpXG5cbiAgICAgIHJldHVybiB0aGlzXG4gICAgfVxuICAgIC8vIFBsYXkgcGF1c2VkIGN1cnJlbnQgYW5pbWF0aW9uXG4gICwgcGxheTogZnVuY3Rpb24oKSB7XG4gICAgICBpZiAodGhpcy5meClcbiAgICAgICAgdGhpcy5meC5wbGF5KClcblxuICAgICAgcmV0dXJuIHRoaXNcbiAgICB9XG4gICAgLy8gU2V0L0dldCB0aGUgc3BlZWQgb2YgdGhlIGFuaW1hdGlvbnNcbiAgLCBzcGVlZDogZnVuY3Rpb24oc3BlZWQpIHtcbiAgICAgIGlmICh0aGlzLmZ4KVxuICAgICAgICBpZiAoc3BlZWQgPT0gbnVsbClcbiAgICAgICAgICByZXR1cm4gdGhpcy5meC5zcGVlZCgpXG4gICAgICAgIGVsc2VcbiAgICAgICAgICB0aGlzLmZ4LnNwZWVkKHNwZWVkKVxuXG4gICAgICByZXR1cm4gdGhpc1xuICAgIH1cbiAgfVxuXG59KVxuXG4vLyBNb3JwaE9iaiBpcyB1c2VkIHdoZW5ldmVyIG5vIG1vcnBoYWJsZSBvYmplY3QgaXMgZ2l2ZW5cblNWRy5Nb3JwaE9iaiA9IFNWRy5pbnZlbnQoe1xuXG4gIGNyZWF0ZTogZnVuY3Rpb24oZnJvbSwgdG8pe1xuICAgIC8vIHByZXBhcmUgY29sb3IgZm9yIG1vcnBoaW5nXG4gICAgaWYoU1ZHLkNvbG9yLmlzQ29sb3IodG8pKSByZXR1cm4gbmV3IFNWRy5Db2xvcihmcm9tKS5tb3JwaCh0bylcbiAgICAvLyBwcmVwYXJlIG51bWJlciBmb3IgbW9ycGhpbmdcbiAgICBpZihTVkcucmVnZXgubnVtYmVyQW5kVW5pdC50ZXN0KHRvKSkgcmV0dXJuIG5ldyBTVkcuTnVtYmVyKGZyb20pLm1vcnBoKHRvKVxuXG4gICAgLy8gcHJlcGFyZSBmb3IgcGxhaW4gbW9ycGhpbmdcbiAgICB0aGlzLnZhbHVlID0gMFxuICAgIHRoaXMuZGVzdGluYXRpb24gPSB0b1xuICB9XG5cbiwgZXh0ZW5kOiB7XG4gICAgYXQ6IGZ1bmN0aW9uKHBvcywgcmVhbCl7XG4gICAgICByZXR1cm4gcmVhbCA8IDEgPyB0aGlzLnZhbHVlIDogdGhpcy5kZXN0aW5hdGlvblxuICAgIH0sXG5cbiAgICB2YWx1ZU9mOiBmdW5jdGlvbigpe1xuICAgICAgcmV0dXJuIHRoaXMudmFsdWVcbiAgICB9XG4gIH1cblxufSlcblxuU1ZHLmV4dGVuZChTVkcuRlgsIHtcbiAgLy8gQWRkIGFuaW1hdGFibGUgYXR0cmlidXRlc1xuICBhdHRyOiBmdW5jdGlvbihhLCB2LCByZWxhdGl2ZSkge1xuICAgIC8vIGFwcGx5IGF0dHJpYnV0ZXMgaW5kaXZpZHVhbGx5XG4gICAgaWYgKHR5cGVvZiBhID09ICdvYmplY3QnKSB7XG4gICAgICBmb3IgKHZhciBrZXkgaW4gYSlcbiAgICAgICAgdGhpcy5hdHRyKGtleSwgYVtrZXldKVxuXG4gICAgfSBlbHNlIHtcbiAgICAgIC8vIHRoZSBNb3JwaE9iaiB0YWtlcyBjYXJlIGFib3V0IHRoZSByaWdodCBmdW5jdGlvbiB1c2VkXG4gICAgICB0aGlzLmFkZChhLCBuZXcgU1ZHLk1vcnBoT2JqKG51bGwsIHYpLCAnYXR0cnMnKVxuICAgIH1cblxuICAgIHJldHVybiB0aGlzXG4gIH1cbiAgLy8gQWRkIGFuaW1hdGFibGUgc3R5bGVzXG4sIHN0eWxlOiBmdW5jdGlvbihzLCB2KSB7XG4gICAgaWYgKHR5cGVvZiBzID09ICdvYmplY3QnKVxuICAgICAgZm9yICh2YXIga2V5IGluIHMpXG4gICAgICAgIHRoaXMuc3R5bGUoa2V5LCBzW2tleV0pXG5cbiAgICBlbHNlXG4gICAgICB0aGlzLmFkZChzLCBuZXcgU1ZHLk1vcnBoT2JqKG51bGwsIHYpLCAnc3R5bGVzJylcblxuICAgIHJldHVybiB0aGlzXG4gIH1cbiAgLy8gQW5pbWF0YWJsZSB4LWF4aXNcbiwgeDogZnVuY3Rpb24oeCwgcmVsYXRpdmUpIHtcbiAgICBpZih0aGlzLnRhcmdldCgpIGluc3RhbmNlb2YgU1ZHLkcpe1xuICAgICAgdGhpcy50cmFuc2Zvcm0oe3g6eH0sIHJlbGF0aXZlKVxuICAgICAgcmV0dXJuIHRoaXNcbiAgICB9XG5cbiAgICB2YXIgbnVtID0gbmV3IFNWRy5OdW1iZXIoKS5tb3JwaCh4KVxuICAgIG51bS5yZWxhdGl2ZSA9IHJlbGF0aXZlXG4gICAgcmV0dXJuIHRoaXMuYWRkKCd4JywgbnVtKVxuICB9XG4gIC8vIEFuaW1hdGFibGUgeS1heGlzXG4sIHk6IGZ1bmN0aW9uKHksIHJlbGF0aXZlKSB7XG4gICAgaWYodGhpcy50YXJnZXQoKSBpbnN0YW5jZW9mIFNWRy5HKXtcbiAgICAgIHRoaXMudHJhbnNmb3JtKHt5Onl9LCByZWxhdGl2ZSlcbiAgICAgIHJldHVybiB0aGlzXG4gICAgfVxuXG4gICAgdmFyIG51bSA9IG5ldyBTVkcuTnVtYmVyKCkubW9ycGgoeSlcbiAgICBudW0ucmVsYXRpdmUgPSByZWxhdGl2ZVxuICAgIHJldHVybiB0aGlzLmFkZCgneScsIG51bSlcbiAgfVxuICAvLyBBbmltYXRhYmxlIGNlbnRlciB4LWF4aXNcbiwgY3g6IGZ1bmN0aW9uKHgpIHtcbiAgICByZXR1cm4gdGhpcy5hZGQoJ2N4JywgbmV3IFNWRy5OdW1iZXIoKS5tb3JwaCh4KSlcbiAgfVxuICAvLyBBbmltYXRhYmxlIGNlbnRlciB5LWF4aXNcbiwgY3k6IGZ1bmN0aW9uKHkpIHtcbiAgICByZXR1cm4gdGhpcy5hZGQoJ2N5JywgbmV3IFNWRy5OdW1iZXIoKS5tb3JwaCh5KSlcbiAgfVxuICAvLyBBZGQgYW5pbWF0YWJsZSBtb3ZlXG4sIG1vdmU6IGZ1bmN0aW9uKHgsIHkpIHtcbiAgICByZXR1cm4gdGhpcy54KHgpLnkoeSlcbiAgfVxuICAvLyBBZGQgYW5pbWF0YWJsZSBjZW50ZXJcbiwgY2VudGVyOiBmdW5jdGlvbih4LCB5KSB7XG4gICAgcmV0dXJuIHRoaXMuY3goeCkuY3koeSlcbiAgfVxuICAvLyBBZGQgYW5pbWF0YWJsZSBzaXplXG4sIHNpemU6IGZ1bmN0aW9uKHdpZHRoLCBoZWlnaHQpIHtcbiAgICBpZiAodGhpcy50YXJnZXQoKSBpbnN0YW5jZW9mIFNWRy5UZXh0KSB7XG4gICAgICAvLyBhbmltYXRlIGZvbnQgc2l6ZSBmb3IgVGV4dCBlbGVtZW50c1xuICAgICAgdGhpcy5hdHRyKCdmb250LXNpemUnLCB3aWR0aClcblxuICAgIH0gZWxzZSB7XG4gICAgICAvLyBhbmltYXRlIGJib3ggYmFzZWQgc2l6ZSBmb3IgYWxsIG90aGVyIGVsZW1lbnRzXG4gICAgICB2YXIgYm94XG5cbiAgICAgIGlmKCF3aWR0aCB8fCAhaGVpZ2h0KXtcbiAgICAgICAgYm94ID0gdGhpcy50YXJnZXQoKS5iYm94KClcbiAgICAgIH1cblxuICAgICAgaWYoIXdpZHRoKXtcbiAgICAgICAgd2lkdGggPSBib3gud2lkdGggLyBib3guaGVpZ2h0ICAqIGhlaWdodFxuICAgICAgfVxuXG4gICAgICBpZighaGVpZ2h0KXtcbiAgICAgICAgaGVpZ2h0ID0gYm94LmhlaWdodCAvIGJveC53aWR0aCAgKiB3aWR0aFxuICAgICAgfVxuXG4gICAgICB0aGlzLmFkZCgnd2lkdGgnICwgbmV3IFNWRy5OdW1iZXIoKS5tb3JwaCh3aWR0aCkpXG4gICAgICAgICAgLmFkZCgnaGVpZ2h0JywgbmV3IFNWRy5OdW1iZXIoKS5tb3JwaChoZWlnaHQpKVxuXG4gICAgfVxuXG4gICAgcmV0dXJuIHRoaXNcbiAgfVxuICAvLyBBZGQgYW5pbWF0YWJsZSBwbG90XG4sIHBsb3Q6IGZ1bmN0aW9uKHApIHtcbiAgICByZXR1cm4gdGhpcy5hZGQoJ3Bsb3QnLCB0aGlzLnRhcmdldCgpLmFycmF5KCkubW9ycGgocCkpXG4gIH1cbiAgLy8gQWRkIGxlYWRpbmcgbWV0aG9kXG4sIGxlYWRpbmc6IGZ1bmN0aW9uKHZhbHVlKSB7XG4gICAgcmV0dXJuIHRoaXMudGFyZ2V0KCkubGVhZGluZyA/XG4gICAgICB0aGlzLmFkZCgnbGVhZGluZycsIG5ldyBTVkcuTnVtYmVyKCkubW9ycGgodmFsdWUpKSA6XG4gICAgICB0aGlzXG4gIH1cbiAgLy8gQWRkIGFuaW1hdGFibGUgdmlld2JveFxuLCB2aWV3Ym94OiBmdW5jdGlvbih4LCB5LCB3aWR0aCwgaGVpZ2h0KSB7XG4gICAgaWYgKHRoaXMudGFyZ2V0KCkgaW5zdGFuY2VvZiBTVkcuQ29udGFpbmVyKSB7XG4gICAgICB0aGlzLmFkZCgndmlld2JveCcsIG5ldyBTVkcuVmlld0JveCh4LCB5LCB3aWR0aCwgaGVpZ2h0KSlcbiAgICB9XG5cbiAgICByZXR1cm4gdGhpc1xuICB9XG4sIHVwZGF0ZTogZnVuY3Rpb24obykge1xuICAgIGlmICh0aGlzLnRhcmdldCgpIGluc3RhbmNlb2YgU1ZHLlN0b3ApIHtcbiAgICAgIGlmICh0eXBlb2YgbyA9PSAnbnVtYmVyJyB8fCBvIGluc3RhbmNlb2YgU1ZHLk51bWJlcikge1xuICAgICAgICByZXR1cm4gdGhpcy51cGRhdGUoe1xuICAgICAgICAgIG9mZnNldDogIGFyZ3VtZW50c1swXVxuICAgICAgICAsIGNvbG9yOiAgIGFyZ3VtZW50c1sxXVxuICAgICAgICAsIG9wYWNpdHk6IGFyZ3VtZW50c1syXVxuICAgICAgICB9KVxuICAgICAgfVxuXG4gICAgICBpZiAoby5vcGFjaXR5ICE9IG51bGwpIHRoaXMuYXR0cignc3RvcC1vcGFjaXR5Jywgby5vcGFjaXR5KVxuICAgICAgaWYgKG8uY29sb3IgICAhPSBudWxsKSB0aGlzLmF0dHIoJ3N0b3AtY29sb3InLCBvLmNvbG9yKVxuICAgICAgaWYgKG8ub2Zmc2V0ICAhPSBudWxsKSB0aGlzLmF0dHIoJ29mZnNldCcsIG8ub2Zmc2V0KVxuICAgIH1cblxuICAgIHJldHVybiB0aGlzXG4gIH1cbn0pXG5cblNWRy5CQm94ID0gU1ZHLmludmVudCh7XG4gIC8vIEluaXRpYWxpemVcbiAgY3JlYXRlOiBmdW5jdGlvbihlbGVtZW50KSB7XG4gICAgLy8gZ2V0IHZhbHVlcyBpZiBlbGVtZW50IGlzIGdpdmVuXG4gICAgaWYgKGVsZW1lbnQpIHtcbiAgICAgIHZhciBib3hcblxuICAgICAgLy8geWVzIHRoaXMgaXMgdWdseSwgYnV0IEZpcmVmb3ggY2FuIGJlIGEgYml0Y2ggd2hlbiBpdCBjb21lcyB0byBlbGVtZW50cyB0aGF0IGFyZSBub3QgeWV0IHJlbmRlcmVkXG4gICAgICB0cnkge1xuXG4gICAgICAgIC8vIHRoZSBlbGVtZW50IGlzIE5PVCBpbiB0aGUgZG9tLCB0aHJvdyBlcnJvclxuICAgICAgICBpZighZG9jdW1lbnQuZG9jdW1lbnRFbGVtZW50LmNvbnRhaW5zKGVsZW1lbnQubm9kZSkpIHRocm93IG5ldyBFeGNlcHRpb24oJ0VsZW1lbnQgbm90IGluIHRoZSBkb20nKVxuXG4gICAgICAgIC8vIGZpbmQgbmF0aXZlIGJib3hcbiAgICAgICAgYm94ID0gZWxlbWVudC5ub2RlLmdldEJCb3goKVxuICAgICAgfSBjYXRjaChlKSB7XG4gICAgICAgIGlmKGVsZW1lbnQgaW5zdGFuY2VvZiBTVkcuU2hhcGUpe1xuICAgICAgICAgIHZhciBjbG9uZSA9IGVsZW1lbnQuY2xvbmUoU1ZHLnBhcnNlci5kcmF3KS5zaG93KClcbiAgICAgICAgICBib3ggPSBjbG9uZS5iYm94KClcbiAgICAgICAgICBjbG9uZS5yZW1vdmUoKVxuICAgICAgICB9ZWxzZXtcbiAgICAgICAgICBib3ggPSB7XG4gICAgICAgICAgICB4OiAgICAgIGVsZW1lbnQubm9kZS5jbGllbnRMZWZ0XG4gICAgICAgICAgLCB5OiAgICAgIGVsZW1lbnQubm9kZS5jbGllbnRUb3BcbiAgICAgICAgICAsIHdpZHRoOiAgZWxlbWVudC5ub2RlLmNsaWVudFdpZHRoXG4gICAgICAgICAgLCBoZWlnaHQ6IGVsZW1lbnQubm9kZS5jbGllbnRIZWlnaHRcbiAgICAgICAgICB9XG4gICAgICAgIH1cbiAgICAgIH1cblxuICAgICAgLy8gcGxhaW4geCBhbmQgeVxuICAgICAgdGhpcy54ID0gYm94LnhcbiAgICAgIHRoaXMueSA9IGJveC55XG5cbiAgICAgIC8vIHBsYWluIHdpZHRoIGFuZCBoZWlnaHRcbiAgICAgIHRoaXMud2lkdGggID0gYm94LndpZHRoXG4gICAgICB0aGlzLmhlaWdodCA9IGJveC5oZWlnaHRcbiAgICB9XG5cbiAgICAvLyBhZGQgY2VudGVyLCByaWdodCBhbmQgYm90dG9tXG4gICAgZnVsbEJveCh0aGlzKVxuICB9XG5cbiAgLy8gRGVmaW5lIFBhcmVudFxuLCBwYXJlbnQ6IFNWRy5FbGVtZW50XG5cbiAgLy8gQ29uc3RydWN0b3JcbiwgY29uc3RydWN0OiB7XG4gICAgLy8gR2V0IGJvdW5kaW5nIGJveFxuICAgIGJib3g6IGZ1bmN0aW9uKCkge1xuICAgICAgcmV0dXJuIG5ldyBTVkcuQkJveCh0aGlzKVxuICAgIH1cbiAgfVxuXG59KVxuXG5TVkcuVEJveCA9IFNWRy5pbnZlbnQoe1xuICAvLyBJbml0aWFsaXplXG4gIGNyZWF0ZTogZnVuY3Rpb24oZWxlbWVudCkge1xuICAgIC8vIGdldCB2YWx1ZXMgaWYgZWxlbWVudCBpcyBnaXZlblxuICAgIGlmIChlbGVtZW50KSB7XG4gICAgICB2YXIgdCAgID0gZWxlbWVudC5jdG0oKS5leHRyYWN0KClcbiAgICAgICAgLCBib3ggPSBlbGVtZW50LmJib3goKVxuXG4gICAgICAvLyB3aWR0aCBhbmQgaGVpZ2h0IGluY2x1ZGluZyB0cmFuc2Zvcm1hdGlvbnNcbiAgICAgIHRoaXMud2lkdGggID0gYm94LndpZHRoICAqIHQuc2NhbGVYXG4gICAgICB0aGlzLmhlaWdodCA9IGJveC5oZWlnaHQgKiB0LnNjYWxlWVxuXG4gICAgICAvLyB4IGFuZCB5IGluY2x1ZGluZyB0cmFuc2Zvcm1hdGlvbnNcbiAgICAgIHRoaXMueCA9IGJveC54ICsgdC54XG4gICAgICB0aGlzLnkgPSBib3gueSArIHQueVxuICAgIH1cblxuICAgIC8vIGFkZCBjZW50ZXIsIHJpZ2h0IGFuZCBib3R0b21cbiAgICBmdWxsQm94KHRoaXMpXG4gIH1cblxuICAvLyBEZWZpbmUgUGFyZW50XG4sIHBhcmVudDogU1ZHLkVsZW1lbnRcblxuICAvLyBDb25zdHJ1Y3RvclxuLCBjb25zdHJ1Y3Q6IHtcbiAgICAvLyBHZXQgdHJhbnNmb3JtZWQgYm91bmRpbmcgYm94XG4gICAgdGJveDogZnVuY3Rpb24oKSB7XG4gICAgICByZXR1cm4gbmV3IFNWRy5UQm94KHRoaXMpXG4gICAgfVxuICB9XG5cbn0pXG5cblxuU1ZHLlJCb3ggPSBTVkcuaW52ZW50KHtcbiAgLy8gSW5pdGlhbGl6ZVxuICBjcmVhdGU6IGZ1bmN0aW9uKGVsZW1lbnQpIHtcbiAgICBpZiAoZWxlbWVudCkge1xuICAgICAgdmFyIGUgICAgPSBlbGVtZW50LmRvYygpLnBhcmVudCgpXG4gICAgICAgICwgYm94ICA9IGVsZW1lbnQubm9kZS5nZXRCb3VuZGluZ0NsaWVudFJlY3QoKVxuICAgICAgICAsIHpvb20gPSAxXG5cbiAgICAgIC8vIGdldCBzY3JlZW4gb2Zmc2V0XG4gICAgICB0aGlzLnggPSBib3gubGVmdFxuICAgICAgdGhpcy55ID0gYm94LnRvcFxuXG4gICAgICAvLyBzdWJ0cmFjdCBwYXJlbnQgb2Zmc2V0XG4gICAgICB0aGlzLnggLT0gZS5vZmZzZXRMZWZ0XG4gICAgICB0aGlzLnkgLT0gZS5vZmZzZXRUb3BcblxuICAgICAgd2hpbGUgKGUgPSBlLm9mZnNldFBhcmVudCkge1xuICAgICAgICB0aGlzLnggLT0gZS5vZmZzZXRMZWZ0XG4gICAgICAgIHRoaXMueSAtPSBlLm9mZnNldFRvcFxuICAgICAgfVxuXG4gICAgICAvLyBjYWxjdWxhdGUgY3VtdWxhdGl2ZSB6b29tIGZyb20gc3ZnIGRvY3VtZW50c1xuICAgICAgZSA9IGVsZW1lbnRcbiAgICAgIHdoaWxlIChlLnBhcmVudCAmJiAoZSA9IGUucGFyZW50KCkpKSB7XG4gICAgICAgIGlmIChlLnZpZXdib3gpIHtcbiAgICAgICAgICB6b29tICo9IGUudmlld2JveCgpLnpvb21cbiAgICAgICAgICB0aGlzLnggLT0gZS54KCkgfHwgMFxuICAgICAgICAgIHRoaXMueSAtPSBlLnkoKSB8fCAwXG4gICAgICAgIH1cbiAgICAgIH1cblxuICAgICAgLy8gcmVjYWxjdWxhdGUgdmlld2JveCBkaXN0b3J0aW9uXG4gICAgICB0aGlzLndpZHRoICA9IGJveC53aWR0aCAgLz0gem9vbVxuICAgICAgdGhpcy5oZWlnaHQgPSBib3guaGVpZ2h0IC89IHpvb21cbiAgICB9XG5cbiAgICAvLyBhZGQgY2VudGVyLCByaWdodCBhbmQgYm90dG9tXG4gICAgZnVsbEJveCh0aGlzKVxuXG4gICAgLy8gb2Zmc2V0IGJ5IHdpbmRvdyBzY3JvbGwgcG9zaXRpb24sIGJlY2F1c2UgZ2V0Qm91bmRpbmdDbGllbnRSZWN0IGNoYW5nZXMgd2hlbiB3aW5kb3cgaXMgc2Nyb2xsZWRcbiAgICB0aGlzLnggKz0gd2luZG93LnBhZ2VYT2Zmc2V0XG4gICAgdGhpcy55ICs9IHdpbmRvdy5wYWdlWU9mZnNldFxuICB9XG5cbiAgLy8gZGVmaW5lIFBhcmVudFxuLCBwYXJlbnQ6IFNWRy5FbGVtZW50XG5cbiAgLy8gQ29uc3RydWN0b3JcbiwgY29uc3RydWN0OiB7XG4gICAgLy8gR2V0IHJlY3QgYm94XG4gICAgcmJveDogZnVuY3Rpb24oKSB7XG4gICAgICByZXR1cm4gbmV3IFNWRy5SQm94KHRoaXMpXG4gICAgfVxuICB9XG5cbn0pXG5cbi8vIEFkZCB1bml2ZXJzYWwgbWVyZ2UgbWV0aG9kXG47W1NWRy5CQm94LCBTVkcuVEJveCwgU1ZHLlJCb3hdLmZvckVhY2goZnVuY3Rpb24oYykge1xuXG4gIFNWRy5leHRlbmQoYywge1xuICAgIC8vIE1lcmdlIHJlY3QgYm94IHdpdGggYW5vdGhlciwgcmV0dXJuIGEgbmV3IGluc3RhbmNlXG4gICAgbWVyZ2U6IGZ1bmN0aW9uKGJveCkge1xuICAgICAgdmFyIGIgPSBuZXcgYygpXG5cbiAgICAgIC8vIG1lcmdlIGJveGVzXG4gICAgICBiLnggICAgICA9IE1hdGgubWluKHRoaXMueCwgYm94LngpXG4gICAgICBiLnkgICAgICA9IE1hdGgubWluKHRoaXMueSwgYm94LnkpXG4gICAgICBiLndpZHRoICA9IE1hdGgubWF4KHRoaXMueCArIHRoaXMud2lkdGgsICBib3gueCArIGJveC53aWR0aCkgIC0gYi54XG4gICAgICBiLmhlaWdodCA9IE1hdGgubWF4KHRoaXMueSArIHRoaXMuaGVpZ2h0LCBib3gueSArIGJveC5oZWlnaHQpIC0gYi55XG5cbiAgICAgIHJldHVybiBmdWxsQm94KGIpXG4gICAgfVxuXG4gIH0pXG5cbn0pXG5cblNWRy5NYXRyaXggPSBTVkcuaW52ZW50KHtcbiAgLy8gSW5pdGlhbGl6ZVxuICBjcmVhdGU6IGZ1bmN0aW9uKHNvdXJjZSkge1xuICAgIHZhciBpLCBiYXNlID0gYXJyYXlUb01hdHJpeChbMSwgMCwgMCwgMSwgMCwgMF0pXG5cbiAgICAvLyBlbnN1cmUgc291cmNlIGFzIG9iamVjdFxuICAgIHNvdXJjZSA9IHNvdXJjZSBpbnN0YW5jZW9mIFNWRy5FbGVtZW50ID9cbiAgICAgIHNvdXJjZS5tYXRyaXhpZnkoKSA6XG4gICAgdHlwZW9mIHNvdXJjZSA9PT0gJ3N0cmluZycgP1xuICAgICAgc3RyaW5nVG9NYXRyaXgoc291cmNlKSA6XG4gICAgYXJndW1lbnRzLmxlbmd0aCA9PSA2ID9cbiAgICAgIGFycmF5VG9NYXRyaXgoW10uc2xpY2UuY2FsbChhcmd1bWVudHMpKSA6XG4gICAgdHlwZW9mIHNvdXJjZSA9PT0gJ29iamVjdCcgP1xuICAgICAgc291cmNlIDogYmFzZVxuXG4gICAgLy8gbWVyZ2Ugc291cmNlXG4gICAgZm9yIChpID0gYWJjZGVmLmxlbmd0aCAtIDE7IGkgPj0gMDsgLS1pKVxuICAgICAgdGhpc1thYmNkZWZbaV1dID0gc291cmNlICYmIHR5cGVvZiBzb3VyY2VbYWJjZGVmW2ldXSA9PT0gJ251bWJlcicgP1xuICAgICAgICBzb3VyY2VbYWJjZGVmW2ldXSA6IGJhc2VbYWJjZGVmW2ldXVxuICB9XG5cbiAgLy8gQWRkIG1ldGhvZHNcbiwgZXh0ZW5kOiB7XG4gICAgLy8gRXh0cmFjdCBpbmRpdmlkdWFsIHRyYW5zZm9ybWF0aW9uc1xuICAgIGV4dHJhY3Q6IGZ1bmN0aW9uKCkge1xuICAgICAgLy8gZmluZCBkZWx0YSB0cmFuc2Zvcm0gcG9pbnRzXG4gICAgICB2YXIgcHggICAgPSBkZWx0YVRyYW5zZm9ybVBvaW50KHRoaXMsIDAsIDEpXG4gICAgICAgICwgcHkgICAgPSBkZWx0YVRyYW5zZm9ybVBvaW50KHRoaXMsIDEsIDApXG4gICAgICAgICwgc2tld1ggPSAxODAgLyBNYXRoLlBJICogTWF0aC5hdGFuMihweC55LCBweC54KSAtIDkwXG5cbiAgICAgIHJldHVybiB7XG4gICAgICAgIC8vIHRyYW5zbGF0aW9uXG4gICAgICAgIHg6ICAgICAgICB0aGlzLmVcbiAgICAgICwgeTogICAgICAgIHRoaXMuZlxuICAgICAgLCB0cmFuc2Zvcm1lZFg6KHRoaXMuZSAqIE1hdGguY29zKHNrZXdYICogTWF0aC5QSSAvIDE4MCkgKyB0aGlzLmYgKiBNYXRoLnNpbihza2V3WCAqIE1hdGguUEkgLyAxODApKSAvIE1hdGguc3FydCh0aGlzLmEgKiB0aGlzLmEgKyB0aGlzLmIgKiB0aGlzLmIpXG4gICAgICAsIHRyYW5zZm9ybWVkWToodGhpcy5mICogTWF0aC5jb3Moc2tld1ggKiBNYXRoLlBJIC8gMTgwKSArIHRoaXMuZSAqIE1hdGguc2luKC1za2V3WCAqIE1hdGguUEkgLyAxODApKSAvIE1hdGguc3FydCh0aGlzLmMgKiB0aGlzLmMgKyB0aGlzLmQgKiB0aGlzLmQpXG4gICAgICAgIC8vIHNrZXdcbiAgICAgICwgc2tld1g6ICAgIC1za2V3WFxuICAgICAgLCBza2V3WTogICAgMTgwIC8gTWF0aC5QSSAqIE1hdGguYXRhbjIocHkueSwgcHkueClcbiAgICAgICAgLy8gc2NhbGVcbiAgICAgICwgc2NhbGVYOiAgIE1hdGguc3FydCh0aGlzLmEgKiB0aGlzLmEgKyB0aGlzLmIgKiB0aGlzLmIpXG4gICAgICAsIHNjYWxlWTogICBNYXRoLnNxcnQodGhpcy5jICogdGhpcy5jICsgdGhpcy5kICogdGhpcy5kKVxuICAgICAgICAvLyByb3RhdGlvblxuICAgICAgLCByb3RhdGlvbjogc2tld1hcbiAgICAgICwgYTogdGhpcy5hXG4gICAgICAsIGI6IHRoaXMuYlxuICAgICAgLCBjOiB0aGlzLmNcbiAgICAgICwgZDogdGhpcy5kXG4gICAgICAsIGU6IHRoaXMuZVxuICAgICAgLCBmOiB0aGlzLmZcbiAgICAgICwgbWF0cml4OiBuZXcgU1ZHLk1hdHJpeCh0aGlzKVxuICAgICAgfVxuICAgIH1cbiAgICAvLyBDbG9uZSBtYXRyaXhcbiAgLCBjbG9uZTogZnVuY3Rpb24oKSB7XG4gICAgICByZXR1cm4gbmV3IFNWRy5NYXRyaXgodGhpcylcbiAgICB9XG4gICAgLy8gTW9ycGggb25lIG1hdHJpeCBpbnRvIGFub3RoZXJcbiAgLCBtb3JwaDogZnVuY3Rpb24obWF0cml4KSB7XG4gICAgICAvLyBzdG9yZSBuZXcgZGVzdGluYXRpb25cbiAgICAgIHRoaXMuZGVzdGluYXRpb24gPSBuZXcgU1ZHLk1hdHJpeChtYXRyaXgpXG5cbiAgICAgIHJldHVybiB0aGlzXG4gICAgfVxuICAgIC8vIEdldCBtb3JwaGVkIG1hdHJpeCBhdCBhIGdpdmVuIHBvc2l0aW9uXG4gICwgYXQ6IGZ1bmN0aW9uKHBvcykge1xuICAgICAgLy8gbWFrZSBzdXJlIGEgZGVzdGluYXRpb24gaXMgZGVmaW5lZFxuICAgICAgaWYgKCF0aGlzLmRlc3RpbmF0aW9uKSByZXR1cm4gdGhpc1xuXG4gICAgICAvLyBjYWxjdWxhdGUgbW9ycGhlZCBtYXRyaXggYXQgYSBnaXZlbiBwb3NpdGlvblxuICAgICAgdmFyIG1hdHJpeCA9IG5ldyBTVkcuTWF0cml4KHtcbiAgICAgICAgYTogdGhpcy5hICsgKHRoaXMuZGVzdGluYXRpb24uYSAtIHRoaXMuYSkgKiBwb3NcbiAgICAgICwgYjogdGhpcy5iICsgKHRoaXMuZGVzdGluYXRpb24uYiAtIHRoaXMuYikgKiBwb3NcbiAgICAgICwgYzogdGhpcy5jICsgKHRoaXMuZGVzdGluYXRpb24uYyAtIHRoaXMuYykgKiBwb3NcbiAgICAgICwgZDogdGhpcy5kICsgKHRoaXMuZGVzdGluYXRpb24uZCAtIHRoaXMuZCkgKiBwb3NcbiAgICAgICwgZTogdGhpcy5lICsgKHRoaXMuZGVzdGluYXRpb24uZSAtIHRoaXMuZSkgKiBwb3NcbiAgICAgICwgZjogdGhpcy5mICsgKHRoaXMuZGVzdGluYXRpb24uZiAtIHRoaXMuZikgKiBwb3NcbiAgICAgIH0pXG5cbiAgICAgIC8vIHByb2Nlc3MgcGFyYW1ldHJpYyByb3RhdGlvbiBpZiBwcmVzZW50XG4gICAgICBpZiAodGhpcy5wYXJhbSAmJiB0aGlzLnBhcmFtLnRvKSB7XG4gICAgICAgIC8vIGNhbGN1bGF0ZSBjdXJyZW50IHBhcmFtZXRyaWMgcG9zaXRpb25cbiAgICAgICAgdmFyIHBhcmFtID0ge1xuICAgICAgICAgIHJvdGF0aW9uOiB0aGlzLnBhcmFtLmZyb20ucm90YXRpb24gKyAodGhpcy5wYXJhbS50by5yb3RhdGlvbiAtIHRoaXMucGFyYW0uZnJvbS5yb3RhdGlvbikgKiBwb3NcbiAgICAgICAgLCBjeDogICAgICAgdGhpcy5wYXJhbS5mcm9tLmN4XG4gICAgICAgICwgY3k6ICAgICAgIHRoaXMucGFyYW0uZnJvbS5jeVxuICAgICAgICB9XG5cbiAgICAgICAgLy8gcm90YXRlIG1hdHJpeFxuICAgICAgICBtYXRyaXggPSBtYXRyaXgucm90YXRlKFxuICAgICAgICAgICh0aGlzLnBhcmFtLnRvLnJvdGF0aW9uIC0gdGhpcy5wYXJhbS5mcm9tLnJvdGF0aW9uICogMikgKiBwb3NcbiAgICAgICAgLCBwYXJhbS5jeFxuICAgICAgICAsIHBhcmFtLmN5XG4gICAgICAgIClcblxuICAgICAgICAvLyBzdG9yZSBjdXJyZW50IHBhcmFtZXRyaWMgdmFsdWVzXG4gICAgICAgIG1hdHJpeC5wYXJhbSA9IHBhcmFtXG4gICAgICB9XG5cbiAgICAgIHJldHVybiBtYXRyaXhcbiAgICB9XG4gICAgLy8gTXVsdGlwbGllcyBieSBnaXZlbiBtYXRyaXhcbiAgLCBtdWx0aXBseTogZnVuY3Rpb24obWF0cml4KSB7XG4gICAgICByZXR1cm4gbmV3IFNWRy5NYXRyaXgodGhpcy5uYXRpdmUoKS5tdWx0aXBseShwYXJzZU1hdHJpeChtYXRyaXgpLm5hdGl2ZSgpKSlcbiAgICB9XG4gICAgLy8gSW52ZXJzZXMgbWF0cml4XG4gICwgaW52ZXJzZTogZnVuY3Rpb24oKSB7XG4gICAgICByZXR1cm4gbmV3IFNWRy5NYXRyaXgodGhpcy5uYXRpdmUoKS5pbnZlcnNlKCkpXG4gICAgfVxuICAgIC8vIFRyYW5zbGF0ZSBtYXRyaXhcbiAgLCB0cmFuc2xhdGU6IGZ1bmN0aW9uKHgsIHkpIHtcbiAgICAgIHJldHVybiBuZXcgU1ZHLk1hdHJpeCh0aGlzLm5hdGl2ZSgpLnRyYW5zbGF0ZSh4IHx8IDAsIHkgfHwgMCkpXG4gICAgfVxuICAgIC8vIFNjYWxlIG1hdHJpeFxuICAsIHNjYWxlOiBmdW5jdGlvbih4LCB5LCBjeCwgY3kpIHtcbiAgICAgIC8vIHN1cHBvcnQgdW5pZm9ybWFsIHNjYWxlXG4gICAgICBpZiAoYXJndW1lbnRzLmxlbmd0aCA9PSAxKSB7XG4gICAgICAgIHkgPSB4XG4gICAgICB9IGVsc2UgaWYgKGFyZ3VtZW50cy5sZW5ndGggPT0gMykge1xuICAgICAgICBjeSA9IGN4XG4gICAgICAgIGN4ID0geVxuICAgICAgICB5ID0geFxuICAgICAgfVxuXG4gICAgICByZXR1cm4gdGhpcy5hcm91bmQoY3gsIGN5LCBuZXcgU1ZHLk1hdHJpeCh4LCAwLCAwLCB5LCAwLCAwKSlcbiAgICB9XG4gICAgLy8gUm90YXRlIG1hdHJpeFxuICAsIHJvdGF0ZTogZnVuY3Rpb24ociwgY3gsIGN5KSB7XG4gICAgICAvLyBjb252ZXJ0IGRlZ3JlZXMgdG8gcmFkaWFuc1xuICAgICAgciA9IFNWRy51dGlscy5yYWRpYW5zKHIpXG5cbiAgICAgIHJldHVybiB0aGlzLmFyb3VuZChjeCwgY3ksIG5ldyBTVkcuTWF0cml4KE1hdGguY29zKHIpLCBNYXRoLnNpbihyKSwgLU1hdGguc2luKHIpLCBNYXRoLmNvcyhyKSwgMCwgMCkpXG4gICAgfVxuICAgIC8vIEZsaXAgbWF0cml4IG9uIHggb3IgeSwgYXQgYSBnaXZlbiBvZmZzZXRcbiAgLCBmbGlwOiBmdW5jdGlvbihhLCBvKSB7XG4gICAgICByZXR1cm4gYSA9PSAneCcgPyB0aGlzLnNjYWxlKC0xLCAxLCBvLCAwKSA6IHRoaXMuc2NhbGUoMSwgLTEsIDAsIG8pXG4gICAgfVxuICAgIC8vIFNrZXdcbiAgLCBza2V3OiBmdW5jdGlvbih4LCB5LCBjeCwgY3kpIHtcbiAgICAgIC8vIHN1cHBvcnQgdW5pZm9ybWFsIHNrZXdcbiAgICAgIGlmIChhcmd1bWVudHMubGVuZ3RoID09IDEpIHtcbiAgICAgICAgeSA9IHhcbiAgICAgIH0gZWxzZSBpZiAoYXJndW1lbnRzLmxlbmd0aCA9PSAzKSB7XG4gICAgICAgIGN5ID0gY3hcbiAgICAgICAgY3ggPSB5XG4gICAgICAgIHkgPSB4XG4gICAgICB9XG5cbiAgICAgIC8vIGNvbnZlcnQgZGVncmVlcyB0byByYWRpYW5zXG4gICAgICB4ID0gU1ZHLnV0aWxzLnJhZGlhbnMoeClcbiAgICAgIHkgPSBTVkcudXRpbHMucmFkaWFucyh5KVxuXG4gICAgICByZXR1cm4gdGhpcy5hcm91bmQoY3gsIGN5LCBuZXcgU1ZHLk1hdHJpeCgxLCBNYXRoLnRhbih5KSwgTWF0aC50YW4oeCksIDEsIDAsIDApKVxuICAgIH1cbiAgICAvLyBTa2V3WFxuICAsIHNrZXdYOiBmdW5jdGlvbih4LCBjeCwgY3kpIHtcbiAgICAgIHJldHVybiB0aGlzLnNrZXcoeCwgMCwgY3gsIGN5KVxuICAgIH1cbiAgICAvLyBTa2V3WVxuICAsIHNrZXdZOiBmdW5jdGlvbih5LCBjeCwgY3kpIHtcbiAgICAgIHJldHVybiB0aGlzLnNrZXcoMCwgeSwgY3gsIGN5KVxuICAgIH1cbiAgICAvLyBUcmFuc2Zvcm0gYXJvdW5kIGEgY2VudGVyIHBvaW50XG4gICwgYXJvdW5kOiBmdW5jdGlvbihjeCwgY3ksIG1hdHJpeCkge1xuICAgICAgcmV0dXJuIHRoaXNcbiAgICAgICAgLm11bHRpcGx5KG5ldyBTVkcuTWF0cml4KDEsIDAsIDAsIDEsIGN4IHx8IDAsIGN5IHx8IDApKVxuICAgICAgICAubXVsdGlwbHkobWF0cml4KVxuICAgICAgICAubXVsdGlwbHkobmV3IFNWRy5NYXRyaXgoMSwgMCwgMCwgMSwgLWN4IHx8IDAsIC1jeSB8fCAwKSlcbiAgICB9XG4gICAgLy8gQ29udmVydCB0byBuYXRpdmUgU1ZHTWF0cml4XG4gICwgbmF0aXZlOiBmdW5jdGlvbigpIHtcbiAgICAgIC8vIGNyZWF0ZSBuZXcgbWF0cml4XG4gICAgICB2YXIgbWF0cml4ID0gU1ZHLnBhcnNlci5uYXRpdmUuY3JlYXRlU1ZHTWF0cml4KClcblxuICAgICAgLy8gdXBkYXRlIHdpdGggY3VycmVudCB2YWx1ZXNcbiAgICAgIGZvciAodmFyIGkgPSBhYmNkZWYubGVuZ3RoIC0gMTsgaSA+PSAwOyBpLS0pXG4gICAgICAgIG1hdHJpeFthYmNkZWZbaV1dID0gdGhpc1thYmNkZWZbaV1dXG5cbiAgICAgIHJldHVybiBtYXRyaXhcbiAgICB9XG4gICAgLy8gQ29udmVydCBtYXRyaXggdG8gc3RyaW5nXG4gICwgdG9TdHJpbmc6IGZ1bmN0aW9uKCkge1xuICAgICAgcmV0dXJuICdtYXRyaXgoJyArIHRoaXMuYSArICcsJyArIHRoaXMuYiArICcsJyArIHRoaXMuYyArICcsJyArIHRoaXMuZCArICcsJyArIHRoaXMuZSArICcsJyArIHRoaXMuZiArICcpJ1xuICAgIH1cbiAgfVxuXG4gIC8vIERlZmluZSBwYXJlbnRcbiwgcGFyZW50OiBTVkcuRWxlbWVudFxuXG4gIC8vIEFkZCBwYXJlbnQgbWV0aG9kXG4sIGNvbnN0cnVjdDoge1xuICAgIC8vIEdldCBjdXJyZW50IG1hdHJpeFxuICAgIGN0bTogZnVuY3Rpb24oKSB7XG4gICAgICByZXR1cm4gbmV3IFNWRy5NYXRyaXgodGhpcy5ub2RlLmdldENUTSgpKVxuICAgIH0sXG4gICAgLy8gR2V0IGN1cnJlbnQgc2NyZWVuIG1hdHJpeFxuICAgIHNjcmVlbkNUTTogZnVuY3Rpb24oKSB7XG4gICAgICByZXR1cm4gbmV3IFNWRy5NYXRyaXgodGhpcy5ub2RlLmdldFNjcmVlbkNUTSgpKVxuICAgIH1cblxuICB9XG5cbn0pXG5cblNWRy5Qb2ludCA9IFNWRy5pbnZlbnQoe1xuICAvLyBJbml0aWFsaXplXG4gIGNyZWF0ZTogZnVuY3Rpb24oeCx5KSB7XG4gICAgdmFyIGksIHNvdXJjZSwgYmFzZSA9IHt4OjAsIHk6MH1cblxuICAgIC8vIGVuc3VyZSBzb3VyY2UgYXMgb2JqZWN0XG4gICAgc291cmNlID0gQXJyYXkuaXNBcnJheSh4KSA/XG4gICAgICB7eDp4WzBdLCB5OnhbMV19IDpcbiAgICB0eXBlb2YgeCA9PT0gJ29iamVjdCcgP1xuICAgICAge3g6eC54LCB5OngueX0gOlxuICAgIHggIT0gbnVsbCA/XG4gICAgICB7eDp4LCB5Oih5ICE9IG51bGwgPyB5IDogeCl9IDogYmFzZSAvLyBJZiB5IGhhcyBubyB2YWx1ZSwgdGhlbiB4IGlzIHVzZWQgaGFzIGl0cyB2YWx1ZVxuXG4gICAgLy8gbWVyZ2Ugc291cmNlXG4gICAgdGhpcy54ID0gc291cmNlLnhcbiAgICB0aGlzLnkgPSBzb3VyY2UueVxuICB9XG5cbiAgLy8gQWRkIG1ldGhvZHNcbiwgZXh0ZW5kOiB7XG4gICAgLy8gQ2xvbmUgcG9pbnRcbiAgICBjbG9uZTogZnVuY3Rpb24oKSB7XG4gICAgICByZXR1cm4gbmV3IFNWRy5Qb2ludCh0aGlzKVxuICAgIH1cbiAgICAvLyBNb3JwaCBvbmUgcG9pbnQgaW50byBhbm90aGVyXG4gICwgbW9ycGg6IGZ1bmN0aW9uKHgsIHkpIHtcbiAgICAgIC8vIHN0b3JlIG5ldyBkZXN0aW5hdGlvblxuICAgICAgdGhpcy5kZXN0aW5hdGlvbiA9IG5ldyBTVkcuUG9pbnQoeCwgeSlcblxuICAgICAgcmV0dXJuIHRoaXNcbiAgICB9XG4gICAgLy8gR2V0IG1vcnBoZWQgcG9pbnQgYXQgYSBnaXZlbiBwb3NpdGlvblxuICAsIGF0OiBmdW5jdGlvbihwb3MpIHtcbiAgICAgIC8vIG1ha2Ugc3VyZSBhIGRlc3RpbmF0aW9uIGlzIGRlZmluZWRcbiAgICAgIGlmICghdGhpcy5kZXN0aW5hdGlvbikgcmV0dXJuIHRoaXNcblxuICAgICAgLy8gY2FsY3VsYXRlIG1vcnBoZWQgbWF0cml4IGF0IGEgZ2l2ZW4gcG9zaXRpb25cbiAgICAgIHZhciBwb2ludCA9IG5ldyBTVkcuUG9pbnQoe1xuICAgICAgICB4OiB0aGlzLnggKyAodGhpcy5kZXN0aW5hdGlvbi54IC0gdGhpcy54KSAqIHBvc1xuICAgICAgLCB5OiB0aGlzLnkgKyAodGhpcy5kZXN0aW5hdGlvbi55IC0gdGhpcy55KSAqIHBvc1xuICAgICAgfSlcblxuICAgICAgcmV0dXJuIHBvaW50XG4gICAgfVxuICAgIC8vIENvbnZlcnQgdG8gbmF0aXZlIFNWR1BvaW50XG4gICwgbmF0aXZlOiBmdW5jdGlvbigpIHtcbiAgICAgIC8vIGNyZWF0ZSBuZXcgcG9pbnRcbiAgICAgIHZhciBwb2ludCA9IFNWRy5wYXJzZXIubmF0aXZlLmNyZWF0ZVNWR1BvaW50KClcblxuICAgICAgLy8gdXBkYXRlIHdpdGggY3VycmVudCB2YWx1ZXNcbiAgICAgIHBvaW50LnggPSB0aGlzLnhcbiAgICAgIHBvaW50LnkgPSB0aGlzLnlcblxuICAgICAgcmV0dXJuIHBvaW50XG4gICAgfVxuICAgIC8vIHRyYW5zZm9ybSBwb2ludCB3aXRoIG1hdHJpeFxuICAsIHRyYW5zZm9ybTogZnVuY3Rpb24obWF0cml4KSB7XG4gICAgICByZXR1cm4gbmV3IFNWRy5Qb2ludCh0aGlzLm5hdGl2ZSgpLm1hdHJpeFRyYW5zZm9ybShtYXRyaXgubmF0aXZlKCkpKVxuICAgIH1cblxuICB9XG5cbn0pXG5cblNWRy5leHRlbmQoU1ZHLkVsZW1lbnQsIHtcblxuICAvLyBHZXQgcG9pbnRcbiAgcG9pbnQ6IGZ1bmN0aW9uKHgsIHkpIHtcbiAgICByZXR1cm4gbmV3IFNWRy5Qb2ludCh4LHkpLnRyYW5zZm9ybSh0aGlzLnNjcmVlbkNUTSgpLmludmVyc2UoKSk7XG4gIH1cblxufSlcblxuU1ZHLmV4dGVuZChTVkcuRWxlbWVudCwge1xuICAvLyBTZXQgc3ZnIGVsZW1lbnQgYXR0cmlidXRlXG4gIGF0dHI6IGZ1bmN0aW9uKGEsIHYsIG4pIHtcbiAgICAvLyBhY3QgYXMgZnVsbCBnZXR0ZXJcbiAgICBpZiAoYSA9PSBudWxsKSB7XG4gICAgICAvLyBnZXQgYW4gb2JqZWN0IG9mIGF0dHJpYnV0ZXNcbiAgICAgIGEgPSB7fVxuICAgICAgdiA9IHRoaXMubm9kZS5hdHRyaWJ1dGVzXG4gICAgICBmb3IgKG4gPSB2Lmxlbmd0aCAtIDE7IG4gPj0gMDsgbi0tKVxuICAgICAgICBhW3Zbbl0ubm9kZU5hbWVdID0gU1ZHLnJlZ2V4LmlzTnVtYmVyLnRlc3QodltuXS5ub2RlVmFsdWUpID8gcGFyc2VGbG9hdCh2W25dLm5vZGVWYWx1ZSkgOiB2W25dLm5vZGVWYWx1ZVxuXG4gICAgICByZXR1cm4gYVxuXG4gICAgfSBlbHNlIGlmICh0eXBlb2YgYSA9PSAnb2JqZWN0Jykge1xuICAgICAgLy8gYXBwbHkgZXZlcnkgYXR0cmlidXRlIGluZGl2aWR1YWxseSBpZiBhbiBvYmplY3QgaXMgcGFzc2VkXG4gICAgICBmb3IgKHYgaW4gYSkgdGhpcy5hdHRyKHYsIGFbdl0pXG5cbiAgICB9IGVsc2UgaWYgKHYgPT09IG51bGwpIHtcbiAgICAgICAgLy8gcmVtb3ZlIHZhbHVlXG4gICAgICAgIHRoaXMubm9kZS5yZW1vdmVBdHRyaWJ1dGUoYSlcblxuICAgIH0gZWxzZSBpZiAodiA9PSBudWxsKSB7XG4gICAgICAvLyBhY3QgYXMgYSBnZXR0ZXIgaWYgdGhlIGZpcnN0IGFuZCBvbmx5IGFyZ3VtZW50IGlzIG5vdCBhbiBvYmplY3RcbiAgICAgIHYgPSB0aGlzLm5vZGUuZ2V0QXR0cmlidXRlKGEpXG4gICAgICByZXR1cm4gdiA9PSBudWxsID9cbiAgICAgICAgU1ZHLmRlZmF1bHRzLmF0dHJzW2FdIDpcbiAgICAgIFNWRy5yZWdleC5pc051bWJlci50ZXN0KHYpID9cbiAgICAgICAgcGFyc2VGbG9hdCh2KSA6IHZcblxuICAgIH0gZWxzZSB7XG4gICAgICAvLyBCVUcgRklYOiBzb21lIGJyb3dzZXJzIHdpbGwgcmVuZGVyIGEgc3Ryb2tlIGlmIGEgY29sb3IgaXMgZ2l2ZW4gZXZlbiB0aG91Z2ggc3Ryb2tlIHdpZHRoIGlzIDBcbiAgICAgIGlmIChhID09ICdzdHJva2Utd2lkdGgnKVxuICAgICAgICB0aGlzLmF0dHIoJ3N0cm9rZScsIHBhcnNlRmxvYXQodikgPiAwID8gdGhpcy5fc3Ryb2tlIDogbnVsbClcbiAgICAgIGVsc2UgaWYgKGEgPT0gJ3N0cm9rZScpXG4gICAgICAgIHRoaXMuX3N0cm9rZSA9IHZcblxuICAgICAgLy8gY29udmVydCBpbWFnZSBmaWxsIGFuZCBzdHJva2UgdG8gcGF0dGVybnNcbiAgICAgIGlmIChhID09ICdmaWxsJyB8fCBhID09ICdzdHJva2UnKSB7XG4gICAgICAgIGlmIChTVkcucmVnZXguaXNJbWFnZS50ZXN0KHYpKVxuICAgICAgICAgIHYgPSB0aGlzLmRvYygpLmRlZnMoKS5pbWFnZSh2LCAwLCAwKVxuXG4gICAgICAgIGlmICh2IGluc3RhbmNlb2YgU1ZHLkltYWdlKVxuICAgICAgICAgIHYgPSB0aGlzLmRvYygpLmRlZnMoKS5wYXR0ZXJuKDAsIDAsIGZ1bmN0aW9uKCkge1xuICAgICAgICAgICAgdGhpcy5hZGQodilcbiAgICAgICAgICB9KVxuICAgICAgfVxuXG4gICAgICAvLyBlbnN1cmUgY29ycmVjdCBudW1lcmljIHZhbHVlcyAoYWxzbyBhY2NlcHRzIE5hTiBhbmQgSW5maW5pdHkpXG4gICAgICBpZiAodHlwZW9mIHYgPT09ICdudW1iZXInKVxuICAgICAgICB2ID0gbmV3IFNWRy5OdW1iZXIodilcblxuICAgICAgLy8gZW5zdXJlIGZ1bGwgaGV4IGNvbG9yXG4gICAgICBlbHNlIGlmIChTVkcuQ29sb3IuaXNDb2xvcih2KSlcbiAgICAgICAgdiA9IG5ldyBTVkcuQ29sb3IodilcblxuICAgICAgLy8gcGFyc2UgYXJyYXkgdmFsdWVzXG4gICAgICBlbHNlIGlmIChBcnJheS5pc0FycmF5KHYpKVxuICAgICAgICB2ID0gbmV3IFNWRy5BcnJheSh2KVxuXG4gICAgICAvLyBzdG9yZSBwYXJhbWV0cmljIHRyYW5zZm9ybWF0aW9uIHZhbHVlcyBsb2NhbGx5XG4gICAgICBlbHNlIGlmICh2IGluc3RhbmNlb2YgU1ZHLk1hdHJpeCAmJiB2LnBhcmFtKVxuICAgICAgICB0aGlzLnBhcmFtID0gdi5wYXJhbVxuXG4gICAgICAvLyBpZiB0aGUgcGFzc2VkIGF0dHJpYnV0ZSBpcyBsZWFkaW5nLi4uXG4gICAgICBpZiAoYSA9PSAnbGVhZGluZycpIHtcbiAgICAgICAgLy8gLi4uIGNhbGwgdGhlIGxlYWRpbmcgbWV0aG9kIGluc3RlYWRcbiAgICAgICAgaWYgKHRoaXMubGVhZGluZylcbiAgICAgICAgICB0aGlzLmxlYWRpbmcodilcbiAgICAgIH0gZWxzZSB7XG4gICAgICAgIC8vIHNldCBnaXZlbiBhdHRyaWJ1dGUgb24gbm9kZVxuICAgICAgICB0eXBlb2YgbiA9PT0gJ3N0cmluZycgP1xuICAgICAgICAgIHRoaXMubm9kZS5zZXRBdHRyaWJ1dGVOUyhuLCBhLCB2LnRvU3RyaW5nKCkpIDpcbiAgICAgICAgICB0aGlzLm5vZGUuc2V0QXR0cmlidXRlKGEsIHYudG9TdHJpbmcoKSlcbiAgICAgIH1cblxuICAgICAgLy8gcmVidWlsZCBpZiByZXF1aXJlZFxuICAgICAgaWYgKHRoaXMucmVidWlsZCAmJiAoYSA9PSAnZm9udC1zaXplJyB8fCBhID09ICd4JykpXG4gICAgICAgIHRoaXMucmVidWlsZChhLCB2KVxuICAgIH1cblxuICAgIHJldHVybiB0aGlzXG4gIH1cbn0pXG5TVkcuZXh0ZW5kKFNWRy5FbGVtZW50LCB7XG4gIC8vIEFkZCB0cmFuc2Zvcm1hdGlvbnNcbiAgdHJhbnNmb3JtOiBmdW5jdGlvbihvLCByZWxhdGl2ZSkge1xuICAgIC8vIGdldCB0YXJnZXQgaW4gY2FzZSBvZiB0aGUgZnggbW9kdWxlLCBvdGhlcndpc2UgcmVmZXJlbmNlIHRoaXNcbiAgICB2YXIgdGFyZ2V0ID0gdGhpc1xuICAgICAgLCBtYXRyaXhcblxuICAgIC8vIGFjdCBhcyBhIGdldHRlclxuICAgIGlmICh0eXBlb2YgbyAhPT0gJ29iamVjdCcpIHtcbiAgICAgIC8vIGdldCBjdXJyZW50IG1hdHJpeFxuICAgICAgbWF0cml4ID0gbmV3IFNWRy5NYXRyaXgodGFyZ2V0KS5leHRyYWN0KClcblxuICAgICAgcmV0dXJuIHR5cGVvZiBvID09PSAnc3RyaW5nJyA/IG1hdHJpeFtvXSA6IG1hdHJpeFxuICAgIH1cblxuICAgIC8vIGdldCBjdXJyZW50IG1hdHJpeFxuICAgIG1hdHJpeCA9IG5ldyBTVkcuTWF0cml4KHRhcmdldClcblxuICAgIC8vIGVuc3VyZSByZWxhdGl2ZSBmbGFnXG4gICAgcmVsYXRpdmUgPSAhIXJlbGF0aXZlIHx8ICEhby5yZWxhdGl2ZVxuXG4gICAgLy8gYWN0IG9uIG1hdHJpeFxuICAgIGlmIChvLmEgIT0gbnVsbCkge1xuICAgICAgbWF0cml4ID0gcmVsYXRpdmUgP1xuICAgICAgICAvLyByZWxhdGl2ZVxuICAgICAgICBtYXRyaXgubXVsdGlwbHkobmV3IFNWRy5NYXRyaXgobykpIDpcbiAgICAgICAgLy8gYWJzb2x1dGVcbiAgICAgICAgbmV3IFNWRy5NYXRyaXgobylcblxuICAgIC8vIGFjdCBvbiByb3RhdGlvblxuICAgIH0gZWxzZSBpZiAoby5yb3RhdGlvbiAhPSBudWxsKSB7XG4gICAgICAvLyBlbnN1cmUgY2VudHJlIHBvaW50XG4gICAgICBlbnN1cmVDZW50cmUobywgdGFyZ2V0KVxuXG4gICAgICAvLyBhcHBseSB0cmFuc2Zvcm1hdGlvblxuICAgICAgbWF0cml4ID0gcmVsYXRpdmUgP1xuICAgICAgICAvLyByZWxhdGl2ZVxuICAgICAgICBtYXRyaXgucm90YXRlKG8ucm90YXRpb24sIG8uY3gsIG8uY3kpIDpcbiAgICAgICAgLy8gYWJzb2x1dGVcbiAgICAgICAgbWF0cml4LnJvdGF0ZShvLnJvdGF0aW9uIC0gbWF0cml4LmV4dHJhY3QoKS5yb3RhdGlvbiwgby5jeCwgby5jeSlcblxuICAgIC8vIGFjdCBvbiBzY2FsZVxuICAgIH0gZWxzZSBpZiAoby5zY2FsZSAhPSBudWxsIHx8IG8uc2NhbGVYICE9IG51bGwgfHwgby5zY2FsZVkgIT0gbnVsbCkge1xuICAgICAgLy8gZW5zdXJlIGNlbnRyZSBwb2ludFxuICAgICAgZW5zdXJlQ2VudHJlKG8sIHRhcmdldClcblxuICAgICAgLy8gZW5zdXJlIHNjYWxlIHZhbHVlcyBvbiBib3RoIGF4ZXNcbiAgICAgIG8uc2NhbGVYID0gby5zY2FsZSAhPSBudWxsID8gby5zY2FsZSA6IG8uc2NhbGVYICE9IG51bGwgPyBvLnNjYWxlWCA6IDFcbiAgICAgIG8uc2NhbGVZID0gby5zY2FsZSAhPSBudWxsID8gby5zY2FsZSA6IG8uc2NhbGVZICE9IG51bGwgPyBvLnNjYWxlWSA6IDFcblxuICAgICAgaWYgKCFyZWxhdGl2ZSkge1xuICAgICAgICAvLyBhYnNvbHV0ZTsgbXVsdGlwbHkgaW52ZXJzZWQgdmFsdWVzXG4gICAgICAgIHZhciBlID0gbWF0cml4LmV4dHJhY3QoKVxuICAgICAgICBvLnNjYWxlWCA9IG8uc2NhbGVYICogMSAvIGUuc2NhbGVYXG4gICAgICAgIG8uc2NhbGVZID0gby5zY2FsZVkgKiAxIC8gZS5zY2FsZVlcbiAgICAgIH1cblxuICAgICAgbWF0cml4ID0gbWF0cml4LnNjYWxlKG8uc2NhbGVYLCBvLnNjYWxlWSwgby5jeCwgby5jeSlcblxuICAgIC8vIGFjdCBvbiBza2V3XG4gICAgfSBlbHNlIGlmIChvLnNrZXcgIT0gbnVsbCB8fCBvLnNrZXdYICE9IG51bGwgfHwgby5za2V3WSAhPSBudWxsKSB7XG4gICAgICAvLyBlbnN1cmUgY2VudHJlIHBvaW50XG4gICAgICBlbnN1cmVDZW50cmUobywgdGFyZ2V0KVxuXG4gICAgICAvLyBlbnN1cmUgc2tldyB2YWx1ZXMgb24gYm90aCBheGVzXG4gICAgICBvLnNrZXdYID0gby5za2V3ICE9IG51bGwgPyBvLnNrZXcgOiBvLnNrZXdYICE9IG51bGwgPyBvLnNrZXdYIDogMFxuICAgICAgby5za2V3WSA9IG8uc2tldyAhPSBudWxsID8gby5za2V3IDogby5za2V3WSAhPSBudWxsID8gby5za2V3WSA6IDBcblxuICAgICAgaWYgKCFyZWxhdGl2ZSkge1xuICAgICAgICAvLyBhYnNvbHV0ZTsgcmVzZXQgc2tldyB2YWx1ZXNcbiAgICAgICAgdmFyIGUgPSBtYXRyaXguZXh0cmFjdCgpXG4gICAgICAgIG1hdHJpeCA9IG1hdHJpeC5tdWx0aXBseShuZXcgU1ZHLk1hdHJpeCgpLnNrZXcoZS5za2V3WCwgZS5za2V3WSwgby5jeCwgby5jeSkuaW52ZXJzZSgpKVxuICAgICAgfVxuXG4gICAgICBtYXRyaXggPSBtYXRyaXguc2tldyhvLnNrZXdYLCBvLnNrZXdZLCBvLmN4LCBvLmN5KVxuXG4gICAgLy8gYWN0IG9uIGZsaXBcbiAgICB9IGVsc2UgaWYgKG8uZmxpcCkge1xuICAgICAgbWF0cml4ID0gbWF0cml4LmZsaXAoXG4gICAgICAgIG8uZmxpcFxuICAgICAgLCBvLm9mZnNldCA9PSBudWxsID8gdGFyZ2V0LmJib3goKVsnYycgKyBvLmZsaXBdIDogby5vZmZzZXRcbiAgICAgIClcblxuICAgIC8vIGFjdCBvbiB0cmFuc2xhdGVcbiAgICB9IGVsc2UgaWYgKG8ueCAhPSBudWxsIHx8IG8ueSAhPSBudWxsKSB7XG4gICAgICBpZiAocmVsYXRpdmUpIHtcbiAgICAgICAgLy8gcmVsYXRpdmVcbiAgICAgICAgbWF0cml4ID0gbWF0cml4LnRyYW5zbGF0ZShvLngsIG8ueSlcbiAgICAgIH0gZWxzZSB7XG4gICAgICAgIC8vIGFic29sdXRlXG4gICAgICAgIGlmIChvLnggIT0gbnVsbCkgbWF0cml4LmUgPSBvLnhcbiAgICAgICAgaWYgKG8ueSAhPSBudWxsKSBtYXRyaXguZiA9IG8ueVxuICAgICAgfVxuICAgIH1cblxuICAgIHJldHVybiB0aGlzLmF0dHIoJ3RyYW5zZm9ybScsIG1hdHJpeClcbiAgfVxufSlcblxuU1ZHLmV4dGVuZChTVkcuRlgsIHtcbiAgdHJhbnNmb3JtOiBmdW5jdGlvbihvLCByZWxhdGl2ZSkge1xuICAgIC8vIGdldCB0YXJnZXQgaW4gY2FzZSBvZiB0aGUgZnggbW9kdWxlLCBvdGhlcndpc2UgcmVmZXJlbmNlIHRoaXNcbiAgICB2YXIgdGFyZ2V0ID0gdGhpcy50YXJnZXQoKVxuICAgICAgLCBtYXRyaXhcblxuICAgIC8vIGFjdCBhcyBhIGdldHRlclxuICAgIGlmICh0eXBlb2YgbyAhPT0gJ29iamVjdCcpIHtcbiAgICAgIC8vIGdldCBjdXJyZW50IG1hdHJpeFxuICAgICAgbWF0cml4ID0gbmV3IFNWRy5NYXRyaXgodGFyZ2V0KS5leHRyYWN0KClcblxuICAgICAgcmV0dXJuIHR5cGVvZiBvID09PSAnc3RyaW5nJyA/IG1hdHJpeFtvXSA6IG1hdHJpeFxuICAgIH1cblxuICAgIC8vIGVuc3VyZSByZWxhdGl2ZSBmbGFnXG4gICAgcmVsYXRpdmUgPSAhIXJlbGF0aXZlIHx8ICEhby5yZWxhdGl2ZVxuXG4gICAgLy8gYWN0IG9uIG1hdHJpeFxuICAgIGlmIChvLmEgIT0gbnVsbCkge1xuICAgICAgbWF0cml4ID0gbmV3IFNWRy5NYXRyaXgobylcblxuICAgIC8vIGFjdCBvbiByb3RhdGlvblxuICAgIH0gZWxzZSBpZiAoby5yb3RhdGlvbiAhPSBudWxsKSB7XG4gICAgICAvLyBlbnN1cmUgY2VudHJlIHBvaW50XG4gICAgICBlbnN1cmVDZW50cmUobywgdGFyZ2V0KVxuXG4gICAgICAvLyBhcHBseSB0cmFuc2Zvcm1hdGlvblxuICAgICAgbWF0cml4ID0gbmV3IFNWRy5Sb3RhdGUoby5yb3RhdGlvbiwgby5jeCwgby5jeSlcblxuICAgIC8vIGFjdCBvbiBzY2FsZVxuICAgIH0gZWxzZSBpZiAoby5zY2FsZSAhPSBudWxsIHx8IG8uc2NhbGVYICE9IG51bGwgfHwgby5zY2FsZVkgIT0gbnVsbCkge1xuICAgICAgLy8gZW5zdXJlIGNlbnRyZSBwb2ludFxuICAgICAgZW5zdXJlQ2VudHJlKG8sIHRhcmdldClcblxuICAgICAgLy8gZW5zdXJlIHNjYWxlIHZhbHVlcyBvbiBib3RoIGF4ZXNcbiAgICAgIG8uc2NhbGVYID0gby5zY2FsZSAhPSBudWxsID8gby5zY2FsZSA6IG8uc2NhbGVYICE9IG51bGwgPyBvLnNjYWxlWCA6IDFcbiAgICAgIG8uc2NhbGVZID0gby5zY2FsZSAhPSBudWxsID8gby5zY2FsZSA6IG8uc2NhbGVZICE9IG51bGwgPyBvLnNjYWxlWSA6IDFcblxuICAgICAgbWF0cml4ID0gbmV3IFNWRy5TY2FsZShvLnNjYWxlWCwgby5zY2FsZVksIG8uY3gsIG8uY3kpXG5cbiAgICAvLyBhY3Qgb24gc2tld1xuICAgIH0gZWxzZSBpZiAoby5za2V3WCAhPSBudWxsIHx8IG8uc2tld1kgIT0gbnVsbCkge1xuICAgICAgLy8gZW5zdXJlIGNlbnRyZSBwb2ludFxuICAgICAgZW5zdXJlQ2VudHJlKG8sIHRhcmdldClcblxuICAgICAgLy8gZW5zdXJlIHNrZXcgdmFsdWVzIG9uIGJvdGggYXhlc1xuICAgICAgby5za2V3WCA9IG8uc2tld1ggIT0gbnVsbCA/IG8uc2tld1ggOiAwXG4gICAgICBvLnNrZXdZID0gby5za2V3WSAhPSBudWxsID8gby5za2V3WSA6IDBcblxuICAgICAgbWF0cml4ID0gbmV3IFNWRy5Ta2V3KG8uc2tld1gsIG8uc2tld1ksIG8uY3gsIG8uY3kpXG5cbiAgICAvLyBhY3Qgb24gZmxpcFxuICAgIH0gZWxzZSBpZiAoby5mbGlwKSB7XG4gICAgICBtYXRyaXggPSBuZXcgU1ZHLk1hdHJpeCgpLm1vcnBoKG5ldyBTVkcuTWF0cml4KCkuZmxpcChcbiAgICAgICAgby5mbGlwXG4gICAgICAsIG8ub2Zmc2V0ID09IG51bGwgPyB0YXJnZXQuYmJveCgpWydjJyArIG8uZmxpcF0gOiBvLm9mZnNldFxuICAgICAgKSlcblxuICAgIC8vIGFjdCBvbiB0cmFuc2xhdGVcbiAgICB9IGVsc2UgaWYgKG8ueCAhPSBudWxsIHx8IG8ueSAhPSBudWxsKSB7XG4gICAgICBtYXRyaXggPSBuZXcgU1ZHLlRyYW5zbGF0ZShvLngsIG8ueSlcbiAgICB9XG5cbiAgICBpZighbWF0cml4KSByZXR1cm4gdGhpc1xuXG4gICAgbWF0cml4LnJlbGF0aXZlID0gcmVsYXRpdmVcblxuICAgIHRoaXMubGFzdCgpLnRyYW5zZm9ybXMucHVzaChtYXRyaXgpXG5cbiAgICBzZXRUaW1lb3V0KGZ1bmN0aW9uKCl7dGhpcy5zdGFydCgpfS5iaW5kKHRoaXMpLCAwKVxuXG4gICAgcmV0dXJuIHRoaXNcbiAgfVxufSlcblxuU1ZHLmV4dGVuZChTVkcuRWxlbWVudCwge1xuICAvLyBSZXNldCBhbGwgdHJhbnNmb3JtYXRpb25zXG4gIHVudHJhbnNmb3JtOiBmdW5jdGlvbigpIHtcbiAgICByZXR1cm4gdGhpcy5hdHRyKCd0cmFuc2Zvcm0nLCBudWxsKVxuICB9LFxuICAvLyBtZXJnZSB0aGUgd2hvbGUgdHJhbnNmb3JtYXRpb24gY2hhaW4gaW50byBvbmUgbWF0cml4IGFuZCByZXR1cm5zIGl0XG4gIG1hdHJpeGlmeTogZnVuY3Rpb24oKSB7XG5cbiAgICB2YXIgbWF0cml4ID0gKHRoaXMuYXR0cigndHJhbnNmb3JtJykgfHwgJycpXG4gICAgICAvLyBzcGxpdCB0cmFuc2Zvcm1hdGlvbnNcbiAgICAgIC5zcGxpdCgvXFwpXFxzKiw/XFxzKi8pLnNsaWNlKDAsLTEpLm1hcChmdW5jdGlvbihzdHIpe1xuICAgICAgICAvLyBnZW5lcmF0ZSBrZXkgPT4gdmFsdWUgcGFpcnNcbiAgICAgICAgdmFyIGt2ID0gc3RyLnRyaW0oKS5zcGxpdCgnKCcpXG4gICAgICAgIHJldHVybiBba3ZbMF0sIGt2WzFdLnNwbGl0KFNWRy5yZWdleC5tYXRyaXhFbGVtZW50cykubWFwKGZ1bmN0aW9uKHN0cil7IHJldHVybiBwYXJzZUZsb2F0KHN0cikgfSldXG4gICAgICB9KVxuICAgICAgLy8gY2FsY3VsYXRlIGV2ZXJ5IHRyYW5zZm9ybWF0aW9uIGludG8gb25lIG1hdHJpeFxuICAgICAgLnJlZHVjZShmdW5jdGlvbihtYXRyaXgsIHRyYW5zZm9ybSl7XG5cbiAgICAgICAgaWYodHJhbnNmb3JtWzBdID09ICdtYXRyaXgnKSByZXR1cm4gbWF0cml4Lm11bHRpcGx5KGFycmF5VG9NYXRyaXgodHJhbnNmb3JtWzFdKSlcbiAgICAgICAgcmV0dXJuIG1hdHJpeFt0cmFuc2Zvcm1bMF1dLmFwcGx5KG1hdHJpeCwgdHJhbnNmb3JtWzFdKVxuXG4gICAgICB9LCBuZXcgU1ZHLk1hdHJpeCgpKVxuXG4gICAgcmV0dXJuIG1hdHJpeFxuICB9LFxuICAvLyBhZGQgYW4gZWxlbWVudCB0byBhbm90aGVyIHBhcmVudCB3aXRob3V0IGNoYW5naW5nIHRoZSB2aXN1YWwgcmVwcmVzZW50YXRpb24gb24gdGhlIHNjcmVlblxuICB0b1BhcmVudDogZnVuY3Rpb24ocGFyZW50KSB7XG4gICAgaWYodGhpcyA9PSBwYXJlbnQpIHJldHVybiB0aGlzXG4gICAgdmFyIGN0bSA9IHRoaXMuc2NyZWVuQ1RNKClcbiAgICB2YXIgdGVtcCA9IHBhcmVudC5yZWN0KDEsMSlcbiAgICB2YXIgcEN0bSA9IHRlbXAuc2NyZWVuQ1RNKCkuaW52ZXJzZSgpXG4gICAgdGVtcC5yZW1vdmUoKVxuXG4gICAgdGhpcy5hZGRUbyhwYXJlbnQpLnVudHJhbnNmb3JtKCkudHJhbnNmb3JtKHBDdG0ubXVsdGlwbHkoY3RtKSlcblxuICAgIHJldHVybiB0aGlzXG4gIH0sXG4gIC8vIHNhbWUgYXMgYWJvdmUgd2l0aCBwYXJlbnQgZXF1YWxzIHJvb3Qtc3ZnXG4gIHRvRG9jOiBmdW5jdGlvbigpIHtcbiAgICByZXR1cm4gdGhpcy50b1BhcmVudCh0aGlzLmRvYygpKVxuICB9XG5cbn0pXG5cblNWRy5UcmFuc2Zvcm1hdGlvbiA9IFNWRy5pbnZlbnQoe1xuXG4gIGNyZWF0ZTogZnVuY3Rpb24oc291cmNlLCBpbnZlcnNlZCl7XG5cbiAgICBpZihhcmd1bWVudHMubGVuZ3RoID4gMSAmJiB0eXBlb2YgaW52ZXJzZWQgIT0gJ2Jvb2xlYW4nKXtcbiAgICAgIHJldHVybiB0aGlzLmNyZWF0ZShbXS5zbGljZS5jYWxsKGFyZ3VtZW50cykpXG4gICAgfVxuXG4gICAgaWYodHlwZW9mIHNvdXJjZSA9PSAnb2JqZWN0Jyl7XG4gICAgICBmb3IodmFyIGkgPSAwLCBsZW4gPSB0aGlzLmFyZ3VtZW50cy5sZW5ndGg7IGkgPCBsZW47ICsraSl7XG4gICAgICAgIHRoaXNbdGhpcy5hcmd1bWVudHNbaV1dID0gc291cmNlW3RoaXMuYXJndW1lbnRzW2ldXVxuICAgICAgfVxuICAgIH1cblxuICAgIGlmKEFycmF5LmlzQXJyYXkoc291cmNlKSl7XG4gICAgICBmb3IodmFyIGkgPSAwLCBsZW4gPSB0aGlzLmFyZ3VtZW50cy5sZW5ndGg7IGkgPCBsZW47ICsraSl7XG4gICAgICAgIHRoaXNbdGhpcy5hcmd1bWVudHNbaV1dID0gc291cmNlW2ldXG4gICAgICB9XG4gICAgfVxuXG4gICAgdGhpcy5pbnZlcnNlZCA9IGZhbHNlXG5cbiAgICBpZihpbnZlcnNlZCA9PT0gdHJ1ZSl7XG4gICAgICB0aGlzLmludmVyc2VkID0gdHJ1ZVxuICAgIH1cblxuICB9XG5cbiwgZXh0ZW5kOiB7XG5cbiAgICBhdDogZnVuY3Rpb24ocG9zKXtcblxuICAgICAgdmFyIHBhcmFtcyA9IFtdXG5cbiAgICAgIGZvcih2YXIgaSA9IDAsIGxlbiA9IHRoaXMuYXJndW1lbnRzLmxlbmd0aDsgaSA8IGxlbjsgKytpKXtcbiAgICAgICAgcGFyYW1zLnB1c2godGhpc1t0aGlzLmFyZ3VtZW50c1tpXV0pXG4gICAgICB9XG5cbiAgICAgIHZhciBtID0gdGhpcy5fdW5kbyB8fCBuZXcgU1ZHLk1hdHJpeCgpXG5cbiAgICAgIG0gPSBuZXcgU1ZHLk1hdHJpeCgpLm1vcnBoKFNWRy5NYXRyaXgucHJvdG90eXBlW3RoaXMubWV0aG9kXS5hcHBseShtLCBwYXJhbXMpKS5hdChwb3MpXG5cbiAgICAgIHJldHVybiB0aGlzLmludmVyc2VkID8gbS5pbnZlcnNlKCkgOiBtXG5cbiAgICB9XG5cbiAgLCB1bmRvOiBmdW5jdGlvbihvKXtcbiAgICAgIGZvcih2YXIgaSA9IDAsIGxlbiA9IHRoaXMuYXJndW1lbnRzLmxlbmd0aDsgaSA8IGxlbjsgKytpKXtcbiAgICAgICAgb1t0aGlzLmFyZ3VtZW50c1tpXV0gPSB0eXBlb2YgdGhpc1t0aGlzLmFyZ3VtZW50c1tpXV0gPT0gJ3VuZGVmaW5lZCcgPyAwIDogb1t0aGlzLmFyZ3VtZW50c1tpXV1cbiAgICAgIH1cblxuICAgICAgLy8gVGhlIG1ldGhvZCBTVkcuTWF0cml4LmV4dHJhY3Qgd2hpY2ggd2FzIHVzZWQgYmVmb3JlIGNhbGxpbmcgdGhpc1xuICAgICAgLy8gbWV0aG9kIHRvIG9idGFpbiBhIHZhbHVlIGZvciB0aGUgcGFyYW1ldGVyIG8gZG9lc24ndCByZXR1cm4gYSBjeCBhbmRcbiAgICAgIC8vIGEgY3kgc28gd2UgdXNlIHRoZSBvbmVzIHRoYXQgd2VyZSBwcm92aWRlZCB0byB0aGlzIG9iamVjdCBhdCBpdHMgY3JlYXRpb25cbiAgICAgIG8uY3ggPSB0aGlzLmN4XG4gICAgICBvLmN5ID0gdGhpcy5jeVxuXG4gICAgICB0aGlzLl91bmRvID0gbmV3IFNWR1tjYXBpdGFsaXplKHRoaXMubWV0aG9kKV0obywgdHJ1ZSkuYXQoMSlcblxuICAgICAgcmV0dXJuIHRoaXNcbiAgICB9XG5cbiAgfVxuXG59KVxuXG5TVkcuVHJhbnNsYXRlID0gU1ZHLmludmVudCh7XG5cbiAgcGFyZW50OiBTVkcuTWF0cml4XG4sIGluaGVyaXQ6IFNWRy5UcmFuc2Zvcm1hdGlvblxuXG4sIGNyZWF0ZTogZnVuY3Rpb24oc291cmNlLCBpbnZlcnNlZCl7XG4gICAgaWYodHlwZW9mIHNvdXJjZSA9PSAnb2JqZWN0JykgdGhpcy5jb25zdHJ1Y3Rvci5jYWxsKHRoaXMsIHNvdXJjZSwgaW52ZXJzZWQpXG4gICAgZWxzZSB0aGlzLmNvbnN0cnVjdG9yLmNhbGwodGhpcywgW10uc2xpY2UuY2FsbChhcmd1bWVudHMpKVxuICB9XG5cbiwgZXh0ZW5kOiB7XG4gICAgYXJndW1lbnRzOiBbJ3RyYW5zZm9ybWVkWCcsICd0cmFuc2Zvcm1lZFknXVxuICAsIG1ldGhvZDogJ3RyYW5zbGF0ZSdcbiAgfVxuXG59KVxuXG5TVkcuUm90YXRlID0gU1ZHLmludmVudCh7XG5cbiAgcGFyZW50OiBTVkcuTWF0cml4XG4sIGluaGVyaXQ6IFNWRy5UcmFuc2Zvcm1hdGlvblxuXG4sIGNyZWF0ZTogZnVuY3Rpb24oc291cmNlLCBpbnZlcnNlZCl7XG4gICAgaWYodHlwZW9mIHNvdXJjZSA9PSAnb2JqZWN0JykgdGhpcy5jb25zdHJ1Y3Rvci5jYWxsKHRoaXMsIHNvdXJjZSwgaW52ZXJzZWQpXG4gICAgZWxzZSB0aGlzLmNvbnN0cnVjdG9yLmNhbGwodGhpcywgW10uc2xpY2UuY2FsbChhcmd1bWVudHMpKVxuICB9XG5cbiwgZXh0ZW5kOiB7XG4gICAgYXJndW1lbnRzOiBbJ3JvdGF0aW9uJywgJ2N4JywgJ2N5J11cbiAgLCBtZXRob2Q6ICdyb3RhdGUnXG4gICwgYXQ6IGZ1bmN0aW9uKHBvcyl7XG4gICAgICB2YXIgbSA9IG5ldyBTVkcuTWF0cml4KCkucm90YXRlKG5ldyBTVkcuTnVtYmVyKCkubW9ycGgodGhpcy5yb3RhdGlvbiAtICh0aGlzLl91bmRvID8gdGhpcy5fdW5kby5yb3RhdGlvbiA6IDApKS5hdChwb3MpLCB0aGlzLmN4LCB0aGlzLmN5KVxuICAgICAgcmV0dXJuIHRoaXMuaW52ZXJzZWQgPyBtLmludmVyc2UoKSA6IG1cbiAgICB9XG4gICwgdW5kbzogZnVuY3Rpb24obyl7XG4gICAgICB0aGlzLl91bmRvID0gb1xuICAgIH1cbiAgfVxuXG59KVxuXG5TVkcuU2NhbGUgPSBTVkcuaW52ZW50KHtcblxuICBwYXJlbnQ6IFNWRy5NYXRyaXhcbiwgaW5oZXJpdDogU1ZHLlRyYW5zZm9ybWF0aW9uXG5cbiwgY3JlYXRlOiBmdW5jdGlvbihzb3VyY2UsIGludmVyc2VkKXtcbiAgICBpZih0eXBlb2Ygc291cmNlID09ICdvYmplY3QnKSB0aGlzLmNvbnN0cnVjdG9yLmNhbGwodGhpcywgc291cmNlLCBpbnZlcnNlZClcbiAgICBlbHNlIHRoaXMuY29uc3RydWN0b3IuY2FsbCh0aGlzLCBbXS5zbGljZS5jYWxsKGFyZ3VtZW50cykpXG4gIH1cblxuLCBleHRlbmQ6IHtcbiAgICBhcmd1bWVudHM6IFsnc2NhbGVYJywgJ3NjYWxlWScsICdjeCcsICdjeSddXG4gICwgbWV0aG9kOiAnc2NhbGUnXG4gIH1cblxufSlcblxuU1ZHLlNrZXcgPSBTVkcuaW52ZW50KHtcblxuICBwYXJlbnQ6IFNWRy5NYXRyaXhcbiwgaW5oZXJpdDogU1ZHLlRyYW5zZm9ybWF0aW9uXG5cbiwgY3JlYXRlOiBmdW5jdGlvbihzb3VyY2UsIGludmVyc2VkKXtcbiAgICBpZih0eXBlb2Ygc291cmNlID09ICdvYmplY3QnKSB0aGlzLmNvbnN0cnVjdG9yLmNhbGwodGhpcywgc291cmNlLCBpbnZlcnNlZClcbiAgICBlbHNlIHRoaXMuY29uc3RydWN0b3IuY2FsbCh0aGlzLCBbXS5zbGljZS5jYWxsKGFyZ3VtZW50cykpXG4gIH1cblxuLCBleHRlbmQ6IHtcbiAgICBhcmd1bWVudHM6IFsnc2tld1gnLCAnc2tld1knLCAnY3gnLCAnY3knXVxuICAsIG1ldGhvZDogJ3NrZXcnXG4gIH1cblxufSlcblxuU1ZHLmV4dGVuZChTVkcuRWxlbWVudCwge1xuICAvLyBEeW5hbWljIHN0eWxlIGdlbmVyYXRvclxuICBzdHlsZTogZnVuY3Rpb24ocywgdikge1xuICAgIGlmIChhcmd1bWVudHMubGVuZ3RoID09IDApIHtcbiAgICAgIC8vIGdldCBmdWxsIHN0eWxlXG4gICAgICByZXR1cm4gdGhpcy5ub2RlLnN0eWxlLmNzc1RleHQgfHwgJydcblxuICAgIH0gZWxzZSBpZiAoYXJndW1lbnRzLmxlbmd0aCA8IDIpIHtcbiAgICAgIC8vIGFwcGx5IGV2ZXJ5IHN0eWxlIGluZGl2aWR1YWxseSBpZiBhbiBvYmplY3QgaXMgcGFzc2VkXG4gICAgICBpZiAodHlwZW9mIHMgPT0gJ29iamVjdCcpIHtcbiAgICAgICAgZm9yICh2IGluIHMpIHRoaXMuc3R5bGUodiwgc1t2XSlcblxuICAgICAgfSBlbHNlIGlmIChTVkcucmVnZXguaXNDc3MudGVzdChzKSkge1xuICAgICAgICAvLyBwYXJzZSBjc3Mgc3RyaW5nXG4gICAgICAgIHMgPSBzLnNwbGl0KCc7JylcblxuICAgICAgICAvLyBhcHBseSBldmVyeSBkZWZpbml0aW9uIGluZGl2aWR1YWxseVxuICAgICAgICBmb3IgKHZhciBpID0gMDsgaSA8IHMubGVuZ3RoOyBpKyspIHtcbiAgICAgICAgICB2ID0gc1tpXS5zcGxpdCgnOicpXG4gICAgICAgICAgdGhpcy5zdHlsZSh2WzBdLnJlcGxhY2UoL1xccysvZywgJycpLCB2WzFdKVxuICAgICAgICB9XG4gICAgICB9IGVsc2Uge1xuICAgICAgICAvLyBhY3QgYXMgYSBnZXR0ZXIgaWYgdGhlIGZpcnN0IGFuZCBvbmx5IGFyZ3VtZW50IGlzIG5vdCBhbiBvYmplY3RcbiAgICAgICAgcmV0dXJuIHRoaXMubm9kZS5zdHlsZVtjYW1lbENhc2UocyldXG4gICAgICB9XG5cbiAgICB9IGVsc2Uge1xuICAgICAgdGhpcy5ub2RlLnN0eWxlW2NhbWVsQ2FzZShzKV0gPSB2ID09PSBudWxsIHx8IFNWRy5yZWdleC5pc0JsYW5rLnRlc3QodikgPyAnJyA6IHZcbiAgICB9XG5cbiAgICByZXR1cm4gdGhpc1xuICB9XG59KVxuU1ZHLlBhcmVudCA9IFNWRy5pbnZlbnQoe1xuICAvLyBJbml0aWFsaXplIG5vZGVcbiAgY3JlYXRlOiBmdW5jdGlvbihlbGVtZW50KSB7XG4gICAgdGhpcy5jb25zdHJ1Y3Rvci5jYWxsKHRoaXMsIGVsZW1lbnQpXG4gIH1cblxuICAvLyBJbmhlcml0IGZyb21cbiwgaW5oZXJpdDogU1ZHLkVsZW1lbnRcblxuICAvLyBBZGQgY2xhc3MgbWV0aG9kc1xuLCBleHRlbmQ6IHtcbiAgICAvLyBSZXR1cm5zIGFsbCBjaGlsZCBlbGVtZW50c1xuICAgIGNoaWxkcmVuOiBmdW5jdGlvbigpIHtcbiAgICAgIHJldHVybiBTVkcudXRpbHMubWFwKFNWRy51dGlscy5maWx0ZXJTVkdFbGVtZW50cyh0aGlzLm5vZGUuY2hpbGROb2RlcyksIGZ1bmN0aW9uKG5vZGUpIHtcbiAgICAgICAgcmV0dXJuIFNWRy5hZG9wdChub2RlKVxuICAgICAgfSlcbiAgICB9XG4gICAgLy8gQWRkIGdpdmVuIGVsZW1lbnQgYXQgYSBwb3NpdGlvblxuICAsIGFkZDogZnVuY3Rpb24oZWxlbWVudCwgaSkge1xuICAgICAgaWYgKGkgPT0gbnVsbClcbiAgICAgICAgdGhpcy5ub2RlLmFwcGVuZENoaWxkKGVsZW1lbnQubm9kZSlcbiAgICAgIGVsc2UgaWYgKGVsZW1lbnQubm9kZSAhPSB0aGlzLm5vZGUuY2hpbGROb2Rlc1tpXSlcbiAgICAgICAgdGhpcy5ub2RlLmluc2VydEJlZm9yZShlbGVtZW50Lm5vZGUsIHRoaXMubm9kZS5jaGlsZE5vZGVzW2ldKVxuXG4gICAgICByZXR1cm4gdGhpc1xuICAgIH1cbiAgICAvLyBCYXNpY2FsbHkgZG9lcyB0aGUgc2FtZSBhcyBgYWRkKClgIGJ1dCByZXR1cm5zIHRoZSBhZGRlZCBlbGVtZW50IGluc3RlYWRcbiAgLCBwdXQ6IGZ1bmN0aW9uKGVsZW1lbnQsIGkpIHtcbiAgICAgIHRoaXMuYWRkKGVsZW1lbnQsIGkpXG4gICAgICByZXR1cm4gZWxlbWVudFxuICAgIH1cbiAgICAvLyBDaGVja3MgaWYgdGhlIGdpdmVuIGVsZW1lbnQgaXMgYSBjaGlsZFxuICAsIGhhczogZnVuY3Rpb24oZWxlbWVudCkge1xuICAgICAgcmV0dXJuIHRoaXMuaW5kZXgoZWxlbWVudCkgPj0gMFxuICAgIH1cbiAgICAvLyBHZXRzIGluZGV4IG9mIGdpdmVuIGVsZW1lbnRcbiAgLCBpbmRleDogZnVuY3Rpb24oZWxlbWVudCkge1xuICAgICAgcmV0dXJuIFtdLnNsaWNlLmNhbGwodGhpcy5ub2RlLmNoaWxkTm9kZXMpLmluZGV4T2YoZWxlbWVudC5ub2RlKVxuICAgIH1cbiAgICAvLyBHZXQgYSBlbGVtZW50IGF0IHRoZSBnaXZlbiBpbmRleFxuICAsIGdldDogZnVuY3Rpb24oaSkge1xuICAgICAgcmV0dXJuIFNWRy5hZG9wdCh0aGlzLm5vZGUuY2hpbGROb2Rlc1tpXSlcbiAgICB9XG4gICAgLy8gR2V0IGZpcnN0IGNoaWxkXG4gICwgZmlyc3Q6IGZ1bmN0aW9uKCkge1xuICAgICAgcmV0dXJuIHRoaXMuZ2V0KDApXG4gICAgfVxuICAgIC8vIEdldCB0aGUgbGFzdCBjaGlsZFxuICAsIGxhc3Q6IGZ1bmN0aW9uKCkge1xuICAgICAgcmV0dXJuIHRoaXMuZ2V0KHRoaXMubm9kZS5jaGlsZE5vZGVzLmxlbmd0aCAtIDEpXG4gICAgfVxuICAgIC8vIEl0ZXJhdGVzIG92ZXIgYWxsIGNoaWxkcmVuIGFuZCBpbnZva2VzIGEgZ2l2ZW4gYmxvY2tcbiAgLCBlYWNoOiBmdW5jdGlvbihibG9jaywgZGVlcCkge1xuICAgICAgdmFyIGksIGlsXG4gICAgICAgICwgY2hpbGRyZW4gPSB0aGlzLmNoaWxkcmVuKClcblxuICAgICAgZm9yIChpID0gMCwgaWwgPSBjaGlsZHJlbi5sZW5ndGg7IGkgPCBpbDsgaSsrKSB7XG4gICAgICAgIGlmIChjaGlsZHJlbltpXSBpbnN0YW5jZW9mIFNWRy5FbGVtZW50KVxuICAgICAgICAgIGJsb2NrLmFwcGx5KGNoaWxkcmVuW2ldLCBbaSwgY2hpbGRyZW5dKVxuXG4gICAgICAgIGlmIChkZWVwICYmIChjaGlsZHJlbltpXSBpbnN0YW5jZW9mIFNWRy5Db250YWluZXIpKVxuICAgICAgICAgIGNoaWxkcmVuW2ldLmVhY2goYmxvY2ssIGRlZXApXG4gICAgICB9XG5cbiAgICAgIHJldHVybiB0aGlzXG4gICAgfVxuICAgIC8vIFJlbW92ZSBhIGdpdmVuIGNoaWxkXG4gICwgcmVtb3ZlRWxlbWVudDogZnVuY3Rpb24oZWxlbWVudCkge1xuICAgICAgdGhpcy5ub2RlLnJlbW92ZUNoaWxkKGVsZW1lbnQubm9kZSlcblxuICAgICAgcmV0dXJuIHRoaXNcbiAgICB9XG4gICAgLy8gUmVtb3ZlIGFsbCBlbGVtZW50cyBpbiB0aGlzIGNvbnRhaW5lclxuICAsIGNsZWFyOiBmdW5jdGlvbigpIHtcbiAgICAgIC8vIHJlbW92ZSBjaGlsZHJlblxuICAgICAgd2hpbGUodGhpcy5ub2RlLmhhc0NoaWxkTm9kZXMoKSlcbiAgICAgICAgdGhpcy5ub2RlLnJlbW92ZUNoaWxkKHRoaXMubm9kZS5sYXN0Q2hpbGQpXG5cbiAgICAgIC8vIHJlbW92ZSBkZWZzIHJlZmVyZW5jZVxuICAgICAgZGVsZXRlIHRoaXMuX2RlZnNcblxuICAgICAgcmV0dXJuIHRoaXNcbiAgICB9XG4gICwgLy8gR2V0IGRlZnNcbiAgICBkZWZzOiBmdW5jdGlvbigpIHtcbiAgICAgIHJldHVybiB0aGlzLmRvYygpLmRlZnMoKVxuICAgIH1cbiAgfVxuXG59KVxuXG5TVkcuZXh0ZW5kKFNWRy5QYXJlbnQsIHtcblxuICB1bmdyb3VwOiBmdW5jdGlvbihwYXJlbnQsIGRlcHRoKSB7XG4gICAgaWYoZGVwdGggPT09IDAgfHwgdGhpcyBpbnN0YW5jZW9mIFNWRy5EZWZzKSByZXR1cm4gdGhpc1xuXG4gICAgcGFyZW50ID0gcGFyZW50IHx8ICh0aGlzIGluc3RhbmNlb2YgU1ZHLkRvYyA/IHRoaXMgOiB0aGlzLnBhcmVudChTVkcuUGFyZW50KSlcbiAgICBkZXB0aCA9IGRlcHRoIHx8IEluZmluaXR5XG5cbiAgICB0aGlzLmVhY2goZnVuY3Rpb24oKXtcbiAgICAgIGlmKHRoaXMgaW5zdGFuY2VvZiBTVkcuRGVmcykgcmV0dXJuIHRoaXNcbiAgICAgIGlmKHRoaXMgaW5zdGFuY2VvZiBTVkcuUGFyZW50KSByZXR1cm4gdGhpcy51bmdyb3VwKHBhcmVudCwgZGVwdGgtMSlcbiAgICAgIHJldHVybiB0aGlzLnRvUGFyZW50KHBhcmVudClcbiAgICB9KVxuXG4gICAgdGhpcy5ub2RlLmZpcnN0Q2hpbGQgfHwgdGhpcy5yZW1vdmUoKVxuXG4gICAgcmV0dXJuIHRoaXNcbiAgfSxcblxuICBmbGF0dGVuOiBmdW5jdGlvbihwYXJlbnQsIGRlcHRoKSB7XG4gICAgcmV0dXJuIHRoaXMudW5ncm91cChwYXJlbnQsIGRlcHRoKVxuICB9XG5cbn0pXG5TVkcuQ29udGFpbmVyID0gU1ZHLmludmVudCh7XG4gIC8vIEluaXRpYWxpemUgbm9kZVxuICBjcmVhdGU6IGZ1bmN0aW9uKGVsZW1lbnQpIHtcbiAgICB0aGlzLmNvbnN0cnVjdG9yLmNhbGwodGhpcywgZWxlbWVudClcbiAgfVxuXG4gIC8vIEluaGVyaXQgZnJvbVxuLCBpbmhlcml0OiBTVkcuUGFyZW50XG5cbn0pXG5cblNWRy5WaWV3Qm94ID0gU1ZHLmludmVudCh7XG5cbiAgY3JlYXRlOiBmdW5jdGlvbihzb3VyY2UpIHtcbiAgICB2YXIgaSwgYmFzZSA9IFswLCAwLCAwLCAwXVxuXG4gICAgdmFyIHgsIHksIHdpZHRoLCBoZWlnaHQsIGJveCwgdmlldywgd2UsIGhlXG4gICAgICAsIHdtICAgPSAxIC8vIHdpZHRoIG11bHRpcGxpZXJcbiAgICAgICwgaG0gICA9IDEgLy8gaGVpZ2h0IG11bHRpcGxpZXJcbiAgICAgICwgcmVnICA9IC9bKy1dPyg/OlxcZCsoPzpcXC5cXGQqKT98XFwuXFxkKykoPzplWystXT9cXGQrKT8vZ2lcblxuICAgIGlmKHNvdXJjZSBpbnN0YW5jZW9mIFNWRy5FbGVtZW50KXtcblxuICAgICAgd2UgPSBzb3VyY2VcbiAgICAgIGhlID0gc291cmNlXG4gICAgICB2aWV3ID0gKHNvdXJjZS5hdHRyKCd2aWV3Qm94JykgfHwgJycpLm1hdGNoKHJlZylcbiAgICAgIGJveCA9IHNvdXJjZS5iYm94XG5cbiAgICAgIC8vIGdldCBkaW1lbnNpb25zIG9mIGN1cnJlbnQgbm9kZVxuICAgICAgd2lkdGggID0gbmV3IFNWRy5OdW1iZXIoc291cmNlLndpZHRoKCkpXG4gICAgICBoZWlnaHQgPSBuZXcgU1ZHLk51bWJlcihzb3VyY2UuaGVpZ2h0KCkpXG5cbiAgICAgIC8vIGZpbmQgbmVhcmVzdCBub24tcGVyY2VudHVhbCBkaW1lbnNpb25zXG4gICAgICB3aGlsZSAod2lkdGgudW5pdCA9PSAnJScpIHtcbiAgICAgICAgd20gKj0gd2lkdGgudmFsdWVcbiAgICAgICAgd2lkdGggPSBuZXcgU1ZHLk51bWJlcih3ZSBpbnN0YW5jZW9mIFNWRy5Eb2MgPyB3ZS5wYXJlbnQoKS5vZmZzZXRXaWR0aCA6IHdlLnBhcmVudCgpLndpZHRoKCkpXG4gICAgICAgIHdlID0gd2UucGFyZW50KClcbiAgICAgIH1cbiAgICAgIHdoaWxlIChoZWlnaHQudW5pdCA9PSAnJScpIHtcbiAgICAgICAgaG0gKj0gaGVpZ2h0LnZhbHVlXG4gICAgICAgIGhlaWdodCA9IG5ldyBTVkcuTnVtYmVyKGhlIGluc3RhbmNlb2YgU1ZHLkRvYyA/IGhlLnBhcmVudCgpLm9mZnNldEhlaWdodCA6IGhlLnBhcmVudCgpLmhlaWdodCgpKVxuICAgICAgICBoZSA9IGhlLnBhcmVudCgpXG4gICAgICB9XG5cbiAgICAgIC8vIGVuc3VyZSBkZWZhdWx0c1xuICAgICAgdGhpcy54ICAgICAgPSAwXG4gICAgICB0aGlzLnkgICAgICA9IDBcbiAgICAgIHRoaXMud2lkdGggID0gd2lkdGggICogd21cbiAgICAgIHRoaXMuaGVpZ2h0ID0gaGVpZ2h0ICogaG1cbiAgICAgIHRoaXMuem9vbSAgID0gMVxuXG4gICAgICBpZiAodmlldykge1xuICAgICAgICAvLyBnZXQgd2lkdGggYW5kIGhlaWdodCBmcm9tIHZpZXdib3hcbiAgICAgICAgeCAgICAgID0gcGFyc2VGbG9hdCh2aWV3WzBdKVxuICAgICAgICB5ICAgICAgPSBwYXJzZUZsb2F0KHZpZXdbMV0pXG4gICAgICAgIHdpZHRoICA9IHBhcnNlRmxvYXQodmlld1syXSlcbiAgICAgICAgaGVpZ2h0ID0gcGFyc2VGbG9hdCh2aWV3WzNdKVxuXG4gICAgICAgIC8vIGNhbGN1bGF0ZSB6b29tIGFjY29yaW5nIHRvIHZpZXdib3hcbiAgICAgICAgdGhpcy56b29tID0gKCh0aGlzLndpZHRoIC8gdGhpcy5oZWlnaHQpID4gKHdpZHRoIC8gaGVpZ2h0KSkgP1xuICAgICAgICAgIHRoaXMuaGVpZ2h0IC8gaGVpZ2h0IDpcbiAgICAgICAgICB0aGlzLndpZHRoICAvIHdpZHRoXG5cbiAgICAgICAgLy8gY2FsY3VsYXRlIHJlYWwgcGl4ZWwgZGltZW5zaW9ucyBvbiBwYXJlbnQgU1ZHLkRvYyBlbGVtZW50XG4gICAgICAgIHRoaXMueCAgICAgID0geFxuICAgICAgICB0aGlzLnkgICAgICA9IHlcbiAgICAgICAgdGhpcy53aWR0aCAgPSB3aWR0aFxuICAgICAgICB0aGlzLmhlaWdodCA9IGhlaWdodFxuXG4gICAgICB9XG5cbiAgICB9ZWxzZXtcblxuICAgICAgLy8gZW5zdXJlIHNvdXJjZSBhcyBvYmplY3RcbiAgICAgIHNvdXJjZSA9IHR5cGVvZiBzb3VyY2UgPT09ICdzdHJpbmcnID9cbiAgICAgICAgc291cmNlLm1hdGNoKHJlZykubWFwKGZ1bmN0aW9uKGVsKXsgcmV0dXJuIHBhcnNlRmxvYXQoZWwpIH0pIDpcbiAgICAgIEFycmF5LmlzQXJyYXkoc291cmNlKSA/XG4gICAgICAgIHNvdXJjZSA6XG4gICAgICB0eXBlb2Ygc291cmNlID09ICdvYmplY3QnID9cbiAgICAgICAgW3NvdXJjZS54LCBzb3VyY2UueSwgc291cmNlLndpZHRoLCBzb3VyY2UuaGVpZ2h0XSA6XG4gICAgICBhcmd1bWVudHMubGVuZ3RoID09IDQgP1xuICAgICAgICBbXS5zbGljZS5jYWxsKGFyZ3VtZW50cykgOlxuICAgICAgICBiYXNlXG5cbiAgICAgIHRoaXMueCA9IHNvdXJjZVswXVxuICAgICAgdGhpcy55ID0gc291cmNlWzFdXG4gICAgICB0aGlzLndpZHRoID0gc291cmNlWzJdXG4gICAgICB0aGlzLmhlaWdodCA9IHNvdXJjZVszXVxuICAgIH1cblxuXG4gIH1cblxuLCBleHRlbmQ6IHtcblxuICAgIHRvU3RyaW5nOiBmdW5jdGlvbigpIHtcbiAgICAgIHJldHVybiB0aGlzLnggKyAnICcgKyB0aGlzLnkgKyAnICcgKyB0aGlzLndpZHRoICsgJyAnICsgdGhpcy5oZWlnaHRcbiAgICB9XG4gICwgbW9ycGg6IGZ1bmN0aW9uKHYpe1xuXG4gICAgICB2YXIgdiA9IGFyZ3VtZW50cy5sZW5ndGggPT0gMSA/XG4gICAgICAgIFt2LngsIHYueSwgdi53aWR0aCwgdi5oZWlnaHRdIDpcbiAgICAgICAgW10uc2xpY2UuY2FsbChhcmd1bWVudHMpXG5cbiAgICAgIHRoaXMuZGVzdGluYXRpb24gPSBuZXcgU1ZHLlZpZXdCb3godilcblxuICAgICAgcmV0dXJuIHRoaXNcblxuICAgIH1cblxuICAsIGF0OiBmdW5jdGlvbihwb3MpIHtcblxuICAgIGlmKCF0aGlzLmRlc3RpbmF0aW9uKSByZXR1cm4gdGhpc1xuXG4gICAgcmV0dXJuIG5ldyBTVkcuVmlld0JveChbXG4gICAgICAgIHRoaXMueCArICh0aGlzLmRlc3RpbmF0aW9uLnggLSB0aGlzLngpICogcG9zXG4gICAgICAsIHRoaXMueSArICh0aGlzLmRlc3RpbmF0aW9uLnkgLSB0aGlzLnkpICogcG9zXG4gICAgICAsIHRoaXMud2lkdGggKyAodGhpcy5kZXN0aW5hdGlvbi53aWR0aCAtIHRoaXMud2lkdGgpICogcG9zXG4gICAgICAsIHRoaXMuaGVpZ2h0ICsgKHRoaXMuZGVzdGluYXRpb24uaGVpZ2h0IC0gdGhpcy5oZWlnaHQpICogcG9zXG4gICAgXSlcblxuICAgIH1cblxuICB9XG5cbiAgLy8gRGVmaW5lIHBhcmVudFxuLCBwYXJlbnQ6IFNWRy5Db250YWluZXJcblxuICAvLyBBZGQgcGFyZW50IG1ldGhvZFxuLCBjb25zdHJ1Y3Q6IHtcblxuICAgIC8vIGdldC9zZXQgdmlld2JveFxuICAgIHZpZXdib3g6IGZ1bmN0aW9uKHYpIHtcbiAgICAgIGlmIChhcmd1bWVudHMubGVuZ3RoID09IDApXG4gICAgICAgIC8vIGFjdCBhcyBhIGdldHRlciBpZiB0aGVyZSBhcmUgbm8gYXJndW1lbnRzXG4gICAgICAgIHJldHVybiBuZXcgU1ZHLlZpZXdCb3godGhpcylcblxuICAgICAgLy8gb3RoZXJ3aXNlIGFjdCBhcyBhIHNldHRlclxuICAgICAgdiA9IGFyZ3VtZW50cy5sZW5ndGggPT0gMSA/XG4gICAgICAgIFt2LngsIHYueSwgdi53aWR0aCwgdi5oZWlnaHRdIDpcbiAgICAgICAgW10uc2xpY2UuY2FsbChhcmd1bWVudHMpXG5cbiAgICAgIHJldHVybiB0aGlzLmF0dHIoJ3ZpZXdCb3gnLCB2KVxuICAgIH1cblxuICB9XG5cbn0pXG4vLyBBZGQgZXZlbnRzIHRvIGVsZW1lbnRzXG47WyAgJ2NsaWNrJ1xuICAsICdkYmxjbGljaydcbiAgLCAnbW91c2Vkb3duJ1xuICAsICdtb3VzZXVwJ1xuICAsICdtb3VzZW92ZXInXG4gICwgJ21vdXNlb3V0J1xuICAsICdtb3VzZW1vdmUnXG4gIC8vICwgJ21vdXNlZW50ZXInIC0+IG5vdCBzdXBwb3J0ZWQgYnkgSUVcbiAgLy8gLCAnbW91c2VsZWF2ZScgLT4gbm90IHN1cHBvcnRlZCBieSBJRVxuICAsICd0b3VjaHN0YXJ0J1xuICAsICd0b3VjaG1vdmUnXG4gICwgJ3RvdWNobGVhdmUnXG4gICwgJ3RvdWNoZW5kJ1xuICAsICd0b3VjaGNhbmNlbCcgXS5mb3JFYWNoKGZ1bmN0aW9uKGV2ZW50KSB7XG5cbiAgLy8gYWRkIGV2ZW50IHRvIFNWRy5FbGVtZW50XG4gIFNWRy5FbGVtZW50LnByb3RvdHlwZVtldmVudF0gPSBmdW5jdGlvbihmKSB7XG4gICAgdmFyIHNlbGYgPSB0aGlzXG5cbiAgICAvLyBiaW5kIGV2ZW50IHRvIGVsZW1lbnQgcmF0aGVyIHRoYW4gZWxlbWVudCBub2RlXG4gICAgdGhpcy5ub2RlWydvbicgKyBldmVudF0gPSB0eXBlb2YgZiA9PSAnZnVuY3Rpb24nID9cbiAgICAgIGZ1bmN0aW9uKCkgeyByZXR1cm4gZi5hcHBseShzZWxmLCBhcmd1bWVudHMpIH0gOiBudWxsXG5cbiAgICByZXR1cm4gdGhpc1xuICB9XG5cbn0pXG5cbi8vIEluaXRpYWxpemUgbGlzdGVuZXJzIHN0YWNrXG5TVkcubGlzdGVuZXJzID0gW11cblNWRy5oYW5kbGVyTWFwID0gW11cblNWRy5saXN0ZW5lcklkID0gMFxuXG4vLyBBZGQgZXZlbnQgYmluZGVyIGluIHRoZSBTVkcgbmFtZXNwYWNlXG5TVkcub24gPSBmdW5jdGlvbihub2RlLCBldmVudCwgbGlzdGVuZXIsIGJpbmRpbmcpIHtcbiAgLy8gY3JlYXRlIGxpc3RlbmVyLCBnZXQgb2JqZWN0LWluZGV4XG4gIHZhciBsICAgICA9IGxpc3RlbmVyLmJpbmQoYmluZGluZyB8fCBub2RlLmluc3RhbmNlIHx8IG5vZGUpXG4gICAgLCBpbmRleCA9IChTVkcuaGFuZGxlck1hcC5pbmRleE9mKG5vZGUpICsgMSB8fCBTVkcuaGFuZGxlck1hcC5wdXNoKG5vZGUpKSAtIDFcbiAgICAsIGV2ICAgID0gZXZlbnQuc3BsaXQoJy4nKVswXVxuICAgICwgbnMgICAgPSBldmVudC5zcGxpdCgnLicpWzFdIHx8ICcqJ1xuXG5cbiAgLy8gZW5zdXJlIHZhbGlkIG9iamVjdFxuICBTVkcubGlzdGVuZXJzW2luZGV4XSAgICAgICAgID0gU1ZHLmxpc3RlbmVyc1tpbmRleF0gICAgICAgICB8fCB7fVxuICBTVkcubGlzdGVuZXJzW2luZGV4XVtldl0gICAgID0gU1ZHLmxpc3RlbmVyc1tpbmRleF1bZXZdICAgICB8fCB7fVxuICBTVkcubGlzdGVuZXJzW2luZGV4XVtldl1bbnNdID0gU1ZHLmxpc3RlbmVyc1tpbmRleF1bZXZdW25zXSB8fCB7fVxuXG4gIGlmKCFsaXN0ZW5lci5fc3ZnanNMaXN0ZW5lcklkKVxuICAgIGxpc3RlbmVyLl9zdmdqc0xpc3RlbmVySWQgPSArK1NWRy5saXN0ZW5lcklkXG5cbiAgLy8gcmVmZXJlbmNlIGxpc3RlbmVyXG4gIFNWRy5saXN0ZW5lcnNbaW5kZXhdW2V2XVtuc11bbGlzdGVuZXIuX3N2Z2pzTGlzdGVuZXJJZF0gPSBsXG5cbiAgLy8gYWRkIGxpc3RlbmVyXG4gIG5vZGUuYWRkRXZlbnRMaXN0ZW5lcihldiwgbCwgZmFsc2UpXG59XG5cbi8vIEFkZCBldmVudCB1bmJpbmRlciBpbiB0aGUgU1ZHIG5hbWVzcGFjZVxuU1ZHLm9mZiA9IGZ1bmN0aW9uKG5vZGUsIGV2ZW50LCBsaXN0ZW5lcikge1xuICB2YXIgaW5kZXggPSBTVkcuaGFuZGxlck1hcC5pbmRleE9mKG5vZGUpXG4gICAgLCBldiAgICA9IGV2ZW50ICYmIGV2ZW50LnNwbGl0KCcuJylbMF1cbiAgICAsIG5zICAgID0gZXZlbnQgJiYgZXZlbnQuc3BsaXQoJy4nKVsxXVxuXG4gIGlmKGluZGV4ID09IC0xKSByZXR1cm5cblxuICBpZiAobGlzdGVuZXIpIHtcbiAgICBpZih0eXBlb2YgbGlzdGVuZXIgPT0gJ2Z1bmN0aW9uJykgbGlzdGVuZXIgPSBsaXN0ZW5lci5fc3ZnanNMaXN0ZW5lcklkXG4gICAgaWYoIWxpc3RlbmVyKSByZXR1cm5cblxuICAgIC8vIHJlbW92ZSBsaXN0ZW5lciByZWZlcmVuY2VcbiAgICBpZiAoU1ZHLmxpc3RlbmVyc1tpbmRleF1bZXZdICYmIFNWRy5saXN0ZW5lcnNbaW5kZXhdW2V2XVtucyB8fCAnKiddKSB7XG4gICAgICAvLyByZW1vdmUgbGlzdGVuZXJcbiAgICAgIG5vZGUucmVtb3ZlRXZlbnRMaXN0ZW5lcihldiwgU1ZHLmxpc3RlbmVyc1tpbmRleF1bZXZdW25zIHx8ICcqJ11bbGlzdGVuZXJdLCBmYWxzZSlcblxuICAgICAgZGVsZXRlIFNWRy5saXN0ZW5lcnNbaW5kZXhdW2V2XVtucyB8fCAnKiddW2xpc3RlbmVyXVxuICAgIH1cblxuICB9IGVsc2UgaWYgKG5zICYmIGV2KSB7XG4gICAgLy8gcmVtb3ZlIGFsbCBsaXN0ZW5lcnMgZm9yIGEgbmFtZXNwYWNlZCBldmVudFxuICAgIGlmIChTVkcubGlzdGVuZXJzW2luZGV4XVtldl0gJiYgU1ZHLmxpc3RlbmVyc1tpbmRleF1bZXZdW25zXSkge1xuICAgICAgZm9yIChsaXN0ZW5lciBpbiBTVkcubGlzdGVuZXJzW2luZGV4XVtldl1bbnNdKVxuICAgICAgICBTVkcub2ZmKG5vZGUsIFtldiwgbnNdLmpvaW4oJy4nKSwgbGlzdGVuZXIpXG5cbiAgICAgIGRlbGV0ZSBTVkcubGlzdGVuZXJzW2luZGV4XVtldl1bbnNdXG4gICAgfVxuXG4gIH0gZWxzZSBpZiAobnMpe1xuICAgIC8vIHJlbW92ZSBhbGwgbGlzdGVuZXJzIGZvciBhIHNwZWNpZmljIG5hbWVzcGFjZVxuICAgIGZvcihldmVudCBpbiBTVkcubGlzdGVuZXJzW2luZGV4XSl7XG4gICAgICAgIGZvcihuYW1lc3BhY2UgaW4gU1ZHLmxpc3RlbmVyc1tpbmRleF1bZXZlbnRdKXtcbiAgICAgICAgICAgIGlmKG5zID09PSBuYW1lc3BhY2Upe1xuICAgICAgICAgICAgICAgIFNWRy5vZmYobm9kZSwgW2V2ZW50LCBuc10uam9pbignLicpKVxuICAgICAgICAgICAgfVxuICAgICAgICB9XG4gICAgfVxuXG4gIH0gZWxzZSBpZiAoZXYpIHtcbiAgICAvLyByZW1vdmUgYWxsIGxpc3RlbmVycyBmb3IgdGhlIGV2ZW50XG4gICAgaWYgKFNWRy5saXN0ZW5lcnNbaW5kZXhdW2V2XSkge1xuICAgICAgZm9yIChuYW1lc3BhY2UgaW4gU1ZHLmxpc3RlbmVyc1tpbmRleF1bZXZdKVxuICAgICAgICBTVkcub2ZmKG5vZGUsIFtldiwgbmFtZXNwYWNlXS5qb2luKCcuJykpXG5cbiAgICAgIGRlbGV0ZSBTVkcubGlzdGVuZXJzW2luZGV4XVtldl1cbiAgICB9XG5cbiAgfSBlbHNlIHtcbiAgICAvLyByZW1vdmUgYWxsIGxpc3RlbmVycyBvbiBhIGdpdmVuIG5vZGVcbiAgICBmb3IgKGV2ZW50IGluIFNWRy5saXN0ZW5lcnNbaW5kZXhdKVxuICAgICAgU1ZHLm9mZihub2RlLCBldmVudClcblxuICAgIGRlbGV0ZSBTVkcubGlzdGVuZXJzW2luZGV4XVxuXG4gIH1cbn1cblxuLy9cblNWRy5leHRlbmQoU1ZHLkVsZW1lbnQsIHtcbiAgLy8gQmluZCBnaXZlbiBldmVudCB0byBsaXN0ZW5lclxuICBvbjogZnVuY3Rpb24oZXZlbnQsIGxpc3RlbmVyLCBiaW5kaW5nKSB7XG4gICAgU1ZHLm9uKHRoaXMubm9kZSwgZXZlbnQsIGxpc3RlbmVyLCBiaW5kaW5nKVxuXG4gICAgcmV0dXJuIHRoaXNcbiAgfVxuICAvLyBVbmJpbmQgZXZlbnQgZnJvbSBsaXN0ZW5lclxuLCBvZmY6IGZ1bmN0aW9uKGV2ZW50LCBsaXN0ZW5lcikge1xuICAgIFNWRy5vZmYodGhpcy5ub2RlLCBldmVudCwgbGlzdGVuZXIpXG5cbiAgICByZXR1cm4gdGhpc1xuICB9XG4gIC8vIEZpcmUgZ2l2ZW4gZXZlbnRcbiwgZmlyZTogZnVuY3Rpb24oZXZlbnQsIGRhdGEpIHtcblxuICAgIC8vIERpc3BhdGNoIGV2ZW50XG4gICAgaWYoZXZlbnQgaW5zdGFuY2VvZiBFdmVudCl7XG4gICAgICAgIHRoaXMubm9kZS5kaXNwYXRjaEV2ZW50KGV2ZW50KVxuICAgIH1lbHNle1xuICAgICAgICB0aGlzLm5vZGUuZGlzcGF0Y2hFdmVudChuZXcgQ3VzdG9tRXZlbnQoZXZlbnQsIHtkZXRhaWw6ZGF0YX0pKVxuICAgIH1cblxuICAgIHJldHVybiB0aGlzXG4gIH1cbn0pXG5cblNWRy5EZWZzID0gU1ZHLmludmVudCh7XG4gIC8vIEluaXRpYWxpemUgbm9kZVxuICBjcmVhdGU6ICdkZWZzJ1xuXG4gIC8vIEluaGVyaXQgZnJvbVxuLCBpbmhlcml0OiBTVkcuQ29udGFpbmVyXG5cbn0pXG5TVkcuRyA9IFNWRy5pbnZlbnQoe1xuICAvLyBJbml0aWFsaXplIG5vZGVcbiAgY3JlYXRlOiAnZydcblxuICAvLyBJbmhlcml0IGZyb21cbiwgaW5oZXJpdDogU1ZHLkNvbnRhaW5lclxuXG4gIC8vIEFkZCBjbGFzcyBtZXRob2RzXG4sIGV4dGVuZDoge1xuICAgIC8vIE1vdmUgb3ZlciB4LWF4aXNcbiAgICB4OiBmdW5jdGlvbih4KSB7XG4gICAgICByZXR1cm4geCA9PSBudWxsID8gdGhpcy50cmFuc2Zvcm0oJ3gnKSA6IHRoaXMudHJhbnNmb3JtKHsgeDogeCAtIHRoaXMueCgpIH0sIHRydWUpXG4gICAgfVxuICAgIC8vIE1vdmUgb3ZlciB5LWF4aXNcbiAgLCB5OiBmdW5jdGlvbih5KSB7XG4gICAgICByZXR1cm4geSA9PSBudWxsID8gdGhpcy50cmFuc2Zvcm0oJ3knKSA6IHRoaXMudHJhbnNmb3JtKHsgeTogeSAtIHRoaXMueSgpIH0sIHRydWUpXG4gICAgfVxuICAgIC8vIE1vdmUgYnkgY2VudGVyIG92ZXIgeC1heGlzXG4gICwgY3g6IGZ1bmN0aW9uKHgpIHtcbiAgICAgIHJldHVybiB4ID09IG51bGwgPyB0aGlzLmdib3goKS5jeCA6IHRoaXMueCh4IC0gdGhpcy5nYm94KCkud2lkdGggLyAyKVxuICAgIH1cbiAgICAvLyBNb3ZlIGJ5IGNlbnRlciBvdmVyIHktYXhpc1xuICAsIGN5OiBmdW5jdGlvbih5KSB7XG4gICAgICByZXR1cm4geSA9PSBudWxsID8gdGhpcy5nYm94KCkuY3kgOiB0aGlzLnkoeSAtIHRoaXMuZ2JveCgpLmhlaWdodCAvIDIpXG4gICAgfVxuICAsIGdib3g6IGZ1bmN0aW9uKCkge1xuXG4gICAgICB2YXIgYmJveCAgPSB0aGlzLmJib3goKVxuICAgICAgICAsIHRyYW5zID0gdGhpcy50cmFuc2Zvcm0oKVxuXG4gICAgICBiYm94LnggICs9IHRyYW5zLnhcbiAgICAgIGJib3gueDIgKz0gdHJhbnMueFxuICAgICAgYmJveC5jeCArPSB0cmFucy54XG5cbiAgICAgIGJib3gueSAgKz0gdHJhbnMueVxuICAgICAgYmJveC55MiArPSB0cmFucy55XG4gICAgICBiYm94LmN5ICs9IHRyYW5zLnlcblxuICAgICAgcmV0dXJuIGJib3hcbiAgICB9XG4gIH1cblxuICAvLyBBZGQgcGFyZW50IG1ldGhvZFxuLCBjb25zdHJ1Y3Q6IHtcbiAgICAvLyBDcmVhdGUgYSBncm91cCBlbGVtZW50XG4gICAgZ3JvdXA6IGZ1bmN0aW9uKCkge1xuICAgICAgcmV0dXJuIHRoaXMucHV0KG5ldyBTVkcuRylcbiAgICB9XG4gIH1cbn0pXG5cbi8vICMjIyBUaGlzIG1vZHVsZSBhZGRzIGJhY2t3YXJkIC8gZm9yd2FyZCBmdW5jdGlvbmFsaXR5IHRvIGVsZW1lbnRzLlxuXG4vL1xuU1ZHLmV4dGVuZChTVkcuRWxlbWVudCwge1xuICAvLyBHZXQgYWxsIHNpYmxpbmdzLCBpbmNsdWRpbmcgbXlzZWxmXG4gIHNpYmxpbmdzOiBmdW5jdGlvbigpIHtcbiAgICByZXR1cm4gdGhpcy5wYXJlbnQoKS5jaGlsZHJlbigpXG4gIH1cbiAgLy8gR2V0IHRoZSBjdXJlbnQgcG9zaXRpb24gc2libGluZ3NcbiwgcG9zaXRpb246IGZ1bmN0aW9uKCkge1xuICAgIHJldHVybiB0aGlzLnBhcmVudCgpLmluZGV4KHRoaXMpXG4gIH1cbiAgLy8gR2V0IHRoZSBuZXh0IGVsZW1lbnQgKHdpbGwgcmV0dXJuIG51bGwgaWYgdGhlcmUgaXMgbm9uZSlcbiwgbmV4dDogZnVuY3Rpb24oKSB7XG4gICAgcmV0dXJuIHRoaXMuc2libGluZ3MoKVt0aGlzLnBvc2l0aW9uKCkgKyAxXVxuICB9XG4gIC8vIEdldCB0aGUgbmV4dCBlbGVtZW50ICh3aWxsIHJldHVybiBudWxsIGlmIHRoZXJlIGlzIG5vbmUpXG4sIHByZXZpb3VzOiBmdW5jdGlvbigpIHtcbiAgICByZXR1cm4gdGhpcy5zaWJsaW5ncygpW3RoaXMucG9zaXRpb24oKSAtIDFdXG4gIH1cbiAgLy8gU2VuZCBnaXZlbiBlbGVtZW50IG9uZSBzdGVwIGZvcndhcmRcbiwgZm9yd2FyZDogZnVuY3Rpb24oKSB7XG4gICAgdmFyIGkgPSB0aGlzLnBvc2l0aW9uKCkgKyAxXG4gICAgICAsIHAgPSB0aGlzLnBhcmVudCgpXG5cbiAgICAvLyBtb3ZlIG5vZGUgb25lIHN0ZXAgZm9yd2FyZFxuICAgIHAucmVtb3ZlRWxlbWVudCh0aGlzKS5hZGQodGhpcywgaSlcblxuICAgIC8vIG1ha2Ugc3VyZSBkZWZzIG5vZGUgaXMgYWx3YXlzIGF0IHRoZSB0b3BcbiAgICBpZiAocCBpbnN0YW5jZW9mIFNWRy5Eb2MpXG4gICAgICBwLm5vZGUuYXBwZW5kQ2hpbGQocC5kZWZzKCkubm9kZSlcblxuICAgIHJldHVybiB0aGlzXG4gIH1cbiAgLy8gU2VuZCBnaXZlbiBlbGVtZW50IG9uZSBzdGVwIGJhY2t3YXJkXG4sIGJhY2t3YXJkOiBmdW5jdGlvbigpIHtcbiAgICB2YXIgaSA9IHRoaXMucG9zaXRpb24oKVxuXG4gICAgaWYgKGkgPiAwKVxuICAgICAgdGhpcy5wYXJlbnQoKS5yZW1vdmVFbGVtZW50KHRoaXMpLmFkZCh0aGlzLCBpIC0gMSlcblxuICAgIHJldHVybiB0aGlzXG4gIH1cbiAgLy8gU2VuZCBnaXZlbiBlbGVtZW50IGFsbCB0aGUgd2F5IHRvIHRoZSBmcm9udFxuLCBmcm9udDogZnVuY3Rpb24oKSB7XG4gICAgdmFyIHAgPSB0aGlzLnBhcmVudCgpXG5cbiAgICAvLyBNb3ZlIG5vZGUgZm9yd2FyZFxuICAgIHAubm9kZS5hcHBlbmRDaGlsZCh0aGlzLm5vZGUpXG5cbiAgICAvLyBNYWtlIHN1cmUgZGVmcyBub2RlIGlzIGFsd2F5cyBhdCB0aGUgdG9wXG4gICAgaWYgKHAgaW5zdGFuY2VvZiBTVkcuRG9jKVxuICAgICAgcC5ub2RlLmFwcGVuZENoaWxkKHAuZGVmcygpLm5vZGUpXG5cbiAgICByZXR1cm4gdGhpc1xuICB9XG4gIC8vIFNlbmQgZ2l2ZW4gZWxlbWVudCBhbGwgdGhlIHdheSB0byB0aGUgYmFja1xuLCBiYWNrOiBmdW5jdGlvbigpIHtcbiAgICBpZiAodGhpcy5wb3NpdGlvbigpID4gMClcbiAgICAgIHRoaXMucGFyZW50KCkucmVtb3ZlRWxlbWVudCh0aGlzKS5hZGQodGhpcywgMClcblxuICAgIHJldHVybiB0aGlzXG4gIH1cbiAgLy8gSW5zZXJ0cyBhIGdpdmVuIGVsZW1lbnQgYmVmb3JlIHRoZSB0YXJnZXRlZCBlbGVtZW50XG4sIGJlZm9yZTogZnVuY3Rpb24oZWxlbWVudCkge1xuICAgIGVsZW1lbnQucmVtb3ZlKClcblxuICAgIHZhciBpID0gdGhpcy5wb3NpdGlvbigpXG5cbiAgICB0aGlzLnBhcmVudCgpLmFkZChlbGVtZW50LCBpKVxuXG4gICAgcmV0dXJuIHRoaXNcbiAgfVxuICAvLyBJbnN0ZXJzIGEgZ2l2ZW4gZWxlbWVudCBhZnRlciB0aGUgdGFyZ2V0ZWQgZWxlbWVudFxuLCBhZnRlcjogZnVuY3Rpb24oZWxlbWVudCkge1xuICAgIGVsZW1lbnQucmVtb3ZlKClcblxuICAgIHZhciBpID0gdGhpcy5wb3NpdGlvbigpXG5cbiAgICB0aGlzLnBhcmVudCgpLmFkZChlbGVtZW50LCBpICsgMSlcblxuICAgIHJldHVybiB0aGlzXG4gIH1cblxufSlcblNWRy5NYXNrID0gU1ZHLmludmVudCh7XG4gIC8vIEluaXRpYWxpemUgbm9kZVxuICBjcmVhdGU6IGZ1bmN0aW9uKCkge1xuICAgIHRoaXMuY29uc3RydWN0b3IuY2FsbCh0aGlzLCBTVkcuY3JlYXRlKCdtYXNrJykpXG5cbiAgICAvLyBrZWVwIHJlZmVyZW5jZXMgdG8gbWFza2VkIGVsZW1lbnRzXG4gICAgdGhpcy50YXJnZXRzID0gW11cbiAgfVxuXG4gIC8vIEluaGVyaXQgZnJvbVxuLCBpbmhlcml0OiBTVkcuQ29udGFpbmVyXG5cbiAgLy8gQWRkIGNsYXNzIG1ldGhvZHNcbiwgZXh0ZW5kOiB7XG4gICAgLy8gVW5tYXNrIGFsbCBtYXNrZWQgZWxlbWVudHMgYW5kIHJlbW92ZSBpdHNlbGZcbiAgICByZW1vdmU6IGZ1bmN0aW9uKCkge1xuICAgICAgLy8gdW5tYXNrIGFsbCB0YXJnZXRzXG4gICAgICBmb3IgKHZhciBpID0gdGhpcy50YXJnZXRzLmxlbmd0aCAtIDE7IGkgPj0gMDsgaS0tKVxuICAgICAgICBpZiAodGhpcy50YXJnZXRzW2ldKVxuICAgICAgICAgIHRoaXMudGFyZ2V0c1tpXS51bm1hc2soKVxuICAgICAgdGhpcy50YXJnZXRzID0gW11cblxuICAgICAgLy8gcmVtb3ZlIG1hc2sgZnJvbSBwYXJlbnRcbiAgICAgIHRoaXMucGFyZW50KCkucmVtb3ZlRWxlbWVudCh0aGlzKVxuXG4gICAgICByZXR1cm4gdGhpc1xuICAgIH1cbiAgfVxuXG4gIC8vIEFkZCBwYXJlbnQgbWV0aG9kXG4sIGNvbnN0cnVjdDoge1xuICAgIC8vIENyZWF0ZSBtYXNraW5nIGVsZW1lbnRcbiAgICBtYXNrOiBmdW5jdGlvbigpIHtcbiAgICAgIHJldHVybiB0aGlzLmRlZnMoKS5wdXQobmV3IFNWRy5NYXNrKVxuICAgIH1cbiAgfVxufSlcblxuXG5TVkcuZXh0ZW5kKFNWRy5FbGVtZW50LCB7XG4gIC8vIERpc3RyaWJ1dGUgbWFzayB0byBzdmcgZWxlbWVudFxuICBtYXNrV2l0aDogZnVuY3Rpb24oZWxlbWVudCkge1xuICAgIC8vIHVzZSBnaXZlbiBtYXNrIG9yIGNyZWF0ZSBhIG5ldyBvbmVcbiAgICB0aGlzLm1hc2tlciA9IGVsZW1lbnQgaW5zdGFuY2VvZiBTVkcuTWFzayA/IGVsZW1lbnQgOiB0aGlzLnBhcmVudCgpLm1hc2soKS5hZGQoZWxlbWVudClcblxuICAgIC8vIHN0b3JlIHJldmVyZW5jZSBvbiBzZWxmIGluIG1hc2tcbiAgICB0aGlzLm1hc2tlci50YXJnZXRzLnB1c2godGhpcylcblxuICAgIC8vIGFwcGx5IG1hc2tcbiAgICByZXR1cm4gdGhpcy5hdHRyKCdtYXNrJywgJ3VybChcIiMnICsgdGhpcy5tYXNrZXIuYXR0cignaWQnKSArICdcIiknKVxuICB9XG4gIC8vIFVubWFzayBlbGVtZW50XG4sIHVubWFzazogZnVuY3Rpb24oKSB7XG4gICAgZGVsZXRlIHRoaXMubWFza2VyXG4gICAgcmV0dXJuIHRoaXMuYXR0cignbWFzaycsIG51bGwpXG4gIH1cblxufSlcblxuU1ZHLkNsaXBQYXRoID0gU1ZHLmludmVudCh7XG4gIC8vIEluaXRpYWxpemUgbm9kZVxuICBjcmVhdGU6IGZ1bmN0aW9uKCkge1xuICAgIHRoaXMuY29uc3RydWN0b3IuY2FsbCh0aGlzLCBTVkcuY3JlYXRlKCdjbGlwUGF0aCcpKVxuXG4gICAgLy8ga2VlcCByZWZlcmVuY2VzIHRvIGNsaXBwZWQgZWxlbWVudHNcbiAgICB0aGlzLnRhcmdldHMgPSBbXVxuICB9XG5cbiAgLy8gSW5oZXJpdCBmcm9tXG4sIGluaGVyaXQ6IFNWRy5Db250YWluZXJcblxuICAvLyBBZGQgY2xhc3MgbWV0aG9kc1xuLCBleHRlbmQ6IHtcbiAgICAvLyBVbmNsaXAgYWxsIGNsaXBwZWQgZWxlbWVudHMgYW5kIHJlbW92ZSBpdHNlbGZcbiAgICByZW1vdmU6IGZ1bmN0aW9uKCkge1xuICAgICAgLy8gdW5jbGlwIGFsbCB0YXJnZXRzXG4gICAgICBmb3IgKHZhciBpID0gdGhpcy50YXJnZXRzLmxlbmd0aCAtIDE7IGkgPj0gMDsgaS0tKVxuICAgICAgICBpZiAodGhpcy50YXJnZXRzW2ldKVxuICAgICAgICAgIHRoaXMudGFyZ2V0c1tpXS51bmNsaXAoKVxuICAgICAgdGhpcy50YXJnZXRzID0gW11cblxuICAgICAgLy8gcmVtb3ZlIGNsaXBQYXRoIGZyb20gcGFyZW50XG4gICAgICB0aGlzLnBhcmVudCgpLnJlbW92ZUVsZW1lbnQodGhpcylcblxuICAgICAgcmV0dXJuIHRoaXNcbiAgICB9XG4gIH1cblxuICAvLyBBZGQgcGFyZW50IG1ldGhvZFxuLCBjb25zdHJ1Y3Q6IHtcbiAgICAvLyBDcmVhdGUgY2xpcHBpbmcgZWxlbWVudFxuICAgIGNsaXA6IGZ1bmN0aW9uKCkge1xuICAgICAgcmV0dXJuIHRoaXMuZGVmcygpLnB1dChuZXcgU1ZHLkNsaXBQYXRoKVxuICAgIH1cbiAgfVxufSlcblxuLy9cblNWRy5leHRlbmQoU1ZHLkVsZW1lbnQsIHtcbiAgLy8gRGlzdHJpYnV0ZSBjbGlwUGF0aCB0byBzdmcgZWxlbWVudFxuICBjbGlwV2l0aDogZnVuY3Rpb24oZWxlbWVudCkge1xuICAgIC8vIHVzZSBnaXZlbiBjbGlwIG9yIGNyZWF0ZSBhIG5ldyBvbmVcbiAgICB0aGlzLmNsaXBwZXIgPSBlbGVtZW50IGluc3RhbmNlb2YgU1ZHLkNsaXBQYXRoID8gZWxlbWVudCA6IHRoaXMucGFyZW50KCkuY2xpcCgpLmFkZChlbGVtZW50KVxuXG4gICAgLy8gc3RvcmUgcmV2ZXJlbmNlIG9uIHNlbGYgaW4gbWFza1xuICAgIHRoaXMuY2xpcHBlci50YXJnZXRzLnB1c2godGhpcylcblxuICAgIC8vIGFwcGx5IG1hc2tcbiAgICByZXR1cm4gdGhpcy5hdHRyKCdjbGlwLXBhdGgnLCAndXJsKFwiIycgKyB0aGlzLmNsaXBwZXIuYXR0cignaWQnKSArICdcIiknKVxuICB9XG4gIC8vIFVuY2xpcCBlbGVtZW50XG4sIHVuY2xpcDogZnVuY3Rpb24oKSB7XG4gICAgZGVsZXRlIHRoaXMuY2xpcHBlclxuICAgIHJldHVybiB0aGlzLmF0dHIoJ2NsaXAtcGF0aCcsIG51bGwpXG4gIH1cblxufSlcblNWRy5HcmFkaWVudCA9IFNWRy5pbnZlbnQoe1xuICAvLyBJbml0aWFsaXplIG5vZGVcbiAgY3JlYXRlOiBmdW5jdGlvbih0eXBlKSB7XG4gICAgdGhpcy5jb25zdHJ1Y3Rvci5jYWxsKHRoaXMsIFNWRy5jcmVhdGUodHlwZSArICdHcmFkaWVudCcpKVxuXG4gICAgLy8gc3RvcmUgdHlwZVxuICAgIHRoaXMudHlwZSA9IHR5cGVcbiAgfVxuXG4gIC8vIEluaGVyaXQgZnJvbVxuLCBpbmhlcml0OiBTVkcuQ29udGFpbmVyXG5cbiAgLy8gQWRkIGNsYXNzIG1ldGhvZHNcbiwgZXh0ZW5kOiB7XG4gICAgLy8gQWRkIGEgY29sb3Igc3RvcFxuICAgIGF0OiBmdW5jdGlvbihvZmZzZXQsIGNvbG9yLCBvcGFjaXR5KSB7XG4gICAgICByZXR1cm4gdGhpcy5wdXQobmV3IFNWRy5TdG9wKS51cGRhdGUob2Zmc2V0LCBjb2xvciwgb3BhY2l0eSlcbiAgICB9XG4gICAgLy8gVXBkYXRlIGdyYWRpZW50XG4gICwgdXBkYXRlOiBmdW5jdGlvbihibG9jaykge1xuICAgICAgLy8gcmVtb3ZlIGFsbCBzdG9wc1xuICAgICAgdGhpcy5jbGVhcigpXG5cbiAgICAgIC8vIGludm9rZSBwYXNzZWQgYmxvY2tcbiAgICAgIGlmICh0eXBlb2YgYmxvY2sgPT0gJ2Z1bmN0aW9uJylcbiAgICAgICAgYmxvY2suY2FsbCh0aGlzLCB0aGlzKVxuXG4gICAgICByZXR1cm4gdGhpc1xuICAgIH1cbiAgICAvLyBSZXR1cm4gdGhlIGZpbGwgaWRcbiAgLCBmaWxsOiBmdW5jdGlvbigpIHtcbiAgICAgIHJldHVybiAndXJsKCMnICsgdGhpcy5pZCgpICsgJyknXG4gICAgfVxuICAgIC8vIEFsaWFzIHN0cmluZyBjb252ZXJ0aW9uIHRvIGZpbGxcbiAgLCB0b1N0cmluZzogZnVuY3Rpb24oKSB7XG4gICAgICByZXR1cm4gdGhpcy5maWxsKClcbiAgICB9XG4gICAgLy8gY3VzdG9tIGF0dHIgdG8gaGFuZGxlIHRyYW5zZm9ybVxuICAsIGF0dHI6IGZ1bmN0aW9uKGEsIGIsIGMpIHtcbiAgICAgIGlmKGEgPT0gJ3RyYW5zZm9ybScpIGEgPSAnZ3JhZGllbnRUcmFuc2Zvcm0nXG4gICAgICByZXR1cm4gU1ZHLkNvbnRhaW5lci5wcm90b3R5cGUuYXR0ci5jYWxsKHRoaXMsIGEsIGIsIGMpXG4gICAgfVxuICB9XG5cbiAgLy8gQWRkIHBhcmVudCBtZXRob2RcbiwgY29uc3RydWN0OiB7XG4gICAgLy8gQ3JlYXRlIGdyYWRpZW50IGVsZW1lbnQgaW4gZGVmc1xuICAgIGdyYWRpZW50OiBmdW5jdGlvbih0eXBlLCBibG9jaykge1xuICAgICAgcmV0dXJuIHRoaXMuZGVmcygpLmdyYWRpZW50KHR5cGUsIGJsb2NrKVxuICAgIH1cbiAgfVxufSlcblxuLy8gQWRkIGFuaW1hdGFibGUgbWV0aG9kcyB0byBib3RoIGdyYWRpZW50IGFuZCBmeCBtb2R1bGVcblNWRy5leHRlbmQoU1ZHLkdyYWRpZW50LCBTVkcuRlgsIHtcbiAgLy8gRnJvbSBwb3NpdGlvblxuICBmcm9tOiBmdW5jdGlvbih4LCB5KSB7XG4gICAgcmV0dXJuICh0aGlzLl90YXJnZXQgfHwgdGhpcykudHlwZSA9PSAncmFkaWFsJyA/XG4gICAgICB0aGlzLmF0dHIoeyBmeDogbmV3IFNWRy5OdW1iZXIoeCksIGZ5OiBuZXcgU1ZHLk51bWJlcih5KSB9KSA6XG4gICAgICB0aGlzLmF0dHIoeyB4MTogbmV3IFNWRy5OdW1iZXIoeCksIHkxOiBuZXcgU1ZHLk51bWJlcih5KSB9KVxuICB9XG4gIC8vIFRvIHBvc2l0aW9uXG4sIHRvOiBmdW5jdGlvbih4LCB5KSB7XG4gICAgcmV0dXJuICh0aGlzLl90YXJnZXQgfHwgdGhpcykudHlwZSA9PSAncmFkaWFsJyA/XG4gICAgICB0aGlzLmF0dHIoeyBjeDogbmV3IFNWRy5OdW1iZXIoeCksIGN5OiBuZXcgU1ZHLk51bWJlcih5KSB9KSA6XG4gICAgICB0aGlzLmF0dHIoeyB4MjogbmV3IFNWRy5OdW1iZXIoeCksIHkyOiBuZXcgU1ZHLk51bWJlcih5KSB9KVxuICB9XG59KVxuXG4vLyBCYXNlIGdyYWRpZW50IGdlbmVyYXRpb25cblNWRy5leHRlbmQoU1ZHLkRlZnMsIHtcbiAgLy8gZGVmaW5lIGdyYWRpZW50XG4gIGdyYWRpZW50OiBmdW5jdGlvbih0eXBlLCBibG9jaykge1xuICAgIHJldHVybiB0aGlzLnB1dChuZXcgU1ZHLkdyYWRpZW50KHR5cGUpKS51cGRhdGUoYmxvY2spXG4gIH1cblxufSlcblxuU1ZHLlN0b3AgPSBTVkcuaW52ZW50KHtcbiAgLy8gSW5pdGlhbGl6ZSBub2RlXG4gIGNyZWF0ZTogJ3N0b3AnXG5cbiAgLy8gSW5oZXJpdCBmcm9tXG4sIGluaGVyaXQ6IFNWRy5FbGVtZW50XG5cbiAgLy8gQWRkIGNsYXNzIG1ldGhvZHNcbiwgZXh0ZW5kOiB7XG4gICAgLy8gYWRkIGNvbG9yIHN0b3BzXG4gICAgdXBkYXRlOiBmdW5jdGlvbihvKSB7XG4gICAgICBpZiAodHlwZW9mIG8gPT0gJ251bWJlcicgfHwgbyBpbnN0YW5jZW9mIFNWRy5OdW1iZXIpIHtcbiAgICAgICAgbyA9IHtcbiAgICAgICAgICBvZmZzZXQ6ICBhcmd1bWVudHNbMF1cbiAgICAgICAgLCBjb2xvcjogICBhcmd1bWVudHNbMV1cbiAgICAgICAgLCBvcGFjaXR5OiBhcmd1bWVudHNbMl1cbiAgICAgICAgfVxuICAgICAgfVxuXG4gICAgICAvLyBzZXQgYXR0cmlidXRlc1xuICAgICAgaWYgKG8ub3BhY2l0eSAhPSBudWxsKSB0aGlzLmF0dHIoJ3N0b3Atb3BhY2l0eScsIG8ub3BhY2l0eSlcbiAgICAgIGlmIChvLmNvbG9yICAgIT0gbnVsbCkgdGhpcy5hdHRyKCdzdG9wLWNvbG9yJywgby5jb2xvcilcbiAgICAgIGlmIChvLm9mZnNldCAgIT0gbnVsbCkgdGhpcy5hdHRyKCdvZmZzZXQnLCBuZXcgU1ZHLk51bWJlcihvLm9mZnNldCkpXG5cbiAgICAgIHJldHVybiB0aGlzXG4gICAgfVxuICB9XG5cbn0pXG5cblNWRy5QYXR0ZXJuID0gU1ZHLmludmVudCh7XG4gIC8vIEluaXRpYWxpemUgbm9kZVxuICBjcmVhdGU6ICdwYXR0ZXJuJ1xuXG4gIC8vIEluaGVyaXQgZnJvbVxuLCBpbmhlcml0OiBTVkcuQ29udGFpbmVyXG5cbiAgLy8gQWRkIGNsYXNzIG1ldGhvZHNcbiwgZXh0ZW5kOiB7XG4gICAgLy8gUmV0dXJuIHRoZSBmaWxsIGlkXG4gICAgZmlsbDogZnVuY3Rpb24oKSB7XG4gICAgICByZXR1cm4gJ3VybCgjJyArIHRoaXMuaWQoKSArICcpJ1xuICAgIH1cbiAgICAvLyBVcGRhdGUgcGF0dGVybiBieSByZWJ1aWxkaW5nXG4gICwgdXBkYXRlOiBmdW5jdGlvbihibG9jaykge1xuICAgICAgLy8gcmVtb3ZlIGNvbnRlbnRcbiAgICAgIHRoaXMuY2xlYXIoKVxuXG4gICAgICAvLyBpbnZva2UgcGFzc2VkIGJsb2NrXG4gICAgICBpZiAodHlwZW9mIGJsb2NrID09ICdmdW5jdGlvbicpXG4gICAgICAgIGJsb2NrLmNhbGwodGhpcywgdGhpcylcblxuICAgICAgcmV0dXJuIHRoaXNcbiAgICB9XG4gICAgLy8gQWxpYXMgc3RyaW5nIGNvbnZlcnRpb24gdG8gZmlsbFxuICAsIHRvU3RyaW5nOiBmdW5jdGlvbigpIHtcbiAgICAgIHJldHVybiB0aGlzLmZpbGwoKVxuICAgIH1cbiAgICAvLyBjdXN0b20gYXR0ciB0byBoYW5kbGUgdHJhbnNmb3JtXG4gICwgYXR0cjogZnVuY3Rpb24oYSwgYiwgYykge1xuICAgICAgaWYoYSA9PSAndHJhbnNmb3JtJykgYSA9ICdwYXR0ZXJuVHJhbnNmb3JtJ1xuICAgICAgcmV0dXJuIFNWRy5Db250YWluZXIucHJvdG90eXBlLmF0dHIuY2FsbCh0aGlzLCBhLCBiLCBjKVxuICAgIH1cblxuICB9XG5cbiAgLy8gQWRkIHBhcmVudCBtZXRob2RcbiwgY29uc3RydWN0OiB7XG4gICAgLy8gQ3JlYXRlIHBhdHRlcm4gZWxlbWVudCBpbiBkZWZzXG4gICAgcGF0dGVybjogZnVuY3Rpb24od2lkdGgsIGhlaWdodCwgYmxvY2spIHtcbiAgICAgIHJldHVybiB0aGlzLmRlZnMoKS5wYXR0ZXJuKHdpZHRoLCBoZWlnaHQsIGJsb2NrKVxuICAgIH1cbiAgfVxufSlcblxuU1ZHLmV4dGVuZChTVkcuRGVmcywge1xuICAvLyBEZWZpbmUgZ3JhZGllbnRcbiAgcGF0dGVybjogZnVuY3Rpb24od2lkdGgsIGhlaWdodCwgYmxvY2spIHtcbiAgICByZXR1cm4gdGhpcy5wdXQobmV3IFNWRy5QYXR0ZXJuKS51cGRhdGUoYmxvY2spLmF0dHIoe1xuICAgICAgeDogICAgICAgICAgICAwXG4gICAgLCB5OiAgICAgICAgICAgIDBcbiAgICAsIHdpZHRoOiAgICAgICAgd2lkdGhcbiAgICAsIGhlaWdodDogICAgICAgaGVpZ2h0XG4gICAgLCBwYXR0ZXJuVW5pdHM6ICd1c2VyU3BhY2VPblVzZSdcbiAgICB9KVxuICB9XG5cbn0pXG5TVkcuRG9jID0gU1ZHLmludmVudCh7XG4gIC8vIEluaXRpYWxpemUgbm9kZVxuICBjcmVhdGU6IGZ1bmN0aW9uKGVsZW1lbnQpIHtcbiAgICBpZiAoZWxlbWVudCkge1xuICAgICAgLy8gZW5zdXJlIHRoZSBwcmVzZW5jZSBvZiBhIGRvbSBlbGVtZW50XG4gICAgICBlbGVtZW50ID0gdHlwZW9mIGVsZW1lbnQgPT0gJ3N0cmluZycgP1xuICAgICAgICBkb2N1bWVudC5nZXRFbGVtZW50QnlJZChlbGVtZW50KSA6XG4gICAgICAgIGVsZW1lbnRcblxuICAgICAgLy8gSWYgdGhlIHRhcmdldCBpcyBhbiBzdmcgZWxlbWVudCwgdXNlIHRoYXQgZWxlbWVudCBhcyB0aGUgbWFpbiB3cmFwcGVyLlxuICAgICAgLy8gVGhpcyBhbGxvd3Mgc3ZnLmpzIHRvIHdvcmsgd2l0aCBzdmcgZG9jdW1lbnRzIGFzIHdlbGwuXG4gICAgICBpZiAoZWxlbWVudC5ub2RlTmFtZSA9PSAnc3ZnJykge1xuICAgICAgICB0aGlzLmNvbnN0cnVjdG9yLmNhbGwodGhpcywgZWxlbWVudClcbiAgICAgIH0gZWxzZSB7XG4gICAgICAgIHRoaXMuY29uc3RydWN0b3IuY2FsbCh0aGlzLCBTVkcuY3JlYXRlKCdzdmcnKSlcbiAgICAgICAgZWxlbWVudC5hcHBlbmRDaGlsZCh0aGlzLm5vZGUpXG4gICAgICAgIHRoaXMuc2l6ZSgnMTAwJScsICcxMDAlJylcbiAgICAgIH1cblxuICAgICAgLy8gc2V0IHN2ZyBlbGVtZW50IGF0dHJpYnV0ZXMgYW5kIGVuc3VyZSBkZWZzIG5vZGVcbiAgICAgIHRoaXMubmFtZXNwYWNlKCkuZGVmcygpXG4gICAgfVxuICB9XG5cbiAgLy8gSW5oZXJpdCBmcm9tXG4sIGluaGVyaXQ6IFNWRy5Db250YWluZXJcblxuICAvLyBBZGQgY2xhc3MgbWV0aG9kc1xuLCBleHRlbmQ6IHtcbiAgICAvLyBBZGQgbmFtZXNwYWNlc1xuICAgIG5hbWVzcGFjZTogZnVuY3Rpb24oKSB7XG4gICAgICByZXR1cm4gdGhpc1xuICAgICAgICAuYXR0cih7IHhtbG5zOiBTVkcubnMsIHZlcnNpb246ICcxLjEnIH0pXG4gICAgICAgIC5hdHRyKCd4bWxuczp4bGluaycsIFNWRy54bGluaywgU1ZHLnhtbG5zKVxuICAgICAgICAuYXR0cigneG1sbnM6c3ZnanMnLCBTVkcuc3ZnanMsIFNWRy54bWxucylcbiAgICB9XG4gICAgLy8gQ3JlYXRlcyBhbmQgcmV0dXJucyBkZWZzIGVsZW1lbnRcbiAgLCBkZWZzOiBmdW5jdGlvbigpIHtcbiAgICAgIGlmICghdGhpcy5fZGVmcykge1xuICAgICAgICB2YXIgZGVmc1xuXG4gICAgICAgIC8vIEZpbmQgb3IgY3JlYXRlIGEgZGVmcyBlbGVtZW50IGluIHRoaXMgaW5zdGFuY2VcbiAgICAgICAgaWYgKGRlZnMgPSB0aGlzLm5vZGUuZ2V0RWxlbWVudHNCeVRhZ05hbWUoJ2RlZnMnKVswXSlcbiAgICAgICAgICB0aGlzLl9kZWZzID0gU1ZHLmFkb3B0KGRlZnMpXG4gICAgICAgIGVsc2VcbiAgICAgICAgICB0aGlzLl9kZWZzID0gbmV3IFNWRy5EZWZzXG5cbiAgICAgICAgLy8gTWFrZSBzdXJlIHRoZSBkZWZzIG5vZGUgaXMgYXQgdGhlIGVuZCBvZiB0aGUgc3RhY2tcbiAgICAgICAgdGhpcy5ub2RlLmFwcGVuZENoaWxkKHRoaXMuX2RlZnMubm9kZSlcbiAgICAgIH1cblxuICAgICAgcmV0dXJuIHRoaXMuX2RlZnNcbiAgICB9XG4gICAgLy8gY3VzdG9tIHBhcmVudCBtZXRob2RcbiAgLCBwYXJlbnQ6IGZ1bmN0aW9uKCkge1xuICAgICAgcmV0dXJuIHRoaXMubm9kZS5wYXJlbnROb2RlLm5vZGVOYW1lID09ICcjZG9jdW1lbnQnID8gbnVsbCA6IHRoaXMubm9kZS5wYXJlbnROb2RlXG4gICAgfVxuICAgIC8vIEZpeCBmb3IgcG9zc2libGUgc3ViLXBpeGVsIG9mZnNldC4gU2VlOlxuICAgIC8vIGh0dHBzOi8vYnVnemlsbGEubW96aWxsYS5vcmcvc2hvd19idWcuY2dpP2lkPTYwODgxMlxuICAsIHNwb2Y6IGZ1bmN0aW9uKHNwb2YpIHtcbiAgICAgIHZhciBwb3MgPSB0aGlzLm5vZGUuZ2V0U2NyZWVuQ1RNKClcblxuICAgICAgaWYgKHBvcylcbiAgICAgICAgdGhpc1xuICAgICAgICAgIC5zdHlsZSgnbGVmdCcsICgtcG9zLmUgJSAxKSArICdweCcpXG4gICAgICAgICAgLnN0eWxlKCd0b3AnLCAgKC1wb3MuZiAlIDEpICsgJ3B4JylcblxuICAgICAgcmV0dXJuIHRoaXNcbiAgICB9XG5cbiAgICAgIC8vIFJlbW92ZXMgdGhlIGRvYyBmcm9tIHRoZSBET01cbiAgLCByZW1vdmU6IGZ1bmN0aW9uKCkge1xuICAgICAgaWYodGhpcy5wYXJlbnQoKSkge1xuICAgICAgICB0aGlzLnBhcmVudCgpLnJlbW92ZUNoaWxkKHRoaXMubm9kZSk7XG4gICAgICB9XG5cbiAgICAgIHJldHVybiB0aGlzO1xuICAgIH1cbiAgfVxuXG59KVxuXG5TVkcuU2hhcGUgPSBTVkcuaW52ZW50KHtcbiAgLy8gSW5pdGlhbGl6ZSBub2RlXG4gIGNyZWF0ZTogZnVuY3Rpb24oZWxlbWVudCkge1xuICAgIHRoaXMuY29uc3RydWN0b3IuY2FsbCh0aGlzLCBlbGVtZW50KVxuICB9XG5cbiAgLy8gSW5oZXJpdCBmcm9tXG4sIGluaGVyaXQ6IFNWRy5FbGVtZW50XG5cbn0pXG5cblNWRy5CYXJlID0gU1ZHLmludmVudCh7XG4gIC8vIEluaXRpYWxpemVcbiAgY3JlYXRlOiBmdW5jdGlvbihlbGVtZW50LCBpbmhlcml0KSB7XG4gICAgLy8gY29uc3RydWN0IGVsZW1lbnRcbiAgICB0aGlzLmNvbnN0cnVjdG9yLmNhbGwodGhpcywgU1ZHLmNyZWF0ZShlbGVtZW50KSlcblxuICAgIC8vIGluaGVyaXQgY3VzdG9tIG1ldGhvZHNcbiAgICBpZiAoaW5oZXJpdClcbiAgICAgIGZvciAodmFyIG1ldGhvZCBpbiBpbmhlcml0LnByb3RvdHlwZSlcbiAgICAgICAgaWYgKHR5cGVvZiBpbmhlcml0LnByb3RvdHlwZVttZXRob2RdID09PSAnZnVuY3Rpb24nKVxuICAgICAgICAgIHRoaXNbbWV0aG9kXSA9IGluaGVyaXQucHJvdG90eXBlW21ldGhvZF1cbiAgfVxuXG4gIC8vIEluaGVyaXQgZnJvbVxuLCBpbmhlcml0OiBTVkcuRWxlbWVudFxuXG4gIC8vIEFkZCBtZXRob2RzXG4sIGV4dGVuZDoge1xuICAgIC8vIEluc2VydCBzb21lIHBsYWluIHRleHRcbiAgICB3b3JkczogZnVuY3Rpb24odGV4dCkge1xuICAgICAgLy8gcmVtb3ZlIGNvbnRlbnRzXG4gICAgICB3aGlsZSAodGhpcy5ub2RlLmhhc0NoaWxkTm9kZXMoKSlcbiAgICAgICAgdGhpcy5ub2RlLnJlbW92ZUNoaWxkKHRoaXMubm9kZS5sYXN0Q2hpbGQpXG5cbiAgICAgIC8vIGNyZWF0ZSB0ZXh0IG5vZGVcbiAgICAgIHRoaXMubm9kZS5hcHBlbmRDaGlsZChkb2N1bWVudC5jcmVhdGVUZXh0Tm9kZSh0ZXh0KSlcblxuICAgICAgcmV0dXJuIHRoaXNcbiAgICB9XG4gIH1cbn0pXG5cblxuU1ZHLmV4dGVuZChTVkcuUGFyZW50LCB7XG4gIC8vIENyZWF0ZSBhbiBlbGVtZW50IHRoYXQgaXMgbm90IGRlc2NyaWJlZCBieSBTVkcuanNcbiAgZWxlbWVudDogZnVuY3Rpb24oZWxlbWVudCwgaW5oZXJpdCkge1xuICAgIHJldHVybiB0aGlzLnB1dChuZXcgU1ZHLkJhcmUoZWxlbWVudCwgaW5oZXJpdCkpXG4gIH1cbiAgLy8gQWRkIHN5bWJvbCBlbGVtZW50XG4sIHN5bWJvbDogZnVuY3Rpb24oKSB7XG4gICAgcmV0dXJuIHRoaXMuZGVmcygpLmVsZW1lbnQoJ3N5bWJvbCcsIFNWRy5Db250YWluZXIpXG4gIH1cblxufSlcblNWRy5Vc2UgPSBTVkcuaW52ZW50KHtcbiAgLy8gSW5pdGlhbGl6ZSBub2RlXG4gIGNyZWF0ZTogJ3VzZSdcblxuICAvLyBJbmhlcml0IGZyb21cbiwgaW5oZXJpdDogU1ZHLlNoYXBlXG5cbiAgLy8gQWRkIGNsYXNzIG1ldGhvZHNcbiwgZXh0ZW5kOiB7XG4gICAgLy8gVXNlIGVsZW1lbnQgYXMgYSByZWZlcmVuY2VcbiAgICBlbGVtZW50OiBmdW5jdGlvbihlbGVtZW50LCBmaWxlKSB7XG4gICAgICAvLyBTZXQgbGluZWQgZWxlbWVudFxuICAgICAgcmV0dXJuIHRoaXMuYXR0cignaHJlZicsIChmaWxlIHx8ICcnKSArICcjJyArIGVsZW1lbnQsIFNWRy54bGluaylcbiAgICB9XG4gIH1cblxuICAvLyBBZGQgcGFyZW50IG1ldGhvZFxuLCBjb25zdHJ1Y3Q6IHtcbiAgICAvLyBDcmVhdGUgYSB1c2UgZWxlbWVudFxuICAgIHVzZTogZnVuY3Rpb24oZWxlbWVudCwgZmlsZSkge1xuICAgICAgcmV0dXJuIHRoaXMucHV0KG5ldyBTVkcuVXNlKS5lbGVtZW50KGVsZW1lbnQsIGZpbGUpXG4gICAgfVxuICB9XG59KVxuU1ZHLlJlY3QgPSBTVkcuaW52ZW50KHtcbiAgLy8gSW5pdGlhbGl6ZSBub2RlXG4gIGNyZWF0ZTogJ3JlY3QnXG5cbiAgLy8gSW5oZXJpdCBmcm9tXG4sIGluaGVyaXQ6IFNWRy5TaGFwZVxuXG4gIC8vIEFkZCBwYXJlbnQgbWV0aG9kXG4sIGNvbnN0cnVjdDoge1xuICAgIC8vIENyZWF0ZSBhIHJlY3QgZWxlbWVudFxuICAgIHJlY3Q6IGZ1bmN0aW9uKHdpZHRoLCBoZWlnaHQpIHtcbiAgICAgIHJldHVybiB0aGlzLnB1dChuZXcgU1ZHLlJlY3QoKSkuc2l6ZSh3aWR0aCwgaGVpZ2h0KVxuICAgIH1cbiAgfVxufSlcblNWRy5DaXJjbGUgPSBTVkcuaW52ZW50KHtcbiAgLy8gSW5pdGlhbGl6ZSBub2RlXG4gIGNyZWF0ZTogJ2NpcmNsZSdcblxuICAvLyBJbmhlcml0IGZyb21cbiwgaW5oZXJpdDogU1ZHLlNoYXBlXG5cbiAgLy8gQWRkIHBhcmVudCBtZXRob2RcbiwgY29uc3RydWN0OiB7XG4gICAgLy8gQ3JlYXRlIGNpcmNsZSBlbGVtZW50LCBiYXNlZCBvbiBlbGxpcHNlXG4gICAgY2lyY2xlOiBmdW5jdGlvbihzaXplKSB7XG4gICAgICByZXR1cm4gdGhpcy5wdXQobmV3IFNWRy5DaXJjbGUpLnJ4KG5ldyBTVkcuTnVtYmVyKHNpemUpLmRpdmlkZSgyKSkubW92ZSgwLCAwKVxuICAgIH1cbiAgfVxufSlcblxuU1ZHLmV4dGVuZChTVkcuQ2lyY2xlLCBTVkcuRlgsIHtcbiAgLy8gUmFkaXVzIHggdmFsdWVcbiAgcng6IGZ1bmN0aW9uKHJ4KSB7XG4gICAgcmV0dXJuIHRoaXMuYXR0cigncicsIHJ4KVxuICB9XG4gIC8vIEFsaWFzIHJhZGl1cyB4IHZhbHVlXG4sIHJ5OiBmdW5jdGlvbihyeSkge1xuICAgIHJldHVybiB0aGlzLnJ4KHJ5KVxuICB9XG59KVxuXG5TVkcuRWxsaXBzZSA9IFNWRy5pbnZlbnQoe1xuICAvLyBJbml0aWFsaXplIG5vZGVcbiAgY3JlYXRlOiAnZWxsaXBzZSdcblxuICAvLyBJbmhlcml0IGZyb21cbiwgaW5oZXJpdDogU1ZHLlNoYXBlXG5cbiAgLy8gQWRkIHBhcmVudCBtZXRob2RcbiwgY29uc3RydWN0OiB7XG4gICAgLy8gQ3JlYXRlIGFuIGVsbGlwc2VcbiAgICBlbGxpcHNlOiBmdW5jdGlvbih3aWR0aCwgaGVpZ2h0KSB7XG4gICAgICByZXR1cm4gdGhpcy5wdXQobmV3IFNWRy5FbGxpcHNlKS5zaXplKHdpZHRoLCBoZWlnaHQpLm1vdmUoMCwgMClcbiAgICB9XG4gIH1cbn0pXG5cblNWRy5leHRlbmQoU1ZHLkVsbGlwc2UsIFNWRy5SZWN0LCBTVkcuRlgsIHtcbiAgLy8gUmFkaXVzIHggdmFsdWVcbiAgcng6IGZ1bmN0aW9uKHJ4KSB7XG4gICAgcmV0dXJuIHRoaXMuYXR0cigncngnLCByeClcbiAgfVxuICAvLyBSYWRpdXMgeSB2YWx1ZVxuLCByeTogZnVuY3Rpb24ocnkpIHtcbiAgICByZXR1cm4gdGhpcy5hdHRyKCdyeScsIHJ5KVxuICB9XG59KVxuXG4vLyBBZGQgY29tbW9uIG1ldGhvZFxuU1ZHLmV4dGVuZChTVkcuQ2lyY2xlLCBTVkcuRWxsaXBzZSwge1xuICAgIC8vIE1vdmUgb3ZlciB4LWF4aXNcbiAgICB4OiBmdW5jdGlvbih4KSB7XG4gICAgICByZXR1cm4geCA9PSBudWxsID8gdGhpcy5jeCgpIC0gdGhpcy5yeCgpIDogdGhpcy5jeCh4ICsgdGhpcy5yeCgpKVxuICAgIH1cbiAgICAvLyBNb3ZlIG92ZXIgeS1heGlzXG4gICwgeTogZnVuY3Rpb24oeSkge1xuICAgICAgcmV0dXJuIHkgPT0gbnVsbCA/IHRoaXMuY3koKSAtIHRoaXMucnkoKSA6IHRoaXMuY3koeSArIHRoaXMucnkoKSlcbiAgICB9XG4gICAgLy8gTW92ZSBieSBjZW50ZXIgb3ZlciB4LWF4aXNcbiAgLCBjeDogZnVuY3Rpb24oeCkge1xuICAgICAgcmV0dXJuIHggPT0gbnVsbCA/IHRoaXMuYXR0cignY3gnKSA6IHRoaXMuYXR0cignY3gnLCB4KVxuICAgIH1cbiAgICAvLyBNb3ZlIGJ5IGNlbnRlciBvdmVyIHktYXhpc1xuICAsIGN5OiBmdW5jdGlvbih5KSB7XG4gICAgICByZXR1cm4geSA9PSBudWxsID8gdGhpcy5hdHRyKCdjeScpIDogdGhpcy5hdHRyKCdjeScsIHkpXG4gICAgfVxuICAgIC8vIFNldCB3aWR0aCBvZiBlbGVtZW50XG4gICwgd2lkdGg6IGZ1bmN0aW9uKHdpZHRoKSB7XG4gICAgICByZXR1cm4gd2lkdGggPT0gbnVsbCA/IHRoaXMucngoKSAqIDIgOiB0aGlzLnJ4KG5ldyBTVkcuTnVtYmVyKHdpZHRoKS5kaXZpZGUoMikpXG4gICAgfVxuICAgIC8vIFNldCBoZWlnaHQgb2YgZWxlbWVudFxuICAsIGhlaWdodDogZnVuY3Rpb24oaGVpZ2h0KSB7XG4gICAgICByZXR1cm4gaGVpZ2h0ID09IG51bGwgPyB0aGlzLnJ5KCkgKiAyIDogdGhpcy5yeShuZXcgU1ZHLk51bWJlcihoZWlnaHQpLmRpdmlkZSgyKSlcbiAgICB9XG4gICAgLy8gQ3VzdG9tIHNpemUgZnVuY3Rpb25cbiAgLCBzaXplOiBmdW5jdGlvbih3aWR0aCwgaGVpZ2h0KSB7XG4gICAgICB2YXIgcCA9IHByb3BvcnRpb25hbFNpemUodGhpcywgd2lkdGgsIGhlaWdodClcblxuICAgICAgcmV0dXJuIHRoaXNcbiAgICAgICAgLnJ4KG5ldyBTVkcuTnVtYmVyKHAud2lkdGgpLmRpdmlkZSgyKSlcbiAgICAgICAgLnJ5KG5ldyBTVkcuTnVtYmVyKHAuaGVpZ2h0KS5kaXZpZGUoMikpXG4gICAgfVxufSlcblNWRy5MaW5lID0gU1ZHLmludmVudCh7XG4gIC8vIEluaXRpYWxpemUgbm9kZVxuICBjcmVhdGU6ICdsaW5lJ1xuXG4gIC8vIEluaGVyaXQgZnJvbVxuLCBpbmhlcml0OiBTVkcuU2hhcGVcblxuICAvLyBBZGQgY2xhc3MgbWV0aG9kc1xuLCBleHRlbmQ6IHtcbiAgICAvLyBHZXQgYXJyYXlcbiAgICBhcnJheTogZnVuY3Rpb24oKSB7XG4gICAgICByZXR1cm4gbmV3IFNWRy5Qb2ludEFycmF5KFtcbiAgICAgICAgWyB0aGlzLmF0dHIoJ3gxJyksIHRoaXMuYXR0cigneTEnKSBdXG4gICAgICAsIFsgdGhpcy5hdHRyKCd4MicpLCB0aGlzLmF0dHIoJ3kyJykgXVxuICAgICAgXSlcbiAgICB9XG4gICAgLy8gT3ZlcndyaXRlIG5hdGl2ZSBwbG90KCkgbWV0aG9kXG4gICwgcGxvdDogZnVuY3Rpb24oeDEsIHkxLCB4MiwgeTIpIHtcbiAgICAgIGlmICh0eXBlb2YgeTEgIT09ICd1bmRlZmluZWQnKVxuICAgICAgICB4MSA9IHsgeDE6IHgxLCB5MTogeTEsIHgyOiB4MiwgeTI6IHkyIH1cbiAgICAgIGVsc2VcbiAgICAgICAgeDEgPSBuZXcgU1ZHLlBvaW50QXJyYXkoeDEpLnRvTGluZSgpXG5cbiAgICAgIHJldHVybiB0aGlzLmF0dHIoeDEpXG4gICAgfVxuICAgIC8vIE1vdmUgYnkgbGVmdCB0b3AgY29ybmVyXG4gICwgbW92ZTogZnVuY3Rpb24oeCwgeSkge1xuICAgICAgcmV0dXJuIHRoaXMuYXR0cih0aGlzLmFycmF5KCkubW92ZSh4LCB5KS50b0xpbmUoKSlcbiAgICB9XG4gICAgLy8gU2V0IGVsZW1lbnQgc2l6ZSB0byBnaXZlbiB3aWR0aCBhbmQgaGVpZ2h0XG4gICwgc2l6ZTogZnVuY3Rpb24od2lkdGgsIGhlaWdodCkge1xuICAgICAgdmFyIHAgPSBwcm9wb3J0aW9uYWxTaXplKHRoaXMsIHdpZHRoLCBoZWlnaHQpXG5cbiAgICAgIHJldHVybiB0aGlzLmF0dHIodGhpcy5hcnJheSgpLnNpemUocC53aWR0aCwgcC5oZWlnaHQpLnRvTGluZSgpKVxuICAgIH1cbiAgfVxuXG4gIC8vIEFkZCBwYXJlbnQgbWV0aG9kXG4sIGNvbnN0cnVjdDoge1xuICAgIC8vIENyZWF0ZSBhIGxpbmUgZWxlbWVudFxuICAgIGxpbmU6IGZ1bmN0aW9uKHgxLCB5MSwgeDIsIHkyKSB7XG4gICAgICByZXR1cm4gdGhpcy5wdXQobmV3IFNWRy5MaW5lKS5wbG90KHgxLCB5MSwgeDIsIHkyKVxuICAgIH1cbiAgfVxufSlcblxuU1ZHLlBvbHlsaW5lID0gU1ZHLmludmVudCh7XG4gIC8vIEluaXRpYWxpemUgbm9kZVxuICBjcmVhdGU6ICdwb2x5bGluZSdcblxuICAvLyBJbmhlcml0IGZyb21cbiwgaW5oZXJpdDogU1ZHLlNoYXBlXG5cbiAgLy8gQWRkIHBhcmVudCBtZXRob2RcbiwgY29uc3RydWN0OiB7XG4gICAgLy8gQ3JlYXRlIGEgd3JhcHBlZCBwb2x5bGluZSBlbGVtZW50XG4gICAgcG9seWxpbmU6IGZ1bmN0aW9uKHApIHtcbiAgICAgIHJldHVybiB0aGlzLnB1dChuZXcgU1ZHLlBvbHlsaW5lKS5wbG90KHApXG4gICAgfVxuICB9XG59KVxuXG5TVkcuUG9seWdvbiA9IFNWRy5pbnZlbnQoe1xuICAvLyBJbml0aWFsaXplIG5vZGVcbiAgY3JlYXRlOiAncG9seWdvbidcblxuICAvLyBJbmhlcml0IGZyb21cbiwgaW5oZXJpdDogU1ZHLlNoYXBlXG5cbiAgLy8gQWRkIHBhcmVudCBtZXRob2RcbiwgY29uc3RydWN0OiB7XG4gICAgLy8gQ3JlYXRlIGEgd3JhcHBlZCBwb2x5Z29uIGVsZW1lbnRcbiAgICBwb2x5Z29uOiBmdW5jdGlvbihwKSB7XG4gICAgICByZXR1cm4gdGhpcy5wdXQobmV3IFNWRy5Qb2x5Z29uKS5wbG90KHApXG4gICAgfVxuICB9XG59KVxuXG4vLyBBZGQgcG9seWdvbi1zcGVjaWZpYyBmdW5jdGlvbnNcblNWRy5leHRlbmQoU1ZHLlBvbHlsaW5lLCBTVkcuUG9seWdvbiwge1xuICAvLyBHZXQgYXJyYXlcbiAgYXJyYXk6IGZ1bmN0aW9uKCkge1xuICAgIHJldHVybiB0aGlzLl9hcnJheSB8fCAodGhpcy5fYXJyYXkgPSBuZXcgU1ZHLlBvaW50QXJyYXkodGhpcy5hdHRyKCdwb2ludHMnKSkpXG4gIH1cbiAgLy8gUGxvdCBuZXcgcGF0aFxuLCBwbG90OiBmdW5jdGlvbihwKSB7XG4gICAgcmV0dXJuIHRoaXMuYXR0cigncG9pbnRzJywgKHRoaXMuX2FycmF5ID0gbmV3IFNWRy5Qb2ludEFycmF5KHApKSlcbiAgfVxuICAvLyBNb3ZlIGJ5IGxlZnQgdG9wIGNvcm5lclxuLCBtb3ZlOiBmdW5jdGlvbih4LCB5KSB7XG4gICAgcmV0dXJuIHRoaXMuYXR0cigncG9pbnRzJywgdGhpcy5hcnJheSgpLm1vdmUoeCwgeSkpXG4gIH1cbiAgLy8gU2V0IGVsZW1lbnQgc2l6ZSB0byBnaXZlbiB3aWR0aCBhbmQgaGVpZ2h0XG4sIHNpemU6IGZ1bmN0aW9uKHdpZHRoLCBoZWlnaHQpIHtcbiAgICB2YXIgcCA9IHByb3BvcnRpb25hbFNpemUodGhpcywgd2lkdGgsIGhlaWdodClcblxuICAgIHJldHVybiB0aGlzLmF0dHIoJ3BvaW50cycsIHRoaXMuYXJyYXkoKS5zaXplKHAud2lkdGgsIHAuaGVpZ2h0KSlcbiAgfVxuXG59KVxuLy8gdW5pZnkgYWxsIHBvaW50IHRvIHBvaW50IGVsZW1lbnRzXG5TVkcuZXh0ZW5kKFNWRy5MaW5lLCBTVkcuUG9seWxpbmUsIFNWRy5Qb2x5Z29uLCB7XG4gIC8vIERlZmluZSBtb3JwaGFibGUgYXJyYXlcbiAgbW9ycGhBcnJheTogIFNWRy5Qb2ludEFycmF5XG4gIC8vIE1vdmUgYnkgbGVmdCB0b3AgY29ybmVyIG92ZXIgeC1heGlzXG4sIHg6IGZ1bmN0aW9uKHgpIHtcbiAgICByZXR1cm4geCA9PSBudWxsID8gdGhpcy5iYm94KCkueCA6IHRoaXMubW92ZSh4LCB0aGlzLmJib3goKS55KVxuICB9XG4gIC8vIE1vdmUgYnkgbGVmdCB0b3AgY29ybmVyIG92ZXIgeS1heGlzXG4sIHk6IGZ1bmN0aW9uKHkpIHtcbiAgICByZXR1cm4geSA9PSBudWxsID8gdGhpcy5iYm94KCkueSA6IHRoaXMubW92ZSh0aGlzLmJib3goKS54LCB5KVxuICB9XG4gIC8vIFNldCB3aWR0aCBvZiBlbGVtZW50XG4sIHdpZHRoOiBmdW5jdGlvbih3aWR0aCkge1xuICAgIHZhciBiID0gdGhpcy5iYm94KClcblxuICAgIHJldHVybiB3aWR0aCA9PSBudWxsID8gYi53aWR0aCA6IHRoaXMuc2l6ZSh3aWR0aCwgYi5oZWlnaHQpXG4gIH1cbiAgLy8gU2V0IGhlaWdodCBvZiBlbGVtZW50XG4sIGhlaWdodDogZnVuY3Rpb24oaGVpZ2h0KSB7XG4gICAgdmFyIGIgPSB0aGlzLmJib3goKVxuXG4gICAgcmV0dXJuIGhlaWdodCA9PSBudWxsID8gYi5oZWlnaHQgOiB0aGlzLnNpemUoYi53aWR0aCwgaGVpZ2h0KVxuICB9XG59KVxuU1ZHLlBhdGggPSBTVkcuaW52ZW50KHtcbiAgLy8gSW5pdGlhbGl6ZSBub2RlXG4gIGNyZWF0ZTogJ3BhdGgnXG5cbiAgLy8gSW5oZXJpdCBmcm9tXG4sIGluaGVyaXQ6IFNWRy5TaGFwZVxuXG4gIC8vIEFkZCBjbGFzcyBtZXRob2RzXG4sIGV4dGVuZDoge1xuICAgIC8vIERlZmluZSBtb3JwaGFibGUgYXJyYXlcbiAgICBtb3JwaEFycmF5OiAgU1ZHLlBhdGhBcnJheVxuICAgIC8vIEdldCBhcnJheVxuICAsIGFycmF5OiBmdW5jdGlvbigpIHtcbiAgICAgIHJldHVybiB0aGlzLl9hcnJheSB8fCAodGhpcy5fYXJyYXkgPSBuZXcgU1ZHLlBhdGhBcnJheSh0aGlzLmF0dHIoJ2QnKSkpXG4gICAgfVxuICAgIC8vIFBsb3QgbmV3IHBvbHkgcG9pbnRzXG4gICwgcGxvdDogZnVuY3Rpb24ocCkge1xuICAgICAgcmV0dXJuIHRoaXMuYXR0cignZCcsICh0aGlzLl9hcnJheSA9IG5ldyBTVkcuUGF0aEFycmF5KHApKSlcbiAgICB9XG4gICAgLy8gTW92ZSBieSBsZWZ0IHRvcCBjb3JuZXJcbiAgLCBtb3ZlOiBmdW5jdGlvbih4LCB5KSB7XG4gICAgICByZXR1cm4gdGhpcy5hdHRyKCdkJywgdGhpcy5hcnJheSgpLm1vdmUoeCwgeSkpXG4gICAgfVxuICAgIC8vIE1vdmUgYnkgbGVmdCB0b3AgY29ybmVyIG92ZXIgeC1heGlzXG4gICwgeDogZnVuY3Rpb24oeCkge1xuICAgICAgcmV0dXJuIHggPT0gbnVsbCA/IHRoaXMuYmJveCgpLnggOiB0aGlzLm1vdmUoeCwgdGhpcy5iYm94KCkueSlcbiAgICB9XG4gICAgLy8gTW92ZSBieSBsZWZ0IHRvcCBjb3JuZXIgb3ZlciB5LWF4aXNcbiAgLCB5OiBmdW5jdGlvbih5KSB7XG4gICAgICByZXR1cm4geSA9PSBudWxsID8gdGhpcy5iYm94KCkueSA6IHRoaXMubW92ZSh0aGlzLmJib3goKS54LCB5KVxuICAgIH1cbiAgICAvLyBTZXQgZWxlbWVudCBzaXplIHRvIGdpdmVuIHdpZHRoIGFuZCBoZWlnaHRcbiAgLCBzaXplOiBmdW5jdGlvbih3aWR0aCwgaGVpZ2h0KSB7XG4gICAgICB2YXIgcCA9IHByb3BvcnRpb25hbFNpemUodGhpcywgd2lkdGgsIGhlaWdodClcblxuICAgICAgcmV0dXJuIHRoaXMuYXR0cignZCcsIHRoaXMuYXJyYXkoKS5zaXplKHAud2lkdGgsIHAuaGVpZ2h0KSlcbiAgICB9XG4gICAgLy8gU2V0IHdpZHRoIG9mIGVsZW1lbnRcbiAgLCB3aWR0aDogZnVuY3Rpb24od2lkdGgpIHtcbiAgICAgIHJldHVybiB3aWR0aCA9PSBudWxsID8gdGhpcy5iYm94KCkud2lkdGggOiB0aGlzLnNpemUod2lkdGgsIHRoaXMuYmJveCgpLmhlaWdodClcbiAgICB9XG4gICAgLy8gU2V0IGhlaWdodCBvZiBlbGVtZW50XG4gICwgaGVpZ2h0OiBmdW5jdGlvbihoZWlnaHQpIHtcbiAgICAgIHJldHVybiBoZWlnaHQgPT0gbnVsbCA/IHRoaXMuYmJveCgpLmhlaWdodCA6IHRoaXMuc2l6ZSh0aGlzLmJib3goKS53aWR0aCwgaGVpZ2h0KVxuICAgIH1cblxuICB9XG5cbiAgLy8gQWRkIHBhcmVudCBtZXRob2RcbiwgY29uc3RydWN0OiB7XG4gICAgLy8gQ3JlYXRlIGEgd3JhcHBlZCBwYXRoIGVsZW1lbnRcbiAgICBwYXRoOiBmdW5jdGlvbihkKSB7XG4gICAgICByZXR1cm4gdGhpcy5wdXQobmV3IFNWRy5QYXRoKS5wbG90KGQpXG4gICAgfVxuICB9XG59KVxuU1ZHLkltYWdlID0gU1ZHLmludmVudCh7XG4gIC8vIEluaXRpYWxpemUgbm9kZVxuICBjcmVhdGU6ICdpbWFnZSdcblxuICAvLyBJbmhlcml0IGZyb21cbiwgaW5oZXJpdDogU1ZHLlNoYXBlXG5cbiAgLy8gQWRkIGNsYXNzIG1ldGhvZHNcbiwgZXh0ZW5kOiB7XG4gICAgLy8gKHJlKWxvYWQgaW1hZ2VcbiAgICBsb2FkOiBmdW5jdGlvbih1cmwpIHtcbiAgICAgIGlmICghdXJsKSByZXR1cm4gdGhpc1xuXG4gICAgICB2YXIgc2VsZiA9IHRoaXNcbiAgICAgICAgLCBpbWcgID0gZG9jdW1lbnQuY3JlYXRlRWxlbWVudCgnaW1nJylcblxuICAgICAgLy8gcHJlbG9hZCBpbWFnZVxuICAgICAgaW1nLm9ubG9hZCA9IGZ1bmN0aW9uKCkge1xuICAgICAgICB2YXIgcCA9IHNlbGYucGFyZW50KFNWRy5QYXR0ZXJuKVxuXG4gICAgICAgIGlmKHAgPT09IG51bGwpIHJldHVyblxuXG4gICAgICAgIC8vIGVuc3VyZSBpbWFnZSBzaXplXG4gICAgICAgIGlmIChzZWxmLndpZHRoKCkgPT0gMCAmJiBzZWxmLmhlaWdodCgpID09IDApXG4gICAgICAgICAgc2VsZi5zaXplKGltZy53aWR0aCwgaW1nLmhlaWdodClcblxuICAgICAgICAvLyBlbnN1cmUgcGF0dGVybiBzaXplIGlmIG5vdCBzZXRcbiAgICAgICAgaWYgKHAgJiYgcC53aWR0aCgpID09IDAgJiYgcC5oZWlnaHQoKSA9PSAwKVxuICAgICAgICAgIHAuc2l6ZShzZWxmLndpZHRoKCksIHNlbGYuaGVpZ2h0KCkpXG5cbiAgICAgICAgLy8gY2FsbGJhY2tcbiAgICAgICAgaWYgKHR5cGVvZiBzZWxmLl9sb2FkZWQgPT09ICdmdW5jdGlvbicpXG4gICAgICAgICAgc2VsZi5fbG9hZGVkLmNhbGwoc2VsZiwge1xuICAgICAgICAgICAgd2lkdGg6ICBpbWcud2lkdGhcbiAgICAgICAgICAsIGhlaWdodDogaW1nLmhlaWdodFxuICAgICAgICAgICwgcmF0aW86ICBpbWcud2lkdGggLyBpbWcuaGVpZ2h0XG4gICAgICAgICAgLCB1cmw6ICAgIHVybFxuICAgICAgICAgIH0pXG4gICAgICB9XG5cbiAgICAgIGltZy5vbmVycm9yID0gZnVuY3Rpb24oZSl7XG4gICAgICAgIGlmICh0eXBlb2Ygc2VsZi5fZXJyb3IgPT09ICdmdW5jdGlvbicpe1xuICAgICAgICAgICAgc2VsZi5fZXJyb3IuY2FsbChzZWxmLCBlKVxuICAgICAgICB9XG4gICAgICB9XG5cbiAgICAgIHJldHVybiB0aGlzLmF0dHIoJ2hyZWYnLCAoaW1nLnNyYyA9IHRoaXMuc3JjID0gdXJsKSwgU1ZHLnhsaW5rKVxuICAgIH1cbiAgICAvLyBBZGQgbG9hZGVkIGNhbGxiYWNrXG4gICwgbG9hZGVkOiBmdW5jdGlvbihsb2FkZWQpIHtcbiAgICAgIHRoaXMuX2xvYWRlZCA9IGxvYWRlZFxuICAgICAgcmV0dXJuIHRoaXNcbiAgICB9XG5cbiAgLCBlcnJvcjogZnVuY3Rpb24oZXJyb3IpIHtcbiAgICAgIHRoaXMuX2Vycm9yID0gZXJyb3JcbiAgICAgIHJldHVybiB0aGlzXG4gICAgfVxuICB9XG5cbiAgLy8gQWRkIHBhcmVudCBtZXRob2RcbiwgY29uc3RydWN0OiB7XG4gICAgLy8gY3JlYXRlIGltYWdlIGVsZW1lbnQsIGxvYWQgaW1hZ2UgYW5kIHNldCBpdHMgc2l6ZVxuICAgIGltYWdlOiBmdW5jdGlvbihzb3VyY2UsIHdpZHRoLCBoZWlnaHQpIHtcbiAgICAgIHJldHVybiB0aGlzLnB1dChuZXcgU1ZHLkltYWdlKS5sb2FkKHNvdXJjZSkuc2l6ZSh3aWR0aCB8fCAwLCBoZWlnaHQgfHwgd2lkdGggfHwgMClcbiAgICB9XG4gIH1cblxufSlcblNWRy5UZXh0ID0gU1ZHLmludmVudCh7XG4gIC8vIEluaXRpYWxpemUgbm9kZVxuICBjcmVhdGU6IGZ1bmN0aW9uKCkge1xuICAgIHRoaXMuY29uc3RydWN0b3IuY2FsbCh0aGlzLCBTVkcuY3JlYXRlKCd0ZXh0JykpXG5cbiAgICB0aGlzLmRvbS5sZWFkaW5nID0gbmV3IFNWRy5OdW1iZXIoMS4zKSAgICAvLyBzdG9yZSBsZWFkaW5nIHZhbHVlIGZvciByZWJ1aWxkaW5nXG4gICAgdGhpcy5fcmVidWlsZCA9IHRydWUgICAgICAgICAgICAgICAgICAgICAgLy8gZW5hYmxlIGF1dG9tYXRpYyB1cGRhdGluZyBvZiBkeSB2YWx1ZXNcbiAgICB0aGlzLl9idWlsZCAgID0gZmFsc2UgICAgICAgICAgICAgICAgICAgICAvLyBkaXNhYmxlIGJ1aWxkIG1vZGUgZm9yIGFkZGluZyBtdWx0aXBsZSBsaW5lc1xuXG4gICAgLy8gc2V0IGRlZmF1bHQgZm9udFxuICAgIHRoaXMuYXR0cignZm9udC1mYW1pbHknLCBTVkcuZGVmYXVsdHMuYXR0cnNbJ2ZvbnQtZmFtaWx5J10pXG4gIH1cblxuICAvLyBJbmhlcml0IGZyb21cbiwgaW5oZXJpdDogU1ZHLlNoYXBlXG5cbiAgLy8gQWRkIGNsYXNzIG1ldGhvZHNcbiwgZXh0ZW5kOiB7XG4gICAgLy8gTW92ZSBvdmVyIHgtYXhpc1xuICAgIHg6IGZ1bmN0aW9uKHgpIHtcbiAgICAgIC8vIGFjdCBhcyBnZXR0ZXJcbiAgICAgIGlmICh4ID09IG51bGwpXG4gICAgICAgIHJldHVybiB0aGlzLmF0dHIoJ3gnKVxuXG4gICAgICAvLyBtb3ZlIGxpbmVzIGFzIHdlbGwgaWYgbm8gdGV4dFBhdGggaXMgcHJlc2VudFxuICAgICAgaWYgKCF0aGlzLnRleHRQYXRoKVxuICAgICAgICB0aGlzLmxpbmVzKCkuZWFjaChmdW5jdGlvbigpIHsgaWYgKHRoaXMuZG9tLm5ld0xpbmVkKSB0aGlzLngoeCkgfSlcblxuICAgICAgcmV0dXJuIHRoaXMuYXR0cigneCcsIHgpXG4gICAgfVxuICAgIC8vIE1vdmUgb3ZlciB5LWF4aXNcbiAgLCB5OiBmdW5jdGlvbih5KSB7XG4gICAgICB2YXIgb3kgPSB0aGlzLmF0dHIoJ3knKVxuICAgICAgICAsIG8gID0gdHlwZW9mIG95ID09PSAnbnVtYmVyJyA/IG95IC0gdGhpcy5iYm94KCkueSA6IDBcblxuICAgICAgLy8gYWN0IGFzIGdldHRlclxuICAgICAgaWYgKHkgPT0gbnVsbClcbiAgICAgICAgcmV0dXJuIHR5cGVvZiBveSA9PT0gJ251bWJlcicgPyBveSAtIG8gOiBveVxuXG4gICAgICByZXR1cm4gdGhpcy5hdHRyKCd5JywgdHlwZW9mIHkgPT09ICdudW1iZXInID8geSArIG8gOiB5KVxuICAgIH1cbiAgICAvLyBNb3ZlIGNlbnRlciBvdmVyIHgtYXhpc1xuICAsIGN4OiBmdW5jdGlvbih4KSB7XG4gICAgICByZXR1cm4geCA9PSBudWxsID8gdGhpcy5iYm94KCkuY3ggOiB0aGlzLngoeCAtIHRoaXMuYmJveCgpLndpZHRoIC8gMilcbiAgICB9XG4gICAgLy8gTW92ZSBjZW50ZXIgb3ZlciB5LWF4aXNcbiAgLCBjeTogZnVuY3Rpb24oeSkge1xuICAgICAgcmV0dXJuIHkgPT0gbnVsbCA/IHRoaXMuYmJveCgpLmN5IDogdGhpcy55KHkgLSB0aGlzLmJib3goKS5oZWlnaHQgLyAyKVxuICAgIH1cbiAgICAvLyBTZXQgdGhlIHRleHQgY29udGVudFxuICAsIHRleHQ6IGZ1bmN0aW9uKHRleHQpIHtcbiAgICAgIC8vIGFjdCBhcyBnZXR0ZXJcbiAgICAgIGlmICh0eXBlb2YgdGV4dCA9PT0gJ3VuZGVmaW5lZCcpe1xuICAgICAgICB2YXIgdGV4dCA9ICcnXG4gICAgICAgIHZhciBjaGlsZHJlbiA9IHRoaXMubm9kZS5jaGlsZE5vZGVzXG4gICAgICAgIGZvcih2YXIgaSA9IDAsIGxlbiA9IGNoaWxkcmVuLmxlbmd0aDsgaSA8IGxlbjsgKytpKXtcblxuICAgICAgICAgIC8vIGFkZCBuZXdsaW5lIGlmIGl0cyBub3QgdGhlIGZpcnN0IGNoaWxkIGFuZCBuZXdMaW5lZCBpcyBzZXQgdG8gdHJ1ZVxuICAgICAgICAgIGlmKGkgIT0gMCAmJiBjaGlsZHJlbltpXS5ub2RlVHlwZSAhPSAzICYmIFNWRy5hZG9wdChjaGlsZHJlbltpXSkuZG9tLm5ld0xpbmVkID09IHRydWUpe1xuICAgICAgICAgICAgdGV4dCArPSAnXFxuJ1xuICAgICAgICAgIH1cblxuICAgICAgICAgIC8vIGFkZCBjb250ZW50IG9mIHRoaXMgbm9kZVxuICAgICAgICAgIHRleHQgKz0gY2hpbGRyZW5baV0udGV4dENvbnRlbnRcbiAgICAgICAgfVxuXG4gICAgICAgIHJldHVybiB0ZXh0XG4gICAgICB9XG5cbiAgICAgIC8vIHJlbW92ZSBleGlzdGluZyBjb250ZW50XG4gICAgICB0aGlzLmNsZWFyKCkuYnVpbGQodHJ1ZSlcblxuICAgICAgaWYgKHR5cGVvZiB0ZXh0ID09PSAnZnVuY3Rpb24nKSB7XG4gICAgICAgIC8vIGNhbGwgYmxvY2tcbiAgICAgICAgdGV4dC5jYWxsKHRoaXMsIHRoaXMpXG5cbiAgICAgIH0gZWxzZSB7XG4gICAgICAgIC8vIHN0b3JlIHRleHQgYW5kIG1ha2Ugc3VyZSB0ZXh0IGlzIG5vdCBibGFua1xuICAgICAgICB0ZXh0ID0gdGV4dC5zcGxpdCgnXFxuJylcblxuICAgICAgICAvLyBidWlsZCBuZXcgbGluZXNcbiAgICAgICAgZm9yICh2YXIgaSA9IDAsIGlsID0gdGV4dC5sZW5ndGg7IGkgPCBpbDsgaSsrKVxuICAgICAgICAgIHRoaXMudHNwYW4odGV4dFtpXSkubmV3TGluZSgpXG4gICAgICB9XG5cbiAgICAgIC8vIGRpc2FibGUgYnVpbGQgbW9kZSBhbmQgcmVidWlsZCBsaW5lc1xuICAgICAgcmV0dXJuIHRoaXMuYnVpbGQoZmFsc2UpLnJlYnVpbGQoKVxuICAgIH1cbiAgICAvLyBTZXQgZm9udCBzaXplXG4gICwgc2l6ZTogZnVuY3Rpb24oc2l6ZSkge1xuICAgICAgcmV0dXJuIHRoaXMuYXR0cignZm9udC1zaXplJywgc2l6ZSkucmVidWlsZCgpXG4gICAgfVxuICAgIC8vIFNldCAvIGdldCBsZWFkaW5nXG4gICwgbGVhZGluZzogZnVuY3Rpb24odmFsdWUpIHtcbiAgICAgIC8vIGFjdCBhcyBnZXR0ZXJcbiAgICAgIGlmICh2YWx1ZSA9PSBudWxsKVxuICAgICAgICByZXR1cm4gdGhpcy5kb20ubGVhZGluZ1xuXG4gICAgICAvLyBhY3QgYXMgc2V0dGVyXG4gICAgICB0aGlzLmRvbS5sZWFkaW5nID0gbmV3IFNWRy5OdW1iZXIodmFsdWUpXG5cbiAgICAgIHJldHVybiB0aGlzLnJlYnVpbGQoKVxuICAgIH1cbiAgICAvLyBHZXQgYWxsIHRoZSBmaXJzdCBsZXZlbCBsaW5lc1xuICAsIGxpbmVzOiBmdW5jdGlvbigpIHtcbiAgICAgIHZhciBub2RlID0gKHRoaXMudGV4dFBhdGggJiYgdGhpcy50ZXh0UGF0aCgpIHx8IHRoaXMpLm5vZGVcblxuICAgICAgLy8gZmlsdGVyIHRzcGFucyBhbmQgbWFwIHRoZW0gdG8gU1ZHLmpzIGluc3RhbmNlc1xuICAgICAgdmFyIGxpbmVzID0gU1ZHLnV0aWxzLm1hcChTVkcudXRpbHMuZmlsdGVyU1ZHRWxlbWVudHMobm9kZS5jaGlsZE5vZGVzKSwgZnVuY3Rpb24oZWwpe1xuICAgICAgICByZXR1cm4gU1ZHLmFkb3B0KGVsKVxuICAgICAgfSlcblxuICAgICAgLy8gcmV0dXJuIGFuIGluc3RhbmNlIG9mIFNWRy5zZXRcbiAgICAgIHJldHVybiBuZXcgU1ZHLlNldChsaW5lcylcbiAgICB9XG4gICAgLy8gUmVidWlsZCBhcHBlYXJhbmNlIHR5cGVcbiAgLCByZWJ1aWxkOiBmdW5jdGlvbihyZWJ1aWxkKSB7XG4gICAgICAvLyBzdG9yZSBuZXcgcmVidWlsZCBmbGFnIGlmIGdpdmVuXG4gICAgICBpZiAodHlwZW9mIHJlYnVpbGQgPT0gJ2Jvb2xlYW4nKVxuICAgICAgICB0aGlzLl9yZWJ1aWxkID0gcmVidWlsZFxuXG4gICAgICAvLyBkZWZpbmUgcG9zaXRpb24gb2YgYWxsIGxpbmVzXG4gICAgICBpZiAodGhpcy5fcmVidWlsZCkge1xuICAgICAgICB2YXIgc2VsZiA9IHRoaXNcbiAgICAgICAgICAsIGJsYW5rTGluZU9mZnNldCA9IDBcbiAgICAgICAgICAsIGR5ID0gdGhpcy5kb20ubGVhZGluZyAqIG5ldyBTVkcuTnVtYmVyKHRoaXMuYXR0cignZm9udC1zaXplJykpXG5cbiAgICAgICAgdGhpcy5saW5lcygpLmVhY2goZnVuY3Rpb24oKSB7XG4gICAgICAgICAgaWYgKHRoaXMuZG9tLm5ld0xpbmVkKSB7XG4gICAgICAgICAgICBpZiAoIXRoaXMudGV4dFBhdGgpXG4gICAgICAgICAgICAgIHRoaXMuYXR0cigneCcsIHNlbGYuYXR0cigneCcpKVxuXG4gICAgICAgICAgICBpZih0aGlzLnRleHQoKSA9PSAnXFxuJykge1xuICAgICAgICAgICAgICBibGFua0xpbmVPZmZzZXQgKz0gZHlcbiAgICAgICAgICAgIH1lbHNle1xuICAgICAgICAgICAgICB0aGlzLmF0dHIoJ2R5JywgZHkgKyBibGFua0xpbmVPZmZzZXQpXG4gICAgICAgICAgICAgIGJsYW5rTGluZU9mZnNldCA9IDBcbiAgICAgICAgICAgIH1cbiAgICAgICAgICB9XG4gICAgICAgIH0pXG5cbiAgICAgICAgdGhpcy5maXJlKCdyZWJ1aWxkJylcbiAgICAgIH1cblxuICAgICAgcmV0dXJuIHRoaXNcbiAgICB9XG4gICAgLy8gRW5hYmxlIC8gZGlzYWJsZSBidWlsZCBtb2RlXG4gICwgYnVpbGQ6IGZ1bmN0aW9uKGJ1aWxkKSB7XG4gICAgICB0aGlzLl9idWlsZCA9ICEhYnVpbGRcbiAgICAgIHJldHVybiB0aGlzXG4gICAgfVxuICAgIC8vIG92ZXJ3cml0ZSBtZXRob2QgZnJvbSBwYXJlbnQgdG8gc2V0IGRhdGEgcHJvcGVybHlcbiAgLCBzZXREYXRhOiBmdW5jdGlvbihvKXtcbiAgICAgIHRoaXMuZG9tID0gb1xuICAgICAgdGhpcy5kb20ubGVhZGluZyA9IG5ldyBTVkcuTnVtYmVyKG8ubGVhZGluZyB8fCAxLjMpXG4gICAgICByZXR1cm4gdGhpc1xuICAgIH1cbiAgfVxuXG4gIC8vIEFkZCBwYXJlbnQgbWV0aG9kXG4sIGNvbnN0cnVjdDoge1xuICAgIC8vIENyZWF0ZSB0ZXh0IGVsZW1lbnRcbiAgICB0ZXh0OiBmdW5jdGlvbih0ZXh0KSB7XG4gICAgICByZXR1cm4gdGhpcy5wdXQobmV3IFNWRy5UZXh0KS50ZXh0KHRleHQpXG4gICAgfVxuICAgIC8vIENyZWF0ZSBwbGFpbiB0ZXh0IGVsZW1lbnRcbiAgLCBwbGFpbjogZnVuY3Rpb24odGV4dCkge1xuICAgICAgcmV0dXJuIHRoaXMucHV0KG5ldyBTVkcuVGV4dCkucGxhaW4odGV4dClcbiAgICB9XG4gIH1cblxufSlcblxuU1ZHLlRzcGFuID0gU1ZHLmludmVudCh7XG4gIC8vIEluaXRpYWxpemUgbm9kZVxuICBjcmVhdGU6ICd0c3BhbidcblxuICAvLyBJbmhlcml0IGZyb21cbiwgaW5oZXJpdDogU1ZHLlNoYXBlXG5cbiAgLy8gQWRkIGNsYXNzIG1ldGhvZHNcbiwgZXh0ZW5kOiB7XG4gICAgLy8gU2V0IHRleHQgY29udGVudFxuICAgIHRleHQ6IGZ1bmN0aW9uKHRleHQpIHtcbiAgICAgIGlmKHRleHQgPT0gbnVsbCkgcmV0dXJuIHRoaXMubm9kZS50ZXh0Q29udGVudCArICh0aGlzLmRvbS5uZXdMaW5lZCA/ICdcXG4nIDogJycpXG5cbiAgICAgIHR5cGVvZiB0ZXh0ID09PSAnZnVuY3Rpb24nID8gdGV4dC5jYWxsKHRoaXMsIHRoaXMpIDogdGhpcy5wbGFpbih0ZXh0KVxuXG4gICAgICByZXR1cm4gdGhpc1xuICAgIH1cbiAgICAvLyBTaG9ydGN1dCBkeFxuICAsIGR4OiBmdW5jdGlvbihkeCkge1xuICAgICAgcmV0dXJuIHRoaXMuYXR0cignZHgnLCBkeClcbiAgICB9XG4gICAgLy8gU2hvcnRjdXQgZHlcbiAgLCBkeTogZnVuY3Rpb24oZHkpIHtcbiAgICAgIHJldHVybiB0aGlzLmF0dHIoJ2R5JywgZHkpXG4gICAgfVxuICAgIC8vIENyZWF0ZSBuZXcgbGluZVxuICAsIG5ld0xpbmU6IGZ1bmN0aW9uKCkge1xuICAgICAgLy8gZmV0Y2ggdGV4dCBwYXJlbnRcbiAgICAgIHZhciB0ID0gdGhpcy5wYXJlbnQoU1ZHLlRleHQpXG5cbiAgICAgIC8vIG1hcmsgbmV3IGxpbmVcbiAgICAgIHRoaXMuZG9tLm5ld0xpbmVkID0gdHJ1ZVxuXG4gICAgICAvLyBhcHBseSBuZXcgaHnCoW5cbiAgICAgIHJldHVybiB0aGlzLmR5KHQuZG9tLmxlYWRpbmcgKiB0LmF0dHIoJ2ZvbnQtc2l6ZScpKS5hdHRyKCd4JywgdC54KCkpXG4gICAgfVxuICB9XG5cbn0pXG5cblNWRy5leHRlbmQoU1ZHLlRleHQsIFNWRy5Uc3Bhbiwge1xuICAvLyBDcmVhdGUgcGxhaW4gdGV4dCBub2RlXG4gIHBsYWluOiBmdW5jdGlvbih0ZXh0KSB7XG4gICAgLy8gY2xlYXIgaWYgYnVpbGQgbW9kZSBpcyBkaXNhYmxlZFxuICAgIGlmICh0aGlzLl9idWlsZCA9PT0gZmFsc2UpXG4gICAgICB0aGlzLmNsZWFyKClcblxuICAgIC8vIGNyZWF0ZSB0ZXh0IG5vZGVcbiAgICB0aGlzLm5vZGUuYXBwZW5kQ2hpbGQoZG9jdW1lbnQuY3JlYXRlVGV4dE5vZGUodGV4dCkpXG5cbiAgICByZXR1cm4gdGhpc1xuICB9XG4gIC8vIENyZWF0ZSBhIHRzcGFuXG4sIHRzcGFuOiBmdW5jdGlvbih0ZXh0KSB7XG4gICAgdmFyIG5vZGUgID0gKHRoaXMudGV4dFBhdGggJiYgdGhpcy50ZXh0UGF0aCgpIHx8IHRoaXMpLm5vZGVcbiAgICAgICwgdHNwYW4gPSBuZXcgU1ZHLlRzcGFuXG5cbiAgICAvLyBjbGVhciBpZiBidWlsZCBtb2RlIGlzIGRpc2FibGVkXG4gICAgaWYgKHRoaXMuX2J1aWxkID09PSBmYWxzZSlcbiAgICAgIHRoaXMuY2xlYXIoKVxuXG4gICAgLy8gYWRkIG5ldyB0c3BhblxuICAgIG5vZGUuYXBwZW5kQ2hpbGQodHNwYW4ubm9kZSlcblxuICAgIHJldHVybiB0c3Bhbi50ZXh0KHRleHQpXG4gIH1cbiAgLy8gQ2xlYXIgYWxsIGxpbmVzXG4sIGNsZWFyOiBmdW5jdGlvbigpIHtcbiAgICB2YXIgbm9kZSA9ICh0aGlzLnRleHRQYXRoICYmIHRoaXMudGV4dFBhdGgoKSB8fCB0aGlzKS5ub2RlXG5cbiAgICAvLyByZW1vdmUgZXhpc3RpbmcgY2hpbGQgbm9kZXNcbiAgICB3aGlsZSAobm9kZS5oYXNDaGlsZE5vZGVzKCkpXG4gICAgICBub2RlLnJlbW92ZUNoaWxkKG5vZGUubGFzdENoaWxkKVxuXG4gICAgcmV0dXJuIHRoaXNcbiAgfVxuICAvLyBHZXQgbGVuZ3RoIG9mIHRleHQgZWxlbWVudFxuLCBsZW5ndGg6IGZ1bmN0aW9uKCkge1xuICAgIHJldHVybiB0aGlzLm5vZGUuZ2V0Q29tcHV0ZWRUZXh0TGVuZ3RoKClcbiAgfVxufSlcblxuU1ZHLlRleHRQYXRoID0gU1ZHLmludmVudCh7XG4gIC8vIEluaXRpYWxpemUgbm9kZVxuICBjcmVhdGU6ICd0ZXh0UGF0aCdcblxuICAvLyBJbmhlcml0IGZyb21cbiwgaW5oZXJpdDogU1ZHLlBhcmVudFxuXG4gIC8vIERlZmluZSBwYXJlbnQgY2xhc3NcbiwgcGFyZW50OiBTVkcuVGV4dFxuXG4gIC8vIEFkZCBwYXJlbnQgbWV0aG9kXG4sIGNvbnN0cnVjdDoge1xuICAgIC8vIENyZWF0ZSBwYXRoIGZvciB0ZXh0IHRvIHJ1biBvblxuICAgIHBhdGg6IGZ1bmN0aW9uKGQpIHtcbiAgICAgIC8vIGNyZWF0ZSB0ZXh0UGF0aCBlbGVtZW50XG4gICAgICB2YXIgcGF0aCAgPSBuZXcgU1ZHLlRleHRQYXRoXG4gICAgICAgICwgdHJhY2sgPSB0aGlzLmRvYygpLmRlZnMoKS5wYXRoKGQpXG5cbiAgICAgIC8vIG1vdmUgbGluZXMgdG8gdGV4dHBhdGhcbiAgICAgIHdoaWxlICh0aGlzLm5vZGUuaGFzQ2hpbGROb2RlcygpKVxuICAgICAgICBwYXRoLm5vZGUuYXBwZW5kQ2hpbGQodGhpcy5ub2RlLmZpcnN0Q2hpbGQpXG5cbiAgICAgIC8vIGFkZCB0ZXh0UGF0aCBlbGVtZW50IGFzIGNoaWxkIG5vZGVcbiAgICAgIHRoaXMubm9kZS5hcHBlbmRDaGlsZChwYXRoLm5vZGUpXG5cbiAgICAgIC8vIGxpbmsgdGV4dFBhdGggdG8gcGF0aCBhbmQgYWRkIGNvbnRlbnRcbiAgICAgIHBhdGguYXR0cignaHJlZicsICcjJyArIHRyYWNrLCBTVkcueGxpbmspXG5cbiAgICAgIHJldHVybiB0aGlzXG4gICAgfVxuICAgIC8vIFBsb3QgcGF0aCBpZiBhbnlcbiAgLCBwbG90OiBmdW5jdGlvbihkKSB7XG4gICAgICB2YXIgdHJhY2sgPSB0aGlzLnRyYWNrKClcblxuICAgICAgaWYgKHRyYWNrKVxuICAgICAgICB0cmFjay5wbG90KGQpXG5cbiAgICAgIHJldHVybiB0aGlzXG4gICAgfVxuICAgIC8vIEdldCB0aGUgcGF0aCB0cmFjayBlbGVtZW50XG4gICwgdHJhY2s6IGZ1bmN0aW9uKCkge1xuICAgICAgdmFyIHBhdGggPSB0aGlzLnRleHRQYXRoKClcblxuICAgICAgaWYgKHBhdGgpXG4gICAgICAgIHJldHVybiBwYXRoLnJlZmVyZW5jZSgnaHJlZicpXG4gICAgfVxuICAgIC8vIEdldCB0aGUgdGV4dFBhdGggY2hpbGRcbiAgLCB0ZXh0UGF0aDogZnVuY3Rpb24oKSB7XG4gICAgICBpZiAodGhpcy5ub2RlLmZpcnN0Q2hpbGQgJiYgdGhpcy5ub2RlLmZpcnN0Q2hpbGQubm9kZU5hbWUgPT0gJ3RleHRQYXRoJylcbiAgICAgICAgcmV0dXJuIFNWRy5hZG9wdCh0aGlzLm5vZGUuZmlyc3RDaGlsZClcbiAgICB9XG4gIH1cbn0pXG5TVkcuTmVzdGVkID0gU1ZHLmludmVudCh7XG4gIC8vIEluaXRpYWxpemUgbm9kZVxuICBjcmVhdGU6IGZ1bmN0aW9uKCkge1xuICAgIHRoaXMuY29uc3RydWN0b3IuY2FsbCh0aGlzLCBTVkcuY3JlYXRlKCdzdmcnKSlcblxuICAgIHRoaXMuc3R5bGUoJ292ZXJmbG93JywgJ3Zpc2libGUnKVxuICB9XG5cbiAgLy8gSW5oZXJpdCBmcm9tXG4sIGluaGVyaXQ6IFNWRy5Db250YWluZXJcblxuICAvLyBBZGQgcGFyZW50IG1ldGhvZFxuLCBjb25zdHJ1Y3Q6IHtcbiAgICAvLyBDcmVhdGUgbmVzdGVkIHN2ZyBkb2N1bWVudFxuICAgIG5lc3RlZDogZnVuY3Rpb24oKSB7XG4gICAgICByZXR1cm4gdGhpcy5wdXQobmV3IFNWRy5OZXN0ZWQpXG4gICAgfVxuICB9XG59KVxuU1ZHLkEgPSBTVkcuaW52ZW50KHtcbiAgLy8gSW5pdGlhbGl6ZSBub2RlXG4gIGNyZWF0ZTogJ2EnXG5cbiAgLy8gSW5oZXJpdCBmcm9tXG4sIGluaGVyaXQ6IFNWRy5Db250YWluZXJcblxuICAvLyBBZGQgY2xhc3MgbWV0aG9kc1xuLCBleHRlbmQ6IHtcbiAgICAvLyBMaW5rIHVybFxuICAgIHRvOiBmdW5jdGlvbih1cmwpIHtcbiAgICAgIHJldHVybiB0aGlzLmF0dHIoJ2hyZWYnLCB1cmwsIFNWRy54bGluaylcbiAgICB9XG4gICAgLy8gTGluayBzaG93IGF0dHJpYnV0ZVxuICAsIHNob3c6IGZ1bmN0aW9uKHRhcmdldCkge1xuICAgICAgcmV0dXJuIHRoaXMuYXR0cignc2hvdycsIHRhcmdldCwgU1ZHLnhsaW5rKVxuICAgIH1cbiAgICAvLyBMaW5rIHRhcmdldCBhdHRyaWJ1dGVcbiAgLCB0YXJnZXQ6IGZ1bmN0aW9uKHRhcmdldCkge1xuICAgICAgcmV0dXJuIHRoaXMuYXR0cigndGFyZ2V0JywgdGFyZ2V0KVxuICAgIH1cbiAgfVxuXG4gIC8vIEFkZCBwYXJlbnQgbWV0aG9kXG4sIGNvbnN0cnVjdDoge1xuICAgIC8vIENyZWF0ZSBhIGh5cGVybGluayBlbGVtZW50XG4gICAgbGluazogZnVuY3Rpb24odXJsKSB7XG4gICAgICByZXR1cm4gdGhpcy5wdXQobmV3IFNWRy5BKS50byh1cmwpXG4gICAgfVxuICB9XG59KVxuXG5TVkcuZXh0ZW5kKFNWRy5FbGVtZW50LCB7XG4gIC8vIENyZWF0ZSBhIGh5cGVybGluayBlbGVtZW50XG4gIGxpbmtUbzogZnVuY3Rpb24odXJsKSB7XG4gICAgdmFyIGxpbmsgPSBuZXcgU1ZHLkFcblxuICAgIGlmICh0eXBlb2YgdXJsID09ICdmdW5jdGlvbicpXG4gICAgICB1cmwuY2FsbChsaW5rLCBsaW5rKVxuICAgIGVsc2VcbiAgICAgIGxpbmsudG8odXJsKVxuXG4gICAgcmV0dXJuIHRoaXMucGFyZW50KCkucHV0KGxpbmspLnB1dCh0aGlzKVxuICB9XG5cbn0pXG5TVkcuTWFya2VyID0gU1ZHLmludmVudCh7XG4gIC8vIEluaXRpYWxpemUgbm9kZVxuICBjcmVhdGU6ICdtYXJrZXInXG5cbiAgLy8gSW5oZXJpdCBmcm9tXG4sIGluaGVyaXQ6IFNWRy5Db250YWluZXJcblxuICAvLyBBZGQgY2xhc3MgbWV0aG9kc1xuLCBleHRlbmQ6IHtcbiAgICAvLyBTZXQgd2lkdGggb2YgZWxlbWVudFxuICAgIHdpZHRoOiBmdW5jdGlvbih3aWR0aCkge1xuICAgICAgcmV0dXJuIHRoaXMuYXR0cignbWFya2VyV2lkdGgnLCB3aWR0aClcbiAgICB9XG4gICAgLy8gU2V0IGhlaWdodCBvZiBlbGVtZW50XG4gICwgaGVpZ2h0OiBmdW5jdGlvbihoZWlnaHQpIHtcbiAgICAgIHJldHVybiB0aGlzLmF0dHIoJ21hcmtlckhlaWdodCcsIGhlaWdodClcbiAgICB9XG4gICAgLy8gU2V0IG1hcmtlciByZWZYIGFuZCByZWZZXG4gICwgcmVmOiBmdW5jdGlvbih4LCB5KSB7XG4gICAgICByZXR1cm4gdGhpcy5hdHRyKCdyZWZYJywgeCkuYXR0cigncmVmWScsIHkpXG4gICAgfVxuICAgIC8vIFVwZGF0ZSBtYXJrZXJcbiAgLCB1cGRhdGU6IGZ1bmN0aW9uKGJsb2NrKSB7XG4gICAgICAvLyByZW1vdmUgYWxsIGNvbnRlbnRcbiAgICAgIHRoaXMuY2xlYXIoKVxuXG4gICAgICAvLyBpbnZva2UgcGFzc2VkIGJsb2NrXG4gICAgICBpZiAodHlwZW9mIGJsb2NrID09ICdmdW5jdGlvbicpXG4gICAgICAgIGJsb2NrLmNhbGwodGhpcywgdGhpcylcblxuICAgICAgcmV0dXJuIHRoaXNcbiAgICB9XG4gICAgLy8gUmV0dXJuIHRoZSBmaWxsIGlkXG4gICwgdG9TdHJpbmc6IGZ1bmN0aW9uKCkge1xuICAgICAgcmV0dXJuICd1cmwoIycgKyB0aGlzLmlkKCkgKyAnKSdcbiAgICB9XG4gIH1cblxuICAvLyBBZGQgcGFyZW50IG1ldGhvZFxuLCBjb25zdHJ1Y3Q6IHtcbiAgICBtYXJrZXI6IGZ1bmN0aW9uKHdpZHRoLCBoZWlnaHQsIGJsb2NrKSB7XG4gICAgICAvLyBDcmVhdGUgbWFya2VyIGVsZW1lbnQgaW4gZGVmc1xuICAgICAgcmV0dXJuIHRoaXMuZGVmcygpLm1hcmtlcih3aWR0aCwgaGVpZ2h0LCBibG9jaylcbiAgICB9XG4gIH1cblxufSlcblxuU1ZHLmV4dGVuZChTVkcuRGVmcywge1xuICAvLyBDcmVhdGUgbWFya2VyXG4gIG1hcmtlcjogZnVuY3Rpb24od2lkdGgsIGhlaWdodCwgYmxvY2spIHtcbiAgICAvLyBTZXQgZGVmYXVsdCB2aWV3Ym94IHRvIG1hdGNoIHRoZSB3aWR0aCBhbmQgaGVpZ2h0LCBzZXQgcmVmIHRvIGN4IGFuZCBjeSBhbmQgc2V0IG9yaWVudCB0byBhdXRvXG4gICAgcmV0dXJuIHRoaXMucHV0KG5ldyBTVkcuTWFya2VyKVxuICAgICAgLnNpemUod2lkdGgsIGhlaWdodClcbiAgICAgIC5yZWYod2lkdGggLyAyLCBoZWlnaHQgLyAyKVxuICAgICAgLnZpZXdib3goMCwgMCwgd2lkdGgsIGhlaWdodClcbiAgICAgIC5hdHRyKCdvcmllbnQnLCAnYXV0bycpXG4gICAgICAudXBkYXRlKGJsb2NrKVxuICB9XG5cbn0pXG5cblNWRy5leHRlbmQoU1ZHLkxpbmUsIFNWRy5Qb2x5bGluZSwgU1ZHLlBvbHlnb24sIFNWRy5QYXRoLCB7XG4gIC8vIENyZWF0ZSBhbmQgYXR0YWNoIG1hcmtlcnNcbiAgbWFya2VyOiBmdW5jdGlvbihtYXJrZXIsIHdpZHRoLCBoZWlnaHQsIGJsb2NrKSB7XG4gICAgdmFyIGF0dHIgPSBbJ21hcmtlciddXG5cbiAgICAvLyBCdWlsZCBhdHRyaWJ1dGUgbmFtZVxuICAgIGlmIChtYXJrZXIgIT0gJ2FsbCcpIGF0dHIucHVzaChtYXJrZXIpXG4gICAgYXR0ciA9IGF0dHIuam9pbignLScpXG5cbiAgICAvLyBTZXQgbWFya2VyIGF0dHJpYnV0ZVxuICAgIG1hcmtlciA9IGFyZ3VtZW50c1sxXSBpbnN0YW5jZW9mIFNWRy5NYXJrZXIgP1xuICAgICAgYXJndW1lbnRzWzFdIDpcbiAgICAgIHRoaXMuZG9jKCkubWFya2VyKHdpZHRoLCBoZWlnaHQsIGJsb2NrKVxuXG4gICAgcmV0dXJuIHRoaXMuYXR0cihhdHRyLCBtYXJrZXIpXG4gIH1cblxufSlcbi8vIERlZmluZSBsaXN0IG9mIGF2YWlsYWJsZSBhdHRyaWJ1dGVzIGZvciBzdHJva2UgYW5kIGZpbGxcbnZhciBzdWdhciA9IHtcbiAgc3Ryb2tlOiBbJ2NvbG9yJywgJ3dpZHRoJywgJ29wYWNpdHknLCAnbGluZWNhcCcsICdsaW5lam9pbicsICdtaXRlcmxpbWl0JywgJ2Rhc2hhcnJheScsICdkYXNob2Zmc2V0J11cbiwgZmlsbDogICBbJ2NvbG9yJywgJ29wYWNpdHknLCAncnVsZSddXG4sIHByZWZpeDogZnVuY3Rpb24odCwgYSkge1xuICAgIHJldHVybiBhID09ICdjb2xvcicgPyB0IDogdCArICctJyArIGFcbiAgfVxufVxuXG4vLyBBZGQgc3VnYXIgZm9yIGZpbGwgYW5kIHN0cm9rZVxuO1snZmlsbCcsICdzdHJva2UnXS5mb3JFYWNoKGZ1bmN0aW9uKG0pIHtcbiAgdmFyIGksIGV4dGVuc2lvbiA9IHt9XG5cbiAgZXh0ZW5zaW9uW21dID0gZnVuY3Rpb24obykge1xuICAgIGlmICh0eXBlb2YgbyA9PSAndW5kZWZpbmVkJylcbiAgICAgIHJldHVybiB0aGlzXG4gICAgaWYgKHR5cGVvZiBvID09ICdzdHJpbmcnIHx8IFNWRy5Db2xvci5pc1JnYihvKSB8fCAobyAmJiB0eXBlb2Ygby5maWxsID09PSAnZnVuY3Rpb24nKSlcbiAgICAgIHRoaXMuYXR0cihtLCBvKVxuXG4gICAgZWxzZVxuICAgICAgLy8gc2V0IGFsbCBhdHRyaWJ1dGVzIGZyb20gc3VnYXIuZmlsbCBhbmQgc3VnYXIuc3Ryb2tlIGxpc3RcbiAgICAgIGZvciAoaSA9IHN1Z2FyW21dLmxlbmd0aCAtIDE7IGkgPj0gMDsgaS0tKVxuICAgICAgICBpZiAob1tzdWdhclttXVtpXV0gIT0gbnVsbClcbiAgICAgICAgICB0aGlzLmF0dHIoc3VnYXIucHJlZml4KG0sIHN1Z2FyW21dW2ldKSwgb1tzdWdhclttXVtpXV0pXG5cbiAgICByZXR1cm4gdGhpc1xuICB9XG5cbiAgU1ZHLmV4dGVuZChTVkcuRWxlbWVudCwgU1ZHLkZYLCBleHRlbnNpb24pXG5cbn0pXG5cblNWRy5leHRlbmQoU1ZHLkVsZW1lbnQsIFNWRy5GWCwge1xuICAvLyBNYXAgcm90YXRpb24gdG8gdHJhbnNmb3JtXG4gIHJvdGF0ZTogZnVuY3Rpb24oZCwgY3gsIGN5KSB7XG4gICAgcmV0dXJuIHRoaXMudHJhbnNmb3JtKHsgcm90YXRpb246IGQsIGN4OiBjeCwgY3k6IGN5IH0pXG4gIH1cbiAgLy8gTWFwIHNrZXcgdG8gdHJhbnNmb3JtXG4sIHNrZXc6IGZ1bmN0aW9uKHgsIHksIGN4LCBjeSkge1xuICAgIHJldHVybiBhcmd1bWVudHMubGVuZ3RoID09IDEgIHx8IGFyZ3VtZW50cy5sZW5ndGggPT0gMyA/XG4gICAgICB0aGlzLnRyYW5zZm9ybSh7IHNrZXc6IHgsIGN4OiB5LCBjeTogY3ggfSkgOlxuICAgICAgdGhpcy50cmFuc2Zvcm0oeyBza2V3WDogeCwgc2tld1k6IHksIGN4OiBjeCwgY3k6IGN5IH0pXG4gIH1cbiAgLy8gTWFwIHNjYWxlIHRvIHRyYW5zZm9ybVxuLCBzY2FsZTogZnVuY3Rpb24oeCwgeSwgY3gsIGN5KSB7XG4gICAgcmV0dXJuIGFyZ3VtZW50cy5sZW5ndGggPT0gMSAgfHwgYXJndW1lbnRzLmxlbmd0aCA9PSAzID9cbiAgICAgIHRoaXMudHJhbnNmb3JtKHsgc2NhbGU6IHgsIGN4OiB5LCBjeTogY3ggfSkgOlxuICAgICAgdGhpcy50cmFuc2Zvcm0oeyBzY2FsZVg6IHgsIHNjYWxlWTogeSwgY3g6IGN4LCBjeTogY3kgfSlcbiAgfVxuICAvLyBNYXAgdHJhbnNsYXRlIHRvIHRyYW5zZm9ybVxuLCB0cmFuc2xhdGU6IGZ1bmN0aW9uKHgsIHkpIHtcbiAgICByZXR1cm4gdGhpcy50cmFuc2Zvcm0oeyB4OiB4LCB5OiB5IH0pXG4gIH1cbiAgLy8gTWFwIGZsaXAgdG8gdHJhbnNmb3JtXG4sIGZsaXA6IGZ1bmN0aW9uKGEsIG8pIHtcbiAgICByZXR1cm4gdGhpcy50cmFuc2Zvcm0oeyBmbGlwOiBhLCBvZmZzZXQ6IG8gfSlcbiAgfVxuICAvLyBNYXAgbWF0cml4IHRvIHRyYW5zZm9ybVxuLCBtYXRyaXg6IGZ1bmN0aW9uKG0pIHtcbiAgICByZXR1cm4gdGhpcy5hdHRyKCd0cmFuc2Zvcm0nLCBuZXcgU1ZHLk1hdHJpeChtKSlcbiAgfVxuICAvLyBPcGFjaXR5XG4sIG9wYWNpdHk6IGZ1bmN0aW9uKHZhbHVlKSB7XG4gICAgcmV0dXJuIHRoaXMuYXR0cignb3BhY2l0eScsIHZhbHVlKVxuICB9XG4gIC8vIFJlbGF0aXZlIG1vdmUgb3ZlciB4IGF4aXNcbiwgZHg6IGZ1bmN0aW9uKHgpIHtcbiAgICByZXR1cm4gdGhpcy54KCh0aGlzIGluc3RhbmNlb2YgU1ZHLkZYID8gMCA6IHRoaXMueCgpKSArIHgsIHRydWUpXG4gIH1cbiAgLy8gUmVsYXRpdmUgbW92ZSBvdmVyIHkgYXhpc1xuLCBkeTogZnVuY3Rpb24oeSkge1xuICAgIHJldHVybiB0aGlzLnkoKHRoaXMgaW5zdGFuY2VvZiBTVkcuRlggPyAwIDogdGhpcy55KCkpICsgeSwgdHJ1ZSlcbiAgfVxuICAvLyBSZWxhdGl2ZSBtb3ZlIG92ZXIgeCBhbmQgeSBheGVzXG4sIGRtb3ZlOiBmdW5jdGlvbih4LCB5KSB7XG4gICAgcmV0dXJuIHRoaXMuZHgoeCkuZHkoeSlcbiAgfVxufSlcblxuU1ZHLmV4dGVuZChTVkcuUmVjdCwgU1ZHLkVsbGlwc2UsIFNWRy5DaXJjbGUsIFNWRy5HcmFkaWVudCwgU1ZHLkZYLCB7XG4gIC8vIEFkZCB4IGFuZCB5IHJhZGl1c1xuICByYWRpdXM6IGZ1bmN0aW9uKHgsIHkpIHtcbiAgICB2YXIgdHlwZSA9ICh0aGlzLl90YXJnZXQgfHwgdGhpcykudHlwZTtcbiAgICByZXR1cm4gdHlwZSA9PSAncmFkaWFsJyB8fCB0eXBlID09ICdjaXJjbGUnID9cbiAgICAgIHRoaXMuYXR0cigncicsIG5ldyBTVkcuTnVtYmVyKHgpKSA6XG4gICAgICB0aGlzLnJ4KHgpLnJ5KHkgPT0gbnVsbCA/IHggOiB5KVxuICB9XG59KVxuXG5TVkcuZXh0ZW5kKFNWRy5QYXRoLCB7XG4gIC8vIEdldCBwYXRoIGxlbmd0aFxuICBsZW5ndGg6IGZ1bmN0aW9uKCkge1xuICAgIHJldHVybiB0aGlzLm5vZGUuZ2V0VG90YWxMZW5ndGgoKVxuICB9XG4gIC8vIEdldCBwb2ludCBhdCBsZW5ndGhcbiwgcG9pbnRBdDogZnVuY3Rpb24obGVuZ3RoKSB7XG4gICAgcmV0dXJuIHRoaXMubm9kZS5nZXRQb2ludEF0TGVuZ3RoKGxlbmd0aClcbiAgfVxufSlcblxuU1ZHLmV4dGVuZChTVkcuUGFyZW50LCBTVkcuVGV4dCwgU1ZHLkZYLCB7XG4gIC8vIFNldCBmb250XG4gIGZvbnQ6IGZ1bmN0aW9uKG8pIHtcbiAgICBmb3IgKHZhciBrIGluIG8pXG4gICAgICBrID09ICdsZWFkaW5nJyA/XG4gICAgICAgIHRoaXMubGVhZGluZyhvW2tdKSA6XG4gICAgICBrID09ICdhbmNob3InID9cbiAgICAgICAgdGhpcy5hdHRyKCd0ZXh0LWFuY2hvcicsIG9ba10pIDpcbiAgICAgIGsgPT0gJ3NpemUnIHx8IGsgPT0gJ2ZhbWlseScgfHwgayA9PSAnd2VpZ2h0JyB8fCBrID09ICdzdHJldGNoJyB8fCBrID09ICd2YXJpYW50JyB8fCBrID09ICdzdHlsZScgP1xuICAgICAgICB0aGlzLmF0dHIoJ2ZvbnQtJysgaywgb1trXSkgOlxuICAgICAgICB0aGlzLmF0dHIoaywgb1trXSlcblxuICAgIHJldHVybiB0aGlzXG4gIH1cbn0pXG5cblNWRy5TZXQgPSBTVkcuaW52ZW50KHtcbiAgLy8gSW5pdGlhbGl6ZVxuICBjcmVhdGU6IGZ1bmN0aW9uKG1lbWJlcnMpIHtcbiAgICAvLyBTZXQgaW5pdGlhbCBzdGF0ZVxuICAgIEFycmF5LmlzQXJyYXkobWVtYmVycykgPyB0aGlzLm1lbWJlcnMgPSBtZW1iZXJzIDogdGhpcy5jbGVhcigpXG4gIH1cblxuICAvLyBBZGQgY2xhc3MgbWV0aG9kc1xuLCBleHRlbmQ6IHtcbiAgICAvLyBBZGQgZWxlbWVudCB0byBzZXRcbiAgICBhZGQ6IGZ1bmN0aW9uKCkge1xuICAgICAgdmFyIGksIGlsLCBlbGVtZW50cyA9IFtdLnNsaWNlLmNhbGwoYXJndW1lbnRzKVxuXG4gICAgICBmb3IgKGkgPSAwLCBpbCA9IGVsZW1lbnRzLmxlbmd0aDsgaSA8IGlsOyBpKyspXG4gICAgICAgIHRoaXMubWVtYmVycy5wdXNoKGVsZW1lbnRzW2ldKVxuXG4gICAgICByZXR1cm4gdGhpc1xuICAgIH1cbiAgICAvLyBSZW1vdmUgZWxlbWVudCBmcm9tIHNldFxuICAsIHJlbW92ZTogZnVuY3Rpb24oZWxlbWVudCkge1xuICAgICAgdmFyIGkgPSB0aGlzLmluZGV4KGVsZW1lbnQpXG5cbiAgICAgIC8vIHJlbW92ZSBnaXZlbiBjaGlsZFxuICAgICAgaWYgKGkgPiAtMSlcbiAgICAgICAgdGhpcy5tZW1iZXJzLnNwbGljZShpLCAxKVxuXG4gICAgICByZXR1cm4gdGhpc1xuICAgIH1cbiAgICAvLyBJdGVyYXRlIG92ZXIgYWxsIG1lbWJlcnNcbiAgLCBlYWNoOiBmdW5jdGlvbihibG9jaykge1xuICAgICAgZm9yICh2YXIgaSA9IDAsIGlsID0gdGhpcy5tZW1iZXJzLmxlbmd0aDsgaSA8IGlsOyBpKyspXG4gICAgICAgIGJsb2NrLmFwcGx5KHRoaXMubWVtYmVyc1tpXSwgW2ksIHRoaXMubWVtYmVyc10pXG5cbiAgICAgIHJldHVybiB0aGlzXG4gICAgfVxuICAgIC8vIFJlc3RvcmUgdG8gZGVmYXVsdHNcbiAgLCBjbGVhcjogZnVuY3Rpb24oKSB7XG4gICAgICAvLyBpbml0aWFsaXplIHN0b3JlXG4gICAgICB0aGlzLm1lbWJlcnMgPSBbXVxuXG4gICAgICByZXR1cm4gdGhpc1xuICAgIH1cbiAgICAvLyBHZXQgdGhlIGxlbmd0aCBvZiBhIHNldFxuICAsIGxlbmd0aDogZnVuY3Rpb24oKSB7XG4gICAgICByZXR1cm4gdGhpcy5tZW1iZXJzLmxlbmd0aFxuICAgIH1cbiAgICAvLyBDaGVja3MgaWYgYSBnaXZlbiBlbGVtZW50IGlzIHByZXNlbnQgaW4gc2V0XG4gICwgaGFzOiBmdW5jdGlvbihlbGVtZW50KSB7XG4gICAgICByZXR1cm4gdGhpcy5pbmRleChlbGVtZW50KSA+PSAwXG4gICAgfVxuICAgIC8vIHJldHVucyBpbmRleCBvZiBnaXZlbiBlbGVtZW50IGluIHNldFxuICAsIGluZGV4OiBmdW5jdGlvbihlbGVtZW50KSB7XG4gICAgICByZXR1cm4gdGhpcy5tZW1iZXJzLmluZGV4T2YoZWxlbWVudClcbiAgICB9XG4gICAgLy8gR2V0IG1lbWJlciBhdCBnaXZlbiBpbmRleFxuICAsIGdldDogZnVuY3Rpb24oaSkge1xuICAgICAgcmV0dXJuIHRoaXMubWVtYmVyc1tpXVxuICAgIH1cbiAgICAvLyBHZXQgZmlyc3QgbWVtYmVyXG4gICwgZmlyc3Q6IGZ1bmN0aW9uKCkge1xuICAgICAgcmV0dXJuIHRoaXMuZ2V0KDApXG4gICAgfVxuICAgIC8vIEdldCBsYXN0IG1lbWJlclxuICAsIGxhc3Q6IGZ1bmN0aW9uKCkge1xuICAgICAgcmV0dXJuIHRoaXMuZ2V0KHRoaXMubWVtYmVycy5sZW5ndGggLSAxKVxuICAgIH1cbiAgICAvLyBEZWZhdWx0IHZhbHVlXG4gICwgdmFsdWVPZjogZnVuY3Rpb24oKSB7XG4gICAgICByZXR1cm4gdGhpcy5tZW1iZXJzXG4gICAgfVxuICAgIC8vIEdldCB0aGUgYm91bmRpbmcgYm94IG9mIGFsbCBtZW1iZXJzIGluY2x1ZGVkIG9yIGVtcHR5IGJveCBpZiBzZXQgaGFzIG5vIGl0ZW1zXG4gICwgYmJveDogZnVuY3Rpb24oKXtcbiAgICAgIHZhciBib3ggPSBuZXcgU1ZHLkJCb3goKVxuXG4gICAgICAvLyByZXR1cm4gYW4gZW1wdHkgYm94IG9mIHRoZXJlIGFyZSBubyBtZW1iZXJzXG4gICAgICBpZiAodGhpcy5tZW1iZXJzLmxlbmd0aCA9PSAwKVxuICAgICAgICByZXR1cm4gYm94XG5cbiAgICAgIC8vIGdldCB0aGUgZmlyc3QgcmJveCBhbmQgdXBkYXRlIHRoZSB0YXJnZXQgYmJveFxuICAgICAgdmFyIHJib3ggPSB0aGlzLm1lbWJlcnNbMF0ucmJveCgpXG4gICAgICBib3gueCAgICAgID0gcmJveC54XG4gICAgICBib3gueSAgICAgID0gcmJveC55XG4gICAgICBib3gud2lkdGggID0gcmJveC53aWR0aFxuICAgICAgYm94LmhlaWdodCA9IHJib3guaGVpZ2h0XG5cbiAgICAgIHRoaXMuZWFjaChmdW5jdGlvbigpIHtcbiAgICAgICAgLy8gdXNlciByYm94IGZvciBjb3JyZWN0IHBvc2l0aW9uIGFuZCB2aXN1YWwgcmVwcmVzZW50YXRpb25cbiAgICAgICAgYm94ID0gYm94Lm1lcmdlKHRoaXMucmJveCgpKVxuICAgICAgfSlcblxuICAgICAgcmV0dXJuIGJveFxuICAgIH1cbiAgfVxuXG4gIC8vIEFkZCBwYXJlbnQgbWV0aG9kXG4sIGNvbnN0cnVjdDoge1xuICAgIC8vIENyZWF0ZSBhIG5ldyBzZXRcbiAgICBzZXQ6IGZ1bmN0aW9uKG1lbWJlcnMpIHtcbiAgICAgIHJldHVybiBuZXcgU1ZHLlNldChtZW1iZXJzKVxuICAgIH1cbiAgfVxufSlcblxuU1ZHLkZYLlNldCA9IFNWRy5pbnZlbnQoe1xuICAvLyBJbml0aWFsaXplIG5vZGVcbiAgY3JlYXRlOiBmdW5jdGlvbihzZXQpIHtcbiAgICAvLyBzdG9yZSByZWZlcmVuY2UgdG8gc2V0XG4gICAgdGhpcy5zZXQgPSBzZXRcbiAgfVxuXG59KVxuXG4vLyBBbGlhcyBtZXRob2RzXG5TVkcuU2V0LmluaGVyaXQgPSBmdW5jdGlvbigpIHtcbiAgdmFyIG1cbiAgICAsIG1ldGhvZHMgPSBbXVxuXG4gIC8vIGdhdGhlciBzaGFwZSBtZXRob2RzXG4gIGZvcih2YXIgbSBpbiBTVkcuU2hhcGUucHJvdG90eXBlKVxuICAgIGlmICh0eXBlb2YgU1ZHLlNoYXBlLnByb3RvdHlwZVttXSA9PSAnZnVuY3Rpb24nICYmIHR5cGVvZiBTVkcuU2V0LnByb3RvdHlwZVttXSAhPSAnZnVuY3Rpb24nKVxuICAgICAgbWV0aG9kcy5wdXNoKG0pXG5cbiAgLy8gYXBwbHkgc2hhcGUgYWxpYXNzZXNcbiAgbWV0aG9kcy5mb3JFYWNoKGZ1bmN0aW9uKG1ldGhvZCkge1xuICAgIFNWRy5TZXQucHJvdG90eXBlW21ldGhvZF0gPSBmdW5jdGlvbigpIHtcbiAgICAgIGZvciAodmFyIGkgPSAwLCBpbCA9IHRoaXMubWVtYmVycy5sZW5ndGg7IGkgPCBpbDsgaSsrKVxuICAgICAgICBpZiAodGhpcy5tZW1iZXJzW2ldICYmIHR5cGVvZiB0aGlzLm1lbWJlcnNbaV1bbWV0aG9kXSA9PSAnZnVuY3Rpb24nKVxuICAgICAgICAgIHRoaXMubWVtYmVyc1tpXVttZXRob2RdLmFwcGx5KHRoaXMubWVtYmVyc1tpXSwgYXJndW1lbnRzKVxuXG4gICAgICByZXR1cm4gbWV0aG9kID09ICdhbmltYXRlJyA/ICh0aGlzLmZ4IHx8ICh0aGlzLmZ4ID0gbmV3IFNWRy5GWC5TZXQodGhpcykpKSA6IHRoaXNcbiAgICB9XG4gIH0pXG5cbiAgLy8gY2xlYXIgbWV0aG9kcyBmb3IgdGhlIG5leHQgcm91bmRcbiAgbWV0aG9kcyA9IFtdXG5cbiAgLy8gZ2F0aGVyIGZ4IG1ldGhvZHNcbiAgZm9yKHZhciBtIGluIFNWRy5GWC5wcm90b3R5cGUpXG4gICAgaWYgKHR5cGVvZiBTVkcuRlgucHJvdG90eXBlW21dID09ICdmdW5jdGlvbicgJiYgdHlwZW9mIFNWRy5GWC5TZXQucHJvdG90eXBlW21dICE9ICdmdW5jdGlvbicpXG4gICAgICBtZXRob2RzLnB1c2gobSlcblxuICAvLyBhcHBseSBmeCBhbGlhc3Nlc1xuICBtZXRob2RzLmZvckVhY2goZnVuY3Rpb24obWV0aG9kKSB7XG4gICAgU1ZHLkZYLlNldC5wcm90b3R5cGVbbWV0aG9kXSA9IGZ1bmN0aW9uKCkge1xuICAgICAgZm9yICh2YXIgaSA9IDAsIGlsID0gdGhpcy5zZXQubWVtYmVycy5sZW5ndGg7IGkgPCBpbDsgaSsrKVxuICAgICAgICB0aGlzLnNldC5tZW1iZXJzW2ldLmZ4W21ldGhvZF0uYXBwbHkodGhpcy5zZXQubWVtYmVyc1tpXS5meCwgYXJndW1lbnRzKVxuXG4gICAgICByZXR1cm4gdGhpc1xuICAgIH1cbiAgfSlcbn1cblxuXG5cblxuU1ZHLmV4dGVuZChTVkcuRWxlbWVudCwge1xuICAvLyBTdG9yZSBkYXRhIHZhbHVlcyBvbiBzdmcgbm9kZXNcbiAgZGF0YTogZnVuY3Rpb24oYSwgdiwgcikge1xuICAgIGlmICh0eXBlb2YgYSA9PSAnb2JqZWN0Jykge1xuICAgICAgZm9yICh2IGluIGEpXG4gICAgICAgIHRoaXMuZGF0YSh2LCBhW3ZdKVxuXG4gICAgfSBlbHNlIGlmIChhcmd1bWVudHMubGVuZ3RoIDwgMikge1xuICAgICAgdHJ5IHtcbiAgICAgICAgcmV0dXJuIEpTT04ucGFyc2UodGhpcy5hdHRyKCdkYXRhLScgKyBhKSlcbiAgICAgIH0gY2F0Y2goZSkge1xuICAgICAgICByZXR1cm4gdGhpcy5hdHRyKCdkYXRhLScgKyBhKVxuICAgICAgfVxuXG4gICAgfSBlbHNlIHtcbiAgICAgIHRoaXMuYXR0cihcbiAgICAgICAgJ2RhdGEtJyArIGFcbiAgICAgICwgdiA9PT0gbnVsbCA/XG4gICAgICAgICAgbnVsbCA6XG4gICAgICAgIHIgPT09IHRydWUgfHwgdHlwZW9mIHYgPT09ICdzdHJpbmcnIHx8IHR5cGVvZiB2ID09PSAnbnVtYmVyJyA/XG4gICAgICAgICAgdiA6XG4gICAgICAgICAgSlNPTi5zdHJpbmdpZnkodilcbiAgICAgIClcbiAgICB9XG5cbiAgICByZXR1cm4gdGhpc1xuICB9XG59KVxuU1ZHLmV4dGVuZChTVkcuRWxlbWVudCwge1xuICAvLyBSZW1lbWJlciBhcmJpdHJhcnkgZGF0YVxuICByZW1lbWJlcjogZnVuY3Rpb24oaywgdikge1xuICAgIC8vIHJlbWVtYmVyIGV2ZXJ5IGl0ZW0gaW4gYW4gb2JqZWN0IGluZGl2aWR1YWxseVxuICAgIGlmICh0eXBlb2YgYXJndW1lbnRzWzBdID09ICdvYmplY3QnKVxuICAgICAgZm9yICh2YXIgdiBpbiBrKVxuICAgICAgICB0aGlzLnJlbWVtYmVyKHYsIGtbdl0pXG5cbiAgICAvLyByZXRyaWV2ZSBtZW1vcnlcbiAgICBlbHNlIGlmIChhcmd1bWVudHMubGVuZ3RoID09IDEpXG4gICAgICByZXR1cm4gdGhpcy5tZW1vcnkoKVtrXVxuXG4gICAgLy8gc3RvcmUgbWVtb3J5XG4gICAgZWxzZVxuICAgICAgdGhpcy5tZW1vcnkoKVtrXSA9IHZcblxuICAgIHJldHVybiB0aGlzXG4gIH1cblxuICAvLyBFcmFzZSBhIGdpdmVuIG1lbW9yeVxuLCBmb3JnZXQ6IGZ1bmN0aW9uKCkge1xuICAgIGlmIChhcmd1bWVudHMubGVuZ3RoID09IDApXG4gICAgICB0aGlzLl9tZW1vcnkgPSB7fVxuICAgIGVsc2VcbiAgICAgIGZvciAodmFyIGkgPSBhcmd1bWVudHMubGVuZ3RoIC0gMTsgaSA+PSAwOyBpLS0pXG4gICAgICAgIGRlbGV0ZSB0aGlzLm1lbW9yeSgpW2FyZ3VtZW50c1tpXV1cblxuICAgIHJldHVybiB0aGlzXG4gIH1cblxuICAvLyBJbml0aWFsaXplIG9yIHJldHVybiBsb2NhbCBtZW1vcnkgb2JqZWN0XG4sIG1lbW9yeTogZnVuY3Rpb24oKSB7XG4gICAgcmV0dXJuIHRoaXMuX21lbW9yeSB8fCAodGhpcy5fbWVtb3J5ID0ge30pXG4gIH1cblxufSlcbi8vIE1ldGhvZCBmb3IgZ2V0dGluZyBhbiBlbGVtZW50IGJ5IGlkXG5TVkcuZ2V0ID0gZnVuY3Rpb24oaWQpIHtcbiAgdmFyIG5vZGUgPSBkb2N1bWVudC5nZXRFbGVtZW50QnlJZChpZEZyb21SZWZlcmVuY2UoaWQpIHx8IGlkKVxuICByZXR1cm4gU1ZHLmFkb3B0KG5vZGUpXG59XG5cbi8vIFNlbGVjdCBlbGVtZW50cyBieSBxdWVyeSBzdHJpbmdcblNWRy5zZWxlY3QgPSBmdW5jdGlvbihxdWVyeSwgcGFyZW50KSB7XG4gIHJldHVybiBuZXcgU1ZHLlNldChcbiAgICBTVkcudXRpbHMubWFwKChwYXJlbnQgfHwgZG9jdW1lbnQpLnF1ZXJ5U2VsZWN0b3JBbGwocXVlcnkpLCBmdW5jdGlvbihub2RlKSB7XG4gICAgICByZXR1cm4gU1ZHLmFkb3B0KG5vZGUpXG4gICAgfSlcbiAgKVxufVxuXG5TVkcuZXh0ZW5kKFNWRy5QYXJlbnQsIHtcbiAgLy8gU2NvcGVkIHNlbGVjdCBtZXRob2RcbiAgc2VsZWN0OiBmdW5jdGlvbihxdWVyeSkge1xuICAgIHJldHVybiBTVkcuc2VsZWN0KHF1ZXJ5LCB0aGlzLm5vZGUpXG4gIH1cblxufSlcbmZ1bmN0aW9uIGlzKGVsLCBvYmope1xuICByZXR1cm4gZWwgaW5zdGFuY2VvZiBvYmpcbn1cblxuLy8gdGVzdHMgaWYgYSBnaXZlbiBzZWxlY3RvciBtYXRjaGVzIGFuIGVsZW1lbnRcbmZ1bmN0aW9uIG1hdGNoZXMoZWwsIHNlbGVjdG9yKSB7XG4gIHJldHVybiAoZWwubWF0Y2hlcyB8fCBlbC5tYXRjaGVzU2VsZWN0b3IgfHwgZWwubXNNYXRjaGVzU2VsZWN0b3IgfHwgZWwubW96TWF0Y2hlc1NlbGVjdG9yIHx8IGVsLndlYmtpdE1hdGNoZXNTZWxlY3RvciB8fCBlbC5vTWF0Y2hlc1NlbGVjdG9yKS5jYWxsKGVsLCBzZWxlY3Rvcik7XG59XG5cbi8vIENvbnZlcnQgZGFzaC1zZXBhcmF0ZWQtc3RyaW5nIHRvIGNhbWVsQ2FzZVxuZnVuY3Rpb24gY2FtZWxDYXNlKHMpIHtcbiAgcmV0dXJuIHMudG9Mb3dlckNhc2UoKS5yZXBsYWNlKC8tKC4pL2csIGZ1bmN0aW9uKG0sIGcpIHtcbiAgICByZXR1cm4gZy50b1VwcGVyQ2FzZSgpXG4gIH0pXG59XG5cbi8vIENhcGl0YWxpemUgZmlyc3QgbGV0dGVyIG9mIGEgc3RyaW5nXG5mdW5jdGlvbiBjYXBpdGFsaXplKHMpIHtcbiAgcmV0dXJuIHMuY2hhckF0KDApLnRvVXBwZXJDYXNlKCkgKyBzLnNsaWNlKDEpXG59XG5cbi8vIEVuc3VyZSB0byBzaXgtYmFzZWQgaGV4XG5mdW5jdGlvbiBmdWxsSGV4KGhleCkge1xuICByZXR1cm4gaGV4Lmxlbmd0aCA9PSA0ID9cbiAgICBbICcjJyxcbiAgICAgIGhleC5zdWJzdHJpbmcoMSwgMiksIGhleC5zdWJzdHJpbmcoMSwgMilcbiAgICAsIGhleC5zdWJzdHJpbmcoMiwgMyksIGhleC5zdWJzdHJpbmcoMiwgMylcbiAgICAsIGhleC5zdWJzdHJpbmcoMywgNCksIGhleC5zdWJzdHJpbmcoMywgNClcbiAgICBdLmpvaW4oJycpIDogaGV4XG59XG5cbi8vIENvbXBvbmVudCB0byBoZXggdmFsdWVcbmZ1bmN0aW9uIGNvbXBUb0hleChjb21wKSB7XG4gIHZhciBoZXggPSBjb21wLnRvU3RyaW5nKDE2KVxuICByZXR1cm4gaGV4Lmxlbmd0aCA9PSAxID8gJzAnICsgaGV4IDogaGV4XG59XG5cbi8vIENhbGN1bGF0ZSBwcm9wb3J0aW9uYWwgd2lkdGggYW5kIGhlaWdodCB2YWx1ZXMgd2hlbiBuZWNlc3NhcnlcbmZ1bmN0aW9uIHByb3BvcnRpb25hbFNpemUoZWxlbWVudCwgd2lkdGgsIGhlaWdodCkge1xuICBpZiAod2lkdGggPT0gbnVsbCB8fCBoZWlnaHQgPT0gbnVsbCkge1xuICAgIHZhciBib3ggPSBlbGVtZW50LmJib3goKVxuXG4gICAgaWYgKHdpZHRoID09IG51bGwpXG4gICAgICB3aWR0aCA9IGJveC53aWR0aCAvIGJveC5oZWlnaHQgKiBoZWlnaHRcbiAgICBlbHNlIGlmIChoZWlnaHQgPT0gbnVsbClcbiAgICAgIGhlaWdodCA9IGJveC5oZWlnaHQgLyBib3gud2lkdGggKiB3aWR0aFxuICB9XG5cbiAgcmV0dXJuIHtcbiAgICB3aWR0aDogIHdpZHRoXG4gICwgaGVpZ2h0OiBoZWlnaHRcbiAgfVxufVxuXG4vLyBEZWx0YSB0cmFuc2Zvcm0gcG9pbnRcbmZ1bmN0aW9uIGRlbHRhVHJhbnNmb3JtUG9pbnQobWF0cml4LCB4LCB5KSB7XG4gIHJldHVybiB7XG4gICAgeDogeCAqIG1hdHJpeC5hICsgeSAqIG1hdHJpeC5jICsgMFxuICAsIHk6IHggKiBtYXRyaXguYiArIHkgKiBtYXRyaXguZCArIDBcbiAgfVxufVxuXG4vLyBNYXAgbWF0cml4IGFycmF5IHRvIG9iamVjdFxuZnVuY3Rpb24gYXJyYXlUb01hdHJpeChhKSB7XG4gIHJldHVybiB7IGE6IGFbMF0sIGI6IGFbMV0sIGM6IGFbMl0sIGQ6IGFbM10sIGU6IGFbNF0sIGY6IGFbNV0gfVxufVxuXG4vLyBQYXJzZSBtYXRyaXggaWYgcmVxdWlyZWRcbmZ1bmN0aW9uIHBhcnNlTWF0cml4KG1hdHJpeCkge1xuICBpZiAoIShtYXRyaXggaW5zdGFuY2VvZiBTVkcuTWF0cml4KSlcbiAgICBtYXRyaXggPSBuZXcgU1ZHLk1hdHJpeChtYXRyaXgpXG5cbiAgcmV0dXJuIG1hdHJpeFxufVxuXG4vLyBBZGQgY2VudHJlIHBvaW50IHRvIHRyYW5zZm9ybSBvYmplY3RcbmZ1bmN0aW9uIGVuc3VyZUNlbnRyZShvLCB0YXJnZXQpIHtcbiAgby5jeCA9IG8uY3ggPT0gbnVsbCA/IHRhcmdldC5iYm94KCkuY3ggOiBvLmN4XG4gIG8uY3kgPSBvLmN5ID09IG51bGwgPyB0YXJnZXQuYmJveCgpLmN5IDogby5jeVxufVxuXG4vLyBDb252ZXJ0IHN0cmluZyB0byBtYXRyaXhcbmZ1bmN0aW9uIHN0cmluZ1RvTWF0cml4KHNvdXJjZSkge1xuICAvLyByZW1vdmUgbWF0cml4IHdyYXBwZXIgYW5kIHNwbGl0IHRvIGluZGl2aWR1YWwgbnVtYmVyc1xuICBzb3VyY2UgPSBzb3VyY2VcbiAgICAucmVwbGFjZShTVkcucmVnZXgud2hpdGVzcGFjZSwgJycpXG4gICAgLnJlcGxhY2UoU1ZHLnJlZ2V4Lm1hdHJpeCwgJycpXG4gICAgLnNwbGl0KFNWRy5yZWdleC5tYXRyaXhFbGVtZW50cylcblxuICAvLyBjb252ZXJ0IHN0cmluZyB2YWx1ZXMgdG8gZmxvYXRzIGFuZCBjb252ZXJ0IHRvIGEgbWF0cml4LWZvcm1hdHRlZCBvYmplY3RcbiAgcmV0dXJuIGFycmF5VG9NYXRyaXgoXG4gICAgU1ZHLnV0aWxzLm1hcChzb3VyY2UsIGZ1bmN0aW9uKG4pIHtcbiAgICAgIHJldHVybiBwYXJzZUZsb2F0KG4pXG4gICAgfSlcbiAgKVxufVxuXG4vLyBDYWxjdWxhdGUgcG9zaXRpb24gYWNjb3JkaW5nIHRvIGZyb20gYW5kIHRvXG5mdW5jdGlvbiBhdChvLCBwb3MpIHtcbiAgLy8gbnVtYmVyIHJlY2FsY3VsYXRpb24gKGRvbid0IGJvdGhlciBjb252ZXJ0aW5nIHRvIFNWRy5OdW1iZXIgZm9yIHBlcmZvcm1hbmNlIHJlYXNvbnMpXG4gIHJldHVybiB0eXBlb2Ygby5mcm9tID09ICdudW1iZXInID9cbiAgICBvLmZyb20gKyAoby50byAtIG8uZnJvbSkgKiBwb3MgOlxuXG4gIC8vIGluc3RhbmNlIHJlY2FsY3VsYXRpb25cbiAgbyBpbnN0YW5jZW9mIFNWRy5Db2xvciB8fCBvIGluc3RhbmNlb2YgU1ZHLk51bWJlciB8fCBvIGluc3RhbmNlb2YgU1ZHLk1hdHJpeCA/IG8uYXQocG9zKSA6XG5cbiAgLy8gZm9yIGFsbCBvdGhlciB2YWx1ZXMgd2FpdCB1bnRpbCBwb3MgaGFzIHJlYWNoZWQgMSB0byByZXR1cm4gdGhlIGZpbmFsIHZhbHVlXG4gIHBvcyA8IDEgPyBvLmZyb20gOiBvLnRvXG59XG5cbi8vIFBhdGhBcnJheSBIZWxwZXJzXG5mdW5jdGlvbiBhcnJheVRvU3RyaW5nKGEpIHtcbiAgZm9yICh2YXIgaSA9IDAsIGlsID0gYS5sZW5ndGgsIHMgPSAnJzsgaSA8IGlsOyBpKyspIHtcbiAgICBzICs9IGFbaV1bMF1cblxuICAgIGlmIChhW2ldWzFdICE9IG51bGwpIHtcbiAgICAgIHMgKz0gYVtpXVsxXVxuXG4gICAgICBpZiAoYVtpXVsyXSAhPSBudWxsKSB7XG4gICAgICAgIHMgKz0gJyAnXG4gICAgICAgIHMgKz0gYVtpXVsyXVxuXG4gICAgICAgIGlmIChhW2ldWzNdICE9IG51bGwpIHtcbiAgICAgICAgICBzICs9ICcgJ1xuICAgICAgICAgIHMgKz0gYVtpXVszXVxuICAgICAgICAgIHMgKz0gJyAnXG4gICAgICAgICAgcyArPSBhW2ldWzRdXG5cbiAgICAgICAgICBpZiAoYVtpXVs1XSAhPSBudWxsKSB7XG4gICAgICAgICAgICBzICs9ICcgJ1xuICAgICAgICAgICAgcyArPSBhW2ldWzVdXG4gICAgICAgICAgICBzICs9ICcgJ1xuICAgICAgICAgICAgcyArPSBhW2ldWzZdXG5cbiAgICAgICAgICAgIGlmIChhW2ldWzddICE9IG51bGwpIHtcbiAgICAgICAgICAgICAgcyArPSAnICdcbiAgICAgICAgICAgICAgcyArPSBhW2ldWzddXG4gICAgICAgICAgICB9XG4gICAgICAgICAgfVxuICAgICAgICB9XG4gICAgICB9XG4gICAgfVxuICB9XG5cbiAgcmV0dXJuIHMgKyAnICdcbn1cblxuLy8gRGVlcCBuZXcgaWQgYXNzaWdubWVudFxuZnVuY3Rpb24gYXNzaWduTmV3SWQobm9kZSkge1xuICAvLyBkbyB0aGUgc2FtZSBmb3IgU1ZHIGNoaWxkIG5vZGVzIGFzIHdlbGxcbiAgZm9yICh2YXIgaSA9IG5vZGUuY2hpbGROb2Rlcy5sZW5ndGggLSAxOyBpID49IDA7IGktLSlcbiAgICBpZiAobm9kZS5jaGlsZE5vZGVzW2ldIGluc3RhbmNlb2YgU1ZHRWxlbWVudClcbiAgICAgIGFzc2lnbk5ld0lkKG5vZGUuY2hpbGROb2Rlc1tpXSlcblxuICByZXR1cm4gU1ZHLmFkb3B0KG5vZGUpLmlkKFNWRy5laWQobm9kZS5ub2RlTmFtZSkpXG59XG5cbi8vIEFkZCBtb3JlIGJvdW5kaW5nIGJveCBwcm9wZXJ0aWVzXG5mdW5jdGlvbiBmdWxsQm94KGIpIHtcbiAgaWYgKGIueCA9PSBudWxsKSB7XG4gICAgYi54ICAgICAgPSAwXG4gICAgYi55ICAgICAgPSAwXG4gICAgYi53aWR0aCAgPSAwXG4gICAgYi5oZWlnaHQgPSAwXG4gIH1cblxuICBiLncgID0gYi53aWR0aFxuICBiLmggID0gYi5oZWlnaHRcbiAgYi54MiA9IGIueCArIGIud2lkdGhcbiAgYi55MiA9IGIueSArIGIuaGVpZ2h0XG4gIGIuY3ggPSBiLnggKyBiLndpZHRoIC8gMlxuICBiLmN5ID0gYi55ICsgYi5oZWlnaHQgLyAyXG5cbiAgcmV0dXJuIGJcbn1cblxuLy8gR2V0IGlkIGZyb20gcmVmZXJlbmNlIHN0cmluZ1xuZnVuY3Rpb24gaWRGcm9tUmVmZXJlbmNlKHVybCkge1xuICB2YXIgbSA9IHVybC50b1N0cmluZygpLm1hdGNoKFNWRy5yZWdleC5yZWZlcmVuY2UpXG5cbiAgaWYgKG0pIHJldHVybiBtWzFdXG59XG5cbi8vIENyZWF0ZSBtYXRyaXggYXJyYXkgZm9yIGxvb3BpbmdcbnZhciBhYmNkZWYgPSAnYWJjZGVmJy5zcGxpdCgnJylcbi8vIEFkZCBDdXN0b21FdmVudCB0byBJRTkgYW5kIElFMTBcbmlmICh0eXBlb2YgQ3VzdG9tRXZlbnQgIT09ICdmdW5jdGlvbicpIHtcbiAgLy8gQ29kZSBmcm9tOiBodHRwczovL2RldmVsb3Blci5tb3ppbGxhLm9yZy9lbi1VUy9kb2NzL1dlYi9BUEkvQ3VzdG9tRXZlbnRcbiAgdmFyIEN1c3RvbUV2ZW50ID0gZnVuY3Rpb24oZXZlbnQsIG9wdGlvbnMpIHtcbiAgICBvcHRpb25zID0gb3B0aW9ucyB8fCB7IGJ1YmJsZXM6IGZhbHNlLCBjYW5jZWxhYmxlOiBmYWxzZSwgZGV0YWlsOiB1bmRlZmluZWQgfVxuICAgIHZhciBlID0gZG9jdW1lbnQuY3JlYXRlRXZlbnQoJ0N1c3RvbUV2ZW50JylcbiAgICBlLmluaXRDdXN0b21FdmVudChldmVudCwgb3B0aW9ucy5idWJibGVzLCBvcHRpb25zLmNhbmNlbGFibGUsIG9wdGlvbnMuZGV0YWlsKVxuICAgIHJldHVybiBlXG4gIH1cblxuICBDdXN0b21FdmVudC5wcm90b3R5cGUgPSB3aW5kb3cuRXZlbnQucHJvdG90eXBlXG5cbiAgd2luZG93LkN1c3RvbUV2ZW50ID0gQ3VzdG9tRXZlbnRcbn1cblxuLy8gcmVxdWVzdEFuaW1hdGlvbkZyYW1lIC8gY2FuY2VsQW5pbWF0aW9uRnJhbWUgUG9seWZpbGwgd2l0aCBmYWxsYmFjayBiYXNlZCBvbiBQYXVsIElyaXNoXG4oZnVuY3Rpb24odykge1xuICB2YXIgbGFzdFRpbWUgPSAwXG4gIHZhciB2ZW5kb3JzID0gWydtb3onLCAnd2Via2l0J11cblxuICBmb3IodmFyIHggPSAwOyB4IDwgdmVuZG9ycy5sZW5ndGggJiYgIXdpbmRvdy5yZXF1ZXN0QW5pbWF0aW9uRnJhbWU7ICsreCkge1xuICAgIHcucmVxdWVzdEFuaW1hdGlvbkZyYW1lID0gd1t2ZW5kb3JzW3hdICsgJ1JlcXVlc3RBbmltYXRpb25GcmFtZSddXG4gICAgdy5jYW5jZWxBbmltYXRpb25GcmFtZSAgPSB3W3ZlbmRvcnNbeF0gKyAnQ2FuY2VsQW5pbWF0aW9uRnJhbWUnXSB8fFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgd1t2ZW5kb3JzW3hdICsgJ0NhbmNlbFJlcXVlc3RBbmltYXRpb25GcmFtZSddXG4gIH1cblxuICB3LnJlcXVlc3RBbmltYXRpb25GcmFtZSA9IHcucmVxdWVzdEFuaW1hdGlvbkZyYW1lIHx8XG4gICAgZnVuY3Rpb24oY2FsbGJhY2spIHtcbiAgICAgIHZhciBjdXJyVGltZSA9IG5ldyBEYXRlKCkuZ2V0VGltZSgpXG4gICAgICB2YXIgdGltZVRvQ2FsbCA9IE1hdGgubWF4KDAsIDE2IC0gKGN1cnJUaW1lIC0gbGFzdFRpbWUpKVxuXG4gICAgICB2YXIgaWQgPSB3LnNldFRpbWVvdXQoZnVuY3Rpb24oKSB7XG4gICAgICAgIGNhbGxiYWNrKGN1cnJUaW1lICsgdGltZVRvQ2FsbClcbiAgICAgIH0sIHRpbWVUb0NhbGwpXG5cbiAgICAgIGxhc3RUaW1lID0gY3VyclRpbWUgKyB0aW1lVG9DYWxsXG4gICAgICByZXR1cm4gaWRcbiAgICB9XG5cbiAgdy5jYW5jZWxBbmltYXRpb25GcmFtZSA9IHcuY2FuY2VsQW5pbWF0aW9uRnJhbWUgfHwgdy5jbGVhclRpbWVvdXQ7XG5cbn0od2luZG93KSlcblxucmV0dXJuIFNWR1xuXG59KSk7IiwiaW1wb3J0IFNWRyBmcm9tICdzdmcuanMnO1xuaW1wb3J0IGZpdEN1cnZlIGZyb20gJ2ZpdC1jdXJ2ZSc7XG5cbmZ1bmN0aW9uIEdseXBoc0NvbnRyb2wocGFubmVsKSB7XG5cblx0dGhpcy5zdGFydCA9IGZ1bmN0aW9uKCBwb2ludCApIHtcblx0XHRyYXdQb2ludERhdGEucHVzaCggcG9pbnQgKTtcbiAgICAgICAgcGFpbnRpbmdQb2x5TGluZSA9IHBhbm5lbC5wb2x5bGluZSgpLmZpbGwoJ25vbmUnKS5zdHJva2UoeyB3aWR0aDogMSB9KVxuXG5cdH07XG5cblx0dGhpcy51cGRhdGUgPSBmdW5jdGlvbiggcG9pbnQgKSB7XG5cdFx0cmF3UG9pbnREYXRhLnB1c2goIHBvaW50ICk7XG4gICAgICAgIHVwZGF0ZUxpbmVzKCBwYWludGluZ1BvbHlMaW5lLCByYXdQb2ludERhdGEpO1xuXHR9O1xuXHRcblx0dGhpcy5lbmQgPSBmdW5jdGlvbigpIHtcblx0XHRsZXQgc21vb3RoQml6ZXIgPSBmaXRDdXJ2ZSggcmF3UG9pbnREYXRhLCBlcnJvciApO1xuICAgICAgICBsZXQgcGF0aFN0cmluZyA9IGZpdHRlZEN1cnZlRGF0YVRvUGF0aFN0cmluZyhzbW9vdGhCaXplcik7XG4gICAgICAgIHBhbm5lbC5wYXRoKCBwYXRoU3RyaW5nICkuZmlsbCgnbm9uZScpLnN0cm9rZSh7IHdpZHRoOiAzIH0pLnN0cm9rZSgnI2YwNicpO1xuICAgICAgICByYXdQb2ludERhdGEubGVuZ3RoID0gMDtcblx0fTtcbn1cblxuZXhwb3J0IGRlZmF1bHQgR2x5cGhzQ29udHJvbDsiLCJpbXBvcnQgU1ZHIGZyb20gJ3N2Zy5qcyc7XG5pbXBvcnQgZml0Q3VydmUgZnJvbSAnZml0LWN1cnZlJztcblxuY29uc3QgZXJyb3IgPSAxMDA7XG5cbmZ1bmN0aW9uIFBhaW50Q29udHJvbChwYW5uZWwpIHtcblx0bGV0IHJhd1BvaW50RGF0YSA9IFtdO1xuXHRsZXQgcGFpbnRpbmdQb2x5TGluZSA9IHVuZGVmaW5lZDtcblxuXHR0aGlzLnN0YXJ0ID0gZnVuY3Rpb24oIHBvaW50ICkge1xuXHRcdHJhd1BvaW50RGF0YS5wdXNoKCBwb2ludCApO1xuICAgICAgICBwYWludGluZ1BvbHlMaW5lID0gcGFubmVsLnBvbHlsaW5lKCkuZmlsbCgnbm9uZScpLnN0cm9rZSh7IHdpZHRoOiAxIH0pXG5cblx0fTtcblx0dGhpcy51cGRhdGUgPSBmdW5jdGlvbiggcG9pbnQgKSB7XG5cdFx0cmF3UG9pbnREYXRhLnB1c2goIHBvaW50ICk7XG4gICAgICAgIHVwZGF0ZUxpbmVzKCBwYWludGluZ1BvbHlMaW5lLCByYXdQb2ludERhdGEpO1xuXHR9O1xuXG5cdHRoaXMuZW5kID0gZnVuY3Rpb24oKSB7XG5cdFx0bGV0IHNtb290aEJpemVyID0gZml0Q3VydmUoIHJhd1BvaW50RGF0YSwgZXJyb3IgKTtcbiAgICAgICAgbGV0IHBhdGhTdHJpbmcgPSBmaXR0ZWRDdXJ2ZURhdGFUb1BhdGhTdHJpbmcoc21vb3RoQml6ZXIpO1xuICAgICAgICBwYW5uZWwucGF0aCggcGF0aFN0cmluZyApLmZpbGwoJ25vbmUnKS5zdHJva2UoeyB3aWR0aDogMyB9KS5zdHJva2UoJyNmMDYnKTtcbiAgICAgICAgcmF3UG9pbnREYXRhLmxlbmd0aCA9IDA7XG5cdH07XG59XG5cbmV4cG9ydCBkZWZhdWx0IFBhaW50Q29udHJvbDtcblxuZnVuY3Rpb24gdXBkYXRlTGluZXMocGFpbnRpbmdQb2x5TGluZSwgcmF3UG9pbnREYXRhKSB7XG5cdHBhaW50aW5nUG9seUxpbmUucGxvdCggcmF3UG9pbnREYXRhICk7XG59XG5mdW5jdGlvbiBmaXR0ZWRDdXJ2ZURhdGFUb1BhdGhTdHJpbmcoZml0dGVkTGluZURhdGEpIHtcbiAgICB2YXIgc3RyID0gXCJcIjtcbiAgICBmaXR0ZWRMaW5lRGF0YS5tYXAoZnVuY3Rpb24gKGJlemllciwgaSkge1xuICAgICAgICBpZiAoaSA9PSAwKSB7XG4gICAgICAgICAgICBzdHIgKz0gXCJNIFwiICsgYmV6aWVyWzBdWzBdICsgXCIgXCIgKyBiZXppZXJbMF1bMV07XG4gICAgICAgIH1cbiAgICAgICAgc3RyICs9IFwiQyBcIiArIGJlemllclsxXVswXSArIFwiIFwiICsgYmV6aWVyWzFdWzFdICsgXCIsIFwiICtcbiAgICAgICAgICAgIGJlemllclsyXVswXSArIFwiIFwiICsgYmV6aWVyWzJdWzFdICsgXCIsIFwiICtcbiAgICAgICAgICAgIGJlemllclszXVswXSArIFwiIFwiICsgYmV6aWVyWzNdWzFdICsgXCIgXCI7XG4gICAgfSk7XG5cbiAgICByZXR1cm4gc3RyO1xufSIsImltcG9ydCBHbHlwaHNDb250cm9sIGZyb20gJy4vQ29udHJvbHMvR2x5cGhzQ29udHJvbCc7XG5pbXBvcnQgUGFpbnRDb250cm9sIGZyb20gJy4vQ29udHJvbHMvUGFpbnRDb250cm9sJztcbmltcG9ydCBTVkcgZnJvbSAnc3ZnLmpzJztcblxudmFyIGRyYXcgPSBTVkcoJ2RyYXdpbmcnKS5zaXplKDMwMCwgMzAwKTtcblxuXG5zZXRDb250cm9sKGRyYXcpO1xuXG5cblxuZnVuY3Rpb24gc2V0Q29udHJvbChfY29udGFpbmVyKSB7XG5cdGxldCBpc01vdXNlRG93biA9IGZhbHNlO1xuXHRsZXQgcG9seWxpbmU7XG5cblx0bGV0IGN1cnJuZXRDb250cm9sID0gbmV3IFBhaW50Q29udHJvbChkcmF3KTtcblxuICAgIF9jb250YWluZXIub24oJ21vdXNlZG93bicsIGZ1bmN0aW9uIChlKSB7XG4gICAgICAgIGNvbnN0IHBvaW50ID0gW1xuICAgICAgICBcdGUuY2xpZW50WCxcbiAgICAgICAgXHRlLmNsaWVudFlcbiAgICAgICAgXVxuICAgICAgICBpc01vdXNlRG93biA9IHRydWU7XG4gICAgICAgIGN1cnJuZXRDb250cm9sLnN0YXJ0KHBvaW50KTtcblxuICAgIH0pO1xuICAgIF9jb250YWluZXIub24oJ21vdXNldXAnLCBmdW5jdGlvbiAoKSB7XG4gICAgICAgIGlzTW91c2VEb3duID0gZmFsc2U7XG4gICAgICAgIGN1cnJuZXRDb250cm9sLmVuZCgpO1xuICAgIH0pO1xuICAgIF9jb250YWluZXIub24oJ21vdXNlbW92ZScsIGZ1bmN0aW9uIChlKSB7XG4gICAgICAgIHZhciB4ID0gZS5vZmZzZXRYO1xuICAgICAgICB2YXIgeSA9IGUub2Zmc2V0WTtcbiAgICAgICAgaWYgKGlzTW91c2VEb3duKSB7XG4gICAgICAgICAgICBjdXJybmV0Q29udHJvbC51cGRhdGUoW3gsIHldKTtcbiAgICAgICAgfVxuICAgIH0pO1xufVxuXG4iXX0=
