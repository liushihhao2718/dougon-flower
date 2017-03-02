(function e(t,n,r){function s(o,u){if(!n[o]){if(!t[o]){var a=typeof require=="function"&&require;if(!u&&a)return a(o,!0);if(i)return i(o,!0);var f=new Error("Cannot find module '"+o+"'");throw f.code="MODULE_NOT_FOUND",f}var l=n[o]={exports:{}};t[o][0].call(l.exports,function(e){var n=t[o][1][e];return s(n?n:e)},l,l.exports,e,t,n,r)}return n[o].exports}var i=typeof require=="function"&&require;for(var o=0;o<r.length;o++)s(r[o]);return s})({1:[function(require,module,exports){
module.exports = require('./lib/bezier');

},{"./lib/bezier":2}],2:[function(require,module,exports){
/**
  A javascript Bezier curve library by Pomax.

  Based on http://pomax.github.io/bezierinfo

  This code is MIT licensed.
**/
(function() {
  "use strict";

  // math-inlining.
  var abs = Math.abs,
      min = Math.min,
      max = Math.max,
      acos = Math.acos,
      sqrt = Math.sqrt,
      pi = Math.PI,
      // a zero coordinate, which is surprisingly useful
      ZERO = {x:0,y:0,z:0};

  // quite needed
  var utils = require('./utils.js');

  // not quite needed, but eventually this'll be useful...
  var PolyBezier = require('./poly-bezier.js');

  /**
   * Bezier curve constructor. The constructor argument can be one of three things:
   *
   * 1. array/4 of {x:..., y:..., z:...}, z optional
   * 2. numerical array/8 ordered x1,y1,x2,y2,x3,y3,x4,y4
   * 3. numerical array/12 ordered x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4
   *
   */
  var Bezier = function(coords) {
    var args = (coords && coords.forEach) ? coords : [].slice.call(arguments);
    var coordlen = false;
    if(typeof args[0] === "object") {
      coordlen = args.length;
      var newargs = [];
      args.forEach(function(point) {
        ['x','y','z'].forEach(function(d) {
          if(typeof point[d] !== "undefined") {
            newargs.push(point[d]);
          }
        });
      });
      args = newargs;
    }
    var higher = false;
    var len = args.length;
    if (coordlen) {
      if(coordlen>4) {
        if (arguments.length !== 1) {
          throw new Error("Only new Bezier(point[]) is accepted for 4th and higher order curves");
        }
        higher = true;
      }
    } else {
      if(len!==6 && len!==8 && len!==9 && len!==12) {
        if (arguments.length !== 1) {
          throw new Error("Only new Bezier(point[]) is accepted for 4th and higher order curves");
        }
      }
    }
    var _3d = (!higher && (len === 9 || len === 12)) || (coords && coords[0] && typeof coords[0].z !== "undefined");
    this._3d = _3d;
    var points = [];
    for(var idx=0, step=(_3d ? 3 : 2); idx<len; idx+=step) {
      var point = {
        x: args[idx],
        y: args[idx+1]
      };
      if(_3d) { point.z = args[idx+2] };
      points.push(point);
    }
    this.order = points.length - 1;
    this.points = points;
    var dims = ['x','y'];
    if(_3d) dims.push('z');
    this.dims = dims;
    this.dimlen = dims.length;

    (function(curve) {
      var order = curve.order;
      var points = curve.points;
      var a = utils.align(points, {p1:points[0], p2:points[order]});
      for(var i=0; i<a.length; i++) {
        if(abs(a[i].y) > 0.0001) {
          curve._linear = false;
          return;
        }
      }
      curve._linear = true;
    }(this));

    this._t1 = 0;
    this._t2 = 1;
    this.update();
  };

  Bezier.fromSVG = function(svgString) {
    var list = svgString.match(/[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?/g).map(parseFloat);
    var relative = /[cq]/.test(svgString);
    if(!relative) return new Bezier(list);
    list = list.map(function(v,i) {
      return i < 2 ? v : v + list[i % 2];
    });
    return new Bezier(list);
  };

  function getABC(n,S,B,E,t) {
    if(typeof t === "undefined") { t = 0.5; }
    var u = utils.projectionratio(t,n),
        um = 1-u,
        C = {
          x: u*S.x + um*E.x,
          y: u*S.y + um*E.y
        },
        s = utils.abcratio(t,n),
        A = {
          x: B.x + (B.x-C.x)/s,
          y: B.y + (B.y-C.y)/s
        };
    return { A:A, B:B, C:C };
  }

  Bezier.quadraticFromPoints = function(p1,p2,p3, t) {
    if(typeof t === "undefined") { t = 0.5; }
    // shortcuts, although they're really dumb
    if(t===0) { return new Bezier(p2,p2,p3); }
    if(t===1) { return new Bezier(p1,p2,p2); }
    // real fitting.
    var abc = getABC(2,p1,p2,p3,t);
    return new Bezier(p1, abc.A, p3);
  };

  Bezier.cubicFromPoints = function(S,B,E, t,d1) {
    if(typeof t === "undefined") { t = 0.5; }
    var abc = getABC(3,S,B,E,t);
    if(typeof d1 === "undefined") { d1 = utils.dist(B,abc.C); }
    var d2 = d1 * (1-t)/t;

    var selen = utils.dist(S,E),
        lx = (E.x-S.x)/selen,
        ly = (E.y-S.y)/selen,
        bx1 = d1 * lx,
        by1 = d1 * ly,
        bx2 = d2 * lx,
        by2 = d2 * ly;
    // derivation of new hull coordinates
    var e1  = { x: B.x - bx1, y: B.y - by1 },
        e2  = { x: B.x + bx2, y: B.y + by2 },
        A = abc.A,
        v1  = { x: A.x + (e1.x-A.x)/(1-t), y: A.y + (e1.y-A.y)/(1-t) },
        v2  = { x: A.x + (e2.x-A.x)/(t), y: A.y + (e2.y-A.y)/(t) },
        nc1 = { x: S.x + (v1.x-S.x)/(t), y: S.y + (v1.y-S.y)/(t) },
        nc2 = { x: E.x + (v2.x-E.x)/(1-t), y: E.y + (v2.y-E.y)/(1-t) };
    // ...done
    return new Bezier(S,nc1,nc2,E);
  };

  var getUtils = function() {
    return utils;
  };

  Bezier.getUtils = getUtils;

  Bezier.prototype = {
    getUtils: getUtils,
    valueOf: function() {
      return this.toString();
    },
    toString: function() {
      return utils.pointsToString(this.points);
    },
    toSVG: function(relative) {
      if(this._3d) return false;
      var p = this.points,
          x = p[0].x,
          y = p[0].y,
          s = ["M", x, y, (this.order===2 ? "Q":"C")];
      for(var i=1, last=p.length; i<last; i++) {
        s.push(p[i].x);
        s.push(p[i].y);
      }
      return s.join(" ");
    },
    update: function() {
      // one-time compute derivative coordinates
      this.dpoints = [];
      for(var p=this.points, d=p.length, c=d-1; d>1; d--, c--) {
        var list = [];
        for(var j=0, dpt; j<c; j++) {
          dpt = {
            x: c * (p[j+1].x - p[j].x),
            y: c * (p[j+1].y - p[j].y)
          };
          if(this._3d) {
            dpt.z = c * (p[j+1].z - p[j].z);
          }
          list.push(dpt);
        }
        this.dpoints.push(list);
        p = list;
      };
      this.computedirection();
    },
    computedirection: function() {
      var points = this.points;
      var angle = utils.angle(points[0], points[this.order], points[1]);
      this.clockwise = angle > 0;
    },
    length: function() {
      return utils.length(this.derivative.bind(this));
    },
    _lut: [],
    getLUT: function(steps) {
      steps = steps || 100;
      if (this._lut.length === steps) { return this._lut; }
      this._lut = [];
      for(var t=0; t<=steps; t++) {
        this._lut.push(this.compute(t/steps));
      }
      return this._lut;
    },
    on: function(point, error) {
      error = error || 5;
      var lut = this.getLUT(), hits = [], c, t=0;
      for(var i=0; i<lut.length; i++) {
        c = lut[i];
        if (utils.dist(c,point) < error) {
          hits.push(c)
          t += i / lut.length;
        }
      }
      if(!hits.length) return false;
      return t /= hits.length;
    },
    project: function(point) {
      // step 1: coarse check
      var LUT = this.getLUT(), l = LUT.length-1,
          closest = utils.closest(LUT, point),
          mdist = closest.mdist,
          mpos = closest.mpos;
      if (mpos===0 || mpos===l) {
        var t = mpos/l, pt = this.compute(t);
        pt.t = t;
        pt.d = mdist;
        return pt;
      }

      // step 2: fine check
      var ft, t, p, d,
          t1 = (mpos-1)/l,
          t2 = (mpos+1)/l,
          step = 0.1/l;
      mdist += 1;
      for(t=t1,ft=t; t<t2+step; t+=step) {
        p = this.compute(t);
        d = utils.dist(point, p);
        if (d<mdist) {
          mdist = d;
          ft = t;
        }
      }
      p = this.compute(ft);
      p.t = ft;
      p.d = mdist;
      return p;
    },
    get: function(t) {
      return this.compute(t);
    },
    point: function(idx) {
      return this.points[idx];
    },
    compute: function(t) {
      // shortcuts
      if(t===0) { return this.points[0]; }
      if(t===1) { return this.points[this.order]; }

      var p = this.points;
      var mt = 1-t;

      // linear?
      if(this.order===1) {
        ret = {
          x: mt*p[0].x + t*p[1].x,
          y: mt*p[0].y + t*p[1].y
        };
        if (this._3d) { ret.z = mt*p[0].z + t*p[1].z; }
        return ret;
      }

      // quadratic/cubic curve?
      if(this.order<4) {
        var mt2 = mt*mt,
            t2 = t*t,
            a,b,c,d = 0;
        if(this.order===2) {
          p = [p[0], p[1], p[2], ZERO];
          a = mt2;
          b = mt*t*2;
          c = t2;
        }
        else if(this.order===3) {
          a = mt2*mt;
          b = mt2*t*3;
          c = mt*t2*3;
          d = t*t2;
        }
        var ret = {
          x: a*p[0].x + b*p[1].x + c*p[2].x + d*p[3].x,
          y: a*p[0].y + b*p[1].y + c*p[2].y + d*p[3].y
        };
        if(this._3d) {
          ret.z = a*p[0].z + b*p[1].z + c*p[2].z + d*p[3].z;
        }
        return ret;
      }

      // higher order curves: use de Casteljau's computation
      var dCpts = JSON.parse(JSON.stringify(this.points));
      while(dCpts.length > 1) {
        for (var i=0; i<dCpts.length-1; i++) {
          dCpts[i] = {
            x: dCpts[i].x + (dCpts[i+1].x - dCpts[i].x) * t,
            y: dCpts[i].y + (dCpts[i+1].y - dCpts[i].y) * t
          };
          if (typeof dCpts[i].z !== "undefined") {
            dCpts[i] = dCpts[i].z + (dCpts[i+1].z - dCpts[i].z) * t
          }
        }
        dCpts.splice(dCpts.length-1, 1);
      }
      return dCpts[0];
    },
    raise: function() {
      var p = this.points, np = [p[0]], i, k=p.length, pi, pim;
      for (var i=1; i<k; i++) {
        pi = p[i];
        pim = p[i-1];
        np[i] = {
          x: (k-i)/k * pi.x + i/k * pim.x,
          y: (k-i)/k * pi.y + i/k * pim.y
        };
      }
      np[k] = p[k-1];
      return new Bezier(np);
    },
    derivative: function(t) {
      var mt = 1-t,
          a,b,c=0,
          p = this.dpoints[0];
      if(this.order===2) { p = [p[0], p[1], ZERO]; a = mt; b = t; }
      if(this.order===3) { a = mt*mt; b = mt*t*2; c = t*t; }
      var ret = {
        x: a*p[0].x + b*p[1].x + c*p[2].x,
        y: a*p[0].y + b*p[1].y + c*p[2].y
      };
      if(this._3d) {
        ret.z = a*p[0].z + b*p[1].z + c*p[2].z;
      }
      return ret;
    },
    inflections: function() {
      return utils.inflections(this.points);
    },
    normal: function(t) {
      return this._3d ? this.__normal3(t) : this.__normal2(t);
    },
    __normal2: function(t) {
      var d = this.derivative(t);
      var q = sqrt(d.x*d.x + d.y*d.y)
      return { x: -d.y/q, y: d.x/q };
    },
    __normal3: function(t) {
      // see http://stackoverflow.com/questions/25453159
      var r1 = this.derivative(t),
          r2 = this.derivative(t+0.01),
          q1 = sqrt(r1.x*r1.x + r1.y*r1.y + r1.z*r1.z),
          q2 = sqrt(r2.x*r2.x + r2.y*r2.y + r2.z*r2.z);
      r1.x /= q1; r1.y /= q1; r1.z /= q1;
      r2.x /= q2; r2.y /= q2; r2.z /= q2;
      // cross product
      var c = {
        x: r2.y*r1.z - r2.z*r1.y,
        y: r2.z*r1.x - r2.x*r1.z,
        z: r2.x*r1.y - r2.y*r1.x
      };
      var m = sqrt(c.x*c.x + c.y*c.y + c.z*c.z);
      c.x /= m; c.y /= m; c.z /= m;
      // rotation matrix
      var R = [   c.x*c.x,   c.x*c.y-c.z, c.x*c.z+c.y,
                c.x*c.y+c.z,   c.y*c.y,   c.y*c.z-c.x,
                c.x*c.z-c.y, c.y*c.z+c.x,   c.z*c.z    ];
      // normal vector:
      var n = {
        x: R[0] * r1.x + R[1] * r1.y + R[2] * r1.z,
        y: R[3] * r1.x + R[4] * r1.y + R[5] * r1.z,
        z: R[6] * r1.x + R[7] * r1.y + R[8] * r1.z
      };
      return n;
    },
    hull: function(t) {
      var p = this.points,
          _p = [],
          pt,
          q = [],
          idx = 0,
          i=0,
          l=0;
      q[idx++] = p[0];
      q[idx++] = p[1];
      q[idx++] = p[2];
      if(this.order === 3) { q[idx++] = p[3]; }
      // we lerp between all points at each iteration, until we have 1 point left.
      while(p.length>1) {
        _p = [];
        for(i=0, l=p.length-1; i<l; i++) {
          pt = utils.lerp(t,p[i],p[i+1]);
          q[idx++] = pt;
          _p.push(pt);
        }
        p = _p;
      }
      return q;
    },
    split: function(t1, t2) {
      // shortcuts
      if(t1===0 && !!t2) { return this.split(t2).left; }
      if(t2===1) { return this.split(t1).right; }

      // no shortcut: use "de Casteljau" iteration.
      var q = this.hull(t1);
      var result = {
        left: this.order === 2 ? new Bezier([q[0],q[3],q[5]]) : new Bezier([q[0],q[4],q[7],q[9]]),
        right: this.order === 2 ? new Bezier([q[5],q[4],q[2]]) : new Bezier([q[9],q[8],q[6],q[3]]),
        span: q
      };

      // make sure we bind _t1/_t2 information!
      result.left._t1  = utils.map(0,  0,1, this._t1,this._t2);
      result.left._t2  = utils.map(t1, 0,1, this._t1,this._t2);
      result.right._t1 = utils.map(t1, 0,1, this._t1,this._t2);
      result.right._t2 = utils.map(1,  0,1, this._t1,this._t2);

      // if we have no t2, we're done
      if(!t2) { return result; }

      // if we have a t2, split again:
      t2 = utils.map(t2,t1,1,0,1);
      var subsplit = result.right.split(t2);
      return subsplit.left;
    },
    extrema: function() {
      var dims = this.dims,
          result={},
          roots=[],
          p, mfn;
      dims.forEach(function(dim) {
        mfn = function(v) { return v[dim]; };
        p = this.dpoints[0].map(mfn);
        result[dim] = utils.droots(p);
        if(this.order === 3) {
          p = this.dpoints[1].map(mfn);
          result[dim] = result[dim].concat(utils.droots(p));
        }
        result[dim] = result[dim].filter(function(t) { return (t>=0 && t<=1); });
        roots = roots.concat(result[dim].sort());
      }.bind(this));
      roots = roots.sort().filter(function(v,idx) { return (roots.indexOf(v) === idx); });
      result.values = roots;
      return result;
    },
    bbox: function() {
      var extrema = this.extrema(), result = {};
      this.dims.forEach(function(d) {
        result[d] = utils.getminmax(this, d, extrema[d]);
      }.bind(this));
      return result;
    },
    overlaps: function(curve) {
      var lbbox = this.bbox(),
          tbbox = curve.bbox();
      return utils.bboxoverlap(lbbox,tbbox);
    },
    offset: function(t, d) {
      if(typeof d !== "undefined") {
        var c = this.get(t);
        var n = this.normal(t);
        var ret = {
          c: c,
          n: n,
          x: c.x + n.x * d,
          y: c.y + n.y * d
        };
        if(this._3d) {
          ret.z = c.z + n.z * d;
        };
        return ret;
      }
      if(this._linear) {
        var nv = this.normal(0);
        var coords = this.points.map(function(p) {
          var ret = {
            x: p.x + t * nv.x,
            y: p.y + t * nv.y
          };
          if(p.z && n.z) { ret.z = p.z + t * nv.z; }
          return ret;
        });
        return [new Bezier(coords)];
      }
      var reduced = this.reduce();
      return reduced.map(function(s) {
        return s.scale(t);
      });
    },
    simple: function() {
      if(this.order===3) {
        var a1 = utils.angle(this.points[0], this.points[3], this.points[1]);
        var a2 = utils.angle(this.points[0], this.points[3], this.points[2]);
        if(a1>0 && a2<0 || a1<0 && a2>0) return false;
      }
      var n1 = this.normal(0);
      var n2 = this.normal(1);
      var s = n1.x*n2.x + n1.y*n2.y;
      if(this._3d) { s += n1.z*n2.z; }
      var angle = abs(acos(s));
      return angle < pi/3;
    },
    reduce: function() {
      var i, t1=0, t2=0, step=0.01, segment, pass1=[], pass2=[];
      // first pass: split on extrema
      var extrema = this.extrema().values;
      if(extrema.indexOf(0)===-1) { extrema = [0].concat(extrema); }
      if(extrema.indexOf(1)===-1) { extrema.push(1); }

      for(t1=extrema[0], i=1; i<extrema.length; i++) {
        t2 = extrema[i];
        segment = this.split(t1,t2);
        segment._t1 = t1;
        segment._t2 = t2;
        pass1.push(segment);
        t1 = t2;
      }

      // second pass: further reduce these segments to simple segments
      pass1.forEach(function(p1) {
        t1=0;
        t2=0;
        while(t2 <= 1) {
          for(t2=t1+step; t2<=1+step; t2+=step) {
            segment = p1.split(t1,t2);
            if(!segment.simple()) {
              t2 -= step;
              if(abs(t1-t2)<step) {
                // we can never form a reduction
                return [];
              }
              segment = p1.split(t1,t2);
              segment._t1 = utils.map(t1,0,1,p1._t1,p1._t2);
              segment._t2 = utils.map(t2,0,1,p1._t1,p1._t2);
              pass2.push(segment);
              t1 = t2;
              break;
            }
          }
        }
        if(t1<1) {
          segment = p1.split(t1,1);
          segment._t1 = utils.map(t1,0,1,p1._t1,p1._t2);
          segment._t2 = p1._t2;
          pass2.push(segment);
        }
      });
      return pass2;
    },
    scale: function(d) {
      var order = this.order;
      var distanceFn = false
      if(typeof d === "function") { distanceFn = d; }
      if(distanceFn && order === 2) { return this.raise().scale(distanceFn); }

      // TODO: add special handling for degenerate (=linear) curves.
      var clockwise = this.clockwise;
      var r1 = distanceFn ? distanceFn(0) : d;
      var r2 = distanceFn ? distanceFn(1) : d;
      var v = [ this.offset(0,10), this.offset(1,10) ];
      var o = utils.lli4(v[0], v[0].c, v[1], v[1].c);
      if(!o) { throw new Error("cannot scale this curve. Try reducing it first."); }
      // move all points by distance 'd' wrt the origin 'o'
      var points=this.points, np=[];

      // move end points by fixed distance along normal.
      [0,1].forEach(function(t) {
        var p = np[t*order] = utils.copy(points[t*order]);
        p.x += (t?r2:r1) * v[t].n.x;
        p.y += (t?r2:r1) * v[t].n.y;
      }.bind(this));

      if (!distanceFn) {
        // move control points to lie on the intersection of the offset
        // derivative vector, and the origin-through-control vector
        [0,1].forEach(function(t) {
          if(this.order===2 && !!t) return;
          var p = np[t*order];
          var d = this.derivative(t);
          var p2 = { x: p.x + d.x, y: p.y + d.y };
          np[t+1] = utils.lli4(p, p2, o, points[t+1]);
        }.bind(this));
        return new Bezier(np);
      }

      // move control points by "however much necessary to
      // ensure the correct tangent to endpoint".
      [0,1].forEach(function(t) {
        if(this.order===2 && !!t) return;
        var p = points[t+1];
        var ov = {
          x: p.x - o.x,
          y: p.y - o.y
        };
        var rc = distanceFn ? distanceFn((t+1)/order) : d;
        if(distanceFn && !clockwise) rc = -rc;
        var m = sqrt(ov.x*ov.x + ov.y*ov.y);
        ov.x /= m;
        ov.y /= m;
        np[t+1] = {
          x: p.x + rc*ov.x,
          y: p.y + rc*ov.y
        }
      }.bind(this));
      return new Bezier(np);
    },
    outline: function(d1, d2, d3, d4) {
      d2 = (typeof d2 === "undefined") ? d1 : d2;
      var reduced = this.reduce(),
          len = reduced.length,
          fcurves = [],
          bcurves = [],
          p,
          alen = 0,
          tlen = this.length();

      var graduated = (typeof d3 !== "undefined" && typeof d4 !== "undefined");

      function linearDistanceFunction(s,e, tlen,alen,slen) {
        return function (v) {
          var f1 = alen/tlen, f2 = (alen+slen)/tlen, d = e-s;
          return utils.map(v, 0,1, s+f1*d, s+f2*d);
        };
      };

      // form curve oulines
      reduced.forEach(function(segment) {
        slen = segment.length();
        if (graduated) {
          fcurves.push(segment.scale(  linearDistanceFunction( d1, d3, tlen,alen,slen)  ));
          bcurves.push(segment.scale(  linearDistanceFunction(-d2,-d4, tlen,alen,slen)  ));
        } else {
          fcurves.push(segment.scale( d1));
          bcurves.push(segment.scale(-d2));
        }
        alen += slen;
      });

      // reverse the "return" outline
      bcurves = bcurves.map(function(s) {
        p = s.points;
        if(p[3]) { s.points = [p[3],p[2],p[1],p[0]]; }
        else { s.points = [p[2],p[1],p[0]]; }
        return s;
      }).reverse();

      // form the endcaps as lines
      var fs = fcurves[0].points[0],
          fe = fcurves[len-1].points[fcurves[len-1].points.length-1],
          bs = bcurves[len-1].points[bcurves[len-1].points.length-1],
          be = bcurves[0].points[0],
          ls = utils.makeline(bs,fs),
          le = utils.makeline(fe,be),
          segments = [ls].concat(fcurves).concat([le]).concat(bcurves),
          slen = segments.length;

      return new PolyBezier(segments);
    },
    outlineshapes: function(d1, d2, curveIntersectionThreshold) {
      d2 = d2 || d1;
      var outline = this.outline(d1,d2).curves;
      var shapes = [];
      for(var i=1, len=outline.length; i < len/2; i++) {
        var shape = utils.makeshape(outline[i], outline[len-i], curveIntersectionThreshold);
        shape.startcap.virtual = (i > 1);
        shape.endcap.virtual = (i < len/2-1);
        shapes.push(shape);
      }
      return shapes;
    },
    intersects: function(curve, curveIntersectionThreshold) {
      if(!curve) return this.selfintersects(curveIntersectionThreshold);
      if(curve.p1 && curve.p2) {
        return this.lineIntersects(curve);
      }
      if(curve instanceof Bezier) { curve = curve.reduce(); }
      return this.curveintersects(this.reduce(), curve, curveIntersectionThreshold);
    },
    lineIntersects: function(line) {
      var mx = min(line.p1.x, line.p2.x),
          my = min(line.p1.y, line.p2.y),
          MX = max(line.p1.x, line.p2.x),
          MY = max(line.p1.y, line.p2.y),
          self=this;
      return utils.roots(this.points, line).filter(function(t) {
        var p = self.get(t);
        return utils.between(p.x, mx, MX) && utils.between(p.y, my, MY);
      });
    },
    selfintersects: function(curveIntersectionThreshold) {
      var reduced = this.reduce();
      // "simple" curves cannot intersect with their direct
      // neighbour, so for each segment X we check whether
      // it intersects [0:x-2][x+2:last].
      var i,len=reduced.length-2,results=[],result,left,right;
      for(i=0; i<len; i++) {
        left = reduced.slice(i,i+1);
        right = reduced.slice(i+2);
        result = this.curveintersects(left, right, curveIntersectionThreshold);
        results = results.concat( result );
      }
      return results;
    },
    curveintersects: function(c1, c2, curveIntersectionThreshold) {
      var pairs = [];
      // step 1: pair off any overlapping segments
      c1.forEach(function(l) {
        c2.forEach(function(r) {
          if(l.overlaps(r)) {
            pairs.push({ left: l, right: r });
          }
        });
      });
      // step 2: for each pairing, run through the convergence algorithm.
      var intersections = [];
      pairs.forEach(function(pair) {
        var result = utils.pairiteration(pair.left, pair.right, curveIntersectionThreshold);
        if(result.length > 0) {
          intersections = intersections.concat(result);
        }
      });
      return intersections;
    },
    arcs: function(errorThreshold) {
      errorThreshold = errorThreshold || 0.5;
      var circles = [];
      return this._iterate(errorThreshold, circles);
    },
    _error: function(pc, np1, s, e) {
      var q = (e - s) / 4,
          c1 = this.get(s + q),
          c2 = this.get(e - q),
          ref = utils.dist(pc, np1),
          d1  = utils.dist(pc, c1),
          d2  = utils.dist(pc, c2);
      return abs(d1-ref) + abs(d2-ref);
    },
    _iterate: function(errorThreshold, circles) {
      var s = 0, e = 1, safety;
      // we do a binary search to find the "good `t` closest to no-longer-good"
      do {
        safety=0;

        // step 1: start with the maximum possible arc
        e = 1;

        // points:
        var np1 = this.get(s), np2, np3, arc, prev_arc;

        // booleans:
        var curr_good = false, prev_good = false, done;

        // numbers:
        var m = e, prev_e = 1, step = 0;

        // step 2: find the best possible arc
        do {
          prev_good = curr_good;
          prev_arc = arc;
          m = (s + e)/2;
          step++;

          np2 = this.get(m);
          np3 = this.get(e);

          arc = utils.getccenter(np1, np2, np3);
          
          //also save the t values
          arc.interval = {
            start: s,
            end: e
          };

          var error = this._error(arc, np1, s, e);
          curr_good = (error <= errorThreshold);

          done = prev_good && !curr_good;
          if(!done) prev_e = e;

          // this arc is fine: we can move 'e' up to see if we can find a wider arc
          if(curr_good) {
            // if e is already at max, then we're done for this arc.
            if (e >= 1) {
              prev_e = 1;
              prev_arc = arc;
              break;
            }
            // if not, move it up by half the iteration distance
            e = e + (e-s)/2;
          }

          // this is a bad arc: we need to move 'e' down to find a good arc
          else {
            e = m;
          }
        }
        while(!done && safety++<100);

        if(safety>=100) {
          console.error("arc abstraction somehow failed...");
          break;
        }

        // console.log("[F] arc found", s, prev_e, prev_arc.x, prev_arc.y, prev_arc.s, prev_arc.e);

        prev_arc = (prev_arc ? prev_arc : arc);
        circles.push(prev_arc);
        s = prev_e;
      }
      while(e < 1);
      return circles;
    }
  };

  module.exports = Bezier;

}());

},{"./poly-bezier.js":3,"./utils.js":4}],3:[function(require,module,exports){
(function() {
  "use strict";

  var utils = require('./utils.js');

  /**
   * Poly Bezier
   * @param {[type]} curves [description]
   */
  var PolyBezier = function(curves) {
    this.curves = [];
    this._3d = false;
    if(!!curves) {
      this.curves = curves;
      this._3d = this.curves[0]._3d;
    }
  }

  PolyBezier.prototype = {
    valueOf: function() {
      return this.toString();
    },
    toString: function() {
      return "[" + this.curves.map(function(curve) {
        return utils.pointsToString(curve.points);
      }).join(", ") + "]";
    },
    addCurve: function(curve) {
      this.curves.push(curve);
      this._3d = this._3d || curve._3d;
    },
    length: function() {
      return this.curves.map(function(v) { return v.length(); }).reduce(function(a,b) { return a+b; });
    },
    curve: function(idx) {
      return this.curves[idx];
    },
    bbox: function() {
      var c = this.curves;
      var bbox = c[0].bbox();
      for(var i=1; i<c.length; i++) {
        utils.expandbox(bbox, c[i].bbox());
      }
      return bbox;
    },
    offset: function(d) {
      var offset = [];
      this.curves.forEach(function(v) {
        offset = offset.concat(v.offset(d));
      });
      return new PolyBezier(offset);
    }
  };

  module.exports = PolyBezier;
}());

},{"./utils.js":4}],4:[function(require,module,exports){
(function() {
  "use strict";

  // math-inlining.
  var abs = Math.abs,
      cos = Math.cos,
      sin = Math.sin,
      acos = Math.acos,
      atan2 = Math.atan2,
      sqrt = Math.sqrt,
      pow = Math.pow,
      // cube root function yielding real roots
      crt = function(v) { return (v<0) ? -pow(-v,1/3) : pow(v,1/3); },
      // trig constants
      pi = Math.PI,
      tau = 2*pi,
      quart = pi/2,
      // float precision significant decimal
      epsilon = 0.000001,
      // extremas used in bbox calculation and similar algorithms
      nMax = Number.MAX_SAFE_INTEGER,
      nMin = Number.MIN_SAFE_INTEGER;

  // Bezier utility functions
  var utils = {
    // Legendre-Gauss abscissae with n=24 (x_i values, defined at i=n as the roots of the nth order Legendre polynomial Pn(x))
    Tvalues: [
      -0.0640568928626056260850430826247450385909,
       0.0640568928626056260850430826247450385909,
      -0.1911188674736163091586398207570696318404,
       0.1911188674736163091586398207570696318404,
      -0.3150426796961633743867932913198102407864,
       0.3150426796961633743867932913198102407864,
      -0.4337935076260451384870842319133497124524,
       0.4337935076260451384870842319133497124524,
      -0.5454214713888395356583756172183723700107,
       0.5454214713888395356583756172183723700107,
      -0.6480936519369755692524957869107476266696,
       0.6480936519369755692524957869107476266696,
      -0.7401241915785543642438281030999784255232,
       0.7401241915785543642438281030999784255232,
      -0.8200019859739029219539498726697452080761,
       0.8200019859739029219539498726697452080761,
      -0.8864155270044010342131543419821967550873,
       0.8864155270044010342131543419821967550873,
      -0.9382745520027327585236490017087214496548,
       0.9382745520027327585236490017087214496548,
      -0.9747285559713094981983919930081690617411,
       0.9747285559713094981983919930081690617411,
      -0.9951872199970213601799974097007368118745,
       0.9951872199970213601799974097007368118745
    ],

    // Legendre-Gauss weights with n=24 (w_i values, defined by a function linked to in the Bezier primer article)
    Cvalues: [
      0.1279381953467521569740561652246953718517,
      0.1279381953467521569740561652246953718517,
      0.1258374563468282961213753825111836887264,
      0.1258374563468282961213753825111836887264,
      0.1216704729278033912044631534762624256070,
      0.1216704729278033912044631534762624256070,
      0.1155056680537256013533444839067835598622,
      0.1155056680537256013533444839067835598622,
      0.1074442701159656347825773424466062227946,
      0.1074442701159656347825773424466062227946,
      0.0976186521041138882698806644642471544279,
      0.0976186521041138882698806644642471544279,
      0.0861901615319532759171852029837426671850,
      0.0861901615319532759171852029837426671850,
      0.0733464814110803057340336152531165181193,
      0.0733464814110803057340336152531165181193,
      0.0592985849154367807463677585001085845412,
      0.0592985849154367807463677585001085845412,
      0.0442774388174198061686027482113382288593,
      0.0442774388174198061686027482113382288593,
      0.0285313886289336631813078159518782864491,
      0.0285313886289336631813078159518782864491,
      0.0123412297999871995468056670700372915759,
      0.0123412297999871995468056670700372915759
    ],

    arcfn: function(t, derivativeFn) {
      var d = derivativeFn(t);
      var l = d.x*d.x + d.y*d.y;
      if(typeof d.z !== "undefined") {
        l += d.z*d.z;
      }
      return sqrt(l);
    },

    between: function(v, m, M) {
      return (m <= v && v <= M) || utils.approximately(v, m) || utils.approximately(v, M);
    },

    approximately: function(a,b,precision) {
      return abs(a-b) <= (precision || epsilon);
    },

    length: function(derivativeFn) {
      var z=0.5,sum=0,len=utils.Tvalues.length,i,t;
      for(i=0; i<len; i++) {
        t = z * utils.Tvalues[i] + z;
        sum += utils.Cvalues[i] * utils.arcfn(t,derivativeFn);
      }
      return z * sum;
    },

    map: function(v, ds,de, ts,te) {
      var d1 = de-ds, d2 = te-ts, v2 =  v-ds, r = v2/d1;
      return ts + d2*r;
    },

    lerp: function(r, v1, v2) {
      var ret = {
        x: v1.x + r*(v2.x-v1.x),
        y: v1.y + r*(v2.y-v1.y)
      };
      if(!!v1.z && !!v2.z) {
        ret.z =  v1.z + r*(v2.z-v1.z);
      }
      return ret;
    },

    pointToString: function(p) {
      var s = p.x+"/"+p.y;
      if(typeof p.z !== "undefined") {
        s += "/"+p.z;
      }
      return s;
    },

    pointsToString: function(points) {
      return "[" + points.map(utils.pointToString).join(", ") + "]";
    },

    copy: function(obj) {
      return JSON.parse(JSON.stringify(obj));
    },

    angle: function(o,v1,v2) {
      var dx1 = v1.x - o.x,
          dy1 = v1.y - o.y,
          dx2 = v2.x - o.x,
          dy2 = v2.y - o.y,
          cross = dx1*dy2 - dy1*dx2,
          dot = dx1*dx2 + dy1*dy2;
      return atan2(cross, dot);
    },

    // round as string, to avoid rounding errors
    round: function(v, d) {
      var s = '' + v;
      var pos = s.indexOf(".");
      return parseFloat(s.substring(0,pos+1+d));
    },

    dist: function(p1, p2) {
      var dx = p1.x - p2.x,
          dy = p1.y - p2.y;
      return sqrt(dx*dx+dy*dy);
    },

    closest: function(LUT, point) {
      var mdist = pow(2,63), mpos, d;
      LUT.forEach(function(p, idx) {
        d = utils.dist(point, p);
        if (d<mdist) {
          mdist = d;
          mpos = idx;
        }
      });
      return { mdist:mdist, mpos:mpos };
    },

    abcratio: function(t, n) {
      // see ratio(t) note on http://pomax.github.io/bezierinfo/#abc
      if (n!==2 && n!==3) {
        return false;
      }
      if (typeof t === "undefined") {
        t = 0.5;
      } else if (t===0 || t===1) {
        return t;
      }
      var bottom = pow(t,n) + pow(1-t,n), top = bottom - 1;
      return abs(top/bottom);
    },

    projectionratio: function(t, n) {
      // see u(t) note on http://pomax.github.io/bezierinfo/#abc
      if (n!==2 && n!==3) {
        return false;
      }
      if (typeof t === "undefined") {
        t = 0.5;
      } else if (t===0 || t===1) {
        return t;
      }
      var top = pow(1-t, n), bottom = pow(t,n) + top;
      return top/bottom;
    },

    lli8: function(x1,y1,x2,y2,x3,y3,x4,y4) {
      var nx=(x1*y2-y1*x2)*(x3-x4)-(x1-x2)*(x3*y4-y3*x4),
          ny=(x1*y2-y1*x2)*(y3-y4)-(y1-y2)*(x3*y4-y3*x4),
          d=(x1-x2)*(y3-y4)-(y1-y2)*(x3-x4);
      if(d==0) { return false; }
      return { x: nx/d, y: ny/d };
    },

    lli4: function(p1,p2,p3,p4) {
      var x1 = p1.x, y1 = p1.y,
          x2 = p2.x, y2 = p2.y,
          x3 = p3.x, y3 = p3.y,
          x4 = p4.x, y4 = p4.y;
      return utils.lli8(x1,y1,x2,y2,x3,y3,x4,y4);
    },

    lli: function(v1, v2) {
      return utils.lli4(v1,v1.c,v2,v2.c);
    },

    makeline: function(p1,p2) {
      var Bezier = require('./bezier');
      var x1 = p1.x, y1 = p1.y, x2 = p2.x, y2 = p2.y, dx = (x2-x1)/3, dy = (y2-y1)/3;
      return new Bezier(x1, y1, x1+dx, y1+dy, x1+2*dx, y1+2*dy, x2, y2);
    },

    findbbox: function(sections) {
      var mx=nMax,my=nMax,MX=nMin,MY=nMin;
      sections.forEach(function(s) {
        var bbox = s.bbox();
        if(mx > bbox.x.min) mx = bbox.x.min;
        if(my > bbox.y.min) my = bbox.y.min;
        if(MX < bbox.x.max) MX = bbox.x.max;
        if(MY < bbox.y.max) MY = bbox.y.max;
      });
      return {
        x: { min: mx, mid:(mx+MX)/2, max: MX, size:MX-mx },
        y: { min: my, mid:(my+MY)/2, max: MY, size:MY-my }
      }
    },

    shapeintersections: function(s1, bbox1, s2, bbox2, curveIntersectionThreshold) {
      if(!utils.bboxoverlap(bbox1, bbox2)) return [];
      var intersections = [];
      var a1 = [s1.startcap, s1.forward, s1.back, s1.endcap];
      var a2 = [s2.startcap, s2.forward, s2.back, s2.endcap];
      a1.forEach(function(l1) {
        if(l1.virtual) return;
        a2.forEach(function(l2) {
          if(l2.virtual) return;
          var iss = l1.intersects(l2, curveIntersectionThreshold);
          if(iss.length>0) {
            iss.c1 = l1;
            iss.c2 = l2;
            iss.s1 = s1;
            iss.s2 = s2;
            intersections.push(iss);
          }
        });
      });
      return intersections;
    },

    makeshape: function(forward, back, curveIntersectionThreshold) {
      var bpl = back.points.length;
      var fpl = forward.points.length;
      var start  = utils.makeline(back.points[bpl-1], forward.points[0]);
      var end    = utils.makeline(forward.points[fpl-1], back.points[0]);
      var shape  = {
        startcap: start,
        forward: forward,
        back: back,
        endcap: end,
        bbox: utils.findbbox([start, forward, back, end])
      };
      var self = utils;
      shape.intersections = function(s2) {
        return self.shapeintersections(shape,shape.bbox,s2,s2.bbox, curveIntersectionThreshold);
      };
      return shape;
    },

    getminmax: function(curve, d, list) {
      if(!list) return { min:0, max:0 };
      var min=nMax, max=nMin,t,c;
      if(list.indexOf(0)===-1) { list = [0].concat(list); }
      if(list.indexOf(1)===-1) { list.push(1); }
      for(var i=0,len=list.length; i<len; i++) {
        t = list[i];
        c = curve.get(t);
        if(c[d] < min) { min = c[d]; }
        if(c[d] > max) { max = c[d]; }
      }
      return { min:min, mid:(min+max)/2, max:max, size:max-min };
    },

    align: function(points, line) {
      var tx = line.p1.x,
          ty = line.p1.y,
          a = -atan2(line.p2.y-ty, line.p2.x-tx),
          d = function(v) {
            return {
              x: (v.x-tx)*cos(a) - (v.y-ty)*sin(a),
              y: (v.x-tx)*sin(a) + (v.y-ty)*cos(a)
            };
          };
      return points.map(d);
    },

    roots: function(points, line) {
      line = line || {p1:{x:0,y:0},p2:{x:1,y:0}};
      var order = points.length - 1;
      var p = utils.align(points, line);
      var reduce = function(t) { return 0<=t && t <=1; };

      if (order === 2) {
        var a = p[0].y,
            b = p[1].y,
            c = p[2].y,
            d = a - 2*b + c;
        if(d!==0) {
          var m1 = -sqrt(b*b-a*c),
              m2 = -a+b,
              v1 = -( m1+m2)/d,
              v2 = -(-m1+m2)/d;
          return [v1, v2].filter(reduce);
        }
        else if(b!==c && d===0) {
          return [ (2*b-c)/2*(b-c) ].filter(reduce);
        }
        return [];
      }

      // see http://www.trans4mind.com/personal_development/mathematics/polynomials/cubicAlgebra.htm
      var pa = p[0].y,
          pb = p[1].y,
          pc = p[2].y,
          pd = p[3].y,
          d = (-pa + 3*pb - 3*pc + pd),
          a = (3*pa - 6*pb + 3*pc) / d,
          b = (-3*pa + 3*pb) / d,
          c = pa / d,
          p = (3*b - a*a)/3,
          p3 = p/3,
          q = (2*a*a*a - 9*a*b + 27*c)/27,
          q2 = q/2,
          discriminant = q2*q2 + p3*p3*p3,
          u1,v1,x1,x2,x3;
       if (discriminant < 0) {
        var mp3 = -p/3,
            mp33 = mp3*mp3*mp3,
            r = sqrt( mp33 ),
            t = -q/(2*r),
            cosphi = t<-1 ? -1 : t>1 ? 1 : t,
            phi = acos(cosphi),
            crtr = crt(r),
            t1 = 2*crtr;
        x1 = t1 * cos(phi/3) - a/3;
        x2 = t1 * cos((phi+tau)/3) - a/3;
        x3 = t1 * cos((phi+2*tau)/3) - a/3;
        return [x1, x2, x3].filter(reduce);
      } else if(discriminant === 0) {
        u1 = q2 < 0 ? crt(-q2) : -crt(q2);
        x1 = 2*u1-a/3;
        x2 = -u1 - a/3;
        return [x1,x2].filter(reduce);
      } else {
        var sd = sqrt(discriminant);
        u1 = crt(-q2+sd);
        v1 = crt(q2+sd);
        return [u1-v1-a/3].filter(reduce);;
      }
    },

    droots: function(p) {
      // quadratic roots are easy
      if(p.length === 3) {
        var a = p[0],
            b = p[1],
            c = p[2],
            d = a - 2*b + c;
        if(d!==0) {
          var m1 = -sqrt(b*b-a*c),
              m2 = -a+b,
              v1 = -( m1+m2)/d,
              v2 = -(-m1+m2)/d;
          return [v1, v2];
        }
        else if(b!==c && d===0) {
          return [(2*b-c)/(2*(b-c))];
        }
        return [];
      }

      // linear roots are even easier
      if(p.length === 2) {
        var a = p[0], b = p[1];
        if(a!==b) {
          return [a/(a-b)];
        }
        return [];
      }
    },

    inflections: function(points) {
      if (points.length<4) return [];

      // FIXME: TODO: add in inflection abstraction for quartic+ curves?

      var p = utils.align(points, { p1: points[0], p2: points.slice(-1)[0] }),
          a = p[2].x * p[1].y,
          b = p[3].x * p[1].y,
          c = p[1].x * p[2].y,
          d = p[3].x * p[2].y,
          v1 = 18 * (-3*a + 2*b + 3*c - d),
          v2 = 18 * (3*a - b - 3*c),
          v3 = 18 * (c - a);

      if (utils.approximately(v1,0)){
        if(!utils.approximately(v2,0)){
          var t = -v3/v2;
          if (0 <= t && t <= 1)
             return [t];
        }
        return [];
      }

      var trm = v2*v2 - 4*v1*v3,
          sq = Math.sqrt(trm),
          d = 2 * v1;

      if (utils.approximately(d,0)) return [];

      return [(sq-v2)/d, -(v2+sq)/d].filter(function(r) {
        return (0 <= r && r <= 1);
      });
    },

    bboxoverlap: function(b1,b2) {
      var dims=['x','y'],len=dims.length,i,dim,l,t,d
      for(i=0; i<len; i++) {
        dim = dims[i];
        l = b1[dim].mid;
        t = b2[dim].mid;
        d = (b1[dim].size + b2[dim].size)/2;
        if(abs(l-t) >= d) return false;
      }
      return true;
    },

    expandbox: function(bbox, _bbox) {
      if(_bbox.x.min < bbox.x.min) { bbox.x.min = _bbox.x.min; }
      if(_bbox.y.min < bbox.y.min) { bbox.y.min = _bbox.y.min; }
      if(_bbox.z && _bbox.z.min < bbox.z.min) { bbox.z.min = _bbox.z.min; }
      if(_bbox.x.max > bbox.x.max) { bbox.x.max = _bbox.x.max; }
      if(_bbox.y.max > bbox.y.max) { bbox.y.max = _bbox.y.max; }
      if(_bbox.z && _bbox.z.max > bbox.z.max) { bbox.z.max = _bbox.z.max; }
      bbox.x.mid = (bbox.x.min + bbox.x.max)/2;
      bbox.y.mid = (bbox.y.min + bbox.y.max)/2;
      if(bbox.z) { bbox.z.mid = (bbox.z.min + bbox.z.max)/2; }
      bbox.x.size = bbox.x.max - bbox.x.min;
      bbox.y.size = bbox.y.max - bbox.y.min;
      if(bbox.z) { bbox.z.size = bbox.z.max - bbox.z.min; }
    },

    pairiteration: function(c1, c2, curveIntersectionThreshold) {
      var c1b = c1.bbox(),
          c2b = c2.bbox(),
          r = 100000,
          threshold = curveIntersectionThreshold || 0.5;
      if(c1b.x.size + c1b.y.size < threshold && c2b.x.size + c2b.y.size < threshold) {
        return [ ((r * (c1._t1+c1._t2)/2)|0)/r + "/" + ((r * (c2._t1+c2._t2)/2)|0)/r ];
      }
      var cc1 = c1.split(0.5),
          cc2 = c2.split(0.5),
          pairs = [
            {left: cc1.left, right: cc2.left },
            {left: cc1.left, right: cc2.right },
            {left: cc1.right, right: cc2.right },
            {left: cc1.right, right: cc2.left }];
      pairs = pairs.filter(function(pair) {
        return utils.bboxoverlap(pair.left.bbox(),pair.right.bbox());
      });
      var results = [];
      if(pairs.length === 0) return results;
      pairs.forEach(function(pair) {
        results = results.concat(
          utils.pairiteration(pair.left, pair.right, threshold)
        );
      })
      results = results.filter(function(v,i) {
        return results.indexOf(v) === i;
      });
      return results;
    },

    getccenter: function(p1,p2,p3) {
      var dx1 = (p2.x - p1.x),
          dy1 = (p2.y - p1.y),
          dx2 = (p3.x - p2.x),
          dy2 = (p3.y - p2.y);
      var dx1p = dx1 * cos(quart) - dy1 * sin(quart),
          dy1p = dx1 * sin(quart) + dy1 * cos(quart),
          dx2p = dx2 * cos(quart) - dy2 * sin(quart),
          dy2p = dx2 * sin(quart) + dy2 * cos(quart);
      // chord midpoints
      var mx1 = (p1.x + p2.x)/2,
          my1 = (p1.y + p2.y)/2,
          mx2 = (p2.x + p3.x)/2,
          my2 = (p2.y + p3.y)/2;
      // midpoint offsets
      var mx1n = mx1 + dx1p,
          my1n = my1 + dy1p,
          mx2n = mx2 + dx2p,
          my2n = my2 + dy2p;
      // intersection of these lines:
      var arc = utils.lli8(mx1,my1,mx1n,my1n, mx2,my2,mx2n,my2n),
          r = utils.dist(arc,p1),
          // arc start/end values, over mid point:
          s = atan2(p1.y - arc.y, p1.x - arc.x),
          m = atan2(p2.y - arc.y, p2.x - arc.x),
          e = atan2(p3.y - arc.y, p3.x - arc.x),
          _;
      // determine arc direction (cw/ccw correction)
      if (s<e) {
        // if s<m<e, arc(s, e)
        // if m<s<e, arc(e, s + tau)
        // if s<e<m, arc(e, s + tau)
        if (s>m || m>e) { s += tau; }
        if (s>e) { _=e; e=s; s=_; }
      } else {
        // if e<m<s, arc(e, s)
        // if m<e<s, arc(s, e + tau)
        // if e<s<m, arc(s, e + tau)
        if (e<m && m<s) { _=e; e=s; s=_; } else { e += tau; }
      }
      // assign and done.
      arc.s = s;
      arc.e = e;
      arc.r = r;
      return arc;
    }
  };

  module.exports = utils;
}());

},{"./bezier":2}],5:[function(require,module,exports){
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

},{}],6:[function(require,module,exports){
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
},{}],7:[function(require,module,exports){
'use strict';

Object.defineProperty(exports, "__esModule", {
	value: true
});

var _fitCurve = require('fit-curve');

var _fitCurve2 = _interopRequireDefault(_fitCurve);

var _LevelCurve = require('../model/LevelCurve');

var _LevelCurve2 = _interopRequireDefault(_LevelCurve);

var _UIManagement = require('../model/UIManagement');

var UI = _interopRequireWildcard(_UIManagement);

function _interopRequireWildcard(obj) { if (obj && obj.__esModule) { return obj; } else { var newObj = {}; if (obj != null) { for (var key in obj) { if (Object.prototype.hasOwnProperty.call(obj, key)) newObj[key] = obj[key]; } } newObj.default = obj; return newObj; } }

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
		if (smoothBizer.length == 0) {
			clearRawData();
			return;
		}
		var pathString = fittedCurveToPathString(smoothBizer);

		drawOnPannel(pannel, pathString);

		var lvCurve = new _LevelCurve2.default(smoothBizer, 1, UI.state.levelCurve);
		lvCurve.drawOn(pannel);

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

},{"../model/LevelCurve":12,"../model/UIManagement":14,"fit-curve":5}],8:[function(require,module,exports){
'use strict';

var _typeof2 = typeof Symbol === "function" && typeof Symbol.iterator === "symbol" ? function (obj) { return typeof obj; } : function (obj) { return obj && typeof Symbol === "function" && obj.constructor === Symbol && obj !== Symbol.prototype ? "symbol" : typeof obj; };

(function webpackUniversalModuleDefinition(root, factory) {
	if ((typeof exports === 'undefined' ? 'undefined' : _typeof2(exports)) === 'object' && (typeof module === 'undefined' ? 'undefined' : _typeof2(module)) === 'object') module.exports = factory();else if (typeof define === 'function' && define.amd) define([], factory);else if ((typeof exports === 'undefined' ? 'undefined' : _typeof2(exports)) === 'object') exports["dat"] = factory();else root["dat"] = factory();
})(undefined, function () {
	return (/******/function (modules) {
			// webpackBootstrap
			/******/ // The module cache
			/******/var installedModules = {};
			/******/
			/******/ // The require function
			/******/function __webpack_require__(moduleId) {
				/******/
				/******/ // Check if module is in cache
				/******/if (installedModules[moduleId])
					/******/return installedModules[moduleId].exports;
				/******/
				/******/ // Create a new module (and put it into the cache)
				/******/var module = installedModules[moduleId] = {
					/******/exports: {},
					/******/id: moduleId,
					/******/loaded: false
					/******/ };
				/******/
				/******/ // Execute the module function
				/******/modules[moduleId].call(module.exports, module, module.exports, __webpack_require__);
				/******/
				/******/ // Flag the module as loaded
				/******/module.loaded = true;
				/******/
				/******/ // Return the exports of the module
				/******/return module.exports;
				/******/
			}
			/******/
			/******/
			/******/ // expose the modules object (__webpack_modules__)
			/******/__webpack_require__.m = modules;
			/******/
			/******/ // expose the module cache
			/******/__webpack_require__.c = installedModules;
			/******/
			/******/ // __webpack_public_path__
			/******/__webpack_require__.p = "";
			/******/
			/******/ // Load entry module and return exports
			/******/return __webpack_require__(0);
			/******/
		}(
		/************************************************************************/
		/******/[
		/* 0 */
		/***/function (module, exports, __webpack_require__) {

			'use strict';

			var _index = __webpack_require__(1);

			var _index2 = _interopRequireDefault(_index);

			function _interopRequireDefault(obj) {
				return obj && obj.__esModule ? obj : { default: obj };
			}

			module.exports = _index2.default; /**
                                      * dat-gui JavaScript Controller Library
                                      * http://code.google.com/p/dat-gui
                                      *
                                      * Copyright 2011 Data Arts Team, Google Creative Lab
                                      *
                                      * Licensed under the Apache License, Version 2.0 (the "License");
                                      * you may not use this file except in compliance with the License.
                                      * You may obtain a copy of the License at
                                      *
                                      * http://www.apache.org/licenses/LICENSE-2.0
                                      */

			/***/
		},
		/* 1 */
		/***/function (module, exports, __webpack_require__) {

			'use strict';

			exports.__esModule = true;

			var _Color = __webpack_require__(2);

			var _Color2 = _interopRequireDefault(_Color);

			var _math = __webpack_require__(6);

			var _math2 = _interopRequireDefault(_math);

			var _interpret = __webpack_require__(3);

			var _interpret2 = _interopRequireDefault(_interpret);

			var _Controller = __webpack_require__(7);

			var _Controller2 = _interopRequireDefault(_Controller);

			var _BooleanController = __webpack_require__(8);

			var _BooleanController2 = _interopRequireDefault(_BooleanController);

			var _OptionController = __webpack_require__(10);

			var _OptionController2 = _interopRequireDefault(_OptionController);

			var _StringController = __webpack_require__(11);

			var _StringController2 = _interopRequireDefault(_StringController);

			var _NumberController = __webpack_require__(12);

			var _NumberController2 = _interopRequireDefault(_NumberController);

			var _NumberControllerBox = __webpack_require__(13);

			var _NumberControllerBox2 = _interopRequireDefault(_NumberControllerBox);

			var _NumberControllerSlider = __webpack_require__(14);

			var _NumberControllerSlider2 = _interopRequireDefault(_NumberControllerSlider);

			var _FunctionController = __webpack_require__(15);

			var _FunctionController2 = _interopRequireDefault(_FunctionController);

			var _ColorController = __webpack_require__(16);

			var _ColorController2 = _interopRequireDefault(_ColorController);

			var _dom = __webpack_require__(9);

			var _dom2 = _interopRequireDefault(_dom);

			var _GUI = __webpack_require__(17);

			var _GUI2 = _interopRequireDefault(_GUI);

			function _interopRequireDefault(obj) {
				return obj && obj.__esModule ? obj : { default: obj };
			}

			/**
    * dat-gui JavaScript Controller Library
    * http://code.google.com/p/dat-gui
    *
    * Copyright 2011 Data Arts Team, Google Creative Lab
    *
    * Licensed under the Apache License, Version 2.0 (the "License");
    * you may not use this file except in compliance with the License.
    * You may obtain a copy of the License at
    *
    * http://www.apache.org/licenses/LICENSE-2.0
    */

			exports.default = {
				color: {
					Color: _Color2.default,
					math: _math2.default,
					interpret: _interpret2.default
				},

				controllers: {
					Controller: _Controller2.default,
					BooleanController: _BooleanController2.default,
					OptionController: _OptionController2.default,
					StringController: _StringController2.default,
					NumberController: _NumberController2.default,
					NumberControllerBox: _NumberControllerBox2.default,
					NumberControllerSlider: _NumberControllerSlider2.default,
					FunctionController: _FunctionController2.default,
					ColorController: _ColorController2.default
				},

				dom: {
					dom: _dom2.default
				},

				gui: {
					GUI: _GUI2.default
				},

				GUI: _GUI2.default
			};

			/***/
		},
		/* 2 */
		/***/function (module, exports, __webpack_require__) {

			'use strict';

			exports.__esModule = true;

			var _interpret = __webpack_require__(3);

			var _interpret2 = _interopRequireDefault(_interpret);

			var _math = __webpack_require__(6);

			var _math2 = _interopRequireDefault(_math);

			var _toString = __webpack_require__(4);

			var _toString2 = _interopRequireDefault(_toString);

			var _common = __webpack_require__(5);

			var _common2 = _interopRequireDefault(_common);

			function _interopRequireDefault(obj) {
				return obj && obj.__esModule ? obj : { default: obj };
			}

			function _classCallCheck(instance, Constructor) {
				if (!(instance instanceof Constructor)) {
					throw new TypeError("Cannot call a class as a function");
				}
			} /**
      * dat-gui JavaScript Controller Library
      * http://code.google.com/p/dat-gui
      *
      * Copyright 2011 Data Arts Team, Google Creative Lab
      *
      * Licensed under the Apache License, Version 2.0 (the "License");
      * you may not use this file except in compliance with the License.
      * You may obtain a copy of the License at
      *
      * http://www.apache.org/licenses/LICENSE-2.0
      */

			var Color = function () {
				function Color() {
					_classCallCheck(this, Color);

					this.__state = _interpret2.default.apply(this, arguments);

					if (this.__state === false) {
						throw new Error('Failed to interpret color arguments');
					}

					this.__state.a = this.__state.a || 1;
				}

				Color.prototype.toString = function toString() {
					return (0, _toString2.default)(this);
				};

				Color.prototype.toHexString = function toHexString() {
					return (0, _toString2.default)(this, true);
				};

				Color.prototype.toOriginal = function toOriginal() {
					return this.__state.conversion.write(this);
				};

				return Color;
			}();

			function defineRGBComponent(target, component, componentHexIndex) {
				Object.defineProperty(target, component, {
					get: function get() {
						if (this.__state.space === 'RGB') {
							return this.__state[component];
						}

						Color.recalculateRGB(this, component, componentHexIndex);

						return this.__state[component];
					},

					set: function set(v) {
						if (this.__state.space !== 'RGB') {
							Color.recalculateRGB(this, component, componentHexIndex);
							this.__state.space = 'RGB';
						}

						this.__state[component] = v;
					}
				});
			}

			function defineHSVComponent(target, component) {
				Object.defineProperty(target, component, {
					get: function get() {
						if (this.__state.space === 'HSV') {
							return this.__state[component];
						}

						Color.recalculateHSV(this);

						return this.__state[component];
					},

					set: function set(v) {
						if (this.__state.space !== 'HSV') {
							Color.recalculateHSV(this);
							this.__state.space = 'HSV';
						}

						this.__state[component] = v;
					}
				});
			}

			Color.recalculateRGB = function (color, component, componentHexIndex) {
				if (color.__state.space === 'HEX') {
					color.__state[component] = _math2.default.component_from_hex(color.__state.hex, componentHexIndex);
				} else if (color.__state.space === 'HSV') {
					_common2.default.extend(color.__state, _math2.default.hsv_to_rgb(color.__state.h, color.__state.s, color.__state.v));
				} else {
					throw new Error('Corrupted color state');
				}
			};

			Color.recalculateHSV = function (color) {
				var result = _math2.default.rgb_to_hsv(color.r, color.g, color.b);

				_common2.default.extend(color.__state, {
					s: result.s,
					v: result.v
				});

				if (!_common2.default.isNaN(result.h)) {
					color.__state.h = result.h;
				} else if (_common2.default.isUndefined(color.__state.h)) {
					color.__state.h = 0;
				}
			};

			Color.COMPONENTS = ['r', 'g', 'b', 'h', 's', 'v', 'hex', 'a'];

			defineRGBComponent(Color.prototype, 'r', 2);
			defineRGBComponent(Color.prototype, 'g', 1);
			defineRGBComponent(Color.prototype, 'b', 0);

			defineHSVComponent(Color.prototype, 'h');
			defineHSVComponent(Color.prototype, 's');
			defineHSVComponent(Color.prototype, 'v');

			Object.defineProperty(Color.prototype, 'a', {
				get: function get() {
					return this.__state.a;
				},

				set: function set(v) {
					this.__state.a = v;
				}
			});

			Object.defineProperty(Color.prototype, 'hex', {
				get: function get() {
					if (!this.__state.space !== 'HEX') {
						this.__state.hex = _math2.default.rgb_to_hex(this.r, this.g, this.b);
					}

					return this.__state.hex;
				},

				set: function set(v) {
					this.__state.space = 'HEX';
					this.__state.hex = v;
				}
			});

			exports.default = Color;

			/***/
		},
		/* 3 */
		/***/function (module, exports, __webpack_require__) {

			'use strict';

			exports.__esModule = true;

			var _toString = __webpack_require__(4);

			var _toString2 = _interopRequireDefault(_toString);

			var _common = __webpack_require__(5);

			var _common2 = _interopRequireDefault(_common);

			function _interopRequireDefault(obj) {
				return obj && obj.__esModule ? obj : { default: obj };
			}

			/**
    * dat-gui JavaScript Controller Library
    * http://code.google.com/p/dat-gui
    *
    * Copyright 2011 Data Arts Team, Google Creative Lab
    *
    * Licensed under the Apache License, Version 2.0 (the "License");
    * you may not use this file except in compliance with the License.
    * You may obtain a copy of the License at
    *
    * http://www.apache.org/licenses/LICENSE-2.0
    */

			var INTERPRETATIONS = [
			// Strings
			{
				litmus: _common2.default.isString,
				conversions: {
					THREE_CHAR_HEX: {
						read: function read(original) {
							var test = original.match(/^#([A-F0-9])([A-F0-9])([A-F0-9])$/i);
							if (test === null) {
								return false;
							}

							return {
								space: 'HEX',
								hex: parseInt('0x' + test[1].toString() + test[1].toString() + test[2].toString() + test[2].toString() + test[3].toString() + test[3].toString(), 0)
							};
						},

						write: _toString2.default
					},

					SIX_CHAR_HEX: {
						read: function read(original) {
							var test = original.match(/^#([A-F0-9]{6})$/i);
							if (test === null) {
								return false;
							}

							return {
								space: 'HEX',
								hex: parseInt('0x' + test[1].toString(), 0)
							};
						},

						write: _toString2.default
					},

					CSS_RGB: {
						read: function read(original) {
							var test = original.match(/^rgb\(\s*(.+)\s*,\s*(.+)\s*,\s*(.+)\s*\)/);
							if (test === null) {
								return false;
							}

							return {
								space: 'RGB',
								r: parseFloat(test[1]),
								g: parseFloat(test[2]),
								b: parseFloat(test[3])
							};
						},

						write: _toString2.default
					},

					CSS_RGBA: {
						read: function read(original) {
							var test = original.match(/^rgba\(\s*(.+)\s*,\s*(.+)\s*,\s*(.+)\s*,\s*(.+)\s*\)/);
							if (test === null) {
								return false;
							}

							return {
								space: 'RGB',
								r: parseFloat(test[1]),
								g: parseFloat(test[2]),
								b: parseFloat(test[3]),
								a: parseFloat(test[4])
							};
						},

						write: _toString2.default
					}
				}
			},

			// Numbers
			{
				litmus: _common2.default.isNumber,

				conversions: {

					HEX: {
						read: function read(original) {
							return {
								space: 'HEX',
								hex: original,
								conversionName: 'HEX'
							};
						},

						write: function write(color) {
							return color.hex;
						}
					}

				}

			},

			// Arrays
			{
				litmus: _common2.default.isArray,
				conversions: {
					RGB_ARRAY: {
						read: function read(original) {
							if (original.length !== 3) {
								return false;
							}

							return {
								space: 'RGB',
								r: original[0],
								g: original[1],
								b: original[2]
							};
						},

						write: function write(color) {
							return [color.r, color.g, color.b];
						}
					},

					RGBA_ARRAY: {
						read: function read(original) {
							if (original.length !== 4) return false;
							return {
								space: 'RGB',
								r: original[0],
								g: original[1],
								b: original[2],
								a: original[3]
							};
						},

						write: function write(color) {
							return [color.r, color.g, color.b, color.a];
						}
					}
				}
			},

			// Objects
			{
				litmus: _common2.default.isObject,
				conversions: {

					RGBA_OBJ: {
						read: function read(original) {
							if (_common2.default.isNumber(original.r) && _common2.default.isNumber(original.g) && _common2.default.isNumber(original.b) && _common2.default.isNumber(original.a)) {
								return {
									space: 'RGB',
									r: original.r,
									g: original.g,
									b: original.b,
									a: original.a
								};
							}
							return false;
						},

						write: function write(color) {
							return {
								r: color.r,
								g: color.g,
								b: color.b,
								a: color.a
							};
						}
					},

					RGB_OBJ: {
						read: function read(original) {
							if (_common2.default.isNumber(original.r) && _common2.default.isNumber(original.g) && _common2.default.isNumber(original.b)) {
								return {
									space: 'RGB',
									r: original.r,
									g: original.g,
									b: original.b
								};
							}
							return false;
						},

						write: function write(color) {
							return {
								r: color.r,
								g: color.g,
								b: color.b
							};
						}
					},

					HSVA_OBJ: {
						read: function read(original) {
							if (_common2.default.isNumber(original.h) && _common2.default.isNumber(original.s) && _common2.default.isNumber(original.v) && _common2.default.isNumber(original.a)) {
								return {
									space: 'HSV',
									h: original.h,
									s: original.s,
									v: original.v,
									a: original.a
								};
							}
							return false;
						},

						write: function write(color) {
							return {
								h: color.h,
								s: color.s,
								v: color.v,
								a: color.a
							};
						}
					},

					HSV_OBJ: {
						read: function read(original) {
							if (_common2.default.isNumber(original.h) && _common2.default.isNumber(original.s) && _common2.default.isNumber(original.v)) {
								return {
									space: 'HSV',
									h: original.h,
									s: original.s,
									v: original.v
								};
							}
							return false;
						},

						write: function write(color) {
							return {
								h: color.h,
								s: color.s,
								v: color.v
							};
						}
					}
				}
			}];

			var result = void 0;
			var toReturn = void 0;

			var interpret = function interpret() {
				toReturn = false;

				var original = arguments.length > 1 ? _common2.default.toArray(arguments) : arguments[0];
				_common2.default.each(INTERPRETATIONS, function (family) {
					if (family.litmus(original)) {
						_common2.default.each(family.conversions, function (conversion, conversionName) {
							result = conversion.read(original);

							if (toReturn === false && result !== false) {
								toReturn = result;
								result.conversionName = conversionName;
								result.conversion = conversion;
								return _common2.default.BREAK;
							}
						});

						return _common2.default.BREAK;
					}
				});

				return toReturn;
			};

			exports.default = interpret;

			/***/
		},
		/* 4 */
		/***/function (module, exports) {

			'use strict';

			exports.__esModule = true;

			exports.default = function (color, forceCSSHex) {
				var colorFormat = color.__state.conversionName.toString();

				var r = Math.round(color.r);
				var g = Math.round(color.g);
				var b = Math.round(color.b);
				var a = color.a;
				var h = Math.round(color.h);
				var s = color.s.toFixed(1);
				var v = color.v.toFixed(1);

				if (forceCSSHex || colorFormat === 'THREE_CHAR_HEX' || colorFormat === 'SIX_CHAR_HEX') {
					var str = color.hex.toString(16);
					while (str.length < 6) {
						str = '0' + str;
					}
					return '#' + str;
				} else if (colorFormat === 'CSS_RGB') {
					return 'rgb(' + r + ',' + g + ',' + b + ')';
				} else if (colorFormat === 'CSS_RGBA') {
					return 'rgba(' + r + ',' + g + ',' + b + ',' + a + ')';
				} else if (colorFormat === 'HEX') {
					return '0x' + color.hex.toString(16);
				} else if (colorFormat === 'RGB_ARRAY') {
					return '[' + r + ',' + g + ',' + b + ']';
				} else if (colorFormat === 'RGBA_ARRAY') {
					return '[' + r + ',' + g + ',' + b + ',' + a + ']';
				} else if (colorFormat === 'RGB_OBJ') {
					return '{r:' + r + ',g:' + g + ',b:' + b + '}';
				} else if (colorFormat === 'RGBA_OBJ') {
					return '{r:' + r + ',g:' + g + ',b:' + b + ',a:' + a + '}';
				} else if (colorFormat === 'HSV_OBJ') {
					return '{h:' + h + ',s:' + s + ',v:' + v + '}';
				} else if (colorFormat === 'HSVA_OBJ') {
					return '{h:' + h + ',s:' + s + ',v:' + v + ',a:' + a + '}';
				}

				return 'unknown format';
			};

			/***/
		},
		/* 5 */
		/***/function (module, exports) {

			'use strict';

			exports.__esModule = true;
			/**
    * dat-gui JavaScript Controller Library
    * http://code.google.com/p/dat-gui
    *
    * Copyright 2011 Data Arts Team, Google Creative Lab
    *
    * Licensed under the Apache License, Version 2.0 (the "License");
    * you may not use this file except in compliance with the License.
    * You may obtain a copy of the License at
    *
    * http://www.apache.org/licenses/LICENSE-2.0
    */

			var ARR_EACH = Array.prototype.forEach;
			var ARR_SLICE = Array.prototype.slice;

			/**
    * Band-aid methods for things that should be a lot easier in JavaScript.
    * Implementation and structure inspired by underscore.js
    * http://documentcloud.github.com/underscore/
    */

			var Common = {
				BREAK: {},

				extend: function extend(target) {
					this.each(ARR_SLICE.call(arguments, 1), function (obj) {
						var keys = this.isObject(obj) ? Object.keys(obj) : [];
						keys.forEach(function (key) {
							if (!this.isUndefined(obj[key])) {
								target[key] = obj[key];
							}
						}.bind(this));
					}, this);

					return target;
				},

				defaults: function defaults(target) {
					this.each(ARR_SLICE.call(arguments, 1), function (obj) {
						var keys = this.isObject(obj) ? Object.keys(obj) : [];
						keys.forEach(function (key) {
							if (this.isUndefined(target[key])) {
								target[key] = obj[key];
							}
						}.bind(this));
					}, this);

					return target;
				},

				compose: function compose() {
					var toCall = ARR_SLICE.call(arguments);
					return function () {
						var args = ARR_SLICE.call(arguments);
						for (var i = toCall.length - 1; i >= 0; i--) {
							args = [toCall[i].apply(this, args)];
						}
						return args[0];
					};
				},

				each: function each(obj, itr, scope) {
					if (!obj) {
						return;
					}

					if (ARR_EACH && obj.forEach && obj.forEach === ARR_EACH) {
						obj.forEach(itr, scope);
					} else if (obj.length === obj.length + 0) {
						// Is number but not NaN
						var key = void 0;
						var l = void 0;
						for (key = 0, l = obj.length; key < l; key++) {
							if (key in obj && itr.call(scope, obj[key], key) === this.BREAK) {
								return;
							}
						}
					} else {
						for (var _key in obj) {
							if (itr.call(scope, obj[_key], _key) === this.BREAK) {
								return;
							}
						}
					}
				},

				defer: function defer(fnc) {
					setTimeout(fnc, 0);
				},

				// call the function immediately, but wait until threshold passes to allow it to be called again
				debounce: function debounce(func, threshold) {
					var timeout = void 0;

					return function () {
						var obj = this;
						var args = arguments;
						function delayed() {
							timeout = null;
						}

						var allowCall = !timeout;

						clearTimeout(timeout);
						timeout = setTimeout(delayed, threshold);

						if (allowCall) {
							func.apply(obj, args);
						}
					};
				},

				toArray: function toArray(obj) {
					if (obj.toArray) return obj.toArray();
					return ARR_SLICE.call(obj);
				},

				isUndefined: function isUndefined(obj) {
					return obj === undefined;
				},

				isNull: function isNull(obj) {
					return obj === null;
				},

				isNaN: function (_isNaN) {
					function isNaN(_x) {
						return _isNaN.apply(this, arguments);
					}

					isNaN.toString = function () {
						return _isNaN.toString();
					};

					return isNaN;
				}(function (obj) {
					return isNaN(obj);
				}),

				isArray: Array.isArray || function (obj) {
					return obj.constructor === Array;
				},

				isObject: function isObject(obj) {
					return obj === Object(obj);
				},

				isNumber: function isNumber(obj) {
					return obj === obj + 0;
				},

				isString: function isString(obj) {
					return obj === obj + '';
				},

				isBoolean: function isBoolean(obj) {
					return obj === false || obj === true;
				},

				isFunction: function isFunction(obj) {
					return Object.prototype.toString.call(obj) === '[object Function]';
				}

			};

			exports.default = Common;

			/***/
		},
		/* 6 */
		/***/function (module, exports) {

			"use strict";

			exports.__esModule = true;
			/**
    * dat-gui JavaScript Controller Library
    * http://code.google.com/p/dat-gui
    *
    * Copyright 2011 Data Arts Team, Google Creative Lab
    *
    * Licensed under the Apache License, Version 2.0 (the "License");
    * you may not use this file except in compliance with the License.
    * You may obtain a copy of the License at
    *
    * http://www.apache.org/licenses/LICENSE-2.0
    */

			var tmpComponent = void 0;

			var ColorMath = {
				hsv_to_rgb: function hsv_to_rgb(h, s, v) {
					var hi = Math.floor(h / 60) % 6;

					var f = h / 60 - Math.floor(h / 60);
					var p = v * (1.0 - s);
					var q = v * (1.0 - f * s);
					var t = v * (1.0 - (1.0 - f) * s);

					var c = [[v, t, p], [q, v, p], [p, v, t], [p, q, v], [t, p, v], [v, p, q]][hi];

					return {
						r: c[0] * 255,
						g: c[1] * 255,
						b: c[2] * 255
					};
				},

				rgb_to_hsv: function rgb_to_hsv(r, g, b) {
					var min = Math.min(r, g, b);
					var max = Math.max(r, g, b);
					var delta = max - min;
					var h = void 0;
					var s = void 0;

					if (max !== 0) {
						s = delta / max;
					} else {
						return {
							h: NaN,
							s: 0,
							v: 0
						};
					}

					if (r === max) {
						h = (g - b) / delta;
					} else if (g === max) {
						h = 2 + (b - r) / delta;
					} else {
						h = 4 + (r - g) / delta;
					}
					h /= 6;
					if (h < 0) {
						h += 1;
					}

					return {
						h: h * 360,
						s: s,
						v: max / 255
					};
				},

				rgb_to_hex: function rgb_to_hex(r, g, b) {
					var hex = this.hex_with_component(0, 2, r);
					hex = this.hex_with_component(hex, 1, g);
					hex = this.hex_with_component(hex, 0, b);
					return hex;
				},

				component_from_hex: function component_from_hex(hex, componentIndex) {
					return hex >> componentIndex * 8 & 0xFF;
				},

				hex_with_component: function hex_with_component(hex, componentIndex, value) {
					return value << (tmpComponent = componentIndex * 8) | hex & ~(0xFF << tmpComponent);
				}
			};

			exports.default = ColorMath;

			/***/
		},
		/* 7 */
		/***/function (module, exports) {

			'use strict';

			exports.__esModule = true;

			function _classCallCheck(instance, Constructor) {
				if (!(instance instanceof Constructor)) {
					throw new TypeError("Cannot call a class as a function");
				}
			}

			/**
    * dat-gui JavaScript Controller Library
    * http://code.google.com/p/dat-gui
    *
    * Copyright 2011 Data Arts Team, Google Creative Lab
    *
    * Licensed under the Apache License, Version 2.0 (the "License");
    * you may not use this file except in compliance with the License.
    * You may obtain a copy of the License at
    *
    * http://www.apache.org/licenses/LICENSE-2.0
    */

			/**
    * @class An "abstract" class that represents a given property of an object.
    *
    * @param {Object} object The object to be manipulated
    * @param {string} property The name of the property to be manipulated
    *
    * @member dat.controllers
    */
			var Controller = function () {
				function Controller(object, property) {
					_classCallCheck(this, Controller);

					this.initialValue = object[property];

					/**
      * Those who extend this class will put their DOM elements in here.
      * @type {DOMElement}
      */
					this.domElement = document.createElement('div');

					/**
      * The object to manipulate
      * @type {Object}
      */
					this.object = object;

					/**
      * The name of the property to manipulate
      * @type {String}
      */
					this.property = property;

					/**
      * The function to be called on change.
      * @type {Function}
      * @ignore
      */
					this.__onChange = undefined;

					/**
      * The function to be called on finishing change.
      * @type {Function}
      * @ignore
      */
					this.__onFinishChange = undefined;
				}

				/**
     * Specify that a function fire every time someone changes the value with
     * this Controller.
     *
     * @param {Function} fnc This function will be called whenever the value
     * is modified via this Controller.
     * @returns {Controller} this
     */

				Controller.prototype.onChange = function onChange(fnc) {
					this.__onChange = fnc;
					return this;
				};

				/**
     * Specify that a function fire every time someone "finishes" changing
     * the value wih this Controller. Useful for values that change
     * incrementally like numbers or strings.
     *
     * @param {Function} fnc This function will be called whenever
     * someone "finishes" changing the value via this Controller.
     * @returns {Controller} this
     */

				Controller.prototype.onFinishChange = function onFinishChange(fnc) {
					this.__onFinishChange = fnc;
					return this;
				};

				/**
     * Change the value of <code>object[property]</code>
     *
     * @param {Object} newValue The new value of <code>object[property]</code>
     */

				Controller.prototype.setValue = function setValue(newValue) {
					this.object[this.property] = newValue;
					if (this.__onChange) {
						this.__onChange.call(this, newValue);
					}

					this.updateDisplay();
					return this;
				};

				/**
     * Gets the value of <code>object[property]</code>
     *
     * @returns {Object} The current value of <code>object[property]</code>
     */

				Controller.prototype.getValue = function getValue() {
					return this.object[this.property];
				};

				/**
     * Refreshes the visual display of a Controller in order to keep sync
     * with the object's current value.
     * @returns {Controller} this
     */

				Controller.prototype.updateDisplay = function updateDisplay() {
					return this;
				};

				/**
     * @returns {Boolean} true if the value has deviated from initialValue
     */

				Controller.prototype.isModified = function isModified() {
					return this.initialValue !== this.getValue();
				};

				return Controller;
			}();

			exports.default = Controller;

			/***/
		},
		/* 8 */
		/***/function (module, exports, __webpack_require__) {

			'use strict';

			exports.__esModule = true;

			var _Controller2 = __webpack_require__(7);

			var _Controller3 = _interopRequireDefault(_Controller2);

			var _dom = __webpack_require__(9);

			var _dom2 = _interopRequireDefault(_dom);

			function _interopRequireDefault(obj) {
				return obj && obj.__esModule ? obj : { default: obj };
			}

			function _classCallCheck(instance, Constructor) {
				if (!(instance instanceof Constructor)) {
					throw new TypeError("Cannot call a class as a function");
				}
			}

			function _possibleConstructorReturn(self, call) {
				if (!self) {
					throw new ReferenceError("this hasn't been initialised - super() hasn't been called");
				}return call && ((typeof call === 'undefined' ? 'undefined' : _typeof2(call)) === "object" || typeof call === "function") ? call : self;
			}

			function _inherits(subClass, superClass) {
				if (typeof superClass !== "function" && superClass !== null) {
					throw new TypeError("Super expression must either be null or a function, not " + (typeof superClass === 'undefined' ? 'undefined' : _typeof2(superClass)));
				}subClass.prototype = Object.create(superClass && superClass.prototype, { constructor: { value: subClass, enumerable: false, writable: true, configurable: true } });if (superClass) Object.setPrototypeOf ? Object.setPrototypeOf(subClass, superClass) : subClass.__proto__ = superClass;
			} /**
      * dat-gui JavaScript Controller Library
      * http://code.google.com/p/dat-gui
      *
      * Copyright 2011 Data Arts Team, Google Creative Lab
      *
      * Licensed under the Apache License, Version 2.0 (the "License");
      * you may not use this file except in compliance with the License.
      * You may obtain a copy of the License at
      *
      * http://www.apache.org/licenses/LICENSE-2.0
      */

			/**
    * @class Provides a checkbox input to alter the boolean property of an object.
    * @extends dat.controllers.Controller
    *
    * @param {Object} object The object to be manipulated
    * @param {string} property The name of the property to be manipulated
    *
    * @member dat.controllers
    */
			var BooleanController = function (_Controller) {
				_inherits(BooleanController, _Controller);

				function BooleanController(object, property) {
					_classCallCheck(this, BooleanController);

					var _this2 = _possibleConstructorReturn(this, _Controller.call(this, object, property));

					var _this = _this2;
					_this2.__prev = _this2.getValue();

					_this2.__checkbox = document.createElement('input');
					_this2.__checkbox.setAttribute('type', 'checkbox');

					function onChange() {
						_this.setValue(!_this.__prev);
					}

					_dom2.default.bind(_this2.__checkbox, 'change', onChange, false);

					_this2.domElement.appendChild(_this2.__checkbox);

					// Match original value
					_this2.updateDisplay();
					return _this2;
				}

				BooleanController.prototype.setValue = function setValue(v) {
					var toReturn = _Controller.prototype.setValue.call(this, v);
					if (this.__onFinishChange) {
						this.__onFinishChange.call(this, this.getValue());
					}
					this.__prev = this.getValue();
					return toReturn;
				};

				BooleanController.prototype.updateDisplay = function updateDisplay() {
					if (this.getValue() === true) {
						this.__checkbox.setAttribute('checked', 'checked');
						this.__checkbox.checked = true;
					} else {
						this.__checkbox.checked = false;
					}

					return _Controller.prototype.updateDisplay.call(this);
				};

				return BooleanController;
			}(_Controller3.default);

			exports.default = BooleanController;

			/***/
		},
		/* 9 */
		/***/function (module, exports, __webpack_require__) {

			'use strict';

			exports.__esModule = true;

			var _common = __webpack_require__(5);

			var _common2 = _interopRequireDefault(_common);

			function _interopRequireDefault(obj) {
				return obj && obj.__esModule ? obj : { default: obj };
			}

			var EVENT_MAP = {
				HTMLEvents: ['change'],
				MouseEvents: ['click', 'mousemove', 'mousedown', 'mouseup', 'mouseover'],
				KeyboardEvents: ['keydown']
			}; /**
       * dat-gui JavaScript Controller Library
       * http://code.google.com/p/dat-gui
       *
       * Copyright 2011 Data Arts Team, Google Creative Lab
       *
       * Licensed under the Apache License, Version 2.0 (the "License");
       * you may not use this file except in compliance with the License.
       * You may obtain a copy of the License at
       *
       * http://www.apache.org/licenses/LICENSE-2.0
       */

			var EVENT_MAP_INV = {};
			_common2.default.each(EVENT_MAP, function (v, k) {
				_common2.default.each(v, function (e) {
					EVENT_MAP_INV[e] = k;
				});
			});

			var CSS_VALUE_PIXELS = /(\d+(\.\d+)?)px/;

			function cssValueToPixels(val) {
				if (val === '0' || _common2.default.isUndefined(val)) {
					return 0;
				}

				var match = val.match(CSS_VALUE_PIXELS);

				if (!_common2.default.isNull(match)) {
					return parseFloat(match[1]);
				}

				// TODO ...ems? %?

				return 0;
			}

			/**
    * @namespace
    * @member dat.dom
    */
			var dom = {

				/**
     *
     * @param elem
     * @param selectable
     */
				makeSelectable: function makeSelectable(elem, selectable) {
					if (elem === undefined || elem.style === undefined) return;

					elem.onselectstart = selectable ? function () {
						return false;
					} : function () {};

					elem.style.MozUserSelect = selectable ? 'auto' : 'none';
					elem.style.KhtmlUserSelect = selectable ? 'auto' : 'none';
					elem.unselectable = selectable ? 'on' : 'off';
				},

				/**
     *
     * @param elem
     * @param horizontal
     * @param vert
     */
				makeFullscreen: function makeFullscreen(elem, hor, vert) {
					var vertical = vert;
					var horizontal = hor;

					if (_common2.default.isUndefined(horizontal)) {
						horizontal = true;
					}

					if (_common2.default.isUndefined(vertical)) {
						vertical = true;
					}

					elem.style.position = 'absolute';

					if (horizontal) {
						elem.style.left = 0;
						elem.style.right = 0;
					}
					if (vertical) {
						elem.style.top = 0;
						elem.style.bottom = 0;
					}
				},

				/**
     *
     * @param elem
     * @param eventType
     * @param params
     */
				fakeEvent: function fakeEvent(elem, eventType, pars, aux) {
					var params = pars || {};
					var className = EVENT_MAP_INV[eventType];
					if (!className) {
						throw new Error('Event type ' + eventType + ' not supported.');
					}
					var evt = document.createEvent(className);
					switch (className) {
						case 'MouseEvents':
							{
								var clientX = params.x || params.clientX || 0;
								var clientY = params.y || params.clientY || 0;
								evt.initMouseEvent(eventType, params.bubbles || false, params.cancelable || true, window, params.clickCount || 1, 0, // screen X
								0, // screen Y
								clientX, // client X
								clientY, // client Y
								false, false, false, false, 0, null);
								break;
							}
						case 'KeyboardEvents':
							{
								var init = evt.initKeyboardEvent || evt.initKeyEvent; // webkit || moz
								_common2.default.defaults(params, {
									cancelable: true,
									ctrlKey: false,
									altKey: false,
									shiftKey: false,
									metaKey: false,
									keyCode: undefined,
									charCode: undefined
								});
								init(eventType, params.bubbles || false, params.cancelable, window, params.ctrlKey, params.altKey, params.shiftKey, params.metaKey, params.keyCode, params.charCode);
								break;
							}
						default:
							{
								evt.initEvent(eventType, params.bubbles || false, params.cancelable || true);
								break;
							}
					}
					_common2.default.defaults(evt, aux);
					elem.dispatchEvent(evt);
				},

				/**
     *
     * @param elem
     * @param event
     * @param func
     * @param bool
     */
				bind: function bind(elem, event, func, newBool) {
					var bool = newBool || false;
					if (elem.addEventListener) {
						elem.addEventListener(event, func, bool);
					} else if (elem.attachEvent) {
						elem.attachEvent('on' + event, func);
					}
					return dom;
				},

				/**
     *
     * @param elem
     * @param event
     * @param func
     * @param bool
     */
				unbind: function unbind(elem, event, func, newBool) {
					var bool = newBool || false;
					if (elem.removeEventListener) {
						elem.removeEventListener(event, func, bool);
					} else if (elem.detachEvent) {
						elem.detachEvent('on' + event, func);
					}
					return dom;
				},

				/**
     *
     * @param elem
     * @param className
     */
				addClass: function addClass(elem, className) {
					if (elem.className === undefined) {
						elem.className = className;
					} else if (elem.className !== className) {
						var classes = elem.className.split(/ +/);
						if (classes.indexOf(className) === -1) {
							classes.push(className);
							elem.className = classes.join(' ').replace(/^\s+/, '').replace(/\s+$/, '');
						}
					}
					return dom;
				},

				/**
     *
     * @param elem
     * @param className
     */
				removeClass: function removeClass(elem, className) {
					if (className) {
						if (elem.className === className) {
							elem.removeAttribute('class');
						} else {
							var classes = elem.className.split(/ +/);
							var index = classes.indexOf(className);
							if (index !== -1) {
								classes.splice(index, 1);
								elem.className = classes.join(' ');
							}
						}
					} else {
						elem.className = undefined;
					}
					return dom;
				},

				hasClass: function hasClass(elem, className) {
					return new RegExp('(?:^|\\s+)' + className + '(?:\\s+|$)').test(elem.className) || false;
				},

				/**
     *
     * @param elem
     */
				getWidth: function getWidth(elem) {
					var style = getComputedStyle(elem);

					return cssValueToPixels(style['border-left-width']) + cssValueToPixels(style['border-right-width']) + cssValueToPixels(style['padding-left']) + cssValueToPixels(style['padding-right']) + cssValueToPixels(style.width);
				},

				/**
     *
     * @param elem
     */
				getHeight: function getHeight(elem) {
					var style = getComputedStyle(elem);

					return cssValueToPixels(style['border-top-width']) + cssValueToPixels(style['border-bottom-width']) + cssValueToPixels(style['padding-top']) + cssValueToPixels(style['padding-bottom']) + cssValueToPixels(style.height);
				},

				/**
     *
     * @param el
     */
				getOffset: function getOffset(el) {
					var elem = el;
					var offset = { left: 0, top: 0 };
					if (elem.offsetParent) {
						do {
							offset.left += elem.offsetLeft;
							offset.top += elem.offsetTop;
							elem = elem.offsetParent;
						} while (elem);
					}
					return offset;
				},

				// http://stackoverflow.com/posts/2684561/revisions
				/**
     *
     * @param elem
     */
				isActive: function isActive(elem) {
					return elem === document.activeElement && (elem.type || elem.href);
				}

			};

			exports.default = dom;

			/***/
		},
		/* 10 */
		/***/function (module, exports, __webpack_require__) {

			'use strict';

			exports.__esModule = true;

			var _Controller2 = __webpack_require__(7);

			var _Controller3 = _interopRequireDefault(_Controller2);

			var _dom = __webpack_require__(9);

			var _dom2 = _interopRequireDefault(_dom);

			var _common = __webpack_require__(5);

			var _common2 = _interopRequireDefault(_common);

			function _interopRequireDefault(obj) {
				return obj && obj.__esModule ? obj : { default: obj };
			}

			function _classCallCheck(instance, Constructor) {
				if (!(instance instanceof Constructor)) {
					throw new TypeError("Cannot call a class as a function");
				}
			}

			function _possibleConstructorReturn(self, call) {
				if (!self) {
					throw new ReferenceError("this hasn't been initialised - super() hasn't been called");
				}return call && ((typeof call === 'undefined' ? 'undefined' : _typeof2(call)) === "object" || typeof call === "function") ? call : self;
			}

			function _inherits(subClass, superClass) {
				if (typeof superClass !== "function" && superClass !== null) {
					throw new TypeError("Super expression must either be null or a function, not " + (typeof superClass === 'undefined' ? 'undefined' : _typeof2(superClass)));
				}subClass.prototype = Object.create(superClass && superClass.prototype, { constructor: { value: subClass, enumerable: false, writable: true, configurable: true } });if (superClass) Object.setPrototypeOf ? Object.setPrototypeOf(subClass, superClass) : subClass.__proto__ = superClass;
			} /**
      * dat-gui JavaScript Controller Library
      * http://code.google.com/p/dat-gui
      *
      * Copyright 2011 Data Arts Team, Google Creative Lab
      *
      * Licensed under the Apache License, Version 2.0 (the "License");
      * you may not use this file except in compliance with the License.
      * You may obtain a copy of the License at
      *
      * http://www.apache.org/licenses/LICENSE-2.0
      */

			/**
    * @class Provides a select input to alter the property of an object, using a
    * list of accepted values.
    *
    * @extends dat.controllers.Controller
    *
    * @param {Object} object The object to be manipulated
    * @param {string} property The name of the property to be manipulated
    * @param {Object|string[]} options A map of labels to acceptable values, or
    * a list of acceptable string values.
    *
    * @member dat.controllers
    */
			var OptionController = function (_Controller) {
				_inherits(OptionController, _Controller);

				function OptionController(object, property, opts) {
					_classCallCheck(this, OptionController);

					var _this2 = _possibleConstructorReturn(this, _Controller.call(this, object, property));

					var options = opts;

					var _this = _this2;

					/**
      * The drop down menu
      * @ignore
      */
					_this2.__select = document.createElement('select');

					if (_common2.default.isArray(options)) {
						(function () {
							var map = {};
							_common2.default.each(options, function (element) {
								map[element] = element;
							});
							options = map;
						})();
					}

					_common2.default.each(options, function (value, key) {
						var opt = document.createElement('option');
						opt.innerHTML = key;
						opt.setAttribute('value', value);
						_this.__select.appendChild(opt);
					});

					// Acknowledge original value
					_this2.updateDisplay();

					_dom2.default.bind(_this2.__select, 'change', function () {
						var desiredValue = this.options[this.selectedIndex].value;
						_this.setValue(desiredValue);
					});

					_this2.domElement.appendChild(_this2.__select);
					return _this2;
				}

				OptionController.prototype.setValue = function setValue(v) {
					var toReturn = _Controller.prototype.setValue.call(this, v);

					if (this.__onFinishChange) {
						this.__onFinishChange.call(this, this.getValue());
					}
					return toReturn;
				};

				OptionController.prototype.updateDisplay = function updateDisplay() {
					if (_dom2.default.isActive(this.__select)) return this; // prevent number from updating if user is trying to manually update
					this.__select.value = this.getValue();
					return _Controller.prototype.updateDisplay.call(this);
				};

				return OptionController;
			}(_Controller3.default);

			exports.default = OptionController;

			/***/
		},
		/* 11 */
		/***/function (module, exports, __webpack_require__) {

			'use strict';

			exports.__esModule = true;

			var _Controller2 = __webpack_require__(7);

			var _Controller3 = _interopRequireDefault(_Controller2);

			var _dom = __webpack_require__(9);

			var _dom2 = _interopRequireDefault(_dom);

			function _interopRequireDefault(obj) {
				return obj && obj.__esModule ? obj : { default: obj };
			}

			function _classCallCheck(instance, Constructor) {
				if (!(instance instanceof Constructor)) {
					throw new TypeError("Cannot call a class as a function");
				}
			}

			function _possibleConstructorReturn(self, call) {
				if (!self) {
					throw new ReferenceError("this hasn't been initialised - super() hasn't been called");
				}return call && ((typeof call === 'undefined' ? 'undefined' : _typeof2(call)) === "object" || typeof call === "function") ? call : self;
			}

			function _inherits(subClass, superClass) {
				if (typeof superClass !== "function" && superClass !== null) {
					throw new TypeError("Super expression must either be null or a function, not " + (typeof superClass === 'undefined' ? 'undefined' : _typeof2(superClass)));
				}subClass.prototype = Object.create(superClass && superClass.prototype, { constructor: { value: subClass, enumerable: false, writable: true, configurable: true } });if (superClass) Object.setPrototypeOf ? Object.setPrototypeOf(subClass, superClass) : subClass.__proto__ = superClass;
			} /**
      * dat-gui JavaScript Controller Library
      * http://code.google.com/p/dat-gui
      *
      * Copyright 2011 Data Arts Team, Google Creative Lab
      *
      * Licensed under the Apache License, Version 2.0 (the "License");
      * you may not use this file except in compliance with the License.
      * You may obtain a copy of the License at
      *
      * http://www.apache.org/licenses/LICENSE-2.0
      */

			/**
    * @class Provides a text input to alter the string property of an object.
    *
    * @extends dat.controllers.Controller
    *
    * @param {Object} object The object to be manipulated
    * @param {string} property The name of the property to be manipulated
    *
    * @member dat.controllers
    */
			var StringController = function (_Controller) {
				_inherits(StringController, _Controller);

				function StringController(object, property) {
					_classCallCheck(this, StringController);

					var _this2 = _possibleConstructorReturn(this, _Controller.call(this, object, property));

					var _this = _this2;

					function onChange() {
						_this.setValue(_this.__input.value);
					}

					function onBlur() {
						if (_this.__onFinishChange) {
							_this.__onFinishChange.call(_this, _this.getValue());
						}
					}

					_this2.__input = document.createElement('input');
					_this2.__input.setAttribute('type', 'text');

					_dom2.default.bind(_this2.__input, 'keyup', onChange);
					_dom2.default.bind(_this2.__input, 'change', onChange);
					_dom2.default.bind(_this2.__input, 'blur', onBlur);
					_dom2.default.bind(_this2.__input, 'keydown', function (e) {
						if (e.keyCode === 13) {
							this.blur();
						}
					});

					_this2.updateDisplay();

					_this2.domElement.appendChild(_this2.__input);
					return _this2;
				}

				StringController.prototype.updateDisplay = function updateDisplay() {
					// Stops the caret from moving on account of:
					// keyup -> setValue -> updateDisplay
					if (!_dom2.default.isActive(this.__input)) {
						this.__input.value = this.getValue();
					}
					return _Controller.prototype.updateDisplay.call(this);
				};

				return StringController;
			}(_Controller3.default);

			exports.default = StringController;

			/***/
		},
		/* 12 */
		/***/function (module, exports, __webpack_require__) {

			'use strict';

			exports.__esModule = true;

			var _Controller2 = __webpack_require__(7);

			var _Controller3 = _interopRequireDefault(_Controller2);

			var _common = __webpack_require__(5);

			var _common2 = _interopRequireDefault(_common);

			function _interopRequireDefault(obj) {
				return obj && obj.__esModule ? obj : { default: obj };
			}

			function _classCallCheck(instance, Constructor) {
				if (!(instance instanceof Constructor)) {
					throw new TypeError("Cannot call a class as a function");
				}
			}

			function _possibleConstructorReturn(self, call) {
				if (!self) {
					throw new ReferenceError("this hasn't been initialised - super() hasn't been called");
				}return call && ((typeof call === 'undefined' ? 'undefined' : _typeof2(call)) === "object" || typeof call === "function") ? call : self;
			}

			function _inherits(subClass, superClass) {
				if (typeof superClass !== "function" && superClass !== null) {
					throw new TypeError("Super expression must either be null or a function, not " + (typeof superClass === 'undefined' ? 'undefined' : _typeof2(superClass)));
				}subClass.prototype = Object.create(superClass && superClass.prototype, { constructor: { value: subClass, enumerable: false, writable: true, configurable: true } });if (superClass) Object.setPrototypeOf ? Object.setPrototypeOf(subClass, superClass) : subClass.__proto__ = superClass;
			} /**
      * dat-gui JavaScript Controller Library
      * http://code.google.com/p/dat-gui
      *
      * Copyright 2011 Data Arts Team, Google Creative Lab
      *
      * Licensed under the Apache License, Version 2.0 (the "License");
      * you may not use this file except in compliance with the License.
      * You may obtain a copy of the License at
      *
      * http://www.apache.org/licenses/LICENSE-2.0
      */

			function numDecimals(x) {
				var _x = x.toString();
				if (_x.indexOf('.') > -1) {
					return _x.length - _x.indexOf('.') - 1;
				}

				return 0;
			}

			/**
    * @class Represents a given property of an object that is a number.
    *
    * @extends dat.controllers.Controller
    *
    * @param {Object} object The object to be manipulated
    * @param {string} property The name of the property to be manipulated
    * @param {Object} [params] Optional parameters
    * @param {Number} [params.min] Minimum allowed value
    * @param {Number} [params.max] Maximum allowed value
    * @param {Number} [params.step] Increment by which to change value
    *
    * @member dat.controllers
    */

			var NumberController = function (_Controller) {
				_inherits(NumberController, _Controller);

				function NumberController(object, property, params) {
					_classCallCheck(this, NumberController);

					var _this = _possibleConstructorReturn(this, _Controller.call(this, object, property));

					var _params = params || {};

					_this.__min = _params.min;
					_this.__max = _params.max;
					_this.__step = _params.step;

					if (_common2.default.isUndefined(_this.__step)) {
						if (_this.initialValue === 0) {
							_this.__impliedStep = 1; // What are we, psychics?
						} else {
							// Hey Doug, check this out.
							_this.__impliedStep = Math.pow(10, Math.floor(Math.log(Math.abs(_this.initialValue)) / Math.LN10)) / 10;
						}
					} else {
						_this.__impliedStep = _this.__step;
					}

					_this.__precision = numDecimals(_this.__impliedStep);
					return _this;
				}

				NumberController.prototype.setValue = function setValue(v) {
					var _v = v;

					if (this.__min !== undefined && _v < this.__min) {
						_v = this.__min;
					} else if (this.__max !== undefined && _v > this.__max) {
						_v = this.__max;
					}

					if (this.__step !== undefined && _v % this.__step !== 0) {
						_v = Math.round(_v / this.__step) * this.__step;
					}

					return _Controller.prototype.setValue.call(this, _v);
				};

				/**
     * Specify a minimum value for <code>object[property]</code>.
     *
     * @param {Number} minValue The minimum value for
     * <code>object[property]</code>
     * @returns {dat.controllers.NumberController} this
     */

				NumberController.prototype.min = function min(v) {
					this.__min = v;
					return this;
				};

				/**
     * Specify a maximum value for <code>object[property]</code>.
     *
     * @param {Number} maxValue The maximum value for
     * <code>object[property]</code>
     * @returns {dat.controllers.NumberController} this
     */

				NumberController.prototype.max = function max(v) {
					this.__max = v;
					return this;
				};

				/**
     * Specify a step value that dat.controllers.NumberController
     * increments by.
     *
     * @param {Number} stepValue The step value for
     * dat.controllers.NumberController
     * @default if minimum and maximum specified increment is 1% of the
     * difference otherwise stepValue is 1
     * @returns {dat.controllers.NumberController} this
     */

				NumberController.prototype.step = function step(v) {
					this.__step = v;
					this.__impliedStep = v;
					this.__precision = numDecimals(v);
					return this;
				};

				return NumberController;
			}(_Controller3.default);

			exports.default = NumberController;

			/***/
		},
		/* 13 */
		/***/function (module, exports, __webpack_require__) {

			'use strict';

			exports.__esModule = true;

			var _NumberController2 = __webpack_require__(12);

			var _NumberController3 = _interopRequireDefault(_NumberController2);

			var _dom = __webpack_require__(9);

			var _dom2 = _interopRequireDefault(_dom);

			var _common = __webpack_require__(5);

			var _common2 = _interopRequireDefault(_common);

			function _interopRequireDefault(obj) {
				return obj && obj.__esModule ? obj : { default: obj };
			}

			function _classCallCheck(instance, Constructor) {
				if (!(instance instanceof Constructor)) {
					throw new TypeError("Cannot call a class as a function");
				}
			}

			function _possibleConstructorReturn(self, call) {
				if (!self) {
					throw new ReferenceError("this hasn't been initialised - super() hasn't been called");
				}return call && ((typeof call === 'undefined' ? 'undefined' : _typeof2(call)) === "object" || typeof call === "function") ? call : self;
			}

			function _inherits(subClass, superClass) {
				if (typeof superClass !== "function" && superClass !== null) {
					throw new TypeError("Super expression must either be null or a function, not " + (typeof superClass === 'undefined' ? 'undefined' : _typeof2(superClass)));
				}subClass.prototype = Object.create(superClass && superClass.prototype, { constructor: { value: subClass, enumerable: false, writable: true, configurable: true } });if (superClass) Object.setPrototypeOf ? Object.setPrototypeOf(subClass, superClass) : subClass.__proto__ = superClass;
			} /**
      * dat-gui JavaScript Controller Library
      * http://code.google.com/p/dat-gui
      *
      * Copyright 2011 Data Arts Team, Google Creative Lab
      *
      * Licensed under the Apache License, Version 2.0 (the "License");
      * you may not use this file except in compliance with the License.
      * You may obtain a copy of the License at
      *
      * http://www.apache.org/licenses/LICENSE-2.0
      */

			function roundToDecimal(value, decimals) {
				var tenTo = Math.pow(10, decimals);
				return Math.round(value * tenTo) / tenTo;
			}

			/**
    * @class Represents a given property of an object that is a number and
    * provides an input element with which to manipulate it.
    *
    * @extends dat.controllers.Controller
    * @extends dat.controllers.NumberController
    *
    * @param {Object} object The object to be manipulated
    * @param {string} property The name of the property to be manipulated
    * @param {Object} [params] Optional parameters
    * @param {Number} [params.min] Minimum allowed value
    * @param {Number} [params.max] Maximum allowed value
    * @param {Number} [params.step] Increment by which to change value
    *
    * @member dat.controllers
    */

			var NumberControllerBox = function (_NumberController) {
				_inherits(NumberControllerBox, _NumberController);

				function NumberControllerBox(object, property, params) {
					_classCallCheck(this, NumberControllerBox);

					var _this2 = _possibleConstructorReturn(this, _NumberController.call(this, object, property, params));

					_this2.__truncationSuspended = false;

					var _this = _this2;

					/**
      * {Number} Previous mouse y position
      * @ignore
      */
					var prevY = void 0;

					function onChange() {
						var attempted = parseFloat(_this.__input.value);
						if (!_common2.default.isNaN(attempted)) {
							_this.setValue(attempted);
						}
					}

					function onFinish() {
						if (_this.__onFinishChange) {
							_this.__onFinishChange.call(_this, _this.getValue());
						}
					}

					function onBlur() {
						onFinish();
					}

					function onMouseDrag(e) {
						var diff = prevY - e.clientY;
						_this.setValue(_this.getValue() + diff * _this.__impliedStep);

						prevY = e.clientY;
					}

					function onMouseUp() {
						_dom2.default.unbind(window, 'mousemove', onMouseDrag);
						_dom2.default.unbind(window, 'mouseup', onMouseUp);
						onFinish();
					}

					function onMouseDown(e) {
						_dom2.default.bind(window, 'mousemove', onMouseDrag);
						_dom2.default.bind(window, 'mouseup', onMouseUp);
						prevY = e.clientY;
					}

					_this2.__input = document.createElement('input');
					_this2.__input.setAttribute('type', 'text');

					// Makes it so manually specified values are not truncated.

					_dom2.default.bind(_this2.__input, 'change', onChange);
					_dom2.default.bind(_this2.__input, 'blur', onBlur);
					_dom2.default.bind(_this2.__input, 'mousedown', onMouseDown);
					_dom2.default.bind(_this2.__input, 'keydown', function (e) {
						// When pressing enter, you can be as precise as you want.
						if (e.keyCode === 13) {
							_this.__truncationSuspended = true;
							this.blur();
							_this.__truncationSuspended = false;
							onFinish();
						}
					});

					_this2.updateDisplay();

					_this2.domElement.appendChild(_this2.__input);
					return _this2;
				}

				NumberControllerBox.prototype.updateDisplay = function updateDisplay() {
					this.__input.value = this.__truncationSuspended ? this.getValue() : roundToDecimal(this.getValue(), this.__precision);
					return _NumberController.prototype.updateDisplay.call(this);
				};

				return NumberControllerBox;
			}(_NumberController3.default);

			exports.default = NumberControllerBox;

			/***/
		},
		/* 14 */
		/***/function (module, exports, __webpack_require__) {

			'use strict';

			exports.__esModule = true;

			var _NumberController2 = __webpack_require__(12);

			var _NumberController3 = _interopRequireDefault(_NumberController2);

			var _dom = __webpack_require__(9);

			var _dom2 = _interopRequireDefault(_dom);

			function _interopRequireDefault(obj) {
				return obj && obj.__esModule ? obj : { default: obj };
			}

			function _classCallCheck(instance, Constructor) {
				if (!(instance instanceof Constructor)) {
					throw new TypeError("Cannot call a class as a function");
				}
			}

			function _possibleConstructorReturn(self, call) {
				if (!self) {
					throw new ReferenceError("this hasn't been initialised - super() hasn't been called");
				}return call && ((typeof call === 'undefined' ? 'undefined' : _typeof2(call)) === "object" || typeof call === "function") ? call : self;
			}

			function _inherits(subClass, superClass) {
				if (typeof superClass !== "function" && superClass !== null) {
					throw new TypeError("Super expression must either be null or a function, not " + (typeof superClass === 'undefined' ? 'undefined' : _typeof2(superClass)));
				}subClass.prototype = Object.create(superClass && superClass.prototype, { constructor: { value: subClass, enumerable: false, writable: true, configurable: true } });if (superClass) Object.setPrototypeOf ? Object.setPrototypeOf(subClass, superClass) : subClass.__proto__ = superClass;
			} /**
      * dat-gui JavaScript Controller Library
      * http://code.google.com/p/dat-gui
      *
      * Copyright 2011 Data Arts Team, Google Creative Lab
      *
      * Licensed under the Apache License, Version 2.0 (the "License");
      * you may not use this file except in compliance with the License.
      * You may obtain a copy of the License at
      *
      * http://www.apache.org/licenses/LICENSE-2.0
      */

			function map(v, i1, i2, o1, o2) {
				return o1 + (o2 - o1) * ((v - i1) / (i2 - i1));
			}

			/**
    * @class Represents a given property of an object that is a number, contains
    * a minimum and maximum, and provides a slider element with which to
    * manipulate it. It should be noted that the slider element is made up of
    * <code>&lt;div&gt;</code> tags, <strong>not</strong> the html5
    * <code>&lt;slider&gt;</code> element.
    *
    * @extends dat.controllers.Controller
    * @extends dat.controllers.NumberController
    *
    * @param {Object} object The object to be manipulated
    * @param {string} property The name of the property to be manipulated
    * @param {Number} minValue Minimum allowed value
    * @param {Number} maxValue Maximum allowed value
    * @param {Number} stepValue Increment by which to change value
    *
    * @member dat.controllers
    */

			var NumberControllerSlider = function (_NumberController) {
				_inherits(NumberControllerSlider, _NumberController);

				function NumberControllerSlider(object, property, min, max, step) {
					_classCallCheck(this, NumberControllerSlider);

					var _this2 = _possibleConstructorReturn(this, _NumberController.call(this, object, property, { min: min, max: max, step: step }));

					var _this = _this2;

					_this2.__background = document.createElement('div');
					_this2.__foreground = document.createElement('div');

					_dom2.default.bind(_this2.__background, 'mousedown', onMouseDown);

					_dom2.default.addClass(_this2.__background, 'slider');
					_dom2.default.addClass(_this2.__foreground, 'slider-fg');

					function onMouseDown(e) {
						document.activeElement.blur();

						_dom2.default.bind(window, 'mousemove', onMouseDrag);
						_dom2.default.bind(window, 'mouseup', onMouseUp);

						onMouseDrag(e);
					}

					function onMouseDrag(e) {
						e.preventDefault();

						var bgRect = _this.__background.getBoundingClientRect();

						_this.setValue(map(e.clientX, bgRect.left, bgRect.right, _this.__min, _this.__max));

						return false;
					}

					function onMouseUp() {
						_dom2.default.unbind(window, 'mousemove', onMouseDrag);
						_dom2.default.unbind(window, 'mouseup', onMouseUp);
						if (_this.__onFinishChange) {
							_this.__onFinishChange.call(_this, _this.getValue());
						}
					}

					_this2.updateDisplay();

					_this2.__background.appendChild(_this2.__foreground);
					_this2.domElement.appendChild(_this2.__background);
					return _this2;
				}

				NumberControllerSlider.prototype.updateDisplay = function updateDisplay() {
					var pct = (this.getValue() - this.__min) / (this.__max - this.__min);
					this.__foreground.style.width = pct * 100 + '%';
					return _NumberController.prototype.updateDisplay.call(this);
				};

				return NumberControllerSlider;
			}(_NumberController3.default);

			exports.default = NumberControllerSlider;

			/***/
		},
		/* 15 */
		/***/function (module, exports, __webpack_require__) {

			'use strict';

			exports.__esModule = true;

			var _Controller2 = __webpack_require__(7);

			var _Controller3 = _interopRequireDefault(_Controller2);

			var _dom = __webpack_require__(9);

			var _dom2 = _interopRequireDefault(_dom);

			function _interopRequireDefault(obj) {
				return obj && obj.__esModule ? obj : { default: obj };
			}

			function _classCallCheck(instance, Constructor) {
				if (!(instance instanceof Constructor)) {
					throw new TypeError("Cannot call a class as a function");
				}
			}

			function _possibleConstructorReturn(self, call) {
				if (!self) {
					throw new ReferenceError("this hasn't been initialised - super() hasn't been called");
				}return call && ((typeof call === 'undefined' ? 'undefined' : _typeof2(call)) === "object" || typeof call === "function") ? call : self;
			}

			function _inherits(subClass, superClass) {
				if (typeof superClass !== "function" && superClass !== null) {
					throw new TypeError("Super expression must either be null or a function, not " + (typeof superClass === 'undefined' ? 'undefined' : _typeof2(superClass)));
				}subClass.prototype = Object.create(superClass && superClass.prototype, { constructor: { value: subClass, enumerable: false, writable: true, configurable: true } });if (superClass) Object.setPrototypeOf ? Object.setPrototypeOf(subClass, superClass) : subClass.__proto__ = superClass;
			} /**
      * dat-gui JavaScript Controller Library
      * http://code.google.com/p/dat-gui
      *
      * Copyright 2011 Data Arts Team, Google Creative Lab
      *
      * Licensed under the Apache License, Version 2.0 (the "License");
      * you may not use this file except in compliance with the License.
      * You may obtain a copy of the License at
      *
      * http://www.apache.org/licenses/LICENSE-2.0
      */

			/**
    * @class Provides a GUI interface to fire a specified method, a property of an object.
    *
    * @extends dat.controllers.Controller
    *
    * @param {Object} object The object to be manipulated
    * @param {string} property The name of the property to be manipulated
    *
    * @member dat.controllers
    */
			var FunctionController = function (_Controller) {
				_inherits(FunctionController, _Controller);

				function FunctionController(object, property, text) {
					_classCallCheck(this, FunctionController);

					var _this2 = _possibleConstructorReturn(this, _Controller.call(this, object, property));

					var _this = _this2;

					_this2.__button = document.createElement('div');
					_this2.__button.innerHTML = text === undefined ? 'Fire' : text;

					_dom2.default.bind(_this2.__button, 'click', function (e) {
						e.preventDefault();
						_this.fire();
						return false;
					});

					_dom2.default.addClass(_this2.__button, 'button');

					_this2.domElement.appendChild(_this2.__button);
					return _this2;
				}

				FunctionController.prototype.fire = function fire() {
					if (this.__onChange) {
						this.__onChange.call(this);
					}
					this.getValue().call(this.object);
					if (this.__onFinishChange) {
						this.__onFinishChange.call(this, this.getValue());
					}
				};

				return FunctionController;
			}(_Controller3.default);

			exports.default = FunctionController;

			/***/
		},
		/* 16 */
		/***/function (module, exports, __webpack_require__) {

			'use strict';

			exports.__esModule = true;

			var _Controller2 = __webpack_require__(7);

			var _Controller3 = _interopRequireDefault(_Controller2);

			var _dom = __webpack_require__(9);

			var _dom2 = _interopRequireDefault(_dom);

			var _Color = __webpack_require__(2);

			var _Color2 = _interopRequireDefault(_Color);

			var _interpret = __webpack_require__(3);

			var _interpret2 = _interopRequireDefault(_interpret);

			var _common = __webpack_require__(5);

			var _common2 = _interopRequireDefault(_common);

			function _interopRequireDefault(obj) {
				return obj && obj.__esModule ? obj : { default: obj };
			}

			function _classCallCheck(instance, Constructor) {
				if (!(instance instanceof Constructor)) {
					throw new TypeError("Cannot call a class as a function");
				}
			}

			function _possibleConstructorReturn(self, call) {
				if (!self) {
					throw new ReferenceError("this hasn't been initialised - super() hasn't been called");
				}return call && ((typeof call === 'undefined' ? 'undefined' : _typeof2(call)) === "object" || typeof call === "function") ? call : self;
			}

			function _inherits(subClass, superClass) {
				if (typeof superClass !== "function" && superClass !== null) {
					throw new TypeError("Super expression must either be null or a function, not " + (typeof superClass === 'undefined' ? 'undefined' : _typeof2(superClass)));
				}subClass.prototype = Object.create(superClass && superClass.prototype, { constructor: { value: subClass, enumerable: false, writable: true, configurable: true } });if (superClass) Object.setPrototypeOf ? Object.setPrototypeOf(subClass, superClass) : subClass.__proto__ = superClass;
			} /**
      * dat-gui JavaScript Controller Library
      * http://code.google.com/p/dat-gui
      *
      * Copyright 2011 Data Arts Team, Google Creative Lab
      *
      * Licensed under the Apache License, Version 2.0 (the "License");
      * you may not use this file except in compliance with the License.
      * You may obtain a copy of the License at
      *
      * http://www.apache.org/licenses/LICENSE-2.0
      */

			var ColorController = function (_Controller) {
				_inherits(ColorController, _Controller);

				function ColorController(object, property) {
					_classCallCheck(this, ColorController);

					var _this2 = _possibleConstructorReturn(this, _Controller.call(this, object, property));

					_this2.__color = new _Color2.default(_this2.getValue());
					_this2.__temp = new _Color2.default(0);

					var _this = _this2;

					_this2.domElement = document.createElement('div');

					_dom2.default.makeSelectable(_this2.domElement, false);

					_this2.__selector = document.createElement('div');
					_this2.__selector.className = 'selector';

					_this2.__saturation_field = document.createElement('div');
					_this2.__saturation_field.className = 'saturation-field';

					_this2.__field_knob = document.createElement('div');
					_this2.__field_knob.className = 'field-knob';
					_this2.__field_knob_border = '2px solid ';

					_this2.__hue_knob = document.createElement('div');
					_this2.__hue_knob.className = 'hue-knob';

					_this2.__hue_field = document.createElement('div');
					_this2.__hue_field.className = 'hue-field';

					_this2.__input = document.createElement('input');
					_this2.__input.type = 'text';
					_this2.__input_textShadow = '0 1px 1px ';

					_dom2.default.bind(_this2.__input, 'keydown', function (e) {
						if (e.keyCode === 13) {
							// on enter
							onBlur.call(this);
						}
					});

					_dom2.default.bind(_this2.__input, 'blur', onBlur);

					_dom2.default.bind(_this2.__selector, 'mousedown', function () /* e */{
						_dom2.default.addClass(this, 'drag').bind(window, 'mouseup', function () /* e */{
							_dom2.default.removeClass(_this.__selector, 'drag');
						});
					});

					var valueField = document.createElement('div');

					_common2.default.extend(_this2.__selector.style, {
						width: '122px',
						height: '102px',
						padding: '3px',
						backgroundColor: '#222',
						boxShadow: '0px 1px 3px rgba(0,0,0,0.3)'
					});

					_common2.default.extend(_this2.__field_knob.style, {
						position: 'absolute',
						width: '12px',
						height: '12px',
						border: _this2.__field_knob_border + (_this2.__color.v < 0.5 ? '#fff' : '#000'),
						boxShadow: '0px 1px 3px rgba(0,0,0,0.5)',
						borderRadius: '12px',
						zIndex: 1
					});

					_common2.default.extend(_this2.__hue_knob.style, {
						position: 'absolute',
						width: '15px',
						height: '2px',
						borderRight: '4px solid #fff',
						zIndex: 1
					});

					_common2.default.extend(_this2.__saturation_field.style, {
						width: '100px',
						height: '100px',
						border: '1px solid #555',
						marginRight: '3px',
						display: 'inline-block',
						cursor: 'pointer'
					});

					_common2.default.extend(valueField.style, {
						width: '100%',
						height: '100%',
						background: 'none'
					});

					linearGradient(valueField, 'top', 'rgba(0,0,0,0)', '#000');

					_common2.default.extend(_this2.__hue_field.style, {
						width: '15px',
						height: '100px',
						border: '1px solid #555',
						cursor: 'ns-resize',
						position: 'absolute',
						top: '3px',
						right: '3px'
					});

					hueGradient(_this2.__hue_field);

					_common2.default.extend(_this2.__input.style, {
						outline: 'none',
						//      width: '120px',
						textAlign: 'center',
						//      padding: '4px',
						//      marginBottom: '6px',
						color: '#fff',
						border: 0,
						fontWeight: 'bold',
						textShadow: _this2.__input_textShadow + 'rgba(0,0,0,0.7)'
					});

					_dom2.default.bind(_this2.__saturation_field, 'mousedown', fieldDown);
					_dom2.default.bind(_this2.__field_knob, 'mousedown', fieldDown);

					_dom2.default.bind(_this2.__hue_field, 'mousedown', function (e) {
						setH(e);
						_dom2.default.bind(window, 'mousemove', setH);
						_dom2.default.bind(window, 'mouseup', fieldUpH);
					});

					function fieldDown(e) {
						setSV(e);
						// document.body.style.cursor = 'none';
						_dom2.default.bind(window, 'mousemove', setSV);
						_dom2.default.bind(window, 'mouseup', fieldUpSV);
					}

					function fieldUpSV() {
						_dom2.default.unbind(window, 'mousemove', setSV);
						_dom2.default.unbind(window, 'mouseup', fieldUpSV);
						// document.body.style.cursor = 'default';
						onFinish();
					}

					function onBlur() {
						var i = (0, _interpret2.default)(this.value);
						if (i !== false) {
							_this.__color.__state = i;
							_this.setValue(_this.__color.toOriginal());
						} else {
							this.value = _this.__color.toString();
						}
					}

					function fieldUpH() {
						_dom2.default.unbind(window, 'mousemove', setH);
						_dom2.default.unbind(window, 'mouseup', fieldUpH);
						onFinish();
					}

					function onFinish() {
						if (_this.__onFinishChange) {
							_this.__onFinishChange.call(_this, _this.__color.toOriginal());
						}
					}

					_this2.__saturation_field.appendChild(valueField);
					_this2.__selector.appendChild(_this2.__field_knob);
					_this2.__selector.appendChild(_this2.__saturation_field);
					_this2.__selector.appendChild(_this2.__hue_field);
					_this2.__hue_field.appendChild(_this2.__hue_knob);

					_this2.domElement.appendChild(_this2.__input);
					_this2.domElement.appendChild(_this2.__selector);

					_this2.updateDisplay();

					function setSV(e) {
						e.preventDefault();

						var fieldRect = _this.__saturation_field.getBoundingClientRect();
						var s = (e.clientX - fieldRect.left) / (fieldRect.right - fieldRect.left);
						var v = 1 - (e.clientY - fieldRect.top) / (fieldRect.bottom - fieldRect.top);

						if (v > 1) {
							v = 1;
						} else if (v < 0) {
							v = 0;
						}

						if (s > 1) {
							s = 1;
						} else if (s < 0) {
							s = 0;
						}

						_this.__color.v = v;
						_this.__color.s = s;

						_this.setValue(_this.__color.toOriginal());

						return false;
					}

					function setH(e) {
						e.preventDefault();

						var fieldRect = _this.__hue_field.getBoundingClientRect();
						var h = 1 - (e.clientY - fieldRect.top) / (fieldRect.bottom - fieldRect.top);

						if (h > 1) {
							h = 1;
						} else if (h < 0) {
							h = 0;
						}

						_this.__color.h = h * 360;

						_this.setValue(_this.__color.toOriginal());

						return false;
					}
					return _this2;
				}

				ColorController.prototype.updateDisplay = function updateDisplay() {
					var i = (0, _interpret2.default)(this.getValue());

					if (i !== false) {
						var mismatch = false;

						// Check for mismatch on the interpreted value.

						_common2.default.each(_Color2.default.COMPONENTS, function (component) {
							if (!_common2.default.isUndefined(i[component]) && !_common2.default.isUndefined(this.__color.__state[component]) && i[component] !== this.__color.__state[component]) {
								mismatch = true;
								return {}; // break
							}
						}, this);

						// If nothing diverges, we keep our previous values
						// for statefulness, otherwise we recalculate fresh
						if (mismatch) {
							_common2.default.extend(this.__color.__state, i);
						}
					}

					_common2.default.extend(this.__temp.__state, this.__color.__state);

					this.__temp.a = 1;

					var flip = this.__color.v < 0.5 || this.__color.s > 0.5 ? 255 : 0;
					var _flip = 255 - flip;

					_common2.default.extend(this.__field_knob.style, {
						marginLeft: 100 * this.__color.s - 7 + 'px',
						marginTop: 100 * (1 - this.__color.v) - 7 + 'px',
						backgroundColor: this.__temp.toHexString(),
						border: this.__field_knob_border + 'rgb(' + flip + ',' + flip + ',' + flip + ')'
					});

					this.__hue_knob.style.marginTop = (1 - this.__color.h / 360) * 100 + 'px';

					this.__temp.s = 1;
					this.__temp.v = 1;

					linearGradient(this.__saturation_field, 'left', '#fff', this.__temp.toHexString());

					this.__input.value = this.__color.toString();

					_common2.default.extend(this.__input.style, {
						backgroundColor: this.__color.toHexString(),
						color: 'rgb(' + flip + ',' + flip + ',' + flip + ')',
						textShadow: this.__input_textShadow + 'rgba(' + _flip + ',' + _flip + ',' + _flip + ',.7)'
					});
				};

				return ColorController;
			}(_Controller3.default);

			var vendors = ['-moz-', '-o-', '-webkit-', '-ms-', ''];

			function linearGradient(elem, x, a, b) {
				elem.style.background = '';
				_common2.default.each(vendors, function (vendor) {
					elem.style.cssText += 'background: ' + vendor + 'linear-gradient(' + x + ', ' + a + ' 0%, ' + b + ' 100%); ';
				});
			}

			function hueGradient(elem) {
				elem.style.background = '';
				elem.style.cssText += 'background: -moz-linear-gradient(top,  #ff0000 0%, #ff00ff 17%, #0000ff 34%, #00ffff 50%, #00ff00 67%, #ffff00 84%, #ff0000 100%);';
				elem.style.cssText += 'background: -webkit-linear-gradient(top,  #ff0000 0%,#ff00ff 17%,#0000ff 34%,#00ffff 50%,#00ff00 67%,#ffff00 84%,#ff0000 100%);';
				elem.style.cssText += 'background: -o-linear-gradient(top,  #ff0000 0%,#ff00ff 17%,#0000ff 34%,#00ffff 50%,#00ff00 67%,#ffff00 84%,#ff0000 100%);';
				elem.style.cssText += 'background: -ms-linear-gradient(top,  #ff0000 0%,#ff00ff 17%,#0000ff 34%,#00ffff 50%,#00ff00 67%,#ffff00 84%,#ff0000 100%);';
				elem.style.cssText += 'background: linear-gradient(top,  #ff0000 0%,#ff00ff 17%,#0000ff 34%,#00ffff 50%,#00ff00 67%,#ffff00 84%,#ff0000 100%);';
			}

			exports.default = ColorController;

			/***/
		},
		/* 17 */
		/***/function (module, exports, __webpack_require__) {

			'use strict';

			var _typeof = typeof Symbol === "function" && _typeof2(Symbol.iterator) === "symbol" ? function (obj) {
				return typeof obj === 'undefined' ? 'undefined' : _typeof2(obj);
			} : function (obj) {
				return obj && typeof Symbol === "function" && obj.constructor === Symbol ? "symbol" : typeof obj === 'undefined' ? 'undefined' : _typeof2(obj);
			}; /**
       * dat-gui JavaScript Controller Library
       * http://code.google.com/p/dat-gui
       *
       * Copyright 2011 Data Arts Team, Google Creative Lab
       *
       * Licensed under the Apache License, Version 2.0 (the "License");
       * you may not use this file except in compliance with the License.
       * You may obtain a copy of the License at
       *
       * http://www.apache.org/licenses/LICENSE-2.0
       */

			var _css = __webpack_require__(18);

			var _css2 = _interopRequireDefault(_css);

			var _saveDialogue = __webpack_require__(19);

			var _saveDialogue2 = _interopRequireDefault(_saveDialogue);

			var _ControllerFactory = __webpack_require__(20);

			var _ControllerFactory2 = _interopRequireDefault(_ControllerFactory);

			var _Controller = __webpack_require__(7);

			var _Controller2 = _interopRequireDefault(_Controller);

			var _BooleanController = __webpack_require__(8);

			var _BooleanController2 = _interopRequireDefault(_BooleanController);

			var _FunctionController = __webpack_require__(15);

			var _FunctionController2 = _interopRequireDefault(_FunctionController);

			var _NumberControllerBox = __webpack_require__(13);

			var _NumberControllerBox2 = _interopRequireDefault(_NumberControllerBox);

			var _NumberControllerSlider = __webpack_require__(14);

			var _NumberControllerSlider2 = _interopRequireDefault(_NumberControllerSlider);

			var _ColorController = __webpack_require__(16);

			var _ColorController2 = _interopRequireDefault(_ColorController);

			var _requestAnimationFrame = __webpack_require__(21);

			var _requestAnimationFrame2 = _interopRequireDefault(_requestAnimationFrame);

			var _CenteredDiv = __webpack_require__(22);

			var _CenteredDiv2 = _interopRequireDefault(_CenteredDiv);

			var _dom = __webpack_require__(9);

			var _dom2 = _interopRequireDefault(_dom);

			var _common = __webpack_require__(5);

			var _common2 = _interopRequireDefault(_common);

			var _style = __webpack_require__(23);

			var _style2 = _interopRequireDefault(_style);

			function _interopRequireDefault(obj) {
				return obj && obj.__esModule ? obj : { default: obj };
			}

			// CSS to embed in build

			_css2.default.inject(_style2.default);

			/** Outer-most className for GUI's */
			var CSS_NAMESPACE = 'dg';

			var HIDE_KEY_CODE = 72;

			/** The only value shared between the JS and SCSS. Use caution. */
			var CLOSE_BUTTON_HEIGHT = 20;

			var DEFAULT_DEFAULT_PRESET_NAME = 'Default';

			var SUPPORTS_LOCAL_STORAGE = function () {
				try {
					return 'localStorage' in window && window.localStorage !== null;
				} catch (e) {
					return false;
				}
			}();

			var SAVE_DIALOGUE = void 0;

			/** Have we yet to create an autoPlace GUI? */
			var autoPlaceVirgin = true;

			/** Fixed position div that auto place GUI's go inside */
			var autoPlaceContainer = void 0;

			/** Are we hiding the GUI's ? */
			var hide = false;

			/** GUI's which should be hidden */
			var hideableGuis = [];

			/**
    * A lightweight controller library for JavaScript. It allows you to easily
    * manipulate variables and fire functions on the fly.
    * @class
    *
    * @member dat.gui
    *
    * @param {Object} [params]
    * @param {String} [params.name] The name of this GUI.
    * @param {Object} [params.load] JSON object representing the saved state of
    * this GUI.
    * @param {Boolean} [params.auto=true]
    * @param {dat.gui.GUI} [params.parent] The GUI I'm nested in.
    * @param {Boolean} [params.closed] If true, starts closed
    */
			var GUI = function GUI(pars) {
				var _this = this;

				var params = pars || {};

				/**
     * Outermost DOM Element
     * @type DOMElement
     */
				this.domElement = document.createElement('div');
				this.__ul = document.createElement('ul');
				this.domElement.appendChild(this.__ul);

				_dom2.default.addClass(this.domElement, CSS_NAMESPACE);

				/**
     * Nested GUI's by name
     * @ignore
     */
				this.__folders = {};

				this.__controllers = [];

				/**
     * List of objects I'm remembering for save, only used in top level GUI
     * @ignore
     */
				this.__rememberedObjects = [];

				/**
     * Maps the index of remembered objects to a map of controllers, only used
     * in top level GUI.
     *
     * @private
     * @ignore
     *
     * @example
     * [
     *  {
       *    propertyName: Controller,
       *    anotherPropertyName: Controller
       *  },
     *  {
       *    propertyName: Controller
       *  }
     * ]
     */
				this.__rememberedObjectIndecesToControllers = [];

				this.__listening = [];

				// Default parameters
				params = _common2.default.defaults(params, {
					autoPlace: true,
					width: GUI.DEFAULT_WIDTH
				});

				params = _common2.default.defaults(params, {
					resizable: params.autoPlace,
					hideable: params.autoPlace
				});

				if (!_common2.default.isUndefined(params.load)) {
					// Explicit preset
					if (params.preset) {
						params.load.preset = params.preset;
					}
				} else {
					params.load = { preset: DEFAULT_DEFAULT_PRESET_NAME };
				}

				if (_common2.default.isUndefined(params.parent) && params.hideable) {
					hideableGuis.push(this);
				}

				// Only root level GUI's are resizable.
				params.resizable = _common2.default.isUndefined(params.parent) && params.resizable;

				if (params.autoPlace && _common2.default.isUndefined(params.scrollable)) {
					params.scrollable = true;
				}
				//    params.scrollable = common.isUndefined(params.parent) && params.scrollable === true;

				// Not part of params because I don't want people passing this in via
				// constructor. Should be a 'remembered' value.
				var useLocalStorage = SUPPORTS_LOCAL_STORAGE && localStorage.getItem(getLocalStorageHash(this, 'isLocal')) === 'true';

				var saveToLocalStorage = void 0;

				Object.defineProperties(this,
				/** @lends dat.gui.GUI.prototype */
				{
					/**
      * The parent <code>GUI</code>
      * @type dat.gui.GUI
      */
					parent: {
						get: function get() {
							return params.parent;
						}
					},

					scrollable: {
						get: function get() {
							return params.scrollable;
						}
					},

					/**
      * Handles <code>GUI</code>'s element placement for you
      * @type Boolean
      */
					autoPlace: {
						get: function get() {
							return params.autoPlace;
						}
					},

					/**
      * The identifier for a set of saved values
      * @type String
      */
					preset: {
						get: function get() {
							if (_this.parent) {
								return _this.getRoot().preset;
							}

							return params.load.preset;
						},

						set: function set(v) {
							if (_this.parent) {
								_this.getRoot().preset = v;
							} else {
								params.load.preset = v;
							}
							setPresetSelectIndex(this);
							_this.revert();
						}
					},

					/**
      * The width of <code>GUI</code> element
      * @type Number
      */
					width: {
						get: function get() {
							return params.width;
						},
						set: function set(v) {
							params.width = v;
							setWidth(_this, v);
						}
					},

					/**
      * The name of <code>GUI</code>. Used for folders. i.e
      * a folder's name
      * @type String
      */
					name: {
						get: function get() {
							return params.name;
						},
						set: function set(v) {
							// TODO Check for collisions among sibling folders
							params.name = v;
							if (titleRowName) {
								titleRowName.innerHTML = params.name;
							}
						}
					},

					/**
      * Whether the <code>GUI</code> is collapsed or not
      * @type Boolean
      */
					closed: {
						get: function get() {
							return params.closed;
						},
						set: function set(v) {
							params.closed = v;
							if (params.closed) {
								_dom2.default.addClass(_this.__ul, GUI.CLASS_CLOSED);
							} else {
								_dom2.default.removeClass(_this.__ul, GUI.CLASS_CLOSED);
							}
							// For browsers that aren't going to respect the CSS transition,
							// Lets just check our height against the window height right off
							// the bat.
							this.onResize();

							if (_this.__closeButton) {
								_this.__closeButton.innerHTML = v ? GUI.TEXT_OPEN : GUI.TEXT_CLOSED;
							}
						}
					},

					/**
      * Contains all presets
      * @type Object
      */
					load: {
						get: function get() {
							return params.load;
						}
					},

					/**
      * Determines whether or not to use <a href="https://developer.mozilla.org/en/DOM/Storage#localStorage">localStorage</a> as the means for
      * <code>remember</code>ing
      * @type Boolean
      */
					useLocalStorage: {

						get: function get() {
							return useLocalStorage;
						},
						set: function set(bool) {
							if (SUPPORTS_LOCAL_STORAGE) {
								useLocalStorage = bool;
								if (bool) {
									_dom2.default.bind(window, 'unload', saveToLocalStorage);
								} else {
									_dom2.default.unbind(window, 'unload', saveToLocalStorage);
								}
								localStorage.setItem(getLocalStorageHash(_this, 'isLocal'), bool);
							}
						}
					}
				});

				// Are we a root level GUI?
				if (_common2.default.isUndefined(params.parent)) {
					params.closed = false;

					_dom2.default.addClass(this.domElement, GUI.CLASS_MAIN);
					_dom2.default.makeSelectable(this.domElement, false);

					// Are we supposed to be loading locally?
					if (SUPPORTS_LOCAL_STORAGE) {
						if (useLocalStorage) {
							_this.useLocalStorage = true;

							var savedGui = localStorage.getItem(getLocalStorageHash(this, 'gui'));

							if (savedGui) {
								params.load = JSON.parse(savedGui);
							}
						}
					}

					this.__closeButton = document.createElement('div');
					this.__closeButton.innerHTML = GUI.TEXT_CLOSED;
					_dom2.default.addClass(this.__closeButton, GUI.CLASS_CLOSE_BUTTON);
					this.domElement.appendChild(this.__closeButton);

					_dom2.default.bind(this.__closeButton, 'click', function () {
						_this.closed = !_this.closed;
					});
					// Oh, you're a nested GUI!
				} else {
					if (params.closed === undefined) {
						params.closed = true;
					}

					var _titleRowName = document.createTextNode(params.name);
					_dom2.default.addClass(_titleRowName, 'controller-name');

					var titleRow = addRow(_this, _titleRowName);

					var onClickTitle = function onClickTitle(e) {
						e.preventDefault();
						_this.closed = !_this.closed;
						return false;
					};

					_dom2.default.addClass(this.__ul, GUI.CLASS_CLOSED);

					_dom2.default.addClass(titleRow, 'title');
					_dom2.default.bind(titleRow, 'click', onClickTitle);

					if (!params.closed) {
						this.closed = false;
					}
				}

				if (params.autoPlace) {
					if (_common2.default.isUndefined(params.parent)) {
						if (autoPlaceVirgin) {
							autoPlaceContainer = document.createElement('div');
							_dom2.default.addClass(autoPlaceContainer, CSS_NAMESPACE);
							_dom2.default.addClass(autoPlaceContainer, GUI.CLASS_AUTO_PLACE_CONTAINER);
							document.body.appendChild(autoPlaceContainer);
							autoPlaceVirgin = false;
						}

						// Put it in the dom for you.
						autoPlaceContainer.appendChild(this.domElement);

						// Apply the auto styles
						_dom2.default.addClass(this.domElement, GUI.CLASS_AUTO_PLACE);
					}

					// Make it not elastic.
					if (!this.parent) {
						setWidth(_this, params.width);
					}
				}

				this.__resizeHandler = function () {
					_this.onResizeDebounced();
				};

				_dom2.default.bind(window, 'resize', this.__resizeHandler);
				_dom2.default.bind(this.__ul, 'webkitTransitionEnd', this.__resizeHandler);
				_dom2.default.bind(this.__ul, 'transitionend', this.__resizeHandler);
				_dom2.default.bind(this.__ul, 'oTransitionEnd', this.__resizeHandler);
				this.onResize();

				if (params.resizable) {
					addResizeHandle(this);
				}

				saveToLocalStorage = function saveToLocalStorage() {
					if (SUPPORTS_LOCAL_STORAGE && localStorage.getItem(getLocalStorageHash(_this, 'isLocal')) === 'true') {
						localStorage.setItem(getLocalStorageHash(_this, 'gui'), JSON.stringify(_this.getSaveObject()));
					}
				};

				// expose this method publicly
				this.saveToLocalStorageIfPossible = saveToLocalStorage;

				function resetWidth() {
					var root = _this.getRoot();
					root.width += 1;
					_common2.default.defer(function () {
						root.width -= 1;
					});
				}

				if (!params.parent) {
					resetWidth();
				}
			};

			GUI.toggleHide = function () {
				hide = !hide;
				_common2.default.each(hideableGuis, function (gui) {
					gui.domElement.style.display = hide ? 'none' : '';
				});
			};

			GUI.CLASS_AUTO_PLACE = 'a';
			GUI.CLASS_AUTO_PLACE_CONTAINER = 'ac';
			GUI.CLASS_MAIN = 'main';
			GUI.CLASS_CONTROLLER_ROW = 'cr';
			GUI.CLASS_TOO_TALL = 'taller-than-window';
			GUI.CLASS_CLOSED = 'closed';
			GUI.CLASS_CLOSE_BUTTON = 'close-button';
			GUI.CLASS_DRAG = 'drag';

			GUI.DEFAULT_WIDTH = 245;
			GUI.TEXT_CLOSED = 'Close Controls';
			GUI.TEXT_OPEN = 'Open Controls';

			GUI._keydownHandler = function (e) {
				if (document.activeElement.type !== 'text' && (e.which === HIDE_KEY_CODE || e.keyCode === HIDE_KEY_CODE)) {
					GUI.toggleHide();
				}
			};
			_dom2.default.bind(window, 'keydown', GUI._keydownHandler, false);

			_common2.default.extend(GUI.prototype,

			/** @lends dat.gui.GUI */
			{

				/**
     * @param object
     * @param property
     * @returns {dat.controllers.Controller} The new controller that was added.
     * @instance
     */
				add: function add(object, property) {
					return _add(this, object, property, {
						factoryArgs: Array.prototype.slice.call(arguments, 2)
					});
				},

				/**
     * @param object
     * @param property
     * @returns {dat.controllers.ColorController} The new controller that was added.
     * @instance
     */
				addColor: function addColor(object, property) {
					return _add(this, object, property, {
						color: true
					});
				},

				/**
     * @param controller
     * @instance
     */
				remove: function remove(controller) {
					// TODO listening?
					this.__ul.removeChild(controller.__li);
					this.__controllers.splice(this.__controllers.indexOf(controller), 1);
					var _this = this;
					_common2.default.defer(function () {
						_this.onResize();
					});
				},

				destroy: function destroy() {
					if (this.autoPlace) {
						autoPlaceContainer.removeChild(this.domElement);
					}

					_dom2.default.unbind(window, 'keydown', GUI._keydownHandler, false);
					_dom2.default.unbind(window, 'resize', this.__resizeHandler);

					if (this.saveToLocalStorageIfPossible) {
						_dom2.default.unbind(window, 'unload', this.saveToLocalStorageIfPossible);
					}
				},

				/**
     * @param name
     * @returns {dat.gui.GUI} The new folder.
     * @throws {Error} if this GUI already has a folder by the specified
     * name
     * @instance
     */
				addFolder: function addFolder(name) {
					// We have to prevent collisions on names in order to have a key
					// by which to remember saved values
					if (this.__folders[name] !== undefined) {
						throw new Error('You already have a folder in this GUI by the' + ' name "' + name + '"');
					}

					var newGuiParams = { name: name, parent: this };

					// We need to pass down the autoPlace trait so that we can
					// attach event listeners to open/close folder actions to
					// ensure that a scrollbar appears if the window is too short.
					newGuiParams.autoPlace = this.autoPlace;

					// Do we have saved appearance data for this folder?
					if (this.load && // Anything loaded?
					this.load.folders && // Was my parent a dead-end?
					this.load.folders[name]) {
						// Did daddy remember me?
						// Start me closed if I was closed
						newGuiParams.closed = this.load.folders[name].closed;

						// Pass down the loaded data
						newGuiParams.load = this.load.folders[name];
					}

					var gui = new GUI(newGuiParams);
					this.__folders[name] = gui;

					var li = addRow(this, gui.domElement);
					_dom2.default.addClass(li, 'folder');
					return gui;
				},

				open: function open() {
					this.closed = false;
				},

				close: function close() {
					this.closed = true;
				},

				onResize: function onResize() {
					// we debounce this function to prevent performance issues when rotating on tablet/mobile
					var root = this.getRoot();
					if (root.scrollable) {
						var top = _dom2.default.getOffset(root.__ul).top;
						var h = 0;

						_common2.default.each(root.__ul.childNodes, function (node) {
							if (!(root.autoPlace && node === root.__save_row)) {
								h += _dom2.default.getHeight(node);
							}
						});

						if (window.innerHeight - top - CLOSE_BUTTON_HEIGHT < h) {
							_dom2.default.addClass(root.domElement, GUI.CLASS_TOO_TALL);
							root.__ul.style.height = window.innerHeight - top - CLOSE_BUTTON_HEIGHT + 'px';
						} else {
							_dom2.default.removeClass(root.domElement, GUI.CLASS_TOO_TALL);
							root.__ul.style.height = 'auto';
						}
					}

					if (root.__resize_handle) {
						_common2.default.defer(function () {
							root.__resize_handle.style.height = root.__ul.offsetHeight + 'px';
						});
					}

					if (root.__closeButton) {
						root.__closeButton.style.width = root.width + 'px';
					}
				},

				onResizeDebounced: _common2.default.debounce(function () {
					this.onResize();
				}, 200),

				/**
     * Mark objects for saving. The order of these objects cannot change as
     * the GUI grows. When remembering new objects, append them to the end
     * of the list.
     *
     * @param {Object...} objects
     * @throws {Error} if not called on a top level GUI.
     * @instance
     */
				remember: function remember() {
					if (_common2.default.isUndefined(SAVE_DIALOGUE)) {
						SAVE_DIALOGUE = new _CenteredDiv2.default();
						SAVE_DIALOGUE.domElement.innerHTML = _saveDialogue2.default;
					}

					if (this.parent) {
						throw new Error('You can only call remember on a top level GUI.');
					}

					var _this = this;

					_common2.default.each(Array.prototype.slice.call(arguments), function (object) {
						if (_this.__rememberedObjects.length === 0) {
							addSaveMenu(_this);
						}
						if (_this.__rememberedObjects.indexOf(object) === -1) {
							_this.__rememberedObjects.push(object);
						}
					});

					if (this.autoPlace) {
						// Set save row width
						setWidth(this, this.width);
					}
				},

				/**
     * @returns {dat.gui.GUI} the topmost parent GUI of a nested GUI.
     * @instance
     */
				getRoot: function getRoot() {
					var gui = this;
					while (gui.parent) {
						gui = gui.parent;
					}
					return gui;
				},

				/**
     * @returns {Object} a JSON object representing the current state of
     * this GUI as well as its remembered properties.
     * @instance
     */
				getSaveObject: function getSaveObject() {
					var toReturn = this.load;
					toReturn.closed = this.closed;

					// Am I remembering any values?
					if (this.__rememberedObjects.length > 0) {
						toReturn.preset = this.preset;

						if (!toReturn.remembered) {
							toReturn.remembered = {};
						}

						toReturn.remembered[this.preset] = getCurrentPreset(this);
					}

					toReturn.folders = {};
					_common2.default.each(this.__folders, function (element, key) {
						toReturn.folders[key] = element.getSaveObject();
					});

					return toReturn;
				},

				save: function save() {
					if (!this.load.remembered) {
						this.load.remembered = {};
					}

					this.load.remembered[this.preset] = getCurrentPreset(this);
					markPresetModified(this, false);
					this.saveToLocalStorageIfPossible();
				},

				saveAs: function saveAs(presetName) {
					if (!this.load.remembered) {
						// Retain default values upon first save
						this.load.remembered = {};
						this.load.remembered[DEFAULT_DEFAULT_PRESET_NAME] = getCurrentPreset(this, true);
					}

					this.load.remembered[presetName] = getCurrentPreset(this);
					this.preset = presetName;
					addPresetOption(this, presetName, true);
					this.saveToLocalStorageIfPossible();
				},

				revert: function revert(gui) {
					_common2.default.each(this.__controllers, function (controller) {
						// Make revert work on Default.
						if (!this.getRoot().load.remembered) {
							controller.setValue(controller.initialValue);
						} else {
							recallSavedValue(gui || this.getRoot(), controller);
						}

						// fire onFinishChange callback
						if (controller.__onFinishChange) {
							controller.__onFinishChange.call(controller, controller.getValue());
						}
					}, this);

					_common2.default.each(this.__folders, function (folder) {
						folder.revert(folder);
					});

					if (!gui) {
						markPresetModified(this.getRoot(), false);
					}
				},

				listen: function listen(controller) {
					var init = this.__listening.length === 0;
					this.__listening.push(controller);
					if (init) {
						updateDisplays(this.__listening);
					}
				},

				updateDisplay: function updateDisplay() {
					_common2.default.each(this.__controllers, function (controller) {
						controller.updateDisplay();
					});
					_common2.default.each(this.__folders, function (folder) {
						folder.updateDisplay();
					});
				}
			});

			/**
    * Add a row to the end of the GUI or before another row.
    *
    * @param gui
    * @param [newDom] If specified, inserts the dom content in the new row
    * @param [liBefore] If specified, places the new row before another row
    */
			function addRow(gui, newDom, liBefore) {
				var li = document.createElement('li');
				if (newDom) {
					li.appendChild(newDom);
				}

				if (liBefore) {
					gui.__ul.insertBefore(li, liBefore);
				} else {
					gui.__ul.appendChild(li);
				}
				gui.onResize();
				return li;
			}

			function markPresetModified(gui, modified) {
				var opt = gui.__preset_select[gui.__preset_select.selectedIndex];

				// console.log('mark', modified, opt);
				if (modified) {
					opt.innerHTML = opt.value + '*';
				} else {
					opt.innerHTML = opt.value;
				}
			}

			function augmentController(gui, li, controller) {
				controller.__li = li;
				controller.__gui = gui;

				_common2.default.extend(controller, {
					options: function options(_options) {
						if (arguments.length > 1) {
							var nextSibling = controller.__li.nextElementSibling;
							controller.remove();

							return _add(gui, controller.object, controller.property, {
								before: nextSibling,
								factoryArgs: [_common2.default.toArray(arguments)]
							});
						}

						if (_common2.default.isArray(_options) || _common2.default.isObject(_options)) {
							var _nextSibling = controller.__li.nextElementSibling;
							controller.remove();

							return _add(gui, controller.object, controller.property, {
								before: _nextSibling,
								factoryArgs: [_options]
							});
						}
					},

					name: function name(v) {
						controller.__li.firstElementChild.firstElementChild.innerHTML = v;
						return controller;
					},

					listen: function listen() {
						controller.__gui.listen(controller);
						return controller;
					},

					remove: function remove() {
						controller.__gui.remove(controller);
						return controller;
					}
				});

				// All sliders should be accompanied by a box.
				if (controller instanceof _NumberControllerSlider2.default) {
					(function () {
						var box = new _NumberControllerBox2.default(controller.object, controller.property, { min: controller.__min, max: controller.__max, step: controller.__step });

						_common2.default.each(['updateDisplay', 'onChange', 'onFinishChange', 'step'], function (method) {
							var pc = controller[method];
							var pb = box[method];
							controller[method] = box[method] = function () {
								var args = Array.prototype.slice.call(arguments);
								pb.apply(box, args);
								return pc.apply(controller, args);
							};
						});

						_dom2.default.addClass(li, 'has-slider');
						controller.domElement.insertBefore(box.domElement, controller.domElement.firstElementChild);
					})();
				} else if (controller instanceof _NumberControllerBox2.default) {
					var r = function r(returned) {
						// Have we defined both boundaries?
						if (_common2.default.isNumber(controller.__min) && _common2.default.isNumber(controller.__max)) {
							// Well, then lets just replace this with a slider.

							// lets remember if the old controller had a specific name or was listening
							var oldName = controller.__li.firstElementChild.firstElementChild.innerHTML;
							var wasListening = controller.__gui.__listening.indexOf(controller) > -1;

							controller.remove();
							var newController = _add(gui, controller.object, controller.property, {
								before: controller.__li.nextElementSibling,
								factoryArgs: [controller.__min, controller.__max, controller.__step]
							});

							newController.name(oldName);
							if (wasListening) newController.listen();

							return newController;
						}

						return returned;
					};

					controller.min = _common2.default.compose(r, controller.min);
					controller.max = _common2.default.compose(r, controller.max);
				} else if (controller instanceof _BooleanController2.default) {
					_dom2.default.bind(li, 'click', function () {
						_dom2.default.fakeEvent(controller.__checkbox, 'click');
					});

					_dom2.default.bind(controller.__checkbox, 'click', function (e) {
						e.stopPropagation(); // Prevents double-toggle
					});
				} else if (controller instanceof _FunctionController2.default) {
					_dom2.default.bind(li, 'click', function () {
						_dom2.default.fakeEvent(controller.__button, 'click');
					});

					_dom2.default.bind(li, 'mouseover', function () {
						_dom2.default.addClass(controller.__button, 'hover');
					});

					_dom2.default.bind(li, 'mouseout', function () {
						_dom2.default.removeClass(controller.__button, 'hover');
					});
				} else if (controller instanceof _ColorController2.default) {
					_dom2.default.addClass(li, 'color');
					controller.updateDisplay = _common2.default.compose(function (val) {
						li.style.borderLeftColor = controller.__color.toString();
						return val;
					}, controller.updateDisplay);

					controller.updateDisplay();
				}

				controller.setValue = _common2.default.compose(function (val) {
					if (gui.getRoot().__preset_select && controller.isModified()) {
						markPresetModified(gui.getRoot(), true);
					}

					return val;
				}, controller.setValue);
			}

			function recallSavedValue(gui, controller) {
				// Find the topmost GUI, that's where remembered objects live.
				var root = gui.getRoot();

				// Does the object we're controlling match anything we've been told to
				// remember?
				var matchedIndex = root.__rememberedObjects.indexOf(controller.object);

				// Why yes, it does!
				if (matchedIndex !== -1) {
					// Let me fetch a map of controllers for thcommon.isObject.
					var controllerMap = root.__rememberedObjectIndecesToControllers[matchedIndex];

					// Ohp, I believe this is the first controller we've created for this
					// object. Lets make the map fresh.
					if (controllerMap === undefined) {
						controllerMap = {};
						root.__rememberedObjectIndecesToControllers[matchedIndex] = controllerMap;
					}

					// Keep track of this controller
					controllerMap[controller.property] = controller;

					// Okay, now have we saved any values for this controller?
					if (root.load && root.load.remembered) {
						var presetMap = root.load.remembered;

						// Which preset are we trying to load?
						var preset = void 0;

						if (presetMap[gui.preset]) {
							preset = presetMap[gui.preset];
						} else if (presetMap[DEFAULT_DEFAULT_PRESET_NAME]) {
							// Uhh, you can have the default instead?
							preset = presetMap[DEFAULT_DEFAULT_PRESET_NAME];
						} else {
							// Nada.
							return;
						}

						// Did the loaded object remember thcommon.isObject? &&  Did we remember this particular property?
						if (preset[matchedIndex] && preset[matchedIndex][controller.property] !== undefined) {
							// We did remember something for this guy ...
							var value = preset[matchedIndex][controller.property];

							// And that's what it is.
							controller.initialValue = value;
							controller.setValue(value);
						}
					}
				}
			}

			function _add(gui, object, property, params) {
				if (object[property] === undefined) {
					throw new Error('Object "' + object + '" has no property "' + property + '"');
				}

				var controller = void 0;

				if (params.color) {
					controller = new _ColorController2.default(object, property);
				} else {
					var factoryArgs = [object, property].concat(params.factoryArgs);
					controller = _ControllerFactory2.default.apply(gui, factoryArgs);
				}

				if (params.before instanceof _Controller2.default) {
					params.before = params.before.__li;
				}

				recallSavedValue(gui, controller);

				_dom2.default.addClass(controller.domElement, 'c');

				var name = document.createElement('span');
				_dom2.default.addClass(name, 'property-name');
				name.innerHTML = controller.property;

				var container = document.createElement('div');
				container.appendChild(name);
				container.appendChild(controller.domElement);

				var li = addRow(gui, container, params.before);

				_dom2.default.addClass(li, GUI.CLASS_CONTROLLER_ROW);
				if (controller instanceof _ColorController2.default) {
					_dom2.default.addClass(li, 'color');
				} else {
					_dom2.default.addClass(li, _typeof(controller.getValue()));
				}

				augmentController(gui, li, controller);

				gui.__controllers.push(controller);

				return controller;
			}

			function getLocalStorageHash(gui, key) {
				// TODO how does this deal with multiple GUI's?
				return document.location.href + '.' + key;
			}

			function addPresetOption(gui, name, setSelected) {
				var opt = document.createElement('option');
				opt.innerHTML = name;
				opt.value = name;
				gui.__preset_select.appendChild(opt);
				if (setSelected) {
					gui.__preset_select.selectedIndex = gui.__preset_select.length - 1;
				}
			}

			function showHideExplain(gui, explain) {
				explain.style.display = gui.useLocalStorage ? 'block' : 'none';
			}

			function addSaveMenu(gui) {
				var div = gui.__save_row = document.createElement('li');

				_dom2.default.addClass(gui.domElement, 'has-save');

				gui.__ul.insertBefore(div, gui.__ul.firstChild);

				_dom2.default.addClass(div, 'save-row');

				var gears = document.createElement('span');
				gears.innerHTML = '&nbsp;';
				_dom2.default.addClass(gears, 'button gears');

				// TODO replace with FunctionController
				var button = document.createElement('span');
				button.innerHTML = 'Save';
				_dom2.default.addClass(button, 'button');
				_dom2.default.addClass(button, 'save');

				var button2 = document.createElement('span');
				button2.innerHTML = 'New';
				_dom2.default.addClass(button2, 'button');
				_dom2.default.addClass(button2, 'save-as');

				var button3 = document.createElement('span');
				button3.innerHTML = 'Revert';
				_dom2.default.addClass(button3, 'button');
				_dom2.default.addClass(button3, 'revert');

				var select = gui.__preset_select = document.createElement('select');

				if (gui.load && gui.load.remembered) {
					_common2.default.each(gui.load.remembered, function (value, key) {
						addPresetOption(gui, key, key === gui.preset);
					});
				} else {
					addPresetOption(gui, DEFAULT_DEFAULT_PRESET_NAME, false);
				}

				_dom2.default.bind(select, 'change', function () {
					for (var index = 0; index < gui.__preset_select.length; index++) {
						gui.__preset_select[index].innerHTML = gui.__preset_select[index].value;
					}

					gui.preset = this.value;
				});

				div.appendChild(select);
				div.appendChild(gears);
				div.appendChild(button);
				div.appendChild(button2);
				div.appendChild(button3);

				if (SUPPORTS_LOCAL_STORAGE) {
					(function () {
						var explain = document.getElementById('dg-local-explain');
						var localStorageCheckBox = document.getElementById('dg-local-storage');
						var saveLocally = document.getElementById('dg-save-locally');

						saveLocally.style.display = 'block';

						if (localStorage.getItem(getLocalStorageHash(gui, 'isLocal')) === 'true') {
							localStorageCheckBox.setAttribute('checked', 'checked');
						}

						showHideExplain(gui, explain);

						// TODO: Use a boolean controller, fool!
						_dom2.default.bind(localStorageCheckBox, 'change', function () {
							gui.useLocalStorage = !gui.useLocalStorage;
							showHideExplain(gui, explain);
						});
					})();
				}

				var newConstructorTextArea = document.getElementById('dg-new-constructor');

				_dom2.default.bind(newConstructorTextArea, 'keydown', function (e) {
					if (e.metaKey && (e.which === 67 || e.keyCode === 67)) {
						SAVE_DIALOGUE.hide();
					}
				});

				_dom2.default.bind(gears, 'click', function () {
					newConstructorTextArea.innerHTML = JSON.stringify(gui.getSaveObject(), undefined, 2);
					SAVE_DIALOGUE.show();
					newConstructorTextArea.focus();
					newConstructorTextArea.select();
				});

				_dom2.default.bind(button, 'click', function () {
					gui.save();
				});

				_dom2.default.bind(button2, 'click', function () {
					var presetName = prompt('Enter a new preset name.');
					if (presetName) {
						gui.saveAs(presetName);
					}
				});

				_dom2.default.bind(button3, 'click', function () {
					gui.revert();
				});

				// div.appendChild(button2);
			}

			function addResizeHandle(gui) {
				var pmouseX = void 0;

				gui.__resize_handle = document.createElement('div');

				_common2.default.extend(gui.__resize_handle.style, {

					width: '6px',
					marginLeft: '-3px',
					height: '200px',
					cursor: 'ew-resize',
					position: 'absolute'
					// border: '1px solid blue'

				});

				function drag(e) {
					e.preventDefault();

					gui.width += pmouseX - e.clientX;
					gui.onResize();
					pmouseX = e.clientX;

					return false;
				}

				function dragStop() {
					_dom2.default.removeClass(gui.__closeButton, GUI.CLASS_DRAG);
					_dom2.default.unbind(window, 'mousemove', drag);
					_dom2.default.unbind(window, 'mouseup', dragStop);
				}

				function dragStart(e) {
					e.preventDefault();

					pmouseX = e.clientX;

					_dom2.default.addClass(gui.__closeButton, GUI.CLASS_DRAG);
					_dom2.default.bind(window, 'mousemove', drag);
					_dom2.default.bind(window, 'mouseup', dragStop);

					return false;
				}

				_dom2.default.bind(gui.__resize_handle, 'mousedown', dragStart);
				_dom2.default.bind(gui.__closeButton, 'mousedown', dragStart);

				gui.domElement.insertBefore(gui.__resize_handle, gui.domElement.firstElementChild);
			}

			function setWidth(gui, w) {
				gui.domElement.style.width = w + 'px';
				// Auto placed save-rows are position fixed, so we have to
				// set the width manually if we want it to bleed to the edge
				if (gui.__save_row && gui.autoPlace) {
					gui.__save_row.style.width = w + 'px';
				}
				if (gui.__closeButton) {
					gui.__closeButton.style.width = w + 'px';
				}
			}

			function getCurrentPreset(gui, useInitialValues) {
				var toReturn = {};

				// For each object I'm remembering
				_common2.default.each(gui.__rememberedObjects, function (val, index) {
					var savedValues = {};

					// The controllers I've made for thcommon.isObject by property
					var controllerMap = gui.__rememberedObjectIndecesToControllers[index];

					// Remember each value for each property
					_common2.default.each(controllerMap, function (controller, property) {
						savedValues[property] = useInitialValues ? controller.initialValue : controller.getValue();
					});

					// Save the values for thcommon.isObject
					toReturn[index] = savedValues;
				});

				return toReturn;
			}

			function setPresetSelectIndex(gui) {
				for (var index = 0; index < gui.__preset_select.length; index++) {
					if (gui.__preset_select[index].value === gui.preset) {
						gui.__preset_select.selectedIndex = index;
					}
				}
			}

			function updateDisplays(controllerArray) {
				if (controllerArray.length !== 0) {
					_requestAnimationFrame2.default.call(window, function () {
						updateDisplays(controllerArray);
					});
				}

				_common2.default.each(controllerArray, function (c) {
					c.updateDisplay();
				});
			}

			module.exports = GUI;

			/***/
		},
		/* 18 */
		/***/function (module, exports) {

			'use strict';

			/**
    * dat-gui JavaScript Controller Library
    * http://code.google.com/p/dat-gui
    *
    * Copyright 2011 Data Arts Team, Google Creative Lab
    *
    * Licensed under the Apache License, Version 2.0 (the "License");
    * you may not use this file except in compliance with the License.
    * You may obtain a copy of the License at
    *
    * http://www.apache.org/licenses/LICENSE-2.0
    */

			module.exports = {
				load: function load(url, indoc) {
					var doc = indoc || document;
					var link = doc.createElement('link');
					link.type = 'text/css';
					link.rel = 'stylesheet';
					link.href = url;
					doc.getElementsByTagName('head')[0].appendChild(link);
				},

				inject: function inject(css, indoc) {
					var doc = indoc || document;
					var injected = document.createElement('style');
					injected.type = 'text/css';
					injected.innerHTML = css;
					var head = doc.getElementsByTagName('head')[0];
					try {
						head.appendChild(injected);
					} catch (e) {// Unable to inject CSS, probably because of a Content Security Policy
					}
				}
			};

			/***/
		},
		/* 19 */
		/***/function (module, exports) {

			module.exports = "<div id=\"dg-save\" class=\"dg dialogue\">\n\n  Here's the new load parameter for your <code>GUI</code>'s constructor:\n\n  <textarea id=\"dg-new-constructor\"></textarea>\n\n  <div id=\"dg-save-locally\">\n\n    <input id=\"dg-local-storage\" type=\"checkbox\"/> Automatically save\n    values to <code>localStorage</code> on exit.\n\n    <div id=\"dg-local-explain\">The values saved to <code>localStorage</code> will\n      override those passed to <code>dat.GUI</code>'s constructor. This makes it\n      easier to work incrementally, but <code>localStorage</code> is fragile,\n      and your friends may not see the same values you do.\n\n    </div>\n\n  </div>\n\n</div>";

			/***/
		},
		/* 20 */
		/***/function (module, exports, __webpack_require__) {

			'use strict';

			exports.__esModule = true;

			var _OptionController = __webpack_require__(10);

			var _OptionController2 = _interopRequireDefault(_OptionController);

			var _NumberControllerBox = __webpack_require__(13);

			var _NumberControllerBox2 = _interopRequireDefault(_NumberControllerBox);

			var _NumberControllerSlider = __webpack_require__(14);

			var _NumberControllerSlider2 = _interopRequireDefault(_NumberControllerSlider);

			var _StringController = __webpack_require__(11);

			var _StringController2 = _interopRequireDefault(_StringController);

			var _FunctionController = __webpack_require__(15);

			var _FunctionController2 = _interopRequireDefault(_FunctionController);

			var _BooleanController = __webpack_require__(8);

			var _BooleanController2 = _interopRequireDefault(_BooleanController);

			var _common = __webpack_require__(5);

			var _common2 = _interopRequireDefault(_common);

			function _interopRequireDefault(obj) {
				return obj && obj.__esModule ? obj : { default: obj };
			}

			var ControllerFactory = function ControllerFactory(object, property) {
				var initialValue = object[property];

				// Providing options?
				if (_common2.default.isArray(arguments[2]) || _common2.default.isObject(arguments[2])) {
					return new _OptionController2.default(object, property, arguments[2]);
				}

				// Providing a map?
				if (_common2.default.isNumber(initialValue)) {
					// Has min and max? (slider)
					if (_common2.default.isNumber(arguments[2]) && _common2.default.isNumber(arguments[3])) {
						// has step?
						if (_common2.default.isNumber(arguments[4])) {
							return new _NumberControllerSlider2.default(object, property, arguments[2], arguments[3], arguments[4]);
						}

						return new _NumberControllerSlider2.default(object, property, arguments[2], arguments[3]);
					}

					// number box
					if (_common2.default.isNumber(arguments[4])) {
						// has step
						return new _NumberControllerBox2.default(object, property, { min: arguments[2], max: arguments[3], step: arguments[4] });
					}
					return new _NumberControllerBox2.default(object, property, { min: arguments[2], max: arguments[3] });
				}

				if (_common2.default.isString(initialValue)) {
					return new _StringController2.default(object, property);
				}

				if (_common2.default.isFunction(initialValue)) {
					return new _FunctionController2.default(object, property, '');
				}

				if (_common2.default.isBoolean(initialValue)) {
					return new _BooleanController2.default(object, property);
				}

				return null;
			}; /**
       * dat-gui JavaScript Controller Library
       * http://code.google.com/p/dat-gui
       *
       * Copyright 2011 Data Arts Team, Google Creative Lab
       *
       * Licensed under the Apache License, Version 2.0 (the "License");
       * you may not use this file except in compliance with the License.
       * You may obtain a copy of the License at
       *
       * http://www.apache.org/licenses/LICENSE-2.0
       */

			exports.default = ControllerFactory;

			/***/
		},
		/* 21 */
		/***/function (module, exports) {

			"use strict";

			exports.__esModule = true;
			/**
    * dat-gui JavaScript Controller Library
    * http://code.google.com/p/dat-gui
    *
    * Copyright 2011 Data Arts Team, Google Creative Lab
    *
    * Licensed under the Apache License, Version 2.0 (the "License");
    * you may not use this file except in compliance with the License.
    * You may obtain a copy of the License at
    *
    * http://www.apache.org/licenses/LICENSE-2.0
    */

			function requestAnimationFrame(callback) {
				setTimeout(callback, 1000 / 60);
			}

			exports.default = window.requestAnimationFrame || window.webkitRequestAnimationFrame || window.mozRequestAnimationFrame || window.oRequestAnimationFrame || window.msRequestAnimationFrame || requestAnimationFrame;

			/***/
		},
		/* 22 */
		/***/function (module, exports, __webpack_require__) {

			'use strict';

			exports.__esModule = true;

			var _dom = __webpack_require__(9);

			var _dom2 = _interopRequireDefault(_dom);

			var _common = __webpack_require__(5);

			var _common2 = _interopRequireDefault(_common);

			function _interopRequireDefault(obj) {
				return obj && obj.__esModule ? obj : { default: obj };
			}

			function _classCallCheck(instance, Constructor) {
				if (!(instance instanceof Constructor)) {
					throw new TypeError("Cannot call a class as a function");
				}
			} /**
      * dat-gui JavaScript Controller Library
      * http://code.google.com/p/dat-gui
      *
      * Copyright 2011 Data Arts Team, Google Creative Lab
      *
      * Licensed under the Apache License, Version 2.0 (the "License");
      * you may not use this file except in compliance with the License.
      * You may obtain a copy of the License at
      *
      * http://www.apache.org/licenses/LICENSE-2.0
      */

			var CenteredDiv = function () {
				function CenteredDiv() {
					_classCallCheck(this, CenteredDiv);

					this.backgroundElement = document.createElement('div');
					_common2.default.extend(this.backgroundElement.style, {
						backgroundColor: 'rgba(0,0,0,0.8)',
						top: 0,
						left: 0,
						display: 'none',
						zIndex: '1000',
						opacity: 0,
						WebkitTransition: 'opacity 0.2s linear',
						transition: 'opacity 0.2s linear'
					});

					_dom2.default.makeFullscreen(this.backgroundElement);
					this.backgroundElement.style.position = 'fixed';

					this.domElement = document.createElement('div');
					_common2.default.extend(this.domElement.style, {
						position: 'fixed',
						display: 'none',
						zIndex: '1001',
						opacity: 0,
						WebkitTransition: '-webkit-transform 0.2s ease-out, opacity 0.2s linear',
						transition: 'transform 0.2s ease-out, opacity 0.2s linear'
					});

					document.body.appendChild(this.backgroundElement);
					document.body.appendChild(this.domElement);

					var _this = this;
					_dom2.default.bind(this.backgroundElement, 'click', function () {
						_this.hide();
					});
				}

				CenteredDiv.prototype.show = function show() {
					var _this = this;

					this.backgroundElement.style.display = 'block';

					this.domElement.style.display = 'block';
					this.domElement.style.opacity = 0;
					//    this.domElement.style.top = '52%';
					this.domElement.style.webkitTransform = 'scale(1.1)';

					this.layout();

					_common2.default.defer(function () {
						_this.backgroundElement.style.opacity = 1;
						_this.domElement.style.opacity = 1;
						_this.domElement.style.webkitTransform = 'scale(1)';
					});
				};

				/**
     * Hide centered div
     */

				CenteredDiv.prototype.hide = function hide() {
					var _this = this;

					var hide = function hide() {
						_this.domElement.style.display = 'none';
						_this.backgroundElement.style.display = 'none';

						_dom2.default.unbind(_this.domElement, 'webkitTransitionEnd', hide);
						_dom2.default.unbind(_this.domElement, 'transitionend', hide);
						_dom2.default.unbind(_this.domElement, 'oTransitionEnd', hide);
					};

					_dom2.default.bind(this.domElement, 'webkitTransitionEnd', hide);
					_dom2.default.bind(this.domElement, 'transitionend', hide);
					_dom2.default.bind(this.domElement, 'oTransitionEnd', hide);

					this.backgroundElement.style.opacity = 0;
					//    this.domElement.style.top = '48%';
					this.domElement.style.opacity = 0;
					this.domElement.style.webkitTransform = 'scale(1.1)';
				};

				CenteredDiv.prototype.layout = function layout() {
					this.domElement.style.left = window.innerWidth / 2 - _dom2.default.getWidth(this.domElement) / 2 + 'px';
					this.domElement.style.top = window.innerHeight / 2 - _dom2.default.getHeight(this.domElement) / 2 + 'px';
				};

				return CenteredDiv;
			}();

			exports.default = CenteredDiv;

			/***/
		},
		/* 23 */
		/***/function (module, exports, __webpack_require__) {

			exports = module.exports = __webpack_require__(24)();
			// imports


			// module
			exports.push([module.id, ".dg {\n  /** Clear list styles */\n  /* Auto-place container */\n  /* Auto-placed GUI's */\n  /* Line items that don't contain folders. */\n  /** Folder names */\n  /** Hides closed items */\n  /** Controller row */\n  /** Name-half (left) */\n  /** Controller-half (right) */\n  /** Controller placement */\n  /** Shorter number boxes when slider is present. */\n  /** Ensure the entire boolean and function row shows a hand */ }\n  .dg ul {\n    list-style: none;\n    margin: 0;\n    padding: 0;\n    width: 100%;\n    clear: both; }\n  .dg.ac {\n    position: fixed;\n    top: 0;\n    left: 0;\n    right: 0;\n    height: 0;\n    z-index: 0; }\n  .dg:not(.ac) .main {\n    /** Exclude mains in ac so that we don't hide close button */\n    overflow: hidden; }\n  .dg.main {\n    -webkit-transition: opacity 0.1s linear;\n    -o-transition: opacity 0.1s linear;\n    -moz-transition: opacity 0.1s linear;\n    transition: opacity 0.1s linear; }\n    .dg.main.taller-than-window {\n      overflow-y: auto; }\n      .dg.main.taller-than-window .close-button {\n        opacity: 1;\n        /* TODO, these are style notes */\n        margin-top: -1px;\n        border-top: 1px solid #2c2c2c; }\n    .dg.main ul.closed .close-button {\n      opacity: 1 !important; }\n    .dg.main:hover .close-button,\n    .dg.main .close-button.drag {\n      opacity: 1; }\n    .dg.main .close-button {\n      /*opacity: 0;*/\n      -webkit-transition: opacity 0.1s linear;\n      -o-transition: opacity 0.1s linear;\n      -moz-transition: opacity 0.1s linear;\n      transition: opacity 0.1s linear;\n      border: 0;\n      position: absolute;\n      line-height: 19px;\n      height: 20px;\n      /* TODO, these are style notes */\n      cursor: pointer;\n      text-align: center;\n      background-color: #000; }\n      .dg.main .close-button:hover {\n        background-color: #111; }\n  .dg.a {\n    float: right;\n    margin-right: 15px;\n    overflow-x: hidden; }\n    .dg.a.has-save > ul {\n      margin-top: 27px; }\n      .dg.a.has-save > ul.closed {\n        margin-top: 0; }\n    .dg.a .save-row {\n      position: fixed;\n      top: 0;\n      z-index: 1002; }\n  .dg li {\n    -webkit-transition: height 0.1s ease-out;\n    -o-transition: height 0.1s ease-out;\n    -moz-transition: height 0.1s ease-out;\n    transition: height 0.1s ease-out; }\n  .dg li:not(.folder) {\n    cursor: auto;\n    height: 27px;\n    line-height: 27px;\n    overflow: hidden;\n    padding: 0 4px 0 5px; }\n  .dg li.folder {\n    padding: 0;\n    border-left: 4px solid transparent; }\n  .dg li.title {\n    cursor: pointer;\n    margin-left: -4px; }\n  .dg .closed li:not(.title),\n  .dg .closed ul li,\n  .dg .closed ul li > * {\n    height: 0;\n    overflow: hidden;\n    border: 0; }\n  .dg .cr {\n    clear: both;\n    padding-left: 3px;\n    height: 27px; }\n  .dg .property-name {\n    cursor: default;\n    float: left;\n    clear: left;\n    width: 40%;\n    overflow: hidden;\n    text-overflow: ellipsis; }\n  .dg .c {\n    float: left;\n    width: 60%; }\n  .dg .c input[type=text] {\n    border: 0;\n    margin-top: 4px;\n    padding: 3px;\n    width: 100%;\n    float: right; }\n  .dg .has-slider input[type=text] {\n    width: 30%;\n    /*display: none;*/\n    margin-left: 0; }\n  .dg .slider {\n    float: left;\n    width: 66%;\n    margin-left: -5px;\n    margin-right: 0;\n    height: 19px;\n    margin-top: 4px; }\n  .dg .slider-fg {\n    height: 100%; }\n  .dg .c input[type=checkbox] {\n    margin-top: 9px; }\n  .dg .c select {\n    margin-top: 5px; }\n  .dg .cr.function,\n  .dg .cr.function .property-name,\n  .dg .cr.function *,\n  .dg .cr.boolean,\n  .dg .cr.boolean * {\n    cursor: pointer; }\n  .dg .selector {\n    display: none;\n    position: absolute;\n    margin-left: -9px;\n    margin-top: 23px;\n    z-index: 10; }\n  .dg .c:hover .selector,\n  .dg .selector.drag {\n    display: block; }\n  .dg li.save-row {\n    padding: 0; }\n    .dg li.save-row .button {\n      display: inline-block;\n      padding: 0px 6px; }\n  .dg.dialogue {\n    background-color: #222;\n    width: 460px;\n    padding: 15px;\n    font-size: 13px;\n    line-height: 15px; }\n\n/* TODO Separate style and structure */\n#dg-new-constructor {\n  padding: 10px;\n  color: #222;\n  font-family: Monaco, monospace;\n  font-size: 10px;\n  border: 0;\n  resize: none;\n  box-shadow: inset 1px 1px 1px #888;\n  word-wrap: break-word;\n  margin: 12px 0;\n  display: block;\n  width: 440px;\n  overflow-y: scroll;\n  height: 100px;\n  position: relative; }\n\n#dg-local-explain {\n  display: none;\n  font-size: 11px;\n  line-height: 17px;\n  border-radius: 3px;\n  background-color: #333;\n  padding: 8px;\n  margin-top: 10px; }\n  #dg-local-explain code {\n    font-size: 10px; }\n\n#dat-gui-save-locally {\n  display: none; }\n\n/** Main type */\n.dg {\n  color: #eee;\n  font: 11px 'Lucida Grande', sans-serif;\n  text-shadow: 0 -1px 0 #111;\n  /** Auto place */\n  /* Controller row, <li> */\n  /** Controllers */ }\n  .dg.main {\n    /** Scrollbar */ }\n    .dg.main::-webkit-scrollbar {\n      width: 5px;\n      background: #1a1a1a; }\n    .dg.main::-webkit-scrollbar-corner {\n      height: 0;\n      display: none; }\n    .dg.main::-webkit-scrollbar-thumb {\n      border-radius: 5px;\n      background: #676767; }\n  .dg li:not(.folder) {\n    background: #1a1a1a;\n    border-bottom: 1px solid #2c2c2c; }\n  .dg li.save-row {\n    line-height: 25px;\n    background: #dad5cb;\n    border: 0; }\n    .dg li.save-row select {\n      margin-left: 5px;\n      width: 108px; }\n    .dg li.save-row .button {\n      margin-left: 5px;\n      margin-top: 1px;\n      border-radius: 2px;\n      font-size: 9px;\n      line-height: 7px;\n      padding: 4px 4px 5px 4px;\n      background: #c5bdad;\n      color: #fff;\n      text-shadow: 0 1px 0 #b0a58f;\n      box-shadow: 0 -1px 0 #b0a58f;\n      cursor: pointer; }\n      .dg li.save-row .button.gears {\n        background: #c5bdad url(data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAAsAAAANCAYAAAB/9ZQ7AAAAGXRFWHRTb2Z0d2FyZQBBZG9iZSBJbWFnZVJlYWR5ccllPAAAAQJJREFUeNpiYKAU/P//PwGIC/ApCABiBSAW+I8AClAcgKxQ4T9hoMAEUrxx2QSGN6+egDX+/vWT4e7N82AMYoPAx/evwWoYoSYbACX2s7KxCxzcsezDh3evFoDEBYTEEqycggWAzA9AuUSQQgeYPa9fPv6/YWm/Acx5IPb7ty/fw+QZblw67vDs8R0YHyQhgObx+yAJkBqmG5dPPDh1aPOGR/eugW0G4vlIoTIfyFcA+QekhhHJhPdQxbiAIguMBTQZrPD7108M6roWYDFQiIAAv6Aow/1bFwXgis+f2LUAynwoIaNcz8XNx3Dl7MEJUDGQpx9gtQ8YCueB+D26OECAAQDadt7e46D42QAAAABJRU5ErkJggg==) 2px 1px no-repeat;\n        height: 7px;\n        width: 8px; }\n      .dg li.save-row .button:hover {\n        background-color: #bab19e;\n        box-shadow: 0 -1px 0 #b0a58f; }\n  .dg li.folder {\n    border-bottom: 0; }\n  .dg li.title {\n    padding-left: 16px;\n    background: #000 url(data:image/gif;base64,R0lGODlhBQAFAJEAAP////Pz8////////yH5BAEAAAIALAAAAAAFAAUAAAIIlI+hKgFxoCgAOw==) 6px 10px no-repeat;\n    cursor: pointer;\n    border-bottom: 1px solid rgba(255, 255, 255, 0.2); }\n  .dg .closed li.title {\n    background-image: url(data:image/gif;base64,R0lGODlhBQAFAJEAAP////Pz8////////yH5BAEAAAIALAAAAAAFAAUAAAIIlGIWqMCbWAEAOw==); }\n  .dg .cr.boolean {\n    border-left: 3px solid #806787; }\n  .dg .cr.color {\n    border-left: 3px solid; }\n  .dg .cr.function {\n    border-left: 3px solid #e61d5f; }\n  .dg .cr.number {\n    border-left: 3px solid #2FA1D6; }\n    .dg .cr.number input[type=text] {\n      color: #2FA1D6; }\n  .dg .cr.string {\n    border-left: 3px solid #1ed36f; }\n    .dg .cr.string input[type=text] {\n      color: #1ed36f; }\n  .dg .cr.function:hover, .dg .cr.boolean:hover {\n    background: #111; }\n  .dg .c input[type=text] {\n    background: #303030;\n    outline: none; }\n    .dg .c input[type=text]:hover {\n      background: #3c3c3c; }\n    .dg .c input[type=text]:focus {\n      background: #494949;\n      color: #fff; }\n  .dg .c .slider {\n    background: #303030;\n    cursor: ew-resize; }\n  .dg .c .slider-fg {\n    background: #2FA1D6;\n    max-width: 100%; }\n  .dg .c .slider:hover {\n    background: #3c3c3c; }\n    .dg .c .slider:hover .slider-fg {\n      background: #44abda; }\n", ""]);

			// exports


			/***/
		},
		/* 24 */
		/***/function (module, exports) {

			/*
   	MIT License http://www.opensource.org/licenses/mit-license.php
   	Author Tobias Koppers @sokra
   */
			// css base code, injected by the css-loader
			module.exports = function () {
				var list = [];

				// return the list of modules as css string
				list.toString = function toString() {
					var result = [];
					for (var i = 0; i < this.length; i++) {
						var item = this[i];
						if (item[2]) {
							result.push("@media " + item[2] + "{" + item[1] + "}");
						} else {
							result.push(item[1]);
						}
					}
					return result.join("");
				};

				// import a list of modules into the list
				list.i = function (modules, mediaQuery) {
					if (typeof modules === "string") modules = [[null, modules, ""]];
					var alreadyImportedModules = {};
					for (var i = 0; i < this.length; i++) {
						var id = this[i][0];
						if (typeof id === "number") alreadyImportedModules[id] = true;
					}
					for (i = 0; i < modules.length; i++) {
						var item = modules[i];
						// skip already imported module
						// this implementation is not 100% perfect for weird media query combinations
						//  when a module is imported multiple times with different media queries.
						//  I hope this will never occur (Hey this way we have smaller bundles)
						if (typeof item[0] !== "number" || !alreadyImportedModules[item[0]]) {
							if (mediaQuery && !item[2]) {
								item[2] = mediaQuery;
							} else if (mediaQuery) {
								item[2] = "(" + item[2] + ") and (" + mediaQuery + ")";
							}
							list.push(item);
						}
					}
				};
				return list;
			};

			/***/
		}
		/******/])
	);
});
;


},{}],9:[function(require,module,exports){
'use strict';

var _typeof = typeof Symbol === "function" && typeof Symbol.iterator === "symbol" ? function (obj) { return typeof obj; } : function (obj) { return obj && typeof Symbol === "function" && obj.constructor === Symbol && obj !== Symbol.prototype ? "symbol" : typeof obj; };

//download.js v4.2, by dandavis; 2008-2016. [CCBY2] see http://danml.com/download.html for tests/usage
// v1 landed a FF+Chrome compat way of downloading strings to local un-named files, upgraded to use a hidden frame and optional mime
// v2 added named files via a[download], msSaveBlob, IE (10+) support, and window.URL support for larger+faster saves than dataURLs
// v3 added dataURL and Blob Input, bind-toggle arity, and legacy dataURL fallback was improved with force-download mime and base64 support. 3.1 improved safari handling.
// v4 adds AMD/UMD, commonJS, and plain browser support
// v4.1 adds url download capability via solo URL argument (same domain/CORS only)
// v4.2 adds semantic variable names, long (over 2MB) dataURL support, and hidden by default temp anchors
// https://github.com/rndme/download

(function (root, factory) {
	if (typeof define === 'function' && define.amd) {
		// AMD. Register as an anonymous module.
		define([], factory);
	} else if ((typeof exports === 'undefined' ? 'undefined' : _typeof(exports)) === 'object') {
		// Node. Does not work with strict CommonJS, but
		// only CommonJS-like environments that support module.exports,
		// like Node.
		module.exports = factory();
	} else {
		// Browser globals (root is window)
		root.download = factory();
	}
})(undefined, function () {

	return function download(data, strFileName, strMimeType) {

		var self = window,
		    // this script is only for browsers anyway...
		defaultMime = "application/octet-stream",
		    // this default mime also triggers iframe downloads
		mimeType = strMimeType || defaultMime,
		    payload = data,
		    url = !strFileName && !strMimeType && payload,
		    anchor = document.createElement("a"),
		    toString = function toString(a) {
			return String(a);
		},
		    myBlob = self.Blob || self.MozBlob || self.WebKitBlob || toString,
		    fileName = strFileName || "download",
		    blob,
		    reader;
		myBlob = myBlob.call ? myBlob.bind(self) : Blob;

		if (String(this) === "true") {
			//reverse arguments, allowing download.bind(true, "text/xml", "export.xml") to act as a callback
			payload = [payload, mimeType];
			mimeType = payload[0];
			payload = payload[1];
		}

		if (url && url.length < 2048) {
			// if no filename and no mime, assume a url was passed as the only argument
			fileName = url.split("/").pop().split("?")[0];
			anchor.href = url; // assign href prop to temp anchor
			if (anchor.href.indexOf(url) !== -1) {
				// if the browser determines that it's a potentially valid url path:
				var ajax = new XMLHttpRequest();
				ajax.open("GET", url, true);
				ajax.responseType = 'blob';
				ajax.onload = function (e) {
					download(e.target.response, fileName, defaultMime);
				};
				setTimeout(function () {
					ajax.send();
				}, 0); // allows setting custom ajax headers using the return:
				return ajax;
			} // end if valid url?
		} // end if url?


		//go ahead and download dataURLs right away
		if (/^data\:[\w+\-]+\/[\w+\-]+[,;]/.test(payload)) {

			if (payload.length > 1024 * 1024 * 1.999 && myBlob !== toString) {
				payload = dataUrlToBlob(payload);
				mimeType = payload.type || defaultMime;
			} else {
				return navigator.msSaveBlob ? // IE10 can't do a[download], only Blobs:
				navigator.msSaveBlob(dataUrlToBlob(payload), fileName) : saver(payload); // everyone else can save dataURLs un-processed
			}
		} //end if dataURL passed?

		blob = payload instanceof myBlob ? payload : new myBlob([payload], { type: mimeType });

		function dataUrlToBlob(strUrl) {
			var parts = strUrl.split(/[:;,]/),
			    type = parts[1],
			    decoder = parts[2] == "base64" ? atob : decodeURIComponent,
			    binData = decoder(parts.pop()),
			    mx = binData.length,
			    i = 0,
			    uiArr = new Uint8Array(mx);

			for (i; i < mx; ++i) {
				uiArr[i] = binData.charCodeAt(i);
			}return new myBlob([uiArr], { type: type });
		}

		function saver(url, winMode) {

			if ('download' in anchor) {
				//html5 A[download]
				anchor.href = url;
				anchor.setAttribute("download", fileName);
				anchor.className = "download-js-link";
				anchor.innerHTML = "downloading...";
				anchor.style.display = "none";
				document.body.appendChild(anchor);
				setTimeout(function () {
					anchor.click();
					document.body.removeChild(anchor);
					if (winMode === true) {
						setTimeout(function () {
							self.URL.revokeObjectURL(anchor.href);
						}, 250);
					}
				}, 66);
				return true;
			}

			// handle non-a[download] safari as best we can:
			if (/(Version)\/(\d+)\.(\d+)(?:\.(\d+))?.*Safari\//.test(navigator.userAgent)) {
				url = url.replace(/^data:([\w\/\-\+]+)/, defaultMime);
				if (!window.open(url)) {
					// popup blocked, offer direct download:
					if (confirm("Displaying New Document\n\nUse Save As... to download, then click back to return to this page.")) {
						location.href = url;
					}
				}
				return true;
			}

			//do iframe dataURL download (old ch+FF):
			var f = document.createElement("iframe");
			document.body.appendChild(f);

			if (!winMode) {
				// force a mime that will download:
				url = "data:" + url.replace(/^data:([\w\/\-\+]+)/, defaultMime);
			}
			f.src = url;
			setTimeout(function () {
				document.body.removeChild(f);
			}, 333);
		} //end saver


		if (navigator.msSaveBlob) {
			// IE10+ : (has Blob, but not a[download] or URL)
			return navigator.msSaveBlob(blob, fileName);
		}

		if (self.URL) {
			// simple fast and modern way using Blob and URL:
			saver(self.URL.createObjectURL(blob), true);
		} else {
			// handle non-Blob()+non-URL browsers:
			if (typeof blob === "string" || blob.constructor === toString) {
				try {
					return saver("data:" + mimeType + ";base64," + self.btoa(blob));
				} catch (y) {
					return saver("data:" + mimeType + "," + encodeURIComponent(blob));
				}
			}

			// Blob but not URL support:
			reader = new FileReader();
			reader.onload = function (e) {
				saver(this.result);
			};
			reader.readAsDataURL(blob);
		}
		return true;
	}; /* end download() */
});

},{}],10:[function(require,module,exports){
'use strict';

var _PaintControl = require('./Controls/PaintControl');

var _PaintControl2 = _interopRequireDefault(_PaintControl);

var _svg = require('svg.js');

var _svg2 = _interopRequireDefault(_svg);

var _UIManagement = require('./model/UIManagement');

var UI = _interopRequireWildcard(_UIManagement);

function _interopRequireWildcard(obj) { if (obj && obj.__esModule) { return obj; } else { var newObj = {}; if (obj != null) { for (var key in obj) { if (Object.prototype.hasOwnProperty.call(obj, key)) newObj[key] = obj[key]; } } newObj.default = obj; return newObj; } }

function _interopRequireDefault(obj) { return obj && obj.__esModule ? obj : { default: obj }; }

(function () {
	var draw = (0, _svg2.default)('drawing').size(1000, 1000);
	UI.setGUI();
	setControl(draw);
})();

function setControl(_container) {
	var isMouseDown = false;
	var currnetControl = new _PaintControl2.default(_container);
	var top = _container.node.getBoundingClientRect().top;
	var left = _container.node.getBoundingClientRect().left;

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

},{"./Controls/PaintControl":7,"./model/UIManagement":14,"svg.js":6}],11:[function(require,module,exports){
"use strict";

Object.defineProperty(exports, "__esModule", {
	value: true
});
exports.default = {
	selectedCurve: []
};

},{}],12:[function(require,module,exports){
'use strict';

Object.defineProperty(exports, "__esModule", {
	value: true
});

var _createClass = function () { function defineProperties(target, props) { for (var i = 0; i < props.length; i++) { var descriptor = props[i]; descriptor.enumerable = descriptor.enumerable || false; descriptor.configurable = true; if ("value" in descriptor) descriptor.writable = true; Object.defineProperty(target, descriptor.key, descriptor); } } return function (Constructor, protoProps, staticProps) { if (protoProps) defineProperties(Constructor.prototype, protoProps); if (staticProps) defineProperties(Constructor, staticProps); return Constructor; }; }();

var _bezierJs = require('bezier-js');

var _bezierJs2 = _interopRequireDefault(_bezierJs);

var _MagneticCurve = require('../model/MagneticCurve');

var _MagneticCurve2 = _interopRequireDefault(_MagneticCurve);

var _CurveManagement = require('./CurveManagement');

var _CurveManagement2 = _interopRequireDefault(_CurveManagement);

var _UIManagement = require('../model/UIManagement');

var UI = _interopRequireWildcard(_UIManagement);

function _interopRequireWildcard(obj) { if (obj && obj.__esModule) { return obj; } else { var newObj = {}; if (obj != null) { for (var key in obj) { if (Object.prototype.hasOwnProperty.call(obj, key)) newObj[key] = obj[key]; } } newObj.default = obj; return newObj; } }

function _interopRequireDefault(obj) { return obj && obj.__esModule ? obj : { default: obj }; }

function _classCallCheck(instance, Constructor) { if (!(instance instanceof Constructor)) { throw new TypeError("Cannot call a class as a function"); } }

var LevelCurve = function () {
	/**
 	@param levelParam : array
 	length : number
 	alpha : number
 	branches : number
 */
	function LevelCurve(basePath, trunkWidth, levelParam) {
		_classCallCheck(this, LevelCurve);

		this.basePath = basePath;
		this.trunkWidth = trunkWidth;
		this.levelParam = levelParam;
		this.curveGroup = undefined;
	}

	_createClass(LevelCurve, [{
		key: 'drawLevelCurve',
		value: function drawLevelCurve(beziers, level) {
			var _this = this;

			if (!beziers) return;
			var sign = 1;

			var Bs = beziers.map(function (b) {
				return new _bezierJs2.default(b[0][0], b[0][1], b[1][0], b[1][1], b[2][0], b[2][1], b[3][0], b[3][1]);
			});

			var totalLength = Bs.reduce(function (length, B) {
				return B.length() + length;
			}, 0);

			this.branchPosition(level).forEach(function (i) {
				if (totalLength === 0) return;

				var bezierIndex = 0;

				var pos = totalLength * i;
				while (pos >= Bs[bezierIndex].length()) {
					pos -= Bs[bezierIndex].length();
					bezierIndex++;
				}

				var posOnSinglebezier = pos / Bs[bezierIndex].length();

				_this.drawAt(posOnSinglebezier, Bs[bezierIndex], sign, level);

				sign *= -1;
			});
		}
	}, {
		key: 'drawAt',
		value: function drawAt(t, b, sign, level) {
			var start = b.get(t);
			var v = b.derivative(t);
			var mag = new _MagneticCurve2.default({
				startX: start.x,
				startY: start.y,
				vx: v.x,
				vy: v.y,
				T: this.levelParam[level].length,
				alpha: this.levelParam[level].alpha,
				sign: sign
			});
			mag.drawOn(this.curveGroup);
			if (level < this.levelParam.length - 1) this.drawLevelCurve(mag.points, level + 1);
		}
	}, {
		key: 'drawOn',
		value: function drawOn(pannel) {
			this.pannel = pannel;
			this.curveGroup = pannel.group();
			pannel.add(this.curveGroup);
			_CurveManagement2.default[this.curveGroup.node.id] = this;
			this.curveGroup.on('click', function (e) {
				console.log('clecked');

				var curve_id = e.target.parentElement.id;
				var lvCurve = _CurveManagement2.default[curve_id];
				_CurveManagement2.default.selectedCurve.length = 0;
				_CurveManagement2.default.selectedCurve.push(lvCurve);
			});

			this.drawLevelCurve(this.basePath, 0);
			this.drawStem();
		}
	}, {
		key: 'drawStem',
		value: function drawStem() {
			var _this2 = this;

			var basePath = this.basePath.map(function (c) {
				return new _bezierJs2.default(c[0][0], c[0][1], c[1][0], c[1][1], c[2][0], c[2][1], c[3][0], c[3][1]);
			});

			var totalLength = basePath.reduce(function (length, B) {
				return B.length() + length;
			}, 0);

			var l = 0;
			var w = UI.state.trunkWidth;

			var outline = basePath.map(function (b) {
				var d1 = w * l / totalLength;
				l += b.length();
				var d2 = w * l / totalLength;
				if (d1 === 0) d1 = 1;
				return b.outline(d1, d1, d2, d2);
			});
			outline.forEach(function (b) {
				var pathString = fittedCurveToPathString(b);
				drawOnPannel(_this2.curveGroup, pathString);
			});
		}
	}, {
		key: 'redraw',
		value: function redraw() {
			if (this.pannel === undefined) {
				console.error('can not redraw!');
				return;
			}
			this.curveGroup.remove();
			this.drawOn(this.pannel);
		}
	}, {
		key: 'branchPosition',
		value: function branchPosition(level) {
			var branches = [];
			var branch = this.levelParam[level].branches;

			var pos = 1 / (branch + 1);
			for (var i = 1; i <= branch; i++) {
				branches.push(pos * i);
			}

			return branches;
		}
	}]);

	return LevelCurve;
}();

exports.default = LevelCurve;

function fittedCurveToPathString(fittedLineData) {
	var str = '';
	//bezier : [ [c0], [c1], [c2], [c3] ]
	fittedLineData.curves.map(function (bezier, i) {
		if (i == 0) {
			str += 'M ' + bezier.points[0].x + ' ' + bezier.points[0].y;
		}

		str += 'C ' + bezier.points[1].x + ' ' + bezier.points[1].y + ', ' + bezier.points[2].x + ' ' + bezier.points[2].y + ', ' + bezier.points[3].x + ' ' + bezier.points[3].y + ' ';
	});

	return str;
}
function drawOnPannel(pannel, pathString) {
	pannel.path(pathString).fill('#bec9b3').stroke({ width: 0 });
}

},{"../model/MagneticCurve":13,"../model/UIManagement":14,"./CurveManagement":11,"bezier-js":1}],13:[function(require,module,exports){
'use strict';

Object.defineProperty(exports, "__esModule", {
	value: true
});

var _slicedToArray = function () { function sliceIterator(arr, i) { var _arr = []; var _n = true; var _d = false; var _e = undefined; try { for (var _i = arr[Symbol.iterator](), _s; !(_n = (_s = _i.next()).done); _n = true) { _arr.push(_s.value); if (i && _arr.length === i) break; } } catch (err) { _d = true; _e = err; } finally { try { if (!_n && _i["return"]) _i["return"](); } finally { if (_d) throw _e; } } return _arr; } return function (arr, i) { if (Array.isArray(arr)) { return arr; } else if (Symbol.iterator in Object(arr)) { return sliceIterator(arr, i); } else { throw new TypeError("Invalid attempt to destructure non-iterable instance"); } }; }();

var _createClass = function () { function defineProperties(target, props) { for (var i = 0; i < props.length; i++) { var descriptor = props[i]; descriptor.enumerable = descriptor.enumerable || false; descriptor.configurable = true; if ("value" in descriptor) descriptor.writable = true; Object.defineProperty(target, descriptor.key, descriptor); } } return function (Constructor, protoProps, staticProps) { if (protoProps) defineProperties(Constructor.prototype, protoProps); if (staticProps) defineProperties(Constructor, staticProps); return Constructor; }; }();

var _fitCurve = require('fit-curve');

var _fitCurve2 = _interopRequireDefault(_fitCurve);

function _interopRequireDefault(obj) { return obj && obj.__esModule ? obj : { default: obj }; }

function _classCallCheck(instance, Constructor) { if (!(instance instanceof Constructor)) { throw new TypeError("Cannot call a class as a function"); } }

var error = 1;

var MagneticCurve = function () {
	function MagneticCurve(param) {
		_classCallCheck(this, MagneticCurve);

		/*
  startX, startY : 初始點
  vx, vy : 初始速度
  T : 總點數
  alpha : 等角螺線參數
  sign : 電荷正負
  */
		this.param = param;
	}

	_createClass(MagneticCurve, [{
		key: 'makeCurve',
		value: function makeCurve() {
			var points = [];

			var sign = this.param.sign || 1;
			var x = this.param.startX;
			var y = this.param.startY;

			var _normalize = normalize([this.param.vx, this.param.vy]),
			    _normalize2 = _slicedToArray(_normalize, 2),
			    vx = _normalize2[0],
			    vy = _normalize2[1];

			var T = this.param.T;
			var t = 0;
			var dt = 1;

			while (t < T) {
				points.push([x, y]);
				var q = sign * Math.pow(T - t, -1 * this.param.alpha);
				x += vx * dt;
				y += vy * dt;

				var ax = -1 * vy * q;
				var ay = vx * q;

				vx += ax * dt;
				vy += ay * dt;

				t += dt;
			}

			return points;
		}
	}, {
		key: 'drawOn',
		value: function drawOn(pannel) {
			var mag = this.makeCurve();
			var smoothBizer = (0, _fitCurve2.default)(mag, error);
			this.points = smoothBizer;

			var pathString = fittedCurveToPathString(smoothBizer);

			pannel.path(pathString).fill('none').stroke({ width: 3 }).stroke('#f00');
		}
	}]);

	return MagneticCurve;
}();

exports.default = MagneticCurve;

function fittedCurveToPathString(fittedLineData) {
	var str = '';
	//bezier : [ [c0], [c1], [c2], [c3] ]
	fittedLineData.forEach(function (bezier, i) {
		if (i == 0) {
			str += 'M ' + bezier[0][0] + ' ' + bezier[0][1];
		}

		str += 'C ' + bezier[1][0] + ' ' + bezier[1][1] + ', ' + bezier[2][0] + ' ' + bezier[2][1] + ', ' + bezier[3][0] + ' ' + bezier[3][1] + ' ';
	});

	return str;
}

function normalize(vector) {
	var x = vector[0];
	var y = vector[1];
	var length = Math.sqrt(Math.pow(x, 2) + Math.pow(y, 2));

	return [x / length, y / length];
}

},{"fit-curve":5}],14:[function(require,module,exports){
'use strict';

Object.defineProperty(exports, "__esModule", {
	value: true
});
exports.state = undefined;
exports.setGUI = setGUI;

var _dat = require('../lib/dat.gui');

var dat = _interopRequireWildcard(_dat);

var _download2 = require('../lib/download');

var _download3 = _interopRequireDefault(_download2);

var _CurveManagement = require('./CurveManagement');

var _CurveManagement2 = _interopRequireDefault(_CurveManagement);

function _interopRequireDefault(obj) { return obj && obj.__esModule ? obj : { default: obj }; }

function _interopRequireWildcard(obj) { if (obj && obj.__esModule) { return obj; } else { var newObj = {}; if (obj != null) { for (var key in obj) { if (Object.prototype.hasOwnProperty.call(obj, key)) newObj[key] = obj[key]; } } newObj.default = obj; return newObj; } }

var gui = void 0,
    folders = [];
var controls = [];
var state = exports.state = {
	trunkWidth: 15,
	intersect: false,
	levelCurve: [{
		length: 100,
		alpha: 0.9,
		branches: 5
	}, {
		length: 20,
		alpha: 0.8,
		branches: 5
	}]
};

var features = {
	download: function download() {
		var svg = document.getElementsByTagName('svg')[0];
		(0, _download3.default)(svg.outerHTML, 'file.svg', 'text/plain');
	}
};

function setGUI() {
	gui = new dat.GUI();

	var c = gui.add(state, 'trunkWidth', 1, 20);
	controls.push(c);
	//gui.add(state, 'intersect');

	levelFolder(0);
	levelFolder(1);
	setOnChange(controls);

	gui.add(features, 'download');
}

function levelFolder(index) {
	var folder = gui.addFolder('Level ' + index);
	controls.push(folder.add(state.levelCurve[index], 'length'));
	controls.push(folder.add(state.levelCurve[index], 'alpha'));
	controls.push(folder.add(state.levelCurve[index], 'branches').step(1));
	folders.push(folder);
}

function setOnChange(controls) {
	controls.forEach(function (c) {
		c.onChange(function () {
			if (_CurveManagement2.default.selectedCurve.length === 1) {
				_CurveManagement2.default.selectedCurve[0].redraw();
			}
		});
	});
}

},{"../lib/dat.gui":8,"../lib/download":9,"./CurveManagement":11}]},{},[10])