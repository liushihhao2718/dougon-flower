import fitCurve from 'fit-curve';
import Bezier from 'bezier-js';
let glm = require('glm-js');
import _ from 'lodash';
export class BezierSpline {
	static makeByPoints(points, error=50){
		let controlPoints = fitCurve(points, error);//[ [c0,c1,c2,c3], ... ]
		return new BezierSpline(controlPoints);
	}
	constructor(controlPoints){
		this.controlPoints = controlPoints;//[ [c0,c1,c2,c3], ... ]
		this.beziers = this.controlPoints.map(b =>	
			new Bezier(
				b[0][0], b[0][1],
				b[1][0], b[1][1],
				b[2][0], b[2][1],
				b[3][0], b[3][1]
			)
		);
		this.lengths = this.beziers.map(b=>b.length());
		this.length = this.lengths.reduce( (a,b)=>a+b , 0);
		this.sample = function(density){
			let burgeons = [];
			let pos = 1 / (density + 1);
			for (let i = 1; i <= density; i++) {
				burgeons.push(pos * i);
			}
			return burgeons.map(x=>this.test(x));
		};
		
		this.colliders = this.controlPoints.reduce((acc, val)=>{
			let collider = makeCollider(val);
			return acc.concat( collider );
		},[]);
	}
	paramColliders(p=0.6){
		return this.controlPoints.reduce((acc, val)=>{
			let collider = makeCollider(val, p);
			return acc.concat( collider );
		},[]);
	}
	outline(beginWidth, endWidth){
		let l = 0;
		let w = (endWidth-beginWidth)/this.length;
		let outline = this.beziers.map(b => {
			let d1 = beginWidth + l * w;
			l += b.length();
			let d2 = beginWidth + l * w;
			if(d1 === 0) d1 = 1;
			return b.outline(d1, d1, d2, d2);
		});

		let f = [];
		let r = [];
		let head=outline[0].curves[0];
		let tail;
		const lastPath = outline[outline.length-1];
		tail = lastPath.curves[lastPath.curves.length/2];

		outline.forEach(b => {
			const length = Math.floor(b.curves.length/2)- 1;
			f = f.concat(b.curves.slice(1, length+1));
			r = r.concat(b.curves.slice(length+2).reverse());
		});

		return [].concat([head],f,[tail], r.reverse());
	}
	test(x){
		
		let {bezierIndex, posOnSinglebezier} = this.sampleAt(x);
		let {point, direction} = this.toRay(this.beziers[bezierIndex], posOnSinglebezier);
		return {point, direction};

	}
	sampleAt(index){

		if (index === 1) return {bezierIndex:this.beziers.length-1, posOnSinglebezier:1};
		let Bs = this.beziers;
		let totalLength = this.length;

		if( totalLength === 0) return;

		let bezierIndex = 0;
		let pos = totalLength * index;
		while( pos >= Bs[bezierIndex].length() ) {
			pos -= Bs[bezierIndex].length();
			bezierIndex++;
		}
		let posOnSinglebezier = pos / Bs[bezierIndex].length();

		return {bezierIndex, posOnSinglebezier};
	}
	toRay(bezierAtIndex, posOnSinglebezier){
		let point =  bezierAtIndex.get( posOnSinglebezier );
		let direction = bezierAtIndex.derivative( posOnSinglebezier );

		return {point, direction};
	}
	segmentRange(i, j){
		if(i === j) return;
		let start = this.sampleAt(i);
		let end = this.sampleAt(j);

		let segments = this.beziers.slice(start.bezierIndex,end.bezierIndex+1);
		if (segments.length == 1) {
			segments[0] = segments[0].split(start.posOnSinglebezier, end.posOnSinglebezier);
		}
		else{
			segments[0] = segments[0].split(start.posOnSinglebezier, 1);
			segments[segments.length-1] = segments[segments.length-1].split(0, end.posOnSinglebezier);
		}
		return segments;
	}

	range(i,j){
		let segments = this.segmentRange(i,j);
		let controlPoints = segments.map(b=>b.points.map(p=>[p.x,p.y]));
		return new BezierSpline(controlPoints);
	}
	svgString() {
		var str = '';
		//bezier : [ [c0], [c1], [c2], [c3] ]
		this.controlPoints.map(function (bezier, i) {
			if (i == 0) {
				str += 'M ' + bezier[0][0] + ' ' + bezier[0][1];
			}

			str += 'C ' + bezier[1][0] + ' ' + bezier[1][1] + ', ' +
			bezier[2][0] + ' ' + bezier[2][1] + ', ' +
			bezier[3][0] + ' ' + bezier[3][1] + ' ';	
						
		});

		return str;
	}
}
//this.controlPoints.map(a=>a.map(b=>b.map(x=>Math.round(x)).join(',')).join('    ')).join(' : ')
function makeCollider(controlPoint, long){
	let c0 = glm.vec2(controlPoint[0][0], controlPoint[0][1]);
	let c1 = glm.vec2(controlPoint[1][0], controlPoint[1][1]);
	let c2 = glm.vec2(controlPoint[2][0], controlPoint[2][1]);
	let c3 = glm.vec2(controlPoint[3][0], controlPoint[3][1]);

	let reflect_C1 = reflecttPoint(c0,c1,c3 );
	let reflect_C2 = reflecttPoint(c0,c2,c3 );

	let collider = controlPoint.concat([reflect_C2, reflect_C1]);
	[collider[1],collider[5]]=longer(collider[1],collider[5], long);
	[collider[2],collider[4]]=longer(collider[2],collider[4], long);
	return collider;
}
function longer(p0,p1, long = 0.6){
	let v0 = glm.vec2(p0[0], p0[1]);
	let v1 = glm.vec2(p1[0], p1[1]);
	let a = glm.mul(glm.sub(v0,v1),long);
	let b = glm.mul(glm.add(v0,v1),0.5);
	return [p0,p1]=[glm.add(b,a).elements,glm.sub(b,a).elements];
}
function reflecttPoint(p0,p1,p2){
	let a = glm.sub(p1,p0);
	let b = glm.normalize(glm.sub(p2,p0));
	let proj_a_b = glm.mul( b , glm.dot(a, b) );
	let reflect = glm.sub( glm.mul(proj_a_b,2),a );
	let reflect_p = glm.add(p0, reflect);
	return [ reflect_p.x, reflect_p.y ];
}