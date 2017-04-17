import fitCurve from 'fit-curve';
import Bezier from 'bezier-js';
let glm = require('glm-js');

export class BezierSpline {
	constructor(points, error=50){
		this.points = points;
		this.controlPoints = fitCurve(points, error);//[ [c0,c1,c2,c3], ... ]
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
			let Bs = this.beziers;

			let totalLength = this.length;
			let burgeons = [];

			let pos = 1 / (density + 1);
			for (let i = 1; i <= density; i++) {
				burgeons.push(pos * i);
			}
			return burgeons.map(i => {
				if( totalLength === 0) return;

				let bezierIndex = 0;

				let pos = totalLength * i;
				while( pos >= Bs[bezierIndex].length() ) {
					pos -= Bs[bezierIndex].length();
					bezierIndex++;
				}
				let bezierAtIndex = Bs[bezierIndex];
				let posOnSinglebezier = pos / bezierAtIndex.length();

				let bezier = Bs[ bezierIndex ];
				let point = bezier.get( posOnSinglebezier );
				let direction = bezier.derivative( posOnSinglebezier );

				return {point, direction};
			});
		};

		this.colliders = this.controlPoints.reduce((acc, val)=>{
			let collider = makeCollider(val);
			return acc.concat( collider );
		},[]);
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
function makeCollider(controlPoint){
	let c0 = glm.vec2(controlPoint[0][0], controlPoint[0][1]);
	let c1 = glm.vec2(controlPoint[1][0], controlPoint[1][1]);
	let c2 = glm.vec2(controlPoint[2][0], controlPoint[2][1]);
	let c3 = glm.vec2(controlPoint[3][0], controlPoint[3][1]);

	let proj_C1 = reflecttPoint(c0,c1,c3 );
	let proj_C2 = reflecttPoint(c0,c2,c3 );

	let collider = controlPoint.concat([proj_C2, proj_C1]);
	return collider;
}
function reflecttPoint(p0,p1,p2){
	let a = glm.sub(p1,p0);
	let b = glm.normalize(glm.sub(p2,p0));
	let proj_a_b = glm.mul( b , glm.dot(a, b) );
	let reflect = glm.sub( glm.mul(proj_a_b,2),a );
	let reflect_p = glm.add(p0, reflect);
	return [ reflect_p.x, reflect_p.y ];
}