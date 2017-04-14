import fitCurve from 'fit-curve';
import Bezier from 'bezier-js';

export class BezierSpline {
	constructor(points, error){
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

		this.colliders = this.controlPoints.reduce((acc, val)=>acc.concat(val),[]);
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