import fitCurve from 'fit-curve';
import Bezier from 'bezier-js';
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
			burgeons.reverse();
			return burgeons.map(x=>this.test(x));
		};
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