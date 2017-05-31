import fitCurve from 'fit-curve';
import Bezier from 'bezier-js';

const error = 20;
// let timeout = 20;
export default class MagneticCurve {
	constructor(param) {
		/*
		startX, startY : 初始點
		vx, vy : 初始速度
		T : 總點數
		alpha : 等角螺線參數
		sign : 電荷正負
		level : curve lv
		*/
		this.param = param;
	}

	sample(number){
		let Bs = this.controlPoints.map(b =>{
			return new Bezier(
				b[0][0], b[0][1],
				b[1][0], b[1][1],
				b[2][0], b[2][1],
				b[3][0], b[3][1]
			);
		});
		let totalLength = Bs.reduce( (length, B) =>{
			return B.length() + length;
		}, 0);
		let branches = [];
		let pos = 1 / (number + 1);
		for (let i = 1; i <= number; i++) {
			branches.push(pos * i);
		}
		let sample = [];
		branches.forEach( i => {
			if( totalLength === 0) return;

			let bezierIndex = 0;

			let pos = totalLength * i;
			while( pos >= Bs[bezierIndex].length() ) {
				pos -= Bs[bezierIndex].length();
				bezierIndex++;
			}
			let bezierAtIndex = Bs[bezierIndex];
			let posOnSinglebezier = pos / bezierAtIndex.length();

			sample.push( bezierAtIndex.compute(posOnSinglebezier) );// [x, y]
		});
		return sample;
	}
	makeCurve() {
		const points = [];
		
		let sign = this.param.sign || 1;
		let x = this.param.startX;
		let y = this.param.startY;

		let [vx, vy] = normalize( [this.param.vx, this.param.vy] );

		let T = this.param.T;
		let t = 0;
		let dt = T/7;

		while( t < T ) {
			points.push( [x, y] );
			let q = sign * Math.pow( (T-t), -1 * this.param.alpha );
			x += vx * dt;
			y += vy * dt;

			let ax = -1 * vy * q;
			let ay = vx * q;

			vx += ax * dt;
			vy += ay * dt;

			t += dt;
		}

		return points;
	}
	getCurve() {
		if( !this._curve ){
			this._curve = this.makeCurve();
		} 
		return this._curve;
	}
	get controlPoints(){
		if(!this._controlPoints)
			this._controlPoints = fitCurve( this._curve, error );		
		return this._controlPoints;
	}

	makeBBox(points){
		let minX = Number.MAX_VALUE;
		let minY = Number.MAX_VALUE;
		let maxX = Number.MIN_VALUE;
		let maxY = Number.MIN_VALUE;
		points.forEach(p => {
			let x = p[0];
			let y = p[1];

			if(x < minX) minX = x;
			if(y < minY) minY = y;
			if(x > maxX) maxX = x;
			if(y > maxY) maxY = y;
		});

		return {
			x: minX,
			y: minY,
			width: maxX - minX,
			height: maxY - minY
		};
	}
	bbox(){
		if(!this._bbox) this._bbox = this.makeBBox( this.getCurve() );
		return this._bbox;
	}
}
function fittedCurveToPathString(fittedLineData) {
	var str = '';
	//bezier : [ [c0], [c1], [c2], [c3] ]
	fittedLineData.forEach(function (bezier, i) {
		if (i == 0) {
			str += 'M ' + bezier[0][0] + ' ' + bezier[0][1];
		}

		str += 'C ' + bezier[1][0] + ' ' + bezier[1][1] + ', ' +
		bezier[2][0] + ' ' + bezier[2][1] + ', ' +
		bezier[3][0] + ' ' + bezier[3][1] + ' ';	
					
	});

	return str;
}

function normalize(vector) {
	let x = vector[0];
	let y = vector[1];
	let length = Math.sqrt(
		Math.pow(x, 2) + Math.pow(y, 2)
	);

	return [x/length, y/length];
}