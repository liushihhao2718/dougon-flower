import fitCurve from 'fit-curve';

const error = 1;


export default class MagneticCurve {
	constructor(param) {
		/*
		startX, startY : 初始點
		vx, vy : 初始速度
		T : 總點數
		alpha : 等角螺線參數
		sign : 電荷正負
		*/
		this.param = param;
	}

	makeCurve() {
		const points = [];
		
		let sign = this.param.sign || 1;
		let x = this.param.startX;
		let y = this.param.startY;

		let [vx, vy] = normalize( [this.param.vx, this.param.vy] );

		let T = this.param.T;
		let t = 0;
		let dt = 1;

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

	drawOn(pannel){
		let mag = this.makeCurve();
		let smoothBizer = fitCurve( mag, error );
		this.points = smoothBizer;

		let pathString = fittedCurveToPathString(smoothBizer);

		// pannel.path(pathString).fill('none').stroke({ width: 3 }).stroke('#f00');
		pannel.path(pathString).fill('none').stroke({ width: 10 }).stroke('#CED5D0');
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