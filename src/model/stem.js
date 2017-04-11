import fitCurve from 'fit-curve';
import MagneticCurve from '../model/MagneticCurve';

const error = 10;
// let timeout = 20;
export class Burgeon{
	constructor(x,y,direction){
		this.x=x;
		this.y=y;
		this.direction = direction;
		this.parent = parent;
	}

	/*
		@param length : number
		the length of generated Stem

		@param density : number
		number of Stem to generated (*both side total)

		@return : Stem[]
	*/
	germinate(length,density){
		
	}
}

export class Stem {
	constructor(startX,startY,vx,vy,alpha,sign,length) {
		this.curve = new MagneticCurve({
			startX,
			startY,
			vx,
			vy,
			T: length,
			alpha,
			sign
		});

	}
	/*
		@param segment : number 
		number of branch
	*/
	*Burgeons(segment){
		//currently generate from root
		//inverse the iterator order to generate from tail
		let sample = this.curve.sample(segment);
		for (var i = sample.length - 1; i >= 0; i--) {
			yield sample[i];
		}
	}
	
	get boundingBox(){
		delete this.boundingBox;
		return this.boundingBox = makeBBox( this.getCurve() );
	}
}

function makeBBox(points){
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
function branchPosition(level){
	let branches = [];
	let branch = this.levelParam[level].branches;

	let pos = 1 / (branch + 1);
	for (let i = 1; i <= branch; i++) {
		branches.push(pos * i);
	}

	return branches;
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