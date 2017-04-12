import MagneticCurve from '../model/MagneticCurve';
import * as Collision from '../model/Collision';
import * as UI from '../model/UIManagement';
import CurveManagement from './CurveManagement';
let shortid = require('shortid');

export class Floral{
	constructor(basePath, flowerType='海石榴華', flowerRotation=45){
		this.id = shortid.generate();
		this.basePath = basePath;
		this.flowerType = flowerType;
		this.flowerRotation = flowerRotation;

		let points = this.basePath.points;
		this.flowerPosition = {
			x: points[points.length-1][0],
			y: points[points.length-1][1]
		};
	}
	burgeons(amount){
		let samples = this.basePath.sample(amount);
		let burgeons = [];
		samples.forEach(s => {
			const {posOnSinglebezier, bezierIndex} = s;
			let bezier = this.basePath[ bezierIndex ];
			let point = bezier.get( posOnSinglebezier );
			let direction = bezier.derivative( posOnSinglebezier );

			burgeons.push( new Burgeon(point[0], point[1], direction) );
		});
	}
}
export class Burgeon{
	constructor(x,y,direction, parent){
		this.x=x;
		this.y=y;
		this.direction = direction;
		this.parent = parent;
	}

	germinate(level, max, sign){
		if(level === max) return;

		let stem = new Stem(this.x, this.y, 
			this.direction.x, this.direction.y,
			level,
			sign
		);

		if( Collision.testCollision(stem.boundingBox, CurveManagement.leafCollisionScene, this.parent )){
			return this.germinate(level+1 , max, sign);
		}else{
			return stem;
		}
	}
}

export class Stem {
	constructor(startX,startY,vx,vy,level,sign) {
		const levelParam = UI.state.levelCurve[level];

		this.curve = new MagneticCurve({
			startX,
			startY,
			vx,
			vy,
			T: levelParam.length,
			alpha:levelParam.alpha,
			sign
		});
	}
	/*
		@param segment : number 
		number of branch
	*/
	burgeons(segment){
		//currently generate from root
		//inverse the iterator order to generate from tail
		let sample = this.curve.sample(segment);
		// for (var i = sample.length - 1; i >= 0; i--) {
		// 	yield sample[i];
		// }
		return sample;
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
