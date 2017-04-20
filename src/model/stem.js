import MagneticCurve from './MagneticCurve';
import { BezierSpline } from './Spline';
import shortid from 'shortid';
import nextType from './MarkovLeaf';

export class Floral{
	constructor(basePath, r,trunkHead,trunkTail,flowerType='海石榴華', aspect = '正面',flowerRotation=30){
		this.id = shortid.generate();
		this.curve = basePath;
		this.flowerType = flowerType;
		this.flowerRotation = flowerRotation;
		this.colliders = undefined;
		let points = this.curve.points;
		this.flowerPosition = {
			x: points[points.length-1][0],
			y: points[points.length-1][1],
			r
		};
		this.trunkHead = trunkHead;
		this.trunkTail = trunkTail;
		this.aspect = aspect;
	}
	burgeons(amount){
		if (this.curve.length === 0) return;
		let samples = this.curve.sample(amount);
		return samples.map(s => {
			const {point, direction} = s;
			return (new Burgeon(point.x, point.y, direction, 0, this) );
		});
	}
}

export class Burgeon{
	constructor(x,y,direction,type,parent){
		this.x=x;
		this.y=y;
		this.direction = direction;
		this.parent = parent;
		this.type = type;
	}

	germinate( length,alpha,sign){
		let leaf = new Leaf(this.x, this.y, 
			this.direction.x, this.direction.y,
			length,
			alpha,
			sign,
			this.type
		);

		return leaf;
	}
}

export class Leaf {
	constructor(startX,startY,vx,vy,length, alpha, sign,type) {
		let mag = new MagneticCurve({
			startX,
			startY,
			vx,
			vy,
			T: length,
			alpha,
			sign
		});
		this.startX = startX;
		this.startY = startY;
		const lastOne = mag.getCurve().length -1;
		this.endX = mag.getCurve()[lastOne][0];
		this.endY = mag.getCurve()[lastOne][1];
		this.type = type;
		this.curve = new BezierSpline(mag.getCurve() );
		this.colliders = this.curve.colliders;
		this.sign = sign;

		this.mag = mag;
	}
	/*
		@param amount : number 
		number of branch
	*/
	burgeons(amount){
		let samples = this.curve.sample(amount);
		return samples.map(s => {
			const {point, direction} = s;
			return (new Burgeon(point.x, point.y, direction, nextType(this.type), this) );
		});
	}
}