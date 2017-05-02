import MagneticCurve from './MagneticCurve';
import { BezierSpline } from './Spline';
import shortid from 'shortid';
import nextType from './MarkovLeaf';
import * as Drawer from './Drawer';
import _ from 'lodash';
import {leafType} from '../images/LeafImage';

export class Floral{
	constructor(basePath, r,trunkHead,trunkTail,flowerType='海石榴華', aspect = '正面',flowerRotation=30){
		this.id = shortid.generate();
		this.curve = (aspect === '正面') ?basePath:basePath.range(0, 0.7);
		this.flowerType = flowerType;
		this.flowerRotation = flowerRotation;
		this.colliders = undefined;
		let point = _.last(this.curve.controlPoints)[3];
		this.flowerPosition = {
			x: point[0],
			y: point[1],
			r
		};
		this.trunkHead = trunkHead;
		this.trunkTail = trunkTail;
		this.aspect = aspect;
	}
	burgeons(amount){
		if (this.curve.length === 0) return;
		let samples = this.curve.sample(amount);

		let type = {
			name: leafType.leaf,
			order: 0
		};

		return samples.map(s => {
			const {point, direction} = s;
			return (new Burgeon(point.x, point.y, direction, type, this) );
		});
	}
	draw(){
		Drawer.drawStem(this);
		Drawer.drawFlower(this);
		Drawer.drawBasePath(this.curve.svgString());
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
		let mag = new MagneticCurve({
			startX:this.x, 
			startY:this.y, 
			vx:this.direction.x,
			vy: this.direction.y,
			T: length,
			alpha,
			sign
		});
		

		return Leaf.makeLeafByPoint(mag.getCurve(), sign, this.type);
	}
}

export class Leaf {
	static makeLeafByPoint(points, sign, type){
		const lastOne = points.length -1;
		let startX = points[0][0];
		let startY = points[0][1];
		let endX = points[lastOne][0];
		let endY = points[lastOne][1];

		return new Leaf(startX,startY,
			endX,endY,
			BezierSpline.makeByPoints(points),
			sign,
			type
		);
	}

	constructor(startX,startY,endX,endY,curve, sign,type) {
		this.startX = startX;
		this.startY = startY;
		this.endX = endX;
		this.endY = endY;
		this.type = type;
		this.curve = curve;
		this.colliders = this.curve.colliders;
		this.sign = sign;
	}
	burgeons(amount){
		let samples = this.curve.sample(amount);
		return samples.map(s => {
			const {point, direction} = s;
			return (new Burgeon(point.x, point.y, direction, nextType(this.type), this) );
		});
	}
}