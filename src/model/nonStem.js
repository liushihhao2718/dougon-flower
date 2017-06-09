import shortid from 'shortid';
import * as Drawer from './Drawer';
import {Leaf} from './stem';
import _ from 'lodash';

export class Flower{
	constructor(x,y, r,flowerType='海石榴華',flowerRotation=30){
		this.id = shortid.generate();
		this.flowerType = flowerType;
		this.flowerRotation = flowerRotation;
		this.flowerPosition = { x,y,r };
	}
	burgeons(){
		return undefined;
	}
	draw(){
		Drawer.drawFlower(this);
	}
}

export class LeafBranch{
	constructor(spline){
		this.spline = spline;
	}
	static segmentedIndex(){
		let segment = [3,3,3,4,4,6,15];
		segment.reverse();
		let total = _.sum(segment);
		let normalizedSegment = _.map(segment, x=>x/total);
		partialSum(normalizedSegment);
		normalizedSegment.unshift(0);
		return normalizedSegment;
	}
	get colliders(){
		if(!this._colliders){
			this._colliders = this.computeSegmentLeaf().map(x=>x.colliders);
		}
		return this._colliders;
	}
	computeSegmentLeaf(){
		if(!this.leafs){
			this.leafs = [];

			let segment = LeafBranch.segmentedIndex();
			for (var i = 0; i < segment.length-1; i++) {
				let start = segment[i];
				let end= segment[i+1];

				if(i !== 0)	start -= 0.01;

				let beziers = this.spline.segmentRange(start, end);
				let points = _.flatMap(beziers, b=> b.getLUT(100)).map(p=>[p.x, p.y]);
				let type = 6 + i;
				const sign = -1;
				this.leafs.push(Leaf.makeLeafByPoint(points, sign, type));
			}
			this.leafs.reverse();
		}
		return this.leafs;
	}

	draw(){
		this.computeSegmentLeaf().forEach(leaf=>{
			Drawer.drawLeaf(leaf);
		});
		Drawer.drawPolygon(this.colliders );
	}
	burgeons(){
		return undefined;
	}
}

function partialSum(numArray) {
	let partialSum = 0;
	for (var i = 0; i < numArray.length; i++) {
		partialSum += numArray[i];
		numArray[i] = partialSum;
	}
}