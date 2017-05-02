import shortid from 'shortid';
import * as Drawer from './Drawer';
import {Leaf} from './stem';
import _ from 'lodash';
import {leafType} from '../images/LeafImage';

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
	constructor(spline, colliderWidth = 90){
		this.spline = spline;
		let beziers =  spline.outline(colliderWidth, 10);
		this.colliders = _.flatMap(beziers,b => b.getLUT(10))
			.map(p=>[p.x, p.y]);
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

	computeSegmentLeaf(){
		if(!this._leaf) this.leafs = [];

		let segment = LeafBranch.segmentedIndex();
		for (var i = 0; i < segment.length-1; i++) {
			let start = segment[i];
			let end= segment[i+1];

			if(i !== 0)	start -= 0.01;

			let beziers = this.spline.segmentRange(start, end);
			let points = _.flatMap(beziers, b=> b.getLUT(100)).map(p=>[p.x, p.y]);
			let type = {
				name :leafType.leafBranch,
				order: i
			};
			this.leafs.push(Leaf.makeLeafByPoint(points, 0, type));
		}
		this.leafs.reverse();
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