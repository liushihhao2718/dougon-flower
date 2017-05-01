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
	segmentedIndex(){
		let segment = [2,2,2,3,4,6,10];
		segment.reverse();
		let total = _.sum(segment);
		let normalizedSegment = _.map(segment, x=>x/total);
		partialSum(normalizedSegment);
		normalizedSegment.unshift(0);
		return normalizedSegment;
	}

	computeSegmentLeaf(){
		let leafs = [];
		let segment = this.segmentedIndex();
		for (var i = 0; i < segment.length-1; i++) {
			let start = segment[i], end= segment[i+1];
			let beziers = this.spline.segmentRange(start, end);
			let points = _.flatMap(beziers, b=> b.getLUT(100)).map(p=>[p.x, p.y]);
			leafs.push(Leaf.makeLeafByPoint(points, 1, 1));
		}
		return leafs;
	}

	draw(){
		this.computeSegmentLeaf().forEach( leaf=>Drawer.drawLeaf(leaf) );
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

	//partialSum = 0;
	//return numArray.map(x=>partialSum+=x);
}