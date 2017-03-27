import Bezier from 'bezier-js';
import MagneticCurve from '../model/MagneticCurve';
import CurveManagement from './CurveManagement';
import * as UI from '../model/UIManagement';
import * as Collision from '../model/Collision';
import LeafImage from '../images/LeafImage';
import nextType from './MarkovLeaf';

let shortid = require('shortid');
let flowerString = require('../images/海石榴心.svg');

export default class LevelCurve {
	/**
		@param levelParam : array
		length : number
		alpha : number
		branches : number
	*/
	constructor(basePath, trunkWidth, levelParam, _scene){
		this.id = shortid.generate();
		this.basePath = basePath;
		this.trunkWidth = trunkWidth;
		this.levelParam = levelParam;

		this.leafLayer = undefined;
		this.stemLayer = undefined;
		this.flowerLayer = undefined;
		this.debugCurveLayer = undefined;

		this.scene = _scene;
		this.mags = [];
		this.leafQueue = [];
	}
	drawLevelCurve(beziers, level, type, parent){
		if( !beziers ) return;

		let Bs = beziers.map(b =>{
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
		let tempMags = [];
		this.branchPosition(level).forEach((i) => {
			if( totalLength === 0) return;

			let bezierIndex = 0;

			let pos = totalLength * i;
			while( pos >= Bs[bezierIndex].length() ) {
				pos -= Bs[bezierIndex].length();
				bezierIndex++;
			}
			let bezierAtIndex = Bs[bezierIndex];
			let posOnSinglebezier = pos / bezierAtIndex.length();

			tempMags.push({posOnSinglebezier, bezierAtIndex, level, parent, type});
		});
		tempMags.reverse();
		this.mags = this.mags.concat(tempMags);
	}
	drawAt(t, b, sign, level, type, parent){
		let start = b.get(t);
		let v = b.derivative(t);

		// do collision
		let mag1 = new MagneticCurve({
			startX: start.x,
			startY: start.y,
			vx: v.x,
			vy: v.y,
			T: this.levelParam[level].length,
			alpha: this.levelParam[level].alpha,
			sign:1
		});
		let mag2 = new MagneticCurve({
			startX: start.x,
			startY: start.y,
			vx: v.x,
			vy: v.y,
			T: this.levelParam[level].length,
			alpha: this.levelParam[level].alpha,
			sign: -1
		});
		//draw mag1 ignore mag2
		if( Collision.testCollision(mag1.bbox(), this.scene, parent, mag2.bbox() ) ){
			//retry
			if (level < this.levelParam.length-1 ) this.drawAt(t, b, sign, level+1, type, parent);
			return;
		}
		else
		{	
			let curve = mag1.getCurve();
			const x1 = curve[0][0];
			const y1 = curve[0][1];
			const x2 = curve[curve.length-1][0];
			const y2 = curve[curve.length-1][1];
			sign = 1;
			this.leafQueue.push({x1, y1, x2, y2, sign, type});
			mag1.drawOn(this.debugCurveLayer, level);

			if (level < this.levelParam.length-1 ) this.drawLevelCurve(mag1.points, level, nextType(type), mag1.bbox() );
		}
		//draw mag2 ignore mag1
		if( Collision.testCollision(mag2.bbox(), this.scene, parent, mag1.bbox() ) ){
			if (level < this.levelParam.length-1 ) this.drawAt(t, b, -1*sign, level+1, type, parent);
			return;
		}
		else
		{		
			let curve = mag2.getCurve();
			const x1 = curve[0][0];
			const y1 = curve[0][1];
			const x2 = curve[curve.length-1][0];
			const y2 = curve[curve.length-1][1];
			sign = -1*sign;
			this.leafQueue.push({x1, y1, x2, y2, sign, type});

			mag2.drawOn(this.debugCurveLayer, level);

			if (level < this.levelParam.length-1 ) this.drawLevelCurve(mag2.points, level, nextType(type), mag2.bbox() );
		}
	}
	drawOn( pannel ){
		this.pannel = pannel;
		let { leafLayer, stemLayer, flowerLayer, debugCurveLayer } = CurveManagement.layer;
		this.leafLayer = leafLayer;
		this.stemLayer = stemLayer;
		this.flowerLayer = flowerLayer;
		this.debugCurveLayer = debugCurveLayer;

		CurveManagement[this.id] = this;

		this.drawLevelCurve(this.basePath, 0, 0);
		while(this.mags.length > 0){
			let { posOnSinglebezier, bezierAtIndex, level, parent, type } = this.mags[0];
			this.drawAt( posOnSinglebezier, bezierAtIndex,  1, level, type, parent );
			this.mags.shift();
		}

		this.leafQueue.reverse();
		this.leafQueue.forEach(	({x1, y1, x2, y2, sign, type}) => this.drawLeaf(x1, y1, x2, y2, sign, type) );
		this.leafQueue.length = 0;
		this.drawStem();
		this.drawFlower();
	}
	drawStem(){
		let stem = this.stemLayer.group();
		stem.addClass('clickable');
		stem.data({ id:this.id });
		stem.click(() =>{
			CurveManagement.selectedCurve.push(this);
		});
		this.drawOutLine( stem, UI.state.trunkHead,UI.state.trunkTail,'#7B5A62');
		this.drawOutLine( stem, UI.state.trunkHead-1,UI.state.trunkTail-1,'#F9F2F4');
		this.drawOutLine( stem, UI.state.trunkHead/1.111,UI.state.trunkTail/1.111, '#CED5D0');
		this.drawOutLine( stem, UI.state.trunkHead/2,UI.state.trunkTail/2, '#9FB9A8');
		this.drawOutLine( stem, UI.state.trunkHead/8,UI.state.trunkTail/8, '#7C8168');
	}
	drawOutLine(pannel, beginWidth, endWidth, color, _basePath = this.basePath){
		let basePath = _basePath.map(c =>{
			return new Bezier(
				c[0][0], c[0][1],
				c[1][0], c[1][1],
				c[2][0], c[2][1],
				c[3][0], c[3][1]
			);
		});
		

		let totalLength = basePath.reduce( (length, Bezier) =>{
			return Bezier.length() + length;
		}, 0);

		let l = 0;
		let w = (endWidth-beginWidth)/totalLength;
		let outline = basePath.map(b => {
			let d1 = beginWidth + l * w;
			l += b.length();
			let d2 = beginWidth + l * w;
			if(d1 === 0) d1 = 1;
			return b.outline(d1, d1, d2, d2);
		});
		outline.forEach(b => {
			let pathString = fittedCurveToPathString(b);
			drawOnPannel( pannel, pathString, color );
		});
	}
	drawFlower(){
		let blackCircle = {
			cx: this.basePath[this.basePath.length-1][3][0],
			cy: this.basePath[this.basePath.length-1][3][1],
			r:  UI.state.trunkTail*2 
		};
	
		let g = this.flowerLayer.group();
		g.addClass('clickable');
		g.data({ id:this.id });
		g.click(()=>{
			CurveManagement.selectedCurve.push(this);
		});
		let flower = g.svg(flowerString);
		const boundingCircle = flower.node.children[0].children[1].children[0].children[0];
		const cr = boundingCircle.getAttribute('r');

		const rate = blackCircle.r * 2/cr;
				
		flower.transform({
			scale: rate,
			cx: cr,
			cy: cr,
		}).transform({
			x: blackCircle.cx,
			y: blackCircle.cy,
		}).transform({
			rotation: 30,
			cx: cr,
			cy: cr,
		});
	}
	drawLeaf(x1, y1, x2, y2, sign, type){
		// let num = Math.floor(Math.random() * ((6-2)+1) + 1);
		let leafString = LeafImage[type];

		// let leafString = LeafImage[type];		
		let g = this.leafLayer.group();
		// g.hide();
		let leaf = g.svg(leafString);
		const direct = leaf.node.children[0].getElementById('direct');
		const direct_x1 = direct.getAttribute('x1');
		const direct_y1 = direct.getAttribute('y1');
		const direct_x2 = direct.getAttribute('x2');
		const direct_y2 = direct.getAttribute('y2');

		const skeletonLength = distance(x1, y1, x2, y2);
		const directLength = distance(direct_x1, direct_y1, direct_x2, direct_y2);

		const toDeg = 180/Math.PI;

		const redLineAngle = Math.atan2( direct_y2 - direct_y1, direct_x2-direct_x1 )* toDeg;
		const leafCurveAngle = Math.atan2( y2 - y1, x2 - x1)* toDeg;
		const roateAngle = (leafCurveAngle - redLineAngle );

		if(sign > 0) 
			leaf.flip('y').transform({
				scale: skeletonLength / directLength
			}).transform({
				x: x1,
				y: y1,
			}).transform({
				rotation: +roateAngle +180-(leafCurveAngle*2),
				cx: direct_x1,
				cy: direct_y1,
			});

		else
			leaf.transform({
				scale: skeletonLength / directLength
			}).transform({
				x: x1,
				y: y1,
			}).transform({
				rotation: roateAngle ,
				cx: direct_x1,
				cy: direct_y1,
			});
	}
	redraw() {
		if( this.pannel === undefined ) {
			console.error('can not redraw!');
			return;
		}
		this.scene.length = 0;
		this.mags.length = 0;
		this.leafQueue.length = 0;

		this.flowerLayer.clear();
		this.stemLayer.clear();
		this.leafLayer.clear();
		this.debugCurveLayer.clear();
		this.drawOn(this.pannel);
	}
	branchPosition(level){
		let branches = [];
		let branch = this.levelParam[level].branches;

		let pos = 1 / (branch + 1);
		for (let i = 1; i <= branch; i++) {
			branches.push(pos * i);
		}

		return branches;
	}
}
function fittedCurveToPathString(fittedLineData) {
	var str = '';
	//bezier : [ [c0], [c1], [c2], [c3] ]
	fittedLineData.curves.map(function (bezier, i) {
		if (i == 0) {
			str += 'M ' + bezier.points[0].x + ' ' + bezier.points[0].y;
		}

		str += 'C ' + bezier.points[1].x + ' ' + bezier.points[1].y + ', ' +
		bezier.points[2].x + ' ' + bezier.points[2].y + ', ' +
		bezier.points[3].x + ' ' + bezier.points[3].y + ' ';	
					
	});

	return str;
}
function drawOnPannel(pannel, pathString, color){
	pannel.path( pathString )
	.fill(color)
	.stroke({width: 0});
}
function distance(x1, y1, x2, y2) {
	const a = x1 - x2;
	const b = y1 - y2;

	return Math.sqrt( a*a + b*b );
}
