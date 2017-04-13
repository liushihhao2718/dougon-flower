import CurveManagement from './CurveManagement';
import * as UI from './UIManagement';
import flowerString from '../images/海石榴心.svg';
import LeafImage from '../images/LeafImage';

export function	draw(){
	CurveManagement[this.id] = this;

	this.makeGrowthPoints(this.basePath, 0, 0);
	let sign = 1;
	while(this.mags.length > 0){
		let { posOnSinglebezier, bezierAtIndex, level, parent, type } = this.mags[0];
		this.drawAt( posOnSinglebezier, bezierAtIndex,  sign, level, type, parent );
		sign *= -1;

		this.mags.shift();
	}
	
	this.leafQueue.reverse();
	this.leafQueue.forEach(	({x1, y1, x2, y2, sign, type}) => this.drawLeaf(x1, y1, x2, y2, sign, type) );
	this.leafQueue.length = 0;
	this.drawStem();
	this.drawFlower();
}
export function	drawStem(floral){
	let stem = CurveManagement.layer.stemLayer.group();
	stem.addClass('clickable');
	stem.data({ id:floral.id });
	stem.click(() =>{
		CurveManagement.selectedCurve.push(floral);
	});
	drawOutLine( stem, UI.state.trunkHead,UI.state.trunkTail,'#7B5A62', floral.curve);
	drawOutLine( stem, UI.state.trunkHead-1,UI.state.trunkTail-1,'#F9F2F4',floral.curve);
	drawOutLine( stem, UI.state.trunkHead/1.111,UI.state.trunkTail/1.111, '#CED5D0',floral.curve);
	drawOutLine( stem, UI.state.trunkHead/2,UI.state.trunkTail/2, '#9FB9A8',floral.curve);
	drawOutLine( stem, UI.state.trunkHead/8,UI.state.trunkTail/8, '#7C8168',floral.curve);
}
function drawOutLine( layer, beginWidth, endWidth, color, basePath){
	let totalLength = basePath.length;

	let l = 0;
	let w = (endWidth-beginWidth)/totalLength;
	let outline = basePath.beziers.map(b => {
		let d1 = beginWidth + l * w;
		l += b.length();
		let d2 = beginWidth + l * w;
		if(d1 === 0) d1 = 1;
		return b.outline(d1, d1, d2, d2);
	});

	let f = [];
	let r = [];
	let head=outline[0].curves[0];
	let tail;
	const lastPath = outline[outline.length-1];
	tail = lastPath.curves[lastPath.curves.length/2];

	outline.forEach(b => {
		const length = Math.floor(b.curves.length/2)- 1;
		f = f.concat(b.curves.slice(1, length+1));
		r = r.concat(b.curves.slice(length+2).reverse());
	});

	let fullPath = [].concat([head],f,[tail], r.reverse());
	

	let pathString = svgPathString(fullPath);
	drawOnPannel( layer, pathString, color );
}
function drawOnPannel(pannel, pathString, color){
	pannel.path( pathString )
	.fill(color)
	.stroke({width: 0});
}
function svgPathString(fittedLineData) {
	var str = '';
	//bezier : [ [c0], [c1], [c2], [c3] ]
	fittedLineData.forEach(function (bezier, i) {
		if (i == 0) {
			str += 'M ' + bezier.points[0].x + ' ' + bezier.points[0].y;
		}

		str += 'C ' + bezier.points[1].x + ' ' + bezier.points[1].y + ', ' +
		bezier.points[2].x + ' ' + bezier.points[2].y + ', ' +
		bezier.points[3].x + ' ' + bezier.points[3].y + ' ';	
					
	});

	return str;
}
export function	drawFlower(floral){
	let blackCircle = {
		cx: floral.flowerPosition.x,
		cy: floral.flowerPosition.y,
		r:  UI.state.flowerSize 
	};

	let g = CurveManagement.layer.flowerLayer.group();
	g.addClass('clickable');
	g.data({ id:floral.id });
	g.click(()=>{
		CurveManagement.selectedCurve.push(floral);
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
export function	drawLeaf(leaf){
	let x1 = leaf.startX;
	let y1 = leaf.startY;
	let x2 = leaf.endX;
	let y2 = leaf.endY;
	let sign = leaf.sign;
	let type = leaf.type;

	let leafString = LeafImage[type];
	let g = CurveManagement.layer.leafLayer.group();
	// g.hide();
	let leafSVG = g.svg(leafString);
	const direct = leafSVG.node.children[0].getElementById('direct');
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
		leafSVG.flip('y').transform({
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
		leafSVG.transform({
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

function distance(x1, y1, x2, y2) {
	const a = x1 - x2;
	const b = y1 - y2;

	return Math.sqrt( a*a + b*b );
}