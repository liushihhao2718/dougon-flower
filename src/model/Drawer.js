import CurveManagement from './CurveManagement';
import {flowerString, LeafImage, cap} from '../images/LeafImage';
import SVG from 'svg.js';
const colorMap = require('../color/colorHex.json');

let symbols;

export function initSvgSymbol(){
	const panel = CurveManagement.panel;
	let flowerSymbol = panel.symbol().svg(flowerString);

	let leafSymbol = LeafImage.map(x => {
		let s = panel.symbol();
		return s.group().svg(x);
	});
	let capSymbol = panel.symbol().svg(cap);

	symbols = {flowerSymbol, leafSymbol, capSymbol};
}

function handleClick(floral) {
	if(CurveManagement.selectedCurve.includes(floral)){
		const index = CurveManagement.selectedCurve.indexOf(floral);
		CurveManagement.selectedCurve.splice(index, 1);
	}
	else{
		CurveManagement.selectedCurve.push(floral);
		// console.log(`${floral.id} clicked`);
	}
	CurveManagement.drawHint();
}
export function drawCap(floral) {
	let width = floral.trunkTail / 2;

	let lastBezier = floral.curve.beziers[floral.curve.beziers.length-1];
	let normal = lastBezier.normal(1);//{x, y}
	let lastPoint = lastBezier.get(1);//{x, y}
	let x1 = lastPoint.x - normal.x * width;
	let y1 = lastPoint.y - normal.y * width;
	let x2 = lastPoint.x + normal.x * width;
	let y2 = lastPoint.y + normal.y * width;
	let capString = cap;
	let g = CurveManagement.layer.stemLayer.group();
	g.click(()=> handleClick(floral));
	let capSVG = g.svg(capString);
	let direct = capSVG.node.children[0].getElementById('direct');
	const direct_x1 = direct.getAttribute('x1');
	const direct_y1 = direct.getAttribute('y1');
	const direct_x2 = direct.getAttribute('x2');
	const direct_y2 = direct.getAttribute('y2');

	const rate = 2 * width / distance(direct_x1, direct_y1, direct_x2, direct_y2);
	const boundingCircle = capSVG.node.children[0].getElementById('boundingCircle');
	const cx = Number(boundingCircle.getAttribute('cx'));
	const cy = Number(boundingCircle.getAttribute('cy'));
	const r = Number(boundingCircle.getAttribute('r'));

	const toDeg = 180/Math.PI;

	const redLineAngle = Math.atan2( direct_y2 - direct_y1, direct_x2-direct_x1 )* toDeg;
	const leafCurveAngle = Math.atan2( y2 - y1, x2 - x1)* toDeg;
	const roateAngle = (leafCurveAngle - redLineAngle );

	let matrix = new SVG.Matrix()
		.translate(x1, y1)
		.rotate(roateAngle)
		.scale(rate)
		.translate(-direct_x1, -direct_y1);
	capSVG.attr( 'transform',matrix.toString() );

	let pos = multiplyMatrixAndPoint(matrix.toString().split(/[^\-\d.]+/).filter(x=>x !== ''), [cx,cy,1]);
	floral.flowerPosition.x = pos[0];
	floral.flowerPosition.y = pos[1];
	floral.flowerPosition.r = rate * r;


	setFloralCollider(floral, capSVG);
}
let worker = new Worker('js/OutlineWorker.js');

export function	drawStem(floral){
	let stem = CurveManagement.layer.stemLayer.group();
	stem.addClass('clickable');
	stem.data({ id:floral.id });
	stem.click(() =>{
		CurveManagement.selectedCurve.push(floral);
	});

	let trunkHead = floral.trunkHead;
	let trunkTail = floral.trunkTail;
	if(floral.aspect == '側面') trunkTail*=0.5;

	let message ={
		basePath : floral.curve.controlPoints,
		trunk: [
			[trunkHead, trunkTail],
			[trunkHead/1.111, trunkTail/1.111],
			[trunkHead/2, trunkTail/2],
			[trunkHead/8, trunkTail/8]
		]
	};

	const color = ['赭筆描道', '綠華','二綠','大綠'];
	worker.postMessage(message);
	worker.onmessage = function(e){
		e.data.forEach((s,i)=>{
			drawOnPannel( stem, s, colorMap[color[i]]);
		});
	};
}

function drawOnPannel(pannel, pathString, color){
	pannel.path( pathString )
		.fill(color)
		.stroke({width: 0});
}

export function	drawFlower(floral){

	if (floral.aspect === '側面') {
		drawCap(floral);
		return;
	}
	let blackCircle = {
		cx: floral.flowerPosition.x,
		cy: floral.flowerPosition.y,
		r:  floral.flowerPosition.r
	};

	let g = CurveManagement.layer.flowerLayer.group();
	g.addClass('clickable');
	g.data({ id:floral.id });
	let flower = g.svg(flowerString);
	const boundingCircle = flower.node.children[0].children[1].children[0].children[0];

	g.click(()=> handleClick(floral));
	const cr = boundingCircle.getAttribute('r');
	const cx = boundingCircle.getAttribute('r');
	const cy = boundingCircle.getAttribute('r');

	const rate = blackCircle.r / cr;

	let matrix = new SVG.Matrix()
		.translate(blackCircle.cx, blackCircle.cy)
		.rotate(floral.flowerRotation)
		.scale(rate)
		.translate(-cx, -cy);
	flower.attr( 'transform',matrix.toString() );

	setFloralCollider(floral, flower);
	
}
function setFloralCollider(floral, svg) {
	let matrix = svg.node.getAttribute('transform').split(/[^\-\d.]+/).filter(x=>x !== '');

	const collider = svg.select('#collider').members[0].node;

	let colliderPointsInSVG = Array.from(collider.points);
	let colliderPointsInWorld = colliderPointsInSVG.map( p => 
		multiplyMatrixAndPoint(matrix, [p.x, p.y, 1]).slice(0,2)
	);
	drawPolygon(colliderPointsInWorld);
	floral.colliders = [colliderPointsInWorld];
}
export function drawMagneticCurve(leaf, level){
	let level_color = ['orange', 'green', 'blue', 'black'];
	CurveManagement.layer.debugCurveLayer
		.path(leaf.curve.svgString() )
		.fill('none').stroke({ width: 3 })
		.stroke(level_color[level]);
}
export function	drawLeaf(leaf){
	let g = CurveManagement.layer.leafLayer.group();
	g.use(symbols.leafSymbol[leaf.type]);
	g.attr({'transform': leaf.transformString()});
}

export 	function drawBasePath(pathString){
	CurveManagement.layer.debugCurveLayer
		.path( pathString )
		.fill('none')
		.stroke({ width: 3 })
		.stroke('#f06');
}
export function drawPolygon(polygon){
	CurveManagement.layer.debugCurveLayer
		.polygon()
		.plot(polygon)
		.fill('none')
		.stroke({ width: 3,color:'red' });
}

function distance(x1, y1, x2, y2) {
	const a = x1 - x2;
	const b = y1 - y2;

	return Math.sqrt( a*a + b*b );
}
function multiplyMatrixAndPoint(matrix, point) {
  
	//Give a simple variable name to each part of the matrix, a column and row number
	var c0r0 = matrix[ 0], c1r0 = matrix[ 2], c2r0 = matrix[ 4];
	var c0r1 = matrix[ 1], c1r1 = matrix[ 3], c2r1 = matrix[ 5];
	var c0r2 = 0, c1r2 = 0, c2r2 = 1;

	//Now set some simple names for the point
	var x = point[0];
	var y = point[1];
	var z = point[2];

	//Multiply the point against each part of the 1st column, then add together
	var resultX = (x * c0r0) + (y * c1r0) + (z * c2r0);

	//Multiply the point against each part of the 2nd column, then add together
	var resultY = (x * c0r1) + (y * c1r1) + (z * c2r1);

	//Multiply the point against each part of the 3rd column, then add together
	var resultZ = (x * c0r2) + (y * c1r2) + (z * c2r2);

	return [resultX, resultY, resultZ];
}