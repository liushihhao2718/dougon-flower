import CurveManagement from './CurveManagement';
import {flowerString, LeafImage, cap, leafType, LeafBranch} from '../images/LeafImage';
import {Vector2} from 'three';
function handleClick(floral) {
	if(CurveManagement.selectedCurve.includes(floral)){
		const index = CurveManagement.selectedCurve.indexOf(floral);
		CurveManagement.selectedCurve.splice(index, 1);
	}
	else{
		CurveManagement.selectedCurve.push(floral);
		console.log(`${floral.id} clicked`);
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
	
	capSVG.transform({
		scale: rate
	}).transform({
		x: x1,
		y: y1,
	}).transform({
		rotation: roateAngle ,
		cx: direct_x1,
		cy: direct_y1,
	});


	let matrix = capSVG.node.getAttribute('transform').split(/[^\-\d.]+/).filter(x=>x !== '');

	let pos = multiplyMatrixAndPoint(matrix, [cx,cy,1]);
	floral.flowerPosition.x = pos[0];
	floral.flowerPosition.y = pos[1];
	floral.flowerPosition.r = rate * r;


	setFloralCollider(floral, capSVG);
}

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

	drawOutLine( stem, trunkHead, trunkTail,'#7B5A62', floral.curve);
	drawOutLine( stem, trunkHead/1.111, trunkTail/1.111, '#CED5D0',floral.curve);
	drawOutLine( stem, trunkHead/2, trunkTail/2, '#9FB9A8',floral.curve);
	drawOutLine( stem, trunkHead/8, trunkTail/8, '#7C8168',floral.curve);
}
function drawOutLine( layer, beginWidth, endWidth, color, basePath){
	let fullPath = basePath.outline(beginWidth, endWidth);

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

	const rate = blackCircle.r/cr;
			
	flower.transform({
		scale: rate
	});
	flower.transform({
		x: blackCircle.cx,
		y: blackCircle.cy,
	});
	flower.transform({
		rotation: floral.flowerRotation,
		cx: cr,
		cy: cr,
	});

	setFloralCollider(floral, flower);
	
}
function setFloralCollider(floral, svg) {
	let matrix = svg.node.getAttribute('transform').split(/[^\-\d.]+/).filter(x=>x !== '');

	const collider = svg.select('#collider').members[0].node;

	let colliderPointsInSVG = collider.attributes.points.value.split(' ').map(s=>s.split(',').map(n=>Number(n)) );
	let colliderPointsInWorld = colliderPointsInSVG.map( ([px, py]) => multiplyMatrixAndPoint(matrix, [px, py, 1]).slice(0,2));
	
	floral.colliders = colliderPointsInWorld;
}
export function drawMagneticCurve(leaf, level){
	let level_color = ['orange', 'green', 'blue', 'black'];
	CurveManagement.layer.debugCurveLayer.path(leaf.curve.svgString() ).fill('none').stroke({ width: 3 }).stroke(level_color[level]);
}
export function	drawLeaf(leaf){
	let x1 = leaf.startX;
	let y1 = leaf.startY;
	let x2 = leaf.endX;
	let y2 = leaf.endY;
	let sign = leaf.sign;
	let type = leaf.type;

	let leafString = '';
	switch(leaf.type.name){
	case leafType.leaf:
		leafString = LeafImage[type.order];
		break;
	case leafType.leafBranch:
		leafString = LeafBranch[type.order];
		break;	
	}
	let g = CurveManagement.layer.leafLayer.group();

	let leafSVG = g.svg(leafString);
	const direct = leafSVG.node.children[0].getElementById('direct');
	const direct_x1 = Number(direct.getAttribute('x1'));
	const direct_y1 = Number(direct.getAttribute('y1'));
	const direct_x2 = Number(direct.getAttribute('x2'));
	const direct_y2 = Number(direct.getAttribute('y2'));

	const skeletonLength = distance(x1, y1, x2, y2);
	const directLength = distance(direct_x1, direct_y1, direct_x2, direct_y2);

	const toDeg = 180/Math.PI;

	const redLineAngle = Math.atan2( direct_y2 - direct_y1, direct_x2-direct_x1 )* toDeg;
	const leafCurveAngle = Math.atan2( y2 - y1, x2 - x1)* toDeg;
	const roateAngle = (leafCurveAngle - redLineAngle );
	const rate = skeletonLength / directLength;
	sign = -1;
	if(sign > 0) 
		leafSVG.flip('y').transform({
			scale: rate
		}).transform({
			x: x1,
			y: y1,
		}).transform({
			rotation: roateAngle +180-(leafCurveAngle*2),
			cx: direct_x1,
			cy: direct_y1,
		});

	else
		leafSVG.transform({
			scale: rate
		}).transform({
			x: x1,
			y: y1,
		}).transform({
			rotation: roateAngle ,
			cx: direct_x1,
			cy: direct_y1,
		});

		
	// setFloralCollider(leaf, leafSVG);

	// let matrix = leafSVG.node.getAttribute('transform').split(/[^\-\d.]+/).filter(x=>x !== '');

	// let [px, py] = collider.attributes.points.value.split(' ')[0].split(',').map(n=>Number(n));
	// let [px,py] = multiplyMatrixAndPoint(matrix, [direct_x2, direct_y2, 1]).slice(0,2);
	// let p2 = new Vector2(direct_x2, direct_y2);
	// let origin = new Vector2(direct_x1, direct_y1);
	// let fixPoint = new Vector2(x1, y1);
	// p2.sub(origin).multiplyScalar(rate).add(origin).add(fixPoint).rotateAround(fixPoint, roateAngle);
	
	// console.log(close(p2.x,x2) && close(p2.y, y2));
}
function close(a, b) {
	return Math.abs(a - b) < 0.001;
}
export 	function drawBasePath(pathString){
	CurveManagement.layer.debugCurveLayer.path( pathString ).fill('none').stroke({ width: 3 }).stroke('#f06');
}
export function drawPolygon(polygon){
	CurveManagement.layer.debugCurveLayer.polygon().plot(polygon).fill('none').stroke({ width: 1,color:'red' });
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