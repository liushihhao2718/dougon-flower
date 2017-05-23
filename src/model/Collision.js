import * as UI from '../model/UIManagement';
import {dougonBoundingNodes} from '../images/dougonBounding';
let inside = require('point-in-polygon');
import _ from 'lodash';

function testCollision(test_polygon, polygons, ...ignore){
	const dougon = dougonBoundingNodes[ UI.state.bound ];
	const all_inside = true;
	if(!polygonInsidePolygon(test_polygon, dougon, all_inside)) return true;
	return !!polygons.find(poly=> !ignore.includes(poly)
		&& polygonCollision(test_polygon, poly)
	);
}

function aabbCollision(rect1, rect2){

	return (
		rect1.x < rect2.x + rect2.width &&
		rect1.x + rect1.width > rect2.x &&
		rect1.y < rect2.y + rect2.height &&
		rect1.height + rect1.y > rect2.y
	);
}

function polygonCollision(polygon1, polygon2){
	return polygonInsidePolygon(polygon1, polygon2);
}
function bezierIntersects(polyBezier, others){
	let flag = false;
	for (var i = 0; i < others.length; i++) {
		for (var j = 0; j < polyBezier.length; j++) {
			if( others[i].intersects(polyBezier[j]).length > 0 ){
				flag = true;
				break;
			}
		}
	}
	return flag;
}
function polygonInsidePolygon(poly,fixedPoly, all_inside){
	return all_inside?
		poly.every(p => inside(p, fixedPoly))
		:!!poly.find(p => inside(p, fixedPoly));
}

function circlePolygonCollision(circle, polygon) {
	return !!polygon.find(p => distance(circle.x, circle.y, p[0], p[1]) < circle.r);
}

function distance(x1, y1, x2, y2) {
	const a = x1 - x2;
	const b = y1 - y2;

	return Math.sqrt( a*a + b*b );
}

export default class ColliderCollection{
	constructor(polygons=[]){
		this.colliders = polygons.map(x => (
			{collider:x, bbox:makeBbox(x)}
		));
	}
	test(polygon, ...ignore){
		let bbox = makeBbox(polygon);
		const dougon = dougonBoundingNodes[ UI.state.bound ];

		if(!polygonInsidePolygon(polygon, dougon, true)) return true;

		return !!this.colliders
		.find(
			x=>!ignore.includes(x.collider)
			&& aabbCollision(x.bbox,bbox)
			&& polygonCollision(x.collider,polygon)
		);
	}
	add(polygon){
		this.colliders.push({collider:polygon, bbox:makeBbox(polygon)});
	}
}

function makeBbox(points){
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