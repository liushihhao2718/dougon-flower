import * as UI from '../model/UIManagement';
import {dougonBoundingNodes} from '../images/dougonBounding';
let inside = require('point-in-polygon');

function aabbCollision(rect1, rect2){

	return (
		rect1.x < rect2.x + rect2.width &&
		rect1.x + rect1.width > rect2.x &&
		rect1.y < rect2.y + rect2.height &&
		rect1.height + rect1.y > rect2.y
	);
}

function polygonCollision(polygon1, polygon2){
	return  !!polygon1.find(p => inside(p, polygon2)) ||
			!!polygon2.find(p => inside(p, polygon1));
}

function polygonInsidePolygon(poly,fixedPoly){
	return poly.every(p => inside(p, fixedPoly));
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

		if(!polygonInsidePolygon(polygon, dougon)) return true;

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