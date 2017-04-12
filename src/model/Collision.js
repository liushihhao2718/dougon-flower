import * as UI from '../model/UIManagement';
import {dougonBoundingNodes} from '../images/dougonBounding';
let inside = require('point-in-polygon');


export function testCollision(test_bbox, bboxes, ...ignore){
	let isCollision = false;

	// if(!insideBound(test_bbox)) return true;
	if(!rectInsidePolygon(test_bbox)) return true;

	for (var i = 0; i < bboxes.length; i++) {
		if (ignore.includes( bboxes[i] )) continue; 
		isCollision = aabbCollision(test_bbox, bboxes[i]);
		if(isCollision)break;
	}
	if( !isCollision  && !bboxes.includes(test_bbox) ) {
		bboxes.push(test_bbox);
	} 
	
	return isCollision;
}

export function aabbCollision(rect1, rect2){

	return (
		rect1.x < rect2.x + rect2.width &&
		rect1.x + rect1.width > rect2.x &&
		rect1.y < rect2.y + rect2.height &&
		rect1.height + rect1.y > rect2.y
   );
}

export function bezierIntersects(polyBezier, others){
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

function rectInsidePolygon(rect){
	const polygon = dougonBoundingNodes[ UI.state.bound ];
	let rectPoints = [
		[rect.x, rect.y],
		[rect.x+rect.width, rect.y],
		[rect.x, rect.y + rect.height],
		[rect.x + rect.width, rect.y + rect.height]
	];
	return rectPoints.every(p => inside(p, polygon));
}