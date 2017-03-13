let bboxes = [];

export function testCollision(test_bbox, ignore){
	let isCollision = false;

	if (bboxes.length === 0){
		bboxes.push(test_bbox);
		return false;
	}

	for (var i = 0; i < bboxes.length; i++) {
		if (ignore === bboxes[i]) return; 
		isCollision = aabbCollision(test_bbox, bboxes[i] );
	}
	if( !isCollision  && !bboxes.includes(test_bbox) ) {
		bboxes.push(test_bbox);
	}
	return isCollision;
}
export function aabbCollision(rect1, rect2){
	const padding = 0;
	return (
		rect1.x <= rect2.x + rect2.width + padding&&
		rect1.x + rect1.width+padding >= rect2.x &&
		rect1.y <= rect2.y + rect2.height+padding &&
		rect1.height + rect1.y+padding >= rect2.y
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