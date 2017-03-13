let bboxes = [];

export function testCollision(test_bbox, ignore){
	let collisionBox = [];
	let isCollision = false;

	if (bboxes.length === 0){
		bboxes.push(test_bbox);
		return false;
	}

	bboxes.forEach(bbox => {
		if (ignore === bbox) return; 
		isCollision = aabbCollision(test_bbox, bbox);
	});

	if( !bboxes.includes(test_bbox) ) {
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

function bezierIntersects(b1, b2){
	
}