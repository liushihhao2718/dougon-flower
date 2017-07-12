export function scale(point, sx, sy){
	let [x,y]=point;
	return [x*sx,y*sy];
}
export function rotate(point, r) {
	const [x, y] = point;
	const c = Math.cos(r);
	const s = Math.sin(r);
	return [c*x - s*y, s*x + c*y];
}
export function translate(point, dx, dy){
	const [x, y] = point;
	return [x+dx, y+dy];
}
export function distance(x1, y1, x2, y2) {
	const a = x1 - x2;
	const b = y1 - y2;

	return Math.sqrt( a*a + b*b );
}
export function multiplyMatrixAndPoint(matrix, point) {
  
	//Give a simple variable name to each part of the matrix, a column and row number
	var c0r0 = matrix.a, c1r0 = matrix.c, c2r0 = matrix.e;
	var c0r1 = matrix.b, c1r1 = matrix.d, c2r1 = matrix.f;
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