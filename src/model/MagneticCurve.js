export default class MagneticCurve {
	static makeCurve(param) {
		/*
		startX, startY : 初始點
		vx, vy : 初始速度
		T : 總點數
		alpha : 等角螺線參數
		*/
		const points = [];

		let x = param.startX;
		let y = param.startY;

		let vx = param.vx;
		let vy = param.vy;

		let T = param.T;
		let t = 0;
		let dt = 1;

		while( t < T ) {
			points.push( [x, y] );
			const q = Math.pow( (T-t), this.param.alpha );
			x += vx * dt;
			y += vy * dt;

			const ax = -1 * vy * q;
			const ay = vx * q;

			vx += ax * dt;
			vy += ay * dt;
		}

		return points;
	}
}