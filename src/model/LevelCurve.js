import Bezier from 'bezier-js';
import MagneticCurve from '../model/MagneticCurve';

export default class LevelCurve {
	/**
		@param levelParam : array
		length : number
		alpha : number
		branches : number
	*/
	constructor(basePath, trunkWidth, levelParam){
		this.basePath = basePath;
		this.trunkWidth = trunkWidth;
		this.levelParam = levelParam;
	}
	drawLevelCurve(beziers, level){
		beziers.forEach(b =>{
			const temp_b = new Bezier(
				b[0][0], b[0][1],
				b[1][0], b[1][1],
				b[2][0], b[2][1],
				b[3][0], b[3][1]
			);
			this.branchPosition(level).forEach((i, index) => {
				let sign = index%2 == 0 ? 1 : -1;
				this.drawAt(i, temp_b, sign, level);
			});
		});
	}
	drawAt(t, b, sign, level){
		let start = b.get(t);
		let v = b.derivative(t);
		let mag = new MagneticCurve({
			startX: start.x,
			startY: start.y,
			vx: v.x,
			vy: v.y,
			T: this.levelParam[level].length,
			alpha: this.levelParam[level].alpha,
			sign: sign
		});
		mag.drawOn(this.pannel);
		if (level < this.levelParam.length-1 ) this.drawLevelCurve(mag.points, level+1);
	}

	drawOn( pannel ){
		this.pannel = pannel;
		this.drawLevelCurve(this.basePath, 0);
	}

	branchPosition(level){
		let branches = [];
		let branch = this.levelParam[level].branches;

		let pos = 1 / (branch + 1);
		for (let i = 1; i <= branch; i++) {
			branches.push(pos * i);
		}

		return branches;
	}
}

