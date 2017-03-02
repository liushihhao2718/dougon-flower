import Bezier from 'bezier-js';
import MagneticCurve from '../model/MagneticCurve';
import CurveManagement from './CurveManagement';

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
		this.curveGroup = undefined;
	}
	drawLevelCurve(beziers, level){
		if( !beziers ) return;
		let sign = 1;

		let Bs = beziers.map(b =>{
			return new Bezier(
				b[0][0], b[0][1],
				b[1][0], b[1][1],
				b[2][0], b[2][1],
				b[3][0], b[3][1]
			);
		});

		let totalLength = Bs.reduce( (length, B) =>{
			return B.length() + length;
		}, 0);

		this.branchPosition(level).forEach((i) => {
			if( totalLength === 0) return;

			let bezierIndex = 0;

			let pos = totalLength * i;
			while( pos >= Bs[bezierIndex].length() ) {
				pos -= Bs[bezierIndex].length();
				bezierIndex++;
			}

			let posOnSinglebezier = pos / Bs[bezierIndex].length();

			this.drawAt(posOnSinglebezier, Bs[bezierIndex], sign, level);

			sign *= -1;
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
		mag.drawOn(this.curveGroup);
		if (level < this.levelParam.length-1 ) this.drawLevelCurve(mag.points, level+1);
	}

	drawOn( pannel ){
		this.pannel = pannel;
		this.curveGroup = pannel.group();
		pannel.add(this.curveGroup);
		CurveManagement[this.curveGroup.node.id] = this;
		this.curveGroup.on('click', e => {
			console.log('clecked');

			let curve_id = e.target.parentElement.id;
			let lvCurve = CurveManagement[ curve_id ];
			CurveManagement.selectedCurve.length = 0;
			CurveManagement.selectedCurve.push(lvCurve);

		});

		this.drawLevelCurve(this.basePath, 0);
	}
	redraw() {
		if( this.pannel === undefined ) {
			console.error('can not redraw!');
			return;
		}
		this.curveGroup.remove();
		this.drawOn(this.pannel);
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

