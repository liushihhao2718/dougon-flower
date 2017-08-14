import {Flower} from '../model/nonStem';
import CurveManagement from '../model/CurveManagement';
import {angle} from '../model/UtilMath';

class FloralControl{
	constructor(pannel){
		this.startPoint = [0,0];
		this.endPoint = [0,0];
		this.pannel = pannel;
		this.hintLine = pannel.line(0,0,0,0).fill('none').stroke('none').width(5);
		this.hintCircle = pannel.circle(0).fill('none').stroke('none').width(5);
	}

	start(point) {
		this.startPoint = point;
		this.update(point);
		this.showHint();
	}
	update(point){
		this.endPoint = point;
		this.updateHint();
	}
	end(){
		this.hideHint();
		let rotation = angle(this.startPoint[0], this.startPoint[1], this.endPoint[0], this.endPoint[1]);
		CurveManagement.floralScene.push(new Flower(this.startPoint[0],this.startPoint[1], this.radius,'海石榴華',rotation));
		CurveManagement.draw();
	}
	showHint(){
		this.hintLine.stroke('orange').width(5);
		this.hintCircle.stroke('orange').width(5);
	}
	hideHint(){
		this.hintLine.stroke('none');
		this.hintCircle.stroke('none');
	}
	updateHint(){
		this.hintLine.plot(this.startPoint[0],this.startPoint[1],this.endPoint[0],this.endPoint[1]);
		this.hintCircle.radius(this.radius).cx(this.startPoint[0]).cy(this.startPoint[1]);
	}

	get radius(){
		let dx = this.startPoint[0]-this.endPoint[0];
		let dy = this.startPoint[1]-this.endPoint[1];
		return Math.sqrt(dx*dx+dy*dy);
	}
}

export default FloralControl;