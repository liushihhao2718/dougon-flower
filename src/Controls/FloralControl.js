import {BezierSpline} from '../model/Spline';
import {Floral} from '../model/stem';
import CurveManagement from '../model/CurveManagement';
import * as UI from '../model/UIManagement';

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
		let aspect = UI.state.aspect;
		let rotation = angle(this.startPoint[0], this.startPoint[1], this.endPoint[0], this.endPoint[1]);
		CurveManagement.floralScene.push(new Floral(new BezierSpline([this.startPoint],0.2), this.radius,0, UI.state.trunkTail,'海石榴華', aspect ,rotation));
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

function angle(x1,y1, x2, y2) {
	const toDeg = 180/Math.PI;
	return 	 Math.atan2( y2 - y1, x2 - x1)* toDeg +90;
}