import {BezierSpline} from '../model/Spline';
import {Floral} from '../model/stem';
import CurveManagement from '../model/CurveManagement';
import * as UI from '../model/UIManagement';

const error = 50;

function PaintControl(pannel) {
	let rawPointData = [];
	let paintingPolyLine = undefined;	

	this.start = function( point ) {
		rawPointData.push( point );
		paintingPolyLine = pannel.polyline().fill('none').stroke({ width: 1 });
	};
	this.update = function( point ) {
		rawPointData.push( point );
		updatePolyLineSVG( paintingPolyLine, rawPointData);
	};

	this.end = function() {
		let smoothBizer = BezierSpline.makeByPoints( rawPointData, error );
		if(smoothBizer.length == 0) {
			clearRawData();
			return;
		}
		let aspect = UI.state.aspect;
		
		const pLast = rawPointData[rawPointData.length-1];
		const pNear = rawPointData[rawPointData.length-6];
		const direction = {
			x : pLast[0] - pNear[0],
			y : pLast[1] - pNear[1]
		};
		const toDeg = 180/Math.PI;

		const redLineAngle = Math.atan2( direction.y, direction.x )* toDeg;
		const flowerAngle = -90;
		const rotateAngle = -1*(flowerAngle - redLineAngle );

		CurveManagement.floralScene.push( new Floral(smoothBizer,UI.state.flowerSize, UI.state.trunkHead, UI.state.trunkTail,'海石榴華', aspect, rotateAngle) );
		CurveManagement.draw();
		UI.changeColor(UI.state.color);
		clearRawData();
	};

	function updatePolyLineSVG(paintingPolyLine, rawPointData) {
		paintingPolyLine.plot( rawPointData );
	}
	
	function clearRawData(){
		rawPointData = [];
		paintingPolyLine.remove();
	}	
}

export default PaintControl;