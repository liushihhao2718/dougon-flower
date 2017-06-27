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
		updateLines( paintingPolyLine, rawPointData);
	};

	this.end = function() {
		let smoothBizer = BezierSpline.makeByPoints( rawPointData, error );
		if(smoothBizer.length == 0) {
			clearRawData();
			return;
		}
		let aspect = UI.state.aspect;
		CurveManagement.floralScene.push( new Floral(smoothBizer,UI.state.flowerSize, UI.state.trunkHead, UI.state.trunkTail,'海石榴華', aspect, UI.state.rotation) );
		
		CurveManagement.draw();
		clearRawData();
	};

	function updateLines(paintingPolyLine, rawPointData) {
		paintingPolyLine.plot( rawPointData );
	}
	
	function clearRawData(){
		rawPointData = [];
		paintingPolyLine.remove();
	}	
}

export default PaintControl;