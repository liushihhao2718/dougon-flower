import {BezierSpline} from '../model/Spline';
import {Floral} from '../model/stem';
import CurveManagement from '../model/CurveManagement';
import * as UI from '../model/UIManagement';
const error = 100;

function SkeletonControl(pannel) {
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
		let smoothBizer = new BezierSpline( rawPointData, error );
		if(smoothBizer.length == 0) {
			clearRawData();
			return;
		}
		
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

export default SkeletonControl;

