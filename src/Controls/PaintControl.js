// import LevelCurve from '../model/LevelCurve';
import {BezierSpline} from '../model/Spline';
import {Floral} from '../model/stem';
import CurveManagement from '../model/CurveManagement';
import {drawFlower, drawStem, drawLeaf} from '../model/Drawer';
const error = 100;

function PaintControl(pannel) {
	let rawPointData = [];
	let paintingPolyLine = undefined;	

	this.start = function( point ) {
				// clearRawData();

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
		let pathString = smoothBizer.svgString();

		drawBasePath(pathString);

		// let floral = new LevelCurve( smoothBizer, 1, UI.state.levelCurve );
		// floral.draw();
		let floral = new Floral(smoothBizer);
		CurveManagement.floralScene.push(floral);
		drawFlower(floral);
		drawStem(floral);


		let stems = CurveManagement.growBranches();

		stems.forEach(s => drawLeaf(s) );
		clearRawData();
	};

	function updateLines(paintingPolyLine, rawPointData) {
		paintingPolyLine.plot( rawPointData );
	}
	
	function drawBasePath(pathString){
		CurveManagement.layer.debugCurveLayer.path( pathString ).fill('none').stroke({ width: 3 }).stroke('#f06');
	}
	function clearRawData(){
		rawPointData = [];
		paintingPolyLine.remove();
	}	
}

export default PaintControl;

