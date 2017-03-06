import fitCurve from 'fit-curve';
import LevelCurve from '../model/LevelCurve';
import * as UI from '../model/UIManagement';

const error = 500;

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
		let smoothBizer = fitCurve( rawPointData, error );
		if(smoothBizer.length == 0) {
			clearRawData();
			return;
		}
		let pathString = fittedCurveToPathString(smoothBizer);

		drawOnPannel(pannel, pathString);

		let lvCurve = new LevelCurve( smoothBizer, 1, UI.state.levelCurve);
		lvCurve.drawOn(pannel);

		clearRawData();
	};

	function updateLines(paintingPolyLine, rawPointData) {
		paintingPolyLine.plot( rawPointData );
	}
	function fittedCurveToPathString(fittedLineData) {
		var str = '';
		//bezier : [ [c0], [c1], [c2], [c3] ]
		fittedLineData.map(function (bezier, i) {
			if (i == 0) {
				str += 'M ' + bezier[0][0] + ' ' + bezier[0][1];
			}

			str += 'C ' + bezier[1][0] + ' ' + bezier[1][1] + ', ' +
			bezier[2][0] + ' ' + bezier[2][1] + ', ' +
			bezier[3][0] + ' ' + bezier[3][1] + ' ';	
						
		});

		return str;
	}
	function drawOnPannel(pannel, pathString){
		pannel.path( pathString ).fill('none').stroke({ width: 3 }).stroke('#f06');
	}
	function clearRawData(){
		rawPointData.length = 0;
		paintingPolyLine.remove();
	}	
}

export default PaintControl;

