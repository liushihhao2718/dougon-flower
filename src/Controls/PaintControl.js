import fitCurve from 'fit-curve';
import MagneticCurve from '../model/MagneticCurve';
import Bezier from 'bezier-js';


const error = 100;

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
		let pathString = fittedCurveToPathString(smoothBizer);


		drawLevelCurve(smoothBizer);
		// draw magnetic curve
		

		drawOnPannel(pannel, pathString);
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
	function drawLevelCurve(beziers){
		beziers.forEach(b =>{
			const temp_b = new Bezier(
				b[0][0], b[0][1],
				b[1][0], b[1][1],
				b[2][0], b[2][1],
				b[3][0], b[3][1]
			);
			drawAt(0.2, temp_b, 1);
			drawAt(0.4, temp_b, -1);
			drawAt(0.6, temp_b, 1);
			drawAt(0.8, temp_b, -1);
		});
	}
	function drawAt(t, b, sign){
		let start = b.get(t);
		let v = b.derivative(t);
		let mag = new MagneticCurve({
			startX: start.x,
			startY: start.y,
			vx: v.x,
			vy: v.y,
			T: 50,
			alpha: 0.65,
			sign: sign
		});
		mag.drawOn(pannel);
	}
}

export default PaintControl;

