import SVG from 'svg.js';
import fitCurve from 'fit-curve';

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
		let pathString = fittedCurveDataToPathString(smoothBizer);
		pannel.path( pathString ).fill('none').stroke({ width: 3 }).stroke('#f06');
		rawPointData.length = 0;
	};
}

export default PaintControl;

function updateLines(paintingPolyLine, rawPointData) {
	paintingPolyLine.plot( rawPointData );
}
function fittedCurveDataToPathString(fittedLineData) {
	var str = '';
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