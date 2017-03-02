import LevelCurve from '../model/LevelCurve';
import * as UI from '../model/UIManagement';

function SelectControl(pannel) {
	this.start = function( point ) {
		rawPointData.push( point );
		paintingPolyLine = pannel.polyline().fill('none').stroke({ width: 1 });

	};
	this.update = function( point ) {
		
	};

	this.end = function() {
		
	};
}

export default SelectControl;

