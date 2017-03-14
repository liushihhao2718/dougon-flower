import LevelCurve from '../model/LevelCurve';
import * as UI from '../model/UIManagement';
let rect;
function BoundControl(pannel) {
	this.start = function( point ) {
		rect = pannel.rect().fill('#524B61');
		rect.x(point[0]).y(point[1]);
	};
	this.update = function( point ) {
		rect.size(point[0] - rect.x(), point[1] - rect.y());
	};

	this.end = function() {
		UI.state.bound = {
			x: rect.x(),
			y: rect.y(),
			width: rect.width(),
			height: rect.height()
		};
	};
}

export default BoundControl;

