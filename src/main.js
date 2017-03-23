import PaintControl from './Controls/PaintControl';
import BoundControl from './Controls/BoundControl';
import SVG from 'svg.js';
import * as UI from './model/UIManagement';
import CurveManagement from './model/CurveManagement';

(function(){
	let draw = SVG('drawing').size(1000, 1000);
	UI.setGUI();
	setSVGLayer(draw);
	setControl(draw);

})();

function setSVGLayer(pannel) {
	let drawingLayer = pannel.group();
	let leafLayer = pannel.group();
	let debugCurveLayer = pannel.group();
	let stemLayer = pannel.group();
	let flowerLayer = pannel.group();
	CurveManagement.layer = {drawingLayer, leafLayer, stemLayer, flowerLayer, debugCurveLayer};

	// order is important
	pannel.add( drawingLayer );
	pannel.add( leafLayer );
	pannel.add( debugCurveLayer );
	pannel.add( stemLayer );
	pannel.add( flowerLayer );
}

function setControl(_container) {
	let isMouseDown = false;
	let tools = {
		'paint':new PaintControl(_container),
		'bound':new BoundControl(_container),
		'select': undefined
	};

	const top = _container.node.getBoundingClientRect().top;
	const left = _container.node.getBoundingClientRect().left;

	_container.on('mousedown', function (e) {
		let currnetControl = tools[UI.state.tool];

		const point = [
			e.clientX - top,
			e.clientY - left
		];
		isMouseDown = true;
		currnetControl.start(point);

	});
	_container.on('mouseup', function () {
		let currnetControl = tools[UI.state.tool];

		isMouseDown = false;
		currnetControl.end();
	});
	_container.on('mousemove', function (e) {
		let currnetControl = tools[UI.state.tool];

		var x = e.offsetX;
		var y = e.offsetY;
		if (isMouseDown) {
			currnetControl.update([x, y]);
		}
	});
}

