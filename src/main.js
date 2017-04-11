import PaintControl from './Controls/PaintControl';
import BoundControl from './Controls/BoundControl';
import SVG from 'svg.js';
import * as UI from './model/UIManagement';
import CurveManagement from './model/CurveManagement';

(function(){
	let draw = SVG('drawing').size(1000, 1000);
	setSVGLayer(draw);
	setControl(draw);
	UI.setGUI();
})();

function setSVGLayer(panel) {
	let drawingLayer = panel.group();
	let leafLayer = panel.group();
	let stemLayer = panel.group();
	let flowerLayer = panel.group();
	let debugCurveLayer = panel.group().hide();

	CurveManagement.layer = {drawingLayer, leafLayer, stemLayer, flowerLayer, debugCurveLayer};
}

function setControl(_container) {
	let isMouseDown = false;
	let tools = {
		paint : new PaintControl(_container),
		bound : new BoundControl(_container),
		select : undefined
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

