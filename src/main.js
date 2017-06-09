import Control from './Controls/Control';
import SVG from 'svg.js';
import * as UI from './model/UIManagement';
import CurveManagement from './model/CurveManagement';

(function(){
	let draw = SVG('drawing').size(1500, 1000);
	setSVGLayer(draw);
	setControl(draw);
	UI.setGUI();
})();

function setSVGLayer(panel) {
	var image = panel.image('target.jpg');
	image.width(1371).height(436).x(60.5).y(63.5);
	//order is very importent!
	let dougonLayer = panel.group();
	let leafLayer = panel.group();
	let stemLayer = panel.group();
	let flowerLayer = panel.group();
	let debugCurveLayer = panel.group().hide();
	let hintLayer = panel.group();

	CurveManagement.panel = panel;
	CurveManagement.initSvgSymbol();
	CurveManagement.layer = {dougonLayer, leafLayer, stemLayer, flowerLayer, debugCurveLayer, hintLayer};
}

function setControl(_container) {
	let isMouseDown = false;
	let tools = {
		paint : new Control.PaintControl(_container),
		flower : new Control.FloralControl(_container),
		select : new Control.SelectControl(_container),
		skeleton: new Control.SkeletonControl(_container)
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

