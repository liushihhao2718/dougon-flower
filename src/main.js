import Control from './Controls/Control';
import SVG from 'svg.js';
import * as UI from './model/UIManagement';
import CurveManagement from './model/CurveManagement';

(function(){
	let draw = SVG('drawing').size('100%', '100%').viewbox(0,0,1900, 2500);
	setSVGLayer(draw);
	setControl(draw);
	UI.setGUI();
})();

function setSVGLayer(panel) {
	// var image = panel.image('target.jpg');
	// image.width(1371).height(436).x(60.5).y(63.5);

	//order is very importent!
	let dougonLayer = panel.group();
	let leafLayer = panel.group();
	let stemLayer = panel.group();
	let flowerLayer = panel.group();
	let debugCurveLayer = panel.group().hide();
	let hintLayer = panel.group();

	CurveManagement.panel = panel;
	CurveManagement.initSvgSymbol();
	CurveManagement.layer = { dougonLayer, leafLayer, stemLayer, flowerLayer, debugCurveLayer, hintLayer};
}

function setControl(_container) {
	let isMouseDown = false;
	let tools = {
		paint : new Control.PaintControl(_container),
		flower : new Control.FloralControl(_container),
		select : new Control.SelectControl(_container),
		skeleton: new Control.SkeletonControl(_container)
	};
		
	_container.on('mousedown', function (e) {
		let currnetControl = tools[UI.state.tool];
		const point = getSvgCtmPoint(e, _container);

		isMouseDown = true;
		currnetControl.start(point);
	});
	_container.on('mouseup', function () {
		let currnetControl = tools[UI.state.tool];

		isMouseDown = false;
		currnetControl.end();
	});
	_container.on('mousemove', e => {
		let currnetControl = tools[UI.state.tool];
		
		const [x,y] = getSvgCtmPoint(e, _container);
		if (isMouseDown) {
			currnetControl.update([x, y]);
		}
	});
}

function getSvgCtmPoint(e, svg){
	const p = svg.node.createSVGPoint();
	p.x = e.clientX;
	p.y = e.clientY;
	const ctm = svg.node.getScreenCTM();
	const inverse = ctm.inverse();
	const transforomP = p.matrixTransform(inverse);
	return [transforomP.x, transforomP.y];
}