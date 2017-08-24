import Control from './Controls/Control';
import SVG from 'svg.js';
import * as UI from './model/UIManagement';
import CurveManagement from './model/CurveManagement';
import {changeColorMap} from './color/colorHex';

(function(){
	let draw = SVG('drawing').size(1900, 2500);
	setSVGLayer(draw);
	setControl(draw);
	// UI.setGUI();
	let bound = (new URL(window.location)).searchParams.get('bound');

	UI.setBounding(bound);
	exportFunctionToHtml();
	UI.setFrontFlowerState();
})();

function setSVGLayer(panel) {
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

function exportFunctionToHtml() {
	window['changeColorByIndex']= UI.changeColorByIndex;
	window['setSideFlowerState'] = UI.setSideFlowerState;
	window['setFrontFlowerState'] = UI.setFrontFlowerState;
	window['setLeafState'] = UI.setLeafState;
	window['download'] = UI.features.download;
}

function readSingleFile(e) {
	var file = e.target.files[0];
	if (!file) return;

	var reader = new FileReader();
	reader.onload = function(e) {
		let contents = e.target.result;
		// displayContents(contents);
		const colorMap = JSON.parse(contents);
		changeColorMap(colorMap);
		CurveManagement.draw();
		UI.changeColor(UI.state.color);
	};
	reader.readAsText(file);
}

function displayContents(contents) {
	var element = document.getElementById('file-content');
	element.textContent = contents;
}

document.getElementById('file-input')
  .addEventListener('change', readSingleFile, false);