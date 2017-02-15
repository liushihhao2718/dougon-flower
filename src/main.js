import PaintControl from './Controls/PaintControl';
import SVG from 'svg.js';
import * as dat from './lib/dat.gui';


(function(){
	let draw = SVG('drawing').size(1000, 1000);
	let gui = new dat.GUI();
	
	setControl(draw);

})();

function setControl(_container) {
	let isMouseDown = false;
	let currnetControl = new PaintControl(_container);
	const top = _container.node.getBoundingClientRect().top;
	const left = _container.node.getBoundingClientRect().left;

	_container.on('mousedown', function (e) {
		const point = [
			e.clientX - top,
			e.clientY - left
		];
		isMouseDown = true;
		currnetControl.start(point);

	});
	_container.on('mouseup', function () {
		isMouseDown = false;
		currnetControl.end();
	});
	_container.on('mousemove', function (e) {
		var x = e.offsetX;
		var y = e.offsetY;
		if (isMouseDown) {
			currnetControl.update([x, y]);
		}
	});
}

