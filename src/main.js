import GlyphsControl from './Controls/GlyphsControl';
import PaintControl from './Controls/PaintControl';
import SVG from 'svg.js';

var draw = SVG('drawing').size(300, 300);


setControl(draw);

function setControl(_container) {
	let isMouseDown = false;
	let currnetControl = new PaintControl(draw);

	_container.on('mousedown', function (e) {
		const point = [
			e.clientX,
			e.clientY
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

