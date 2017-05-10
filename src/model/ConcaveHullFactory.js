import SVG from 'svg.js';
import concaveman from 'concaveman';
import flowerString from '../images/海石榴心_v3.svg';
import _ from 'lodash';


var draw = SVG('drawing').size(1000, 1000);

let flower = draw.svg(flowerString);
let points = getPathPoints(flower).filter(x=>typeof x[0] === 'number' && typeof x[1] === 'number');
let polygon = concaveman(points, Infinity, 1);
flower.hide();
draw.polygon().plot(polygon).fill('none').stroke({ width: 3 }).stroke('red');
console.log(polygon);


function getPathPoints(svg) {
	let commands = svg.select('path').members.map(p=>{
		return	draw.path(p.node.getAttribute('d')).array().value;
	});

	return _(commands).flatten().map(cmd=>
		_(cmd).filter(x => typeof x === 'number').chunk(2).value()
	).flatten().value();
}