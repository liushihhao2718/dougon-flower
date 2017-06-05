import SVG from 'svg.js';
import concaveman from 'concaveman';
import _ from 'lodash';
import {LeafImage, LeafBranch} from '../src/images/LeafImage';


function getPathPoints(svg, draw) {
	let commands = svg.select('path').members.map(p=>{
		return	draw.path(p.node.getAttribute('d')).array().value;
	});

	return _(commands).flatten().map(cmd=>
		_(cmd).filter(x => typeof x === 'number').chunk(2).value()
	).flatten().value();
}

function concaveSVG(flowerString, concave, length){
	let draw = SVG('drawing').size(300, 300);

	let flower = draw.svg(flowerString);
	let points = getPathPoints(flower, draw).filter(x=>typeof x[0] === 'number' && typeof x[1] === 'number');
	let polygon = concaveman(points, concave, length);
	draw.polygon().plot(polygon).fill('none').stroke({ width: 3 }).stroke('red');

	return polygon;
}

function 正面(){
	let flowerString = require('../src/images/海石榴心_v3.svg');
	concaveSVG(flowerString, 1, 23);
}

正面();

function 側面(){
	let flowerString = require('../src/images/sideFlower_v6.svg');
	let draw = SVG('drawing').size(300, 300);

	let flower = draw.svg(flowerString);
	let points = getPathPoints(flower, draw).filter(x=>typeof x[0] === 'number' && typeof x[1] === 'number');
	let polygon = concaveman(points, 1.6, 55);
	SVG.adopt(flower.select('#SvgjsG1684').members[0].node).polygon().plot(polygon).fill('none').stroke({ width: 3 }).stroke('red');
}

側面();

function leaf(){
	let json = [];
	_.range(6).forEach(i =>{
		const polygon = concaveSVG(LeafImage[i], 2, 45);
		json.push(polygon);
	});
	console.log(JSON.stringify(json));
}
leaf();

function leafBranch(){
	let json = [];
	_.range(7).forEach(i =>{
		const polygon = concaveSVG(LeafBranch[i], 2, 45);
		json.push(polygon);
	});
	console.log(JSON.stringify(json));
}
leafBranch();