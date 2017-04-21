import * as dat from '../lib/dat.gui';
import download from '../lib/download';
import CurveManagement from './CurveManagement';
import {dougonBounding, dougonBoundingParam} from '../images/dougonBounding';

const styleMap = require('../color/StyleMap.json');
const colorMap = require('../color/colorHex.json');

let gui, folders = [];
let controls = [];
let leafControls = [];
export let state = {
	flowerSize: 100,
	trunkHead: 5,
	trunkTail: 30,
	levelCurve :[
		{
			length: 100,
			alpha: 0.9,
			branches: 3
		},
		{
			length: 62,
			alpha: 0.9,
			branches: 5
		},
		{
			length: 38,
			alpha: 0.9,
			branches: 2
		}
	],
	bound: '撩擔方',
	tool:'paint',
	show:{
	},
	color:'鋪地捲成',
	aspect: '正面'
};

let features = {
	download : function(){
		let svg = document.getElementsByTagName('svg')[0];
		download(svg.outerHTML, 'file.svg', 'text/plain');
	},
	toggleLayer: function(layer){
		layer.visible() ? layer.hide() : layer.show();
	}
};

export function setGUI(){
	gui = new dat.GUI();
	let bound =	gui.add(state, 'bound', Object.keys(dougonBounding));

	let c0 = gui.add(state, 'tool', ['paint', 'flower', 'select']);
	let c1 = gui.add(state, 'trunkHead', 1, 20);
	let c2 = gui.add(state, 'trunkTail', 5, 40);
	let c3 = gui.add(state, 'flowerSize', 10, 200);

	// control/s.push(c0);
	// controls.push(c1);
	// controls.push(c2);
	c1.onChange(head =>{
		for(let floral of CurveManagement.selectedCurve){
			floral.trunkHead = head;	
		}
		
		CurveManagement.draw();
		changeColor( state.color );
	});
	c2.onChange(tail =>{
		for(let floral of CurveManagement.selectedCurve){
			floral.trunkTail = tail;	
		}
		
		CurveManagement.draw();
		changeColor( state.color );
	});
	c3.onChange(r =>{
		for(let floral of CurveManagement.selectedCurve){
			floral.flowerPosition.r = r;	
		}
		
		CurveManagement.draw();
		changeColor( state.color );
	});

	levelFolder(0);
	levelFolder(1);
	levelFolder(2);

	setOnChange(controls, leafControls);

	let colorControl = gui.add(state, 'color', ['青緣紅地', '綠緣青地', '鋪地捲成']);
	colorControl.onChange(value => changeColor(value));

	gui.add(state, 'aspect', ['正面', '側面']);

	gui.add(features, 'download');
	bound.onChange(value => setBounding(value) );
	setBounding(state.bound);


	let folder = gui.addFolder('Layer');
	Object.keys(CurveManagement.layer).forEach(key => {
		let layer = CurveManagement.layer[key];
		state.show[key] = layer.visible();
		let control = folder.add(state.show, key);
		control.onChange(()=>{
			features.toggleLayer(layer);
		});
	});
}

function levelFolder(index){
	let folder = gui.addFolder(`Level ${index}`);
	leafControls.push( folder.add(state.levelCurve[index], 'length',10,500) );
	leafControls.push( folder.add(state.levelCurve[index], 'alpha') );
	leafControls.push( folder.add(state.levelCurve[index], 'branches').step(1) );
	folders.push(folder);
}

function setOnChange(controls, leafControls){
	controls.forEach( c => {
		c.onChange( () => {
			CurveManagement.redrawFlower();
			changeColor( state.color );
		});
	});
	leafControls.forEach( c => {
		c.onChange( () => {
			CurveManagement.draw();
			changeColor( state.color );
		});
	});
}

function setBounding(value){
	let svgString = dougonBounding[value];
	CurveManagement.clearAllLayer();
	CurveManagement.clearScene();

	if(value === '撩擔方'){
		CurveManagement.layer.dougonLayer.svg(svgString);
	}else
		CurveManagement.layer.dougonLayer.path(svgString).addClass('background');

	let backgroundColorKey = styleMap['五彩遍裝'][state.color].background;
	state.flowerSize = dougonBoundingParam[value].flowerSize;
	state.trunkHead = dougonBoundingParam[value].trunkHead;
	state.trunkTail = dougonBoundingParam[value].trunkTail;
	gui.__controllers[2].updateDisplay();
	gui.__controllers[3].updateDisplay();
	gui.__controllers[4].updateDisplay();
	changeBGColor('background',colorMap[backgroundColorKey]);
}

function changeColor(value) {
	let colorTags = styleMap['五彩遍裝'][value];
	Object.keys( colorTags ).forEach(key =>{
		changeBGColor( key, colorMap[colorTags[key]] );
	});
}
function changeBGColor(class_name,new_color) {
	const reg = /[A-Z]6/;
	const cols = document.getElementsByClassName(class_name);

	if( reg.test(class_name) ){
		for(const graph of cols) {
			graph.style.stroke = new_color;
		}
	}
	else{
		for(const graph of cols) {
			graph.style.fill = new_color;
		}
	}	
}
