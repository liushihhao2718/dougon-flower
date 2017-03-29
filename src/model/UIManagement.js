import * as dat from '../lib/dat.gui';
import download from '../lib/download';
import CurveManagement from './CurveManagement';
import dougonBounding from '../images/dougonBounding';

let gui, folders = [];
let controls = [];
export let state = {
	trunkHead: 5,
	trunkTail: 30,
	intersect: false,
	levelCurve :[
		{
			length: 100,
			alpha: 0.9,
			branches: 10
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
		},
		// {
		// 	length: 25,
		// 	alpha: 0.9,
		// 	branches: 5
		// }
	],
	bound: '令栱',
	inside: [],
	tool:'paint',
	'show':{
		'debugCurveLayer':true	
	}
};

let features = {
	download : function(){
		let svg = document.getElementsByTagName('svg')[0];
		download(svg.outerHTML, 'file.svg', 'text/plain');
	},
	toggleLayer: function(){
		let layer = CurveManagement.layer.debugCurveLayer;
		layer.visible() ? layer.hide() : layer.show();
	}
};

export function setGUI(){
	gui = new dat.GUI();
	let c0 = gui.add(state, 'tool', ['paint', 'bound', 'select']);
	let bound =	gui.add(state, 'bound', Object.keys(dougonBounding));
	bound.onChange(value => setBounding(value) );
	setBounding(state.bound);
	let c1 = gui.add(state, 'trunkHead', 1, 20);
	let c2 = gui.add(state, 'trunkTail', 20, 40);
	
	controls.push(c0);
	controls.push(c1);
	controls.push(c2);
	//gui.add(state, 'intersect');


	levelFolder(0);
	levelFolder(1);
	levelFolder(2);
	// levelFolder(3);

	setOnChange(controls);

	gui.add(features, 'download');
	gui.add(features, 'toggleLayer');

}

function levelFolder(index){
	let folder = gui.addFolder(`Level ${index}`);
	controls.push( folder.add(state.levelCurve[index], 'length',10,500) );
	controls.push( folder.add(state.levelCurve[index], 'alpha') );
	controls.push( folder.add(state.levelCurve[index], 'branches').step(1) );
	folders.push(folder);
}

function setOnChange(controls){
	controls.forEach( c => {
		c.onChange( () => {
			if( CurveManagement.selectedCurve.length === 1 ){
				CurveManagement.selectedCurve[0].redraw();
			}
		});
	});
}

function setBounding(value){
	let svgString = dougonBounding[value];
	CurveManagement.layer.drawingLayer.clear();
	CurveManagement.layer.drawingLayer.path(svgString).fill('#524B61');
}
