import * as dat from '../lib/dat.gui';
import download from '../lib/download';
import CurveManagement from './CurveManagement';
import {dougonBounding, dougonBoundingParam} from '../images/dougonBounding';

const styleMap = require('../color/StyleMap.json');
const colorMap = require('../color/colorHex.json');

import {Floral} from './stem';
import {BezierSpline} from './Spline';

let gui, folders = [];
let controls = [];
let leafControls = [];
export let state = {
	flowerSize: 150,
	trunkHead: 1,
	trunkTail: 1,
	density: 1,
	levelCurve :[
		{
			length: 80,
			alpha: 0.9,
			branches: 5
		},
		{
			length: 62,
			alpha: 0.9,
			branches: 5
		},
		{
			length: 40,
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
	},
	test:function() {
		const error = 50;
		let rawPointData = [[204,433],[204,431],[208,421],[212,406],[223,392],[233,376],[242,361],[249,353],[258,342],[265,335],[271,327],[278,321],[285,315],[291,312],[299,306],[305,304],[308,302],[312,300],[318,297],[328,294],[355,286],[384,277],[422,273],[451,273],[474,271],[501,270],[524,268],[548,268],[570,268],[591,268],[604,268],[614,268],[625,271],[639,274],[653,278],[665,284],[673,286],[677,289],[683,293],[689,296],[696,300],[704,304],[713,310],[728,316],[741,321],[752,326],[758,329],[762,332],[765,333],[768,335],[770,337],[774,339],[779,342],[785,345],[792,349],[800,354],[807,356],[812,357],[816,359],[823,362],[832,366],[845,368],[859,373],[870,375],[877,377],[880,378],[883,379],[886,380],[888,382],[894,383],[900,385],[908,388],[919,390],[925,392],[932,392],[937,392],[944,392],[949,392],[956,392],[966,392],[974,392],[980,392],[985,392],[991,392],[999,392],[1009,392],[1014,392],[1021,392],[1027,392],[1036,392],[1045,392],[1052,392],[1061,392],[1075,392],[1099,392],[1120,392],[1132,392],[1139,392],[1140,391],[1145,390],[1152,387],[1160,384],[1167,378],[1180,375],[1189,372],[1199,367],[1220,360],[1236,355],[1248,350],[1253,346],[1257,341],[1260,336],[1262,333],[1263,330],[1264,328],[1265,325],[1265,323],[1267,320],[1268,318],[1268,316],[1269,316],[1269,315],[1269,314],[1269,313],[1269,312],[1270,311],[1270,310],[1270,309]];		

		state.tools.paint.start(rawPointData[0]);
		rawPointData.slice(1, rawPointData.length-1).forEach(p=>{
			state.tools.paint.update(p);	
		});
		state.tools.paint.end();
		// let aspect = state.aspect;
		// CurveManagement.floralScene.push( new Floral(smoothBizer,state.flowerSize, state.trunkHead, state.trunkTail,'海石榴華', aspect) );
		
		// CurveManagement.draw();
	}
};

export function setGUI(){
	gui = new dat.GUI();
	let bound =	gui.add(state, 'bound', Object.keys(dougonBounding));

	let c0 = gui.add(state, 'tool', ['paint', 'flower', 'select', 'skeleton']);
	let c1 = gui.add(state, 'trunkHead', 1, 20);
	let c2 = gui.add(state, 'trunkTail', 1, 40);
	let c3 = gui.add(state, 'flowerSize', 10, 200);
	let c4 = gui.add(state, 'density', 0.5, 2);
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
	c4.onChange(()=>{
		CurveManagement.draw();
		changeColor(state.color);
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

	gui.add(features, 'test');

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
	
	CurveManagement.layer.dougonLayer.svg(svgString);
	
	state.flowerSize = dougonBoundingParam[value].flowerSize;
	state.trunkHead = dougonBoundingParam[value].trunkHead;
	state.trunkTail = dougonBoundingParam[value].trunkTail;
	gui.__controllers[2].updateDisplay();
	gui.__controllers[3].updateDisplay();
	gui.__controllers[4].updateDisplay();

	changeColor( state.color );
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
