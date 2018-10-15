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
	trunkHead: 5,
	trunkTail: 10,
	density: 1,
	rotation: 45,
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
	bound: '圖樣',
	tool:'paint',
	show:{
		image:false,
		dougon:true
	},
	color:'鋪地捲成',
	aspect: '正面',
	leafBranchType: 'test'
};

let features = {
	download : function(){
		let svg = document.getElementsByTagName('svg')[0].cloneNode(true);
		let image = svg.getElementsByTagName('image')[0];
		if(image) image.remove();

		let bounding = document.getElementsByTagName('svg')[0].children[1].children[0].children[1];
		let bbox = bounding.getBBox();
		let {x, y} = convertCoords(bbox.x, bbox.y, document.getElementsByTagName('svg')[0], bounding);
		svg.setAttribute('viewBox', `${x} ${y} ${bbox.width} ${bbox.height}`);
		svg.setAttribute('width', bbox.width);
		svg.setAttribute('height', bbox.height);
		let string = svg.outerHTML;

		download(string, 'file.svg', 'text/plain');
	},
	toggleLayer: function(layer){
		layer.visible() ? layer.hide() : layer.show();
	},
	test:function() {
		const error = 50;
		let rawPointData = [[204,433],[204,431],[208,421],[212,406],[223,392],[233,376],[242,361],[249,353],[258,342],[265,335],[271,327],[278,321],[285,315],[291,312],[299,306],[305,304],[308,302],[312,300],[318,297],[328,294],[355,286],[384,277],[422,273],[451,273],[474,271],[501,270],[524,268],[548,268],[570,268],[591,268],[604,268],[614,268],[625,271],[639,274],[653,278],[665,284],[673,286],[677,289],[683,293],[689,296],[696,300],[704,304],[713,310],[728,316],[741,321],[752,326],[758,329],[762,332],[765,333],[768,335],[770,337],[774,339],[779,342],[785,345],[792,349],[800,354],[807,356],[812,357],[816,359],[823,362],[832,366],[845,368],[859,373],[870,375],[877,377],[880,378],[883,379],[886,380],[888,382],[894,383],[900,385],[908,388],[919,390],[925,392],[932,392],[937,392],[944,392],[949,392],[956,392],[966,392],[974,392],[980,392],[985,392],[991,392],[999,392],[1009,392],[1014,392],[1021,392],[1027,392],[1036,392],[1045,392],[1052,392],[1061,392],[1075,392],[1099,392],[1120,392],[1132,392],[1139,392],[1140,391],[1145,390],[1152,387],[1160,384],[1167,378],[1180,375],[1189,372],[1199,367],[1220,360],[1236,355],[1248,350],[1253,346],[1257,341],[1260,336],[1262,333],[1263,330],[1264,328],[1265,325],[1265,323],[1267,320],[1268,318],[1268,316],[1269,316],[1269,315],[1269,314],[1269,313],[1269,312],[1270,311],[1270,310],[1270,309]];
		let smoothBizer = BezierSpline.makeByPoints( rawPointData, error );
		
		let aspect = state.aspect;
		CurveManagement.floralScene.push( new Floral(smoothBizer,state.flowerSize, state.trunkHead, state.trunkTail,'海石榴華', aspect) );
		
		CurveManagement.draw();
	},
	toggleImage: function(){
		let image = document.getElementsByTagName('image')[0];
		if (image)
			image.style.visibility = state.show.image? 'visible':'hidden';
	},
	toggleDougon: function(){
		Array.from(CurveManagement.layer.dougonLayer.node.children[0].children)
		.filter(x => x.nodeName != 'image')
		.forEach(x => {
			x.style.display = '';
			x.style.visibility = state.show.dougon? 'visible':'hidden';
		});
	},
	whiteColor : whiteColor
};

export function setGUI(){
	gui = new dat.GUI();
	let bound =	gui.add(state, 'bound', Object.keys(dougonBounding));

	let c0 = gui.add(state, 'tool', ['paint', 'flower', 'select', 'skeleton']);
	let c1 = gui.add(state, 'trunkHead', 1, 20);
	let c2 = gui.add(state, 'trunkTail', 1, 40);
	let c3 = gui.add(state, 'flowerSize', 10, 200);
	let c4 = gui.add(state, 'density', 0.5, 2);
	let c5 = gui.add(state, 'rotation', -180, 180);	

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
	c5.onChange( rotation =>{
		for(let floral of CurveManagement.selectedCurve){
			floral.flowerRotation = rotation;	
		}
		
		CurveManagement.draw();
		changeColor( state.color );
	});
	levelFolder(0);
	levelFolder(1);
	levelFolder(2);

	setOnChange(controls, leafControls);

	let colorControl = gui.add(state, 'color', ['青緣紅地','青緣紅地-地色紅粉','綠緣青地','綠緣青地-地色為白', '鋪地捲成', '素描']);
	colorControl.onChange(value => changeColor(value));


	gui.add(state, 'leafBranchType', ['big', 'test']);
	gui.add(state, 'aspect', ['正面', '側面']);

	gui.add(features, 'download');
	bound.onChange(value => setBounding(value) );
	setBounding(state.bound);

	gui.add(features, 'test');
	gui.add(features, 'whiteColor');
	let folder = gui.addFolder('Layer');

	// folder.add(state.show, 'image')
	// .onChange(()=>{
	// 	features.toggleImage();
	// });

	folder.add(state.show, 'dougon')
	.onChange(()=>{
		features.toggleDougon();
	});	
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
	leafControls.push( folder.add(state.levelCurve[index], 'length',1,90) );
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
	
	features.toggleImage();

	state.flowerSize = dougonBoundingParam[value].flowerSize;
	state.trunkHead = dougonBoundingParam[value].trunkHead;
	state.trunkTail = dougonBoundingParam[value].trunkTail;
	state.levelCurve[0].length = dougonBoundingParam[value].level1Size || 80;

	gui.__controllers[2].updateDisplay();
	gui.__controllers[3].updateDisplay();
	gui.__controllers[4].updateDisplay();

	changeColor( state.color );
}
function convertCoords(x,y, svgDoc, elem) {

	let offset = svgDoc.getBoundingClientRect();

	let matrix = elem.getScreenCTM();

	return {
		x:(matrix.a * x) + (matrix.c * y) + matrix.e - offset.left,
		y:(matrix.b * x) + (matrix.d * y) + matrix.f - offset.top
	};
}

export function changeColor(value) {
	let colorTags = styleMap['五彩遍裝'][value];
	Object.keys( colorTags ).forEach(key =>{
		changeBGColor( key, colorMap[colorTags[key]] );
	});
}
function changeBGColor(class_name,new_color) {
	const reg = /[A-Z]6/;
	const cols = Array.from(document.getElementsByClassName(class_name));

	if( reg.test(class_name) ){
		for(let graph of cols) {
			graph.style.stroke = new_color;
		}
	}
	else{
		for(let graph of cols) {
			graph.style.fill = new_color;
		}
	}	
}

function whiteColor() {
	const path = Array.from(document.getElementsByTagName('path'));
	for(let tag of path) {
		tag.style.stroke = 'black';
		tag.style.fill = 'white';
		tag.style.strokeWidth = '1px';
	}
}