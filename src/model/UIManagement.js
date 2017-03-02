import * as dat from '../lib/dat.gui';
import download from '../lib/download';
import CurveManagement from './CurveManagement';

let gui, folders = [];
let controls = [];
export let state = {
	trunkWidth: 15,
	intersect: false,
	levelCurve :[
		{
			length: 100,
			alpha: 0.9,
			branches: 5
		},
		{
			length: 20,
			alpha: 0.8,
			branches: 5
		}
	]
};

let features = {
	download : function(){
		let svg = document.getElementsByTagName('svg')[0];
		download(svg.outerHTML, 'file.svg', 'text/plain');
	}
};

export function setGUI(){
	gui = new dat.GUI();

	gui.add(state, 'trunkWidth', 1, 20);
	gui.add(state, 'intersect');

	levelFolder(0);
	levelFolder(1);
	setOnChange(controls);

	gui.add(features, 'download');

}

function levelFolder(index){
	let folder = gui.addFolder(`Level ${index}`);
	controls.push( folder.add(state.levelCurve[index], 'length') );
	controls.push( folder.add(state.levelCurve[index], 'alpha') );
	controls.push( folder.add(state.levelCurve[index], 'branches') );
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