import * as dat from '../lib/dat.gui';

let gui, folders = [];
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

export function setGUI(){
	gui = new dat.GUI();

	gui.add(state, 'trunkWidth', 1, 20);
	gui.add(state, 'intersect');

	let f0 = gui.addFolder('Level 0');
	f0.add(state.levelCurve[0], 'length');
	f0.add(state.levelCurve[0], 'alpha');
	f0.add(state.levelCurve[0], 'branches');
	folders.push(f0);

	let f1 = gui.addFolder('Level 1');
	f1.add(state.levelCurve[1], 'length');
	f1.add(state.levelCurve[1], 'alpha');
	f1.add(state.levelCurve[1], 'branches');
	folders.push(f1);
}
