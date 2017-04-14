import * as UI from './UIManagement';
import * as Collision from './Collision';
import * as Drawer from '../model/Drawer';

export default {
	selectedCurve :[],
	scene :[],
	floralScene:[],
	layer:{},
	growBranches(){
		let leafDrawingQueue=[];
		let leafCollisionScene=[];
		const amount = UI.state.levelCurve[0].branches;
		let groupedBurgeons = this.floralScene.map( f=>f.burgeons(amount) );
		let burgeons = flatten( groupedBurgeons );
		let sign = 1;
		for(let levelParam of UI.state.levelCurve){
			let nextLevelBurgeons = [];

			while(burgeons.length > 0){

				let burgeon = burgeons[0];
				burgeons.shift();

				let leaf = burgeon.germinate(levelParam.length,levelParam.alpha, sign=-1*sign);

				if( Collision.testCollision(leaf.colliders, leafCollisionScene, burgeon.parent.colliders )){
					nextLevelBurgeons.push(burgeon);
				}
				else{
					leafDrawingQueue.push(leaf);
					leafCollisionScene.push(leaf.colliders);
					burgeons = burgeons.concat( leaf.burgeons( levelParam.branches ) );
				}
			}
		}
		

		return leafDrawingQueue;
	},
	clearAllLayer(){
		Object.values(this.layer).forEach(layer=> layer.clear());
	},
	clearDrawing(){
		this.layer.flowerLayer.clear();
		this.layer.stemLayer.clear();
		this.layer.leafLayer.clear();
		this.layer.debugCurveLayer.clear();
		this.layer.hintLayer.clear();
	},
	draw(){
		this.clearDrawing();
		this.floralScene.forEach(floral=>{

			Drawer.drawFlower(floral);
			Drawer.drawStem(floral);
			Drawer.drawBasePath(floral.curve.svgString());
			
			let stems = this.growBranches().reverse();

			stems.forEach(s => Drawer.drawLeaf(s) );
			stems.forEach(s => Drawer.drawPolygon(s.colliders) );
		});
	},
	redrawFlower(){
		this.layer.flowerLayer.clear();
		this.floralScene.forEach(floral=>{

			Drawer.drawFlower(floral);
			Drawer.drawStem(floral);
			Drawer.drawBasePath(floral.curve.svgString());
		});
	},
	drawHint(){
		this.layer.hintLayer.clear();
		const position = this.selectedCurve.flowerPosition;
		this.layer.hintLayer.circle(position.r * 2)
			.cx(position.x).cy(position.y)
			.fill('none').stroke({color:'red'});
	},
	clearScene(){
		this.floralScene.length = [];
	}
};

function flatten(arr){
	return arr.reduce((acc, val) => acc.concat( Array.isArray(val) ? flatten(val) : val),[]);
}