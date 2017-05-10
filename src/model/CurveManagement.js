import * as UI from './UIManagement';
import * as Collision from './Collision';
import * as Drawer from '../model/Drawer';

export default {
	selectedCurve :[],
	scene :[],
	floralScene:[],
	layer:{},
	debug_burgeons:[],
	growBranches(){
		let leafDrawingQueue=[];
		let flowerPosition = this.floralScene.map(f=>f.flowerPosition).filter(f => f!==undefined);
		let collisionScene= this.floralScene.map(f=>f.colliders).filter(f => f!==undefined);
		const amount = UI.state.levelCurve[0].branches;
		let groupedBurgeons = this.floralScene.map( f=>f.burgeons(amount) ).filter(f => f!==undefined);
		let burgeons = flatten( groupedBurgeons );
		let sign = 1;
		UI.state.levelCurve.forEach((levelParam, index)=>{

			let nextLevelBurgeons = [];

			while(burgeons.length > 0){

				let burgeon = burgeons[0];
				burgeons.shift();

				let leaf = burgeon.germinate(levelParam.length,levelParam.alpha, sign=-1*sign);

				if(Collision.testCollision(leaf.colliders, collisionScene, flowerPosition, burgeon.parent.colliders)){
					nextLevelBurgeons.push(burgeon);
				}
				else{
					leafDrawingQueue.push(leaf);
					collisionScene.push(leaf.colliders);
					this.debug_burgeons.push({leaf,level:index});
					burgeons = burgeons.concat( leaf.burgeons( levelParam.branches ) );
				}
			}
			burgeons = burgeons.concat( nextLevelBurgeons );
		});
		

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
		
		this.floralScene.forEach(floral=>floral.draw());
		
		let stems = this.growBranches().reverse();
		stems.forEach(s => Drawer.drawLeaf(s) );
		stems.forEach(s => Drawer.drawPolygon(s.colliders) );
		
		this.debug_burgeons.forEach(({leaf, level}) => Drawer.drawMagneticCurve(leaf, level));
		this.drawHint();
	},
	drawHint(){
		this.layer.hintLayer.clear();

		for(let floral of this.selectedCurve){
			const position = floral.flowerPosition;

			this.layer.hintLayer.circle(position.r * 2)
				.cx(position.x).cy(position.y)
				.fill('none').stroke({color:'orange', width:10});
		}
		
	},
	clearScene(){
		this.floralScene.length = [];
	}
};

function flatten(arr){
	return arr.reduce((acc, val) => acc.concat( Array.isArray(val) ? flatten(val) : val),[]);
}

