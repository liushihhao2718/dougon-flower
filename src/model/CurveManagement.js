import * as UI from './UIManagement';
import * as Collision from './Collision';

export default {
	selectedCurve :[],
	scene :[],
	leafCollisionScene:[],
	floralScene:[],
	layer:{},
	leafDrawingQueue:[],
	growBranches(){
		const amount = UI.state.levelCurve[0].branches;
		let groupedBurgeons = this.floralScene.map( f=>f.burgeons(amount) );
		let burgeons = flatten( groupedBurgeons );
		let sign = 1;
		for(let levelParam of UI.state.levelCurve){
			let nextLevelBurgeons = [];

			while(burgeons.length > 0){

				let burgeon = burgeons[0];
				burgeons.shift();

				let leaf = burgeon.germinate(levelParam.length,levelParam.curve, sign=-1*sign);

				if( Collision.testCollision(leaf.colliders, this.leafCollisionScene, burgeon.parent.colliders )){
					nextLevelBurgeons.push(burgeon);
				}
				else{
					this.leafDrawingQueue.push(leaf);
					this.leafCollisionScene.push(leaf.colliders);
					burgeons = burgeons.concat( leaf.burgeons() );
				}
			}
		}
		

		return this.leafDrawingQueue;
	},
	clearAllLayer(){
		Object.values(this.layer).forEach(layer=> layer.clear());
	},
	clearDrawing(){
		this.layer.flowerLayer.clear();
		this.layer.stemLayer.clear();
		this.layer.leafLayer.clear();
		this.layer.debugCurveLayer.clear();
	},

	//sudo only c.a.
	// growBranches2(){
	// 	let stems = [this.stem];
	// 	let burgeons = stems to burgeons;
	// 	for(let levelParam of UI.state.levelCurve){
	// 		this.scene.newLevel();
	// 		let next_level_burgeons = [];
	// 		let sign = 1;

	// 		while (!burgeons.empty()){
	// 			let bergeon = burgeons.deque();
	// 			let stem = bergeon.germinate(levelParam.length,sign=-sign);
	// 			if(stem collides with stage)
	// 				next_level_burgeons.push(bergeon);
	// 			else{
	// 				this.scene.addToScene(stem);
	// 				burgeons.enqueue_all(stem.burgeons(levelParam.amount));
	// 			}
	// 		}
	// 	}
	// }
};

function flatten(arr){
	return arr.reduce((acc, val) => acc.concat( Array.isArray(val) ? flatten(val) : val),[]);
}
 function noCollision(argument) {
 	// body...
 }