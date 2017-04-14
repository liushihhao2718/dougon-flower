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

				let leaf = burgeon.germinate(levelParam.length,levelParam.alpha, sign=-1*sign);

				if( Collision.testCollision(leaf.colliders, this.leafCollisionScene, burgeon.parent.colliders )){
					nextLevelBurgeons.push(burgeon);
				}
				else{
					this.leafDrawingQueue.push(leaf);
					this.leafCollisionScene.push(leaf.colliders);
					burgeons = burgeons.concat( leaf.burgeons( levelParam.branches ) );
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
	}
};

function flatten(arr){
	return arr.reduce((acc, val) => acc.concat( Array.isArray(val) ? flatten(val) : val),[]);
}