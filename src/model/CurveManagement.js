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
		const level = 0;
		const max = UI.state.levelCurve.length - 1;

		while(burgeons.length > 0){
			let burgeon = burgeons[0];


			let leaf = burgeon.germinate(sign);

			if( Collision.testCollision(leaf.boundingBox, this.leafCollisionScene, leaf.parent ))


			
			


			burgeons.shift();
			if(leaf === undefined)continue;

			let segment = UI.state.levelCurve[leaf.level];
			burgeons = burgeons.concat( leaf.burgeons() );
			this.leafDrawingQueue.push( leaf );

			sign *= -1;
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
	growBranches2(){
		const amount = UI.state.levelCurve[0].branches;
		let stems = [this.stem];
		let burgeons;
		for(let levelParam of UI.state.levelCurve){
			this.scene.newLevel();
			let new_burgeons = [];
			let new_stems = [];
			let sign = 1;
			burgeons = burgeons.concat(stems.burgeons(levelParam.amount).flatten())
			for(let bergeon of burgeons){
				let stem = bergeon.germinate(levelParam.length,sign);
				if(stem collides with stage)
					new_burgeons.push(bergeon);
				else{
					AddToScene(stem);
					new_stems.push(stem);
				}
			}
			stems = new_stems;
		}
	}
};

function flatten(arr){
	return arr.reduce((acc, val) => acc.concat( Array.isArray(val) ? flatten(val) : val),[]);
}
 function noCollision(argument) {
 	// body...
 }