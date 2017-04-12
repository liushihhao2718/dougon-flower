export default {
	selectedCurve :[],
	scene :[],
	leafCollisionScene:[],
	floralScene:[],
	layer:{},
	leafDrawingQueue:[],
	growBranches(){
		let groupedBurgeons = this.floralScene.map( f=>f.burgeons() );
		let burgeons = flatten( groupedBurgeons );
		let sign = 1;
		while(burgeons.length > 0){
			let burgeon = burgeons[0];
			const level = 0;
			let stem = burgeon.germinate(level, sign);

			burgeons.shift();
			burgeons = burgeons.concat( stem.burgeons() );
			this.leafDrawingQueue.push(stem.burgeons() );

			sign *= -1;
		}
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
	arr.reduce((acc, val) => acc.concat( Array.isArray(val) ? flatten(val) : val),[]);
}
