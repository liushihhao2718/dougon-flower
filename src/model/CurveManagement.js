export default {
	selectedCurve :[],
	scene :[],
	leafCollisionScene:[],
	floralScene:[],
	layer:{},
	leafDrawingQueue:[],
	growBranches(){
		// this.floralScene.map()
	},
	clearLayer(){
		Object.values(this.layer).forEach(layer=> layer.clear());
	}
};