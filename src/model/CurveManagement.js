import * as UI from './UIManagement';
import ColliderCollection from './Collision';
import * as Drawer from '../model/Drawer';
export default {
	panel : undefined,
	selectedCurve :[],
	scene :[],
	floralScene:[],
	layer:{},
	debug_burgeons:[],
	growBranches(){
		let leafDrawingQueue=[];
		let collisionScene= new ColliderCollection([]);
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

				if(collisionScene.test(leaf.colliders, burgeon.parent.colliders)){
					nextLevelBurgeons.push(burgeon);
				}
				else{
					leafDrawingQueue.push(leaf);
					collisionScene.add(leaf.colliders);
					this.debug_burgeons.push({leaf,level:index});
					burgeons = burgeons.concat( leaf.burgeons( levelParam.branches ) );
				}
			}
			burgeons = nextLevelBurgeons;
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
		this.debug_burgeons = [];
	},
	draw(){
		this.clearDrawing();
		
		this.floralScene.forEach(floral=>floral.draw());
		
		let stems = this.growBranches().reverse();
		stems.forEach(s => Drawer.drawLeaf(s) );
		// stems.forEach(s => Drawer.drawPolygon(s.colliders) );
		stems.forEach(s =>	{
			let bbox = makeBbox( s.colliders);
			Drawer.drawPolygon(
			[
				[bbox.x, bbox.y],
				[bbox.x+bbox.width, bbox.y],
				[bbox.x+bbox.width, bbox.y+bbox.height],
				[bbox.x, bbox.y+bbox.height],
			]);
		});
		// this.debug_burgeons.forEach(({leaf, level}) => {
		// 	Drawer.drawMagneticCurve(leaf, level);

			
		// });

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
	},
	initSvgSymbol(){
		Drawer.initSvgSymbol();
	}
};

function flatten(arr){
	return arr.reduce((acc, val) => acc.concat( Array.isArray(val) ? flatten(val) : val),[]);
}

function makeBbox(points){
	let minX = Number.MAX_VALUE;
	let minY = Number.MAX_VALUE;
	let maxX = Number.MIN_VALUE;
	let maxY = Number.MIN_VALUE;
	points.forEach(p => {
		let x = p[0];
		let y = p[1];

		if(x < minX) minX = x;
		if(y < minY) minY = y;
		if(x > maxX) maxX = x;
		if(y > maxY) maxY = y;
	});

	return {
		x: minX,
		y: minY,
		width: maxX-minX,
		height: maxY-minY
	};
}

