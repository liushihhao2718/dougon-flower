export let dougonBounding = {
	令栱:require('./令栱.svg'),
	撩擔方:require('./撩擔方.svg')
};

export let dougonBoundingNodes = {
	令栱 : [[19.205,211.473],[704.384,211.473],[760.09,185.256],[789.952,139.683],[790.342,138.473],[942.385,138.473],[942.385,225.115],[923.24,288.927],[879.538,340.298],[805.837,373.799],[703.597,395.708],[21.866,400.334]],
	撩擔方 :[[80.5,82.5],[1412.5,82.5],[1412.5,482.5],[80.5,482.5]]
};

export let dougonBoundingParam = {
	令栱:{
		flowerSize: 30,
		trunkHead: 5,
		trunkTail: 10
	},
	撩擔方:{
		flowerSize: 90,
		trunkHead: 5,
		trunkTail: 30
	}
};

/*
	這兩個看似無用的function很有用，可以將以拉的svg description轉成陣列貼在 dougonBoundingNodes
*/

function toArray(pathDescription) {
	return pathDescription.split(/[^0-9,.]/).filter(x=>x.length !=0).map(x=>x.split(/\,/));
}

function printArray(array) {
	return '['+array.map(p=>'['+p.join(',')+']').join(',')+']';
}