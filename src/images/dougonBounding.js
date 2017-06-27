export let dougonBounding = {
	令栱:require('./令栱.svg'),
	撩擔方:require('./撩擔方.svg'),
	栱眼壁:require('./new_00224.svg')
};

export let dougonBoundingNodes = {
	令栱 : [[19.205,211.473],[704.384,211.473],[760.09,185.256],[789.952,139.683],[790.342,138.473],[942.385,138.473],[942.385,225.115],[923.24,288.927],[879.538,340.298],[805.837,373.799],[703.597,395.708],[21.866,400.334]],
	撩擔方 :[[80.5,82.5],[1412.5,82.5],[1412.5,482.5],[80.5,482.5]],
	栱眼壁 : [[637.333,750.333],[727.333,687.667],[745.333,615.667],[804,575.667],[825.244,545.541],[834.577,489.541],[873.333,488.333],[884,549.667],[924.667,589.667],[970,615.667],[975.334,675],[1006,709.667],[1074,749],[1074,774.862],[637.333,774.862]]
};

export let dougonBoundingParam = {
	令栱:{
		flowerSize: 30,
		trunkHead: 5,
		trunkTail: 10
	},
	撩擔方:{
		flowerSize: 150,
		trunkHead: 5,
		trunkTail: 10
	},
	栱眼壁:{
		flowerSize: 60,
		trunkHead: 3,
		trunkTail: 5
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