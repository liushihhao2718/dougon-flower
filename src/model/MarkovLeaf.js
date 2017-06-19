const markov = [
	[4,4,6,14,17],
	[0,1,1,1,1],
	[0,0,2,2,2],
	[5,5,6,10,11],
	[0,0,3,4,5]
];

/*
leafType:
@property name symbol
@property order number
*/
export default function nextType(currentType) {

return range(0, 4);



	const row = markov[currentType];
	const max = row[row.length-1];
	const random = range(0, max);



	for (let i = 0; i < row.length; i++) {
		if( random < row[i] ) return i;
	}
	return row.length - 1;
}

function range(min, max){
	return Math.floor(Math.random() * ((max - min)+1) + min);
}