/*
	e
	basePath : Spline
	trunk { beginWidth, endWidth } : {number,number}
*/
importScripts('./bezier.js');
onmessage = function (e) {
	let beziers = e.data.basePath.map(b =>	
		new Bezier(
			b[0][0], b[0][1],
			b[1][0], b[1][1],
			b[2][0], b[2][1],
			b[3][0], b[3][1]
		)
	);
	const length = beziers.map(b=>b.length()).reduce( (a,b)=>a+b , 0);

	const pathes = e.data.trunk.map(([beginWidth, endWidth])=>{
		const fullPath = beziersOutline(beziers,length, beginWidth, endWidth);
		const pathString = svgPathString(fullPath);
		return pathString;
	});
	
	postMessage(pathes);
};

function svgPathString(fittedLineData) {
	var str = '';
	//bezier : [ [c0], [c1], [c2], [c3] ]
	fittedLineData.forEach(function (bezier, i) {
		if (i == 0) {
			str += 'M ' + bezier.points[0].x + ' ' + bezier.points[0].y;
		}

		str += 'C ' + bezier.points[1].x + ' ' + bezier.points[1].y + ', ' +
		bezier.points[2].x + ' ' + bezier.points[2].y + ', ' +
		bezier.points[3].x + ' ' + bezier.points[3].y + ' ';	
					
	});

	return str;
}

function beziersOutline(beziers,length, beginWidth, endWidth){
	let l = 0;
	let w = (endWidth-beginWidth)/length;
	let outline = beziers.map(b => {
		let d1 = beginWidth + l * w;
		l += b.length();
		let d2 = beginWidth + l * w;
		if(d1 === 0) d1 = 1;
		return b.outline(d1, d1, d2, d2);
	});

	let f = [];
	let r = [];
	let head=outline[0].curves[0];
	let tail;
	const lastPath = outline[outline.length-1];
	tail = lastPath.curves[lastPath.curves.length/2];

	outline.forEach(b => {
		const length = Math.floor(b.curves.length/2)- 1;
		f = f.concat(b.curves.slice(1, length+1));
		r = r.concat(b.curves.slice(length+2).reverse());
	});

	return [].concat([head],f,[tail], r.reverse());
}