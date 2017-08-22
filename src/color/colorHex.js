let colorMap = require('./colorHex.json');
window['colorMap'] = colorMap;

export function changeColorMap(newColorMap){
	window['colorMap'] = newColorMap;
}