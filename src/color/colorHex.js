let colorMap = require('./colorHex.json');
export default colorMap;

export function changeColorMap(newColorMap){
	colorMap = newColorMap;
}