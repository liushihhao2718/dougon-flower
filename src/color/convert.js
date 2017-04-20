/*

把正倫色碼表轉換成json
 
第一步 用下面網址把csv換換成json
	http://www.convertcsv.com/csv-to-json.htm
	選csv to keyed json
 	赭筆描道     123  90  98 => { 赭筆描道: { r:123, g:90, b:98} } 

第二步 把r,g,b 三個屬性轉成HEX型式
	{ r:123, g:90, b:98} => #7b5a62

	其實寫成 rgb(123,90,98) css也是接受的
	但是根據這篇stackoverflow上的討論
	http://stackoverflow.com/questions/11620719/efficiency-of-color-selection-in-html-rgb-vs-hex-vs-name
	hex比較快一點，既然可以事先做好，那我們就把它轉換一下

*/

let csv = `name r g b
赭筆描道     123  90  98

粉暈         249 242 244

青華         194 206 242

三青         151 174 207

二青         125 132 181

大青          82  75  97

深墨壓心      75  71  65

綠華         206 213 208

三綠         157 178 176

二綠         159 185 168

大綠         124 129 104

深色草汁壓心  68  71  82

粉暈罩藤黃   252 233 201

朱華合粉罩藤黃 236 187 173

深朱壓心     175  75  64

朱華粉       233 142 115

三朱         228  99  80

二朱         208 112  83

深朱         179  81  66

深色紫礦壓心 113  61  69

綠豆褐       201 205 173

合綠         150 175 143

土黃         140 146 126

雌黃         251 204 155

黃丹(丹)     219 161 149

紫粉(銀朱)   188  99  97

合朱         194  96  84

胭脂         211  97  98

土朱         142  82  72

紫檀         120  75  71

貼真金       230 171 151

`;
let jsonArray = [
 {
   "name": "赭筆描道",
   "r": 123,
   "g": 90,
   "b": 98
 },
 {
   "name": "粉暈",
   "r": 249,
   "g": 242,
   "b": 244
 },
 {
   "name": "青華",
   "r": 194,
   "g": 206,
   "b": 242
 },
 {
   "name": "三青",
   "r": 151,
   "g": 174,
   "b": 207
 },
 {
   "name": "二青",
   "r": 125,
   "g": 132,
   "b": 181
 },
 {
   "name": "大青",
   "r": 82,
   "g": 75,
   "b": 97
 },
 {
   "name": "深墨壓心",
   "r": 75,
   "g": 71,
   "b": 65
 },
 {
   "name": "綠華",
   "r": 206,
   "g": 213,
   "b": 208
 },
 {
   "name": "三綠",
   "r": 157,
   "g": 178,
   "b": 176
 },
 {
   "name": "二綠",
   "r": 159,
   "g": 185,
   "b": 168
 },
 {
   "name": "大綠",
   "r": 124,
   "g": 129,
   "b": 104
 },
 {
   "name": "深色草汁壓心",
   "r": 68,
   "g": 71,
   "b": 82
 },
 {
   "name": "粉暈罩藤黃",
   "r": 252,
   "g": 233,
   "b": 201
 },
 {
   "name": "朱華合粉罩藤黃",
   "r": 236,
   "g": 187,
   "b": 173
 },
 {
   "name": "深朱壓心",
   "r": 175,
   "g": 75,
   "b": 64
 },
 {
   "name": "朱華粉",
   "r": 233,
   "g": 142,
   "b": 115
 },
 {
   "name": "三朱",
   "r": 228,
   "g": 99,
   "b": 80
 },
 {
   "name": "二朱",
   "r": 208,
   "g": 112,
   "b": 83
 },
 {
   "name": "深朱",
   "r": 179,
   "g": 81,
   "b": 66
 },
 {
   "name": "深色紫礦壓心",
   "r": 113,
   "g": 61,
   "b": 69
 },
 {
   "name": "綠豆褐",
   "r": 201,
   "g": 205,
   "b": 173
 },
 {
   "name": "合綠",
   "r": 150,
   "g": 175,
   "b": 143
 },
 {
   "name": "土黃",
   "r": 140,
   "g": 146,
   "b": 126
 },
 {
   "name": "雌黃",
   "r": 251,
   "g": 204,
   "b": 155
 },
 {
   "name": "黃丹(丹)",
   "r": 219,
   "g": 161,
   "b": 149
 },
 {
   "name": "紫粉(銀朱)",
   "r": 188,
   "g": 99,
   "b": 97
 },
 {
   "name": "合朱",
   "r": 194,
   "g": 96,
   "b": 84
 },
 {
   "name": "胭脂",
   "r": 211,
   "g": 97,
   "b": 98
 },
 {
   "name": "土朱",
   "r": 142,
   "g": 82,
   "b": 72
 },
 {
   "name": "紫檀",
   "r": 120,
   "g": 75,
   "b": 71
 },
 {
   "name": "貼真金",
   "r": 230,
   "g": 171,
   "b": 151
 }
];

function rgbToHex(r, g, b) {
    return "#" + ((1 << 24) + (r << 16) + (g << 8) + b).toString(16).slice(1);
}

function convertJsonToColorHex(jsonObject){
	let obj = {};
	jsonObject.forEach(c => { 
		obj[c.name] = rgbToHex(c.r, c.g, c.b);
	});
	return obj;
}

function prettyPrintJSON(jsonObject){
	const nextLine = '\n';
	const space = ' ';
	console.log( JSON.stringify(jsonObject, nextLine, space) );
}


//程式從這裡執行
(function main(){
	let hexColorJson = convertJsonToColorHex( jsonArray );
	prettyPrintJSON( hexColorJson );
})();