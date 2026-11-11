// =====================================================================
// Chapter 4 · Script 01c · Annual greenness (Landsat NDVI 90th percentile)
// Project: riversand2024 · Oct 2026
// ---------------------------------------------------------------------
// For hydrological year YYYY (1 May YYYY-1 to 30 Apr YYYY), the 90th
// percentile of clear-sky NDVI. Cropland that is bare in the dry season is
// green in the monsoon or in the boro season; excavated ground, spoil,
// stockpiles and crusher yards stay below about 0.35 all year. This band
// separates seasonal fallow from extraction surfaces in Scripts 03-04.
// Output: ch4/ls_ndvi_p90/ls_ndvi_p90_YYYY (int16 x 1e4), corridor only.
// =====================================================================
var ROOT = 'projects/riversand2024/assets/ch4/';
var START = 1988, END = 2025;
var frame = ee.Geometry.Rectangle([88.0, 24.75, 92.55, 26.70], null, false);
var corridor = ee.Image(ROOT + 'corridor_mask').gt(0);

function clear(img) {
  var qa = img.select('QA_PIXEL');
  return img.updateMask(qa.bitwiseAnd(parseInt('111110', 2)).eq(0));
}
function ndviTM(img) {
  return clear(img).select(['SR_B4', 'SR_B3']).multiply(0.0000275).add(-0.2)
    .normalizedDifference(['SR_B4', 'SR_B3']).rename('NDVI').toFloat();
}
function ndviOLI(img) {
  return clear(img).select(['SR_B5', 'SR_B4']).multiply(0.0000275).add(-0.2)
    .normalizedDifference(['SR_B5', 'SR_B4']).rename('NDVI').toFloat();
}
var NDVI = ee.ImageCollection('LANDSAT/LT05/C02/T1_L2').map(ndviTM)
  .merge(ee.ImageCollection('LANDSAT/LE07/C02/T1_L2').filterDate('1999-01-01', '2021-01-01').map(ndviTM))
  .merge(ee.ImageCollection('LANDSAT/LC08/C02/T1_L2').map(ndviOLI))
  .merge(ee.ImageCollection('LANDSAT/LC09/C02/T1_L2').map(ndviOLI))
  .filterBounds(frame);

for (var y = START; y <= END; y++) {
  var col = NDVI.filterDate(ee.Date.fromYMD(y - 1, 5, 1), ee.Date.fromYMD(y, 5, 1));
  var img = col.reduce(ee.Reducer.percentile([90])).rename('NDVI_P90')
    .multiply(10000).round().toInt16().updateMask(corridor)
    .set({season: y, 'system:time_start': ee.Date.fromYMD(y, 1, 15).millis(),
          n_scenes: col.size(), scale_factor: 0.0001});
  Export.image.toAsset({image: img, description: 'ch4_ls_ndvi_p90_' + y,
    assetId: ROOT + 'ls_ndvi_p90/ls_ndvi_p90_' + y, region: frame, scale: 30,
    crs: 'EPSG:32646', maxPixels: 1e11});
}
