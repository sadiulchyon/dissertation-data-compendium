// =====================================================================
// Chapter 4 · Script 02 · Annual dry-season Sentinel-2 composites, 2017-2025
// Project: riversand2024 · Oct 2026
// ---------------------------------------------------------------------
// Dry season YYYY = 1 Nov (YYYY-1) to 30 Apr (YYYY). Surface reflectance
// (COPERNICUS/S2_SR_HARMONIZED) with Cloud Score+ masking (cs_cdf >= 0.60).
// SR coverage of Bangladesh begins in 2017, so season 2017 is partial.
// Inputs : ch4/corridor_mask (Script 01a)
// Outputs: ch4/s2_dry/s2_dry_YYYY (10 m, int16 x 1e4, NOBS = clear obs)
// Band names match the Landsat composites (Script 01b).
// =====================================================================
var ROOT = 'projects/riversand2024/assets/ch4/';
var START = 2017, END = 2025;
var CRS = 'EPSG:32646', SCALE = 10, CS_MIN = 0.60;
var frame = ee.Geometry.Rectangle([88.0, 24.75, 92.55, 26.70], null, false);
var corridor = ee.Image(ROOT + 'corridor_mask').gt(0);

var NAMES = ['BLUE', 'GREEN', 'RED', 'NIR', 'SWIR1', 'SWIR2'];
var S2 = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED').filterBounds(frame)
  .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 60))
  .linkCollection(ee.ImageCollection('GOOGLE/CLOUD_SCORE_PLUS/V1/S2_HARMONIZED'), ['cs_cdf'])
  .map(function (img) {
    return ee.Image(img.updateMask(img.select('cs_cdf').gte(CS_MIN))
      .select(['B2', 'B3', 'B4', 'B8', 'B11', 'B12'], NAMES).divide(10000)
      .copyProperties(img, ['system:time_start']));
  });

function indices(x) {
  var ndvi = x.normalizedDifference(['NIR', 'RED']).rename('NDVI');
  var mndwi = x.normalizedDifference(['GREEN', 'SWIR1']).rename('MNDWI');
  var ndbi = x.normalizedDifference(['SWIR1', 'NIR']).rename('NDBI');
  var bsi = x.expression('((S+R)-(N+B))/((S+R)+(N+B))',
    {S: x.select('SWIR1'), R: x.select('RED'), N: x.select('NIR'), B: x.select('BLUE')}).rename('BSI');
  return x.addBands([ndvi, mndwi, ndbi, bsi]);
}

function drySeason(y) {
  var s = ee.Date.fromYMD(y - 1, 11, 1), e = ee.Date.fromYMD(y, 5, 1);
  var col = S2.filterDate(s, e).map(indices);
  var med = col.select(NAMES.concat(['NDVI', 'MNDWI', 'NDBI', 'BSI'])).median();
  var p = col.select(['MNDWI', 'NDVI']).reduce(ee.Reducer.percentile([20, 80]));
  var nobs = col.select('GREEN').count().rename('NOBS');
  return med.multiply(10000).round().toInt16()
    .addBands(p.multiply(10000).round().toInt16())
    .addBands(nobs.min(255).toUint8())
    .updateMask(corridor)
    .set({season: y, 'system:time_start': ee.Date.fromYMD(y, 1, 15).millis(),
          n_scenes: col.size(), scale_factor: 0.0001, sensor: 'Sentinel-2 MSI',
          window: s.format('YYYY-MM-dd').cat('/').cat(e.format('YYYY-MM-dd'))});
}

for (var y = START; y <= END; y++) {
  Export.image.toAsset({image: drySeason(y), description: 'ch4_s2_dry_' + y,
    assetId: ROOT + 's2_dry/s2_dry_' + y, region: frame, scale: SCALE, crs: CRS,
    maxPixels: 1e12});
}

Map.centerObject(frame, 8);
Map.addLayer(drySeason(2020).select(['RED', 'GREEN', 'BLUE']), {min: 0, max: 2500}, 'S2 dry 2020');
