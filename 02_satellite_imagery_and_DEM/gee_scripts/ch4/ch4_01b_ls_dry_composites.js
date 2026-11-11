// =====================================================================
// Chapter 4 · Script 01b · Annual dry-season Landsat composites, 1988-2025
// Project: riversand2024 · Oct 2026
// ---------------------------------------------------------------------
// Dry season YYYY = 1 Nov (YYYY-1) to 30 Apr (YYYY), the low-flow window
// used in Chapters 5 and 6, when bars, pits and stockpiles are exposed.
// Inputs : projects/riversand2024/assets/ch4/corridor_mask (Script 01a)
// Outputs: ch4/ls_dry/ls_dry_YYYY (int16, x 1e4; NOBS = clear observations)
//          ch4/ls_dry_availability (table, km2 by season and sector)
// Landsat TM/ETM+ reflectance is transformed to OLI equivalents with the
// OLS coefficients of Roy et al. (2016). ETM+ is used through 2020 only.
// =====================================================================
var ROOT = 'projects/riversand2024/assets/ch4/';
var START = 1988, END = 2025;
var CRS = 'EPSG:32646', SCALE = 30;
var frame = ee.Geometry.Rectangle([88.0, 24.75, 92.55, 26.70], null, false);
var breaks = [88.0, 89.3, 89.92, 90.3, 90.95, 91.6, 92.55];
var corridorMask = ee.Image(ROOT + 'corridor_mask');
var corridor = corridorMask.gt(0);

// ---------------------------------------------------------------------
// 2. Landsat Collection 2 Level 2 surface reflectance, harmonized
// ---------------------------------------------------------------------
function prep(img) {
  var qa = img.select('QA_PIXEL');
  // bits 1 dilated cloud, 2 cirrus, 3 cloud, 4 shadow, 5 snow
  var clear = qa.bitwiseAnd(parseInt('111110', 2)).eq(0)
    .and(img.select('QA_RADSAT').eq(0));
  return ee.Image(img.updateMask(clear)
    .select('SR_B.').multiply(0.0000275).add(-0.2)
    .copyProperties(img, ['system:time_start', 'SPACECRAFT_ID']));
}
var NAMES = ['BLUE', 'GREEN', 'RED', 'NIR', 'SWIR1', 'SWIR2'];
function tm(img) {   // TM / ETM+ -> OLI equivalent (Roy et al., 2016, OLS)
  var x = prep(img).select(['SR_B1', 'SR_B2', 'SR_B3', 'SR_B4', 'SR_B5', 'SR_B7'], NAMES);
  var slope = ee.Image.constant([0.8474, 0.8483, 0.9047, 0.8462, 0.8937, 0.9071]);
  var icpt = ee.Image.constant([0.0003, 0.0088, 0.0061, 0.0412, 0.0254, 0.0172]);
  return ee.Image(x.multiply(slope).add(icpt).rename(NAMES).toFloat()
    .copyProperties(img, ['system:time_start', 'SPACECRAFT_ID']));
}
function oli(img) {
  return ee.Image(prep(img).select(['SR_B2', 'SR_B3', 'SR_B4', 'SR_B5', 'SR_B6', 'SR_B7'], NAMES).toFloat()
    .copyProperties(img, ['system:time_start', 'SPACECRAFT_ID']));
}
function indices(x) {
  var ndvi = x.normalizedDifference(['NIR', 'RED']).rename('NDVI');
  var mndwi = x.normalizedDifference(['GREEN', 'SWIR1']).rename('MNDWI');
  var ndbi = x.normalizedDifference(['SWIR1', 'NIR']).rename('NDBI');
  var bsi = x.expression('((S+R)-(N+B))/((S+R)+(N+B))',
    {S: x.select('SWIR1'), R: x.select('RED'), N: x.select('NIR'), B: x.select('BLUE')}).rename('BSI');
  return x.addBands([ndvi, mndwi, ndbi, bsi]);
}

var L5 = ee.ImageCollection('LANDSAT/LT05/C02/T1_L2').map(tm);
var L7 = ee.ImageCollection('LANDSAT/LE07/C02/T1_L2')
  .filterDate('1999-01-01', '2021-01-01').map(tm);   // excludes orbit-drift years
var L8 = ee.ImageCollection('LANDSAT/LC08/C02/T1_L2').map(oli);
var L9 = ee.ImageCollection('LANDSAT/LC09/C02/T1_L2').map(oli);
var LS = L5.merge(L7).merge(L8).merge(L9).filterBounds(frame);

function drySeason(y) {
  var s = ee.Date.fromYMD(y - 1, 11, 1), e = ee.Date.fromYMD(y, 5, 1);
  var col = LS.filterDate(s, e).map(indices);
  var med = col.select(NAMES.concat(['NDVI', 'MNDWI', 'NDBI', 'BSI'])).median();
  var p = col.select(['MNDWI', 'NDVI']).reduce(ee.Reducer.percentile([20, 80]));
  var nobs = col.select('GREEN').count().rename('NOBS');
  var img = med.multiply(10000).round().toInt16()
    .addBands(p.multiply(10000).round().toInt16())
    .addBands(nobs.toUint8())
    .updateMask(corridor)
    .set({season: y, 'system:time_start': ee.Date.fromYMD(y, 1, 15).millis(),
          n_scenes: col.size(), scale_factor: 0.0001,
          window: s.format('YYYY-MM-dd').cat('/').cat(e.format('YYYY-MM-dd'))});
  return img;
}

// ---------------------------------------------------------------------
// 3. Availability: corridor area with >= 1 and >= 3 clear observations
// ---------------------------------------------------------------------
function availability(y) {
  var s = ee.Date.fromYMD(y - 1, 11, 1), e = ee.Date.fromYMD(y, 5, 1);
  var col = LS.filterDate(s, e);
  var nobs = col.select('GREEN').count().unmask(0);
  var a = ee.Image.pixelArea().divide(1e6);
  var img = ee.Image.cat([a.rename('area_km2'), a.multiply(nobs.gte(1)).rename('ge1_km2'),
    a.multiply(nobs.gte(3)).rename('ge3_km2'), corridorMask]).updateMask(corridor);
  var d = img.reduceRegion({reducer: ee.Reducer.sum().repeat(3).group(3, 'sector'),
    geometry: frame, scale: 90, crs: CRS, maxPixels: 1e10, tileScale: 4});
  var groups = ee.List(d.get('groups'));
  return ee.FeatureCollection(groups.map(function (g) {
    g = ee.Dictionary(g); var v = ee.List(g.get('sum'));
    return ee.Feature(ee.Geometry.Point([90.5, 25.5]), {season: y, sector: g.get('sector'),
      area_km2: v.get(0), ge1_km2: v.get(1), ge3_km2: v.get(2),
      n_scenes: col.size(),
      n_L5: col.filter(ee.Filter.eq('SPACECRAFT_ID', 'LANDSAT_5')).size(),
      n_L7: col.filter(ee.Filter.eq('SPACECRAFT_ID', 'LANDSAT_7')).size(),
      n_L8L9: col.filter(ee.Filter.inList('SPACECRAFT_ID', ['LANDSAT_8', 'LANDSAT_9'])).size()});
  }));
}

// ---------------------------------------------------------------------
// 4. Exports (start from the Tasks tab, or by the run harness)
// ---------------------------------------------------------------------
var years = [];
for (var yy = START; yy <= END; yy++) years.push(yy);

Export.table.toAsset({collection: ee.FeatureCollection(years.map(availability)).flatten(),
  description: 'ch4_ls_dry_availability', assetId: ROOT + 'ls_dry_availability'});

years.forEach(function (y) {
  Export.image.toAsset({image: drySeason(y), description: 'ch4_ls_dry_' + y,
    assetId: ROOT + 'ls_dry/ls_dry_' + y, region: frame, scale: SCALE, crs: CRS,
    maxPixels: 1e11});
});

// Display
Map.centerObject(frame, 8);
Map.addLayer(drySeason(2010).select(['RED', 'GREEN', 'BLUE']), {min: 0, max: 2500}, 'dry 2010');
