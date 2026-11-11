// =============================================================================
// Someshwari River spatial-resolution sensitivity test
// Sentinel-2 native 10 m versus Sentinel-2 reflectance aggregated to 30 m
// Same 65 mechanized-era dry-season dates used in Figure 4
// Version 2
// -----------------------------------------------------------------------------
//
// PURPOSE
//   Quantify how the 10 m versus 30 m spatial-resolution contrast affects the
//   mapped exposed-bar area. This is a spatial-resolution sensitivity test,
//   not a complete Landsat versus Sentinel-2 cross-sensor calibration.
//
// METHOD
//   1. Select the same pre-screened Sentinel-2 SR scene date used in Figure 4.
//   2. Apply the SCL cloud/shadow mask.
//   3. Compute MNDWI at native resolution; bar = MNDWI <= 0.
//   4. Aggregate B3 and B11 reflectance to 30 m by pixel mean, then recompute
//      MNDWI with the same threshold.
//   5. Export per-scene 10 m and 30 m bar area and their differences.
//
// IMPORTANT CAVEAT
//   This test isolates spatial-resolution and mixed-pixel effects only. It does
//   not quantify differences in spectral response, radiometric calibration,
//   atmospheric correction, or acquisition geometry between sensor systems.
// =============================================================================


// -----------------------------------------------------------------------------
// 1. Configuration
// -----------------------------------------------------------------------------
var AOI_ASSET = 'projects/riversand2024/assets/someshwari/aoi';
var MNDWI_THRESHOLD = 0.0;
var NATIVE_SCALE_M = 10;
var TARGET_SCALE_M = 30;


var aoi = ee.FeatureCollection(AOI_ASSET).geometry();
var aoiAreaKm2 = ee.Number(aoi.area()).divide(1e6);


// -----------------------------------------------------------------------------
// 2. Mechanized-era dates used in Figure 4
// -----------------------------------------------------------------------------
var DATES = [
  '2017-02-11','2018-03-03','2018-03-08','2019-01-02','2019-01-07',
  '2019-01-12','2019-01-17','2019-01-22','2019-02-06','2019-03-03',
  '2019-03-08','2020-01-17','2020-01-22','2020-02-06','2020-02-11',
  '2020-03-02','2020-03-17','2020-03-27','2020-04-01','2020-04-06',
  '2021-01-01','2021-01-06','2021-02-05','2021-02-10','2021-02-15',
  '2021-03-17','2021-03-22','2021-03-27','2021-04-06','2022-01-16',
  '2022-01-31','2022-02-15','2022-03-02','2022-03-07','2022-03-17',
  '2023-01-11','2023-01-21','2023-01-26','2023-02-05','2023-02-10',
  '2023-02-15','2023-02-20','2023-02-25','2023-03-02','2023-03-07',
  '2023-03-12','2023-04-11','2024-01-26','2024-02-10','2024-03-01',
  '2024-03-06','2024-03-11','2024-03-16','2024-04-15','2025-01-10',
  '2025-01-20','2025-01-25','2025-02-04','2025-02-09','2025-02-14',
  '2025-02-19','2025-02-24','2025-03-01','2025-03-06','2025-03-26'
];


// -----------------------------------------------------------------------------
// 3. Helper functions
// -----------------------------------------------------------------------------
function maskS2Scl(img) {
  var scl = img.select('SCL');
  var clear = scl.neq(3)       // cloud shadow
    .and(scl.neq(8))           // medium-probability cloud
    .and(scl.neq(9))           // high-probability cloud
    .and(scl.neq(10));         // cirrus
  return img.updateMask(clear);
}


function areaKm2(binaryImage, bandName, geometry, scale) {
  var areaM2 = binaryImage.multiply(ee.Image.pixelArea()).reduceRegion({
    reducer: ee.Reducer.sum(),
    geometry: geometry,
    scale: scale,
    maxPixels: 1e9,
    bestEffort: false,
    tileScale: 4
  }).get(bandName);
  return ee.Number(areaM2).divide(1e6);
}


function processScene(dateString) {
  var start = ee.Date(dateString);
  var end = start.advance(1, 'day');


  // The Figure 4 workflow used the first chronologically sorted image for each
  // pre-screened date. scene_count and system_index are exported for audit.
  var collection = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
    .filterBounds(aoi)
    .filterDate(start, end)
    .sort('system:time_start');


  var sceneCount = collection.size();
  var image = ee.Image(collection.first());
  var masked = maskS2Scl(image).clip(aoi);


  var projection10 = masked.select('B3').projection();
  var projection30 = projection10.atScale(TARGET_SCALE_M);


  // Native 10 m classification.
  var mndwi10 = masked.normalizedDifference(['B3', 'B11']).rename('mndwi');
  var water10 = mndwi10.gt(MNDWI_THRESHOLD).rename('water');
  var bar10 = mndwi10.lte(MNDWI_THRESHOLD).rename('bar');


  // Aggregate reflectance first, then calculate MNDWI at 30 m.
  var b3_30 = masked.select('B3')
    .reduceResolution({reducer: ee.Reducer.mean(), maxPixels: 65535})
    .reproject(projection30);
  var b11_30 = masked.select('B11')
    .reduceResolution({reducer: ee.Reducer.mean(), maxPixels: 65535})
    .reproject(projection30);


  var mndwi30 = b3_30.subtract(b11_30)
    .divide(b3_30.add(b11_30))
    .rename('mndwi');
  var water30 = mndwi30.gt(MNDWI_THRESHOLD).rename('water');
  var bar30 = mndwi30.lte(MNDWI_THRESHOLD).rename('bar');


  var bar10Km2 = areaKm2(bar10, 'bar', aoi, NATIVE_SCALE_M);
  var water10Km2 = areaKm2(water10, 'water', aoi, NATIVE_SCALE_M);
  var bar30Km2 = areaKm2(bar30, 'bar', aoi, TARGET_SCALE_M);
  var water30Km2 = areaKm2(water30, 'water', aoi, TARGET_SCALE_M);


  var differenceKm2 = bar30Km2.subtract(bar10Km2);
  var differencePercent = differenceKm2.divide(bar10Km2).multiply(100);


  var coverageFraction = masked.select('B3').mask().reduceRegion({
    reducer: ee.Reducer.mean(),
    geometry: aoi,
    scale: NATIVE_SCALE_M,
    maxPixels: 1e9,
    bestEffort: false,
    tileScale: 4
  }).get('B3');


  return ee.Feature(null, {
    image_date: dateString,
    scene_count: sceneCount,
    system_index: image.get('system:index'),
    system_time_start: image.get('system:time_start'),
    cloudy_pixel_percentage: image.get('CLOUDY_PIXEL_PERCENTAGE'),
    aoi_coverage_percentage: ee.Number(coverageFraction).multiply(100),
    bar_km2_10m: bar10Km2,
    water_km2_10m: water10Km2,
    bar_km2_30m: bar30Km2,
    water_km2_30m: water30Km2,
    difference_km2_30m_minus_10m: differenceKm2,
    difference_percentage: differencePercent,
    mndwi_threshold: MNDWI_THRESHOLD,
    aoi_area_km2: aoiAreaKm2,
    aoi_asset: AOI_ASSET,
    run_date: ee.Date(Date.now()).format('YYYY-MM-dd')
  });
}


// -----------------------------------------------------------------------------
// 4. Run, inspect, and export
// -----------------------------------------------------------------------------
var results = ee.FeatureCollection(DATES.map(processScene));


print('Expected rows:', DATES.length);
print('Rows produced:', results.size());
print('Rows requiring scene-count review:',
      results.filter(ee.Filter.neq('scene_count', 1)));
print('Mean difference, percent:', results.aggregate_mean('difference_percentage'));
print('Median difference, percent:',
      results.reduceColumns(ee.Reducer.median(), ['difference_percentage']));
print('First output row:', results.first());


Export.table.toDrive({
  collection: results,
  description: 'Someshwari_S2_10m_vs_30m_resolution_sensitivity',
  folder: 'GEE_Exports',
  fileNamePrefix: 'Someshwari_S2_10m_vs_30m_resolution_sensitivity',
  fileFormat: 'CSV',
  selectors: [
    'image_date', 'scene_count', 'system_index', 'system_time_start',
    'cloudy_pixel_percentage', 'aoi_coverage_percentage',
    'bar_km2_10m', 'water_km2_10m',
    'bar_km2_30m', 'water_km2_30m',
    'difference_km2_30m_minus_10m', 'difference_percentage',
    'mndwi_threshold', 'aoi_area_km2', 'aoi_asset', 'run_date'
  ]
});


// =============================================================================
// End
// =============================================================================