// =============================================================================
// Someshwari River spatial-resolution sensitivity visualization
// Sentinel-2 native 10 m versus Sentinel-2 reflectance aggregated to 30 m
// Version 2
// -----------------------------------------------------------------------------
// Produces:
//   1. A representative map showing native, 30 m, and changed classifications.
//   2. A scatter plot with an actual 1:1 reference line.
//   3. A time-series comparison.
//   4. A histogram of per-scene percentage differences.
//
// This visualizes spatial-resolution effects only. It is not a full
// Landsat versus Sentinel-2 cross-sensor calibration.
// =============================================================================


// -----------------------------------------------------------------------------
// 1. Configuration
// -----------------------------------------------------------------------------
var AOI_ASSET = 'projects/riversand2024/assets/someshwari/aoi';
var MNDWI_THRESHOLD = 0.0;
var NATIVE_SCALE_M = 10;
var TARGET_SCALE_M = 30;
var REPRESENTATIVE_DATE = '2023-03-02';


var aoi = ee.FeatureCollection(AOI_ASSET).geometry();


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
// 3. Processing functions
// -----------------------------------------------------------------------------
function maskS2Scl(img) {
  var scl = img.select('SCL');
  var clear = scl.neq(3)
    .and(scl.neq(8))
    .and(scl.neq(9))
    .and(scl.neq(10));
  return img.updateMask(clear);
}


function areaKm2(binaryImage, bandName, geometry, scale) {
  return ee.Number(binaryImage.multiply(ee.Image.pixelArea()).reduceRegion({
    reducer: ee.Reducer.sum(),
    geometry: geometry,
    scale: scale,
    maxPixels: 1e9,
    bestEffort: false,
    tileScale: 4
  }).get(bandName)).divide(1e6);
}


function classifyAtTwoResolutions(dateString) {
  var start = ee.Date(dateString);
  var image = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
    .filterBounds(aoi)
    .filterDate(start, start.advance(1, 'day'))
    .sort('system:time_start')
    .first();


  image = maskS2Scl(ee.Image(image)).clip(aoi);


  var projection10 = image.select('B3').projection();
  var projection30 = projection10.atScale(TARGET_SCALE_M);


  var mndwi10 = image.normalizedDifference(['B3', 'B11']).rename('mndwi');
  var bar10 = mndwi10.lte(MNDWI_THRESHOLD).rename('bar_10m');
  var water10 = mndwi10.gt(MNDWI_THRESHOLD).rename('water_10m');


  var b3_30 = image.select('B3')
    .reduceResolution({reducer: ee.Reducer.mean(), maxPixels: 65535})
    .reproject(projection30);
  var b11_30 = image.select('B11')
    .reduceResolution({reducer: ee.Reducer.mean(), maxPixels: 65535})
    .reproject(projection30);


  var mndwi30 = b3_30.subtract(b11_30)
    .divide(b3_30.add(b11_30))
    .rename('mndwi');
  var bar30 = mndwi30.lte(MNDWI_THRESHOLD).rename('bar_30m');
  var water30 = mndwi30.gt(MNDWI_THRESHOLD).rename('water_30m');


  return {
    image: image,
    bar10: bar10,
    water10: water10,
    bar30: bar30,
    water30: water30,
    bar10Km2: areaKm2(bar10, 'bar_10m', aoi, NATIVE_SCALE_M),
    bar30Km2: areaKm2(bar30, 'bar_30m', aoi, TARGET_SCALE_M)
  };
}


function makeFeature(dateString) {
  var result = classifyAtTwoResolutions(dateString);
  var differenceKm2 = result.bar30Km2.subtract(result.bar10Km2);
  return ee.Feature(null, {
    image_date: dateString,
    date_millis: ee.Date(dateString).millis(),
    bar_10m_km2: result.bar10Km2,
    bar_30m_km2: result.bar30Km2,
    one_to_one: result.bar10Km2,
    difference_km2: differenceKm2,
    difference_percentage: differenceKm2.divide(result.bar10Km2).multiply(100)
  });
}


var allScenes = ee.FeatureCollection(DATES.map(makeFeature));


// -----------------------------------------------------------------------------
// 4. Representative map
// -----------------------------------------------------------------------------
var representative = classifyAtTwoResolutions(REPRESENTATIVE_DATE);


Map.centerObject(aoi, 13);
Map.setOptions('SATELLITE');
Map.addLayer(
  representative.image,
  {bands: ['B4', 'B3', 'B2'], min: 0, max: 3000},
  'RGB ' + REPRESENTATIVE_DATE,
  true
);
Map.addLayer(
  representative.bar10.selfMask(),
  {palette: ['#e8b957'], opacity: 0.70},
  'Exposed bar, native 10 m',
  true
);
Map.addLayer(
  representative.bar30.selfMask(),
  {palette: ['#c72485'], opacity: 0.70},
  'Exposed bar, aggregated 30 m',
  false
);


var bar30At10 = representative.bar30
  .reproject(representative.bar10.projection())
  .rename('bar_30m_at_10m_grid');
var changedMask = representative.bar10.neq(bar30At10).selfMask();
var changeDirection = representative.bar10.subtract(bar30At10)
  .updateMask(changedMask);


Map.addLayer(
  changeDirection,
  {min: -1, max: 1, palette: ['#1a9850', '#ffffbf', '#d73027']},
  'Changed class: red lost bar, green gained bar',
  true
);


var aoiOutline = ee.Image().byte().paint({
  featureCollection: ee.FeatureCollection(AOI_ASSET),
  color: 1,
  width: 2
});
Map.addLayer(aoiOutline, {palette: ['white']}, 'AOI outline', false);


var legend = ui.Panel({
  style: {position: 'bottom-right', padding: '8px 12px', backgroundColor: 'white'}
});
legend.add(ui.Label(
  'Spatial-resolution sensitivity: ' + REPRESENTATIVE_DATE,
  {fontWeight: 'bold', fontSize: '12px'}
));
[
  ['#e8b957', 'Exposed bar at native 10 m'],
  ['#c72485', 'Exposed bar after 30 m aggregation'],
  ['#d73027', 'Bar at 10 m, water at 30 m'],
  ['#1a9850', 'Water at 10 m, bar at 30 m']
].forEach(function(item) {
  legend.add(ui.Panel([
    ui.Label('', {
      backgroundColor: item[0],
      padding: '8px',
      margin: '2px 6px 2px 0'
    }),
    ui.Label(item[1], {fontSize: '11px', margin: '2px 0'})
  ], ui.Panel.Layout.flow('horizontal')));
});
Map.add(legend);


// -----------------------------------------------------------------------------
// 5. Scatter plot with a true 1:1 line
// -----------------------------------------------------------------------------
var scatterData = allScenes.sort('bar_10m_km2');
var scatter = ui.Chart.feature.byFeature({
  features: scatterData,
  xProperty: 'bar_10m_km2',
  yProperties: ['bar_30m_km2', 'one_to_one']
}).setChartType('LineChart').setOptions({
  title: 'Exposed-bar area: aggregated 30 m versus native 10 m',
  hAxis: {title: 'Native 10 m bar area (km2)'},
  vAxis: {title: 'Aggregated 30 m bar area (km2)'},
  series: {
    0: {lineWidth: 0, pointSize: 5},
    1: {lineWidth: 1.5, pointSize: 0, lineDashStyle: [4, 4]}
  },
  legend: {position: 'top'}
});
print(scatter);


print('Pearson correlation:', allScenes.reduceColumns({
  reducer: ee.Reducer.pearsonsCorrelation(),
  selectors: ['bar_10m_km2', 'bar_30m_km2']
}));


// -----------------------------------------------------------------------------
// 6. Time series
// -----------------------------------------------------------------------------
var timeSeries = ui.Chart.feature.byFeature({
  features: allScenes.sort('date_millis'),
  xProperty: 'date_millis',
  yProperties: ['bar_10m_km2', 'bar_30m_km2']
}).setChartType('LineChart').setOptions({
  title: 'Exposed-bar area at 10 m and 30 m, 2017-2025',
  hAxis: {title: 'Date', format: 'yyyy', gridlines: {count: 9}},
  vAxis: {title: 'Exposed-bar area (km2)'},
  series: {
    0: {pointSize: 4, lineWidth: 1.2},
    1: {pointSize: 4, lineWidth: 1.2}
  },
  legend: {position: 'top'}
});
print(timeSeries);


// -----------------------------------------------------------------------------
// 7. Histogram and numerical summary
// -----------------------------------------------------------------------------
var histogram = ui.Chart.feature.histogram({
  features: allScenes,
  property: 'difference_percentage',
  minBucketWidth: 0.5
}).setChartType('ColumnChart').setOptions({
  title: 'Per-scene difference: (30 m - 10 m) / 10 m',
  hAxis: {title: 'Difference (%)'},
  vAxis: {title: 'Number of scenes'},
  legend: {position: 'none'}
});
print(histogram);


print('Mean difference, percent:',
      allScenes.aggregate_mean('difference_percentage'));
print('Median difference, percent:',
      allScenes.reduceColumns(ee.Reducer.median(), ['difference_percentage']));
print('Representative-scene areas, km2:', ee.Dictionary({
  native_10m: representative.bar10Km2,
  aggregated_30m: representative.bar30Km2,
  difference: representative.bar30Km2.subtract(representative.bar10Km2)
}));


// =============================================================================
// End
// =============================================================================