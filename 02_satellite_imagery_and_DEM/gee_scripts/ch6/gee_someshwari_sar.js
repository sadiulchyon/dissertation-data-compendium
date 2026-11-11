// ============================================================
// SAND MINING DETECTION - SOMESHWARI RIVER
// Section 1: Setup & Data Availability Check
// ============================================================

// Load AOI
var aoi = ee.FeatureCollection('projects/riversand2024/assets/someshwari/aoi');

// Visualize AOI boundary
Map.centerObject(aoi, 13);
Map.setOptions('SATELLITE');

var empty = ee.Image().byte();
var aoiOutline = empty.paint({featureCollection: aoi, color: 1, width: 2});
Map.addLayer(aoiOutline, {palette: 'red'}, 'Study Area');

// ============================================================
// CHECK SENTINEL-1 AVAILABILITY
// ============================================================

// Check DESCENDING pass
var s1_desc = ee.ImageCollection('COPERNICUS/S1_GRD_FLOAT')
  .filterBounds(aoi)
  .filterDate('2015-01-01', '2026-01-01')
  .filter(ee.Filter.eq('orbitProperties_pass', 'DESCENDING'));

// Check ASCENDING pass
var s1_asc = ee.ImageCollection('COPERNICUS/S1_GRD_FLOAT')
  .filterBounds(aoi)
  .filterDate('2015-01-01', '2026-01-01')
  .filter(ee.Filter.eq('orbitProperties_pass', 'ASCENDING'));

print('DESCENDING images available:', s1_desc.size());
print('ASCENDING images available:', s1_asc.size());

// Print first and last image dates for each pass
print('DESCENDING - first image date:', 
  ee.Image(s1_desc.sort('system:time_start').first())
    .date().format('YYYY-MM-dd'));
print('DESCENDING - last image date:', 
  ee.Image(s1_desc.sort('system:time_start', false).first())
    .date().format('YYYY-MM-dd'));

print('ASCENDING - first image date:', 
  ee.Image(s1_asc.sort('system:time_start').first())
    .date().format('YYYY-MM-dd'));
print('ASCENDING - last image date:', 
  ee.Image(s1_asc.sort('system:time_start', false).first())
    .date().format('YYYY-MM-dd'));

// Check which polarizations are available
print('DESCENDING polarizations (first image):', 
  ee.Image(s1_desc.first()).bandNames());
print('ASCENDING polarizations (first image)', 
  ee.Image(s1_asc.first()).bandNames());

// ============================================================
// SAND MINING DETECTION - SOMESHWARI RIVER
// Section 2: Build SAR Time Series
// ============================================================

var aoi = ee.FeatureCollection('projects/riversand2024/assets/someshwari/aoi');
Map.centerObject(aoi, 13);
Map.setOptions('SATELLITE');

var empty = ee.Image().byte();
var aoiOutline = empty.paint({featureCollection: aoi, color: 1, width: 2});
Map.addLayer(aoiOutline, {palette: 'red'}, 'Study Area');

// ============================================================
// STEP 2.1: Incidence angle mask
// ============================================================
function maskAngle(image) {
  var angleMask = image.select('angle');
  return image.updateMask(angleMask.gte(31).and(angleMask.lte(45)));
}

var sarColl = ee.ImageCollection('COPERNICUS/S1_GRD_FLOAT')
  .filterBounds(aoi)
  .filterDate('2015-01-01', '2026-01-01')
  .filter(ee.Filter.eq('orbitProperties_pass', 'DESCENDING'))
  .filter(ee.Filter.eq('instrumentMode', 'IW'))
  .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VH'))
  .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VV'))
  .map(maskAngle)
  .select(['VV', 'VH']);

print('Total images after angle filter:', sarColl.size());

// ============================================================
// STEP 2.3: Build daily mosaics (one image per acquisition date)
// ============================================================
function getDates(dd) {
  return ee.Date(dd).format('YYYY-MM-dd');
}

function mosaicByDate(dateStr) {
  dateStr = ee.Date(dateStr);
  var filtered = sarColl.filterDate(dateStr, dateStr.advance(1, 'day'));
  return filtered.mosaic()
    .clip(aoi)
    .set({
      'system:time_start': dateStr.millis(),
      'dateYMD': dateStr.format('YYYY-MM-dd')
    });
}

var datesMosaic = ee.List(sarColl.aggregate_array('system:time_start'))
  .map(getDates)
  .distinct();

print('Total unique acquisition dates:', datesMosaic.size());

var mosaicList = datesMosaic.map(mosaicByDate);
var mosaicColl = ee.ImageCollection(mosaicList);

print('Mosaic collection:', mosaicColl);

// ============================================================
// STEP 2.4: Visualize one early and one recent image
// ============================================================
var sarVis = {
  bands: ['VV', 'VH', 'VV'],
  min: [-18, -23, 3],
  max: [-4, -11, 15]
};

// Early image
var earlyImage = mosaicColl
  .filter(ee.Filter.eq('dateYMD', '2015-01-19'))
  .first().log10().multiply(10.0);

// Recent image - filter for a 2025 date
var recentImage = mosaicColl
  .sort('system:time_start', false)
  .first().log10().multiply(10.0);

Map.addLayer(earlyImage, sarVis, 'SAR - Early (2015-01-19)');
Map.addLayer(recentImage, sarVis, 'SAR - Most Recent');

// ============================================================
// STEP 2.5: Quick coverage check per year
// ============================================================
var yearlyCount = ee.List.sequence(2015, 2025).map(function(y) {
  var count = mosaicColl
    .filter(ee.Filter.calendarRange(y, y, 'year'))
    .size();
  return ee.List([y, count]);
});
print('Images per year:', yearlyCount);

// ============================================================
// SAND MINING DETECTION - SOMESHWARI RIVER
// Section 3: Omnibus Q Change Detection (Three Periods)
// ============================================================

// --- Paste your full Section 2 code first, then add this below ---

var omb = require('users/mortcanty/changedetection:omnibus');
var significance = 0.001;
var median = true;

// ============================================================
// STEP 3.1: Split into three periods
// ============================================================
var mosaicColl_A = mosaicColl.filterDate('2015-01-01', '2020-01-01');
var mosaicColl_B = mosaicColl.filterDate('2020-01-01', '2023-01-01');
var mosaicColl_C = mosaicColl.filterDate('2023-01-01', '2026-01-01');

print('Period A (2015-2019):', mosaicColl_A.size());
print('Period B (2020-2022):', mosaicColl_B.size());
print('Period C (2023-2025):', mosaicColl_C.size());

var mosaicList_A = mosaicColl_A.toList(mosaicColl_A.size());
var mosaicList_B = mosaicColl_B.toList(mosaicColl_B.size());
var mosaicList_C = mosaicColl_C.toList(mosaicColl_C.size());

// ============================================================
// STEP 3.2: Run Omnibus on all three periods
// ============================================================
var result_A = ee.Dictionary(omb.omnibus(mosaicList_A, significance, median));
var result_B = ee.Dictionary(omb.omnibus(mosaicList_B, significance, median));
var result_C = ee.Dictionary(omb.omnibus(mosaicList_C, significance, median));

print('Period A result:', result_A);
print('Period B result:', result_B);
print('Period C result:', result_C);

// Extract maps
var smap_A = ee.Image(result_A.get('smap')).byte();
var cmap_A = ee.Image(result_A.get('cmap')).byte();
var fmap_A = ee.Image(result_A.get('fmap')).byte();

var smap_B = ee.Image(result_B.get('smap')).byte();
var cmap_B = ee.Image(result_B.get('cmap')).byte();
var fmap_B = ee.Image(result_B.get('fmap')).byte();

var smap_C = ee.Image(result_C.get('smap')).byte();
var cmap_C = ee.Image(result_C.get('cmap')).byte();
var fmap_C = ee.Image(result_C.get('fmap')).byte();

// ============================================================
// STEP 3.3: Visualize — fixed numbers in vis params
// ============================================================
var jet = ['black', 'blue', 'cyan', 'yellow', 'red'];

var vis_A = {min: 0, max: 106, palette: jet};
var vis_B = {min: 0, max: 80,  palette: jet};
var vis_C = {min: 0, max: 80,  palette: jet};

Map.addLayer(smap_A, vis_A, 'smap_A');
Map.addLayer(cmap_A, vis_A, 'cmap_A');
Map.addLayer(fmap_A, vis_A, 'fmap_A');

Map.addLayer(smap_B, vis_B, 'smap_B');
Map.addLayer(cmap_B, vis_B, 'cmap_B');
Map.addLayer(fmap_B, vis_B, 'fmap_B');

Map.addLayer(smap_C, vis_C, 'smap_C');
Map.addLayer(cmap_C, vis_C, 'cmap_C');
Map.addLayer(fmap_C, vis_C, 'fmap_C');

// ============================================================
// STEP 3.4: Export all three smaps as assets
// ============================================================
Export.image.toAsset({
  image: smap_A,
  description: 'smap_2015_2019',
  assetId: 'projects/riversand2024/assets/someshwari/smap_2015_2019',
  region: aoi, scale: 10, maxPixels: 1e13
});

Export.image.toAsset({
  image: smap_B,
  description: 'smap_2020_2022',
  assetId: 'projects/riversand2024/assets/someshwari/smap_2020_2022',
  region: aoi, scale: 10, maxPixels: 1e13
});

Export.image.toAsset({
  image: smap_C,
  description: 'smap_2023_2025',
  assetId: 'projects/riversand2024/assets/someshwari/smap_2023_2025',
  region: aoi, scale: 10, maxPixels: 1e13
});

// Export dates CSVs
var exportDates = function(coll, desc) {
  Export.table.toDrive({
    collection: ee.FeatureCollection(
      ee.List(coll.aggregate_array('dateYMD')).map(function(d) {
        return ee.Feature(null, {date: d});
      })),
    folder: 'someshwari_sar',
    description: desc,
    fileFormat: 'CSV'
  });
};

exportDates(mosaicColl_A, 'dates_2015_2019');
exportDates(mosaicColl_B, 'dates_2020_2022');
exportDates(mosaicColl_C, 'dates_2023_2025');

// ============================================================
// SAND MINING DETECTION - SOMESHWARI RIVER
// Section 4: Filter & Final Mining Change Map
// ============================================================

var aoi = ee.FeatureCollection('projects/riversand2024/assets/someshwari/aoi');
Map.centerObject(aoi, 13);
Map.setOptions('SATELLITE');

var aoiOutline = ee.Image().byte().paint({featureCollection: aoi, color: 1, width: 2});
Map.addLayer(aoiOutline, {palette: 'red'}, 'AOI');

// ============================================================
// STEP 4.1: Load the three smap assets
// ============================================================
var smap_A = ee.Image('projects/riversand2024/assets/someshwari/smap_2015_2019');
var smap_B = ee.Image('projects/riversand2024/assets/someshwari/smap_2020_2022');
var smap_C = ee.Image('projects/riversand2024/assets/someshwari/smap_2023_2025');

// ============================================================
// STEP 4.2: Load filter layers
// ============================================================

// SRTM elevation and slope
var srtm  = ee.Image('USGS/SRTMGL1_003').clip(aoi);
var slope = ee.Terrain.slope(srtm);

// JRC water — include areas near water (buffer 100m around max extent)
var waterJRC = ee.Image('JRC/GSW1_3/GlobalSurfaceWater')
  .select('max_extent').clip(aoi);

// Convert water to vector and buffer 100m, then back to image
var waterBuffer = waterJRC.selfMask()
  .reduceToVectors({
    geometry: aoi,
    scale: 30,
    geometryType: 'polygon',
    maxPixels: 1e10
  })
  .geometry()
  .buffer(100);

var nearWater = ee.Image(1).clip(waterBuffer);

// ============================================================
// STEP 4.3: Apply filters for sand mining context
// ============================================================
function filterSandMining(smap) {
  return smap
    .updateMask(srtm.lt(150))        // low elevation only
    .updateMask(slope.lt(5))          // flat terrain (sandbars/floodplains)
    .updateMask(nearWater.eq(1))      // within 100m of water
    .selfMask();
}

var filtered_A = filterSandMining(smap_A);
var filtered_B = filterSandMining(smap_B);
var filtered_C = filterSandMining(smap_C);

// ============================================================
// STEP 4.4: Remove small patches (min 500 m²)
// ============================================================
function removeSmalPatches(image) {
  var pixelCount = image.gt(0).connectedPixelCount(200);
  var minPixelCount = ee.Image(500).divide(ee.Image.pixelArea());
  return image.updateMask(pixelCount.gte(minPixelCount));
}

var final_A = removeSmalPatches(filtered_A);
var final_B = removeSmalPatches(filtered_B);
var final_C = removeSmalPatches(filtered_C);

// ============================================================
// STEP 4.5: Visualize
// ============================================================
var jet = ['black', 'blue', 'cyan', 'yellow', 'red'];

Map.addLayer(final_A, {min: 0, max: 106, palette: jet}, 'Mining Changes 2015-2019');
Map.addLayer(final_B, {min: 0, max: 80,  palette: jet}, 'Mining Changes 2020-2022');
Map.addLayer(final_C, {min: 0, max: 80,  palette: jet}, 'Mining Changes 2023-2025');

// ============================================================
// STEP 4.6: Calculate mined area per period (sq meters and acres)
// ============================================================
function calcArea(image, label) {
  var area = image.gt(0)
    .multiply(ee.Image.pixelArea())
    .reduceRegion({
      reducer: ee.Reducer.sum(),
      geometry: aoi,
      scale: 10,
      maxPixels: 1e10
    });
  var sqm  = ee.Number(area.get('constant'));
  var acres = sqm.divide(4046.86);
  print(label + ' area (sq m):', sqm);
  print(label + ' area (acres):', acres);
}

calcArea(final_A, '2015-2019');
calcArea(final_B, '2020-2022');
calcArea(final_C, '2023-2025');

// ============================================================
// STEP 4.7: Export final filtered maps
// ============================================================
Export.image.toDrive({
  image: final_A, description: 'mining_2015_2019',
  folder: 'someshwari_sar', region: aoi, scale: 10
});
Export.image.toDrive({
  image: final_B, description: 'mining_2020_2022',
  folder: 'someshwari_sar', region: aoi, scale: 10
});
Export.image.toDrive({
  image: final_C, description: 'mining_2023_2025',
  folder: 'someshwari_sar', region: aoi, scale: 10
});