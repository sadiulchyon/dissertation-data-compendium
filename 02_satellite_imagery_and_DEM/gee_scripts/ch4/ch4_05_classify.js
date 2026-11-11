// =====================================================================
// Chapter 4 · Script 05 · Random forest classification of every dry season
// Project: riversand2024 · Oct 2026
// ---------------------------------------------------------------------
// One random forest (300 trees, sqrt of predictors per split) is trained on
// ch4/training_ls (Script 04), pooled across seasons and sensors, with 30%
// of samples held out for a first internal check. The model is applied to
// each dry season 1988-2025 (composites from ch4/ls_dry; annual NDVI p90
// computed on the fly). Output per season: band class (1 water, 2 natural
// sediment, 3 excavated or disturbed ground, 4 vegetation and cropland,
// 5 built and other) and pmax (highest class probability x 100).
// Outputs: ch4/ls_class/ls_class_YYYY; ch4/rf_holdout (confusion table);
//          ch4/rf_importance (variable importance).
// The holdout check is not the accuracy assessment: training and holdout
// points share sources and seasons. Accuracy is assessed with an
// independent stratified sample (Script 06; Olofsson et al., 2014).
// =====================================================================
var ROOT = 'projects/riversand2024/assets/ch4/';
var CRS = 'EPSG:32646', SCALE = 30;
var frame = ee.Geometry.Rectangle([88.0, 24.75, 92.55, 26.70], null, false);
var corridorMask = ee.Image(ROOT + 'corridor_mask');
var corridor = corridorMask.gt(0);
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


function ndviP90(y) {   // Script 01c: hydrological year May (y-1) to April (y)
  var col = LS.filterDate(ee.Date.fromYMD(y - 1, 5, 1), ee.Date.fromYMD(y, 5, 1))
    .map(function (i) { return i.normalizedDifference(['NIR', 'RED']).rename('NDVI'); });
  return col.reduce(ee.Reducer.percentile([90])).rename('NDVI_P90').multiply(10000).round().toInt16();
}

var pre2000 = ee.ImageCollection('JRC/GSW1_4/YearlyHistory')
  .filter(ee.Filter.calendarRange(1984, 1999, 'year'))
  .map(function (i) { return i.gte(2); }).max().unmask(0).rename('PRE2000_WATER');
var distPre = pre2000.selfMask().fastDistanceTransform(256).sqrt().multiply(30)
  .rename('DIST_PRE2000').unmask(7680);
function predictors(y) {
  var dry = ee.Image(ROOT + 'ls_dry/ls_dry_' + y);
  var tex = dry.select('NIR').reduceNeighborhood(ee.Reducer.stdDev(), ee.Kernel.square(1)).rename('NIR_SD3');
  return dry.addBands(ndviP90(y)).addBands(tex).addBands(pre2000).addBands(distPre)
    .addBands(corridorMask.rename('sector'));
}
var BANDS = ['BLUE', 'GREEN', 'RED', 'NIR', 'SWIR1', 'SWIR2', 'NDVI', 'MNDWI', 'NDBI', 'BSI',
  'MNDWI_p20', 'MNDWI_p80', 'NDVI_p20', 'NDVI_p80', 'NDVI_P90', 'NIR_SD3', 'PRE2000_WATER',
  'DIST_PRE2000', 'sector'];

var samples = ee.FeatureCollection(ROOT + 'training_ls')
  .filter(ee.Filter.inList('class', [1, 2, 3, 4, 5]))
  .filter(ee.Filter.notNull(BANDS))
  .randomColumn('split', 42);
var trainSet = samples.filter(ee.Filter.lt('split', 0.7));
var testSet = samples.filter(ee.Filter.gte('split', 0.7));
var rf = ee.Classifier.smileRandomForest({numberOfTrees: 300, minLeafPopulation: 1, bagFraction: 0.632, seed: 42})
  .train({features: trainSet, classProperty: 'class', inputProperties: BANDS});

var cm = testSet.classify(rf).errorMatrix('class', 'classification');
Export.table.toAsset({collection: ee.FeatureCollection([ee.Feature(ee.Geometry.Point([90.5, 25.5]), {
    matrix: cm.array().toList(), order: cm.order(), overall: cm.accuracy(), kappa: cm.kappa(),
    users: cm.consumersAccuracy().toList(), producers: cm.producersAccuracy().toList(),
    n_train: trainSet.size(), n_test: testSet.size()})]),
  description: 'ch4_rf_holdout', assetId: ROOT + 'rf_holdout'});
var imp = ee.Dictionary(ee.Dictionary(rf.explain()).get('importance'));
Export.table.toAsset({collection: ee.FeatureCollection(imp.keys().map(function (k) {
    return ee.Feature(ee.Geometry.Point([90.5, 25.5]), {band: k, importance: imp.get(k)}); })),
  description: 'ch4_rf_importance', assetId: ROOT + 'rf_importance'});

var rfp = rf.setOutputMode('MULTIPROBABILITY');
function classifySeason(y) {
  var x = predictors(y).select(BANDS);
  var probs = x.classify(rfp);
  var cls = probs.arrayArgmax().arrayGet(0).add(1).rename('class').toUint8();
  var pmax = probs.arrayReduce(ee.Reducer.max(), [0]).arrayGet(0).multiply(100).round().rename('pmax').toUint8();
  return cls.addBands(pmax).updateMask(corridor)
    .set({season: y, 'system:time_start': ee.Date.fromYMD(y, 1, 15).millis(), model: 'rf300_v1'});
}
for (var y = 1988; y <= 2025; y++) {
  Export.image.toAsset({image: classifySeason(y), description: 'ch4_ls_class_' + y,
    assetId: ROOT + 'ls_class/ls_class_' + y, region: frame, scale: SCALE, crs: CRS,
    maxPixels: 1e11, pyramidingPolicy: {'class': 'mode', 'pmax': 'mean'}});
}
