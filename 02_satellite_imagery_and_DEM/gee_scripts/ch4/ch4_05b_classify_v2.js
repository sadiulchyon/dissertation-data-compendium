// =====================================================================
// Chapter 4 · Script 05b · Random forest classification, version 2
// Project: riversand2024 · Oct 2026
// ---------------------------------------------------------------------
// Version 1 (Script 05) under-mapped excavated ground: at its own class 3
// training points it returned class 3 at only 22 of 37, confusing crushing
// and stockpile yards with built land and pit margins with sediment. After
// the post-classification rules, sectors S1 to S5 held almost no class 3 in
// any season, although labelled imagery shows pit fields and sorting yards
// near Tetulia (S1) and pits on the Someshwari (S4).
// Changes in version 2:
//   (a) built-land samples from the stable rules are dropped within 1.5 km
//       of inventory records of extraction, where GHSL built-up includes
//       crushing yards;
//   (b) each labelled candidate contributes its 3 x 3 neighbourhood (9
//       cells) in its own season; neighbours of class 3 points are kept
//       only if bare (MNDWI < 0) and never green (NDVI p90 < 0.35);
//   (c) stable-rule samples are capped at 250 per class.
// Check, 20 Oct 2026: with one third of labelled candidates (and all their
// neighbours) and 30% of the other samples held out, overall accuracy 0.94;
// class 3 recall 0.88 (50 of 57) and precision 0.86 (50 of 58); 10 of 11
// held-out class 3 candidates were mapped as class 3. Held-out samples are
// not independent reference data (Script 07 gives the accuracy assessment).
// Output: ch4/ls_class_v2b/ls_class_v2b_YYYY (class, pmax).
// =====================================================================
var ROOT = 'projects/riversand2024/assets/ch4/';
var CRS = 'EPSG:32646', SCALE = 30;
var frame = ee.Geometry.Rectangle([88.0, 24.75, 92.55, 26.70], null, false);
var cmask = ee.Image(ROOT + 'corridor_mask');
var corridor = cmask.gt(0);

// Reflectance = DN x 0.0000275 - 0.2 must be applied before the index: an
// NDVI of raw Collection 2 DN is compressed toward 0 by the offset. (The
// first submission of v2 omitted it, which made every cell look never green;
// it was cancelled on 20 Oct, 2:45 AM CDT.) TM/ETM+ NDVI is taken from the
// harmonized bands as in Script 01c/04 (Roy et al., 2016 coefficients).
function qa(i) {
  return i.updateMask(i.select('QA_PIXEL').bitwiseAnd(62).eq(0).and(i.select('QA_RADSAT').eq(0)));
}
function ndT(i) {
  var x = qa(i).select(['SR_B4', 'SR_B3']).multiply(0.0000275).add(-0.2);
  var nir = x.select('SR_B4').multiply(0.8462).add(0.0412), red = x.select('SR_B3').multiply(0.9047).add(0.0061);
  return nir.subtract(red).divide(nir.add(red)).rename('N');
}
function ndO(i) {
  return qa(i).select(['SR_B5', 'SR_B4']).multiply(0.0000275).add(-0.2).normalizedDifference(['SR_B5', 'SR_B4']).rename('N');
}
function p90(y) {   // annual NDVI p90, hydrological year May (y-1) to April (y)
  var s = ee.Date.fromYMD(y - 1, 5, 1), e = ee.Date.fromYMD(y, 5, 1);
  var c = ee.ImageCollection('LANDSAT/LT05/C02/T1_L2').filterBounds(frame).filterDate(s, e).map(ndT);
  if (y >= 1999 && y <= 2021) c = c.merge(ee.ImageCollection('LANDSAT/LE07/C02/T1_L2').filterBounds(frame)
    .filterDate(s, e).filterDate('1999-01-01', '2021-01-01').map(ndT));
  if (y >= 2013) c = c.merge(ee.ImageCollection('LANDSAT/LC08/C02/T1_L2').filterBounds(frame).filterDate(s, e).map(ndO));
  if (y >= 2022) c = c.merge(ee.ImageCollection('LANDSAT/LC09/C02/T1_L2').filterBounds(frame).filterDate(s, e).map(ndO));
  return c.select('N').reduce(ee.Reducer.percentile([90])).rename('NDVI_P90').multiply(10000).round().toInt16();
}
var pre = ee.ImageCollection('JRC/GSW1_4/YearlyHistory').filter(ee.Filter.calendarRange(1984, 1999, 'year'))
  .map(function (i) { return i.gte(2); }).max().unmask(0).rename('PRE2000_WATER');
var dist = pre.selfMask().fastDistanceTransform(256).sqrt().multiply(30).rename('DIST_PRE2000').unmask(7680);
function pred(y) {
  var d = ee.Image(ROOT + 'ls_dry/ls_dry_' + y);
  var t = d.select('NIR').reduceNeighborhood(ee.Reducer.stdDev(), ee.Kernel.square(1)).rename('NIR_SD3');
  return d.addBands(p90(y)).addBands(t).addBands(pre).addBands(dist).addBands(cmask.rename('sector'));
}
var BANDS = ['BLUE', 'GREEN', 'RED', 'NIR', 'SWIR1', 'SWIR2', 'NDVI', 'MNDWI', 'NDBI', 'BSI', 'MNDWI_p20',
  'MNDWI_p80', 'NDVI_p20', 'NDVI_p80', 'NDVI_P90', 'NIR_SD3', 'PRE2000_WATER', 'DIST_PRE2000', 'sector'];

var inv = ee.FeatureCollection(ROOT + 'inventory_411').filter(ee.Filter.inList('label_group', ['Sand Mining', 'Other extraction']));
var zone = inv.map(function (f) { return f.buffer(1500, 50); }).union(50).geometry();
var base = ee.FeatureCollection(ROOT + 'training_ls').filter(ee.Filter.notNull(BANDS))
  .map(function (f) { return f.set('class', ee.Number(f.get('class')).round().toInt()); });
var stable = ee.FeatureCollection([1, 2, 4, 5].map(function (c) {
  var s = base.filter(ee.Filter.eq('source', 'stable_rule')).filter(ee.Filter.eq('class', c));
  if (c === 5) s = s.filter(ee.Filter.bounds(zone).not());
  return s.randomColumn('r', 7).sort('r').limit(250);
})).flatten();
var rv = ee.FeatureCollection(ROOT + 'class3_review').filter(ee.Filter.inList('confidence', ['H', 'M']));
var offs = [];
[-1, 0, 1].forEach(function (dx) { [-1, 0, 1].forEach(function (dy) { offs.push([dx * 0.0003, dy * 0.00027]); }); });
var SEASONS = [2005, 2010, 2011, 2012, 2013, 2014, 2015, 2017, 2018, 2019, 2020, 2023, 2024];
var aug = ee.FeatureCollection(SEASONS.map(function (y) {
  var pts = rv.filter(ee.Filter.eq('season', y));
  var nb = ee.FeatureCollection(offs.map(function (o) {
    return pts.map(function (f) {
      var c = f.geometry().coordinates();
      return ee.Feature(ee.Geometry.Point([ee.Number(c.get(0)).add(o[0]), ee.Number(c.get(1)).add(o[1])]),
        {id: f.get('id'), label: f.get('label'), centre: (o[0] === 0 && o[1] === 0) ? 1 : 0});
    });
  })).flatten();
  return pred(y).select(BANDS).sampleRegions({collection: nb, properties: ['id', 'label', 'centre'], scale: 30})
    .map(function (f) { return f.set({'class': ee.Number(f.get('label')).toInt(), season: y, source: 'review_aug'}); });
})).flatten().filter(ee.Filter.or(ee.Filter.eq('centre', 1), ee.Filter.neq('class', 3),
  ee.Filter.and(ee.Filter.lt('MNDWI', 0), ee.Filter.lt('NDVI_P90', 3500))));
var pits = base.filter(ee.Filter.eq('source', 'lubha_pits'));
var samples = stable.merge(aug).merge(pits);

var rf = ee.Classifier.smileRandomForest({numberOfTrees: 300, seed: 42})
  .train({features: samples, classProperty: 'class', inputProperties: BANDS}).setOutputMode('MULTIPROBABILITY');
function classifySeason(y) {
  var probs = pred(y).select(BANDS).classify(rf);
  var cls = probs.arrayArgmax().arrayGet(0).add(1).rename('class').toUint8();
  var pmax = probs.arrayReduce(ee.Reducer.max(), [0]).arrayGet(0).multiply(100).round().rename('pmax').toUint8();
  return cls.addBands(pmax).updateMask(corridor)
    .set({season: y, 'system:time_start': ee.Date.fromYMD(y, 1, 15).millis(), model: 'rf300_v2b'});
}
for (var y = 1988; y <= 2025; y++) {
  Export.image.toAsset({image: classifySeason(y), description: 'ch4_ls_class_v2b_' + y,
    assetId: ROOT + 'ls_class_v2b/ls_class_v2b_' + y, region: frame, scale: SCALE, crs: CRS,
    maxPixels: 1e11, pyramidingPolicy: {'class': 'mode', 'pmax': 'mean'}});
}
