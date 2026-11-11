// =====================================================================
// Chapter 4 · Script 04 · Training samples for the annual classification
// Project: riversand2024 · Oct 2026
// ---------------------------------------------------------------------
// Classes: 1 water; 2 natural bare sediment (bars, channel bed);
//          3 excavated or disturbed ground (pits, spoil, stockpiles, crusher
//            and sorting yards); 4 vegetation and cropland; 5 built and other.
// Sources:
//   classes 1, 2, 4, 5  stable-pixel rules (below), stratified by sector and
//                       sampled in eight training seasons;
//   class 3             (a) candidates labelled in Google Earth (Script 03,
//                       review CSV uploaded as ch4/class3_review), and
//                       (b) Lubha pit polygons 2016-2024 of >= 2,700 m2
//                       (three Landsat pixels; 133 pits), sampled at the centroid.
// Each sample carries the predictor values of its own season (Script 01b
// composite + Script 01c NDVI_P90 + context bands), so one classifier is
// trained across seasons and sensors.
// Output: ch4/training_ls (FeatureCollection with 'class', 'season', 'source').
// =====================================================================
var ROOT = 'projects/riversand2024/assets/ch4/';
var SEASONS = [1990, 1995, 2000, 2005, 2010, 2015, 2020, 2024];
var N_PER = 25;   // per class, sector and season
var corridorMask = ee.Image(ROOT + 'corridor_mask');
var inv = ee.FeatureCollection(ROOT + 'inventory_411');
var polys = ee.FeatureCollection(ROOT + 'study_polygons');
var ext = inv.filter(ee.Filter.inList('label_group', ['Sand Mining', 'Other extraction']));
var keepOut = ee.Image(0).paint(ext.map(function (f) { return f.buffer(1500, 50); }), 1)
  .paint(polys, 1).eq(1);

var pre2000 = ee.ImageCollection('JRC/GSW1_4/YearlyHistory')
  .filter(ee.Filter.calendarRange(1984, 1999, 'year'))
  .map(function (i) { return i.gte(2); }).max().unmask(0).rename('PRE2000_WATER');
var distPre = pre2000.selfMask().fastDistanceTransform(256).sqrt().multiply(30)
  .rename('DIST_PRE2000').unmask(7680);

function predictors(y) {
  var dry = ee.Image(ROOT + 'ls_dry/ls_dry_' + y);
  var p90 = ee.Image(ROOT + 'ls_ndvi_p90/ls_ndvi_p90_' + y);
  var tex = dry.select('NIR').reduceNeighborhood(ee.Reducer.stdDev(), ee.Kernel.square(1))
    .rename('NIR_SD3');
  return dry.addBands(p90).addBands(tex).addBands(pre2000).addBands(distPre)
    .addBands(corridorMask.rename('sector'));
}
function ghslBuilt(y) {
  var e = Math.min(2025, Math.max(1990, Math.round(y / 5) * 5));
  return ee.Image('JRC/GHSL/P2023A/GHS_BUILT_S/' + e).select('built_surface').divide(1e4);
}

function stable(y) {
  var x = predictors(y);
  var water = x.select('MNDWI_p20').gt(1500).and(x.select('NDVI_P90').lt(2000));
  var bare = pre2000.eq(1).and(x.select('BSI').gt(0)).and(x.select('MNDWI').lt(-500))
    .and(x.select('NDVI_P90').lt(3000)).and(keepOut.not());
  var veg = x.select('NDVI_P90').gt(5500).and(x.select('NDVI').gt(3000));
  var built = ghslBuilt(y).gt(0.5).and(x.select('NDVI_P90').lt(4000));
  var cls = ee.Image(0).where(water, 1).where(bare, 2).where(veg, 4).where(built, 5)
    .updateMask(corridorMask.gt(0)).selfMask().rename('class');
  var strata = cls.multiply(10).add(corridorMask).rename('stratum').toInt();
  var s = x.addBands(cls).addBands(strata).stratifiedSample({numPoints: N_PER,
    classBand: 'stratum', region: ee.Geometry.Rectangle([88.0, 24.75, 92.55, 26.70], null, false),
    scale: 30, seed: y, geometries: true, tileScale: 16});
  return s.map(function (f) { return f.set({season: y, source: 'stable_rule'}); });
}

// Class 3 (a): reviewed candidates (uploaded CSV: lon, lat, season, label,
// confidence). Each point is labelled from a dated high-resolution image
// (Google Earth historical imagery or Esri World Imagery Wayback) and is
// sampled in the dry season of that image date, so a candidate drawn for
// epoch 2020 but seen on a Dec 2018 image trains season 2019. Candidates that
// proved to be water, natural sediment, vegetation or built land train those
// classes. Only H and M confidence labels are used.
var review = ee.FeatureCollection(ROOT + 'class3_review')
  .filter(ee.Filter.inList('confidence', ['H', 'M']));
var REVIEW_SEASONS = [2005, 2010, 2011, 2012, 2013, 2014, 2015, 2017, 2018, 2019, 2020, 2023, 2024];
function fromReview(y) {
  var x = predictors(y);
  var pts = review.filter(ee.Filter.eq('season', y));
  return x.sampleRegions({collection: pts, properties: ['label', 'id'], scale: 30, geometries: true})
    .map(function (f) { return f.set({class: ee.Number(f.get('label')).toInt(),
      season: y, source: 'image_review'}); });
}
// Class 3 (b): Lubha pits >= 2,700 m2; pit year YYYY -> dry season YYYY
var pits = ee.FeatureCollection(ROOT + 'lubha_pits_2012_2024')
  .filter(ee.Filter.eq('kind', 'pit_polygon'))
  .filter(ee.Filter.gte('year', 2016))
  .map(function (f) { return f.set('a', f.geometry().area(1)); })
  .filter(ee.Filter.gte('a', 2700));
function fromPits(y) {
  var x = predictors(y);
  var pts = pits.filter(ee.Filter.eq('year', y))
    .map(function (f) { return ee.Feature(f.geometry().centroid(1), {zone: f.get('zone')}); });
  return x.sampleRegions({collection: pts, properties: ['zone'], scale: 30, geometries: true})
    .map(function (f) { return f.set({class: 3, season: y, source: 'lubha_pits'}); });
}

var train = ee.FeatureCollection(SEASONS.map(stable)).flatten()
  .merge(ee.FeatureCollection(REVIEW_SEASONS.map(fromReview)).flatten())
  .merge(ee.FeatureCollection([2016, 2017, 2018, 2019, 2020, 2021, 2022, 2023, 2024].map(fromPits)).flatten());

Export.table.toAsset({collection: train, description: 'ch4_training_ls',
  assetId: ROOT + 'training_ls'});
