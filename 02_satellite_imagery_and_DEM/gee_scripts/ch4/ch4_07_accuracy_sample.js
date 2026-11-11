// =====================================================================
// Chapter 4 · Script 07 · Stratified reference sample for the accuracy
// assessment (Olofsson et al., 2014)
// Project: riversand2024 · Oct 2026 (version v1b, 20 Oct 2026)
// ---------------------------------------------------------------------
// Assessed seasons: 2012 (Google Earth imagery of Nov 2010 to Mar 2012 is
// widespread) and 2023 (Maxar imagery of Nov 2022 to Apr 2023 in Esri
// Wayback). The quantity H1 needs is the area of excavated ground, so the
// assessment is of class 3 against all other classes, with three strata:
//   A  mapped excavated ground in the season (v1b rules, as Script 06);
//   B  not mapped, but within 1 km of excavated ground mapped in any season
//      1988-2025 (first post-processed v1 run, ch4/ls_c3_final) or of an
//      inventory record of extraction, where omission is most likely;
//   C  the rest of the corridor.
// Allocation per season: A 150, B 150, C 100 points (400; 800 in all).
// Points closer than 300 m to a training sample are excluded. Each point is
// labelled into the five classes from sub-meter imagery dated within the
// season, without reference to the map (labeller page, blind set).
// Output: ch4/acc_sample_v1b (season, stratum, W, corridor_m2, map_c3,
// rf_class, lon, lat).
// =====================================================================
var ROOT = 'projects/riversand2024/assets/ch4/';
var frame = ee.Geometry.Rectangle([88.0, 24.75, 92.55, 26.70], null, false);
var cm = ee.Image(ROOT + 'corridor_mask');
var corridor = cm.gt(0);
var inv = ee.FeatureCollection(ROOT + 'inventory_411')
  .filter(ee.Filter.inList('label_group', ['Sand Mining', 'Other extraction']));
var train = ee.FeatureCollection(ROOT + 'training_ls');

function r1(y) {
  var c = ee.Image(ROOT + 'ls_class/ls_class_' + y);
  var d = ee.Image(ROOT + 'ls_dry/ls_dry_' + y);
  return c.select('class').eq(3).and(d.select('MNDWI').lt(-500)).and(d.select('MNDWI_p80').lt(0)).unmask(0);
}
function finalC3(y) {
  var s = r1(y - 1).add(r1(y)).add(r1(y + 1));
  var keep = r1(y).and(s.gte(2));
  return keep.and(keep.selfMask().connectedPixelCount(10, true).gte(3).unmask(0));
}
var ever = ee.ImageCollection(ROOT + 'ls_c3_final').max().unmask(0);
var near = ever.selfMask().fastDistanceTransform(64).sqrt().multiply(30).lte(1000).unmask(0)
  .or(ee.Image(0).paint(inv.map(function (f) { return f.buffer(1000, 50); }), 1).eq(1));
var trainMask = ee.Image(0).paint(train.map(function (f) { return f.buffer(300, 50); }), 1).eq(0);

var ALLOC = [150, 150, 100];
function sampleSeason(y) {
  var c3 = finalC3(y);
  var stratum = ee.Image(3).where(near, 2).where(c3.eq(1), 1).updateMask(corridor).rename('stratum').toInt();
  var area = ee.Image.pixelArea().addBands(stratum)
    .reduceRegion({reducer: ee.Reducer.sum().group(1, 'stratum'), geometry: frame, scale: 30,
      crs: 'EPSG:32646', maxPixels: 1e11, tileScale: 8});
  var groups = ee.List(area.get('groups'));
  var total = ee.Number(groups.iterate(function (g, s) { return ee.Number(s).add(ee.Dictionary(g).get('sum')); }, 0));
  var W = ee.Dictionary(groups.iterate(function (g, d) {
    g = ee.Dictionary(g);
    return ee.Dictionary(d).set(ee.Number(g.get('stratum')).format('%d'), ee.Number(g.get('sum')).divide(total));
  }, ee.Dictionary({})));
  var img = stratum.addBands(c3.rename('map_c3'))
    .addBands(ee.Image(ROOT + 'ls_class/ls_class_' + y).select('class').rename('rf_class'));
  var pts = img.updateMask(trainMask).stratifiedSample({numPoints: 0, classBand: 'stratum',
    classValues: [1, 2, 3], classPoints: ALLOC, region: frame, scale: 30,
    seed: 1000 + y, geometries: true, tileScale: 16});
  return pts.map(function (f) {
    var xy = f.geometry().coordinates();
    return f.set({season: y, W: W.get(ee.Number(f.get('stratum')).format('%d')),
      corridor_m2: total, lon: xy.get(0), lat: xy.get(1)});
  });
}
var sample = sampleSeason(2012).merge(sampleSeason(2023));
Export.table.toAsset({collection: sample, description: 'ch4_acc_sample_v1b', assetId: ROOT + 'acc_sample_v1b'});
