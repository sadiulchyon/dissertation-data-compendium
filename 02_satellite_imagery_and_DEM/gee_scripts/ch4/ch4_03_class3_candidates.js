// =====================================================================
// Chapter 4 · Script 03 · Candidate extraction-surface points for review
// Project: riversand2024 · Oct 2026
// ---------------------------------------------------------------------
// Candidates for class 3 (excavated or disturbed ground) are drawn, for five
// review seasons (2005, 2010, 2015, 2020, 2024), from pixels that are
//   (a) bare in the dry-season median (MNDWI < 0),
//   (b) not green at any time of the hydrological year (Landsat NDVI p90 < 0.35),
//       which excludes dry-season fallow and boro cropland,
//   (c) never water in JRC GSW yearly history 1984-1999 (outside the
//       pre-2000 channel belt, so natural bars are excluded), and
//   (d) within the corridor and within 2 km of an inventory record of
//       extraction, a named site or river, or an uncertain label, or inside
//       the Jaflong and Bholaganj study polygons.
// 30 points per season (stratified random, seed = season). Candidates are
// labelled in Google Earth historical imagery (KML/CSV in chapter-4/working/
// training/) before they enter training; rejected points are kept as
// examples of their true class.
// =====================================================================
var ROOT = 'projects/riversand2024/assets/ch4/';
var corridor = ee.Image(ROOT + 'corridor_mask').gt(0);
var inv = ee.FeatureCollection(ROOT + 'inventory_411');
var polys = ee.FeatureCollection(ROOT + 'study_polygons');
var ext = inv.filter(ee.Filter.inList('label_group',
  ['Sand Mining', 'Other extraction', 'Named feature', 'River/floodplain', 'Uncertain label']));
var zone = ext.map(function (f) { return f.buffer(2000, 50); }).union(50).geometry()
  .union(polys.geometry(), 50);

var pre2000 = ee.ImageCollection('JRC/GSW1_4/YearlyHistory')
  .filter(ee.Filter.calendarRange(1984, 1999, 'year'))
  .map(function (i) { return i.gte(2); }).max().unmask(0);

var epochs = [2005, 2010, 2015, 2020, 2024];
var cands = ee.FeatureCollection(epochs.map(function (y) {
  var dry = ee.Image(ROOT + 'ls_dry/ls_dry_' + y);
  var p90 = ee.Image(ROOT + 'ls_ndvi_p90/ls_ndvi_p90_' + y);
  var c = p90.lt(3500).and(dry.select('MNDWI').lt(0)).and(pre2000.eq(0)).and(corridor)
    .clip(zone).rename('c').selfMask();
  return c.stratifiedSample({numPoints: 30, classBand: 'c', region: zone, scale: 30,
    seed: y, geometries: true, tileScale: 8}).map(function (f) { return f.set('epoch', y); });
})).flatten();

Export.table.toAsset({collection: cands, description: 'ch4_class3_candidates',
  assetId: ROOT + 'class3_candidates'});
// The 150 review points of 19 Oct 2026 were drawn with the same rules from
// on-the-fly composites before the composite assets were complete.
