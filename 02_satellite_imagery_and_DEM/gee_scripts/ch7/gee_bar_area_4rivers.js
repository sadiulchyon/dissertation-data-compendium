// Chapter 6B: dry-season exposed-bar and water area in the channel belt at each gauge reach, 1988 to 2026,
// from Landsat 5/7/8/9 (Collection 2, Level 2) and Sentinel-2 (S2_SR_HARMONIZED), following Section 6A.3.2.
// Channel belt: within BELT m of the low-flow channel traced from the gauge through JRC GSW occurrence >= 10% for
// REACH m (as in gee_inundation_4rivers.js), intersected with the JRC maximum water extent (1984 to 2021); Bangladesh only.
// Per scene (1 Dec to 15 Apr, cloud-masked, >= 95% of the belt observed): water = MNDWI > 0; bar = MNDWI <= 0 and
// NDVI < 0.25 and red reflectance > 0.12 (bright bare sediment); other = the rest. Areas in km2.
// Run in the Code Editor (project riversand2024, 3 Nov 2026).
var R = {
  Jadukata:   {gauge: [91.2508, 25.1891], reach: 8000, belt: 700},
  PiyainSpill:{gauge: [92.0167, 25.1568], reach: 5000, belt: 400},
  PiyainMain: {gauge: [91.9670, 25.1624], reach: 5000, belt: 400},
  SariGowain: {gauge: [92.1173, 25.0939], reach: 8000, belt: 250},
  Lubha:      {gauge: [92.3085, 25.0515], reach: 6000, belt: 300}
};
var bd = ee.FeatureCollection('USDOS/LSIB_SIMPLE/2017').filter(ee.Filter.eq('country_na', 'Bangladesh')).geometry();
var gsw = ee.Image('JRC/GSW1_4/GlobalSurfaceWater');
var occ = gsw.select('occurrence').unmask(0);
var maxext = gsw.select('max_extent').unmask(0);
function beltOf(name) {
  var r = R[name];
  var g = ee.Geometry.Point(r.gauge);
  var box = g.buffer(r.reach + 3000).bounds();
  var seeds = ee.Image().byte().paint(ee.FeatureCollection([ee.Feature(g.buffer(300))]), 1).unmask(0).and(occ.gte(10));
  var chanCand = occ.gte(10).selfMask().clip(box);
  var cum = chanCand.multiply(0).add(1).cumulativeCost({source: seeds.selfMask(), maxDistance: r.reach, geodeticDistance: true});
  var river = cum.mask().and(cum.gte(0)).unmask(0).gt(0).clip(box);
  var px = ee.Image.pixelArea().sqrt();
  var dRiver = river.fastDistanceTransform(256, 'pixels', 'squared_euclidean').sqrt().multiply(px);
  var belt = dRiver.lte(r.belt).and(maxext.eq(1)).clip(box).clip(bd).selfMask().rename('belt');
  return {box: box, belt: belt, geom: box};
}
// ---- sensors
function maskL(img) {  // Landsat C2 L2
  var qa = img.select('QA_PIXEL');
  var clear = qa.bitwiseAnd(1 << 1).eq(0).and(qa.bitwiseAnd(1 << 3).eq(0)).and(qa.bitwiseAnd(1 << 4).eq(0)).and(qa.bitwiseAnd(1 << 2).eq(0));
  return img.updateMask(clear);
}
function prepL57(img) {
  var sr = img.select(['SR_B2', 'SR_B3', 'SR_B4', 'SR_B5']).multiply(0.0000275).add(-0.2).rename(['green', 'red', 'nir', 'swir']);
  return maskL(img.addBands(sr, null, true)).set('sensor', img.get('SPACECRAFT_ID')).set('cloud', img.get('CLOUD_COVER'));
}
function prepL89(img) {
  var sr = img.select(['SR_B3', 'SR_B4', 'SR_B5', 'SR_B6']).multiply(0.0000275).add(-0.2).rename(['green', 'red', 'nir', 'swir']);
  return maskL(img.addBands(sr, null, true)).set('sensor', img.get('SPACECRAFT_ID')).set('cloud', img.get('CLOUD_COVER'));
}
function prepS2(img) {
  var scl = img.select('SCL');
  var clear = scl.neq(3).and(scl.neq(8)).and(scl.neq(9)).and(scl.neq(10)).and(scl.neq(11)).and(scl.neq(1));
  var sr = img.select(['B3', 'B4', 'B8', 'B11']).multiply(0.0001).rename(['green', 'red', 'nir', 'swir']);
  return img.addBands(sr, null, true).updateMask(clear).set('sensor', 'S2').set('cloud', img.get('CLOUDY_PIXEL_PERCENTAGE'));
}
function drySeason(col) {
  return col.filter(ee.Filter.or(ee.Filter.calendarRange(12, 12, 'month'), ee.Filter.calendarRange(1, 3, 'month'), ee.Filter.and(ee.Filter.calendarRange(1, 15, 'day_of_month'), ee.Filter.calendarRange(4, 4, 'month'))));
}
function classify(img) {
  var mndwi = img.normalizedDifference(['green', 'swir']);
  var ndvi = img.normalizedDifference(['nir', 'red']);
  var water = mndwi.gt(0);
  var bar = mndwi.lte(0).and(ndvi.lt(0.25)).and(img.select('red').gt(0.12));
  var other = water.not().and(bar.not());
  return img.addBands(ee.Image.cat(water.rename('water'), bar.rename('bar'), other.rename('other'), img.select('green').mask().rename('valid')).toFloat());
}
function run(name) {
  var B = beltOf(name); var geom = B.geom; var r = R[name];
  var box = B.box;
  var l5 = ee.ImageCollection('LANDSAT/LT05/C02/T1_L2').filterBounds(box).map(prepL57);
  var l7 = ee.ImageCollection('LANDSAT/LE07/C02/T1_L2').filterBounds(box).filterDate('1999-01-01', '2003-05-31').map(prepL57); // before the SLC failure
  var l8 = ee.ImageCollection('LANDSAT/LC08/C02/T1_L2').filterBounds(box).map(prepL89);
  var l9 = ee.ImageCollection('LANDSAT/LC09/C02/T1_L2').filterBounds(box).map(prepL89);
  var s2 = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED').filterBounds(box).map(prepS2);
  var col = drySeason(l5.merge(l7).merge(l8).merge(l9).merge(s2).filter(ee.Filter.lt('cloud', 40)));
  var beltArea = ee.Image.pixelArea().updateMask(B.belt).reduceRegion({reducer: ee.Reducer.sum(), geometry: geom, scale: 30, maxPixels: 1e9}).get('area');
  var areaImg = ee.Image.pixelArea();
  var out = col.map(function (img) {
    var c = classify(img); var sc = ee.Number(ee.Algorithms.If(ee.String(img.get('sensor')).equals('S2'), 20, 30));
    var stats = ee.Image.cat(areaImg.updateMask(c.select('water')).rename('water'), areaImg.updateMask(c.select('bar')).rename('bar'),
                             areaImg.updateMask(c.select('other')).rename('other'), areaImg.updateMask(c.select('valid')).rename('valid'))
      .updateMask(B.belt).reduceRegion({reducer: ee.Reducer.sum(), geometry: geom, scale: sc, maxPixels: 1e9});
    return ee.Feature(null, {river: name, date: img.date().format('YYYY-MM-dd'), sensor: img.get('sensor'), cloud: img.get('cloud'),
      belt_km2: ee.Number(beltArea).divide(1e6), valid_frac: ee.Number(stats.get('valid')).divide(beltArea),
      water_km2: ee.Number(stats.get('water')).divide(1e6), bar_km2: ee.Number(stats.get('bar')).divide(1e6), other_km2: ee.Number(stats.get('other')).divide(1e6)});
  });
  return ee.FeatureCollection(out).filter(ee.Filter.gte('valid_frac', 0.95));
}
// 3 Nov 2026: first batch exports (c. 02:10 CST, to Drive folder ch6B_bar_area on sandboat025) ran with a bug (cloud mask not
// applied to the reflectance bands); corrected c. 02:50 and recomputed per river with getDownloadURL, read from the Code Editor
// page (window.__csvs2), packed and transcribed to bar_area_<river>.csv. The Drive exports are superseded.
Object.keys(R).forEach(function (n) {
  Export.table.toDrive({collection: run(n), description: 'ch6B_bar_' + n, folder: 'ch6B_bar_area', fileNamePrefix: 'bar_area_' + n, fileFormat: 'CSV',
    selectors: ['river', 'date', 'sensor', 'cloud', 'belt_km2', 'valid_frac', 'water_km2', 'bar_km2', 'other_km2']});
});
Map.centerObject(ee.Geometry.Point(R.PiyainSpill.gauge), 12);
Object.keys(R).forEach(function (n) { Map.addLayer(beltOf(n).belt, {palette: ['ff8800']}, 'belt ' + n); });
