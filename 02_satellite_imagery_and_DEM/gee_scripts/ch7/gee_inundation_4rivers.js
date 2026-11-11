// Chapter 7, Section 7.3.6: Sentinel-1 floodplain inundation for the Jadukata, Piyain, Sari-Gowain and Lubha rivers.
// Code Editor version of the browser run of 30 Oct 2026 (same zone rule, track, months and thresholds).
// Output: one table per river with, per image date, zone area, observed (valid) zone area and water area at
// VV < -15, -16 and -17 dB (m2). Export to Drive and join to qdaily.csv with analyse.py.

var R = {
  Jadukata:   {box: [91.19, 25.10, 91.32, 25.215], gauges: [[91.2508, 25.1891]], reach: 7000},
  Piyain:     {box: [91.92, 25.11, 92.05, 25.20],  gauges: [[92.0167, 25.1568], [91.967, 25.1624]], reach: 7000},
  SariGowain: {box: [91.955, 25.05, 92.20, 25.17], gauges: [[92.1173, 25.0939], [91.9809, 25.0886]], reach: 9000},
  Lubha:      {box: [92.28, 25.036, 92.37, 25.085], gauges: [[92.3085, 25.0515]], reach: 6000}
};

function zoneOf(name) {
  var r = R[name];
  var box = ee.Geometry.Rectangle(r.box, null, false);
  var bd = ee.FeatureCollection('USDOS/LSIB_SIMPLE/2017').filter(ee.Filter.eq('country_na', 'Bangladesh')).geometry();
  var occ = ee.Image('JRC/GSW1_4/GlobalSurfaceWater').select('occurrence').unmask(0);
  var water30 = occ.gte(30);
  var chanCand = occ.gte(10).selfMask().clip(box);
  var seeds = ee.Image().byte().paint(ee.FeatureCollection(r.gauges.map(function (p) {
    return ee.Feature(ee.Geometry.Point(p).buffer(250)); })), 1).unmask(0).and(occ.gte(10));
  // channel traced from the gauges through pixels with >= 10% water occurrence, limited to r.reach metres
  var cum = chanCand.multiply(0).add(1).cumulativeCost({source: seeds.selfMask(), maxDistance: r.reach, geodeticDistance: true});
  var river = cum.mask().and(cum.gte(0)).unmask(0).gt(0).clip(box);
  var px = ee.Image.pixelArea().sqrt();
  var dRiver = river.fastDistanceTransform(256, 'pixels', 'squared_euclidean').sqrt().multiply(px);
  var dWater = river.or(water30).fastDistanceTransform(64, 'pixels', 'squared_euclidean').sqrt().multiply(px);
  var dem = ee.ImageCollection('JAXA/ALOS/AW3D30/V3_2').select('DSM').mosaic().setDefaultProjection('EPSG:4326', null, 30);
  var slope = ee.Terrain.slope(dem);
  var zone = dRiver.lte(2000).and(dWater.gt(200)).and(slope.lt(5)).and(dem.lte(20))
    .clip(box).clip(bd).selfMask().rename('zone');
  return {box: box, zone: zone};
}

function inundation(name) {
  var Z = zoneOf(name);
  var s1 = ee.ImageCollection('COPERNICUS/S1_GRD').filterBounds(Z.box)
    .filter(ee.Filter.eq('instrumentMode', 'IW'))
    .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VV'))
    .filter(ee.Filter.eq('orbitProperties_pass', 'DESCENDING'))
    .filter(ee.Filter.eq('relativeOrbitNumber_start', 77))
    .filterDate('2015-01-01', '2026-01-01')
    .filter(ee.Filter.calendarRange(6, 9, 'month'))
    .select('VV');
  var dates = ee.List(s1.aggregate_array('system:time_start'))
    .map(function (t) { return ee.Date(t).format('YYYY-MM-dd'); }).distinct();
  var area = ee.Image.pixelArea();
  var rows = dates.map(function (d) {
    var d0 = ee.Date.parse('YYYY-MM-dd', d);
    var vv = s1.filterDate(d0, d0.advance(1, 'day')).mosaic();
    var valid = vv.mask().reduce(ee.Reducer.min()).gt(0);
    var z = Z.zone.unmask(0).gt(0);
    var img = ee.Image.cat([
      area.updateMask(z).rename('zone'),
      area.updateMask(z.and(valid)).rename('valid'),
      area.updateMask(z.and(vv.lt(-15))).rename('w15'),
      area.updateMask(z.and(vv.lt(-16))).rename('w16'),
      area.updateMask(z.and(vv.lt(-17))).rename('w17')]);
    var s = img.reduceRegion({reducer: ee.Reducer.sum(), geometry: Z.box, scale: 10, maxPixels: 1e10, tileScale: 4});
    return ee.Feature(null, s).set({river: name, date: d});
  });
  return ee.FeatureCollection(rows);
}

Object.keys(R).forEach(function (n) {
  Export.table.toDrive({collection: inundation(n), description: 'inun_D77_' + n, fileFormat: 'CSV',
    selectors: ['river', 'date', 'zone', 'valid', 'w15', 'w16', 'w17']});
  Map.addLayer(zoneOf(n).zone, {palette: ['ffb300']}, 'zone ' + n);
});
Map.setCenter(91.8, 25.12, 10);
