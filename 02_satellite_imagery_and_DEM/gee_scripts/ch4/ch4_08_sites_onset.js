// =====================================================================
// DRAFT, not yet run (project in restricted mode until 1 Nov 2026).
// Chapter 4 · Script 08 · Extraction sites, onset and persistence
// Project: riversand2024 · 21 Oct 2026
// ---------------------------------------------------------------------
// Input: ch4/ls_c3_final_v1b (excavated ground, 1988-2025, Script 06 v1b).
// Per cell: number of seasons mapped (n), first and last season mapped.
// Sites: cells mapped in at least 2 seasons, grown by 300 m and joined into
// polygons at 90 m (8-connected); each site gets its annual mapped area
// (km2) for every season, its centroid and its sector.
// Outputs: table ch4/sites_v1b (one feature per site, fields a1988..a2025,
// n_seasons_any, first, last, peak_year, sector, lon, lat, site_km2).
// Mapped areas are before error adjustment (Script 07 sample).
// =====================================================================
var ROOT = 'projects/riversand2024/assets/ch4/';
var frame = ee.Geometry.Rectangle([88.0, 24.75, 92.55, 26.70], null, false);
var cm = ee.Image(ROOT + 'corridor_mask');
var START = 1988, END = 2025;

var imgs = [];
for (var y = START; y <= END; y++) {
  imgs.push(ee.Image(ROOT + 'ls_c3_final_v1b/ls_c3_final_v1b_' + y).unmask(0).rename('a' + y));
}
var stack = ee.Image.cat(imgs);
var n = stack.reduce(ee.Reducer.sum()).rename('n');
var yr = ee.ImageCollection(imgs.map(function (im, k) {
  return im.multiply(START + k).selfMask().rename('y').toInt16();
}));
var first = yr.min().rename('first'), last = yr.max().rename('last');

var ever = n.gte(2);
var grown = ever.selfMask().focalMax(300, 'circle', 'meters').reproject('EPSG:32646', null, 90);
var sites = grown.updateMask(cm.gt(0)).toInt().reduceToVectors({geometry: frame, scale: 90, crs: 'EPSG:32646',
  geometryType: 'polygon', eightConnected: true, maxPixels: 1e10, tileScale: 8});

var areaStack = stack.multiply(ee.Image.pixelArea().divide(1e6));
var withArea = areaStack.reduceRegions({collection: sites, reducer: ee.Reducer.sum(), scale: 30,
  crs: 'EPSG:32646', tileScale: 8});
var withSector = cm.rename('sector').reduceRegions({collection: withArea, reducer: ee.Reducer.mode(),
  scale: 90, crs: 'EPSG:32646', tileScale: 8});
var out = withSector.map(function (f) {
  var c = f.geometry().centroid(10).coordinates();
  return f.set({lon: c.get(0), lat: c.get(1), site_km2: f.geometry().area(10).divide(1e6)}).setGeometry(null);
});
Export.table.toAsset({collection: out, description: 'ch4_sites_v1b', assetId: ROOT + 'sites_v1b'});

// Rasters for the maps (RGB encodes first season - 1980, seasons mapped, water %).
var water = ee.Image('JRC/GSW1_4/GlobalSurfaceWater').select('occurrence').unmask(0);
var rgb = first.subtract(1980).unmask(0).updateMask(ever.or(1)).multiply(ever)
  .addBands(n.multiply(ever)).addBands(water).toUint8().rename(['r', 'g', 'b']);
