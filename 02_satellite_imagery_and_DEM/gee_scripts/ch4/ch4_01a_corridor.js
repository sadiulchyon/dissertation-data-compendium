// =====================================================================
// Chapter 4 · Script 01a · Northern-fringe river corridor
// Project: riversand2024 · Oct 2026
// ---------------------------------------------------------------------
// Corridor = Bangladesh land, north of 24.75 N, within 30 km of the India
// border, within 1 km of surface water mapped by JRC GSW (1984-2021), and
// within a river-size-dependent distance of a HydroSHEDS river with an
// upstream area >= 50 km2 (1.5 km if < 1e4 km2; 3 km if < 2e5 km2; 8 km
// for the Brahmaputra-Jamuna and Teesta). Band 'sector' = 1-6, west to east.
// Output: projects/riversand2024/assets/ch4/corridor_mask (30 m, EPSG:32646)
// =====================================================================
var ROOT = 'projects/riversand2024/assets/ch4/';
var frame = ee.Geometry.Rectangle([88.0, 24.75, 92.55, 26.70], null, false);
var lsib = ee.FeatureCollection('USDOS/LSIB_SIMPLE/2017');
var bgd = lsib.filter(ee.Filter.eq('country_na', 'Bangladesh'));
var ind = lsib.filter(ee.Filter.eq('country_na', 'India')).filterBounds(frame.buffer(40000));

var inBgd = ee.Image(0).paint(bgd, 1).eq(1);
var nearBorder = ind.distance(30000, 50).lte(30000).unmask(0);
var nearWater = ee.Image('JRC/GSW1_4/GlobalSurfaceWater').select('max_extent').unmask(0)
  .eq(1).focalMax({radius: 1000, units: 'meters'});
var ffr = ee.FeatureCollection('WWF/HydroSHEDS/v1/FreeFlowingRivers').filterBounds(frame);
var nearRiver = ffr.filter(ee.Filter.rangeContains('UPLAND_SKM', 50, 1e4)).distance(1500, 30).lte(1500).unmask(0)
  .or(ffr.filter(ee.Filter.rangeContains('UPLAND_SKM', 1e4, 2e5)).distance(3000, 30).lte(3000).unmask(0))
  .or(ffr.filter(ee.Filter.gte('UPLAND_SKM', 2e5)).distance(8000, 50).lte(8000).unmask(0));

var lonlat = ee.Image.pixelLonLat();
var lon = lonlat.select('longitude');
var breaks = [88.0, 89.3, 89.92, 90.3, 90.95, 91.6, 92.55];
var sector = ee.Image(0);
for (var i = 0; i < 6; i++) {
  sector = sector.where(lon.gte(breaks[i]).and(lon.lt(breaks[i + 1])), i + 1);
}
var corridor = inBgd.and(nearBorder).and(nearWater).and(nearRiver)
  .and(lonlat.select('latitude').gte(24.75));
var corridorMask = sector.updateMask(corridor).rename('sector').toByte();

Export.image.toAsset({image: corridorMask, description: 'ch4_corridor_mask',
  assetId: ROOT + 'corridor_mask', region: frame, scale: 30, crs: 'EPSG:32646',
  maxPixels: 1e11, pyramidingPolicy: {sector: 'mode'}});

Map.centerObject(frame, 8);
Map.addLayer(corridorMask.randomVisualizer(), {}, 'corridor by sector');
