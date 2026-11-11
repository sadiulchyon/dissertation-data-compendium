// =====================================================================
// Chapter 4 · Script 06 · Post-classification rules and annual class 3 area
// Project: riversand2024 · Oct 2026 (version v1b, 20 Oct 2026)
// ---------------------------------------------------------------------
// Input: the version 1 random forest maps ch4/ls_class (Script 05), complete
// for 1988-2025. (A version 2 classification, Script 05b, was stopped on
// 20 Oct: it was no better at class 3 in held-out checks, and Earth Engine
// ran it at one or two tasks at a time.)
// A cell is kept as excavated ground in season y when:
//   R1  the RF class is 3, the dry-season median MNDWI is below -0.05 and the
//       80th percentile of MNDWI is below 0 (bare in at least 80% of clear
//       observations, which removes channel margins);
//   R2  R1 holds in at least two of seasons y-1, y and y+1 (persistence;
//       seasons 1988 and 2025 use the two available seasons);
//   R3  the cell belongs to a connected group (8-neighbour) of at least
//       3 cells (0.27 ha) that pass R2 (5 cells in the first run, which
//       dropped small genuine patches).
// Outputs: ch4/ls_c3_final_v1b/ls_c3_final_v1b_YYYY (1 = excavated ground)
// and the table ch4/c3_area_by_sector_v1b (km2 by season and sector).
// Mapped areas are biased by classification error; the accuracy assessment
// (Script 07) gives error-adjusted areas.
// =====================================================================
var ROOT = 'projects/riversand2024/assets/ch4/';
var frame = ee.Geometry.Rectangle([88.0, 24.75, 92.55, 26.70], null, false);
var cm = ee.Image(ROOT + 'corridor_mask');
var START = 1988, END = 2025, MMU = 3;

function r1(y) {
  var c = ee.Image(ROOT + 'ls_class/ls_class_' + y);
  var d = ee.Image(ROOT + 'ls_dry/ls_dry_' + y);
  return c.select('class').eq(3).and(d.select('MNDWI').lt(-500)).and(d.select('MNDWI_p80').lt(0)).unmask(0);
}
function finalC3(y) {
  var lo = Math.max(START, y - 1), hi = Math.min(END, y + 1);
  var s = ee.Image(0);
  for (var k = lo; k <= hi; k++) s = s.add(r1(k));
  var keep = r1(y).and(s.gte(2));
  var mmu = keep.selfMask().connectedPixelCount(10, true).gte(MMU);
  return keep.and(mmu.unmask(0)).updateMask(cm.gt(0)).rename('c3').toUint8()
    .set({season: y, 'system:time_start': ee.Date.fromYMD(y, 1, 15).millis(), rules: 'R1-R3 v1b (MMU 3)'});
}

var rows = [];
for (var y = START; y <= END; y++) {
  var img = finalC3(y);
  Export.image.toAsset({image: img, description: 'ch4_ls_c3_final_v1b_' + y,
    assetId: ROOT + 'ls_c3_final_v1b/ls_c3_final_v1b_' + y, region: frame, scale: 30, crs: 'EPSG:32646',
    maxPixels: 1e11, pyramidingPolicy: {c3: 'max'}});
  var a = ee.Image.pixelArea().divide(1e6).updateMask(img.eq(1)).addBands(cm.rename('s'));
  var d = a.reduceRegion({reducer: ee.Reducer.sum().group(1, 's'), geometry: frame, scale: 30,
    crs: 'EPSG:32646', maxPixels: 1e11, tileScale: 8});
  rows.push(ee.FeatureCollection(ee.List(d.get('groups')).map(function (g) {
    g = ee.Dictionary(g);
    return ee.Feature(ee.Geometry.Point([90.5, 25.5]), {season: y, sector: g.get('s'), c3_km2: g.get('sum')});
  })));
}
Export.table.toAsset({collection: ee.FeatureCollection(rows).flatten(), description: 'ch4_c3_area_by_sector_v1b',
  assetId: ROOT + 'c3_area_by_sector_v1b'});
