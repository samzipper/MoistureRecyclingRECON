// RUN SCRIPT ON GOOGLE EARTH ENGINE:
// https://code.earthengine.google.com/c72b0e45d2fe46c4f8106aaa644c8a2a
// THEN COPY OUTPUT FILE FROM Google Drive TO THIS REPO

// DataPrep_02_ExtractLANID
// This script extracts LANID irrigation status for polygons created using DataPrep_01_CreatePolygonsERA5.R
//

/***************************************
 * INPUTS
 ***************************************/

var grid = ee.FeatureCollection(
  'projects/ee-samzipper/assets/ERA5_grid_CONUS'
);

var lanid = ee.Image(
  'users/xyhuwmir4/LANID/LANID_v1_rse'
);

var years = [2008, 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017];

/***************************************
 * FUNCTION: SUMMARIZE YEAR
 ***************************************/

function summarizeYear(year) {

  // LANID band name
  var yy = String(year).slice(2);
  var lanidBand = 'irMap' + yy;

  // CDL image for same year
  var cdl = ee.ImageCollection('USDA/NASS/CDL')
    .filter(ee.Filter.calendarRange(year, year, 'year'))
    .first()
    .select('cropland');

  var cdlProj = cdl.projection();

  // CDL class 111 = Water
  // Everything else is considered land
  var landArea = cdl
  .neq(111)
  .multiply(ee.Image.pixelArea())
  .rename('land_area_m2');

  // LANID irrigated pixels
  var irrArea = lanid
    .select(lanidBand)
    .unmask(0)
    .eq(1)
    .multiply(ee.Image.pixelArea())
    .rename('irr_area_m2');

  var img = ee.Image.cat([
    irrArea,
    landArea
  ]);

  var stats = img.reduceRegions({
    collection: grid,
    reducer: ee.Reducer.sum(),
    crs: cdlProj,
    scale: cdlProj.nominalScale(),
    tileScale: 4
  });

  return stats.map(function(f) {

    var irrArea = ee.Number(f.get('irr_area_m2'));
    var landArea = ee.Number(f.get('land_area_m2'));

    return ee.Feature(null, {
      cellid: f.get('cellid'),
      year: year,
      irr_area_m2: irrArea,
      land_area_m2: landArea,
      dry_area_m2: landArea.subtract(irrArea),
      irr_prc: irrArea.divide(landArea)
    });
  });
}


/***************************************
 * BUILD OUTPUT TABLE
 ***************************************/

var allYears = ee.FeatureCollection(
  years.map(summarizeYear)
).flatten();

print('Rows:', allYears.size());
print('Example:', allYears.first());


/***************************************
 * EXPORT
 ***************************************/

Export.table.toDrive({
  collection: allYears,
  description: 'ERA5-CONUS_LANID_PolygonCounts_Annual',
  folder: 'GEEexports',
  fileFormat: 'CSV'
});