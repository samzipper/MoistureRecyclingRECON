// GDE_01_ExportMonthlyETandP.js
// Run in the Google Earth Engine Code Editor to prepare monthly climate inputs.
// Export monthly totals, monthly climatologies, and coverage layers on source grids.
// GDE masking and ET/P calculations are performed locally in later steps.\
//
// Author: Sam + Copilot
// Date: 11/9/2026

// 1. Set analysis period, tile batch, and export options
// Set tileNames to null for all tiles; otherwise batch only the named subset.
// Batches use sorted names and a zero-based start; null batchSize selects all remaining.
// Drive exports accept folder names, not nested paths.
var config = {
  startYear: 2008,
  endYear: 2017,
  driveFolder: 'GlobalGDEs',
  tileAsset: 'projects/ee-samzipper/assets/GDE_Tiles_RohdeEtAl',
  tileNames: ['n35w110'],
  batchStart: 0,
  batchSize: 20,
  exportYearlyStacks: false,
  exportClimatologies: true,
  exportTranspiration: false,
  exportCoverage: false,  //set false for production runs
  showDiagnostics: false,
  fileDimensions: 8192,
  noData: -9999
};
var selectedTiles = ee.FeatureCollection(config.tileAsset).select(['name']);
if (config.tileNames !== null) {
  selectedTiles = selectedTiles.filter(ee.Filter.inList('name', config.tileNames));
}
selectedTiles = selectedTiles.sort('name');
var etBands = config.exportTranspiration ? ['ET', 'T'] : ['ET'];

// 2. Load the PML and TerraClimate products used in the GDE model code
// PML components are daily rates (mm/day); TerraClimate pr is a monthly total (mm).
// Land ET = transpiration (Ec) + soil evaporation (Es) + interception evaporation (Ei).
var pmlId = 'CAS/IGSNRR/PML/V2_v017';
var precipitationId = 'IDAHO_EPSCOR/TERRACLIMATE';
var startDate = ee.Date.fromYMD(config.startYear, 1, 1);
var endDate = ee.Date.fromYMD(config.endYear + 1, 1, 1);
var years = ee.List.sequence(config.startYear, config.endYear);
var months = ee.List.sequence(1, 12);
var yearCount = config.endYear - config.startYear + 1;
// Include preceding composites in case their intervals overlap the analysis start.
var pml = ee.ImageCollection(pmlId)
  .filterDate(startDate.advance(-8, 'day'), endDate)
  .select(['Ec', 'Es', 'Ei'])
  .sort('system:time_start');
var precipitation = ee.ImageCollection(precipitationId)
  .filterDate(startDate, endDate)
  .select('pr')
  .sort('system:time_start');
// Retain each product's original grid; do not resample ET and P to a common grid.
var etProjection = ee.Image(pml.first()).select('Ec').projection();
var pProjection = ee.Image(precipitation.first()).projection();

if (config.showDiagnostics) {
  print('PML image', pml.first());
  print('PML grid', etProjection);
  print('Precipitation grid', pProjection);
  print('Precipitation images', precipitation.size());
  print('Expected precipitation images', yearCount * 12);
  Map.centerObject(selectedTiles, 6);
  Map.addLayer(selectedTiles, {}, 'Selected GDE tiles');
  print('Selected tile names', selectedTiles.aggregate_array('name'));
}

// 3. Define calendar helpers and check the assumed input schedule
// These client-side helpers also allow calendar arithmetic to be tested locally.
// PML intervals last eight days, except for the shortened final interval of a year.
function intervalEndMillis(startMillis) {
  var compositeStart = new Date(startMillis);
  var nextYear = Date.UTC(compositeStart.getUTCFullYear() + 1, 0, 1);
  return Math.min(startMillis + 8 * 86400000, nextYear);
}

// Count only the part of a composite interval falling within a calendar month.
function overlapDays(startMillis, endMillis, monthStartMillis, monthEndMillis) {
  return Math.max(0, Math.min(endMillis, monthEndMillis) -
    Math.max(startMillis, monthStartMillis)) / 86400000;
}

// Stop exports if PML dates differ from the expected January 1-based schedule,
// or if TerraClimate does not contain exactly one image for every requested month.
function validateCalendar(pmlTimes, precipitationTimes) {
  var expectedPml = [];
  var expectedPrecipitation = [];
  for (var year = config.startYear; year <= config.endYear; year++) {
    var yearEnd = Date.UTC(year + 1, 0, 1);
    for (var time = Date.UTC(year, 0, 1); time < yearEnd; time += 8 * 86400000) {
      expectedPml.push(time);
    }
    for (var month = 0; month < 12; month++) {
      expectedPrecipitation.push(Date.UTC(year, month, 1));
    }
  }
  var selectedPmlTimes = pmlTimes.filter(function(time) {
    return time >= Date.UTC(config.startYear, 0, 1);
  });
  if (JSON.stringify(selectedPmlTimes) !== JSON.stringify(expectedPml)) {
    throw new Error('PML timestamps are missing, duplicated, or do not follow ' +
      'the assumed January 1-based 8-day calendar. No exports created.');
  }
  if (JSON.stringify(precipitationTimes) !== JSON.stringify(expectedPrecipitation)) {
    throw new Error('Expected exactly one TerraClimate image per calendar ' +
      'month. No exports created.');
  }
}

// 4. Prepare ET rates and optional transpiration for each PML composite
// Require all three ET components to be present and nonnegative at each pixel.
// ET_water is excluded; these outputs represent land ET, not open-water evaporation.
function prepareComposite(image) {
  var compositeStart = ee.Date(image.get('system:time_start'));
  var nextYear = ee.Date.fromYMD(ee.Number(compositeStart.get('year')).add(1), 1, 1);
  var compositeEnd = compositeStart.advance(8, 'day').millis().min(nextYear.millis());
  var components = image.select(['Ec', 'Es', 'Ei']);
  var valid = components.mask().reduce(ee.Reducer.min())
    .and(components.gte(0).reduce(ee.Reducer.min()));
  var rates = components.reduce(ee.Reducer.sum()).rename('ET');
  if (config.exportTranspiration) {
    rates = rates.addBands(image.select('Ec').rename('T'));
  }
  rates = rates.updateMask(valid);
  return rates.set({
    'system:time_start': compositeStart.millis(),
    interval_end: compositeEnd
  });
}

var composites = pml.map(prepareComposite);

// 5. Integrate PML rates into monthly depths and retain monthly precipitation
function monthlyTotals(year, month) {
  year = ee.Number(year);
  month = ee.Number(month);
  var monthStart = ee.Date.fromYMD(year, month, 1);
  var monthEnd = monthStart.advance(1, 'month');
  var expectedDays = monthEnd.difference(monthStart, 'day');
  // Select composites overlapping this month, including those starting last month.
  var overlapping = composites
    .filter(ee.Filter.lt('system:time_start', monthEnd.millis()))
    .filter(ee.Filter.gt('interval_end', monthStart.millis()));
  // Weight daily rates by overlap days, assuming a constant rate within each interval.
  var integrated = overlapping.map(function(image) {
    var days = ee.Number(image.get('interval_end')).min(monthEnd.millis())
      .subtract(ee.Number(image.get('system:time_start')).max(monthStart.millis()))
      .divide(86400000).max(0);
    var depths = image.select(etBands).multiply(days);
    var coverage = image.select('ET').multiply(0).add(days)
      .unmask(0).rename('ET_valid_days');
    return depths.addBands(coverage).toFloat();
  }).sum().setDefaultProjection(etProjection);
  // Mask incomplete ET months rather than treating missing days as zero evaporation.
  var etDays = integrated.select('ET_valid_days');
  var et = integrated.select(etBands).updateMask(etDays.eq(expectedDays));
  // Precipitation is already accumulated by month; no daily integration is needed.
  var monthlyP = precipitation.filterDate(monthStart, monthEnd).map(function(image) {
    return image.updateMask(image.gte(0));
  });
  var pCount = monthlyP.count().unmask(0).rename('P_valid_months')
    .setDefaultProjection(pProjection);
  var p = monthlyP.sum().rename('P').updateMask(pCount.eq(1))
    .setDefaultProjection(pProjection);
  return et.addBands(etDays).addBands(p).addBands(pCount).set({
    'system:time_start': monthStart.millis(),
    year: year,
    month: month,
    expected_days: expectedDays,
    'system:index': monthStart.format('YYYY_MM')
  });
}

// Build 120 monthly images for 2008-2017, with dates and coverage bands attached.
var monthly = ee.ImageCollection.fromImages(years.map(function(year) {
  return months.map(function(month) {
    return monthlyTotals(year, month);
  });
}).flatten()).sort('system:time_start');

// Interactive previews are disabled for production; only diagnostics combine tile geometries.
if (config.showDiagnostics) {
  var previewRegion = selectedTiles.geometry();
  var previewStart = ee.Date.fromYMD(config.startYear, 7, 1);
  var july = ee.Image(monthly.filterDate(previewStart, previewStart.advance(1, 'month')).first());
  print('July bands and properties', july);
  var julyVis = {
    min: 0,
    max: 400,
    palette: ['#b2182b', '#ef8a62', '#f7f7f7', '#67a9cf', '#2166ac']
  };
  Map.addLayer(july.select('ET').clip(previewRegion),
    julyVis, 'July ET');
  Map.addLayer(july.select('P').clip(previewRegion),
    julyVis, 'July P');

  // Use the same monthly-depth color scale for the two July preview layers.
  var julyLegend = ui.Panel({
    style: {position: 'bottom-left', padding: '8px 12px'}
  });
  julyLegend.add(ui.Label('July ET / P (mm)', {
    fontWeight: 'bold', margin: '0 0 6px 0'
  }));
  julyLegend.add(ui.Thumbnail({
    image: ee.Image.pixelLonLat().select('longitude').multiply(julyVis.max),
    params: {
      bbox: [0, 0, 1, 0.1],
      dimensions: '255x16',
      min: julyVis.min,
      max: julyVis.max,
      palette: julyVis.palette,
      format: 'png'
    },
    style: {width: '255px', height: '16px', margin: '0'}
  }));
  julyLegend.add(ui.Panel({
    widgets: [0, 100, 200, 300, 400].map(function(value, index) {
      return ui.Label(String(value), {
        width: '51px', margin: '4px 0 0 0',
        textAlign: index === 0 ? 'left' : (index === 4 ? 'right' : 'center')
      });
    }),
    layout: ui.Panel.Layout.flow('horizontal')
  }));
  Map.add(julyLegend);
}

// 6. Average each calendar month across years to create 12 monthly climatologies
// Require complete monthly totals in every year; export valid-year counts separately.
function climatology(month) {
  month = ee.Number(month);
  var matchingMonths = monthly.filter(ee.Filter.eq('month', month));
  var etCount = matchingMonths.select('ET').count().unmask(0)
    .rename('ET_valid_years').setDefaultProjection(etProjection);
  var pCount = matchingMonths.select('P').count().unmask(0)
    .rename('P_valid_years').setDefaultProjection(pProjection);
  var etMean = matchingMonths.select(etBands).mean()
    .updateMask(etCount.eq(yearCount)).setDefaultProjection(etProjection);
  var pMean = matchingMonths.select('P').mean()
    .updateMask(pCount.eq(yearCount)).setDefaultProjection(pProjection);
  return etMean.addBands(etCount).addBands(pMean).addBands(pCount).set({
    month: month,
    'system:index': month.format('%02d')
  });
}

var monthlyClimatology = ee.ImageCollection.fromImages(months.map(climatology));

// 7. Prepare each enabled native-grid stack once, before the tile export loop
function prepareStack(collection, band) {
  var stack = collection.select(band).toBands();
  var names = collection.aggregate_array('system:index').map(function(index) {
    return ee.String(band).cat('_').cat(ee.String(index));
  });
  var sourceProjection = ee.Image(collection.first()).select(band).projection();
  return stack.rename(names).toFloat().setDefaultProjection(sourceProjection);
}

function prepareExportStacks(metadata) {
  var period = config.startYear + '_' + config.endYear;
  var variables = [
    {band: 'ET', coverage: 'ET_valid_days', validYears: 'ET_valid_years',
      projection: metadata.etProjection},
    {band: 'P', coverage: 'P_valid_months', validYears: 'P_valid_years',
      projection: metadata.pProjection}
  ];
  if (config.exportTranspiration) {
    variables.push({band: 'T', projection: metadata.etProjection});
  }
  var outputs = [];
  function addStack(collection, band, prefix, projection, record) {
    outputs.push({
      image: prepareStack(collection, band),
      prefix: prefix,
      projection: projection,
      record: record
    });
  }
  variables.forEach(function(variable) {
    if (config.exportYearlyStacks) {
      for (var year = config.startYear; year <= config.endYear; year++) {
        var annualMonths = monthly.filter(ee.Filter.eq('year', year));
        var prefix = 'GDE_' + variable.band + '_monthly_' + year;
        addStack(annualMonths, variable.band, prefix, variable.projection, {
          prefix: prefix, band: variable.band, year: year, product: 'monthly',
          units: 'mm/month', expected_bands: 12,
          crs: variable.projection.crs,
          transform: JSON.stringify(variable.projection.transform)
        });
        if (variable.coverage && config.exportCoverage) {
          addStack(annualMonths, variable.coverage,
            'GDE_' + variable.coverage + '_monthly_' + year, variable.projection, null);
        }
      }
    }
    if (config.exportClimatologies) {
      var climatologyPrefix = 'GDE_' + variable.band + '_monthly_climatology_' + period;
      addStack(monthlyClimatology, variable.band, climatologyPrefix, variable.projection, {
        prefix: climatologyPrefix, band: variable.band, product: 'monthly_climatology',
        units: 'mm/month', expected_bands: 12,
        crs: variable.projection.crs,
        transform: JSON.stringify(variable.projection.transform)
      });
      if (variable.validYears && config.exportCoverage) {
        addStack(monthlyClimatology, variable.validYears,
          'GDE_' + variable.validYears + '_' + period, variable.projection, null);
      }
    }
  });
  return outputs;
}

// Only clipping and task creation are tile-specific; retain source grids and filenames.
function exportStack(output, tile) {
  var tileRegion = ee.Geometry(tile.geometry);
  var tilePrefix = output.prefix + '_' + tile.properties.name;
  Export.image.toDrive({
    image: output.image.clip(tileRegion).unmask(config.noData),
    description: tilePrefix,
    folder: config.driveFolder,
    fileNamePrefix: tilePrefix,
    region: tileRegion,
    crs: output.projection.crs,
    crsTransform: output.projection.transform,
    maxPixels: 1e13,
    shardSize: 256,
    fileDimensions: config.fileDimensions,
    skipEmptyTiles: true,
    fileFormat: 'GeoTIFF',
    formatOptions: {cloudOptimized: true, noData: config.noData}
  });
}

// 8. Validate tile names and choose a deterministic batch before creating any tasks
function isIntegerNumber(value) {
  return typeof value === 'number' && isFinite(value) && Math.floor(value) === value;
}

function selectTileBatch(tiles) {
  if (!isIntegerNumber(config.batchStart) || config.batchStart < 0 ||
      (config.batchSize !== null &&
        (!isIntegerNumber(config.batchSize) || config.batchSize <= 0))) {
    throw new Error('batchStart must be a nonnegative integer and batchSize ' +
      'must be a positive integer or null. No export tasks created.');
  }
  var tileNames = tiles.map(function(tile) {
    return tile.properties && tile.properties.name;
  });
  if (tileNames.some(function(name, index) {
    return typeof name !== 'string' || name.length === 0 || tileNames.indexOf(name) !== index;
  })) {
    throw new Error('Tile names must be nonempty, unique strings. No export tasks created.');
  }
  if (config.tileNames !== null && (tiles.length !== config.tileNames.length ||
      config.tileNames.some(function(name) { return tileNames.indexOf(name) === -1; }))) {
    throw new Error('Each requested tile name must match exactly one feature. ' +
      'Requested: ' + config.tileNames.join(', ') + '; found: ' + tileNames.join(', ') +
      '. No export tasks created.');
  }
  var sortedTiles = tiles.slice().sort(function(first, second) {
    var firstName = first.properties.name;
    var secondName = second.properties.name;
    return firstName < secondName ? -1 : (firstName > secondName ? 1 : 0);
  });
  var batchEnd = config.batchSize === null ? tiles.length : config.batchStart + config.batchSize;
  var batch = sortedTiles.slice(config.batchStart, batchEnd);
  if (batch.length === 0) {
    throw new Error('No tiles in this batch. Check batchStart, batchSize, and tileNames. ' +
      'No export tasks created.');
  }
  return batch;
}

function createExports(metadata) {
  validateCalendar(metadata.pmlTimes, metadata.precipitationTimes);
  [metadata.etProjection, metadata.pProjection].forEach(function(projection) {
    if (!projection || typeof projection.crs !== 'string' ||
        !Array.isArray(projection.transform) || projection.transform.length !== 6 ||
        !projection.transform.every(function(value) {
          return typeof value === 'number' && isFinite(value);
        })) {
      throw new Error('Source band CRS or affine transform is missing or invalid. ' +
        'No export tasks created.');
    }
  });
  if (!metadata.tiles || !Array.isArray(metadata.tiles.features)) {
    throw new Error('Tile metadata must contain a GeoJSON features array. ' +
      'No export tasks created.');
  }
  var batch = selectTileBatch(metadata.tiles.features);
  var outputs = prepareExportStacks(metadata);
  print('Export batch', {
    period: config.startYear + '_' + config.endYear,
    matched_tiles: metadata.tiles.features.length,
    batch_start: config.batchStart,
    batch_end_exclusive: config.batchStart + batch.length,
    tile_names: batch.map(function(tile) { return tile.properties.name; }),
    raster_tasks: batch.length * outputs.length,
    manifest_tasks: batch.length
  });
  batch.forEach(function(tile) {
    createTileExports(outputs, tile);
  });
}

// Create separate raster files and a manifest for each selected tile.
function createTileExports(outputs, tile) {
  var tileName = tile.properties.name;
  var period = config.startYear + '_' + config.endYear;
  var manifest = [];
  outputs.forEach(function(output) {
    exportStack(output, tile);
    if (output.record) {
      manifest.push(ee.Feature(null, output.record)
        .set('prefix', output.prefix + '_' + tileName));
    }
  });
  // Record source products, units, grids, and processing assumptions in a CSV manifest.
  var provenance = ee.FeatureCollection(manifest).map(function(feature) {
    return feature.set({
      pml_collection: pmlId,
      precipitation_collection: precipitationId,
      start_date: config.startYear + '-01-01',
      end_date_exclusive: (config.endYear + 1) + '-01-01',
      et_definition: 'Ec + Es + Ei; ET_water excluded; T = Ec',
      temporal_integration: 'Daily rate times calendar-month overlap days; year-end truncated',
      missing_data_policy: 'Complete months and all years required; no gap filling',
      coverage_units: 'ET_valid_days: days; P_valid_months: count; valid_years: count',
      resampling: 'None; each variable exported on its source grid',
      no_data: config.noData,
      tile_name: tileName,
      tile_asset: config.tileAsset,
      region_geojson: JSON.stringify(tile.geometry),
      generated_utc: new Date().toISOString()
    });
  });
  Export.table.toDrive({
    collection: provenance,
    description: 'GDE_monthly_climate_manifest_' + period + '_' + tileName,
    folder: config.driveFolder,
    fileFormat: 'CSV'
  });
}

// 9. Request input metadata and create tasks when it becomes available
// Start the generated export tasks manually from the Earth Engine Tasks tab.
print('PML collection is deprecated; no newer version is substituted.');
print('ET is land ET (Ec + Es + Ei), not Rohde et al.\'s unweighted Ec predictor.');
print('Tile selection: ' + (config.tileNames === null ? 'all asset tiles' : config.tileNames.join(', ')));
print('Drive export folder: ' + config.driveFolder);
if (config.showDiagnostics) {
  print('Monthly image count', monthly.size());
  print('Climatology image count', monthlyClimatology.size());
}
ee.Dictionary({
  pmlTimes: pml.aggregate_array('system:time_start'),
  precipitationTimes: precipitation.aggregate_array('system:time_start')
}).evaluate(function(metadata, error) {
  if (error) {
    print('Input metadata could not be loaded. No export tasks created.', error);
    return;
  }
  var etBand = ee.Image(pml.first()).select('Ec').getInfo().bands[0];
  var pBand = ee.Image(precipitation.first()).select('pr').getInfo().bands[0];
  metadata.etProjection = {crs: etBand.crs, transform: etBand.crs_transform};
  metadata.pProjection = {crs: pBand.crs, transform: pBand.crs_transform};
  // Evaluate the collection directly to obtain GeoJSON, not a nested collection value.
  selectedTiles.evaluate(function(tileCollection, tileError) {
    if (tileError) {
      print('Tile boundaries could not be loaded. No export tasks created.', tileError);
      return;
    }
    metadata.tiles = tileCollection;
    createExports(metadata);
  });
});
