// =====================================================================
// TREES EVI extraction - LATEST version (2026-07-17)
// Google Earth Engine (JavaScript, code.earthengine.google.com)
//
// Extracts annual mean & median EVI over a 1-ha footprint around each farm,
// per country x chunk x year, and exports to Drive. One script Run loops
// over ALL countries (no per-country editing).
//
// RUN MODE (edit RUN_MODE below): 'control' or 'treated'. Everything else
// (asset, country field, export folder, de-dup behaviour) follows from it.
//
// CLOUD MASKING: Cloud Score+ (GOOGLE/CLOUD_SCORE_PLUS/V1/S2_HARMONIZED),
// which is reliable across 2017-2024 (QA60 was zeroed for much of
// 2022-early 2024). Per-pixel mask on the cs_cdf band at CLEAR_THRESHOLD.
//
// EVI CLAMPING: exports BOTH the raw EVI reductions (evi_mean, evi_median)
// and clamped-to-[-1,1] versions (evi_mean_clamp, evi_median_clamp), so the
// sensitivity to blow-up pixels can be checked in Stata.
//
// FOOTPRINT: 1 ha (radius ~56.4 m) for BOTH treated and control, matching the
// historical treated-farm extraction.
//
// GEOMETRY: dropped from the export (setGeometry(null)) - not used downstream
// and it bloats the CSVs.
// =====================================================================

// =======================================================
// RUN MODE - edit here only
// =======================================================
var RUN_MODE = 'control';   // 'control'  or  'treated'

var MODES = {
    control: {
        asset:        'projects/ee-claudiorizzi93/assets/random_control_farms_50km',
        countryField: 'country_iso_trees',   // control asset country field
        exportFolder: 'gee_evi_contr_farms_50km',
        bufferLabel:  '50km',                 // change with the ring (40km/60km)
        dedupGrid:    true                    // de-dup near-coincident control points
    },
    treated: {
        asset:        'projects/ee-claudiorizzi93/assets/farm_coordinates_all',
        countryField: 'iso_trees',            // treated asset country field
        exportFolder: 'gee_evi_results',
        bufferLabel:  'treated',
        dedupGrid:    false                   // keep every treated farm (see note)
    }
};
var CONFIG = MODES[RUN_MODE];

// =======================================================
// Other configuration
// =======================================================

// Countries to process in ONE run (CMR dropped: very few observations).
var COUNTRIES = ['KEN', 'SEN', 'TZA', 'UGA', 'MLI', 'GIN'];

var N_CHUNKS = 8;
var years = [2017, 2018, 2019, 2020, 2021, 2022, 2023, 2024];

// Scene-level cloudiness pre-filter (optional with Cloud Score+; drops only
// near-fully-cloudy scenes early).
var CLOUDY_THRESHOLD = 80;

// Cloud Score+ per-pixel masking
var CS_PLUS = ee.ImageCollection('GOOGLE/CLOUD_SCORE_PLUS/V1/S2_HARMONIZED');
var CS_BAND = 'cs_cdf';       // 'cs' (instantaneous) or 'cs_cdf' (recommended)
var CLEAR_THRESHOLD = 0.60;   // keep pixels with cs_cdf >= 0.60

// Buffer size: 1 ha
var bufferAreaHa = 1.0;
var bufferAreaM2 = 10000;                               // 1 ha
var bufferRadiusM = Math.sqrt(bufferAreaM2 / Math.PI);  // ~56.4 m

var farms_raw = ee.FeatureCollection(CONFIG.asset);

print('EVI extraction | mode:', RUN_MODE, '| label:', CONFIG.bufferLabel);
print('Countries:', COUNTRIES.length, '| tasks to queue:',
    COUNTRIES.length * N_CHUNKS * years.length);
print('Footprint (ha):', bufferAreaHa, '| radius (m):', bufferRadiusM);
print('Cloud mask: Cloud Score+', CS_BAND, '>=', CLEAR_THRESHOLD);
print('Export folder:', CONFIG.exportFolder);

// =======================================================
// Country loop
// =======================================================

COUNTRIES.forEach(function (COUNTRY) {

    // ---- Filter to country ----
    var farms_country = farms_raw.filter(ee.Filter.eq(CONFIG.countryField, COUNTRY));

    // ---- De-dup near-coincident points (country-level, BEFORE chunking) ----
    // Done once per country (not per chunk) so it is global within a country.
    // Controls only: two random control points within ~111 m are effectively
    // the same location. Treated farms are distinct households and are NOT
    // de-duped (a coarse grid would drop legitimately close farms).
    if (CONFIG.dedupGrid) {
        farms_country = dedupByGrid(farms_country);
    }

    // ---- Assign deterministic chunks ----
    farms_country = farms_country.randomColumn('rand_chunk', 42);
    var farms_chunked = farms_country.map(function (f) {
        var chunk = ee.Number(f.get('rand_chunk')).multiply(N_CHUNKS).int();
        return f.set('chunk_id', chunk);
    });

    // ---- Chunk loop ----
    for (var chunkIdx = 0; chunkIdx < N_CHUNKS; chunkIdx++) {

        (function (chunk) {

            var farms_chunk = farms_chunked.filter(ee.Filter.eq('chunk_id', chunk));

            // 1-ha polygons around each farm
            var farms_buffer = farms_chunk.map(function (f) {
                var buffered = f.geometry().buffer(bufferRadiusM);
                return ee.Feature(buffered, f.toDictionary().set('buffer_area_ha', bufferAreaHa));
            });
            var farmsAOI = farms_buffer.geometry();

            years.forEach(function (year) {

                var s2 = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
                    .filterDate(year + '-01-01', (year + 1) + '-01-01')
                    .filterBounds(farmsAOI)
                    .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', CLOUDY_THRESHOLD))
                    .linkCollection(CS_PLUS, [CS_BAND]);   // attach Cloud Score+ band

                var eviCol = s2.map(maskCloudScorePlus).map(addEVI)
                    .select(['EVI', 'EVI_clamp']);

                // Annual reductions: raw and clamped
                var annualStack = eviCol.select('EVI').mean().rename('evi_mean')
                    .addBands(eviCol.select('EVI').median().rename('evi_median'))
                    .addBands(eviCol.select('EVI_clamp').mean().rename('evi_mean_clamp'))
                    .addBands(eviCol.select('EVI_clamp').median().rename('evi_median_clamp'));

                var stats = annualStack.reduceRegions({
                    collection: farms_buffer,
                    reducer: ee.Reducer.mean(),
                    scale: 10,
                    tileScale: 8
                });

                // Add metadata AND drop geometry (not exported to the CSV)
                stats = stats.map(function (f) {
                    return f.set({
                        'country_iso_trees': COUNTRY,
                        'year': year,
                        'buffer_area_ha': bufferAreaHa,
                        'chunk_id': chunk,
                        'n_chunks': N_CHUNKS
                    }).setGeometry(null);
                });

                Export.table.toDrive({
                    collection: stats,
                    description: COUNTRY + 'EVI_1ha_' + CONFIG.bufferLabel +
                        '_year_' + year + '_chunk_' + chunk,
                    fileFormat: 'CSV',
                    folder: CONFIG.exportFolder
                });
            });

        })(chunkIdx);
    }
});

// =======================================================
// (Visualization block removed - see previous version if a preview is
//  needed; a full country loop would stack too many map layers.)
// =======================================================

// =======================================================
// Helper Functions (declarations are hoisted; safe at the end)
// =======================================================

// Cloud Score+ per-pixel mask, then scale reflectance to 0-1.
// cs_cdf is 0-1 already, so it is read BEFORE the divide.
function maskCloudScorePlus(image) {
    return image.updateMask(image.select(CS_BAND).gte(CLEAR_THRESHOLD)).divide(10000);
}

// EVI = 2.5 * ((NIR - RED) / (NIR + 6*RED - 7.5*BLUE + 1)); bands B8/B4/B2.
// Adds both the raw EVI and a clamped-to-[-1,1] copy.
function addEVI(image) {
    var evi = image.expression(
        '2.5 * ((NIR - RED) / (NIR + 6 * RED - 7.5 * BLUE + 1))', {
        'NIR': image.select('B8'),
        'RED': image.select('B4'),
        'BLUE': image.select('B2')
    }).rename('EVI');
    var eviClamp = evi.clamp(-1, 1).rename('EVI_clamp');
    return image.addBands(evi).addBands(eviClamp);
}

// De-dup features onto a ~0.001 deg (~111 m) grid; keeps the first per cell.
function dedupByGrid(fc) {
    var withGrid = fc.map(function (f) {
        var coords = f.geometry().coordinates();
        var lon = ee.Number(coords.get(0));
        var lat = ee.Number(coords.get(1));
        var cellLon = lon.multiply(1000).round().divide(1000);
        var cellLat = lat.multiply(1000).round().divide(1000);
        var cellId = cellLat.format('%.3f').cat('_').cat(cellLon.format('%.3f'));
        return f.set('grid_id', cellId);
    });
    return withGrid.distinct(['grid_id']);
}
