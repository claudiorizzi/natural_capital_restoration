// =====================================================================
// TREES control-farm selection - 50 km MAIN SPEC, LATEST version (2026-07-17)
// Google Earth Engine (JavaScript, code.earthengine.google.com)
//
// Purpose: regenerate the MAIN 1-50 km control-farm draw with the fixes
// reviewed on 2026-07-17. Single classic ring (no distance bands) - the
// spillover-gradient bands design lives in gee_control_farms_LATEST.txt.
//
// CHANGES vs the previous 50 km script:
//  (1) ROBUST SAMPLING: numPixels 1 -> 50 with an explicit seed; exported
//      point is a seeded-random pick among the candidates. numPixels:1
//      frequently returned EMPTY where valid land cover was sparse, so
//      seeds in sparse (browner) areas silently got no control -
//      systematic selection toward greener areas.
//  (2) EXCLUSION: controls must be >= 1 km from ALL treated farms of the
//      country (was ~3 x buffer radius, and only the control's own seed
//      farm had a 1 km minimum). Cross-country exclusion ignored for now.
//  (3) BUFFER CONSTANTS: 1 ha = 10,000 m2 (old script used 40,000 m2 =
//      4 ha while the comment claimed 1 ha / 56.4 m). Metadata here, but
//      CHECK the EVI-extraction script for the same 40,000 constant -
//      there it sets the control EVI footprint.
//  (4) LAND COVER: ESA WorldCover 2021 -> Dynamic World 2017 mode
//      composite. Screening on a 2021 snapshot conditioned selection on
//      post-baseline dynamics; 2017 is the panel baseline.
//  (5) EXPORTS: folder 'gee_control_farms_50km_LATEST' (kept separate
//      from the bands folder); control_method label built from the
//      actual parameters; 'no_control' seed farms exported too (type
//      column) so failed seeds are visible in the data.
//  (6) All countries loop in ONE script run (no more editing COUNTRY
//      between runs). Export tasks are still started MANUALLY from the
//      Tasks tab - process one country per session if preferred by
//      setting COUNTRIES to a single entry (e.g. ['KEN'] -> 8 tasks).
//
// INGESTION NOTE (Stata): exports contain BOTH type=="control" and
//   type=="no_control" rows. The join code must add:
//       keep if type=="control"
//   (keep the no_control rows for diagnostics - they are the failed seeds).
// =====================================================================

// =======================================================
// Parameters
// =======================================================

// Countries to process in this run
var COUNTRIES = ['KEN', 'SEN', 'TZA', 'UGA', 'MLI', 'GIN', 'CMR'];

// Chunking: process a fraction of farms per export to avoid memory issues.
var N_CHUNKS = 8;

// Asset with all treated farm coordinates
var farms_raw = ee.FeatureCollection(
    'projects/ee-claudiorizzi93/assets/farm_coordinates_all'
);

// Distance configuration: MAIN SPEC single ring, 1-50 km
var MIN_DIST_M = 1000;   // 1 km
var MAX_DIST_M = 50000;  // 50 km

// Minimum distance of a control from ANY treated farm of the country (m)
var EXCLUSION_M = 1000;

// Buffer size metadata: 1 ha (correct math; used by the EVI extraction step)
var bufferAreaHa = 1.0;
var bufferAreaM2 = bufferAreaHa * 10000;               // 1 ha = 10,000 m2
var bufferRadiusM = Math.sqrt(bufferAreaM2 / Math.PI); // ~56.4 m

// Sampling parameters (fix #1)
var N_CANDIDATES = 50;   // candidate pixels drawn per seed farm
var SAMPLE_SEED = 42;    // seed for the pixel sampler (reproducible)
var PICK_SEED = 7;       // seed for the random pick among candidates
var SAMPLE_SCALE = 30;   // sampling scale (m); DW native is 10 m

// Land Cover: Dynamic World 2017 mode composite (fix #4)
// Valid classes for controls: 2 = grass, 4 = crops, 5 = shrub_and_scrub
var VALID_CLASSES = [2, 4, 5];
var dw2017 = ee.ImageCollection('GOOGLE/DYNAMICWORLD/V1')
    .filterDate('2017-01-01', '2018-01-01')
    .select('label')
    .mode();

// control_method label built from the ACTUAL parameters (fix #5)
var methodLabel = 'ring_' + (MIN_DIST_M / 1000) + 'km_' +
    (MAX_DIST_M / 1000) + 'km_dw2017';

// Export folder (fix #5)
var EXPORT_FOLDER = 'gee_control_farms_50km_LATEST';

print('Generating 50 km control points');
print('Ring:', MIN_DIST_M, '-', MAX_DIST_M, 'm | exclusion:', EXCLUSION_M, 'm from all treated farms');
print('Tasks to be queued:', COUNTRIES.length * N_CHUNKS);

// =======================================================
// Loop: countries x chunks
// =======================================================

COUNTRIES.forEach(function (COUNTRY) {

    var farms_country = farms_raw.filter(ee.Filter.eq('iso_trees', COUNTRY));

    // Stable chunk assignment (seeded; a farm keeps its chunk across runs)
    var farms_chunked = farms_country.randomColumn('rand_chunk', 42)
        .map(function (f) {
            var chunk = ee.Number(f.get('rand_chunk')).multiply(N_CHUNKS).int();
            return f.set('chunk_id', chunk);
        });

    // Exclusion mask: >= EXCLUSION_M from ALL treated farms of this country (fix #2)
    var treatmentBuffer = farms_country.map(function (f) {
        return f.buffer(EXCLUSION_M);
    });
    var treatmentMask = ee.Image(0).paint(treatmentBuffer, 1);
    var notTreatment = treatmentMask.not();

    // Valid land cover AND outside the exclusion zone
    var lcMask = dw2017.eq(VALID_CLASSES[0])
        .or(dw2017.eq(VALID_CLASSES[1]))
        .or(dw2017.eq(VALID_CLASSES[2]))
        .and(notTreatment);
    var validLC = dw2017.updateMask(lcMask);

    for (var i = 0; i < N_CHUNKS; i++) {

        // Local scope to capture the loop variable 'i'
        (function (chunkIdx) {

            var farms_current_chunk = farms_chunked
                .filter(ee.Filter.eq('chunk_id', chunkIdx));

            var processFarm = function (farm) {
                farm = ee.Feature(farm);
                var geom = farm.geometry();

                // "Donut" geometry (annulus), 1-50 km
                var outer = geom.buffer(MAX_DIST_M, 100);
                var inner = geom.buffer(MIN_DIST_M, 100);
                var ring = outer.difference(inner, 100);

                // Robust sampling (fix #1): draw N_CANDIDATES pixels, then
                // take a seeded-random pick among them.
                var sample = validLC.sample({
                    region: ring,
                    scale: SAMPLE_SCALE,
                    numPixels: N_CANDIDATES,
                    seed: SAMPLE_SEED,
                    geometries: true,
                    dropNulls: true
                });

                var pick = ee.Feature(
                    sample.randomColumn('r', PICK_SEED).sort('r').first()
                );

                return ee.Algorithms.If(
                    sample.size().gt(0),
                    pick.set({
                        'original_farm_id': farm.get('forestgardenid'),
                        'country_iso_trees': COUNTRY,
                        'type': 'control',
                        'chunk_id': chunkIdx,
                        'buffer_area_ha': bufferAreaHa,
                        'control_method': methodLabel
                    }),
                    // failed seeds are EXPORTED now (fix #5), not dropped
                    farm.set({
                        'original_farm_id': farm.get('forestgardenid'),
                        'country_iso_trees': COUNTRY,
                        'type': 'no_control',
                        'chunk_id': chunkIdx,
                        'control_method': methodLabel
                    })
                );
            };

            var controls_all = farms_current_chunk.map(processFarm);

            // Export BOTH control and no_control rows (filter in Stata:
            // keep if type=="control").
            Export.table.toDrive({
                collection: controls_all,
                description: COUNTRY + '_ctrl_1km_50km_chunk_' + chunkIdx + '_LATEST',
                fileFormat: 'CSV',
                folder: EXPORT_FOLDER
            });

        })(i);
    }
});
