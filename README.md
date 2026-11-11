# Dissertation data compendium (public release v1.0)

Open data of the PhD dissertation on sand and gravel extraction on the rivers of northern and northeast Bangladesh, by MD Sadiul Alam Chyon (Department of Geography and GIS, University of Illinois Urbana-Champaign; supervisor James L. Best).

The thesis reconstructs where and when sand and gravel extraction developed along the northern fringe of Bangladesh since the late 1980s (Chapter 4), dates extraction and recovery on the Lubha River (Chapter 5), traces the impact chain of sand mining across channel, aquifer and floodplain on the Someshwari River (Chapter 6), tests whether the same chain recurs on the Jadukata, Piyain, Sari-Gowain and Lubha rivers (Chapter 7), and examines the governance of extraction on the Someshwari (Chapter 8).

## Contents

| File or folder | Contents |
| --- | --- |
| `Dissertation_Data_Guide_public_v1.0.docx` | Description of every group of data, the stations and sections, imagery, field campaigns, documentary sources and withheld data |
| `Dissertation_Data_Catalog_public_v1.0.xlsx` | One row per dataset (62), stations and record coverage, cross-sections, imagery and scenes, field campaigns, press and lease sources, code and archives, withheld data, MD5 manifest of every file |
| `Dissertation_Data_Tables_public_v1.0.xlsx` | The open datasets as tidy tables |
| `01_hydrometric_records_BWDB_BMD/` | Station metadata and coordinates; data request lists |
| `02_satellite_imagery_and_DEM/` | Scene lists, Google Earth Engine scripts (Chapters 4, 6, 7), Copernicus DEM GLO-30 tiles |
| `03_extraction_mapping_and_derived_results/` | Derived results by chapter: Ch2 hydroclimatic baseline; Ch4 site set (337 sites), January Landsat series, footprint, inventory, training and accuracy samples, screens; Ch5 Lubha pits, volumes, barracks and 2020 pit maps; Ch6 Someshwari bar area; Ch7 bar area, section metrics, stage at the same discharge, wells, Sentinel-1 inundation; Ch9 synthesis series |
| `05_field_photos_and_drone/` | Indices of field and drone photographs (time, GPS position) |
| `06_press_reports_and_official_documents/` | Bangladesh Sand Watch records (555 at 10 November 2026; 489 at 2 November 2026, the count used in the thesis), press and official sources cited in the chapters, governance events |
| `07_leases_gazettes_land_records/` | District balumahal lease lists and tenders, Bangladesh Gazette quarry areas (2013, 2015, 2016), national land-use map mining class, mouza gazetteer, original documents |
| `08_GIS_layers/` | Areas of interest, Chapter 4 site windows and inventory, rivers and boundaries |
| `09_national_indicators_and_literature/` | National demand series of Chapter 4; values digitized for Figure 1.1a |

## Withheld from this release

- Raw hydrometric records of the Bangladesh Water Development Board (BWDB) and the Bangladesh Meteorological Department (BMD): available from the agencies on request. The catalog lists every station and section with its period and number of records.
- Identifiable field data (transcripts, field forms, files with personal details), protected under UIUC IRB protocol IRB23-0168.
- The de-identified Someshwari household survey and focus group and key informant aggregates, to be added after acceptance of the related journal article.
- Licensed imagery (PlanetScope, SkySat, TripleSat, SPOT, WorldView); scenes are listed in the catalog.
- Field photographs and drone media; third-party series with restrictive terms.

## Conventions

Water level and bed elevation in m PWD (Public Works Datum, c. 0.46 m below mean sea level); discharge in m3 s-1; rainfall in mm; groundwater as depth below ground (m); areas in ha or km2. Chapter numbers follow the final structure of the thesis: 1 Introduction; 2 Study region; 3 Research design and methods; 4 Northern fringe; 5 Lubha River; 6 Someshwari River; 7 Four further rivers; 8 Governance (social); 9 Synthesis.

## Related work

- Chyon, M. S. A., Best, J., & Haq, S. M. A. (2026). From extraction to recovery: Geomorphic and social transformations of gravel mining, Lubha River, Northeast Bangladesh. *Cell Reports Sustainability*, 3(1), 100605. https://doi.org/10.1016/j.crsus.2025.100605
- Lubha River code: sadiulchyon/gravel-mining-in-sylhet v3 (2025). https://doi.org/10.5281/zenodo.17692625
- Bangladesh Sand Watch archive: https://github.com/sadiulchyon/bangladesh-sand-watch
- Northern fringe extraction atlas: https://sadiulchyon.github.io/northern-fringe-extraction-atlas/

## Licence

The author's data and documents are released under the Creative Commons Attribution 4.0 International licence (CC BY 4.0; see `LICENSE`). Third-party files keep their own terms: Copernicus DEM tiles (Copernicus licence, attribution required), HydroRIVERS clip (CC BY 4.0; Lehner and Grill, 2013), ISCGM Global Map waterways (attribution), and government lease and gazette documents (public records). Press records give bibliographic details, the author's summaries and links; article texts belong to the outlets.

## How to cite

Please cite the Zenodo record of this release (see `CITATION.cff`).
