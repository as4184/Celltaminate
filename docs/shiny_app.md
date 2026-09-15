# Shiny application

The Shiny application provides an interactive interface for Celltaminate.

Launch it with:

```bash
bin/celltaminate app
```

The application allows users to:

- upload one or more Kraken2 or KrakenUniq reports;
- choose genus- or species-level analysis;
- adjust score thresholds;
- view sample-level and cohort-level summaries;
- inspect interactive taxon tables;
- review abundance and score visualizations;
- inspect taxonomy and curated contextual information for selected taxa.

The bundled kitome/background and clinical-pathogen panels are loaded automatically. The reference background table is loaded automatically when `data/reference/refined_cell.lines.tsv` is present.
