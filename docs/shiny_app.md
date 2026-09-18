# Shiny application

Launch the application with:

```bash
bin/celltaminate app
```

The application allows users to upload Kraken2 or KrakenUniq reports, choose genus- or species-level analysis, adjust the prioritization threshold, inspect binary prioritized/not-prioritized results, and review abundance, score, taxonomy, and curated contextual information.

The bundled kitome/background and clinical-pathogen panels are loaded automatically. The reference background table is loaded automatically when `data/reference/refined_cell.lines.tsv` is present.
