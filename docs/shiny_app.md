# Shiny app

The Shiny app provides a point-and-click interface for users who do not want to run command-line analysis directly.

Launch it with:

```bash
bin/celltaminate app
```

The app allows users to:

- upload one or more Kraken2 or KrakenUniq reports;
- choose genus-level or species-level scoring;
- adjust score thresholds and scoring strictness;
- view cohort-level and sample-level summaries;
- inspect interactive taxon tables;
- review abundance plots and score distributions;
- select taxa for contextual summaries with taxonomy and curated remarks when available.

For deployment, the app can be run on a local workstation, Shiny Server, Posit Connect, or a containerized environment.
