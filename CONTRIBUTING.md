# Contributing

Contributions should be focused, reproducible, and easy to review.

## Development workflow

1. Create a new branch for each change.
2. Keep changes small and documented.
3. Add or update examples when changing user-facing behavior.
4. Do not commit patient data, protected health information, local filesystem paths, or large benchmark outputs.
5. Run a smoke test before opening a pull request.

## Code style

- Keep the command-line interface stable unless a change is necessary.
- Prefer explicit argument names and clear error messages.
- Preserve backward compatibility with Kraken2 6-column reports and Kraken2/KrakenUniq 8-column reports.
- Document new output files in `docs/output_files.md`.
