# Repository migration

The QSAR module was moved from `dokuyucubatu-beep/S.-aureus` on 2026-11-11.
All tracked files under `qsar/` were copied byte-for-byte; the model, data, split,
predictions and evidence manifests were not regenerated or changed.

Historical commit references in the module refer to the original repository:
- Protocol/data freeze: https://github.com/dokuyucubatu-beep/S.-aureus/commit/74bd292
- Completed QSAR release: https://github.com/dokuyucubatu-beep/S.-aureus/commit/1dabb40

The old repository retains its docking experiment. Removing the QSAR files from
its current branch does not erase historical commits.
