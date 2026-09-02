# Changelog

## v0.5.0

- Switched from `TimeseriesTools` to `TimeseriesBase` dependency.
- Dropped `Dierckx`; removed `stimulustrace` and `interpmatch` interpolation helpers from `Behaviour.jl`.
- Removed `Log𝑓` dim export (kept `Chan`, `Unit`, `Depth`).
- `Session(::AbstractSession)` identity constructor added.
- Fixed unit metrics for Visual Coding sessions.
- Bumped compat: `AllenSDK 0.4`, `DimensionalData 0.30`, `HDF5 0.17`, `Colors 0.13`, `GeometryBasics 0.5`, `JSON 1`, `julia 1.10`.

## v0.4.0

- Bumped `NWBStream` to 0.2; updated for newer Python packages (AllenSDK pipeline).

## v0.3.1

- Replaced legacy `Ti` dim with `𝑡` throughout.
- Windows fix for S3 data path parsing (`splitpath` → manual split).

## v0.3.0

- Migrated `SparseDimArray` → `SparseToolsArray` (built on `TimeseriesTools` 0.5 `ToolsArray`).
- Added custom dims `Chan`, `Unit`, `Depth`, `Log𝑓`.
- `AbstractDimArray` signatures updated to `AbstractToolsArray`.
- Removed in-repo docs build.
- LFP: layer calculation, depth method choices, brain-visualisation depth update (`inbrain` now generally 0).
- Visual Behavior updates incl. theta-score comparison.

## v0.2.0

- Moved cache/manifest storage from package dir to a `Scratch.jl` scratchspace.
- Performance improvements: preload cache, save depth data.
- Improved depth calculations with selectable methods.
- Offline-mode fixes for behavior cache.
- Moved streamline download to `__init__`.
- Filter noisy AllenSDK Python warnings on init.
- Fixed `loaddataframe` export name.
- Bumped `AllenSDK 0.2`, `DimensionalData 0.27`, `TimeseriesTools 0.4`.

## v0.1.0

- Initial release: extracted from `AllenNeuropixels.jl`.
- Session abstractions (`AbstractSession`, `Session`, `NWBSession`, `HybridSession`).
- Ecephys / Brain Observatory / MouseConnectivity / ReferenceSpace caches.
- LFP, spike band, behaviour, and visual-behaviour data access.
- Removed `PyCall` in favour of `PythonCall`.
