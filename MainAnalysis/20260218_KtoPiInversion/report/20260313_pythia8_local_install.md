# Local PYTHIA8 install

- Exact external version used here: `PYTHIA 8.317`
- Upstream release tarball: `pythia8317.tgz`
- Official source URL: `https://pythia.org/download/pythia83/pythia8317.tgz`
- Upstream banner reported by the built binary:
  - `This is PYTHIA version 8.317`
  - `Last date of change: 20 Jan 2026`
- Tune used by the current local `e^+e^-` Z-pole test:
  - no explicit tune override was set in our wrapper or example
  - therefore the run uses the PYTHIA8 defaults
  - relevant `e^+e^-` default from the local XML documentation:
    - `Tune:ee = 7`
    - settings file: `tunes/Monash2013-ee.cmnd`
    - description: `Monash 2013 ee tune`
- For completeness, the global hadron-collision default in the same release is:
  - `Tune:pp = 14`
  - settings file: `tunes/Monash2013.cmnd`
  - this is not the controlling tune for the present `e^+e^-` run
- Install prefix:
  - `external/pythia8317-src/install`
- Local config helper:
  - `external/pythia8317-src/install/bin/pythia8-config`

## Verified local test

- Example source:
  - `tools/pythia8_zpole_example.cc`
- Wrapper script:
  - `tools/run_pythia8_zpole_example.sh`

Run:

```bash
tools/run_pythia8_zpole_example.sh 1000
```

This generates hadronic `e+e-` events at the `Z` pole with:

- `Beams:idA = -11`
- `Beams:idB = 11`
- `Beams:eCM = 91.2`
- `WeakSingleBoson:ffbar2gmZ = on`
- `23:onIfAny = 1 2 3 4 5`

## Quick validation

The local test completed successfully with a 200-event smoke run and reported:

- accepted events: `200`
- mean charged multiplicity: `20.68`
