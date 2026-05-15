# Plot Final Results

This note documents the commands to reproduce the final K/pi plots from the current workflow.

## 1) Run full workflow (MC closure + nominal MC/data)

```bash
./analysis.sh
```

Outputs:
- `output/KtoPi-MC-Gen-Closure.root`
- `output/KtoPi-MC-Reco-Closure.root`
- `output/KtoPi-MC-Gen-Nominal.root`
- `output/KtoPi-MC-Reco-Nominal.root`
- `output/KtoPi-Data-Reco-Nominal.root`

## 2) Plot final K/pi comparisons

```bash
root -l -b -q plotClosure.C
root -l -b -q plotResultCor.C
root -l -b -q plotClosureVsData.C
```

Main figures:
- `KtoPiClosure_Overlay.pdf`
- `KtoPiRatio_Overlay.pdf`
- `KtoPiClosureVsData_Overlay.pdf`

## 3) Optional closure diagnostics by species

```bash
root -l -b -q plotClosureSpectra.C
root -l -b -q plotClosureByNtagSummary.C
```

Additional figures:
- `KPtClosure_Overlay.pdf`
- `PiPtClosure_Overlay.pdf`
- `PPtClosure_Overlay.pdf`
- `Closure_ByNtag_Summary.pdf`
