# scripts-08-pipelineB-spm

## Files:
- `Step01_MainSPM.m`: script to run the SPM pipeline (1st level analysis, single and multi run)
- `Step02_ExtractTCfromROI.m`: script to extract the time course from a ROI based on spherical subject-specific coordinates and the activation map of each functional run.
- `Step03_CalculateERAs.m`: script to calculate the event-related averages (ERAs) from the time course of each ROI. Generates the figures.
- `Step04_ExtractTvalues.m`: script to extract the t-values from the activation map of each functional run based on spherical subject-specific coordinates.
