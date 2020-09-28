# notes
- for all analyses, the **origin** of the coordinate systems is: top left corner of most rostral image
- one indexing is used for all coordinate system, such that the corner pixel is either (1,1,1) (for pixels), (2,2,2) (for microns), etc...
- 3D tensors have dimension order: AP, DV, ML

# histo todo
- function to plot cell over outlines
  - should cells be in pixels, or real units?
- turn pipeline into function
  - save
    - imgs: 3d original // 3d warped over ccf // 2d projection original // 2d projection over ccf
    - table: unit_id, orig_pix, orig_mm, ccf_pix, ccf_mm, session, ddi, shank
- incorporate 'final section' into analysis
- analyze all brains!

## files needed for reconstruction
- intrabrain per-cell coordinates
  - ephys/ephysHistoData/ephysHistoData.mat (single file with all brains)
  - Reconstruction/probe.mat: contains locations of good channels (don't use this)
- aligned histo tiff files
  - mouse_name/Neurotrace_final.tif
- tracing
  - Reconstruction/brain.mat: locations of ones for tracing for all brain regions

# qz notes
- `cer18` has all nuclei traced front to back - use this for testing
- section thickness is in the `histology` in the sessionInfo
- res for all images is 2um/pixel

# rick qs
- when will ephysHistoData_old be ported to ephysHistoData?: yes
- what is 'corner' of cellLocations, and are AP units correct?
- make sure 0 indexing no mistakes when scaling from pix to mm
- WE NEED TO DEFINE A CORNER!
