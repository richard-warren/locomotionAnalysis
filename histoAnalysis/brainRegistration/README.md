# notes
- for all analyses, the **origin** of the coordinate systems is: top left corner of most rostral image
- one indexing is used for all coordinate system, such that the corner pixel is either (1,1,1) (for pixels), (2,2,2) (for microns), etc...
- 3D tensors have dimension order: AP, DV, ML

# histo todo
- tform plan
  - crop and resize images
  - find matrices that implement cropping and resizing, and save
  - find left and right matrices to perform subsequent steps


- [ ] make everything 25 micron resolution?
- [ ] try tform with coordinate args
- [ ] incorporate 'final section' into analysis
- [ ] apply tform to tracked cells
- [ ] figure out pipeline
- [ ] analyze all brains!

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
