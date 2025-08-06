
# bugs/todos

- [ ] **unable to set reconstruction reference by clicking**
- [ ] **resume on non-critical exceptions**
- [ ] **catch empty intervals, but only if needed**, else ignore
- [ ] **check if items are cyclic** - dry run option
- [ ] **test orientation of images of different formats, reconstructed or not**
- [ ] **stitch multiple input files**
- [ ] **compare reconstruction to pre-reconstructed files**
- [ ] when opening files after one another, intervals/flags are not displayed correctly
- [ ] changing detectors on empty intervals causes errors because detect is triggered
- [ ] catch "item does not exist"
- [ ] validate jsons
- [ ] Exception when first loading reconstructed, then raw


# features to come

- [ ] preprocessors
    - [ ] resample (time)
    - [ ] average breaths
    - [ ] diff (convert to flow)
    - [ ] rotate
    - [ ] normalize
    - [ ] filter (over time)
    - [ ] resample time
- [ ] base items
    - [ ] FRIC
    - [ ] optical flow
    - [ ] phase delay
    - [ ] min/max delay
    - [ ] correlation with sum signal
    - [ ] point-wise time constants?
    - [ ] keep dimensions-option for EELI
    - [ ] all thorax ROI
    - [ ] all thorax ROI hardcoded
- [ ] operations
    - [ ] ApEn (approximate entropy)
    - [ ] correlation
    - [ ] flatten dimension
    - [ ] trim last (e.g. for breath_times)
    - [ ] custom non-time regressions
    - [ ] custom non-time fits?
    - [ ] normalize std 1 mean 0.
- [ ] detector
    - [ ] manual flag placement/manipulation, keep flag positions when switching to manual
    - [ ] restrict detector selection to cyclic/maneuver detectors
- [ ] data handler
    - [ ] different reconstruction algorithms, generic representation for GUI
    - [ ] sentec raw support
    - [ ] timple reconstructed support
    - [ ] timple raw support
- [ ] backend
    - [ ] calculation multithreading
- [ ] GUI
    - [ ] drag-and-drop open
    - [ ] change plot options black-on-white for screenshots
    - [ ] **item preview** graph view / item editor
    - [ ] automatic all-file-interval
- [ ] exporter
    - [ ] export pdf report
    - [ ] export 2d as npy, 2d as graphic
- [ ] extras: visualization, interpretation, aggregation; FV-loops histograms https://pftforum.com/review/tutorials/spirometry-tutorials/assessing-flow-volume-loops/
