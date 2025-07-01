# FLAIR Test Plan 

## Objectives
- **Validate Accuracy:**  
  Ensure the FLAIR pipeline produces correct and high-confidence isoform results for a variety of RNA sequencing datasets. This includes verifying that known transcripts are identified and false positives are minimized.
- **Validate Performance:**  
  Verify that the pipeline runs efficiently within acceptable time and memory limits on different dataset sizes. Performance targets (runtime, memory) are established to support use in both development and production environments.
- **Regression Testing:**  
  Detect any regressions in results or performance when introducing new versions of FLAIR. The test suite should confirm that updates do not break existing functionality or expected outcomes.
- **Quality Control (QC):**  
  Provide confidence in the pipeline’s QC steps (e.g., read alignment quality, splice-site correction, isoform filtering) by using test cases that reflect both typical and edge-case scenarios.
- **Documentation and Reproducibility:**  
  Define a clear process for running tests and interpreting metrics, so developers and users can reproduce validation results and understand the quality of pipeline outputs.

---

## Scope
- **In Scope:**  
  - All major stages of the FLAIR pipeline (alignment, correction, collapse, quantification) are covered by test cases.  
  - Both accuracy (e.g. isoform identification correctness) and performance (runtime/memory) aspects are included.  
  - Covers multiple data types (PacBio Hi-Fi, simulated reads from LRGASP) and use cases (single-sample and multi-sample analysis).

- **Out of Scope:**  
  - Upstream data generation or wet-lab processes.  
  - Visualization tools or downstream analysis beyond FLAIR’s outputs.  
  - Internal correctness testing of third-party tools (e.g. minimap2), limited to their integration with FLAIR.

---

## Test Strategy
- **Test Types:**  
  Emphasis on integration testing using realistic datasets. Each test case runs the FLAIR CLI (or its modules) on input data, then validates output via automated comparisons to expected results or known benchmarks. Unit tests of individual functions are minimal—priority is on end-to-end and module-level tests.

- **Accuracy Verification:**  
  - For datasets with known ground truth (e.g. simulated data), use precision and recall of isoform identification.  
  - For real biological datasets, use expected outcomes (known isoform annotations, SQANTI classifications, TSS/TTS evidence) to qualitatively assess correctness.

- **Performance Testing:**  
  Track runtime and memory usage for each test. Include large-scale scenarios to ensure the pipeline meets performance criteria. Each test case will have pass/fail criteria (e.g. `< X` hours, `< Y` GB RAM).

- **Regression & Version Comparison:**  
  Run the suite in parallel against a baseline version (e.g. previous release). Flag differences in outputs (isoform counts, detected isoforms). Small differences may be acceptable if justified; significant deviations require review.

- **Test Automation:**  
  Automate tests where possible (e.g. CI pipelines). Provide scripts to run all cases, collect results, and save logs for auditing and comparison.

---

## Test Environment
- **Hardware/Platform:**  
  Testing in an HPC environment with sufficient resources.

- **Software Setup:**  
  - Install FLAIR via CLI (pip or Conda) at specific versions under test.  
  - All dependencies (Python, minimap2, samtools, etc.) in a controlled environment.  
  - Reference genome FASTA and annotation GTF are pre-indexed and available.

- **Data Availability:**  
  - Test inputs (FASTA/Q reads, reference files) stored in an accessible location.  
  - Some tests use small example data packaged with the suite; others require downloading public data (scripts or instructions provided).

- **Isolation:**  
  - Each test runs in its own directory to avoid interference.  
  - Temporary files cleaned up or stored per test.  
  - Cluster job scripts request consistent resources for comparable metrics.

- **Configuration:**  
  Default pipeline parameters are used unless testing a specific option. Ensures tests measure default behavior for reproducibility.

---

## Schedule and Execution
- **Development Phase:**  
  During new CLI development, add tests for new features or bug fixes. Developers run tests locally before merging.

- **Continuous Integration (CI):**  
  Run the suite (or core subset) on every pull request/merge to main. Immediate feedback on functionality or performance breaks.

- **Pre-release Validation:**  
  Before a new FLAIR release, execute the full suite in a clean environment. Address any failures or regressions prior to release.

- **Periodic Regression Runs:**  
  On a regular schedule (e.g. nightly/weekly), run critical accuracy and heavy-load performance tests against the latest development version.

- **Maintenance:**  
  When dependencies or reference data change, re-run relevant tests. A re-validation plan defines the scope needed for significant updates.

---

## Metrics and Criteria
- **Accuracy Metrics:**  
  - Simulated tests: require ≥ 95% isoform precision and recall.  
  - Biological tests: key isoforms or splicing events must be detected.  
  - Failure: metrics outside expected ranges or missing features.

- **Performance Metrics:**  
  - Defined runtime and memory thresholds per scenario (eg small tests in minutes; whole transcriptome in hours).  
  - Exceeding limits results in test failure.

- **Result Comparison:**  
  - Automated QC of output metrics and comparison to previous tests
  - Tolerance for minor differences (e.g. 1–2 extra novel isoforms).  
  - Significant discrepancies flagged.

- **Documentation of Outcomes:**  
  - Each test case documents expected outcomes and pass criteria.  
  - Reports indicate “Pass” or “Fail” for accuracy and performance.  
  - Failures trigger investigation and bug fixes before release.

---

## Roles and Responsibilities
- **Pipeline Developers:**  
  - Write and update tests for new features and bug fixes.  
  - Define expected outcomes and run the suite during development.  
  - Fix any regressions uncovered by tests.

- **QC/Testing Lead:**  
  - Own the test strategy and suite maintenance.  
  - Curate scenarios, review results, and update expected outputs.  
  - Decide when to add new datasets for coverage.

- **End Users (Optional UAT):**  
  - Experienced FLAIR users may run the suite as User Acceptance Testing.  
  - Provide feedback on real-world use cases.  
  - Not responsible for writing tests but help ensure practical coverage.

- **Release Manager:**  
  - Verify that all tests have passed before approving a new release.  
  - Coordinate with developers and QC lead to resolve issues.  
  - No release is approved until the suite gives a green light.
