# proteomics-batch-effect-correction-benchmarking
# v1.0
This Docker image provides a pre-configured R environment for the proteomics batch effect correction benchmarking project. It includes all necessary system dependencies and R packages to ensure consistency and reproducibility in your analysis environment.

## Quick Start
To use this Docker image, please ensure you have Docker installed.

### Pull the Docker Image:
run the following command (which should be completed within a few minutes):
`docker pull qiaochuchen/proteomics-batch-effect-correction-benchmarking:v2.0`

### Run the Docker Container:
run the following command:
`docker run --rm proteomics-batch-effect-correction-benchmarking:v2.0 test.R`

## Technical Details

### Base Environment
Base Image: rocker/r-ver:4.4.1 (R 4.4.1 environment based on Debian)

### Working Directory
The working directory inside the image is set to /proteomics-batch-effect-correction-benchmarking. All your project files will be copied to this directory.

### Dependencies
The image includes various system libraries required to run R packages, such as libssl-dev, libcurl4-openssl-dev, libcairo2-dev, etc. These dependencies have been optimized to resolve common R package compilation issues.

### R Package Management (renv)
The project uses renv for R package dependency management. The renv.lock file locks all R package versions, ensuring a highly reproducible environment. When you build the image, renv::restore() automatically installs the packages specified in renv.lock.

gPCA: The gPCA package is pre-downloaded locally and included in the image to address its potential installation issues.

fgsea: For fgsea package compilation issues, the C++14 standard is enforced (ENV CXX14FLAGS="-O2 -Wall -pedantic -std=c++14").

### System Compatibility
Our code leverages multi-core processing for performance optimization. We strongly recommend running this Docker container on a Linux or macOS environment, as the R parallelization functions (mclapply) rely on the fork mechanism, which is not supported natively on Windows.

If you are using Docker Desktop on Windows, you may experience significant performance degradation when running multi-core tasks, as the container's multi-core functionality is limited by the underlying Windows Subsystem for Linux (WSL).

### Performance Optimization
Parallelism Limit: By default, the number of OpenBLAS threads in the container is limited to 8 (ENV OPENBLAS_NUM_THREADS=8) to optimize computational resource usage. You can adjust this value as needed.
