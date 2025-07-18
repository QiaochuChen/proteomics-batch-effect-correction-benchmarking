# proteomics-batch-effect-correction-benchmarking
# v1.0
This Docker image provides a pre-configured R environment for the proteomics batch effect correction benchmarking project. It includes all necessary system dependencies and R packages to ensure consistency and reproducibility in your analysis environment.

## Quick Start
To use this Docker image, please ensure you have Docker installed.

### Build the Docker Image:
Open your terminal in the root directory of your proteomics-batch-effect-correction-benchmarking project (i.e., the directory containing Dockerfile and renv.lock), then run the following command:

`docker build -t proteomics-batch-effect-correction-benchmarking .`

This command builds a Docker image named proteomics-batch-effect-correction-benchmarking based on the Dockerfile. The build process may take some time as it downloads the base image, installs system dependencies, and R packages.

### Run the Docker Container:
Once the build is complete, you can run the container using the following command:

`docker run --rm proteomics-batch-effect-correction-benchmarking:v1.0 test.R`


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

### Performance Optimization
Parallelism Limit: By default, the number of OpenBLAS threads in the container is limited to 8 (ENV OPENBLAS_NUM_THREADS=8) to optimize computational resource usage. You can adjust this value as needed.


### Project Structure
To ensure the image builds and runs correctly, your local project directory (proteomics-batch-effect-correction-benchmarking) should contain the following key files and directories:

proteomics-batch-effect-correction-benchmarking/
├── Dockerfile
├── renv.lock
├── .Rprofile
├── renv/
│   └── activate.R
├── utils/
│   └── gPCA_1.0.tar.gz  # Ensure this path is correct and corresponds to COPY ./utils/gPCA_1.0.tar.gz . in Dockerfile
├── test.R             # Your default R script
└── # Other R scripts, data, analysis files, etc...
