# Regulation-of-Extramedullary-Erythropoiesis-by-CD47-and-THBS1 (CCBR-1072)

This code accompanies the paper entitled: Differential regulation by CD47 and thrombospondin-1 of extramedullary erythropoiesis in mouse spleen (In Press).

To reproduce these results, follow these steps:

1.  Clone this GitHub repo (i.e. the page you are on):
    * ```git clone https://github.com/NIDAP-Community/Regulation-of-Extramedullary-Erythropoiesis-by-CD47-and-THBS1.git```

2.  The input files for this pipeline will be available upon request. Please reach out to the authors before continue to following steps

3.  Install docker and build the docker container:
    * Navigate to the cloned repository directory. 
    * Move to the ./Docker_file/ directory of this repo

4.  Build the container:
    * ```docker build --tag nidap-r3 .```

5.  Navigate to the cloned repository directory, Run the conainer by mounting the ./src/ directory of the repo to /tmp/ in the container:
    * ```docker run -ti -v $(pwd)/src:/mnt nidap-r3```
    
6.  Run the following code.
    * ```cd /mnt```
    * ```bash run_pipeline.sh```
