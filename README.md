HOW TO:
- before compiling, specify number of variables n in metropolis_data.cpp
- compile using metropolis_compiler batch script
- run main_pipeline.py python script
-- in main_pipeline.py: specify number of communities, nodes per community and hyper-edge degree
-- script automatically generates inter- and intra-community connection probabilities p and q and checks for consistency
-- in addition, the hgsbm algorithm calculates the realized value of the mixing parameter mu

