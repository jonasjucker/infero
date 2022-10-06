# Build
```bash
# infero
git clone git@github.com:jonasjucker/infero.git
cd infero
git checkout ecrad
. dev/euler_env.sh # or dev/daint_env.sh for daint
cd ..
mkdir infero_build 
cd infero_build
export INFERO_HOME=$(pwd); . ../infero/dev/env.sh
./../infero/dev/1_install_deps.sh
./../infero/dev/2_install_infero.sh
./../infero/dev/3_run_tests.sh

# ecrad ML
cd ../infero/ecrad
./euler_compile.sh # or daint_compile.sh for daint
```

# Run
```bash
cd infero/ecrad
ln -s /path/to/icon/grid icon_grid.nc
ln -s /path/to/input-data input_data.nc

# Euler
. ../dev/euler/env.sh
bsub -R "rusage[mem=30000]" ./ecrad_ml.sh path/to/regress-radiation/saved_models

# Daint
./ecrad_ml.sh path/to/regress-radiation/saved_models
