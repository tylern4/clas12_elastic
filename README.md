# clas12_elastic



Analysis program for elastic events in the clas12 detector.



##### To Build

```shell
git clone --recurse-submodules https://github.com/tylern4/clas12_elastic.git;
mkdir -p clas12_elastic/build;
cd clas12_elastic/build;
cmake ..;
make;
```



##### To Run

Program takes a beam energy and number of threads as environment variables. Input files are converted using dst2root as part of the [hipo_tools](https://github.com/JeffersonLab/hipo_tools) package.

```shell
NUM_THREADS=8 BEAM_E=10.6041 ./clas12_elastic output.root input_files_*.root
```
