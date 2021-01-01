## Introduction

## External dependecies
- libigl: https://libigl.github.io/

## Compilation
Starting in this directory, issue:

    mkdir build
    cd build
    cmake .. -DCMAKE_BUILD_TYPE=Release

If you are using Mac or Linux, then issue:

    make

## Execution

Once built, you can execute the assignment from inside the `build/` using 

    ./cloth-simulation

While running, you can activate or de-activate the collision sphere by pressing `c`. 
Toggle wind effect using `w`.
`f` can be used to toggle fixing the bottom of the cloth.
`d` can be used to drop the cloth completely.

## Material parameters

`sterch_stffness` and `bend_stiffness` can be changed in `assignment_setup.h`

## Video

A video of the results can be found here: https://youtu.be/E0S3BfkH8nM


