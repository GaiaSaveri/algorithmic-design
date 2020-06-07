## Homework for the course “Algorithmic Design”

**Author: Gaia Saveri**

**AA 2019/2020**

### Structure of the repository

---

This repository is structured as follows:

* the folder `Strassen` contains the solution to the first homework. The required codes are inside the sub-folder `strassen`, while the file `hw_12-03.pdf` contains a short report explaining and commenting the implemented codes.

* the folder `Heaps`  contains the solution to the second (sub-folder `homework1`) and the third homework (sub-folder `homework2`). The required codes are inside the (sub-)sub-folders `AD_binheaps`, while the files `hw_17-03.pdf` and `hw_24-03.pdf` contains a short report explaining and commenting the implemented codes, as well as the solution to theoretical exercises.

* the folder `Sorting` contains the solution to the fourth (sub-folder `homework1`) and the fifth homework (sub-folder `homework2`). The required codes are inside the (sub-)sub-folders `sorting`, while the files `hw_31-03.pdf` and `hw_2-04.pdf` contains a short report explaining and commenting the implemented codes, as well as the solution to theoretical exercises.

* the folder `Dijkstra` contains the solution to the sixth homework. The required codes are inside the subfolder `Ad_dijkstra`, while the file `hw_30-04.pdf` contains a small report commenting the implemented codes.

Inside each of the previously described folders there is a file named `request.pdf`, which is the original text of the given homework. When present, the sub-folder `data` contains the plots included in the report (all plots are realized using the package `matplotlib` of `Python3`).

### Compilation

---

In order to compile the codes one needs to produce a Makefile:

```
cmake -G "Unix Makefiles" CMakeLists.pdf
```

and then compile just typing `make`.

### Timings

---

All the timings have been taken in second. 
