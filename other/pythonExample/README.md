# Working from Python Example


### :warning: ALPHA: Package under construction (pre-release) :warning:
Until I add tests, I make no guarantees as to the suitability of this code.

While working from R meets the workflow of many people if you have a large
amount of existing python code, or simply feel more productive working from
python, then using these network comparison tools can be difficult.

In the example file in this directory, we have included a simple way of calling the NetEMD code from python. 

It should be noted that it is not fast (as data has to be passed between the
runtimes), nor is it memory efficient as there are duplicate versions of data.
However, for networks with a relatively small of nodes, it is a simple and
relatively pain-free way of experimenting with this code. 

To run this code the following python packages are required:

* networkx - A python network library
* numpy -  Python matrix and array library
* rpy2 -  A library which allows R functions to be called from inside of python.

Note: If you are using python 2.7, then you cannot install the current version
of rpy2 as it is python 3 only, however the final version that is python 2.7
compatible will be fine.
