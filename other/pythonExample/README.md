# Working from Python Example


### :warning: ALPHA: Package under construction (pre-release) :warning:
Until I add tests, I make no guarantees as to the suitability of this code.

While working from R meets the workflow of many people if you have a large
amount of existing python code, or simply feel more productive working from
python, then using these network comparison tools can be difficult.

In the example file in this directory we have included a simple way of calling the NetEMD code from python. 

It should be noted that it is not fast (as data has to be passed between the
runtimes), nor is it memory efficient as there is duplicate versions of data.
However, for networks with relatively small numbers of edges it is a simple and
relatively painfree way of experimenting with this code. 

