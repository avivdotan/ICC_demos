# Introduction to Computation and Cognition (19.1.0006) - TA's Demos

TA's demos for the ["Introduction to Computation and Cognition" course](https://moodle2.bgu.ac.il/moodle/course/view.php?id=31934) @ BGU. 

All notebooks are written in [julia](https://julialang.org/) using [Pluto.jl](https://github.com/fonsp/Pluto.jl). \
It's best to run the notebooks using [Google Chrome](https://www.google.com/chrome/) (some plots might be rendered incorrectly in other browsers). 

## Run a notebook online

<!--  - Copy the notebook's URL and paste it [here](http://pluto-on-binder.glitch.me/). 
- Open the output link. -->
- Copy the notebook's name and paste it [here](https://five-orange-canopy.glitch.me/). 
- Click Run!. 
- Wait, this will take some time. 
- Play with the notebook!

**Note:** Running Pluto.jl notebooks via [Binder](https://mybinder.org/) is ***slow*** (especially at loading). Some notebooks *might even fail to load* due to timeouts. Try running *locally* if possible. 


## Run a notebook locally

***If you encounter any problems in the process, either run the notebooks online or google your error and deal with it yourself; we <ins>cannot</ins> provide any technical support.***

### Install [julia](https://julialang.org/) and [Pluto.jl](https://github.com/fonsp/Pluto.jl)

Follow the instructions [here](https://www.youtube.com/watch?v=OOjKEgbt8AI) or [here](https://github.com/fonsp/Pluto.jl#Installation). 

### Run a notebook

- Open a julia REPL ("command window") and type: 
    ```julia
    import Pluto; Pluto.run()
    ```
- Choose a notebook (either a URL or a local path). 
- Play with the notebook!

**Note:** Each notebook will automatically install any missing packages and compile them *on the fly*. This means The **first run of each notebook will be slow**, especially the first run of the first notebook using *[Plots.jl](http://docs.juliaplots.org/latest/)* (all notebooks here use it). 
