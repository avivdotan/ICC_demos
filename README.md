# ICC_demos

TA's demos for the ["Introduction to Computation and Cognition" course](https://moodle2.bgu.ac.il/moodle/course/view.php?id=31934) @ BGU. 

All notebooks are written in [julia](https://julialang.org/) using [Pluto.jl](https://github.com/fonsp/Pluto.jl). 

## Run a notebook online

- Copy the notebook's URL and paste it [here](http://pluto-on-binder.glitch.me/). 
- Open the output link. 
- Go have a coffee, this will take some time. 
- Play with the notebook!

**Note:** Running Pluto notebooks via binder is ***slow*** (especially at loading). Some notebooks *might even fail to load* due to timeouts. Try running *locally* if possible. 


## Run a notebook locally

### Install [julia](https://julialang.org/) and [Pluto.jl](https://github.com/fonsp/Pluto.jl)

Follow the instructions [here](https://www.youtube.com/watch?v=OOjKEgbt8AI) or [here](https://github.com/fonsp/Pluto.jl#Installation). 

### Run a notebook

- Open a julia REPL ("command window") and type: 

    ```julia
    import Pluto; Pluto.run()
    ```

- Choose a notebook (either a URL or a local path). 
- Play with the notebook!

**Note:** Each notebook will automatically install any missing packages and compile them *on the fly*. This means The **first runs wil be slow**, especially the first run of the first notebook using *Plots.jl* (all notebooks here use it). 
