(ths repo was cloned from a previous [gist](https://gist.github.com/fsau/78f506ddc8434d4b896aaeb1c3e49f3b))

# Arbitrary Particle Simulator

## What is this?

A small pair of programs to simulate physical point-like particles and plot them. The number of dimensions, forces and particles is arbitrary, limited only by memory and processing requirements.
The `sim.c` is responsible for generating the data from the initial conditions and the `myplot.m` is responsible for normalizing the data in time (using splines) and plotting it in real-time.

## How to use?

The C code in this project comes with a 4-particle system. Taking that as example, it shouldn't be hard to adapt it to a different system. New forces can be easily created and added to a system. Note that you can add internal and external forces (every force is calculated with every pair, including itself, that it can be used as an external force).

The modifiers are _generalized forces_ which can alter directly the particle properties as speed, position, etc. That's useful for collisions and similar processes which are easily described analytically but aren't very easy to implement only with forces.

If a particular force or modifier requires some particle property (for example, charge), one can easily add it to the Particle struct and include it on the system.

The `myplot.m` is a Octave script that creates an animated real-time plot of the system, so it's important to keep the time scale on a visible range (not too small neither too large), although it can be easily adapted.

### Running

To run the simulator, just put every files on the same folder and run `octave -qf ./myplot.m` there. It will compile the code, run it and plot it automatically. To print GIF you need Imagemagick and Gifsicle installed on your system.
