Asteroid Orbit Example for CSc 620

Code Author:
    Kat Volk
    Postdoctoral Researcher
    The University of Arizona
    Department of Planetary Sciences/
    Lunar and Planetary Laboratory
    1629 E University Blvd
    Tucson, AZ 85721
    http://katvolk.com/wp/

Context
    A simulation code can inject particles in various locations in the
    solar system and do an in-depth simulation of that particles orbit.
    Scientists are trying to determine which locations lead to particles
    that are stably in mean motion resonance with Neptune.
    
    See the technique-summary.pdf for more information.

Checking Resonance

    The libration-check.pl file is a perl script that reads in information
    about various test particles and checks their orbital resonance.
    The output is some graphs in postscript format.
    
Running the code

    You will have to be on a system that has perl and gnuplot installed.
    
        chmod u+x libration-check.pl
        time ./libration-check.pl < inputs.txt


Resources

    
