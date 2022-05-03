# Abl1 Kinase aMD Scripts with PMFMe Calculator + Modi/Dunbrack Clustering for humans

This has improved, it's a simple Python archive of code with Java and Python processing scripts w/ two examples of protonated/deprotonated Asp381,
It wasn't too special from the start but Wonpil Im had some pretty good basic scripts you can have some fun with all your
force switching and all the other gear you might call a "proper CHARMM MD simulation" with PME and LJ6-12 and all that jazz.
Now it's quite better with some examples of boosting both PeriodicTorsionForce and CMAPTorsionForce as group 2 and a PMF script that preps up this
Dunbrack clustering assignment to feed a trajectory (like perhaps every 100 ps as implemented right now) and get a "time series".
It is nice running me on OpenCL with an NVIDIA RTX 2080 Ti or something cool like that, for some unknown reason it's no different than CUDA.

The modifications here do the aMD sampling over dihedrals (PeriodicTorsionForce and/or CMAPTorsionForce) with an AMDForceGroupLangevinVRORVIntegrator. 
Nothing fancy intended, just a complete Boltzmann sampling of x values for the configurational partition function at the given temperature.

^^^^^ Evolving storyline! ^_^

Running it - a restart XML file prod_0.rst (or prodD_0.rst, prodCMAP_0.rst) is provided for the jobs in production/openmm folder:
From the CHARMM-GUI OpenMM scripts, we added a special VRORV Integrator that also includes aMD protocol for a force group (e.g.)
dihedral boosting and saved all the collective variables (two distances, 10 dihedrals) every 10 steps (20 fs) in order to
average out noise from the PMF.  Does it work worse vs. with doing this every 1-2 steps?   It seems a little, but we only have a little
data to compare anymore to an old implementation with NAMD where it seems to factor out the variation of the "noise" in deltaV
at coordinate instances that determine the mean values and the PMF with a calculation:

< A > = <A * exp(beta * deltaV)> / <exp(beta * deltaV)>
  
PMF(coordiate) = -kT ln ( sum(exp(beta*deltaV)) )

- summation is done over all configurations falling into a certain PMF bin, now they don't have a weight of just "1" but ones that
  were lower energy got the biggest boost so that's why a large deltaV increases the weight, they are the one with lowest true "U_total,dih".    
- Depending on your parameterization, there may be no points where deltaV = 0 and thus the weight is the lowest 1 because that 
would have to be super high in energy space such that it was over the limit we put on the energy with the parameter E, so it gets no boost at all.  

Yes, you probably could just set this E to 150% of the <U_dihedral,total> from the simple MD run to never see this weight 1, but alpha should be large too, 
not too large though, one suggestion was 300 kcal/mol (1255.2 kJ/mol so that it is not totally smoothing it (flatlining U*_dihedral,total) and 
making the variations in deltaV even bigger (approaching the original variation in U_dihedral,total).   alpha = 0.2 (E - < V >) was suggested and the values
here are extremely close to this estimation formula from looking at results in the citation below.

So this is a essentially a very interesting example of histogram reweighting and can perhaps even get you 2D PMFs if you
are brave enough to get it running long enough and spin it up for 7 days over multiple runs out to a month or two.  Said who?  Well Andrew McCammon
published a lot on this with tons of references, he is this old guy who has done tons with many great scientists, who at one point worked upstairs from me
but I only really introduced myself once as "the guy who worked in Wilma Olson's lab" and I said she was doing well and left.
[Useful reference: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3115733/ it's more about alanine dipeptide than NAMD, without much description of their code]

But there is this concept of ideas propagating without direct contact (quantum mechanics or COVID-19?) so this had worked well on another system before, where 
it was fun to see how easy it is to dock drugs after you've massively "shook up" a crystal structure that is way too compact because of those
crystal packing forces (totally not "fake news" but YMMV!).  So you can throw away a 'hellovalot' (100 ns?) of aMD for equilibration and then just
start sampling away!  And someone else told me you need a drug docking expert for this, but a lot of blind algorithms just are killin' it automatically.


### PBS/SLURM/LSF Errata

Remember to change the LSF scripts if you are just taking this off the grid entirely to run freely on your desktop system with either
an NVIDIA card or an AMD card - factors like # of OpenCL/CUDA cores on the card and clock speeds and plenty of other things probably determine
the overall.  For that it's just a couple lines in the run_lilac.csh files
  
```csh

  set platform    = OpenCL       # one way to try it, vroom, will work on many hardware platforms
  set platform    = CUDA         # the other way to try it, can be just the same in most cases

  set cnt = 1                    # to resume with the prod_0.rst file

```


It is also machine-dependent so my observations on an LSF node where the OpenCL version on the NVIDIA RTX 2080 Ti matches CUDA
should be checked on other GPUs that seemed to run into issues with CUDA and other allocations - such as TeslaV100, in many instances performance is just equal.  
It's a great card but recently (2022) slightly out of date.
Maybe the RTX2080 is even better - these are all marketing tactics where other vendors have decent offerings too, you might be surprised to
find 3 separate OpenCL devices on your machine/laptop from vendors including all of Intel, AMD, and NVIDIA.  DON'T BELIEVE THE HYPE! :fire:

## Running Scripts

### Compiling PMFMe, and then running it

```bash
  javac PMFMe.java
  java PMFMe fileName.out
```

PMFMe outputs the JSON files (dihedrals.json and distances.json) for the following calculation.

### Running Modi/Dunbrack Clustering Assignments for time trajectories:

```bash
  python3 dunbrack_cluster.py dihedrals.json distances.json
```

### Simulation Protocols

Example Protocol:
- Energy minimization
- Short 125ps "protein equilibration" with constraints (400 bb/40 sc, in standard units kJ/(mol*nm^2) - your basic CHARMM-type protocol)
- 20 ns unconstrained "regular MD" (no constraints, Langevin 2fs, to calculate <U_dihedral,total> which is around 2675 kcal/mol or 11218 kJ/mol)
- "Infinite" aMD, throwing away, e.g., 100 ns in equilibration (get me away from that crystal!), for final statistics.
- Production data (4 simulations) was run with 500 ns aMD warmup and 3000 ns+ trial runs.  The pocket is too hydrophobic to open without a drug.

