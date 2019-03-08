Guide to the data structures ‘massdata’ and ‘redmaela’

The data and problem setup are in matlab structure formats, and have basically the following fields:
    
    finalSC - sequence cluster ids (this is just 0:40)
    finalfreqSC - other final sequence cluster frequencies (not needed)
    G - matrix encoding the genome data. Rows are isolates, cols are loci. 1 if locus present, 0 if absent
    taxon - names of the rows (isolates, taxa)
    serotype - serotype of each row
    sc - sequence cluster of each row
    vt - 0 or 1; whether the type was in some particular vaccine (not needed) 
    dev- ‘deviation’ numbers (not needed)
    ics - initial conditions
    seronames - names of the serotypes in the dataset, order is important. 
    locusids - names of loci
    SerotypeToStrain - matrix encoding which serotype each strain is
    SeroNames - redundant
    strainInv - old invasiveness numbers
    sigma - parameter for strength of NFDS
    antigenvals - something from Nick about each antige n’s effectiveness (or 0 if not an antigen) 
    AntiNames - names of the antigens; order is important
    AllNames - not used
    AntiINDEX - locations of the antigens in the order they are listed in AntiNames
    weighteddev - info from weights of the NFDS effect by loci
    locusweights - locus weights 
    locusfreq - equilibrium locus frequencies
    DR_Inv_combo - old combined score
    DRscore - DR score of each genotype
    Invasiveness - structure with invasiveness numbers from the meta-nalaysis
    note - note about the new invasiveness 
    v - vaccine efficacy parameter
    m - migration (small parameter)


File list:

example_script.m  -- script showing how to set up the model and run the Bayesian optimization in matlab

fitness_for_bo.m -- takes the Bayesian optimization variable and evaluates the fitness. This function converts some formats and then calls the ODE solver (runODEmodel.m)

runODEmodel.m -- takes in a number of problem-specific parameters and runs the ODE. It returns either just the fitness or a structure with the solution to the ODE. This calls getInvasiveness at the end 

getInvErrorbars.m -- info about how I got invasiveness error bars (not needed to run the code)

getInvasiveness.m -- get the invasivness of a simulated population

makeallvars.m -- creates the 'optimizable variables' needed by bayesopt

myxconstraint -- constraint function allowing me to constrain the number of serotypes in the Bayesian optimisation

