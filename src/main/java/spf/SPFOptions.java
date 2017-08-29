package spf;

import java.util.Random;

import briefj.opt.Option;

public class SPFOptions {

    @Option 
    public boolean verbose = true;

    @Option(gloss = "Number of particles stored in memory.")
    public int numberOfConcreteParticles = DEFAULT_N_CONCRETE_PARTICLES;

    @Option(gloss = "Maximum number of implicit particles, represented via random number generation replay. Only costs CPU, not memory.")
    public int maxNumberOfVirtualParticles = 1000;

    @Option(gloss = "Virtual particles will be used until that relative effective sampling size is reached (or maxNumberOfVirtualParticles is reached)")
    public double targetedRelativeESS = 0.5;
    
    @Option
    public Random resamplingRandom = new Random(1);

    @Option
    public Random mainRandom = new Random(1);

    @Option
    public ResamplingScheme resamplingScheme = ResamplingScheme.MULTINOMIAL;
    
    @Option
    public boolean storeParticleWeights = true;

    public static final int DEFAULT_N_CONCRETE_PARTICLES = 1000;
    
}
