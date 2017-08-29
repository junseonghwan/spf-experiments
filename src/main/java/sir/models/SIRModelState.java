package sir.models;

/**
 * Implements a simple state of SIR model to be used for SMC.
 * S(t), I(t), R(t): susceptible, infected, and recovered population
 * N_{SI}(t): count of the number of individuals who have transitioned from S to I
 * N_{IR}(t): count of the number of individuals who have transitioned from I to R
 * B: count of the number of birth
 * D_S(t), D_I(t), D_R(t): count of the deaths from S, I, R respectively
 * 
 * @author Seong-Hwan Jun (s2jun.uw@gmail.com)
 *
 */
public class SIRModelState 
{
	private int S, I, R; // the number of individuals in S, I, R
	private int NSI, NIR; // the count of individuals from S->I and I->R
	private int B, DS, DI, DR; // the count of birth, death from S, I, R
	
	public static SIRModelState initialize(int S0, int I0, int R0)
	{
		SIRModelState state = new SIRModelState();
		state.S = S0;
		state.I = I0;
		state.R = R0;

		state.NSI = 0;
		state.NIR = 0;
		
		state.B = 0;
		state.DS = 0;
		state.DI = 0;
		state.DR = 0;

		return state;
	}
	
	public static SIRModelState copyState(SIRModelState src)
	{
		SIRModelState state = new SIRModelState();
		state.S = src.S;
		state.I = src.I;
		state.R = src.R;

		state.NSI = src.NSI;
		state.NIR = src.NIR;
		
		state.B = src.B;
		state.DS = src.DS;
		state.DI = src.DI;
		state.DR = src.DR;

		return state;
	}

}
