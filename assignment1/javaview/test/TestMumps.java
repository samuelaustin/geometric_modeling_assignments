package test;

import dev6.numeric.PnMumpsSolver;
import jv.object.PsDebug;
import jv.vecmath.PdVector;
import jvx.numeric.PnSparseMatrix;


/**
* Tests whether MUMPS solver is working
*/
public class TestMumps {
	
	
	public static void main(String args[]) {
		if (!PnMumpsSolver.isAvailable()) {
			System.out.println("Mumps not available. Check VM arguments and dll paths.");
			return;
		}
		
		int n=10, k=5;
		PnSparseMatrix M = new PnSparseMatrix(n, n);
		for (int i=0; i<M.getNumRows(); i++) {
			//Big entry on diagonal
			M.setEntry(i, i, Math.random() * 5000 + 5000);
			//k small entries per row at random position, while keeping symmetry
			for (int count=0; count<k; count++) {
				int j = (int) (Math.random() * n);
				if (j != i) {
					M.setEntry(i, j, Math.random()*20 - 10);
					M.setEntry(j, i, Math.random()*20 - 10);
				}
			}
		}
		long startTime = System.currentTimeMillis();
		long factorization = -1;
		try {
			factorization = PnMumpsSolver.factor(M, PnMumpsSolver.Type.SYMMETRIC_POSITIVE_DEFINITE);
		} catch (Exception e) {
			System.out.println("Could not factorize matrix.");
			return;
		}
		System.out.println("Time for factorization: "+(System.currentTimeMillis()-startTime)+" ms");
		PdVector allOne = new PdVector(n);
		allOne.setConstant(1.d);
		PdVector solution = new PdVector(n);
		startTime = System.currentTimeMillis();
		try {
			PnMumpsSolver.solve(factorization, solution, allOne);
		} catch (Exception e) {
			System.out.println("Could not solve linear system.");
			return;
		}
		System.out.println("Time for one solve: "+(System.currentTimeMillis()-startTime)+" ms");
		
	}
}