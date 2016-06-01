package workshop;

import java.awt.Color;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Random;
import java.util.Vector;

import jv.geom.PgBndPolygon;
import jv.geom.PgElementSet;
import jv.geom.PgPolygonSet;
import jv.geom.PgVectorField;
import jv.geom.PuCleanMesh;
import jv.number.PdColor;
import jv.object.PsConfig;
import jv.object.PsDebug;
import jv.object.PsObject;
import jv.project.PgGeometry;
import jv.vecmath.PdMatrix;
import jv.vecmath.PdVector;
import jv.vecmath.PiVector;
import jv.vecmath.PuMath;
import jv.viewer.PvDisplay;
import jv.project.PvGeometryIf;
import jvx.numeric.PnMatrix;
import jvx.numeric.PnSparseMatrix;
import jvx.project.PjWorkshop;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.LUDecomposition;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.SingularValueDecomposition;

import dev6.numeric.PnMumpsSolver;
/**
 *  Workshop for surface registration
 */

public class Registration extends PjWorkshop {
	
	/** First surface to be registered. */	
	PgElementSet	m_surfP;	
	/** Second surface to be registered. */
	PgElementSet	m_surfQ;	
	
	
	/** Constructor */
	public Registration() {
		super("Surface Registration");
		if (getClass() == Registration.class) {
			init();
		}
	}
	
	/** Initialization */
	public void init() {
		super.init();
	}
	
	
	/** Set two Geometries. */
	public void setGeometries(PgElementSet surfP, PgElementSet surfQ) {
		m_surfP = surfP;
		m_surfQ = surfQ;
	}
	
	public void iterativeClosestPoint(boolean pointtoplane){
		if(m_surfP == null || m_surfQ == null)
		{
			System.out.println("Surfaces have not been set.");
			return;
		}
		
		// Choose random amount of vertices from P.
		int n = Math.min(m_surfP.getNumVertices(), 1000);
		Random r = new Random(System.currentTimeMillis());
		PdVector[] vectors;
		boolean[] used = new boolean[m_surfP.getNumVertices()];
		ArrayList<PdVector> allVectors = new ArrayList<PdVector>(Arrays.asList(m_surfP.getVertices()));
		if(n == m_surfP.getNumVertices()){
			vectors = m_surfP.getVertices();
		} else {
			vectors = new PdVector[n];
			for(int i = 0; i < n; i++){
				int next = r.nextInt(allVectors.size());
				vectors[i] = allVectors.remove(next);
			}
		}
		
		// Calculate min. distance between this set and points of Q.
		
		if(pointtoplane){
			HashMap<PdVector, HashMap<PdVector, Integer>> corresponding = GetClosestPointToPlaneDistances(vectors, m_surfQ);
			pointToPlane(corresponding, m_surfQ);
		}
		else {
			// Calculate min. distance between this set and points of Q.
			HashMap<PdVector, PdVector> corresponding = GetClosestPointToPointDistances(vectors, m_surfQ);
			pointToPoint(corresponding);
		}
		/*
		// Compute centroid of P and Q.
		PdVector pCentroid = computeCentroid(corresponding.keySet());
		PdVector qCentroid = computeCentroid(corresponding.values());
		
		// Compute optimal rotation and translation.
		PdMatrix m = new PdMatrix(3);
		for(PdVector vector : corresponding.keySet()){
			PdMatrix a = new PdMatrix(3);
			a.adjoint(PdVector.subNew(vector, pCentroid), PdVector.subNew(corresponding.get(vector), qCentroid));
			m.add(a);
		}
		double n2 = corresponding.keySet().size();
		PdMatrix rotation = divide(m, n2);
		
		// Compute SVD.
		Array2DRowRealMatrix realmatrix = new Array2DRowRealMatrix(rotation.m_data);
		SingularValueDecomposition svd = new SingularValueDecomposition(realmatrix);
		
		double[][] inter = ones(3);
		
		RealMatrix uvt = svd.getU().multiply(svd.getVT());
		PdMatrix vut = new PdMatrix(svd.getV().multiply(svd.getUT()).getData());
		inter[2][2] = PnMatrix.determinant(vut.m_data, 3);
		RealMatrix rOptInter = svd.getV().multiply(new Array2DRowRealMatrix(inter).multiply(svd.getUT()));
		PdMatrix rOpt = new PdMatrix(rOptInter.getData());
		PdVector tOpt = PdVector.subNew(qCentroid, matrixMult(rOpt, pCentroid));
		
		rotateAndTranslate(m_surfP, rOpt, tOpt);*/
	}
	
	public void pointToPlane(HashMap<PdVector, HashMap<PdVector, Integer>> corresponding, PgElementSet surface) {
		// Get all vertex normals.
		surface.makeVertexNormals();
		PdVector[] normalsOfQ = surface.getVertexNormals();
		
		// Construct A.
		RealMatrix A = constructA(corresponding, normalsOfQ);
		RealMatrix AT = A.transpose();
		
		// Calculate (A^T*A)^-1.
		RealMatrix ATA = AT.multiply(A);
		RealMatrix ATAinv = new LUDecomposition(ATA).getSolver().getInverse();
		
		// Construct b.
		RealMatrix b = constructB(corresponding, normalsOfQ);
		
		// Calculate x.
		RealMatrix x = ATAinv.multiply(AT).multiply(b);
		
		// Calculate R and t.
		double r1 = x.getEntry(0, 0);
		double r2 = x.getEntry(1, 0);
		double r3 = x.getEntry(2, 0);
		PdMatrix R = new PdMatrix(3);
		R.setEntry(1, 0, r3);
		R.setEntry(2, 0, -r2);
		R.setEntry(0, 1, -r3);
		R.setEntry(2, 1, r1);
		R.setEntry(0, 2, r2);
		R.setEntry(1, 2, -r1);
		
		// Compute SVD.
		Array2DRowRealMatrix realmatrix = new Array2DRowRealMatrix(R.m_data);
		SingularValueDecomposition svd = new SingularValueDecomposition(realmatrix);
				
		double[][] inter = ones(3);

		PdMatrix uvt = new PdMatrix(svd.getU().multiply(svd.getVT()).getData());
		inter[2][2] = PnMatrix.determinant(uvt.m_data, 3);
		RealMatrix rOptInter = svd.getU().multiply(new Array2DRowRealMatrix(inter).multiply(svd.getVT()));
		
		R = new PdMatrix(rOptInter.getData());
		
		double[] t_entries = {x.getEntry(3, 0), x.getEntry(4, 0), x.getEntry(5, 0)};
		PdVector t = new PdVector(t_entries);
		
		// Update the mesh.
		rotateAndTranslate(m_surfP, R, t);
	}
	
	public void pointToPoint(HashMap<PdVector, PdVector> corresponding) {
		// Compute centroid of P and Q.
		PdVector pCentroid = computeCentroid(corresponding.keySet());
		PdVector qCentroid = computeCentroid(corresponding.values());
				
		// Compute optimal rotation and translation.
		PdMatrix m = new PdMatrix(3);
		for(PdVector vector : corresponding.keySet()){
			PdMatrix a = new PdMatrix(3);
			a.adjoint(PdVector.subNew(vector, pCentroid), PdVector.subNew(corresponding.get(vector), qCentroid));
			m.add(a);
		}
		double n2 = corresponding.keySet().size();
		PdMatrix rotation = divide(m, n2);
		
		// Compute SVD.
		Array2DRowRealMatrix realmatrix = new Array2DRowRealMatrix(rotation.m_data);
		SingularValueDecomposition svd = new SingularValueDecomposition(realmatrix);
				
		double[][] inter = ones(3);
				
		RealMatrix uvt = svd.getU().multiply(svd.getVT());
		PdMatrix vut = new PdMatrix(svd.getV().multiply(svd.getUT()).getData());
		inter[2][2] = PnMatrix.determinant(vut.m_data, 3);
		RealMatrix rOptInter = svd.getV().multiply(new Array2DRowRealMatrix(inter).multiply(svd.getUT()));
		PdMatrix rOpt = new PdMatrix(rOptInter.getData());
		PdVector tOpt = PdVector.subNew(qCentroid, matrixMult(rOpt, pCentroid));
				
		rotateAndTranslate(m_surfP, rOpt, tOpt);
	}
	
	private PdVector computeCentroid(Collection<PdVector> vectors){
		PdVector centroid = new PdVector(0.0, 0.0, 0.0);
		for(PdVector v : vectors){
			centroid.add(v);
		}
		double divisor = vectors.size();
		double[] data = centroid.m_data;
		data[0] = data[0]/divisor;
		data[1] = data[1]/divisor;
		data[2] = data[2]/divisor;
		return new PdVector(data);
	}

	private HashMap<PdVector, PdVector> GetClosestPointToPointDistances(PdVector[] vectors, PgElementSet surface)
	{
		// Calculate min. distance between this set and points of Q.
		HashMap<PdVector, Double> results = new HashMap<PdVector, Double>();
		HashMap<PdVector, PdVector> corresponding = new HashMap<PdVector, PdVector>();
		double meanDist = 0.0;
		for(int j = 0; j < vectors.length; j++){
			double dist = Double.MAX_VALUE;
			PdVector opt = new PdVector();
			PdVector[] vectorsOfQ = surface.getVertices();
			for(int j2 = 0; j2 < vectorsOfQ.length; j2++){
				if(PdVector.dist(vectors[j], vectorsOfQ[j2]) < dist){
					dist = PdVector.dist(vectors[j], vectorsOfQ[j2]);
					opt = vectorsOfQ[j2];
				}
			}
			meanDist += dist;
			results.put(vectors[j], dist);
			corresponding.put(vectors[j], opt);
		}
		meanDist = meanDist/vectors.length;
		
		// Remove all vectors with distance greater than k*meanDist.
		double k = 2.0;
		for(PdVector v : results.keySet()) {
			if(results.get(v) > k*meanDist){
				corresponding.remove(v);
			}
		}

		return corresponding;
	}

	private HashMap<PdVector, HashMap<PdVector, Integer>> GetClosestPointToPlaneDistances(PdVector[] vectors, PgElementSet surface)
	{
		// Calculate min. distance between this set and points of Q.
		HashMap<PdVector, Double> results = new HashMap<PdVector, Double>();
		HashMap<PdVector, HashMap<PdVector, Integer>> corresponding = new HashMap<PdVector, HashMap<PdVector, Integer>>();
		double meanDist = 0.0;

		PdVector[] vectorsOfQ = surface.getVertices();
		//surface.makeVertexNormals();
		//PdVector[] normalsOfQ = surface.getVertexNormals();

		for(int j = 0; j < vectors.length; j++){
			double minDist = Double.MAX_VALUE;
			PdVector closestPoint = null;
			int indexOfQ = -1;
			for(int j2 = 0; j2 < vectorsOfQ.length; j2++){
				//PdVector point = PdVector.subNew(vectors[j], vectorsOfQ[j2]);
				//PdVector projectionOntoNormal = GetProjectionPOntoQ(point, normalsOfQ[j2]);
				double dist = PdVector.dist(vectors[j], vectorsOfQ[j2]);
				if(dist < minDist){
					minDist = dist;
					closestPoint = vectorsOfQ[j2];
					indexOfQ = j2;
				}
			}
			meanDist += minDist;
			results.put(vectors[j], minDist);
			HashMap<PdVector, Integer> map = new HashMap<PdVector, Integer>();
			map.put(closestPoint, indexOfQ);
			corresponding.put(vectors[j], map);
		}
		meanDist = meanDist/vectors.length;

		double k = 5.0;
		for(PdVector v : results.keySet()) {
			if(results.get(v) > k*meanDist){
				corresponding.remove(v);
			}
		}

		return corresponding;
	}

	private PdMatrix divide(PdMatrix input, double divisor){
		double[][] data = input.m_data;
		for(int i = 0; i < data.length; i++){
			for(int j = 0; j < data[0].length; j++){
				data[i][j] = data[i][j]/divisor;
			}
		}
		return new PdMatrix(data);
	}
	
	private double[][] ones(int dim){
		double[][] res = new double[dim][dim];
		for(int i = 0; i < dim; i++){
			for(int j = 0; j < dim; j++){
				if(i == j){
					res[i][j] = 1.0;
				}
				else {
					res[i][j] = 0.0;
				}
			}
		}
		return res;
	}
	
	private PdVector matrixMult(PdMatrix m, PdVector v){
		double[] results = new double[3];
		for(int i = 0; i < 3; i++){
			results[i] = (m.getEntry(0,i)*v.getEntry(0)) + (m.getEntry(1,i)*v.getEntry(1)) + (m.getEntry(2,i)*v.getEntry(2));
		}
		return new PdVector(results);
	}
	
	/*private PdVector GetProjectionPOntoQ(PdVector p, PdVector q)
	{
		PdVector result = new PdVector(q.m_data);
		result.multScalar(p.dot(q) / q.dot(q));
		return result;
	}*/
	
	private void rotateAndTranslate(PgElementSet set, PdMatrix r, PdVector t)
	{
		System.out.println("Rotating and translating");
		System.out.println(r);
		System.out.println(t);
		for(int i = 0; i < set.getNumVertices(); i++){
			PdVector old = set.getVertex(i);
			set.setVertex(i, PdVector.addNew(matrixMult(r, old), t));
		}
		m_surfP.update(m_surfP);
	}
	
	private RealMatrix constructA(HashMap<PdVector, HashMap<PdVector, Integer>> corresponding, PdVector[] normals){
		// Constructing A as shown here: http://resources.mpi-inf.mpg.de/deformableShapeMatching/EG2011_Tutorial/slides/2.1%20Rigid%20ICP.pdfS
		int n = corresponding.keySet().size();
		RealMatrix res = new Array2DRowRealMatrix(n, 6);
		int i = 0;
		for(PdVector pi : corresponding.keySet()) {
			HashMap<PdVector, Integer> results = corresponding.get(pi);
			PdVector qi = results.keySet().iterator().next();
			int indexOfNorm = results.get(qi);
			PdVector crossRes = PdVector.crossNew(pi, normals[indexOfNorm]);
			res.setEntry(i,0, crossRes.getEntry(0));
			res.setEntry(i,1, crossRes.getEntry(1));
			res.setEntry(i,2, crossRes.getEntry(2));
			res.setEntry(i,3, pi.getEntry(0));
			res.setEntry(i,4, pi.getEntry(1));
			res.setEntry(i,5, pi.getEntry(2));
			i++;
		}
		
		return res;
	}
	
	private RealMatrix constructB(HashMap<PdVector, HashMap<PdVector, Integer>> corresponding, PdVector[] normals){
		// Constructing b as shown here: http://resources.mpi-inf.mpg.de/deformableShapeMatching/EG2011_Tutorial/slides/2.1%20Rigid%20ICP.pdfS
		int n = corresponding.keySet().size();
		RealMatrix res = new Array2DRowRealMatrix(n, 1);
		int i = 0;
		for(PdVector pi : corresponding.keySet()) {
			HashMap<PdVector, Integer> results = corresponding.get(pi);
			PdVector qi = results.keySet().iterator().next();
			int indexOfNorm = results.get(qi);
			PdVector norm = normals[indexOfNorm];
			
			PdVector diff = PdVector.subNew(pi, qi);
			diff.multScalar(-1.0);
			
			double entry = diff.dot(norm);
			res.setEntry(i,0, entry);
			i++;
		}
		
		return res;
	}
	
	// ========================================================================================================================
	// =====================================                                ===================================================
	// =====================================        ASSIGNMENT 2            ===================================================
	// =====================================                                ===================================================
	// ========================================================================================================================
	
	public PnSparseMatrix computeGradientMatrix(){
		// Init PnSparseMatrix.
		int amtRows = m_surfP.getNumElements() * 3;
		int amtCols = m_surfP.getNumVertices();
		PnSparseMatrix res = new PnSparseMatrix(amtRows, amtCols, 3);
		
		int amtTriangles = m_surfP.getNumElements();
		for(int i = 0; i < amtTriangles; i++){
			PdMatrix g = computeTriangleMatrix(i, m_surfP.getElementNormals()[i]);
			int[] indices = m_surfP.getElement(i).getEntries();
			for(int j = 0; j < 3; j++){
				res.setEntry(3*i, indices[j], g.getEntry(0, j));
				res.setEntry(3*i + 1, indices[j], g.getEntry(1, j));
				res.setEntry(3*i + 2, indices[j], g.getEntry(2, j));
			}
		}
		
		return res;
	}
	
	private PdMatrix computeTriangleMatrix(int faceIndex, PdVector normal){
		int[] indicesOfVertices = m_surfP.getElement(faceIndex).getEntries();
		// Get p1, p2 and p3.
		PdVector p1 = m_surfP.getVertex(indicesOfVertices[0]);
		PdVector p2 = m_surfP.getVertex(indicesOfVertices[1]);
		PdVector p3 = m_surfP.getVertex(indicesOfVertices[2]);
		
		// Calculate e1, e2 and e3.
		PdVector e1 = PdVector.subNew(p3, p2);
		PdVector e2 = PdVector.subNew(p1, p3);
		PdVector e3 = PdVector.subNew(p2, p1);
		
		// Calculate area of the triangle.
		double area = calculateArea(faceIndex);
		
		// Construct the matrix.
		double multiplier = 1.0/(2.0 * area);
		PdMatrix res = new PdMatrix(3);
		
		res.setColumn(0, PdVector.crossNew(e1, normal));
		res.setColumn(1, PdVector.crossNew(e2, normal));
		res.setColumn(2, PdVector.crossNew(e3, normal));
		
		res.multScalar(multiplier);
		
		return res;
	}
	
	private double calculateArea(int faceIndex) {
		int[] indicesOfVertices = m_surfP.getElement(faceIndex).getEntries();
		// Get p1, p2 and p3.
		PdVector p1 = m_surfP.getVertex(indicesOfVertices[0]);
		PdVector p2 = m_surfP.getVertex(indicesOfVertices[1]);
		PdVector p3 = m_surfP.getVertex(indicesOfVertices[2]);
		
		// Calculate e1, e2 and e3.
		PdVector e1 = PdVector.subNew(p3, p2);
		PdVector e2 = PdVector.subNew(p1, p3);
		PdVector e3 = PdVector.subNew(p2, p1);
		
		// Calculate area of triangle, using Heron's formula.
		double a = e1.length();
		double b = e2.length();
		double c = e3.length();
		double s = (a+b+c)/2.0;
		double area = Math.sqrt(s*(s - a)*(s - b)*(s - c));
		return area;
	}
	
	public void transform(PdMatrix A) throws Exception {
		m_surfP.makeElementNormals();
		
		// Calculate G and G^T.
		PnSparseMatrix G = computeGradientMatrix();
		PnSparseMatrix GT = G.transposeNew();
		
		// Calculate Mv.
		int amtRows = m_surfP.getNumElements() * 3;
		PnSparseMatrix Mv = new PnSparseMatrix(amtRows, amtRows, 1);
		int faceIndex = 0;
		for(int i = 0; i < amtRows; i++){
			double area = calculateArea(faceIndex);
			Mv.setEntry(i, i, area);
			i++;
			Mv.setEntry(i, i, area);
			i++;
			Mv.setEntry(i, i, area);
			i++;
			faceIndex++;
		}
		
		// Calculate g-tilde.
		PdMatrix gTilde = new PdMatrix(amtRows, 3);
		for(int i = 0; i < amtRows; i = i + 3){
			int index = i/3;
			PdMatrix subG = computeTriangleMatrix(index, m_surfP.getElementNormal(index));
			PdMatrix subRes;
			if(m_surfP.getElement(index).hasTag(PsObject.IS_SELECTED)){
				subRes = new PdMatrix(A.m_data);
				subRes.mult(subRes, subG);
				subRes.transpose();
			}
			else {
				subRes = subG;
			}
			for(int j = 0; j < 3; j++){
				gTilde.setEntry(i, j, subRes.getEntry(i, j));
				gTilde.setEntry(i + 1, j, subRes.getEntry(i + 1, j));
				gTilde.setEntry(i + 2, j, subRes.getEntry(i + 2, j));
			}
		}
		
		// Compute G^TMvG.
		PnSparseMatrix interMatrix = PnSparseMatrix.multMatrices(GT, Mv, null);
		PnSparseMatrix matrix = PnSparseMatrix.multMatrices(interMatrix, G, null);
		
		// Compute G^TMvg-tilde. This has three results, called the x-, y-, and z-vectors.
		PdVector gTilde_x = gTilde.getColumn(0);
		PdVector gTilde_y = gTilde.getColumn(1);
		PdVector gTilde_z = gTilde.getColumn(2);
		
		PdVector b_x = PnSparseMatrix.rightMultVector(interMatrix, gTilde_x, null);
		PdVector b_y = PnSparseMatrix.rightMultVector(interMatrix, gTilde_y, null);
		PdVector b_z = PnSparseMatrix.rightMultVector(interMatrix, gTilde_z, null);
		
		// Prepare the x for the equation Ax = b.
		PdMatrix x = new PdMatrix(m_surfP.getNumVertices(), 3);
		
		PdVector x_x = new PdVector(m_surfP.getNumVertices());
		PdVector x_y = new PdVector(m_surfP.getNumVertices());
		PdVector x_z = new PdVector(m_surfP.getNumVertices());
		
		// Solve system.
		long factorization = PnMumpsSolver.factor(matrix, PnMumpsSolver.Type.GENERAL_SYMMETRIC);
		PnMumpsSolver.solve(factorization, x_x, b_x);
		PnMumpsSolver.solve(factorization, x_y, b_y);
		PnMumpsSolver.solve(factorization, x_z, b_z);
			
		x.setColumn(0, x_x);
		x.setColumn(1, x_y);
		x.setColumn(2, x_z);
			
		x.transpose();
		
		replaceVertices(x);
	}
	
	private void replaceVertices(PdMatrix x) {
		for(int i = 0; i < x.getNumCols(); i++) {
			PdVector newX = new PdVector(new double[]{x.getEntry(0, i), x.getEntry(1, i), x.getEntry(2, i)});
			m_surfP.setVertex(i, newX);
		}
	}
}