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
import jvx.project.PjWorkshop;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.LUDecomposition;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.SingularValueDecomposition;
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
		int n = Math.min(m_surfP.getNumVertices(), 200);
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
		
		// Calculate r and t.
		double[] r_entries = {x.getEntry(0, 0), x.getEntry(1, 0), x.getEntry(2, 0)};
		PdVector r = new PdVector(r_entries);
		
		double[] t_entries = {x.getEntry(3, 0), x.getEntry(4, 0), x.getEntry(5, 0)};
		PdVector t = new PdVector(t_entries);
		
		// Update the mesh.
		rotateAndTranslatePointToPlane(m_surfP, r, t);
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
	
	private void rotateAndTranslatePointToPlane(PgElementSet set, PdVector r, PdVector t)
	{
		System.out.println("Rotating and translating");
		
		for(int i = 0; i < set.getNumVertices(); i++){
			PdVector old = set.getVertex(i);
			PdVector trans = PdVector.addNew(old, PdVector.crossNew(old, r));
			trans = PdVector.addNew(trans, t);
			set.setVertex(i, trans);
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
}
