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
			HashMap<PdVector, PdVector> corresponding = GetClosestPointToPlaneDistances(vectors, m_surfQ);
		}
		else {
			HashMap<PdVector, PdVector> corresponding = GetClosestPointToPointDistances(vectors, m_surfQ);
		}
		
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
		HashSet<PdVector> entrySet = new HashSet<PdVector>();
		for(PdVector v : results.keySet()) {
			if(results.get(v) > k*meanDist){
				corresponding.remove(v);
			}
		}

		return corresponding;
	}

	private HashMap<PdVector, PdVector> GetClosestPointToPlaneDistances(PdVector[] vectors, PgElementSet surface)
	{
		// Calculate min. distance between this set and points of Q.
		HashMap<PdVector, Double> results = new HashMap<PdVector, Double>();
		HashMap<PdVector, PdVector> corresponding = new HashMap<PdVector, PdVector>();
		double meanDist = 0.0;

		PdVector[] vectorsOfQ = surface.getVertices();
		PdVector[] normalsOfQ = surface.getElementNormals();

		for(int j = 0; j < vectors.length; j++){
			double minDist = Double.MAX_VALUE;
			PdVector closestPoint = null;
			for(int j2 = 0; j2 < vectorsOfQ.length; j2++){
				PdVector point = PdVector.subNew(vectors[j], vectorsOfQ[j2]);
				PdVector projectionOntoNormal = GetProjectionPOntoQ(point, normalsOfQ[j2]);
				double dist = projectionOntoNormal.dot(projectionOntoNormal);
				if(dist < minDist){
					minDist = PdVector.dist(vectors[j], vectorsOfQ[j2]);
					closestPoint = vectorsOfQ[j2];
				}
			}
			meanDist += minDist;
			results.put(vectors[j], minDist);
			corresponding.put(vectors[j], closestPoint);
		}
		meanDist = meanDist/vectors.length;

		double k = 2.0;
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
			results[i] = (m.getEntry(i, 0)*v.getEntry(0)) + (m.getEntry(i, 1)*v.getEntry(1)) + (m.getEntry(i, 2)*v.getEntry(2));
		}
		return new PdVector(results);
	}
	
	private PdVector GetProjectionPOntoQ(PdVector p, PdVector q)
	{
		PdVector result = (PdVector) q.clone();
		result.multScalar(p.dot(q) / q.dot(q));
		return result;
	}
	
	private void rotateAndTranslate(PgElementSet set, PdMatrix r, PdVector t)
	{
		for(int i = 0; i < set.getNumVertices(); i++){
			set.setVertex(i, PdVector.addNew(matrixMult(r, set.getElement(i)), t));
		}
	}
}
