package workshop;

import java.awt.Color;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashMap;
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
import jv.vecmath.PdVector;
import jv.vecmath.PiVector;
import jv.vecmath.PuMath;
import jv.viewer.PvDisplay;
import jv.project.PvGeometryIf;

import jvx.project.PjWorkshop;

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
	
	public void iterativeClosestPoint(){
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
		HashMap<PdVector, Double> results = new HashMap<PdVector, Double>();
		HashMap<PdVector, PdVector> corresponding = new HashMap<PdVector, PdVector>();
		double meanDist = 0.0;
		for(int j = 0; j < n; j++){
			double dist = Double.MAX_VALUE;
			PdVector[] vectorsOfQ = m_surfQ.getVertices();
			for(int j2 = 0; j2 < vectorsOfQ.length; j2++){
				if(PdVector.dist(vectors[j], vectorsOfQ[j2]) < dist){
					dist = PdVector.dist(vectors[j], vectorsOfQ[j2]);
				}
			}
			meanDist += dist;
			results.put(vectors[j], dist);
		}
		meanDist = meanDist/n;
		
		// Remove all vectors with distance greater than k*meanDist.
		double k = 2.0;
		for(PdVector v : results.keySet()) {
			if(results.get(v) > k*meanDist){
				results.remove(v);
			}
		}
		
		// Compute centroid of P and Q.
		PdVector pCentroid = computeCentroid(results.keySet());
		PdVector qCentroid = computeCentroid(Arrays.asList(m_surfQ.getVertices()));
		
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
	
}
