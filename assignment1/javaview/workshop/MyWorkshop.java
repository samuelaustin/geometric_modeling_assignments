package workshop;

import java.awt.Color;
import java.util.Random;

import jv.geom.PgEdgeStar;
import jv.geom.PgElementSet;
import jv.project.PgGeometry;
import jv.vecmath.PdVector;
import jv.vecmath.PiVector;
import jvx.project.PjWorkshop;

public class MyWorkshop extends PjWorkshop {

	PgElementSet m_geom;
	PgElementSet m_geomSave;
	double[] shapeRegularities;
	
	public MyWorkshop() {
		super("My Workshop");
		init();
	}
	
	@Override
	public void setGeometry(PgGeometry geom) {
		super.setGeometry(geom);
		m_geom 		= (PgElementSet)super.m_geom;
		m_geomSave 	= (PgElementSet)super.m_geomSave;
	}
	
	public void init() {		
		super.init();
	}
	
	public void makeRandomElementColors() {
		//assure that the color array is allocated
		m_geom.assureElementColors();
		
		Random rand = new Random();
		Color randomColor;
		
		int noe = m_geom.getNumElements();
		for(int i=0; i<noe; i++){
			randomColor = Color.getHSBColor(rand.nextFloat(), 1.0f, 1.0f);//new Color(rand.nextFloat(), rand.nextFloat(), rand.nextFloat());
			m_geom.setElementColor(i, randomColor);
		}
		m_geom.showElementColorFromVertices(false);
		m_geom.showElementColors(true);	
		m_geom.showSmoothElementColors(false);
	}
	
	public void makeRandomVertexColors() {
		//assure that the color array is allocated
		m_geom.assureVertexColors();
		
		Random rand = new Random();
		Color randomColor;
		
		int nov = m_geom.getNumVertices();
		for(int i=0; i<nov; i++){
			randomColor = Color.getHSBColor(rand.nextFloat(), 1.0f, 1.0f);
			m_geom.setVertexColor(i, randomColor);
		}
		
		m_geom.showElementColors(true);	
		m_geom.showVertexColors(true);
		m_geom.showElementColorFromVertices(true);	
		m_geom.showSmoothElementColors(true);
	}
	
	
	public void setXOff(double xOff) {
		int nov = m_geom.getNumVertices();
		PdVector v = new PdVector(3);
		// the double array is v.m_data 
		for (int i=0; i<nov; i++) {
			v.copyArray(m_geomSave.getVertex(i));
			v.setEntry(0, v.getEntry(0)+xOff);
			m_geom.setVertex(i, v);
		}
	}
	
	public int calcGenus(){
		int nov = m_geom.getNumVertices();
		int noe = m_geom.getNumEdges();
		int nof = m_geom.getNumElements();
		int var = nov - noe + nof;
		return 1 - (var/2);
		
	}
	
	public double calcVolume(){
		int numOfElems = m_geom.getNumElements();
		double sum = 0.0;
		for(int i = 0; i < numOfElems; i++){
			PiVector indices = m_geom.getElement(i);
			PdVector x = m_geom.getVertex(indices.m_data[0]);
			PdVector y = m_geom.getVertex(indices.m_data[1]);
			PdVector z = m_geom.getVertex(indices.m_data[2]);
			sum = sum + (dotProduct(x, crossProduct(y, z))/6.0);
		}
		return Math.abs(sum);
	}
	
	public double[] calcShapeReg(){
		double minShapeReg = Double.MAX_VALUE;
		double maxShapeReg = Double.MIN_VALUE;
		double sum;
		int numOfElems = m_geom.getNumElements();
		double[] regs = new double[numOfElems];
		for(int i = 0; i < numOfElems; i++){
			PiVector triangle = m_geom.getElement(i);
			PdVector vect1 = m_geom.getVertex(triangle.m_data[0]);
			PdVector vect2 = m_geom.getVertex(triangle.m_data[1]);
			PdVector vect3 = m_geom.getVertex(triangle.m_data[2]);
			double area = PdVector.area(vect1, vect2, vect3);
			double a = distance(vect1, vect2);
			double b = distance(vect2, vect3);
			double c = distance(vect1, vect3);
			double divisor = 0.5 * (a + b + c);
			double radiusInner = area/divisor;
			double abc = a*b*c;
			double divisor2 = Math.sqrt((a+b+c)*(b+c-a)*(c+a-b)*(a+b-c));
			double radiusOuter = abc/divisor2;
			double reg = radiusInner/radiusOuter;
			sum += reg;
			if(reg < minShapeReg){
				minShapeReg = reg;
			}
			if(reg > maxShapeReg){
				maxShapeReg = reg;
			}
			regs[i] = reg;
		}
		double mean = sum/numOfElems;
		double stddev;
		for(int j = 0; j < numOfElems; j++){
			stddev += Math.pow(regs[j] - mean, 2);
		}
		stddev = Math.sqrt(stddev/numOfElems);
		shapeRegularities = regs;
		double[] results = new double[4];
		results[0] = minShapeReg;
		results[1] = maxShapeReg;
		results[2] = mean;
		results[3] = stddev;
		return results;
	}
	
	public double[] calcValence(){
		// Valence of vertex: amount of adjacent edges.
		int amtEdges = m_geom.getNumEdges();
		int[] valences = new int[m_geom.getNumVertices()];
		double sum = 0;
		for(int i = 0; i < amtEdges; i++){
			PgEdgeStar vector = m_geom.getEdgeStar(i);
			valences[vector.getVertexInd(0)]++;
			valences[vector.getVertexInd(1)]++;
			sum = sum + 2.0;
		}
		double avgValence = sum/m_geom.getNumVertices();
		double stddev = 0.0;
		int min_valence = Integer.MAX_VALUE;
		int max_valence = Integer.MIN_VALUE;
		for(int j = 0; j < valences.length; j++){
			stddev += Math.pow(valences[j] - avgValence, 2);
			if(valences[j] < min_valence){
				min_valence = valences[j];
			}
			if(valences[j] > max_valence){
				max_valence = valences[j];
			}
		}
		stddev = Math.sqrt(stddev/m_geom.getNumVertices());
		double[] results = new double[4];
		results[0] = avgValence;
		results[1] = min_valence;
		results[2] = max_valence;
		results[3] = stddev;
		return results;
	}
	
	private double dotProduct(PdVector x, PdVector y){
		return ((x.getEntry(0) * y.getEntry(0)) + (x.getEntry(1) * y.getEntry(1)) + (x.getEntry(2) * y.getEntry(2)));
	}
	
	private PdVector crossProduct(PdVector x, PdVector y){
		PdVector res = new PdVector(3);
		double cross_x = ((x.getEntry(1) * y.getEntry(2)) - (x.getEntry(2) * y.getEntry(1)));
		double cross_y = ((x.getEntry(2) * y.getEntry(0)) - (x.getEntry(0) * y.getEntry(2)));
		double cross_z = ((x.getEntry(0) * y.getEntry(1)) - (x.getEntry(1) * y.getEntry(0)));
		res.setEntry(0, cross_x);
		res.setEntry(1, cross_y);
		res.setEntry(2, cross_z);
		return res;
	}
	
	private double distance(PdVector x, PdVector y){
		return Math.sqrt((
				Math.pow((x.getEntry(0) - y.getEntry(0),2)) + 
				Math.pow((x.getEntry(1) - y.getEntry(1),2)) +
				Math.pow((x.getEntry(2) - y.getEntry(2),2))
				));
	}
}
