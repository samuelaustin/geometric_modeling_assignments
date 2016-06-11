package workshop;

import java.awt.Color;
import java.util.List;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Random;

import jv.geom.PgEdgeStar;
import jv.geom.PgElementSet;
import jv.project.PgGeometry;
import jv.vecmath.PdVector;
import jv.vecmath.PiVector;
import jvx.project.PjWorkshop;
import jvx.numeric.PnMatrix;
import jvx.numeric.PnSparseMatrix;
import jv.vecmath.PdMatrix;
import jvx.numeric.PnConjugateGradientMatrix;

public class Smoothing extends PjWorkshop {

	PgElementSet m_geom;
	PgElementSet m_geomSave;

	PnSparseMatrix G;
	PnSparseMatrix GT;
	PnSparseMatrix Mv;
	PnSparseMatrix M;
	PdMatrix X;
	
	HashMap<Integer, List<Integer>> neighbours = new HashMap<Integer, List<Integer>>();
	
	public Smoothing() {
		super("Smoothing");
		init();
	}
	
	@Override
	public void setGeometry(PgGeometry geom) {
		super.setGeometry(geom);
		m_geom 		= (PgElementSet)super.m_geom;
		m_geomSave 	= (PgElementSet)super.m_geomSave;
		addAllNeighbours();
	}
	
	public void init() {		
		super.init();
	}

	private void initMatrices(PgElementSet surface)
	{
		// Calculate G and G^T.
		surface.makeElementNormals();
		G = computeGradientMatrix(surface);
		GT = G.transposeNew();

		// Calculate Mv.
		int amtRows = surface.getNumElements() * 3;
		Mv = new PnSparseMatrix(amtRows, amtRows, 1);
		int faceIndex = 0;
		for(int i = 0; i < amtRows; i++){
			double area = calculateArea(surface, faceIndex);
			Mv.setEntry(i, i, area);
			i++;
			Mv.setEntry(i, i, area);
			i++;
			Mv.setEntry(i, i, area);
			faceIndex++;
		}
		
		M = new PnSparseMatrix(surface.getNumVertices(), surface.getNumVertices(), 1);
		// calculate M: diagonal matrix with entry (i,i) = 1/3 of the sum of the areas of all
		// triangles adjacent to vertex i.
		for(int i = 0; i < surface.getNumElements(); i++) {
			double area = calculateArea(surface, i)/3.0;
			int[] vertices = surface.getElement(i).m_data;
			M.setEntry(vertices[0], vertices[0], M.getEntry(vertices[0], vertices[0]) + area);
			M.setEntry(vertices[1], vertices[1], M.getEntry(vertices[1], vertices[1]) + area);
			M.setEntry(vertices[2], vertices[2], M.getEntry(vertices[2], vertices[2]) + area);
		}

		// Calculate g-tilde.
		X = new PdMatrix(surface.getNumVertices(), 3);
		for(int vIndex = 0; vIndex < surface.getNumVertices(); vIndex++){
			X.setRow(vIndex, surface.getVertex(vIndex));
		}
		
		System.out.println("Matrices created.");
	}

	private double calculateArea(PgElementSet surface, int faceIndex) {
		int[] indicesOfVertices = surface.getElement(faceIndex).getEntries();
		// Get p1, p2 and p3.
		PdVector p1 = surface.getVertex(indicesOfVertices[0]);
		PdVector p2 = surface.getVertex(indicesOfVertices[1]);
		PdVector p3 = surface.getVertex(indicesOfVertices[2]);
		
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

	private PdMatrix computeTriangleMatrix(PgElementSet surface, int faceIndex, PdVector normal){
		int[] indicesOfVertices = surface.getElement(faceIndex).getEntries();
		// Get p1, p2 and p3.
		PdVector p1 = surface.getVertex(indicesOfVertices[0]);
		PdVector p2 = surface.getVertex(indicesOfVertices[1]);
		PdVector p3 = surface.getVertex(indicesOfVertices[2]);
		
		// Calculate e1, e2 and e3.
		PdVector e1 = PdVector.subNew(p3, p2);
		PdVector e2 = PdVector.subNew(p1, p3);
		PdVector e3 = PdVector.subNew(p2, p1);
		
		// Calculate area of the triangle.
		double area = calculateArea(surface, faceIndex);
		
		// Construct the matrix.
		double multiplier = 1.0/(2.0 * area);
		PdMatrix res = new PdMatrix(3);
		
		res.setColumn(0, PdVector.crossNew(normal, e1));
		res.setColumn(1, PdVector.crossNew(normal, e2));
		res.setColumn(2, PdVector.crossNew(normal, e3));
		
		res.multScalar(multiplier);
		
		return res;
	}

	public PnSparseMatrix computeGradientMatrix(PgElementSet surface){
		if(surface == null)
		{
			System.out.println("Surface is null (computeGradientMatrix)");
		}
		int amtRows = surface.getNumElements() * 3;
		int amtCols = surface.getNumVertices();
		PnSparseMatrix res = new PnSparseMatrix(amtRows, amtCols, 3);
		
		int amtTriangles = surface.getNumElements();
		for(int i = 0; i < amtTriangles; i++){
			PdMatrix g = computeTriangleMatrix(surface, i, surface.getElementNormals()[i]);
			
			int[] indices = surface.getElement(i).getEntries();
			for(int j = 0; j < 3; j++){
				res.setEntry(3*i, indices[j], g.getEntry(0, j));
				res.setEntry(3*i + 1, indices[j], g.getEntry(1, j));
				res.setEntry(3*i + 2, indices[j], g.getEntry(2, j));
			}
		}
		
		return res;
	}

	public void iteratedAveraging(double timeStep) {
		int numVertices = m_geom.getNumVertices();
		
		PdVector[] newVertices = new PdVector[numVertices];
		for(int i = 0; i < numVertices; i++){
			newVertices[i] = PdVector.addNew(m_geom.getVertex(i), average(timeStep, i));
		}
		for(int j = 0; j < numVertices; j++){
			m_geom.setVertex(j, newVertices[j]);
		}
		m_geom.update(m_geom);
	}
	
	private PdVector average(double t, int index) {
		List<Integer> list = neighbours.get(index);
		double divisor = list.size();
		PdVector res = new PdVector(new double[]{0.0, 0.0, 0.0});
		
		// Calculate average of neighbours of vertex at index.
		for(int i = 0; i < list.size(); i++){
			res.add(m_geom.getVertex(list.get(i)));
		}
		res.multScalar(1.0/divisor);
		
		res.sub(m_geom.getVertex(index));
		res.multScalar(t);
		
		return res;
	}
	
	public void explicitMeanCurvatureFlow(double timeStep) {
		// Reconstruct matrices.
		initMatrices(m_geom);
				
		// Construct L = M^-1*G^T*Mv*G.
		PnSparseMatrix MInverse = new PnSparseMatrix(M.getNumRows(), M.getNumCols());
		for(int i = 0; i < M.getNumRows(); i++){
			MInverse.setEntry(i, i, 1.0/M.getEntry(i, i));
		}
				
		PnSparseMatrix inter1 = PnSparseMatrix.multMatrices(MInverse, GT, null);
		PnSparseMatrix inter2 = PnSparseMatrix.multMatrices(inter1, Mv, null);
		PnSparseMatrix L = PnSparseMatrix.multMatrices(inter2, G, null);
		
		PdVector X_x = new PdVector();
		PdVector X_y = new PdVector();
		PdVector X_z = new PdVector();
		
		for(int rowIndex = 0; rowIndex < X.getNumRows(); rowIndex++) {
			X_x.addEntry(X.getEntry(rowIndex,0));
			X_y.addEntry(X.getEntry(rowIndex,1));
			X_z.addEntry(X.getEntry(rowIndex,2));
		}
		
		L.multScalar(timeStep);
		
		// Compute -timeStep*L*x.
		PdVector L_x = PnSparseMatrix.rightMultVector(L, X_x, null);
		PdVector L_y = PnSparseMatrix.rightMultVector(L, X_y, null);
		PdVector L_z = PnSparseMatrix.rightMultVector(L, X_z, null);
		
		// Replace all vectors x with x - timeStep*L*x.
		for(int i = 0; i < m_geom.getNumVertices(); i++) {
			PdVector res = new PdVector(new double[]{X_x.getEntry(i) - L_x.getEntry(i),
					X_y.getEntry(i) - L_y.getEntry(i),
					X_z.getEntry(i) - L_z.getEntry(i)});
			m_geom.setVertex(i, res);
		}
		m_geom.update(m_geom);
		System.out.println("Done");
	}
	
	public void implicitMeanCurvatureFlow(double timeStep) {
		// Reconstruct matrices.
		initMatrices(m_geom);
		
		// Construct S = G^T*Mv*G.
		PnSparseMatrix interMatrix = PnSparseMatrix.multMatrices(GT, Mv, null);
		PnSparseMatrix S = PnSparseMatrix.multMatrices(interMatrix, G, null);
		
		// Compute M + tS.
		S.multScalar(timeStep);
		S.add(M);
		
		// Set up X-vectors.
		PdVector X_x = new PdVector();
		PdVector X_y = new PdVector();
		PdVector X_z = new PdVector();
		
		for(int rowIndex = 0; rowIndex < m_geom.getNumVertices(); rowIndex++) {
			X_x.addEntry(X.getEntry(rowIndex,0));
			X_y.addEntry(X.getEntry(rowIndex,1));
			X_z.addEntry(X.getEntry(rowIndex,2));
		}
		
		// Set up b for Ax = b.
		PdVector b_x = PnSparseMatrix.rightMultVector(M, X_x, null);
		PdVector b_y = PnSparseMatrix.rightMultVector(M, X_y, null);
		PdVector b_z = PnSparseMatrix.rightMultVector(M, X_z, null);
				
		PdVector x_x = new PdVector(m_geom.getNumVertices());
		PdVector x_y = new PdVector(m_geom.getNumVertices());
		PdVector x_z = new PdVector(m_geom.getNumVertices());
				
		// Solve system.
		PnConjugateGradientMatrix conjGradMatrix = new PnConjugateGradientMatrix();
		conjGradMatrix.solve(S, x_x, b_x);
		conjGradMatrix.solve(S, x_y, b_y);
		conjGradMatrix.solve(S, x_z, b_z);
		
		// Prepare the x for the equation Ax = b.
		PdMatrix x = new PdMatrix(m_geom.getNumVertices(), 3);
		
		x.setColumn(0, x_x);
		x.setColumn(1, x_y);
		x.setColumn(2, x_z);
		
		replaceVertices(m_geom, x);
		System.out.println("Done");
	}

	private void replaceVertices(PgElementSet surface, PdMatrix x) {
		for(int i = 0; i < x.getNumRows(); i++) {
			PdVector old = surface.getVertex(i);
			PdVector newX = new PdVector(new double[]{x.getEntry(i, 0), x.getEntry(i, 1), x.getEntry(i, 2)});

			surface.setVertex(i, newX);
		}
		surface.update(surface);
	}
	
	private void addAllNeighbours() {
		int amtFaces = m_geom.getNumElements();
		for(int i = 0; i < amtFaces; i++) {
			int[] data = m_geom.getElement(i).m_data;
			addNeighbour(data[0], data[1]);
			addNeighbour(data[1], data[0]);
			addNeighbour(data[1], data[2]);
			addNeighbour(data[2], data[1]);
			addNeighbour(data[2], data[0]);
			addNeighbour(data[0], data[2]);
		}
	}
	
	private void addNeighbour(int index, int neighbour) {
		List<Integer> list = neighbours.get(index);
		if(list == null) {
			list = new ArrayList<Integer>();
		}
		if(!list.contains(neighbour)) {
			list.add(neighbour);
			neighbours.put(index, list);
		}
	}
}
