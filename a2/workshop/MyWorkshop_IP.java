package workshop;

import java.awt.Button;
import java.awt.FlowLayout;
import java.awt.Panel;
import java.awt.TextField;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;
import java.awt.geom.Area;
import java.util.*;

import org.omg.Messaging.SYNC_WITH_TRANSPORT;

import jv.geom.PgBndPolygon;
import jv.geom.PgEdgeStar;
import jv.number.PuDouble;
import jv.object.PsConfig;
import jv.object.PsDialog;
import jv.object.PsUpdateIf;
import jv.vecmath.PiVector;
import jvx.project.PjWorkshop_IP;

import jv.vecmath.PdVector;
import jv.geom.PgElementSet;
import jv.vecmath.PdMatrix;

public class MyWorkshop_IP extends PjWorkshop_IP implements ActionListener {

	protected Button m_bMakeRandomElementColors;
	protected Button m_bMakeRandomVertexColors;
	protected PuDouble m_xOff;
	protected TextField	m_textField;

	// task 1.1 Compute G sparse matrix
	protected Button m_bCalculateGMatrix;
	protected TextField m_GMatrixResult;
	protected PdMatrix m_GMatrix;

	// task 1.2 Compute combinatorial Laplacian matrix L
	protected Button m_bCalculateLMatrix;
	protected TextField m_LMatrixResult;
	protected PdMatrix m_LMatrix; // Combinatorial Laplacian matrix

	// task 1.3 Compute Mass Matrix M
	protected Button m_bCalculateMMatrix;
	protected TextField m_MMatrixResult;
	protected PdMatrix m_MMatrix; // Mass matrix

	// task 1.4 Compute Cotangent matrix S
	protected Button m_bCalculateSMatrix;
	protected TextField m_SMatrixResult;
	protected PdMatrix m_SMatrix; // Cotangent matrix
	
	MyWorkshop m_ws;
	
	public MyWorkshop_IP() {
		super();
		if(getClass() == MyWorkshop_IP.class)
			init();
	}
	
	public void init() {
		super.init();
		setTitle("Mesh and Surface Analysis");
	}
	
	public String getNotice() {
		return "This workshop can be used to calculate perform mesh and surface analysis to retrieve properties such as genus, volume and amount of components.";
	}
	
	public void setParent(PsUpdateIf parent) {
		super.setParent(parent);
		m_ws = (MyWorkshop)parent;
	
		addSubTitle("Properties");

		m_textField = addTextField("Number of Triangles:", 20);
		m_textField.setEditable(false);
		m_textField.setText(String.valueOf(m_ws.m_geom.getNumElements()));

		m_xOff = new PuDouble("X Offset");
		m_xOff.setDefBounds(-10,10,0.1,1);
		m_xOff.addUpdateListener(this);
		m_xOff.init();
		add(m_xOff.getInfoPanel());

		// task 1.1 Compute G sparse matrix
        m_bCalculateGMatrix = new Button("Calculate G Matrix");
        m_bCalculateGMatrix.addActionListener(this);
        m_GMatrixResult = new TextField("G Matrix: ", 20);
        m_GMatrixResult.setEditable(false);
        Panel panel2 = new Panel(new FlowLayout(FlowLayout.CENTER));
        panel2.add(m_bCalculateGMatrix);
        panel2.add(m_GMatrixResult);
        add(panel2);

		// task 1.2 Compute combinatorial Laplacian matrix L
		m_bCalculateLMatrix = new Button("Calculate L Matrix");
		m_bCalculateLMatrix.addActionListener(this);
		m_LMatrixResult = new TextField("L Matrix: ", 20);
		m_LMatrixResult.setEditable(false);
		Panel panel3 = new Panel(new FlowLayout(FlowLayout.CENTER));
		panel3.add(m_bCalculateLMatrix);
		panel3.add(m_LMatrixResult);
		add(panel3);

		// task 1.3 Compute Mass Matrix M
		m_bCalculateMMatrix = new Button("Calculate M Matrix");
		m_bCalculateMMatrix.addActionListener(this);
		m_MMatrixResult = new TextField("M Matrix: ", 20);
		m_MMatrixResult.setEditable(false);
		Panel panel4 = new Panel(new FlowLayout(FlowLayout.CENTER));
		panel4.add(m_bCalculateMMatrix);
		panel4.add(m_MMatrixResult);
		add(panel4);

		// task 1.4 Compute Cotangent matrix S
		m_bCalculateSMatrix = new Button("Calculate S Matrix");
		m_bCalculateSMatrix.addActionListener(this);
		m_SMatrixResult = new TextField("S Matrix: ", 20);
		m_SMatrixResult.setEditable(false);
		Panel panel5 = new Panel(new FlowLayout(FlowLayout.CENTER));
		panel5.add(m_bCalculateSMatrix);
		panel5.add(m_SMatrixResult);
		add(panel5);

		
		validate();
	}
	
	
	public boolean update(Object event) {
		if (event == m_xOff) {
			m_ws.setXOff(m_xOff.getValue());
			m_ws.m_geom.update(m_ws.m_geom);
			return true;
		} else
			return super.update(event);
	}
	
	/**
	 * Handle action events fired by buttons etc.
	 */
	public void actionPerformed(ActionEvent event) {
		Object source = event.getSource();
		if (source == m_bMakeRandomElementColors) {
			m_ws.makeRandomElementColors();
			m_ws.m_geom.update(m_ws.m_geom);
			return;
		}
		else if (source == m_bMakeRandomVertexColors) {
			m_ws.makeRandomVertexColors();
			m_ws.m_geom.update(m_ws.m_geom);
			return;
		}

		// task 1.1 calculate G	
		else if (source == m_bCalculateGMatrix) {
			m_GMatrix = calculateGMatrix(m_ws.m_geom);
            testGMatrix();
            return;
        }

		// task 1.2 calculate L
		else if (source == m_bCalculateLMatrix) {
			m_LMatrix = calculateLMatrix(m_ws.m_geom);
			testLMatrix();
			return;
		}

		// task 1.3 calculate M
		else if (source == m_bCalculateMMatrix) {
			m_MMatrix = calculateMassMatrix(m_ws.m_geom);
			testMassMatrix();
			return;
		}

		// task 1.4 calculate S
		else if (source == m_bCalculateSMatrix) {
			// m_SMatrix = calculateSMatrix(m_ws.m_geom);
			testSMatrix();
			return;
		}

		else
			super.actionPerformed(event);
	}

	// task 1.1 G matrix
	private PdMatrix calculateGMatrix(PgElementSet mesh) {
		m_GMatrixResult.setText("Calculating...");

		int nVertices = mesh.getNumVertices();
		int nTriangles = mesh.getNumElements();

		PdMatrix G = new PdMatrix(3 * nTriangles, nVertices);

		// Loop over all triangles in the mesh
		for (int i = 0; i < mesh.getNumElements(); i++) {
			PiVector triangle = mesh.getElement(i);

			// Print the triangle
			// System.out.println("Triangle " + i + ": " + triangle.getEntry(0) + ", " + triangle.getEntry(1) + ", "
			// 		+ triangle.getEntry(2));

			// Get the vertices of the triangle
			PdVector[] vertices = new PdVector[3];
			for (int j = 0; j < 3; j++) {
				vertices[j] = mesh.getVertex(triangle.getEntry(j));
			}

			// List of indices of the vertices of the triangle
			List<Integer> vertexIndices = new ArrayList<Integer>();
			for (int j = 0; j < 3; j++) {
				vertexIndices.add(triangle.getEntry(j));
			}
			
			// Calculate the gradient of the function on this triangle
			PdMatrix gradient = calculateTriangleGradient(vertices, vertexIndices);

			// Add the gradient to the G matrix (vertial stacking), gradient is 3x3 matrix
			for (int j = 0; j < 3; j++) {
				int ind = vertexIndices.get(j);
				for (int k = 0; k < 3; k++) {
					// set entry in gradientMatrix_true at correct row and column position
					G.setEntry(i * 3 + k, ind, gradient.getEntry(k, j));
				}
			}
		}
		// print G matrix
		// System.out.println("G matrix: " + G.toString());
		m_GMatrixResult.setText("Calculated!");

		return G;
	}

	// task 1.1 G matrix (Gradient of a function on a triangle)
	private PdMatrix calculateTriangleGradient(PdVector[] vertices, List<Integer> vertexIndices) {
		// Make vertex and index pairs
		Map<PdVector, Integer> vertexIndexMap = new HashMap<PdVector, Integer>();
		for (int i = 0; i < 3; i++) {
			vertexIndexMap.put(vertices[i], vertexIndices.get(i));
		}

		// Order the vertices by index, in ascending order
		List<Map.Entry<PdVector, Integer>> entryList = new ArrayList<>(vertexIndexMap.entrySet());
		entryList.sort(Map.Entry.comparingByValue());

		// Get the vertices in the correct order
		PdVector[] orderedVertices = new PdVector[3];
		for (int i = 0; i < 3; i++) {
			orderedVertices[i] = entryList.get(i).getKey();
		}

		// Calculate the area of the triangle
		double area = calculateTriangleArea(orderedVertices);

		// Calculate the gradient as a 3x3 matrix
		PdMatrix gradientMatrix = new PdMatrix(3, 3);

		// Calculate the triangle normal, the normal is consistent with the order of the vertices
		PdVector triangleNormal = PdVector.crossNew(PdVector.subNew(orderedVertices[1], orderedVertices[0]), PdVector.subNew(orderedVertices[2], orderedVertices[0]));
		triangleNormal.normalize();

		// Calculate the gradient of the function on the triangle
		for (int i = 0; i < 3; i++) {
			// Edge vector, following the order of the vertices
			PdVector edgeVector = PdVector.subNew(orderedVertices[(i + 2) % 3], orderedVertices[(i + 1) % 3]);

			// Calculate the gradient of the function on the triangle
			PdVector gradient = PdVector.crossNew(triangleNormal, edgeVector);

			// Add the gradient to the gradient matrix
			gradientMatrix.setColumn(i, gradient);
		}
		gradientMatrix.multScalar(1 / (2 * area));
		return gradientMatrix;
	}

	// task 1.1 G matrix (Area of a triangle)
	private Double calculateTriangleArea(PdVector[] vertices) {
		PdVector crossProduct = PdVector.crossNew(PdVector.subNew(vertices[1], vertices[0]), PdVector.subNew(vertices[2], vertices[0]));
		double area = crossProduct.length() / 2;

		return area;
	}

	// task 1.1 G matrix (Test G using a generated triangle mesh)
	private void testGMatrix() {
		// Create a 2 triangle mesh
		PgElementSet mesh = new PgElementSet(3);

		// Add vertices
		mesh.addVertex(new PdVector(0, 0, 0)); // 0
		mesh.addVertex(new PdVector(1, 0, 0)); // 1
		mesh.addVertex(new PdVector(0, 1, 0)); // 2
		mesh.addVertex(new PdVector(0, 0, 1)); // 3

		// Add triangles
		mesh.addElement(new PiVector(0, 1, 2)); // 0
		mesh.addElement(new PiVector(0, 1, 3)); // 1

		List<List<Integer>> vertexIndices =  new ArrayList<List<Integer>>();
		vertexIndices.add(Arrays.asList(0, 1, 2));
		vertexIndices.add(Arrays.asList(0, 1, 3));

		// Manually calculate the gradient matrix
		PdMatrix gradientMatrix_true = new PdMatrix(6, 4);

		// Triangle 0
		double area0 = calculateTriangleArea(new PdVector[] { mesh.getVertex(0), mesh.getVertex(1), mesh.getVertex(2) });
		PdVector e1 = PdVector.subNew(mesh.getVertex(2), mesh.getVertex(1));
		PdVector e2 = PdVector.subNew(mesh.getVertex(0), mesh.getVertex(2));
		PdVector e3 = PdVector.subNew(mesh.getVertex(1), mesh.getVertex(0));
		PdVector n0 = PdVector.crossNew(e1, e2);
		n0.normalize();


		PdVector g0_0 = PdVector.crossNew(n0, e1);
		PdVector g0_1 = PdVector.crossNew(n0, e2);
		PdVector g0_2 = PdVector.crossNew(n0, e3);
		g0_0.multScalar(1 / (2 * area0));
		g0_1.multScalar(1 / (2 * area0));
		g0_2.multScalar(1 / (2 * area0));

		PdMatrix g0 = new PdMatrix(3, 3);
		g0.setColumn(0, g0_0);
		g0.setColumn(1, g0_1);
		g0.setColumn(2, g0_2);
		// System.out.println("g0: " + g0.toString());

		// Triangle 1
		double area1 = calculateTriangleArea(new PdVector[] { mesh.getVertex(0), mesh.getVertex(1), mesh.getVertex(3) });
		PdVector e4 = PdVector.subNew(mesh.getVertex(3), mesh.getVertex(1));
		PdVector e5 = PdVector.subNew(mesh.getVertex(0), mesh.getVertex(3));
		PdVector e6 = PdVector.subNew(mesh.getVertex(1), mesh.getVertex(0));

		PdVector n1 = PdVector.crossNew(e4, e5);
		n1.normalize();

		PdVector g1_0 = PdVector.crossNew(n1, e4);
		PdVector g1_1 = PdVector.crossNew(n1, e5);
		PdVector g1_2 = PdVector.crossNew(n1, e6);
		g1_0.multScalar(1 / (2 * area1));
		g1_1.multScalar(1 / (2 * area1));
		g1_2.multScalar(1 / (2 * area1));

		PdMatrix g1 = new PdMatrix(3, 3);
		g1.setColumn(0, g1_0);
		g1.setColumn(1, g1_1);
		g1.setColumn(2, g1_2);

		List<PdMatrix> gMatrices = new ArrayList<PdMatrix>();
		gMatrices.add(g0);
		gMatrices.add(g1);

		// Add the gradient to the G matrix (vertial stacking), at the correct column (vertex index)
		for (int i = 0; i < 2; i++) {
			List<Integer> vertexIndexList = vertexIndices.get(i);
			PdMatrix currentGMatrix = gMatrices.get(i); // current gradient matrix
			for (int j = 0; j < 3; j++) {
				int ind = vertexIndexList.get(j);
				for (int k = 0; k < 3; k++) {
					// set entry in gradientMatrix_true at correct row and column position
					gradientMatrix_true.setEntry(i * 3 + k, ind, currentGMatrix.getEntry(k, j));
				}
			}
		}

		// print gradient matrix
		// System.out.println("Gradient_true matrix: " + gradientMatrix_true.toString());

		// Calculate the gradient matrix using the function
		PdMatrix gradientMatrix_calculated = calculateGMatrix(mesh);
		// System.out.println("Gradient_calculated matrix: " + gradientMatrix_calculated.toString());

		// Compare the two gradient matrices
		for (int i = 0; i < 6; i++) {
			for (int j = 0; j < 4; j++) {
				if (Math.abs(gradientMatrix_true.getEntry(i, j) - gradientMatrix_calculated.getEntry(i, j)) > 0.0001) {
					m_GMatrixResult.setText("Failed");
					return;
				}
			}
		}
		m_GMatrixResult.setText("Success");
	}

	// task 1.2 L matrix (Calculate the Combinatorial Laplacian matrix)
	private PdMatrix calculateLMatrix(PgElementSet mesh) {
		m_LMatrixResult.setText("Calculating...");
		int nVertices = mesh.getNumVertices();
		PdMatrix L = new PdMatrix(nVertices, nVertices);
		List<Set<Integer>> neighbors = new ArrayList<>(nVertices);

		// Initialize neighbors sets
		for (int i = 0; i < nVertices; i++) {
			neighbors.add(new HashSet<>());
		}

		// Populate neighbors sets
		for (int i = 0; i < mesh.getNumElements(); i++) {
			PiVector triangle = mesh.getElement(i);
			for (int j = 0; j < 3; j++) {
				neighbors.get(triangle.getEntry(j)).add(triangle.getEntry((j+1)%3));
				neighbors.get(triangle.getEntry(j)).add(triangle.getEntry((j+2)%3));
			}
		}

		// Populate Laplacian matrix
		for (int i = 0; i < nVertices; i++) {
			L.setEntry(i, i, neighbors.get(i).size());
			for (Integer j : neighbors.get(i)) {
				L.setEntry(i, j, -1);
			}
		}
		m_LMatrixResult.setText("Calculated!");
		return L;
	}


	// task 1.2 L matrix (Test L using a generated triangle mesh)
	private void testLMatrix() {
		// Create a triangle mesh
		PgElementSet mesh = new PgElementSet(3);

		mesh.addVertex(new PdVector(0, 0, 0)); // p1
		mesh.addVertex(new PdVector(1, 0, 0)); // p2
		mesh.addVertex(new PdVector(0, 1, 0)); // p3

		mesh.addElement(new PiVector(0, 1, 2)); // a triangle

		// Manually calculate the L matrix
		PdMatrix L_true = new PdMatrix(3, 3);
		L_true.setEntry(0, 0, 2);
		L_true.setEntry(0, 1, -1);
		L_true.setEntry(0, 2, -1);
		L_true.setEntry(1, 0, -1);
		L_true.setEntry(1, 1, 2);
		L_true.setEntry(1, 2, -1);
		L_true.setEntry(2, 0, -1);
		L_true.setEntry(2, 1, -1);
		L_true.setEntry(2, 2, 2);
	

		// System.out.println("L_true matrix: " + L_true.toString());

		// Calculate the gradient matrix using the function
		PdMatrix L_calculated = calculateLMatrix(mesh);
		// System.out.println("L_calculated matrix: " + L_calculated.toString());

		// Compare the two matrices
		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < 3; j++) {
				if (Math.abs(L_true.getEntry(i, j) - L_calculated.getEntry(i, j)) > 0.0001) {
					m_LMatrixResult.setText("Failed");
					return;
				}
			}
		}
		m_LMatrixResult.setText("Success");
	}

	// task 1.3 Compute Mass Matrix M
	private PdMatrix calculateMassMatrix(PgElementSet mesh) {
		m_MMatrixResult.setText("Calculating...");

		int nVertices = mesh.getNumVertices();
		PdMatrix M = new PdMatrix(nVertices, nVertices);
		List<Double> areas = new ArrayList<>(nVertices);

		// Initialize areas list
		for (int i = 0; i < nVertices; i++) {
			areas.add(0.0);
		}

		// Populate areas list
		for (int i = 0; i < mesh.getNumElements(); i++) {
			PiVector triangle = mesh.getElement(i);
			PdVector p1 = mesh.getVertex(triangle.getEntry(0));
			PdVector p2 = mesh.getVertex(triangle.getEntry(1));
			PdVector p3 = mesh.getVertex(triangle.getEntry(2));
			PdVector[] vertices = {p1, p2, p3}; 
			double area = calculateTriangleArea(vertices);
			areas.set(triangle.getEntry(0), areas.get(triangle.getEntry(0)) + area);
			areas.set(triangle.getEntry(1), areas.get(triangle.getEntry(1)) + area);
			areas.set(triangle.getEntry(2), areas.get(triangle.getEntry(2)) + area);
		}

		// Populate Mass matrix
		for (int i = 0; i < nVertices; i++) {
			M.setEntry(i, i, areas.get(i) / 3);
		}

		m_MMatrixResult.setText("Calculated!");
		return M;
	}

	// task 1.3 Compute Mass Matrix M (Test M using a generated triangle mesh)
	private void testMassMatrix() {
		// Create a triangle mesh
		PgElementSet mesh = new PgElementSet(3);

		mesh.addVertex(new PdVector(0, 0, 0)); // p1
		mesh.addVertex(new PdVector(1, 0, 0)); // p2
		mesh.addVertex(new PdVector(0, 1, 0)); // p3

		mesh.addElement(new PiVector(0, 1, 2)); // a triangle

		// Manually calculate the M matrix
		PdMatrix M_true = new PdMatrix(3, 3);
		M_true.setEntry(0, 0, 0.5 / 3);
		M_true.setEntry(0, 1, 0);
		M_true.setEntry(0, 2, 0);
		M_true.setEntry(1, 0, 0);
		M_true.setEntry(1, 1,  0.5 / 3);
		M_true.setEntry(1, 2, 0);
		M_true.setEntry(2, 0, 0);
		M_true.setEntry(2, 1, 0);
		M_true.setEntry(2, 2,  0.5 / 3);
	

		// System.out.println("M_true matrix: " + M_true.toString());

		// Calculate the gradient matrix using the function
		PdMatrix M_calculated = calculateMassMatrix(mesh);
		// System.out.println("M_calculated matrix: " + M_calculated.toString());

		// Compare the two matrices
		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < 3; j++) {
				if (Math.abs(M_true.getEntry(i, j) - M_calculated.getEntry(i, j)) > 0.0001) {
					m_MMatrixResult.setText("Failed");
					return;
				}
			}
		}
		m_MMatrixResult.setText("Success");
	}

	// task 1.4 Compute the cotangent matrix S
	private PdMatrix calculateSMatrix(PgElementSet mesh) {
		m_SMatrixResult.setText("Calculating...");
		PdMatrix M = calculateMassMatrix(mesh);
		PdMatrix G = calculateGMatrix(mesh);
		PdMatrix G_copy = new PdMatrix(G.getNumRows(), G.getNumCols());

		// calculate the transpose of G
		for (int i = 0; i < G.getNumRows(); i++) {
			for (int j = 0; j < G.getNumCols(); j++) {
				G_copy.setEntry(j, i, G.getEntry(i, j));
			}
		}

		G_copy.rightMult(M);
		G_copy.rightMult(G);

		m_SMatrixResult.setText("Calculated!");
		return G_copy;
	}

	// task 1.4 Compute the cotangent matrix S (Test S using a generated triangle mesh)
	private void testSMatrix() {
		// Create a triangle mesh
		PgElementSet mesh = new PgElementSet(3);

		mesh.addVertex(new PdVector(0, 0, 0)); // p1
		mesh.addVertex(new PdVector(1, 0, 0)); // p2
		mesh.addVertex(new PdVector(0, 1, 0)); // p3

		mesh.addElement(new PiVector(0, 1, 2)); // a triangle

		// Manually calculate the S matrix
		PdMatrix S_true = new PdMatrix(3, 3);
		S_true.setEntry(0, 0, 1.0/3);
		S_true.setEntry(0, 1, -0.5/3);
		S_true.setEntry(0, 2, -0.5/3);
		S_true.setEntry(1, 0, -0.5/3);
		S_true.setEntry(1, 1, 0.5/3);
		S_true.setEntry(1, 2, 0);
		S_true.setEntry(2, 0, -0.5/3);
		S_true.setEntry(2, 1, 0);
		S_true.setEntry(2, 2, 0.5/3);


		// System.out.println("S_true matrix: " + S_true.toString());

		// Calculate the gradient matrix using the function
		PdMatrix S_calculated = calculateSMatrix(mesh);
		// System.out.println("S_calculated matrix: " + S_calculated.toString());

		// Compare the two matrices
		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < 3; j++) {
				if (Math.abs(S_true.getEntry(i, j) - S_calculated.getEntry(i, j)) > 0.0001) {
					m_SMatrixResult.setText("Failed");
					return;
				}
			}
		}
		m_SMatrixResult.setText("Success");
	}
	

	/**
	 * Get information which bottom buttons a dialog should create
	 * when showing this info panel.
	 */
	protected int getDialogButtons()		{
		return PsDialog.BUTTON_OK;
	}
}
