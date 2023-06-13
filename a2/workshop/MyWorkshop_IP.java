package workshop;

import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;
import java.awt.geom.Area;
import java.util.*;
import java.util.List;

import javax.swing.text.html.parser.Element;

//import org.omg.Messaging.SYNC_WITH_TRANSPORT;

import dev6.numeric.PnMumpsSolver;
import jv.geom.PgBndPolygon;
import jv.geom.PgEdgeStar;
import jv.number.PdVector_IP;
import jv.number.PuDouble;
import jv.object.PsConfig;
import jv.object.PsDialog;
import jv.object.PsObject;
import jv.object.PsUpdateIf;
import jv.vecmath.PiVector;
import jvx.numeric.PnSparseMatrix;
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
	protected PnSparseMatrix m_GMatrix;

	// task 1.2 Compute combinatorial Laplacian matrix L
	protected Button m_bCalculateLMatrix;
	protected TextField m_LMatrixResult;
	protected PnSparseMatrix m_LMatrix; // Combinatorial Laplacian matrix

	// task 1.3 Compute Mass Matrix M
	protected Button m_bCalculateMMatrix;
	protected TextField m_MMatrixResult;
	protected PnSparseMatrix m_MMatrix; // Mass matrix

	// task 1.4 Compute Cotangent matrix S
	protected Button m_bCalculateSMatrix;
	protected TextField m_SMatrixResult;
	protected PnSparseMatrix m_SMatrix; // Cotangent matrix

	// task 2 Input Matrix A
	protected PdVector_IP m_vRow1;
	protected PdVector_IP m_vRow2;
	protected PdVector_IP m_vRow3;
	protected PdMatrix m_inputMatrix = new PdMatrix(3, 3);
	protected Button m_bComputeSurface;
	protected Button m_bResetSurface;
	
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


		// task 2 Vector fields to specify 3x3 matrix
		m_vRow1 = new PdVector_IP();
		m_vRow1.setLayout(new FlowLayout(FlowLayout.CENTER));
		m_vRow1.setVector(new PdVector(2, 0, 0));
		m_vRow1.setParent(this);
		add(m_vRow1);

		m_vRow2 = new PdVector_IP();
		m_vRow2.setLayout(new FlowLayout(FlowLayout.CENTER));
		m_vRow2.setVector(new PdVector(0, 2, 0));
		m_vRow2.setParent(this);
		add(m_vRow2);

		m_vRow3 = new PdVector_IP();
		m_vRow3.setLayout(new FlowLayout(FlowLayout.CENTER));
		m_vRow3.setVector(new PdVector(0, 0, 2));
		m_vRow3.setParent(this);
		add(m_vRow3);
		m_inputMatrix.setRow(0, m_vRow1.getVector());
		m_inputMatrix.setRow(1, m_vRow2.getVector());
		m_inputMatrix.setRow(2, m_vRow3.getVector());

		// task 2 compute surface deformation
		Panel computePanel = new Panel(new BorderLayout());
		m_bComputeSurface = new Button("Compute deformed surface");
		m_bComputeSurface.addActionListener(this);
		computePanel.add(m_bComputeSurface, BorderLayout.CENTER);
		add(computePanel);

		// task 2 reset mesh to original
		Panel resetPanel = new Panel(new BorderLayout());
		m_bResetSurface = new Button("Reset");
		m_bResetSurface.addActionListener(this);
		resetPanel.add(m_bResetSurface, BorderLayout.CENTER);
		add(resetPanel);
		
		validate();
	}
	
	
	public boolean update(Object event) {
		if (event == m_xOff) {
			m_ws.setXOff(m_xOff.getValue());
			m_ws.m_geom.update(m_ws.m_geom);
			return true;
		} else if (event == m_vRow1) {
			System.out.println("row1 " + m_vRow1.getVector().toShortString());
			m_inputMatrix.setRow(0, m_vRow1.getVector());
			return true;
		} else if (event == m_vRow2) {
			System.out.println("row2 " + m_vRow2.getVector().toShortString());
			m_inputMatrix.setRow(1, m_vRow2.getVector());
			return true;
		} else if (event == m_vRow3) {
			System.out.println("row3 " + m_vRow3.getVector().toShortString());
			m_inputMatrix.setRow(2, m_vRow3.getVector());
			return true;
		} else return super.update(event);
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
			m_GMatrix = calculateGmatrix(m_ws.m_geom);
            testGmatrix();
            return;
        }

		// task 1.2 calculate L
		else if (source == m_bCalculateLMatrix) {
			m_LMatrix = calculateLmatrix(m_ws.m_geom);
			testLmatrix();
			return;
		}

		// task 1.3 calculate M
		else if (source == m_bCalculateMMatrix) {
			m_MMatrix = calculateMVmatrix(m_ws.m_geom);
			testMVmatrix();
			return;
		}

		// task 1.4 calculate S
		else if (source == m_bCalculateSMatrix) {
			m_SMatrix = calculateSmatrix(m_ws.m_geom);
			testSmatrix();
			return;
		}

		// task 2
		else if (source == m_bComputeSurface) {
			System.out.println("Compute deformed surface");
			computeSurface();
		}

		else if (source == m_bResetSurface) {
			System.out.println("Reset");
			m_ws.m_geom.setVertices(m_ws.m_geomSave.getVertices().clone());
			m_ws.m_geom.update(m_ws.m_geom);
		}

		else
			super.actionPerformed(event);
	}

	// task 1.1 G matrix
	private PnSparseMatrix calculateGmatrix(PgElementSet mesh) {
		m_GMatrixResult.setText("Calculating...");

		int nVertices = mesh.getNumVertices();
		int nTriangles = mesh.getNumElements();

		PnSparseMatrix G = new PnSparseMatrix(3 * nTriangles, nVertices);

		for (int i = 0; i < nTriangles; i++) {
			PiVector triangle = mesh.getElement(i);

			List<Integer> vertexIndices = new ArrayList<Integer>();
			PdVector[] vertices = new PdVector[3];
			for (int j = 0; j < 3; j++) {
				vertices[j] = mesh.getVertex(triangle.getEntry(j));
				vertexIndices.add(triangle.getEntry(j));
			}
			
			PdMatrix gradient = calculateTriangleGradient(vertices, vertexIndices);

			// Add the gradient to the G matrix (vertical stacking), gradient is 3x3 matrix
			for (int j = 0; j < 3; j++) {
				int ind = vertexIndices.get(j);
				for (int k = 0; k < 3; k++) {
					G.setEntry(i * 3 + k, ind, gradient.getEntry(k, j));
				}
			}
		}
		m_GMatrixResult.setText("Calculated!");

		return G;
	}

	// task 1.1 G matrix (Gradient of a function on a triangle)
	private PdMatrix calculateTriangleGradient(PdVector[] vertices, List<Integer> vertexIndices) {
		double area = calculateTriangleArea(vertices);
		PdMatrix gradientMatrix = new PdMatrix(3, 3);

		// Calculate the triangle normal, the normal is consistent with the order of the vertices
		PdVector triangleNormal = PdVector.crossNew(PdVector.subNew(vertices[1], vertices[0]), PdVector.subNew(vertices[2], vertices[0]));
		triangleNormal.normalize();

		// Calculate the gradient of the function on the triangle
		for (int i = 0; i < 3; i++) {
			// Edge vector, following the order of the vertices
			PdVector edgeVector = PdVector.subNew(vertices[(i + 2) % 3], vertices[(i + 1) % 3]);

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
	private void testGmatrix() {
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

		// Calculate the gradient matrix using the function
		PnSparseMatrix gradientMatrix_calculated = calculateGmatrix(mesh);
		// System.out.println("Gradient_calculated matrix: " + gradientMatrix_calculated.toString());

		PnSparseMatrix gradientMatrix_true_ = new PnSparseMatrix(gradientMatrix_true);
		// Compare the two gradient matrices
		for (int i = 0; i < 6; i++) {
			for (int j = 0; j < 4; j++) {
				if (Math.abs(gradientMatrix_true_.getEntry(i, j) - gradientMatrix_calculated.getEntry(i, j)) > 0.0001) {
					m_GMatrixResult.setText("Failed");
					return;
				}
			}
		}
		m_GMatrixResult.setText("Success");
	}

	// task 1.2 L matrix (Calculate the Combinatorial Laplacian matrix)
	private PnSparseMatrix calculateLmatrix(PgElementSet mesh) {
		m_LMatrixResult.setText("Calculating...");
		int nVertices = mesh.getNumVertices();
		PnSparseMatrix L = new PnSparseMatrix(nVertices, nVertices);
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
	private void testLmatrix() {
		// Create a triangle mesh
		PgElementSet mesh = new PgElementSet(3);

		mesh.addVertex(new PdVector(0, 0, 0)); // p1
		mesh.addVertex(new PdVector(1, 0, 0)); // p2
		mesh.addVertex(new PdVector(0, 1, 0)); // p3

		mesh.addElement(new PiVector(0, 1, 2)); // a triangle

		// Manually calculate the L matrix
		PnSparseMatrix L_true = new PnSparseMatrix(3, 3);
		L_true.setEntry(0, 0, 2);
		L_true.setEntry(0, 1, -1);
		L_true.setEntry(0, 2, -1);
		L_true.setEntry(1, 0, -1);
		L_true.setEntry(1, 1, 2);
		L_true.setEntry(1, 2, -1);
		L_true.setEntry(2, 0, -1);
		L_true.setEntry(2, 1, -1);
		L_true.setEntry(2, 2, 2);
	
		PnSparseMatrix L_calculated = calculateLmatrix(mesh);

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
	private PnSparseMatrix calculateMVmatrix(PgElementSet mesh) {
		m_MMatrixResult.setText("Calculating...");

		int nTriangles = mesh.getNumElements();
		PnSparseMatrix M = new PnSparseMatrix(3*nTriangles, 3*nTriangles);

		// Populate areas list
		for (int i = 0; i < mesh.getNumElements(); i++) {
			PiVector triangle = mesh.getElement(i);
			PdVector p1 = mesh.getVertex(triangle.getEntry(0));
			PdVector p2 = mesh.getVertex(triangle.getEntry(1));
			PdVector p3 = mesh.getVertex(triangle.getEntry(2));
			PdVector[] vertices = {p1, p2, p3};
			double area = calculateTriangleArea(vertices);
			int index = 3*i;
			M.setEntry(index, index, area);
			M.setEntry(index+1, index+1, area);
			M.setEntry(index+2, index+2, area);
		}

		m_MMatrixResult.setText("Calculated!");
		return M;
	}

	// task 1.3 Compute Mass Matrix M (Test M using a generated triangle mesh)
	private void testMVmatrix() {
		// Create a triangle mesh
		PgElementSet mesh = new PgElementSet(3);

		mesh.addVertex(new PdVector(0, 0, 0)); // p1
		mesh.addVertex(new PdVector(1, 0, 0)); // p2
		mesh.addVertex(new PdVector(0, 1, 0)); // p3

		mesh.addElement(new PiVector(0, 1, 2)); // a triangle

		// Manually calculate the M matrix
		PnSparseMatrix M_true = new PnSparseMatrix(3, 3);
		M_true.setEntry(0, 0, 0.5);
		M_true.setEntry(0, 1, 0);
		M_true.setEntry(0, 2, 0);
		M_true.setEntry(1, 0, 0);
		M_true.setEntry(1, 1,  0.5);
		M_true.setEntry(1, 2, 0);
		M_true.setEntry(2, 0, 0);
		M_true.setEntry(2, 1, 0);
		M_true.setEntry(2, 2,  0.5);
	

		PnSparseMatrix M_calculated = calculateMVmatrix(mesh);

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
	private PnSparseMatrix calculateSmatrix(PgElementSet mesh) {
		m_SMatrixResult.setText("Calculating...");
		m_GMatrix = calculateGmatrix(mesh);
		PnSparseMatrix M = calculateMVmatrix(mesh);
		PnSparseMatrix G_T = m_GMatrix.transposeNew();

		PnSparseMatrix temp = PnSparseMatrix.multMatrices(G_T, M, null);
		PnSparseMatrix result = PnSparseMatrix.multMatrices(temp, m_GMatrix, null);

		return result;
	}

	// task 1.4 Compute the cotangent matrix S (Test S using a generated triangle mesh)
	private void testSmatrix() {
		PgElementSet mesh = new PgElementSet(3);

		mesh.addVertex(new PdVector(0, 0, 0)); // p1
		mesh.addVertex(new PdVector(1, 0, 0)); // p2
		mesh.addVertex(new PdVector(0, 1, 0)); // p3

		mesh.addElement(new PiVector(0, 1, 2)); // a triangle

		// Manually calculate the S matrix
		PnSparseMatrix S_true = new PnSparseMatrix(3, 3);
		S_true.setEntry(0, 0, 1.0);
		S_true.setEntry(0, 1, -0.5);
		S_true.setEntry(0, 2, -0.5);
		S_true.setEntry(1, 0, -0.5);
		S_true.setEntry(1, 1, 0.5);
		S_true.setEntry(1, 2, 0);
		S_true.setEntry(2, 0, -0.5);
		S_true.setEntry(2, 1, 0);
		S_true.setEntry(2, 2, 0.5);

		PnSparseMatrix S_calculated = calculateSmatrix(mesh);

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

	// task 2 Shape editing using differential coordinates
	private void computeSurface() {
		PgElementSet mesh = m_ws.m_geom; // get mesh
		m_GMatrix = calculateGmatrix(mesh); // init G matrix
		PdMatrix grads = getGradients(); // get gradients with input matrix applied to selected vertices

		PnSparseMatrix G_T = m_GMatrix.transposeNew(); // G transpose
		PnSparseMatrix M = calculateMVmatrix(mesh); // get Mv matrix
		PnSparseMatrix G_TM = PnSparseMatrix.multMatrices(G_T, M, null); // get right matrix
		PnSparseMatrix S = calculateSmatrix(mesh); // get left matrix

		// Init output vertices and get the right gradient vectors
		PdVector[] g = new PdVector[3];
		PdVector[] v = new PdVector[3];
		for (int i = 0; i < 3; i++) {
			g[i] = PnSparseMatrix.rightMultVector(G_TM, grads.getRow(i), null);
			v[i] = new PdVector(mesh.getNumVertices());
		}

		PdVector[] vertices = new PdVector[mesh.getNumVertices()];
		PdMatrix verticesM = new PdMatrix(mesh.getNumVertices(), 3);

		// solve the linear systems using factorization to get vectors v
		try {
			S.validate();
			long factorization = dev6.numeric.PnMumpsSolver.factor(S, PnMumpsSolver.Type.GENERAL_SYMMETRIC);
			for (int i = 0; i < 3; i++) {
				dev6.numeric.PnMumpsSolver.solve(factorization, v[i], g[i]);
				verticesM.setColumn(i, v[i]);
			}

		} catch (Exception e) {
			e.printStackTrace();
			System.out.println("solving did not succeed");
		}

		PdVector old_centroid = mesh.getCenterOfGravity(); // save old centroid

		// update mesh
		for (int i = 0; i < mesh.getNumVertices(); i++) {
			vertices[i] = verticesM.getRow(i);
		}

		mesh.setVertices(vertices);
		mesh.update(mesh);

		PdVector centroid_diff = PdVector.subNew(mesh.getCenterOfGravity(), old_centroid); // get the difference between new and old centroid

		// subtract it from all vertices to keep the barycentre constant
		for (int i = 0; i < mesh.getNumVertices(); i++) {
			vertices[i].sub(centroid_diff);
		}

		// update mesh with correct vertices
		mesh.setVertices(vertices);
		mesh.update(mesh);
	}

	private PdMatrix getGradients() {
		PgElementSet mesh = m_ws.m_geom;
		PdMatrix vertices = new PdMatrix(3, mesh.getNumVertices());
		PdMatrix grads = new PdMatrix(3, m_GMatrix.getNumRows());

		// save vertices to matrix
		for (int i = 0; i < mesh.getNumVertices(); i++) {
			vertices.setColumn(i, mesh.getVertex(i));
		}

		// multiply each coordinate vector of vertex matrix by G matrix to get gradient vectors
		for (int i = 0; i < 3; i++) {
			grads.setRow(i, PnSparseMatrix.rightMultVector(m_GMatrix, vertices.getRow(i),null));
		}

		// for each selected triangle, multiply the corresponding 3x3 gradient by the 3x3 input matrix
		PdMatrix elementM = new PdMatrix(3, 3);
		for (int i = 0; i < mesh.getNumElements(); i++) {
			if (mesh.getElement(i).hasTag(PsObject.IS_SELECTED)) {
				for (int j = 0; j < 3; j++) {
					elementM.setColumn(j, grads.getColumn(i*3+j));
				}
				elementM.leftMult(m_inputMatrix);
				for (int j = 0; j < 3; j++) {
					grads.setColumn(i*3+j, elementM.getColumn(j));
				}
			}
		}

		return grads;
	}

	/**
	 * Get information which bottom buttons a dialog should create
	 * when showing this info panel.
	 */
	protected int getDialogButtons()		{
		return PsDialog.BUTTON_OK;
	}
}
