package workshop;

import java.awt.Button;
import java.awt.FlowLayout;
import java.awt.Panel;
import java.awt.TextField;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;
import java.util.*;

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
		// task 1.1 calculate genus	
		else if (source == m_bCalculateGMatrix) {
            testGMatrix();
            return;
        }
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

			// Add the gradient to the G matrix
			for (int j = 0; j < 3; j++) {
				for (int k = 0; k < nVertices; k++) {
					G.setEntry(i * 3 + j, k, gradient.getEntry(j, k));
				}
			}
		}
		// print G matrix
		// System.out.println("G matrix: " + G.toString());

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
		// Create a triangle mesh
		PgElementSet mesh = new PgElementSet(3);

		mesh.addVertex(new PdVector(0, 0, 0)); // p1
		mesh.addVertex(new PdVector(1, 0, 0)); // p2
		mesh.addVertex(new PdVector(0, 1, 0)); // p3

		mesh.addElement(new PiVector(0, 1, 2));

		// Manually calculate the gradient matrix of the function on the triangle
		PdVector normal = new PdVector(0, 0, 1);
		double area = 0.5;

		PdVector e1 = new PdVector(-1,1,0); // 3-2
		PdVector e2 = new PdVector(0,-1,0); // 1-3
		PdVector e3 = new PdVector(1,0,0); // 2-1

		PdMatrix gradientMatrix_true = new PdMatrix(3, 3);
		gradientMatrix_true.setColumn(0, PdVector.crossNew(normal, e1));
		gradientMatrix_true.setColumn(1, PdVector.crossNew(normal, e2));
		gradientMatrix_true.setColumn(2, PdVector.crossNew(normal, e3));

		// print gradient matrix
		System.out.println("Gradient_true matrix: " + gradientMatrix_true.toString());

		// Calculate the gradient matrix using the function
		PdMatrix gradientMatrix_calculated = calculateGMatrix(mesh);
		System.out.println("Gradient_calculated matrix: " + gradientMatrix_calculated.toString());

		if (gradientMatrix_true.toString().equals(gradientMatrix_calculated.toString())) {
			m_GMatrixResult.setText("Success");
		} else {
			m_GMatrixResult.setText("Fail");
		}
	}
	

	/**
	 * Get information which bottom buttons a dialog should create
	 * when showing this info panel.
	 */
	protected int getDialogButtons()		{
		return PsDialog.BUTTON_OK;
	}
}
