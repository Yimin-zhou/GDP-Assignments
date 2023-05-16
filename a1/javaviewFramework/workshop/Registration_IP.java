package workshop;

import java.awt.BorderLayout;
import java.awt.Button;
import java.awt.Dimension;
import java.awt.GridLayout;
import java.awt.Label;
import java.awt.List;
import java.awt.Panel;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Vector;

import Jama.Matrix;
import Jama.SingularValueDecomposition;
import jv.geom.PgElementSet;
import jv.object.PsConfig;
import jv.object.PsDialog;
import jv.object.PsUpdateIf;
import jv.objectGui.PsList;
import jv.project.PgGeometryIf;
import jv.project.PvGeometryIf;
import jv.viewer.PvDisplay;
import jvx.project.PjWorkshop_IP;

import jv.vecmath.PdVector;
import jv.vecmath.PdMatrix;
import java.util.Random;



/**
 * Info Panel of Workshop for surface registration
 *
 */
public class Registration_IP extends PjWorkshop_IP implements ActionListener{
	protected	List			m_listActive;
	protected	List			m_listPassive;
	protected	Vector			m_geomList;
	protected	Registration	m_registration;
	protected   Button			m_bSetSurfaces;

	protected Button m_bStartICP;
	// task 2.1 step 1
	private static final int NUM_ICP_POINTS = 50;
	private Random random = new Random();


	/** Constructor */
	public Registration_IP () {
		super();
		if (getClass() == Registration_IP.class)
			init();
	}

	/**
	 * Informational text on the usage of the dialog.
	 * This notice will be displayed if this info panel is shown in a dialog.
	 * The text is split at line breaks into individual lines on the dialog.
	 */
	public String getNotice() {
		return "This text should explain what the workshop is about and how to use it.";
	}

	/** Assign a parent object. */
	public void setParent(PsUpdateIf parent) {
		super.setParent(parent);
		m_registration = (Registration)parent;

		addSubTitle("Select Surfaces to be Registered");

		Panel pGeometries = new Panel();
		pGeometries.setLayout(new GridLayout(1, 2));

		Panel Passive = new Panel();
		Passive.setLayout(new BorderLayout());
		Panel Active = new Panel();
		Active.setLayout(new BorderLayout());
		Label ActiveLabel = new Label("Surface P");
		Active.add(ActiveLabel, BorderLayout.NORTH);
		m_listActive = new PsList(5, true);
		Active.add(m_listActive, BorderLayout.CENTER);
		pGeometries.add(Active);
		Label PassiveLabel = new Label("Surface Q");
		Passive.add(PassiveLabel, BorderLayout.NORTH);
		m_listPassive = new PsList(5, true);
		Passive.add(m_listPassive, BorderLayout.CENTER);
		pGeometries.add(Passive);
		add(pGeometries);

		Panel pSetSurfaces = new Panel(new BorderLayout());
		m_bSetSurfaces = new Button("Set selected surfaces");
		m_bSetSurfaces.addActionListener(this);
		pSetSurfaces.add(m_bSetSurfaces, BorderLayout.CENTER);
		add(pSetSurfaces);

		// task 2
		m_bStartICP = new Button("Start ICP");
		m_bStartICP.addActionListener(this);
		Panel pStartICP = new Panel(new BorderLayout());
		pStartICP.add(m_bStartICP, BorderLayout.CENTER);
		add(pStartICP);

		updateGeomList();
		validate();
	}

	/** Initialisation */
	public void init() {
		super.init();
		setTitle("Surface Registration");

	}

	/** Set the list of geometries in the lists to the current state of the display. */
	public void updateGeomList() {
		Vector displays = m_registration.getGeometry().getDisplayList();
		int numDisplays = displays.size();
		m_geomList = new Vector();
		for (int i=0; i<numDisplays; i++) {
			PvDisplay disp =((PvDisplay)displays.elementAt(i));
			PgGeometryIf[] geomList = disp.getGeometries();
			int numGeom = geomList.length;
			for (int j=0; j<numGeom; j++) {
				if (!m_geomList.contains(geomList[j])) {
					//Take just PgElementSets from the list.
					if (geomList[j].getType() == PvGeometryIf.GEOM_ELEMENT_SET)
						m_geomList.addElement(geomList[j]);
				}
			}
		}
		int nog = m_geomList.size();
		m_listActive.removeAll();
		m_listPassive.removeAll();
		for (int i=0; i<nog; i++) {
			String name = ((PgGeometryIf)m_geomList.elementAt(i)).getName();
			m_listPassive.add(name);
			m_listActive.add(name);
		}
	}
	/**
	 * Handle action events fired by buttons etc.
	 */
	public void actionPerformed(ActionEvent event) {
		Object source = event.getSource();
		if (source == m_bSetSurfaces) {
			m_registration.setGeometries((PgElementSet)m_geomList.elementAt(m_listActive.getSelectedIndex()),
					(PgElementSet)m_geomList.elementAt(m_listPassive.getSelectedIndex()));
			return;
		}
		else if (source == m_bStartICP) {
			performICP();
			return;
		}
	}

	// task 2
	public void performICP() {
		PgElementSet p = (PgElementSet) m_geomList.elementAt(m_listActive.getSelectedIndex());
		PgElementSet q = (PgElementSet) m_geomList.elementAt(m_listPassive.getSelectedIndex());

		int[] selectedIndices = selectRandomVertices(p, NUM_ICP_POINTS);
		ArrayList<PdVector> pVertices = new ArrayList<>();
		ArrayList<PdVector> qVertices = new ArrayList<>();
		ArrayList<Double> distances = new ArrayList<>();

		for (int i : selectedIndices) {
			PdVector pi = p.getVertex(i);
			PdVector qi = findClosestVertex(pi, q);

			// Store vertices and distance
			pVertices.add(pi);
			qVertices.add(qi);
			distances.add(pi.dist(qi));
		}

		// Compute the median distance
		Collections.sort(distances);
		double medianDistance;
		if (distances.size() % 2 == 0) {
			medianDistance = (distances.get(distances.size()/2 - 1) + distances.get(distances.size()/2)) / 2.0;
		} else {
			medianDistance = distances.get(distances.size()/2);
		}

		double k = 1.5;
		// Determine the threshold for removing pairs
		double threshold = k * medianDistance;

		// Create new ArrayLists for filtered vertices
		ArrayList<PdVector> newPVertices = new ArrayList<>();
		ArrayList<PdVector> newQVertices = new ArrayList<>();

		// Filter pairs with a distance larger than the threshold
		for (int i = 0; i < distances.size(); i++) {
			if (distances.get(i) <= threshold) {
				newPVertices.add(pVertices.get(i));
				newQVertices.add(qVertices.get(i));
			}
		}

		// Replace the original vertex lists with the filtered ones
		pVertices = newPVertices;
		qVertices = newQVertices;

		// task 2.3
		// Compute centroids of P and Q points
		PdVector centroidP = new PdVector(3);
		PdVector centroidQ = new PdVector(3);
		for (int i = 0; i < pVertices.size(); i++) {
			centroidP.add(pVertices.get(i));
			centroidQ.add(qVertices.get(i));
		}
		centroidP.multScalar(1.0 / pVertices.size());
		centroidQ.multScalar(1.0 / qVertices.size());

		// Build matrix H
		PdMatrix M = new PdMatrix(3, 3);
		for (int i = 0; i < pVertices.size(); i++) {
			PdVector pCentroid = PdVector.subNew(pVertices.get(i), centroidP);
			PdVector qCentroid = PdVector.subNew(qVertices.get(i), centroidQ);
			PdMatrix pqT = outerProductNew(pCentroid, qCentroid);
			M.add(pqT);
		}

		// Compute SVD of H
		Matrix jamaM = new Matrix(M.getEntries());
		SingularValueDecomposition SVD = jamaM.svd();
		PdMatrix U = new PdMatrix(SVD.getU().getArray());
		PdMatrix V = new PdMatrix(SVD.getV().getArray());
		PdMatrix s = new PdMatrix(SVD.getS().getArray());
		// M.svd(U, s, V); // TODO calculate SVD

		// Compute rotation matrix R
		PdMatrix R = new PdMatrix(3, 3);
		R.leftMult(U);
		R.rightMult(V);

		if (R.det() < 0) {
			PdMatrix diag = new PdMatrix(new double[][]{{1, 0, 0}, {0, 1, 0}, {0, 0, -1}});
			U.leftMult(diag);
			R = U;
			R.rightMult(V);
		}

		// Compute translation vector t
		PdVector t = PdVector.subNew(centroidQ, R.leftMultMatrix(centroidP, centroidP)); // TODO calculate t

		// Apply the optimal rigid transformation to P
		for (PdVector vertex : pVertices) {
			vertex.leftMultMatrix(R);
			vertex.add(t);
		}

		// Call this to update the geometry of P
		p.update(p);
	}



	private int[] selectRandomVertices(PgElementSet mesh, int numVertices) {
		int[] indices = new int[numVertices];
		for (int i = 0; i < numVertices; i++) {
			indices[i] = random.nextInt(mesh.getNumVertices());
		}
		return indices;
	}

	private PdVector findClosestVertex(PdVector point, PgElementSet mesh) {
		PdVector closest = null;
		double closestDistance = Double.MAX_VALUE;

		// brute force search
		for (int i = 0; i < mesh.getNumVertices(); i++) {
			PdVector vertex = mesh.getVertex(i);
			double distance = point.dist(vertex);
			if (distance < closestDistance) {
				closest = vertex;
				closestDistance = distance;
			}
		}

		return closest;
	}
	public static PdMatrix outerProductNew(PdVector a, PdVector b) {
		PdMatrix result = new PdMatrix(3, 3);
		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < 3; j++) {
				result.setEntry(i, j, a.getEntry(i) * b.getEntry(j));
			}
		}
		return result;
	}

	/**
	 * Get information which bottom buttons a dialog should create
	 * when showing this info panel.
	 */
	protected int getDialogButtons()		{
		return PsDialog.BUTTON_OK;
	}
}